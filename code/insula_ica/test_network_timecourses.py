"""Checks for the replacement's data direction, confounds, contrasts and GLS."""
import os
import json
from types import SimpleNamespace
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import nibabel as nib
import numpy as np

import run_network_timecourse_analysis as n
import remove_withdrawn_ica_projections as withdrawal


class NetworkTimecourseTests(unittest.TestCase):
    def test_cleanup_removes_only_withdrawn_derivatives_and_preserves_source_maps(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d); record = root/'code/insula_ica/records/linux2-completed-ICA-exploration-test'
            record.mkdir(parents=True)
            output = root/'derivatives/insula_ica/completed-ICA-exploration-test'; output.mkdir(parents=True)
            (record/'state.json').write_text(json.dumps(dict(output_root=str(output))))
            for parent in (record, output):
                for name in withdrawal.NAMES:
                    path = parent/name; path.parent.mkdir(exist_ok=True); path.write_text('withdrawn')
            source = root/'derivatives/fsl/source.feat/stats/cope10.nii.gz'; source.parent.mkdir(parents=True); source.write_bytes(b'preserve')
            marker = output/'runs/run1/complete.json'; marker.parent.mkdir(parents=True)
            marker.write_text(json.dumps(dict(result=dict(effects=[1, 2], diagnostics=['keep']), signature=dict(cope='old', model='keep'))))
            retained = output/'runs/run1/joint_bold_bilateral_connectivity.nii.gz'; retained.write_bytes(b'keep map')
            with patch.object(withdrawal.core, 'ROOT', root), patch.object(withdrawal.subprocess, 'check_output', return_value=''):
                result = withdrawal.clean()
            self.assertEqual(len(result['removed']), 12)
            self.assertEqual(source.read_bytes(), b'preserve'); self.assertEqual(retained.read_bytes(), b'keep map')
            cleaned = json.loads(marker.read_text()); self.assertNotIn('effects', cleaned['result'])
            self.assertEqual(cleaned['result']['diagnostics'], ['keep']); self.assertEqual(cleaned['signature'], dict(model='keep'))

    def test_stage1_preserves_task_signal_and_joint_model_handles_shared_nuisance(self):
        rng = np.random.default_rng(404); t = 180
        maps = rng.normal(size=(100, 5)); w = n.exploration.weights(maps)
        task = np.sin(np.arange(t)/10); nuisance = .7*task+rng.normal(size=t)
        x = np.column_stack([task, nuisance]); x -= x.mean(0)
        beta = rng.normal(size=(2, 5)); tc = x@beta
        y = maps@tc.T+5000+rng.normal(size=(100, 1))
        recovered = n.task_timecourses(w, y)
        np.testing.assert_allclose(recovered, tc, atol=1e-10)
        np.testing.assert_allclose(np.linalg.pinv(x)@recovered, beta, atol=1e-10)
        # Task-cleaned residuals would have destroyed precisely this response.
        residual = recovered-x@np.linalg.pinv(x)@recovered
        self.assertLess(np.linalg.norm(residual), 1e-9)

    def test_condition_means_are_fitted_contrasts_not_fabricated_from_difference(self):
        rng = np.random.default_rng(1); design = rng.normal(size=(150, 13)); settings = {'fmri(evs_real)': '10'}
        for i, title in enumerate(('C_def', 'C_rec', 'F_def', 'F_rec', 'S_def', 'S_rec'), 4): settings[f'fmri(evtitle{i})'] = title
        saved = np.zeros((18, 13)); saved[9, :9] = [0, 0, 0, -1, 1, -1, 1, -1, 1]
        cs, labels = n.contrasts_for(dict(task='trust', contrast_index=10), settings, design, saved)
        np.testing.assert_allclose(cs[0], 3*(cs[1]-cs[2]))
        self.assertEqual(labels[1], 'betrayal_mean_CFS')
        settings['fmri(evtitle4)'] = 'wrong'
        with self.assertRaisesRegex(ValueError, 'EV order'): n.contrasts_for(dict(task='trust', contrast_index=10), settings, design, saved)

    def test_group_uses_participants_and_all_primary_tests(self):
        rng = np.random.default_rng(20); rows = []
        subjects = [str(k) for k in range(132)]
        for mode in n.g.inputs.MODES:
            for roi in n.g.core.ROIS:
                for task, contrast in [('trust', 'betrayal_minus_reciprocity_sum'), ('ugr', 'offer_pmod_as_stored')]:
                    for k in range(1, 6):
                        for j, subject in enumerate(subjects):
                            value = rng.normal()
                            for run, delta in ([(1, 0)] if j%2 else [(1, -2), (2, 2)]):
                                rows.append(dict(mode=mode, roi=roi, task=task, contrast=contrast, component=k,
                                                 subject=subject, run=run, effect=value+delta))
        participant, group = n.summarize(rows, subjects)
        self.assertEqual(len(group), 60); self.assertEqual(len(participant), 60*132)
        for r in group: self.assertAlmostEqual(r['p_bonferroni_60'], min(1, 60*r['p']))
        with self.assertRaisesRegex(ValueError, 'Missing participant'): n.summarize(rows[:-2], subjects)

    def test_real_film_estimates_each_network_separately_without_spatial_smoothing(self):
        film = Path(os.environ.get('FSLDIR', '/Users/tug87422/fsl'))/'bin/film_gls'
        if not film.is_file(): self.skipTest('FSL is not installed')
        rng = np.random.default_rng(5); t = 170
        x = rng.normal(size=(t, 4)); x -= x.mean(0)
        noise = rng.normal(size=(t, 5))
        for i in range(1, t): noise[i] += .55*noise[i-1]
        tc = x@rng.normal(size=(4, 5))+noise; tc -= tc.mean(0)
        contrasts = np.eye(4)[:2]
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            effect, variance, ols = n.fit_timecourses(tc, x, contrasts, root/'all', 1.710, film)
            self.assertEqual(effect.shape, (2, 5)); self.assertTrue((variance > 0).all())
            np.testing.assert_allclose(ols, contrasts@np.linalg.pinv(x)@tc)
            # Independent single-network reference: no cross-network AC pooling.
            from subprocess import run
            folder = root/'single'; folder.mkdir()
            image = nib.Nifti1Image(tc[:, 2].reshape(1, 1, 1, t).astype('float32'), np.eye(4)); nib.save(image, folder/'tc.nii.gz')
            result = run([str(film), '--in='+str(folder/'tc.nii.gz'), '--rn='+str(folder/'stats'),
                          '--pd='+str(root/'all/design.mat'), '--con='+str(root/'all/design.con'), '--thr=-1000000000'], capture_output=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            for k in range(2): self.assertAlmostEqual(effect[k, 2], nib.load(folder/f'stats/cope{k+1}.nii.gz').get_fdata().item(), places=5)
            command = (root/'all/command.json').read_text()
            self.assertNotIn('--sa', command); self.assertNotIn('--noest', command)

    def test_three_run_pilot_reads_4d_and_resumes_without_any_source_task_map(self):
        film = Path(os.environ.get('FSLDIR', '/Users/tug87422/fsl'))/'bin/film_gls'
        if not film.is_file(): self.skipTest('FSL is not installed')
        rng = np.random.default_rng(8); shape = (6, 6, 6); t = 150
        with tempfile.TemporaryDirectory() as d:
            root = Path(d); source = root/'source'; (source/'masks').mkdir(parents=True)
            reference = nib.Nifti1Image(np.ones(shape, dtype='float32'), np.eye(4))
            completed = {}; rows = []; designs = {}; identities = {}
            for roi in n.g.core.ROIS:
                mask = np.ones(shape, dtype='uint8')
                if roi == 'left': mask[3:] = 0
                if roi == 'right': mask[:3] = 0
                nib.save(nib.Nifti1Image(mask, reference.affine), source/'masks'/(roi+'.nii.gz'))
                for mode in n.g.inputs.MODES:
                    data = rng.normal(size=(*shape, 5)).astype('float32'); data[mask == 0] = 0
                    path = source/f'{mode}_{roi}.nii.gz'; nib.save(nib.Nifti1Image(data, reference.affine), path)
                    completed[f'{mode}/{roi}/dim-5'] = dict(oriented_maps=str(path))
            for j, (task, tr) in enumerate([('trust', 1.615), ('ugr', 1.615), ('ugr', 1.710)]):
                feat = root/f'feat{j}'; feat.mkdir(); nib.save(reference, feat/'mean_func.nii.gz')
                design = rng.normal(size=(t, 13)); design -= design.mean(0); designs[str(feat)] = design
                n.vest(feat/'design.mat', design)
                con = np.zeros((18, 13)); con[9 if task == 'trust' else 10, :9] = ([0, 0, 0, -1, 1, -1, 1, -1, 1] if task == 'trust' else [0, 1, 0, 1, 0, 1, 0, 1, 0])
                n.vest(feat/'design.con', con)
                settings = {'fmri(evs_real)': '10', 'fmri(temphp_yn)': '0', 'fmri(ndelete)': '0'}
                for i, title in enumerate(('C_def', 'C_rec', 'F_def', 'F_rec', 'S_def', 'S_rec'), 4): settings[f'fmri(evtitle{i})'] = title
                (feat/'design.fsf').write_text('\n'.join(f'set {k} "{v}"' for k, v in settings.items()))
                y = (1000+rng.normal(size=(*shape, t))).astype('float32')
                image = nib.Nifti1Image(y, reference.affine); image.header.set_zooms((1, 1, 1, tr)); nib.save(image, feat/'filtered_func_data.nii.gz')
                (feat/'confoundevs.txt').write_text('synthetic source inspector is tested separately')
                row = dict(key=f'run{j}', subject='synthetic', task=task, tr=tr, run=j+1, nvol=t,
                           active_feat=str(feat), filtered_data=str(feat/'filtered_func_data.nii.gz'),
                           confounds=str(feat/'confoundevs.txt'), contrast_index=10 if task == 'trust' else 11)
                rows.append(row); identities[str(feat)] = [n.g.core.identity(feat/'filtered_func_data.nii.gz', False)]
            exploration_record = root/'exploration'; exploration_record.mkdir()
            (exploration_record/'source_snapshot.json').write_text(json.dumps(dict(source_record=str(source))))
            record = root/'record'; record.mkdir(); (record/'errors').mkdir(); output = root/'output'; output.mkdir()
            args = SimpleNamespace(exploration_record=exploration_record, publish=False, pilot_only=True, jobs=1)
            state = dict(runs={})
            inspect = lambda row, ref: (designs[row['active_feat']], identities[row['active_feat']])
            with patch.object(n.exploration, 'snapshot', return_value=(source, rows, completed)), patch.object(n.g.inputs, 'inspect_source', side_effect=inspect), patch.dict(os.environ, FSLDIR=str(film.parents[1])):
                n.run(args, record, output, state)
                self.assertTrue(state['pilot_passed']); self.assertEqual(len(state['runs']), 3)
                path = output/'run0/joint_bold_bilateral/stage1_timecourses.tsv'; before = path.stat().st_mtime_ns
                n.run(args, record, output, state)
                self.assertEqual(path.stat().st_mtime_ns, before)
                self.assertFalse(any((root/'feat0').glob('stats/*')))


if __name__ == '__main__': unittest.main()
