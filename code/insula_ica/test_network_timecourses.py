"""Checks for the replacement's data direction, standardization, OLS and the five requested bars."""
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

    def test_bar_contrasts_match_saved_EV_order_and_average_after_standardization(self):
        rng = np.random.default_rng(1); design = rng.normal(size=(150, 13)); settings = {'fmri(evs_real)': '10'}
        for i, title in enumerate(('C_def', 'C_rec', 'F_def', 'F_rec', 'S_def', 'S_rec'), 4): settings[f'fmri(evtitle{i})'] = title
        saved = np.zeros((18, 13)); saved[9, :9] = [0, 0, 0, -1, 1, -1, 1, -1, 1]
        cs, labels = n.contrasts_for(dict(task='trust', contrast_index=10), settings, design, saved)
        self.assertEqual(labels, list(n.BAR_LABELS[2:]))
        np.testing.assert_allclose(cs.sum(0), saved[9])
        beta = np.arange(13); np.testing.assert_allclose(cs@beta, [1, 1, 1])
        settings['fmri(evtitle4)'] = 'wrong'
        with self.assertRaisesRegex(ValueError, 'EV order'): n.contrasts_for(dict(task='trust', contrast_index=10), settings, design, saved)
        titles = ('nonsocial_high_constant', 'nonsocial_high_pmod', 'nonsocial_low_constant', 'nonsocial_low_pmod',
                  'social_high_constant', 'social_high_pmod', 'social_low_constant', 'social_low_pmod')
        for i, title in enumerate(titles, 1): settings[f'fmri(evtitle{i})'] = title
        saved[10, :9] = [0, 1, 0, 1, 0, 1, 0, 1, 0]
        cs, labels = n.contrasts_for(dict(task='ugr', contrast_index=11), settings, design, saved)
        self.assertEqual(labels[:2], list(n.BAR_LABELS[:2])); np.testing.assert_allclose(cs@beta, [2, 6, 1, 3, 5, 7])

    def test_standardizes_all_model_columns_and_response_before_joint_OLS(self):
        rng = np.random.default_rng(5); t = 170
        x = rng.normal(size=(t, 7))*np.arange(1, 8)+np.arange(7)*10
        x[:, 6] = 25; x[:, 5] = 0
        tc = rng.normal(size=(t, 5))+x[:, :2]@rng.normal(size=(2, 5))
        cs = np.eye(7)[:3]
        with tempfile.TemporaryDirectory() as d:
            folder = Path(d)/'fit'; effects, beta = n.fit_timecourses(tc, x, cs, folder)
            zx = np.zeros_like(x); zx[:, :5] = (x[:, :5]-x[:, :5].mean(0))/x[:, :5].std(0, ddof=1)
            zy = (tc-tc.mean(0))/tc.std(0, ddof=1)
            expected = np.linalg.lstsq(zx, zy, rcond=1e-10)[0]
            np.testing.assert_allclose(beta, expected, atol=1e-12)
            np.testing.assert_allclose(effects, cs@expected, atol=1e-12)
            np.testing.assert_allclose(np.load(folder/'standardized_design.npy')[:, :5].std(0, ddof=1), 1)
            np.testing.assert_allclose(np.load(folder/'standardized_timecourses.npy').std(0, ddof=1), 1)
            self.assertEqual(json.loads((folder/'standardization.json').read_text())['constant_design_columns_zeroed'], [5, 6])
            bad = np.eye(7)[6:]
            with self.assertRaisesRegex(ValueError, 'not estimable'): n.fit_timecourses(tc, x, bad, Path(d)/'bad')

    def test_group_uses_participants_and_exactly_five_bars_per_component(self):
        rng = np.random.default_rng(20); rows = []
        subjects = [str(k) for k in range(132)]
        for mode in n.g.inputs.MODES:
            for roi in n.g.core.ROIS:
                for label in n.BAR_LABELS:
                    task = 'ugr' if label.startswith('ugr') else 'trust'
                    for k in range(1, 6):
                        for j, subject in enumerate(subjects):
                            value = rng.normal()
                            for run, delta in ([(1, 0)] if j%2 else [(1, -2), (2, 2)]):
                                rows.append(dict(mode=mode, roi=roi, task=task, contrast=label, component=k,
                                                 subject=subject, run=run, effect=value+delta))
        participant, group = n.summarize(rows, subjects)
        self.assertEqual(len(group), 150); self.assertEqual(len(participant), 150*132)
        self.assertTrue(all(r['plotted'] for r in group))
        first = group[0]; estimates = [r['effect'] for r in participant if all(r[f] == first[f] for f in ('mode', 'roi', 'task', 'component', 'contrast'))]
        self.assertAlmostEqual(first['effect'], np.mean(estimates))
        with self.assertRaisesRegex(ValueError, 'Missing participant'): n.summarize(rows[:-2], subjects)

    def test_three_run_pilot_reads_4d_and_resumes_without_any_source_task_map(self):
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
                settings = {'fmri(evs_real)': '10', 'fmri(temphp_yn)': '0', 'fmri(ndelete)': '0', 'feat_files(1)': str(feat/'input_to_FEAT.nii.gz')}
                ugr_titles = ('nonsocial_high_constant', 'nonsocial_high_pmod', 'nonsocial_low_constant', 'nonsocial_low_pmod', 'social_high_constant', 'social_high_pmod', 'social_low_constant', 'social_low_pmod')
                for i, title in enumerate(('C_def', 'C_rec', 'F_def', 'F_rec', 'S_def', 'S_rec'), 4): settings[f'fmri(evtitle{i})'] = title
                if task == 'ugr':
                    for i, title in enumerate(ugr_titles, 1): settings[f'fmri(evtitle{i})'] = title
                (feat/'design.fsf').write_text('\n'.join(f'set {k} "{v}"' for k, v in settings.items()))
                y = (1000+rng.normal(size=(*shape, t))).astype('float32')
                image = nib.Nifti1Image(y, reference.affine); image.header.set_zooms((1, 1, 1, tr)); nib.save(image, feat/'input_to_FEAT.nii.gz')
                # A different processed series must not be used accidentally.
                nib.save(nib.Nifti1Image(np.zeros_like(y), reference.affine), feat/'filtered_func_data.nii.gz')
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
            with patch.object(n.exploration, 'snapshot', return_value=(source, rows, completed)), patch.object(n.g.inputs, 'inspect_source', side_effect=inspect):
                n.run(args, record, output, state)
                self.assertTrue(state['pilot_passed']); self.assertEqual(len(state['runs']), 3)
                mask = nib.load(source/'masks/bilateral.nii.gz').get_fdata() > 0
                maps = nib.load(source/'joint_bold_bilateral.nii.gz').get_fdata()[mask]
                data = nib.load(root/'feat0/input_to_FEAT.nii.gz').get_fdata()[mask]
                expected = n.task_timecourses(n.exploration.weights(maps), data)
                np.testing.assert_allclose(np.loadtxt(output/'run0/joint_bold_bilateral/stage1_timecourses.tsv'), expected, atol=1e-10)
                path = output/'run0/joint_bold_bilateral/stage1_timecourses.tsv'; before = path.stat().st_mtime_ns
                n.run(args, record, output, state)
                self.assertEqual(path.stat().st_mtime_ns, before)
                self.assertFalse(any((root/'feat0').glob('stats/*')))


if __name__ == '__main__': unittest.main()
