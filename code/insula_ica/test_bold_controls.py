"""Controls must be independent and comparisons must capture temporal changes."""
from pathlib import Path
import argparse
import contextlib
import io
import json
import sys
import tempfile
import unittest
from unittest.mock import patch

import nibabel as nib
import numpy as np

from compare_bold_controls import frame_metrics, series_metrics
from run_bold_controls import ARMS, SEED, command, select_runs
import run_bold_controls as controls
from recover_bold_controls import completed_case, recover


class BoldControlsTests(unittest.TestCase):
    def make_images(self, root, a, b, affine=None):
        affine = np.eye(4) if affine is None else affine
        for name, values in [('a', a), ('b', b)]:
            nib.save(nib.Nifti1Image(values, affine), root / (name + '.nii.gz'))
        nib.save(nib.Nifti1Image(np.ones(a.shape[:3]), affine), root / 'mask.nii.gz')
        return root / 'a.nii.gz', root / 'b.nii.gz', root / 'mask.nii.gz'

    def test_identical_series_preserves_precision_and_returns_exact_equality(self):
        with tempfile.TemporaryDirectory() as folder:
            a = np.random.default_rng(42).normal(size=(2, 2, 2, 20)) + 1000
            paths = self.make_images(Path(folder), a, a.copy())
            m = series_metrics(*paths)
            self.assertTrue(m['exact_numeric_equality'])
            self.assertEqual(m['relative_rms_difference'], 0.)
            self.assertEqual(m['voxel_demeaned_relative_rms'], 0.)
            self.assertAlmostEqual(m['voxel_temporal_correlation_percentiles']['median'], 1.)
            paths = self.make_images(Path(folder), a, a + 1e-8)
            self.assertFalse(series_metrics(*paths)['exact_numeric_equality'])

    def test_uniform_scale_and_temporal_change_are_distinguished(self):
        with tempfile.TemporaryDirectory() as folder:
            a = np.random.default_rng(2).normal(size=(2, 2, 2, 30)) + 100
            m = series_metrics(*self.make_images(Path(folder), a, a*2))
            self.assertAlmostEqual(m['relative_rms_difference'], 1.)
            self.assertAlmostEqual(m['global_scale'], 2.)
            self.assertAlmostEqual(m['relative_rms_after_global_scale'], 0.)
            self.assertAlmostEqual(m['voxel_demeaned_relative_rms'], 1.)
            y = a.mean(axis=3, keepdims=True) + 2*(a-a.mean(axis=3, keepdims=True))
            m = series_metrics(*self.make_images(Path(folder), a, y))
            self.assertLess(m['relative_rms_difference'], .02)
            self.assertAlmostEqual(m['voxel_demeaned_relative_rms'], 1.)

    def test_grid_mismatch_is_not_silently_resampled(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            a = np.ones((2, 2, 2, 5))
            paths = self.make_images(root, a, a)
            affine = np.eye(4); affine[0, 3] = 2
            nib.save(nib.Nifti1Image(a, affine), paths[1])
            self.assertFalse(series_metrics(*paths)['compared'])

    def test_historical_frame_uses_requested_index(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            frame = np.arange(1., 9).reshape(2, 2, 2)
            a = np.stack([frame*2, frame, frame*5], axis=3)
            _, bold, mask = self.make_images(root, a, a)
            nib.save(nib.Nifti1Image(frame, np.eye(4)), root / 'ref.nii.gz')
            self.assertEqual(frame_metrics(root / 'ref.nii.gz', bold, mask, 1)['relative_rms_difference'], 0.)

    def test_repeat_commands_share_controls_but_not_work_or_output(self):
        args = ('apptainer', Path('/tools/fmriprep-24.1.1.simg'), Path('/upstream'))
        a = command(*args, Path('/a'), Path('/work-a'), 'sub-10402', 14, 28000)
        b = command(*args, Path('/b'), Path('/work-b'), 'sub-10402', 14, 28000)
        for cmd in (a, b):
            self.assertEqual(cmd[cmd.index('--random-seed')+1], str(SEED))
            self.assertEqual(cmd[cmd.index('--omp-nthreads')+1], '1')
            self.assertIn('--skull-strip-fixed-seed', cmd)
            self.assertIn('/upstream/bids:/bids:ro', cmd)
            self.assertIn('/ZPOOL/data/tools/templateflow:/opt/templateflow:ro', cmd)
        self.assertNotEqual(a, b)
        self.assertNotIn('/work-a:/scratch', b)
        self.assertEqual(ARMS, {'v24-a': '24.1.1', 'v24-b': '24.1.1', 'v25-a': '25.2.5'})
        ordinary = command(*args, Path('/ordinary'), Path('/work-ordinary'), 'sub-10402', 14, 28000, controlled=False)
        self.assertNotIn('--random-seed', ordinary)
        self.assertNotIn('--skull-strip-fixed-seed', ordinary)
        self.assertEqual(ordinary[ordinary.index('--omp-nthreads')+1], '8')

    def test_extra_subject_is_selected_without_outcome_information(self):
        rows = [dict(subject=s, run=r, warpkit=dict(ready=True)) for s in
                ('sub-10369', 'sub-10402', 'sub-10478', 'sub-10486') for r in (1, 2)]
        plan = dict(commands=[dict(subject='sub-10369'), dict(subject='sub-10402')], all_runs=rows)
        selected = select_runs(plan)
        self.assertEqual({r['subject'] for r in selected}, {'sub-10402', 'sub-10478'})

    def test_real_case_collection_after_successful_subprocess(self):
        # Exercise the actual launcher success branch, rather than mocking it
        # away as the experiment-wide fixture does. This catches the undefined
        # output-directory variable that affected the first linux2 collection.
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            case, record = root / 'case', root / 'record'
            record.mkdir()
            key = 'sub-10402_ses-01_task-trust_run-1'
            name = key + '_space-MNI152NLin6Asym_desc-preproc_bold.nii.gz'
            bold = case / 'fmriprep/sub-10402/ses-01/func' / name
            bold.parent.mkdir(parents=True)
            nib.save(nib.Nifti1Image(np.ones((2, 2, 2, 3)), np.eye(4)), bold)
            row = dict(key=key, subject='sub-10402', new_bold=str(root / name), nvol=3)
            item = dict(name='fixture', case=str(case), controlled=True, argv=[sys.executable, '-c', 'pass'], runs=[row])
            status = dict(cases={})
            with contextlib.redirect_stdout(io.StringIO()):
                result = controls.launch_case(item, {}, record, status)
            self.assertNotIn('error', result)
            self.assertEqual(result['returncode'], 0)
            self.assertEqual(result['bold_outputs'][key], str(bold))
            self.assertEqual(json.loads((record / 'status.json').read_text())['cases']['fixture']['returncode'], 0)
            # A zero exit or an output file alone is not sufficient for recovery
            # when the saved fMRIPrep success marker is absent.
            with self.assertRaisesRegex(ValueError, 'success marker absent'):
                completed_case(item, record)

    def test_end_to_end_control_report_with_fake_preprocessing_preserves_anchors(self):
        with tempfile.TemporaryDirectory() as folder:
            base = Path(folder)
            rebuild, upstream, record = [base / s for s in ('rebuild', 'upstream', 'record')]
            for path in (rebuild, upstream, record): path.mkdir()
            data = np.random.default_rng(8).normal(size=(2, 2, 2, 6)) + 1000
            rows, anchors = [], {}
            for subject in ('sub-10402', 'sub-10478'):
                for run in (1, 2):
                    key = f'{subject}_ses-01_task-trust_run-{run}'
                    old = base / (key + '.feat'); old.mkdir()
                    nib.save(nib.Nifti1Image(data[..., 3], np.eye(4)), old / 'example_func.nii.gz')
                    nib.save(nib.Nifti1Image(np.ones(data.shape[:3]), np.eye(4)), old / 'mask.nii.gz')
                    anchors.update({str(p): controls.sha(p) for p in old.iterdir()})
                    rows.append(dict(subject=subject, key=key, nvol=6, run=run, warpkit=dict(ready=True),
                        original_feat=str(old), new_bold=str(base / (key + '_space-MNI152NLin6Asym_desc-preproc_bold.nii.gz'))))
            (rebuild / 'plan-pilot.json').write_text(json.dumps(dict(all_runs=rows,
                commands=[dict(subject='sub-10369'), dict(subject='sub-10402')])))
            (rebuild / 'container.json').write_text('{"sha256":"fixture"}')
            for row in rows:
                if row['subject'] == 'sub-10402':
                    nib.save(nib.Nifti1Image(data, np.eye(4)), row['new_bold'])
            args = argparse.Namespace(rebuild_root=rebuild, upstream_root=upstream, tools_root=base / 'tools',
                jobs=8, nprocs=14, mem_mb=28000, execute=True)
            inventory = [dict(path=str(args.tools_root / f'fmriprep-{v}.simg'), size=1,
                             sha256='fixture', package_versions=dict(fmriprep=v)) for v in ('24.1.1', '25.2.5')]

            def fake_launch(item, env, record, status):
                output = {}
                for row in item['runs']:
                    path = Path(item['case']) / 'fmriprep' / row['subject'] / 'ses-01/func' / Path(row['new_bold']).name
                    path.parent.mkdir(parents=True, exist_ok=True)
                    nib.save(nib.Nifti1Image(data * (1.01 if item['arm'] == 'v25-a' else 1.), np.eye(4)), path)
                    output[row['key']] = str(path)
                logs = Path(item['case']) / 'fmriprep' / item['subject'] / 'log/fixture'; logs.mkdir(parents=True)
                fixed = 'true' if item['controlled'] else 'false'
                omp = 1 if item['controlled'] else 8
                (logs / 'fmriprep.toml').write_text(f'[seeds]\nmaster={SEED}\n[workflow]\nskull_strip_fixed_seed={fixed}\n[nipype]\nomp_nthreads={omp}\n')
                result = dict(name=item['name'], bold_outputs=output, returncode=0)
                status['cases'][item['name']] = result
                return result

            status = dict(state='preflight', cases={})
            with patch.object(controls, 'ROOT', base), patch.object(controls, 'WORK_BASE', base / 'scratch'), \
                 patch.object(controls.shutil, 'which', return_value='fake-apptainer'), \
                 patch.object(controls, 'image_inventory', return_value=inventory), \
                 patch.object(controls, 'input_snapshot', return_value={}), \
                 patch.object(controls, 'verify_snapshot', return_value=[]), \
                 patch.object(controls.subprocess, 'check_output', return_value='fixture-commit'), \
                 patch.object(controls, 'launch_case', side_effect=fake_launch), \
                 contextlib.redirect_stdout(io.StringIO()):
                controls.run(args, record, status)
            self.assertEqual(status['state'], 'complete')
            self.assertEqual(len(status['cases']), 8)
            self.assertTrue(all(controls.sha(p) == h for p, h in anchors.items()))
            results = json.loads((record / 'results.json').read_text())
            for row in results.values():
                self.assertTrue(row['full_series']['v24_repeat']['exact_numeric_equality'])
                self.assertAlmostEqual(row['full_series']['v25_vs_v24']['relative_rms_difference'], .01)
            self.assertTrue((record / 'report.md').is_file())
            self.assertEqual(len(json.loads((record / 'plan.json').read_text())['cases']), 8)
            # Recollect the saved fixture without any preprocessing calls and
            # preserve every byte of the original experiment record.
            (record / 'input-sha256.json').write_text('{}\n')
            for item in json.loads((record / 'plan.json').read_text())['cases']:
                (record / (item['name'] + '.fmriprep-log.txt')).write_text('fMRIPrep finished successfully!\n')
            original = {str(p): controls.sha(p) for p in record.rglob('*') if p.is_file()}
            with patch('recover_bold_controls.subprocess.check_output', return_value='fixture-commit'), \
                 contextlib.redirect_stdout(io.StringIO()):
                recovered = recover(record)
            self.assertEqual(json.loads((recovered / 'status.json').read_text())['state'], 'complete')
            self.assertTrue(all(controls.sha(p) == h for p, h in original.items()))
            self.assertTrue((recovered / 'results.json').is_file())
            self.assertTrue((recovered / 'bold-output-sha256.json').is_file())


if __name__ == '__main__':
    unittest.main()
