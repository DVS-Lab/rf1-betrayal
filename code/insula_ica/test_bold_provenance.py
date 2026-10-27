"""Provenance checks must preserve unknowns and leave source data untouched."""
import contextlib
import io
import json
from pathlib import Path
import tarfile
import tempfile
import unittest

import nibabel as nib
import numpy as np

from diagnose_bold_provenance import collect, frame_comparison, frozen_check, setting_differences, task_scope
from rebuild_trust import sha


class BoldProvenanceTests(unittest.TestCase):
    def test_missing_setting_is_distinct_from_explicit_null(self):
        diff = setting_differences({'workflow': {}}, {'workflow': {'use_bbr': None}})
        self.assertFalse(diff['workflow.use_bbr']['historical_present'])
        self.assertTrue(diff['workflow.use_bbr']['rebuilt_present'])

    def test_compare_unanticipated_workflow_keys_and_seeds(self):
        old = {'workflow': {'new_combination_option': 'a'}, 'seeds': {'ants': 1}}
        new = {'workflow': {'new_combination_option': 'b'}, 'seeds': {'ants': 2}}
        self.assertEqual(set(setting_differences(old, new)), {'workflow.new_combination_option', 'seeds.ants'})

    def test_config_task_coverage_is_not_guessed(self):
        self.assertEqual(task_scope({}), 'task_not_recorded')
        self.assertEqual(task_scope({'execution': {'task_id': 'ugr'}}), 'other_task')
        self.assertEqual(task_scope({'execution': {'bids_filters': {'bold': {'task': ['trust']}}}}), 'trust_named')

    def test_frozen_json_detects_content_change_even_with_same_size(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / 'metadata.json'
            path.write_text('{"EchoTime": 1}')
            expected = {'sha256': sha(path)}
            self.assertTrue(frozen_check(path, expected)['matches_rebuild_snapshot'])
            path.write_text('{"EchoTime": 2}')
            self.assertFalse(frozen_check(path, expected)['matches_rebuild_snapshot'])
            self.assertFalse(frozen_check(Path(folder) / 'missing', expected)['matches_rebuild_snapshot'])

    def test_echo_frame_check_uses_requested_frame_and_does_not_resample(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            values = np.arange(1, 9, dtype=float).reshape(2, 2, 2)
            nib.save(nib.Nifti1Image(values, np.eye(4)), root / 'ref.nii.gz')
            series = np.stack([values * 9, values, values * 3], axis=3)
            nib.save(nib.Nifti1Image(series, np.eye(4)), root / 'bold.nii.gz')
            match = frame_comparison(root / 'ref.nii.gz', root / 'bold.nii.gz', np.ones_like(values), 1)
            self.assertEqual(match['relative_rms_difference'], 0.)
            affine = np.eye(4); affine[0, 3] = 10
            nib.save(nib.Nifti1Image(series, affine), root / 'shifted.nii.gz')
            mismatch = frame_comparison(root / 'ref.nii.gz', root / 'shifted.nii.gz', np.ones_like(values), 1)
            self.assertFalse(mismatch['compared'])

    def test_complete_bundle_preserves_sources_and_includes_prior_report(self):
        with tempfile.TemporaryDirectory() as folder:
            base = Path(folder)
            rebuild, upstream, output, old = [base / name for name in ('rebuild', 'upstream', 'output', 'old.feat')]
            for path in (rebuild, upstream, output, old): path.mkdir()
            key = 'sub-10001_ses-01_task-trust_run-1'
            func = rebuild / 'fmriprep-24/sub-10001/ses-01/func'
            func.mkdir(parents=True)
            values = np.arange(1, 9, dtype=float).reshape(2, 2, 2)
            nib.save(nib.Nifti1Image(values, np.eye(4)), old / 'example_func.nii.gz')
            nib.save(nib.Nifti1Image(np.ones_like(values), np.eye(4)), old / 'mask.nii.gz')
            bold = func / (key + '_space-MNI152NLin6Asym_desc-preproc_bold.nii.gz')
            nib.save(nib.Nifti1Image(np.stack([values]*4, axis=3), np.eye(4)), bold)
            row = dict(key=key, subject='sub-10001', nvol=4, original_feat=str(old), new_bold=str(bold))
            (rebuild / 'plan-pilot.json').write_text(json.dumps(dict(all_runs=[row], commands=[dict(runs=[key])])))
            (rebuild / 'container.json').write_text('{"version": "24.1.1"}')
            diagnostics = rebuild / 'diagnostics'; diagnostics.mkdir()
            (diagnostics / 'trust-preprocessing-evidence-20261027.json').write_text('{"preserved": true}')
            for label, folder in [('historical', upstream / 'derivatives/fmriprep-24'), ('rebuilt', rebuild / 'fmriprep-24')]:
                logs = folder / 'sub-10001/log/2026'; logs.mkdir(parents=True)
                option = 'loglin' if label == 'historical' else 'curvefit'
                (logs / 'fmriprep.toml').write_text('[workflow]\nme_t2s_fit_method = "' + option + '"\n[execution]\ntask_id = "trust"\n')
                (folder / 'sub-10001.html').write_text('<p>fMRIPrep 24.1.1</p>')
            before = {str(p): sha(p) for p in base.rglob('*') if p.is_file()}
            with contextlib.redirect_stdout(io.StringIO()): collect(rebuild, upstream, output)
            self.assertTrue(all(sha(p) == h for p, h in before.items()))
            result = json.loads((output / 'evidence.json').read_text())
            diff = result['subjects']['sub-10001']['config_comparisons'][0]['differences']
            self.assertEqual(diff['workflow.me_t2s_fit_method']['historical'], 'loglin')
            self.assertEqual(result['subjects']['sub-10001']['historical']['subject_report_fmriprep_versions'], ['24.1.1'])
            self.assertEqual(json.loads((output / 'previous-preprocessing-evidence.json').read_text()), {'preserved': True})
            self.assertEqual(result['runs'][key]['rebuilt_frame_candidates'][0]['relative_rms_difference'], 0.)
            self.assertEqual(result['errors'], [])
            for line in (output / 'SHA256SUMS').read_text().splitlines():
                digest, name = line.split('  ', 1)
                self.assertEqual(digest, sha(output / name))
            with tarfile.open(output.with_suffix('.tar.gz')) as archive:
                names = archive.getnames()
                self.assertIn('output/evidence.json', names)
                self.assertFalse(any(name.endswith(('.nii', '.nii.gz')) for name in names))


if __name__ == '__main__':
    unittest.main()
