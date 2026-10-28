"""Safety checks for costly Trust preprocessing and historical-fit recovery."""
from pathlib import Path
import json
import tempfile
import unittest

import nibabel as nib
import numpy as np

import rebuild_trust as r
from test_insula_ica import vest


class TrustRebuildTests(unittest.TestCase):
    def test_container_migration_accepts_changed_notices_but_preserves_legacy_record(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            image = root / '24.sif'
            image.write_bytes(b'container fixture')
            first = 'Proxy warning\nfMRIPrep 24.1.1\nA newer version is available: 25.2.6\n'
            second = 'Network check failed\nfMRIPrep 24.1.1\nA newer version is available: 25.2.7\n'
            record = root / 'container.json'
            r.save_json(record, dict(path=str(image), sha256=r.sha(image), version=first))
            legacy = record.read_bytes()
            r.record_container(root, image, second)
            self.assertEqual((root / 'container.legacy.json').read_bytes(), legacy)
            stable = record.read_bytes()
            r.record_container(root, image, first)
            self.assertEqual(record.read_bytes(), stable)
            self.assertEqual(json.loads(stable)['version'], '24.1.1')
            self.assertEqual(len((root / 'container_observations.jsonl').read_text().splitlines()), 2)

    def test_container_digest_change_still_blocks_without_overwriting_record(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            image = root / '24.sif'
            image.write_bytes(b'original')
            r.record_container(root, image, 'fMRIPrep 24.1.1\n')
            original = (root / 'container.json').read_bytes()
            image.write_bytes(b'changed image')
            with self.assertRaisesRegex(ValueError, 'SHA256 changed'):
                r.record_container(root, image, 'fMRIPrep 24.1.1\n')
            self.assertEqual((root / 'container.json').read_bytes(), original)

    def test_version_is_the_application_line_not_the_new_release_notice(self):
        self.assertEqual(r.parsed_version('fMRIPrep v24.1.1\nLatest release: 25.2.7\n'), '24.1.1')
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            with self.assertRaisesRegex(ValueError, 'not fMRIPrep'):
                r.record_container(root, root / 'missing', 'fMRIPrep 25.2.5\nOld release: 24.1.1\n')
        with self.assertRaisesRegex(ValueError, 'unique'):
            r.parsed_version('fMRIPrep 24.1.1\nfMRIPrep 25.2.5\n')

    def test_production_uses_eight_subjects_and_pilot_retains_two(self):
        self.assertEqual(r.subject_jobs(False, None), 8)
        self.assertEqual(r.subject_jobs(True, None), 2)
        self.assertEqual(r.subject_jobs(False, 4), 4)
        with self.assertRaisesRegex(ValueError, 'concurrent'):
            r.subject_jobs(False, 0)

    def warpkit_fixture(self, folder):
        upstream = Path(folder)
        sub = 'sub-10001'
        func = upstream / 'bids' / sub / 'ses-01/func'
        fmap = func.parent / 'fmap'
        func.mkdir(parents=True)
        fmap.mkdir()
        raw = []
        for echo in range(1, 5):
            path = func / f'{sub}_ses-01_task-trust_run-1_echo-{echo}_part-mag_bold.nii.gz'
            raw.append(path)
            path.touch()
            Path(str(path)[:-7] + '.json').write_text('{}')
        stem = fmap / (sub + '_ses-01_acq-trust_run-1')
        for suffix in ('_fieldmap.nii.gz', '_magnitude.nii.gz'):
            nib.save(nib.Nifti1Image(np.ones((2, 2, 2)), np.eye(4)), Path(str(stem) + suffix))
        Path(str(stem) + '_magnitude.json').write_text('{}')
        metadata = Path(str(stem) + '_fieldmap.json')
        metadata.write_text(json.dumps(dict(Units='Hz', IntendedFor=[p.relative_to(upstream / 'bids' / sub).as_posix() for p in raw])))
        return upstream, sub, raw, metadata

    def test_reuses_bids_warpkit_products_without_derivative_work_files(self):
        with tempfile.TemporaryDirectory() as folder:
            upstream, sub, raw, metadata = self.warpkit_fixture(folder)
            result = r.warpkit_inputs(upstream, sub, 1, raw)
            self.assertTrue(result['ready'])
            self.assertFalse(result['retained_4d_exists'])
            self.assertEqual(len(result['files']), 4)
            metadata.write_text(json.dumps(dict(Units='Hz', IntendedFor=[])))
            result = r.warpkit_inputs(upstream, sub, 1, raw)
            self.assertFalse(result['ready'])
            self.assertIn('not associated', result['problems'][0])

    def test_warpkit_checks_units_and_b0_identifier_links(self):
        with tempfile.TemporaryDirectory() as folder:
            upstream, sub, raw, metadata = self.warpkit_fixture(folder)
            for p in raw:
                Path(str(p)[:-7] + '.json').write_text(json.dumps(dict(B0FieldSource='trust1')))
            metadata.write_text(json.dumps(dict(Units='Hz', B0FieldIdentifier='trust1')))
            self.assertTrue(r.warpkit_inputs(upstream, sub, 1, raw)['ready'])
            metadata.write_text(json.dumps(dict(Units='rad/s', B0FieldIdentifier='trust1')))
            result = r.warpkit_inputs(upstream, sub, 1, raw)
            self.assertFalse(result['ready'])
            self.assertIn('Units', result['problems'][0])

    def test_filter_keeps_historical_anatomy_and_fieldmaps_but_only_accepted_trust(self):
        selected = r.filter_for([{'run': 2}])
        self.assertEqual(selected['bold']['task'], 'trust')
        self.assertEqual(selected['bold']['run'], ['2'])
        self.assertNotIn('run', selected['sbref'])
        self.assertEqual(selected['t1w']['session'], '01')
        self.assertEqual(selected['fmap'], {'datatype': 'fmap'})
        self.assertEqual(selected['bold']['part'], [None, 'mag'])

    def test_filter_entities_are_compatible_with_24_series_string_deserializer(self):
        # 24.1.1 applies `'Query' in value` to every non-enum scalar after
        # converting JSON null to Query.NONE. Integers therefore fail even
        # though PyBIDS itself would accept numeric runs.
        selected = json.loads(json.dumps(r.filter_for([{'run': 2}, {'run': 1}, {'run': '2'}])))
        for query in selected.values():
            for value in query.values():
                for leaf in value if isinstance(value, list) else [value]:
                    if leaf is not None:
                        self.assertIsInstance(leaf, str)
        self.assertEqual(selected['bold']['run'], ['1', '2'])

    def test_command_does_not_write_upstream_or_run_other_tasks(self):
        cmd = r.command('apptainer', Path('/images/24.sif'), Path('/upstream'),
                        Path('/isolated'), Path('/scratch/sub-10001'), 'sub-10001', 14, 28000)
        self.assertIn('/upstream/bids:/bids:ro', cmd)
        self.assertEqual(cmd[cmd.index('--task-id') + 1], 'trust')
        self.assertIn('/rebuild/fmriprep-24', cmd)
        self.assertIn('--fs-no-reconall', cmd)
        self.assertIn('--me-output-echos', cmd)
        self.assertEqual(cmd[cmd.index('--output-spaces') + 1:cmd.index('--fs-no-reconall')],
                         ['T1w', 'MNI152NLin6Asym'])
        self.assertNotIn('--overwrite', cmd)
        self.assertNotIn('tedana', cmd)

    def test_fsf_changes_only_three_locations(self):
        text = ('# untouched\nset feat_files(1) "/old/bold"\n'
                'set confoundev_files(1) "/old/conf"\nset fmri(outputdir) "/old/out"\n'
                'set fmri(prewhiten_yn) 1\nset fmri(smooth) 5\n')
        updated = r.recovery_fsf(text, '/new/bold', '/new/conf', '/new/out')
        self.assertEqual(updated.replace('/new/', '/old/'), text)
        with self.assertRaisesRegex(ValueError, 'ambiguous'):
            r.recovery_fsf(text + 'set feat_files(1) "/another"\n', 'b', 'c', 'd')
        with self.assertRaisesRegex(ValueError, 'Tcl'):
            r.recovery_fsf(text, '/bad/$path', '/c', '/o')

    def test_preserved_nuisance_matrix_cannot_be_silently_replaced(self):
        with tempfile.TemporaryDirectory() as folder:
            source, dest = Path(folder) / 'source', Path(folder) / 'saved'
            source.write_text('1 2\n3 4\n')
            r.preserve(source, dest)
            self.assertEqual(dest.read_bytes(), source.read_bytes())
            source.write_text('5 6\n7 8\n')
            with self.assertRaisesRegex(ValueError, 'changed'):
                r.preserve(source, dest)
            self.assertEqual(dest.read_text(), '1 2\n3 4\n')

    def test_input_snapshot_rejects_changed_raw_identity(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / 'inputs.json'
            r.freeze_json(path, {'input': {'size': 500, 'mtime_ns': 1}})
            with self.assertRaisesRegex(ValueError, 'changed'):
                r.freeze_json(path, {'input': {'size': 500, 'mtime_ns': 2}})

    def test_historical_statistic_change_fails_comparison_and_is_reported(self):
        with tempfile.TemporaryDirectory() as folder:
            old, new = Path(folder) / 'old', Path(folder) / 'new'
            for tree in (old, new):
                (tree / 'stats').mkdir(parents=True)
                vest(tree / 'design.mat', np.ones((12, 1)))
                vest(tree / 'design.con', np.ones((1, 1)))
                for name in ('cope1', 'zstat1'):
                    nib.save(nib.Nifti1Image(np.ones((2, 2, 2)), np.eye(4)), tree / 'stats' / (name + '.nii.gz'))
            nib.save(nib.Nifti1Image(np.ones((2, 2, 2, 12)), np.eye(4)), new / 'stats/res4d.nii.gz')
            run = dict(key='synthetic', original_feat=str(old), nvol=12)
            r.compare_fit(run, new)
            nib.save(nib.Nifti1Image(np.full((2, 2, 2), 2.), np.eye(4)), new / 'stats/cope1.nii.gz')
            with self.assertRaisesRegex(ValueError, 'COPE/ZSTAT'):
                r.compare_fit(run, new)
            report = json.loads((new / 'historical_comparison.json').read_text())
            self.assertFalse(next(x['matched'] for x in report['comparisons'] if x['map'] == 'cope1.nii.gz'))

    def test_comparison_uses_model_contrasts_and_reports_extra_historical_maps(self):
        with tempfile.TemporaryDirectory() as folder:
            old, new = Path(folder) / 'old', Path(folder) / 'new'
            for tree in (old, new):
                (tree / 'stats').mkdir(parents=True)
                vest(tree / 'design.mat', np.ones((12, 1)))
                vest(tree / 'design.con', np.ones((1, 1)))
                for name in ('cope1', 'zstat1'):
                    nib.save(nib.Nifti1Image(np.ones((2, 2, 2)), np.eye(4)), tree / 'stats' / (name + '.nii.gz'))
            nib.save(nib.Nifti1Image(np.ones((2, 2, 2, 12)), np.eye(4)), new / 'stats/res4d.nii.gz')
            nib.save(nib.Nifti1Image(np.full((2, 2, 2), 999.), np.eye(4)), old / 'stats/zstat20.nii.gz')
            run = dict(key='synthetic', original_feat=str(old), nvol=12)
            r.compare_fit(run, new)
            report = json.loads((new / 'historical_comparison.json').read_text())
            self.assertEqual(report['ncontrasts'], 1)
            self.assertEqual(report['outside_original_model'], ['zstat20.nii.gz'])
            self.assertEqual([x['map'] for x in report['comparisons']], ['cope1.nii.gz', 'zstat1.nii.gz'])
            self.assertTrue(all(x['matched'] for x in report['comparisons']))

    def test_missing_model_defined_map_is_never_silently_skipped(self):
        with tempfile.TemporaryDirectory() as folder:
            old, new = Path(folder) / 'old', Path(folder) / 'new'
            for tree in (old, new):
                (tree / 'stats').mkdir(parents=True)
                vest(tree / 'design.mat', np.ones((12, 1)))
                vest(tree / 'design.con', np.ones((2, 1)))
                for name in ('cope1', 'zstat1', 'cope2', 'zstat2'):
                    if tree == new and name == 'zstat2':
                        continue
                    nib.save(nib.Nifti1Image(np.ones((2, 2, 2)), np.eye(4)), tree / 'stats' / (name + '.nii.gz'))
            with self.assertRaisesRegex(ValueError, 'missing model-defined maps: zstat2'):
                r.compare_fit(dict(key='synthetic', original_feat=str(old)), new)
            report = json.loads((new / 'historical_comparison.json').read_text())
            self.assertEqual(report['missing_reconstructed'], ['zstat2.nii.gz'])
            self.assertEqual(report['comparisons'], [])


if __name__ == '__main__':
    unittest.main()
