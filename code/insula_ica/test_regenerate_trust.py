from pathlib import Path
import json
import tempfile
import unittest
import os
from unittest.mock import patch

import nibabel as nib
import numpy as np

import regenerate_trust as g
import rebuild_trust as r
from test_insula_ica import vest


class RegenerateTrustTests(unittest.TestCase):
    def test_nonidentical_statistics_allowed_but_design_changes_are_not(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); old, new = root / 'old', root / 'new'
            for tree, value in ((old, 1.), (new, 2.)):
                (tree / 'stats').mkdir(parents=True)
                vest(tree / 'design.mat', np.ones((12, 1)))
                vest(tree / 'design.con', np.ones((1, 1)))
                for name in ('cope1', 'zstat1'):
                    nib.save(nib.Nifti1Image(np.full((2, 2, 2), value), np.eye(4)), tree / 'stats' / (name + '.nii.gz'))
            nib.save(nib.Nifti1Image(np.ones((2, 2, 2, 12)), np.eye(4)), new / 'stats/res4d.nii.gz')
            row = dict(key='fixture', original_feat=str(old), nvol=12)
            r.compare_fit(row, new, strict_historical=False)
            self.assertFalse(all(x['matched'] for x in json.loads((new / 'historical_comparison.json').read_text())['comparisons']))
            with self.assertRaisesRegex(ValueError, 'COPE/ZSTAT'):
                r.compare_fit(row, new)
            vest(new / 'design.con', np.full((1, 1), 2.))
            with self.assertRaisesRegex(ValueError, 'design/contrasts'):
                r.compare_fit(row, new, strict_historical=False)

    def test_exact_nuisance_values_remain_required(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); old, new = root / 'old', root / 'new'
            old.mkdir(); new.mkdir()
            values = np.arange(20.).reshape(10, 2)
            for path in (old / 'confoundevs.txt', new / 'confoundevs.txt', root / 'saved.tsv'):
                np.savetxt(path, values)
            row = dict(key='fixture', confounds=str(root / 'saved.tsv'), original_feat=str(old),
                       confounds_sha256=g.core.sha(root / 'saved.tsv'))
            self.assertTrue(g.exact_confounds(row, new)['exact_numeric_equality'])
            values[0, 0] += 1e-9
            np.savetxt(new / 'confoundevs.txt', values)
            with self.assertRaisesRegex(ValueError, 'nuisance values differ'):
                g.exact_confounds(row, new)

    def test_production_uses_24_with_recorded_seed_and_normal_thread_resources(self):
        cmd = g.production_command('apptainer', Path('/tools/fmriprep-24.1.1.simg'), Path('/upstream'),
                                   Path('/output'), Path('/work'), 'sub-10478')
        self.assertIn('/tools/fmriprep-24.1.1.simg', cmd)
        self.assertEqual(cmd[cmd.index('--nthreads')+1], '14')
        self.assertEqual(cmd[cmd.index('--omp-nthreads')+1], '8')
        self.assertEqual(cmd[cmd.index('--random-seed')+1], str(g.SEED))
        self.assertIn('--skull-strip-fixed-seed', cmd)
        self.assertIn('/upstream/bids:/bids:ro', cmd)

    def test_completed_control_is_linked_without_reprocessing_or_overwriting(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); subject = 'sub-10478'; key = subject + '_ses-01_task-trust_run-1'
            output = root / 'derivatives/insula_ica/rebuild'
            source = root / 'derivatives/insula_ica/controls/case/fmriprep' / subject
            func = source / 'ses-01/func'; func.mkdir(parents=True)
            name = key + '_part-mag_space-MNI152NLin6Asym_desc-preproc_bold.nii.gz'
            bold = func / name
            image = nib.Nifti1Image(np.ones((2, 2, 2, 12)), np.eye(4)); image.header.set_zooms((1, 1, 1, 2))
            nib.save(image, bold)
            old = root / 'old'; old.mkdir()
            nib.save(nib.Nifti1Image(np.ones((2, 2, 2)), np.eye(4)), old / 'mean_func.nii.gz')
            record = root / 'code/insula_ica/records/linux2-bold-controls-20261027T113051Z-b1efa8'
            recovered = record / 'recovered-20261027T153848Z-31fc23'; recovered.mkdir(parents=True)
            (recovered / 'status.json').write_text(json.dumps(dict(cases={subject+'-v24-a': dict(returncode=0, bold_outputs={key: str(bold)})})))
            (recovered / 'bold-output-sha256.json').write_text(json.dumps({str(bold): dict(sha256=g.core.sha(bold))}))
            (record / 'containers.json').write_text('[{"package_versions":{"fmriprep":"24.1.1"},"sha256":"image-sha"}]')
            row = dict(key=key, subject=subject, new_bold=str(output / 'fmriprep-24' / subject / 'ses-01/func' / name),
                       original_feat=str(old), nvol=12, tr=2.)
            before = g.core.sha(bold)
            with patch.object(g.core, 'ROOT', root):
                result = g.reuse_control(subject, [row], output, 'image-sha')
                self.assertEqual(result['source'], 'completed_v24_a_control')
                self.assertTrue((output / 'fmriprep-24' / subject).is_symlink())
                self.assertEqual(g.core.sha(row['new_bold']), before)
                with self.assertRaisesRegex(ValueError, 'Existing incomplete'):
                    g.reuse_control(subject, [row], output, 'image-sha')
            self.assertEqual(g.core.sha(bold), before)

    def test_qc_flags_low_coverage_without_excluding_and_rejects_nonfinite_data(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); old, target = root / 'old', root / 'feat_reconstruction/fixture.feat'
            old.mkdir(); (target / 'stats').mkdir(parents=True)
            values = np.arange(96., dtype=float).reshape(2, 2, 2, 12) + 100
            residuals = values.copy(); residuals[0, 0, 0] = 0
            bold = root / 'bold.nii.gz'
            nib.save(nib.Nifti1Image(values, np.eye(4)), bold)
            nib.save(nib.Nifti1Image(residuals, np.eye(4)), target / 'stats/res4d.nii.gz')
            for directory in (old, target):
                nib.save(nib.Nifti1Image(np.ones((2, 2, 2)), np.eye(4)), directory / 'mask.nii.gz')
                np.savetxt(directory / 'confoundevs.txt', np.arange(12.).reshape(12, 1))
            nib.save(nib.Nifti1Image(values[..., 6], np.eye(4)), old / 'example_func.nii.gz')
            saved = root / 'saved.tsv'; np.savetxt(saved, np.arange(12.).reshape(12, 1))
            (target / 'historical_comparison.json').write_text('{"comparisons":[{"matched":false}]}')
            row = dict(key='fixture', subject='sub-10001', new_bold=str(bold), original_feat=str(old),
                       nvol=12, confounds=str(saved), confounds_sha256=g.core.sha(saved))
            left = np.zeros((2, 2, 2), bool); left[0] = True
            masks = dict(left=left, right=~left, bilateral=np.ones((2, 2, 2), bool))
            with patch.dict(os.environ, {'FSLDIR': str(root)}), patch.object(g.core, 'anatomical_masks', return_value=(masks, {}, 1)):
                qc = g.qc_run(row, root)
                self.assertIn('left_insula_coverage_below_0.8', qc['review_flags'])
                self.assertEqual(qc['constant_mask_voxels'], 1)
                self.assertTrue(qc['confounds']['exact_numeric_equality'])
                self.assertEqual(qc['historical_statistics_matched'], 0)
                residuals[1, 1, 1, 0] = np.nan
                nib.save(nib.Nifti1Image(residuals, np.eye(4)), target / 'stats/res4d.nii.gz')
                with self.assertRaisesRegex(ValueError, 'Nonfinite'):
                    g.qc_run(row, root)


if __name__ == '__main__':
    unittest.main()
