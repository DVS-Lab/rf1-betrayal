"""Evidence checks must distinguish missing provenance from verified identity."""
from pathlib import Path
import tempfile
import unittest

import nibabel as nib
import numpy as np

from inspect_trust_preprocessing import checkpoint, matrix_match, provenance, tedana_check


class PreprocessingEvidenceTests(unittest.TestCase):
    def test_confounds_require_numeric_identity_not_just_matching_dimensions(self):
        with tempfile.TemporaryDirectory() as folder:
            a, b = Path(folder) / 'a', Path(folder) / 'b'
            np.savetxt(a, [[1., 2.], [3., 4.]])
            np.savetxt(b, [[1., 2.], [3., 5.]])
            self.assertFalse(matrix_match(a, b)['equal'])
            self.assertIsNone(matrix_match(a, Path(folder) / 'absent')['equal'])

    def test_tedana_rejected_timecourses_must_match_saved_suffix_in_original_order(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder);prefix = 'sub-10001_ses-01_task-trust_run-1'
            (root / (prefix + '_desc-tedana_metrics.tsv')).write_text('Component\tclassification\nICA_000\taccepted\nICA_001\trejected\nICA_002\trejected\n')
            mixing = np.arange(12).reshape(4, 3).astype(float)
            np.savetxt(root / (prefix + '_desc-ICA_mixing.tsv'), mixing, delimiter='\t',
                       header='ICA_000\tICA_001\tICA_002', comments='')
            saved = root / 'saved'
            np.savetxt(saved, np.column_stack([np.ones(4), mixing[:, [1, 2]]]))
            self.assertTrue(tedana_check(root, prefix, saved)['rejected_suffix_matches'])
            np.savetxt(saved, np.column_stack([np.ones(4), mixing[:, [2, 1]]]))
            self.assertFalse(tedana_check(root, prefix, saved)['rejected_suffix_matches'])

    def test_absent_metadata_and_tedana_are_unknown_not_confirmation(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            check = tedana_check(root, 'missing', root / 'missing')
            self.assertFalse(check['available'])
            self.assertIsNone(check['rejected_suffix_matches'])
            self.assertEqual(provenance(root, 'sub-10001')['execution_configs'], [])

    def test_example_checkpoint_separates_uniform_scaling_from_spatial_change(self):
        with tempfile.TemporaryDirectory() as folder:
            old, new = Path(folder) / 'old', Path(folder) / 'new'
            old.mkdir();new.mkdir()
            values = np.arange(1, 9, dtype=float).reshape(2, 2, 2)
            nib.save(nib.Nifti1Image(values, np.eye(4)), old / 'example_func.nii.gz')
            nib.save(nib.Nifti1Image(2*values, np.eye(4)), new / 'example_func.nii.gz')
            result = checkpoint(old, new, 'example_func.nii.gz', np.ones((2, 2, 2), bool))
            self.assertAlmostEqual(result['relative_rms_difference'], 1.)
            self.assertAlmostEqual(result['new_on_old_global_scale'], 2.)
            self.assertAlmostEqual(result['relative_rms_after_global_scale'], 0.)


if __name__ == '__main__':
    unittest.main()
