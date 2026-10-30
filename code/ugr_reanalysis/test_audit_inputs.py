import unittest
import numpy as np
from audit_inputs import contrast_diagnostics


class ContrastAuditTests(unittest.TestCase):
    def test_original_contrast_units_survive_column_scaling(self):
        rng = np.random.default_rng(2); a, b = rng.normal(size=(2, 90))
        design = np.column_stack([a, 7*b, 2*a])
        result = contrast_diagnostics(design, np.array([[1, 0, 2], [1, 0, 0]]))
        self.assertTrue(all(not r['requires_review'] for r in result if r['contrast'] == 1))
        self.assertTrue(all(r['requires_review'] for r in result if r['contrast'] == 2))

    def test_constant_column_requires_zero_contrast_weight(self):
        rng = np.random.default_rng(3); design = np.column_stack([rng.normal(size=(90, 2)), np.ones(90)])
        result = contrast_diagnostics(design, np.array([[1, -1, 0], [0, 0, 1]]))
        self.assertTrue(all(not r['requires_review'] for r in result if r['contrast'] == 1))
        self.assertTrue(all(r['weights_on_constant_columns'] for r in result if r['contrast'] == 2))

    def test_full_rank_model_all_contrasts_estimable(self):
        design = np.random.default_rng(4).normal(size=(90, 12))
        self.assertTrue(all(not r['requires_review'] for r in contrast_diagnostics(design, np.eye(12))))


if __name__ == '__main__': unittest.main()
