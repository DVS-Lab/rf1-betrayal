"""Pooling is a common raw-unit slope, standardized only after combining EVs."""
from pathlib import Path
import tempfile
import unittest

import numpy as np

import pooled_network_tasks as p


def settings():
    titles = ('nonsocial_high_constant', 'nonsocial_high_pmod', 'nonsocial_low_constant',
              'nonsocial_low_pmod', 'social_high_constant', 'social_high_pmod',
              'social_low_constant', 'social_low_pmod')
    return {f'fmri(evtitle{i})':title for i, title in enumerate(titles, 1)}


class PooledModelTests(unittest.TestCase):
    def test_pools_raw_pmods_preserves_four_constants_and_all_other_columns(self):
        rng = np.random.default_rng(22); x = rng.normal(size=(150, 17)); before = x.copy()
        pooled, cs, mapping = p.pool_design(x, settings())
        self.assertEqual(pooled.shape, (150, 15)); np.testing.assert_array_equal(x, before)
        np.testing.assert_array_equal(pooled[:, 1], x[:, 1]+x[:, 3])
        np.testing.assert_array_equal(pooled[:, 4], x[:, 5]+x[:, 7])
        np.testing.assert_array_equal(pooled[:, [0, 2, 3, 5]], x[:, [0, 2, 4, 6]])
        np.testing.assert_array_equal(pooled[:, 6:], x[:, 8:])
        np.testing.assert_array_equal(cs@np.arange(15), [1, 4])
        self.assertEqual(mapping[4]['original_columns'], [5, 7])

    def test_constant_high_condition_pmod_is_identifiable_when_pooled_with_varying_low(self):
        rng = np.random.default_rng(33); x = rng.normal(size=(180, 17))
        x[:, 5] = -1.458333*x[:, 4]  # exact alias from fixed-amplitude high trials
        self.assertEqual(np.linalg.matrix_rank(x), 16)
        pooled, cs, _ = p.pool_design(x, settings())
        self.assertEqual(np.linalg.matrix_rank(pooled), 15)
        zy = pooled@rng.normal(size=(15, 5))+rng.normal(size=(180, 5))
        with tempfile.TemporaryDirectory() as d:
            effect, beta = p.n.fit_timecourses(zy, pooled, cs, Path(d)/'fit')
            zx, _, _, _ = p.n.zscore_columns(pooled); zresponse, _, _, _ = p.n.zscore_columns(zy)
            reference = np.linalg.lstsq(zx, zresponse, rcond=1e-10)[0]
            np.testing.assert_allclose(beta, reference, atol=1e-12)
            np.testing.assert_allclose(effect, cs@reference, atol=1e-12)

    def test_standardize_after_pooling_is_distinct_from_pooling_standardized_EVs(self):
        rng = np.random.default_rng(44); x = rng.normal(size=(150, 17)); x[:, 1] *= 10
        pooled, _, _ = p.pool_design(x, settings())
        zpooled, _, _, _ = p.n.zscore_columns(pooled)
        zx, _, _, _ = p.n.zscore_columns(x)
        self.assertGreater(np.linalg.norm(zpooled[:, 1]-(zx[:, 1]+zx[:, 3])/2), 1)
        expected = (x[:, 1]+x[:, 3]); expected = (expected-expected.mean())/expected.std(ddof=1)
        np.testing.assert_allclose(zpooled[:, 1], expected)


if __name__ == '__main__': unittest.main()
