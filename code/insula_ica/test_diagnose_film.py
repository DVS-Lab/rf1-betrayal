"""Tests for interpreting agreement in the same-input FILM control."""
import unittest
import numpy as np
from diagnose_film_replay import agreement


class FilmDiagnosticTests(unittest.TestCase):
    def test_perfect_spatial_correlation_does_not_hide_amplitude_change(self):
        a = np.array([-2., -1., 1., 2.])
        result = agreement(a, 2*a)
        self.assertAlmostEqual(result['spatial_r'], 1.)
        self.assertAlmostEqual(result['relative_rms_difference'], 1.)

    def test_identical_maps_pass_and_constant_zero_maps_are_explicit(self):
        a = np.array([-2., -1., 1., 2.])
        self.assertEqual(agreement(a, a)['relative_rms_difference'], 0.)
        result = agreement(np.zeros(4), np.zeros(4))
        self.assertIsNone(result['spatial_r'])
        self.assertIsNone(result['relative_rms_difference'])
        self.assertEqual(result['max_absolute_difference'], 0.)


if __name__ == '__main__':
    unittest.main()
