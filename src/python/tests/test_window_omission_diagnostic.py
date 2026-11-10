"""Default singular-pair omission, adjoints and optional cutoff diagnostics."""
import unittest
import warnings
import numpy as np
from frozen_sigma_tangent import CachedFrozenSigmaResponse, real_inner
from test_frozen_sigma_tangent import complex_model


class WindowDiagnosticTests(unittest.TestCase):
    def make(self, cutoff=None, gap=0.0):
        m = complex_model()
        e = m.eps.copy()
        e[0, 4] = e[0, 3] + gap
        return CachedFrozenSigmaResponse(
            m.meta, e, m.raw, m.apply_x(m.raw, m.X), m.kw, m.W,
            diagnostic_gap_cutoff=cutoff, cache_mb=2,
        )

    def test_default_warns_and_retains_regular_pairs(self):
        with self.assertWarnsRegex(RuntimeWarning, 'omitted 1 singular.*1 k points'):
            model = self.make()
        self.assertFalse(model.diagnostic_only)
        report = model.window_gap_report()
        self.assertEqual(report['omitted_pairs'], 1)
        self.assertEqual(report['affected_kpoints'], 1)
        self.assertAlmostEqual(report['affected_kpoint_weight'], .3)
        for k, denominator in enumerate(model.window_denominators()):
            exact = model.eps[k, model.W][None, :] - model.eps[k, model.Q][:, None]
            regular = abs(exact) > 1e-10
            np.testing.assert_array_equal(denominator[regular], exact[regular])
            np.testing.assert_array_equal(1/denominator[~regular], 0)

    def test_default_matches_tested_diagnostic_and_transpose(self):
        with self.assertWarns(RuntimeWarning):
            default = self.make()
        diagnostic = self.make(1e-10)
        rng = np.random.default_rng(377)
        shape = (2, 5, 5)
        h = rng.normal(size=shape) + 1j*rng.normal(size=shape)
        h = (h+h.swapaxes(1, 2).conj())/2
        y = rng.normal(size=shape) + 1j*rng.normal(size=shape)
        y = (y+y.swapaxes(1, 2).conj())/2
        for fixed_n in (False, True):
            actual, _ = default.response(h, fixed_number=fixed_n)
            expected, _ = diagnostic.response(h, fixed_number=fixed_n)
            np.testing.assert_allclose(actual, expected, atol=1e-14)
            np.testing.assert_allclose(default.adjoint(y, fixed_number=fixed_n),
                                       diagnostic.adjoint(y, fixed_number=fixed_n), atol=1e-14)
            self.assertAlmostEqual(real_inner(y, actual),
                real_inner(default.adjoint(y, fixed_number=fixed_n), h), places=11)
        np.testing.assert_array_equal(default.central_density(), diagnostic.central_density())
        np.testing.assert_array_equal(default.u, diagnostic.u)
        np.testing.assert_array_equal(default.X, diagnostic.X)
        np.testing.assert_array_equal(default.kw, diagnostic.kw)
        self.assertTrue(diagnostic.diagnostic_only)

    def test_regular_gap_is_unchanged_and_silent(self):
        with warnings.catch_warnings(record=True) as emitted:
            warnings.simplefilter('always')
            regular = self.make(gap=1e-6)
        self.assertFalse(emitted)
        self.assertEqual(regular.window_gap_report()['omitted_pairs'], 0)
        self.assertIsNone(regular.window_gap_warning())
        for k, denominator in enumerate(regular.window_denominators()):
            np.testing.assert_array_equal(denominator,
                regular.eps[k, regular.W][None, :] - regular.eps[k, regular.Q][:, None])

    def test_threshold_and_diagnostic_plateau(self):
        with self.assertWarns(RuntimeWarning):
            self.make(gap=5e-11)
        a, b = self.make(1e-12), self.make(1e-8)
        for x, y in zip(a.window_denominators(), b.window_denominators()):
            np.testing.assert_array_equal(x, y)
        for cutoff in (0, -1e-10, float('nan'), float('inf')):
            with self.assertRaisesRegex(ValueError, 'finite and positive'):
                self.make(cutoff)


if __name__ == '__main__':
    unittest.main()
