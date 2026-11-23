"""Shear at fixed fractional sites and the full-frequency scalar adjoint check."""
import unittest
from types import SimpleNamespace
from unittest.mock import patch
import numpy as np

from lapw_soc_response import SOCMatching, SOCSite
from lapw_total_stress_driver import solve_charge_adjoint
from response_parallel import SERIAL
from stress_response_adjoint import PairedOperator, ScreenedResponse


class StrainChecks(unittest.TestCase):
    def test_general_matching_in_local_frame_with_lo_and_gamma(self):
        # Regenerate the native matching coefficients after deforming all
        # reciprocal vectors. Include q=0 and a nontrivial orbital rotation.
        rng = np.random.default_rng(197)
        frac = rng.integers(-3, 4, (24, 3)).astype(float)
        frac[0] = 0
        coeff = rng.normal(size=(24, 6)) + 1j*rng.normal(size=(24, 6))
        coeff[-5:, :2] *= 9e-7
        b = np.array([[.8, .11, .03], [0., .95, -.27], [0., 0., 1.1]])
        axes, _ = np.linalg.qr(rng.normal(size=(3, 3)))
        ell = 2
        ctx = SimpleNamespace(rmt=2.35, lapw=True, primitive=rng.normal(size=(1, 4, 3)),
                              lo_coefficients=rng.normal(size=(3, 1)))
        for legacy in (False, True):
            def build(e):
                f = np.eye(3) + e
                br = np.linalg.solve(f.T, b)
                geo = SimpleNamespace(br1=br, nlo=5, volume=153*np.linalg.det(f),
                    channels={(0, ell): SimpleNamespace(ctx=ctx, slots=[0, 1, 2])},
                    lo_indices={(0, ell): [np.arange(5)]}, legacy_lo_filter=legacy)
                site = SOCSite(0, axes @ br, np.array([.13, -.2, .45]), {})
                return SOCMatching(geo, site, frac, coeff, ell)
            central = build(np.zeros((3,3)))
            for i,j in zip(*np.triu_indices(3)):
                e = np.zeros((3,3))
                e[i,j] = e[j,i] = 1.
                actual = central.action(None, None, np.zeros_like(ctx.primitive[-1]),
                                        np.zeros_like(ctx.lo_coefficients), strain=e)
                dh = 1e-6
                numeric = (build(dh*e).z - build(-dh*e).z)/(2*dh)
                np.testing.assert_allclose(actual, numeric, atol=1e-8, rtol=2e-7)

    def test_full_operator_residual_recorded_without_frequency_compression(self):
        jacobian = np.array([[0.3, 0.1], [0.0, -0.2]])
        source = np.array([1.0, -0.5])
        op = PairedOperator(2, 2, lambda x: jacobian @ x, lambda x: jacobian.T @ x)
        branch = SimpleNamespace(records=[0], model=SimpleNamespace())
        # Either compression is disabled, or no k point benefits from it.
        for factor_tolerance in (0., 1e-7):
            with self.subTest(factor_tolerance=factor_tolerance), patch(
                "lapw_frequency_modes.prepare_frequency_modes", return_value={}
            ):
                (solution, report), stats = solve_charge_adjoint(
                    ScreenedResponse(op), source, [branch], SERIAL, 1., factor_tolerance,
                    lambda _: None,
                )
                residual = np.linalg.norm(source - (solution - jacobian.T @ solution))
                self.assertAlmostEqual(residual, report["residual_norm"], places=14)
                self.assertLessEqual(residual, report["threshold"])
                self.assertEqual(stats["compressed_branch_kpoints"], 0)
                self.assertTrue(report["full_frequency_check"])
                self.assertFalse(report["full_frequency_fallback"])

    def test_exact_residual_and_restoration_for_two_spin_channels(self):
        grad = np.array([1., -1.5])
        exact = np.diag([.23, -.1])
        for altered in (False, True):
            branches = [SimpleNamespace(records=[0], model=SimpleNamespace(generalized_adjoint_modes=None))
                        for _ in range(2)]
            def jt(x):
                fitted = all(b.model.generalized_adjoint_modes for b in branches)
                return np.zeros_like(x) if fitted and altered else exact @ x
            op = PairedOperator(2, 2, lambda x: exact @ x, jt)
            table = SimpleNamespace(signs=[1], error=1e-8, width=5)
            with patch('lapw_frequency_modes.prepare_frequency_modes', return_value={0:table}):
                (lam, report), stats = solve_charge_adjoint(ScreenedResponse(op), grad, branches,
                                                             SERIAL, 2, 1e-7, lambda line: None)
            self.assertEqual(stats["compressed_branch_kpoints"], 2)
            self.assertTrue(report["full_frequency_check"])
            self.assertEqual(report["full_frequency_fallback"], altered)
            self.assertTrue(all(b.model.generalized_adjoint_modes is None for b in branches))
            np.testing.assert_allclose(lam, np.linalg.solve(np.eye(2)-exact, grad), atol=1e-12)
