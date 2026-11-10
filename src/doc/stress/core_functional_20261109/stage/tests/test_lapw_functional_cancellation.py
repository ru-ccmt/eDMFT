"""Core units/degeneracy, spatial conventions and incomplete-cancellation guards."""

from types import SimpleNamespace
import unittest
import numpy as np
from lapw_functional_cancellation import (
    SpatialPotentialPairing,
    CoreEigenvalueFunctional,
    stationary_potential_readout,
)
from lapw_core_response import (
    CorePotentialExtension,
    MatchedCoreState,
    AtomCoreResponse,
    PeriodicCoreResponse,
)
from lapw_hartree_periodic import PeriodicHartreeResponse, native_charge_weights
from test_lapw_hartree_periodic import fixture


def setup():
    d, s, _ = fixture()
    layout = PeriodicHartreeResponse(d, s, (8, 8, 8), inversion=False)
    return d, s, layout


def core_fixture(rela=False, exterior=False):
    r = np.geomspace(1e-5, 1.3 if exterior else 20.0, 1201)
    z = 5.0
    occ = 2.0
    ext = CorePotentialExtension(
        r, np.full(len(r), -2 * z), extra_points=400 if exterior else 0
    )
    gamma = np.sqrt(1 - (z / 137.0359895) ** 2) if rela else 1.0
    energy = -z * z / (1 + gamma)
    if exterior:
        # The native two-point exterior slope is not the exact Coulomb slope.
        # Solve its actual matching condition, rather than assigning a
        # continuum Coulomb eigenvalue to a different exterior potential.
        from lapw_core_response import _CoreTrajectory, _origin, _outer
        from scipy.optimize import root_scalar

        rr = ext.r
        v = ext.central
        c = 137.0359895 if rela else 1e30
        match = int(np.argmin(abs(rr - 0.3)))
        end = len(rr) - 1

        def residual(e):
            yo = _CoreTrajectory(
                rr,
                v,
                e,
                -1,
                c,
                np.arange(5),
                np.arange(5, match + 1),
                _origin(rr, z, -1, c, e, v[0]),
            ).values
            ii = np.arange(end, end - 5, -1)
            yi = _CoreTrajectory(
                rr,
                v,
                e,
                -1,
                c,
                ii,
                np.arange(end - 5, match - 1, -1),
                _outer(rr[ii], c, e),
            ).values
            return yo[match, 1] / yo[match, 0] - yi[match, 1] / yi[match, 0]

        energy = root_scalar(residual, x0=energy, x1=energy + 0.001, xtol=1e-12).root
    st = MatchedCoreState(
        ext.r,
        ext.central,
        energy,
        z=z,
        kappa=-1,
        match=int(np.argmin(abs(ext.r - 0.3))),
        end=len(ext.r) - 1,
        occupation=occ,
        relativistic=rela,
        output_points=len(r),
    )
    layout = SimpleNamespace(
        integrations=[SimpleNamespace(r=r)],
        size=len(r),
        keys=((0, 0, 0),),
        slices=[slice(0, len(r))],
    )
    core = PeriodicCoreResponse(layout, {0: AtomCoreResponse(ext, [st])})
    return core, r, z


class FunctionalCancellationTests(unittest.TestCase):
    def test_spatial_constant_potential_counts_mt_and_ir_once(self):
        d, s, h = setup()
        p = SpatialPotentialPairing(h, s.mult)
        density = np.zeros(h.size)
        v = np.zeros(h.size)
        r = h.integrations[0].r
        density[h.slices[0]] = 4 * np.pi * r * r * 0.03
        density[h.pw_offset] = 0.03
        v[h.slices[0]] = r
        constant = np.zeros(h.stars.nstar, complex)
        constant[0] = 1
        v[h.pw_offset :] = h.reanalysis.action(constant).view(float)
        # Native radial CHARGE integral plus exact analytic interstitial volume.
        mt = native_charge_weights(r) @ density[h.slices[0]]
        exact = mt + 0.03 * (h.volume - 4 * np.pi * r[-1] ** 3 / 3)
        self.assertAlmostEqual(p.value(density, v), exact, places=12)

    def test_spatial_complex_pairing_and_bilinear_derivatives(self):
        _, s, h = setup()
        p = SpatialPotentialPairing(h, s.mult)
        rng = np.random.default_rng(41)
        n, v, dn, dv = rng.normal(size=(4, h.size))
        step = 0.003
        numeric = (
            p.value(n + step * dn, v + step * dv)
            - p.value(n - step * dn, v - step * dv)
        ) / (2 * step)
        analytic = p.gradient(n) @ dv + p.gradient(v) @ dn
        np.testing.assert_allclose(numeric, analytic, rtol=1e-11, atol=1e-9)
        expected = (
            h.volume
            * np.vdot(
                n[h.pw_offset :].view(complex), v[h.pw_offset :].view(complex)
            ).real
        )
        self.assertAlmostEqual(
            float(p.gradient(n)[h.pw_offset :] @ v[h.pw_offset :]), expected, places=11
        )
        self.assertLess(p.operator.check_transpose(samples=2), 1e-13)

    def test_spatial_rejects_coordinate_reordering(self):
        d, s, h = setup()
        p = SpatialPotentialPairing(h, s.mult)
        d.rad = dict(reversed(list(d.rad.items())))
        with self.assertRaisesRegex(ValueError, "coordinates differ"):
            p.density_vector(d)

    def test_core_units_occupations_multiplicity_and_fixed_shift(self):
        core, r, z = core_fixture()
        g = CoreEigenvalueFunctional(core, [3], [0.7])
        self.assertAlmostEqual(g.value_ry, 3 * 2 * (-z * z + 0.7), places=12)
        self.assertEqual(g.core_number, 6)
        # dV_Ry/dZ=-2/r, and hence d(rV_Ry)/dZ=-2.
        self.assertAlmostEqual(
            g.action(np.full(len(r), -2.0)), 3 * 2 * (-2 * z), places=5
        )
        np.testing.assert_array_equal(
            g.potential_gradient,
            CoreEigenvalueFunctional(core, [3], [-0.8]).potential_gradient,
        )

    def test_relativistic_core_coulomb_scalar(self):
        core, r, z = core_fixture(True)
        g = CoreEigenvalueFunctional(core, [1], [0])
        c = 137.0359895
        gamma = np.sqrt(1 - (z / c) ** 2)
        self.assertAlmostEqual(g.value_ry, -4 * z * z / (1 + gamma), places=12)
        # Includes the fixed-Z origin parameter limitation at finite R0.
        np.testing.assert_allclose(
            g.action(np.full(len(r), -2.0)), -4 * z / gamma, rtol=3e-6
        )

    def test_core_exterior_is_retained(self):
        core, r, z = core_fixture(exterior=True)
        g = CoreEigenvalueFunctional(core, [1], [0])
        d = r * np.sin(r)
        self.assertGreater(abs(g.exterior_gradient @ d), 1e-5)
        self.assertAlmostEqual(
            g.action(d),
            float((g.interior_gradient + g.exterior_gradient) @ d),
            places=14,
        )
        self.assertEqual(np.count_nonzero(g.exterior_gradient[:-2]), 0)

    def test_empty_core_and_bad_weights(self):
        _, s, h = setup()
        core = PeriodicCoreResponse(h, {}, empty_atoms=[0])
        g = CoreEigenvalueFunctional(core, [1], [0])
        self.assertEqual(g.value_ry, 0)
        np.testing.assert_array_equal(g.potential_gradient, np.zeros(h.size))
        for mult in [[0], [1.5], [np.nan]]:
            with self.assertRaises(ValueError):
                CoreEigenvalueFunctional(core, mult, [0])
            with self.assertRaises(ValueError):
                SpatialPotentialPairing(h, mult)

    def test_stationary_identity_uses_input_and_never_claims_full_stress(self):
        from lapw_potential_core_state import save_potential_core_state
        from test_lapw_potential_core_state import inputs
        from pathlib import Path
        import tempfile

        with tempfile.TemporaryDirectory() as tmp:
            pc = save_potential_core_state(Path(tmp) / "state.npz", **inputs())
            pairing = SpatialPotentialPairing(pc.hartree, pc.structure.mult)
            spatial = pairing.gradient(pc.input_density.vector())
            result = stationary_potential_readout(pc, spatial)
            np.testing.assert_array_equal(result["potential"], spatial * 0)
            np.testing.assert_array_equal(result["charge"], spatial * 0)
            for key in [
                "stress_applied",
                "baseline_subtracted",
                "full_klein_functional",
                "native_density_functional_validated",
            ]:
                self.assertFalse(result[key])


if __name__ == "__main__":
    unittest.main()
