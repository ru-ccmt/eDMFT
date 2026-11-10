"""Explicit spatial/normalization geometry against nonlinear endpoint maps."""

from copy import deepcopy
from types import SimpleNamespace
import unittest
import numpy as np
from lapw_charge_response import NativeChargeNormalization, native_charge_covectors
from lapw_strain_charge import (
    charge_integration_geometry,
    material_source_from_raw,
    spatial_matching_source,
)
from lapw_density_response import (
    MuffinTinDensitySum,
    SpatialValenceDensity,
    InterstitialDensity,
    FourierStars,
)
from wien_response_density import Density
from test_lapw_charge_response import fixture
from test_lapw_density_response import fixture as density_fixture


class StrainChargeTests(unittest.TestCase):
    def test_exact_ir_integration_covector_geometry(self):
        density, s = fixture()
        e = np.diag([0.4, -0.2, 0.7])
        h = 2e-6
        for inversion in [False, True]:
            actual = charge_integration_geometry(density, s, e, inversion=inversion)
            ends = []
            for sign in [-1, 1]:
                sp = deepcopy(s)
                for name, rate in zip(["a", "b", "c"], np.diag(e)):
                    setattr(sp, name, getattr(s, name) * (1 + sign * h * rate))
                ends.append(
                    native_charge_covectors(density, sp, inversion=inversion)[1]
                )
            np.testing.assert_allclose(
                actual, (ends[1] - ends[0]) / (2 * h), rtol=3e-8, atol=2e-8
            )

    def test_two_normalizations_and_output_inverse_transport(self):
        density, s = fixture()
        s.timat = np.eye(3, dtype=int)[None]
        s.tau = np.zeros((1, 3))
        n = NativeChargeNormalization(density, s, 7.8, inversion=False)
        # Deliberately different input and raw/normalized output. Using input
        # charge in inverse transport must fail this nonlinear reference.
        input_density = density.copy().scale(0.87)
        state = SimpleNamespace(
            normalization=n,
            potential_core=SimpleNamespace(input_density=input_density, structure=s),
        )
        rng = np.random.default_rng(317)
        d = rng.normal(size=n.size) * 0.0002
        e = np.diag([0.9, 1.1, 1.0])
        h = 1e-6
        actual, parts = material_source_from_raw(state, d, e)
        ends = []
        for t in [-h, h]:
            sp = deepcopy(s)
            for name, rate in zip(["a", "b", "c"], np.diag(e)):
                setattr(sp, name, getattr(s, name) * (1 + t * rate))
            raw = density.from_vector(n.base + t * d)
            physical = NativeChargeNormalization(
                raw, sp, n.target, inversion=False
            ).central
            pulled = density.from_vector(physical).transport(1 + t, True)
            ends.append(
                NativeChargeNormalization(pulled, s, n.target, inversion=False).central
            )
        np.testing.assert_allclose(
            actual, (ends[1] - ends[0]) / (2 * h), rtol=3e-7, atol=1e-9
        )
        self.assertLess(abs((n.mt + n.ir) @ actual), 2e-14)
        self.assertGreater(np.linalg.norm(parts["post_normalization_change"]), 1e-4)

    def test_no_normalization_policy_still_transports_output(self):
        density, s = fixture()
        s.timat = np.eye(3, dtype=int)[None]
        s.tau = np.zeros((1, 3))
        n = NativeChargeNormalization(density, s, 7.8, inversion=False, enabled=False)
        state = SimpleNamespace(
            normalization=n,
            potential_core=SimpleNamespace(input_density=density, structure=s),
        )
        h = 1e-6
        d = np.linspace(-0.001, 0.001, n.size)
        actual, parts = material_source_from_raw(state, d, np.eye(3))
        expected = (
            density.from_vector(n.base + h * d).transport(1 + h, True).vector()
            - density.from_vector(n.base - h * d).transport(1 - h, True).vector()
        ) / (2 * h)
        np.testing.assert_allclose(actual, expected, rtol=1e-6, atol=1e-9)
        np.testing.assert_array_equal(parts["normalization_geometry"], 0)
        np.testing.assert_allclose(parts["post_normalization_change"], 0, atol=2e-15)

    def test_spatial_matching_and_volume_with_complex_offdiagonal_rho(self):
        rng, phi, z, rho, hs = density_fixture()
        template = Density.__new__(Density)
        template.rad = {(0, *h): np.zeros(7) for h in hs}
        template.pw = {(0, 0, 0): 0j, (1, 0, 0): 0j, (-1, 0, 0): 0j}
        template.grids = [np.linspace(-5, 1, 7)]
        template.mult = np.array([1])
        stars = FourierStars(
            (7, 7, 7),
            np.array(list(template.pw)),
            np.eye(3, dtype=int)[None],
            np.zeros((1, 3)),
            complex_phases=True,
        )
        labels = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]])
        c = rng.normal(size=(3, 4)) + 1j * rng.normal(size=(3, 4))
        dz = rng.normal(size=z.shape) + 1j * rng.normal(size=z.shape)

        def endpoint(t):
            atom = MuffinTinDensitySum(
                [(0, z + t * dz, 1.0)], rho[None], phi, hs, 0.0001
            )
            ir = InterstitialDensity(c, labels, 33 * (1 + t) ** 3, stars)
            return SpatialValenceDensity(template, {0: atom}, [ir])

        central = endpoint(0)
        actual = spatial_matching_source(central, {0: [dz]}, 3.0)
        h = 2e-6
        expected = (endpoint(h).value() - endpoint(-h).value()) / (2 * h)
        np.testing.assert_allclose(actual, expected, rtol=1e-7, atol=2e-7)


if __name__ == "__main__":
    unittest.main()
