"""XC scalar derivatives: physical limits, discrete fields and vacuum policy."""

import unittest
from pathlib import Path
import tempfile
import numpy as np
from lapw_xc_response import unpolarized_xc, MuffinTinXCResponse, InterstitialXCResponse
from lapw_functional_cancellation import PreparedXCEnergy
from lapw_potential_core_state import save_potential_core_state
from test_lapw_potential_core_state import inputs
from test_lapw_xc_response import mt, ir


class XCEnergyTests(unittest.TestCase):
    def test_uniform_gas_energy_potential_relation(self):
        n = np.geomspace(1e-4, 1e3, 60)
        for functional in ["LDA", "PBE"]:

            def e(x):
                return unpolarized_xc(
                    x,
                    0,
                    0,
                    0,
                    functional=functional,
                    charge_cut=1e-8,
                    quantity="energy_per_electron",
                )

            derivative = e(n) + n * np.imag(e(n + 1j * 1e-22)) / 1e-22
            v = unpolarized_xc(n, 0, 0, 0, functional=functional, charge_cut=1e-8)
            np.testing.assert_allclose(derivative, v, rtol=1e-14, atol=1e-14)

    def test_mt_energy_field_finite_difference_and_transpose(self):
        for functional in ["LDA", "PBE"]:
            base, n = mt(functional)

            def prepare(c):
                extra = (
                    {}
                    if functional == "LDA"
                    else dict(
                        sine=base.sine.ravel(),
                        cosine=base.cosine.ravel(),
                        angular_derivatives=base.angular[1:],
                    )
                )
                return MuffinTinXCResponse(
                    base.r,
                    base.ells,
                    c,
                    base.Y,
                    base.weights,
                    functional=functional,
                    quantity="energy_per_electron",
                    **extra
                )

            e = prepare(n)
            d = n * np.array([0.2, -0.3])
            h = 0.0003
            fd = (prepare(n + h * d).central - prepare(n - h * d).central) / (2 * h)
            np.testing.assert_allclose(e.action(d), fd, rtol=2e-7, atol=2e-9)
            self.assertLess(e.operator.check_transpose(samples=2), 1e-12)

    def test_ir_energy_field_finite_difference_and_transpose(self):
        for functional in ["LDA", "PBE"]:
            base, n = ir(functional)

            def prepare(c):
                return InterstitialXCResponse(
                    base.reanalysis,
                    c,
                    [8.0, 9.0, 10.0],
                    functional=functional,
                    quantity="energy_per_electron",
                )

            e = prepare(n)
            d = np.array([0.01, -0.02, 0.005, 0.002], complex)
            h = 0.0003
            np.testing.assert_allclose(
                e.action(d),
                (prepare(n + h * d).central - prepare(n - h * d).central) / (2 * h),
                rtol=1e-7,
                atol=1e-10,
            )
            self.assertLess(e.operator.check_transpose(samples=2), 1e-12)

    def test_energy_does_not_use_potential_vacuum_average(self):
        base, n = ir("LDA", rho_min=(0.16, 0.14))
        e = InterstitialXCResponse(
            base.reanalysis,
            n,
            [8.0, 9.0, 10.0],
            functional="LDA",
            rho_min=(0.16, 0.14),
            quantity="energy_per_electron",
        )
        self.assertGreater(e.average_count, 0)
        x = np.random.default_rng(3).normal(size=e.shape)
        np.testing.assert_array_equal(e.average(x), x)
        np.testing.assert_array_equal(e.average(x, True), x)
        self.assertGreater(np.linalg.norm(base.average(x) - x), 1)

    def test_complete_xc_density_term_derivative(self):
        with tempfile.TemporaryDirectory() as tmp:
            data = inputs()
            n = data["density"]
            d = n.vector() * 0.1
            h = 0.0001

            def prepare(t):
                args = data | dict(density=n.from_vector(n.vector() + t * d))
                return PreparedXCEnergy(
                    save_potential_core_state(Path(tmp) / "central.npz", **args)
                )

            e, p, m = prepare(0), prepare(h), prepare(-h)
            self.assertAlmostEqual(
                e.energy_gradient @ d, (p.value_ry - m.value_ry) / (2 * h), places=8
            )
            self.assertAlmostEqual(
                e.density_term_gradient @ d,
                (p.density_term_value_ry - m.density_term_value_ry) / (2 * h),
                places=8,
            )

    def test_invalid_observable_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "quantity"):
            unpolarized_xc(
                1, 0, 0, 0, functional="LDA", charge_cut=1e-8, quantity="unknown"
            )


if __name__ == "__main__":
    unittest.main()
