import unittest
from unittest.mock import patch
from types import SimpleNamespace
import numpy as np
from wien_response_density import Density
from wien_response_operators import (
    charge_covectors,
    InterstitialNormalization,
    hydrostatic_transport_source,
)


def density_fixture():
    d = Density.__new__(Density)
    d.grids = [np.linspace(-4, 0.4, 41)]
    d.mult = [2]
    r = np.exp(d.grids[0])
    d.rad = {(0, 0, 0): 0.3 * r * r * np.exp(-r), (0, 2, 0): 0.02 * r * r}
    d.pw = {(0, 0, 0): 0.03 + 0j, (1, 0, 0): 0.002 + 0.001j, (0, 1, 1): 0.001 - 0.002j}
    structure = SimpleNamespace(
        lattyp="P",
        alpha=90.0,
        beta=90.0,
        gamma=90.0,
        a=5.0,
        b=6.0,
        c=7.0,
        nat=1,
        mult=[2],
        rmt=[float(r[-1])],
        pos=np.array([[0.1, 0.2, 0.3], [0.6, 0.7, 0.8]]),
        ReadStruct=lambda path: None,
    )
    return d, structure


class DensityOperatorTests(unittest.TestCase):
    def test_normalization_tangent_and_transpose_off_fixed_point(self):
        density, structure = density_fixture()
        with patch("wien_response_operators.Struct", return_value=structure):
            normal = InterstitialNormalization(density, "unused", 7.8)
        self.assertGreater(abs(normal.alpha - 1), 0.1)
        self.assertLess(normal.operator.check_transpose(), 1e-14)
        rng = np.random.default_rng(98)
        direction = rng.normal(size=len(normal.base)) * 0.01
        h = 1e-4
        fd = (
            normal.endpoint(normal.base + h * direction)
            - normal.endpoint(normal.base - h * direction)
        ) / (2 * h)
        np.testing.assert_allclose(fd, normal.operator.matvec(direction), atol=1e-11)
        self.assertAlmostEqual(
            (normal.mt + normal.ir) @ normal.operator.matvec(direction), 0, places=13
        )
        self.assertAlmostEqual(
            (normal.mt + normal.ir) @ normal.endpoint(normal.base), 7.8, places=13
        )

    def test_transport_derivative_matches_same_spline_endpoints(self):
        density, _ = density_fixture()
        for pullback in [False, True]:
            answer = hydrostatic_transport_source(density, pullback=pullback).vector()
            errors = []
            for h in [2e-4, 1e-4]:
                fd = (
                    density.transport(1 + h, pullback).vector()
                    - density.transport(1 - h, pullback).vector()
                ) / (2 * h)
                errors.append(np.linalg.norm(fd - answer))
            self.assertLess(errors[-1], 1e-7)
            self.assertGreater(errors[0] / errors[1], 3.8)


if __name__ == "__main__":
    unittest.main()
