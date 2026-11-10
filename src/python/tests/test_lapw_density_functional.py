"""Independent scalar finite differences including ionic-boundary partners."""

from pathlib import Path
from copy import deepcopy
import tempfile, unittest
import numpy as np
from lapw_density_functional import PreparedDensityFunctional
from lapw_potential_core_state import save_potential_core_state
from test_lapw_potential_core_state import inputs


class DensityFunctionalTests(unittest.TestCase):
    def prepare(self, root, args):
        return PreparedDensityFunctional(
            save_potential_core_state(root / "state.npz", **args)
        )

    def test_density_gradient_of_complete_discrete_scalar(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            args = inputs()
            args["density"].pw.update({(1, 0, 0): 0.0001 + 0j, (-1, 0, 0): 0.0001 + 0j})
            a = self.prepare(root, args)
            r = np.exp(args["density"].grids[0])
            d = args["density"].zeros()
            d.rad[0, 0, 0] = r * r * np.sin(r) * 0.01
            d.pw[0, 0, 0] = 0.003 + 0j
            d.pw[1, 0, 0] = 0.001 + 0j
            d.pw[-1, 0, 0] = 0.001 + 0j
            h = 2e-4
            values = []
            for sign in [-1, 1]:
                b = deepcopy(args)
                b["density"] = args["density"].from_vector(
                    args["density"].vector() + sign * h * d.vector()
                )
                values.append(self.prepare(root, b).value_ry)
            fd = (values[1] - values[0]) / (2 * h)
            np.testing.assert_allclose(a.action(d.vector()), fd, rtol=3e-9, atol=2e-9)
            # Omitting the nuclear boundary covector would fail this test.
            missing = a.state.hartree.pullback(a.boundary_gradient) @ d.vector()
            self.assertGreater(abs(missing), 1e-4)

    def test_explicit_geometry_with_fixed_radial_density(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            args = inputs()
            a = self.prepare(root, args)
            e = np.diag([0.7, -0.4, 0.9])
            h = 2e-5
            values = []
            for sign in [-1, 1]:
                b = deepcopy(args)
                for name, factor in zip(["a", "b", "c"], np.diag(e)):
                    setattr(
                        b["structure"],
                        name,
                        getattr(args["structure"], name) * (1 + sign * h * factor),
                    )
                values.append(self.prepare(root, b).value_ry)
            np.testing.assert_allclose(
                a.geometry(e)[0],
                (values[1] - values[0]) / (2 * h),
                rtol=2e-8,
                atol=2e-8,
            )

    def test_reject_nonfinite_density_direction(self):
        with tempfile.TemporaryDirectory() as tmp:
            a = self.prepare(Path(tmp), inputs())
            d = np.zeros_like(a.density)
            d[0] = np.nan
            with self.assertRaises(ValueError):
                a.action(d)


if __name__ == "__main__":
    unittest.main()
