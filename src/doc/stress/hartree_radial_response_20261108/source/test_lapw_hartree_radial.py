"""Discrete native Coulomb quadrature, physical units, and paired reverses.

Native controls use EDMFT_HARTREE_REFERENCE_LIBRARY, built from the user's
WIEN2k CoulInt/Spline and the small checkpoint reference wrapper. They do not
ship or substitute a second native implementation in this repository.
"""

import ctypes
import os
from pathlib import Path
import unittest
import numpy as np
from lapw_hartree_radial import NativeCoulombIntegral, RadialHartreeDirichlet


class RadialHartreeTests(unittest.TestCase):
    def test_integral_of_cubic_in_log_coordinate(self):
        r = np.geomspace(0.02, 2.0, 71)
        q = NativeCoulombIntegral(r)
        x = np.arange(1, len(r) + 1) * q.dx
        f = 1 + 0.4 * x - 0.3 * x * x + 0.2 * x**3
        primitive = x + 0.2 * x * x - 0.1 * x**3 + 0.05 * x**4
        np.testing.assert_allclose(
            q.action(f / r), primitive - primitive[0], atol=2e-13
        )
        np.testing.assert_allclose(
            q.action(f / r, tail=True), primitive[-1] - primitive, atol=2e-13
        )

    def test_batched_prefix_and_tail_transposes(self):
        rng = np.random.default_rng(871)
        for n in [4, 5, 97, 781]:
            q = NativeCoulombIntegral(np.geomspace(1e-4, 2.2, n))
            x, y = rng.normal(size=(2, n, 3))
            for tail in [False, True]:
                ax, aty = q.action(x, tail=tail), q.pullback(y, tail=tail)
                scale = max(
                    np.linalg.norm(y) * np.linalg.norm(ax),
                    np.linalg.norm(x) * np.linalg.norm(aty),
                )
                self.assertLess(abs(np.vdot(y, ax) - np.vdot(aty, x)) / scale, 2e-14)
            prefix = q.action(x)
            np.testing.assert_allclose(
                q.action(x, tail=True), prefix[-1] - prefix, atol=5e-14
            )

    def test_dirichlet_units_boundary_and_transpose(self):
        r = np.geomspace(1e-4, 2.0, 401)
        integration = NativeCoulombIntegral(r)
        for l in [0, 1, 2, 4, 8]:
            q = RadialHartreeDirichlet(integration, l)
            self.assertLess(q.operator.check_transpose(), 2e-13)
            boundary = 0.31
            c = r ** (l + 2)
            self.assertAlmostEqual(
                q.action(c, boundary)[-1], boundary * (r[-1] if l == 0 else 1)
            )
            profile = boundary * r if l == 0 else boundary * (r / r[-1]) ** l
            np.testing.assert_allclose(
                q.action(np.zeros_like(r), boundary), profile, atol=1e-15
            )

    def test_uniform_multipole_continuum_convergence(self):
        # Independent Poisson solution for rho_L(r)=r^l. For l=0 the
        # input is CLMSUM c00=r^2, and the output is r*V_spherical.
        errors = {}
        for n in [401, 801]:
            r = np.geomspace(1e-5, 2.0, n)
            integration = NativeCoulombIntegral(r)
            for l in [0, 1, 2, 4, 8]:
                q = RadialHartreeDirichlet(integration, l)
                exact = (
                    r * (r[-1] ** 2 - r * r) / 3
                    if l == 0
                    else 4 * np.pi / (2 * l + 3) * r**l * (r[-1] ** 2 - r * r)
                )
                errors[n, l] = np.linalg.norm(
                    q.action(r ** (l + 2), 0) - exact
                ) / np.linalg.norm(exact)
        for l in [0, 1, 2, 4, 8]:
            self.assertGreater(errors[401, l] / errors[801, l], 16)
            self.assertLess(errors[801, l], 3e-5)

    def test_grid_and_real_input_validation(self):
        for r in [
            np.ones(10),
            np.arange(10),
            np.geomspace(0.01, 1, 3),
            np.geomspace(0.01, 1, 10) ** np.linspace(1, 2, 10),
        ]:
            with self.assertRaises(ValueError):
                NativeCoulombIntegral(r)
        q = NativeCoulombIntegral(np.geomspace(0.01, 1, 51))
        with self.assertRaises(ValueError):
            q.action(np.ones(51, dtype=complex))
        with self.assertRaises(ValueError):
            q.pullback(np.ones(50))
        with self.assertRaises(ValueError):
            RadialHartreeDirichlet(q, -1)
        with self.assertRaises(ValueError):
            RadialHartreeDirichlet(q, 0).action(np.ones(51), np.nan)


LIB = os.environ.get("EDMFT_HARTREE_REFERENCE_LIBRARY", "")


@unittest.skipUnless(
    LIB and Path(LIB).is_file(),
    "set EDMFT_HARTREE_REFERENCE_LIBRARY for native controls",
)
class NativeHartreeTests(unittest.TestCase):
    def test_integrator_against_compiled_native_coulint(self):
        lib = ctypes.CDLL(LIB)
        ptr = ctypes.POINTER(ctypes.c_double)
        for n in [4, 5, 101, 781]:
            r = np.geomspace(1e-4, 2.2, n)
            q = NativeCoulombIntegral(r)
            for c in [r * r * np.sin(r), np.random.default_rng(n).normal(size=n)]:
                expected = c.copy()
                lib.coulint_(
                    r.ctypes.data_as(ptr),
                    ctypes.byref(ctypes.c_double(q.dx)),
                    ctypes.byref(ctypes.c_int(n)),
                    expected.ctypes.data_as(ptr),
                )
                np.testing.assert_allclose(
                    q.action(c), expected, atol=5e-14, rtol=2e-13
                )

    def test_dirichlet_against_native_radial_assembly(self):
        lib = ctypes.CDLL(LIB)
        ptr = ctypes.POINTER(ctypes.c_double)
        for n in [101, 401, 781]:
            r = np.geomspace(1e-4, 2.2, n)
            integration = NativeCoulombIntegral(r)
            for l in [0, 1, 2, 4, 8]:
                q = RadialHartreeDirichlet(integration, l)
                c = r ** (l + 2) * (1 + 0.1 * np.cos(r))
                boundary = 0.314
                expected = np.empty(n)
                lib.radial_dirichlet_ref_(
                    ctypes.byref(ctypes.c_int(n)),
                    r.ctypes.data_as(ptr),
                    ctypes.byref(ctypes.c_double(integration.dx)),
                    ctypes.byref(ctypes.c_int(l)),
                    c.ctypes.data_as(ptr),
                    ctypes.byref(ctypes.c_double(boundary)),
                    expected.ctypes.data_as(ptr),
                )
                error = np.linalg.norm(
                    q.action(c, boundary) - expected
                ) / np.linalg.norm(expected)
                self.assertLess(error, 5e-13)


if __name__ == "__main__":
    unittest.main()
