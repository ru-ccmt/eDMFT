"""Independent nonlinear geometry, native-coordinate and transport controls."""

from copy import deepcopy
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace as NS
import unittest
import numpy as np
from scipy.interpolate import CubicSpline
from lapw_hartree_periodic import (
    NativeStarTable,
    NativeHartreeReanalysis,
    PeriodicHartreeResponse,
)
from lapw_strain_potential import (
    PotentialFile,
    read_potential_file,
    NativePotentialFileLayout,
    unequal_centered_derivative,
    hydrostatic_mt_transport,
    strain_mask_derivative,
    masked_potential_strain_source,
)
from lapw_hydrostatic_potential import hydrostatic_coulomb_geometry
from test_lapw_hartree_periodic import fixture
from wien_response_density import Density


class ExplicitStrainPotentialTests(unittest.TestCase):
    def setUp(self):
        self.d, self.s, self.reps = fixture()
        self.h = PeriodicHartreeResponse(self.d, self.s, (12, 14, 16), inversion=False)
        self.b = np.diag(2 * np.pi / np.array([self.s.a, self.s.b, self.s.c]))

    def test_file_coordinate_selection_conjugation_and_labels(self):
        h = self.h
        sph = PotentialFile({(0, 0, 0): np.arange(301.0)}, {})
        ns = PotentialFile(
            {
                key: np.full(301, i + 0.5)
                for i, key in enumerate(h.keys)
                if key[1:] != (0, 0)
            },
            {tuple(g): complex(i, -i / 3) for i, g in enumerate(self.reps)},
        )
        layout = NativePotentialFileLayout(h, sph, ns, complex_interstitial=True)
        field = np.zeros(h.size)
        for key, sl in zip(h.keys, h.slices):
            field[sl] = (sph.radial | ns.radial)[key]
        # Internal REAN4 is the conjugate of VNS file data.
        field[h.pw_offset :] = (
            np.array(list(ns.interstitial.values())).conj().view(float)
        )
        np.testing.assert_array_equal(layout.selector.matvec(field), layout.central)
        # Dict insertion/file-row order must not replace explicit G labels.
        reversed_ns = PotentialFile(
            dict(reversed(list(ns.radial.items()))),
            dict(reversed(list(ns.interstitial.items()))),
        )
        np.testing.assert_array_equal(layout.pack(sph, reversed_ns), layout.central)
        native = np.arange(layout.coordinates.size, dtype=float)
        np.testing.assert_array_equal(
            layout.covector(layout.selector.rmatvec(native)), native
        )
        bad = PotentialFile(
            ns.radial, {g: v for g, v in ns.interstitial.items() if g != (0, 0, 0)}
        )
        with self.assertRaises(ValueError):
            layout.pack(sph, bad)

    def test_parser_rejects_truncation_and_duplicate_channels(self):
        header = "   ATOMNUMBER = 1\n VLM(R) FOR L 0 M= 0\n"
        values = "   " + "".join(f"{v:19.12E}" for v in [1.0, 2.0, 3.0, 4.0]) + "\n"
        pw = (
            " 1 NUMBER OF PW\n"
            + "   "
            + "".join(f"{v:5d}" for v in [0, 0, 0])
            + f"{2.:19.12E}{-3.:19.12E}\n"
        )
        with TemporaryDirectory() as temp:
            p = Path(temp) / "potential"
            p.write_text(header + values + pw)
            x = read_potential_file(p, [4])
            np.testing.assert_array_equal(x.radial[0, 0, 0], [1, 2, 3, 4])
            self.assertEqual(x.interstitial[0, 0, 0], 2 - 3j)
            p.write_text(header + values + header + values)
            with self.assertRaises(ValueError):
                read_potential_file(p, [4])
            p.write_text(header + values + pw)
            with self.assertRaises(ValueError):
                read_potential_file(p, [5])

    def test_unequal_steps_are_exact_for_quadratic(self):
        a = np.array([2.0, 3.0, 4.0])
        b = np.array([1.7, -0.4, 2.0])
        c = np.array([3.0, -9.0, 0.1])
        hm, hp = -0.00021, 0.00019
        actual = unequal_centered_derivative(
            a + b * hm + c * hm**2, a, a + b * hp + c * hp**2, hm, hp
        )
        np.testing.assert_allclose(actual, b, atol=2e-12, rtol=0)
        with self.assertRaises(ValueError):
            unequal_centered_derivative(a, a, a, 0.1, 0.2)

    def test_potential_material_transport_mixed_radial_units(self):
        x = np.linspace(-4, 1, 101)
        radial = {(0, 0, 0): -8 + 3 * x + x**2, (0, 2, 0): 1 + 2 * x + x**3 / 10}
        expected = hydrostatic_mt_transport(radial, [x])
        step = 1e-5
        for key, v in radial.items():
            endpoints = []
            for sign in [-1, 1]:
                lam = 1 + sign * step
                endpoints.append(
                    CubicSpline(x, v)(x + np.log(lam))
                    / (lam if key[1:] == (0, 0) else 1)
                )
            np.testing.assert_allclose(
                (endpoints[1] - endpoints[0]) / (2 * step),
                expected[key],
                rtol=2e-9,
                atol=2e-9,
            )
        # The point nucleus is constant in rV at fixed r, but is not constant
        # when pulled back to material coordinates. No zero-at-RMT assumption.
        np.testing.assert_allclose(
            hydrostatic_mt_transport({(0, 0, 0): np.full(101, -8.0)}, [x])[0, 0, 0], 8.0
        )

    def reanalysis(self, t, e):
        f = np.eye(3) + t * e
        return NativeHartreeReanalysis(
            self.h.stars,
            self.h.reanalysis.shape,
            np.linalg.inv(f).T @ self.b,
            self.s.pos,
            self.s.rmt,
            self.h.volume * np.linalg.det(f),
        )

    def test_mask_hydrostatic_and_shear_against_full_geometry(self):
        for e in [
            np.eye(3),
            np.array([[0.2, 0.3, 0], [0.3, -0.1, 0.1], [0, 0.1, 0.4]]),
        ]:
            exact = strain_mask_derivative(
                self.h.reanalysis, self.b, self.s.pos, self.s.rmt, self.h.volume, e
            )
            errors = []
            for step in [4e-4, 2e-4]:
                fd = (
                    self.reanalysis(step, e).mask_fourier
                    - self.reanalysis(-step, e).mask_fourier
                ) / (2 * step)
                errors.append(np.linalg.norm(fd - exact) / np.linalg.norm(exact))
            self.assertLess(errors[-1], 2e-6)
            self.assertGreater(errors[0] / errors[1], 3.8)
        expected = 4 * np.pi * self.s.rmt[0] ** 3 / self.h.volume
        exact = strain_mask_derivative(
            self.h.reanalysis, self.b, self.s.pos, self.s.rmt, self.h.volume, np.eye(3)
        )
        self.assertAlmostEqual(exact[0, 0, 0].real, expected, places=15)

    def test_reanalysis_surface_source_for_nonconstant_potential(self):
        rng = np.random.default_rng(928)
        potential = rng.normal(size=len(self.reps)) + 1j * rng.normal(
            size=len(self.reps)
        )
        e = np.array([[1, 0.1, 0], [0.1, 0.7, 0], [0, 0, -0.3]])
        mask = strain_mask_derivative(
            self.h.reanalysis, self.b, self.s.pos, self.s.rmt, self.h.volume, e
        )
        actual = masked_potential_strain_source(self.h.reanalysis, potential, mask)
        step = 1e-5
        expected = (
            self.reanalysis(step, e).action(potential)
            - self.reanalysis(-step, e).action(potential)
        ) / (2 * step)
        np.testing.assert_allclose(actual, expected, rtol=2e-8, atol=2e-10)

    def test_coulomb_geometry_includes_nucleus_and_boundary(self):
        rng = np.random.default_rng(172)
        d = Density.__new__(Density)
        d.__dict__.update(self.d.__dict__)
        d.mult = [1]
        d.rad = {
            key: rng.normal(size=301) * np.exp(2 * self.d.grids[0])
            for key in self.h.keys
        }
        d.pw = {tuple(g): rng.normal() + 1j * rng.normal() for g in self.reps}
        s = deepcopy(self.s)
        s.Znuc = np.array([3.0])
        central = NS(hartree=self.h, structure=s, input_density=d)
        source, _, _ = hydrostatic_coulomb_geometry(central)

        # Nonlinear reference: rebuild native Poisson matrices at endpoints;
        # add the point-nucleus multipole and local Dirichlet potential.
        def value(t):
            ss = deepcopy(s)
            ss.a *= 1 + t
            ss.b *= 1 + t
            ss.c *= 1 + t
            h = PeriodicHartreeResponse(d, ss, (12, 14, 16), inversion=False)
            z = np.zeros(len(h.keys))
            z[0] = 3 / np.sqrt(4 * np.pi)
            ion_pw = h.A @ (h.S @ z)
            bound = (h.B @ ion_pw).real
            ion_mt = [op.boundary * bb for op, bb in zip(h.radial, bound)]
            r = h.integrations[0].r
            ion_mt[0] -= 6 * (1 - r / ss.rmt[0])
            return h.action(d.vector()) + np.concatenate(
                ion_mt + [h.reanalysis.action(ion_pw).view(float)]
            )

        errors = []
        for step in [4e-4, 2e-4]:
            fd = (value(step) - value(-step)) / (2 * step)
            errors.append(np.linalg.norm(fd - source) / np.linalg.norm(source))
        self.assertLess(errors[-1], 2e-6)
        self.assertGreater(errors[0] / errors[1], 3.8)


if __name__ == "__main__":
    unittest.main()
