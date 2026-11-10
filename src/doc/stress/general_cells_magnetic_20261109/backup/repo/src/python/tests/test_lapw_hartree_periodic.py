"""Independent Poisson controls, native comparisons and complete transposes."""

import ctypes
import os
from pathlib import Path
from types import SimpleNamespace
import unittest
import numpy as np
from scipy.special import sph_harm_y
from lapw_hartree_periodic import (
    native_charge_weights,
    native_local_axes,
    NativeStarTable,
    NativeHartreeReanalysis,
    PeriodicHartreeResponse,
)


def fixture():
    """A P1 cell with all real harmonics through l=2 and complex PW data."""
    r = np.geomspace(1e-5, 1.3, 301)
    keys = [(0, 0, 0)]
    for l in [1, 2]:
        keys += [(0, l, 0)]
        for m in range(1, l + 1):
            keys += [(0, l, m), (0, -l, m)]
    rep = np.array(
        [
            (0, 0, 0),
            (1, 0, 0),
            (-1, 0, 0),
            (0, 1, 0),
            (0, -1, 0),
            (1, 1, 0),
            (-1, -1, 0),
        ]
    )

    def density():
        return SimpleNamespace(
            rad={k: np.zeros(len(r)) for k in keys},
            pw={tuple(g): 0j for g in rep},
            grids=[np.log(r)],
            zeros=density,
        )

    s = SimpleNamespace(
        nat=1,
        lattyp="P",
        a=6.0,
        b=7.0,
        c=8.0,
        alpha=90.0,
        beta=90.0,
        gamma=90.0,
        iatnr=[-1],
        jrj=[len(r)],
        r0=[r[0]],
        rmt=[r[-1]],
        mult=[1],
        pos=np.array([[0.17, 0.31, 0.41]]),
        rotloc=np.eye(3)[None],
        timat=np.eye(3, dtype=int)[None],
        tau=np.zeros((1, 3)),
    )
    return density(), s, rep


class PeriodicHartreeTests(unittest.TestCase):
    def test_charge_continuum_convergence_and_no_origin(self):
        errors = []
        for n in [101, 201, 401]:
            r = np.geomspace(1e-5, 2.0, n)
            w = native_charge_weights(r)
            # int r^2 dr, with the native lower limit r0, not zero.
            errors.append(abs(w @ r**2 - (r[-1] ** 3 - r[0] ** 3) / 3))
        self.assertGreater(errors[0] / errors[1], 10)
        self.assertGreater(errors[1] / errors[2], 10)

    def test_native_rounded_axis_repair(self):
        axes = np.array(
            [[0.7071068, 0.7071068, 0], [-0.7071068, 0.7071068, 0], [0, 0, 1]]
        )
        original = axes.copy()
        actual = native_local_axes(axes)
        np.testing.assert_array_equal(axes, original)
        np.testing.assert_allclose(actual.T @ actual, np.eye(3), atol=3e-16)
        arbitrary = np.eye(3) * 0.998
        np.testing.assert_array_equal(native_local_axes(arbitrary), arbitrary)

    def test_nonsymmorphic_star_average(self):
        reps = np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0]])
        # Half-cell translation extinguishes the odd star; repeated rotations
        # must average phases, not create two copies of the same G member.
        stars = NativeStarTable(
            reps,
            np.repeat(np.eye(3, dtype=int)[None], 2, axis=0),
            [[0, 0, 0], [0.5, 0, 0]],
            inversion=False,
        )
        np.testing.assert_array_equal(stars.sizes, [1, 1, 1])
        np.testing.assert_allclose(stars.phases, [1, 0, 1], atol=3e-16)

    def test_fft_constant_potential_mask_and_transpose(self):
        d, s, rep = fixture()
        b = np.diag(2 * np.pi / np.array([s.a, s.b, s.c]))
        stars = NativeStarTable(rep, s.timat, s.tau, inversion=False)
        q = NativeHartreeReanalysis(
            stars, (12, 14, 16), b, s.pos, s.rmt, s.a * s.b * s.c
        )
        v = np.zeros(len(rep), complex)
        v[0] = 0.43
        masked = q.action(v)
        # Constant potential times theta_I has analytical Fourier coefficients.
        expected = 0.43 * q.mask_fourier.conj()[q.indices]
        np.testing.assert_allclose(masked, expected, atol=2e-16)
        self.assertLess(q.operator.check_transpose(samples=3), 1e-14)

    def test_real_harmonic_multipole_conversion(self):
        d, s, rep = fixture()
        q = PeriodicHartreeResponse(d, s, (12, 14, 16), inversion=False)
        # A physical angular density assembled in either real or complex
        # harmonics must agree. This independently checks both sine signs.
        rng = np.random.default_rng(31)
        moments = rng.normal(size=len(q.keys))
        complex_m = q.S @ moments
        lm = sorted(set((abs(l), mm) for a, l, m in q.keys for mm in {m, -m}))
        theta = np.array([0.2, 1.4, 2.1])
        phi = np.array([0.5, -0.7, 2.7])
        complex_field = sum(
            c * sph_harm_y(l, m, theta, phi) for c, (l, m) in zip(complex_m, lm)
        )
        real_field = np.zeros(3, complex)
        for c, (a, L, m) in zip(moments, q.keys):
            l = abs(L)
            yp = sph_harm_y(l, m, theta, phi)
            ym = sph_harm_y(l, -m, theta, phi)
            h = (
                yp
                if m == 0
                else (
                    ((-1) ** m * yp + ym) / np.sqrt(2)
                    if L > 0
                    else -1j * ((-1) ** m * yp - ym) / np.sqrt(2)
                )
            )
            real_field += c * h
        np.testing.assert_allclose(complex_field, real_field, atol=1e-15)
        np.testing.assert_allclose(complex_field.imag, 0, atol=1e-15)

    def test_complete_transpose_linearity_and_boundary(self):
        d, s, rep = fixture()
        q = PeriodicHartreeResponse(d, s, (12, 14, 16), inversion=False)
        self.assertLess(q.operator.check_transpose(samples=3), 1e-13)
        rng = np.random.default_rng(71)
        x = rng.normal(size=q.size)
        y = rng.normal(size=q.size)
        np.testing.assert_allclose(
            q.action(2 * x - 0.3 * y),
            2 * q.action(x) - 0.3 * q.action(y),
            atol=2e-12,
            rtol=1e-12,
        )
        np.testing.assert_array_equal(q.action(np.zeros(q.size)), np.zeros(q.size))
        out = q.action(x)
        _, pot, b = q.intermediates(x)
        self.assertEqual(pot[0], 0)
        for k, sl, bd in zip(q.keys, q.slices, b):
            self.assertAlmostEqual(
                out[sl][-1], bd * (s.rmt[k[0]] if k[1] == 0 else 1), places=12
            )

    def test_native_potential_coordinate_selection(self):
        d, s, rep = fixture()
        q = PeriodicHartreeResponse(d, s, (12, 14, 16), inversion=False)
        for complex_values in [True, False]:
            ns = {0: [(-2, 1), (1, 0)]}
            selector, coords = q.native_potential_output(
                [0], ns, rep[:3], complex_interstitial=complex_values
            )
            self.assertLess(selector.check_transpose(), 1e-14)
            rng = np.random.default_rng(29)
            x = rng.normal(size=q.size)
            spherical, nonspherical, pw = coords.unpack(selector.action(x))
            lookup = dict(zip(q.keys, q.slices))
            np.testing.assert_array_equal(spherical, x[lookup[0, 0, 0]])
            np.testing.assert_array_equal(
                nonspherical[0], np.column_stack([x[lookup[0, *lm]] for lm in ns[0]])
            )
            expected = x[q.pw_offset :].view(complex)[:3]
            np.testing.assert_array_equal(
                pw, expected.conj() if complex_values else expected.real
            )
            composed = selector @ q.operator
            self.assertLess(composed.check_transpose(), 1e-13)

    def test_complex_vns_file_reconstructs_same_real_potential(self):
        # SETFFT/forward FFT uses exp(-i G.r); LAPW1's FILE coefficients
        # use exp(+i G.r). A sine term detects the required conjugation;
        # inversion-symmetric cosine-only tests cannot detect this error.
        d, s, rep = fixture()
        q = PeriodicHartreeResponse(d, s, (12, 14, 16), inversion=False)
        selector, coords = q.native_potential_output(
            [0], {0: [(1, 0)]}, rep[:3], complex_interstitial=True
        )
        x = np.zeros(q.size)
        z = np.array([0.4, 0.7 + 0.3j, 0.7 - 0.3j])
        x[q.pw_offset : q.pw_offset + 6] = z.view(float)
        _, _, file_coeff = coords.unpack(selector.matvec(x))
        theta = np.linspace(0, 2 * np.pi, 31)
        wave = np.array([np.ones(len(theta)), np.exp(1j * theta), np.exp(-1j * theta)])
        physical = (file_coeff @ wave).real
        np.testing.assert_allclose(
            physical, 0.4 + 1.4 * np.cos(theta) + 0.6 * np.sin(theta), atol=1e-14
        )

    def test_guards_reject_unsupported_or_inconsistent_inputs(self):
        d, s, rep = fixture()
        with self.assertRaises(ValueError):
            NativeStarTable(rep[[1, 0]], s.timat, s.tau, inversion=False)
        with self.assertRaises(ValueError):
            NativeStarTable(
                np.concatenate([rep, rep[:1]]), s.timat, s.tau, inversion=False
            )
        with self.assertRaises(ValueError):
            PeriodicHartreeResponse(d, s, (2, 14, 16), inversion=False)
        with self.assertRaises(ValueError):
            PeriodicHartreeResponse(d, s, (12, 14, 16), ncon=1, inversion=False)
        s.lattyp = "F"
        with self.assertRaises(NotImplementedError):
            PeriodicHartreeResponse(d, s, (12, 14, 16), inversion=False)
        s.lattyp = "P"
        s.iatnr = [1]
        with self.assertRaises(NotImplementedError):
            PeriodicHartreeResponse(d, s, (12, 14, 16), inversion=False)
        s.iatnr = [-1]
        s.rmt = [1.4]
        with self.assertRaises(ValueError):
            PeriodicHartreeResponse(d, s, (12, 14, 16), inversion=False)


CHARGE = os.environ.get("EDMFT_HARTREE_CHARGE_LIBRARY", "")


@unittest.skipUnless(
    CHARGE and Path(CHARGE).is_file(),
    "set EDMFT_HARTREE_CHARGE_LIBRARY for native CHARGE control",
)
class NativeChargeTests(unittest.TestCase):
    def test_default_charge_against_compiled_native(self):
        lib = ctypes.CDLL(CHARGE)
        ptr = ctypes.POINTER(ctypes.c_double)
        for n in [4, 5, 6, 7, 8, 9, 101, 781]:
            r = np.geomspace(1e-5, 2.1, n)
            dx = np.log(r[-1] / r[0]) / (n - 1)
            for field in [r**2, np.random.default_rng(n).normal(size=n)]:
                out = ctypes.c_double()
                lib.charge_(
                    r.ctypes.data_as(ptr),
                    ctypes.byref(ctypes.c_double(dx)),
                    field.ctypes.data_as(ptr),
                    ctypes.byref(ctypes.c_int(1)),
                    ctypes.byref(ctypes.c_int(n)),
                    ctypes.byref(out),
                )
                self.assertAlmostEqual(
                    native_charge_weights(r) @ field, out.value, places=12
                )


if __name__ == "__main__":
    unittest.main()
