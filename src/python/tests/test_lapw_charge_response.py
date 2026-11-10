"""Native MIXER reference, charge closure, finite differences and adjoints."""

import ctypes, os, unittest
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from wien_response_density import Density
from lapw_charge_response import (
    native_charge_covectors,
    NativeChargeNormalization,
    native_mixer_inversion,
)


def fixture():
    d = Density.__new__(Density)
    r = np.geomspace(1e-4, 1.3, 201)
    d.grids = [np.log(r)]
    d.mult = [2]
    d.rad = {(0, 0, 0): 2 * r * r * np.exp(-r), (0, 2, 0): 0.01 * r**4 * np.exp(-r)}
    d.pw = {
        (0, 0, 0): 0.03 + 0j,
        (1, 0, 0): 0.002 + 0.001j,
        (-1, 0, 0): 0.002 - 0.001j,
        (0, 1, 1): 0.001 - 0.002j,
        (0, -1, -1): 0.001 + 0.002j,
    }
    s = SimpleNamespace(
        lattyp="P",
        alpha=90.0,
        beta=90.0,
        gamma=90.0,
        a=5.0,
        b=6.0,
        c=7.0,
        nat=1,
        mult=np.array([2]),
        jrj=np.array([len(r)]),
        r0=np.array([r[0]]),
        rmt=np.array([r[-1]]),
        pos=np.array([[0.1, 0.2, 0.3], [0.6, 0.7, 0.8]]),
    )
    return d, s


def native_endpoint(d, s, target, inversion):
    lib = ctypes.CDLL(os.environ["EDMFT_MIXER_NORMALIZATION_LIBRARY"])
    nrad = max(s.jrj)
    r = np.zeros((nrad, s.nat), order="F")
    rad = r.copy(order="F")
    for a in range(s.nat):
        r[: s.jrj[a], a] = np.exp(d.grids[a])
        rad[: s.jrj[a], a] = d.rad[a, 0, 0]
    pos = np.zeros((3, s.nat * 48), order="F")
    pos[:, : len(s.pos)] = s.pos.T
    arrays = [
        r,
        rad,
        np.array(s.jrj, np.int32),
        np.array(s.mult, np.int32),
        np.array(s.rmt, float),
        pos,
        np.array([s.a, s.b, s.c]),
        np.array(list(d.pw), np.int32).T.copy(order="F"),
        np.array(list(d.pw.values()), complex),
        np.array([target]),
        np.zeros(s.nat + 1),
    ]
    pointers = [a.ctypes.data_as(ctypes.c_void_p) for a in arrays]
    lib.native_normalization(
        ctypes.c_int(nrad),
        ctypes.c_int(s.nat),
        ctypes.c_int(len(d.pw)),
        *pointers[:10],
        ctypes.c_int(inversion),
        pointers[10]
    )
    out = d.copy()
    out.pw = dict(zip(out.pw, arrays[8]))
    return out, arrays[-1]


class ChargeResponseTests(unittest.TestCase):
    def test_native_reality_requires_inversion_at_origin(self):
        s = SimpleNamespace(
            timat=np.array([np.eye(3, dtype=int), -np.eye(3, dtype=int)]),
            tau=np.zeros((2, 3)),
        )
        self.assertTrue(native_mixer_inversion(s))
        s.tau[1, 0] = 1e-4
        self.assertFalse(native_mixer_inversion(s))
        s.tau[1, 0] = 1e-6
        self.assertTrue(native_mixer_inversion(s))
        s.tau = s.tau[:1]
        with self.assertRaises(ValueError):
            native_mixer_inversion(s)

    def test_fixed_rmt_geometry_source_matches_changing_cell(self):
        from copy import deepcopy

        d, s = fixture()
        op = NativeChargeNormalization(d, s, 7.8, inversion=False)
        h = 1e-5
        sp = deepcopy(s)
        sm = deepcopy(s)
        for name in ["a", "b", "c"]:
            setattr(sp, name, getattr(s, name) * (1 + h))
            setattr(sm, name, getattr(s, name) * (1 - h))
        mp, ip = native_charge_covectors(d, sp, inversion=False)
        mm, im = native_charge_covectors(d, sm, inversion=False)
        explicit = op.explicit_source(
            mt_geometry=((mp - mm) / (2 * h)) @ op.base,
            ir_geometry=((ip - im) / (2 * h)) @ op.base,
        )
        fd = (
            NativeChargeNormalization(d, sp, 7.8, inversion=False).central
            - NativeChargeNormalization(d, sm, 7.8, inversion=False).central
        ) / (2 * h)
        np.testing.assert_allclose(explicit, fd, rtol=2e-8, atol=1e-10)

    def test_constant_density_integrates_to_cell_volume(self):
        errors = []
        for n in [201, 401, 801]:
            d, s = fixture()
            r = np.geomspace(s.r0[0], s.rmt[0], n)
            s.jrj[0] = n
            d.grids = [np.log(r)]
            d.rad = {(0, 0, 0): 4 * np.pi * r * r * 0.03}
            for g in d.pw:
                d.pw[g] = 0.03 if g == (0, 0, 0) else 0j
            mt, ir = native_charge_covectors(d, s, inversion=False)
            errors.append(abs((mt + ir) @ d.vector() - 0.03 * s.a * s.b * s.c))
        # CHARGE is a finite quadrature of r^3 on log(r), not an exact
        # continuum integral. Verify its fourth-order approach to volume.
        self.assertGreater(errors[0] / errors[1], 14)
        self.assertGreater(errors[1] / errors[2], 14)
        self.assertLess(errors[2], 1e-8)

    def test_charge_closure_and_transpose_away_from_fixed_point(self):
        d, s = fixture()
        op = NativeChargeNormalization(d, s, 7.8, inversion=False)
        self.assertGreater(abs(op.alpha - 1), 0.1)
        self.assertLess(op.operator.check_transpose(samples=3), 1e-13)
        rng = np.random.default_rng(12)
        x = rng.normal(size=op.size) * 0.01
        self.assertAlmostEqual((op.mt + op.ir) @ op.action(x), 0, places=12)
        self.assertAlmostEqual((op.mt + op.ir) @ op.central, 7.8, places=12)
        h = 1e-4
        fd = (op.endpoint(op.base + h * x) - op.endpoint(op.base - h * x)) / (2 * h)
        np.testing.assert_allclose(op.action(x), fd, rtol=1e-8, atol=1e-10)

    def test_second_normalization_is_redundant_after_first(self):
        d, s = fixture()
        a = NativeChargeNormalization(d, s, 7.8, inversion=False)
        b = NativeChargeNormalization(d.from_vector(a.central), s, 7.8, inversion=False)
        x = np.random.default_rng(3).normal(size=a.size)
        np.testing.assert_allclose(
            b.action(a.action(x)), a.action(x), atol=1e-13, rtol=1e-13
        )

    def test_disabled_policy_and_target_source(self):
        d, s = fixture()
        op = NativeChargeNormalization(d, s, 7.8, inversion=False, enabled=False)
        x = np.random.default_rng(2).normal(size=op.size)
        np.testing.assert_array_equal(op.action(x), x)
        np.testing.assert_array_equal(op.pullback(x), x)
        np.testing.assert_array_equal(op.explicit_source(target_derivative=1), x * 0)
        active = NativeChargeNormalization(d, s, 7.8, inversion=False)
        y = active.explicit_source(target_derivative=0.2)
        self.assertAlmostEqual((active.mt + active.ir) @ y, 0.2, places=14)

    def test_guards(self):
        d, s = fixture()
        with self.assertRaises(ValueError):
            NativeChargeNormalization(d, s, float("nan"), inversion=False)
        with self.assertRaises(ValueError):
            NativeChargeNormalization(d, s, 7.8, inversion="false")
        with self.assertRaises(ValueError):
            NativeChargeNormalization(d, s, -1, inversion=False)
        s.a = 0.0
        with self.assertRaises(ValueError):
            native_charge_covectors(d, s, inversion=False)

    @unittest.skipUnless(
        os.environ.get("EDMFT_MIXER_NORMALIZATION_LIBRARY"),
        "native mixer reference required",
    )
    def test_native_mixer_endpoint_and_directional_derivative(self):
        d, s = fixture()
        rng = np.random.default_rng(49)
        x = rng.normal(size=len(d.vector())) * 0.001
        for inv in [False, True]:
            op = NativeChargeNormalization(d, s, 7.8, inversion=inv)
            native, q = native_endpoint(d, s, 7.8, inv)
            self.assertAlmostEqual(op.M, q[:-1] @ s.mult, places=12)
            self.assertAlmostEqual(op.I, q[-1], places=12)
            np.testing.assert_allclose(
                op.central, native.vector(), rtol=5e-14, atol=1e-14
            )
            h = 1e-4
            p, _ = native_endpoint(d.from_vector(d.vector() + h * x), s, 7.8, inv)
            m, _ = native_endpoint(d.from_vector(d.vector() - h * x), s, 7.8, inv)
            np.testing.assert_allclose(
                op.action(x), (p.vector() - m.vector()) / (2 * h), rtol=2e-8, atol=2e-11
            )
