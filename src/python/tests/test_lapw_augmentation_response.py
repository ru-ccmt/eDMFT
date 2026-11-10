"""Angular response controls with independent coefficient and Dyson endpoints."""

import unittest
import numpy as np
from scipy.special import spherical_jn
from lapw_augmentation_response import (
    PreparedAugmentationChannel,
    AugmentationContractions,
    RadialAngularLayout,
    harmonics,
    real_harmonic_gaunt,
    real_pair,
)
from frozen_sigma_response import herm
from test_frozen_sigma_tangent import complex_model
from test_generalized_sigma_response import endpoint


def fixture():
    rng = np.random.default_rng(542)
    q = rng.normal(size=(21, 3))
    ql = np.linalg.norm(q, axis=1)
    phase = np.exp(1j * rng.normal(size=21))
    c = rng.normal(size=(21, 5)) + 1j * rng.normal(size=(21, 5))
    groups = [np.array([13, 14, 15, 16]), np.array([17, 18, 19, 20])]
    channel = PreparedAugmentationChannel(2, 1.9, 400.0, q, ql, phase, c, 13, groups)
    boundary = rng.normal(size=(4, 4))
    boundary[0, 0] = 1.3
    lo = rng.normal(size=(4, 2))
    return rng, channel, boundary, lo, (q, ql, phase, c, groups)


def literal_coefficients(boundary, lo, lapw, data):
    q, ql, phase, c, groups = data
    r = 1.9
    ell = 2
    geom = 4 * np.pi / np.sqrt(400.0) * 1j**ell * harmonics(ell, q).conj() * phase
    coefficients = np.zeros((5, 4, 21), complex)
    for m in range(5):
        for g in range(13):
            j = spherical_jn(ell, ql[g] * r)
            dj = ql[g] * spherical_jn(ell, ql[g] * r, derivative=True)
            if lapw:
                coefficients[m, 0, g] = (
                    r * (boundary[0, 1] * dj - boundary[2, 1] * j) * geom[m, g]
                )
                coefficients[m, 1, g] = (
                    r * (boundary[2, 0] * j - boundary[0, 0] * dj) * geom[m, g]
                )
            else:
                coefficients[m, 0, g] = r * j / boundary[0, 0] * geom[m, g]
        for slot, group in enumerate(groups):
            for g in group:
                coefficients[m, :, g] = r * r * lo[:, slot] * geom[m, g]
    return coefficients @ c


class AugmentationTests(unittest.TestCase):
    def test_prepared_sums_match_full_native_coefficients(self):
        rng, ch, b, lo, data = fixture()
        for lapw in (False, True):
            np.testing.assert_allclose(
                ch.value(b, lo, lapw=lapw),
                literal_coefficients(b, lo, lapw, data),
                rtol=2e-14,
                atol=2e-14,
            )
        # The cache contains no array whose size scales with the original PW count.
        self.assertEqual(ch.fj.shape, (5, 5))
        self.assertEqual(ch.local.shape, (2, 5, 5))

    def test_prepared_derivative_and_reverse(self):
        rng, ch, b, lo, data = fixture()
        db = rng.normal(size=b.shape)
        dl = rng.normal(size=lo.shape)
        bar = rng.normal(size=(5, 4, 5)) + 1j * rng.normal(size=(5, 4, 5))
        for lapw in (False, True):
            actual = ch.action(b, db, dl, lapw=lapw)
            h = 1e-5
            fd = (
                literal_coefficients(b + h * db, lo + h * dl, lapw, data)
                - literal_coefficients(b - h * db, lo - h * dl, lapw, data)
            ) / (2 * h)
            np.testing.assert_allclose(actual, fd, atol=5e-10, rtol=2e-9)
            bb, lb = ch.pullback(b, bar, lapw=lapw)
            self.assertAlmostEqual(
                real_pair(bar, actual), real_pair(bb, db) + real_pair(lb, dl), places=12
            )

    def test_no_local_orbitals(self):
        rng, ch, b, lo, data = fixture()
        q, ql, phase, c, groups = data
        pure = PreparedAugmentationChannel(2, 1.9, 400.0, q, ql, phase, c, 13, [])
        boundary = b[:, :2]
        empty = np.empty((2, 0))
        actual = pure.value(boundary, empty, lapw=True)
        self.assertEqual(actual.shape, (5, 2, 5))
        self.assertTrue(np.iscomplexobj(actual))
        expected = literal_coefficients(b, np.zeros_like(lo), True, data)[:, :2]
        np.testing.assert_allclose(actual, expected, atol=2e-14)

    def test_gaunt_spherical_limit_and_sine_hermiticity(self):
        np.testing.assert_allclose(
            real_harmonic_gaunt(4, 0, 0), np.eye(25) / np.sqrt(4 * np.pi), atol=6e-15
        )
        cosine = real_harmonic_gaunt(4, 2, 0)
        self.assertAlmostEqual(cosine[0, 6].real, 1 / np.sqrt(4 * np.pi), places=14)
        sine = real_harmonic_gaunt(4, -3, 2)
        np.testing.assert_allclose(sine, sine.conj().T, atol=2e-15)
        self.assertGreater(np.linalg.norm(sine - sine.T), 0.5)

    def test_radial_angular_transposes(self):
        rng = np.random.default_rng(32)
        layout = RadialAngularLayout({0: 3, 1: 2, 2: 4}, [(2, 0), (-3, 2), (4, 4)])
        pairs = {
            (l, lp): rng.normal(size=(n, np_, 3))
            for l, n in layout.counts.items()
            for lp, np_ in layout.counts.items()
        }
        cov = rng.normal(size=(layout.size, layout.size)) + 1j * rng.normal(
            size=(layout.size, layout.size)
        )
        pull = layout.potential_pullback(cov)
        self.assertAlmostEqual(
            real_pair(cov, layout.potential(pairs)),
            sum(real_pair(pull[k], a) for k, a in pairs.items()),
            places=12,
        )
        metrics = {l: rng.normal(size=(n, n)) for l, n in layout.counts.items()}
        pull = layout.overlap_pullback(cov)
        self.assertAlmostEqual(
            real_pair(cov, layout.overlap(metrics)),
            sum(real_pair(pull[l], a) for l, a in metrics.items()),
            places=12,
        )
        q = rng.normal(size=4)
        t = rng.normal(size=(5, 3)) + 1j * rng.normal(size=(5, 3))
        dq = rng.normal(size=q.shape)
        dt = rng.normal(size=t.shape) + 1j * rng.normal(size=t.shape)
        bar = rng.normal(size=(layout.size, 3)) + 1j * rng.normal(size=(layout.size, 3))
        qb, tb = layout.projector_pullback(2, q, t, bar)
        dp = layout.projector(2, dq, t) + layout.projector(2, q, dt)
        self.assertAlmostEqual(
            real_pair(bar, dp), real_pair(qb, dq) + real_pair(tb, dt), places=12
        )

    def test_joint_contractions_nonlinear_and_reverse(self):
        rng = np.random.default_rng(33)
        z = rng.normal(size=(14, 5)) + 1j * rng.normal(size=(14, 5))
        k = herm(rng.normal(size=(14, 14)) + 1j * rng.normal(size=(14, 14)))
        s = np.eye(14)
        p = rng.normal(size=(14, 3)) + 1j * rng.normal(size=(14, 3))
        model = AugmentationContractions(z, k, s, p)
        dz = 0.1 * (rng.normal(size=z.shape) + 1j * rng.normal(size=z.shape))
        dk = 0.1 * herm(rng.normal(size=k.shape) + 1j * rng.normal(size=k.shape))
        ds = 0.1 * herm(rng.normal(size=s.shape) + 1j * rng.normal(size=s.shape))
        dp = 0.1 * (rng.normal(size=p.shape) + 1j * rng.normal(size=p.shape))
        actual = model.action(dz, dk, ds, dp)
        errors = []
        for h in (1e-3, 5e-4):
            ends = [
                AugmentationContractions(
                    z + sign * h * dz,
                    k + sign * h * dk,
                    s + sign * h * ds,
                    p + sign * h * dp,
                ).value()
                for sign in (-1, 1)
            ]
            fd = [(a - b) / (2 * h) for a, b in zip(ends[1], ends[0])]
            errors.append(sum(np.linalg.norm(a - b) for a, b in zip(fd, actual)))
        self.assertGreater(errors[0] / errors[1], 3.8)
        bars = [
            herm(rng.normal(size=(5, 5)) + 1j * rng.normal(size=(5, 5)))
            for _ in range(2)
        ] + [rng.normal(size=(5, 3)) + 1j * rng.normal(size=(5, 3))]
        left = sum(real_pair(a, b) for a, b in zip(bars, actual))
        right = sum(
            real_pair(a, b) for a, b in zip(model.pullback(*bars), (dz, dk, ds, dp))
        )
        self.assertAlmostEqual(left, right, places=11)

    def test_through_normalized_window_dmft_response(self):
        rng = np.random.default_rng(74)
        electronic = complex_model()
        nk, nb = electronic.eps.shape
        prepared = []
        directions = []
        h = []
        o = []
        raw = []
        for ik in range(nk):
            z = np.linalg.qr(rng.normal(size=(9, nb)) + 1j * rng.normal(size=(9, nb)))[
                0
            ]
            k = z @ np.diag(electronic.eps[ik]) @ z.conj().T + 1.3 * (
                np.eye(9) - z @ z.conj().T
            )
            p = z @ electronic.raw[ik, :, :, 0]
            s = np.eye(9)
            a = AugmentationContractions(z, k, s, p)
            prepared.append(a)
            dz = 0.02 * (rng.normal(size=z.shape) + 1j * rng.normal(size=z.shape))
            dk = 0.04 * herm(rng.normal(size=k.shape) + 1j * rng.normal(size=k.shape))
            ds = np.zeros_like(s)
            dp = np.zeros_like(p)
            directions.append((dz, dk, ds, dp))
            dh, do, dphi = a.action(dz, dk, ds, dp)
            h.append(dh)
            o.append(do)
            raw.append(dphi[:, :, None])
        h, o, raw = map(np.array, (h, o, raw))
        actual = electronic.generalized_response(h, o, raw, fixed_number=False)[0]
        errors = []
        for step in (2e-3, 1e-3):
            values = []
            for sign in (-1, 1):
                hs = []
                os = []
                ps = []
                for a, d in zip(prepared, directions):
                    args = [
                        base + sign * step * delta
                        for base, delta in zip((a.z, a.k, a.s, a.p), d)
                    ]
                    hh, oo, pp = AugmentationContractions(*args).value()
                    hs.append(hh)
                    os.append(oo)
                    ps.append(pp[:, :, None])
                dh = np.array(hs) - np.array([np.diag(e) for e in electronic.eps])
                do = np.array(os) - np.eye(nb)
                dr = np.array(ps) - electronic.raw
                values.append(endpoint(electronic, dh, do, dr, 1.0)[0])
            errors.append(np.linalg.norm((values[1] - values[0]) / (2 * step) - actual))
        self.assertLess(errors[-1], 1e-7)
        self.assertGreater(errors[0] / errors[1], 3.7)

    def test_native_matching_and_gaunt_inputs_are_used(self):
        rng, ch, b, lo, data = fixture()
        q, ql, phase, c, groups = data
        matching = np.array(
            [
                spherical_jn(2, ql[:13] * 1.9),
                ql[:13] * spherical_jn(2, ql[:13] * 1.9, derivative=True),
            ]
        )
        native = PreparedAugmentationChannel(
            2, 1.9, 400.0, q, ql, phase, c, 13, groups, matching_values=2 * matching
        )
        expected = 2 * ch.value(b, np.zeros_like(lo), lapw=True)
        np.testing.assert_allclose(
            native.value(b, np.zeros_like(lo), lapw=True), expected, atol=2e-14
        )
        g = real_harmonic_gaunt(2, 2, 0)
        base = RadialAngularLayout({0: 2, 1: 2, 2: 2}, [(2, 0)])
        native = RadialAngularLayout(
            {0: 2, 1: 2, 2: 2}, [(2, 0)], gaunt_matrices=np.array([3 * g])
        )
        pairs = {
            (l, lp): rng.normal(size=(2, 2, 1))
            for l in base.counts
            for lp in base.counts
        }
        np.testing.assert_allclose(
            native.potential(pairs), 3 * base.potential(pairs), atol=2e-14
        )

    def test_invalid_and_missing_inputs_rejected(self):
        rng, ch, b, lo, data = fixture()
        with self.assertRaises(ValueError):
            ch.value(b[:2], lo, lapw=True)
        b[0, 0] = 0
        with self.assertRaises(ValueError):
            ch.value(b, lo, lapw=False)
        layout = RadialAngularLayout({0: 2, 1: 2}, [(2, 0)])
        with self.assertRaises(ValueError):
            layout.potential({})
        with self.assertRaises(ValueError):
            real_harmonic_gaunt(4, -3, 0)


if __name__ == "__main__":
    unittest.main()
