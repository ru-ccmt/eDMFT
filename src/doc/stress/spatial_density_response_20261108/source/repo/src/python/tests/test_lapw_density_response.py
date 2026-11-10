"""Independent spatial/quadrature controls for forward and reverse density maps."""

import unittest
import numpy as np
from numpy.polynomial.legendre import leggauss
from lapw_density_response import (
    MuffinTinDensity,
    MuffinTinDensitySum,
    SpatialValenceDensity,
    InterstitialDensity,
    FourierStars,
    herm,
)
from lapw_augmentation_response import harmonics


def pair(a, b):
    return np.vdot(a, b).real


def fixture():
    rng = np.random.default_rng(731)
    phi = {l: rng.normal(size=(7, 2, 2 + (l % 2))) for l in range(3)}
    n = sum((2 * l + 1) * p.shape[2] for l, p in phi.items())
    z = rng.normal(size=(n, 4)) + 1j * rng.normal(size=(n, 4))
    rho = herm(rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4)))
    hs = [(0, 0), (2, 0), (2, 1), (-2, 1), (-3, 2), (4, 3)]
    return rng, phi, z, rho, hs


def angular_oracle(phi, z, rho, hs, small, clmsum):
    x, w = leggauss(13)
    az = 2 * np.pi * np.arange(27) / 27
    points = (
        np.array(
            [
                np.sqrt(1 - x * x)[:, None] * np.cos(az),
                np.sqrt(1 - x * x)[:, None] * np.sin(az),
                np.broadcast_to(x[:, None], (13, 27)),
            ]
        )
        .reshape(3, -1)
        .T
    )
    weights = np.repeat(w, 27) * 2 * np.pi / 27
    psi = np.zeros((7, 2, len(points), 4), complex)
    offset = 0
    for l, p in phi.items():
        size = (2 * l + 1) * p.shape[2]
        zz = z[offset : offset + size].reshape(2 * l + 1, p.shape[2], 4)
        offset += size
        psi += np.einsum("mp,rci,min->rcpn", harmonics(l, points), p, zz)
    density = sum(
        weight * np.einsum("rpn,nm,rpm->rp", psi[:, c], rho, psi[:, c].conj()).real
        for c, weight in enumerate([1.0, small])
    )
    result = []
    for L, m in hs:
        y = harmonics(abs(L), points)
        k = y[abs(L) + m]
        if m:
            k = (
                ((-1) ** m * k + y[abs(L) - m]) / np.sqrt(2)
                if L > 0
                else -1j * ((-1) ** m * k - y[abs(L) - m]) / np.sqrt(2)
            )
        scale = np.sqrt(4 * np.pi) if clmsum and (L, m) == (0, 0) else 1.0
        result.append(scale * np.einsum("rp,p->r", density, weights * k.real))
    return np.array(result)


class DensityTests(unittest.TestCase):
    def test_k_and_site_sum_before_radial_contraction(self):
        rng, phi, z, rho, hs = fixture()
        rhos = np.array([rho, 2 * rho])
        zs = [z, 0.7 * z.conj(), 1.2 * z]
        entries = [(0, zs[0], 0.5), (0, zs[1], 0.5), (1, zs[2], 1.0)]
        summed = MuffinTinDensitySum(entries, rhos, phi, hs, 0.0001)
        expected = sum(
            w * angular_oracle(phi, zz, rhos[k], hs, 0.0001, True)
            for k, zz, w in entries
        )
        np.testing.assert_allclose(summed.value(), expected, atol=1e-11)
        dr = np.array(
            [
                herm(rng.normal(size=rho.shape) + 1j * rng.normal(size=rho.shape))
                for _ in rhos
            ]
        )
        dz = [rng.normal(size=z.shape) + 1j * rng.normal(size=z.shape) for _ in entries]
        dp = {l: rng.normal(size=p.shape) for l, p in phi.items()}
        y = rng.normal(size=expected.shape)
        derivative = summed.action(dr, dz, dp)
        expected = sum(
            w * MuffinTinDensity(zz, rhos[k], phi, hs, 0.0001).action(dr[k], dz[i], dp)
            for i, (k, zz, w) in enumerate(entries)
        )
        np.testing.assert_allclose(derivative, expected, atol=1e-11)
        rb, zb, pb = summed.pullback(y)
        self.assertAlmostEqual(
            pair(y, derivative),
            pair(rb, dr)
            + sum(pair(a, b) for a, b in zip(zb, dz))
            + sum(pair(pb[l], dp[l]) for l in phi),
            places=10,
        )
        light = summed.pullback(y, include_basis=False)
        np.testing.assert_allclose(light[0], rb, atol=1e-12)
        self.assertEqual(light[1], [])
        self.assertIsNone(light[2])
        np.testing.assert_allclose(
            summed.action(None, dz, dp) + summed.action(dr), derivative, atol=1e-11
        )

    def test_complete_density_coordinates_and_chain_operator(self):
        from wien_response_density import Density

        rng, phi, z, rho, hs = fixture()
        rhos = np.array([rho, 1.3 * rho])
        entries = [(0, z, 1.0), (1, 0.7 * z.conj(), 1.0)]
        mt = MuffinTinDensitySum(entries, rhos, phi, hs, 0.0001)
        reps = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [-1, 0, 0]])
        stars = FourierStars(
            (9, 9, 9), reps, np.eye(3, dtype=int)[None], np.zeros((1, 3))
        )
        labels = np.array([[0, 0, 0], [1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, 0, 1]])
        coeff = [rng.normal(size=(5, 4)) + 1j * rng.normal(size=(5, 4)) for _ in rhos]
        ir = [InterstitialDensity(c, labels, 123.0, stars) for c in coeff]
        template = Density.__new__(Density)
        template.rad = {(0, *h): np.zeros(7) for h in hs}
        template.pw = {tuple(g): 0j for g in reps}
        full = SpatialValenceDensity(template, {0: mt}, ir)
        actual = template.from_vector(full.value())
        np.testing.assert_allclose(
            np.array(list(actual.rad.values())), mt.value(), atol=1e-13
        )
        np.testing.assert_allclose(
            np.array(list(actual.pw.values())),
            sum(m.value(r) for m, r in zip(ir, rhos)),
            atol=1e-13,
        )
        dr = np.array(
            [
                herm(rng.normal(size=rho.shape) + 1j * rng.normal(size=rho.shape))
                for _ in rhos
            ]
        )
        dz = {
            0: [
                rng.normal(size=z.shape) + 1j * rng.normal(size=z.shape)
                for _ in entries
            ]
        }
        dp = {0: {l: rng.normal(size=p.shape) for l, p in phi.items()}}
        y = rng.normal(size=len(template.vector()))
        derivative = full.action(dr, dz, dp)
        rb, zb, pb = full.pullback(y)
        self.assertAlmostEqual(
            y @ derivative,
            pair(rb, dr)
            + sum(pair(a, b) for a, b in zip(zb[0], dz[0]))
            + sum(pair(pb[0][l], dp[0][l]) for l in phi),
            places=10,
        )
        op, coords = full.fixed_basis_operator()
        np.testing.assert_allclose(
            op.matvec(coords.pack(dr)), full.action(dr), atol=1e-13
        )
        np.testing.assert_allclose(coords.unpack(op.rmatvec(y)), rb, atol=1e-11)
        np.testing.assert_allclose(
            full.action(dr) + full.action(None, dz, dp), derivative, atol=1e-11
        )
        with self.assertRaises(ValueError):
            full.action(dr, matching_dots={})

    def test_mt_value_against_spatial_angular_integration(self):
        _, phi, z, rho, hs = fixture()
        for small in [0.0, 1 / 137.0359895**2]:
            for clmsum in [False, True]:
                m = MuffinTinDensity(z, rho, phi, hs, small, clmsum=clmsum)
                np.testing.assert_allclose(
                    m.value(),
                    angular_oracle(phi, z, rho, hs, small, clmsum),
                    atol=5e-12,
                    rtol=3e-13,
                )

    def test_mt_joint_derivative_and_reverse(self):
        rng, phi, z, rho, hs = fixture()
        small = 1 / 137.0359895**2
        m = MuffinTinDensity(z, rho, phi, hs, small)
        dr = herm(rng.normal(size=rho.shape) + 1j * rng.normal(size=rho.shape))
        dz = rng.normal(size=z.shape) + 1j * rng.normal(size=z.shape)
        dp = {l: rng.normal(size=p.shape) for l, p in phi.items()}
        y = rng.normal(size=(len(hs), 7))
        h = 2e-5
        ends = [
            angular_oracle(
                {l: p + s * h * dp[l] for l, p in phi.items()},
                z + s * h * dz,
                rho + s * h * dr,
                hs,
                small,
                True,
            )
            for s in [-1, 1]
        ]
        actual = m.action(dr, dz, dp)
        np.testing.assert_allclose(
            actual, (ends[1] - ends[0]) / (2 * h), atol=2e-6, rtol=3e-8
        )
        rb, zb, pb = m.pullback(y)
        self.assertAlmostEqual(
            pair(y, actual),
            pair(rb, dr) + pair(zb, dz) + sum(pair(pb[l], dp[l]) for l in phi),
            places=10,
        )

    def test_mt_basis_terms_are_separate_and_required(self):
        _, phi, z, rho, hs = fixture()
        m = MuffinTinDensity(z, rho, phi, hs, 0.0)
        zero = np.zeros_like(rho)
        np.testing.assert_allclose(m.action(zero, z_dot=z), 2 * m.value(), atol=1e-12)
        np.testing.assert_allclose(
            m.action(zero, primitive_dot=phi), 2 * m.value(), atol=1e-12
        )
        np.testing.assert_allclose(m.action(rho), m.value(), atol=1e-12)
        # Arbitrary complex band gauge must change rho contragrediently.
        q = np.linalg.qr(z[:4])[0]
        transformed = MuffinTinDensity(z @ q, q.conj().T @ rho @ q, phi, hs, 0.0)
        np.testing.assert_allclose(m.value(), transformed.value(), atol=1e-12)

    def test_stars_with_duplicate_operations_and_complex_phases(self):
        shape = (9, 11, 13)
        reps = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0]])
        rot = np.array(
            [np.eye(3, dtype=int), np.eye(3, dtype=int), -np.eye(3, dtype=int)]
        )
        tau = np.array([[0.0, 0, 0], [0.25, 0, 0], [0, 0, 0]])
        rng = np.random.default_rng(97)
        f = rng.normal(size=shape) + 1j * rng.normal(size=shape)
        y = rng.normal(size=3) + 1j * rng.normal(size=3)
        for phases in [True, False]:
            stars = FourierStars(shape, reps, rot, tau, complex_phases=phases)
            result = stars.action(f)
            self.assertEqual(result[0], f[0, 0, 0])
            factor = (1 - 1j) / 2 if phases else 0.5
            self.assertAlmostEqual(
                abs(result[1] - (factor * f[1, 0, 0] + f[-1, 0, 0])), 0, places=14
            )
            self.assertAlmostEqual(
                pair(y, result), pair(stars.pullback(y), f), places=13
            )

    def test_interstitial_pair_sum_derivative_reverse_and_cache(self):
        rng = np.random.default_rng(18)
        shape = (9, 11, 13)
        g = np.array([[0, 0, 0], [1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, 0, 1]])
        reps = np.array(
            [
                [i, j, k]
                for i in range(-2, 3)
                for j in range(-2, 3)
                for k in range(-2, 3)
            ]
        )
        stars = FourierStars(shape, reps, np.eye(3, dtype=int)[None], np.zeros((1, 3)))
        c = rng.normal(size=(5, 3)) + 1j * rng.normal(size=(5, 3))
        rho = herm(rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3)))
        dr = herm(rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3)))
        dc = rng.normal(size=c.shape) + 1j * rng.normal(size=c.shape)

        def exact(c, rho):
            d = c @ rho @ c.conj().T
            return np.array(
                [
                    sum(
                        d[i, j]
                        for i in range(len(g))
                        for j in range(len(g))
                        if np.array_equal(g[i] - g[j], h)
                    )
                    / 321.0
                    for h in reps
                ]
            )

        y = rng.normal(size=len(reps)) + 1j * rng.normal(size=len(reps))
        h = 1e-5
        for cache in [True, False]:
            m = InterstitialDensity(c, g, 321.0, stars, cache_wavefunctions=cache)
            np.testing.assert_allclose(m.value(rho), exact(c, rho), atol=3e-17)
            act = m.action(rho, dr, dc)
            fd = (exact(c + h * dc, rho + h * dr) - exact(c - h * dc, rho - h * dr)) / (
                2 * h
            )
            np.testing.assert_allclose(act, fd, atol=2e-11, rtol=2e-9)
            rb, cb = m.pullback(y, rho)
            self.assertAlmostEqual(pair(y, act), pair(rb, dr) + pair(cb, dc), places=13)
            np.testing.assert_allclose(rb, m.pullback(y), atol=1e-15)

    def test_reject_inconsistent_inputs(self):
        _, phi, z, rho, hs = fixture()
        with self.assertRaises(ValueError):
            MuffinTinDensity(z, rho + 1j * np.eye(4), phi, hs, 0.0)
        m = MuffinTinDensity(z, rho, phi, hs, 0.0)
        with self.assertRaises(ValueError):
            m.action(rho, primitive_dot={})
        stars = FourierStars(
            (5, 5, 5),
            np.zeros((1, 3), int),
            np.eye(3, dtype=int)[None],
            np.zeros((1, 3)),
        )
        with self.assertRaises(ValueError):
            InterstitialDensity(
                np.ones((2, 1)), np.array([[-2, 0, 0], [2, 0, 0]]), 1.0, stars
            )


if __name__ == "__main__":
    unittest.main()
