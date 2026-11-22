"""Compare complex band-space atomic derivatives against independent displaced geometries."""
import unittest
import numpy as np
from lapw_geometry_response import interstitial_geometry
from lapw_hartree_periodic import NativeStarTable
from lapw_strain_field import harmonic_average_direction


class AtomicGeometryTests(unittest.TestCase):
    def test_moving_interstitial(self):
        rng = np.random.default_rng(215)
        # More rows than fit in one work block, with an oblique reciprocal basis.
        labels = rng.choice(1000, 295, replace=False)
        labels = np.array(np.unravel_index(labels, (10, 10, 10))).T - 5
        reciprocal = np.array([[1., .2, 0.], [0., .8, .12], [.1, 0., .9]])
        q = (labels + [.1, .2, -.15]) @ reciprocal.T
        c = rng.normal(size=(len(q), 6)) + 1j*rng.normal(size=(len(q), 6))
        c /= np.sqrt(len(q))
        p = np.array([[.27, .24, .17], [.15, .71, .36]])
        u = np.array([[.47, -.31, .12], [-.1, .08, .23]])
        dp = u @ reciprocal / (2*np.pi)
        def blocks(pos, **kw):
            return interstitial_geometry(q, labels, c, pos, [1.2, 1.6], 502., np.zeros((3, 3)), 3, **kw)
        actual = blocks(p, displacements=u, central=False)[:2]
        delta = 2e-6
        plus = blocks(p + delta * dp)[2:]
        minus = blocks(p - delta * dp)[2:]
        for a, b, d in zip(actual, plus, minus):
            np.testing.assert_allclose(a, (b-d)/(2*delta), atol=3e-8, rtol=3e-7)

    def test_moving_star_phase(self):
        reps = [[0,0,0], [1,-2,1], [2,0,1]]
        stars = NativeStarTable(np.array(reps), np.eye(3, dtype=int)[None], np.zeros((1,3)), inversion=False)
        b = np.array([[.9,.11,0.], [0.,1.,.2], [.1,0.,.68]])
        u = np.array([.1,-.32,.28])
        pos = np.array([.17,.21,.36])
        dp = b.T @ u / (2*np.pi)
        lm = [(l,m) for l in range(4) for m in range(-l,l+1)]
        h = 1e-6
        def average(site):
            return stars.harmonic_average(b, site, np.eye(3), lm)
        actual = harmonic_average_direction(stars, b, pos, np.eye(3), lm, np.zeros((3,3)), displacement=u)
        np.testing.assert_allclose(actual, (average(pos+h*dp)-average(pos-h*dp))/(2*h), atol=2e-10)


if __name__ == '__main__':
    unittest.main()
