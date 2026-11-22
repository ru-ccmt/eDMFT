import unittest
import numpy as np
from scipy.linalg import eigh
from soc_scalar_window import outside_tangent


class ScalarWindowTest(unittest.TestCase):
    def test_off_diagonal_overlap_and_gauge(self):
        O = np.array([[1.2, 0.17, 0.08], [0.17, 0.9, -0.10], [0.08, -0.10, 1.3]])
        H = np.array([[0.4, 0.3, 0.01], [0.3, -0.8, 0.14], [0.01, 0.14, 0.2]])
        dH = np.array([[0.3, 0.1, 0.08], [0.1, 0.06, -0.3], [0.08, -0.3, -0.1]])
        dO = np.array([[-0.2, 0.13, -0.16], [0.13, 0.1, 0.2], [-0.16, 0.2, 0.22]])
        e, C = eigh(H, O)
        # The retained space contains holes at either end of the spectrum.
        A = C[:, 1:2] * np.exp(0.3j)
        tangent, gap = outside_tangent(H, O, e[1:2], A, dH, dO)
        t = 1.0e-5
        vectors = []
        for x in (-t, t):
            _, B = eigh(H + x * dH, O + x * dO)
            B = B[:, 1:2] * np.exp(0.3j)
            B *= np.exp(-1j * np.angle((A.conj().T @ O @ B).item()))
            vectors.append(B)
        fd = (vectors[1] - vectors[0]) / (2 * t)
        fd -= A @ (A.conj().T @ O @ fd)
        np.testing.assert_allclose(tangent, fd, atol=2e-9)
        self.assertGreater(gap, 0.0)
        np.testing.assert_allclose(A.conj().T @ O @ tangent, 0, atol=1e-14)
        whole, gap = outside_tangent(H, O, e, C, dH, dO)
        np.testing.assert_array_equal(whole, 0)
        self.assertEqual(gap, np.inf)

    def test_retained_degenerate_pair_can_rotate(self):
        H = np.diag([0.1, 0.1, 1.6])
        O = np.eye(3)
        dH = np.array([[0., .1, .3j], [.1, 0., .2], [-.3j, .2, .4]])
        dO = np.zeros_like(O)
        A = O[:, :2]
        U = np.array([[1., 1j], [1j, 1.]])/np.sqrt(2.)
        e = np.array([.1, .1])
        dA, _ = outside_tangent(H, O, e, A, dH, dO)
        rotated, _ = outside_tangent(H, O, e, A @ U, dH, dO)
        np.testing.assert_allclose(rotated, dA @ U, atol=1e-14)

    def test_split_degeneracy(self):
        eye = np.eye(2)
        with self.assertRaisesRegex(ValueError, "degenerate scalar"):
            outside_tangent(eye, eye, np.ones(1), eye[:, :1], eye, eye)
