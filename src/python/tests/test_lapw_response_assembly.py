"""Independent native-triangle controls, including complex coefficient gauges."""

import unittest
import numpy as np

from lapw_response_assembly import NativeBandProjection


class AssemblyTests(unittest.TestCase):
    def test_native_order_and_transpose(self):
        rng = np.random.default_rng(167)
        for complex_case in (False, True):
            c = rng.normal(size=(7, 4))
            angular = rng.normal(size=(7, 7))
            if complex_case:
                c = c + 1j * rng.normal(size=c.shape)
                angular = angular + 1j * rng.normal(size=angular.shape)
            radial = rng.normal(size=(7, 7))  # deliberately nonsymmetric
            native = np.zeros((7, 7), complex)
            for j in range(7):
                native[j, j] = (radial[j, j] * angular[j, j]).real
                for i in range(j):
                    native[i, j] = radial[i, j] * angular[i, j]
                    native[j, i] = native[i, j].conjugate()
            op = NativeBandProjection(angular, c)
            np.testing.assert_allclose(op.project(radial), c.conj().T @ native @ c)
            self.assertLess(op.operator.check_transpose(), 1e-13)
            # Changing a discarded triangle must have no effect.
            np.testing.assert_allclose(op.project(radial), op.project(np.triu(radial)))

    def test_central_band_gauge_covariance(self):
        rng = np.random.default_rng(53)
        c = rng.normal(size=(8, 4)) + 1j * rng.normal(size=(8, 4))
        q, _ = np.linalg.qr(rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4)))
        a = rng.normal(size=(8, 8)) + 1j * rng.normal(size=(8, 8))
        r = rng.normal(size=(8, 8))
        p, rotated = NativeBandProjection(a, c), NativeBandProjection(a, c @ q)
        np.testing.assert_allclose(rotated.project(r), q.conj().T @ p.project(r) @ q)


if __name__ == "__main__":
    unittest.main()
