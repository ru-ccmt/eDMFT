from unittest import TestCase
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from lapw_soc_force_reference import NativeSOCForceMap


class NativeDirectionalTests(TestCase):
    def test_screening_sign_and_reused_direction(self):
        # The charge source and the free-energy gradient need not be collinear.
        j = np.array([[0.16, 0.21], [-0.3, 0.34]])
        source = np.array([0.13, 0.08])
        c = np.array([-0.7, 0.31])
        a = 0.2
        exact = np.linalg.solve(np.eye(2) - j, source)
        obj = NativeSOCForceMap.__new__(NativeSOCForceMap)
        obj.structure = SimpleNamespace(mult=[2])
        obj.n = np.zeros(2)
        calls = []

        def native(pattern, h, q):
            calls.append(q.copy())
            return h * (source + j @ q), np.array([1 + h * (a + c @ q)])

        obj.evaluate = native
        result, q = obj.solve(
            1, np.array([0.0, 0.0, 1.0]), None, 0.001, 1e-9, 20, lambda m: None
        )
        np.testing.assert_allclose(q, exact, atol=1.0e-9)
        self.assertAlmostEqual(
            result["force_component_mRy_Bohr"], -500 * (a + c @ exact)
        )
        self.assertAlmostEqual(result["explicit_force_mRy_per_Bohr"], -500 * a)
        self.assertAlmostEqual(result["screened_force_mRy_per_Bohr"], -500 * c @ exact)
        calls.clear()
        smaller, q2 = obj.solve(
            1, np.array([0.0, 0.0, 1.0]), None, 0.0005, 1e-9, 20, lambda m: None, q
        )
        self.assertEqual(len(calls), 4)
        self.assertTrue(smaller["response"]["converged"])
        np.testing.assert_array_equal(q2, q)
