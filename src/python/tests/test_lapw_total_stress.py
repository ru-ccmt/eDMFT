"""Tensor reconstruction must not turn a measured trace into arbitrary axes."""

from itertools import permutations, product
from types import SimpleNamespace as NS
import unittest
import warnings
from unittest.mock import patch
import numpy as np
from lapw_total_stress import diagonal_invariant_strains, stress_from_contractions


def structure(rotations, cell=(4.0, 4.0, 6.0)):
    return NS(
        lattyp="P",
        alpha=90.0,
        beta=90.0,
        gamma=90.0,
        a=cell[0],
        b=cell[1],
        c=cell[2],
        timat=rotations,
    )


class TensorTests(unittest.TestCase):
    def test_unresolved_window_cannot_be_published(self):
        from lapw_total_stress import PreparedTotalStress
        state = NS(valence=NS(branches=(NS(model=NS(diagnostic_only=True)),)))
        with self.assertRaisesRegex(NotImplementedError, "window-omission"):
            PreparedTotalStress(state, ".", "case", complex_vectors=False,
                                library="unused")

    def charge_context(self, residual, *, root=True, threshold=1e-3, warning=None):
        from lapw_total_stress import PreparedTotalStress
        meta = NS(sigma_inf=np.zeros((1, 1, 1)), sigma=np.zeros((1, 1, 1, 1)),
                  omega=np.ones(1), sym=lambda x: x)
        state = NS(valence=NS(branches=(NS(model=NS(diagnostic_only=False, meta=meta)),)),
                   metadata={'central_fixed_point_residual_norm': residual},
                   potential_core=NS(structure=structure([np.eye(3)])),
                   parallel=NS(root=root))
        with patch('lapw_total_stress.invariant_strains', return_value=np.eye(3)[None]), \
             patch('lapw_total_stress.PreparedKleinFunctional') as functional, \
             patch('lapw_total_stress.SearchedEnergyResponse') as energy:
            result = PreparedTotalStress(state, '.', 'case', complex_vectors=False,
                library='unused', fixed_point_tolerance=threshold, warning=warning)
            functional.assert_called_once_with(state)
            energy.assert_called_once()
            return result

    def test_finite_charge_residual_warns_and_continues(self):
        with self.assertWarnsRegex(RuntimeWarning, 'Central charge-map residual'):
            total = self.charge_context(.00503141)
        self.assertFalse(total.diagnostic_only)
        self.assertFalse(total.charge_fixed_point['within_threshold'])
        self.assertEqual(total.charge_fixed_point['residual_norm'], .00503141)
        self.assertEqual(total.charge_fixed_point['warning_threshold'], .001)
        self.assertIn('not WIEN :DIS:', total.charge_fixed_point['warning'])

    def test_charge_warning_uses_root_callback_and_is_silent_below_threshold(self):
        messages = []
        self.charge_context(.005, warning=messages.append)
        self.assertEqual(len(messages), 1)
        with warnings.catch_warnings(record=True) as emitted:
            warnings.simplefilter('always')
            worker = self.charge_context(.005, root=False, warning=messages.append)
            converged = self.charge_context(.0001, warning=messages.append)
            boundary = self.charge_context(.001, warning=messages.append)
        self.assertFalse(emitted)
        self.assertEqual(len(messages), 1)
        self.assertIsNotNone(worker.charge_fixed_point['warning'])
        self.assertTrue(converged.charge_fixed_point['within_threshold'])
        self.assertTrue(boundary.charge_fixed_point['within_threshold'])
        self.assertIsNone(converged.charge_fixed_point['warning'])

    def test_invalid_charge_residual_or_threshold_still_fails(self):
        for residual in [float('nan'), float('inf'), -1.0]:
            with self.subTest(residual=residual), self.assertRaisesRegex(ValueError, 'finite and nonnegative'):
                self.charge_context(residual)
        for threshold in [0, -1, float('nan'), float('inf')]:
            with self.subTest(threshold=threshold), self.assertRaisesRegex(ValueError, 'positive finite'):
                self.charge_context(.005, threshold=threshold)

    def test_tetragonal_two_sources_reconstruct_three_axes(self):
        r = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
        s = structure([np.linalg.matrix_power(r, k) for k in range(4)])
        e = diagonal_invariant_strains(s)
        self.assertEqual(len(e), 2)
        sigma = np.diag([0.2, 0.2, -0.6])
        volume = 96.0
        d = volume * np.einsum("ij,aij->a", sigma, e)
        np.testing.assert_allclose(
            stress_from_contractions(e, d, volume), sigma, atol=1e-15
        )

    def test_cubic_only_hydro_is_enough(self):
        rotations = [
            np.eye(3)[list(p)] @ np.diag(v)
            for p in permutations(range(3))
            for v in product([-1, 1], repeat=3)
        ]
        e = diagonal_invariant_strains(structure(rotations, (4.0, 4.0, 4.0)))
        self.assertEqual(len(e), 1)
        np.testing.assert_allclose(
            stress_from_contractions(e, [-12.0], 20.0), -0.2 * np.eye(3)
        )

    def test_orthorhombic_has_three_independent_axes(self):
        e = diagonal_invariant_strains(
            structure([np.diag(v) for v in product([-1, 1], repeat=3)], (3.0, 4.0, 5.0))
        )
        self.assertEqual(len(e), 3)
        sigma = np.diag([1.0, 2.0, 3.0])
        d = 60 * np.einsum("ij,aij->a", sigma, e)
        np.testing.assert_allclose(
            stress_from_contractions(e, d, 60.0), sigma, atol=1e-14
        )

    def test_allowed_shear_is_rejected_instead_of_set_to_zero(self):
        for rotations in [[np.eye(3)], [np.eye(3), np.diag([-1, -1, 1])]]:
            with self.assertRaisesRegex(NotImplementedError, "shear"):
                diagonal_invariant_strains(structure(rotations))


if __name__ == "__main__":
    unittest.main()
