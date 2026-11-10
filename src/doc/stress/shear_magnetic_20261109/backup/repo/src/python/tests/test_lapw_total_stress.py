"""Tensor reconstruction must not turn a measured trace into arbitrary axes."""

from itertools import permutations, product
from types import SimpleNamespace as NS
import unittest
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
    def test_unvalidated_magnetic_total_cannot_be_published(self):
        from lapw_total_stress import PreparedTotalStress
        state = NS(valence=NS(branches=(NS(), NS())))
        with self.assertRaisesRegex(NotImplementedError, "magnetic response"):
            PreparedTotalStress(state, ".", "case", complex_vectors=False,
                                library="unused")

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
