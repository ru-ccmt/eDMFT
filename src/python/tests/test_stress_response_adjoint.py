import unittest
import numpy as np
from frozen_sigma_tangent import real_inner
from stress_response_adjoint import (
    PairedOperator,
    HermitianCoordinates,
    electronic_operator,
    ScreenedResponse,
    ResponseConvergenceError,
    tensile_tensor,
    engineering_strain_bases,
    RY_BOHR3_TO_GPA,
)
from test_frozen_sigma_tangent import complex_model


def matrix_operator(a):
    return PairedOperator(a.shape[1], a.shape[0], lambda x: a @ x, lambda y: a.T @ y)


class AdjointTests(unittest.TestCase):
    def test_hermitian_coordinates_are_isometric(self):
        coordinates = HermitianCoordinates((2, 5, 5))
        rng = np.random.default_rng(90)
        x, y = rng.normal(size=(2, coordinates.size))
        a, b = coordinates.unpack(x), coordinates.unpack(y)
        np.testing.assert_allclose(coordinates.pack(a), x, atol=5e-16)
        self.assertAlmostEqual(real_inner(a, b), x @ y, places=13)

    def test_composed_rectangular_transpose(self):
        rng = np.random.default_rng(18)
        a, b = rng.normal(size=(6, 7)), rng.normal(size=(7, 3))
        operator = matrix_operator(a) @ matrix_operator(b)
        self.assertLess(operator.check_transpose(), 1e-14)

    def test_electronic_transpose_coordinates(self):
        operator, _ = electronic_operator(complex_model())
        self.assertLess(operator.check_transpose(), 1e-13)

    def test_one_adjoint_matches_six_forward_solves(self):
        rng = np.random.default_rng(49)
        j = rng.normal(size=(17, 17)) * 0.05
        # A nonnormal map with unstable simple mixing; do not assume SPD.
        j[0, 0] = -7
        j[0, 1] = 1.7
        solver = ScreenedResponse(matrix_operator(j))
        sources = rng.normal(size=(17, 6))
        gradient = rng.normal(size=17)
        explicit = rng.normal(size=6)
        answer, lam, report = solver.contractions(
            gradient, sources, explicit, rtol=1e-12, atol=1e-13
        )
        forwards = np.column_stack(
            [solver.solve(sources[:, nu], rtol=1e-12, atol=1e-13)[0] for nu in range(6)]
        )
        np.testing.assert_allclose(
            answer, explicit + gradient @ forwards, rtol=1e-11, atol=1e-11
        )
        np.testing.assert_allclose(
            lam, np.linalg.solve(np.eye(17) - j.T, gradient), atol=1e-11
        )
        self.assertTrue(report["converged"])

    def test_physical_chain_electrons_and_common_gradient(self):
        model = complex_model()
        chi, coordinates = electronic_operator(model)
        rng = np.random.default_rng(112)
        # Reduced-coordinate diagnostic: a paired potential/reconstruction
        # map. This is not a claimed LAPW spatial Hartree/XC operator.
        b = rng.normal(size=(coordinates.size, 8)) * 0.01
        potential = matrix_operator(b)
        density = matrix_operator(b.T)
        j = density @ chi @ potential
        j.check_transpose()
        c = potential.rmatvec(coordinates.pack(model.projector_potential_gradient()))
        sources = rng.normal(size=(8, 6))
        solver = ScreenedResponse(j)
        answer, _, _ = solver.contractions(c, sources, rtol=1e-12, atol=1e-14)
        forward = [
            solver.solve(sources[:, i], rtol=1e-12, atol=1e-14)[0] for i in range(6)
        ]
        np.testing.assert_allclose(answer, np.array(forward) @ c, atol=1e-13)

    def test_bad_transpose_rejected(self):
        j = np.array([[0.1, 1.2], [-0.2, 0.4]])
        wrong = PairedOperator(2, 2, lambda x: j @ x, lambda y: j @ y)
        with self.assertRaisesRegex(ValueError, "transpose test failed"):
            wrong.check_transpose()

    def test_singular_screening_is_not_reported_as_stress(self):
        solver = ScreenedResponse(matrix_operator(np.eye(2)))
        with self.assertRaises(ResponseConvergenceError):
            solver.contractions(np.ones(2), np.ones((2, 6)))

    def test_shear_convention_and_hydrostatic_trace(self):
        stress = np.array([[1.0, 4.0, 5.0], [4.0, 2.0, 6.0], [5.0, 6.0, 3.0]])
        volume = 71.0
        derivatives = (
            np.einsum("vij,ij->v", engineering_strain_bases(), stress)
            * volume
            / RY_BOHR3_TO_GPA
        )
        np.testing.assert_allclose(tensile_tensor(derivatives, volume), stress)
        self.assertAlmostEqual(
            derivatives[:3].sum() * RY_BOHR3_TO_GPA / (3 * volume), np.trace(stress) / 3
        )


if __name__ == "__main__":
    unittest.main()
