"""Compare searched-energy screening to independent nonlinear equilibria."""

from types import SimpleNamespace as NS
from unittest.mock import patch
import unittest
import numpy as np
from scipy.optimize import root
from stress_response_adjoint import PairedOperator
from lapw_energy_response import SearchedEnergyResponse
from lapw_radial_energy_policy import search_root_terms


def operator(a):
    return PairedOperator(a.shape[1], a.shape[0], lambda x: a @ x, lambda y: a.T @ y)


class EnergyScreeningTests(unittest.TestCase):
    def test_native_cont_sentinel_is_not_a_radial_root(self):
        self.assertEqual(search_root_terms(0.3, -0.204, -520.0), [])
        self.assertEqual(search_root_terms(0.3, -520.0, -520.0), [])
        self.assertEqual(
            search_root_terms(-1.8137, -2.106, -1.521),
            [(-2.106, 2, 0.5), (-1.521, 0, 0.5)],
        )
        with self.assertRaises(NotImplementedError):
            search_root_terms(-0.204, -0.204, -520.0)

    def test_joint_search_screening_scalar_and_adjoint(self):
        # A nonlinear charge map and nonlinear searched energies have their
        # common fixed point at zero. Only their central derivatives are
        # supplied to the response class; the endpoint oracle solves the
        # nonlinear map independently for each finite strain.
        rng = np.random.default_rng(119)
        n = 7
        m = 9
        p = 2
        K = rng.normal(size=(m, n)) * 0.11
        F = rng.normal(size=(n, m)) * 0.09
        A = rng.normal(size=(p, m)) * 0.14
        P = rng.normal(size=(n, p)) * 0.17
        g = rng.normal(size=p)
        c0 = rng.normal(size=n)
        s0 = rng.normal(size=n)
        vexp = rng.normal(size=m)
        a0 = 0.81
        chain = NS(
            potential=operator(K),
            charge_from_potential=operator(F),
            normalization=operator(np.eye(n)),
            density=operator(np.zeros((n, p))),
            electronic=operator(np.eye(p)),
        )
        state = NS(
            chain=chain,
            jacobian=operator(F @ K),
            valence=None,
            potential_core=NS(hartree=NS(keys=[(0, 0, 0)], slices=[slice(0, m)])),
        )
        obs = NS(charge_gradient=c0, lattice=NS(joint_gradient=g))
        details = [dict(atom=0, l=0, slots=[i]) for i in range(p)]

        def direction(v, energies):
            e = energies[0, 0]
            return e, P @ e

        with patch(
            "lapw_energy_response.searched_energy_derivatives",
            return_value=({(0, 0): A}, details),
        ), patch("lapw_energy_response.energy_direction", side_effect=direction):
            model = SearchedEnergyResponse(
                state, obs, ".", "test", complex_vectors=False, library="unused"
            )
        self.assertLess(model.jacobian.check_transpose(), 1e-14)
        s = model.complete_source(s0, vexp)
        q, report = model.screened.solve(s, rtol=1e-13, atol=1e-13)
        explicit = a0 + model.explicit_addition(vexp)
        d = explicit + model.charge_gradient @ q
        values, lam, _ = model.screened.contractions(
            model.charge_gradient,
            s[:, None],
            explicit=[explicit],
            rtol=1e-13,
            atol=1e-13,
        )
        self.assertAlmostEqual(d, values[0], places=12)

        def scalar(t):
            def energies(x):
                return np.sin(A @ (K @ x + t * vexp))

            sol = root(
                lambda x: x - t * s0 - np.tanh(F @ (K @ x)) - P @ energies(x),
                t * q,
                tol=1e-11,
            )
            if np.linalg.norm(sol.fun) > 1e-12:
                raise AssertionError("nonlinear oracle failed")
            return a0 * t + c0 @ sol.x + g @ energies(sol.x) + 0.3 * np.sum(sol.x**2)

        errors = []
        for h in [0.002, 0.001]:
            errors.append(abs((scalar(h) - scalar(-h)) / (2 * h) - d))
        self.assertLess(errors[-1], 1e-6)
        self.assertGreater(errors[0] / errors[1], 3.8)


if __name__ == "__main__":
    unittest.main()
