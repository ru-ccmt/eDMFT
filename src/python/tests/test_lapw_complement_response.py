"""Independent spectral, native-metric and nonlinear complement controls."""

import unittest
from unittest.mock import patch
import numpy as np
from scipy.linalg import block_diag, eigh

from frozen_sigma_response import herm
from frozen_sigma_tangent import CachedFrozenSigmaResponse, real_inner
from lapw_complement_response import (
    ShiftedComplement,
    LAPWComplementBasis,
    EmptyComplementResponse,
    complement_electronic_operator,
)
from test_frozen_sigma_tangent import complex_model
from test_generalized_sigma_response import endpoint


def unitary(rng, n):
    return np.linalg.qr(rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n)))[0]


def fixture():
    retained = complex_model()
    nk, nr = retained.eps.shape
    rng = np.random.default_rng(754)
    nq = 4
    eq = np.array([[3.0, 3.0, 4.0, 6.0], [3.2, 3.7, 4.3, 6.8]])
    rotations = np.array([unitary(rng, nq) for k in range(nk)])
    raw_q = rng.normal(size=(nk, nq, 2, 1)) + 1j * rng.normal(size=(nk, nq, 2, 1))
    raw = np.concatenate((retained.raw, raw_q), axis=1)
    # This fixture has one impurity block, so local/native layouts coincide.
    full = CachedFrozenSigmaResponse(
        retained.meta,
        np.concatenate((retained.eps, eq), axis=1),
        raw,
        retained.apply_x(raw, retained.X),
        retained.kw,
        retained.W,
    )
    solvers = [
        ShiftedComplement(
            herm(u.conj().T @ np.diag(e) @ u),
            er,
            mu=retained.meta.mu,
            temperature=retained.meta.temperature,
        )
        for u, e, er in zip(rotations, eq, retained.eps)
    ]
    rq = [np.einsum("an,nic->aic", u.conj().T, q) for u, q in zip(rotations, raw_q)]
    model = EmptyComplementResponse(retained, solvers, rq)
    transform = np.array([block_diag(np.eye(nr), u) for u in rotations])
    return retained, full, model, transform, rng


class ComplementTests(unittest.TestCase):
    def test_block_coordinates_and_paired_operator(self):
        small, full, m, t, rng = fixture()
        for fixed in [False, True]:
            op, inputs, outputs = complement_electronic_operator(m, fixed_number=fixed)
            self.assertLess(op.check_transpose(samples=2, tolerance=2e-11), 2e-11)
            v = rng.normal(size=inputs.size)
            np.testing.assert_allclose(inputs.pack(*inputs.unpack(v)), v, atol=5e-16)
            y = rng.normal(size=outputs.size)
            rr, q = outputs.unpack(y)
            self.assertAlmostEqual(
                y @ y,
                real_inner(rr, rr) + 2 * sum(real_inner(a, a) for a in q),
                places=11,
            )

    def test_shifted_solve_matches_full_spectral_and_transpose(self):
        small, full, m, t, rng = fixture()
        for solver in m.solvers:
            rhs = rng.normal(size=(solver.nq, solver.nr)) + 1j * rng.normal(
                size=(solver.nq, solver.nr)
            )
            y = rng.normal(size=rhs.shape) + 1j * rng.normal(size=rhs.shape)
            e, u = eigh(solver.h)
            expected = u @ (
                (u.conj().T @ rhs) / (solver.energies[None, :] - e[:, None])
            )
            with patch(
                "lapw_complement_response.hessenberg",
                side_effect=AssertionError("unexpected reduction"),
            ):
                actual = solver.action(rhs)
                backward = solver.action(y)
            np.testing.assert_allclose(actual, expected, atol=2e-14)
            self.assertLess(solver.residual(actual, rhs), 1e-13)
            self.assertAlmostEqual(
                real_inner(y, actual), real_inner(backward, rhs), places=12
            )

    def test_full_generalized_response_both_ensembles_and_reverse(self):
        small, full, m, t, rng = fixture()
        nr = small.eps.shape[1]
        n = full.eps.shape[1]
        h = herm(rng.normal(size=(2, n, n)) + 1j * rng.normal(size=(2, n, n)))
        o = 0.3 * herm(rng.normal(size=h.shape) + 1j * rng.normal(size=h.shape))
        p = rng.normal(size=full.raw.shape) + 1j * rng.normal(size=full.raw.shape)
        y = herm(rng.normal(size=h.shape) + 1j * rng.normal(size=h.shape))
        td = t.swapaxes(1, 2).conj()
        hh, oo, yy = td @ h @ t, td @ o @ t, td @ y @ t
        pp = np.einsum("knm,kmic->knic", td, p)
        for fixed in [False, True]:
            expected, info = full.generalized_response(h, o, p, fixed_number=fixed)
            rr, qr, meta = m.response(
                hh[:, :nr, :nr],
                oo[:, :nr, :nr],
                pp[:, :nr],
                hh[:, nr:, :nr],
                oo[:, nr:, :nr],
                fixed_number=fixed,
            )
            mixed = np.zeros_like(h)
            mixed[:, :nr, :nr] = rr
            mixed[:, nr:, :nr] = qr
            mixed[:, :nr, nr:] = np.array(qr).swapaxes(1, 2).conj()
            np.testing.assert_allclose(t @ mixed @ td, expected, atol=3e-12, rtol=2e-12)
            self.assertAlmostEqual(meta["dmu"], info["dmu"], places=11)
            hb, ob, rb, hq, oq = m.adjoint(
                yy[:, :nr, :nr], yy[:, nr:, :nr], fixed_number=fixed
            )
            forward = real_inner(yy[:, :nr, :nr], rr) + 2 * real_inner(
                yy[:, nr:, :nr], qr
            )
            reverse = (
                real_inner(hb, hh[:, :nr, :nr])
                + real_inner(ob, oo[:, :nr, :nr])
                + real_inner(rb, pp[:, :nr])
                + 2 * real_inner(hq, hh[:, nr:, :nr])
                + 2 * real_inner(oq, oo[:, nr:, :nr])
            )
            self.assertAlmostEqual(forward, reverse, places=10)

    def test_full_nonlinear_endpoint(self):
        small, full, m, t, rng = fixture()
        nr = small.eps.shape[1]
        n = full.eps.shape[1]
        h = 0.1 * herm(rng.normal(size=(2, n, n)) + 1j * rng.normal(size=(2, n, n)))
        o = 0.05 * herm(rng.normal(size=h.shape) + 1j * rng.normal(size=h.shape))
        p = 0.05 * (
            rng.normal(size=full.raw.shape) + 1j * rng.normal(size=full.raw.shape)
        )
        td = t.swapaxes(1, 2).conj()
        hh = td @ h @ t
        oo = td @ o @ t
        pp = np.einsum("knm,kmic->knic", td, p)
        rr, q, _ = m.response(
            hh[:, :nr, :nr],
            oo[:, :nr, :nr],
            pp[:, :nr],
            hh[:, nr:, :nr],
            oo[:, nr:, :nr],
            fixed_number=False,
        )
        mixed = np.zeros_like(h)
        mixed[:, :nr, :nr] = rr
        mixed[:, nr:, :nr] = q
        mixed[:, :nr, nr:] = np.array(q).swapaxes(1, 2).conj()
        analytical = t @ mixed @ td
        errors = []
        for step in [2e-4, 1e-4]:
            plus, *_ = endpoint(full, h, o, p, step)
            minus, *_ = endpoint(full, h, o, p, -step)
            errors.append(np.linalg.norm((plus - minus) / (2 * step) - analytical))
        self.assertLess(errors[-1], 1e-7)
        self.assertGreater(errors[0] / errors[1], 3.5)

    def test_native_nonorthogonal_basis_needs_only_retained_vectors(self):
        rng = np.random.default_rng(792)
        n = 12
        nr = 4
        a = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
        o = herm(a.conj().T @ a + np.eye(n))
        eo, uo = eigh(o)
        c = (uo / np.sqrt(eo)) @ uo.conj().T @ unitary(rng, n)
        e = np.r_[-2.0, -0.2, 0.3, 0.8, np.linspace(3.0, 9.0, n - nr)]
        h = herm(o @ c @ np.diag(e) @ c.conj().T @ o)
        basis = LAPWComplementBasis.prepare(
            h, o, c[:, :nr], e[:nr], mu=0.1, temperature=0.03
        )
        full = np.column_stack((basis.retained, basis.complement))
        np.testing.assert_allclose(full.conj().T @ o @ full, np.eye(n), atol=2e-14)
        np.testing.assert_allclose(
            basis.complement.conj().T @ h @ basis.retained, 0, atol=3e-13
        )
        dh = herm(rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n)))
        do = herm(rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n)))
        hr, hq = basis.project_matrix(dh)
        orr, oq = basis.project_matrix(do)
        tangent = basis.solver.action(hq - oq * e[:nr])
        expected = c[:, nr:] @ (
            (c[:, nr:].conj().T @ (dh @ c[:, :nr] - (do @ c[:, :nr]) * e[:nr]))
            / (e[:nr] - e[nr:, None])
        )
        np.testing.assert_allclose(basis.complement @ tangent, expected, atol=3e-14)
        raw = rng.normal(size=(n, 2, 1)) + 1j * rng.normal(size=(n, 2, 1))
        pr, pq = basis.project_raw(raw)
        np.testing.assert_allclose(
            pr, np.einsum("an,nic->aic", c[:, :nr].conj().T, raw)
        )
        rho = herm(rng.normal(size=(nr, nr)) + 1j * rng.normal(size=(nr, nr)))
        np.testing.assert_allclose(
            np.trace(o @ basis.density(rho, tangent)).real,
            np.trace(rho).real,
            atol=4e-14,
        )

    def test_rejects_occupied_or_degenerate_complement(self):
        for h in [np.diag([0.1, 2.0]), np.diag([1.0, 2.0])]:
            with self.assertRaises(ValueError):
                ShiftedComplement(h, np.array([-0.2, 1.0]), mu=0.1, temperature=0.01)
        with self.assertRaisesRegex(ValueError, "thermally empty"):
            ShiftedComplement(
                np.diag([1.1, 2.0]), np.array([-0.2, 1.0]), mu=1.15, temperature=0.1
            )
        with self.assertRaisesRegex(ValueError, "O-orthonormal"):
            LAPWComplementBasis.prepare(
                np.diag([-1.0, 2.0, 3.0]),
                np.eye(3),
                np.array([[2.0], [0.0], [0.0]]),
                np.array([-1.0]),
                mu=0.0,
                temperature=0.01,
            )
        with self.assertRaisesRegex(ValueError, "inconsistent"):
            LAPWComplementBasis.prepare(
                np.diag([-1.0, 2.0, 3.0]),
                np.eye(3),
                np.array([[1.0], [0.0], [0.0]]),
                np.array([-0.9]),
                mu=0.0,
                temperature=0.01,
            )


if __name__ == "__main__":
    unittest.main()
