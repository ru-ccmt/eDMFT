"""Full moving-basis charge-map tests against independent nonlinear endpoints.

The small complex real-space model is a physical Galerkin discretization, not
a material calculation. Its potential changes both the Hamiltonian operator
and trial functions. The oracle rediagonalizes and renormalizes the projectors.
"""

import unittest
import numpy as np

from frozen_sigma_response import herm
from frozen_sigma_tangent import real_inner
from stress_response_adjoint import PairedOperator
from lapw_response_chain import ChargeResponseChain, generalized_electronic_operator
from test_frozen_sigma_tangent import complex_model
from test_generalized_sigma_response import endpoint


def dense_map(matrix, name):
    return PairedOperator(
        matrix.shape[1],
        matrix.shape[0],
        lambda x: matrix @ x,
        lambda y: matrix.T @ y,
        name,
    )


class SpatialFixture:
    def __init__(self, fixed_number=False):
        self.model = complex_model()
        self.electronic, self.coordinates = generalized_electronic_operator(
            self.model, fixed_number
        )
        self.ng = 12
        nk, nb = self.model.eps.shape
        rng = np.random.default_rng(196)
        self.psi = np.array(
            [
                np.linalg.qr(
                    rng.normal(size=(self.ng, nb)) + 1j * rng.normal(size=(self.ng, nb))
                )[0]
                for k in range(nk)
            ]
        )
        self.operator = np.array(
            [
                b @ np.diag(e) @ b.conj().T + 1.7 * (np.eye(self.ng) - b @ b.conj().T)
                for b, e in zip(self.psi, self.model.eps)
            ]
        )
        self.trial_change = 0.025 * (
            rng.normal(size=(self.ng, nk, self.ng, nb))
            + 1j * rng.normal(size=(self.ng, nk, self.ng, nb))
        )
        self.chi = []
        for k, b in enumerate(self.psi):
            outer = rng.normal(size=(self.ng, 2)) + 1j * rng.normal(size=(self.ng, 2))
            self.chi.append(
                (
                    b @ self.model.raw[k, :, :, 0]
                    + 0.15 * (np.eye(self.ng) - b @ b.conj().T) @ outer
                )[:, :, None]
            )
        self.chi = np.array(self.chi)
        self.rho = self.model.central_density()
        columns = []
        direct = []
        for g in range(self.ng):
            v = np.eye(self.ng)[g]
            db = self.trial_change[g]
            h = herm(
                db.swapaxes(-1, -2).conj() @ self.operator @ self.psi
                + self.psi.swapaxes(-1, -2).conj() @ self.operator @ db
                + self.psi.swapaxes(-1, -2).conj() @ (v[None, :, None] * self.psi)
            )
            o = herm(
                db.swapaxes(-1, -2).conj() @ self.psi
                + self.psi.swapaxes(-1, -2).conj() @ db
            )
            raw = np.einsum("kgn,kgic->knic", db.conj(), self.chi)
            columns.append(self.coordinates.pack(h, o, raw))
            direct.append(
                2 * np.einsum("kgn,knm,kgm->g", db, self.rho, self.psi.conj()).real
            )
        basis = dense_map(np.array(columns).T, "native basis fixture")
        density = PairedOperator(
            self.coordinates.band.size,
            self.ng,
            lambda x: np.einsum(
                "kgn,knm,kgm->g",
                self.psi,
                self.coordinates.band.unpack(x),
                self.psi.conj(),
            ).real,
            lambda y: self.coordinates.band.pack(
                np.einsum("kgn,g,kgm->knm", self.psi.conj(), y, self.psi)
            ),
            "spatial density fixture",
        )
        a = rng.normal(size=(self.ng, self.ng))
        self.kernel = 0.03 * (a @ a.T) / self.ng
        self.chain = ChargeResponseChain(
            dense_map(self.kernel, "Hartree/XC fixture"),
            basis,
            self.electronic,
            density,
            dense_map(np.array(direct).T, "direct basis density"),
            dense_map(np.zeros((self.ng, self.ng)), "no-core fixture"),
            dense_map(np.eye(self.ng), "unit normalization fixture"),
        )

    def nonlinear(self, v):
        b = self.psi + np.einsum("p,pkgn->kgn", v, self.trial_change)
        bd = b.swapaxes(-1, -2).conj()
        h = herm(bd @ (self.operator + np.diag(v)[None]) @ b)
        o = herm(bd @ b)
        raw = np.einsum("kgn,kgic->knic", b.conj(), self.chi)
        rho, _, _ = endpoint(
            self.model,
            h - np.array([np.diag(e) for e in self.model.eps]),
            o - np.eye(h.shape[-1]),
            raw - self.model.raw,
            1.0,
        )
        return np.einsum("kgn,knm,kgm->g", b, rho, b.conj()).real


class ChainTests(unittest.TestCase):
    def test_joint_coordinates_pairing(self):
        f = SpatialFixture()
        rng = np.random.default_rng(31)
        x = rng.normal(size=f.coordinates.size)
        y = rng.normal(size=x.size)
        self.assertAlmostEqual(
            sum(
                real_inner(a, b)
                for a, b in zip(f.coordinates.unpack(x), f.coordinates.unpack(y))
            ),
            x @ y,
            places=12,
        )
        np.testing.assert_allclose(
            f.coordinates.pack(*f.coordinates.unpack(x)), x, atol=1e-15
        )

    def test_complete_spatial_nonlinear_response(self):
        f = SpatialFixture()
        v = np.sin(np.arange(f.ng) + 0.2)
        actual = f.chain.charge_from_potential.matvec(v)
        errors = []
        for h in (2e-4, 1e-4):
            fd = (f.nonlinear(h * v) - f.nonlinear(-h * v)) / (2 * h)
            errors.append(np.linalg.norm(fd - actual))
        self.assertLess(errors[-1], 2e-7)
        self.assertGreater(errors[0] / errors[1], 3.8)
        incomplete = f.chain.valence.matvec(v)
        self.assertGreater(np.linalg.norm(incomplete - actual), 0.05)

    def test_full_chain_transpose_both_ensembles(self):
        for fixed_number in (False, True):
            f = SpatialFixture(fixed_number)
            errors = f.chain.check_transposes(samples=2)
            self.assertLess(max(errors.values()), 1e-12)
            if fixed_number:
                actual = f.chain.charge_from_potential.matvec(np.sin(np.arange(f.ng)))
                self.assertLess(abs(actual.sum()), 1e-13)

    def test_screened_chain_matches_dense_reference_and_one_adjoint(self):
        f = SpatialFixture(True)
        c = f.chain
        rng = np.random.default_rng(4)
        jac = np.column_stack([c.jacobian.matvec(e) for e in np.eye(f.ng)])
        sources = rng.normal(size=(f.ng, 6))
        gradient = rng.normal(size=f.ng)
        direct = rng.normal(size=6)
        expected = direct + gradient @ np.linalg.solve(np.eye(f.ng) - jac, sources)
        actual, _, report = c.screened.contractions(
            gradient, sources, direct, rtol=1e-12, atol=1e-13
        )
        self.assertTrue(report["converged"])
        np.testing.assert_allclose(actual, expected, atol=1e-11)

    def test_explicit_source_includes_each_path_once(self):
        f = SpatialFixture()
        c = f.chain
        rng = np.random.default_rng(51)
        band = rng.normal(size=f.coordinates.size)
        basis = rng.normal(size=f.ng)
        core = rng.normal(size=f.ng)
        pot = rng.normal(size=f.ng)
        actual = c.explicit_source(band, basis, core, pot)
        expected = c.density.matvec(c.electronic.matvec(band + c.basis.matvec(pot)))
        expected += basis + core + c.basis_density.matvec(pot) + c.core.matvec(pot)
        np.testing.assert_allclose(actual, expected, atol=2e-14)

    def test_missing_coordinate_adapter_is_rejected(self):
        f = SpatialFixture()
        c = f.chain
        with self.assertRaises(ValueError):
            ChargeResponseChain(
                c.potential,
                c.basis,
                c.electronic,
                c.density,
                c.basis_density,
                dense_map(np.zeros((f.ng, 1)), "wrong core coordinates"),
                c.normalization,
            )


if __name__ == "__main__":
    unittest.main()
