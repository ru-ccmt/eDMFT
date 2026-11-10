"""Fixed-basis integral, transpose, and ordinary potential scalar controls."""

from types import SimpleNamespace as NS
import unittest
import numpy as np
from scipy.special import expit
from lapw_physical_potential import FixedBasisPhysicalPotential
from lapw_complement_response import (
    ComplementMatrixCoordinates,
    ComplementPerturbationCoordinates,
)
from lapw_spherical_response import PreparedSphericalBand
from lapw_nonspherical_response import PreparedNonspherical
from lapw_potential_response import (
    NativePotentialCoordinates,
    NativeInterstitialPotential,
)
from frozen_sigma_response import herm
from test_lapw_potential_assembly import fixture


class PhysicalPotentialTests(unittest.TestCase):
    def setUp(self):
        rng, ch, p, lo, weight, radius, energy, elo, triangle, _ = fixture()
        nr = 2
        sphere = PreparedSphericalBand(
            ch, p, lo, weight, radius, energy, elo, triangle[:, :nr], lapw=True
        )
        z = ch.value(p[-1], lo, lapw=True).reshape(-1, ch.nband)
        ns = PreparedNonspherical(
            {1: p},
            weight,
            [(2, 0), (-2, 1)],
            rng.normal(size=(len(p), 2)),
            [(z, {1: ch}, {1: lo}, {1: True})],
        )
        labels = np.array([[0, 0, 0], [1, 0, 0]])
        reps = np.array([[0, 0, 0], [1, 0, 0], [-1, 0, 0]])
        ir = NativeInterstitialPotential(
            rng.normal(size=(2, 4)) + 1j * rng.normal(size=(2, 4)),
            labels,
            reps,
            np.eye(3, dtype=int)[None],
            np.zeros((1, 3)),
            complex_values=True,
        )
        self.bands = NS(
            compact=True,
            nr=nr,
            nk=1,
            sizes=[4],
            input=NativePotentialCoordinates(
                len(p), {0: ns.v.shape}, 3, complex_interstitial=True
            ),
            output=ComplementPerturbationCoordinates(
                ComplementMatrixCoordinates((1, nr, nr), [2]), (1, nr, 1, 1)
            ),
            spherical=NS(
                spherical=[(0, 1, 0, sphere)],
                slices={0: slice(0, len(p))},
                input_size=len(p),
            ),
            ns={0: (ns, [0])},
            ir=[ir],
        )
        self.op = FixedBasisPhysicalPotential(self.bands)
        self.rng = rng
        self.sphere = sphere

    def test_spherical_matrix_is_independent_fixed_basis_integral(self):
        op = self.sphere
        dv = self.rng.normal(size=len(op.radius))
        x = self.bands.input.pack(
            dv, {0: np.zeros_like(self.bands.ns[0][0].v)}, np.zeros(3, complex)
        )
        h, o, raw, hq, oq = self.op.output.unpack(self.op.action(x))
        # Independently evaluate the native real-space quadratic integral.
        # No radial tangent, surface kernel, or energy derivative enters.
        integral = np.einsum(
            "rci,rc,rcj,r->ij", op.phi, op.metric, op.phi, dv / op.radius
        )
        expected = sum(z.conj().T @ integral @ z[:, :2] for z in op.z)
        np.testing.assert_allclose(np.concatenate((h[0], hq[0])), expected, atol=2e-11)
        np.testing.assert_array_equal(o, 0)
        np.testing.assert_array_equal(raw, 0)
        np.testing.assert_array_equal(oq, 0)

    def test_all_sectors_cross_blocks_and_transpose(self):
        x = self.rng.normal(size=self.op.input.size)
        y = self.rng.normal(size=self.op.output.size)
        forward = self.op.action(x)
        lhs, rhs = y @ forward, self.op.pullback(y) @ x
        self.assertLess(abs(lhs - rhs) / max(1, abs(lhs), abs(rhs)), 3e-13)

    def test_density_covector_matches_nonlinear_fermi_scalar(self):
        x = 0.001 * self.rng.normal(size=self.op.input.size)
        h, _, _, _, _ = self.op.output.unpack(self.op.action(x))
        center = herm(self.rng.normal(size=(2, 2)) + 1j * self.rng.normal(size=(2, 2)))
        e, u = np.linalg.eigh(center)
        temperature = 0.2
        rho = 2 * (u * expit(-e / temperature)) @ u.conj().T
        analytic = self.op.density_gradient(rho[None]) @ x
        scalar = (
            lambda t: -2
            * temperature
            * np.logaddexp(
                0, -np.linalg.eigvalsh(center + t * h[0]) / temperature
            ).sum()
        )
        step = 1e-5
        fd = (scalar(step) - scalar(-step)) / (2 * step)
        np.testing.assert_allclose(fd, analytic, rtol=2e-7, atol=2e-8)


if __name__ == "__main__":
    unittest.main()
