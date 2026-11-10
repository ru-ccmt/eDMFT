"""Native Fourier phases/triangle, rectangular kernels and coordinate controls."""

import unittest
import numpy as np
from lapw_density_response import herm
from lapw_potential_response import (
    NativeInterstitialPotential,
    NativePotentialCoordinates,
)
from lapw_nonspherical_response import PreparedNonspherical
from lapw_augmentation_response import PreparedAugmentationChannel


class GeneralPotentialTests(unittest.TestCase):
    def test_native_ir_against_explicit_star_and_triangle(self):
        rng = np.random.default_rng(92)
        g = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0]])
        reps = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0], [1, -1, 0]])
        rotations = np.array([np.eye(3, dtype=int), -np.eye(3, dtype=int)])
        translations = np.array([[0.0, 0.0, 0.0], [0.0, 0.25, 0.0]])
        c = rng.normal(size=(4, 5)) + 1j * rng.normal(size=(4, 5))
        for complex_values in [False, True]:
            op = NativeInterstitialPotential(
                c, g, reps, rotations, translations, complex_values=complex_values
            )
            v = rng.normal(size=len(reps))
            if complex_values:
                v = v + 1j * rng.normal(size=len(reps))
            # Independent nested native assembly, including distinct star
            # members, duplicate-operation phase average and upper copying.
            lookup = {}
            for coeff, rep in zip(v, reps):
                members = {}
                for r, t in zip(rotations, translations):
                    angle = 2 * np.pi * (rep @ t)
                    factor = np.exp(-1j * angle) if complex_values else np.cos(angle)
                    members.setdefault(tuple(r @ rep), []).append(factor)
                for key, factors in members.items():
                    lookup[key] = coeff * sum(factors) / len(factors)
            native = np.zeros((4, 4), complex)
            for i in range(4):
                for j in range(i, 4):
                    # Complex TAPEWF exports CONJG(Z); the matrix diagnostic
                    # conjugates HPANEL to match those coefficient vectors.
                    z = lookup.get(tuple(g[i] - g[j]), 0j).conjugate()
                    native[i, j] = z.real if i == j else z
                    native[j, i] = native[i, j].conjugate()
            np.testing.assert_allclose(
                op.action(v), c.conj().T @ native @ c, atol=3e-14
            )
            for nr in [1, 3]:
                rr, q = op.action(v, retained=nr)
                y = herm(rng.normal(size=(nr, nr)) + 1j * rng.normal(size=(nr, nr)))
                yq = rng.normal(size=q.shape) + 1j * rng.normal(size=q.shape)
                vb = op.pullback(y, yq)
                self.assertAlmostEqual(
                    np.vdot(y, rr).real + 2 * np.vdot(yq, q).real,
                    np.vdot(vb, v).real,
                    places=11,
                )
            y = herm(rng.normal(size=(5, 5)) + 1j * rng.normal(size=(5, 5)))
            self.assertAlmostEqual(
                np.vdot(y, op.action(v)).real,
                np.vdot(op.pullback(y), v).real,
                places=11,
            )

    def test_complex_native_vector_convention(self):
        # A real potential with imaginary Fourier coefficients distinguishes
        # the solver eigenvector convention from the exported case.vector one.
        op = NativeInterstitialPotential(
            np.eye(2),
            np.array([[0, 0, 0], [1, 0, 0]]),
            np.array([[0, 0, 0], [-1, 0, 0], [1, 0, 0]]),
            np.eye(3, dtype=int)[None],
            np.zeros((1, 3)),
            complex_values=True,
        )
        np.testing.assert_allclose(
            op.coefficient_matrix(np.array([0.7, 1 + 0.3j, 1 - 0.3j])),
            np.array([[0.7, 1 - 0.3j], [1 + 0.3j, 0.7]]),
        )

    def test_ir_duplicate_phases_are_averaged_not_summed(self):
        g = np.array([[0, 0, 0], [1, 0, 0]])
        rotation = np.tile(np.eye(3, dtype=int), (2, 1, 1))
        tau = np.array([[0, 0, 0], [0.5, 0, 0.0]])
        op = NativeInterstitialPotential(
            np.eye(2),
            g,
            np.array([[0, 0, 0], [-1, 0, 0]]),
            rotation,
            tau,
            complex_values=True,
        )
        np.testing.assert_allclose(
            op.action(np.array([0.7, 1.0 + 0.3j])), 0.7 * np.eye(2), atol=1e-15
        )

    def test_nonspherical_rectangular_potential_and_reverse(self):
        rng = np.random.default_rng(98)
        n = 11
        nb = 7
        phi = {0: rng.normal(size=(n, 4, 2)), 1: rng.normal(size=(n, 4, 2))}
        metric = np.zeros((n, 4))
        metric[:, :2] = rng.uniform(0.1, 0.5, (n, 2))
        q = rng.normal(size=(8, 3))
        qn = np.linalg.norm(q, axis=1)
        coeff = rng.normal(size=(8, nb)) + 1j * rng.normal(size=(8, nb))
        channels = {
            l: PreparedAugmentationChannel(
                l, 2.0, 40.0, q, qn, np.exp(0.2j * q[:, 0]), coeff, 8, []
            )
            for l in phi
        }
        lo = {l: np.zeros((2, 0)) for l in phi}
        lapw = {l: True for l in phi}
        z = np.concatenate(
            [
                channels[l].value(phi[l][-1], lo[l], lapw=True).reshape(-1, nb)
                for l in phi
            ]
        )
        v = rng.normal(size=(n, 2))
        op = PreparedNonspherical(
            phi, metric, [(1, 0), (-1, 1)], v, [(z, channels, lo, lapw)]
        )
        dv = rng.normal(size=v.shape)
        nr = 3
        old = op.action({l: np.zeros_like(p) for l, p in phi.items()}, lo, dv)[0]
        full = op.potential_action(dv)[0]
        rr, qr = op.potential_action(dv, retained=nr)[0]
        np.testing.assert_allclose(full, old, atol=5e-13)
        np.testing.assert_allclose(rr, full[:nr, :nr], atol=5e-13)
        np.testing.assert_allclose(qr, full[nr:, :nr], atol=5e-13)
        y = herm(rng.normal(size=(nr, nr)) + 1j * rng.normal(size=(nr, nr)))
        yq = rng.normal(size=qr.shape) + 1j * rng.normal(size=qr.shape)
        vb = op.potential_pullback([y], [yq])
        self.assertAlmostEqual(
            np.vdot(y, rr).real + 2 * np.vdot(yq, qr).real,
            np.vdot(vb, dv).real,
            places=10,
        )

    def test_native_potential_coordinates_reality_and_roundtrip(self):
        for complex_values in [False, True]:
            coords = NativePotentialCoordinates(
                3, {0: (2, 2), 1: (3, 1)}, 4, complex_interstitial=complex_values
            )
            v = np.random.default_rng(19).normal(size=coords.size)
            np.testing.assert_array_equal(coords.pack(*coords.unpack(v)), v)
            s, ns, pw = coords.unpack(v)
            with self.assertRaises(ValueError):
                coords.pack(s, {0: ns[0]}, pw)
            if not complex_values:
                with self.assertRaises(ValueError):
                    coords.pack(s, ns, pw.astype(complex))

    def test_duplicate_stars_and_lo_rows_are_not_silently_accepted(self):
        g = np.array([[0, 0, 0], [1, 0, 0]])
        with self.assertRaisesRegex(ValueError, "overlapping"):
            NativeInterstitialPotential(
                np.eye(2),
                g,
                np.array([[1, 0, 0], [-1, 0, 0]]),
                np.array([np.eye(3, dtype=int), -np.eye(3, dtype=int)]),
                np.zeros((2, 3)),
                complex_values=True,
            )
        with self.assertRaisesRegex(ValueError, "duplicate PW"):
            NativeInterstitialPotential(
                np.eye(2),
                np.zeros((2, 3), int),
                g,
                np.eye(3, dtype=int)[None],
                np.zeros((1, 3)),
                complex_values=False,
            )


if __name__ == "__main__":
    unittest.main()
