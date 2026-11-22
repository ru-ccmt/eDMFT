"""Differentials of the two retained spaces in a second-variational SOC calculation.

A displacement and an input-density direction use the same maps. The scalar
eigenvectors A are in LAPW coordinates with metric O. The spinors C are in
two copies of that retained scalar space, in LAPWSO order (down, up).
All adjoints use the real Frobenius product Re Tr(X^H Y).
"""

import numpy as np
from scipy.linalg import cholesky, eigh, lstsq, qr, solve_triangular


def hermitian(x):
    return (x + x.conj().T) * 0.5


class RetainedEigenspace:
    """Parallel-transport differential of H A = O A E, including its metric.

    Rotations within an exactly degenerate retained multiplet are a gauge.
    Only retained/excluded differences occur in the denominators. A cutoff
    splitting a multiplet gives a nonunique spectral derivative; the caller
    can inspect singular_pairs to report the omissions.
    """

    def __init__(self, h, o, a, energy, *, gap_tolerance=1e-10):
        self.a = np.asarray(a)
        self.energy = np.asarray(energy)
        if (
            self.a.ndim != 2
            or self.a.shape != (len(h), len(self.energy))
            or not len(self.energy)
        ):
            raise ValueError("invalid retained eigenvector dimensions")
        if not np.isfinite(gap_tolerance) or gap_tolerance < 0:
            raise ValueError("invalid scalar/spinor gap tolerance")
        lower = None if o is None else cholesky(o, lower=True)
        y = self.a if lower is None else lower.conj().T @ self.a
        oa = self.a if o is None else o @ self.a
        residual = h @ self.a - oa * self.energy
        if (
            np.max(abs(y.conj().T @ y - np.eye(len(self.energy)))) > 1e-7
            or np.max(abs(residual)) > 1e-7
        ):
            raise ValueError("eigenvectors, eigenvalues and central H/O disagree")
        if y.shape[0] == len(self.energy):
            # A complete retained space has no cross-window band pairs.
            # Some LAPACK drivers reject an eigensolve of dimension zero.
            self.outside = self.a[:, :0]
            ev = self.energy[:0]
        else:
            q = qr(y, mode="full", check_finite=False)[0][:, len(self.energy) :]
            if lower is not None:
                q = solve_triangular(lower.conj().T, q, lower=False)
            ev, rotation = eigh(hermitian(q.conj().T @ h @ q))
            self.outside = q @ rotation
        gaps = self.energy[None, :] - ev[:, None]
        pairs = abs(gaps) <= gap_tolerance
        self.singular_pairs = int(np.count_nonzero(pairs))
        self.minimum_gap = np.min(abs(gaps), initial=np.inf)
        self.reciprocal_gaps = np.divide(
            1.0, gaps, out=np.zeros_like(gaps), where=~pairs
        )

    def action(self, dh, do=None):
        """Return d(A^H H A) and dA. The first matrix can be off diagonal."""
        a, b, e = self.a, self.outside, self.energy
        ha = dh @ a
        hrr = hermitian(a.conj().T @ ha)
        if do is None:
            return hrr, b @ ((b.conj().T @ ha) * self.reciprocal_gaps)
        oa = do @ a
        orr = hermitian(a.conj().T @ oa)
        # O-orthonormality fixes the Hermitian part of A^H O dA to -orr/2.
        da = b @ ((b.conj().T @ (ha - oa * e)) * self.reciprocal_gaps) - 0.5 * a @ orr
        return hrr - 0.5 * (e[:, None] + e) * orr, da

    def pullback(self, projected_bar, vector_bar, *, metric=True):
        a, b, e = self.a, self.outside, self.energy
        g = hermitian(projected_bar)
        outside_bar = b @ ((b.conj().T @ vector_bar) * self.reciprocal_gaps)
        bh = hermitian((a @ g + outside_bar) @ a.conj().T)
        if not metric:
            return bh
        bo_rr = -0.5 * ((e[:, None] + e) * g + hermitian(a.conj().T @ vector_bar))
        return bh, hermitian((a @ bo_rr - outside_bar * e) @ a.conj().T)


def spinor_coefficients(a, native):
    """Recover both scalar-band amplitudes from the two DMFT2 vector records."""
    if native.ndim != 3 or native.shape[2] != 2 or len(native) != len(a):
        raise ValueError("expected two components in the central LAPW ordering")
    columns = np.column_stack((native[:, :, 1], native[:, :, 0]))
    amplitudes, _, rank, _ = lstsq(a, columns, lapack_driver="gelsy")
    if rank != a.shape[1] or np.max(abs(a @ amplitudes - columns)) > 1e-8:
        raise ValueError("spinors do not lie in the supplied scalar input space")
    c = np.concatenate(np.hsplit(amplitudes, 2))
    if np.max(abs(c.conj().T @ c - np.eye(c.shape[1]))) > 1e-7:
        raise ValueError("spinor components are not jointly orthonormal")
    return c


class SecondVariationalResponse:
    """One k point, from scalar H/O and the SOC operator to physical spinors.

    The output wavefunctions (down,up) are in original LAPW coordinates.
    Their differential includes both scalar and spinor subspace motion. A
    downstream density or local overlap must differentiate both factors.
    """

    def __init__(self, h, o, scalar, energies, soc, spinors, spin_energies):
        self.scalar = RetainedEigenspace(h, o, scalar, energies)
        self.soc = soc
        hs = soc.action(central=True)
        hs[np.diag_indices_from(hs)] += np.tile(energies, 2)
        self.spin = RetainedEigenspace(hs, None, spinors, spin_energies)
        self.spin_blocks = np.split(spinors, 2)
        self.waves = np.stack([scalar @ c for c in self.spin_blocks])

    def action(self, dh, do, radial=None, velocities=None):
        hs, da = self.scalar.action(dh, do)
        hso = self.soc.action(radial, da, velocities)
        n = len(hs)
        hso[:n, :n] += hs
        hso[n:, n:] += hs
        hk, dc = self.spin.action(hso)
        dpsi = [
            da @ c + self.scalar.a @ d
            for c, d in zip(self.spin_blocks, np.split(dc, 2))
        ]
        return hk, np.stack(dpsi)

    def pullback(self, hbar, psibar):
        a = self.scalar.a
        da = sum(y @ c.conj().T for y, c in zip(psibar, self.spin_blocks))
        dc = np.concatenate([a.conj().T @ y for y in psibar])
        hk = self.spin.pullback(hbar, dc, metric=False)
        radial, matching, velocity = self.soc.pullback(hk)
        n = a.shape[1]
        dh, do = self.scalar.pullback(hk[:n, :n] + hk[n:, n:], da + matching)
        return dh, do, radial, velocity
