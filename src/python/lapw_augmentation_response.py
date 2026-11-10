"""Prepared muffin-tin angular assembly for the native response engine.

Plane-wave sums are done at preparation, not inside a response iteration.
The functions here keep the central band gauge fixed. They assemble potential,
overlap and raw-projector blocks; the spherical kinetic/boundary Hamiltonian
and its native triangle convention remain separate contributions.
"""

from functools import lru_cache
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import sph_harm_y, spherical_jn


def real_pair(a, b):
    return float(np.vdot(a, b).real)


def harmonics(ell, q_local):
    """Condon--Shortley Y_lm, ordered m=-l,...,l, for local Cartesian q."""
    q = np.asarray(q_local, float)
    if q.ndim != 2 or q.shape[1] != 3 or not np.all(np.isfinite(q)):
        raise ValueError("invalid local reciprocal vectors")
    length = np.linalg.norm(q, axis=1)
    cosine = np.divide(q[:, 2], length, out=np.ones_like(length), where=length > 0)
    theta, phi = np.arccos(np.clip(cosine, -1, 1)), np.arctan2(q[:, 1], q[:, 0])
    return np.array([sph_harm_y(ell, m, theta, phi) for m in range(-ell, ell + 1)])


class PreparedAugmentationChannel:
    """Cache the fixed geometry sums for one atom/l/k and central band gauge.

    ``q_local`` supplies directions in the atom's native local frame;
    ``q_length`` supplies global Cartesian lengths used for Bessel matching.
    Keeping these separate preserves the native treatment of rounded ROTLOC.
    ``phase`` is the native site phase, chosen consistently with the raw
    projector gauge. ``lo_groups[j]`` identifies all coefficient rows belonging
    to ordinary LO slot j of this atom/l, including equivalent atoms.

    F=sum_G geometry*j_l*C_G, D=sum_G geometry*q*j_l'*C_G,
    L_j=sum_(G in LO_j) geometry*C_G, geometry=4*pi*i^l*Y_lm^*phase/sqrt(Omega).
    Complete LO factors carry the additional native RMT^2 prefactor.
    The prepared object retains only these small band arrays.
    A native caller can supply its existing (j_l, q*j_l') arrays through
    matching_values; this preserves the host's numerical matching convention.
    """

    def __init__(
        self,
        ell,
        rmt,
        volume,
        q_local,
        q_length,
        phase,
        coefficients,
        plane_wave_count,
        lo_groups,
        *,
        matching_values=None,
    ):
        q = np.asarray(q_length, float)
        c, phase = np.asarray(coefficients), np.asarray(phase)
        npw = int(plane_wave_count)
        if ell < 0 or int(ell) != ell or rmt <= 0 or volume <= 0:
            raise ValueError("invalid augmentation geometry")
        if c.ndim != 2 or not 0 <= npw <= c.shape[0] or c.shape[1] < 1:
            raise ValueError("invalid central coefficient dimensions")
        if q.shape != (c.shape[0],) or phase.shape != q.shape or np.any(q < 0):
            raise ValueError("incompatible reciprocal lengths or phases")
        if not all(np.all(np.isfinite(v)) for v in (c, q, phase)):
            raise ValueError("nonfinite augmentation input")
        y = harmonics(ell, q_local)
        if y.shape[1] != len(q):
            raise ValueError("local directions and coefficient rows differ")
        geometry = 4 * np.pi / np.sqrt(volume) * (1j**ell) * y.conj() * phase
        if matching_values is None:
            j = spherical_jn(ell, q[:npw] * rmt)
            dj = q[:npw] * spherical_jn(ell, q[:npw] * rmt, derivative=True)
        else:
            matching = np.asarray(matching_values)
            if (
                matching.shape != (2, npw)
                or np.iscomplexobj(matching)
                or not np.all(np.isfinite(matching))
            ):
                raise ValueError("invalid native matching values")
            j, dj = matching
        self.fj = (geometry[:, :npw] * j) @ c[:npw]
        self.dfj = (geometry[:, :npw] * dj) @ c[:npw]
        seen = set()
        local = []
        for group in lo_groups:
            indices = np.asarray(group)
            if (
                indices.ndim != 1
                or not np.issubdtype(indices.dtype, np.integer)
                or np.any(indices < npw)
                or np.any(indices >= len(q))
            ):
                raise ValueError("invalid native LO coefficient rows")
            if len(set(indices)) != len(indices) or seen.intersection(indices):
                raise ValueError("duplicate LO coefficient rows")
            seen.update(indices)
            local.append(geometry[:, indices] @ c[indices])
        self.local = np.array(local, dtype=complex).reshape(
            (len(local), 2 * ell + 1, c.shape[1])
        )
        self.ell, self.rmt = int(ell), float(rmt)
        self.nlo, self.nband = len(local), c.shape[1]

    def _check(self, boundary, lo):
        if np.shape(boundary) != (4, 2 + self.nlo) or np.shape(lo) != (
            2 + self.nlo,
            self.nlo,
        ):
            raise ValueError("radial primitive/LO dimensions disagree")
        if np.iscomplexobj(boundary) or np.iscomplexobj(lo):
            raise ValueError("radial components must be real")
        if not np.all(np.isfinite(boundary)) or not np.all(np.isfinite(lo)):
            raise ValueError("nonfinite radial components")

    def value(self, boundary, lo, *, lapw):
        """Return Z[m,primitive,band], for primitive boundary (A,B,AP,BP)."""
        self._check(boundary, lo)
        r = self.rmt
        z = r * r * np.einsum("is,smn->min", lo, self.local)
        if lapw:
            z[:, 0] += r * (boundary[0, 1] * self.dfj - boundary[2, 1] * self.fj)
            z[:, 1] += r * (boundary[2, 0] * self.fj - boundary[0, 0] * self.dfj)
        else:
            if boundary[0, 0] == 0:
                raise ValueError("singular APW boundary value")
            z[:, 0] += r * self.fj / boundary[0, 0]
        return z

    def action(self, boundary, boundary_dot, lo_dot, *, lapw):
        self._check(boundary_dot, lo_dot)
        self._check(boundary, lo_dot)
        r = self.rmt
        z = r * r * np.einsum("is,smn->min", lo_dot, self.local)
        if lapw:
            z[:, 0] += r * (
                boundary_dot[0, 1] * self.dfj - boundary_dot[2, 1] * self.fj
            )
            z[:, 1] += r * (
                boundary_dot[2, 0] * self.fj - boundary_dot[0, 0] * self.dfj
            )
        else:
            if boundary[0, 0] == 0:
                raise ValueError("singular APW boundary value")
            z[:, 0] -= r * self.fj * boundary_dot[0, 0] / boundary[0, 0] ** 2
        return z

    def pullback(self, boundary, covector, *, lapw):
        """Real boundary and ABC covectors, to sum before radial reverse solve."""
        self._check(boundary, np.zeros((2 + self.nlo, self.nlo)))
        g = np.asarray(covector)
        if g.shape != (2 * self.ell + 1, 2 + self.nlo, self.nband) or not np.all(
            np.isfinite(g)
        ):
            raise ValueError("invalid augmentation covector")
        r = self.rmt
        gb = np.zeros_like(boundary, dtype=float)
        gm = r * r * np.einsum("min,smn->is", g.conj(), self.local).real
        if lapw:
            gb[0, 1] = r * real_pair(g[:, 0], self.dfj)
            gb[2, 1] = -r * real_pair(g[:, 0], self.fj)
            gb[2, 0] = r * real_pair(g[:, 1], self.fj)
            gb[0, 0] = -r * real_pair(g[:, 1], self.dfj)
        else:
            if boundary[0, 0] == 0:
                raise ValueError("singular APW boundary value")
            gb[0, 0] = -r * real_pair(g[:, 0], self.fj) / boundary[0, 0] ** 2
        return gb, gm


@lru_cache(maxsize=128)
def real_harmonic_gaunt(lmax, signed_l, m):
    """Gaunt matrix for a noncubic WIEN real lattice harmonic.

    Positive signed_l selects cosine, negative selects sine. The caller must
    expand cubic harmonics first. Numerical product quadrature is exact for
    the finite harmonic degrees here; it is performed only during preparation.
    """
    L = abs(signed_l)
    if min(lmax, L, m) < 0 or m > L or (signed_l < 0 and m == 0):
        raise ValueError("invalid real lattice harmonic")
    degree = 2 * lmax + L
    nx, nphi = max(2, degree + 2), max(3, 2 * degree + 3)
    x, weights = leggauss(nx)
    phi = 2 * np.pi * np.arange(nphi) / nphi
    sin_theta = np.sqrt(1 - x * x)
    q = (
        np.array(
            [
                sin_theta[:, None] * np.cos(phi),
                sin_theta[:, None] * np.sin(phi),
                np.broadcast_to(x[:, None], (nx, nphi)),
            ]
        )
        .reshape(3, -1)
        .T
    )
    y = np.concatenate([harmonics(l, q) for l in range(max(lmax, L) + 1)])
    potential = y[L * (L + 1) + m]
    if m:
        negative = y[L * (L + 1) - m]
        potential = (
            (((-1) ** m) * potential + negative) / np.sqrt(2)
            if signed_l >= 0
            else (-1j * (((-1) ** m) * potential - negative) / np.sqrt(2))
        )
    y = y[: (lmax + 1) ** 2]
    result = (
        y.conj() * (np.repeat(weights, nphi) * 2 * np.pi / nphi * potential)
    ) @ y.T
    result.setflags(write=False)
    return result


class AugmentationContractions:
    """Small-space V/O/raw-projector assembly and its exact reverse.

    Z[component,band] describes central band wavefunctions in MT primitives.
    K and S are potential and overlap kernels; P[component,orbital] contains
    fixed-frame local-orbital overlaps. Return (Z^dagger K Z, Z^dagger S Z,
    Z^dagger P). P and the result are RAW: window selection, symmetry and
    inverse-square-root normalization belong to the shared electronic map.
    This potential block does not replace spherical kinetic/surface H.
    """

    def __init__(self, z, potential, overlap, projector):
        self.z = np.array(z, complex, copy=True)
        self.k = np.array(potential, complex, copy=True)
        self.s = np.array(overlap, complex, copy=True)
        self.p = np.array(projector, complex, copy=True)
        if self.z.ndim != 2 or min(self.z.shape) < 1:
            raise ValueError("invalid augmentation band coefficients")
        n = self.z.shape[0]
        if (
            self.k.shape != (n, n)
            or self.s.shape != (n, n)
            or self.p.ndim != 2
            or self.p.shape[0] != n
        ):
            raise ValueError("incompatible augmentation kernels")
        if not all(np.all(np.isfinite(a)) for a in (self.z, self.k, self.s, self.p)):
            raise ValueError("nonfinite augmentation kernels")
        for a in (self.k, self.s):
            if np.linalg.norm(a - a.conj().T) > 1e-11 * max(1.0, np.linalg.norm(a)):
                raise ValueError("augmentation kernel is not Hermitian")
        self.kz, self.sz = self.k @ self.z, self.s @ self.z

    def value(self):
        return (
            self.z.conj().T @ self.kz,
            self.z.conj().T @ self.sz,
            self.z.conj().T @ self.p,
        )

    def action(self, z_dot, k_dot, s_dot, p_dot):
        z_dot, k_dot, s_dot, p_dot = map(np.asarray, (z_dot, k_dot, s_dot, p_dot))
        for value, base in zip(
            (z_dot, k_dot, s_dot, p_dot), (self.z, self.k, self.s, self.p)
        ):
            if np.shape(value) != base.shape or not np.all(np.isfinite(value)):
                raise ValueError("invalid augmentation derivative")
        for a in (k_dot, s_dot):
            if np.linalg.norm(a - a.conj().T) > 1e-11 * max(1.0, np.linalg.norm(a)):
                raise ValueError("augmentation kernel derivative must be Hermitian")
        a, b = z_dot.conj().T @ self.kz, z_dot.conj().T @ self.sz
        return (
            a + a.conj().T + self.z.conj().T @ k_dot @ self.z,
            b + b.conj().T + self.z.conj().T @ s_dot @ self.z,
            z_dot.conj().T @ self.p + self.z.conj().T @ p_dot,
        )

    def pullback(self, v_bar, o_bar, p_bar):
        nb = self.z.shape[1]
        for a in (v_bar, o_bar):
            if np.shape(a) != (nb, nb) or not np.all(np.isfinite(a)):
                raise ValueError("invalid band covector")
            if np.linalg.norm(a - a.conj().T) > 1e-11 * max(1.0, np.linalg.norm(a)):
                raise ValueError("band covector must be Hermitian")
        if np.shape(p_bar) != (nb, self.p.shape[1]) or not np.all(np.isfinite(p_bar)):
            raise ValueError("invalid projector covector")
        z_bar = 2 * self.kz @ v_bar + 2 * self.sz @ o_bar + self.p @ p_bar.conj().T
        return (
            z_bar,
            self.z @ v_bar @ self.z.conj().T,
            self.z @ o_bar @ self.z.conj().T,
            self.z @ p_bar,
        )


class RadialAngularLayout:
    """Linear map between radial integrals and the small angular kernel.

    Components are ordered by l, then m=-l..l, then radial primitive. All
    ordered (l,l') pairs are supplied, including their Hermitian partners.
    The transpose uses the same full-matrix pairing; no extra factor of two
    is inserted by guessing which pairs are independent.
    """

    def __init__(self, radial_counts, harmonics, *, gaunt_matrices=None):
        self.counts = dict(sorted(radial_counts.items()))
        if not self.counts or any(
            int(l) != l or l < 0 or int(n) != n or n < 2 for l, n in self.counts.items()
        ):
            raise ValueError("invalid primitive radial counts")
        self.harmonics = tuple(tuple(h) for h in harmonics)
        self.slices = {}
        size = 0
        for l, n in self.counts.items():
            self.slices[l] = slice(size, size + (2 * l + 1) * n)
            size += (2 * l + 1) * n
        self.size = size
        lmax = max(self.counts)
        # Native exports may provide their exact Gaunt convention, including
        # already expanded cubic harmonics. Never infer a cubic combination
        # from the signed-L convention used for noncubic real harmonics.
        if gaunt_matrices is None:
            angular = np.array(
                [real_harmonic_gaunt(lmax, L, M) for L, M in self.harmonics]
            ).reshape(len(self.harmonics), (lmax + 1) ** 2, (lmax + 1) ** 2)
        else:
            angular = np.asarray(gaunt_matrices)
            if angular.shape != (
                len(self.harmonics),
                (lmax + 1) ** 2,
                (lmax + 1) ** 2,
            ) or not np.all(np.isfinite(angular)):
                raise ValueError("invalid native Gaunt matrices")
        self.gaunt = {}
        for l in self.counts:
            for lp in self.counts:
                self.gaunt[l, lp] = np.array(
                    angular[:, l * l : (l + 1) ** 2, lp * lp : (lp + 1) ** 2], copy=True
                )

    def potential(self, radial_pairs):
        if set(radial_pairs) != set(self.gaunt):
            raise ValueError("missing or extra angular momentum pairs")
        out = np.zeros((self.size, self.size), complex)
        for (l, lp), g in self.gaunt.items():
            value = np.asarray(radial_pairs[l, lp])
            if (
                value.shape != (self.counts[l], self.counts[lp], len(self.harmonics))
                or np.iscomplexobj(value)
                or not np.all(np.isfinite(value))
            ):
                raise ValueError("invalid primitive potential integrals")
            block = np.einsum("qmn,ijq->minj", g, value)
            out[self.slices[l], self.slices[lp]] = block.reshape(
                (2 * l + 1) * self.counts[l], (2 * lp + 1) * self.counts[lp]
            )
        return out

    def potential_pullback(self, covector):
        g = np.asarray(covector)
        if g.shape != (self.size, self.size) or not np.all(np.isfinite(g)):
            raise ValueError("invalid angular potential covector")
        out = {}
        for (l, lp), angular in self.gaunt.items():
            block = g[self.slices[l], self.slices[lp]].reshape(
                2 * l + 1, self.counts[l], 2 * lp + 1, self.counts[lp]
            )
            out[l, lp] = np.einsum("qmn,minj->ijq", angular.conj(), block).real
        return out

    def overlap(self, radial_metrics):
        if set(radial_metrics) != set(self.counts):
            raise ValueError("missing or extra overlap channel")
        out = np.zeros((self.size, self.size), complex)
        for l, n in self.counts.items():
            a = np.asarray(radial_metrics[l])
            if a.shape != (n, n) or np.iscomplexobj(a) or not np.all(np.isfinite(a)):
                raise ValueError("invalid radial overlap")
            out[self.slices[l], self.slices[l]] = np.kron(np.eye(2 * l + 1), a)
        return out

    def overlap_pullback(self, covector):
        g = np.asarray(covector)
        if g.shape != (self.size, self.size) or not np.all(np.isfinite(g)):
            raise ValueError("invalid overlap covector")
        out = {}
        for l, n in self.counts.items():
            block = g[self.slices[l], self.slices[l]].reshape(
                2 * l + 1, n, 2 * l + 1, n
            )
            out[l] = np.einsum("mimj->ij", block).real
        return out

    def projector(self, ell, radial_overlap, local_transform):
        """Raw fixed-frame projector kernel; normalization is applied later."""
        q, t = np.asarray(radial_overlap), np.asarray(local_transform)
        if (
            ell not in self.counts
            or q.shape != (self.counts[ell],)
            or np.iscomplexobj(q)
        ):
            raise ValueError("invalid radial projector overlap")
        if (
            t.ndim != 2
            or t.shape[0] != 2 * ell + 1
            or not np.all(np.isfinite(t))
            or not np.all(np.isfinite(q))
        ):
            raise ValueError("invalid local orbital transformation")
        out = np.zeros((self.size, t.shape[1]), complex)
        out[self.slices[ell]] = np.einsum("i,ma->mia", q, t).reshape(-1, t.shape[1])
        return out

    def projector_pullback(self, ell, radial_overlap, local_transform, covector):
        self.projector(ell, radial_overlap, local_transform)  # validate fixed inputs
        q, t, g = (
            np.asarray(radial_overlap),
            np.asarray(local_transform),
            np.asarray(covector),
        )
        if g.shape != (self.size, t.shape[1]) or not np.all(np.isfinite(g)):
            raise ValueError("invalid raw-projector covector")
        block = g[self.slices[ell]].reshape(2 * ell + 1, len(q), t.shape[1])
        return np.einsum("ma,mia->i", t.conj(), block).real, np.einsum(
            "i,mia->ma", q, block
        )
