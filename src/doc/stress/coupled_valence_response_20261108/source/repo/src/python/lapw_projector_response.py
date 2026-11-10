"""Explicit raw projector-5 radial/matching response in one fixed local frame.

Window selection, orbital symmetry and inverse-square-root normalization are
handled by CachedFrozenSigmaResponse, after all k/site raw derivatives have
been assembled. No density or k weights belong in this raw overlap map.
"""

import numpy as np
from lapw_density_response import checked


class PreparedProjector5:
    """One correlated site/l/k, with supplied native local axes and site phase.

    orbital_transform[m,a] equals CF^T for a single non-SOC orbital block.
    For clusters, supply the corresponding block of the full transformation
    and accumulate its output in the correct impurity columns. Fixed projector
    tails may be supplied in the SAME band/orbital gauge. At fixed geometry
    their derivative with respect to the basis-generating potential is zero.
    This does not supply geometric/moving-boundary derivatives of those tails.
    """

    def __init__(
        self,
        channel,
        primitive,
        lo_coefficients,
        projector_metric,
        projector_radial,
        orbital_transform,
        *,
        lapw,
        tail=None
    ):
        self.channel = channel
        self.lapw = bool(lapw)
        array = np.asarray(primitive)
        if array.ndim != 3:
            raise ValueError("invalid projector primitive dimensions")
        n = len(array)
        m = 2 + channel.nlo
        self.phi = np.array(
            checked(array, (n, 4, m), "projector primitive", real=True),
            dtype=float,
            copy=True,
        )
        self.lo = np.array(
            checked(lo_coefficients, (m, channel.nlo), "LO factors", real=True),
            copy=True,
        )
        w = checked(projector_metric, (n, 4), "projector metric", real=True)
        p = checked(projector_radial, (n, 2), "fixed radial projector", real=True)
        t = np.asarray(orbital_transform)
        if t.ndim != 2 or t.shape[1] < 1:
            raise ValueError("invalid orbital transformation")
        self.transform = np.array(
            checked(t, (2 * channel.ell + 1, t.shape[1]), "orbital transformation"),
            complex,
        )
        self.weighted_projector = w[:, :2] * p
        self.q = np.einsum("rci,rc->i", self.phi[:, :2], self.weighted_projector)
        self.z = channel.value(self.phi[-1], self.lo, lapw=lapw)
        self.shape = (channel.nband, t.shape[1])
        self.tail = (
            np.zeros(self.shape, complex)
            if tail is None
            else np.array(checked(tail, self.shape, "fixed projector tail"), complex)
        )

    def value(self):
        return (
            np.einsum("min,i,ma->na", self.z.conj(), self.q, self.transform) + self.tail
        )

    def action(self, primitive_dot, lo_dot):
        u = checked(
            primitive_dot, self.phi.shape, "projector primitive derivative", real=True
        )
        dm = checked(lo_dot, self.lo.shape, "LO derivative", real=True)
        dz = self.channel.action(self.phi[-1], u[-1], dm, lapw=self.lapw)
        dq = np.einsum("rci,rc->i", u[:, :2], self.weighted_projector)
        return np.einsum("min,i,ma->na", dz.conj(), self.q, self.transform) + np.einsum(
            "min,i,ma->na", self.z.conj(), dq, self.transform
        )

    def pullback(self, covector):
        y = checked(covector, self.shape, "raw-projector covector")
        zb = np.einsum("i,ma,na->min", self.q, self.transform, y.conj())
        qb = np.einsum("na,min,ma->i", y.conj(), self.z.conj(), self.transform).real
        gp = np.zeros_like(self.phi)
        gp[:, :2] = self.weighted_projector[:, :, None] * qb
        gb, gm = self.channel.pullback(self.phi[-1], zb, lapw=self.lapw)
        gp[-1] += gb
        return gp, gm
