"""Small primitive-space nonspherical potential and band response.

Uses the supplied native quadrature and harmonic convention. Both radial
legs, the potential and matching coefficients are differentiated. All ordered
l/l' pairs are kept, so the transpose never guesses a factor of two.
"""

import numpy as np
from lapw_augmentation_response import RadialAngularLayout
from lapw_density_response import checked, herm
from lapw_compact_blocks import (
    checked_columns,
    hermitian_columns,
    right_multiply_columns,
)


class PreparedNonspherical:
    def __init__(self, primitive, metric, harmonics, potential, matching):
        if not primitive:
            raise ValueError("empty nonspherical radial set")
        arrays = {l: np.asarray(p) for l, p in primitive.items()}
        if any(p.ndim != 3 for p in arrays.values()):
            raise ValueError("invalid nonspherical primitive dimensions")
        self.n = next(iter(arrays.values())).shape[0]
        self.phi = {
            l: np.array(
                checked(
                    p, (self.n, 4, p.shape[2]), "nonspherical primitive", real=True
                ),
                dtype=float,
                copy=True,
            )
            for l, p in arrays.items()
        }
        self.metric = np.array(
            checked(metric, (self.n, 4), "radial metric", real=True), copy=True
        )
        self.layout = RadialAngularLayout(
            {l: p.shape[2] for l, p in self.phi.items()}, harmonics
        )
        self.v = np.array(
            checked(
                potential,
                (self.n, len(harmonics)),
                "nonspherical potential (Ry)",
                real=True,
            ),
            copy=True,
        )
        # matching: (Z, [PreparedAugmentationChannel by l], LO factors by l,
        #             LAPW flags by l). Each entry can represent one site/k.
        self.entries = []
        self.pairs = self._integrate(self.phi, self.phi, self.v)
        self.kernel = self.layout.potential(self.pairs)
        for z, channels, lo, lapw in matching:
            if (
                set(channels) != set(self.phi)
                or set(lo) != set(self.phi)
                or set(lapw) != set(self.phi)
            ):
                raise ValueError("missing nonspherical channel")
            z = np.asarray(z)
            if z.ndim != 2:
                raise ValueError("invalid nonspherical matching dimensions")
            z = np.array(
                checked(z, (self.layout.size, z.shape[1]), "nonspherical matching"),
                complex,
                copy=True,
            )
            self.entries.append((z, channels, lo, lapw, self.kernel @ z))
        if not self.entries:
            raise ValueError("no nonspherical site/k entries")

    def _integrate(self, left, right, potential):
        out = {}
        for (l, lp), g in self.layout.gaunt.items():
            if np.max(abs(g), initial=0) < 5e-14:
                out[l, lp] = np.zeros(
                    (
                        self.layout.counts[l],
                        self.layout.counts[lp],
                        len(self.layout.harmonics),
                    )
                )
            else:
                out[l, lp] = np.einsum(
                    "rci,rc,rcj,rq->ijq",
                    left[l],
                    self.metric,
                    right[lp],
                    potential,
                    optimize=True,
                )
        return out

    def value(self):
        return [herm(z.conj().T @ kz) for z, ch, lo, lapw, kz in self.entries]

    def potential_action(self, potential_dot, retained=None):
        """Fixed-radial potential input, optionally ONLY RR and QR band blocks.

        This is the explicit nonspherical delta V leg. Changes of the central
        nonspherical kernel through radial/matching response belong to action().
        The potential has the same Ry/harmonic convention as self.v, without r.
        """
        dv = checked(
            potential_dot, self.v.shape, "nonspherical potential derivative", real=True
        )
        dk = self.layout.potential(self._integrate(self.phi, self.phi, dv))
        out = []
        for z, channels, lo, lapw, kz in self.entries:
            n = z.shape[1]
            if retained is not None and (
                int(retained) != retained or not 0 < retained < n
            ):
                raise ValueError("invalid retained nonspherical band count")
            nr = n if retained is None else int(retained)
            columns = z.conj().T @ (dk @ z[:, :nr])
            out.append(
                herm(columns)
                if retained is None
                else (herm(columns[:nr]), columns[nr:])
            )
        return out

    def potential_pullback(self, covectors, cross_covectors=None):
        """Transpose of potential_action; QR pairing includes its RQ partner.

        Sum all k/site covectors in the small angular space BEFORE the radial
        contractions. With cross_covectors, do not form the unused QQ block.
        """
        if len(covectors) != len(self.entries) or (
            cross_covectors is not None and len(cross_covectors) != len(self.entries)
        ):
            raise ValueError("missing nonspherical potential covectors")
        gk = np.zeros_like(self.kernel)
        for j, (z, channels, lo, lapw, kz) in enumerate(self.entries):
            y = np.asarray(covectors[j])
            if cross_covectors is None:
                checked(y, (z.shape[1],) * 2, "potential band covector", hermitian=True)
                gk += z @ y @ z.conj().T
            else:
                if y.ndim != 2 or not 0 < len(y) < z.shape[1]:
                    raise ValueError("invalid retained potential covector")
                nr = len(y)
                checked(y, (nr, nr), "retained potential covector", hermitian=True)
                q = checked(
                    cross_covectors[j],
                    (z.shape[1] - nr, nr),
                    "cross potential covector",
                )
                zr, zq = z[:, :nr], z[:, nr:]
                cross = (zq @ q) @ zr.conj().T
                gk += (zr @ y) @ zr.conj().T + cross + cross.conj().T
        gv = np.zeros_like(self.v)
        for (l, lp), bar in self.layout.potential_pullback(gk).items():
            if np.max(abs(self.layout.gaunt[l, lp]), initial=0) < 5e-14:
                continue
            gv += np.einsum(
                "rci,rc,rcj,ijq->rq",
                self.phi[l],
                self.metric,
                self.phi[lp],
                bar,
                optimize=True,
            )
        return gv

    def action(self, primitive_dot, lo_dot, potential_dot=None, *, retained=None):
        """Differentiate radial/matching legs; optionally return only RR/QR columns."""
        if set(primitive_dot) != set(self.phi) or set(lo_dot) != set(self.phi):
            raise ValueError("missing nonspherical tangents")
        for l, p in primitive_dot.items():
            checked(p, self.phi[l].shape, "nonspherical radial derivative", real=True)
        dv = (
            np.zeros_like(self.v)
            if potential_dot is None
            else checked(potential_dot, self.v.shape, "potential derivative", real=True)
        )
        a = self._integrate(primitive_dot, self.phi, self.v)
        b = self._integrate(self.phi, primitive_dot, self.v)
        c = self._integrate(self.phi, self.phi, dv)
        dk = self.layout.potential({key: a[key] + b[key] + c[key] for key in a})
        out = []
        for z, channels, lo, lapw, kz in self.entries:
            dz = np.concatenate(
                [
                    channels[l]
                    .action(
                        self.phi[l][-1], primitive_dot[l][-1], lo_dot[l], lapw=lapw[l]
                    )
                    .reshape(-1, z.shape[1])
                    for l in self.layout.counts
                ]
            )
            n = z.shape[1] if retained is None else int(retained)
            if not 0 < n <= z.shape[1] or (retained is not None and n != retained):
                raise ValueError("invalid retained nonspherical count")
            out.append(
                hermitian_columns(
                    dz.conj().T @ kz[:, :n]
                    + kz.conj().T @ dz[:, :n]
                    + z.conj().T @ (dk @ z[:, :n])
                )
            )
        return out

    def pullback(self, covectors):
        if len(covectors) != len(self.entries):
            raise ValueError("missing nonspherical covector")
        gp = {l: np.zeros_like(p) for l, p in self.phi.items()}
        gm = {l: np.zeros_like(self.entries[0][2][l]) for l in self.phi}
        gk = np.zeros_like(self.kernel)
        for y, (z, channels, lo, lapw, kz) in zip(covectors, self.entries):
            checked_columns(y, z.shape[1], "nonspherical band covector")
            zb = 2 * right_multiply_columns(kz, y)
            gk += right_multiply_columns(z, y) @ z.conj().T
            for l, ch in channels.items():
                bar = zb[self.layout.slices[l]].reshape(
                    2 * l + 1, self.layout.counts[l], z.shape[1]
                )
                gb, gl = ch.pullback(self.phi[l][-1], bar, lapw=lapw[l])
                gp[l][-1] += gb
                gm[l] += gl
        pairs = self.layout.potential_pullback(gk)
        gv = np.zeros_like(self.v)
        for (l, lp), bar in pairs.items():
            if np.max(abs(self.layout.gaunt[l, lp]), initial=0) < 5e-14:
                continue
            weighted = np.einsum("rq,ijq->rij", self.v, bar, optimize=True)
            gp[l] += self.metric[:, :, None] * np.einsum(
                "rij,rcj->rci", weighted, self.phi[lp], optimize=True
            )
            gp[lp] += self.metric[:, :, None] * np.einsum(
                "rij,rci->rcj", weighted, self.phi[l], optimize=True
            )
            gv += np.einsum(
                "rci,rc,rcj,ijq->rq",
                self.phi[l],
                self.metric,
                self.phi[lp],
                bar,
                optimize=True,
            )
        return gp, gm, gv
