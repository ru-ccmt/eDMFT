"""Small primitive-space nonspherical potential and band response.

Uses the supplied native quadrature and harmonic convention. Both radial
legs, the potential and matching coefficients are differentiated. All ordered
l/l' pairs are kept, so the transpose never guesses a factor of two.
"""

import numpy as np
from lapw_augmentation_response import RadialAngularLayout
from lapw_density_response import checked, herm


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

    def action(self, primitive_dot, lo_dot, potential_dot=None):
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
            a = dz.conj().T @ kz
            out.append(herm(a + a.conj().T + z.conj().T @ dk @ z))
        return out

    def pullback(self, covectors):
        if len(covectors) != len(self.entries):
            raise ValueError("missing nonspherical covector")
        gp = {l: np.zeros_like(p) for l, p in self.phi.items()}
        gm = {l: np.zeros_like(self.entries[0][2][l]) for l in self.phi}
        gk = np.zeros_like(self.kernel)
        for y, (z, channels, lo, lapw, kz) in zip(covectors, self.entries):
            checked(y, (z.shape[1],) * 2, "nonspherical band covector", hermitian=True)
            zb = 2 * kz @ y
            gk += z @ y @ z.conj().T
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
