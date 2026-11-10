"""Joint H/O/raw-projector derivative for spherical muffin-tin perturbations.

This is the spherical-input restriction of B, not a complete potential map.
The nonspherical CENTRAL Hamiltonian is retained and its moving radial legs
are differentiated. Interstitial geometry and potential are held fixed.
No normalization/window derivative belongs here: E performs that chain once.
"""

import numpy as np
from lapw_density_response import checked
from lapw_response_chain import GeneralizedCoordinates
from stress_response_adjoint import PairedOperator


class SphericalPotentialBands:
    """Assemble shared radial tangents and accumulate covectors before reverse.

    contexts[atom][l] are NativeRadialContext objects owned by the caller.
    spherical entries: (atom,l,k,PreparedSphericalBand), including each site.
    nonspherical entries: (atom,PreparedNonspherical,k_indices_for_entries).
    projectors: (atom,l,k,impurity_block,column_slice,PreparedProjector5).
    raw_shape is the LOCAL impurity-block layout of the electronic response.
    Potentials are packed spherical r*delta V_Ry, one mesh per inequivalent
    atom in sorted order. All radial energies are explicitly held fixed.
    """

    def __init__(
        self,
        contexts,
        spherical,
        nonspherical,
        projectors,
        band_shape,
        raw_shape,
        *,
        energy_policy
    ):
        if energy_policy != "fixed":
            raise NotImplementedError("only fixed radial energies are supported")
        self.contexts = contexts
        self.spherical = list(spherical)
        self.nonspherical = list(nonspherical)
        self.projectors = list(projectors)
        self.coordinates = GeneralizedCoordinates(band_shape, raw_shape)
        self.slices = {}
        offset = 0
        for a in sorted(contexts):
            if not contexts[a]:
                raise ValueError("empty radial context set")
            n = next(iter(contexts[a].values())).n
            if any(c.n != n for c in contexts[a].values()):
                raise ValueError("inconsistent atom grids")
            self.slices[a] = slice(offset, offset + n)
            offset += n
        nk, nb, _ = band_shape
        if raw_shape[:2] != (nk, nb):
            raise ValueError("band/projector dimensions disagree")
        covered = set()
        for a, l, k, op in self.spherical:
            c = contexts[a][l]
            if (
                not 0 <= k < nk
                or op.triangle.shape != (nb, nb)
                or not np.array_equal(op.phi, c.primitive)
            ):
                raise ValueError("inconsistent spherical context")
            covered.add((a, l, k))
        required = {(a, l, k) for a in contexts for l in contexts[a] for k in range(nk)}
        if covered != required:
            raise ValueError("missing spherical atom/l/k sectors")
        if {a for a, op, ks in self.nonspherical} != set(contexts):
            raise ValueError("missing nonspherical atom")
        for a, op, ks in self.nonspherical:
            if len(ks) != len(op.entries) or any(not 0 <= k < nk for k in ks):
                raise ValueError("nonspherical entry/k mismatch")
            for l, p in op.phi.items():
                if not np.array_equal(p, contexts[a][l].primitive):
                    raise ValueError("inconsistent nonspherical primitive")
        for a, l, k, ic, cols, op in self.projectors:
            if (
                not 0 <= k < nk
                or not 0 <= ic < raw_shape[3]
                or op.shape != (nb, len(range(raw_shape[2])[cols]))
            ):
                raise ValueError("projector output layout mismatch")
            if not np.array_equal(op.phi, contexts[a][l].primitive):
                raise ValueError("inconsistent projector primitive")
        self.input_size = offset
        self.operator = PairedOperator(
            offset,
            self.coordinates.size,
            self.action,
            self.pullback,
            "spherical potential to joint H/O/raw projector",
        )

    def pack(self, potentials):
        if set(potentials) != set(self.slices):
            raise ValueError("missing spherical atom potential")
        v = np.empty(self.input_size)
        for a, sl in self.slices.items():
            v[sl] = checked(
                potentials[a], (sl.stop - sl.start,), "spherical potential", real=True
            )
        return v

    def tangents(self, vector):
        v = checked(
            vector, (self.input_size,), "spherical potential coordinates", real=True
        )
        return {
            a: {l: ctx.action(v[self.slices[a]]) for l, ctx in channels.items()}
            for a, channels in self.contexts.items()
        }

    def action(self, vector):
        tangent = self.tangents(vector)
        h = np.zeros(self.coordinates.band.shape, complex)
        o = np.zeros_like(h)
        raw = np.zeros(self.coordinates.raw_shape, complex)
        for a, l, k, op in self.spherical:
            dh, do = op.action(*tangent[a][l])
            h[k] += dh
            o[k] += do
        for a, op, ks in self.nonspherical:
            dp = {l: tangent[a][l][0] for l in op.phi}
            dm = {l: tangent[a][l][1] for l in op.phi}
            for k, dh in zip(ks, op.action(dp, dm)):
                h[k] += dh
        for a, l, k, ic, cols, op in self.projectors:
            raw[k, :, cols, ic] += op.action(*tangent[a][l])
        return self.coordinates.pack(h, o, raw)

    def pullback(self, vector):
        hb, ob, rb = self.coordinates.unpack(vector)
        gp = {
            a: {l: np.zeros_like(c.primitive) for l, c in cs.items()}
            for a, cs in self.contexts.items()
        }
        gm = {
            a: {l: np.zeros_like(c.lo_coefficients) for l, c in cs.items()}
            for a, cs in self.contexts.items()
        }
        for a, l, k, op in self.spherical:
            p, m, eb, elb, db = op.pullback(hb[k], ob[k])
            gp[a][l] += p
            gm[a][l] += m
            # Fixed E and Vphysical=Vbasis: explicit dE and d(Vphys-Vbasis)
            # are zero. Their covectors have no route to this input space.
        for a, op, ks in self.nonspherical:
            p, m, v = op.pullback([hb[k] for k in ks])
            for l in op.phi:
                gp[a][l] += p[l]
                gm[a][l] += m[l]
        for a, l, k, ic, cols, op in self.projectors:
            p, m = op.pullback(rb[k, :, cols, ic])
            gp[a][l] += p
            gm[a][l] += m
        answer = np.zeros(self.input_size)
        for a, channels in self.contexts.items():
            for l, c in channels.items():
                answer[self.slices[a]] += c.pullback(gp[a][l], gm[a][l])[0]
        return answer
