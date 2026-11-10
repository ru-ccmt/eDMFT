"""In-memory moving-basis part of the spatial valence charge response.

This supplies D_spherical: delta(r*V_sph) -> delta n at fixed band density.
The full D = D_spherical S, where S selects/converts the spherical radial
potential from the common KS-potential coordinates. Its transpose must use
S^T as well. This module does not infer units or freeze missing components of
the electronic response. R E B, core and charge normalization are separate.
"""

import numpy as np
from lapw_density_response import checked
from stress_response_adjoint import PairedOperator


class BasisDensityResponse:
    """Share one native radial context per atom/l across all k points/sites.

    spatial: prepared SpatialValenceDensity with weighted central rho;
    contexts[atom][l]: NativeRadialContext;
    matching[atom][entry][l]: PreparedAugmentationChannel, ordered exactly
        like spatial.atoms[atom].entries (equivalent site, k point).

    Input is one r*delta V_spherical_Ry mesh per inequivalent atom, packed in
    sorted atom order. Grid/RMT, coefficient gauge and linearization energies
    are fixed. No physical strain geometry is supplied by this operation.
    Automatic energy searches need a separate differentiated energy policy.
    The caller owns contexts and must keep them open for this object's life.
    """

    def __init__(self, spatial, contexts, matching, *, energy_policy):
        if energy_policy != "fixed":
            raise NotImplementedError(
                "only explicitly fixed radial energies are supported"
            )
        self.spatial, self.contexts, self.matching = spatial, contexts, matching
        if set(contexts) != set(spatial.atoms) or set(matching) != set(spatial.atoms):
            raise ValueError("missing or extra atom contexts")
        self.slices = {}
        offset = 0
        for a in sorted(spatial.atoms):
            atom = spatial.atoms[a]
            kernel = atom.kernel
            self.slices[a] = slice(offset, offset + kernel.nr)
            offset += kernel.nr
            if set(contexts[a]) != set(kernel.phi) or len(matching[a]) != len(
                atom.entries
            ):
                raise ValueError("radial or matching channel layout mismatch")
            for l, ctx in contexts[a].items():
                if (
                    ctx.ell != l
                    or ctx.n != kernel.nr
                    or not np.allclose(
                        ctx.primitive[:, :2], kernel.phi[l], rtol=1e-11, atol=1e-12
                    )
                ):
                    raise ValueError(
                        "central radial context does not match density basis"
                    )
            for channels, (_, z, _) in zip(matching[a], atom.entries):
                if set(channels) != set(contexts[a]):
                    raise ValueError("missing matching channel")
                for l, ch in channels.items():
                    ctx = contexts[a][l]
                    if ch.ell != l or ch.nlo != ctx.nlo or ch.rmt != ctx.rmt:
                        raise ValueError("matching geometry/radial channel mismatch")
                    central = ch.value(
                        ctx.primitive[-1], ctx.lo_coefficients, lapw=ctx.lapw
                    )
                    expected = z[kernel.layout.slices[l]].reshape(central.shape)
                    if not np.allclose(central, expected, rtol=1e-11, atol=1e-12):
                        raise ValueError(
                            "matching coefficients are in a different central gauge"
                        )
        self.input_size = offset
        self.operator = PairedOperator(
            offset,
            len(spatial.template.vector()),
            self.action,
            self.pullback,
            "spherical-potential moving-basis density",
        )

    def pack(self, potentials):
        if set(potentials) != set(self.slices):
            raise ValueError("missing or extra spherical potential")
        out = np.empty(self.input_size)
        for a, sl in self.slices.items():
            out[sl] = checked(
                potentials[a], (sl.stop - sl.start,), "spherical potential", real=True
            )
        return out

    def action(self, vector):
        v = checked(vector, (self.input_size,), "spherical potential vector", real=True)
        dz = {}
        dp = {}
        for a, sl in self.slices.items():
            # Radial/ABC work is independent of k and equivalent site. Prepare
            # once, then reuse its boundary values in the small matching sums.
            tangent = {l: ctx.action(v[sl]) for l, ctx in self.contexts[a].items()}
            dp[a] = {l: p[:, :2] for l, (p, m) in tangent.items()}
            dz[a] = []
            for channels in self.matching[a]:
                block = []
                for l in sorted(channels):
                    ctx = self.contexts[a][l]
                    p, m = tangent[l]
                    ch = channels[l]
                    block.append(
                        ch.action(ctx.primitive[-1], p[-1], m, lapw=ctx.lapw).reshape(
                            -1, ch.nband
                        )
                    )
                dz[a].append(np.concatenate(block))
        return self.spatial.action(matching_dots=dz, primitive_dots=dp)

    def pullback(self, covector):
        # Density-matrix and interstitial covectors are unnecessary for D.
        _, zbar, pbar = self.spatial.pullback(covector, include_density=False)
        result = np.zeros(self.input_size)
        for a, sl in self.slices.items():
            atom = self.spatial.atoms[a]
            gp = {
                l: np.zeros_like(ctx.primitive) for l, ctx in self.contexts[a].items()
            }
            gm = {
                l: np.zeros_like(ctx.lo_coefficients)
                for l, ctx in self.contexts[a].items()
            }
            for l in gp:
                gp[l][:, :2] = pbar[a][l]
            for channels, bar in zip(self.matching[a], zbar[a]):
                for l, ch in channels.items():
                    ctx = self.contexts[a][l]
                    b = bar[atom.kernel.layout.slices[l]].reshape(
                        2 * l + 1, 2 + ctx.nlo, ch.nband
                    )
                    gb, gl = ch.pullback(ctx.primitive[-1], b, lapw=ctx.lapw)
                    gp[l][-1] += gb
                    gm[l] += gl
            # Sum all matching and spatial radial covectors BEFORE the native
            # reverse pass. Discard energy covectors only under fixed E policy.
            for l, ctx in self.contexts[a].items():
                result[sl] += ctx.pullback(gp[l], gm[l])[0]
        return result
