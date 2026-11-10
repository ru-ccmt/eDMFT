"""Small-rank completion for successful native radial-energy searches.

The ordinary response holds radial energies fixed. SELECT's searched energies
are extra scalar inputs e=A V_sph. Compute each independent energy column once:
P_i = N [R E B_energy,i + D_energy,i], g_i = g_lattice.B_energy,i.
Then J_total=J_fixed+P A K and c_total=c_fixed+K^T A^T g. The source and
explicit scalar get the matching P A V_exp and g.A V_exp partners. This
retains their cancellations and adds no independent fitted stress term.
"""

import numpy as np
from response_parallel import parallel_of
from lapw_density_response import checked
from lapw_radial_energy_policy import searched_energy_derivatives
from stress_response_adjoint import PairedOperator, ScreenedResponse


def energy_direction(valence, energies):
    """Joint band and direct spatial derivatives at fixed physical potential."""
    if hasattr(valence, "branches"):
        pieces = [energy_direction(v, energies) for v in valence.branches]
        return valence.bands.output.pack(*(x[0] for x in pieces)), .5 * sum(x[1] for x in pieces)
    op = valence.bands.spherical
    c = op.coordinates
    nk, nr = c.raw_shape[:2]
    t = {}
    pt = {}
    for (a, l), de in energies.items():
        ctx = op.contexts[a][l]
        de = checked(de, (1 + ctx.nlo,), "radial energy direction", real=True)
        t[a, l] = ctx.action(np.zeros(ctx.n), de[0], de[1:])
        ctx = op.projector_contexts[a][l]
        pt[a, l] = ctx.action(np.zeros(ctx.n), de[0], de[1:])
    h = np.zeros((nk, nr, nr), complex)
    o = np.zeros_like(h)
    hq = [np.zeros(shape, complex) for shape in c.matrix.shapes]
    oq = [np.zeros_like(v) for v in hq]
    raw = np.zeros(c.raw_shape, complex)
    for a, l, k, s in op.spherical:
        if (a, l) not in energies:
            continue
        de = energies[a, l]
        dh, do = s.action(*t[a, l], energy_dot=de[0], lo_energies_dot=de[1:])
        h[k] += dh[:nr]
        hq[k] += dh[nr:]
        o[k] += do[:nr]
        oq[k] += do[nr:]
    for a, ns, ks in op.nonspherical:
        if not any(aa == a for aa, l in energies):
            continue
        dp = {
            l: t[a, l][0] if (a, l) in t else np.zeros_like(op.contexts[a][l].primitive)
            for l in ns.phi
        }
        dm = {
            l: (
                t[a, l][1]
                if (a, l) in t
                else np.zeros_like(op.contexts[a][l].lo_coefficients)
            )
            for l in ns.phi
        }
        for k, dh in zip(ks, ns.action(dp, dm, retained=nr)):
            h[k] += dh[:nr]
            hq[k] += dh[nr:]
    for a, l, k, ic, cols, p in op.projectors:
        if (a, l) in pt:
            raw[k, :, cols, ic] += p.action(*pt[a, l])[:nr]
    direct = valence.direct
    dp = {a: None for a in direct.contexts}
    dz = {a: None for a in direct.contexts}
    for a in direct.contexts:
        if not any(aa == a for aa, l in energies):
            continue
        dp[a] = {
            l: (
                pt[a, l][0][:, :2]
                if (a, l) in pt
                else np.zeros_like(ctx.primitive[:, :2])
            )
            for l, ctx in direct.contexts[a].items()
        }
        dz[a] = []
        for channels in direct.matching[a]:
            blocks = []
            for l, ch in sorted(channels.items()):
                ctx = direct.contexts[a][l]
                if (a, l) in pt:
                    p, m = pt[a, l]
                    block = ch.action(ctx.primitive[-1], p[-1], m, lapw=ctx.lapw)
                else:
                    block = np.zeros((2 * l + 1, 2 + ctx.nlo, ch.nband), complex)
                blocks.append(block.reshape(-1, ch.nband))
            dz[a].append(np.concatenate(blocks))
    density = direct.spatial.action(matching_dots=dz, primitive_dots=dp)
    return c.pack(h, o, raw, hq, oq), density


class SearchedEnergyResponse:
    """Attach searched radial energies, holding the IN1 reference EF fixed.

    Fixed RMT means no separate explicit geometry of the radial root. A
    changing IN1 reference EF, a root/fallback switch, HDLO, or moving RMT
    requires a different policy and is not silently supplied here.
    """

    def __init__(
        self,
        state,
        observable,
        directory,
        case,
        *,
        complex_vectors,
        library,
        progress=None
    ):
        parallel = parallel_of(state)
        maps, details = searched_energy_derivatives(
            directory, case, complex_vectors=complex_vectors, library=library
        )
        self.state = state
        self.observable = observable
        self.details = details
        self.A = np.zeros((len(details), state.chain.potential.output_size))
        columns = []
        direct = []
        for i, item in enumerate(details):
            a, l = item["atom"], item["l"]
            slots = item["slots"]
            mat = maps[a, l]
            sl = next(
                sl
                for key, sl in zip(
                    state.potential_core.hartree.keys,
                    state.potential_core.hartree.slices,
                )
                if key == (a, 0, 0)
            )
            self.A[i, sl] = mat[slots[0]]
            de = np.zeros(len(mat))
            de[slots] = 1
            band, density = energy_direction(state.valence, {(a, l): de})
            ch = state.chain
            columns.append(
                ch.normalization.matvec(
                    ch.density.matvec(ch.electronic.matvec(band))
                    + parallel.sum(density)
                )
            )
            direct.append(parallel.sum(observable.lattice.joint_gradient @ band))
            if progress is not None:
                progress(i, item)
        self.P = (
            np.column_stack(columns)
            if columns
            else np.zeros((state.chain.normalization.output_size, 0))
        )
        self.g = np.array(direct)
        self.K = state.chain.potential
        self.J = state.jacobian
        # Factor on the common potential input: K (especially PBE) is
        # expensive and must be applied only once per forward/reverse action.
        fixed = state.chain.charge_from_potential
        self.charge_from_potential = PairedOperator(
            fixed.input_size,
            fixed.output_size,
            lambda v: fixed.matvec(v) + self.P @ (self.A @ v),
            lambda y: fixed.rmatvec(y) + self.A.T @ (self.P.T @ y),
            "charge response including searched radial energies",
        )
        self.jacobian = self.charge_from_potential @ self.K
        self.screened = ScreenedResponse(self.jacobian)
        self.charge_gradient = observable.charge_gradient + self.K.rmatvec(
            self.A.T @ self.g
        )

    def complete_source(self, fixed_source, explicit_potential):
        return fixed_source + self.P @ (self.A @ explicit_potential)

    def explicit_addition(self, explicit_potential):
        return self.g @ (self.A @ explicit_potential)

    def total_energy_directions(self, explicit_potential, density_response):
        return self.A @ (explicit_potential + self.K.matvec(density_response))
