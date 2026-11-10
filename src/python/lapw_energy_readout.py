"""Diagnostic pullbacks to native radial linearization energies.

These are covectors, not a prescription for the strain derivative of an
energy-search algorithm. They expose the energy inputs discarded under the
fixed-energy policy; both native matrix and ATPARN projector/density banks
are retained. All coefficients are Ry per Ry of the corresponding energy.
"""

import numpy as np


def band_energy_gradient(spherical, vector):
    hb, ob, rb, hq, oq = spherical.coordinates.unpack(vector)
    nr = rb.shape[1]
    hc = [np.concatenate((h, q)) for h, q in zip(hb, hq)]
    oc = [np.concatenate((o, q)) for o, q in zip(ob, oq)]
    gp = {
        a: {l: np.zeros_like(c.primitive) for l, c in cs.items()}
        for a, cs in spherical.contexts.items()
    }
    gm = {
        a: {l: np.zeros_like(c.lo_coefficients) for l, c in cs.items()}
        for a, cs in spherical.contexts.items()
    }
    ge = {
        (a, l): np.zeros(1 + c.nlo)
        for a, cs in spherical.contexts.items()
        for l, c in cs.items()
    }
    for a, l, k, op in spherical.spherical:
        p, m, e, el, _ = op.pullback(hc[k], oc[k])
        gp[a][l] += p
        gm[a][l] += m
        ge[a, l] += np.r_[e, el]
    for a, op, ks in spherical.nonspherical:
        p, m, _ = op.pullback([hc[k] for k in ks])
        for l in op.phi:
            gp[a][l] += p[l]
            gm[a][l] += m[l]
    pg = {a: {} for a in spherical.contexts}
    mg = {a: {} for a in spherical.contexts}
    for a, l in {(a, l) for a, l, *_ in spherical.projectors}:
        c = spherical.projector_contexts[a][l]
        pg[a][l] = np.zeros_like(c.primitive)
        mg[a][l] = np.zeros_like(c.lo_coefficients)
    for a, l, k, ic, cols, op in spherical.projectors:
        y = np.zeros(op.shape, complex)
        y[:nr] = rb[k, :, cols, ic]
        p, m = op.pullback(y)
        pg[a][l] += p
        mg[a][l] += m
    for a, cs in spherical.contexts.items():
        for l, c in cs.items():
            _, e, el = c.pullback(gp[a][l], gm[a][l])
            ge[a, l] += np.r_[e, el]
    for a, cs in pg.items():
        for l in cs:
            _, e, el = spherical.projector_contexts[a][l].pullback(pg[a][l], mg[a][l])
            ge[a, l] += np.r_[e, el]
    return ge


def density_energy_gradient(direct, covector):
    _, zbar, pbar = direct.spatial.pullback(covector, include_density=False)
    result = {}
    for a in direct.slices:
        atom = direct.spatial.atoms[a]
        gp = {l: np.zeros_like(c.primitive) for l, c in direct.contexts[a].items()}
        gm = {
            l: np.zeros_like(c.lo_coefficients) for l, c in direct.contexts[a].items()
        }
        for l in gp:
            gp[l][:, :2] = pbar[a][l]
        for channels, bar in zip(direct.matching[a], zbar[a]):
            for l, ch in channels.items():
                c = direct.contexts[a][l]
                b = bar[atom.kernel.layout.slices[l]].reshape(
                    2 * l + 1, 2 + c.nlo, ch.nband
                )
                gb, gl = ch.pullback(c.primitive[-1], b, lapw=c.lapw)
                gp[l][-1] += gb
                gm[l] += gl
        for l, c in direct.contexts[a].items():
            _, e, el = c.pullback(gp[l], gm[l])
            result[a, l] = np.r_[e, el]
    return result


def screened_energy_gradient(state, observable, adjoint):
    """dGamma/dE, including the induced self-consistent density via one lambda.

    E denotes the independent native radial-energy inputs here. The result
    contracts the TOTAL E_nu from any separately specified energy-search
    policy. It must not be added if that policy already entered B and J.
    """
    ch = state.chain
    normalized = ch.normalization.rmatvec(adjoint)
    joint = observable.lattice.joint_gradient + ch.electronic.rmatvec(
        ch.density.rmatvec(normalized)
    )
    ge = band_energy_gradient(state.valence.bands.spherical, joint)
    for key, value in density_energy_gradient(state.valence.direct, normalized).items():
        ge[key] += value
    return ge
