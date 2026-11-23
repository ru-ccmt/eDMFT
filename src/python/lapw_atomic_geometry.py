"""Cartesian position gradient of the LAPW band and valence-density readouts.

At fixed cell, RMT and native field coefficients, a site enters the LAPW
matching sums through exp(i G.R). The projector uses exp(i (k+G).R).
Reverse the band/charge contractions into these phases before selecting a
displacement. Each local k point supplies a gradient for every active site.
No site-resolved dH, dO, or spatial charge derivative needs to be stored.
"""

from collections import defaultdict

import numpy as np
from scipy.special import eval_legendre, spherical_jn

from lapw_augmentation_response import harmonics
from lapw_compact_blocks import right_multiply_columns
from lapw_geometry_response import native_band_geometry_source, prepared_geometry_inputs
from lapw_strain_charge import native_spatial_geometry_source
from lapw_lattice import NativeLattice


def matching_position_gradient(
    ell,
    radius,
    volume,
    local_q,
    lengths,
    phase,
    rates,
    coefficients,
    npw,
    lo_rows,
    boundary,
    lo,
    lapw,
    zbar,
    pw_only=None,
):
    """Transpose the three Cartesian derivatives of an augmentation channel.

    For F = sum_g Y_g phase_g j_g C_g the real pairing is
    Re sum_mn conj(Fbar_mn) dF_mn. Multiplication by i*G when differentiating
    the phase gives -Im of the corresponding row weight, times G. For LAPW
    plane waves both the Bessel value and slope contribute. The additional
    pw_only covector belongs to the boundary value with LO coefficients zero.
    """
    angular = 4 * np.pi * (1j**ell) / np.sqrt(volume) * harmonics(ell, local_q).conj() * phase
    pw = zbar if pw_only is None else zbar + pw_only
    if lapw:
        fbar = radius * (-boundary[2, 1] * pw[:, 0] + boundary[2, 0] * pw[:, 1])
        dbar = radius * (boundary[0, 1] * pw[:, 0] - boundary[0, 0] * pw[:, 1])
    else:
        fbar = radius * pw[:, 0] / boundary[0, 0]
    length = lengths[:npw]
    c = coefficients[:npw].T
    j = spherical_jn(ell, radius * length)
    row = (fbar.conj() @ c) * j
    if lapw:
        slope = length * spherical_jn(ell, radius * length, derivative=True)
        row += (dbar.conj() @ c) * slope
    weight = -np.einsum("mg,mg->g", angular[:, :npw], row).imag
    gradient = weight @ rates[:npw]
    for slot, indices in enumerate(lo_rows):
        # The LO term carries RMT^2 and no Bessel factor.
        local_bar = radius**2 * np.einsum("min,i->mn", zbar, lo[:, slot])
        row = local_bar.conj() @ coefficients[indices].T
        gradient -= np.einsum("mg,mg->g", angular[:, indices], row).imag @ rates[indices]
    return gradient


def plane_wave_position_gradients(q, labels, c, structure, volume, weights, covectors, active):
    """Transpose interstitial O/T and the ordered MT Wronskian together.

    For an explicit or screened readout Y in RR/QR coordinates, move its
    Hermitian completion to the PW basis: B = C Y C†. Every translated
    site then contributes a scalar array times (q_i-q_j) to its gradient.
    Block the i index so scratch is independent of the number of displacements.
    """
    n = len(q)
    result = np.zeros((len(covectors), len(structure.pos), 3))
    if n == 0:
        return result
    length = np.linalg.norm(q, axis=1)
    unit = np.divide(q, length[:, None], out=np.zeros_like(q), where=length[:, None] > 0)
    phase = np.exp(2j * np.pi * (structure.pos @ labels.T))
    orbit = np.repeat(np.arange(structure.nat), structure.mult)
    orbit_data = []
    for a in np.unique(orbit[active]):
        sites = np.flatnonzero(active & (orbit == a))
        r = structure.rmt[a]
        channels = []
        for ell, w in weights[a].items():
            if w:
                channels.append(
                    (
                        ell,
                        4 * np.pi * (2 * ell + 1) * w / volume,
                        spherical_jn(ell, r * length),
                        length * spherical_jn(ell, r * length, derivative=True),
                    )
                )
        orbit_data.append((sites, r, channels))
    index = np.arange(n)
    rows = max(1, min(n, (8 * 1024**2) // (224 * n)))
    ket = c.conj().T
    for start in range(0, n, rows):
        block = slice(start, start + rows)
        metric = q[block] @ q.T
        # Direct differences also resolve the diagonal and very short G-G'.
        diff = np.linalg.norm(q[block, None] - q[None, :], axis=2)
        cosine = np.clip(unit[block] @ unit.T, -1, 1)
        order = np.sign(index[None, :] - index[block, None])
        bars = []
        for hbar, obar, _ in covectors:
            h = (right_multiply_columns(c[block], hbar) @ ket).conj()
            o = (right_multiply_columns(c[block], obar) @ ket).conj()
            bars.append((h, o + metric * h))
        for sites, r, channels in orbit_data:
            x = r * diff
            sphere = (
                4
                * np.pi
                * r**3
                / volume
                * np.divide(spherical_jn(1, x), x, out=np.full_like(x, 1 / 3), where=x != 0)
            )
            surface = np.zeros_like(x)
            for ell, prefactor, j, slope in channels:
                surface += (
                    prefactor
                    * eval_legendre(ell, cosine)
                    * (slope[block, None] * j - j[block, None] * slope)
                )
            surface *= order
            for b, (h, overlap_kinetic) in enumerate(bars):
                # Removed sphere has +i(q_i-q_j); the Wronskian has the
                # opposite sign, with the conjugate of the native upper triangle.
                pair = sphere * overlap_kinetic - surface * h
                for site in sites:
                    site_phase = phase[site, block, None].conj() * phase[site, None, :]
                    w = -(pair * site_phase).imag
                    result[b, site] += w.sum(axis=1) @ q[block] - w.sum(axis=0) @ q
    return result


def atomic_geometry_gradients(
    state, directory, case, modes, fixed_band, screened_band, density_weight
):
    """Local Cartesian gradients of the fixed and induced geometry scalars.

    Matrix covectors already contain the DMFT window, density, spin and k
    weights. The raw valence-density derivative enters only the induced
    readout. Projection onto the symmetry orbits and MPI reduction follow
    outside this routine. All global fields are held in the central gauge.
    """
    s = state.structure
    result = np.zeros((2, len(s.pos), 3))
    if not modes:
        return result
    active = np.any([np.any(pattern != 0, axis=1) for _, _, pattern in modes], axis=0)
    geo = prepared_geometry_inputs(state, directory, case)
    lattice = NativeLattice.from_structure(s)
    nk = len(state.records)
    bands = state.bands.spherical
    spherical, nonspherical, projector = (defaultdict(list) for _ in range(3))
    occurrence = defaultdict(int)
    for a, ell, k, op in bands.spherical:
        key = (a, ell, k)
        site = geo.first[a] + occurrence[key]
        occurrence[key] += 1
        if active[site]:
            spherical[k].append((a, site, ell, op))
    for a, op, ks in bands.nonspherical:
        for i, (k, entry) in enumerate(zip(ks, op.entries)):
            site = geo.first[a] + i // nk
            if active[site]:
                nonspherical[k].append((a, site, op, entry))
    for a, ell, k, ic, cols, op in bands.projectors:
        orbs = geo.inl.cix[ic + 1]
        if len(orbs) != 1:
            raise NotImplementedError("single-site impurity required")
        site = orbs[0][0] - 1
        if active[site]:
            xyz, shift = geo.inl.atoms[site + 1][1:3]
            if shift is not None:
                raise NotImplementedError("shifted projectors not supplied")
            frame = (np.asarray(xyz) if len(xyz) else s.rotloc[a]) @ geo.axes[site]
            projector[k].append((a, site, ell, ic, cols, op, frame))

    spatial = state.direct.spatial
    density = spatial.template.from_vector(density_weight)
    mt = {}
    for a, atom in spatial.atoms.items():
        if np.any(active[geo.first[a] : geo.first[a] + s.mult[a]]):
            ym = np.array([density.rad[(a, *h)] for h in atom.kernel.layout.harmonics])
            mt[a], _ = atom.kernel.component_pullback(ym, include_radial=False)

    unpack = state.bands.output.point_covectors
    for k, (fixed, screened) in enumerate(zip(unpack(fixed_band), unpack(screened_band))):
        record, base = state.records[k], state.bases[k]
        c = getattr(base, "coefficients", None)
        if c is None:
            c = np.column_stack((base.retained, base.complement))
        npw = record.nbasis - record.nlo
        labels = record.labels.T
        q = geo.momenta[k]
        g = lattice.momenta(labels)
        phase = np.exp(2j * np.pi * (s.pos @ labels.T))
        weights = [{} for _ in range(s.nat)]

        def matching(
            a, site, ell, frame, ch, boundary, lo, lapw, bar, *, pw=None, projector_phase=False
        ):
            local_phase = (
                np.exp(2j * np.pi * (geo.fractional[k] @ geo.locations[site]))
                if projector_phase
                else phase[site]
            )
            groups = [npw + geo.lo[a, ell, j] for j in range(ch.nlo)]
            return matching_position_gradient(
                ell,
                ch.rmt,
                lattice.volume,
                q @ frame.T,
                geo.lengths[k],
                local_phase,
                q if projector_phase else g,
                c[:, : bar.shape[-1]],
                npw,
                groups,
                boundary,
                lo,
                lapw,
                bar,
                pw_only=pw,
            )

        for a, site, ell, op in spherical[k]:
            boundary = op.phi[-1]
            z, zpw = op._matching()
            kz, sz = op._products(z)
            value, slope = op._boundary_value(zpw), op._boundary_slope(z)
            if op.lapw:
                weights[a][ell] = 0.5 * op.r**2 * op.w * (op.w - 1)
            for b, (h, o, _) in enumerate((fixed, screened)):
                bar = 2 * (right_multiply_columns(kz, h) + right_multiply_columns(sz, o))
                # The boundary slope has all LO channels; the boundary value
                # uses only the plane-wave leg, including for an LAPW+LO.
                bar += (
                    op.r * boundary[2][None, :, None] * right_multiply_columns(value, h)[:, None, :]
                )
                pw = (
                    op.r * boundary[0][None, :, None] * right_multiply_columns(slope, h)[:, None, :]
                )
                result[b, site] += matching(
                    a, site, ell, np.eye(3), op.channel, boundary, op.lo, op.lapw, bar, pw=pw
                )
        for a, site, op, entry in nonspherical[k]:
            z, channels, lo, lapw, kz = entry
            frame = s.rotloc[a] @ geo.axes[site]
            for b, (h, _, _) in enumerate((fixed, screened)):
                bar = 2 * right_multiply_columns(kz, h)
                offset = 0
                for ell, ch in sorted(channels.items()):
                    size = (2 * ell + 1) * (2 + ch.nlo)
                    channel_bar = bar[offset : offset + size].reshape(2 * ell + 1, 2 + ch.nlo, -1)
                    result[b, site] += matching(
                        a, site, ell, frame, ch, op.phi[ell][-1], lo[ell], lapw[ell], channel_bar
                    )
                    offset += size
        for a, site, ell, ic, cols, op, frame in projector[k]:
            for b, (_, _, raw) in enumerate((fixed, screened)):
                y = raw[:, cols, ic]
                if np.any(y):
                    bar = np.einsum("na,ma,i->min", y.conj(), op.transform, op.q)
                    result[b, site] += matching(
                        a,
                        site,
                        ell,
                        frame,
                        op.channel,
                        op.phi[-1],
                        op.lo,
                        op.lapw,
                        bar,
                        projector_phase=True,
                    )
        for a, density_bar in mt.items():
            atom = spatial.atoms[a]
            for i in range(s.mult[a]):
                site = geo.first[a] + i
                if not active[site]:
                    continue
                entry = i * nk + k
                frame = s.rotloc[a] @ geo.axes[site]
                _, _, w = atom.entries[entry]
                bar = 2 * w * (density_bar @ atom.zrho[entry])
                for ell, ch in state.direct.matching[a][entry].items():
                    ctx = state.direct.contexts[a][ell]
                    channel_bar = bar[atom.kernel.layout.slices[ell]].reshape(
                        2 * ell + 1, 2 + ch.nlo, -1
                    )
                    result[1, site] += matching(
                        a,
                        site,
                        ell,
                        frame,
                        ch,
                        ctx.primitive[-1],
                        ctx.lo_coefficients,
                        ctx.lapw,
                        channel_bar,
                    )
        result += plane_wave_position_gradients(
            q[:npw], labels[:npw], c[:npw], s, lattice.volume, weights, (fixed, screened), active
        )
    return result


def atomic_geometry_readouts(
    state, directory, case, modes, fixed_band, screened_band, density_weight
):
    """Return two rank-local scalars per symmetry-allowed atomic coordinate.

    Each directional source builds a single band/density tangent. The reverse
    kernel instead projects four band matrices (fixed and induced H and O)
    and reuses them across coordinates. For a few independent orbits with
    only one active direction on each, there is little to amortize: the
    directional calculation is faster in the SrVO3 control. Both paths use
    the same analytic derivatives and k weights. Reduce only the final rows.
    """
    result = np.zeros((2, len(modes)))
    if not modes:
        return result
    patterns = np.array([pattern for _, _, pattern in modes])
    distinct_orbits = len({atom for atom, _, _ in modes})
    if len(modes) <= 4 and distinct_orbits == len(modes):
        fixed_cell = np.zeros((3, 3))
        for i, pattern in enumerate(patterns):
            band, _ = native_band_geometry_source(
                state, directory, case, fixed_cell, displacements=pattern
            )
            spatial = native_spatial_geometry_source(
                state, directory, case, fixed_cell, displacements=pattern
            )
            result[:, i] = fixed_band @ band, screened_band @ band + density_weight @ spatial
        return result
    fixed, induced = atomic_geometry_gradients(
        state, directory, case, modes, fixed_band, screened_band, density_weight
    )
    return np.stack([np.einsum("sij,ij->s", patterns, gradient) for gradient in (fixed, induced)])
