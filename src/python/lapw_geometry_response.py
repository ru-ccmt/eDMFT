"""Explicit strain derivatives of native LAPW H/O and raw projector overlaps.

Potential FILE coefficients, radial energies, RMT and local axes are held fixed.
This is the geometric source to add to B*(s_V+K*q), not a second potential
response. Matrices stay in the CENTRAL coefficient/band gauge, retaining both
off-diagonal legs and the native upper-triangle Wronskian convention.
The adapter supports prepared WIEN cell geometries and symmetry-preserving
strains, with separate electronic branches for DMFT-only magnetism. It does not publish a stress or construct the
full direct charge/core/normalization geometry source.
"""

from copy import copy
from functools import lru_cache
from pathlib import Path
import hashlib
from lapw_lattice import NativeLattice, native_site_frames
import numpy as np
from scipy.special import eval_legendre, spherical_jn
from lapw_augmentation_response import harmonics
from lapw_compact_blocks import hermitian_columns
from lapw_density_response import checked
from lapw_native_snapshot import read_basis_matrix
from lapw_strain_potential import read_potential_file
from indmffile import Indmfl


def harmonic_direction(ell, q, dq):
    """Y_lm and its analytic Cartesian directional derivative, including poles.

    For n=q/|q|, dn=(dq-n*n.dq)/|q| and omega=n x dn. The infinitesimal
    rotation gives dY=i*omega.L*Y, evaluated with angular-momentum ladders.
    This avoids dividing by sin(theta) near an axis. q=0 stays at q=0 under
    homogeneous strain and is assigned zero directional derivative.
    """
    q = np.asarray(q, float)
    dq = checked(dq, q.shape, "reciprocal derivative", real=True)
    y = harmonics(ell, q)
    length = np.linalg.norm(q, axis=1)
    n = np.divide(q, length[:, None], out=np.zeros_like(q), where=length[:, None] != 0)
    dn = np.divide(
        dq - n * np.sum(n * dq, axis=1)[:, None],
        length[:, None],
        out=np.zeros_like(q),
        where=length[:, None] != 0,
    )
    omega = np.cross(n, dn)
    dy = np.zeros_like(y)
    for i, m in enumerate(range(-ell, ell + 1)):
        plus = np.sqrt(ell * (ell + 1) - m * (m + 1)) * y[i + 1] if m < ell else 0
        minus = np.sqrt(ell * (ell + 1) - m * (m - 1)) * y[i - 1] if m > -ell else 0
        dy[i] = (
            0.5j * omega[:, 0] * (plus + minus)
            + 0.5 * omega[:, 1] * (plus - minus)
            + 1j * omega[:, 2] * m * y[i]
        )
    return y, dy


def matching_direction(ell, rmt, length, length_dot):
    """j_l(qR), q*j_l'(qR), and their derivative at fixed RMT.

    The spherical-Bessel ODE supplies j'', with analytic q=0 limits. Avoid
    numerical differences in j_l and avoid ratios to j_l at its zeros.
    """
    q = np.asarray(length, float)
    dq = checked(length_dot, q.shape, "length derivative", real=True)
    x = q * rmt
    j = spherical_jn(ell, x)
    jp = spherical_jn(ell, x, derivative=True)
    jpp = np.zeros_like(x)
    nz = x != 0
    jpp[nz] = -2 * jp[nz] / x[nz] + (ell * (ell + 1) / x[nz] ** 2 - 1) * j[nz]
    jpp[~nz] = {0: -1 / 3, 2: 2 / 15}.get(ell, 0.0)
    return np.array([j, q * jp]), np.array([rmt * dq * jp, dq * (jp + x * jpp)])


def augmentation_geometry_direction(
    channel,
    q_local,
    dq_local,
    length,
    length_dot,
    phase,
    coefficients,
    npw,
    lo_groups,
    volume,
    trace,
    phase_rate=None,
):
    """Return a channel whose linear matching arrays are their geometry tangent.

    g=4*pi*i^l*Y_lm* exp(iG.R)/sqrt(Omega). Fractional phases are fixed, so
    dg uses dY and -tr(strain)/2. LO factors carry the same normalization.
    Radial primitives and ABC coefficients are held fixed here.
    """
    c = np.asarray(coefficients)
    phase = np.asarray(phase)
    if phase_rate is None:
        y, dy = harmonic_direction(channel.ell, q_local, dq_local)
    else:
        # A translation changes only the Bloch phase at fixed cell and axes.
        # The matching arrays include LO rows, which carry their native G labels.
        y = harmonics(channel.ell, q_local)
        dy = np.zeros_like(y)
    pref = 4 * np.pi * (1j**channel.ell) / np.sqrt(volume) * phase
    g = pref * y.conj()
    dg = (1j * phase_rate * g if phase_rate is not None
          else pref * (dy.conj() - 0.5 * trace * y.conj()))
    match, dmatch = matching_direction(
        channel.ell, channel.rmt, length[:npw], length_dot[:npw]
    )
    out = copy(channel)
    out.fj = (dg[:, :npw] * match[0] + g[:, :npw] * dmatch[0]) @ c[:npw]
    out.dfj = (dg[:, :npw] * match[1] + g[:, :npw] * dmatch[1]) @ c[:npw]
    out.local = np.array(
        [dg[:, indices] @ c[indices] for indices in lo_groups], complex
    ).reshape(channel.local.shape)
    return out, (g[:, :npw], dg[:, :npw], match, dmatch)


def wronskian_geometry_direction(data, coefficients, retained):
    """Differentiate BEFORE the ordered native upper-triangle copy."""
    g, dg, m, dm = data
    angular = g.conj().T @ g
    da = dg.conj().T @ g + g.conj().T @ dg
    j, d = m
    jp, dp = dm
    w = d[:, None] * j[None, :] - j[:, None] * d[None, :]
    dw = (
        dp[:, None] * j[None, :]
        + d[:, None] * jp[None, :]
        - jp[:, None] * d[None, :]
        - j[:, None] * dp[None, :]
    )
    a = da * w + angular * dw
    upper = np.triu(a, 1)
    native = upper + upper.conj().T
    np.fill_diagonal(native, np.diag(a).real)
    c = coefficients
    return hermitian_columns(c.conj().T @ (native @ c[:, :retained]))


def atomic_wronskian_columns(q, labels, c, structure, velocities, volume, weights, retained):
    """Derivative of the native ordered surface kernel for an internal coordinate.

    At fixed cell, the radial factors and P_l(q_i.q_j) stay fixed. For each
    site only exp[i(G_j-G_i).R] changes. Sum the phase tangents over the orbit
    and the radial/angular kernels over l BEFORE projecting to band columns.
    The lower triangle comes from the conjugate of the native upper triangle:
    j'_i*j_j - j_i*j'_j is antisymmetric, hence the sign(j-i) below.
    """
    n = len(q)
    result = np.zeros((c.shape[1], retained), complex)
    length = np.linalg.norm(q, axis=1)
    unit = np.divide(q, length[:, None], out=np.zeros_like(q), where=length[:, None] > 0)
    atoms = np.repeat(np.arange(structure.nat), structure.mult)
    orbits = []
    for a, r in enumerate(structure.rmt):
        moving = np.flatnonzero((atoms == a) & np.any(velocities != 0, axis=1))
        if not len(moving):
            continue
        channels = []
        x = length * r
        for ell, weight in weights[a].items():
            if weight:
                channels.append((ell, (4*np.pi / volume) * (2*ell+1) * weight,
                                 spherical_jn(ell, x), length * spherical_jn(ell, x, derivative=True)))
        if channels:
            phases = np.exp(2j*np.pi*(structure.pos[moving] @ labels.T))
            rates = velocities[moving] @ q.T
            orbits.append((channels, phases, rates))
    if not orbits:
        return result
    rows = max(1, min(n, (8 * 1024**2) // (128 * n)))
    indices = np.arange(n)
    for start in range(0, n, rows):
        block = slice(start, min(n, start+rows))
        cosine = np.clip(unit[block] @ unit.T, -1, 1)
        tangent = np.zeros(cosine.shape, complex)
        for channels, phases, rates in orbits:
            dr_phase = np.zeros_like(tangent)
            for phase, rate in zip(phases, rates):
                dr_phase += (1j * phase[block, None].conj() * phase[None, :]
                             * (rate[None, :] - rate[block, None]))
            radial = np.zeros_like(cosine)
            for ell, prefactor, j, dj in channels:
                radial += (prefactor * eval_legendre(ell, cosine)
                           * (dj[block, None]*j[None, :] - j[block, None]*dj[None, :]))
            tangent += radial * dr_phase
        tangent *= np.sign(indices[None, :] - indices[block, None])
        result += c[block].conj().T @ (tangent @ c[:, :retained])
    return hermitian_columns(result)


def spherical_geometry_direction(op, channel_dot, triangle_dot):
    """Both matching legs of the small spherical kernels and surface term."""
    dz = channel_dot.value(op.phi[-1], op.lo, lapw=op.lapw)
    dzpw = (dz if channel_dot.nlo == 0 else
            channel_dot.value(op.phi[-1], np.zeros_like(op.lo), lapw=op.lapw))
    z, zpw = op._matching()
    kz, sz = op._products(z)
    p, d = op._boundary_value(zpw), op._boundary_slope(z)
    n = op.retained
    h = np.einsum("min,mip->np", dz.conj(), kz[..., :n], optimize=True) + np.einsum(
        "min,mip->np", kz.conj(), dz[..., :n], optimize=True
    )
    o = np.einsum("min,mip->np", dz.conj(), sz[..., :n], optimize=True) + np.einsum(
        "min,mip->np", sz.conj(), dz[..., :n], optimize=True
    )
    dp = np.einsum("i,min->mn", op.phi[-1, 0] / op.r, dzpw)
    dd = np.einsum("i,min->mn", op.phi[-1, 2] / op.r, dz)
    h += (
        op.r**2
        / 2
        * (
            dp.conj().T @ d[:, :n]
            + p.conj().T @ dd[:, :n]
            + d.conj().T @ dp[:, :n]
            + dd.conj().T @ p[:, :n]
        )
    )
    if op.lapw:
        h += 0.5 * op.r**2 * op.w * (op.w - 1) * triangle_dot
    return hermitian_columns(h), hermitian_columns(o)


def atomic_interstitial_columns(q, labels, c, positions, radii, volume, velocities, retained):
    """Project the moving-sphere overlap and kinetic matrices a few PW rows at a time.

    At fixed cell the kinetic metric stays fixed. The sphere transform in O_ij
    changes by -i(G_i-G_j).u times its central value. Its contribution to the
    interstitial overlap has a second minus sign because the sphere is excluded.
    Blocking bounds the pairwise workspace even for a large APW basis.
    """
    n = len(q)
    h = np.zeros((c.shape[1], retained), complex)
    o = np.zeros_like(h)
    # A row needs several real/complex arrays of length n, independent of nk.
    rows = max(1, min(n, (8 * 1024**2) // (128 * n)))
    moving = [(4*np.pi*r**3 / volume, r, np.exp(2j*np.pi*(labels @ pos)), q @ u)
              for pos, r, u in zip(positions, radii, velocities) if np.any(u)]
    for start in range(0, n, rows):
        block = slice(start, min(n, start + rows))
        length = np.linalg.norm(q[block, None, :] - q[None, :, :], axis=-1)
        overlap = np.zeros(length.shape, complex)
        for scale, r, phase, rate in moving:
            x = length * r
            shape = np.divide(spherical_jn(1, x), x, out=np.full_like(x, 1/3), where=x != 0)
            overlap += (1j * scale * shape * phase[block, None].conj() * phase[None, :]
                        * (rate[block, None] - rate[None, :]))
        bra = c[block].conj().T
        ket = c[:, :retained]
        o += bra @ (overlap @ ket)
        overlap *= q[block] @ q.T
        h += bra @ (overlap @ ket)
    return hermitian_columns(h), hermitian_columns(o)


def interstitial_geometry(
    q, labels, coefficients, positions, radii, volume, strain, retained, *, central=True, displacements=None
):
    """Exact analytic PW overlap and kinetic derivative, no potential derivative.

    Theta is the untruncated sphere-exclusion Fourier matrix, unlike the
    LAPW0 finite-FFT potential mask. T_ij=(q_i.q_j)*Theta_ij in Ry. The native
    masked potential coefficients are fixed, so their matrix derivative is
    zero in this coefficient gauge. LO rows have no interstitial contribution.
    """
    q = np.asarray(q)
    labels = np.asarray(labels)
    c = np.asarray(coefficients)
    if displacements is not None and not central:
        dh, do = atomic_interstitial_columns(q, labels, c, positions, radii, volume, displacements, retained)
        return dh, do, None, None
    diff = q[:, None] - q[None, :]
    norm = np.linalg.norm(diff, axis=-1)
    direction = np.divide(
        np.einsum("...i,ij,...j->...", diff, strain, diff),
        norm**2,
        out=np.zeros_like(norm),
        where=norm != 0,
    )
    theta = np.eye(len(q), dtype=complex)
    dtheta = np.zeros_like(theta)
    for site, (pos, r) in enumerate(zip(positions, radii)):
        x = norm * r
        f = np.divide(spherical_jn(1, x), x, out=np.full_like(x, 1 / 3), where=x != 0)
        phase = np.exp(2j * np.pi * (labels @ pos))
        phase = phase.conj()[:, None] * phase[None, :]
        scale = 4 * np.pi * r**3 / volume
        theta -= scale * f * phase
        if displacements is not None:
            # diff=G_i-G_j. The removed sphere carries exp[-i diff.R].
            dtheta += 1j * scale * f * phase * (diff @ displacements[site])
            continue
        dtheta += (
            scale
            * (np.trace(strain) * f + direction * (spherical_jn(0, x) - 3 * f))
            * phase
        )
    kinetic = q @ q.T
    dt = -2 * (q @ strain @ q.T) * theta + kinetic * dtheta
    project = lambda a: hermitian_columns(c.conj().T @ (a @ c[:, :retained]))
    return (
        project(dt),
        project(dtheta),
        project(kinetic * theta) if central else None,
        project(theta) if central else None,
    )


def prepared_geometry_inputs(state, directory, case):
    """Small immutable central geometry shared by all strain/density sources.

    H/O exports were already read and fingerprinted during valence preparation.
    Keep their three-component k coordinates instead of rereading dense N^2
    files for every strain. Matrices stay distributed with their owning rank.
    This cache belongs to one prepared state, never to a mutable material run.
    """
    if hasattr(state, "_strain_geometry"):
        return state._strain_geometry
    from types import SimpleNamespace
    from response_parallel import parallel_of

    parallel = parallel_of(state.model)
    directory = Path(directory)
    s = state.structure

    channel = getattr(state, "spin_channel", "")
    projector_name = case + ".indmfl" + channel
    has_projector = projector_name in state.fingerprints
    def validate():
        for name in ([projector_name] if has_projector else []) + [case + ".vns", case + ".struct"]:
            if (
                hashlib.sha256((directory / name).read_bytes()).hexdigest()
                != state.fingerprints[name]
            ):
                raise ValueError("changed geometry-source input: " + name)

    parallel.root_call(validate)
    inl = None
    if has_projector:
        inl = Indmfl(str(directory / case), extn="indmfl" + channel)
        inl.read()
    first = np.r_[0, np.cumsum(s.mult)[:-1]]
    axes, locations = native_site_frames(s)
    lo, offset = {}, 0
    for a in range(s.nat):
        for l in range(4):
            for j in range(state.projector_contexts[a][l].nlo):
                n = s.mult[a] * (2 * l + 1)
                lo[a, l, j] = np.arange(offset, offset + n)
                offset += n
    if any(r.nlo != offset for r in state.records):
        raise ValueError("native LO rows differ")
    kcoords = [r.reciprocal_k for r in state.records]
    fractional = [r.labels.T + k for r, k in zip(state.records, kcoords)]
    momenta = [NativeLattice.from_structure(s).momenta(g) for g in fractional]
    state._strain_geometry = SimpleNamespace(
        inl=inl,
        first=first,
        axes=axes,
        locations=locations,
        lo=lo,
        fractional=fractional,
        momenta=momenta,
        lengths=[np.linalg.norm(q, axis=1) for q in momenta],
        vns=read_potential_file(directory / (case + ".vns"), s.jrj),
        central_report=None,
    )
    return state._strain_geometry


def native_band_geometry_source(state, directory, case, strain, *, displacements=None):
    """Assemble explicit joint H/O/raw geometry in central RR/QR coordinates.

    Reuse prepared radial kernels and current native input conventions.
    Symmetry-breaking strains require expanding equivalent-site/orbital and
    potential layouts; this adapter rejects them instead of silently applying
    a ground-state-only symmetry representation. Keep all state contexts open.
    No self-energy, normalization X or window derivative is inserted here:
    the shared generalized electronic/projector map applies those afterwards.
    """
    if hasattr(state, "branches"):
        pieces = [native_band_geometry_source(v, directory, case, strain, displacements=displacements) for v in state.branches]
        return state.bands.output.pack(*(x[0] for x in pieces)), dict(spin_branches=[x[1] for x in pieces])
    s = state.structure
    bands = state.bands.spherical
    if displacements is not None:
        displacements = checked(displacements, s.pos.shape, "site velocities (Bohr/Bohr)", real=True)
        if np.any(strain):
            raise ValueError("atomic translation requires a fixed cell")
    e = checked(strain, (3, 3), "strain direction", real=True)
    if not np.allclose(e, e.T, atol=1e-14, rtol=0):
        raise ValueError("symmetric strain required")
    for rot in NativeLattice.from_structure(s).cartesian_rotations(s.timat):
        if np.linalg.norm(e @ rot - rot @ e) > 1e-10:
            raise NotImplementedError(
                "expand symmetry/equivalent-site layouts before this strain"
            )
    from response_parallel import parallel_of

    parallel = parallel_of(state.model)
    prepared = prepared_geometry_inputs(state, directory, case)
    inl = prepared.inl
    first, axes, locations, lo_indices = (
        prepared.first,
        prepared.axes,
        prepared.locations,
        prepared.lo,
    )
    nr = state.model.eps.shape[1]
    nk = len(state.records)
    trace = np.trace(e)
    lattice = NativeLattice.from_structure(s)
    volume = lattice.volume
    check_central = displacements is None and prepared.central_report is None

    @lru_cache(maxsize=1 if displacements is not None else None)
    def geometry(k):
        rec = state.records[k]
        base = state.bases[k]
        c = getattr(base, "coefficients", None)
        if c is None:
            c = np.column_stack((base.retained, base.complement))
        npw = rec.nbasis - rec.nlo
        frac = prepared.fractional[k]
        q = prepared.momenta[k]
        dq = -q @ e
        length = prepared.lengths[k]
        dl = np.divide(
            np.sum(q * dq, axis=1), length, out=np.zeros_like(length), where=length != 0
        )
        return c, npw, frac, q, dq, length, dl

    def channel_dot(ch, a, l, k, site, frame, projector=False):
        c, npw, frac, q, dq, length, dl = geometry(k)
        phase = np.exp(
            2j
            * np.pi
            * (
                frac @ locations[site]
                if projector
                else state.records[k].labels.T @ s.pos[site]
            )
        )
        phase_rate = None
        if displacements is not None:
            # The native matrix gauge uses G.R; the projector uses (k+G).R.
            # Equivalent-site positions may differ by a constant lattice translation.
            wavevectors = q if projector else lattice.momenta(state.records[k].labels.T)
            phase_rate = wavevectors @ displacements[site]
        groups = [npw + lo_indices[a, l, j] for j in range(ch.nlo)]
        return augmentation_geometry_direction(
            ch,
            q @ frame.T,
            dq @ frame.T,
            length,
            dl,
            phase,
            c,
            npw,
            groups,
            volume,
            trace, phase_rate=phase_rate,
        )

    h = [np.zeros((b.retained.shape[0], nr), complex) for b in state.bases]
    o = [np.zeros_like(a) for a in h]
    hc = [] if displacements is not None else [a.copy() for a in h]
    oc = [] if displacements is not None else [a.copy() for a in h]
    raw = np.zeros(state.model.raw.shape, complex)
    counts = {}
    wronskian_weights = [{} for _ in range(s.nat)]
    for a, l, k, op in bands.spherical:
        key = (a, l, k)
        site = first[a] + counts.get(key, 0)
        counts[key] = counts.get(key, 0) + 1
        if displacements is not None:
            if not np.any(displacements[site]):
                continue
            if op.lapw and l not in wronskian_weights[a]:
                wronskian_weights[a][l] = 0.5 * op.r**2 * op.w * (op.w - 1)
        dot, data = channel_dot(op.channel, a, l, k, site, np.eye(3))
        c, npw, *_ = geometry(k)
        tri = wronskian_geometry_direction(data, c[:npw], nr) if op.lapw and displacements is None else 0
        dh, do = spherical_geometry_direction(op, dot, tri)
        h[k] += dh
        o[k] += do
        if check_central:
            hv, ov = op.value(orbit_sum=True)
            hc[k] += hv
            oc[k] += ov
    for a, op, ks in bands.nonspherical:
        for index, (k, entry) in enumerate(zip(ks, op.entries)):
            z, channels, lo, lapw, kz = entry
            site = first[a] + index // nk
            if displacements is not None and not np.any(displacements[site]):
                continue
            dz = []
            for l, ch in sorted(channels.items()):
                dot, _ = channel_dot(ch, a, l, k, site, s.rotloc[a] @ axes[site])
                dz.append(
                    dot.value(op.phi[l][-1], lo[l], lapw=lapw[l]).reshape(
                        -1, z.shape[1]
                    )
                )
            dz = np.concatenate(dz)
            h[k] += hermitian_columns(
                dz.conj().T @ kz[:, :nr] + kz.conj().T @ dz[:, :nr]
            )
            if check_central:
                hc[k] += hermitian_columns(z.conj().T @ kz[:, :nr])
    for a, l, k, ic, cols, op in bands.projectors:
        orbs = inl.cix[ic + 1]
        if len(orbs) != 1:
            raise NotImplementedError("single-site impurity required")
        site = orbs[0][0] - 1
        if displacements is not None and not np.any(displacements[site]):
            continue
        _, xyz, shift, _ = inl.atoms[site + 1]
        if shift is not None:
            raise NotImplementedError("shifted projectors not supplied")
        xyz = np.array(xyz) if len(xyz) else s.rotloc[a]
        dot, _ = channel_dot(
            op.channel, a, l, k, site, xyz @ axes[site], projector=True
        )
        dz = dot.value(op.phi[-1], op.lo, lapw=op.lapw)
        raw[k, :, cols, ic] += np.einsum("min,i,ma->na", dz.conj(), op.q, op.transform)[
            :nr
        ]
    vns = prepared.vns
    for k in range(nk):
        c, npw, frac, q, dq, length, dl = geometry(k)
        dh, do, t0, o0 = interstitial_geometry(
            q[:npw],
            state.records[k].labels[:, :npw].T,
            c[:npw],
            s.pos,
            np.repeat(s.rmt, s.mult),
            volume,
            e,
            nr,
            central=check_central, displacements=displacements,
        )
        if displacements is not None:
            dh += atomic_wronskian_columns(
                q[:npw], state.records[k].labels[:, :npw].T, c[:npw], s,
                displacements, volume, wronskian_weights, nr,
            )
        h[k] += dh
        o[k] += do
        if check_central:
            hc[k] += t0
            oc[k] += o0
            ir = state.bands.ir[k]
            pw = np.array([vns.interstitial[tuple(g)] for g in ir.representatives])
            if not ir.complex_values:
                pw = pw.real
            rr, qr = ir.action(pw, retained=nr)
            hc[k] += np.concatenate((rr, qr))
    joint = state.bands.output.pack(
        np.array([a[:nr] for a in h]), np.array([a[:nr] for a in o]), raw,
        [a[nr:] for a in h], [a[nr:] for a in o],
    )
    if displacements is not None:
        return joint, dict(scope="atomic translation at fixed native potential coefficients")
    expected_h = [
        np.concatenate((np.diag(ev), np.zeros((len(hh) - nr, nr))))
        for ev, hh in zip(state.model.eps, h)
    ]
    expected_o = [
        np.concatenate((np.eye(nr), np.zeros((len(oo) - nr, nr)))) for oo in o
    ]
    report = dict(
        strain=e.tolist(),
        scope="explicit H/O/raw geometry at fixed native potential coefficients",
        central_H_relative=float(
            np.sqrt(
                parallel.sum(
                    sum(np.linalg.norm(a - b) ** 2 for a, b in zip(hc, expected_h))
                )
                / parallel.sum(sum(np.linalg.norm(b) ** 2 for b in expected_h))
            )
        ),
        central_O_relative=float(
            np.sqrt(
                parallel.sum(
                    sum(np.linalg.norm(a - b) ** 2 for a, b in zip(oc, expected_o))
                )
                / parallel.sum(nk * nr)
            )
        ),
        includes_window_normalization=False,
        stress_applied=False,
    )
    if check_central:
        prepared.central_report = {
            k: report[k] for k in ("central_H_relative", "central_O_relative")
        }
    else:
        report.update(prepared.central_report)
    return joint, report
