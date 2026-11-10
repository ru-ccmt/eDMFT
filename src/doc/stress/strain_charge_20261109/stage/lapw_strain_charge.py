"""Explicit hydrostatic charge source in the native reference's coordinates.

This is one UNMIXED physical update at fixed input material charge and fixed
local Sigma/DC. No impurity vertex, charge iteration or numerical derivative
is used. The source includes native normalization before inverse transport
AND central-cell normalization after it, as in the checked native reference.
The returned source is ready for (I-J)q=s; it is not a stress correction.
"""

import hashlib
from pathlib import Path
import numpy as np
from scipy.special import spherical_jn
from lapw_density_response import checked
from lapw_geometry_response import (
    augmentation_geometry_direction,
    native_band_geometry_source,
)
from lapw_native_snapshot import read_basis_matrix
from lapw_charge_response import NativeChargeNormalization, native_mixer_inversion
from lapw_hydrostatic_potential import hydrostatic_explicit_potential
from wien_response_operators import hydrostatic_transport_source


def charge_integration_geometry(density, structure, strain, *, inversion):
    """Derivative of MIXER's IR integral at fixed native density coefficients.

    With fixed RMT, the MT integration covector does not move. The integrated
    IR mask is Omega*delta_G0-sum_a 4*pi*R_a^3*j1(GR)/(GR)*phase_a.
    Its derivative is tr(e)*Omega at G=0, plus the reciprocal-length term
    sum_a 4*pi*R_a^3*(G.e.G/G^2)*(j0-3*j1/x)*phase_a.
    Density PW entries are star sums; no extra star multiplicity is inserted.
    Returns the IR covector derivative in Density.vector real coordinates.
    """
    e = checked(strain, (3, 3), "strain", real=True)
    s = structure
    if not isinstance(inversion, (bool, np.bool_)):
        raise ValueError("explicit inversion flag required")
    if s.lattyp.strip() != "P" or not np.allclose(
        [s.alpha, s.beta, s.gamma], 90, atol=1e-10, rtol=0
    ):
        raise NotImplementedError("primitive orthogonal charge geometry required")
    if not np.allclose(e, e.T, atol=1e-14, rtol=0):
        raise ValueError("symmetric strain required")
    reps = np.array(list(density.pw))
    cell = np.array([s.a, s.b, s.c])
    q = 2 * np.pi * reps / cell
    length = np.linalg.norm(q, axis=1)
    direction = np.divide(
        np.einsum("ni,ij,nj->n", q, e, q),
        length**2,
        out=np.zeros_like(length),
        where=length != 0,
    )
    mask = np.zeros(len(q), complex)
    mask[np.all(reps == 0, axis=1)] = np.trace(e) * np.prod(cell)
    for pos, r in zip(s.pos, np.repeat(s.rmt, s.mult)):
        x = length * r
        f = np.divide(spherical_jn(1, x), x, out=np.full_like(x, 1 / 3), where=x != 0)
        mask += (
            4
            * np.pi
            * r**3
            * direction
            * (spherical_jn(0, x) - 3 * f)
            * np.exp(-2j * np.pi * (reps @ pos))
        )
    if inversion:
        mask = mask.real.astype(complex)
    out = np.zeros_like(density.vector())
    offset = sum(map(len, density.rad.values()))
    out[offset:] = mask.view(float)
    return out


def spatial_matching_source(spatial, matching_dots, volume_trace):
    """Both MT matching legs and the native PW 1/Omega derivative at fixed rho.

    Fractional coordinates, G labels, star/symmetry phases and rho are fixed.
    The IR Fourier series has no moving mask here: clmval is its unmasked
    continuation. The changing IR integration mask enters normalization later.
    """
    out = spatial.action(matching_dots=matching_dots)
    nrad = sum(map(len, spatial.template.rad.values()))
    out[nrad:] -= volume_trace * spatial.value()[nrad:]
    return out


def native_spatial_geometry_source(valence, directory, case, strain):
    """Explicit valence matching/volume source in the stored harmonic frames."""
    s = valence.structure
    e = checked(strain, (3, 3), "strain", real=True)
    if not np.allclose(e, e.T, atol=1e-14, rtol=0):
        raise ValueError("symmetric strain required")
    if any(np.linalg.norm(e @ rot - rot @ e) > 1e-10 for rot in s.timat):
        raise NotImplementedError("expand symmetry layout before this strain")
    directory = Path(directory)
    direct = valence.direct
    first = np.r_[0, np.cumsum(s.mult)[:-1]]
    axes = []
    for site, a in enumerate(np.repeat(np.arange(s.nat), s.mult)):
        for rot, tau in zip(s.timat, s.tau):
            delta = rot @ s.pos[first[a]] + tau - s.pos[site]
            delta -= np.rint(delta)
            if max(abs(delta)) < 2e-7:
                break
        else:
            raise ValueError("missing equivalent-site frame")
        axes.append(rot.T)
    lo = {}
    offset = 0
    for a in range(s.nat):
        for l in range(4):
            for j in range(valence.projector_contexts[a][l].nlo):
                count = s.mult[a] * (2 * l + 1)
                lo[a, l, j] = np.arange(offset, offset + count)
                offset += count
    nk = len(valence.records)
    geom = []
    for k, (record, base) in enumerate(zip(valence.records, valence.bases)):
        path = directory / f"native.k{k+1:08d}.bin"
        if (
            hashlib.sha256(path.read_bytes()).hexdigest()
            != valence.fingerprints[path.name]
        ):
            raise ValueError("changed native snapshot")
        snapshot = read_basis_matrix(path, b"EDMFTOV1")
        q = (
            2
            * np.pi
            * (record.labels.T + snapshot.reciprocal_k)
            / np.array([s.a, s.b, s.c])
        )
        dq = -q @ e
        length = np.linalg.norm(q, axis=1)
        dl = np.divide(
            np.sum(q * dq, axis=1), length, out=np.zeros_like(length), where=length != 0
        )
        geom.append((base.retained, record.nbasis - record.nlo, q, dq, length, dl))
    zd = {}
    for a, entries in direct.matching.items():
        zd[a] = []
        for index, channels in enumerate(entries):
            k = index % nk
            site = first[a] + index // nk
            c, npw, q, dq, length, dl = geom[k]
            frame = s.rotloc[a] @ axes[site]
            phase = np.exp(2j * np.pi * (valence.records[k].labels.T @ s.pos[site]))
            blocks = []
            for l, ch in sorted(channels.items()):
                ctx = direct.contexts[a][l]
                dot, _ = augmentation_geometry_direction(
                    ch,
                    q @ frame.T,
                    dq @ frame.T,
                    length,
                    dl,
                    phase,
                    c,
                    npw,
                    [npw + lo[a, l, j] for j in range(ch.nlo)],
                    s.a * s.b * s.c,
                    np.trace(e),
                )
                blocks.append(
                    dot.value(
                        ctx.primitive[-1], ctx.lo_coefficients, lapw=ctx.lapw
                    ).reshape(-1, ch.nband)
                )
            zd[a].append(np.concatenate(blocks))
    return spatial_matching_source(direct.spatial, zd, np.trace(e))


def material_source_from_raw(state, raw_source, strain, *, transport_rate=1.0):
    """Normalize in strained cell, pull back, normalize in central cell.

    T_eta input pushforward has derivative tau_n; inverse transport has
    -tau evaluated at the central NORMALIZED OUTPUT, not at input n0.
    The last normalization has alpha=1 at that output. This second operation
    is part of the native reference's material-coordinate map and is generally
    necessary for the explicit source even though it leaves ordinary J intact.
    """
    pc = state.potential_core
    n = state.normalization
    template = pc.input_density
    inversion = native_mixer_inversion(pc.structure)
    ir = charge_integration_geometry(
        template, pc.structure, strain, inversion=inversion
    )
    ng = n.explicit_source(ir_geometry=float(ir @ n.base))
    normalized = n.action(raw_source) + ng
    output = template.from_vector(n.central)
    pullback = (
        transport_rate * hydrostatic_transport_source(output, pullback=True).vector()
    )
    post = NativeChargeNormalization(
        output, pc.structure, n.target, inversion=inversion, enabled=n.enabled
    )
    source = post.action(normalized + pullback)
    return source, dict(
        normalized_native=normalized,
        normalization_geometry=ng,
        inverse_transport=pullback,
        post_normalization_change=source - normalized - pullback,
    )


def hydrostatic_charge_source(
    state, directory, case, *, strain=None, potential_source=None
):
    """Complete explicit hydrostatic source for the prepared native charge map.

    Optional strain/potential_source are an explicit diagnostic override for
    rounded near-hydrostatic native structures; density transport still uses
    trace(strain)/3. The default e=I uses the fully analytic installed s_V.
    A nonhydrostatic override requires its matching potential source and is
    not an implementation of general anisotropic density transport.
    """
    e = np.eye(3) if strain is None else checked(strain, (3, 3), "strain", real=True)
    rate = np.trace(e) / 3
    if potential_source is None:
        if not np.allclose(e, rate * np.eye(3), rtol=0, atol=1e-12):
            raise ValueError(
                "nonhydrostatic diagnostic requires matching potential source"
            )
        potential_source, _ = hydrostatic_explicit_potential(state.potential_core)
        potential_source = rate * potential_source
    potential_source = checked(
        potential_source,
        (state.chain.potential.output_size,),
        "explicit potential source",
        real=True,
    )
    band, report = native_band_geometry_source(state.valence, directory, case, e)
    geometry = native_spatial_geometry_source(state.valence, directory, case, e)
    ch = state.chain
    joint = band + ch.basis.matvec(potential_source)
    electronic = ch.electronic.matvec(joint)
    spatial = ch.density.matvec(electronic)
    radial = ch.basis_density.matvec(potential_source)
    core = ch.core.matvec(potential_source)
    raw = spatial + radial + core + geometry
    source, parts = material_source_from_raw(state, raw, e, transport_rate=rate)
    parts.update(
        band_geometry=band,
        potential_source=potential_source,
        joint=joint,
        electronic_density=electronic,
        valence_electronic=spatial,
        valence_radial=radial,
        core=core,
        valence_geometry=geometry,
        raw_output=raw,
    )
    report.update(
        scope="explicit source of native normalized material charge update",
        geometry_source_includes_window_normalization=False,
        includes_window_normalization=True,
        source_norm=float(np.linalg.norm(source)),
        source_charge=float((state.normalization.mt + state.normalization.ir) @ source),
        transport_rate=float(rate),
        includes_post_transport_normalization=True,
    )
    return source, parts, report
