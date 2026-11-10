"""General symmetry-preserving strain of the native potential map.

All six symmetric Cartesian components are supported, including shear. The
native star/equivalent-atom layout can be reused precisely when e commutes
with its Cartesian point group. Symmetry-forbidden first stress derivatives
are zero; symmetry-breaking response/elastic constants need expanded layouts.
The density-coordinate transport is the independent scalar convention tr(e)/3.
"""

from lapw_lattice import NativeLattice
import numpy as np
from scipy.special import spherical_jn, factorial2
from lapw_density_response import checked
from lapw_geometry_response import harmonic_direction
from lapw_hartree_periodic import native_local_axes
from lapw_strain_potential import strain_mask_derivative, masked_potential_strain_source
from wien_response_operators import hydrostatic_transport_source


def harmonic_average_direction(stars, reciprocal, position, local_axes, lm, strain):
    """Differentiate the native YKAV star sum, retaining phases and star weights."""
    q = stars.members @ reciprocal.T
    local = q @ local_axes.T
    direction = (-q @ strain) @ local_axes.T
    phase = (
        stars.phases
        * np.exp(2j * np.pi * (stars.members @ position))
        / stars.sizes[stars.owners]
    )
    cache = {l: harmonic_direction(l, local, direction)[1] for l in {l for l, m in lm}}
    out = np.empty((len(lm), stars.nstar), complex)
    for i, (l, m) in enumerate(lm):
        values = phase * cache[l][l + m].conj()
        out[i] = np.bincount(
            stars.owners, values.real, minlength=stars.nstar
        ) + 1j * np.bincount(stars.owners, values.imag, minlength=stars.nstar)
        out[i, 0] = 0
    return out


def diagonal_strain_potential(prepared, strain, *, transport_rate=None):
    """Analytic potential source for symmetric cell strain plus scalar transport.

    q'=-e q, Omega'=tr(e)Omega. Let gamma_G=(q.e.q)/q^2. Differentiate the
    Coulomb pseudocharge matrices, including YKAV angular motion, without
    ratios to spherical harmonics or Bessel functions. XC gradients transform
    as -e grad and Hessians as -e Hessian-Hessian e. RMT/radial grids stay fixed.
    A caller that supplies general symmetric e must deliberately choose its
    density-coordinate convention; default trace(e)/3 matches native probes.
    """
    e = checked(strain, (3, 3), "symmetric strain", real=True)
    if not np.allclose(e, e.T, atol=1e-14, rtol=0):
        raise ValueError("symmetric cell strain required")
    h, s = prepared.hartree, prepared.structure
    if any(np.linalg.norm(e @ rot - rot @ e) > 1e-10 for rot in NativeLattice.from_structure(s).cartesian_rotations(s.timat)):
        raise NotImplementedError(
            "expand native symmetry/star layout before this strain"
        )
    trace = np.trace(e)
    rate = trace / 3 if transport_rate is None else float(transport_rate)
    if not np.isfinite(rate):
        raise ValueError("finite transport rate required")
    b = NativeLattice.from_structure(s).reciprocal
    q = h.stars.representatives @ b.T
    norm = np.linalg.norm(q, axis=1)
    gamma = np.divide(
        np.einsum("ni,ij,nj->n", q, e, q),
        norm**2,
        out=np.zeros_like(norm),
        where=norm != 0,
    )
    dA = []
    dD = []
    dB = np.zeros_like(h.B)
    first = np.r_[0, np.cumsum(s.mult)[:-1]]
    for atom in range(s.nat):
        keys = [k for k, key in enumerate(h.keys) if key[0] == atom]
        lm = sorted(
            set(
                (abs(h.keys[k][1]), m)
                for k in keys
                for m in {h.keys[k][2], -h.keys[k][2]}
            )
        )
        lookup = {key: i for i, key in enumerate(lm)}
        args = (b, s.pos[first[atom]], native_local_axes(s.rotloc[atom]), lm)
        y = h.stars.harmonic_average(*args)
        dy = harmonic_average_direction(h.stars, *args, e)
        radius = s.rmt[atom]
        x = norm[1:] * radius
        g = gamma[1:]
        dd = np.zeros((len(lm), h.stars.nstar), complex)
        da = np.zeros((h.stars.nstar, len(lm)), complex)
        for i, (l, m) in enumerate(lm):
            j = spherical_jn(l + 1, x)
            jp = spherical_jn(l + 1, x, derivative=True)
            pref = 4 * np.pi * radius ** (l + 2) * 1j**l / norm[1:]
            dd[i, 1:] = pref * (g * (j - x * jp) * y[i, 1:] + j * dy[i, 1:])
            power = h.ncon - l + 1
            pref = (
                -h.diagonal[1:]
                * (4 * np.pi / h.volume)
                * h.stars.sizes[1:]
                * s.mult[atom]
                * (-1j / radius) ** l
                * factorial2(2 * h.ncon + 3)
                / factorial2(2 * l + 1)
            )
            j = spherical_jn(h.ncon + 1, x)
            jp = spherical_jn(h.ncon + 1, x, derivative=True)
            da[1:, i] = pref * (
                (((2 + power) * g - trace) * j / x**power - g * jp / x ** (power - 1))
                * y[i, 1:].conj()
                + j / x**power * dy[i, 1:].conj()
            )
        for k in keys:
            a, sl, m = h.keys[k]
            l = abs(sl)
            phase = (-1j if sl < 0 else 1) * (-1) ** m
            factor = 4 * np.pi * 1j**l
            if m:
                factor *= np.sqrt(2) / phase
            if l == 0:
                factor /= np.sqrt(4 * np.pi)
            dB[k, 1:] = h.boundary_factors[k] * factor * (
                -g * x * spherical_jn(l, x, derivative=True) * y[lookup[l, m], 1:]
                + spherical_jn(l, x) * dy[lookup[l, m], 1:]
            )
        dA.append(da)
        dD.append(dd)
    dA = np.concatenate(dA, axis=1)
    dD = np.concatenate(dD, axis=0)
    radial, rho = h._split(prepared.input_density.vector())
    moments = np.array([w @ c for w, c in zip(h.moments, radial)])
    z = np.array(
        [
            s.Znuc[a] / np.sqrt(4 * np.pi) if (l, m) == (0, 0) else 0
            for a, l, m in h.keys
        ]
    )
    multipoles = -h.S @ moments + h.D @ rho + h.S @ z
    p = h.diagonal * rho + h.A @ multipoles
    dp = 2 * gamma * h.diagonal * rho + dA @ multipoles + h.A @ (dD @ rho)
    db = (dB @ p + h.B @ dp).real
    mask = strain_mask_derivative(
        h.reanalysis, b, s.pos, np.repeat(s.rmt, s.mult), h.volume, e
    )
    mt = [op.boundary * bound for op, bound in zip(h.radial, db)]
    ir = h.reanalysis.action(dp) + masked_potential_strain_source(h.reanalysis, p, mask)
    coulomb = np.concatenate(mt + [ir.view(float)])
    xc = prepared.xc.ir
    density = np.maximum(xc.expand(rho), 0)
    if xc.functional == "PBE":
        fields = xc.fields(density)
        direction = np.zeros_like(fields)
        direction[..., 1:4] = -fields[..., 1:4] @ e
        pairs = [(0,0),(1,1),(2,2),(1,2),(0,2),(0,1)]
        hessian = np.empty(fields.shape[:-1] + (3,3))
        for i,(a,b) in enumerate(pairs):
            hessian[...,a,b] = hessian[...,b,a] = fields[...,4+i]
        derivative = -e @ hessian - hessian @ e
        for i,(a,b) in enumerate(pairs):
            direction[...,4+i] = derivative[...,a,b]
        point = np.imag(xc._point(fields + 1j * 1e-24 * direction)) / 1e-24
    else:
        point = np.zeros(xc.shape)
    dp = xc.select(xc.average(point))
    pw = h.reanalysis.action(dp) + masked_potential_strain_source(
        h.reanalysis, xc.central_stars, mask
    )
    vxc = np.zeros(h.size)
    vxc[h.pw_offset :] = pw.view(float)
    transport = rate * prepared.potential_operator.matvec(
        hydrostatic_transport_source(prepared.input_density).vector()
    )
    return transport + coulomb + vxc, dict(
        transport=transport,
        coulomb_geometry=coulomb,
        xc_geometry=vxc,
        mask_derivative=mask,
    )

# Historical function name retained for existing importers.
symmetric_strain_potential = diagonal_strain_potential
