"""Analytic hydrostatic LAPW0 source at fixed RMT and fractional ionic sites.

For A_eta=(1+eta)A and input Density.transport(1+eta), compute
    v_exp = K * tau_n + (dH_geometry)n + dV_ion + dV_xc,geometry.
The output is in PreparedPotentialCore field coordinates. There are no native
launches, geometry finite differences, nuclear finite differences or self-
consistent iterations. Frozen branch/cutoff assumptions match the prepared
potential/core state. Only hydrostatic primitive orthogonal cells are supplied.
"""

import numpy as np
from scipy.special import spherical_jn, factorial2
from lapw_hartree_periodic import native_local_axes
from wien_response_operators import hydrostatic_transport_source
from lapw_strain_potential import strain_mask_derivative, masked_potential_strain_source


def hydrostatic_coulomb_geometry(prepared):
    """Coulomb source at fixed NATIVE density coefficients, including nuclei.

    Native notation: Q=-S*m+D*rho+Q_ion; P=diag*rho+A*Q; boundary=Re(B*P).
    At fixed radial grid/RMT, S and m do not move; dQ=dD*rho. Reciprocal lengths
    decrease, volume increases, while Y_lm(local G) and fractional phases stay
    fixed. Differentiate spherical Bessel functions without division by them,
    so their zeros are regular. The point nucleus has Q_ion=S*Z/sqrt(4*pi);
    its local Dirichlet potential -2Z(1-r/R) in rV units is geometry invariant.
    """
    h, s = prepared.hartree, prepared.structure
    b = np.diag(2 * np.pi / np.array([s.a, s.b, s.c]))
    norm = np.linalg.norm(h.stars.representatives @ b.T, axis=1)
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
        y = h.stars.harmonic_average(
            b, s.pos[first[atom]], native_local_axes(s.rotloc[atom]), lm
        )
        radius = s.rmt[atom]
        x = norm[1:] * radius
        dd = np.zeros((len(lm), h.stars.nstar), complex)
        da = np.zeros((h.stars.nstar, len(lm)), complex)
        for i, (l, m) in enumerate(lm):
            dd[i, 1:] = (
                4
                * np.pi
                * radius ** (l + 2)
                * 1j**l
                * y[i, 1:]
                / norm[1:]
                * (spherical_jn(l + 1, x) - x * spherical_jn(l + 1, x, derivative=True))
            )
            power = h.ncon - l + 1
            pref = (
                -h.diagonal[1:]
                * (4 * np.pi / h.volume)
                * h.stars.sizes[1:]
                * s.mult[atom]
                * (-1j / radius) ** l
                * factorial2(2 * h.ncon + 3)
                / factorial2(2 * l + 1)
                * y[i, 1:].conj()
            )
            da[1:, i] = pref * (
                (power - 1) * spherical_jn(h.ncon + 1, x) / x**power
                - spherical_jn(h.ncon + 1, x, derivative=True) / x ** (power - 1)
            )
        for k in keys:
            a, signed_l, m = h.keys[k]
            l = abs(signed_l)
            phase = (-1j if signed_l < 0 else 1) * (-1) ** m
            factor = 4 * np.pi * 1j**l
            if m:
                factor *= np.sqrt(2) / phase
            if l == 0:
                factor /= np.sqrt(4 * np.pi)
            dB[k, 1:] = (
                factor
                * (-x * spherical_jn(l, x, derivative=True))
                * y[lookup[l, m], 1:]
            )
        dA.append(da)
        dD.append(dd)
    dA = np.concatenate(dA, axis=1)
    dD = np.concatenate(dD, axis=0)
    radial, rho = h._split(prepared.input_density.vector())
    m = np.array([w @ c for w, c in zip(h.moments, radial)])
    z = np.array(
        [
            s.Znuc[a] / np.sqrt(4 * np.pi) if (l, mm) == (0, 0) else 0
            for a, l, mm in h.keys
        ]
    )
    multipoles = -h.S @ m + h.D @ rho + h.S @ z
    p = h.diagonal * rho + h.A @ multipoles
    dp = 2 * h.diagonal * rho + dA @ multipoles + h.A @ (dD @ rho)
    db = (dB @ p + h.B @ dp).real
    mask = strain_mask_derivative(
        h.reanalysis, b, s.pos, np.repeat(s.rmt, s.mult), h.volume, np.eye(3)
    )
    mt = [op.boundary * bound for op, bound in zip(h.radial, db)]
    ir = h.reanalysis.action(dp) + masked_potential_strain_source(h.reanalysis, p, mask)
    source = np.concatenate(mt + [ir.view(float)])
    # Central nonlinear Coulomb reconstruction is an independent provenance
    # check against native total potential after adding prepared XC.
    central = []
    boundaries = (h.B @ p).real
    for key, op, c, bb in zip(h.keys, h.radial, radial, boundaries):
        value = op.action(c, bb)
        a, l, mm = key
        if (l, mm) == (0, 0):
            r = h.integrations[a].r
            value = value - 2 * s.Znuc[a] * (1 - r / s.rmt[a])
        central.append(value)
    central = np.concatenate(central + [h.reanalysis.action(p).view(float)])
    return source, central, mask


def hydrostatic_xc_geometry(prepared, mask_derivative):
    """XC source at fixed native density coefficients.

    MT grids and angular axes are fixed, so their explicit geometry source is
    zero. In the IR native quadratic fit, first spatial derivatives scale as
    lambda^-1 and second derivatives as lambda^-2 at fixed density samples.
    The existing analytic XC formulas supply a directional complex-step
    derivative on the central branch. Native vacuum averaging, star truncation
    and moving-mask response follow in exactly that order.
    """
    h, xc = prepared.hartree, prepared.xc
    ir = xc.ir
    _, rho = h._split(prepared.input_density.vector())
    density = np.maximum(ir.expand(rho), 0)
    if ir.functional == "PBE":
        fields = ir.fields(density)
        direction = fields.copy()
        direction[..., 0] = 0
        direction[..., 1:4] *= -1
        direction[..., 4:] *= -2
        point = np.imag(ir._point(fields + 1j * 1e-24 * direction)) / 1e-24
    else:
        point = np.zeros(ir.shape)
    dp = ir.select(ir.average(point))
    pw = h.reanalysis.action(dp) + masked_potential_strain_source(
        h.reanalysis, ir.central_stars, mask_derivative
    )
    source = np.zeros(h.size)
    source[h.pw_offset :] = pw.view(float)
    central = np.zeros(h.size)
    for atom, slices in zip(xc.mt, xc.atom_slices):
        for sl, col in zip(slices, atom.central.T):
            central[sl] = col
    central[h.pw_offset :] = ir.central.view(float)
    return source, central


def hydrostatic_explicit_potential(prepared):
    """Return analytic source and diagnostics; no total-stress publication.

    tau_n is the SAME Density.transport derivative used to define the screened
    density coordinate. K*q must be added separately. In particular do not add
    a further r.grad(V) to this fixed-r source unless transforming all band
    geometry and functional derivatives to material coordinates as well.
    """
    tau = hydrostatic_transport_source(prepared.input_density).vector()
    transported = prepared.potential_operator.matvec(tau)
    coulomb, vc, mask = hydrostatic_coulomb_geometry(prepared)
    xc, vxc = hydrostatic_xc_geometry(prepared, mask)
    return transported + coulomb + xc, dict(
        transport=transported,
        coulomb_geometry=coulomb,
        xc_geometry=xc,
        central_coulomb=vc,
        central_xc=vxc,
        central_potential=vc + vxc,
        mask_derivative=mask,
    )
