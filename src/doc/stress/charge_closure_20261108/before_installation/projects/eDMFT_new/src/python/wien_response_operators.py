"""Discrete charge normalization and transport derivatives for WIEN densities.

These adapters use the exact coefficients and integration convention of the
existing finite-difference response driver. They supply genuine transpose
actions for the charge-assembly part of screening. They do not implement the
Hartree/XC, core-state, radial-basis, or general-strain source operators.
"""

from __future__ import annotations

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.special import spherical_jn

from stress_response_adjoint import PairedOperator
from wstruct import Struct


def charge_covectors(density, structure):
    """MT and interstitial charge covectors in Density.vector() coordinates.

    c_lm(r)=r^2 rho_lm(r). The existing charge assembler integrates the
    not-a-knot cubic spline of r*c_00 on the log-r grid. Differentiating this
    discrete integral, rather than substituting different quadrature weights,
    is necessary for agreement and a correct transpose. Fourier coefficients
    are star sums; theta includes every symmetry-equivalent sphere.
    """
    s = Struct()
    s.ReadStruct(str(structure))
    if s.lattyp.strip() != "P" or not np.allclose([s.alpha, s.beta, s.gamma], 90.0):
        raise ValueError(
            "charge adapter currently requires a primitive orthogonal cell"
        )
    vector = density.vector()
    mt, ir = np.zeros_like(vector), np.zeros_like(vector)
    offset = 0
    integrals = {}
    for (atom, ell, m), values in density.rad.items():
        n = len(values)
        if ell == m == 0:
            grid = density.grids[atom]
            key = tuple(grid)
            if key not in integrals:
                # Batched spline solves, not one nonlinear charge evaluation
                # per coefficient. Temporary O(Nrad^2) storage is released
                # after constructing this small, reusable quadrature vector.
                weights = CubicSpline(grid, np.eye(n), axis=0).integrate(
                    grid[0], grid[-1]
                )
                integrals[key] = weights * np.exp(grid)
            mt[offset : offset + n] = density.mult[atom] * integrals[key]
        offset += n
    cell = np.array([s.a, s.b, s.c])
    volume = float(np.prod(cell))
    for g in density.pw:
        g = np.asarray(g)
        magnitude = np.linalg.norm(2 * np.pi * g / cell)
        theta = complex(np.all(g == 0))
        position = 0
        for atom in range(s.nat):
            x = magnitude * s.rmt[atom]
            sphere = 1.0 if x == 0 else 3 * spherical_jn(1, x) / x
            for _ in range(s.mult[atom]):
                theta -= (
                    4
                    * np.pi
                    * s.rmt[atom] ** 3
                    / (3 * volume)
                    * sphere
                    * np.exp(-2j * np.pi * np.dot(s.pos[position], g))
                )
                position += 1
        ir[offset : offset + 2] = volume * np.array([theta.real, theta.imag])
        offset += 2
    return mt, ir


class InterstitialNormalization:
    """Differentiate the actual interstitial-only normalization, not uniform scaling.

    N(q) leaves all MT coefficients unchanged and multiplies the PW part p
    by alpha=(Ntarget-M(q))/I(q). At q0,
        dN = D dq - p0 ell^T dq,
        D = diag(1_MT, alpha0_PW),
        ell = (m + alpha0*i)/I0.
    Its transpose is D y - ell*(p0^T y). The alpha0 factors remain necessary
    away from an exactly normalized central output. Core density must already
    be included before applying this map.
    """

    def __init__(self, density, structure, target_charge):
        if not np.isfinite(target_charge):
            raise ValueError("target charge must be finite")
        self.base = density.vector().copy()
        self.mt, self.ir = charge_covectors(density, structure)
        self.target = float(target_charge)
        self.nrad = sum(len(v) for v in density.rad.values())
        m, i = self.mt @ self.base, self.ir @ self.base
        if i <= 0:
            raise ValueError("nonpositive central interstitial charge")
        self.alpha = (self.target - m) / i
        self.i0 = i
        self.diagonal = np.ones(len(self.base))
        self.diagonal[self.nrad :] = self.alpha
        self.p0 = self.base.copy()
        self.p0[: self.nrad] = 0
        self.ell = (self.mt + self.alpha * self.ir) / i
        self.operator = PairedOperator(
            len(self.base),
            len(self.base),
            lambda x: self.diagonal * x - self.p0 * np.dot(self.ell, x),
            lambda y: self.diagonal * y - self.ell * np.dot(self.p0, y),
            "interstitial charge normalization",
        )

    def endpoint(self, vector):
        """Independent nonlinear normalization used for directional checks."""
        vector = np.asarray(vector, float)
        if vector.shape != self.base.shape or not np.all(np.isfinite(vector)):
            raise ValueError("invalid normalization endpoint")
        i = self.ir @ vector
        if i <= 0:
            raise ValueError("nonpositive endpoint interstitial charge")
        result = vector.copy()
        result[self.nrad :] *= (self.target - self.mt @ vector) / i
        return result

    def explicit_geometry(self, mt_charge_derivative, ir_charge_derivative):
        """Normalization source from changing integration geometry at fixed q.

        The caller supplies M_nu(q0), I_nu(q0), derived with the same fixed-RMT
        convention. Ignoring these in a general strain source changes total N.
        This helper does not infer them from hydrostatic density transport.
        """
        return (
            -self.p0
            * (mt_charge_derivative + self.alpha * ir_charge_derivative)
            / self.i0
        )


def hydrostatic_transport_source(density, *, pullback=False):
    """Exact derivative of Density.transport at lambda=1, including MT edges.

    Pushforward: dc/deta=-c-dc/dlog(r), dPW/deta=-3 PW.
    Pullback has the opposite signs. This is the existing scalar affine
    transport convention; it is not an anisotropic/shear source generator.
    """
    result = density.copy()
    sign = 1.0 if pullback else -1.0
    for key, values in density.rad.items():
        grid = density.grids[key[0]]
        result.rad[key] = sign * (values + CubicSpline(grid, values)(grid, 1))
    for g, value in density.pw.items():
        result.pw[g] = 3 * sign * value
    return result
