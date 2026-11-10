"""Prepared radial part of LAPW0 Hartree response, with an exact transpose.

This implements the discrete CoulEM (CoulInt) branch, including its native
cubic-spline endpoint conditions and Euler--Maclaurin corrections. It is not
a replacement quadrature. It is also not the whole periodic Hartree operator:
the caller must supply the boundary potential from the periodic multipole/PW
solution. Geometry and RMT are fixed. There is no XC or core response here.
"""

import numpy as np
from scipy.linalg.lapack import dgttrf, dgttrs
from stress_response_adjoint import PairedOperator


class NativeCoulombIntegral:
    """F -> integral from the first radial point to each grid point.

    The input is F(r), integrated with dr, not with dlog(r). Like CoulInt,
    this map returns zero at the first point and adds no origin term. The
    native spline solves a tridiagonal system for one sixth of the second
    derivative, with endpoint RHS derived from four-point differences.
    These are NOT the default SciPy not-a-knot boundary conditions.

    One pivoted tridiagonal LU is prepared. Actions and transposes cost O(Nr)
    per field and accept either (Nr,) or batched (Nr,Nfield) real arrays.
    There is no dense quadrature matrix or coordinate-by-coordinate probing.
    """

    def __init__(self, radius, *, log_step=None):
        r = np.asarray(radius)
        if (
            r.ndim != 1
            or len(r) < 4
            or np.iscomplexobj(r)
            or not np.all(np.isfinite(r))
            or np.any(r <= 0)
            or np.any(np.diff(r) <= 0)
        ):
            raise ValueError("need at least four positive increasing radial points")
        self.r = np.array(r, float, copy=True)
        self.n = len(r)
        self.dx = float(
            np.log(r[-1] / r[0]) / (self.n - 1) if log_step is None else log_step
        )
        if (
            not np.isfinite(self.dx)
            or self.dx <= 0
            or not np.allclose(np.diff(np.log(r)), self.dx, rtol=3e-12, atol=5e-15)
        ):
            raise ValueError("CoulInt requires a uniform logarithmic radial grid")
        # Match the native abscissae, including floating-point spacing.
        x = np.arange(1, self.n + 1) * self.dx
        self.h = np.diff(x)
        diagonal = np.r_[-self.h[0], 2 * (self.h[:-1] + self.h[1:]), -self.h[-1]]
        *self.factors, info = dgttrf(self.h.copy(), diagonal, self.h.copy())
        if info != 0:
            raise ValueError(f"native spline factorization failed: {info}")
        self.left = (
            self.h[0] ** 2
            / (x[3] - x[0])
            * np.array([-1 / (x[2] - x[0]), 1 / (x[3] - x[1])])
        )
        self.right = (
            self.h[-1] ** 2
            / (x[-1] - x[-4])
            * np.array([1 / (x[-2] - x[-4]), -1 / (x[-1] - x[-3])])
        )
        for a in [self.r, self.h, self.left, self.right, *self.factors]:
            a.setflags(write=False)
        self.operator = PairedOperator(
            self.n, self.n, self.action, self.pullback, "native radial Coulomb integral"
        )

    def _columns(self, value):
        a = np.asarray(value)
        if (
            a.ndim not in (1, 2)
            or a.shape[0] != self.n
            or (a.ndim == 2 and a.shape[1] == 0)
            or np.iscomplexobj(a)
            or not np.all(np.isfinite(a))
        ):
            raise ValueError("invalid real radial field")
        return (a[:, None] if a.ndim == 1 else a), a.ndim == 1

    def _solve(self, rhs, *, transpose=False):
        solution, info = dgttrs(
            *self.factors,
            np.array(rhs, float, order="F", copy=True),
            trans="T" if transpose else "N",
            overwrite_b=True,
        )
        if info != 0:
            raise RuntimeError(f"native spline solve failed: {info}")
        return solution

    def _spline(self, f):
        slope = np.diff(f, axis=0) / self.h[:, None]
        rhs = np.zeros_like(f)
        rhs[1:-1] = np.diff(slope, axis=0)
        rhs[0] = self.left @ rhs[1:3]
        rhs[-1] = self.right @ rhs[-3:-1]
        q = self._solve(rhs)
        first = np.empty_like(f)
        first[:-1] = slope - self.h[:, None] * (2 * q[:-1] + q[1:])
        first[-1] = slope[-1] + self.h[-1] * (q[-2] + 2 * q[-1])
        cubic = np.empty_like(f)
        cubic[:-1] = np.diff(q, axis=0) / self.h[:, None]
        cubic[-1] = cubic[-2]
        return first, cubic

    def _spline_pullback(self, first_bar, cubic_bar):
        slope_bar = first_bar[:-1].copy()
        slope_bar[-1] += first_bar[-1]
        qb = np.zeros_like(first_bar)
        hb = self.h[:, None] * first_bar[:-1]
        qb[:-1] -= 2 * hb
        qb[1:] -= hb
        qb[-2] += self.h[-1] * first_bar[-1]
        qb[-1] += 2 * self.h[-1] * first_bar[-1]
        cb = cubic_bar[:-1].copy()
        cb[-1] += cubic_bar[-1]
        cb /= self.h[:, None]
        qb[1:] += cb
        qb[:-1] -= cb
        rb = self._solve(qb, transpose=True)
        rb[1:3] += self.left[:, None] * rb[0]
        rb[-3:-1] += self.right[:, None] * rb[-1]
        slope_bar[1:] += rb[1:-1]
        slope_bar[:-1] -= rb[1:-1]
        slope_bar /= self.h[:, None]
        fb = np.zeros_like(first_bar)
        fb[1:] += slope_bar
        fb[:-1] -= slope_bar
        return fb

    def action(self, field, *, tail=False):
        y, squeeze = self._columns(field)
        f = self.r[:, None] * y
        first, cubic = self._spline(f)
        if tail:
            # I(R)-I(r), without subtracting two large prefix integrals.
            # This matters when high-l radial kernels grow near the origin.
            trap = np.cumsum(f[::-1], axis=0)[::-1] - 0.5 * (f[-1] + f)
            result = self.dx * (
                trap
                - self.dx / 12 * (first[-1] - first)
                + self.dx**3 / 120 * (cubic[-1] - cubic)
            )
            result[-1] = 0
        else:
            trap = np.cumsum(f, axis=0) - 0.5 * (f[0] + f)
            result = self.dx * (
                trap
                - self.dx / 12 * (first - first[0])
                + self.dx**3 / 120 * (cubic - cubic[0])
            )
            result[0] = 0
        return result[:, 0] if squeeze else result

    def pullback(self, covector, *, tail=False):
        y, squeeze = self._columns(covector)
        trap_bar = self.dx * y
        if tail:
            fb = np.cumsum(trap_bar, axis=0) - 0.5 * trap_bar
            endpoint, sign = -1, -1
        else:
            fb = np.cumsum(trap_bar[::-1], axis=0)[::-1] - 0.5 * trap_bar
            endpoint, sign = 0, 1
        fb[endpoint] -= 0.5 * np.sum(trap_bar, axis=0)
        first_bar = -sign * self.dx**2 / 12 * y
        first_bar[endpoint] -= np.sum(first_bar, axis=0)
        cubic_bar = sign * self.dx**4 / 120 * y
        cubic_bar[endpoint] -= np.sum(cubic_bar, axis=0)
        result = self.r[:, None] * (fb + self._spline_pullback(first_bar, cubic_bar))
        return result[:, 0] if squeeze else result


class RadialHartreeDirichlet:
    """One MT harmonic: density and boundary value -> native potential output.

    Input density is c_00=4*pi*r^2*rho_spherical for ell=0, and c_L=r^2*rho_L
    otherwise, matching Density.vector(). The boundary input is delta V_0(R)
    in Ry for ell=0, and delta V_L(R) otherwise. Output is r*delta V_0 for
    ell=0 (vsp), and delta V_L for ell>0 (vns). Signed-l sine/cosine harmonics
    share this radial map; their angular rotations occur outside it.

    Only the CoulEM branch is supplied. Nuclear, background-charge and strain
    geometry sources are excluded. The boundary input is REQUIRED when this
    is connected to the full periodic Hartree operator; zero is appropriate
    only for the isolated Dirichlet control, not as a substitute for screening.
    """

    def __init__(self, integration, ell):
        if not isinstance(integration, NativeCoulombIntegral):
            raise TypeError("explicit native Coulomb integration context required")
        if int(ell) != ell or ell < 0:
            raise ValueError("invalid angular momentum")
        self.integration = integration
        self.ell = int(ell)
        r = integration.r
        try:
            with np.errstate(over="raise", divide="raise", invalid="raise"):
                self.a = r**self.ell
                self.b = (1 - (r / r[-1]) ** (2 * self.ell + 1)) / r ** (self.ell + 1)
        except FloatingPointError as exc:
            raise ValueError(
                "radial powers overflow for this grid/angular momentum"
            ) from exc
        self.scale = (
            2 * r if self.ell == 0 else np.full(len(r), 8 * np.pi / (2 * self.ell + 1))
        )
        self.boundary = r if self.ell == 0 else (r / r[-1]) ** self.ell
        for a in [self.a, self.b, self.scale, self.boundary]:
            a.setflags(write=False)
        self.operator = PairedOperator(
            len(r) + 1,
            len(r),
            self._packed_action,
            self._packed_pullback,
            "native MT radial Hartree with explicit boundary",
        )

    def action(self, density, boundary_value):
        c = np.asarray(density)
        if (
            c.shape != (self.integration.n,)
            or np.iscomplexobj(c)
            or not np.all(np.isfinite(c))
        ):
            raise ValueError("invalid radial density")
        if (
            np.ndim(boundary_value) != 0
            or np.iscomplexobj(boundary_value)
            or not np.isfinite(boundary_value)
        ):
            raise ValueError("explicit real boundary value required")
        volume = self.b * self.integration.action(
            self.a * c
        ) + self.a * self.integration.action(self.b * c, tail=True)
        return self.scale * volume + self.boundary * boundary_value

    def pullback(self, covector):
        y, _ = self.integration._columns(covector)
        if y.shape[1] != 1:
            raise ValueError("one radial covector required")
        y = y[:, 0]
        v = self.scale * y
        back = self.a * self.integration.pullback(
            self.b * v
        ) + self.b * self.integration.pullback(self.a * v, tail=True)
        return back, float(self.boundary @ y)

    def _packed_action(self, vector):
        return self.action(vector[:-1], vector[-1])

    def _packed_pullback(self, vector):
        density, boundary = self.pullback(vector)
        return np.r_[density, boundary]
