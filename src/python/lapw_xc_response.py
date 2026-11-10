"""Prepared nonmagnetic LDA/PBE native XC response and exact transpose.

The optional energy_per_electron quantity differentiates the actual native
scalar-energy field. Its discrete derivative need not equal the separately
implemented potential; lapw_functional_cancellation assembles that difference.

Local partial derivatives are obtained by complex-step differentiation of
analytic formulas on the central active branch, once during preparation.
There are no subtractive finite differences or per-density-coordinate probes.
Repeated actions use the cached local Jacobian, native spline operations,
angular contractions and FFT-grid stencils. Meta-GGA/SOC/spin response and
arbitrary functional switches are deliberately outside this initial scope.
"""

import numpy as np
from scipy.linalg.lapack import dgttrf, dgttrs
from lapw_density_response import checked
from stress_response_adjoint import PairedOperator


def _gcor(rs):
    a = 0.0310907
    a1 = 0.21370
    q0 = -2 * a * (1 + a1 * rs)
    root = np.sqrt(rs)
    q1 = 2 * a * (7.5957 * root + 3.5876 * rs + 1.6382 * rs * root + 0.49294 * rs**2)
    log = np.log1p(1 / q1)
    q3 = a * (7.5957 / root + 2 * 3.5876 + 3 * 1.6382 * root + 4 * 0.49294 * rs)
    return q0 * log, -2 * a * a1 * log - q0 * q3 / (q1 * q1 + q1)


def unpolarized_xc(
    n,
    g,
    lap,
    gdgg,
    *,
    functional,
    charge_cut,
    exchange_gradient=None,
    exchange_gdgg=None,
    quantity="potential"
):
    """Native potential or energy/electron in Ry from local density invariants.

    Branch decisions use real central values. Preserve the built-in native
    low-density fallbacks; replacing them by another library's PBE would
    differentiate a different implemented map. Only standard XC_LDA/PBE.
    """
    if quantity not in ("potential", "energy_per_electron"):
        raise ValueError("invalid XC quantity")
    if functional not in ("LDA", "PBE"):
        raise NotImplementedError("only native LDA/PBE supported")
    n, g, lap, gdgg = np.broadcast_arrays(n, g, lap, gdgg)
    active = n.real > 2e-15
    rho = np.where(active, n, 1.0)
    rs = (0.75 / (np.pi * rho)) ** (1 / 3)
    ec, ecrs = _gcor(rs)
    vc = ec - rs * ecrs / 3
    lda_x = -2 * (3 / np.pi * np.where(rho.real >= 2e-12, rho, 2e-12)) ** (1 / 3)
    lda_e = 0.75 * lda_x
    if functional == "LDA":
        if quantity == "energy_per_electron":
            return np.where(
                active, lda_e + np.where(rho.real > charge_cut, 2 * ec, 0), 0
            )
        return np.where(active, lda_x + np.where(rho.real > charge_cut, 2 * vc, 0), 0)
    kf = (29.6088132032680758565 * rho) ** (1 / 3)
    rs = 1.91915829267751300662 / kf
    ec, ecrs = _gcor(rs)
    vc = ec - rs * ecrs / 3
    gx = g if exchange_gradient is None else exchange_gradient
    ux = gdgg if exchange_gdgg is None else exchange_gdgg
    s = gx / (2 * kf * rho)
    u = ux / (rho**2 * (2 * kf) ** 3)
    v = lap / (rho * (2 * kf) ** 2)
    mu = 0.2195149727645171
    kappa = 0.804
    ratio = mu / kappa
    p0 = 1 + ratio * s * s
    enhance = 1 + kappa - kappa / p0
    fs = 2 * mu / p0**2
    fss = -4 * ratio * s * fs / p0
    ex = -0.738558766382022405884 * rho ** (1 / 3)
    vx = ex * (4 / 3 * enhance - (u - 4 / 3 * s**3) * fss - v * fs)
    cut = (rho.real < 1e-6) & (abs(vx.real - lda_x.real) > 0.4 * abs(lda_x.real))
    vx = np.where(cut, lda_x, vx)
    beta = 0.06672455060314922
    gamma = 0.031090690869654895
    delta = beta / gamma
    ks = np.sqrt(4 * kf / np.pi)
    twoks = 2 * ks
    t = g / (twoks * rho)
    uu = gdgg / (rho**2 * twoks**3)
    vv = lap / (rho * twoks**2)
    b = delta / np.expm1(-ec / gamma)
    t2 = t * t
    t4 = t2 * t2
    q4 = 1 + b * t2
    q5 = 1 + b * t2 + b * b * t4
    h = gamma * np.log1p(delta * q4 * t2 / q5)
    bec = b * b * (delta / b + 1) / beta
    q8 = q5 * q5 + delta * q4 * q5 * t2
    q9 = 1 + 2 * b * t2
    hb = -beta * b * t4 * t2 * (2 + b * t2) / q8
    hrs = -rs / 3 * hb * bec * ecrs
    hbt = (
        2
        * beta
        * t4
        * (q4 * q5 * (2 * delta - 6 * b) - delta * (q5 * q9 + q4 * q9 * q9))
        / q8**2
    )
    hrst = rs / 3 * t2 * hbt * bec * ecrs
    ht = 2 * beta * q9 / q8
    htt = (
        4
        * beta
        * t
        * (
            2 * b / q8
            - q9
            * (2 * b * q5 * q9 + delta * (q4 * q5 + b * t2 * (q4 * q9 + q5)))
            / q8**2
        )
    )
    dvc = h + hrs + hrst + t2 * ht / 6 + 7 * t2 * t * htt / 6 - uu * htt - vv * ht
    pbe_c = vc + np.where(cut, 0, dvc)
    pbe_c = np.where((rho.real < 2e-6) & (abs(dvc.real) > 1.1 * abs(vc.real)), 0, pbe_c)
    result = 2 * (vx + pbe_c)
    if quantity == "energy_per_electron":
        # Native PBEX replaces EX by XCENER in its low-density fallback.
        # That function returns Ry, yet VXCLM2 doubles it alongside the usual
        # Hartree PBEX output. Preserve that implemented branch, as for VX.
        # The large-dvc PBEC safeguard suppresses the potential only, not EC+H.
        result = 2 * (np.where(cut, lda_e, ex * enhance) + ec + np.where(cut, 0, h))
        result = np.where(rho.real <= 2 * charge_cut, lda_e + 2 * ec, result)
        result = np.where(rho.real < 4e-9, lda_e + 2 * ec, result)
        return np.where(active, result, 0)
    # VXCLM2 returns the LDA limit before PBEX/PBEC when BOTH spins are low.
    result = np.where(rho.real <= 2 * charge_cut, lda_x + 2 * vc, result)
    result = np.where(rho.real < 4e-9, lda_x + 2 * vc, result)
    return np.where(active, result, 0)


def local_jacobian(function, fields):
    """Cache local analytic partials; branch tests are held at real inputs."""
    fields = np.asarray(fields, float)
    if not np.all(np.isfinite(fields)):
        raise ValueError("nonfinite XC fields")
    central = np.asarray(function(fields))
    jac = np.empty_like(fields)
    for j in range(fields.shape[-1]):
        h = 1e-24 * np.maximum(1, abs(fields[..., j]))
        x = fields.astype(complex)
        x[..., j] += 1j * h
        jac[..., j] = np.imag(function(x)) / h
    if not np.all(np.isfinite(jac)):
        raise ValueError("nonfinite local XC Jacobian")
    return central, jac


class NativeSplineDerivative:
    """Native SPLINE first derivative on arbitrary increasing abscissae.

    Factor once; O(N) forward/reverse. This is not a not-a-knot spline.
    Only the first derivative is needed here. The endpoint equations match
    the same four-point conditions used by native CoulInt.
    """

    def __init__(self, x):
        x = np.asarray(x, float)
        self.n = len(x)
        if (
            x.ndim != 1
            or len(x) < 4
            or np.any(np.diff(x) <= 0)
            or not np.all(np.isfinite(x))
        ):
            raise ValueError("invalid spline grid")
        self.h = np.diff(x)
        h = self.h
        diag = np.r_[-h[0], 2 * (h[:-1] + h[1:]), -h[-1]]
        *self.factor, info = dgttrf(h.copy(), diag, h.copy())
        if info:
            raise ValueError("spline factorization failed")
        self.left = (
            h[0] ** 2
            / (x[3] - x[0])
            * np.array([-1 / (x[2] - x[0]), 1 / (x[3] - x[1])])
        )
        self.right = (
            h[-1] ** 2
            / (x[-1] - x[-4])
            * np.array([1 / (x[-2] - x[-4]), -1 / (x[-1] - x[-3])])
        )

    def _solve(self, x, transpose=False):
        y, info = dgttrs(
            *self.factor,
            np.array(x, order="F", copy=True),
            trans="T" if transpose else "N"
        )
        if info:
            raise RuntimeError("spline solve failed")
        return y

    def action(self, x):
        x = np.asarray(x)
        single = x.ndim == 1
        x = x[:, None] if single else x
        if x.shape[0] != self.n:
            raise ValueError("spline dimension mismatch")
        slopes = np.diff(x, axis=0) / self.h[:, None]
        diff = np.diff(slopes, axis=0)
        rhs = np.zeros_like(x)
        rhs[1:-1] = diff
        rhs[0] = self.left @ diff[:2]
        rhs[-1] = self.right @ diff[-2:]
        c = self._solve(rhs)
        y = np.empty_like(x)
        y[:-1] = slopes - self.h[:, None] * (c[1:] + 2 * c[:-1])
        y[-1] = slopes[-1] + self.h[-1] * (c[-2] + 2 * c[-1])
        return y[:, 0] if single else y

    def pullback(self, y):
        y = np.asarray(y)
        single = y.ndim == 1
        y = y[:, None] if single else y
        c = np.zeros_like(y)
        slopes = y[:-1].copy()
        slopes[-1] += y[-1]
        c[:-1] -= 2 * self.h[:, None] * y[:-1]
        c[1:] -= self.h[:, None] * y[:-1]
        c[-2] += self.h[-1] * y[-1]
        c[-1] += 2 * self.h[-1] * y[-1]
        rhs = self._solve(c, True)
        diff = rhs[1:-1].copy()
        diff[:2] += self.left[:, None] * rhs[0]
        diff[-2:] += self.right[:, None] * rhs[-1]
        slopes[1:] += diff
        slopes[:-1] -= diff
        x = np.zeros_like(y)
        x[1:] += slopes / self.h[:, None]
        x[:-1] -= slopes / self.h[:, None]
        return x[:, 0] if single else x


class NativeDensityRadialDerivatives:
    """DFUNDR's parity-extended spline and near-origin replacements."""

    def __init__(self, r, ells):
        self.r = np.asarray(r, float)
        self.ell = np.asarray(ells)
        self.n = len(r)
        if (
            self.r.ndim != 1
            or self.n < 7
            or self.r[0] <= 0
            or np.any(np.diff(self.r) <= 0)
            or not np.all(np.isfinite(self.r))
        ):
            raise ValueError("invalid radial density grid")
        if self.ell.ndim != 1 or not np.issubdtype(self.ell.dtype, np.integer):
            raise ValueError("invalid harmonic degrees")
        self.spline = NativeSplineDerivative(np.r_[-self.r[::-1], self.r])
        self.sign = np.where(self.ell % 2 == 0, 1.0, -1.0)
        self.replacements = []
        for targets, source in [
            (range(3), np.arange(3, 6)),
            (range(5), np.arange(5, 7)),
        ]:
            rows = []
            for i in targets:
                row = np.ones(len(source))
                for j in range(len(source)):
                    for k in range(len(source)):
                        if k != j:
                            row[j] *= (r[i] - r[source[k]]) / (
                                r[source[j]] - r[source[k]]
                            )
                rows.append(row)
            self.replacements.append((len(rows), source, np.array(rows)))

    def _one(self, f, sign, which):
        ext = np.concatenate((f[::-1] * sign, f))
        y = self.spline.action(ext)[self.n :]
        n, indices, weights = self.replacements[which]
        y[:n] = weights @ y[indices]
        return y

    def _back(self, y, sign, which):
        y = np.array(y, copy=True)
        n, indices, weights = self.replacements[which]
        y[indices] += weights.T @ y[:n]
        y[:n] = 0
        full = np.zeros((2 * self.n, y.shape[1]))
        full[self.n :] = y
        x = self.spline.pullback(full)
        return x[self.n :] + x[: self.n][::-1] * sign

    def action(self, f):
        d1 = self._one(f, self.sign, 0)
        return d1, self._one(d1, -self.sign, 1)

    def pullback(self, y1, y2):
        return self._back(y1 + self._back(y2, -self.sign, 1), self.sign, 0)


def spherical_invariants(fields, r, sine, cosine):
    """GRANS in its local (r,theta,phi) frame, including the polar branch."""
    n, dr, drr, dt, dtt, dp, dpp, drt, drp, dtp = np.moveaxis(fields, -1, 0)
    pole = abs(sine) < 1e-10
    sn = np.where(pole, 1.0, sine)
    cot = cosine / sn
    gx = dr
    gy = dt / r
    gz = np.where(pole, 0, dp / (r * sn))
    square = gx * gx + gy * gy + gz * gz
    nonzero = square.real >= 1e-40
    mag = np.sqrt(np.where(nonzero, square, 1.0))
    lap = drr + 2 * dr / r + (dtt + np.where(pole, 0, cot * dt + dpp / sn**2)) / r**2
    ax = (
        dr * drr
        + dt * (drt - dt / r) / r**2
        + np.where(pole, 0, dp * (drp - dp / r) / (r * sn) ** 2)
    )
    ay = (
        dr * drt
        + dt * dtt / r**2
        + np.where(pole, 0, dp * (dtp - cot * dp) / (r * sn) ** 2)
    ) / r
    az = np.where(
        pole, 0, (dr * drp + dt * dtp / r**2 + dp * dpp / (r * sn) ** 2) / (r * sn)
    )
    gdgg = np.where(nonzero, (gx * ax + gy * ay + gz * az) / mag, 0)
    return n, np.where(nonzero, mag, 0), lap, gdgg


class MuffinTinXCResponse:
    """One atom: CLMSUM harmonics -> native vsp/vns potential harmonics.

    Angular tables come from the actual native quadrature. This avoids using
    a different Lebedev grid and preserves angular mean subtraction. The
    material loader/export is explicit, not silently inferred here.
    """

    def __init__(
        self,
        r,
        ells,
        central,
        angular,
        weights,
        *,
        functional,
        sine=None,
        cosine=None,
        angular_derivatives=None,
        harmonics=None,
        quantity="potential",
        cubic=False,
    ):
        self.r = np.asarray(r, float)
        self.ells = np.asarray(ells, int)
        self.functional = functional
        self.quantity = quantity
        self.nr = len(r)
        self.nl = len(ells)
        self.shape = (self.nr, self.nl)
        if (
            self.nr < 7
            or self.nl < 1
            or self.r[0] <= 0
            or not np.all(np.isfinite(self.r))
            or np.any(np.diff(self.r) <= 0)
        ):
            raise ValueError("invalid MT radial grid")
        self.harmonics = None if harmonics is None else tuple(map(tuple, harmonics))
        if self.harmonics is not None and (
            len(self.harmonics) != self.nl
            or any(
                l != e or int(l) != l or int(m) != m or not 0 <= m <= abs(l)
                for (l, m), e in zip(self.harmonics, self.ells)
            )
        ):
            raise ValueError("invalid signed harmonic labels")
        center = checked(central, self.shape, "central MT density", real=True)
        self.cubic = bool(cubic)
        self.native_angular = np.array(angular, float)
        self.native_angular_derivatives = angular_derivatives
        self.Y = self.native_angular.copy()
        self.weights = np.asarray(weights, float)
        if not np.all(np.isfinite(self.Y)) or not np.all(np.isfinite(self.weights)):
            raise ValueError("nonfinite angular quadrature")
        if self.Y.shape != (len(weights), self.nl) or self.ells[0] != 0:
            raise ValueError("invalid MT angular layout")
        self.scale = np.ones(self.nl)
        self.scale[0] = 1 / np.sqrt(4 * np.pi)
        projection_y = self.Y
        if self.cubic:
            from lapw_cubic_harmonics import cubic_maps
            if self.harmonics is None:
                raise ValueError("cubic XC requires signed harmonic labels")
            c, d = cubic_maps(self.harmonics)
            # Exported YY has the compressed entries at its front; the tail
            # is workspace and must not be interpreted as ordinary harmonics.
            self.Y = self.native_angular[:, :len(c)] @ c
            projection_y = self.native_angular[:, :len(c)] @ d.T
        yw = projection_y * self.weights[:, None]
        self.project = yw - yw[:, 0, None] * np.sum(yw, axis=0) / np.sum(yw[:, 0])
        self.project[:, 0] += yw[:, 0]
        self.output_scale = np.ones(self.shape)
        self.output_scale[:, 0] = self.r / np.sqrt(4 * np.pi)
        if functional == "PBE":
            if angular_derivatives is None or sine is None or cosine is None:
                raise ValueError("PBE requires native angular derivatives")
            self.angular = [self.Y] + [
                (np.asarray(a, float)[:, :len(c)] @ c if self.cubic else np.asarray(a, float))
                for a in angular_derivatives
            ]
            if len(self.angular) != 6 or any(
                a.shape != self.Y.shape or not np.all(np.isfinite(a))
                for a in self.angular
            ):
                raise ValueError("invalid angular derivative shape")
            self.sine = checked(sine, (len(weights),), "polar sine", real=True)[None, :]
            self.cosine = checked(cosine, (len(weights),), "polar cosine", real=True)[
                None, :
            ]
            if not np.allclose(
                self.sine**2 + self.cosine**2, 1, rtol=1e-10, atol=1e-12
            ):
                raise ValueError("invalid angular polar coordinates")
            self.radial = NativeDensityRadialDerivatives(self.r, abs(self.ells))
        elif functional != "LDA":
            raise NotImplementedError("unsupported XC functional")
        fields = self.fields(center)
        self.central_points, self.jac = local_jacobian(self._point, fields)
        self.central = (self.central_points @ self.project) * self.output_scale
        self.operator = PairedOperator(
            center.size,
            center.size,
            lambda x: self.action(x.reshape(self.shape)).ravel(),
            lambda y: self.pullback(y.reshape(self.shape)).ravel(),
            "MT " + functional + " XC",
        )

    def fields(self, c):
        f = c / self.r[:, None] ** 2 * self.scale
        if self.functional == "LDA":
            return (f @ self.Y.T)[..., None]
        y, dt, dtt, dp, dpp, dtp = self.angular
        d1, d2 = self.radial.action(f)
        return np.stack(
            (
                f @ y.T,
                d1 @ y.T,
                d2 @ y.T,
                f @ dt.T,
                f @ dtt.T,
                f @ dp.T,
                f @ dpp.T,
                d1 @ dt.T,
                d1 @ dp.T,
                f @ dtp.T,
            ),
            axis=-1,
        )

    def _point(self, x):
        args = (
            (x[..., 0], 0, 0, 0)
            if self.functional == "LDA"
            else spherical_invariants(x, self.r[:, None], self.sine, self.cosine)
        )
        return unpolarized_xc(
            *args, functional=self.functional, charge_cut=1e-8, quantity=self.quantity
        )

    def action(self, c):
        x = checked(c, self.shape, "MT density perturbation", real=True)
        return (
            np.sum(self.jac * self.fields(x), axis=-1) @ self.project
        ) * self.output_scale

    def pullback(self, y):
        y = checked(y, self.shape, "MT potential covector", real=True)
        point = (y * self.output_scale) @ self.project.T
        z = np.moveaxis(point[..., None] * self.jac, -1, 0)
        if self.functional == "LDA":
            f = z[0] @ self.Y
        else:
            yy, dt, dtt, dp, dpp, dtp = self.angular
            f = z[0] @ yy + z[3] @ dt + z[4] @ dtt + z[5] @ dp + z[6] @ dpp + z[9] @ dtp
            f += self.radial.pullback(z[1] @ yy + z[7] @ dt + z[8] @ dp, z[2] @ yy)
        return f / self.r[:, None] ** 2 * self.scale


def cartesian_invariants(fields):
    """PWXAD5 invariants from n and the nine quadratic-fit derivatives.

    The stored six Hessian entries are xx,yy,zz,yz,zx,xy (true derivatives,
    not polynomial coefficients). At a vanishing gradient the native branch
    avoids dividing by its norm. Its PBE potential has a finite limit there.
    """
    n, gx, gy, gz, xx, yy, zz, yz, zx, xy = np.moveaxis(fields, -1, 0)
    square = gx * gx + gy * gy + gz * gz
    active = square.real > 1e-28
    mag = np.sqrt(np.where(active, square, 1.0))
    contraction = (
        xx * gx * gx
        + yy * gy * gy
        + zz * gz * gz
        + 2 * (xy * gx * gy + yz * gy * gz + zx * gz * gx)
    )
    gdgg = np.where(
        active, contraction / mag, (gx * xx + gy * yy + gz * zz) / np.sqrt(3)
    )
    return n, np.where(active, mag, 0), xx + yy + zz, gdgg


class InterstitialXCResponse:
    """Native periodic density-star -> masked potential-star XC derivative.

    A quadratic fit on the native lattice-specific neighbours reproduces PWXAD4/5.
    Preparation contracts the local PBE Jacobian with that fit, leaving a
    variable-coefficient stencil for each forward/reverse call.
    Both positivity and low-density averaging follow the central branch.
    GETFF1 selects star representatives BEFORE truncation/re-expansion and
    interstitial masking; replacing it by a member average changes this map.
    """

    def __init__(
        self,
        reanalysis,
        central,
        cell_lengths,
        *,
        functional,
        rho_min=(2e-8, 1e-9),
        quantity="potential",
        lattice_type="P",
    ):
        self.reanalysis = reanalysis
        self.stars = reanalysis.stars
        self.shape = reanalysis.shape
        self.nstar = self.stars.nstar
        self.functional = functional
        self.quantity = quantity
        cell = np.asarray(cell_lengths)
        if cell.shape == (3,):
            if np.any(cell <= 0):
                raise ValueError("cell lengths must be positive")
            cell = np.diag(cell)
        cell = checked(cell, (3, 3), "Cartesian FFT cell", real=True)
        if abs(np.linalg.det(cell)) < 1e-14:
            raise ValueError("singular FFT cell")
        self.rep_indices = tuple((self.stars.representatives % np.array(self.shape)).T)
        center = checked(central, (self.nstar,), "central density stars")
        raw = self.expand(center)
        self.positive = raw > 0
        density = np.maximum(raw, 0)
        if (
            len(rho_min) != 2
            or not np.all(np.isfinite(rho_min))
            or not 0 <= rho_min[1] < rho_min[0]
        ):
            raise ValueError("invalid vacuum averaging thresholds")
        self.average_from = (raw / 2 > rho_min[1]) & (raw / 2 < rho_min[0])
        self.average_to = raw / 2 < rho_min[0]
        self.average_count = int(np.sum(self.average_from))
        # PWXAD4 uses the 120-degree planar metric for H cells. Its
        # neighbour set is independent of the CornUse switch in that case.
        # Cartesian derivatives still use AVEC/FFT-size, never fractional
        # coordinate lengths. Keep native loop ordering for reproducibility.
        extent = 2 if lattice_type.strip() == "H" else 1
        self.offsets = np.array([
            (i, j, k)
            for j in range(-extent, extent + 1)
            for i in range(-extent, extent + 1)
            for k in (-1, 0, 1)
            if lattice_type.strip() != "H" or i*i + j*j - i*j + k*k < 3
        ])
        if functional == "PBE":
            xyz = (self.offsets / np.array(self.shape)) @ cell.T
            x, y, z = xyz.T
            design = np.array(
                [np.ones(len(self.offsets)), x, y, z, x * x, y * y, z * z, y * z, z * x, x * y]
            )
            fit = np.linalg.solve(design @ design.T, design)[1:]
            fit[3:6] *= 2
            # Neighbour-minus-centre form avoids differentiating a fitted
            # constant and is exactly the native convention.
            self.fit = fit
            fields = self.fields(density)
            central_points, jac = local_jacobian(self._point, fields)
            self.stencil = np.einsum("...i,ij->j...", jac[..., 1:], fit)
            self.direct = jac[..., 0] - np.sum(self.stencil, axis=0)
        elif functional == "LDA":
            central_points, jac = local_jacobian(self._point, density[..., None])
            self.stencil = None
            self.direct = jac[..., 0]
        else:
            raise NotImplementedError("unsupported XC functional")
        self.central_points = self.average(central_points)
        self.central_stars = self.select(self.central_points)
        self.central = self.reanalysis.action(self.central_stars)
        self.operator = PairedOperator(
            2 * self.nstar,
            2 * self.nstar,
            lambda x: self.action(np.ascontiguousarray(x).view(complex)).view(float),
            lambda y: self.pullback(np.ascontiguousarray(y).view(complex)).view(float),
            "interstitial " + functional + " XC",
        )

    def expand(self, stars):
        stars = checked(stars, (self.nstar,), "density stars")
        grid = np.zeros(self.shape, complex)
        grid[self.reanalysis.indices] = (
            stars[self.stars.owners] * self.reanalysis.weight
        )
        return np.fft.fftn(grid).real

    def expand_back(self, grid):
        member = (self.reanalysis.points * np.fft.ifftn(grid))[
            self.reanalysis.indices
        ] * self.reanalysis.weight.conj()
        return np.bincount(
            self.stars.owners, member.real, minlength=self.nstar
        ) + 1j * np.bincount(self.stars.owners, member.imag, minlength=self.nstar)

    def select(self, grid):
        return self.stars.sizes * np.fft.ifftn(grid)[self.rep_indices]

    def select_back(self, y):
        grid = np.zeros(self.shape, complex)
        grid[self.rep_indices] = y * self.stars.sizes
        return np.fft.fftn(grid).real / self.reanalysis.points

    def fields(self, density):
        result = np.zeros(self.shape + (10,))
        result[..., 0] = density
        for off, weights in zip(self.offsets, self.fit.T):
            delta = np.roll(density, tuple(-off), axis=(0, 1, 2)) - density
            result[..., 1:] += delta[..., None] * weights
        return result

    def _point(self, fields):
        if self.functional == "LDA":
            return unpolarized_xc(
                fields[..., 0],
                0,
                0,
                0,
                functional="LDA",
                charge_cut=1e-9,
                quantity=self.quantity,
            )
        args = cartesian_invariants(fields)
        # Native spin gradients are half the total gradient. The tiny-gradient
        # branch therefore activates at a factor-four different squared norm
        # in PBEX and PBEC, even in an unpolarized calculation.
        square = np.sum(fields[..., 1:4] ** 2, axis=-1)
        active = square.real > 4e-28
        exg = np.where(active, args[1], 0)
        exu = np.where(
            active,
            args[3],
            np.sum(fields[..., 1:4] * fields[..., 4:7], axis=-1) / np.sqrt(3),
        )
        return unpolarized_xc(
            *args,
            functional="PBE",
            charge_cut=1e-9,
            exchange_gradient=exg,
            exchange_gdgg=exu,
            quantity=self.quantity
        )

    def average(self, x, transpose=False):
        x = np.array(x, copy=True)
        # XCPOT3 averages V_xc in vacuum points; UFFT's energy per electron
        # is left unchanged. The scalar-energy derivative must preserve this.
        if self.quantity == "energy_per_electron":
            return x
        if self.average_count:
            if transpose:
                total = np.sum(x[self.average_to])
                x[self.average_to] = 0
                x[self.average_from] += total / self.average_count
            else:
                x[self.average_to] = np.mean(x[self.average_from])
        return x

    def action(self, x):
        density = self.expand(x) * self.positive
        result = self.direct * density
        if self.stencil is not None:
            for off, coef in zip(self.offsets, self.stencil):
                result += coef * np.roll(density, tuple(-off), axis=(0, 1, 2))
        return self.reanalysis.action(self.select(self.average(result)))

    def pullback(self, y):
        local = self.average(self.select_back(self.reanalysis.pullback(y)), True)
        density = self.direct * local
        if self.stencil is not None:
            for off, coef in zip(self.offsets, self.stencil):
                density += np.roll(coef * local, tuple(off), axis=(0, 1, 2))
        return self.expand_back(density * self.positive)


class PeriodicXCResponse:
    """Assemble prepared MT/IR maps in PeriodicHartreeResponse coordinates.

    All atoms and all density harmonics must be present. The Hartree object
    supplies an explicitly validated common layout and native FFT geometry;
    this constructor neither runs LAPW0 nor changes any production file.
    """

    def __init__(self, hartree, muffin_tins, interstitial):
        self.layout = hartree
        self.size = hartree.size
        self.mt = tuple(muffin_tins)
        self.ir = interstitial
        quantities = {op.quantity for op in self.mt} | {self.ir.quantity}
        if len(quantities) != 1:
            raise ValueError("cannot combine potential and energy XC fields")
        self.quantity = quantities.pop()
        if (
            len(self.mt) != len(hartree.integrations)
            or self.ir.reanalysis is not hartree.reanalysis
        ):
            raise ValueError("XC/Hartree geometry mismatch")
        self.atom_slices = []
        for a, op in enumerate(self.mt):
            rows = [
                (key, sl)
                for key, sl in zip(hartree.keys, hartree.slices)
                if key[0] == a
            ]
            if (
                op.harmonics is None
                or tuple(key[1:] for key, sl in rows) != op.harmonics
            ):
                raise ValueError("XC requires exact signed (l,m) channel ordering")
            if (
                len(rows) != op.nl
                or not np.array_equal([abs(key[1]) for key, sl in rows], abs(op.ells))
                or len(hartree.integrations[a].r) != op.nr
                or not np.allclose(hartree.integrations[a].r, op.r, rtol=1e-12, atol=0)
            ):
                raise ValueError("XC radial layout mismatch")
            self.atom_slices.append([sl for key, sl in rows])
        self.operator = PairedOperator(
            self.size, self.size, self.action, self.pullback, "periodic XC"
        )
        self._hxc_operator = PairedOperator(
            self.size,
            self.size,
            lambda x: hartree.action(x) + self.action(x),
            lambda y: hartree.pullback(y) + self.pullback(y),
            "periodic Hartree + XC",
        )

    @property
    def hxc_operator(self):
        if self.quantity != "potential":
            raise ValueError("an XC energy field is not an Hxc potential response")
        return self._hxc_operator

    def _map(self, x, reverse):
        x = checked(x, (self.size,), "XC field", real=True)
        result = np.zeros(self.size)
        for op, slices in zip(self.mt, self.atom_slices):
            radial = np.array([x[sl] for sl in slices]).T
            out = op.pullback(radial) if reverse else op.action(radial)
            for sl, col in zip(slices, out.T):
                result[sl] = col
        z = np.ascontiguousarray(x[self.layout.pw_offset :]).view(complex)
        result[self.layout.pw_offset :] = (
            self.ir.pullback(z) if reverse else self.ir.action(z)
        ).view(float)
        return result

    def action(self, x):
        return self._map(x, False)

    def pullback(self, x):
        return self._map(x, True)
