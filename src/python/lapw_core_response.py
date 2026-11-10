"""Discrete radial core response with matched inward/outward propagation.

This differentiates the five-point predictor/corrector used by native LCORE,
including the origin series, outer boundary, eigenvalue matching and CHARGE
normalization. A prepared state supplies its native matching/end indices.
The active discrete branch is held fixed; no new eigenvalue search is done
in a response call. The small component is scaled by c for NREL stability.
"""

import numpy as np
from lapw_density_response import checked
from lapw_hartree_periodic import native_charge_weights
from stress_response_adjoint import PairedOperator


class CorePotentialExtension:
    """Native vsp Ry -> shifted/extrapolated Hartree potential and transpose.

    The fixed case.inc shift is subtracted before conversion. Outside RMT the
    a/r+b/r^2 fit matches the last value and two-point slope; the monotonic
    clamp follows the central active branch. No derivative of a fixed shift
    is included. A branch exactly at its kink is reported to the caller.
    """

    def __init__(self, radius, vsp, *, shift_ry=0, extra_points=None):
        r = np.asarray(radius, float)
        self.n = len(r)
        native_charge_weights(r)  # logarithmic-grid validation
        v = checked(vsp, (len(r),), "spherical potential", real=True)
        if extra_points is None:
            extra_points = 140 if len(r) % 2 else 141
        if not np.isfinite(shift_ry):
            raise ValueError("nonfinite core shift")
        if int(extra_points) != extra_points or extra_points < 0:
            raise ValueError("invalid core extension length")
        extra_points = int(extra_points)
        dx = np.log(r[-1] / r[0]) / (len(r) - 1)
        self.r = r[0] * np.exp(dx * np.arange(len(r) + extra_points))
        self.scale = 1 / (2 * r)
        base = (v - shift_ry * r) * self.scale
        R = r[-1]
        dr = R - r[-2]
        tail = self.r[len(r) :]
        # Coefficients of V[-2] and V[-1] in the native exterior fit.
        self.fit = np.column_stack(
            (
                -(R**2) / (dr * tail) + R**3 / (dr * tail**2),
                (2 * R + R**2 / dr) / tail + (-(R**2) - R**3 / dr) / tail**2,
            )
        )
        self.clamped = np.zeros(len(tail), bool)
        self.kink_margin = np.inf
        central = list(base)
        for i, row in enumerate(self.fit):
            candidate = row @ base[-2:]
            previous = central[-1]
            self.kink_margin = min(self.kink_margin, abs(candidate - previous))
            self.clamped[i] = candidate < previous
            central.append(previous if self.clamped[i] else candidate)
        self.central = np.array(central)
        self.operator = PairedOperator(
            self.n, len(self.r), self.action, self.pullback, "core potential extension"
        )

    def action(self, v):
        v = checked(v, (self.n,), "core potential perturbation", real=True) * self.scale
        result = np.empty(len(self.r))
        result[: self.n] = v
        for i, row in enumerate(self.fit):
            result[self.n + i] = (
                result[self.n + i - 1] if self.clamped[i] else row @ v[-2:]
            )
        return result

    def pullback(self, y):
        z = np.array(
            checked(y, (len(self.r),), "core potential covector", real=True), copy=True
        )
        for i in range(len(self.fit) - 1, -1, -1):
            if self.clamped[i]:
                z[self.n + i - 1] += z[self.n + i]
            else:
                z[self.n - 2 : self.n] += self.fit[i] * z[self.n + i]
        return z[: self.n] * self.scale


def _origin(r, z, kappa, c, energy, potential0, tolerance=1e-13):
    """First five native series values/derivatives and d/d(V0-E)."""
    ci2 = 1 / c**2
    gamma = np.sqrt(kappa * kappa - z * z * ci2)
    a = potential0 + z / r[0] - energy
    coeff = (
        np.array([1.0, z / (kappa - gamma)])
        if kappa < 0
        else np.array([z * ci2 / (kappa + gamma), 1.0])
    )
    dc = np.zeros(2)
    p = np.tile(coeff, (5, 1))
    dp = np.zeros_like(p)
    f = gamma * p
    df = np.zeros_like(p)
    for n in range(1, 21):
        t = n + gamma
        den = n * (n + 2 * gamma)
        mat = (
            np.array(
                [
                    [z * ci2 * a, (t - kappa) * (2 - a * ci2)],
                    [(t + kappa) * a, -z * (2 - a * ci2)],
                ]
            )
            / den
        )
        md = np.array([[z * ci2, -(t - kappa) * ci2], [t + kappa, z * ci2]]) / den
        dc, coeff = mat @ dc + md @ coeff, mat @ coeff
        term = r[:5, None] ** n * coeff
        dterm = r[:5, None] ** n * dc
        stop = n > 1 and np.any(
            np.all(abs(term / np.where(p != 0, p, 1)) <= tolerance, axis=1)
        )
        p += term
        dp += dterm
        f += t * term
        df += t * dterm
        if stop:
            break
    else:
        raise ValueError("core origin series did not converge")
    power = r[:5, None] ** gamma
    return p * power, f * power, dp * power, df * power


def _outer(r, c, e):
    decay = np.sqrt(-e * (2 + e / c**2))
    dd = -(1 + e / c**2) / decay
    ratio = -decay / (2 + e / c**2)
    dratio = -dd / (2 + e / c**2) + decay / c**2 / (2 + e / c**2) ** 2
    p = np.exp(-decay * r)
    dp = -dd * r * p
    y = np.column_stack((p, ratio * p))
    dy = np.column_stack((dp, dratio * p + ratio * dp))
    f = -decay * r[:, None] * y
    df = -dd * r[:, None] * y - decay * r[:, None] * dy
    return y, f, dy, df


class _CoreTrajectory:
    """Cached 12-state recurrence; forward and reverse cost O(Nrad)."""

    def __init__(self, r, v, e, kappa, c, initial_indices, indices, initial):
        self.r = r
        self.indices = np.array(indices, int)
        self.initial_indices = np.array(initial_indices, int)
        y, f, dy, df = initial
        self.initial_derivative = np.r_[dy[-1], df.ravel()]
        self.initial_y_derivative = dy
        self.values = np.zeros((len(r), 2))
        self.values[self.initial_indices] = y
        state = np.r_[y[-1], f.ravel()]
        dx = np.log(r[1] / r[0])
        dm = dx / 720 * (1 if self.indices[0] > initial_indices[-1] else -1)
        alpha = 475 / 502
        beta = 27 / 502
        pred = np.column_stack(
            (
                np.eye(2),
                np.kron(np.array([[251, -1274, 2616, -2774, 1901]]) * dm, np.eye(2)),
            )
        )
        base = np.column_stack(
            (np.eye(2), np.kron(np.array([[0, -19, 106, -264, 646]]) * dm, np.eye(2)))
        )
        self.matrices = []
        self.forcing = []
        for i in self.indices:
            t = e - v[i]
            L = np.array([[-kappa, 2 * r[i] + t * r[i] / c**2], [-t * r[i], kappa]])
            Lt = np.array([[0, r[i] / c**2], [-r[i], 0]])
            top = alpha * base + (alpha * dm * 251 * L + beta * np.eye(2)) @ pred
            matrix = np.zeros((12, 12))
            matrix[:2] = top
            matrix[2:10, 4:12] = np.eye(8)
            matrix[10:] = L @ top
            derivative = np.zeros((12, 12))
            derivative[:2] = alpha * dm * 251 * Lt @ pred
            derivative[10:] = Lt @ top + L @ derivative[:2]
            self.forcing.append(derivative @ state)
            self.matrices.append(matrix)
            state = matrix @ state
            self.values[i] = state[:2]
        self.matrices = np.array(self.matrices)
        self.forcing = np.array(self.forcing)

    def action(self, dv, de, initial_dot):
        state = self.initial_derivative * initial_dot
        out = np.zeros_like(self.values)
        out[self.initial_indices] = self.initial_y_derivative * initial_dot
        for i, M, f in zip(self.indices, self.matrices, self.forcing):
            state = M @ state + f * (de - dv[i])
            out[i] = state[:2]
        return out

    def pullback(self, y):
        state = np.zeros(12)
        gv = np.zeros(len(self.r))
        ge = 0.0
        for i, M, f in zip(self.indices[::-1], self.matrices[::-1], self.forcing[::-1]):
            state[:2] += y[i]
            g = f @ state
            ge += g
            gv[i] -= g
            state = M.T @ state
        gi = state @ self.initial_derivative + np.sum(
            y[self.initial_indices] * self.initial_y_derivative
        )
        return gv, ge, gi


class MatchedCoreState:
    """One fixed-occupation native core state, differentiated implicitly.

    energy and potential are Hartree; indices are zero-based, from the native
    central matching algorithm. kappa is the Dirac angular quantum number.
    NREL is c=1e30. No small-component response is discarded in RELA.
    If F=Qtilde_out/P_out-Qtilde_in/P_in at the matching point, stationarity
    of the matching equation gives dE=-F_V[dV]/F_E. Both the origin series
    and asymptotic boundary contribute to F_E. The prepared energy covector
    and energy-driven wave response are reused in every action. Normalizing
    the matched density q gives dn=occ*(dq-q*dI/I)/I, I=CHARGE(q).
    The input operator acts on the FULL extended potential; combine with
    CorePotentialExtension for native r*V inputs. Output is 4*pi*r^2*dn_core.
    """

    def __init__(
        self,
        radius,
        potential,
        energy,
        *,
        z,
        kappa,
        match,
        end,
        occupation,
        relativistic,
        output_points
    ):
        self.r = np.asarray(radius, float)
        self.v = checked(potential, self.r.shape, "core potential", real=True)
        native_charge_weights(self.r)
        if (
            any(int(i) != i for i in (match, end, output_points))
            or not 5 <= match < end - 5
            or not end < len(self.r)
            or not 0 < output_points <= len(self.r)
        ):
            raise ValueError("invalid core matching/output indices")
        if (
            not np.all(np.isfinite([energy, z, kappa, occupation]))
            or energy >= 0
            or z <= 0
            or int(kappa) != kappa
            or kappa == 0
            or occupation < 0
        ):
            raise ValueError("invalid bound core state")
        if not isinstance(relativistic, (bool, np.bool_)):
            raise ValueError("explicit relativistic flag required")
        if relativistic and z >= abs(kappa) * 137.0359895:
            raise ValueError("supercritical point-nucleus core state")
        match, end, output_points = map(int, (match, end, output_points))
        self.e = float(energy)
        self.c = 137.0359895 if relativistic else 1e30
        self.im = int(match)
        self.end = int(end)
        self.nout = int(output_points)
        self.occ = float(occupation)
        self.metric = np.array([1.0, 1 / self.c**2])
        self.out = _CoreTrajectory(
            self.r,
            self.v,
            energy,
            kappa,
            self.c,
            np.arange(5),
            np.arange(5, match + 1),
            _origin(self.r, z, kappa, self.c, energy, self.v[0]),
        )
        ii = np.arange(end, end - 5, -1)
        self.inside = _CoreTrajectory(
            self.r,
            self.v,
            energy,
            kappa,
            self.c,
            ii,
            np.arange(end - 5, match - 1, -1),
            _outer(self.r[ii], self.c, energy),
        )
        yo, yi = self.out.values, self.inside.values
        self.factor = yo[match, 0] / yi[match, 0]
        self.joined = yo.copy()
        self.joined[match:] = self.factor * yi[match:]
        self.weights = np.zeros(len(self.r))
        self.weights[: end + 1] = native_charge_weights(self.r[: end + 1])
        self.unscaled = np.sum(self.joined**2 * self.metric, axis=1)
        self.norm = self.weights @ self.unscaled
        if self.norm <= 0 or not np.isfinite(self.norm):
            raise ValueError("invalid core norm")
        self.density = self.occ * self.unscaled[:output_points] / self.norm
        self.residual = yo[match, 1] / yo[match, 0] - yi[match, 1] / yi[match, 0]
        if not np.isfinite(self.residual) or abs(self.residual) > 1e-7 * max(
            1, abs(yo[match, 1] / yo[match, 0]), abs(yi[match, 1] / yi[match, 0])
        ):
            raise ValueError("central core eigenvalue is not converged")
        self.energy_wave, self.energy_residual = self._fixed_action(
            np.zeros(len(self.r)), 1.0
        )
        gv, ge = self._fixed_pullback(np.zeros_like(self.joined), 1.0)
        if abs(ge - self.energy_residual) > 1e-10 * max(abs(ge), 1):
            raise ValueError("core energy derivative transpose failed")
        if not np.isfinite(ge) or abs(ge) < 1e-20:
            raise ValueError("singular core matching derivative")
        self.energy_gradient = -gv / ge
        self.operator = PairedOperator(
            len(self.r), self.nout, self.action, self.pullback, "matched core response"
        )

    def _fixed_action(self, dv, de):
        a = self.out.action(dv, de, dv[0] - de)
        b = self.inside.action(dv, de, de)
        i = self.im
        yo, yi = self.out.values, self.inside.values
        df = (a[i, 0] - self.factor * b[i, 0]) / yi[i, 0]
        joined = a.copy()
        joined[i:] = self.factor * b[i:] + df * yi[i:]
        residual = (a[i, 1] - yo[i, 1] / yo[i, 0] * a[i, 0]) / yo[i, 0] - (
            b[i, 1] - yi[i, 1] / yi[i, 0] * b[i, 0]
        ) / yi[i, 0]
        return joined, residual

    def _fixed_pullback(self, y, residual_bar):
        i = self.im
        yo, yi = self.out.values, self.inside.values
        ao = y.copy()
        ao[i:] = 0
        bi = np.zeros_like(y)
        bi[i:] = self.factor * y[i:]
        factor_bar = np.sum(y[i:] * yi[i:])
        ao[i, 0] += factor_bar / yi[i, 0]
        bi[i, 0] -= factor_bar * self.factor / yi[i, 0]
        ao[i] += residual_bar * np.array([-yo[i, 1] / yo[i, 0] ** 2, 1 / yo[i, 0]])
        bi[i] -= residual_bar * np.array([-yi[i, 1] / yi[i, 0] ** 2, 1 / yi[i, 0]])
        av, ae, ai = self.out.pullback(ao)
        bv, be, ei = self.inside.pullback(bi)
        av[0] += ai
        return av + bv, ae - ai + be + ei

    def action(self, dv):
        dv = checked(
            dv, self.r.shape, "extended core potential perturbation", real=True
        )
        de = self.energy_gradient @ dv
        wave, _ = self._fixed_action(dv, 0.0)
        wave += de * self.energy_wave
        dd = 2 * np.sum(self.joined * wave * self.metric, axis=1)
        return (
            self.occ
            / self.norm
            * (dd - self.unscaled * (self.weights @ dd) / self.norm)[: self.nout]
        )

    def pullback(self, y):
        y = checked(y, (self.nout,), "core density covector", real=True)
        full = np.zeros(len(self.r))
        full[: self.nout] = y
        weighted = full - self.weights * (full @ self.unscaled) / self.norm
        wave_bar = (
            2 * self.occ / self.norm * weighted[:, None] * self.joined * self.metric
        )
        gv, ge = self._fixed_pullback(wave_bar, 0.0)
        return gv + ge * self.energy_gradient


class AtomCoreResponse:
    """Sum occupied states after one shared native potential extension.

    Each state's density already includes its occupation. The output uses
    CLMCOR/CLMSUM's c00 convention, so no sqrt(4*pi) belongs in this sum.
    Normalization covers the native extended core domain; charge inside RMT
    can change through core leakage. Occupations and case.inc shift are fixed.
    """

    def __init__(self, extension, states):
        self.extension = extension
        self.states = tuple(states)
        self.n = extension.n
        if not self.states:
            raise ValueError("atom core response requires its occupied states")
        for s in self.states:
            if (
                s.nout != self.n
                or s.r.shape != extension.r.shape
                or not np.allclose(s.r, extension.r, rtol=1e-12, atol=0)
                or not np.allclose(s.v, extension.central, rtol=1e-10, atol=1e-10)
            ):
                raise ValueError("core state/extension central data mismatch")
        self.density = sum((s.density for s in self.states), start=np.zeros(self.n))
        self.operator = PairedOperator(
            self.n, self.n, self.action, self.pullback, "summed atom core"
        )

    def action(self, dvsp):
        dv = self.extension.action(dvsp)
        return sum((s.action(dv) for s in self.states), start=np.zeros(self.n))

    def pullback(self, y):
        y = checked(y, (self.n,), "atom core density covector", real=True)
        return self.extension.pullback(
            sum(
                (s.pullback(y) for s in self.states),
                start=np.zeros(len(self.extension.r)),
            )
        )


class PeriodicCoreResponse:
    """Full packed native potential -> spherical core-density response C.

    Atoms without core must be listed explicitly in empty_atoms. No core
    response is assigned to nonspherical or interstitial density coordinates,
    matching LCORE's spherical approximation. This is the unmixed C map; no
    mixer normalization, screening, or derivative of occupations is implied.
    """

    def __init__(self, layout, atoms, *, empty_atoms=()):
        self.layout = layout
        self.size = layout.size
        self.atoms = dict(atoms)
        self.empty_atoms = set(empty_atoms)
        if self.empty_atoms.intersection(self.atoms) or set(
            self.atoms
        ) | self.empty_atoms != set(range(len(layout.integrations))):
            raise ValueError("every atom must have explicit core or no-core metadata")
        lookup = dict(zip(layout.keys, layout.slices))
        self.slices = {}
        for a, op in self.atoms.items():
            sl = lookup[a, 0, 0]
            if sl.stop - sl.start != op.n or not np.allclose(
                layout.integrations[a].r, op.extension.r[: op.n], rtol=1e-12, atol=0
            ):
                raise ValueError("core output radial layout mismatch")
            self.slices[a] = sl
        self.operator = PairedOperator(
            self.size, self.size, self.action, self.pullback, "periodic spherical core"
        )

    def action(self, x):
        x = checked(x, (self.size,), "core potential field", real=True)
        out = np.zeros(self.size)
        for a, op in self.atoms.items():
            out[self.slices[a]] = op.action(x[self.slices[a]])
        return out

    def pullback(self, y):
        y = checked(y, (self.size,), "core density field covector", real=True)
        out = np.zeros(self.size)
        for a, op in self.atoms.items():
            out[self.slices[a]] = op.pullback(y[self.slices[a]])
        return out
