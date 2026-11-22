"""Central derivatives of LAPWSO's nonmagnetic second-variational Hamiltonian.

Radial directions use r*V in Ry and native linearization energies in Ry.
The spin order of all block matrices here is LAPWSO's (down, up). The
DMFT density dump labels its vector components in the opposite order.
"""

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix

from lapw_augmentation_response import PreparedAugmentationChannel, harmonics
from lapw_radial_native import NativeRadialContext


def _radial_gradient(radius):
    """Five-node Lagrange derivative at each native VDERIV radius."""
    n = len(radius)
    offsets = np.clip(np.arange(n) - 2, 0, n - 5)
    columns = offsets[:, None] + np.arange(5)
    # Dimensionless offsets avoid poorly scaled Vandermonde rows near r=0.
    scale = radius[columns[:, -1]] - radius[columns[:, 0]]
    x = (radius[columns] - radius[:, None]) / scale[:, None]
    powers = x[:, None, :] ** np.arange(5)[None, :, None]
    rhs = np.zeros((n, 5, 1))
    rhs[:, 1, 0] = 1 / scale
    weights = np.linalg.solve(powers, rhs)[:, :, 0]
    return csr_matrix(
        (weights.ravel(), columns.ravel(), np.arange(0, 5 * n + 1, 5)), shape=(n, n)
    )


class SOCRadialChannel:
    """Radial matrix and transpose of its discrete derivative for one site/l.

    GARADME uses only the large components in its SOC integral. The small
    components still enter their normalization and the matching coefficients.
    Filtering the first APW+lo primitive is essential: LAPWSO packs genuine
    radial wavefunctions consecutively, whereas the scalar API reserves an
    empty radial slot for that local combination of u and udot.
    """

    def __init__(self, context, potential):
        if context.quadrature != "lapwso" or not context.relativistic:
            raise ValueError("SOC requires native LAPWSO radial functions")
        self.ctx = context
        self.slots = np.arange(2 + context.nlo)
        if not context.lapw and context.nlo:
            self.slots = np.delete(self.slots, 2)
        radius = context.radius
        self.potential = np.asarray(potential, dtype=float) / radius
        self.gradient = _radial_gradient(radius)
        self.twice_c = 2 * 137.0359895
        energies = np.r_[context.energy, context.energy, context.lo_energies][
            self.slots
        ]
        self.denominator = (
            self.twice_c + (energies - self.potential[:, None]) / self.twice_c
        )
        self.waves = context.primitive[:, 0, self.slots] / self.denominator
        # Native CALI sums r*f(r) in log(r), with its own even-grid start.
        # The r here cancels the explicit 1/r in GARADME.
        n = len(radius)
        start = n % 2
        w = np.zeros(n)
        w[start : n - 1 : 2] += 4
        w[start + 1 : n : 2] += 2
        w[-1] -= 1
        if start:
            w[0] += 1
        self.weights = w * (2 * np.log(radius[-1] / radius[0]) / (3 * (n - 1)))
        self.field = self.weights * (self.gradient @ self.potential)
        self.matrix = self.waves.T @ (self.field[:, None] * self.waves)

    def action(self, dv, de=0.0, delo=None):
        """Return dR alongside the boundary and ordinary LO matching tangents."""
        ctx = self.ctx
        if delo is None:
            delo = np.zeros(ctx.nlo)
        du, dabc = ctx.action(dv, de, delo)
        dpotential = dv / ctx.radius
        dden = (np.r_[de, de, delo][self.slots] - dpotential[:, None]) / self.twice_c
        dw = (du[:, 0, self.slots] - self.waves * dden) / self.denominator
        leg = dw.T @ (self.field[:, None] * self.waves)
        dfield = self.weights * (self.gradient @ dpotential)
        dr = leg + leg.T + self.waves.T @ (dfield[:, None] * self.waves)
        return dr, du[-1], dabc

    def pullback(self, matrix_bar, boundary_bar, abc_bar):
        """Transpose to the potential and all radial energy coordinates."""
        ctx = self.ctx
        bar = np.asarray(matrix_bar).real
        wave_bar = self.field[:, None] * (self.waves @ (bar + bar.T))
        du_bar = np.zeros_like(ctx.primitive)
        du_bar[-1] = boundary_bar
        du_bar[:, 0, self.slots] += wave_bar / self.denominator
        denominator_bar = -wave_bar * self.waves / (self.twice_c * self.denominator)
        energy_bar = np.zeros(2 + ctx.nlo)
        energy_bar[self.slots] = denominator_bar.sum(axis=0)
        field_bar = ((self.waves @ bar) * self.waves).sum(axis=1)
        potential_bar = self.gradient.T @ (
            self.weights * field_bar
        ) - denominator_bar.sum(axis=1)
        dv, de, delo = ctx.pullback(du_bar, abc_bar)
        return (
            dv + potential_bar / ctx.radius,
            de + energy_bar[:2].sum(),
            delo + energy_bar[2:],
        )


@dataclass
class SOCSite:
    atom: int
    local_momenta: np.ndarray
    location: np.ndarray
    angular: dict


class SOCGeometry:
    """Native angular-spin matrices, positions, central radial data and LO layout."""

    def __init__(self, path, structure, potentials, stack, library=None):
        self.channels = {}
        self.sites = []
        self.lo_indices = {}
        self.nlo = 0
        self.matrix_error = 0.0
        with Path(path).open("rb") as stream:

            def array(dtype, shape):
                count = int(np.prod(shape))
                values = np.fromfile(stream, dtype=dtype, count=count)
                if len(values) != count:
                    raise ValueError("truncated LAPWSO geometry")
                return values.reshape(shape, order="F")

            magic = stream.read(8)
            if magic not in (b"SOCGEOM2", b"SOCGEOM3"):
                raise ValueError("invalid LAPWSO response geometry")
            nat, nsites, self.lmax, self.nlo = array("i4", (4,))
            policy = array("i4", (1,))[0] if magic == b"SOCGEOM3" else 1
            if policy not in (0, 1):
                raise ValueError("invalid LAPWSO local-orbital matching policy")
            self.legacy_lo_filter = bool(policy)
            if nat != structure.nat or nsites != sum(structure.mult):
                raise ValueError("LAPWSO geometry and structure disagree")
            self.volume = array("f8", (1,))[0]
            self.br1 = array("f8", (3, 3))
            for atom, v in enumerate(potentials):
                mult, enabled, offset = array("i4", (3,))
                if mult != structure.mult[atom]:
                    raise ValueError("LAPWSO site multiplicity mismatch")
                for ell in range(self.lmax + 1):
                    nr, nlo, lapw = array("i4", (3,))
                    e = array("f8", (1,))[0]
                    elo = array("f8", (nlo,))
                    native = array("f8", (nr, nr))
                    groups = []
                    for j in range(nlo):
                        size = mult * (2 * ell + 1)
                        groups.append(np.arange(offset, offset + size))
                        offset += size
                    self.lo_indices[atom, ell] = groups
                    if not ell or not enabled:
                        continue
                    r0, rmt = structure.r0[atom], structure.rmt[atom]
                    ctx = stack.enter_context(
                        NativeRadialContext(
                            v,
                            r0=r0,
                            rmt=rmt,
                            dx=np.log(rmt / r0) / (len(v) - 1),
                            energy=e,
                            lo_energies=elo,
                            ell=ell,
                            charge=structure.Znuc[atom],
                            relativistic=True,
                            lapw=bool(lapw),
                            quadrature="lapwso",
                            library=library,
                        )
                    )
                    channel = SOCRadialChannel(ctx, v)
                    self.matrix_error = max(
                        self.matrix_error, float(np.max(abs(channel.matrix - native)))
                    )
                    if not np.allclose(channel.matrix, native, rtol=3e-8, atol=2e-10):
                        raise ValueError(
                            f"LAPWSO radial reconstruction failed at atom {atom+1}, l={ell}: {self.matrix_error:g}"
                        )
                    self.channels[atom, ell] = channel
                for _ in range(mult):
                    rot, loc = array("f8", (3, 3)), array("f8", (3,))
                    angular = {
                        l: array("c16", (2 * l + 1, 2 * l + 1, 2, 2))
                        for l in range(1, self.lmax + 1)
                    }
                    self.sites.append(
                        SOCSite(atom, rot, loc, angular if enabled else {})
                    )
            if stream.read(1):
                raise ValueError("unexpected data at end of LAPWSO geometry")


class SOCMatching:
    """LAPWSO matching in the convention that produced the parent spinors."""

    @staticmethod
    def component_projection(z, masks):
        # Old spectra used a real-linear projection onto active components.
        # Its real Frobenius transpose is the same projection. New spectra
        # retain the entire complex coefficient, including tangents at zero.
        if masks is None:
            return z
        re, im = masks
        return z.real * re + 1j * z.imag * im

    def __init__(self, geometry, site, frac, coefficients, ell):
        radial = geometry.channels[site.atom, ell]
        self.ctx = radial.ctx
        self.slots = radial.slots
        self.c = coefficients
        self.q = frac @ geometry.br1.T
        npw = len(frac) - geometry.nlo
        phase = np.exp(2j * np.pi * (frac @ site.location))
        local = frac @ site.local_momenta.T
        ch = PreparedAugmentationChannel(
            ell,
            self.ctx.rmt,
            geometry.volume,
            local,
            np.linalg.norm(self.q, axis=1),
            phase,
            coefficients,
            npw,
            [npw + g for g in geometry.lo_indices[site.atom, ell]],
        )
        # Store the linear maps multiplying the plane-wave coefficient tangents.
        # The local harmonic and k+G phase are the native LAPWSO conventions.
        from scipy.special import spherical_jn

        angular = (
            4
            * np.pi
            * (1j**ell)
            / np.sqrt(geometry.volume)
            * harmonics(ell, local).conj()
        )
        q = np.linalg.norm(self.q[:npw], axis=1)
        self.f = angular[:, :npw] * (phase[:npw] * spherical_jn(ell, q * self.ctx.rmt))
        self.df = angular[:, :npw] * (
            phase[:npw] * q * spherical_jn(ell, q * self.ctx.rmt, derivative=True)
        )
        self.npw = npw
        self.locals = []
        for j, group in enumerate(geometry.lo_indices[site.atom, ell]):
            rows = npw + group
            z = phase[rows, None] * coefficients[rows]
            masks = (
                (abs(z.real) > 1e-6, abs(z.imag) > 1e-6)
                if geometry.legacy_lo_filter
                else None
            )
            ch.local[j] = angular[:, rows] @ self.component_projection(z, masks)
            self.locals.append((rows, angular[:, rows], phase[rows], masks))
        self.channel = ch
        self.z = ch.value(
            self.ctx.primitive[-1], self.ctx.lo_coefficients, lapw=self.ctx.lapw
        )[:, self.slots, :]

    def action(self, dc, velocity, boundary_dot, abc_dot):
        from copy import copy

        ctx, ch = self.ctx, self.channel
        dz = ch.action(ctx.primitive[-1], boundary_dot, abc_dot, lapw=ctx.lapw)
        if dc is not None or velocity is not None:
            # Differentiate the phase before applying the parent LO convention.
            x = dc if velocity is None else 1j * (self.q @ velocity)[:, None] * self.c
            if dc is not None and velocity is not None:
                x = x + dc
            dch = copy(ch)
            dch.fj = self.f @ x[: self.npw]
            dch.dfj = self.df @ x[: self.npw]
            dch.local = np.empty_like(ch.local)
            for j, (rows, angular, phase, masks) in enumerate(self.locals):
                v = phase[:, None] * x[rows]
                dch.local[j] = angular @ self.component_projection(v, masks)
            dz += dch.value(ctx.primitive[-1], ctx.lo_coefficients, lapw=ctx.lapw)
        return dz[:, self.slots, :]

    def pullback(self, zbar):
        """Real Frobenius transpose to radial matching, scalar vectors and position."""
        ctx, ch = self.ctx, self.channel
        full = np.zeros((len(zbar), len(ctx.lo_energies) + 2, zbar.shape[-1]), complex)
        full[:, self.slots, :] = zbar
        gb, gabc = ch.pullback(ctx.primitive[-1], full, lapw=ctx.lapw)
        b = ctx.primitive[-1]
        r = ctx.rmt
        if ctx.lapw:
            g_f = r * (-b[2, 1] * full[:, 0] + b[2, 0] * full[:, 1])
            g_df = r * (b[0, 1] * full[:, 0] - b[0, 0] * full[:, 1])
            gc_pw = self.f.conj().T @ g_f + self.df.conj().T @ g_df
        else:
            gc_pw = self.f.conj().T @ (r / b[0, 0] * full[:, 0])
        gc = np.zeros_like(self.c, dtype=complex)
        gc[: self.npw] = gc_pw
        for j, (rows, angular, phase, masks) in enumerate(self.locals):
            g_loc = r * r * np.einsum("i,min->mn", ctx.lo_coefficients[:, j], full)
            v = angular.conj().T @ g_loc
            gc[rows] += phase[:, None].conj() * self.component_projection(v, masks)
        gu = self.q.T @ (gc.conj() * (1j * self.c)).real.sum(axis=1)
        return gb, gabc, gc, gu


class SOCBandHamiltonian:
    """Second-variational SOC addition in the central scalar-band frame at one k.

    Input coefficient directions must use the same ordered LAPW rows as the
    central scalar vector. The scalar-Hamiltonian variation is added by the
    caller. Keeping the radial and scalar-subspace paths separate here makes
    the same actions usable for both an atomic direction and an input charge.
    """

    def __init__(self, geometry, k, labels, coefficients):
        self.geometry = geometry
        self.c = np.asarray(coefficients)
        self.blocks = []
        frac = np.asarray(labels).T + np.asarray(k)
        for site_idx, site in enumerate(geometry.sites):
            for ell, angular in site.angular.items():
                key = (site.atom, ell)
                match = SOCMatching(geometry, site, frac, self.c, ell)
                self.blocks.append((site_idx, key, angular, match))

    @staticmethod
    def _apply(angular, radial, z):
        # The radial integral is symmetric. <m,s|L.S|m',s'> has native
        # ordering angular[m,m',s,s']; z[m',radial,band] is a ket.
        rz = (radial @ z).reshape(len(z), -1)
        return (angular @ rz).reshape(z.shape)

    def action(self, radial_dots=None, dc=None, velocities=None, *, central=False):
        """SOC matrix derivative, or the central SOC matrix when central=True."""
        n = self.c.shape[1]
        h = np.zeros((2 * n, 2 * n), complex)
        for site, key, angular, match in self.blocks:
            channel = self.geometry.channels[key]
            z = match.z
            if not central:
                if radial_dots is None:
                    dr = np.zeros_like(channel.matrix)
                    db = np.zeros_like(channel.ctx.primitive[-1])
                    da = np.zeros_like(channel.ctx.lo_coefficients)
                else:
                    dr, db, da = radial_dots[key]
                dz = match.action(
                    dc, None if velocities is None else velocities[site], db, da
                )
            for i, j in ((0, 0), (0, 1), (1, 1)):
                a = angular[:, :, i, j]
                rz = self._apply(a, channel.matrix, z).reshape(-1, n)
                if central:
                    block = z.reshape(-1, n).conj().T @ rz
                else:
                    rdz = self._apply(a, channel.matrix, dz) + self._apply(a, dr, z)
                    block = dz.reshape(-1, n).conj().T @ rz + z.reshape(
                        -1, n
                    ).conj().T @ rdz.reshape(-1, n)
                h[i * n : (i + 1) * n, j * n : (j + 1) * n] += block
                if i != j:
                    h[j * n : (j + 1) * n, i * n : (i + 1) * n] += block.conj().T
        return h

    def pullback(self, covector):
        """Transpose of action. Sum small radial covectors over k before OUTWIN."""
        n = self.c.shape[1]
        g = np.asarray(covector)
        if g.shape != (2 * n, 2 * n) or not np.allclose(g, g.conj().T):
            raise ValueError(
                "expected Hermitian covector in the central scalar-spin frame"
            )
        radial = {}
        gc = np.zeros_like(self.c, dtype=complex)
        gu = np.zeros((len(self.geometry.sites), 3))
        for site, key, angular, match in self.blocks:
            channel = self.geometry.channels[key]
            z = match.z
            gz = np.zeros_like(z)
            gr = np.zeros_like(channel.matrix)
            for i in range(2):
                for j in range(2):
                    # Re Tr(g^H Z^H R Z). Both matching legs give the same
                    # real variation after summing the Hermitian spin blocks.
                    a = angular[:, :, i, j]
                    block = g[j * n : (j + 1) * n, i * n : (i + 1) * n]
                    v = (z.reshape(-1, n) @ block).reshape(z.shape)
                    gz += 2 * self._apply(a, channel.matrix, v)
                    av = (a @ v.reshape(len(z), -1)).reshape(z.shape)
                    gr += np.einsum("min,mjn->ij", z.conj(), av).real
            gb, ga, dc, du = match.pullback(gz)
            if key not in radial:
                radial[key] = (gr, gb, ga)
            else:
                for total, add in zip(radial[key], (gr, gb, ga)):
                    total += add
            gc += dc
            gu[site] += du
        return radial, gc, gu
