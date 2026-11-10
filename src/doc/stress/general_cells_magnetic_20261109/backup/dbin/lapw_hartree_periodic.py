"""Prepared fixed-geometry LAPW0 periodic Hartree derivative and transpose.

The input is the CLMSUM convention: c00=4*pi*r**2*n_spherical, other
cL=r**2*nL, and Fourier density STAR SUMS. The output is r*dV_spherical
in Ry, dV_L in Ry, and INTERNAL REAN4 interstitial per-member coefficients in Ry.
The file-coordinate adapter conjugates these for complex native VNS inputs.
Ions, RMT, cell, local frames and Fourier cutoffs stay fixed. This is the
charge-to-Hartree part of the response, not an explicit strain derivative.

The initial adapter deliberately rejects cubic-harmonic compression and
non-primitive/non-orthogonal cells. It uses the default CoulEM/CHARGE rules;
LDMTest SIMP/NCEM are not supported. XC, core and screening are separate.
"""

import numpy as np
from scipy.special import spherical_jn, sph_harm_y, factorial2
from lapw_hartree_radial import NativeCoulombIntegral, RadialHartreeDirichlet
from lapw_density_response import checked
from stress_response_adjoint import PairedOperator


def native_local_axes(axes):
    """Apply LAPW0 SymmRot's repair of rounded structure-file entries.

    For example 0.7071068 is replaced by 1/sqrt(2). A general polar
    orthogonalization would be a different native algorithm, so do not use it.
    """
    out = np.array(checked(axes, (3, 3), "local axes", real=True), copy=True)
    for value in [*(1 / np.sqrt(np.arange(2, 7))), np.sqrt(0.75), np.sqrt(2 / 3)]:
        select = abs(abs(out) - value) < 1e-6
        out[select] = np.copysign(value, out[select])
    return out


def native_charge_weights(radius):
    """Default LAPW0 CHARGE whole-grid weights, including dr=r*dlog(r).

    CHARGE uses the extended Simpson endpoint rule for >=8 points. Unlike
    CoulInt, it has no spline solve. No origin correction is added by the
    native default branch. MIXER uses the same default CHARGE rule. Keep it
    distinct from the radial Green-function integration and the legacy spline
    quadrature used by older response-driver charge-normalization tests.
    """
    grid = NativeCoulombIntegral(radius)
    n = grid.n
    if n < 8:
        # Default ConOpt=.true. short-grid branch, useful for small controls.
        w = {
            4: np.array([1, 3, 3, 1]) * 3 / 8,
            5: np.array([1, 4, 2, 4, 1]) / 3,
            6: np.array([17, 59, 44, 44, 59, 17]) / 48,
            7: np.array([35, 113, 97, 86, 97, 113, 35]) / 96,
        }[n]
    else:
        w = np.ones(n)
        w[:4] = np.array([17, 59, 43, 49]) / 48
        w[-4:] = w[:4][::-1]
    return w * grid.dx * grid.r


class NativeStarTable:
    """Distinct STERN members and mean nonsymmorphic phases, prepared once.

    rotations are the native reciprocal matrices (wstruct.timat.transpose).
    The two representations must not be confused: density/POTK store star
    sums, whereas the reanalyzed potential consumed by LAPW1 stores averages.
    """

    def __init__(self, representatives, rotations, translations, *, inversion):
        rep = np.asarray(representatives)
        rot = np.asarray(rotations)
        if (
            rep.ndim != 2
            or rep.shape[1] != 3
            or not len(rep)
            or not np.issubdtype(rep.dtype, np.integer)
        ):
            raise ValueError("invalid Fourier representatives")
        if (
            rot.ndim != 3
            or rot.shape[1:] != (3, 3)
            or not len(rot)
            or not np.issubdtype(rot.dtype, np.integer)
        ):
            raise ValueError("invalid reciprocal rotations")
        tau = checked(translations, (len(rot), 3), "translations", real=True)
        if not isinstance(inversion, (bool, np.bool_)):
            raise ValueError("inversion flag must be explicit")
        if np.any(rep[0]) or np.any(np.all(rep[1:] == 0, axis=1)):
            raise ValueError("unique G=0 must be the first star")
        self.representatives = rep.copy()
        self.inversion = bool(inversion)
        self.nstar = len(rep)
        self.stars = []
        seen = set()
        for g in rep:
            members = {}
            for rotation, shift in zip(rot, tau):
                phase = np.exp(2j * np.pi * (g @ shift))
                if inversion:
                    phase = complex(phase.real)
                members.setdefault(tuple(rotation @ g), []).append(phase)
            if seen.intersection(members):
                raise ValueError("overlapping Fourier stars")
            seen.update(members)
            self.stars.append(
                (
                    np.array(list(members), int),
                    np.array([np.mean(p) for p in members.values()]),
                )
            )
        self.sizes = np.array([len(g) for g, p in self.stars])
        self.members = np.concatenate([g for g, p in self.stars])
        self.phases = np.concatenate([p for g, p in self.stars])
        self.owners = np.repeat(np.arange(self.nstar), self.sizes)
        for a in [
            self.representatives,
            self.sizes,
            self.members,
            self.phases,
            self.owners,
        ]:
            a.setflags(write=False)

    def harmonic_average(self, reciprocal, position, local_axes, lm):
        """YKAV = mean_star[tau_G exp(iG.Ra) Y_lm*(local G)]."""
        q = self.members @ reciprocal.T @ local_axes.T
        norm = np.linalg.norm(q, axis=1)
        theta = np.arccos(
            np.clip(
                np.divide(q[:, 2], norm, out=np.ones_like(norm), where=norm > 0), -1, 1
            )
        )
        phi = np.arctan2(q[:, 1], q[:, 0])
        phase = (
            self.phases
            * np.exp(2j * np.pi * (self.members @ position))
            / self.sizes[self.owners]
        )
        result = np.empty((len(lm), self.nstar), complex)
        for i, (l, m) in enumerate(lm):
            row = phase * sph_harm_y(l, m, theta, phi).conj()
            result[i] = np.bincount(
                self.owners, row.real, minlength=self.nstar
            ) + 1j * np.bincount(self.owners, row.imag, minlength=self.nstar)
            result[i, 0] = 1 / np.sqrt(4 * np.pi) if l == 0 else 0
        return result


class NativeHartreeReanalysis:
    """REAN0/SETFFT/REAN4/GETFFT, including native finite FFT-grid aliasing.

    Input: unmasked potential star sums (POTK). Output: masked per-member
    potential coefficients (vns convention). The step function is represented
    on the FULL native FFT reciprocal box, not just the density star list.
    Taking real fields matches the native discrete extension for redundant
    complex input coordinates and makes the exact transpose real-linear.
    """

    def __init__(self, stars, shape, reciprocal, positions, radii, volume):
        self.stars = stars
        if len(shape) != 3 or any(int(n) != n or n < 4 or n % 2 for n in shape):
            raise ValueError("native reanalysis requires an explicit even FFT grid")
        self.shape = tuple(map(int, shape))
        b = checked(reciprocal, (3, 3), "reciprocal matrix", real=True)
        pos = np.asarray(positions)
        checked(pos, (len(pos), 3), "site positions", real=True)
        radii = checked(radii, (len(pos),), "site RMT", real=True)
        if np.any(radii <= 0) or not np.isfinite(volume) or volume <= 0:
            raise ValueError("invalid sphere radii/volume")
        if np.any(abs(stars.members) >= np.array(shape) // 2):
            raise ValueError("star touches FFT Nyquist plane or is outside grid")
        self.indices = tuple((stars.members % np.array(shape)).T)
        # Native REAN0 stores +iG.R phases. REAN4 conjugates them before FFT.
        axes = [np.fft.fftfreq(n) * n for n in shape]
        labels = np.stack(np.meshgrid(*axes, indexing="ij"), axis=-1)
        norm = np.linalg.norm(labels @ b.T, axis=-1)
        mask = np.zeros(self.shape, complex)
        for p, r in zip(pos, radii):
            x = norm * r
            ratio = np.divide(
                spherical_jn(1, x), x, out=np.full_like(x, 1 / 3), where=x != 0
            )
            phase = np.exp(2j * np.pi * (labels @ p))
            if stars.inversion:
                phase = phase.real
            mask -= 4 * np.pi * r**3 * ratio * phase / volume
        mask[0, 0, 0] += 1
        self.mask_fourier = mask
        self.mask_real = np.fft.fftn(mask.conj()).real
        self.weight = stars.phases / stars.sizes[stars.owners]
        self.points = int(np.prod(shape))
        self.operator = PairedOperator(
            2 * stars.nstar,
            2 * stars.nstar,
            self._forward_vector,
            self._reverse_vector,
            "Hartree FFT reanalysis",
        )
        for a in [self.mask_fourier, self.mask_real, self.weight]:
            a.setflags(write=False)

    def _collect(self, member_values):
        return np.bincount(
            self.stars.owners, member_values.real, minlength=self.stars.nstar
        ) + 1j * np.bincount(
            self.stars.owners, member_values.imag, minlength=self.stars.nstar
        )

    def action(self, potential_star_sums):
        z = checked(potential_star_sums, (self.stars.nstar,), "potential star sums")
        grid = np.zeros(self.shape, complex)
        grid[self.indices] = z[self.stars.owners] * self.weight
        physical = np.fft.fftn(grid).real
        transformed = np.fft.ifftn(physical * self.mask_real)
        return self._collect(transformed[self.indices] * self.weight.conj())

    def pullback(self, covector):
        y = checked(covector, (self.stars.nstar,), "reanalysis covector")
        grid = np.zeros(self.shape, complex)
        grid[self.indices] = y[self.stars.owners] * self.weight
        physical = np.fft.fftn(grid).real / self.points
        transformed = np.fft.ifftn(physical * self.mask_real) * self.points
        return self._collect(transformed[self.indices] * self.weight.conj())

    def _forward_vector(self, x):
        return self.action(np.ascontiguousarray(x).view(complex)).view(float)

    def _reverse_vector(self, x):
        return self.pullback(np.ascontiguousarray(x).view(complex)).view(float)


class PeriodicHartreeResponse:
    """Electron-density derivative through native pseudocharge and radial solve.

    q=-S m + D rho: replace the Fourier continuation's sphere multipoles
    by the actual radial-density multipoles. S converts native real harmonics
    to complex ones; D integrates the PW continuation analytically in spheres.
    p=8*pi*rho/G**2 + A q, p(G=0)=0: periodic Poisson/pseudocharge solution.
    b=Re(B p): boundary data, with spherical b in physical scalar Ry.
    v_MT=Dirichlet(c,b); v_I=REAN4(p).

    Each factor is prepared once. No density-coordinate probes or dense
    whole-density Jacobian are constructed. The reverse follows this exact
    chain and uses Re(vdot) for complex intermediate arrays. The fixed G=0
    convention is the native Coulomb gauge, not a neutrality correction.
    """

    def __init__(self, template, structure, fft_shape, *, ncon=9, inversion):
        s = structure
        if s.lattyp != "P" or not np.allclose(
            [s.alpha, s.beta, s.gamma], 90, atol=1e-10
        ):
            raise NotImplementedError(
                "initial Hartree adapter requires primitive orthogonal cells"
            )
        if any(i >= 0 for i in s.iatnr):
            raise NotImplementedError(
                "compressed cubic harmonics are not yet supported"
            )
        if int(ncon) != ncon or not 0 <= ncon <= 19:
            raise ValueError("invalid pseudocharge smoothing order")
        self.template = template.zeros()
        self.keys = tuple(template.rad)
        self.ncon = int(ncon)
        reps = np.array(list(template.pw), int)
        self.stars = NativeStarTable(
            reps, s.timat.transpose(0, 2, 1), s.tau, inversion=inversion
        )
        reciprocal = np.diag(2 * np.pi / np.array([s.a, s.b, s.c]))
        volume = s.a * s.b * s.c
        self.volume = volume
        norm = np.linalg.norm(reps @ reciprocal.T, axis=1)
        self.diagonal = np.zeros(len(reps))
        self.diagonal[1:] = 8 * np.pi / norm[1:] ** 2
        self.integrations = [NativeCoulombIntegral(np.exp(g)) for g in template.grids]
        if (
            len(self.integrations) != s.nat
            or len(s.pos) != sum(s.mult)
            or any(m < 1 for m in s.mult)
        ):
            raise ValueError("inconsistent structure/density atom layout")
        for a, grid in enumerate(self.integrations):
            if grid.n != s.jrj[a] or not np.allclose(
                [grid.r[0], grid.r[-1]], [s.r0[a], s.rmt[a]], rtol=1e-12, atol=0
            ):
                raise ValueError("density grid differs from structure")
        if any(not 0 <= a < s.nat for a, l, m in self.keys) or any(
            (a, 0, 0) not in self.keys for a in range(s.nat)
        ):
            raise ValueError("missing spherical harmonic or invalid atom index")
        if any(
            np.shape(template.rad[key]) != (self.integrations[key[0]].n,)
            for key in self.keys
        ):
            raise ValueError("invalid radial coefficient length")
        self.radial = [
            RadialHartreeDirichlet(self.integrations[a], abs(l))
            for a, l, m in self.keys
        ]
        self.moments = [
            native_charge_weights(self.integrations[a].r)
            * self.integrations[a].r ** abs(l)
            / (np.sqrt(4 * np.pi) if l == 0 else 1)
            for a, l, m in self.keys
        ]
        self.slices = []
        offset = 0
        for key in self.keys:
            n = len(template.rad[key])
            self.slices.append(slice(offset, offset + n))
            offset += n
        self.pw_offset = offset
        self.size = offset + 2 * len(reps)
        # Same radial ordering as the Density template; spherical output differs
        # in units (r V) and PW differs in normalization, documented above.
        acomp, dcomp, scomp = [], [], []
        boundary = np.zeros((len(self.keys), len(reps)), complex)
        first = np.r_[0, np.cumsum(s.mult)[:-1]]
        for atom in range(s.nat):
            keys = [k for k, key in enumerate(self.keys) if key[0] == atom]
            lm = sorted(
                set(
                    (abs(self.keys[k][1]), m)
                    for k in keys
                    for m in ({self.keys[k][2], -self.keys[k][2]})
                )
            )
            if not lm or any(l > ncon for l, m in lm):
                raise ValueError("missing atom harmonics or l exceeds smoothing order")
            y = self.stars.harmonic_average(
                reciprocal, s.pos[first[atom]], native_local_axes(s.rotloc[atom]), lm
            )
            lm_index = {key: i for i, key in enumerate(lm)}
            conversion = np.zeros((len(lm), len(self.keys)), complex)
            for k in keys:
                a, signed_l, m = self.keys[k]
                l = abs(signed_l)
                if m < 0 or m > l or (m == 0 and signed_l < 0):
                    raise ValueError("invalid native real-harmonic label")
                phase = (-1j if signed_l < 0 else 1) * (-1) ** m
                conversion[lm_index[l, m], k] = phase / np.sqrt(2) if m else 1
                if m:
                    conversion[lm_index[l, -m], k] = (
                        1j if signed_l < 0 else 1
                    ) / np.sqrt(2)
                factor = 4 * np.pi * 1j**l
                if m:
                    factor *= np.sqrt(2) / phase
                if l == 0:
                    factor /= np.sqrt(4 * np.pi)
                boundary[k] = (
                    factor * spherical_jn(l, norm * s.rmt[atom]) * y[lm_index[l, m]]
                )
            correction = np.zeros((len(lm), len(reps)), complex)
            pseudo = np.zeros((len(reps), len(lm)), complex)
            radius = s.rmt[atom]
            x = norm[1:] * radius
            for i, (l, m) in enumerate(lm):
                correction[i, 1:] = (
                    4
                    * np.pi
                    * radius ** (l + 2)
                    * 1j**l
                    * y[i, 1:]
                    * spherical_jn(l + 1, x)
                    / norm[1:]
                )
                if l == 0:
                    correction[i, 0] = 4 * np.pi * radius**3 / (3 * np.sqrt(4 * np.pi))
                mp = factorial2(2 * ncon + 3) / factorial2(2 * l + 1)
                pseudo[1:, i] = (
                    -self.diagonal[1:]
                    * (4 * np.pi / volume)
                    * self.stars.sizes[1:]
                    * s.mult[atom]
                    * (-1j / radius) ** l
                    * mp
                    * spherical_jn(ncon + 1, x)
                    / x ** (ncon - l + 1)
                    * y[i, 1:].conj()
                )
            acomp.append(pseudo)
            dcomp.append(correction)
            scomp.append(conversion)
        self.A = np.concatenate(acomp, axis=1)
        self.D = np.concatenate(dcomp, axis=0)
        self.S = np.concatenate(scomp, axis=0)
        self.B = boundary
        site_radii = np.repeat(s.rmt, s.mult)
        self.reanalysis = NativeHartreeReanalysis(
            self.stars, fft_shape, reciprocal, s.pos, site_radii, volume
        )
        self.operator = PairedOperator(
            self.size, self.size, self.action, self.pullback, "periodic Hartree"
        )
        for a in [self.A, self.D, self.S, self.B, self.diagonal, *self.moments]:
            a.setflags(write=False)

    def _split(self, x):
        x = checked(x, (self.size,), "Hartree field", real=True)
        return [x[sl] for sl in self.slices], np.ascontiguousarray(
            x[self.pw_offset :]
        ).view(complex)

    def intermediates(self, x):
        radial, rho = self._split(x)
        moment = np.array([w @ c for w, c in zip(self.moments, radial)])
        q = -self.S @ moment + self.D @ rho
        potential = self.diagonal * rho + self.A @ q
        b = (self.B @ potential).real
        return radial, potential, b

    def action(self, x):
        radial, potential, boundary = self.intermediates(x)
        mt = [op.action(c, b) for op, c, b in zip(self.radial, radial, boundary)]
        ir = self.reanalysis.action(potential)
        return np.concatenate(mt + [ir.view(float)])

    def pullback(self, y):
        radial_y, ir_y = self._split(y)
        bars = [op.pullback(z) for op, z in zip(self.radial, radial_y)]
        boundary_y = np.array([b for c, b in bars])
        p_y = self.B.conj().T @ boundary_y + self.reanalysis.pullback(ir_y)
        q_y = self.A.conj().T @ p_y
        moment_y = -(self.S.conj().T @ q_y).real
        radial = [c + w * b for (c, _), w, b in zip(bars, self.moments, moment_y)]
        rho_y = self.diagonal * p_y + self.D.conj().T @ q_y
        return np.concatenate(radial + [rho_y.view(float)])

    def native_potential_output(
        self,
        spherical_atoms,
        nonspherical_harmonics,
        representatives,
        *,
        complex_interstitial
    ):
        """Explicit projection/reordering to NativePotentialBands input layout.

        The selected LM/G cutoffs are supplied by the actual potential reader.
        This matters if charge contains more harmonics than the Hamiltonian.
        Missing requested channels fail; no absent potential is silently zeroed.
        The packed input contains INTERNAL REAN4 Fourier coefficients. LAPW0
        writes conjg(cvalue) to a complex VNS file (real part only for a real
        file). NativePotentialBands accepts that FILE convention. Conjugation
        belongs at this adapter boundary exactly once, in both directions.
        Returns (paired selector, NativePotentialCoordinates). Compose selector
        @ self.operator for the physical charge-to-band-input Hartree map.
        """
        from lapw_potential_response import NativePotentialCoordinates

        spherical_atoms = tuple(spherical_atoms)
        if len(set(spherical_atoms)) != len(spherical_atoms) or set(
            spherical_atoms
        ) != set(range(len(self.integrations))):
            raise ValueError("spherical atoms must include every atom exactly once")
        lookup = {key: sl for key, sl in zip(self.keys, self.slices)}
        rows = []
        for a in spherical_atoms:
            sl = lookup[a, 0, 0]
            rows.extend(range(sl.start, sl.stop))
        spherical_size = len(rows)
        shapes = {}
        for a, harmonics in sorted(nonspherical_harmonics.items()):
            harmonics = tuple(map(tuple, harmonics))
            if (
                a not in spherical_atoms
                or not harmonics
                or len(set(harmonics)) != len(harmonics)
                or (0, 0) in harmonics
            ):
                raise ValueError("invalid nonspherical output harmonic list")
            if any((a, l, m) not in lookup for l, m in harmonics):
                raise ValueError(
                    "requested potential harmonic absent in density layout"
                )
            columns = [
                np.arange(lookup[a, l, m].start, lookup[a, l, m].stop)
                for l, m in harmonics
            ]
            shapes[a] = (len(columns[0]), len(columns))
            rows.extend(np.array(columns).T.ravel())
        rep = np.asarray(representatives)
        if (
            rep.ndim != 2
            or rep.shape[1] != 3
            or not len(rep)
            or not np.issubdtype(rep.dtype, np.integer)
        ):
            raise ValueError("invalid output potential representatives")
        lookup_pw = {tuple(g): i for i, g in enumerate(self.stars.representatives)}
        if len(set(map(tuple, rep))) != len(rep) or any(
            tuple(g) not in lookup_pw for g in rep
        ):
            raise ValueError("duplicate or missing output potential representative")
        indices = np.array([lookup_pw[tuple(g)] for g in rep])
        rows.extend(self.pw_offset + 2 * indices)
        if complex_interstitial:
            rows.extend(self.pw_offset + 2 * indices + 1)
        rows = np.array(rows, int)
        signs = np.ones(len(rows))
        if complex_interstitial:
            signs[-len(rep) :] = -1.0
        coordinates = NativePotentialCoordinates(
            spherical_size, shapes, len(rep), complex_interstitial=complex_interstitial
        )

        def forward(x):
            return (
                signs * checked(x, (self.size,), "Hartree potential", real=True)[rows]
            )

        def reverse(y):
            y = checked(y, (coordinates.size,), "potential covector", real=True)
            result = np.zeros(self.size)
            np.add.at(result, rows, signs * y)
            return result

        return (
            PairedOperator(
                self.size,
                coordinates.size,
                forward,
                reverse,
                "Hartree to native potential coordinates",
            ),
            coordinates,
        )
