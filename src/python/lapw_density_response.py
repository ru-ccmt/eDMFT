"""Spatial valence density and exact transposes in a fixed LAPW geometry.

The band density includes its spin/k weight already. No occupation clipping,
extra k weight, or core/charge renormalization is applied here. Muffin-tin
coefficients are r**2 times the harmonic density, not the density itself.
Changing radial functions and matching coefficients is independent of changing
the band density matrix: all three product-rule terms are required.
"""

import numpy as np
from lapw_augmentation_response import RadialAngularLayout


def herm(a):
    return (a + a.conj().T) / 2


def checked(a, shape, name, real=False, hermitian=False):
    a = np.asarray(a)
    if a.shape != shape or not np.all(np.isfinite(a)):
        raise ValueError("invalid " + name)
    if real and np.iscomplexobj(a):
        raise ValueError(name + " must be real")
    if hermitian and np.linalg.norm(a - a.conj().T) > 1e-11 * max(
        1.0, np.linalg.norm(a)
    ):
        raise ValueError(name + " must be Hermitian")
    return a


class MuffinTinDensity:
    """Prepared one-site MT density; sum equivalent sites with 1/multiplicity.

    Z[component,band] uses the RadialAngularLayout component ordering.
    primitive[l][radius,large/small,radial] contains the actual stored OUTWIN
    components (r times the physical radial functions). small_weight is the
    native *density* small-component factor; it need not equal a radial
    integration routine's nonrelativistic placeholder.

    Output is [harmonic,radius]. clmsum=True includes sqrt(4*pi) on ONLY the
    spherical harmonic, matching wien_response_density.Density coordinates.
    Noncubic real harmonics follow the existing angular adapter; callers may
    supply native Gaunt matrices for other harmonic conventions.
    """

    def __init__(
        self,
        z,
        rho,
        primitive,
        harmonics,
        small_weight,
        *,
        clmsum=True,
        gaunt_matrices=None
    ):
        self.phi = {l: np.array(v, float, copy=True) for l, v in primitive.items()}
        if not self.phi or any(np.iscomplexobj(v) for v in primitive.values()):
            raise ValueError("radial primitives must be real")
        counts = {l: v.shape[2] for l, v in self.phi.items() if v.ndim == 3}
        if len(counts) != len(self.phi):
            raise ValueError("invalid radial primitive dimensions")
        self.layout = RadialAngularLayout(
            counts, harmonics, gaunt_matrices=gaunt_matrices
        )
        self.nr = next(iter(self.phi.values())).shape[0]
        for l, v in self.phi.items():
            checked(v, (self.nr, 2, counts[l]), "radial primitive", real=True)
        if not np.isfinite(small_weight) or small_weight < 0:
            raise ValueError("invalid small-component weight")
        self.weights = np.array([1.0, small_weight])
        z = np.asarray(z)
        if z.ndim != 2 or z.shape[1] < 1:
            raise ValueError("invalid matching coefficients")
        self.z = np.array(
            checked(z, (self.layout.size, z.shape[1]), "matching coefficients"), complex
        )
        self.rho = np.array(
            checked(rho, (z.shape[1], z.shape[1]), "band density", hermitian=True),
            complex,
        )
        self.scale = np.array(
            [
                np.sqrt(4 * np.pi) if clmsum and h == (0, 0) else 1.0
                for h in self.layout.harmonics
            ]
        )
        self.zrho = self.z @ self.rho
        self.blocks = self._angular(self.zrho @ self.z.conj().T)

    def _angular(self, matrix):
        out = {}
        for (l, lp), g in self.layout.gaunt.items():
            # Exact selection rules leave numerical quadrature roundoff in
            # forbidden blocks. Avoid radial work on those machine-zero blocks.
            if np.max(abs(g), initial=0) < 5e-14:
                continue
            nl, nlp = self.layout.counts[l], self.layout.counts[lp]
            block = matrix[self.layout.slices[l], self.layout.slices[lp]].reshape(
                2 * l + 1, nl, 2 * lp + 1, nlp
            )
            # Tr(G D), not a transpose of the complex sine Gaunt matrix.
            out[l, lp] = np.einsum("qmn,minj->qij", g.conj(), block).real
        return out

    def _radial(self, blocks, left, right):
        out = np.zeros((len(self.scale), self.nr))
        for (l, lp), a in blocks.items():
            out += np.einsum(
                "qij,rci,rcj,c->qr", a, left[l], right[lp], self.weights, optimize=True
            )
        return self.scale[:, None] * out

    def value(self):
        return self._radial(self.blocks, self.phi, self.phi)

    def action(self, rho_dot, z_dot=None, primitive_dot=None):
        drho = checked(rho_dot, self.rho.shape, "density derivative", hermitian=True)
        dd = self.z @ drho @ self.z.conj().T
        if z_dot is not None:
            dz = checked(z_dot, self.z.shape, "matching derivative")
            a = dz @ self.zrho.conj().T
            dd += a + a.conj().T
        return self.component_action(dd, primitive_dot)

    def component_action(self, matrix_dot, primitive_dot=None):
        """Act on the already accumulated primitive angular density matrix."""
        dd = checked(
            matrix_dot,
            (self.layout.size,) * 2,
            "component density derivative",
            hermitian=True,
        )
        result = self._radial(self._angular(dd), self.phi, self.phi)
        if primitive_dot is not None:
            if set(primitive_dot) != set(self.phi):
                raise ValueError("missing or extra radial derivative")
            dp = {
                l: checked(v, self.phi[l].shape, "radial derivative", real=True)
                for l, v in primitive_dot.items()
            }
            result += self._radial(self.blocks, dp, self.phi)
            result += self._radial(self.blocks, self.phi, dp)
        return result

    def pullback(self, covector):
        """Return (rho_bar, Z_bar, primitive_bar) in real Frobenius pairing."""
        dbar, pbar = self.component_pullback(covector)
        return herm(self.z.conj().T @ dbar @ self.z), 2 * dbar @ self.zrho, pbar

    def component_pullback(self, covector, *, include_radial=True):
        """Transpose after the k/site sum; radial reverse is needed only once."""
        y = checked(
            covector, (len(self.scale), self.nr), "MT density covector", real=True
        )
        y = y * self.scale[:, None]
        dbar = np.zeros((self.layout.size, self.layout.size), complex)
        pbar = (
            {l: np.zeros_like(v) for l, v in self.phi.items()}
            if include_radial
            else None
        )
        for (l, lp), a in self.blocks.items():
            left, right = self.phi[l], self.phi[lp]
            ab = np.einsum(
                "qr,rci,rcj,c->qij", y, left, right, self.weights, optimize=True
            )
            g = self.layout.gaunt[l, lp]
            block = np.einsum("qmn,qij->minj", g, ab)
            dbar[self.layout.slices[l], self.layout.slices[lp]] = block.reshape(
                (2 * l + 1) * left.shape[2], (2 * lp + 1) * right.shape[2]
            )
            if include_radial:
                pbar[l] += np.einsum(
                    "qr,qij,rcj,c->rci", y, a, right, self.weights, optimize=True
                )
                pbar[lp] += np.einsum(
                    "qr,qij,rci,c->rcj", y, a, left, self.weights, optimize=True
                )
        dbar = herm(dbar)
        return dbar, pbar


class MuffinTinDensitySum:
    """Sum k points/equivalent sites BEFORE the expensive radial contractions.

    Each entry is (k_index, Z, site_average). rho[k] already includes its
    spin/k weight; site_average is normally 1/MULT and must not reweight k.
    Every site uses the same representative atom's local harmonic frame.
    Only one radial data set is retained for the atom, not one per k point.
    """

    def __init__(self, entries, rho, primitive, harmonics, small_weight, **options):
        rho = np.asarray(rho)
        if rho.ndim != 3 or rho.shape[1] != rho.shape[2] or rho.shape[0] < 1:
            raise ValueError("invalid k-resolved density")
        self.rho = np.array(
            [checked(r, r.shape, "band density", hermitian=True) for r in rho], complex
        )
        self.entries = []
        self.zrho = []
        size = sum((2 * l + 1) * v.shape[2] for l, v in primitive.items())
        total = np.zeros((size, size), complex)
        for k, z, weight in entries:
            if (
                int(k) != k
                or not 0 <= k < len(rho)
                or not np.isfinite(weight)
                or weight <= 0
            ):
                raise ValueError("invalid site contribution")
            z = np.array(
                checked(z, (size, rho.shape[1]), "site matching coefficients"), complex
            )
            zr = z @ rho[k]
            total += weight * zr @ z.conj().T
            self.entries.append((int(k), z, float(weight)))
            self.zrho.append(zr)
        if not self.entries:
            raise ValueError("empty atom density sum")
        # In the component representation, the matching matrix is identity
        # and the density is the accumulated Z rho Z^dagger. Actions below
        # use component interfaces directly, without multiplying that identity.
        self.kernel = MuffinTinDensity(
            np.eye(size), herm(total), primitive, harmonics, small_weight, **options
        )

    def value(self):
        return self.kernel.value()

    def action(self, rho_dot, z_dots=None, primitive_dot=None):
        drho = (
            None
            if rho_dot is None
            else checked(rho_dot, self.rho.shape, "k-resolved density derivative")
        )
        if drho is not None:
            for d in drho:
                checked(d, d.shape, "band density derivative", hermitian=True)
        if z_dots is not None and len(z_dots) != len(self.entries):
            raise ValueError("missing site matching derivatives")
        dd = np.zeros((self.kernel.layout.size,) * 2, complex)
        for i, ((k, z, w), zr) in enumerate(zip(self.entries, self.zrho)):
            if drho is not None:
                dd += w * z @ drho[k] @ z.conj().T
            if z_dots is not None:
                dz = checked(z_dots[i], z.shape, "site matching derivative")
                a = dz @ zr.conj().T
                dd += w * (a + a.conj().T)
        return self.kernel.component_action(herm(dd), primitive_dot)

    def pullback(self, covector, *, include_basis=True, include_density=True):
        dbar, pbar = self.kernel.component_pullback(
            covector, include_radial=include_basis
        )
        rbar = np.zeros_like(self.rho) if include_density else None
        zbar = []
        for (k, z, w), zr in zip(self.entries, self.zrho):
            if include_density:
                rbar[k] += w * herm(z.conj().T @ dbar @ z)
            if include_basis:
                zbar.append(2 * w * dbar @ zr)
        return rbar, zbar, pbar


class FourierStars:
    """Native STERN star SUM and transpose, including nonsymmorphic phases.

    rotations are the native reciprocal action IZ (transpose of wstruct.timat).
    The sum is over distinct star members; duplicate operations contribute the
    average translation phase for that member. There is NO division by star size.
    Volume normalization is supplied by the density map, not repeated here.
    """

    def __init__(
        self, shape, representatives, rotations, translations, *, complex_phases=True
    ):
        self.shape = tuple(shape)
        if len(self.shape) != 3 or any(int(n) != n or n < 1 for n in self.shape):
            raise ValueError("invalid FFT dimensions")
        rep = np.asarray(representatives)
        rot = np.asarray(rotations)
        tau = np.asarray(translations)
        if (
            rep.ndim != 2
            or rep.shape[1] != 3
            or not np.issubdtype(rep.dtype, np.integer)
        ):
            raise ValueError("invalid reciprocal star representatives")
        if (
            rot.ndim != 3
            or rot.shape[1:] != (3, 3)
            or not np.issubdtype(rot.dtype, np.integer)
        ):
            raise ValueError("invalid reciprocal rotations")
        checked(tau, (len(rot), 3), "translations", real=True)
        if len(rot) == 0:
            raise ValueError("empty symmetry group")
        self.representatives = rep.copy()
        self.members = []
        for g in rep:
            star = {}
            for r, t in zip(rot, tau):
                q = tuple(r @ g)
                angle = 2 * np.pi * (g @ t)
                phase = np.exp(1j * angle) if complex_phases else complex(np.cos(angle))
                star.setdefault(q, []).append(phase)
            indices = []
            factors = []
            for q, phases in star.items():
                if any(abs(v) > n // 2 for v, n in zip(q, self.shape)):
                    raise ValueError("Fourier star outside the FFT grid")
                indices.append(tuple(v % n for v, n in zip(q, self.shape)))
                factors.append(np.mean(phases).conjugate())
            self.members.append((np.array(indices, int), np.array(factors)))

    def action(self, fourier):
        f = checked(fourier, self.shape, "Fourier density")
        return np.array([np.sum(f[tuple(i.T)] * phase) for i, phase in self.members])

    def pullback(self, covector):
        y = checked(covector, (len(self.members),), "star covector")
        out = np.zeros(self.shape, complex)
        for value, (indices, phase) in zip(y, self.members):
            np.add.at(out, tuple(indices.T), value * phase.conj())
        return out


class InterstitialDensity:
    """FFT reconstruction without diagonalizing the signed density matrix.

    coefficients contains ONLY plane-wave rows, never LO rows. Integer G labels,
    cell volume and the FFT grid are fixed. A central periodic wavefunction cache
    is optional; only one k point need reside in memory. Applying rho in
    coefficient space before FFT avoids an O(Ngrid*Nband**2) real-space product.
    The output is the native star-summed full-cell Fourier continuation used
    for interstitial density; no sphere mask is multiplied into this field.
    """

    def __init__(
        self, coefficients, labels, volume, stars, *, cache_wavefunctions=True
    ):
        c = np.asarray(coefficients)
        g = np.asarray(labels)
        if c.ndim != 2 or min(c.shape) < 1:
            raise ValueError("invalid plane-wave coefficients")
        self.c = checked(c, c.shape, "plane-wave coefficients").copy()
        if g.shape != (len(c), 3) or not np.issubdtype(g.dtype, np.integer):
            raise ValueError("invalid plane-wave labels")
        if not np.isfinite(volume) or volume <= 0:
            raise ValueError("invalid cell volume")
        self.volume = volume
        self.stars = stars
        self.shape = stars.shape
        self.ng = int(np.prod(self.shape))
        if any(np.ptp(g[:, j]) * 2 >= self.shape[j] for j in range(3)):
            raise ValueError("FFT grid aliases plane-wave density products")
        self.index = tuple((g % self.shape).T)
        if len(set(map(tuple, g))) != len(g):
            raise ValueError("duplicate plane-wave labels")
        self.psi = self._waves(self.c) if cache_wavefunctions else None

    def _waves(self, c):
        # Transform contiguous spatial grids, one per band. Keeping bands
        # last during the FFT gives long-stride transforms and substantially
        # more memory traffic. Expose the same (..., band) view afterwards;
        # no wavefunction cache or approximation is introduced.
        grid = np.zeros((c.shape[1], *self.shape), complex)
        grid[(slice(None),) + self.index] = c.T
        waves = np.fft.ifftn(grid, axes=(-3, -2, -1)) * self.ng
        return np.moveaxis(waves, 0, -1)

    def _wave_pullback(self, value):
        waves = np.moveaxis(value, -1, 0)
        grid = np.fft.fftn(waves, axes=(-3, -2, -1))
        return grid[(slice(None),) + self.index].T

    def value(self, rho):
        rho = checked(rho, (self.c.shape[1],) * 2, "band density", hermitian=True)
        psi = self.psi if self.psi is not None else self._waves(self.c)
        n = np.einsum("...n,...n->...", self._waves(self.c @ rho), psi.conj()).real
        return self.stars.action(np.fft.fftn(n) / (self.ng * self.volume))

    def action(self, rho, rho_dot, coefficient_dot=None):
        result = self.value(rho_dot)
        if coefficient_dot is not None:
            rho = checked(
                rho, (self.c.shape[1],) * 2, "central density", hermitian=True
            )
            dc = checked(coefficient_dot, self.c.shape, "coefficient derivative")
            n = (
                2
                * np.einsum(
                    "...n,...n->...", self._waves(dc), self._waves(self.c @ rho).conj()
                ).real
            )
            result += self.stars.action(np.fft.fftn(n) / (self.ng * self.volume))
        return result

    def pullback(self, covector, rho=None):
        # FFT/N has adjoint IFFT. The real projection is essential: the
        # independent coordinates of the physical grid density are real.
        v = np.fft.ifftn(self.stars.pullback(covector)).real / self.volume
        psi = self.psi if self.psi is not None else self._waves(self.c)
        vc = self._wave_pullback(v[..., None] * psi)
        rhobar = herm(self.c.conj().T @ vc)
        if rho is None:
            return rhobar
        rho = checked(rho, rhobar.shape, "central density", hermitian=True)
        return rhobar, 2 * vc @ rho


class SpatialValenceDensity:
    """Connect MT and IR operations to the common Density.vector coordinates.

    atoms maps inequivalent atom indices to MuffinTinDensitySum; interstitial
    contains one plane-wave map per k point. The representation of rho is
    fixed across these maps. Potential-induced radial changes contribute in
    MT only: at fixed geometry and central coefficients, interstitial trial
    plane waves do not change. Electronic eigenvector response belongs in
    rho_dot and must not also be passed as a change of plane-wave coefficients.

    Core density and output charge normalization remain separate operations
    in ChargeResponseChain, preventing accidental double normalization here.
    """

    def __init__(self, template, atoms, interstitial):
        self.template = template.zeros()
        self.atoms = dict(atoms)
        if set(self.atoms) != {key[0] for key in template.rad}:
            raise ValueError("missing or extra muffin-tin atom")
        if not self.atoms:
            raise ValueError("no muffin-tin densities")
        self.rho = next(iter(self.atoms.values())).rho
        for atom, m in self.atoms.items():
            if m.rho.shape != self.rho.shape or not np.array_equal(m.rho, self.rho):
                raise ValueError("inconsistent central density matrices")
            hs = {key[1:] for key in template.rad if key[0] == atom}
            if set(m.kernel.layout.harmonics) != hs:
                raise ValueError("density harmonic layout mismatch")
            for h in hs:
                if template.rad[(atom, *h)].shape != (m.kernel.nr,):
                    raise ValueError("density radial grid length mismatch")
        self.interstitial = list(interstitial)
        if len(self.interstitial) != len(self.rho):
            raise ValueError("missing interstitial k point")
        for m in self.interstitial:
            if m.c.shape[1] != self.rho.shape[1] or not np.array_equal(
                m.stars.representatives, np.array(list(template.pw))
            ):
                raise ValueError("interstitial band/star layout mismatch")

    def _assemble(self, radial, fourier):
        out = self.template.zeros()
        for atom, values in radial.items():
            for h, value in zip(self.atoms[atom].kernel.layout.harmonics, values):
                out.rad[(atom, *h)] = value
        for g, v in zip(out.pw, fourier):
            out.pw[g] = v
        return out.vector()

    def value(self):
        return self._assemble(
            {a: m.value() for a, m in self.atoms.items()},
            sum(m.value(r) for m, r in zip(self.interstitial, self.rho)),
        )

    def action(self, rho_dot=None, matching_dots=None, primitive_dots=None):
        for changes in [matching_dots, primitive_dots]:
            if changes is not None and set(changes) != set(self.atoms):
                raise ValueError("missing or extra atom derivatives")
        rad = {
            a: m.action(
                rho_dot,
                None if matching_dots is None else matching_dots[a],
                None if primitive_dots is None else primitive_dots[a],
            )
            for a, m in self.atoms.items()
        }
        ir = np.zeros(len(self.template.pw), complex)
        if rho_dot is not None:
            for m, r in zip(self.interstitial, rho_dot):
                ir += m.value(r)
        return self._assemble(rad, ir)

    def pullback(self, covector, *, include_basis=True, include_density=True):
        v = checked(
            covector,
            self.template.vector().shape,
            "spatial density covector",
            real=True,
        )
        y = self.template.from_vector(v)
        rb = np.zeros_like(self.rho) if include_density else None
        zb = {}
        pb = {}
        for a, m in self.atoms.items():
            ym = np.array([y.rad[(a, *h)] for h in m.kernel.layout.harmonics])
            r, z, p = m.pullback(
                ym, include_basis=include_basis, include_density=include_density
            )
            if include_density:
                rb += r
            if include_basis:
                zb[a] = z
                pb[a] = p
        if include_density:
            yi = np.array(list(y.pw.values()))
            for k, m in enumerate(self.interstitial):
                rb[k] += m.pullback(yi)
        return rb, zb, pb

    def fixed_basis_operator(self):
        """R in ChargeResponseChain; input matrices already contain k/spin weights."""
        from stress_response_adjoint import HermitianCoordinates, PairedOperator

        coordinates = HermitianCoordinates(self.rho.shape)
        operator = PairedOperator(
            coordinates.size,
            len(self.template.vector()),
            lambda v: self.action(coordinates.unpack(v)),
            lambda y: coordinates.pack(self.pullback(y, include_basis=False)[0]),
            "native spatial valence-density reconstruction",
        )
        return operator, coordinates
