"""Differentiate one second-variational k point in the central LAPW frame.

The scalar material maps have retained/retained (RR) and outside/retained
(QR) blocks. On output the two coherent spinor components are traced into
one charge. DMFT2's band density, including the k weight, is 2*w*rho.
The factor one half below appears once for every physical density trace.
"""

from types import SimpleNamespace
from pathlib import Path
import numpy as np
from response_parallel import SERIAL
from lapw_complement_response import LAPWComplementBasis
from lapw_energy_response import energy_direction
from lapw_energy_readout import band_energy_gradient, density_energy_gradient
from lapw_geometry_response import native_band_geometry_source
from lapw_strain_charge import native_spatial_geometry_source
from lapw_soc_subspace import SecondVariationalResponse, hermitian
from lapw_valence_state import _prepare_representation


class SOCValenceAtK:
    """Select the radial, scalar-window and spinor-window paths exactly once."""

    def __init__(
        self,
        directory,
        case,
        ik,
        mat,
        k,
        labels,
        h,
        o,
        a,
        energies,
        soc,
        spinors,
        spin_energies,
        rho,
        projectors,
        pc,
        search,
        fingerprints,
        stack,
        library,
    ):
        self.directory, self.case, self.ik = Path(directory), case, ik
        self.problem = SecondVariationalResponse(
            h, o, a, energies, soc, spinors, spin_energies
        )
        self.project = projectors.prepare(k, labels)
        self.rho = rho
        self.search = search
        self.o = o
        self.scale = 0.5
        blocks = self.problem.spin_blocks
        scalar_rho = self.scale * sum((c @ rho) @ c.conj().T for c in blocks)
        scal = self.problem.scalar
        mix = np.column_stack((a, scal.outside))
        mix.setflags(write=False)
        self.dual = o @ mix
        nr = len(energies)
        basis = LAPWComplementBasis(
            a, scal.outside, SimpleNamespace(nq=len(a) - nr), {}, mix
        )
        model = SimpleNamespace(
            eps=energies[None, :],
            parallel=SERIAL,
            raw=np.zeros((1, nr, 1, 1), complex),
            central_density=lambda: scalar_rho[None, :, :],
        )
        rec = SimpleNamespace(
            kpoint=ik, reciprocal_k=k, nbasis=len(a), nlo=mat.nlo_total, labels=labels
        )
        from copy import copy

        mat = copy(mat)
        mat.cv = [(energies, mix, labels)]
        self.val = _prepare_representation(
            case,
            self.directory,
            mat,
            [rec],
            [k],
            [basis],
            model,
            nr,
            stack,
            complex_vectors=mat.complex,
            library=library,
            cache_wavefunctions=False,
            matrix_quadrature="lapwso",
            fft_grid="fast",
            parallel=SERIAL,
            basis_cache_mb=0,
            disk_matching=True,
            spin_channel="",
            dft=False,
            scalar_spin_frame=True,
            geometry_only=pc is None,
        )
        self.val.fingerprints = fingerprints
        self.projector_contexts = projectors.contexts
        if pc is None:
            return
        self.select, _ = pc.hartree.native_potential_output(
            tuple(self.val.bands.spherical.slices),
            {a: op.layout.harmonics for a, (op, _) in self.val.bands.ns.items()},
            self.val.bands.ir[0].representatives,
            complex_interstitial=self.val.bands.input.complex_interstitial,
        )
        self.slices = {
            a: sl
            for (a, l, m), sl in zip(pc.hartree.keys, pc.hartree.slices)
            if l == m == 0
        }

    def _radial_directions(self, v):
        energies = {
            key: matrix @ v[self.slices[key[0]]] for key, matrix in self.search.items()
        }

        def apply(key, op):
            dv = v[self.slices[key[0]]]
            de = energies.get(key)
            return op.action(dv) if de is None else op.action(dv, de[0], de[1:])

        radial = {
            key: apply(key, ch)
            for key, ch in self.problem.soc.geometry.channels.items()
        }
        projector = {key: apply(key, ch) for key, ch in self.projector_contexts.items()}
        return energies, radial, projector

    def action(self, v, velocity=None, *, include_charge=True):
        """Return dHspin, raw dPhi and charge contribution at fixed band occupations."""
        val, p = self.val, self.problem
        energies, radial, projector = self._radial_directions(v)
        potential = self.select.matvec(v)
        band = val.bands.action(potential)
        spatial = (
            val.direct.operator.matvec(val.bands.spherical_selector.matvec(potential))
            if include_charge
            else None
        )
        if energies:
            part, charge = energy_direction(val, energies)
            band += part
            if include_charge:
                spatial += charge
        if velocity is not None:
            e = np.zeros((3, 3))
            geo, _ = native_band_geometry_source(
                val, self.directory, self.case, e, displacements=velocity
            )
            band += geo
            if include_charge:
                spatial += native_spatial_geometry_source(
                    val, self.directory, self.case, e, displacements=velocity
                )
        hh, oo, _, hq, oq = val.bands.output.unpack(band)
        scal = p.scalar
        dh = hh[0] - 0.5 * (scal.energy[:, None] + scal.energy) * oo[0]
        da = (
            scal.outside @ ((hq[0] - oq[0] * scal.energy) * scal.reciprocal_gaps)
            - 0.5 * scal.a @ oo[0]
        )
        hso = p.soc.action(radial, da, velocity)
        n = len(dh)
        hso[:n, :n] += dh
        hso[n:, n:] += dh
        hk, dc = p.spin.action(hso)
        dpsi = np.stack(
            [da @ c + scal.a @ dd for c, dd in zip(p.spin_blocks, np.split(dc, 2))]
        )
        du = self.project.evaluate(
            p.waves, dpsi=dpsi, radial=projector, velocity=velocity
        )
        if not include_charge:
            return hk, du, None
        # Only RR and QR blocks are needed because the central density has no
        # QQ block. One conjugate QR leg is supplied by spatial.action.
        dr = self.scale * sum(
            (self.dual.conj().T @ x) @ (self.rho @ c.conj().T)
            for x, c in zip(dpsi, p.spin_blocks)
        )
        spatial += val.spatial.action(
            val.spatial.coordinates.pack((dr[:n] + dr[:n].conj().T)[None], [dr[n:]])
        )
        return hk, du, spatial

    def frozen_soc_force_covectors(self, hbar, raw_bar):
        """Derivative of the embedded SOC operator at frozen radial shape.

        In the scalar LAPW frame its covariant matrix is Q^H W Q,
        with Q=A A^H O. Inside S, (dQ) A = A_out Z. The retained
        metric connection -A (A^H dO A)/2 belongs to the usual overlap
        force. Including it here would count part of that force twice.

        hbar includes the physical density and normalized DMFT-window
        gradient. The change of eigenvectors above the last occupied or
        correlated spinor also changes the raw projector. Reverse that
        dependence once into the same second-variational Hamiltonian.
        """
        p = self.problem
        waves_bar, _, _ = self.project.pullback(p.waves, raw_bar)
        dc_bar = np.concatenate([p.scalar.a.conj().T @ x for x in waves_bar])
        spin_bar = p.spin.pullback(hbar, dc_bar, metric=False)
        _, scalar_bar, rigid_bar = p.soc.pullback(spin_bar)
        # Pair once with the QR part of the scalar perturbation for any atom.
        cross_bar = (p.scalar.outside.conj().T @ scalar_bar) * p.scalar.reciprocal_gaps
        return rigid_bar, cross_bar

    def frozen_soc_direction(self, velocity, covectors, geo=None):
        """Return the rigid and scalar-window energy derivatives, in Ry/Bohr."""
        if geo is None:
            geo, _ = native_band_geometry_source(
                self.val,
                self.directory,
                self.case,
                np.zeros((3, 3)),
                displacements=velocity,
            )
        _, _, _, hq, oq = self.val.bands.output.unpack(geo)
        rigid, cross = covectors
        return np.array(
            [
                np.vdot(rigid, velocity).real,
                np.vdot(cross, hq[0] - oq[0] * self.problem.scalar.energy).real,
            ]
        )

    def occupation_charge(self, drho):
        p = self.problem
        rr = self.scale * sum(c @ drho @ c.conj().T for c in p.spin_blocks)
        coords = self.val.spatial.coordinates
        return self.val.spatial.action(
            coords.pack(rr[None], [np.zeros(coords.shapes[0], complex)])
        )

    def charge_readout(self, y):
        """Adjoint of spatial reconstruction to band density and wavefunctions."""
        p, val = self.problem, self.val
        rr, qr = val.spatial.coordinates.unpack(val.spatial.pullback(y))
        nr = rr.shape[-1]
        g = self.scale * sum(c.conj().T @ rr[0] @ c for c in p.spin_blocks)
        ypsi = self.dual[:, :nr] @ rr[0] + self.dual[:, nr:] @ qr[0]
        waves = np.stack(
            [2 * self.scale * ypsi @ (c @ self.rho) for c in p.spin_blocks]
        )
        return hermitian(g), waves

    def pullback(self, hbar, raw_bar, waves_bar, y):
        """Return the central potential covector; position sources use action."""
        p, val = self.problem, self.val
        dw, radial_proj, _ = self.project.pullback(p.waves, raw_bar)
        dw += waves_bar
        da = sum(x @ c.conj().T for x, c in zip(dw, p.spin_blocks))
        dc = np.concatenate([p.scalar.a.conj().T @ x for x in dw])
        hso = p.spin.pullback(hbar, dc, metric=False)
        radial, db, _ = p.soc.pullback(hso)
        da += db
        a, b, e = p.scalar.a, p.scalar.outside, p.scalar.energy
        n = len(e)
        hrr = hermitian(hso[:n, :n] + hso[n:, n:])
        orr = -0.5 * ((e[:, None] + e) * hrr + hermitian(a.conj().T @ da))
        hqr = 0.5 * (b.conj().T @ da) * p.scalar.reciprocal_gaps
        obr = -hqr * e
        # Packed QR covectors include both Hermitian legs.
        packed = val.bands.output.pack(
            hrr[None], orr[None], np.zeros_like(val.model.raw), [hqr], [obr]
        )
        gv = val.bands.pullback(packed)
        has_density = bool(np.any(y))
        if has_density:
            gv += val.bands.spherical_selector.rmatvec(val.direct.operator.rmatvec(y))
        result = self.select.rmatvec(gv)
        ge = {}
        if self.search:
            ge = band_energy_gradient(val.bands.spherical, packed)
            if has_density:
                for key, add in density_energy_gradient(val.direct, y).items():
                    ge[key] += add
        for bank, bars in (
            (p.soc.geometry.channels, radial),
            (self.projector_contexts, radial_proj),
        ):
            for key, bar in bars.items():
                dv, de, dlo = bank[key].pullback(*bar)
                result[self.slices[key[0]]] += dv
                if key in self.search:
                    ge[key] += np.r_[de, dlo]
        for key, matrix in self.search.items():
            result[self.slices[key[0]]] += matrix.T @ ge[key]
        return result
