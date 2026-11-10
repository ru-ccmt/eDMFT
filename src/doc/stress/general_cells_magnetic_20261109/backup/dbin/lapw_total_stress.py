"""Complete native frozen-Sigma stress for symmetry-enforced diagonal tensors.

All functional terms are combined before contraction. One screened adjoint
serves every independent strain. This module consumes prepared central data;
it neither changes the calculation nor appends a partial correction to mixer.

Supported: the primitive orthogonal, nonmagnetic/no-SOC, fixed-N projector-5
scope of the native loaders, fixed RMT and local frames, successful ordinary
radial-energy searches with fixed IN1 reference EF. Symmetry must forbid every
off-diagonal stress component. General shear, magnetic and fixed-mu publication
require their own validation and are deliberately rejected here.
"""

from dataclasses import dataclass
import time
import numpy as np
from scipy.linalg import null_space
from lapw_density_response import checked
from lapw_strain_field import diagonal_strain_potential
from lapw_strain_charge import (
    hydrostatic_charge_source,
    PreparedChargeSourceContraction,
    native_spatial_geometry_source,
)
from lapw_geometry_response import native_band_geometry_source
from lapw_klein_functional import PreparedKleinFunctional
from lapw_energy_response import SearchedEnergyResponse
from stress_response_adjoint import RY_BOHR3_TO_GPA


def diagonal_invariant_strains(structure):
    """A complete strain basis when point symmetry forces stress to be diagonal.

    Determine the invariant subspace of *all six* symmetric tensor components
    first. Testing just diagonal strains would silently discard allowed shear
    stress in monoclinic/triclinic cases. The hydrostatic direction comes first;
    averaged zz, yy, xx directions complete its span without arbitrary signs.
    Fractional symmetry operations are converted to Cartesian coordinates.
    """
    s = structure
    if s.lattyp.strip() != "P" or not np.allclose(
        [s.alpha, s.beta, s.gamma], 90, rtol=0, atol=1e-10
    ):
        raise NotImplementedError("primitive orthogonal cells required")
    cell = np.diag([s.a, s.b, s.c])
    rotations = [cell @ r @ np.linalg.inv(cell) for r in s.timat]
    if not rotations or any(
        not np.allclose(r.T @ r, np.eye(3), atol=1e-10, rtol=0) for r in rotations
    ):
        raise ValueError("invalid Cartesian point symmetry")
    basis = [np.diag(np.eye(3)[i]) for i in range(3)]
    for i, j in [(1, 2), (0, 2), (0, 1)]:
        e = np.zeros((3, 3))
        e[i, j] = e[j, i] = 1 / np.sqrt(2)
        basis.append(e)
    constraints = np.vstack(
        [np.column_stack([(r @ e @ r.T - e).ravel() for e in basis]) for r in rotations]
    )
    invariant = null_space(constraints, rcond=1e-10)
    if invariant.shape[1] == 0 or np.linalg.norm(invariant[3:]) > 1e-9:
        raise NotImplementedError("symmetry allows shear; expanded sources required")
    directions = [np.eye(3)]
    for i in [2, 1, 0]:
        e = basis[i]
        candidate = sum(r @ e @ r.T for r in rotations) / len(rotations)
        columns = np.column_stack([v.ravel() for v in directions])
        residual = (
            candidate.ravel()
            - columns @ np.linalg.lstsq(columns, candidate.ravel(), rcond=None)[0]
        )
        if np.linalg.norm(residual) > 1e-9:
            directions.append(candidate)
    if len(directions) != invariant.shape[1]:
        raise ValueError("incomplete invariant diagonal strain basis")
    return np.array(directions)


def stress_from_contractions(strains, derivatives, volume):
    """Reconstruct sigma from D_nu=Omega Tr(sigma e_nu), in Ry/Bohr^3.

    This is called only with a complete invariant basis from the function
    above. A hydrostatic derivative alone determines a tensor only in cubic
    symmetry. Pressure in the thermodynamic compression convention is -tr/3.
    """
    e = np.asarray(strains)
    d = checked(derivatives, (len(e),), "functional derivatives", real=True)
    if not np.isfinite(volume) or volume <= 0:
        raise ValueError("positive finite volume required")
    gram = np.einsum("aij,bij->ab", e, e)
    if np.linalg.matrix_rank(gram) != len(e):
        raise ValueError("dependent strain directions")
    return np.einsum("a,aij->ij", np.linalg.solve(gram, d / volume), e)


@dataclass
class PreparedStrainSource:
    strain: np.ndarray
    source: np.ndarray
    geometry: np.ndarray
    potential: np.ndarray
    report: dict
    induced: float | None = None


class PreparedTotalStress:
    """One complete functional, radial-search completion and screened adjoint.

    The material density transport is an explicit coordinate choice:
    tau_nu=(tr(e_nu)/3) tau_hydro. It need not affinely deform the physical
    density for anisotropic strain: the solved q supplies the remaining
    response. The identical convention enters charge sources and the scalar,
    so the final derivative is coordinate independent at the fixed point.

    Central states must be freshly exported using the same IN1/search policy
    as the finite-difference reference. The fixed-point check is distinct from
    the native reconstruction check and catches stale-but-internally-consistent
    vector exports. Its norm is in native density coordinates, not electrons.
    """

    def __init__(
        self,
        state,
        directory,
        case,
        *,
        complex_vectors,
        library,
        fixed_point_tolerance=1e-3,
        progress=None,
    ):
        residual = state.metadata["central_fixed_point_residual_norm"]
        if not state.fixed_number:
            raise NotImplementedError("fixed-N total stress is validated here")
        if not np.isfinite(fixed_point_tolerance) or fixed_point_tolerance <= 0:
            raise ValueError("positive fixed-point tolerance required")
        if residual > fixed_point_tolerance:
            raise ValueError(f"central charge map is not converged: {residual:g}")
        self.strains = diagonal_invariant_strains(state.potential_core.structure)
        # Structural symmetry alone does not prove that a noisy frozen Sigma
        # respects it. Its Hermitian real/imaginary parts must each be invariant
        # under the actual orbital/cix map. Otherwise only the tensor projected
        # onto the imposed charge symmetry is determined by these sources.
        meta = state.valence.model.meta
        errors = []
        for sigma in [meta.sigma_inf] + [
            meta.sigma[..., n] for n in range(len(meta.omega))
        ]:
            dagger = sigma.swapaxes(0, 1).conj()
            for part in [(sigma + dagger) / 2, (sigma - dagger) / (2j)]:
                errors.append(
                    float(
                        np.linalg.norm(meta.sym(part) - part)
                        / max(np.linalg.norm(part), 1e-30)
                    )
                )
        self.sigma_symmetry_error = max(errors, default=0.0)
        self.state, self.directory, self.case = state, directory, case
        self.functional = PreparedKleinFunctional(state)
        self.energy = SearchedEnergyResponse(
            state,
            self.functional,
            directory,
            case,
            complex_vectors=complex_vectors,
            library=library,
            progress=progress,
        )

    def prepare_source(self, strain, *, readout=None):
        """Build one explicit source, including radial-search geometry input."""
        e = checked(strain, (3, 3), "strain direction", real=True)
        potential, _ = diagonal_strain_potential(self.state.potential_core, e)
        if readout is not None:
            if readout.state is not self.state:
                raise ValueError("source contraction belongs to a different state")
            band, report = native_band_geometry_source(
                self.state.valence, self.directory, self.case, e
            )
            geometry = native_spatial_geometry_source(
                self.state.valence, self.directory, self.case, e
            )
            induced = readout.contraction(band, geometry, potential, e)
            induced += readout.adjoint @ (self.energy.P @ (self.energy.A @ potential))
            report.update(
                source_contracted_by_adjoint=True,
                radial_energy_policy="differentiated native searches",
            )
            return PreparedStrainSource(
                e, np.empty(0), band, potential, report, float(induced)
            )
        source, parts, report = hydrostatic_charge_source(
            self.state, self.directory, self.case, strain=e, potential_source=potential
        )
        source = self.energy.complete_source(source, potential)
        report.update(
            radial_energy_policy="differentiated native searches",
            scalar_transport_rate=float(np.trace(e) / 3),
        )
        return PreparedStrainSource(
            e, source, parts["band_geometry"], potential, report
        )

    def contract(
        self,
        sources,
        *,
        solved_adjoint=None,
        initial=None,
        rtol=2e-12,
        atol=1e-11,
        restart=30,
        max_cycles=8,
    ):
        """Compute the complete tensor with one adjoint and true-residual check."""
        if len(sources) != len(self.strains) or not np.allclose(
            [x.strain for x in sources], self.strains, atol=1e-12, rtol=0
        ):
            raise ValueError("complete ordered invariant sources required")
        start = time.perf_counter()
        geometry = np.column_stack([x.geometry for x in sources])
        potential = np.column_stack([x.potential for x in sources])
        charge = np.column_stack([x.source for x in sources])
        explicit, pieces = self.functional.explicit(geometry, potential, self.strains)
        explicit += self.energy.explicit_addition(potential)
        contracted = [x.induced is not None for x in sources]
        if any(contracted) and (not all(contracted) or solved_adjoint is None):
            raise ValueError("contracted sources require their common solved adjoint")
        if solved_adjoint is None:
            values, lam, solver = self.energy.screened.contractions(
                self.energy.charge_gradient,
                charge,
                explicit=explicit,
                initial=initial,
                rtol=rtol,
                atol=atol,
                restart=restart,
                max_cycles=max_cycles,
            )
        else:
            lam, solver = solved_adjoint
            if not solver.get("converged"):
                raise ValueError("unconverged source adjoint")
            values = explicit + (
                np.array([x.induced for x in sources])
                if all(contracted)
                else lam @ charge
            )
        sigma = stress_from_contractions(
            self.strains, values, self.state.potential_core.hartree.volume
        )
        report = dict(
            observable="complete native frozen-Sigma Klein tangent",
            stress_Ry_Bohr3=sigma.tolist(),
            stress_GPa=(sigma * RY_BOHR3_TO_GPA).tolist(),
            mean_tensile_stress_GPa=float(np.trace(sigma) * RY_BOHR3_TO_GPA / 3),
            compression_pressure_GPa=float(-np.trace(sigma) * RY_BOHR3_TO_GPA / 3),
            derivatives_Ry=values.tolist(),
            strains=self.strains.tolist(),
            explicit_Ry=explicit.tolist(),
            induced_Ry=(values - explicit).tolist(),
            sources_contracted_by_adjoint=all(contracted),
            fixed_energy_explicit_parts=pieces,
            adjoint=solver,
            central=self.state.metadata,
            search_details=self.energy.details,
            contraction_seconds=time.perf_counter() - start,
            stress_applied=False,
            full_tensor_within_symmetry=True,
            sigma_symmetry_relative_error=self.sigma_symmetry_error,
            tensor_is_symmetry_projection=self.sigma_symmetry_error > 1e-10,
            unrestricted_tensor_validated=False,
            scope="nonmagnetic fixed-N primitive orthogonal; symmetry-enforced diagonal tensor",
        )
        return sigma, lam, report
