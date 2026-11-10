"""Connect prepared native LAPW components to frozen-Sigma screening.

This is the material adapter between the separately validated potential/core
and full-k valence states and ChargeResponseChain. It performs no native
launches or file writes. Preparation is shared by all forward/adjoint solves.
The caller owns the valence state's radial contexts and must keep them open.

No strain source or stress observable is guessed here. A native hydrostatic
source can be supplied to solve_source for an independent screening check.
Analytic geometry sources belong in chain.solve_perturbation. Publication of
a corrected stress tensor is a separate, still unfinished operation.
"""

from dataclasses import dataclass
import numpy as np
from lapw_complement_response import complement_electronic_operator
from lapw_density_response import checked
from lapw_response_chain import ChargeResponseChain


def _same_density_layout(left, right):
    """Equal size is insufficient: ordering, grids and star labels must agree."""
    if (
        tuple(left.rad) != tuple(right.rad)
        or tuple(left.pw) != tuple(right.pw)
        or not np.array_equal(left.mult, right.mult)
        or len(left.grids) != len(right.grids)
        or any(not np.array_equal(x, y) for x, y in zip(left.grids, right.grids))
    ):
        raise ValueError("central density layouts differ")


def _same_structure(left, right):
    if (
        left.nat != right.nat
        or left.lattyp.strip() != right.lattyp.strip()
        or left.mode != right.mode
    ):
        raise ValueError("central structures differ")
    for key in (
        "a",
        "b",
        "c",
        "alpha",
        "beta",
        "gamma",
        "mult",
        "jrj",
        "r0",
        "rmt",
        "pos",
        "timat",
        "tau",
        "rotloc",
        "iatnr",
        "Znuc",
    ):
        if not np.array_equal(getattr(left, key), getattr(right, key)):
            raise ValueError("central structures differ: " + key)


@dataclass(frozen=True)
class ScreenedDensityResult:
    """The density response and INDUCED potential K q, not total strained V.

    Any explicit geometric/ionic potential still needs to be added once.
    Values use the potential/core archive's real native field coordinates.
    """

    density: np.ndarray
    induced_potential: np.ndarray
    report: dict


class NativeScreenedState:
    """One global J=N[REB+D+C]K with shared projector normalization and mu.

    Parameters are independently prepared central states. No per-k screening,
    normalization, or chemical-potential fit occurs. ``raw_core`` must be the
    native unnormalized central LCORE output. Normalization is differentiated
    at the total raw output, never at the input or already normalized density.

    Fixed-N refers to the ELECTRONIC ensemble. MIXER's fixed target and its
    normalization are separate choices, also present when fixed_number=False.
    Radial energies, RMT, local frames and Sigma/DC stay fixed in ordinary J.
    The allowed geometry/spin/projector scope remains that of the loaders.
    Full snapshot provenance is the caller/loader's responsibility; this
    adapter additionally compares the shared structure and spherical-potential
    fingerprints. It is not an automatic exporter or execution-cache loader.
    """

    def __init__(
        self,
        valence,
        potential_core,
        raw_core,
        *,
        target_charge,
        fixed_number,
        normalization_enabled=True,
        closure_tolerance=2e-8,
    ):
        if not isinstance(fixed_number, (bool, np.bool_)):
            raise ValueError("explicit fixed_number policy required")
        if not np.isfinite(closure_tolerance) or closure_tolerance <= 0:
            raise ValueError("positive central closure tolerance required")
        _same_structure(valence.structure, potential_core.structure)
        _same_density_layout(valence.reference, potential_core.input_density)
        if (
            not np.array_equal(raw_core.mult, valence.reference.mult)
            or len(raw_core.grids) != len(valence.reference.grids)
            or any(
                not np.array_equal(x, y)
                for x, y in zip(raw_core.grids, valence.reference.grids)
            )
        ):
            raise ValueError("central core density layout differs")
        self.valence = valence
        self.potential_core = potential_core
        self.fixed_number = bool(fixed_number)
        select, coordinates = potential_core.hartree.native_potential_output(
            tuple(valence.bands.spherical.slices),
            {a: op.layout.harmonics for a, (op, ks) in valence.bands.ns.items()},
            valence.bands.ir[0].representatives,
            complex_interstitial=valence.bands.input.complex_interstitial,
        )
        if (
            coordinates.size != valence.bands.input.size
            or coordinates.nonspherical_shapes
            != valence.bands.input.nonspherical_shapes
            or coordinates.spherical_size != valence.bands.input.spherical_size
            or coordinates.nstar != valence.bands.input.nstar
        ):
            raise ValueError("native potential coordinates differ")
        self.potential_selector = select

        # These physical overlaps cannot be chosen independently even if both
        # prepared objects happen to have the same array dimensions.
        case_keys = [k for k in valence.fingerprints if k.endswith(".struct")]
        if len(case_keys) != 1:
            raise ValueError("valence snapshot needs a unique structure fingerprint")
        case = case_keys[0][:-7]
        for key in ("struct", "vsp"):
            expected = potential_core.metadata["source_fingerprints"].get(key)
            if (
                expected is None
                or valence.fingerprints.get(case + "." + key) != expected
            ):
                raise ValueError(
                    "potential/core and valence snapshot fingerprints differ: " + key
                )

        raw = valence.reference.copy()
        raw.add(raw_core)
        actual_core = raw.vector() - valence.reference.vector()
        expected_core = np.zeros_like(actual_core)
        for a, op in potential_core.core.atoms.items():
            expected_core[potential_core.core.slices[a]] = op.density
        core_error = np.linalg.norm(actual_core - expected_core) / max(
            np.linalg.norm(expected_core), 1.0
        )
        valence_error = np.linalg.norm(
            valence.spatial.value() - valence.reference.vector()
        ) / max(np.linalg.norm(valence.reference.vector()), 1.0)
        if max(core_error, valence_error) > closure_tolerance:
            raise ValueError(
                f"central native density reconstruction failed: core={core_error:g}, valence={valence_error:g}"
            )
        self.normalization = potential_core.normalization(
            raw, target_charge, enabled=normalization_enabled
        )
        electronic, inputs, outputs = complement_electronic_operator(
            valence.complement, fixed_number=self.fixed_number
        )
        if (
            inputs.size != valence.bands.output.size
            or outputs.size != valence.spatial.coordinates.size
        ):
            raise ValueError("compact electronic coordinate mismatch")
        self.chain = ChargeResponseChain(
            potential=potential_core.potential_operator,
            basis=valence.bands.operator @ select,
            electronic=electronic,
            density=valence.spatial.operator,
            basis_density=valence.direct.operator
            @ valence.bands.spherical_selector
            @ select,
            core=potential_core.core.operator,
            normalization=self.normalization.operator,
        )
        self.jacobian = self.chain.jacobian
        self.screened = self.chain.screened
        self.metadata = dict(
            fixed_number=self.fixed_number,
            target_charge=float(target_charge),
            normalization_enabled=bool(normalization_enabled),
            normalization_alpha=float(self.normalization.alpha),
            central_valence_relative_error=float(valence_error),
            central_core_relative_error=float(core_error),
            central_fixed_point_residual_norm=float(
                np.linalg.norm(
                    self.normalization.central - potential_core.input_density.vector()
                )
            ),
            kpoints=len(valence.records),
            retained_bands=valence.model.eps.shape[1],
            matrix_quadrature=valence.matrix_quadrature,
            stress_applied=False,
        )

    def solve_source(self, source, **controls):
        """Solve one explicit normalized source, with a checked true residual.

        The source must use central Density.vector coordinates, including all
        geometry/transport/normalization contributions exactly once. Native
        +/-strain evaluations are a reference source; they do not constitute
        an analytic implementation of the six strain sources.
        This also works away from a fixed point as a local tangent equation.
        The caller must check central convergence before interpreting it as
        the derivative of a self-consistent material calculation.
        """
        source = checked(
            source, (self.jacobian.input_size,), "normalized density source", real=True
        )
        density, report = self.screened.solve(source, **controls)
        potential = self.chain.potential.matvec(density)
        report = dict(
            report,
            **self.metadata,
            source_charge=float(
                (self.normalization.mt + self.normalization.ir) @ source
            ),
            response_charge=float(
                (self.normalization.mt + self.normalization.ir) @ density
            ),
        )
        return ScreenedDensityResult(density, potential, report)

    def observable_gradient(self, density_covector, potential_covector):
        """c = c_density + K^T c_potential for ONE specified scalar observable.

        This is an exact coordinate pullback, not an assumption that the
        observable is energy or stress. Its explicit geometry derivative is
        supplied separately. One adjoint can serve many sources for that same
        scalar; six unrelated observable gradients require separate solves.
        """
        d = checked(
            density_covector, (self.jacobian.input_size,), "density covector", real=True
        )
        v = checked(
            potential_covector,
            (self.chain.potential.output_size,),
            "potential covector",
            real=True,
        )
        return d + self.chain.potential.rmatvec(v)
