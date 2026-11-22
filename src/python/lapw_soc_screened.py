"""Matrix-free, fixed-self-energy charge and Klein derivatives for SOC.

All inputs are central native files. Electronic normalization and the fixed-N
chemical potential couple the k points. Charge and potential grids are shared
across ranks. The large second-variational and LAPW matrices are streamed one
k point per rank, including when transposing the map for a shared adjoint.
"""

from contextlib import ExitStack
from pathlib import Path
import struct
from tempfile import TemporaryFile
import numpy as np

from force_response_modes import internal_modes
from lapw_atomic_force import atomic_charge_integral_derivative
from lapw_density_functional import PreparedDensityFunctional
from lapw_dft_vector import _record, _header
from lapw_functional_cancellation import CoreEigenvalueFunctional
from lapw_potential_core_state import fingerprint_files, load_potential_core_state
from lapw_radial_energy_policy import searched_energy_derivatives
from lapw_soc_electronic import prepare_electronic, operator_readout
from lapw_soc_projector import SpinorProjectors
from lapw_soc_response import SOCGeometry, SOCBandHamiltonian
from lapw_soc_subspace import spinor_coefficients
from lapw_soc_valence import SOCValenceAtK
from lapw_strain_field import diagonal_strain_potential
from lapw_strain_potential import masked_potential_strain_source
from lapw_valence_state import _snapshot_inputs
from response_parallel import SERIAL
from soc_scalar_window import matrices
from stress_response_adjoint import PairedOperator, ScreenedResponse
from wien_response_density import Density


def scalar_inventory(path, nat):
    offsets = []
    with path.open("rb") as f:
        for _ in range(2 * nat):
            _record(f, skip=True)
        while True:
            start = f.tell()
            head = _header(_record(f))
            if head is None:
                break
            offsets.append(start)
            for _ in range(1 + 2 * head[5]):
                if _record(f, skip=True) is None:
                    raise ValueError("truncated scalar LAPW vector")
    return offsets


def scalar_record(path, offset):
    with path.open("rb") as f:
        f.seek(offset)
        x, y, z, _, n, b, _ = _header(_record(f))
        lab = np.frombuffer(_record(f), "<i4").reshape((3, n), order="F")
        a = np.empty((n, b), order="F", dtype=complex)
        e = np.empty(b)
        for j in range(b):
            idx, e[j] = struct.unpack("<id", _record(f))
            data = _record(f)
            if idx != j + 1 or len(data) not in (8 * n, 16 * n):
                raise ValueError("invalid scalar vector band record")
            a[:, j] = np.frombuffer(data, "<c16" if len(data) == 16 * n else "<f8")
    return np.array([x, y, z]), lab, e, a


class SOCScreenedForce:
    def __init__(
        self,
        directory,
        case,
        stack,
        *,
        geometry,
        prefix,
        matrices_prefix,
        potential_core,
        library=None,
        parallel=SERIAL,
        cache_mb=128.0,
        fixed_number=True,
        operator_only=False,
        progress=lambda message: None,
    ):
        d = self.directory = Path(directory)
        self.case, self.parallel, self.progress = case, parallel, progress
        self.library, self.fixed_number = library, fixed_number
        self.stack = stack
        self.single = None
        self.pc = None if operator_only else load_potential_core_state(potential_core)
        from dmft2_step import entries

        scalar_control = entries((d / "lapw1.def").read_text())[5]
        self.mat = _snapshot_inputs(case, d, scalar_control.endswith(".in1c"))
        self.fingerprints = fingerprint_files(
            {
                case + "." + ext: d / (case + "." + ext)
                for ext in ("struct", "vns", "vsp", "clmsum")
            }
        )
        for ext in (() if operator_only else ("struct", "vsp", "clmsum")):
            expected = self.pc.metadata["source_fingerprints"].get(ext)
            if expected != self.fingerprints[case + "." + ext]:
                raise ValueError(
                    f"the potential/core archive and spinors must use the same {case}.{ext}"
                )
        self.geo = SOCGeometry(
            geometry, self.mat.s, [v[0, 0] for v in self.mat.v], stack, library
        )
        self.model, self.records = prepare_electronic(
            prefix,
            parallel=parallel,
            cache_mb=cache_mb,
            require_fresh=not operator_only,
        )
        self.project = SpinorProjectors(
            self.mat, d, case, prefix, self.model.meta, stack, library
        )
        self.search = {}
        if not operator_only:
            self.search, _ = searched_energy_derivatives(
                d,
                case,
                complex_vectors=self.mat.complex,
                library=library,
                scalar_solver="dftKS",
            )
        self.vector = d / (case + ".vector")
        self.offsets = parallel.root_call(
            lambda: scalar_inventory(self.vector, self.mat.s.nat)
        )
        if len(self.offsets) != parallel.sum(len(self.records)):
            raise ValueError("scalar and spinor meshes differ")
        self.matrices = matrices_prefix
        moments = (
            operator_readout(prefix, self.model, self.records)
            if operator_only
            else None
        )
        if moments is None:
            self.density = self.model.central_density()
            self.hbar, self.ubar = self.model.projector_gradients()
            self.hbar = 0.5 * (self.hbar + self.density)
            self.ubar *= 0.5
        else:
            self.density, self.hbar, self.ubar = moments
            self.progress("Reusing the native DMFT2 density and projector moments")
        if operator_only:
            # The frozen-shape operator contracts with occupations at the
            # current charge iterate. A converged charge fixed point is
            # required only for the screened total.
            return
        self.reference = Density(d / (case + ".clmval"), d / (case + ".struct"))
        self.n = len(self.pc.input_density.vector())
        raw = self.reference.copy()
        core = self.pc.core
        zero = np.zeros(self.n)
        for atom, op in core.atoms.items():
            zero[core.slices[atom]] = op.density
        raw = raw.from_vector(raw.vector() + zero)
        inm = (d / (case + ".inm")).read_text().splitlines()[0].split()
        if float(inm[1]) != 0 or inm[2].upper() not in ("YES", "NO"):
            raise ValueError("unsupported mixer normalization")
        self.norm = self.pc.normalization(
            raw,
            np.dot(self.mat.s.mult, self.mat.s.Znuc),
            enabled=inm[2].upper() == "YES",
        )
        central = self.pc.input_density.vector()
        self.charge_residual = np.linalg.norm(self.norm.central - central) / max(
            1.0, np.linalg.norm(central)
        )
        self.progress(
            f"Central SOC charge-map relative residual {self.charge_residual:.6g}"
        )
        if self.charge_residual > 1e-5:
            raise ValueError(
                f"central SOC charge is unconverged (relative residual {self.charge_residual:.6g})"
            )
        self.den = PreparedDensityFunctional(self.pc)
        self.core = CoreEigenvalueFunctional.from_state(self.pc)
        self.minimum_gaps = np.full(2, np.inf)
        self.singular = np.zeros(2, int)

    def points(self):
        if self.single is not None:
            yield 0, self.single
            return
        nr = self.model.eps.shape[1]
        for j, r in enumerate(self.records):
            k, l, e, a = scalar_record(self.vector, self.offsets[r.index - 1])
            es, c, u, labels = r.arrays()
            if not np.array_equal(l, labels) or r.nlo != self.mat.nlo_total:
                raise ValueError(
                    f"mismatched scalar and spinor LAPW rows at k={r.index}"
                )
            h, o, _ = matrices(self.matrices, r.index, k, l, self.mat.s, r.nlo)
            spin = spinor_coefficients(a, c[:, :nr])
            with ExitStack() as stack:
                item = SOCValenceAtK(
                    self.directory,
                    self.case,
                    r.index,
                    self.mat,
                    k,
                    l,
                    h,
                    o,
                    a,
                    e,
                    SOCBandHamiltonian(self.geo, k, l, a),
                    spin,
                    es[:nr],
                    self.density[j],
                    self.project,
                    self.pc,
                    self.search,
                    self.fingerprints,
                    stack,
                    self.library,
                )
                current = item.project.evaluate(item.problem.waves)
                if not np.allclose(current, self.model.raw[j], rtol=2e-9, atol=2e-11):
                    raise ValueError(
                        f"assembled spinor projector differs from the native export at k={r.index}"
                    )
                if len(self.records) == 1:
                    self.single = item
                    self.stack.enter_context(stack.pop_all())
                yield j, item
                del item

    def central_closure(self):
        n = np.zeros(self.n)
        for _, item in self.points():
            n += item.val.spatial.value()
            spaces = (item.problem.scalar, item.problem.spin)
            self.singular += [s.singular_pairs for s in spaces]
            self.minimum_gaps = np.minimum(
                self.minimum_gaps, [s.minimum_gap for s in spaces]
            )
            del item
        n = self.parallel.sum(n)
        target = self.reference.vector()
        error = np.linalg.norm(n - target) / max(1.0, np.linalg.norm(target))
        if error > 2e-8:
            raise ValueError(
                f"central spinor valence density disagrees with DMFT2 ({error:g})"
            )
        return error

    def valence_direction(self, v, velocity):
        dh = np.empty_like(self.hbar)
        du = np.empty_like(self.ubar)
        direct = np.zeros(self.n)
        for j, item in self.points():
            dh[j], du[j], part = item.action(v, velocity)
            direct += part
            del item
        drho, _ = self.model.response(
            dh, fixed_number=self.fixed_number, raw_derivative=du
        )
        for j, item in self.points():
            direct += item.occupation_charge(drho[j])
            del item
        lattice = np.vdot(self.hbar, dh).real + np.vdot(self.ubar, du).real
        return self.parallel.sum(direct), self.parallel.sum(float(lattice))

    def potential_readout(self, y=None):
        """Lattice gradient or physical charge-map transpose in potential coordinates."""
        result = np.zeros(self.pc.hartree.size)
        with ExitStack() as scratch:
            # The normalization and fixed-N constraint require a k reduction
            # between these sweeps. Save the first sweep's wavefunction
            # covector so the expensive spatial pullback is evaluated once.
            # For several k points use rank-private storage, released even on
            # failure. Keep only one k point's large radial operators in RAM.
            has_density = y is not None
            spool = None
            if has_density:
                if len(self.records) > 1:
                    spool = scratch.enter_context(TemporaryFile())
                g = np.empty_like(self.hbar)
                for j, item in self.points():
                    g[j], waves = item.charge_readout(y)
                    if spool is not None:
                        np.save(spool, waves, allow_pickle=False)
                        del waves
                    del item
                if spool is not None:
                    spool.seek(0)
                bh, bu = self.model.adjoint(
                    g, fixed_number=self.fixed_number, return_raw=True
                )
            else:
                bh, bu = self.hbar, self.ubar
                y = np.zeros(self.n)
            for j, item in self.points():
                if spool is not None:
                    waves = np.load(spool, allow_pickle=False)
                elif not has_density:
                    waves = np.zeros_like(item.problem.waves)
                result += item.pullback(bh[j], bu[j], waves, y)
                del item
        return self.parallel.sum(result)

    def frozen_soc_force(self, atoms):
        """Return only the explicit derivative of Q^H W_MT Q for ordinary FOR.

        Work in the current spinor frame at frozen radial shapes. There is no
        density/adjoint solve. The same SOC operator already participates in
        the full FORT energy derivative, whose publication replaces the total.
        """
        modes = internal_modes(self.mat.s, atoms)
        values = np.zeros((len(modes), 2))
        for j, item in self.points():
            covectors = item.frozen_soc_force_covectors(self.hbar[j], self.ubar[j])
            for m, (_, _, vel) in enumerate(modes):
                values[m] += item.frozen_soc_direction(vel, covectors)
            del item
        values = self.parallel.sum(values)
        vectors = np.zeros((self.mat.s.nat, 3))
        rows = []
        for (a, u, _), pair in zip(modes, values):
            rigid, window = -1000.0 * pair / self.mat.s.mult[a - 1]
            component = float(rigid + window)
            vectors[a - 1] += component * u
            rows.append(
                dict(
                    atom=int(a),
                    direction_cartesian=u.tolist(),
                    rigid_force_mRy_per_Bohr=float(rigid),
                    scalar_window_force_mRy_per_Bohr=float(window),
                    force_component_mRy_Bohr=component,
                )
            )
            self.progress(
                f"Atom {a}: explicit SOC rigid {rigid:.9f}, scalar window {window:.9f} mRy/Bohr"
            )
        return dict(
            case=self.case,
            structure=str(self.directory / (self.case + ".struct")),
            observable="frozen-radial SOC operator force",
            source_method="analytic",
            scope="add the rigid and scalar-window SOC operator to the ordinary force",
            application="add_soc_operator",
            selected_atoms=sorted(set(atoms or range(1, self.mat.s.nat + 1))),
            components=rows,
            force_mRy_Bohr=vectors.tolist(),
            converged=True,
        )

    def displacement_source(self, atom, direction, velocities):
        """At fixed native charge coefficients, differentiate the output map and free energy."""
        pc = self.pc
        v, parts = diagonal_strain_potential(
            pc, np.zeros((3, 3)), displacements=velocities
        )
        dv, band = self.valence_direction(v, velocities)
        ir = slice(pc.hartree.pw_offset, None)
        int_d = atomic_charge_integral_derivative(
            pc.input_density, self.mat.s, velocities
        )
        norm_geo = self.norm.explicit_source(
            ir_geometry=float(int_d @ self.norm.base[ir])
        )
        source = self.norm.operator.matvec(dv + pc.core.operator.matvec(v)) + norm_geo
        electro = self.den.boundary_gradient - 0.5 * self.den.pairing.gradient(
            self.den.density
        )
        xc = masked_potential_strain_source(
            pc.hartree.reanalysis,
            self.den.xc.field_response.ir.central_stars,
            parts["mask_derivative"],
        )
        den = electro @ parts["coulomb_geometry"] + pc.hartree.volume * (
            self.den.density[ir] @ (xc.view(float) - parts["xc_geometry"][ir])
        )
        scale = -1000.0 / self.mat.s.mult[atom - 1]
        explicit = scale * (den + self.core.potential_gradient @ v + band)
        row = dict(
            atom=int(atom),
            direction_cartesian=direction.tolist(),
            density_energy_force_mRy_per_Bohr=float(scale * den),
            occupied_core_force_mRy_per_Bohr=float(
                scale * (self.core.potential_gradient @ v)
            ),
            lattice_force_mRy_per_Bohr=float(scale * band),
            explicit_force_mRy_per_Bohr=float(explicit),
        )
        return source, row

    def run(self, atoms, *, rtol=1e-8, max_cycles=20, fixed_charge_only=False):
        closure = self.central_closure()
        self.progress(f"Central spinor valence charge relative error {closure:.6g}")
        directions = internal_modes(self.mat.s, atoms)
        sources, rows = [], []
        for a, u, vel in directions:
            self.progress(f"Computing analytic atomic source on site {a} along {u}")
            source, row = self.displacement_source(a, u, vel)
            sources.append(source)
            rows.append(row)
            self.progress(
                f"Atom {a} fixed-charge force {row['explicit_force_mRy_per_Bohr']:.9f} mRy/Bohr"
            )

        info = dict(
            source_method="analytic",
            case=self.case,
            components=rows,
            selected_atoms=sorted(set(atoms or range(1, self.mat.s.nat + 1))),
            structure=str(self.directory / (self.case + ".struct")),
            central_valence_relative_error=closure,
            central_relative_charge_residual=self.charge_residual,
            omitted_cross_window_pairs=self.parallel.sum(self.singular).tolist(),
            minimum_cross_window_gaps_Ry=[
                float(v) if np.isfinite(v) else None
                for v in self.parallel.reduce(self.minimum_gaps, "min")
            ],
        )
        if fixed_charge_only:
            return dict(
                info,
                observable="fixed-charge SOC Klein derivative",
                diagnostic_only=True,
            )
        pc = self.pc
        K, core, N = pc.potential_operator, pc.core.operator, self.norm.operator
        self.progress("Forming the SOC Klein-functional charge gradient")
        c = self.den.density_gradient + K.rmatvec(
            self.potential_readout() + self.core.potential_gradient
        )

        def jt(y):
            z = N.rmatvec(y)
            return K.rmatvec(self.potential_readout(z) + core.rmatvec(z))

        def forward(y):
            v = K.matvec(y)
            z, _ = self.valence_direction(v, None)
            return N.matvec(z + core.matvec(v))

        J = PairedOperator(self.n, self.n, forward, jt, "SOC charge Jacobian")
        solution, report = ScreenedResponse(J, progress=self.progress).solve(
            c, adjoint=True, rtol=rtol, max_cycles=max_cycles
        )
        if not report["converged"]:
            raise RuntimeError(f"SOC charge adjoint did not converge: {report}")
        forces = np.zeros((self.mat.s.nat, 3))
        for (a, u, _), source, row in zip(directions, sources, rows):
            induced = float(-1000.0 * (solution @ source) / self.mat.s.mult[a - 1])
            component = row["explicit_force_mRy_per_Bohr"] + induced
            row.update(
                screened_force_mRy_per_Bohr=induced, force_component_mRy_Bohr=component
            )
            forces[a - 1] += component * u
            self.progress(
                f"Atom {a}: charge response {induced:.9f}, screened total {component:.9f} mRy/Bohr"
            )
        return dict(
            info,
            converged=True,
            observable="screened frozen-Sigma atomic force",
            scope="complete nonmagnetic spinor force at fixed self-energy",
            force_mRy_Bohr=forces.tolist(),
            adjoint=report,
        )


def main():
    import argparse
    import json
    from dmft2_step import atomic_json

    a = argparse.ArgumentParser(description=__doc__)
    a.add_argument("--snapshot", type=Path, required=True)
    a.add_argument("--case", required=True)
    a.add_argument("--progress-file", type=Path)
    a.add_argument("--output", type=Path, required=True)
    a.add_argument("--library", type=Path)
    a.add_argument("--force-atoms", type=int, nargs="+")
    diagnostic = a.add_mutually_exclusive_group()
    diagnostic.add_argument(
        "--fixed-charge-only",
        action="store_true",
        help="diagnose the explicit derivative without solving an adjoint",
    )
    diagnostic.add_argument(
        "--soc-operator-only",
        action="store_true",
        help="evaluate only the frozen-radial operator term for ordinary FOR",
    )
    a.add_argument("--mpi", action="store_true")
    a.add_argument("--cache-mb", type=float, default=128.0)
    a.add_argument("--rtol", type=float, default=1e-8)
    args = a.parse_args()
    from response_parallel import ResponseParallel, limit_blas_threads

    limit_blas_threads()

    par = ResponseParallel.world() if args.mpi else SERIAL
    d = args.snapshot
    details = json.loads((d / "preparation.json").read_text())
    if not details.get("spinor") or details["case"] != args.case:
        raise ValueError("expected central spinor response preparation")

    def progress(line):
        if par.root:
            print(line, flush=True)
            if args.progress_file:
                from stress_response_summary import progress as record_progress

                record_progress(args.progress_file, line)

    with ExitStack() as stack:
        if "maximum_parent_normalization_change" in details:
            progress(
                "Maximum change from parent SOlapm.dat: "
                f"{details['maximum_parent_normalization_change']:.6g}"
            )
        model = SOCScreenedForce(
            d,
            args.case,
            stack,
            geometry=d / "soc.geometry",
            prefix=d / "response",
            matrices_prefix=d / "native",
            potential_core=d / "potential_core.npz",
            library=args.library,
            cache_mb=args.cache_mb,
            parallel=par,
            fixed_number=details["physical_policy"]["fixed_number"],
            operator_only=args.soc_operator_only,
            progress=progress,
        )
        report = (
            model.frozen_soc_force(args.force_atoms)
            if args.soc_operator_only
            else model.run(
                args.force_atoms,
                rtol=args.rtol,
                fixed_charge_only=args.fixed_charge_only,
            )
        )
        if par.root:
            args.output.mkdir(parents=True, exist_ok=True)
            atomic_json(args.output / "report.json", report)


if __name__ == "__main__":
    main()
