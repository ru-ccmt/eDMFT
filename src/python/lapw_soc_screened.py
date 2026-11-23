"""Matrix-free, fixed-self-energy charge and Klein derivatives for SOC.

All inputs are central native files. Electronic normalization and the fixed-N
chemical potential couple the k points. Charge and potential grids are shared
across ranks. Prepared LAPW/SOC data are reused within this geometry up to a per-rank
memory budget. Remaining k points are streamed in their original order.
"""

from contextlib import ExitStack
from pathlib import Path
import os
import shutil
import struct
import time
from tempfile import TemporaryFile, gettempdir
import tracemalloc
import numpy as np

from force_response_modes import internal_modes
from lapw_atomic_force import atomic_charge_integral_derivative
from lapw_atomic_geometry import atomic_geometry_readouts
from lapw_density_functional import PreparedDensityFunctional
from lapw_dft_vector import _record, _header
from lapw_functional_cancellation import CoreEigenvalueFunctional
from lapw_potential_core_state import fingerprint_files, load_potential_core_state
from lapw_radial_energy_policy import searched_energy_derivatives
from lapw_soc_electronic import prepare_electronic, operator_readout
from lapw_soc_projector import SpinorProjectors
from lapw_soc_response import SOCGeometry, SOCBandHamiltonian
from lapw_soc_subspace import spinor_coefficients
from lapw_soc_valence import SOCValenceAtK, native_band_geometry_source, native_spatial_geometry_source
from lapw_frequency_modes import prepare_frequency_modes
from lapw_strain_field import diagonal_strain_potential
from lapw_strain_potential import masked_potential_strain_source
from lapw_valence_state import _snapshot_inputs, prepare_radial_contexts
from response_parallel import SERIAL
from soc_scalar_window import matrices
from stress_response_adjoint import PairedOperator, ScreenedResponse, ResponseConvergenceError
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


def point_cache_limits(parallel, limit_mb):
    """Limit retained k data separately from resolvents and active work arrays.

    MPI shared-memory groups give the actual number of processes on this node.
    By default use at most 1 GiB per rank and one eighth of node memory overall.
    Matching blocks already use anonymous scratch files; keep room for the
    streamed point and other workers. A zero limit restores streaming.
    """
    if limit_mb is not None and (not np.isfinite(limit_mb) or limit_mb < 0):
        raise ValueError("k-point cache budget must be finite and nonnegative")
    local_size = 1
    if parallel.size > 1:
        from mpi4py import MPI
        shared = parallel.comm.Split_type(MPI.COMM_TYPE_SHARED)
        try:
            local_size = shared.Get_size()
        finally:
            shared.Free()
    if limit_mb is None:
        try:
            memory = os.sysconf("SC_PHYS_PAGES") * os.sysconf("SC_PAGE_SIZE")
        except (ValueError, OSError):
            memory = 0
        meminfo = Path("/proc/meminfo")
        if meminfo.is_file():
            for line in meminfo.read_text().splitlines():
                if line.startswith("MemAvailable:"):
                    memory = min(memory, int(line.split()[1]) * 1024)
                    break
        memory_limit = min(1024**3, max(0, memory) // (8 * local_size))
    else:
        memory_limit = int(limit_mb * 1024**2)
    scratch_limit = shutil.disk_usage(gettempdir()).free // (4 * local_size)
    return memory_limit, scratch_limit


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
        kpoint_cache_mb=None,
        fixed_number=True,
        operator_only=False,
        progress=lambda message: None,
    ):
        d = self.directory = Path(directory)
        self.case, self.parallel, self.progress = case, parallel, progress
        self.start_time = time.monotonic()
        self.library, self.fixed_number = library, fixed_number
        self.stack = stack
        self.cached_points = {}
        stack.callback(self.cached_points.clear)
        self.cache_bytes = self.cache_scratch = 0
        self.cache_limit, self.scratch_limit = (0, 0) if operator_only else point_cache_limits(parallel, kpoint_cache_mb)
        self.admit_points = self.cache_limit > 0
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
        # Keeping several k points must not multiply the live Fortran radial
        # handles. Both quadratures depend only on this geometry and potential.
        in1 = d / (case + (".in1c" if self.mat.complex else ".in1"))
        ltop = int(in1.read_text().splitlines()[1].split()[1])
        self.radial_contexts = prepare_radial_contexts(self.mat, ltop, "lapwso", stack, library)
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

    def _prepare_point(self, j, r, stack):
        nr = self.model.eps.shape[1]
        k, l, e, a = scalar_record(self.vector, self.offsets[r.index - 1])
        es, c, u, labels = r.arrays()
        if not np.array_equal(l, labels) or r.nlo != self.mat.nlo_total:
            raise ValueError(
                f"mismatched scalar and spinor LAPW rows at k={r.index}"
            )
        h, o, _ = matrices(self.matrices, r.index, k, l, self.mat.s, r.nlo)
        spin = spinor_coefficients(a, c[:, :nr])
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
            radial_contexts=self.radial_contexts,
        )
        current = item.project.evaluate(item.problem.waves)
        if not np.allclose(current, self.model.raw[j], rtol=2e-9, atol=2e-11):
            raise ValueError(
                f"assembled spinor projector differs from the native export at k={r.index}"
            )
        return item

    def points(self):
        for j, record in enumerate(self.records):
            if j in self.cached_points:
                yield j, self.cached_points[j]
                continue
            with ExitStack() as resources:
                if self.admit_points:
                    # Measure the allocations that survive construction, including
                    # NumPy buffers. RSS deltas undercount when malloc reuses a
                    # temporary. Shared native radial contexts were built above.
                    # Trace only during admission, never during a response action.
                    started_trace = not tracemalloc.is_tracing()
                    if started_trace:
                        tracemalloc.start()
                    try:
                        before = tracemalloc.get_traced_memory()[0]
                        item = self._prepare_point(j, record, resources)
                        used = tracemalloc.get_traced_memory()[0] - before
                    finally:
                        if started_trace:
                            tracemalloc.stop()
                    # Allow for allocator overhead. A nonpositive reading cannot
                    # certify a size (e.g. unrelated allocations were freed).
                    reserved = max(0, used) * 5 // 4 + 1024**2
                    scratch = item.val.matching_scratch_bytes
                    if (used > 0 and self.cache_bytes + reserved <= self.cache_limit
                            and self.cache_scratch + scratch <= self.scratch_limit):
                        self.stack.enter_context(resources.pop_all())
                        self.cached_points[j] = item
                        self.cache_bytes += reserved
                        self.cache_scratch += scratch
                    else:
                        # Keep the admitted prefix. Evicting on a sequential sweep
                        # would repeatedly rebuild every point on the next pass.
                        self.admit_points = False
                else:
                    item = self._prepare_point(j, record, resources)
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

    def valence_direction(self, v, velocity, *, hydrostatic=False, strain=None):
        dh = np.empty_like(self.hbar)
        du = np.empty_like(self.ubar)
        direct = np.zeros(self.n)
        for j, item in self.points():
            dh[j], du[j], part = item.action(v, velocity, hydrostatic=hydrostatic, strain=strain)
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

    def potential_readout(self, y=None, *, frequency_modes=None, geometry_file=None, with_strain=False):
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
                    g, fixed_number=self.fixed_number, return_raw=True, frequency_modes=frequency_modes
                )
            else:
                bh, bu = self.hbar, self.ubar
                y = np.zeros(self.n)
            if geometry_file is not None:
                geometry_file.seek(0)
            for j, item in self.points():
                if spool is not None:
                    waves = np.load(spool, allow_pickle=False)
                elif not has_density:
                    waves = np.zeros_like(item.problem.waves)
                if geometry_file is None:
                    result += item.pullback(bh[j], bu[j], waves, y)
                else:
                    out = item.pullback(bh[j], bu[j], waves, y, with_geometry=True, with_strain=with_strain)
                    result += out[0]
                    for covector in out[1:]:
                        np.save(geometry_file, covector, allow_pickle=False)
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
        den, norm_geo = self._density_geometry(velocities, parts)
        source = self.norm.operator.matvec(dv + pc.core.operator.matvec(v)) + norm_geo
        return source, self._force_row(atom, direction, v, den, band)

    def _density_geometry(self, velocities, parts):
        pc = self.pc
        ir = slice(pc.hartree.pw_offset, None)
        int_d = atomic_charge_integral_derivative(
            pc.input_density, self.mat.s, velocities
        )
        norm_geo = self.norm.explicit_source(
            ir_geometry=float(int_d @ self.norm.base[ir])
        )
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
        return float(den), norm_geo

    def _force_row(self, atom, direction, potential, den, band):
        scale = -1000.0 / self.mat.s.mult[atom - 1]
        core = self.core.potential_gradient @ potential
        return dict(atom=int(atom), direction_cartesian=direction.tolist(),
                    density_energy_force_mRy_per_Bohr=float(scale * den),
                    occupied_core_force_mRy_per_Bohr=float(scale * core),
                    lattice_force_mRy_per_Bohr=float(scale * band),
                    explicit_force_mRy_per_Bohr=float(scale * (den + core + band)))

    def pressure_source(self):
        return self.strain_source(np.eye(3))

    def strain_source(self, e):
        """Partial charge and Klein derivatives for one symmetric strain at fixed native density.

        q solves (I-J)q=s in this same density representation. Hence the
        total derivative is a+lambda.s with no additional input-charge
        transport. Normalization still sees the changing interstitial measure.
        """
        potential, density, parts, core, normalization = self._strain_fields(e)
        hydro = np.array_equal(e, np.eye(3))
        valence, bands = self.valence_direction(potential, None, hydrostatic=hydro,
                                                strain=None if hydro else e)
        source = self.norm.operator.matvec(valence + self.pc.core.operator.matvec(potential)) + normalization
        return source, self._strain_energy_parts(density, parts, core, bands)

    def _strain_fields(self, e):
        """Geometry-dependent fields in the native fixed charge coordinates."""
        from lapw_charge_response import native_mixer_inversion
        from lapw_strain_charge import charge_integration_geometry
        pc = self.pc
        potential, fields = diagonal_strain_potential(pc, e, transport_rate=0)
        density, parts = self.den.geometry(e, potential_parts=fields)
        core = float(self.core.potential_gradient @ potential)
        integral = charge_integration_geometry(pc.input_density, self.mat.s, e,
                                                inversion=native_mixer_inversion(self.mat.s))
        normalization = self.norm.explicit_source(ir_geometry=float(integral @ self.norm.base))
        return potential, density, parts, core, normalization

    @staticmethod
    def _strain_energy_parts(density, parts, core, bands):
        return dict(electrostatic_Ry=parts["electrostatic"], xc_Ry=parts["xc"],
                    core_Ry=core, bands_Ry=float(bands),
                    fixed_charge_derivative_Ry=float(density + core + bands))

    def _geometry_contractions(self, directions, strains, solution, fixed_file, response_file):
        """Reuse two exact electronic covectors across all strains and displacements.

        The expensive reverse propagation through scalar/spinor bands, the
        normalized projector and charge readout was done for c and for λ.
        A single streamed k sweep now pairs each explicit geometry derivative
        with both covectors. No additional frequency sum is needed per strain.
        """
        fixed_file.seek(0)
        response_file.seek(0)
        z = self.norm.operator.rmatvec(solution)
        forces = np.zeros((2, len(directions)))
        stress = np.zeros((2, len(strains)))
        for _, item in self.points():
            fb, fp = (np.load(fixed_file, allow_pickle=False) for _ in range(2))
            sb, sp = (np.load(response_file, allow_pickle=False) for _ in range(2))
            if len(strains):
                fh, fu = (np.load(fixed_file, allow_pickle=False) for _ in range(2))
                sh, su = (np.load(response_file, allow_pickle=False) for _ in range(2))
            if directions:
                forces += atomic_geometry_readouts(
                    item.val, self.directory, self.case, directions, fb, sb, z)
                # Spinor and correlated-projector translations have already
                # been differentiated in Cartesian site coordinates.
                for i, (_, _, pattern) in enumerate(directions):
                    forces[0, i] += np.vdot(fp, pattern).real
                    forces[1, i] += np.vdot(sp, pattern).real
            for i, e in enumerate(strains):
                band, _ = native_band_geometry_source(item.val, self.directory, self.case, e)
                density = native_spatial_geometry_source(item.val, self.directory, self.case, e)
                hydro = np.array_equal(e, np.eye(3))
                geometry = dict(hydrostatic=hydro, strain=None if hydro else e)
                dh = item.problem.soc.action(**geometry)
                du = item.project.evaluate(item.problem.waves, **geometry)
                stress[0, i] += fb @ band + np.vdot(fh, dh).real + np.vdot(fu, du).real
                stress[1, i] += sb @ band + np.vdot(sh, dh).real + np.vdot(su, du).real + z @ density
        return self.parallel.sum(forces), self.parallel.sum(stress)

    def _strain_rows(self, strains, solution, fixed_v, response_v, sums):
        partials, induced = [], []
        for i, e in enumerate(strains):
            potential, density, parts, core, norm = self._strain_fields(e)
            bands = fixed_v @ potential + sums[0, i]
            partials.append(self._strain_energy_parts(density, parts, core, bands))
            induced.append(float(response_v @ potential + sums[1, i] + solution @ norm))
        return partials, np.array(induced)

    def _force_rows(self, directions, solution, fixed_v, response_v, sums):
        zero = np.zeros((3, 3))
        rows = []
        forces = np.zeros((self.mat.s.nat, 3))
        for i, (a, u, velocity) in enumerate(directions):
            v, parts = diagonal_strain_potential(self.pc, zero, displacements=velocity)
            den, norm_geo = self._density_geometry(velocity, parts)
            row = self._force_row(a, u, v, den, fixed_v @ v + sums[0, i])
            induced = float(-1000.0 / self.mat.s.mult[a - 1]
                            * (response_v @ v + sums[1, i] + solution @ norm_geo))
            total = row["explicit_force_mRy_per_Bohr"] + induced
            row.update(screened_force_mRy_per_Bohr=induced, force_component_mRy_Bohr=total)
            rows.append(row)
            forces[a - 1] += total * u
            self.progress(f"Atom {a}: charge response {induced:.9f}, screened total {total:.9f} mRy/Bohr")
        return forces, rows

    def _screened_solution(self, fixed_file, response_file, resources, total_points,
                           rtol, max_cycles, frequency_rtol, with_strain=False):
        """Solve once for all directions and verify λ against the original discrete Jᵀ."""
        pc = self.pc
        K, core, N = pc.potential_operator, pc.core.operator, self.norm.operator
        self.progress("Forming the SOC Klein-functional charge gradient")
        geometry_options = dict(with_strain=True) if with_strain else {}
        fixed_v = self.potential_readout(geometry_file=fixed_file, **geometry_options)
        c = self.den.density_gradient + K.rmatvec(fixed_v + self.core.potential_gradient)

        def forward(y):
            v = K.matvec(y)
            z, _ = self.valence_direction(v, None)
            return N.matvec(z + core.matvec(v))

        def solve(tables, initial=None):
            def jt(y):
                z = N.rmatvec(y)
                return K.rmatvec(self.potential_readout(z, frequency_modes=tables) + core.rmatvec(z))
            jacobian = PairedOperator(self.n, self.n, forward, jt, "SOC charge Jacobian")
            return ScreenedResponse(jacobian, progress=self.progress).solve(
                c, adjoint=True, rtol=rtol, max_cycles=max_cycles, initial=initial)

        tables = {}
        if frequency_rtol:
            self.progress("Preparing signed Matsubara modes for the screened adjoint")
            tables = prepare_frequency_modes(self.model, frequency_rtol,
                                             self.scratch_limit - self.cache_scratch, resources)
        local_stats = [(self.records[k].index, len(t.signs), t.error) for k, t in tables.items()]
        stats = sorted(entry for rank in self.parallel.gather(local_stats) for entry in rank)
        mode_info = dict(
            factor_rtol=frequency_rtol, total_kpoints=total_points,
            compressed_kpoints=len(stats), kpoint_modes=[(int(i), r) for i, r, _ in stats],
            max_relative_factor_error=max((err for _, _, err in stats), default=0.0))
        if stats:
            ranks = [r for _, r, _ in stats]
            self.progress(f"Using {min(ranks)} to {max(ranks)} signed modes on {len(stats)}/{total_points} k points")
        # All ranks follow the same fallback decision, including those with
        # only full-frequency local points. Verify against the original Jᵀ.
        # During that check also collect the exact force covectors for λ.
        started = time.monotonic()
        attempts = []
        solution = None
        retry = False
        seed = initial = None
        if os.environ.get("EDMFT_RELAX_BLOCK"):
            from structure_relaxation import AdjointSeed
            seed = AdjointSeed(pc, "soc")
            initial = self.parallel.root_call(lambda: seed.load(os.environ.get("EDMFT_RELAX_INITIAL")))
            if initial is not None:
                self.progress("Starting GMRES from the preceding geometry adjoint")
        if stats:
            try:
                solution, report = solve(tables, initial=initial)
                attempts.append(report.copy())
            except ResponseConvergenceError as error:
                attempts.append(error.report.copy())
                retry = True
        if solution is None:
            if retry:
                self.progress("Mode solve did not converge; solving with all Matsubara frequencies")
            try:
                solution, report = solve(None, initial=initial)
            except ResponseConvergenceError:
                if initial is None:
                    raise
                self.progress("Restarting the full-frequency adjoint with a zero initial iterate")
                solution, report = solve(None)
            attempts.append(report.copy())

        def exact_readout(lam):
            z = N.rmatvec(lam)
            dv = self.potential_readout(z, geometry_file=response_file, **geometry_options) + core.rmatvec(z)
            residual = np.linalg.norm(c - (lam - K.rmatvec(dv)))
            return dv, float(self.parallel.reduce(float(residual), "max"))
        response_v = None
        checks = 0
        if response_file is not None or stats:
            response_v, residual = exact_readout(solution)
            checks += 1
            if stats and not retry and not (residual <= report["threshold"]):
                self.progress("Full-frequency residual exceeds target; refining with the full operator")
                solution, report = solve(None, initial=solution)
                attempts.append(report.copy())
                checks += 1
                response_v, residual = exact_readout(solution)
                retry = True
            report["residual_norm"] = residual
            report["converged"] = bool(np.isfinite(residual) and residual <= report["threshold"])
            report["full_frequency_check"] = True
            scale = np.linalg.norm(c) or 1.0
            self.progress(f"Full-frequency adjoint residual {residual / scale:.6e} "
                          f"(scaled target {report['threshold'] / scale:.6e})")
            if not report["converged"]:
                raise ResponseConvergenceError(report)
        report.update(full_frequency_fallback=retry, full_frequency_checks=checks,
                      actions=sum(a["actions"] for a in attempts) + checks,
                      seconds=time.monotonic() - started)
        if retry:
            report["mode_attempt"] = attempts[0]
        if seed is not None:
            self.parallel.root_call(lambda: seed.save(os.environ["EDMFT_RELAX_BLOCK"], solution))
        return solution, report, mode_info, fixed_v, response_v

    def run(self, atoms, *, rtol=1e-8, max_cycles=20, fixed_charge_only=False, pressure=False,
            frequency_rtol=None, tensor=False, atomic_forces=False):
        if frequency_rtol is None:
            frequency_rtol = min(1e-6, 0.01 * np.sqrt(rtol))
        if not np.isfinite(frequency_rtol) or frequency_rtol < 0 or frequency_rtol >= 1:
            raise ValueError("SOC frequency factor tolerance must be in [0,1)")
        closure = self.central_closure()
        self.progress(f"Central spinor valence charge relative error {closure:.6g}")
        cache = dict(
            retained=int(self.parallel.sum(len(self.cached_points))),
            total=int(self.parallel.sum(len(self.records))),
            max_rank_reserved_MiB=self.parallel.reduce(self.cache_bytes, "max") / 1024**2,
            max_rank_limit_MiB=self.parallel.reduce(self.cache_limit, "max") / 1024**2,
            max_rank_scratch_MiB=self.parallel.reduce(self.cache_scratch, "max") / 1024**2,
        )
        self.progress(f"Retaining {cache['retained']}/{cache['total']} prepared SOC k points; "
                      f"maximum reserved memory {cache['max_rank_reserved_MiB']:.1f} MiB/rank, "
                      f"matching scratch {cache['max_rank_scratch_MiB']:.1f} MiB/rank")
        directions = [] if pressure or (tensor and not atomic_forces) else internal_modes(self.mat.s, atoms)
        rows = []
        if fixed_charge_only:
            for a, u, vel in directions:
                _, row = self.displacement_source(a, u, vel)
                rows.append(row)

        info = dict(
            source_method="analytic",
            case=self.case,
            components=rows,
            kpoint_cache=cache,
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
        strains = []
        if pressure or tensor:
            from lapw_total_stress import invariant_strains
            strains = invariant_strains(self.mat.s) if tensor else np.eye(3)[None, :, :]
            if not np.isclose(self.geo.volume, self.pc.hartree.volume, rtol=1e-10, atol=0):
                raise ValueError("LAPWSO and LAPW0 primitive volumes disagree")
            info.update(primitive_volume_Bohr3=float(self.geo.volume),
                        ensemble="fixed_N" if self.fixed_number else "fixed_mu")
            if not atomic_forces:
                del info["selected_atoms"], info["components"]
        # For the hydrostatic-only case one forward source needs no covector
        # files. Multiple strains and joint stress/force readouts amortize
        # the exact transpose and geometry sweep over all the directions.
        share_strain = tensor and (len(strains) > 1 or bool(directions))
        strain_sources, partials = [], []
        if not share_strain:
            for i, e in enumerate(strains):
                self.progress(f"Preparing SOC strain direction {i+1}/{len(strains)}")
                source, partial = self.strain_source(e)
                strain_sources.append(source)
                partials.append(partial)
            if pressure:
                info.update(partials[0])
        if fixed_charge_only:
            return dict(
                info,
                observable="fixed-charge SOC Klein derivative",
                diagnostic_only=True,
            )
        with ExitStack() as resources:
            save_geometry = bool(directions) or share_strain
            fixed_file = resources.enter_context(TemporaryFile()) if save_geometry else None
            response_file = resources.enter_context(TemporaryFile()) if save_geometry else None
            solution, report, mode_info, fixed_v, response_v = self._screened_solution(
                fixed_file, response_file, resources, cache["total"], rtol, max_cycles, frequency_rtol,
                with_strain=share_strain)
            info["frequency_compression"] = mode_info
            if save_geometry:
                self.progress(f"Contracting the shared adjoint with {len(directions)} displacements "
                              f"and {len(strains) if share_strain else 0} strain directions")
                forces_sum, stress_sum = self._geometry_contractions(
                    directions, strains if share_strain else [], solution, fixed_file, response_file)
                if share_strain:
                    partials, induced = self._strain_rows(strains, solution, fixed_v, response_v, stress_sum)
            if directions:
                forces, info["components"] = self._force_rows(
                    directions, solution, fixed_v, response_v, forces_sum)
            elif not pressure:
                forces = np.zeros((self.mat.s.nat, 3))
        if tensor:
            from lapw_total_stress import stress_from_contractions
            from stress_response_adjoint import RY_BOHR3_TO_GPA

            explicit = np.array([p["fixed_charge_derivative_Ry"] for p in partials])
            if not share_strain:
                induced = np.array([solution @ source for source in strain_sources])
            values = explicit + induced
            sigma = stress_from_contractions(strains, values, self.geo.volume)
            info["central"] = dict(kpoints=cache["total"], retained_bands=self.model.eps.shape[1],
                                   fixed_number=self.fixed_number, dmft_spin_channels=1,
                                   normalization_enabled=self.norm.enabled, common_mu_Ry=self.model.meta.mu,
                                   mpi_ranks=self.parallel.size, central_valence_relative_error=closure)
            info["charge_fixed_point"] = dict(warning_threshold=1e-5, within_threshold=True)
            info["timings"] = dict(adjoint_solve_seconds=report["seconds"],
                                   total_seconds=self.parallel.reduce(time.monotonic() - self.start_time, "max"))
            if atomic_forces:
                info["atomic_forces"] = dict(case=self.case, structure=info["structure"],
                    selected_atoms=info.pop("selected_atoms"), components=info.pop("components"),
                    force_mRy_Bohr=forces.tolist(), source_method="analytic",
                    scope="screened frozen-Sigma SOC atomic forces")
            mean_stress = float(np.trace(sigma) * RY_BOHR3_TO_GPA / 3)
            self.progress(f"Mean tensile SOC stress {mean_stress:.9f} GPa")
            return dict(info, converged=True, observable="screened frozen-Sigma SOC stress",
                        stress_Ry_Bohr3=sigma.tolist(), stress_GPa=(sigma*RY_BOHR3_TO_GPA).tolist(),
                        compression_pressure_GPa=-mean_stress, mean_tensile_stress_GPa=mean_stress,
                        strains=strains.tolist(), derivatives_Ry=values.tolist(),
                        explicit_Ry=explicit.tolist(), induced_Ry=induced.tolist(),
                        strain_explicit_contributions=partials, shared_strain_readout=bool(share_strain), adjoint=report,
                        full_tensor_within_symmetry=True, stress_ready=True,
                        scope="nonmagnetic SOC, fixed correlated axes and imposed crystal symmetry")
        if pressure:
            from stress_response_adjoint import RY_BOHR3_TO_GPA
            induced = float(solution @ strain_sources[0])
            derivative = partials[0]["fixed_charge_derivative_Ry"] + induced
            pressure_gpa = -RY_BOHR3_TO_GPA * derivative / (3 * self.geo.volume)
            self.progress(f"Screened SOC pressure {pressure_gpa:.9f} GPa (compression positive)")
            return dict(info, converged=True, observable="screened hydrostatic SOC pressure",
                        charge_response_derivative_Ry=induced, dF_dlog_linear_scale_Ry=derivative,
                        compression_pressure_GPa=float(pressure_gpa), adjoint=report)
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
    quantity = a.add_mutually_exclusive_group()
    a.add_argument("--force-atoms", type=int, nargs="+")
    quantity.add_argument("--stress", action="store_true", help="complete symmetry-allowed SOC tensor")
    a.add_argument("--atomic-forces", action="store_true", help="also compute forces with the stress adjoint")
    quantity.add_argument("--pressure", action="store_true",
                          help="compute -dF/dV for uniform dilation of the input cell")
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
    a.add_argument("--kpoint-cache-mb", type=float, default=None,
                   help="MiB per MPI rank for prepared SOC k points; default automatic, 0 streams")
    a.add_argument("--frequency-rtol", type=float, default=None,
                   help="relative error in each Matsubara factor; zero uses full frequencies")
    a.add_argument("--rtol", type=float, default=1e-8)
    args = a.parse_args()
    if args.atomic_forces and not args.stress:
        a.error("--atomic-forces accompanies --stress")
    if args.stress and (args.soc_operator_only or args.fixed_charge_only):
        a.error("force-only diagnostics require the atomic force path")
    if args.pressure and args.soc_operator_only:
        a.error("--soc-operator-only applies to atomic forces")
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
            kpoint_cache_mb=args.kpoint_cache_mb,
            parallel=par,
            fixed_number=details.get("response_policy", details["physical_policy"])["fixed_number"],
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
                pressure=args.pressure, frequency_rtol=args.frequency_rtol,
                tensor=args.stress, atomic_forces=args.atomic_forces,
            )
        )
        if par.root:
            args.output.mkdir(parents=True, exist_ok=True)
            atomic_json(args.output / "report.json", report)


if __name__ == "__main__":
    main()
