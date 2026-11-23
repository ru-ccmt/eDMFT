"""Opt-in cell and internal-coordinate relaxation for run_dmft.py.

MSR1a updates charge and internal coordinates at fixed cell. On a converged
charge, a fixed-position STRFT measurement supplies the tensile stress and
screened forces. Only the outer loop writes lattice constants. Force units are
mRy/Bohr, stress and the trial elastic moduli are in GPa, lengths in Bohr.
"""

import fcntl
import hashlib
import itertools
import json
import os
from contextlib import contextmanager
from datetime import datetime
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import zipfile

import numpy as np

from dmft2_step import atomic_json
from lapw_lattice import NativeLattice
from wstruct import Struct


DEFAULTS = {
    "relax": False,
    "relax_pressure": 0.0,  # positive external compression, GPa
    "relax_stress_tol": 0.02,  # largest absolute residual principal stress, GPa
    "relax_force_tol": 0.5,  # largest allowed per-atom force norm, mRy/Bohr
    "relax_max_strain": 0.01,  # largest absolute principal log strain per cell step
    "relax_bulk_modulus": 120.0,  # trial isotropic bulk and shear stiffnesses, GPa
    "relax_shear_modulus": 60.0,
    "relax_refresh_bohr": 0.08,  # refresh the force correction after this motion
    "relax_max_steps": 20,  # screened evaluations per invocation
}


def structure(path):
    s = Struct(inAngs=False)
    s.ReadStruct(str(path))
    return s


def options(params):
    opts = {k: params[k] for k in DEFAULTS}
    if (
        type(opts["relax"]) is not bool
        or type(opts["relax_max_steps"]) is not int
        or opts["relax_max_steps"] < 1
    ):
        raise ValueError("relax must be boolean and relax_max_steps a positive integer")
    for key in DEFAULTS.keys() - {"relax", "relax_max_steps"}:
        x = opts[key]
        if isinstance(x, bool) or not isinstance(x, (float, int)) or not np.isfinite(x):
            raise ValueError(f"{key} requires a finite number")
        if key != "relax_pressure" and x <= 0:
            raise ValueError(f"{key} must be positive")
    return opts


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def atomic_write(path, data):
    path = Path(path)
    with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as f:
        tmp = Path(f.name)
    try:
        tmp.write_bytes(data)
        if path.exists():
            shutil.copymode(path, tmp)
        tmp.replace(path)
    finally:
        tmp.unlink(missing_ok=True)


def log(directory, text):
    text = f"[{datetime.now().isoformat(timespec='seconds')}] {text}"
    print("RELAX", text, flush=True)
    with (directory / "history").open("a") as f:
        print(text, file=f)


def compact_number(value, format_spec=".3g"):
    return format(value, format_spec).replace("e-0", "e-").replace("e+0", "e+")


def reset_mixer_history(case, directory, reason):
    """Restart only this case's WIEN multisecant mixer; retain the charge."""
    files = [p for p in Path.cwd().glob(case + ".broyd*") if p.is_file() or p.is_symlink()]
    for file in files:
        file.unlink()
    log(directory, f"{case}.broyd* reset ({len(files)} files): {reason}")


def set_mode(path, value):
    """Change the fixed-width first field; keep the rest of the input verbatim."""
    path = Path(path).resolve()
    data = path.read_bytes()
    m = re.match(rb"(\S+)([ \t]*)(?=\S|\r?\n|$)", data)
    if m is None:
        raise ValueError(f"missing first input field in {path}")
    token = value.encode()
    if m[1] != token:
        width = max(len(m[0]), len(token) + bool(m[2]))
        atomic_write(path, token.ljust(width) + data[m.end() :])


def restore_headers(directory):
    journal = directory / "input_headers.json"
    if journal.exists():
        for filename, first in json.loads(journal.read_text()).items():
            path = Path(filename)
            _, newline, tail = path.read_bytes().partition(b"\n")
            atomic_write(path, first.encode() + newline + tail)
        journal.unlink()


@contextmanager
def managed_headers(case, directory):
    restore_headers(directory)
    modes = {case + ".in0": "STR", case + ".inm": "MSR1"}
    for suffix in (".in2", ".in2c"):
        p = Path(case + suffix)
        if p.is_file() and p.stat().st_size:
            modes[str(p)] = "FOR"
    if len(modes) < 3:
        raise FileNotFoundError("relaxation requires an active case.in2 or case.in2c")
    files = {str(Path(p).resolve()): v for p, v in modes.items()}
    for p in files:
        if Path(p).parent != Path.cwd():
            raise ValueError(f"make a case-local copy of the managed input {p}")
    atomic_json(
        directory / "input_headers.json",
        {p: Path(p).read_bytes().split(b"\n", 1)[0].decode() for p in files},
    )
    try:
        for p, mode in files.items():
            set_mode(p, mode)
        yield
    finally:
        restore_headers(directory)


def agree_on_positions(path):
    """Give all native codes the same geometry after an MSR1a move.

    LAPW1 reads the optional 15-decimal tail; DMFT2 and parts of the response
    read the 8-decimal site lines. The resulting tiny phase differences matter
    in the cancellation of muffin-tin and interstitial matrix elements. Use
    the common precision throughout, just as in the finite-displacement tests.
    """
    path = Path(path)
    lines = path.read_bytes().splitlines(keepends=True)
    i = next((i for i, line in enumerate(lines) if line.strip() == b"Precise positions"), None)
    if i is None:
        return
    pos = structure(path).pos
    tail = lines[i + 1 : i + 1 + len(pos)]
    precise = np.asarray([list(map(float, line.split())) for line in tail])
    if precise.shape != pos.shape or not np.allclose(pos, precise, atol=1.01e-8, rtol=0):
        raise ValueError(
            "inconsistent site and Precise positions records; reconcile the STRUCT geometry"
        )
    if np.array_equal(pos, precise):
        return
    for j, row in enumerate(pos):
        lines[i + 1 + j] = ("".join(f"{x:20.15f}" for x in row) + "\n").encode()
    atomic_write(path, b"".join(lines))


def native_forces(text, s):
    """Sum the same site-local parts as MIXER, including any explicit SOC FVA."""
    from force_response_modes import internal_modes

    force = np.zeros((s.nat, 3))
    seen = [set() for _ in range(s.nat)]
    for line in text.splitlines():
        m = re.match(r"^:(FHF|FVA|FSU|FCO)(\d{3}):", line)
        if m:
            i = int(m[2]) - 1
            numbers = line[17:].split()
            if not 0 <= i < s.nat or len(numbers) < 4:
                raise ValueError("invalid native force row: " + line)
            force[i] += [float(v.replace("D", "E")) for v in numbers[1:4]]
            seen[i].add(m[1])
    projected = np.zeros_like(force)
    for atom, direction, _ in internal_modes(s):
        i = atom - 1
        if not {"FHF", "FVA", "FCO"} <= seen[i]:
            raise ValueError(f"incomplete current force on atom {atom}: {seen[i]}")
        # The response uses global Cartesian forces. Native FGL may use a
        # different crystal frame in an oblique cell; rotate the local components.
        projected[i] += direction * np.dot(direction, np.linalg.solve(s.rotloc[i], force[i]))
    if not np.isfinite(projected).all():
        raise ValueError("nonfinite forces")
    return projected


def correction_rows(delta, s):
    rows = []
    for i, f in enumerate(np.einsum("aij,aj->ai", s.rotloc, delta), 1):
        if np.any(f):
            rows.append(
                f":FVA{i:03d}: {i:3d}.ATOM"
                + "".join(f"{v:18.9f}" for v in (np.linalg.norm(f), *f))
                + " previous screened minus ordinary force\n"
            )
    return "".join(rows)


def stress_residual(s, tensor, pressure):
    stress = np.asarray(tensor, float)
    if (
        stress.shape != (3, 3)
        or not np.isfinite(stress).all()
        or not np.allclose(stress, stress.T, atol=1e-8)
    ):
        raise ValueError("invalid screened stress tensor")
    rotations = NativeLattice.from_structure(s).cartesian_rotations(s.timat)
    residual = stress + pressure * np.eye(3)
    return np.mean([r @ residual @ r.T for r in rotations], axis=0)


def elastic_compliance(opts):
    """Isotropic compliance acting on flattened Cartesian symmetric tensors."""
    hydro = np.outer(np.eye(3).ravel(), np.eye(3).ravel()) / 3
    return hydro / (3 * opts["relax_bulk_modulus"]) + (np.eye(9) - hydro) / (
        2 * opts["relax_shear_modulus"]
    )


def updated_compliance(cell, residual, opts, prior):
    """Inverse BFGS from two measured stresses in the same Cartesian frame.

    Keep all nine tensor entries so the Euclidean dot product is the physical
    double contraction (including both shear entries). Only symmetric,
    crystal-invariant tensors are ever applied. WIEN can rotate the Cartesian
    frame when an oblique metric changes; transport the previous tensors
    using the polar rotation of the actual cell change.
    """
    base = elastic_compliance(opts)
    if (
        not prior
        or prior["controls"] != [opts[k] for k in ELASTIC_CONTROLS]
        or not np.allclose(cell, prior["next_cell"], atol=1e-9, rtol=0)
    ):
        return base, False
    f = cell @ np.linalg.inv(prior["cell"])
    u, stretch, vt = np.linalg.svd(f)
    rotation = u @ vt
    transform = np.kron(rotation, rotation)
    compliance = transform @ prior["compliance"] @ transform.T
    step = ((u * np.log(stretch)) @ u.T).ravel()
    change = residual.ravel() - transform @ np.ravel(prior["stress"])
    curvature = step @ change
    if curvature <= max(1e-12, 0.05 * np.linalg.norm(step) * np.linalg.norm(change)):
        return base, False
    update = np.eye(9) - np.outer(step, change) / curvature
    trial = update @ compliance @ update.T + np.outer(step, step) / curvature
    # A noisy stress difference can imply essentially zero elastic stiffness.
    # Discard such an update; the independently bounded step is always retained.
    factor = np.linalg.cholesky(base)
    scaled = np.linalg.solve(factor, np.linalg.solve(factor, trial).T)
    eig = np.linalg.eigvalsh(scaled)
    if eig[0] < 0.05 or eig[-1] > 20:
        return base, False
    return (trial + trial.T) * 0.5, True


ELASTIC_CONTROLS = ("relax_pressure", "relax_bulk_modulus", "relax_shear_modulus")


def cell_step(residual, opts, compliance=None):
    if compliance is None:
        compliance = elastic_compliance(opts)
    strain = -(compliance @ residual.ravel()).reshape(3, 3)
    strain = (strain + strain.T) * 0.5
    size = max(abs(np.linalg.eigvalsh(strain)))
    return strain * min(1.0, opts["relax_max_strain"] / size) if size else strain


def strained_struct(source, destination, strain):
    """Write the new metric in WIEN's conventional parameters, including R cells.

    LATGEN's R coordinates are primitive rhombohedral, but the STRUCT fields
    are hexagonal lengths. Verify the metric after native fixed-width rounding.
    Keep all atom, local-frame, optional precise-position and symmetry records.
    """
    s = structure(source)
    lattice = NativeLattice.from_structure(s)
    val, vec = np.linalg.eigh(strain)
    cell = ((vec * np.exp(val)) @ vec.T) @ lattice.cell
    old = np.array([s.a, s.b, s.c, s.alpha, s.beta, s.gamma])
    kind = s.lattyp.strip()
    if kind in ("R", "H"):
        lengths = old[:3] * np.linalg.norm(cell, axis=1) / np.linalg.norm(lattice.cell, axis=1)
        params = [*lengths, *old[3:]]
    else:
        lengths = np.linalg.norm(cell, axis=0)
        directions = cell / lengths
        angle = np.rad2deg(np.arccos(np.clip(directions.T @ directions, -1, 1)))
        params = [*lengths, angle[1, 2], angle[0, 2], angle[0, 1]]
    line = "".join(f"{x:10.6f}" for x in params)
    if len(line) != 60:
        raise ValueError("cell exceeds the WIEN STRUCT field widths")
    lines = Path(source).read_bytes().splitlines(keepends=True)
    lines[3] = line.encode() + lines[3][60:]
    Path(destination).write_bytes(b"".join(lines))
    next_s = structure(destination)
    new = NativeLattice.from_structure(next_s)
    if not np.allclose(new.cell.T @ new.cell, cell.T @ cell, rtol=4e-7, atol=3e-6):
        raise ValueError(
            "proposed strain cannot be represented in the current WIEN lattice setting"
        )
    new.cartesian_rotations(next_s.timat)
    return next_s


def spheres_clear(s):
    lattice = NativeLattice.from_structure(s)
    inv = np.linalg.inv(lattice.primitive_cell)
    frac = (inv @ lattice.cell @ s.pos.T).T
    radii = np.repeat(s.rmt, s.mult)
    for i in range(len(frac)):
        for j in range(i + 1):
            cutoff = radii[i] + radii[j]
            dif = frac[i] - frac[j]
            span = cutoff * np.linalg.norm(inv, axis=1) + 1e-9
            ranges = [
                range(int(np.ceil(d - b)), int(np.floor(d + b)) + 1) for d, b in zip(dif, span)
            ]
            for t in itertools.product(*ranges):
                if i == j and t == (0, 0, 0):
                    continue
                if np.linalg.norm(lattice.primitive_cell @ (dif - t)) < cutoff - 1e-8:
                    return False
    return True


def density_seed(struct_file, next_struct, source, output):
    """Keep atom-centred radial coefficients. Correct the interstitial constant
    so the restart has the same electron count in the new native quadrature."""
    from lapw_charge_response import native_charge_covectors, native_mixer_inversion
    from wien_response_density import Density

    old, new = structure(struct_file), structure(next_struct)
    d = Density(source, struct_file)
    v = d.vector()

    def weights(s):
        return native_charge_covectors(d, s, inversion=native_mixer_inversion(s))

    mt, ir = weights(old)
    target = (mt + ir) @ v
    mt, ir = weights(new)
    offset = len(v) - 2 * len(d.pw) + 2 * list(d.pw).index((0, 0, 0))
    if ir[offset] <= 0:
        raise ValueError("nonpositive interstitial volume")
    d.pw[0, 0, 0] += (target - (mt + ir) @ v) / ir[offset]
    # CLMSUM's second record includes the conventional lattice constants.
    line = d.lines[1]
    old_cell = "".join(f"{x:10.6f}" for x in (old.a, old.b, old.c))
    new_cell = "".join(f"{x:10.6f}" for x in (new.a, new.b, new.c))
    if old_cell in line:
        d.lines[1] = line.replace(old_cell, new_cell, 1)
    d.write(output)
    stored = Density(output, next_struct)
    new_charge = float((mt + ir) @ stored.vector())
    if not np.isclose(new_charge, target, atol=5e-8, rtol=0):
        raise ValueError("transported charge differs after writing the native density")
    return float(target), new_charge


def read_json(path, default=None):
    return json.loads(path.read_text()) if path.is_file() else default


class AdjointSeed:
    """Transport the starting iterate only. GMRES still certifies the new solution.

    Muffin-tin coordinates require the same site order, local harmonics, axes
    and radial grids. Native plane waves are keyed by their integer G labels:
    their ordering and number may change as the lattice constants change.
    """

    def __init__(self, pc, backend):
        d, s = pc.input_density, pc.structure
        layout = [backend, s.lattyp.strip(), np.asarray(list(d.rad), int).tolist()]
        layout += [
            np.asarray(getattr(s, key)).tolist()
            for key in ("mult", "Znuc", "jrj", "r0", "rmt", "rotloc", "timat", "tau")
        ]
        self.signature = hashlib.sha256(json.dumps(layout).encode()).hexdigest()
        self.g = np.asarray(list(d.pw), dtype=int).reshape(-1, 3)
        self.radial = sum(len(v) for v in d.rad.values())
        self.size = self.radial + 2 * len(self.g)

    def load(self, path):
        if not path:
            return None
        try:
            with np.load(path, allow_pickle=False) as z:
                if str(z["layout"]) != self.signature:
                    return None
                old, labels = z["adjoint"], z["g"]
                if (
                    labels.ndim != 2
                    or labels.shape[1] != 3
                    or old.shape != (self.radial + 2 * len(labels),)
                    or not np.isfinite(old).all()
                ):
                    return None
                result = np.zeros(self.size)
                result[: self.radial] = old[: self.radial]
                coeffs = dict(zip(map(tuple, labels), old[self.radial :].reshape(-1, 2)))
                for i, g in enumerate(self.g):
                    if tuple(g) in coeffs:
                        result[self.radial + 2 * i : self.radial + 2 * i + 2] = coeffs[tuple(g)]
                return result
        except (OSError, ValueError, KeyError, EOFError, zipfile.BadZipFile):
            return None

    def save(self, directory, value):
        if value.shape != (self.size,) or not np.isfinite(value).all():
            return
        path = Path(directory) / "adjoint_seed.npz"
        with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as f:
            tmp = Path(f.name)
        try:
            with tmp.open("wb") as f:
                np.savez_compressed(f, layout=self.signature, g=self.g, adjoint=value)
            tmp.replace(path)
        finally:
            tmp.unlink(missing_ok=True)


def refresh_projector(directory):
    """Remove an older dmft1 normalization from the native dmft2 input path.

    FOR already traverses the overlap sum to differentiate X. Leaving the
    cache in place makes it use the previous X while STRFT reconstructs X
    from the current eigenvectors. Preserve the input before the native pass.
    """
    cache = Path("SOlapm.dat")
    if not cache.is_file():
        return
    backup = Path(directory) / ("SOlapm.input." + digest(cache) + ".dat")
    if cache.is_symlink():
        atomic_write(backup, cache.read_bytes())
        cache.unlink()
    else:
        cache.replace(backup)
    log(
        Path(directory).parent,
        f"{Path(directory).name} saved {backup.name}; "
        "dmft2 will normalize projectors on the current LAPW states",
    )


def save_measured_normalization(directory, export, channels):
    """Stage the measured X in the same fixed-width impurity frame as dmft1."""
    from frozen_sigma_response import Metadata

    metadata = Metadata.read(str(export / "response.meta"))
    if len(channels) > 1:
        other = Metadata.read(str(export / "response_dn.meta"))
        if not (
            np.array_equal(metadata.dimensions, other.dimensions)
            and np.allclose(metadata.normalization, other.normalization, rtol=0, atol=1e-11)
        ):
            # A single cache cannot represent both spins. In its absence each
            # ordinary native call computes its own X, just as in the measurement.
            return
    rows = []
    for ic, dim in enumerate(metadata.dimensions):
        rows.append(f"{ic+1:12d}\n")
        for row in metadata.normalization[:dim, :dim, ic]:
            if not np.isfinite(row).all():
                raise ValueError("nonfinite projector normalization in response snapshot")
            rows.append("".join(f"{z.real:20.16f}{z.imag:20.16f}" for z in row) + "\n")
    atomic_write(Path(directory) / "SOlapm.measured.dat", "".join(rows).encode())


class ElectronicBlock:
    """Hooks used by one ordinary run_dmft electronic/impurity loop."""

    def __init__(self, case, directory, opts):
        self.case = case
        self.directory = directory
        self.opts = opts
        self.path = Path(case + ".struct")
        agree_on_positions(self.path)
        self.frozen = False
        self.inner_iteration = 0
        self.small_force_count = 0
        self.charge_met = False
        self.last_mode = None
        self.predictor = read_json(directory.parent / "predictor.json")
        self.ordinary = None

    def mix(self, offset, call_mixer):
        s = structure(self.path)
        scf = Path(self.case + ".scf")
        with scf.open("rb") as stream:
            stream.seek(offset)
            current = native_forces(stream.read().decode(), s)
        refresh = False
        if self.predictor:
            p = self.predictor
            old = structure(p["structure"])
            if not (
                old.lattyp == s.lattyp
                and np.allclose(old.rmt, s.rmt, rtol=0, atol=1e-10)
                and np.allclose(old.rotloc, s.rotloc, rtol=0, atol=1e-7)
                and np.array_equal(old.mult, s.mult)
                and np.array_equal(old.timat, s.timat)
                and np.array_equal(old.Znuc, s.Znuc)
            ):
                raise ValueError(
                    "force predictor does not match the current site and symmetry layout"
                )
            lattice = NativeLattice.from_structure(s)
            frac = (s.pos - old.pos) @ lattice.cell.T @ np.linalg.inv(lattice.primitive_cell.T)
            movement = np.linalg.norm((frac - np.rint(frac)) @ lattice.primitive_cell.T, axis=1)
            fmat = lattice.cell @ np.linalg.inv(NativeLattice.from_structure(old).cell)
            u, stretch, vt = np.linalg.svd(fmat)
            refresh = (
                max(movement) >= self.opts["relax_refresh_bohr"]
                or max(abs(np.log(stretch))) > 1.1 * self.opts["relax_max_strain"]
            )
            if not refresh:
                delta = np.asarray(p["correction"]) @ (u @ vt).T
                current += delta
                with scf.open("a") as f:
                    f.write(correction_rows(delta, s))
        peak = float(max(np.linalg.norm(current, axis=1)))
        threshold = self.opts["relax_force_tol"] * (1.0 if self.predictor else 1.5)
        # Hold the geometry until charge equilibration once an ordinary force
        # is small. Unconverged charge oscillations otherwise switch MSR1a
        # on and off, repeatedly discarding useful native mixing history.
        self.small_force_count = self.small_force_count + 1 if peak <= threshold else 0
        self.frozen = (
            self.frozen
            or refresh
            or (self.charge_met and self.small_force_count >= 2)
            or not np.any(current)
        )
        self.hold_for_refresh = refresh
        self.peak = peak
        mode = "MSR1" if self.frozen else "MSR1a"
        if self.last_mode and self.last_mode != mode:
            reset_mixer_history(
                self.case,
                self.directory.parent,
                f"{self.directory.name} mixer changes from {self.last_mode} to {mode}",
            )
        self.last_mode = mode
        self.inner_iteration += 1
        self.force_kind = "FOR+scrn" if self.predictor and not refresh else "FOR"
        set_mode(self.case + ".inm", mode)
        try:
            call_mixer()
            if mode == "MSR1a":
                agree_on_positions(self.path)
        finally:
            set_mode(self.case + ".inm", "MSR1")

    def check_ordinary(self, drho, dene, p):
        """Release an early position hold if the converged charge changed the force."""
        step = f"{self.directory.name}.{self.inner_iteration:03d}"
        # Diff uses 1.0 until two SCF records exist. Read the first residual from MIXER.
        charge = drho
        if self.inner_iteration == 1:
            rows = [
                line
                for line in Path(self.case + ".scfm").read_text().splitlines()
                if line.startswith(":DIS")
            ]
            charge = float(rows[-1].split()[-1]) if rows else None
        self.charge_met = charge is not None and 0 <= charge <= p["cc"]
        charge_text = (
            "rho n/a"
            if charge is None
            else f"rho {compact_number(charge, '.2e'):<9} (<={compact_number(p['cc'])},"
            f"{'yes' if self.charge_met else 'no'})"
        )
        position_mode = "(posi.-held)" if self.frozen else "(int.-relax)"
        log(
            self.directory.parent,
            f"{step} {self.last_mode:<5} {position_mode} "
            f"{self.force_kind:<8} {self.peak:9.5f} mRy/Bohr, {charge_text}",
        )
        if (
            self.predictor
            and self.frozen
            and not self.hold_for_refresh
            and drho <= p["cc"]
            and dene <= p["ec"]
            and self.peak > self.opts["relax_force_tol"]
        ):
            self.frozen = False
            self.small_force_count = 0

    def before_screen(self, drho, dene, occupation_change, p, native_pass):
        if not (
            self.frozen
            and 0 <= drho <= p["cc"]
            and 0 <= dene <= p["ec"]
            and (not p["runIMP"] or occupation_change <= p["nc"])
        ):
            raise RuntimeError(
                "inner relaxation exhausted finish/max_lda_iterations before force "
                "and electronic convergence; resume from the current charge and positions"
            )
        # Measure FOR on exactly the same final LAPW state as STRFT. This also
        # includes the explicit SOC operator if enabled in the ordinary force path.
        log(
            self.directory.parent,
            f"{self.directory.name} electronic convergence reached; "
            f"measuring FOR and STRFT at identical positions. Screened solve: {self.case}.outputdmf2_response",
        )
        self.geometry_hash = digest(self.path)
        path = self.directory / "ordinary.scf"
        native_pass(str(path))
        self.ordinary = native_forces(path.read_text(), structure(self.path))
        for file in (Path(self.case + ".in2"), Path(self.case + ".in2c")):
            if file.is_file() and file.stat().st_size:
                set_mode(file, "STRFT")
        self.previous_id = (
            read_json(Path(self.case + ".stress_response.status.json"), {}) or {}
        ).get("run_id")

    def after_screen(self, p):
        status = read_json(Path(self.case + ".stress_response.status.json"), {})
        if (
            status.get("run_id") == self.previous_id
            or status.get("requested_mode") != "STRFT"
            or status.get("status") != "complete"
            or not status.get("stress_applied")
            or not status.get("force_applied")
            or digest(self.path) != self.geometry_hash
        ):
            raise RuntimeError("missing fresh published fixed-geometry STRFT measurement")
        sf, ff = Path(status["report"]), Path(status["force_report"])
        if digest(sf) != status["report_sha256"] or digest(ff) != status["force_report_sha256"]:
            raise RuntimeError("screened report modified after publication")
        stress, forces = read_json(sf), read_json(ff)
        s = structure(self.path)
        chosen = forces.get("selected_atoms", list(range(1, s.nat + 1)))
        if sorted(chosen) != list(range(1, s.nat + 1)):
            raise ValueError("relax=True requires screened forces on all inequivalent atoms")
        f = np.asarray(forces["force_mRy_Bohr"], float)
        if f.shape != self.ordinary.shape or not np.isfinite(f).all():
            raise ValueError("invalid complete screened force")
        sigma = np.asarray(stress["stress_GPa"], float)
        stress_residual(s, sigma, self.opts["relax_pressure"])
        lines = Path(self.case + ".scfm").read_text().splitlines()
        dis = [float(line.split()[-1]) for line in lines if line.startswith(":DIS")]
        if not dis or not np.isfinite(dis[-1]):
            raise ValueError("missing charge residual in final MIXER output")
        shutil.copyfile(self.path, self.directory / "measured.struct")
        record = dict(
            structure=str(self.directory / "measured.struct"),
            run_id=status["run_id"],
            report=str(sf),
            stress_GPa=sigma.tolist(),
            force_mRy_Bohr=f.tolist(),
            ordinary_force=self.ordinary.tolist(),
            correction=(f - self.ordinary).tolist(),
            charge_distance=dis[-1],
            charge_converged=0 <= dis[-1] <= p["cc"],
        )
        measured_overlap = self.directory / "SOlapm.measured.dat"
        if measured_overlap.exists():
            atomic_write("SOlapm.dat", measured_overlap.read_bytes())
        atomic_json(self.directory / "measurement.json", record)


def commit_cell(directory):
    """Roll a staged structure/charge pair forward after an interrupted update."""
    journal = directory / "cell_commit.json"
    plan = read_json(journal)
    if not plan:
        return
    pending = []
    for item in plan:
        dest, src = Path(item["dest"]), Path(item["src"])
        current = digest(dest)
        if current == item["after"]:
            continue
        if current != item["before"] or digest(src) != item["after"]:
            raise RuntimeError(f"external change during cell update: {dest}")
        pending.append((dest, src))
    for dest, src in pending:
        atomic_write(dest, src.read_bytes())
    journal.unlink()


def advance(case, directory, step, residual, opts):
    src, rho = Path(case + ".struct"), Path(case + ".clmsum")
    old = structure(src)
    cell = NativeLattice.from_structure(old).cell
    compliance, updated = updated_compliance(
        cell, residual, opts, read_json(directory / "elastic.json")
    )
    strain = stress_residual(old, cell_step(residual, opts, compliance), 0.0)
    out, seed = step / "next.struct", step / "next.clmsum"
    for _ in range(20):
        new = strained_struct(src, out, strain)
        if spheres_clear(new):
            break
        strain *= 0.5
    else:
        raise RuntimeError("cell step cannot keep current muffin-tin radii without overlap")
    if out.read_bytes() == src.read_bytes():
        raise RuntimeError("cell step is below STRUCT precision; inspect residual and settings")
    before, after = density_seed(src, out, rho, seed)
    plan = [
        dict(dest=str(p), src=str(q), before=digest(p), after=digest(q))
        for p, q in ((src, out), (rho, seed))
    ]
    atomic_json(
        directory / "elastic.json",
        dict(
            cell=cell.tolist(),
            next_cell=NativeLattice.from_structure(new).cell.tolist(),
            stress=residual.tolist(),
            compliance=compliance.tolist(),
            controls=[opts[k] for k in ELASTIC_CONTROLS],
        ),
    )
    atomic_json(directory / "cell_commit.json", plan)
    commit_cell(directory)
    reset_mixer_history(case, directory, f"{step.name} cell metric changed")
    log(
        directory,
        f"{step.name} cell updated to {new.a:.6f} {new.b:.6f} {new.c:.6f} Bohr; "
        f"restart charge {before:.10f} -> {after:.10f}, "
        + ("BFGS compliance" if updated else "initial elastic compliance"),
    )


def run(case, params, script):
    """Outer fixed-pressure loop; reuse the normal DMFT driver for each block."""
    opts = options(params)
    if params["finish"] < 1:
        raise ValueError("relax=True requires charge self-consistency (finish >= 1)")
    directory = Path(case + ".relax").resolve()
    directory.mkdir(exist_ok=True)
    with (directory / "lock").open("a+") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        commit_cell(directory)
        log(
            directory,
            f"{case} relaxation at Pext={opts['relax_pressure']:.5g} GPa. "
            f"Targets: max screened force {opts['relax_force_tol']} mRy/Bohr, max principal stress residual "
            f"{opts['relax_stress_tol']} GPa, rho <= {compact_number(params['cc'])}.",
        )
        with managed_headers(case, directory):
            # Electron/ion mixing and screened publication must use all allowed
            # sites. Refuse a diagnostic-only or selected-atom response request.
            config = read_json(Path(case + ".stress_response.json"), {})
            atoms = config.get("solver", {}).get("force_atoms")
            if config.get("enabled") is False or (
                atoms is not None
                and sorted(atoms) != list(range(1, structure(case + ".struct").nat + 1))
            ):
                raise ValueError("relax=True requires enabled STRFT response on all atoms")
            if config.get("solver", {}).get("directions") is not None:
                raise ValueError("relax=True requires the full symmetry-allowed stress")
            pending = directory / "pending_step.json"
            previous = read_json(pending)
            if previous:
                # A step limit leaves a measured structure available for inspection.
                # On restart, continue the saved proposal if its inputs are intact.
                if (
                    previous["options"] == opts
                    and previous["geometry"] == digest(case + ".struct")
                    and previous["density"] == digest(case + ".clmsum")
                ):
                    advance(
                        case,
                        directory,
                        Path(previous["step"]),
                        np.asarray(previous["residual"]),
                        opts,
                    )
                pending.unlink()
            for iteration in range(opts["relax_max_steps"]):
                index = 1 + max(
                    (int(p.name) for p in directory.iterdir() if p.is_dir() and p.name.isdigit()),
                    default=0,
                )
                step = directory / f"{index:04d}"
                step.mkdir()
                shutil.copyfile(case + ".struct", step / "initial.struct")
                pred = read_json(directory / "predictor.json")
                if pred and pred.get("case") != case:
                    raise ValueError("force predictor belongs to a different material")
                for file in (Path(case + ".in2"), Path(case + ".in2c")):
                    if file.exists() and file.stat().st_size:
                        set_mode(file, "FOR")
                env = dict(os.environ, EDMFT_RELAX_BLOCK=str(step))
                env.pop("EDMFT_RELAX_INITIAL", None)
                if pred:
                    previous_seed = Path(pred["structure"]).parent / "adjoint_seed.npz"
                    if previous_seed.is_file():
                        env["EDMFT_RELAX_INITIAL"] = str(previous_seed)
                st = structure(step / "initial.struct")
                log(
                    directory,
                    f"{step.name} fixed-cell MSR1a/FOR followed by fixed-position STRFT. "
                    f"Cell (Bohr) {st.a:.6f} {st.b:.6f} {st.c:.6f}; "
                    f"angles (degrees) {st.alpha:.6f} {st.beta:.6f} {st.gamma:.6f}. "
                    "During FOR, MIXER stress marked partial is incomplete. "
                    "The cell step uses the screened STRFT tensor below.",
                )

                atomic_json(
                    directory / "status.json", dict(step=index, phase="running", converged=False)
                )
                try:
                    subprocess.run(
                        [sys.executable, str(Path(script).resolve())],
                        env=env,
                        check=True,
                        pass_fds=(lock.fileno(),),
                    )
                except (subprocess.CalledProcessError, KeyboardInterrupt) as error:
                    log(
                        directory,
                        f"{step.name} electronic worker stopped: {error}. See {case}.dayfile, "
                        f"dmft_info.out and {case}.outputdmf2_response. Cell step not applied.",
                    )
                    atomic_json(
                        directory / "status.json",
                        dict(
                            step=index,
                            phase=(
                                "interrupted" if isinstance(error, KeyboardInterrupt) else "failed"
                            ),
                            converged=False,
                        ),
                    )
                    raise
                result = read_json(step / "measurement.json")
                if not result or digest(result["structure"]) != digest(case + ".struct"):
                    raise RuntimeError("missing screened measurement on the current structure")
                atomic_json(
                    directory / "predictor.json",
                    dict(case=case, structure=result["structure"], correction=result["correction"]),
                )
                s = structure(case + ".struct")
                residual = stress_residual(s, result["stress_GPa"], opts["relax_pressure"])
                stress = float(max(abs(np.linalg.eigvalsh(residual))))
                force = float(max(np.linalg.norm(result["force_mRy_Bohr"], axis=1)))
                converged = bool(
                    result["charge_converged"]
                    and stress <= opts["relax_stress_tol"]
                    and force <= opts["relax_force_tol"]
                )
                atomic_json(
                    directory / "status.json",
                    dict(
                        step=index,
                        phase="converged" if converged else "measured",
                        converged=converged,
                        max_force_mRy_Bohr=force,
                        max_stress_residual_GPa=stress,
                        charge_distance=result["charge_distance"],
                        measurement=str(step / "measurement.json"),
                    ),
                )
                log(
                    directory,
                    f"{step.name} screened max force {force:.6f} mRy/Bohr, "
                    f"max stress residual {stress:.6f} GPa, rho {compact_number(result['charge_distance'], '.2e')}",
                )
                for i, start in enumerate(np.r_[0, np.cumsum(s.mult)[:-1]]):
                    xyz = " ".join(f"{v:.8f}" for v in s.pos[start])
                    force_xyz = " ".join(f"{v:.5f}" for v in result["force_mRy_Bohr"][i])
                    delta_xyz = " ".join(f"{v:.5f}" for v in result["correction"][i])
                    log(
                        directory,
                        f"{step.name} atom {i+1} {s.aname[i]} fractional ({xyz}), screened force "
                        f"({force_xyz}), screened minus FOR ({delta_xyz}) mRy/Bohr",
                    )
                for label, row in zip("xyz", result["stress_GPa"]):
                    log(
                        directory,
                        f"{step.name} tensile stress {label} (GPa) "
                        + " ".join(f"{v:.6f}" for v in row),
                    )
                pressure = -np.trace(result["stress_GPa"]) / 3
                log(directory, f"{step.name} pressure (GPa) {pressure:.6f}")
                if converged:
                    log(directory, "converged on a fresh STRFT measurement")
                    return
                # Build cell secants from internally relaxed structures. When the
                # first STRFT substantially corrects FOR, let MSR1a equilibrate
                # the internal parameters at this cell before changing its metric.
                move_cell = (
                    result["charge_converged"]
                    and force <= opts["relax_force_tol"]
                    and stress > opts["relax_stress_tol"]
                )
                if iteration + 1 == opts["relax_max_steps"]:
                    if move_cell:
                        atomic_json(
                            pending,
                            dict(
                                options=opts,
                                geometry=digest(case + ".struct"),
                                density=digest(case + ".clmsum"),
                                step=str(step),
                                residual=residual.tolist(),
                            ),
                        )
                    log(
                        directory,
                        f"{step.name} reached relax_max_steps. Restart from this measured geometry to continue.",
                    )
                    raise RuntimeError(
                        "relax_max_steps reached; restart to continue from the measured structure"
                    )
                if move_cell:
                    advance(case, directory, step, residual, opts)
                else:
                    log(
                        directory,
                        f"{step.name} retaining cell for another charge/internal-coordinate block",
                    )
