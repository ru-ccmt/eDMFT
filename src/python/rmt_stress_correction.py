#!/usr/bin/env python3
"""Measure total-energy sensitivity to a uniform muffin-tin-radius change.

This is a diagnostic, not a stress correction.  The published LAPW stress
uses fixed spherical muffin-tin boundaries; its frozen-augmentation
approximation explicitly discards the radial linear response.  A uniform
radius derivative therefore probes sensitivity to a different basis
partition and must not be added to the analytic stress.

For ``R -> (1+eta_R) R`` this program evaluates

    dE/deta_R = sum_a R_a dE/dR_a

and reports the stress-like scale

    S_R = -(1/(3 Omega)) sum_a R_a dE/dR_a

only to express that sensitivity in GPa.  ``S_R`` is not a stress-tensor
component and is never combined with the analytic tensor by this utility.

The derivative includes changes of the radial solutions, matching
coefficients, MT/interstitial partition, and WIEN2k's RMT-dependent energy
reference.  It is evaluated in two temporary calculations at
RMT*(1-delta) and RMT*(1+delta), with the plane-wave cutoff held fixed and
both densities converged by ordinary DFT cycles.

The program writes ``case.rmtstress`` plus ``case.rmtstress.json`` in the
reference directory; these historical filenames are retained for
compatibility.  It never edits ``case.struct``, ``case.clmsum``, or
``case.scf``.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict, dataclass
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
from typing import Iterable


RY_PER_BOHR3_TO_GPA = 14710.507641078
FLOAT = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?"
# WIEN2k reads the atom record with a fixed format.  The RMT field occupies
# exactly ten columns after ``RMT=``; expanding a short number in place shifts
# the following Z field and produces an opaque Fortran floating-point error.
RMT_RE = re.compile(r"(\bRMT\s*=)(.{10})(\s+Z:)", re.IGNORECASE)
RKMAX_RE = re.compile(
    r"^(\s*)([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?)(.*)$"
)


class RMTStressError(RuntimeError):
    """A controlled input, execution, or validation failure."""


@dataclass(frozen=True)
class SideResult:
    side: str
    scale: float
    directory: str
    energy_ry: float
    energy_history_ry: tuple[float, ...]
    matrix_sizes: tuple[int, ...]
    energy_records: int
    command_returncode: int
    converged: bool
    charge_distance: float | None


@dataclass(frozen=True)
class RadiusSensitivityResult:
    case: str
    reference_directory: str
    delta: float
    cycles: int
    energy_tolerance_ry: float
    charge_tolerance: float
    volume_bohr3: float
    energy_minus_ry: float
    energy_plus_ry: float
    radial_derivative_ry: float
    stress_like_scale_ry_per_bohr3: float
    stress_like_scale_gpa: float
    stress_like_scale_history_gpa: tuple[float, ...]
    is_stress_correction: bool
    analytic_stress_gpa: tuple[tuple[float, float, float], ...] | None
    analytic_pressure_gpa: float | None
    minus: SideResult
    plus: SideResult


def nonempty(path: Path) -> bool:
    return path.is_file() and path.stat().st_size > 0


def detect_case(directory: Path, requested: str | None) -> str:
    if requested:
        if not nonempty(directory / f"{requested}.struct"):
            raise RMTStressError(f"missing {directory / f'{requested}.struct'}")
        if not nonempty(directory / f"{requested}.clmsum"):
            raise RMTStressError(f"missing {directory / f'{requested}.clmsum'}")
        return requested

    candidates = [
        p.name[: -len(".struct")]
        for p in sorted(directory.glob("*.struct"))
        if nonempty(directory / (p.name[: -len(".struct")] + ".clmsum"))
    ]
    if not candidates:
        raise RMTStressError("no case has both case.struct and case.clmsum")
    if len(candidates) != 1:
        raise RMTStressError(
            "ambiguous cases: " + ", ".join(candidates) + "; use -f CASE"
        )
    return candidates[0]


def numbers(text: str) -> list[float]:
    return [float(x.replace("D", "E").replace("d", "e")) for x in re.findall(FLOAT, text)]


def last_scalar(paths: Iterable[Path], tag: str) -> float:
    value: float | None = None
    for path in paths:
        if not nonempty(path):
            continue
        for line in path.read_text(errors="replace").splitlines():
            if not line.lstrip().startswith(tag):
                continue
            tail = line.split("=", 1)[1] if "=" in line else line[len(tag) :]
            vals = numbers(tail)
            if vals:
                value = vals[0]
    if value is None or not math.isfinite(value):
        joined = ", ".join(str(p) for p in paths)
        raise RMTStressError(f"could not find a finite {tag} in {joined}")
    return value


def parse_stress(path: Path) -> tuple[tuple[float, float, float], ...] | None:
    """Read the last complete total stress tensor in GPa from case.scf."""

    if not nonempty(path):
        return None
    rows: dict[int, tuple[float, float, float]] = {}
    completed: tuple[tuple[float, float, float], ...] | None = None
    for line in path.read_text(errors="replace").splitlines():
        match = re.match(r"\s*:STRESS_GPa00([123]):(.*)$", line)
        if match is None or "total" not in line.lower():
            continue
        vals = numbers(match.group(2))
        if len(vals) < 3:
            continue
        rows[int(match.group(1))] = (vals[0], vals[1], vals[2])
        if len(rows) == 3:
            completed = (rows[1], rows[2], rows[3])
            rows = {}
    return completed


def parse_pressure(path: Path) -> float | None:
    if not nonempty(path):
        return None
    value = None
    for line in path.read_text(errors="replace").splitlines():
        if line.lstrip().startswith(":PRESSURE:") and "total" in line.lower():
            vals = numbers(line.split(":PRESSURE:", 1)[1])
            if vals:
                value = vals[0]
    return value


def scaled_struct(text: str, scale: float) -> str:
    """Scale every inequivalent-atom RMT while preserving the line layout."""

    count = 0

    def replace(match: re.Match[str]) -> str:
        nonlocal count
        count += 1
        old = float(match.group(2).strip().replace("D", "E").replace("d", "e"))
        field = f"{old * scale:10.5f}"
        if len(field) != 10:
            raise RMTStressError(f"scaled RMT does not fit WIEN2k's F10 field: {field}")
        return match.group(1) + field + match.group(3)

    result = RMT_RE.sub(replace, text)
    if count == 0:
        raise RMTStressError("no RMT= fields were found in case.struct")
    return result


def scaled_in1(text: str, scale: float) -> str:
    """Scale RKMAX so changing RMT does not change the plane-wave basis.

    WIEN2k specifies the cutoff as ``RKMAX = RMT_min*Kmax`` on the second
    line of case.in1/in1c.  When every RMT is multiplied by ``scale``, RKMAX
    must be multiplied by the same number to keep Kmax and hence the integer
    G-vector list fixed.  Otherwise the requested radius derivative is mixed
    with a discontinuous plane-wave cutoff derivative.
    """

    lines = text.splitlines(keepends=True)
    if len(lines) < 2:
        raise RMTStressError("case.in1/in1c has fewer than two lines")
    newline = "\n" if lines[1].endswith("\n") else ""
    body = lines[1][:-1] if newline else lines[1]
    match = RKMAX_RE.match(body)
    if match is None:
        raise RMTStressError("could not parse RKMAX on line 2 of case.in1/in1c")
    old = float(match.group(2).replace("D", "E").replace("d", "e"))
    # Five decimals retain substantially more than the precision normally
    # present in WIEN2k input while leaving the rest of the line untouched.
    replacement = f"{old * scale:.5f}"
    lines[1] = match.group(1) + replacement + match.group(3) + newline
    return "".join(lines)


def copy_reference_inputs(source: Path, destination: Path, case: str) -> None:
    """Copy the small calculation state needed for two DFT cycles.

    Large wave-function and diagnostic files are regenerated and would make
    each temporary side hundreds of megabytes.  Everything else at the top
    level is copied so uncommon WIEN2k input extensions remain supported.
    """

    destination.mkdir(parents=True, exist_ok=False)
    skip_starts = (
        f"{case}.almblm",
        f"{case}.vector",
        f"{case}.energy",
        f"{case}.output",
        f"{case}.scf",
        f"{case}.broyd",
        f"{case}.radwf",
        f"{case}.nsh",
        f"{case}.help",
        f"{case}.vsp",
        f"{case}.vns",
        f"{case}.clmval",
        f"{case}.corewf",
        f"{case}.dayfile",
    )
    for item in source.iterdir():
        if not item.is_file():
            continue
        if (
            item.name in {":log", ".processes"}
            or item.name.startswith(skip_starts)
            or item.name.endswith(("_new", ".error", ".def"))
        ):
            continue
        shutil.copy2(item, destination / item.name)

    for required in (f"{case}.struct", f"{case}.clmsum"):
        if not nonempty(destination / required):
            raise RMTStressError(f"temporary copy is missing required {required}")


def select_energy_only_inputs(directory: Path, case: str) -> tuple[str, ...]:
    """Disable force/stress output in temporary radius branches.

    The calibration uses only converged total energies.  Inheriting ``STR``
    from the reference calculation makes every temporary LAPW0/LAPW2 cycle
    evaluate the expensive analytic stress again, without changing the
    stationary density or energy.  Modify only the private temporary copies;
    the reference inputs remain untouched.
    """

    changed: list[str] = []
    for suffix in ("in0", "in2", "in2c", "in2up", "in2dn"):
        path = directory / f"{case}.{suffix}"
        if not nonempty(path):
            continue
        lines = path.read_text().splitlines(keepends=True)
        if not lines:
            continue
        first = lines[0]
        body = first.rstrip("\r\n")
        ending = first[len(body) :]
        match = re.match(r"^(\s*)(STRF?|FOR)(\s+.*)$", body, re.IGNORECASE)
        if match is None:
            continue
        lines[0] = match.group(1) + "TOT" + match.group(3) + ending
        path.write_text("".join(lines))
        changed.append(path.name)
    return tuple(changed)


def run_checked(command: list[str], cwd: Path, log: Path, allow_failure: bool = False) -> int:
    with log.open("a", encoding="utf-8") as stream:
        stream.write("\n$ " + " ".join(command) + "\n")
        stream.flush()
        completed = subprocess.run(
            command,
            cwd=cwd,
            stdout=stream,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
        )
    if completed.returncode and not allow_failure:
        raise RMTStressError(
            f"command failed with status {completed.returncode}: {' '.join(command)}; "
            f"see {log}"
        )
    return completed.returncode


def parse_side(side: str, scale: float, directory: Path, case: str, returncode: int) -> SideResult:
    scf = directory / f"{case}.scf"
    if not nonempty(scf):
        raise RMTStressError(f"{side}: run_lapw did not create {scf}")
    text = scf.read_text(errors="replace")
    energies = [
        float(x.replace("D", "E").replace("d", "e"))
        for x in re.findall(r":ENE\s+:.*?=\s*(" + FLOAT + r")", text)
    ]
    if not energies:
        raise RMTStressError(f"{side}: no :ENE record in {scf}")

    matrix_sizes: list[int] = []
    for line in text.splitlines():
        if line.startswith(":RKM") and "MATRIX SIZE" in line:
            match = re.search(r"MATRIX SIZE\s+(\d+)", line)
            if match:
                matrix_sizes.append(int(match.group(1)))

    charge = None
    for line in text.splitlines():
        if line.startswith(":DIS") and "CHARGE DISTANCE" in line:
            vals = numbers(line.split(")", 1)[-1])
            if vals:
                charge = vals[-1]

    return SideResult(
        side=side,
        scale=scale,
        directory=str(directory),
        energy_ry=energies[-1],
        energy_history_ry=tuple(energies),
        matrix_sizes=tuple(matrix_sizes),
        energy_records=len(energies),
        command_returncode=returncode,
        # run_lapw returns zero when its convergence criterion is satisfied.
        # Reaching the -i limit without convergence returns nonzero.
        converged=returncode == 0,
        charge_distance=charge,
    )


def validate_side_convergence(
    result: SideResult,
    max_cycles: int,
    log: Path,
) -> None:
    """Require convergence, while accepting fewer than ``max_cycles`` cycles.

    ``run_lapw -i N`` interprets N as an upper bound.  A well-converged side
    can therefore contain fewer than N energy records and must not be rejected.
    Conversely, N records accompanied by a nonzero return code mean that the
    iteration limit was reached without convergence and must not be used for a
    quantitative RMT derivative.
    """

    if result.converged:
        return
    charge = (
        "unknown"
        if result.charge_distance is None
        else f"{result.charge_distance:.6e}"
    )
    raise RMTStressError(
        f"{result.side}: run_lapw did not converge within {max_cycles} cycles "
        f"({result.energy_records} completed, final charge distance {charge}); "
        f"increase --cycles and see {log}"
    )


def make_lapw_command(
    run_lapw: str,
    case: str,
    max_cycles: int,
    energy_tolerance: float,
    charge_tolerance: float,
) -> list[str]:
    """Build a run_lapw command that requires both energy and charge tests."""

    return [
        run_lapw,
        "-f",
        case,
        "-i",
        str(max_cycles),
        "-ec",
        f"{energy_tolerance:.12g}",
        "-cc",
        f"{charge_tolerance:.12g}",
    ]


def common_matrix_size(
    minus_sizes: tuple[int, ...],
    plus_sizes: tuple[int, ...],
) -> int | None:
    """Return a common fixed basis size, independent of iteration counts.

    Early convergence can legitimately give the two sides different numbers
    of records.  Requiring the complete sequences to be identical therefore
    produced a false basis-change error.  What matters is that every recorded
    iteration on both sides used one and the same LAPW matrix size.
    """

    if not minus_sizes or not plus_sizes:
        return None
    minus_unique = set(minus_sizes)
    plus_unique = set(plus_sizes)
    if len(minus_unique) != 1 or minus_unique != plus_unique:
        return None
    return next(iter(minus_unique))


def sensitivity_history_gpa(
    minus_energies: tuple[float, ...],
    plus_energies: tuple[float, ...],
    delta: float,
    volume: float,
) -> tuple[float, ...]:
    """Return a final-aligned convergence history for unequal cycle counts."""

    common = min(len(minus_energies), len(plus_energies))
    if common < 1:
        raise RMTStressError("cannot form sensitivity history without both energies")
    minus_tail = minus_energies[-common:]
    plus_tail = plus_energies[-common:]
    return tuple(
        -((energy_plus - energy_minus) / (2.0 * delta))
        / (3.0 * volume)
        * RY_PER_BOHR3_TO_GPA
        for energy_minus, energy_plus in zip(minus_tail, plus_tail)
    )


def calculate_side(
    side: str,
    scale: float,
    source: Path,
    root: Path,
    case: str,
    cycles: int,
    energy_tolerance: float,
    charge_tolerance: float,
    x_command: str,
    run_lapw: str,
) -> SideResult:
    directory = root / side / case
    copy_reference_inputs(source, directory, case)
    select_energy_only_inputs(directory, case)
    struct = directory / f"{case}.struct"
    struct_new = directory / f"{case}.struct_new"
    struct_new.write_text(scaled_struct(struct.read_text(), scale))

    in1_candidates = [
        path for path in (directory / f"{case}.in1", directory / f"{case}.in1c")
        if nonempty(path)
    ]
    if not in1_candidates:
        raise RMTStressError(f"{side}: neither {case}.in1 nor {case}.in1c exists")
    for in1 in in1_candidates:
        in1.write_text(scaled_in1(in1.read_text(), scale))

    log = directory / "rmt_stress_run.log"
    run_checked([x_command, "clminter", "-f", case], directory, log)
    shutil.copy2(struct_new, struct)
    clmsum_new = directory / f"{case}.clmsum_new"
    if not nonempty(clmsum_new):
        raise RMTStressError(f"{side}: clminter did not create {clmsum_new}")
    shutil.copy2(clmsum_new, directory / f"{case}.clmsum")

    # ``cycles`` is a maximum.  run_lapw may stop earlier with return code zero
    # when it converges.  A nonzero return code after reaching the limit is not
    # accepted: an unconverged radial derivative can be wrong by several GPa.
    # run_lapw detects a complex calculation from case.in1c itself.  Unlike
    # the lower-level ``x lapw1 -c`` command, WIEN2k-24's run_lapw has no
    # ``-c`` option; supplying one prints an error and is version-dependent.
    lapw_command = make_lapw_command(
        run_lapw,
        case,
        cycles,
        energy_tolerance,
        charge_tolerance,
    )
    returncode = run_checked(
        lapw_command,
        directory,
        log,
        allow_failure=True,
    )
    result = parse_side(side, scale, directory, case, returncode)
    validate_side_convergence(result, cycles, log)
    return result


def write_outputs(reference: Path, result: RadiusSensitivityResult) -> tuple[Path, Path]:
    json_path = reference / f"{result.case}.rmtstress.json"
    text_path = reference / f"{result.case}.rmtstress"
    json_path.write_text(json.dumps(asdict(result), indent=2, sort_keys=True) + "\n")

    lines = [
        "# Uniform-RMT energy sensitivity (NOT A STRESS CORRECTION)",
        "# Do not add this diagnostic scale to the analytic stress tensor.",
        f":RMT_SENSITIVITY_DELTA: {result.delta:.12e}",
        f":RMT_SENSITIVITY_MAX_CYCLES: {result.cycles}",
        f":RMT_SENSITIVITY_EC_RY: {result.energy_tolerance_ry:.12e}",
        f":RMT_SENSITIVITY_CC: {result.charge_tolerance:.12e}",
        f":RMT_SENSITIVITY_COMPLETED_CYCLES: {result.minus.energy_records} {result.plus.energy_records}",
        f":RMT_SENSITIVITY_DERIV_RY: {result.radial_derivative_ry:.15e}",
        f":RMT_SENSITIVITY_SCALE_RYBOHR3: {result.stress_like_scale_ry_per_bohr3:.15e}",
        f":RMT_SENSITIVITY_SCALE_GPa: {result.stress_like_scale_gpa:.9f}",
        ":RMT_SENSITIVITY_HISTORY_FINAL_ALIGNED_GPa: "
        + " ".join(f"{value:.9f}" for value in result.stress_like_scale_history_gpa),
    ]
    text_path.write_text("\n".join(lines) + "\n")
    return text_path, json_path


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("directory", nargs="?", default=".", help="reference calculation directory")
    parser.add_argument("-f", "--case", help="WIEN2k case name; normally auto-detected")
    parser.add_argument("--delta", type=float, default=1.0e-3, help="fractional RMT step (default: 0.001)")
    parser.add_argument(
        "--cycles",
        type=int,
        default=20,
        help="maximum DFT cycles on each side; convergence is required (default: 20)",
    )
    parser.add_argument(
        "--ec",
        type=float,
        default=1.0e-6,
        help="required run_lapw energy convergence in Ry (default: 1e-6)",
    )
    parser.add_argument(
        "--cc",
        type=float,
        default=1.0e-5,
        help="required run_lapw charge convergence (default: 1e-5)",
    )
    parser.add_argument("--keep", action="store_true", help="keep temporary plus/minus calculations")
    parser.add_argument("--workdir", type=Path, help="parent for a retained work directory")
    parser.add_argument(
        "--allow-basis-change",
        action="store_true",
        help="accept unequal LAPW matrix sizes (normally an error)",
    )
    return parser.parse_args()


def main() -> int:
    options = parse_arguments()
    try:
        reference = Path(options.directory).expanduser().resolve()
        if not reference.is_dir():
            raise RMTStressError(f"not a directory: {reference}")
        if not (0.0 < options.delta < 0.02):
            raise RMTStressError("--delta must be positive and smaller than 0.02")
        if options.cycles < 1:
            raise RMTStressError("--cycles must be positive")
        if not math.isfinite(options.ec) or options.ec <= 0.0:
            raise RMTStressError("--ec must be finite and positive")
        if not math.isfinite(options.cc) or options.cc <= 0.0:
            raise RMTStressError("--cc must be finite and positive")

        case = detect_case(reference, options.case)
        volume = last_scalar(
            (reference / f"{case}.scf", reference / f"{case}.scfm"), ":VOL"
        )
        stress = parse_stress(reference / f"{case}.scf")
        pressure = parse_pressure(reference / f"{case}.scf")

        x_command = shutil.which("x")
        run_lapw = shutil.which("run_lapw")
        if x_command is None or run_lapw is None:
            raise RMTStressError(
                "x and run_lapw must be on PATH; activate the WIEN2k environment first"
            )

        if options.keep:
            parent = (options.workdir or Path(tempfile.gettempdir())).expanduser().resolve()
            parent.mkdir(parents=True, exist_ok=True)
            root = Path(tempfile.mkdtemp(prefix=f"{case}_rmtstress_", dir=parent))
            cleanup = False
        else:
            temporary = tempfile.TemporaryDirectory(prefix=f"{case}_rmtstress_")
            root = Path(temporary.name)
            cleanup = True

        try:
            with ThreadPoolExecutor(max_workers=2) as pool:
                minus_future = pool.submit(
                    calculate_side,
                    "minus",
                    1.0 - options.delta,
                    reference,
                    root,
                    case,
                    options.cycles,
                    options.ec,
                    options.cc,
                    x_command,
                    run_lapw,
                )
                plus_future = pool.submit(
                    calculate_side,
                    "plus",
                    1.0 + options.delta,
                    reference,
                    root,
                    case,
                    options.cycles,
                    options.ec,
                    options.cc,
                    x_command,
                    run_lapw,
                )
                minus = minus_future.result()
                plus = plus_future.result()

            if (
                not options.allow_basis_change
                and common_matrix_size(minus.matrix_sizes, plus.matrix_sizes) is None
            ):
                raise RMTStressError(
                    "LAPW matrix size changed or could not be read "
                    f"(minus={minus.matrix_sizes}, plus={plus.matrix_sizes}); "
                    "reduce --delta or use "
                    "--allow-basis-change only for a deliberate cutoff test"
                )

            derivative = (plus.energy_ry - minus.energy_ry) / (2.0 * options.delta)
            sensitivity_ry = -derivative / (3.0 * volume)
            sensitivity_gpa = sensitivity_ry * RY_PER_BOHR3_TO_GPA
            sensitivity_history = sensitivity_history_gpa(
                minus.energy_history_ry,
                plus.energy_history_ry,
                options.delta,
                volume,
            )

            result = RadiusSensitivityResult(
                case=case,
                reference_directory=str(reference),
                delta=options.delta,
                cycles=options.cycles,
                energy_tolerance_ry=options.ec,
                charge_tolerance=options.cc,
                volume_bohr3=volume,
                energy_minus_ry=minus.energy_ry,
                energy_plus_ry=plus.energy_ry,
                radial_derivative_ry=derivative,
                stress_like_scale_ry_per_bohr3=sensitivity_ry,
                stress_like_scale_gpa=sensitivity_gpa,
                stress_like_scale_history_gpa=sensitivity_history,
                is_stress_correction=False,
                analytic_stress_gpa=stress,
                analytic_pressure_gpa=pressure,
                minus=minus,
                plus=plus,
            )
            text_path, json_path = write_outputs(reference, result)
        finally:
            if cleanup:
                temporary.cleanup()

        print(f"case                              {case}")
        print(f"RMT step                          {options.delta:.9e}")
        print(f"maximum DFT cycles                {options.cycles}")
        print(f"energy/charge tolerances          {options.ec:.3e} Ry / {options.cc:.3e}")
        print(
            "completed DFT cycles              "
            f"R- {minus.energy_records}, R+ {plus.energy_records}"
        )
        print(f"E(R-)                             {minus.energy_ry:.12f} Ry")
        print(f"E(R+)                             {plus.energy_ry:.12f} Ry")
        print(f"sum R dE/dR                      {derivative:.12e} Ry")
        print("WARNING                           diagnostic only; do not add to stress")
        print(f"RMT stress-like sensitivity       {sensitivity_gpa:.9f} GPa")
        print(
            "sensitivity by DFT cycle          "
            + " ".join(f"{value:.9f}" for value in sensitivity_history)
            + " GPa"
        )
        if pressure is not None:
            print(f"analytic pressure                 {pressure:.9f} GPa")
        print(f"wrote                             {text_path}")
        print(f"wrote                             {json_path}")
        if options.keep:
            print(f"kept calculations                 {root}")
        return 0
    except RMTStressError as error:
        print(f"rmt_stress_correction.py: error: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
