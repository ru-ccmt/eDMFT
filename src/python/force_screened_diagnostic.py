#!/usr/bin/env python3
"""Check a frozen-self-energy atomic force with one screened adjoint.

This is an explicit, default-off diagnostic. For each pair of displaced
structures, run ONE unmixed LAPW0/LAPW1/DMFT2/LCORE pass from the same central
input charge. It does not reconverge charge or solve the impurity problem.
The existing central stress-response adjoint is then contracted with the
finite-difference source of the actual normalized charge map.

The full first derivative is a_A + lambda.s_A. Both terms must be present:
for FeSe they are individually much larger than the final force. A second
step size allows removal of the leading O(h^2) difference error. The output
is a force of the frozen-Sigma Klein tangent, not an impurity-vertex result.

Two modes are supported:
  --pair H MINUS_SNAPSHOT PLUS_SNAPSHOT
      Read already prepared unmixed snapshots.
  --prepare-structures --workspace DIR --pair H MINUS.struct PLUS.struct
      Generate each unmixed snapshot inside DIR using the reference response
      run's exact inputs and workers; never edit the production calculations.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace
import gc
import re
import shutil

import numpy as np

from lapw_density_functional import PreparedDensityFunctional
from lapw_functional_cancellation import CoreEigenvalueFunctional
from lapw_lattice import NativeLattice
from lapw_potential_core_state import load_potential_core_state
from wien_response_density import Density


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prepare_endpoint(response_run, case, structure, workspace, label, use_mpi):
    """One isolated, unmixed native pass from the central material charge."""
    from lapw_response_prepare import prepare

    source = response_run / "reference"
    if not source.is_dir():
        raise ValueError(f"central reference inputs missing: {source}")
    structure = structure.resolve()
    if not structure.is_file() or structure.name != case + ".struct":
        raise ValueError(f"expected {case}.struct: {structure}")
    directory = workspace / label
    inputs, snapshot = directory / "input", directory / "snapshot"
    if (snapshot / "preparation.json").is_file():
        if sha256(snapshot / (case + ".struct")) != sha256(structure):
            raise ValueError(f"cached snapshot has a different structure: {snapshot}")
        if sha256(snapshot / (case + ".clmsum")) != sha256(source / (case + ".clmsum")):
            raise ValueError(f"cached snapshot has a different input charge: {snapshot}")
        return snapshot
    if directory.exists():
        raise ValueError(f"incomplete preparation directory exists: {directory}")
    directory.mkdir(parents=True)
    shutil.copytree(source, inputs)
    shutil.copy2(structure, inputs / (case + ".struct"))
    config = json.loads((response_run / "central/preparation.json").read_text())["config"]
    config = dict(config, reference=str(inputs), workspace=str(snapshot))
    if not use_mpi:
        config["native_mpi"] = []
        config["lapw1_mpi"] = []
    (directory / "config.json").write_text(json.dumps(config, indent=2) + "\n")
    prepare(config)
    return snapshot


def last_vsum(path: Path) -> float:
    values = []
    for line in path.read_text().splitlines():
        if line.startswith(":VSUM"):
            match = re.search(r"=\s*([-+\d.Ee]+)", line)
            if match is None:
                raise ValueError(f"malformed VSUM in {path}")
            values.append(float(match.group(1)))
    if not values:
        raise ValueError(f"VSUM missing in {path}")
    return values[-1]


def endpoint(snapshot: Path, case: str, normalization_enabled: bool):
    pc = load_potential_core_state(snapshot / "potential_core.npz")
    structure = snapshot / (case + ".struct")
    raw = Density(snapshot / (case + ".clmval"), structure)
    raw.add(Density(snapshot / (case + ".clmcor"), structure))
    if (tuple(raw.rad), tuple(raw.pw)) != (
        tuple(pc.input_density.rad), tuple(pc.input_density.pw)
    ):
        raise ValueError(f"output and input density layouts differ: {snapshot}")
    target = float(np.dot(pc.structure.mult, pc.structure.Znuc))
    normalized = pc.normalization(
        raw, target, enabled=normalization_enabled
    )
    channels = [snapshot / (case + ".scf2")]
    if (snapshot / (case + ".scf2dn")).is_file():
        channels.append(snapshot / (case + ".scf2dn"))
    vsum = float(np.mean([last_vsum(path) for path in channels]))
    density = float(PreparedDensityFunctional(pc).value_ry)
    core = float(CoreEigenvalueFunctional.from_state(pc).value_ry)
    return dict(
        pc=pc,
        normalized=normalized,
        energy=density + core + vsum,
        parts=dict(density_Ry=density, core_Ry=core, vsum_Ry=vsum),
    )


def compact_endpoint(snapshot, case, normalization_enabled):
    """Read only the values needed when comparing many displaced endpoints.

    The XC, core and Hartree operators keep bound methods of their owners.
    These reference cycles own large NumPy arrays. Python schedules cyclic
    collection by object allocation counts, which do not reflect array sizes.
    Reclaim the preceding native preparation before loading this endpoint,
    then retain only its geometry, charge coefficients and scalar readouts.
    The full endpoint() remains available for diagnostics needing operators.
    """
    gc.collect()
    result = endpoint(snapshot, case, normalization_enabled)
    pc = result["pc"]
    result["pc"] = SimpleNamespace(structure=pc.structure, input_density=pc.input_density)
    norm = result["normalized"]
    result["normalized"] = SimpleNamespace(central=norm.central, mt=norm.mt,
                                           ir=norm.ir, alpha=norm.alpha)
    del pc, norm
    gc.collect()
    return result


def verify_pair(center, minus, plus, h, multiplicity, adjoint):
    """Match fixed inputs, native coordinates and the physical displacement."""
    c, m, p = (x["pc"] for x in (center, minus, plus))
    structures = (m.structure, p.structure)
    for s in structures:
        for name in (
            "lattyp", "mode", "a", "b", "c", "alpha", "beta", "gamma",
            "mult", "jrj", "r0", "rmt", "timat", "tau", "rotloc", "iatnr", "Znuc",
        ):
            if not np.array_equal(getattr(c.structure, name), getattr(s, name)):
                raise ValueError(f"endpoint changes nonposition structural field {name}")
    for x in (m, p):
        if not np.array_equal(c.input_density.vector(), x.input_density.vector()):
            raise ValueError("endpoint input charge differs from central response")
        if (tuple(c.input_density.rad), tuple(c.input_density.pw)) != (
            tuple(x.input_density.rad), tuple(x.input_density.pw)
        ):
            raise ValueError("endpoint charge-coordinate layout differs")
    lattice = NativeLattice.from_structure(c.structure)
    dm = (m.structure.pos - c.structure.pos) @ lattice.cell.T
    dp = (p.structure.pos - c.structure.pos) @ lattice.cell.T
    if not np.allclose(dm + dp, 0.0, atol=2e-5, rtol=0):
        raise ValueError("plus/minus displacements are not opposite")
    radii = np.linalg.norm(dp, axis=1)
    moved = radii > 2e-5
    if np.count_nonzero(moved) != multiplicity:
        raise ValueError("moving-site count differs from multiplicity")
    if not np.allclose(radii[moved], h, atol=2e-5, rtol=0):
        raise ValueError("displacement magnitude differs from supplied step")
    if adjoint.shape != plus["normalized"].central.shape:
        raise ValueError("adjoint and charge-source sizes differ")
    return radii.tolist()


def one_pair(center, minus, plus, h, multiplicity, adjoint):
    displacements = verify_pair(center, minus, plus, h, multiplicity, adjoint)
    source = (plus["normalized"].central - minus["normalized"].central) / (2 * h)
    explicit = (plus["energy"] - minus["energy"]) / (2 * h)
    induced = float(adjoint @ source)
    total = explicit + induced
    return dict(
        step_Bohr=h,
        displaced_site_magnitudes_Bohr=displacements,
        fixed_input_energy_minus_Ry=minus["energy"],
        fixed_input_energy_plus_Ry=plus["energy"],
        energy_parts_minus=minus["parts"],
        energy_parts_plus=plus["parts"],
        explicit_derivative_Ry_per_Bohr=explicit,
        screened_derivative_Ry_per_Bohr=induced,
        total_derivative_Ry_per_Bohr=total,
        explicit_force_mRy_per_Bohr=-1000 * explicit / multiplicity,
        screened_force_mRy_per_Bohr=-1000 * induced / multiplicity,
        full_force_mRy_per_Bohr=-1000 * total / multiplicity,
        source_norm_per_Bohr=float(np.linalg.norm(source)),
        source_charge_per_Bohr=float(
            (minus["normalized"].mt + minus["normalized"].ir) @ source
        ),
        normalization_alpha_minus=float(minus["normalized"].alpha),
        normalization_alpha_plus=float(plus["normalized"].alpha),
    )


def run(args):
    response_run = args.response_run.resolve()
    report = json.loads((response_run / "total_stress/report.json").read_text())
    if report.get("converged") is not True or report["adjoint"].get("converged") is not True:
        raise ValueError("central adjoint is not converged")
    if report.get("diagnostic_only", False):
        raise ValueError("central response has an omitted-window diagnostic")
    if report["central"]["dmft_spin_channels"] != 1:
        raise NotImplementedError("magnetic atomic-force diagnostic needs spin validation")
    if "frozen-Sigma" not in report.get("observable", ""):
        raise NotImplementedError("only frozen-Sigma eDMFT is supported; DFT needs its reference-EF component")
    case = report["case"]
    normalization_enabled = bool(report["central"]["normalization_enabled"])
    adjoint = np.load(response_run / "total_stress/adjoint.npz")["adjoint"]
    central = {"pc": load_potential_core_state(response_run / "central/potential_core.npz")}
    if args.prepare_structures:
        if args.workspace is None:
            raise ValueError("--workspace required with --prepare-structures")
        args.workspace.mkdir(parents=True, exist_ok=True)
    pairs = []
    for index, (step, minus_path, plus_path) in enumerate(args.pair):
        h = float(step)
        if not np.isfinite(h) or h <= 0:
            raise ValueError("positive finite step required")
        minus_path, plus_path = Path(minus_path), Path(plus_path)
        if args.prepare_structures:
            minus_path = prepare_endpoint(
                response_run, case, minus_path, args.workspace,
                f"pair{index}_minus", args.use_mpi,
            )
            plus_path = prepare_endpoint(
                response_run, case, plus_path, args.workspace,
                f"pair{index}_plus", args.use_mpi,
            )
        minus = compact_endpoint(minus_path, case, normalization_enabled)
        plus = compact_endpoint(plus_path, case, normalization_enabled)
        value = one_pair(central, minus, plus, h, args.multiplicity, adjoint)
        value["minus_snapshot"] = str(minus_path.resolve())
        value["plus_snapshot"] = str(plus_path.resolve())
        pairs.append(value)
    result = dict(
        observable="frozen-Sigma Klein atomic-force diagnostic",
        response_run=str(response_run),
        case=case,
        multiplicity=args.multiplicity,
        central_fixed_point_residual_norm=report["central"]["central_fixed_point_residual_norm"],
        adjoint_residual_norm=report["adjoint"]["residual_norm"],
        pairs=pairs,
        note="Finite-difference explicit source at fixed input charge; one shared screened adjoint; no impurity vertex.",
    )
    if len(pairs) == 2:
        (a, b) = sorted(pairs, key=lambda item: item["step_Bohr"], reverse=True)
        h1, h2 = a["step_Bohr"], b["step_Bohr"]
        if np.isclose(h1, h2):
            raise ValueError("two distinct step sizes required")
        # Central differences have leading error proportional to h^2.
        result["extrapolated_force_mRy_per_Bohr"] = (
            h1**2 * b["full_force_mRy_per_Bohr"]
            - h2**2 * a["full_force_mRy_per_Bohr"]
        ) / (h1**2 - h2**2)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--response-run", type=Path, required=True)
    parser.add_argument("--pair", nargs=3, action="append", metavar=("STEP_BOHR", "MINUS", "PLUS"), required=True)
    parser.add_argument("--multiplicity", type=int, required=True)
    parser.add_argument("--prepare-structures", action="store_true")
    parser.add_argument("--workspace", type=Path)
    parser.add_argument("--use-mpi", action="store_true", help="inherit MPI commands from the central preparation")
    parser.add_argument("--output", type=Path, help="optional JSON result file")
    args = parser.parse_args()
    if args.multiplicity < 1 or len(args.pair) > 2:
        parser.error("one or two pairs and positive multiplicity required")
    result = run(args)
    text = json.dumps(result, indent=2, allow_nan=False) + "\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
