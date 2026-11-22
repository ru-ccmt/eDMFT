"""Generate symmetry-preserving internal coordinates and native displaced force sources.

The analytic and displaced FORT methods share the same Cartesian orbit
patterns. These span internal motion compatible with the supplied STRUCT
and its charge representation. A symmetry-breaking displacement requires
electronic inputs prepared in the reduced symmetry.

The displaced method in this module also provides an independent check on
the analytic atomic source. Each endpoint holds the input density
coefficients and impurity self-energy fixed. One native pass supplies both
the explicit Klein derivative and the charge-map source. The two-step
default extrapolates the central differences to zero displacement. All modes
share the central screened adjoint.
"""

from __future__ import annotations

import json
from pathlib import Path
import re

import numpy as np

from lapw_lattice import NativeLattice
from wstruct import Struct
from force_screened_diagnostic import compact_endpoint, one_pair, prepare_endpoint
from lapw_potential_core_state import load_potential_core_state


def check_force_frame(structure, atoms):
    """Supported native mixer conventions for nonmagnetic force publication."""
    # FHF/FVA/FSU/FCO are site-local Cartesian for every WIEN lattice.
    # MIXER itself applies FORCINT to form FGL; using BR2 here would
    # rotate/scale a centred or oblique cell a second time.
    NativeLattice.from_structure(structure)
    for atom in atoms:
        rot = structure.rotloc[atom - 1]
        if not np.allclose(rot @ rot.T, np.eye(3), atol=4e-7, rtol=0):
            raise ValueError(f"nonorthogonal WIEN local rotation on atom {atom}")


def internal_modes(structure, atoms=None):
    """Return (inequivalent atom, Cartesian unit vector, orbit pattern).

    A site-stabilizer Reynolds average projects each Cartesian direction onto
    the allowed displacement subspace. All equivalent sites are then moved by
    the actual crystal operations, with consistency checked for every mapping.
    Thus each parameter moves an entire orbit and keeps the native symmetry.
    """
    lattice = NativeLattice.from_structure(structure)
    if not len(structure.timat) or not any(
        np.array_equal(r, np.eye(3)) and lattice.equivalent_displacement(t)
        for r, t in zip(structure.timat, structure.tau)
    ):
        raise ValueError("force modes require explicit STRUCT symmetry operations including identity")
    rotations = lattice.cartesian_rotations(structure.timat)
    selected = set(range(1, structure.nat + 1)) if atoms is None else set(atoms)
    if not selected or min(selected) < 1 or max(selected) > structure.nat:
        raise ValueError("force_atoms must contain inequivalent 1-based atom indices")
    modes = []
    first = 0
    for atom, mult in enumerate(structure.mult, start=1):
        representative = structure.pos[first]
        positions = structure.pos[first:first + mult]
        if atom in selected:
            mappings = []
            for position in positions:
                matches = [i for i, (r, tau) in enumerate(zip(structure.timat, structure.tau))
                           if lattice.equivalent_displacement(r @ representative + tau - position)]
                if not matches:
                    raise ValueError(f"symmetry does not generate atom {atom}'s orbit")
                mappings.append(matches)
            stabilizer = mappings[0]
            projection = sum(rotations[i] for i in stabilizer) / len(stabilizer)
            directions = []
            for axis in np.eye(3):
                vector = projection @ axis
                for previous in directions:
                    vector -= previous * (previous @ vector)
                norm = np.linalg.norm(vector)
                if norm > 1e-8:
                    directions.append(vector / norm)
            for direction in directions:
                pattern = np.zeros_like(structure.pos, dtype=float)
                for site, matches in enumerate(mappings):
                    images = [rotations[i] @ direction for i in matches]
                    if any(np.linalg.norm(image - images[0]) > 1e-7 for image in images[1:]):
                        raise ValueError(f"inconsistent site-symmetry displacement for atom {atom}")
                    pattern[first + site] = images[0]
                if not np.allclose(np.linalg.norm(pattern[first:first + mult], axis=1), 1,
                                   atol=1e-7, rtol=0):
                    raise ValueError("non-unit symmetry orbit displacement")
                modes.append((atom, direction, pattern))
        first += mult
    return modes


def displaced_struct(source, target, positions):
    """Change only positions, including WIEN's optional precise-position tail.

    Struct.WriteStruct hardcodes RELA, so rewriting the entire file would
    silently change NREL calculations. Preserve all nonposition source bytes.
    """
    source, target = Path(source), Path(target)
    lines = source.read_text().splitlines(keepends=True)
    number = len(positions)
    index = [i for i, line in enumerate(lines)
             if re.search(r"\bX=\s*[-+\d.]", line)]
    if len(index) != number:
        raise ValueError("STRUCT site count does not match parsed positions")
    # Native WIEN readers can consult either the eight-decimal site records
    # or the optional precise-position tail. Make them numerically identical:
    # a 5e-10 fractional mismatch otherwise perturbs large cancelling terms
    # in the force reference at an observable level.
    rounded = np.array([[float(f"{value:.8f}") for value in position]
                        for position in positions])
    for line_no, position in zip(index, rounded):
        line = lines[line_no]
        prefix = line[:line.index("X=")]
        ending = "\n" if line.endswith("\n") else ""
        # A minus sign makes -0.dddddddd eleven characters. The native
        # F10.8 reader truncates its last digit and shifts subsequent axes.
        # Omit the leading zero on negative fractions as WIEN permits.
        fields = [f"{v:.8f}".replace("-0.", "-.").rjust(10) for v in position]
        if any(len(field) != 10 for field in fields):
            raise ValueError("displaced position does not fit a WIEN F10.8 field")
        lines[line_no] = prefix + f"X={fields[0]} Y={fields[1]} Z={fields[2]}" + ending
    precise = next((i for i, line in enumerate(lines) if line.strip() == "Precise positions"), None)
    if precise is not None:
        if len(lines) < precise + 1 + number:
            raise ValueError("truncated STRUCT precise positions")
        for site, position in enumerate(rounded):
            lines[precise + 1 + site] = "".join(f"{value:20.15f}" for value in position) + "\n"
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text("".join(lines))


def resolve_displaced_modes(records, steps, nat):
    """Recover per-atom Cartesian forces for the displacements WIEN actually used.

    At finite STRUCT precision, row i can acquire a component along allowed
    mode j, especially in an oblique cell. At each requested step the measured
    derivatives obey y_i = sum_j B_ij f_j, where B contains the projections
    of the written displacement patterns divided by the nominal step. Solve
    this system (at most 3 by 3) before Richardson extrapolation. The raw
    directional derivatives stay in each record for independent comparisons.
    """
    force = np.zeros((nat, 3))
    columns = ("explicit_force_mRy_per_Bohr", "screened_force_mRy_per_Bohr",
               "full_force_mRy_per_Bohr")
    for atom in sorted({r["atom"] for r in records}):
        site = [r for r in records if r["atom"] == atom]
        at_step = []
        for index, h in enumerate(steps):
            rows = [r["steps"][index] for r in site]
            basis = np.asarray([r["sampled_mode_coefficients"] for r in rows])
            values = np.asarray([[r[key] for key in columns] for r in rows])
            try:
                at_step.append(np.linalg.solve(basis, values))
            except np.linalg.LinAlgError as exc:
                raise ValueError(f"written displacements at atom {atom} do not span "
                                 f"the requested modes at h={h:g} Bohr") from exc
        resolved = at_step[0]
        if len(steps) == 2:
            h1, h2 = steps
            resolved = (h1*h1*at_step[1] - h2*h2*at_step[0]) / (h1*h1 - h2*h2)
        for record, (explicit, induced, total) in zip(site, resolved):
            record["explicit_force_mRy_per_Bohr"] = float(explicit)
            record["screened_force_mRy_per_Bohr"] = float(induced)
            record["force_component_mRy_Bohr"] = float(total)
            force[atom-1] += total * np.asarray(record["direction_cartesian"])
    return force


def run_force_modes(response_run, workspace, *, steps=(0.01, 0.005),
                    atoms=None, use_mpi=True, progress=print):
    """Evaluate all requested invariant forces with one shared adjoint."""
    response_run, workspace = Path(response_run).resolve(), Path(workspace).resolve()
    report = json.loads((response_run / "total_stress/report.json").read_text())
    if not report.get("converged") or not report.get("adjoint", {}).get("converged"):
        raise ValueError("central charge adjoint is not converged")
    if report.get("diagnostic_only") or report["central"]["dmft_spin_channels"] != 1:
        raise NotImplementedError("screened atomic forces require a complete nonmagnetic response")
    if "frozen-Sigma" not in report.get("observable", ""):
        raise NotImplementedError("DFT atomic reference-Fermi-energy source is not validated")
    case = report["case"]
    structure_path = response_run / "reference" / (case + ".struct")
    structure = Struct(inAngs=False)
    structure.ReadStruct(str(structure_path))
    lattice = NativeLattice.from_structure(structure)
    modes = internal_modes(structure, atoms)
    if not modes:
        raise ValueError("no symmetry-allowed internal displacement for selected atoms")
    steps = tuple(float(h) for h in steps)
    if len(steps) not in (1, 2) or any(not np.isfinite(h) or h <= 0 for h in steps):
        raise ValueError("one or two positive finite force steps required")
    if len(steps) == 2 and np.isclose(steps[0], steps[1]):
        raise ValueError("force steps must be distinct")
    adjoint = np.load(response_run / "total_stress/adjoint.npz")["adjoint"]
    central = {"pc": load_potential_core_state(response_run / "central/potential_core.npz")}
    normalize = bool(report["central"]["normalization_enabled"])
    workspace.mkdir(parents=True, exist_ok=False)
    records = []
    for atom, direction, pattern in modes:
        site_patterns = [p for a, _, p in modes if a == atom]
        values = []
        for h in steps:
            label = f"atom{atom:03d}_axis{len([r for r in records if r['atom']==atom])+1}_h{h:g}".replace('.', 'p')
            snapshots = []
            for sign, side in [(-1, "minus"), (1, "plus")]:
                coordinate = structure.pos + sign * h * (pattern @ np.linalg.inv(lattice.cell).T)
                path = workspace / "structures" / label / side / (case + ".struct")
                displaced_struct(structure_path, path, coordinate)
                progress(f"Preparing {case} atom {atom}, h={h:g} Bohr, {side}")
                snapshot = prepare_endpoint(response_run, case, path, workspace,
                                            label + "_" + side, use_mpi)
                snapshots.append(compact_endpoint(snapshot, case, normalize))
            value = one_pair(central, *snapshots, h, int(structure.mult[atom - 1]), adjoint)
            minus, plus = (snapshot["pc"].structure.pos for snapshot in snapshots)
            actual_velocity = ((plus - minus) @ lattice.cell.T) / (2*h)
            # Equivalent sites obey the orbit rotations. Virtual work per
            # atom uses their average projection onto each allowed pattern.
            value["sampled_mode_coefficients"] = (
                np.einsum("sij,ij->s", site_patterns, actual_velocity) / structure.mult[atom-1]
            ).tolist()
            values.append(value)
        records.append(dict(atom=atom, direction_cartesian=direction.tolist(), steps=values))
    force = resolve_displaced_modes(records, steps, structure.nat)
    for record in records:
        progress(f"Screened force on atom {record['atom']} along {record['direction_cartesian']}: "
                 f"{record['force_component_mRy_Bohr']:.9f} mRy/Bohr")
    return dict(observable="symmetry-preserving frozen-Sigma screened atomic force",
                selected_atoms=sorted(set(atoms or range(1, structure.nat + 1))),
                response_run=str(response_run), case=case, structure=str(structure_path),
                steps_Bohr=steps,
                force_mRy_Bohr=force.tolist(), components=records,
                scope="nonmagnetic, fixed self-energy, charge-screened; native one-pass atomic sources")
