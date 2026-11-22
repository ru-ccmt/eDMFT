"""Publish validated screened stress and/or atomic forces at SCF assembly.

WIEN mixer accumulates T=VAL+CORE-DENS and prints sigma=-T. Therefore the
additional VAL rows are -sigma_new-T_old. They REPLACE the old total; adding
sigma_new as another valence term would both double count and reverse its sign.
For forces the added FVA row replaces the native *total* after all FHF/FVA/
FSU/FCO pieces are assembled; it is not merely the induced-force term. Only
the newly assembled iteration is parsed, never prior SCF history.
"""

from pathlib import Path
import hashlib, json, re
import numpy as np


def mixer_components(text):
    """Read this iteration's raw tensors; reporting shares the publication parser."""
    arrays = {name: np.zeros((3, 3)) for name in ["VAL", "CORE", "DENS"]}
    seen = {name: set() for name in arrays}
    pattern = re.compile(r"^:STR_(VAL|CORE|DENS)(00[123]):\s*(.*)$")
    for line in text.splitlines():
        match = pattern.match(line)
        if match is None:
            continue
        name, row, values = match.groups()
        i = int(row) - 1
        fields = values.split()
        if len(fields) < 3:
            raise ValueError("incomplete native stress row")
        value = np.array([float(x.replace("D", "E")) for x in fields[:3]])
        if not np.all(np.isfinite(value)):
            raise ValueError("nonfinite native stress row")
        if name == "DENS":
            arrays[name][i] = value
        else:
            arrays[name][i] += value
        seen[name].add(i)
    if any(rows != {0, 1, 2} for rows in seen.values()):
        raise ValueError(
            "response publication requires complete DENS, VAL and CORE tensors; run lcore before SCF assembly"
        )
    return arrays


def mixer_accumulator(text):
    arrays = mixer_components(text)
    return arrays["VAL"] + arrays["CORE"] - arrays["DENS"]


def replacement_rows(text, sigma, *, reference_id):
    if ":STRRESP_BEGIN:" in text:
        raise ValueError("response already applied to this assembled iteration")
    sigma = np.asarray(sigma, float)
    if (
        sigma.shape != (3, 3)
        or not np.all(np.isfinite(sigma))
        or not np.allclose(sigma, sigma.T, atol=1e-12, rtol=0)
    ):
        raise ValueError("finite symmetric total tensor required")
    delta = -sigma - mixer_accumulator(text)
    rows = [f"\n:STRRESP_BEGIN: {reference_id} complete analytic replacement\n"]
    rows.extend(
        ":STR_VAL%03d:  " % (i + 1) + "".join(f"{v:20.10f}" for v in row) + "\n"
        for i, row in enumerate(delta)
    )
    rows.append(":STRRESP_END: symmetry-constrained frozen-Sigma Klein tangent\n")
    return "".join(rows)


def force_replacement_rows(text, report, *, reference_id):
    """Replace the *total* native force through one additive mixer FVA row.

    WIEN mixer sums FHF, FVA, FSU and FCO components in the site-local frame.
    The screened response evaluates the complete force, so adding its induced
    part to FVA alone would double count uncancelled terms. A displacement
    derivative is in global Cartesian axes. Convert its force to the same
    site-local Cartesian axes as the native components before subtraction.
    """
    if ":FORRESP_BEGIN:" in text:
        raise ValueError("screened force already applied to this assembled iteration")
    from wstruct import Struct
    structure = Struct(inAngs=False)
    structure.ReadStruct(report["structure"])
    force = np.asarray(report["force_mRy_Bohr"], float)
    if force.shape != (structure.nat, 3) or not np.all(np.isfinite(force)):
        raise ValueError("screened force report has invalid dimensions or values")
    active = {int(item["atom"]) for item in report["components"]}
    selected = set(report.get("selected_atoms", active))
    if not active <= selected or any(type(i) is not int or not 1 <= i <= structure.nat for i in selected):
        raise ValueError("invalid atomic force selection")

    from force_response_modes import check_force_frame
    check_force_frame(structure, selected)
    pieces = {atom: {tag: np.zeros(3) for tag in ("FHF", "FVA", "FSU", "FCO")}
              for atom in selected}
    seen = {atom: set() for atom in selected}
    pattern = re.compile(r"^:(FHF|FVA|FSU|FCO)(\d{3}):")
    for line in text.splitlines():
        match = pattern.match(line)
        if match is None:
            continue
        tag, atom = match[1], int(match[2])
        if atom not in selected:
            continue
        fields = line[17:].split()
        if len(fields) < 4:
            raise ValueError(f"incomplete {tag} force on atom {atom}")
        values = np.array([float(value.replace("D", "E")) for value in fields[1:4]])
        if not np.all(np.isfinite(values)):
            raise ValueError(f"nonfinite {tag} force on atom {atom}")
        pieces[atom][tag] += values
        seen[atom].add(tag)
    additive = report.get("application") == "add_soc_operator"
    label = "explicit SOC operator correction" if additive else "complete screened-force replacement"
    rows = [f"\n:FORRESP_BEGIN: {reference_id} {label}\n"]
    for atom in sorted(selected):
        if atom not in active and not seen[atom]:
            # Some native programs omit all force rows for a site fixed by symmetry.
            continue
        if atom in active and not {"FHF", "FVA", "FCO"}.issubset(seen[atom]):
            raise ValueError(f"native force parts missing on atom {atom}: {seen[atom]}")
        # WIEN reads its three STRUCT rows as columns (wstruct does likewise).
        # Mixer prints FCA = inv(rotloc) @ FOR before forming the internal FGL.
        local_force = structure.rotloc[atom-1] @ force[atom-1]
        delta = local_force if additive else local_force - sum(pieces[atom].values())
        rows.append(
            f":FVA{atom:03d}: {atom:3d}.ATOM"
            + "".join(f"{value:18.9f}" for value in (np.linalg.norm(delta), *delta))
            + f" {label}\n"
        )
    rows.append(f":FORRESP_END: {label}\n")
    return "".join(rows)


def response_rows(case, current_iteration, *, directory="."):
    """Return no rows for ordinary/deferred steps; reject stale ready results."""
    directory = Path(directory)
    path = directory / (case + ".stress_response.status.json")
    if not path.is_file():
        return ""
    status = json.loads(path.read_text())
    if status.get("backend") not in {"analytic","native_soc_forward"} or not (status.get("stress_ready") or status.get("force_ready")):
        return ""
    if status.get("status") != "complete" or status.get("converged") is not True:
        raise ValueError("analytic stress is not complete")

    def digest(path):
        return hashlib.sha256(path.read_bytes()).hexdigest()

    for name, expected in status["parent_fingerprints"].items():
        source = directory / name
        if not source.is_file() or digest(source) != expected:
            raise ValueError("stale analytic stress: changed input " + name)
    rows = ""
    if status.get("stress_ready"):
        result = Path(status["report"])
        if digest(result) != status["report_sha256"]:
            raise ValueError("analytic stress report changed")
        report = json.loads(result.read_text())
        if (
            report.get("converged") is not True
            or report.get("full_tensor_within_symmetry") is not True
            or report.get("diagnostic_only", False)
        ):
            raise ValueError("a complete supported tensor is required")
        rows += replacement_rows(
            current_iteration, report["stress_Ry_Bohr3"],
            reference_id=status["reference_id"],
        )
    if status.get("force_ready"):
        path = Path(status["force_report"])
        if digest(path) != status["force_report_sha256"]:
            raise ValueError("screened-force report changed")
        force_report = json.loads(path.read_text())
        if force_report.get("case") != case:
            raise ValueError("screened-force report belongs to a different case")
        rows += force_replacement_rows(
            current_iteration, force_report, reference_id=status["reference_id"]
        )
    return rows


def mark_published(case, correction, scf, *, directory=".", current_iteration=None):
    """Record publication only after the SCF write succeeds, under the run lock."""
    if not correction:
        return
    import fcntl
    from dmft2_step import atomic_json

    directory = Path(directory)
    path = directory / (case + ".stress_response.status.json")
    with (directory / ("." + case + ".dmft2-response.lock")).open("a+") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        status = json.loads(path.read_text())
        markers = []
        if status.get("stress_ready"):
            markers.append(":STRRESP_BEGIN: " + status["reference_id"] + " ")
        if status.get("force_ready"):
            markers.append(":FORRESP_BEGIN: " + status["reference_id"] + " ")
        if status.get("status") != "complete" or any(marker not in correction for marker in markers):
            raise RuntimeError("response state changed during SCF publication")
        status.update(stress_applied=bool(status.get("stress_ready")),
                      force_applied=bool(status.get("force_ready")),
                      published_scf=str(Path(scf).resolve()))
        atomic_json(path, status)
        # Report only after successful SCF publication. Do not read whole SCF
        # history: the caller passes the exact freshly assembled iteration.
        from stress_response_summary import published
        if status.get("stress_ready"):
            parts = mixer_components(current_iteration) if current_iteration is not None else None
            published(case, status, correction, components=parts, directory=directory)
