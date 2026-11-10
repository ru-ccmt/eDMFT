"""Publish a validated complete response tensor at SCF assembly, exactly once.

WIEN mixer accumulates T=VAL+CORE-DENS and prints sigma=-T. Therefore the
additional VAL rows are -sigma_new-T_old. They REPLACE the old total; adding
sigma_new as another valence term would both double count and reverse its sign.
Only the newly assembled iteration is parsed, never prior SCF history.
"""

from pathlib import Path
import hashlib, json, re
import numpy as np


def mixer_accumulator(text):
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


def response_rows(case, current_iteration, *, directory="."):
    """Return no rows for ordinary/deferred steps; reject stale ready results."""
    directory = Path(directory)
    path = directory / (case + ".stress_response.status.json")
    if not path.is_file():
        return ""
    status = json.loads(path.read_text())
    if status.get("backend") != "analytic" or not status.get("stress_ready"):
        return ""
    if status.get("status") != "complete" or status.get("converged") is not True:
        raise ValueError("analytic stress is not complete")

    def digest(path):
        return hashlib.sha256(path.read_bytes()).hexdigest()

    for name, expected in status["parent_fingerprints"].items():
        source = directory / name
        if not source.is_file() or digest(source) != expected:
            raise ValueError("stale analytic stress: changed input " + name)
    result = Path(status["report"])
    if digest(result) != status["report_sha256"]:
        raise ValueError("analytic stress report changed")
    report = json.loads(result.read_text())
    if (
        report.get("converged") is not True
        or report.get("full_tensor_within_symmetry") is not True
    ):
        raise ValueError("a complete supported tensor is required")
    return replacement_rows(
        current_iteration,
        report["stress_Ry_Bohr3"],
        reference_id=status["reference_id"],
    )


def mark_published(case, correction, scf, *, directory="."):
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
        if (
            status.get("status") != "complete"
            or ":STRRESP_BEGIN: " + status["reference_id"] + " " not in correction
        ):
            raise RuntimeError("response state changed during SCF publication")
        status.update(stress_applied=True, published_scf=str(Path(scf).resolve()))
        atomic_json(path, status)
