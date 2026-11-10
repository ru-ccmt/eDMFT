"""Readable response journal shared by x_dmft and run_dmft SCF assembly.

This module only reports existing results. JSON reports and their hashes remain
the authority for stress publication. In particular, induced_Ry is NOT the
correction to the old mixer stress. Both decompositions are printed separately.
Only the launcher and MPI rank zero write this journal, sequentially. Reporting
errors warn without hiding a solver failure or invalidating completed physics.
"""
from datetime import datetime
from functools import wraps
from pathlib import Path
import hashlib
import json
import shutil
import sys

# Same conversion as the response JSON driver (stress_response_adjoint.py).
RY_BOHR3_TO_GPA = 14710.507641078


def best_effort(function):
    @wraps(function)
    def call(*args, **kwargs):
        try:
            return function(*args, **kwargs)
        except Exception as error:
            print(f"Response text summary: {function.__name__} failed: {error}",
                  file=sys.stderr, flush=True)
    return call


def summary_path(case, directory="."):
    return Path(directory) / (case + ".outputdmf2_response")


def _time():
    return datetime.now().astimezone().isoformat(timespec="seconds")


def _run_id(status):
    return status.get("run_id", status.get("reference_id", "unknown"))


def _append(path, text):
    with Path(path).open("a") as stream:
        stream.write(text)
        stream.flush()


@best_effort
def progress(path, message):
    """Append a flushed stage/action message; MPI callers use rank zero only."""
    if path is not None:
        _append(path, f"[{_time()}] {message}\n")


@best_effort
def begin(case, status, *, directory=".", reconstructed=False):
    """Start the latest-run view. Completed runs also keep an archived copy."""
    text = (
        f"DFT+eDMFT STRESS RESPONSE -- {case}\n"
        + "=" * 78 + "\n"
        + f"Run ID: {_run_id(status)}\n"
        + f"Summary created: {_time()}\n"
        + f"Execution: {status.get('phase', 'unknown')}; IN2 mode: "
          f"{status.get('requested_mode', 'unknown')}\n"
        + "Workflow: dmft2 -> isolated response -> lcore -> SCF assembly -> mixer\n"
        + "run_dmft: response is deferred until the final measurement pass.\n"
        + "Stress convention: tensile positive; compression pressure = -trace/3.\n"
        + "The native :PRESSURE: line currently prints trace/3 (mean tensile stress).\n"
        + "This file is informational; status/report JSON controls publication.\n\n"
    )
    if reconstructed:
        text += "Reconstructed from saved results; no calculation was rerun.\n\n"
    summary_path(case, directory).write_text(text)


def _archive(case, status, directory):
    if status.get("workspace"):
        work = Path(status["workspace"])
        if work.is_dir():
            shutil.copy2(summary_path(case, directory), work / "outputdmf2_response.txt")


def _number(value, precision=9):
    return "not recorded" if value is None else f"{float(value):.{precision}g}"


def _stress_table(rows):
    """Physical Cartesian stress contributions, including all allowed shear."""
    import numpy as np
    lines = ["All entries below are GPa; mean = (xx+yy+zz)/3.\n",
             f"{'Contribution':29s}" + "".join(f"{x:>17s}" for x in
             ("xx", "yy", "zz", "mean")) + "\n"]
    for label, values in rows:
        t = np.asarray(values, float)
        elements = [t[0, 0], t[1, 1], t[2, 2], np.trace(t)/3]
        lines.append(f"{label:29s}" + "".join(f"{v:17.6f}" for v in elements) + "\n")
    if any(np.max(abs(np.asarray(t)-np.diag(np.diag(t)))) > 1e-12 for _, t in rows):
        lines.append(f"{'Shear components':29s}" +
                     "".join(f"{x:>17s}" for x in ("yz", "xz", "xy")) + "\n")
        for label, values in rows:
            t = np.asarray(values)
            lines.append(f"{label:29s}" + "".join(f"{v:17.6f}" for v in
                         (t[1, 2], t[0, 2], t[0, 1])) + "\n")
    else:
        lines.append("All shear components in this table are zero.\n")
    return "".join(lines)


def _report_text(report):
    import numpy as np
    sigma = np.asarray(report["stress_Ry_Bohr3"], float) * RY_BOHR3_TO_GPA
    lines = ["\nCOMPLETE RESPONSE RESULT\n", "-" * 78 + "\n",
             _stress_table([("Complete stress", sigma)]),
             f"Compression-positive pressure: {-np.trace(sigma)/3:.9f} GPa\n",
             f"Scope: {report.get('scope', 'see JSON report')}\n"]
    if "sigma_symmetry_relative_error" in report:
        lines.append("Self-energy symmetry relative error: " +
                     _number(report["sigma_symmetry_relative_error"]) + "\n")
    if report.get("tensor_is_symmetry_projection"):
        lines.append("Tensor is projected onto the imposed charge symmetry.\n")

    # The explicit part includes radial-energy-search completion. The older
    # fixed-energy pieces alone do not sum to explicit_Ry; show that difference
    # explicitly rather than leaving an unexplained discrepancy in the table.
    if all(k in report for k in ("strains", "explicit_Ry", "induced_Ry", "derivatives_Ry")):
        e = report["strains"]
        explicit = np.asarray(report["explicit_Ry"])
        lines += ["\nFUNCTIONAL DERIVATIVE BREAKDOWN\n",
                  "D_nu = dK/deta for epsilon = eta*e_nu; entries are Ry per unit eta.\n",
                  "D_nu = a_nu + lambda.s_nu, with (I-J)^T lambda = c.\n",
                  "J: charge-update Jacobian; c: functional gradient in charge coordinates;\n",
                  "s_nu: explicit charge source; a_nu: explicit functional derivative.\n",
                  "This explicit/induced split is NOT the old-stress/net-correction split.\n"]
        for i, strain in enumerate(e):
            label = "hydrostatic" if np.allclose(strain, np.eye(3)) else "general"
            for axis, name in enumerate(("xx", "yy", "zz")):
                if np.allclose(strain, np.diag(np.eye(3)[axis])):
                    label = name
            lines.append(f"Direction {i+1} ({label}): e = {strain}\n")
        lines.append(f"{'Derivative contribution':36s}" +
                     "".join(f"{'direction '+str(i+1):>21s}" for i in range(len(e))) + "\n")
        parts = report.get("fixed_energy_explicit_parts", {})
        labels = {
            "lattice_geometry_Ry": "Lattice geometry",
            "lattice_explicit_potential_Ry": "Lattice explicit potential",
            "core_explicit_Ry": "Core explicit potential",
            "density_explicit_Ry": "Density functional explicit",
        }
        rows = [(labels.get(k, k), np.asarray(v)) for k, v in parts.items()]
        if parts:
            rows.append(("Radial-energy-search completion*",
                         explicit - sum(np.asarray(v) for v in parts.values())))
        rows += [("Explicit subtotal a_nu", explicit),
                 ("Screened induced lambda.s_nu", report["induced_Ry"]),
                 ("TOTAL D_nu", report["derivatives_Ry"])]
        for label, values in rows:
            lines.append(f"{label:36s}" + "".join(f"{v:21.12f}" for v in values) + "\n")
        if parts:
            lines.append("* Calculated as explicit_Ry minus the fixed-energy explicit pieces.\n")
        lines += ["Large individual terms cancel; only consistently combined totals are physical.\n",
                  "D_nu = Omega*sum_ij(sigma_ij*e_nu,ij). For hydrostatic strain,\n",
                  "mean tensile stress = D_hydro/(3*Omega), with conversion to GPa.\n"]

    for gap in report.get("window_gap_treatment", []):
        if gap.get("omitted_pairs", 0):
            channel = gap.get("spin_channel") or "primary"
            lines += [f"\nWARNING [{channel}]: Window cutoff uses band-pair omission.\n",
                      f"Omitted {gap['omitted_pairs']} singular cross-window band pairs on "
                      f"{gap['affected_kpoints']} k points; total k weight "
                      f"{gap['affected_kpoint_weight']:.8g}.\n",
                      f"Threshold: |energy gap| <= {gap['cutoff_Ry']:g} Ry. "
                      "All regular pairs, k points and central projectors are retained.\n"]
    fixed_point = report.get("charge_fixed_point", {})
    if fixed_point.get("warning"):
        lines.append("\nWARNING: " + fixed_point["warning"] + "\n")
    central = report.get("central", {})
    adjoint = report.get("adjoint", {})
    timing = report.get("timings", {})
    lines += ["\nCONVERGENCE AND COMPUTATIONAL WORK\n",
              f"Response solver converged: {report.get('converged', 'not recorded')}\n",
              f"k points: {central.get('kpoints', 'not recorded')}; retained bands: "
              f"{central.get('retained_bands', 'not recorded')}; MPI ranks: "
              f"{report.get('parallel', {}).get('ranks', central.get('mpi_ranks', 'not recorded'))}\n",
              "Ensemble: " + ("fixed total particle number" if
                  report.get('physical_policy', {}).get('fixed_number', central.get('fixed_number', True))
                  else "fixed chemical potential") + "; DMFT spin channels: " +
                  str(central.get('dmft_spin_channels', 'not recorded')) + "\n",
              "Charge normalization: " + str(central.get('normalization_enabled', 'not recorded')) +
                  "; common mu: " + _number(central.get('common_mu_Ry')) + " Ry\n",
              "Central charge fixed-point residual (native coordinates): " +
              _number(central.get("central_fixed_point_residual_norm")) + "\n",
              "Central charge warning threshold / within threshold: " +
              _number(fixed_point.get("warning_threshold")) + " / " +
              str(fixed_point.get("within_threshold", "not recorded")) + "\n",
              "Native valence/core reconstruction relative errors: " +
              _number(central.get("central_valence_relative_error")) + " / " +
              _number(central.get("central_core_relative_error")) + "\n",
              f"Adjoint: converged={adjoint.get('converged', 'not recorded')}, "
              f"operator applications={adjoint.get('actions', 'not recorded')}, "
              f"GMRES info={adjoint.get('gmres_info', 'not recorded')}\n",
              "Adjoint true residual / threshold: " + _number(adjoint.get("residual_norm")) +
              " / " + _number(adjoint.get("threshold")) + "\n"]
    for key, label in [
        ("prepare_potential_core_seconds", "Potential/core setup"),
        ("prepare_valence_seconds", "Electronic response setup"),
        ("prepare_functional_search_seconds", "Functional/radial-search setup"),
        ("adjoint_solve_seconds", "Adjoint solve"),
        ("source_readout_seconds", "Shared source readout"),
        ("total_seconds", "Response driver total"),
    ]:
        if key in timing:
            lines.append(f"{label:37s}{timing[key]:12.3f} s\n")
    if "source_seconds" in timing:
        lines.append("Individual strain-source times: " +
                     ", ".join(f"{v:.3f} s" for v in timing["source_seconds"]) + "\n")
    lines.append("Setup sub-times may overlap; do not sum them to estimate wall time.\n")
    return "".join(lines)


@best_effort
def finish(case, status, *, directory="."):
    path = summary_path(case, directory)
    state = status.get("status", "unknown")
    explanations = {
        "deferred_until_final": "Response deferred until the final run_dmft measurement pass.",
        "disabled": "Response disabled for this step; no response correction is available.",
        "skipped_noncharge": "No response for this non-charge dmft2 step.",
        "failed": "FAILED. This run has no publishable response correction.",
        "complete": "Response calculation complete.",
    }
    progress(path, explanations.get(state, state))
    if state == "failed":
        _append(path, "Failure: " + status.get("error", "unspecified") + "\n")
    if state == "complete" and status.get("report"):
        report = json.loads(Path(status["report"]).read_text())
        if "stress_Ry_Bohr3" in report:
            _append(path, _report_text(report))
        else:
            _append(path, "Legacy response backend: no complete stress tensor supplied.\n")
        if status.get("workspace"):
            work = Path(status["workspace"])
            prep = work / "central/preparation.json"
            if prep.is_file():
                times = json.loads(prep.read_text()).get("timings", [])
                _append(path, "\nIsolated native preparation (before response driver):\n" +
                        "".join(f"  {r['program']:30s}{r['seconds']:12.3f} s\n" for r in times))
            _append(path, "\nDETAILED RECORDS\n" +
                    f"Workspace: {work}\nReport: {status['report']}\n" +
                    f"Progress log: {work/'response.log'}\nNative preparation log: {work/'prepare.log'}\n" +
                    f"Input/executable hashes: {work/'reference.json'}\n" +
                    f"Saved inputs: {work/'reference'}\nNative exports: {work/'central'}\n")
        _append(path, f"Status/SCF publication record: {Path(directory).resolve()/(case+'.stress_response.status.json')}\n")
        if not status.get("stress_applied"):
            progress(path, "NEXT: lcore and SCF assembly. Original stress and net correction "
                     "will be recorded here after successful assembly.")
        else:
            progress(path, "Status JSON records successful SCF publication; see comparison below.")
    _archive(case, status, directory)


@best_effort
def published(case, status, correction, *, components=None, directory="."):
    """Report the actual assembly, after SCF write and status update succeed."""
    import numpy as np
    path = summary_path(case, directory)
    if not path.is_file() or f"Run ID: {_run_id(status)}\n" not in path.read_text():
        begin(case, status, directory=directory, reconstructed=True)
        finish(case, status, directory=directory)
    report = json.loads(Path(status["report"]).read_text())
    final = np.asarray(report["stress_Ry_Bohr3"], float)
    # These are added VAL rows, whose physical stress contribution has the
    # opposite sign. Reconstruct even when an older caller supplied no parts.
    values = [list(map(float, row.split()[1:4])) for row in correction.splitlines()
              if row.startswith(":STR_VAL")]
    delta_val = np.asarray(values)
    if delta_val.shape != (3, 3):
        raise ValueError("three publication correction rows required")
    old = final + delta_val
    rows = []
    if components is not None:
        rows = [("Density contribution (+D)", components["DENS"]),
                ("Valence contribution (-V)", -components["VAL"]),
                ("Core contribution (-C)", -components["CORE"])]
        old = components["DENS"] - components["VAL"] - components["CORE"]
    rows += [("Original total", old), ("Net response correction", final-old),
             ("FINAL TOTAL", final)]
    _append(path, "\nSCF PUBLICATION: ORIGINAL STRESS AND CORRECTION\n" + "-"*78 + "\n" +
            _stress_table([(label, value*RY_BOHR3_TO_GPA) for label,value in rows]) +
            "sigma_old = D-V-C; appended delta_VAL = sigma_old-sigma_new.\n" +
            "The extra STR_VAL rows carry a TOTAL replacement, not a pure valence term.\n" +
            f"Published SCF: {status.get('published_scf', 'not recorded')}\n" +
            f"Reference ID: {status.get('reference_id', 'not recorded')}\n" +
            "The response does not rewrite case.scf2. Mixer reads the appended SCF block.\n")
    progress(path, "SCF publication complete. Tensor is ready for mixer; mixer completion "
             "is not tracked by this journal. Its final output is case.scfm.")
    _archive(case, status, directory)


def _rebuild(case, *, directory="."):
    """Recreate a readable report from an existing result, without rerunning it.

    If already published, use the same reference's last SCF correction block
    and the preceding density/valence/core rows, never sums over SCF history.
    This also works after mixer has updated clmsum: it is historical reporting,
    not permission to reuse that stress on a new charge iterate.
    """
    directory = Path(directory)
    status = json.loads((directory/(case+".stress_response.status.json")).read_text())
    if status.get("report"):
        data = Path(status["report"]).read_bytes()
        if hashlib.sha256(data).hexdigest() != status.get("report_sha256"):
            raise ValueError("saved response report hash mismatch")
    correction = None
    components = None
    if status.get("stress_applied"):
        from stress_response_publish import mixer_components
        text = Path(status["published_scf"]).read_text()
        tag = ":STRRESP_BEGIN: " + status["reference_id"] + " "
        begin_at = text.rfind(tag)
        end_at = text.find(":STRRESP_END:", begin_at)
        if begin_at < 0 or end_at < 0:
            raise ValueError("matching publication block is absent from saved SCF")
        start = text.rfind(":STR_DENS001:", 0, begin_at)
        if start < 0:
            raise ValueError("original density stress missing before publication")
        components = mixer_components(text[start:begin_at])
        correction = text[begin_at:end_at]
    begin(case, status, directory=directory, reconstructed=True)
    finish(case, status, directory=directory)
    if correction is not None:
        published(case, status, correction, components=components, directory=directory)
    return summary_path(case, directory)


def rebuild(case, *, directory="."):
    """Hold the launcher lock while reconstructing a historical text summary."""
    import fcntl
    with (Path(directory)/("."+case+".dmft2-response.lock")).open("a+") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        return _rebuild(case, directory=directory)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Rebuild the readable response summary from saved results.")
    parser.add_argument("case")
    parser.add_argument("--directory", default=".")
    args = parser.parse_args()
    print(rebuild(args.case, directory=args.directory))
