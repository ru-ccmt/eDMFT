#!/usr/bin/env python3
"""Compute a complete supported stress tensor from one prepared native snapshot.

This command uses one central snapshot, never +/- relaxed calculations. It
does not modify the parent SCF or publish a partly supported tensor to mixer.
Example:
  lapw_total_stress_driver.py --snapshot CASE.response_runs/ID/central \
      --case CASE --potential-core potential_core.npz --output total_stress

With --potential-core omitted, the archive is built from the isolated native
LAPW0/LCORE exports. Matrix/vector and dmft2 exports must be present as described
by lapw_valence_state. The same reusable entry point can be called by the dmft2
executor after its native preparation pass.
"""
from pathlib import Path
import argparse
import hashlib
import json
import os
import time
import resource
import numpy as np
from lapw_valence_state import load_native_valence_state
from lapw_potential_core_state import load_potential_core_state, fingerprint_files
from lapw_native_potential_export import prepare_native_potential_core
from lapw_screened_state import NativeScreenedState
from lapw_total_stress import PreparedTotalStress
from wien_response_density import Density


def run_prepared(
    snapshot,
    case,
    output,
    *,
    potential_core=None,
    library=None,
    fixed_point_tolerance=1e-3,
    progress=print,
):
    """Prepare once, form invariant sources, solve one adjoint, save evidence."""
    start = time.perf_counter()
    snapshot = Path(snapshot).resolve()
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=False)
    library = Path(library or Path(__file__).with_name("liblapw_response.so")).resolve()
    if not library.is_file():
        raise ValueError("native radial library missing")
    for suffix in ["inso", "indmfldn", "clmup", "clmdn"]:
        p = snapshot / (case + "." + suffix)
        if p.is_file() and p.stat().st_size:
            raise NotImplementedError("nonmagnetic/no-SOC snapshot required")
    files = {
        name: snapshot / (case + "." + name)
        for name in ["struct", "clmsum", "in0", "vsp", "inc"]
    }
    expected = fingerprint_files(files)
    if potential_core is None:
        potential_core = output / "potential_core.npz"
        prepare_native_potential_core(snapshot, case, potential_core)
    pc = load_potential_core_state(potential_core, expected_fingerprints=expected)
    mixer = (snapshot / (case + ".inm")).read_text().splitlines()[0].split()
    if len(mixer) < 3 or float(mixer[1]) != 0 or mixer[2].upper() != "YES":
        raise NotImplementedError("neutral normalized native charge policy required")
    import re

    unit5 = re.search(
        r"(?m)^\s*5\s*,\s*'([^']+)'", (snapshot / "lapw1.def").read_text()
    )
    if unit5 is None:
        raise ValueError("LAPW1 DEF unit 5 missing")
    complex_vectors = unit5[1].endswith(".in1c")
    core = Density(snapshot / (case + ".clmcor"), snapshot / (case + ".struct"))
    target = float(np.dot(pc.structure.mult, pc.structure.Znuc))
    timings = dict(prepare_potential_core_seconds=time.perf_counter() - start)
    with load_native_valence_state(
        snapshot,
        case,
        complex_vectors=complex_vectors,
        energy_policy="fixed",
        matrix_quadrature="wien2k24",
        library=library,
    ) as v:
        timings["prepare_valence_seconds"] = v.preparation_seconds
        state = NativeScreenedState(
            v, pc, core, target_charge=target, fixed_number=True
        )
        t = time.perf_counter()
        total = PreparedTotalStress(
            state,
            snapshot,
            case,
            complex_vectors=complex_vectors,
            library=library,
            fixed_point_tolerance=fixed_point_tolerance,
            progress=lambda i, x: progress(
                f"Prepared radial energy search {i+1}: atom {x['atom']+1}, l={x['l']}"
            ),
        )
        timings["prepare_functional_search_seconds"] = time.perf_counter() - t
        sources = []
        timings["source_seconds"] = []
        for i, e in enumerate(total.strains):
            t = time.perf_counter()
            x = total.prepare_source(e)
            sources.append(x)
            timings["source_seconds"].append(time.perf_counter() - t)
            np.savez_compressed(
                output / f"source{i}.npz",
                strain=e,
                source=x.source,
                geometry=x.geometry,
                potential=x.potential,
            )
            progress(f"Prepared strain source {i+1}/{len(total.strains)}")
        # Profiling is development-only and default-off. It measures a warm
        # action, not context construction or six repeated adjoint solves.
        from stress_response_adjoint import PairedOperator

        operator = total.energy.jacobian
        action_times = []

        def transpose(x):
            before = time.perf_counter()
            if os.environ.get("EDMFT_RESPONSE_PROFILE") and len(action_times) == 1:
                import cProfile

                profile = cProfile.Profile()
                answer = profile.runcall(operator.rmatvec, x)
                profile.dump_stats(str(output / "adjoint_action.prof"))
            else:
                answer = operator.rmatvec(x)
            action_times.append(time.perf_counter() - before)
            progress(f"Adjoint action {len(action_times)}: {action_times[-1]:.3f} s")
            return answer

        total.energy.screened.jacobian = PairedOperator(
            operator.input_size,
            operator.output_size,
            operator.matvec,
            transpose,
            "timed total adjoint",
        )
        sigma, lam, report = total.contract(sources)
        timings["adjoint_action_seconds"] = action_times
        np.savez_compressed(
            output / "adjoint.npz",
            adjoint=lam,
            charge_gradient=total.energy.charge_gradient,
        )
        if fingerprint_files(files) != expected:
            raise RuntimeError("central snapshot changed during response")
        report.update(
            converged=True,
            case=case,
            snapshot=str(snapshot),
            source_fingerprints=expected,
            valence_fingerprints=v.fingerprints,
            radial_library_sha256=hashlib.sha256(library.read_bytes()).hexdigest(),
            timings=timings | dict(total_seconds=time.perf_counter() - start),
            peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        )
    (output / "report.json").write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n"
    )
    progress("Mean tensile stress: %.9f GPa" % report["mean_tensile_stress_GPa"])
    return report


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--snapshot", required=True, type=Path)
    p.add_argument("--case", required=True)
    p.add_argument("--output", required=True, type=Path)
    p.add_argument("--potential-core", type=Path)
    p.add_argument("--library", type=Path)
    p.add_argument("--fixed-point-tolerance", type=float, default=1e-3)
    a = p.parse_args()
    run_prepared(
        a.snapshot,
        a.case,
        a.output,
        potential_core=a.potential_core,
        library=a.library,
        fixed_point_tolerance=a.fixed_point_tolerance,
        progress=lambda x: print(x, flush=True),
    )


if __name__ == "__main__":
    main()
