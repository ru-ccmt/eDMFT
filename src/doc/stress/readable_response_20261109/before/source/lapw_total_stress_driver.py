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
from contextlib import ExitStack
from response_parallel import SERIAL, ResponseParallel, limit_blas_threads
import argparse
import hashlib
import json
import os
import time
import resource
import sys
import numpy as np
from lapw_valence_state import load_native_valence_state
from lapw_potential_core_state import load_potential_core_state, fingerprint_files
from lapw_native_potential_export import prepare_native_potential_core
from lapw_screened_state import NativeScreenedState
from lapw_total_stress import PreparedTotalStress
from lapw_strain_charge import PreparedChargeSourceContraction
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
    parallel=SERIAL,
    cache_mb=256.0,
    basis_cache_mb=128.0,
    full_sources=False,
    save_geometry=False,
    ensemble="auto",
    diagnostic_gap_cutoff=None,
):
    """Prepare once, form invariant sources, solve one adjoint, save evidence."""
    thread_policy = limit_blas_threads()
    start = time.perf_counter()
    snapshot = Path(snapshot).resolve()
    output = Path(output).resolve()
    parallel.root_call(lambda: output.mkdir(parents=True, exist_ok=False))
    if not parallel.root:
        progress = lambda message: None
    library = Path(library or Path(__file__).with_name("liblapw_response.so")).resolve()
    if not library.is_file():
        raise ValueError("native radial library missing")
    from lapw_response_policy import detect_policy
    from lapw_magnetic_response import combine_valence
    provenance = snapshot / "preparation.json"
    override = json.loads(provenance.read_text()).get("physical_policy", {}) if provenance.is_file() else {}
    if ensemble not in ("auto", "fixed-n", "fixed-mu"):
        raise ValueError("invalid electronic ensemble")
    if ensemble != "auto":
        override["fixed_number"] = ensemble == "fixed-n"
    policy = detect_policy(snapshot, case, override=override)
    files = {
        name: snapshot / (case + "." + name)
        for name in ["struct", "clmsum", "in0", "vsp", "inc"]
    }
    expected = parallel.root_call(lambda: fingerprint_files(files))
    if potential_core is None:
        potential_core = output / "potential_core.npz"
        parallel.root_call(
            lambda: prepare_native_potential_core(snapshot, case, potential_core)
        )
    pc = load_potential_core_state(potential_core, expected_fingerprints=expected)
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
    with ExitStack() as stack:
        branches = [stack.enter_context(load_native_valence_state(
            snapshot, case, complex_vectors=complex_vectors, energy_policy="fixed",
            matrix_quadrature="wien2k24", library=library, parallel=parallel,
            cache_mb=cache_mb, basis_cache_mb=basis_cache_mb, spin_channel=channel,
            diagnostic_gap_cutoff=diagnostic_gap_cutoff,
        )) for channel in policy.channels]
        v = combine_valence(*branches) if policy.magnetic else branches[0]
        timings["prepare_valence_seconds"] = v.preparation_seconds
        state = NativeScreenedState(
            v, pc, core, target_charge=target, fixed_number=policy.fixed_number,
            normalization_enabled=policy.normalization_enabled,
        )
        t = time.perf_counter()
        total = PreparedTotalStress(
            state,
            snapshot,
            case,
            complex_vectors=complex_vectors,
            library=library,
            fixed_point_tolerance=fixed_point_tolerance,
            diagnostic=diagnostic_gap_cutoff is not None,
            progress=lambda i, x: progress(
                f"Prepared radial energy search {i+1}: atom {x['atom']+1}, l={x['l']}"
            ),
        )
        timings["prepare_functional_search_seconds"] = time.perf_counter() - t
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
                profile.dump_stats(
                    str(
                        output
                        / (
                            "adjoint_action.prof"
                            if parallel.size == 1
                            else f"adjoint_action.rank{parallel.rank}.prof"
                        )
                    )
                )
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
        solved = None
        readout = None
        if not full_sources:
            t = time.perf_counter()
            solved = total.energy.screened.solve(
                total.energy.charge_gradient,
                adjoint=True,
                rtol=2e-12,
                atol=1e-11,
                restart=30,
                max_cycles=8,
            )
            timings["adjoint_solve_seconds"] = time.perf_counter() - t
            t = time.perf_counter()
            readout = PreparedChargeSourceContraction(state, solved[0])
            timings["source_readout_seconds"] = time.perf_counter() - t
        sources = []
        timings["source_seconds"] = []
        for i, e in enumerate(total.strains):
            t = time.perf_counter()
            x = total.prepare_source(e, readout=readout)
            sources.append(x)
            timings["source_seconds"].append(time.perf_counter() - t)
            parallel.root_call(
                lambda: np.savez_compressed(
                    output / f"source{i}.npz",
                    strain=e,
                    source=x.source,
                    induced_Ry=np.empty(0) if x.induced is None else x.induced,
                    geometry=(
                        x.geometry
                        if save_geometry and parallel.size == 1
                        else np.empty(0)
                    ),
                    geometry_saved=save_geometry,
                    potential=x.potential,
                    geometry_is_rank_local=parallel.size > 1,
                )
            )
            if save_geometry and parallel.size > 1:
                np.savez_compressed(
                    output / f"geometry{i}.rank{parallel.rank:05d}.npz",
                    kpoints=[r.kpoint for r in v.records],
                    geometry=x.geometry,
                )
            progress(f"Prepared strain source {i+1}/{len(total.strains)}")
        sigma, lam, report = total.contract(sources, solved_adjoint=solved)
        timings["adjoint_action_seconds"] = action_times
        parallel.root_call(
            lambda: np.savez_compressed(
                output / "adjoint.npz",
                adjoint=lam,
                charge_gradient=total.energy.charge_gradient,
            )
        )
        if parallel.root_call(lambda: fingerprint_files(files)) != expected:
            raise RuntimeError("central snapshot changed during response")
        if parallel.fingerprints(v.fingerprint_paths) != v.fingerprints:
            raise RuntimeError("native valence snapshot changed during response")
        report.update(
            converged=True,
            case=case,
            physical_policy=policy.record(),
            snapshot=str(snapshot),
            source_fingerprints=expected,
            valence_fingerprints=v.fingerprints,
            radial_library_sha256=hashlib.sha256(library.read_bytes()).hexdigest(),
            timings=timings | dict(total_seconds=time.perf_counter() - start),
            peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
            * (1 if sys.platform == "darwin" else 1024),
            band_geometry_saved=save_geometry,
        )
    rank_statistics = parallel.gather(
        dict(
            rank=parallel.rank,
            kpoints=[int(r.kpoint) for r in v.records],
            peak_rss_bytes=report["peak_rss_bytes"],
            timings=report["timings"],
            electronic_cache=[b.model.cache_info() for b in branches],
        )
    )
    report["parallel"] = dict(
        ranks=parallel.size,
        distribution="k points",
        replicated_coordinates="real-space charge and potential; no MPI inner-product sum",
        rank_statistics=rank_statistics,
        blas_thread_policy=thread_policy,
    )
    report["peak_rss_bytes"] = max(r["peak_rss_bytes"] for r in rank_statistics)
    report["timings"]["total_seconds"] = max(
        r["timings"]["total_seconds"] for r in rank_statistics
    )
    parallel.root_call(
        lambda: (output / "report.json").write_text(
            json.dumps(report, indent=2, allow_nan=False) + "\n"
        )
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
    p.add_argument(
        "--mpi", action="store_true", help="distribute k points over MPI ranks"
    )
    p.add_argument(
        "--cache-mb", type=float, default=256.0, help="electronic cache budget per rank"
    )
    p.add_argument(
        "--basis-cache-mb",
        type=float,
        default=128.0,
        help="central spherical product cache budget per rank",
    )
    p.add_argument(
        "--full-sources",
        action="store_true",
        help="diagnostic: build every charge source instead of its adjoint contraction",
    )
    p.add_argument(
        "--save-geometry",
        action="store_true",
        help="diagnostic: retain local band-geometry arrays (large)",
    )
    p.add_argument("--ensemble", choices=["auto", "fixed-n", "fixed-mu"], default="auto",
                   help="auto uses preparation provenance; override for legacy exports")
    p.add_argument("--diagnostic-gap-cutoff", type=float,
                   help="omit singular window pairs ONLY in a non-publishable diagnostic")
    a = p.parse_args()
    parallel = ResponseParallel.world() if a.mpi else SERIAL
    run_prepared(
        a.snapshot,
        a.case,
        a.output,
        potential_core=a.potential_core,
        library=a.library,
        fixed_point_tolerance=a.fixed_point_tolerance,
        progress=lambda x: print(x, flush=True),
        parallel=parallel,
        cache_mb=a.cache_mb,
        basis_cache_mb=a.basis_cache_mb,
        full_sources=a.full_sources,
        save_geometry=a.save_geometry,
        ensemble=a.ensemble,
        diagnostic_gap_cutoff=a.diagnostic_gap_cutoff,
    )


if __name__ == "__main__":
    try:
        main()
    except BaseException:
        # Fail all ranks rather than hanging in another rank's next reduction.
        # Library users retain ordinary Python exceptions; the CLI owns MPI.
        import sys, traceback

        traceback.print_exc()
        if "--mpi" in sys.argv:
            from mpi4py import MPI

            MPI.COMM_WORLD.Abort(1)
        raise SystemExit(1)
