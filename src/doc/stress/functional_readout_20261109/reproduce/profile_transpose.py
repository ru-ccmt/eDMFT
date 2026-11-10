"""Profile one actual prepared FeSe charge-Jacobian transpose; no edits/runs."""

from pathlib import Path
import cProfile
import pstats
import numpy as np
from check_material import root
from lapw_valence_state import load_native_valence_state
from lapw_potential_core_state import load_potential_core_state
from lapw_screened_state import NativeScreenedState
from wien_response_density import Density

case = "FeSe"
source = Path("/private/tmp/edmft_full_k_response_20261108") / case / "charge/0"
with load_native_valence_state(
    source,
    case,
    complex_vectors=False,
    energy_policy="fixed",
    matrix_quadrature="wien2k24",
    library="/Users/haule/dbin/liblapw_response.so",
) as v:
    pc = load_potential_core_state(
        "/Users/haule/projects/eDMFT_new/src/doc/stress/charge_closure_20261108/"
        + case
        + "_potential_core.npz"
    )
    core = Density(
        Path("/private/tmp/edmft_xc_core_response_20261108")
        / case
        / "hxc_core/0"
        / (case + ".clmcor"),
        source / (case + ".struct"),
    )
    state = NativeScreenedState(v, pc, core, target_charge=120, fixed_number=True)
    c = np.load(root / "three_k" / case / "covectors.npz")["lattice_adjoint"]
    profiler = cProfile.Profile()
    profiler.enable()
    answer = state.jacobian.rmatvec(c)
    profiler.disable()
    profiler.dump_stats(root / "FeSe_transpose.prof")
    with (root / "FeSe_transpose_profile.txt").open("w") as f:
        pstats.Stats(profiler, stream=f).strip_dirs().sort_stats(
            "cumulative"
        ).print_stats(45)
    print(
        "Profiled actual native transpose; output norm",
        np.linalg.norm(answer),
        flush=True,
    )
