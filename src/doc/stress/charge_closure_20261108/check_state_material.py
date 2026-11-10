"""Roundtrip native central potential/core contexts without scratch dependencies."""

from pathlib import Path
import sys, time, json
import numpy as np

root = Path(__file__).resolve().parent
prev = Path("/private/tmp/edmft_xc_core_response_20261108")
sys.path[:0] = [
    str(root / "stage"),
    str(prev),
    "/Users/haule/projects/eDMFT_new/src/python",
]
from lapw_potential_core_state import (
    save_potential_core_state,
    load_potential_core_state,
    fingerprint_files,
)
from check_xc_periodic import build, read_geometry, read_native
from check_hxc_core import build_core
from check_core import case_inputs, read_states


def relative(a, b):
    return float(np.linalg.norm(a - b) / max(np.linalg.norm(b), 1e-30))


reports = {}
for case in sys.argv[1:] or ["FeSe", "SrVO3"]:
    xc, h, d = build(case)
    core, _ = build_core(case, h)
    folder = prev / case / "hxc_core/0"
    xbase = prev / case / "xc/0"
    s, orbitals, vsp = case_inputs(folder, case)
    functional = "PBE" if case == "FeSe" else "LDA"
    angular = []
    for a in range(s.nat):
        geo, lm = read_geometry(xbase, a, functional)
        assert [tuple(x) for x in lm] == [(l, m) for aa, l, m in d.rad if aa == a]
        names = ["angular", "weights"] + (
            ["sine", "cosine", "angular_derivatives"] if functional == "PBE" else []
        )
        angular.append({key: geo[key] for key in names})
    states = read_states(folder / "core_native.bin")
    i = 0
    core_data = []
    for a, (shift, occ) in enumerate(orbitals):
        local = []
        for occupation in occ:
            state = states[i]
            i += 1
            assert state["nuc"] == 0 and state["z"] == s.Znuc[a]
            assert (state["c"] < 1e10) == (s.mode.strip() == "RELA")
            local.append(
                dict(
                    energy_Ha=float(state["e"]),
                    kappa=int(state["kappa"]),
                    occupation=float(occupation),
                    match=int(state["match"]),
                    end=int(state["end"]),
                )
            )
        core_data.append(dict(vsp=vsp[a], shift_Ry=float(shift), states=local))
    files = {name: xbase / (case + "." + name) for name in ["struct", "clmsum", "in0"]}
    files.update({name: folder / (case + "." + name) for name in ["vsp", "inc"]})
    archive = root / (case + "_potential_core.npz")
    args = dict(
        density=d,
        structure=s,
        fft_shape=read_native(xbase)["shape"],
        inversion=case == "FeSe",
        functional=functional,
        angular_data=angular,
        core_data=core_data,
        source_files=files,
    )
    t = time.perf_counter()
    save_potential_core_state(archive, **args)
    write = time.perf_counter() - t
    t = time.perf_counter()
    loaded = load_potential_core_state(
        archive, expected_fingerprints=fingerprint_files(files)
    )
    read = time.perf_counter() - t
    direction = np.load(
        "/private/tmp/edmft_hartree_periodic_20261108/" + case + "_direction.npz"
    )["density"]
    rng = np.random.default_rng(29)
    potential = rng.normal(size=h.size)
    result = dict(
        archive_bytes=archive.stat().st_size,
        prepare_and_save_seconds=write,
        load_and_prepare_seconds=read,
        hxc_action_relative=relative(
            loaded.potential_operator.matvec(direction),
            xc.hxc_operator.matvec(direction),
        ),
        hxc_reverse_relative=relative(
            loaded.potential_operator.rmatvec(potential),
            xc.hxc_operator.rmatvec(potential),
        ),
        core_action_relative=relative(
            loaded.core.action(potential), core.action(potential)
        ),
        core_reverse_relative=relative(
            loaded.core.pullback(potential), core.pullback(potential)
        ),
        hxc_transpose=loaded.potential_operator.check_transpose(
            samples=1, tolerance=1e-10
        ),
        core_transpose=loaded.core.operator.check_transpose(samples=1, tolerance=1e-10),
    )
    assert (
        result["hxc_action_relative"] < 1e-10 and result["hxc_reverse_relative"] < 1e-10
    )
    assert (
        result["core_action_relative"] < 1e-7 and result["core_reverse_relative"] < 1e-7
    )
    print(case, json.dumps(result, indent=2), flush=True)
    reports[case] = result
(root / "state_material.json").write_text(json.dumps(reports, indent=2) + "\n")
