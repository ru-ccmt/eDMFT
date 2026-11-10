from pathlib import Path
import sys, time, json
import numpy as np

root = Path(__file__).resolve().parent
sys.path[:0] = [
    str(root / "stage"),
    str(root),
    "/Users/haule/projects/eDMFT_new/src/python",
]
from lapw_core_response import *
from check_core import read_states, case_inputs
from check_xc_periodic import build, err


def build_core(case, layout):
    folder = root / case / "hxc_core" / "0"
    s, orbitals, vsp = case_inputs(folder, case)
    states = read_states(folder / "core_native.bin")
    i = 0
    atoms = {}
    empty = []
    for a, (shift, occupations) in enumerate(orbitals):
        if not occupations:
            empty.append(a)
            continue
        r = np.geomspace(s.r0[a], s.rmt[a], s.jrj[a])
        ext = CorePotentialExtension(r, vsp[a], shift_ry=shift)
        ops = []
        for occ in occupations:
            n = states[i]
            if n["nuc"] != 0:
                raise ValueError("only point-nucleus core validated")
            ops.append(
                MatchedCoreState(
                    n["r"],
                    n["v"],
                    n["e"],
                    z=n["z"],
                    kappa=n["kappa"],
                    match=n["match"],
                    end=n["end"],
                    occupation=occ,
                    relativistic=n["c"] < 1e10,
                    output_points=len(r),
                )
            )
            i += 1
        atoms[a] = AtomCoreResponse(ext, ops)
    return PeriodicCoreResponse(layout, atoms, empty_atoms=empty), orbitals


if __name__ == "__main__":
    reports = {}
    for case in sys.argv[1:] or ["FeSe", "SrVO3"]:
        xc, h, density = build(case)
        t = time.perf_counter()
        core, orbitals = build_core(case, h)
        prep = time.perf_counter() - t
        combined = core.operator @ xc.hxc_operator
        d = np.load(
            "/private/tmp/edmft_hartree_periodic_20261108/" + case + "_direction.npz"
        )["density"]
        analytic = combined.matvec(d)
        result = dict(
            core_preparation_seconds=prep,
            transpose=combined.check_transpose(samples=1, tolerance=1e-10),
            steps={},
        )
        for step in [0.01, 0.005]:
            endpoints = [
                read_states(root / case / "hxc_core" / str(q) / "core_native.bin")
                for q in [step, -step]
            ]
            numeric = np.zeros(h.size)
            i = 0
            atoms = []
            for a, (shift, occupations) in enumerate(orbitals):
                sl = core.slices[a]
                val = np.zeros(sl.stop - sl.start)
                for occ in occupations:
                    p, m = endpoints[0][i], endpoints[1][i]
                    i += 1
                    val += (
                        occ
                        * ((p["p"] ** 2 + p["q"] ** 2) - (m["p"] ** 2 + m["q"] ** 2))[
                            : len(val)
                        ]
                        / (2 * step)
                    )
                numeric[sl] = val
                atoms.append(err(analytic[sl], val))
            result["steps"][step] = dict(total=err(analytic, numeric), atoms=atoms)
        result["timing"] = {}
        for label, func in [
            ("core_forward", core.action),
            ("core_reverse", core.pullback),
            ("hxc_core_forward", combined.matvec),
        ]:
            times = []
            for _ in range(4):
                t = time.perf_counter()
                func(d)
                times.append(time.perf_counter() - t)
            result["timing"][label] = float(np.median(times[1:]))
        reports[case] = result
        print(case, json.dumps(result, indent=2), flush=True)
    (root / "hxc_core_material.json").write_text(json.dumps(reports, indent=2) + "\n")
