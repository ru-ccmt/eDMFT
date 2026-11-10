from pathlib import Path
import sys, json, time, re
import numpy as np

root = Path(__file__).resolve().parent
sys.path[:0] = [str(root / "stage"), "/Users/haule/projects/eDMFT_new/src/python"]
from lapw_core_response import MatchedCoreState, CorePotentialExtension
from wien_response_density import FLOAT
from wstruct import Struct


def read_states(path):
    states = []
    with path.open("rb") as f:
        while True:
            h = np.fromfile(f, "i4", 7)
            if not len(h):
                break
            n, l, k, npnt, end, match, nuc = h
            e, z, c, gamma = np.fromfile(f, "f8", 4)
            r, v, p, q = np.fromfile(f, "f8", 4 * npnt).reshape(4, npnt)
            states.append(
                dict(
                    n=n,
                    l=l,
                    kappa=k,
                    r=r,
                    v=v,
                    p=p,
                    q=q,
                    e=e,
                    z=z,
                    c=c,
                    end=end - 1,
                    match=match - 1,
                    nuc=nuc,
                )
            )
    return states


def case_inputs(folder, case):
    s = Struct()
    s.ReadStruct(str(folder / (case + ".struct")))
    lines = (folder / (case + ".inc")).read_text().splitlines()
    i = 0
    orbitals = []
    for a in range(s.nat):
        words = lines[i].split()
        num = int(words[0])
        shift = float(words[1])
        i += 1
        occ = []
        for j in range(num):
            occ.append(float(lines[i].split(",")[2].split()[0]))
            i += 1
        orbitals.append((shift, occ))
    text = (folder / (case + ".vsp")).read_text()
    vsp = []
    for atom in re.finditer(r"ATOMNUMBER\s*=\s*(\d+)", text):
        a = int(atom[1]) - 1
        end = text.find("ATOMNUMBER", atom.end())
        end = len(text) if end < 0 else end
        ch = re.search(r"VLM\(R\) FOR L\s*0\s+M=\s*0", text[atom.end() : end])
        offset = atom.end() + ch.end()
        vsp.append(
            np.array(
                [
                    float(m.group().replace("D", "E"))
                    for m in list(FLOAT.finditer(text[offset:end]))[: s.jrj[a]]
                ]
            )
        )
    return s, orbitals, vsp


if __name__ == "__main__":
    report = {}
    for case in sys.argv[1:] or ["FeSe", "SrVO3"]:
        folder = root / case / "core"
        s, orbitals, vsp = case_inputs(folder / "0", case)
        states = read_states(folder / "0/core_native.bin")
        index = 0
        rows = []
        total_time = 0
        assert len(states) == sum(len(occ) for shift, occ in orbitals)
        endpoints = {
            step: read_states(folder / str(step) / "core_native.bin")
            for step in [0.01, -0.01, 0.005, -0.005]
        }
        for a, (shift, occupations) in enumerate(orbitals):
            r = np.geomspace(s.r0[a], s.rmt[a], s.jrj[a])
            ex = CorePotentialExtension(r, vsp[a], shift_ry=shift)
            dv = 0.07 * r * np.sin(1.3 * r)
            dext = ex.action(dv)
            for occ in occupations:
                native = states[index]
                t = time.perf_counter()
                op = MatchedCoreState(
                    native["r"],
                    native["v"],
                    native["e"],
                    z=native["z"],
                    kappa=native["kappa"],
                    match=native["match"],
                    end=native["end"],
                    occupation=occ,
                    relativistic=native["c"] < 1e10,
                    output_points=len(r),
                )
                prepare = time.perf_counter() - t
                central = occ * (native["p"] ** 2 + native["q"] ** 2)[: len(r)]
                row = dict(
                    atom=a,
                    n=int(native["n"]),
                    kappa=int(native["kappa"]),
                    central_relative=float(
                        np.linalg.norm(op.density - central) / np.linalg.norm(central)
                    ),
                    extension_max=float(np.max(abs(ex.central - native["v"]))),
                    match_residual=float(op.residual),
                    preparation_seconds=prepare,
                )
                t = time.perf_counter()
                analytic = op.action(dext)
                total_time += time.perf_counter() - t
                row["transpose"] = op.operator.check_transpose(
                    samples=1, tolerance=1e-10
                )
                row["steps"] = {}
                for h in [0.01, 0.005]:
                    p, m = endpoints[h][index], endpoints[-h][index]
                    numeric = (
                        occ
                        * ((p["p"] ** 2 + p["q"] ** 2) - (m["p"] ** 2 + m["q"] ** 2))[
                            : len(r)
                        ]
                        / (2 * h)
                    )
                    de = (p["e"] - m["e"]) / (2 * h)
                    row["steps"][h] = dict(
                        density_relative=float(
                            np.linalg.norm(analytic - numeric)
                            / max(np.linalg.norm(numeric), 1e-30)
                        ),
                        energy_absolute=float(abs(op.energy_gradient @ dext - de)),
                        same_indices=bool(
                            p["match"] == m["match"] == native["match"]
                            and p["end"] == m["end"] == native["end"]
                        ),
                    )
                rows.append(row)
                index += 1
        report[case] = dict(states=rows, forward_seconds=total_time)
        print(
            case,
            "states",
            len(rows),
            "central",
            max(v["central_relative"] for v in rows),
            "FD",
            max(v["steps"][0.01]["density_relative"] for v in rows),
            "transpose",
            max(v["transpose"] for v in rows),
            "forward seconds",
            total_time,
            flush=True,
        )
    (root / "core_material.json").write_text(json.dumps(report, indent=2) + "\n")
