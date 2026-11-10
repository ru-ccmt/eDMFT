from pathlib import Path
import sys, json, time
import numpy as np

root = Path(__file__).resolve().parent
sys.path[:0] = [str(root / "stage"), "/Users/haule/projects/eDMFT_new/src/python"]
from lapw_xc_response import *
from wien_response_density import Density


def read_native(folder):
    with (folder / "xc_native.bin").open("rb") as f:
        nat, nkk, *shape = np.fromfile(f, "i4", 5)
        volume = np.fromfile(f, "f8", 1)[0]
        reps = np.fromfile(f, "i4", 3 * nkk).reshape(nkk, 3)
        rad = []
        for a in range(nat):
            nr, nl = np.fromfile(f, "i4", 2)
            lm = np.fromfile(f, "i4", 2 * nl).reshape(nl, 2)
            v = np.fromfile(f, "f8", nr * nl).reshape(nl, nr).T
            rad.append((lm, v))
        fft = np.fromfile(f, "c16", np.prod(shape)).reshape(shape, order="F")
        assert not f.read()
    return dict(shape=shape, volume=volume, reps=reps, radial=rad, fft=fft)


def read_geometry(folder, a, functional):
    with (folder / f"xc_geometry_{a+1:04d}.bin").open("rb") as f:
        nr, nl, ng = np.fromfile(f, "i4", 3)
        lm = np.fromfile(f, "i4", 2 * nl).reshape(nl, 2)
        r = np.fromfile(f, "f8", nr)
        Y = np.fromfile(f, "f8", nl * ng).reshape(ng, nl)
        w = np.fromfile(f, "f8", ng)
        if functional == "PBE":
            sine, cosine = np.fromfile(f, "f8", 2 * ng).reshape(2, ng)
            d = np.fromfile(f, "f8", 5 * nl * ng).reshape(5, ng, nl)
        else:
            sine = cosine = d = None
        assert not f.read()
    return (
        dict(
            r=r,
            ells=lm[:, 0],
            angular=Y,
            weights=w,
            functional=functional,
            sine=sine,
            cosine=cosine,
            angular_derivatives=d,
            harmonics=lm,
        ),
        lm,
    )


def err(x, y):
    return float(np.linalg.norm(x - y) / max(np.linalg.norm(y), 1e-30))


if __name__ == "__main__":
    report = {}
    for case in sys.argv[1:] or ["FeSe", "SrVO3"]:
        functional = "PBE" if case == "FeSe" else "LDA"
        folder = root / case / "xc"
        central = folder / "0"
        density = Density(central / (case + ".clmsum"), central / (case + ".struct"))
        nat = len(density.grids)
        refs = {
            h: read_native(folder / str(h)) for h in [0, 0.01, -0.01, 0.005, -0.005]
        }
        rows = []
        direction = np.load(
            "/private/tmp/edmft_hartree_periodic_20261108/" + case + "_direction.npz"
        )["density"]
        slices = {}
        off = 0
        for key, value in density.rad.items():
            slices[key] = slice(off, off + len(value))
            off += len(value)
        for a in range(nat):
            geom, lm = read_geometry(central, a, functional)
            keys = [(a, int(l), int(m)) for l, m in lm]
            c = np.array([density.rad[key] for key in keys]).T
            t = time.perf_counter()
            op = MuffinTinXCResponse(central=c, **geom)
            prep = time.perf_counter() - t
            dt = np.dtype(
                [
                    ("ir", "i4"),
                    ("k", "i4"),
                    ("n", "f8"),
                    ("g", "f8"),
                    ("lap", "f8"),
                    ("gdgg", "f8"),
                    ("v", "f8"),
                ]
            )
            points = np.fromfile(central / f"xc_points_{a+1:04d}.bin", dt).reshape(
                len(geom["r"]), len(geom["weights"])
            )
            point = unpolarized_xc(
                *[points[key] for key in ["n", "g", "lap", "gdgg"]],
                functional=functional,
                charge_cut=1e-8,
            )
            row = dict(
                atom=a,
                point_formula_relative=err(point, points["v"]),
                point_formula_max=float(np.max(abs(point - points["v"]))),
                preparation_seconds=prep,
            )
            if functional == "PBE":
                invariants = spherical_invariants(
                    op.fields(c), op.r[:, None], op.sine, op.cosine
                )
                row["invariants"] = {
                    key: err(val, points[key])
                    for key, val in zip(["n", "g", "lap", "gdgg"], invariants)
                }
            row["central_relative"] = err(
                op.central, refs[0]["radial"][a][1] * op.output_scale
            )
            row["transpose"] = op.operator.check_transpose(samples=2, tolerance=1e-10)
            dc = np.array([direction[slices[key]] for key in keys]).T
            dxc = op.action(dc)
            row["steps"] = {}
            for h in [0.01, 0.005]:
                fd = (
                    (refs[h]["radial"][a][1] - refs[-h]["radial"][a][1])
                    / (2 * h)
                    * op.output_scale
                )
                row["steps"][h] = dict(
                    relative=err(dxc, fd),
                    absolute=float(np.max(abs(dxc - fd))),
                    magnitude=float(np.max(abs(fd))),
                )
            rows.append(row)
            print(case, json.dumps(row), flush=True)
        report[case] = rows
    (root / "xc_mt_material.json").write_text(json.dumps(report, indent=2) + "\n")
