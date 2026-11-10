from pathlib import Path
import sys, json, time
import numpy as np

root = Path(__file__).resolve().parent
sys.path[:0] = [str(root / "stage"), "/Users/haule/projects/eDMFT_new/src/python"]
from wstruct import Struct
from wien_response_density import Density
from lapw_hartree_periodic import PeriodicHartreeResponse


def read_native(folder):
    with (folder / "hartree_native.bin").open("rb") as f:

        def get(n, dtype):
            x = np.fromfile(f, dtype, count=n)
            if len(x) != n:
                raise ValueError("short export")
            return x

        nat, nkk, ncon, *shape = get(6, "i4")
        volume = get(1, "f8")[0]
        reps = get(3 * nkk, "i4").reshape(nkk, 3)
        inst = get(nkk, "i4")
        norm = get(nkk, "f8")
        rho = get(nkk, "c16")
        pot = get(nkk, "c16")
        radial = {}
        charge = {}
        for a in range(nat):
            nr, nlm = get(2, "i4")
            lm = get(2 * nlm, "i4").reshape(nlm, 2)
            c = get(nr * nlm, "f8").reshape(nlm, nr)
            v = get(nr * nlm, "f8").reshape(nlm, nr)
            for (l, m), c1, v1 in zip(lm, c, v):
                radial[a, l, m] = v1
                charge[a, l, m] = c1
        mask = get(np.prod(shape), "c16").reshape(shape, order="F")
        masked = get(nkk, "c16")
        if f.read():
            raise ValueError("extra export data")
    return dict(
        volume=volume,
        ncon=ncon,
        shape=tuple(shape),
        reps=reps,
        inst=inst,
        norm=norm,
        rho=rho,
        potential=pot,
        radial=radial,
        charge=charge,
        mask=mask,
        masked=masked,
    )


def normerror(x, y):
    return float(np.linalg.norm(x - y) / max(np.linalg.norm(y), 1e-30))


if __name__ == "__main__":
    report = {}
    for case in sys.argv[1:] or ["FeSe", "SrVO3"]:
        folder = root / case / "0"
        template = Density(folder / (case + ".clmsum"), folder / (case + ".struct"))
        s = Struct()
        s.ReadStruct(str(folder / (case + ".struct")))
        native = read_native(folder)
        assert np.array_equal(list(template.pw), native["reps"])
        t = time.perf_counter()
        op = PeriodicHartreeResponse(
            template, s, native["shape"], ncon=native["ncon"], inversion=case == "FeSe"
        )
        prep = time.perf_counter() - t
        assert np.array_equal(op.stars.sizes, native["inst"])
        d = np.load(root / (case + "_direction.npz"))["density"]
        analytic = op.action(d)
        radial, pot, bdy = op.intermediates(d)
        result = {
            "preparation_seconds": prep,
            "mask_relative": normerror(
                op.reanalysis.mask_fourier, native["mask"] / native["volume"]
            ),
            "reanalysis_native": normerror(
                op.reanalysis.action(native["potential"]), native["masked"]
            ),
            "reanalysis_transpose": op.reanalysis.operator.check_transpose(
                samples=2, tolerance=1e-10
            ),
            "transpose": op.operator.check_transpose(samples=2, tolerance=1e-10),
            "steps": {},
        }
        for step in [0.01, 0.005]:
            plus = read_native(root / case / str(step))
            minus = read_native(root / case / str(-step))
            den = 2 * step
            fdpot = (plus["potential"] - minus["potential"]) / den
            fdir = (plus["masked"] - minus["masked"]) / den
            fdmt = []
            fdb = []
            components = []
            for key in op.keys:
                a, l, m = key
                r = op.integrations[a].r
                fd = (plus["radial"][key] - minus["radial"][key]) / den
                if l == 0:
                    fd /= np.sqrt(4 * np.pi)
                fdb.append(fd[-1])
                if l == 0:
                    fd = fd * r
                fdmt.append(fd)
            fd = np.concatenate(fdmt + [fdir.view(float)])
            result["steps"][step] = {
                "potential_star_sums": normerror(pot, fdpot),
                "boundary": normerror(bdy, fdb),
                "radial": normerror(analytic[: op.pw_offset], fd[: op.pw_offset]),
                "masked_PW": normerror(analytic[op.pw_offset :], fd[op.pw_offset :]),
                "all": normerror(analytic, fd),
            }
        result["timing"] = {}
        for name, func in [("forward", op.action), ("reverse", op.pullback)]:
            samples = []
            for i in range(8):
                t = time.perf_counter()
                func(d)
                samples.append(time.perf_counter() - t)
            result["timing"][name] = float(np.median(samples[1:]))
        report[case] = result
        print(case, json.dumps(result, indent=2), flush=True)
    (root / "material_periodic_hartree.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
