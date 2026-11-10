from pathlib import Path
import sys, json, time, importlib.util
import numpy as np

root = Path(__file__).resolve().parent
sys.path[:0] = [
    str(root / "stage"),
    str(root),
    "/Users/haule/projects/eDMFT_new/src/python",
]
from lapw_xc_response import *
from lapw_hartree_periodic import PeriodicHartreeResponse
from wstruct import Struct
from wien_response_density import Density
from check_xc import read_geometry, read_native, err

spec = importlib.util.spec_from_file_location(
    "hartree_reference",
    "/private/tmp/edmft_hartree_periodic_20261108/check_material.py",
)
hr = importlib.util.module_from_spec(spec)
spec.loader.exec_module(hr)


def build(case):
    folder = root / case / "xc"
    base = folder / "0"
    functional = "PBE" if case == "FeSe" else "LDA"
    density = Density(base / (case + ".clmsum"), base / (case + ".struct"))
    s = Struct()
    s.ReadStruct(str(base / (case + ".struct")))
    ref = read_native(base)
    h = PeriodicHartreeResponse(density, s, ref["shape"], inversion=case == "FeSe")
    mt = []
    for a in range(s.nat):
        geo, lm = read_geometry(base, a, functional)
        c = np.array([density.rad[a, int(l), int(m)] for l, m in lm]).T
        mt.append(MuffinTinXCResponse(central=c, **geo))
    ir = InterstitialXCResponse(
        h.reanalysis,
        np.array(list(density.pw.values())),
        [s.a, s.b, s.c],
        functional=functional,
    )
    return PeriodicXCResponse(h, mt, ir), h, density


def native_masked(h, fft):
    grid = np.fft.ifftn(np.fft.fftn(fft).real * h.reanalysis.mask_real)
    member = grid[h.reanalysis.indices] * h.reanalysis.weight.conj()
    return np.bincount(
        h.stars.owners, member.real, minlength=h.stars.nstar
    ) + 1j * np.bincount(h.stars.owners, member.imag, minlength=h.stars.nstar)


if __name__ == "__main__":
    reports = {}
    for case in sys.argv[1:] or ["FeSe", "SrVO3"]:
        t = time.perf_counter()
        xc, h, density = build(case)
        prep = time.perf_counter() - t
        folder = root / case / "xc"
        native = read_native(folder / "0")
        h0 = hr.read_native(
            Path("/private/tmp/edmft_hartree_periodic_20261108") / case / "0"
        )
        total_stars = xc.ir.central_stars + h0["potential"]
        grid = np.zeros(h.reanalysis.shape, complex)
        grid[h.reanalysis.indices] = total_stars[h.stars.owners] * h.reanalysis.weight
        result = dict(
            preparation_seconds=prep,
            central_total_fft_relative=err(grid, native["fft"]),
            ir_transpose=xc.ir.operator.check_transpose(samples=2, tolerance=1e-10),
            xc_transpose=xc.operator.check_transpose(samples=1, tolerance=1e-10),
            hxc_transpose=xc.hxc_operator.check_transpose(samples=1, tolerance=1e-10),
            vacuum_points=xc.ir.average_count,
            steps={},
        )
        d = np.load(
            "/private/tmp/edmft_hartree_periodic_20261108/" + case + "_direction.npz"
        )["density"]
        analytic = xc.hxc_operator.matvec(d)
        for step in [0.01, 0.005]:
            refs = [read_native(folder / str(q)) for q in [step, -step]]
            he = [
                hr.read_native(
                    Path("/private/tmp/edmft_hartree_periodic_20261108") / case / str(q)
                )
                for q in [step, -step]
            ]
            rows = []
            for a, l, m in h.keys:
                lm = refs[0]["radial"][a][0]
                j = np.nonzero(np.all(lm == [l, m], axis=1))[0][0]
                val = (
                    refs[0]["radial"][a][1][:, j]
                    - refs[1]["radial"][a][1][:, j]
                    + he[0]["radial"][a, l, m]
                    - he[1]["radial"][a, l, m]
                ) / (2 * step)
                if l == 0:
                    val *= xc.mt[a].r / np.sqrt(4 * np.pi)
                rows.append(val)
            fd_ir = (
                native_masked(h, refs[0]["fft"]) - native_masked(h, refs[1]["fft"])
            ) / (2 * step)
            numeric = np.concatenate(rows + [fd_ir.view(float)])
            result["steps"][step] = dict(
                radial=err(analytic[: h.pw_offset], numeric[: h.pw_offset]),
                interstitial=err(analytic[h.pw_offset :], numeric[h.pw_offset :]),
                all=err(analytic, numeric),
            )
        result["timing"] = {}
        for name, func in [
            ("xc_forward", xc.action),
            ("xc_reverse", xc.pullback),
            ("hxc_forward", xc.hxc_operator.matvec),
        ]:
            times = []
            for _ in range(4):
                t = time.perf_counter()
                func(d)
                times.append(time.perf_counter() - t)
            result["timing"][name] = float(np.median(times[1:]))
        reports[case] = result
        print(case, json.dumps(result, indent=2), flush=True)
    (root / "xc_periodic_material.json").write_text(
        json.dumps(reports, indent=2) + "\n"
    )
