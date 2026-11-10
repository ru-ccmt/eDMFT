"""Native MIXER normalization after the independently computed XC/core endpoints."""

from pathlib import Path
import sys, json, time
import numpy as np

root = Path(__file__).resolve().parent
previous = Path("/private/tmp/edmft_xc_core_response_20261108")
sys.path[:0] = [
    str(root / "stage"),
    str(root / "stage/tests"),
    "/Users/haule/projects/eDMFT_new/src/python",
    str(previous),
]
from lapw_charge_response import NativeChargeNormalization, native_mixer_inversion
from wien_frozen_sigma_response import integrated_charge
from wien_response_density import Density
from wstruct import Struct
from check_core import case_inputs
from test_lapw_charge_response import native_endpoint


def output(case, step, fixed=False):
    branch = "fixed_n" if fixed and step else "charge"
    folder = previous / case / branch / str(step)
    sfile = folder / (case + ".struct")
    d = Density(folder / (case + ".clmval"), sfile)
    d.add(Density(previous / case / "hxc_core" / str(step) / (case + ".clmcor"), sfile))
    return d


reports = {}
for case in ["FeSe", "SrVO3"]:
    folder = previous / case / "charge/0"
    d = output(case, 0)
    s = Struct()
    s.ReadStruct(str(folder / (case + ".struct")))
    inv = native_mixer_inversion(s)
    _, orbitals, _ = case_inputs(previous / case / "hxc_core/0", case)
    ncore = sum(mult * sum(occ) for mult, (shift, occ) in zip(s.mult, orbitals))
    band_target = json.loads((previous / (case + "_fixed_n_roots.json")).read_text())[
        "target_band_number"
    ]
    target = band_target + ncore
    op = NativeChargeNormalization(d, s, target, inversion=inv)
    native, charges = native_endpoint(d, s, target, inv)
    np.testing.assert_allclose(op.central, native.vector(), rtol=2e-13, atol=2e-13)
    legacy = integrated_charge(
        d, folder / (case + ".struct"), quadrature="legacy_spline"
    )
    result = dict(
        inversion=inv,
        target=target,
        raw_charge=op.M + op.I,
        alpha=op.alpha,
        native_minus_legacy_charge=op.M + op.I - sum(legacy),
        transpose=op.operator.check_transpose(samples=2, tolerance=1e-11),
        branches={},
    )
    for fixed in [False, True]:
        label = "fixed_N" if fixed else "fixed_mu"
        raw = np.load(
            previous
            / (case + "_hxc_valence_core" + ("_fixed_n" if fixed else "") + ".npz")
        )["density"]
        t = time.perf_counter()
        analytic = op.action(raw)
        elapsed = time.perf_counter() - t
        row = dict(
            charge_derivative=float((op.mt + op.ir) @ analytic),
            unscaled_charge_derivative=float((op.mt + op.ir) @ raw),
            forward_seconds=elapsed,
            steps={},
        )
        for h in [0.001, 0.0005]:
            end = []
            for step in [h, -h]:
                endpoint = output(case, step, fixed)
                norm, _ = native_endpoint(endpoint, s, target, inv)
                end.append(norm.vector())
            fd = (end[0] - end[1]) / (2 * h)
            row["steps"][h] = dict(
                relative=float(np.linalg.norm(fd - analytic) / np.linalg.norm(fd)),
                maximum=float(np.max(abs(fd - analytic))),
            )
        result["branches"][label] = row
    reports[case] = result
    print(case, json.dumps(result, indent=2), flush=True)
(root / "normalized_material.json").write_text(json.dumps(reports, indent=2) + "\n")
