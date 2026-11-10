"""Resolve finite-step chemical-potential effects without new native runs."""

from pathlib import Path
import json
import sys
import numpy as np
from check_material import native_scalar, root
from frozen_sigma_response import FrozenSigmaResponse
from lapw_strain_potential import unequal_centered_derivative

case = sys.argv[1]
mesh = sys.argv[2] if len(sys.argv) > 2 else "three_k"
old = Path("/private/tmp/edmft_full_k_response_20261108")
source = (old / mesh if mesh == "full_bz" else old) / case / "charge/0"
m, _ = FrozenSigmaResponse.from_export(source / "response")
central = native_scalar(source, m)
native = Path("/private/tmp/edmft_screened_native_20261108") / mesh / case
report_path = root / mesh / case / "REPORT.json"
analytic = json.loads(report_path.read_text())["explicit_lattice_derivative_Ry"]
rows = []
for step in [0.00025, 0.000125] + ([0.0000625] if mesh == "three_k" else []):
    if step == 0.0000625:
        base = Path("/private/tmp/edmft_strain_charge_20261109/native_extra") / case
        folders = [base / f"source_{label}_6.25e-05" for label in ["m", "p"]]
    else:
        folders = [native / f"source_{label}_{step}" for label in ["m", "p"]]
    records = [json.loads((p / "done.json").read_text()) for p in folders]
    offsets = [r["actual_scale"] - 1 for r in records]
    row = dict(step=step, offsets=offsets)
    for policy in ["fixed_mu", "legendre_fixed_number"]:
        kw = (
            [{}, {}]
            if policy == "fixed_mu"
            else [
                dict(mu=r["mu_Ry"], number=r["request"]["target_band"]) for r in records
            ]
        )
        values = [native_scalar(p, m, **k) for p, k in zip(folders, kw)]
        fd = float(unequal_centered_derivative(values[0], central, values[1], *offsets))
        row[policy] = dict(values_Ry=values, fd_Ry=fd, error_Ry=fd - analytic)
    rows.append(row)
    print(case, mesh, row, flush=True)
report = dict(analytic_Ry=analytic, rows=rows, richardson=[])
for a, b in zip(rows, rows[1:]):
    ratio = (
        (a["offsets"][1] - a["offsets"][0]) / (b["offsets"][1] - b["offsets"][0])
    ) ** 2
    result = dict(steps=[a["step"], b["step"]])
    for policy in ["fixed_mu", "legendre_fixed_number"]:
        extrap = (ratio * b[policy]["fd_Ry"] - a[policy]["fd_Ry"]) / (ratio - 1)
        result[policy] = dict(
            fd_Ry=extrap,
            error_Ry=extrap - analytic,
            relative_error=abs(extrap - analytic) / abs(analytic),
        )
    report["richardson"].append(result)
    print(case, mesh, "RICHARDSON", result, flush=True)
(root / mesh / case / "scalar_refinement.json").write_text(
    json.dumps(report, indent=2) + "\n"
)
