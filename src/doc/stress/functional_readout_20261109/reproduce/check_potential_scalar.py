"""Independent fixed-geometry lattice scalar potential-direction endpoints."""

from pathlib import Path
import json
import sys
import numpy as np
from check_material import native_scalar, root
from frozen_sigma_response import FrozenSigmaResponse
from lapw_potential_core_state import load_potential_core_state

case = sys.argv[1]
old = Path("/private/tmp/edmft_full_k_response_20261108")
source = old / case / "charge/0"
m, _ = FrozenSigmaResponse.from_export(source / "response")
arrays = np.load(root / "three_k" / case / "covectors.npz")
pc = load_potential_core_state(
    "/Users/haule/projects/eDMFT_new/src/doc/stress/charge_closure_20261108/"
    + case
    + "_potential_core.npz"
)
direction = np.load(
    "/private/tmp/edmft_hartree_periodic_20261108/" + case + "_direction.npz"
)["density"]
potential = pc.potential_operator.matvec(direction)
analytic = float(arrays["lattice_potential"] @ potential)
report = dict(analytic_Ry=analytic, steps=[])
for step in [0.001, 0.0005]:
    folders = [source.parent / str(sign * step) for sign in [-1, 1]]
    values = [native_scalar(p, m) for p in folders]
    fd = (values[1] - values[0]) / (2 * step)
    report["steps"].append(
        dict(step=step, values_Ry=values, fd_Ry=fd, error_Ry=fd - analytic)
    )
extrap = (4 * report["steps"][1]["fd_Ry"] - report["steps"][0]["fd_Ry"]) / 3
report["richardson"] = dict(
    fd_Ry=extrap,
    error_Ry=extrap - analytic,
    relative_error=abs(extrap - analytic) / abs(analytic),
)
(root / "three_k" / case / "potential_scalar.json").write_text(
    json.dumps(report, indent=2) + "\n"
)
print(case, json.dumps(report), flush=True)
