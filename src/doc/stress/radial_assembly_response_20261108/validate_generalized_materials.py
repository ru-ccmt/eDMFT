"""Algebraic moving-basis controls on saved full material exports, read-only.

The directions are synthetic matrix perturbations, not physical strain sources.
This validates the metric/projector interface on actual spectra and self-energy;
it does not construct native H/O derivatives or report a material stress.
"""

from pathlib import Path
import argparse
import json
import time

import numpy as np
from frozen_sigma_tangent import CachedFrozenSigmaResponse, real_inner
from frozen_sigma_response import herm


def validate(prefix):
    m, _ = CachedFrozenSigmaResponse.from_export(prefix, cache_mb=1024)
    start = time.monotonic()
    m.prepare()
    setup = time.monotonic() - start
    rng = np.random.default_rng(615)
    shape = (*m.eps.shape, m.eps.shape[1])
    h = 0.02 * herm(rng.normal(size=shape) + 1j * rng.normal(size=shape))
    o = 0.02 * herm(rng.normal(size=shape) + 1j * rng.normal(size=shape))
    y = 0.02 * herm(rng.normal(size=shape) + 1j * rng.normal(size=shape))
    raw = 0.02 * (rng.normal(size=m.raw.shape) + 1j * rng.normal(size=m.raw.shape))
    start = time.monotonic()
    f, info = m.generalized_response(h, o, raw)
    forward_time = time.monotonic() - start
    start = time.monotonic()
    hb, ob, rb = m.generalized_adjoint(y)
    reverse_time = time.monotonic() - start
    left = real_inner(y, f)
    right = real_inner(hb, h) + real_inner(ob, o) + real_inner(rb, raw)
    error = abs(left - right) / max(abs(left), abs(right), 1e-30)
    a = h + 1j * o
    ad = a.swapaxes(-1, -2).conj()
    basis_h = ad * m.eps[:, None, :] + m.eps[:, :, None] * a
    basis_o = ad + a
    basis_raw = np.einsum("knm,kmic->knic", ad, m.raw)
    rho = m.central_density()
    expected = -a @ rho - rho @ ad
    gauges = []
    for fixed_n in (False, True):
        actual, binfo = m.generalized_response(basis_h, basis_o, basis_raw, fixed_n)
        rel = float(np.linalg.norm(actual - expected) / np.linalg.norm(expected))
        gauges.append(
            dict(
                fixed_number=fixed_n,
                relative_error=rel,
                dN=binfo["dN"],
                dmu=binfo["dmu"],
            )
        )
        assert rel < 1e-10
        assert abs(binfo["dN"]) < 1e-10
    assert error < 1e-10
    assert abs(info["dN"]) < 1e-10
    return dict(
        prefix=str(prefix),
        shape=list(shape),
        setup_seconds=setup,
        forward_seconds=forward_time,
        adjoint_seconds=reverse_time,
        transpose_relative=error,
        fixed_number_charge=info["dN"],
        basis_gauge_controls=gauges,
        cache=m.cache_info(),
        physical_strain_sources=False,
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--export", action="append", required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    results = []
    for prefix in args.export:
        result = validate(prefix)
        results.append(result)
        args.output.write_text(json.dumps(results, indent=2) + "\n")
        print(json.dumps(result), flush=True)
