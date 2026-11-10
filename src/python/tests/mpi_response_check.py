"""Run with mpiexec -n 1/2/4 ... mpi_response_check.py; no native fixtures needed.

Compare distributed complex, real-linear symmetry, dX, dmu and generalized
metric responses to an independent SERIAL context containing the entire mesh.
Test both cache and no-cache paths, both ensembles, and collective transposes.
"""

import json
import os
import warnings
import numpy as np
from response_parallel import ResponseParallel
from frozen_sigma_tangent import CachedFrozenSigmaResponse, real_inner
from frozen_sigma_response import herm
from test_frozen_sigma_tangent import complex_model
from stress_response_adjoint import PairedOperator

p = ResponseParallel.world()
try:
    central = complex_model()
    copies = 4
    eps = np.tile(central.eps, (copies, 1))
    if os.environ.get("EDMFT_TEST_WINDOW_DEGENERACY"):
        eps[0, 4] = eps[0, 3]
        warnings.filterwarnings("ignore", message="Window cutoff:", category=RuntimeWarning)
    raw = np.tile(central.raw, (copies, 1, 1, 1))
    weights = np.tile(central.kw, copies) / copies
    normalized = central.apply_x(raw, central.X)
    ids = p.partition(np.arange(len(eps)) + 1)
    rng = np.random.default_rng(1861)
    h = herm(
        rng.normal(size=(*eps.shape, eps.shape[1]))
        + 1j * rng.normal(size=(*eps.shape, eps.shape[1]))
    )
    o = herm(rng.normal(size=h.shape) + 1j * rng.normal(size=h.shape))
    y = herm(rng.normal(size=h.shape) + 1j * rng.normal(size=h.shape))
    r = rng.normal(size=raw.shape) + 1j * rng.normal(size=raw.shape)
    max_error = 0.0

    def same(a, b):
        global max_error
        max_error = max(max_error, float(np.max(abs(a - b))))
        np.testing.assert_allclose(a, b, atol=2e-12, rtol=2e-11)

    for cache in (0, 4):
        allmodel = CachedFrozenSigmaResponse(
            central.meta, eps, raw, normalized, weights, central.W, cache_mb=cache
        )
        local = CachedFrozenSigmaResponse(
            central.meta,
            eps[ids],
            raw[ids],
            normalized[ids],
            weights[ids],
            central.W,
            parallel=p,
            cache_mb=cache,
        )
        assert local.window_gap_report() == allmodel.window_gap_report()
        for fixed in (False, True):
            full, fm = allmodel.generalized_response(h, o, r, fixed_number=fixed)
            part, pm = local.generalized_response(
                h[ids], o[ids], r[ids], fixed_number=fixed
            )
            same(part, full[ids])
            same(np.array(pm["dmu"]), np.array(fm["dmu"]))
            fb = allmodel.generalized_adjoint(y, fixed_number=fixed)
            pb = local.generalized_adjoint(y[ids], fixed_number=fixed)
            for a, b in zip(pb, fb):
                same(a, b[ids])
            lhs = p.sum(real_inner(y[ids], part))
            rhs = p.sum(
                sum(real_inner(a, b) for a, b in zip(pb, (h[ids], o[ids], r[ids])))
            )
            same(np.array(lhs), np.array(rhs))
        for a, b in zip(local.projector_gradients(), allmodel.projector_gradients()):
            same(a, b[ids])
    # Global field -> local k matrix -> global field, including a separate
    # replicated core term. Its transpose must sum only the k contributions.
    b = rng.normal(size=(8, 3))
    a = rng.normal(size=(4, 8))
    core = rng.normal(size=(4, 3))
    bop = PairedOperator(3, len(ids), lambda x: b[ids] @ x, lambda y: b[ids].T @ y)
    aop = PairedOperator(
        len(ids), 4, lambda x: a[:, ids] @ x, lambda y: a[:, ids].T @ y
    )
    direct = PairedOperator(
        3, 4, lambda x: (a[:, ids] @ b[ids]) @ x, lambda y: (a[:, ids] @ b[ids]).T @ y
    )
    p.check_transpose(p.global_to_local(bop), output_local=True)
    p.check_transpose(p.local_to_global(aop), input_local=True)
    chain = p.local_to_global(aop) @ p.global_to_local(bop)
    summed = p.sum_operator(direct)
    for op in (chain, summed):
        op.check_transpose(tolerance=2e-13)
        same(op.matvec(np.ones(3)), a @ b @ np.ones(3))
        same(op.rmatvec(np.ones(4)), b.T @ a.T @ np.ones(4))
    if p.root:
        print(
            json.dumps(
                dict(
                    ranks=p.size,
                    passed=True,
                    max_absolute_error=p.reduce(max_error, "max"),
                )
            ),
            flush=True,
        )
    else:
        p.reduce(max_error, "max")
except BaseException:
    import traceback

    traceback.print_exc()
    p.comm.Abort(1)
