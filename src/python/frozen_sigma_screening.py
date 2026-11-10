#!/usr/bin/env python3
"""Solve a frozen-self-energy linear charge response without forming chi.

In material coordinates q=J*n(F*r), the strain response obeys

    q1 = s + A(q1),       A = chi_frozen * K,
    s  = b + chi_frozen * v_affine.

The supplied action must be a directional derivative of ONE charge-map
evaluation at the central structure.  It must include the chosen chemical-
potential constraint and the derivative of the normalized finite-window
embedding.  An impurity vertex is intentionally absent because Sigma is
fixed.  This is not a solver for the fully interacting DMFT susceptibility.

Anderson acceleration is used instead of plain mixing: the latter diverges
for the FeSe response tested here.  We operate on vectors, not dense charge
susceptibility matrices.  The caller chooses a physically appropriate norm
and must separately check finite-difference step and basis convergence.

The optional CLI calls an external action adapter with {input} and {output}
placeholders.  Arrays are ordinary real .npy vectors.  The adapter must write
A(input), NOT s+A(input).  No shell expansion is performed.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import time
from pathlib import Path

import numpy as np


def solve_response(
    source,
    action,
    *,
    initial=None,
    mixing=0.3,
    history=7,
    absolute_tolerance=0.02,
    relative_tolerance=0.0,
    max_iterations=25,
    callback=None,
):
    """Return a verified iterate and its residual history.

    With x=q1, f=s+A(x)-x, differences in successive residuals give the
    least-squares Anderson update

        gamma = argmin || f - DeltaF gamma ||,
        x_new = x + beta*f - (DeltaX+beta*DeltaF) gamma.

    Nearly dependent history columns are discarded by the SVD tolerance.
    This SVD is only a small history solve, not an inverse of the electronic
    response operator.  Nonconvergence is reported, never silently accepted.
    """
    source = np.asarray(source, dtype=float)
    if source.ndim != 1 or not np.all(np.isfinite(source)):
        raise ValueError("source must be a finite real vector")
    if not (0 < mixing <= 1) or history < 0 or max_iterations < 1:
        raise ValueError("invalid screening iteration parameters")
    if absolute_tolerance < 0 or relative_tolerance < 0:
        raise ValueError("response tolerances must be nonnegative")
    x = source.copy() if initial is None else np.array(initial, dtype=float)
    if x.shape != source.shape or not np.all(np.isfinite(x)):
        raise ValueError("initial response has incompatible dimensions")
    threshold = max(absolute_tolerance, relative_tolerance * np.linalg.norm(source))
    xs, fs, records = [], [], []
    for iteration in range(max_iterations):
        started = time.monotonic()
        ax = np.asarray(action(x.copy()), dtype=float)
        if ax.shape != x.shape or not np.all(np.isfinite(ax)):
            raise ValueError("charge-map action returned an invalid response")
        residual = source + ax - x
        norm = float(np.linalg.norm(residual))
        record = {
            "iteration": iteration,
            "residual_norm": norm,
            "seconds": time.monotonic() - started,
        }
        records.append(record)
        if callback is not None:
            callback(x.copy(), residual.copy(), record)
        if norm <= threshold:
            return x, {"converged": True, "threshold": threshold, "iterations": records}
        # Keep the last returned vector tied to an actually evaluated
        # residual even if the iteration budget is exhausted.
        if iteration == max_iterations - 1:
            break
        xs.append(x.copy())
        fs.append(residual.copy())
        xs, fs = xs[-history - 1 :], fs[-history - 1 :]
        next_x = x + mixing * residual
        if len(xs) > 1:
            dx = np.column_stack(np.diff(np.asarray(xs), axis=0))
            df = np.column_stack(np.diff(np.asarray(fs), axis=0))
            coefficients = np.linalg.lstsq(df, residual, rcond=1e-10)[0]
            next_x -= (dx + mixing * df) @ coefficients
        x = next_x
    return x, {"converged": False, "threshold": threshold, "iterations": records}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--action-command", nargs="+", required=True)
    parser.add_argument("--tolerance", type=float, default=0.02)
    parser.add_argument("--max-iterations", type=int, default=25)
    parser.add_argument("--mixing", type=float, default=0.3)
    args = parser.parse_args(argv)
    work = args.work.resolve()
    work.mkdir(parents=True, exist_ok=True)
    counter = 0

    def action(vector):
        nonlocal counter
        folder = work / f"action_{counter:04d}"
        folder.mkdir()  # Refuse to mix a new solve with stale action files.
        counter += 1
        incoming, outgoing = folder / "input.npy", folder / "output.npy"
        np.save(incoming, vector)
        command = [
            arg.replace("{input}", str(incoming)).replace("{output}", str(outgoing))
            for arg in args.action_command
        ]
        with (folder / "action.log").open("w") as log:
            subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT)
        if not outgoing.exists():
            raise ValueError("response adapter did not produce output.npy")
        return np.load(outgoing, allow_pickle=False)

    def checkpoint(vector, residual, row):
        np.save(work / "response.npy", vector)
        np.save(work / "residual.npy", residual)
        print(json.dumps(row), flush=True)

    vector, report = solve_response(
        np.load(args.source, allow_pickle=False),
        action,
        absolute_tolerance=args.tolerance,
        max_iterations=args.max_iterations,
        mixing=args.mixing,
        callback=checkpoint,
    )
    np.save(work / "response.npy", vector)
    (work / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    return 0 if report["converged"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
