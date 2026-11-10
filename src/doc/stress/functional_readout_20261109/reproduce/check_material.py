"""Strained native lattice scalar and explicit+screened readout diagnostics."""

from pathlib import Path
import json
import sys
import time
import resource

import numpy as np
from scipy.linalg import eigvalsh, lu_factor

root = Path(__file__).resolve().parent
if root.name == "reproduce":
    root = root.parent
sys.path[:0] = [str(root / "stage"), "/Users/haule/projects/eDMFT_new/src/python"]
from frozen_sigma_response import FrozenSigmaResponse, herm
from lapw_valence_state import load_native_valence_state
from lapw_potential_core_state import load_potential_core_state
from lapw_screened_state import NativeScreenedState
from lapw_strain_functional import (
    lattice_functional_gradients,
    PreparedStrainFunctional,
)
from lapw_strain_potential import unequal_centered_derivative
from wien_response_density import Density
from stress_response_adjoint import PairedOperator


def logabsdet(a):
    values = []
    for matrix in a:
        lu, _ = lu_factor(matrix, check_finite=True)
        pivots = abs(np.diag(lu))
        if np.any(pivots == 0) or not np.all(np.isfinite(pivots)):
            raise ValueError("singular scalar-oracle matrix")
        values.append(np.log(pivots).sum())
    return np.array(values)


def native_scalar(folder, central, *, mu=None, number=None):
    """Nonlinear native H and U, fixed CENTRAL mu/Sigma/T/quadrature.

    The spectrum makes this oracle invariant to native band phases and avoids
    needing a common endpoint basis. The independent low-rank log determinant
    is evaluated in correlated space. No response gradients enter this code.
    """
    m, records = FrozenSigmaResponse.from_export(folder / "response")
    for name in [
        "temperature",
        "mu",
        "omega",
        "weights",
        "sigma_inf",
        "sigma",
        "dimensions",
        "sym_real",
        "sym_imag",
    ]:
        np.testing.assert_allclose(
            getattr(m.meta, name), getattr(central.meta, name), rtol=0, atol=1e-13
        )
    np.testing.assert_array_equal(m.W, central.W)
    np.testing.assert_array_equal(m.kw, central.kw)
    assert m.eps.shape == central.eps.shape
    T = central.meta.temperature
    if mu is None:
        mu = central.meta.mu
    elif number is None:
        raise ValueError("a changed mu requires the Legendre particle number")
    total = 0.0
    for k, u in enumerate(m.uc):
        eps = m.eps[k]
        static = np.diag(eps - mu).astype(complex) + u @ central.sinf @ u.conj().T
        eigenvalues = eigvalsh(herm(static))
        value = -T * np.logaddexp(0, -eigenvalues / T).sum()
        for j in range(0, len(central.meta.omega), 64):
            sl = slice(j, j + 64)
            z = mu + 1j * central.meta.omega[sl]
            overlap = np.einsum(
                "ni,wn,nj->wij", u.conj(), 1 / (z[:, None] - eps), u, optimize=True
            )
            ls = logabsdet(np.eye(u.shape[1]) - central.sinf @ overlap)
            ld = logabsdet(np.eye(u.shape[1]) - central.sigma[sl] @ overlap)
            value += 2 * T * np.dot(central.meta.weights[sl], ls - ld)
        total += 2 * central.kw[k] * value
    if number is not None:
        total += (mu - central.meta.mu) * number
    return float(total)


def run(case, mesh="three_k", solve=True):
    old = Path("/private/tmp/edmft_full_k_response_20261108")
    source = (old / mesh if mesh == "full_bz" else old) / case / "charge/0"
    native = Path("/private/tmp/edmft_screened_native_20261108") / mesh / case
    previous = Path("/private/tmp/edmft_strain_charge_20261109") / mesh / case
    output = root / mesh / case
    output.mkdir(parents=True, exist_ok=True)
    parts = np.load(previous / "analytic.npz")
    with load_native_valence_state(
        source,
        case,
        complex_vectors=(case == "SrVO3"),
        energy_policy="fixed",
        matrix_quadrature="wien2k24",
        library="/Users/haule/dbin/liblapw_response.so",
    ) as v:
        print(case, mesh, "prepared", v.preparation_seconds, flush=True)
        pc = load_potential_core_state(
            "/Users/haule/projects/eDMFT_new/src/doc/stress/charge_closure_20261108/"
            + case
            + "_potential_core.npz"
        )
        core = Density(
            Path("/private/tmp/edmft_xc_core_response_20261108")
            / case
            / "hxc_core/0"
            / (case + ".clmcor"),
            source / (case + ".struct"),
        )
        state = NativeScreenedState(
            v,
            pc,
            core,
            target_charge={"SrVO3": 85, "FeSe": 120}[case],
            fixed_number=True,
        )
        start = time.monotonic()
        prepared = lattice_functional_gradients(v.complement)
        gradients, coords = prepared
        report = dict(
            case=case,
            mesh=mesh,
            kpoints=len(v.records),
            central=state.metadata,
            valence_preparation_seconds=v.preparation_seconds,
            gradient_seconds=time.monotonic() - start,
            scope="regularized lattice scalar, not full Klein or total stress",
            sources=v.fingerprints,
            components={},
        )
        joint = parts["joint"]
        analytic = float(gradients["lattice"] @ joint)
        report["explicit_lattice_derivative_Ry"] = analytic
        center = native_scalar(source, v.model)
        report["native_steps"] = {}
        for step in [0.00025, 0.000125]:
            folders = [native / f"source_{label}_{step}" for label in ["m", "p"]]
            offsets = [
                json.loads((p / "done.json").read_text())["actual_scale"] - 1
                for p in folders
            ]
            values = [native_scalar(p, v.model) for p in folders]
            fd = float(
                unequal_centered_derivative(values[0], center, values[1], *offsets)
            )
            row = dict(
                offsets=offsets,
                values_Ry=values,
                fd_Ry=fd,
                error_Ry=fd - analytic,
                relative_error=abs(fd - analytic) / abs(analytic),
            )
            report["native_steps"][str(step)] = row
            print(case, mesh, "native lattice", step, row, flush=True)
        steps = list(report["native_steps"].values())
        ratio = (
            (steps[0]["offsets"][1] - steps[0]["offsets"][0])
            / (steps[1]["offsets"][1] - steps[1]["offsets"][0])
        ) ** 2
        extrap = (ratio * steps[1]["fd_Ry"] - steps[0]["fd_Ry"]) / (ratio - 1)
        report["richardson"] = dict(
            fd_Ry=extrap,
            error_Ry=extrap - analytic,
            relative_error=abs(extrap - analytic) / abs(analytic),
        )
        geo, pot = parts["band_geometry"][:, None], parts["potential_source"][:, None]
        q = np.load(previous / "analytic_density.npy")[:, None]
        s = np.load(previous / "analytic_source.npy")[:, None]
        residual = s[:, 0] - q[:, 0] + state.jacobian.matvec(q[:, 0])
        report["forward_residual_norm"] = float(np.linalg.norm(residual))
        saved = {}
        for component in ["projector", "lattice"]:
            obs = PreparedStrainFunctional(
                state, component=component, gradients=prepared
            )
            values, row = obs.forward_contractions(geo, pot, q)
            row["forward_Ry"] = values.tolist()
            row["explicit_joint_Ry"] = float(obs.joint_gradient @ joint)
            explicit = (
                row["geometry_derivatives"][0]
                + row["explicit_potential_derivatives"][0]
            )
            assert abs(explicit - row["explicit_joint_Ry"]) < 2e-9
            saved.update(
                {
                    component + "_" + key: value
                    for key, value in dict(
                        joint=obs.joint_gradient,
                        potential=obs.potential_gradient,
                        charge=obs.charge_gradient,
                    ).items()
                }
            )
            if solve:
                initial = None
                if mesh == "full_bz":
                    # Charge coordinates are identical; only the tangent map
                    # and readout differ. This is just an initial guess, never
                    # a substitute for the final full-mesh true-residual check.
                    warm = root / "three_k" / case / "covectors.npz"
                    initial = np.load(warm)[component + "_adjoint"]
                    row["adjoint_initial"] = str(warm)
                if mesh == "three_k" and component == "projector":
                    initial = np.load(
                        Path("/private/tmp/edmft_projector_functional_20261108")
                        / case
                        / "projector_covectors.npz"
                    )["adjoint"]
                base_j = state.screened.jacobian
                action_seconds = []

                def transpose(x):
                    began = time.monotonic()
                    answer = base_j.rmatvec(x)
                    action_seconds.append(time.monotonic() - began)
                    print(
                        case,
                        mesh,
                        component,
                        "adjoint action",
                        len(action_seconds),
                        "seconds",
                        action_seconds[-1],
                        flush=True,
                    )
                    return answer

                state.screened.jacobian = PairedOperator(
                    base_j.input_size,
                    base_j.output_size,
                    base_j.matvec,
                    transpose,
                    "timed exact native Jacobian",
                )
                try:
                    actual, lam, adjoint = obs.adjoint_contractions(
                        geo,
                        pot,
                        s,
                        initial=initial,
                        rtol=2e-9,
                        atol=1e-11,
                        restart=30,
                        max_cycles=20,
                    )
                finally:
                    state.screened.jacobian = base_j
                row["adjoint_action_seconds"] = action_seconds
                bound = float(
                    np.linalg.norm(lam) * np.linalg.norm(residual)
                    + np.linalg.norm(q) * adjoint["adjoint"]["residual_norm"]
                )
                row.update(
                    adjoint=adjoint,
                    adjoint_Ry=actual.tolist(),
                    forward_adjoint_error_Ry=float(abs(actual[0] - values[0])),
                    residual_bound_Ry=bound,
                )
                assert abs(actual[0] - values[0]) < 1.001 * bound + 2e-10
                saved[component + "_adjoint"] = lam
            report["components"][component] = row
            print(case, mesh, component, json.dumps(row), flush=True)
        direct = {}
        for key in [
            "geometry_derivatives",
            "explicit_potential_derivatives",
            "induced_derivatives",
            "forward_Ry",
        ]:
            direct[key] = (
                np.array(report["components"]["lattice"][key])
                - np.array(report["components"]["projector"][key])
            ).tolist()
        report["components"]["direct"] = direct
        report["peak_rss_bytes"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        np.savez_compressed(output / "covectors.npz", **saved)
        (output / "REPORT.json").write_text(json.dumps(report, indent=2) + "\n")
        print(case, mesh, "DONE", report["richardson"], flush=True)


if __name__ == "__main__":
    run(
        sys.argv[1],
        sys.argv[2] if len(sys.argv) > 2 else "three_k",
        "--no-solve" not in sys.argv,
    )
