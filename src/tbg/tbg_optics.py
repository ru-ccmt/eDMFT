#!/usr/bin/env python
"""Optical conductivity prototype for the TBG DMFT data.

This first version evaluates the real-axis bubble in the spectral form

    sigma_ab(Omega) = pi/Omega sum_k int dw [f(w)-f(w+Omega)]
                      Tr[v_a A(w+Omega) v_b A(w)]

where A = i (G^R - G^A)/(2*pi).  The self-energy is diagonal and scalar in
the local four-orbital space, but becomes low-rank in the retained band basis:

    Sigma_band(k,w) = sum_site sigma_site(w) U_site(k)^H U_site(k).

The production path uses this low-rank structure through Woodbury, with an
optional SVD compression of each projector.  A full-inversion path is kept for
debugging conventions on small k/frequency windows.
"""

from __future__ import annotations

import argparse
import os
import re
from dataclasses import dataclass
from typing import Iterable

import numpy as np


try:
    from numba import njit
except Exception:  # pragma: no cover - numba is optional
    njit = None


PI = np.pi
E2_OVER_H_SI = 3.87404586493e-5  # Siemens


@dataclass
class Channel:
    name: str
    projector_key: str
    sigma_path: str
    sigma_column: int = 0
    zero: bool = False


def fermi(x: np.ndarray, beta: float) -> np.ndarray:
    bx = np.clip(beta * x, -700.0, 700.0)
    return 1.0 / (np.exp(bx) + 1.0)


def minus_fermi_derivative(x: np.ndarray, beta: float) -> np.ndarray:
    bx = np.clip(beta * x, -700.0, 700.0)
    f = 1.0 / (np.exp(bx) + 1.0)
    return beta * f * (1.0 - f)


def read_header_value(path: str, key: str) -> float | None:
    try:
        with open(path, "r") as fh:
            first = fh.readline()
    except OSError:
        return None
    if not first.startswith("#"):
        return None
    pattern = re.compile(rf"(?:^|\s){re.escape(key)}=([^\s,]+)")
    match = pattern.search(first[1:])
    if match is None:
        return None
    try:
        return float(match.group(1))
    except ValueError:
        return None


def read_sigma(path: str, sigma_column: int = 0) -> tuple[np.ndarray, np.ndarray]:
    data = np.loadtxt(path, comments="#")
    if data.ndim != 2 or data.shape[1] < 3:
        raise ValueError(f"{path}: expected columns w ReSigma ImSigma [Re Im ...]")
    pair0 = 1 + 2 * sigma_column
    pair1 = pair0 + 1
    if pair1 >= data.shape[1]:
        npairs = (data.shape[1] - 1) // 2
        raise ValueError(f"{path}: sigma_column={sigma_column} but only {npairs} complex columns exist")
    w = data[:, 0].astype(float)
    sig = data[:, pair0] + 1j * data[:, pair1]
    order = np.argsort(w)
    return w[order], sig[order]


def read_sigma_temperature(path: str) -> float | None:
    return read_header_value(path, "T")


def read_beta(workdir: str, channels: list[Channel], explicit_beta: float | None) -> float:
    if explicit_beta is not None:
        return explicit_beta
    for ch in channels:
        if ch.zero:
            continue
        spath = ch.sigma_path if os.path.isabs(ch.sigma_path) else os.path.join(workdir, ch.sigma_path)
        temp = read_sigma_temperature(spath)
        if temp is not None and temp > 0:
            return 1.0 / temp
    raise ValueError("Could not determine temperature from Sig.out header; pass --beta explicitly")


def read_top_level_param(workdir: str, name: str) -> float | None:
    for fname in ("tbg_params.dat", "tbg_params.py"):
        path = os.path.join(workdir, fname)
        if not os.path.exists(path):
            continue
        pattern = re.compile(rf"^\s*{re.escape(name)}\s*=\s*([-+0-9.eEdD]+)")
        with open(path, "r") as fh:
            for line in fh:
                match = pattern.match(line)
                if match:
                    return float(match.group(1).replace("D", "E").replace("d", "e"))
    return None


def first_sigma_path(workdir: str, channels: list[Channel]) -> str | None:
    for ch in channels:
        if ch.zero:
            continue
        return ch.sigma_path if os.path.isabs(ch.sigma_path) else os.path.join(workdir, ch.sigma_path)
    return None


def compute_formula_mu(workdir: str, channels: list[Channel]) -> tuple[float, dict[str, float]]:
    spath = first_sigma_path(workdir, channels)
    if spath is None:
        raise ValueError("Can not compute formula mu without a nonzero self-energy channel")
    sig_mu = read_header_value(spath, "mu")
    u0 = read_top_level_param(workdir, "U0")
    nf = read_top_level_param(workdir, "Nf")
    missing = [name for name, value in (("mu", sig_mu), ("U0", u0), ("Nf", nf)) if value is None]
    if missing:
        raise ValueError("Can not compute formula mu; missing " + ", ".join(missing))
    ef_formula = sig_mu - u0 * (nf - 0.5)
    return ef_formula, {"sig_mu": sig_mu, "U0": u0, "Nf": nf}


def resolve_mu(workdir: str, channels: list[Channel], explicit_mu: float | None, mu_source: str, warn_tol: float) -> float:
    if explicit_mu is not None:
        return explicit_mu

    ef_path = os.path.join(workdir, "EF.dat")
    ef_mu = None
    if os.path.exists(ef_path):
        vals = np.loadtxt(ef_path, comments="#", ndmin=1)
        ef_mu = float(np.ravel(vals)[0])

    formula_mu = None
    formula_info = None
    try:
        formula_mu, formula_info = compute_formula_mu(workdir, channels)
    except ValueError as exc:
        if mu_source == "formula":
            raise
        print(f"Warning: could not check EF.dat against formula: {exc}", flush=True)

    if ef_mu is not None and formula_mu is not None:
        diff = ef_mu - formula_mu
        msg = (
            f"EF check: EF.dat={ef_mu:.12g} formula=mu-U0*(Nf-0.5)={formula_mu:.12g} "
            f"diff={diff:.3e} (mu={formula_info['sig_mu']:.12g}, U0={formula_info['U0']:.12g}, Nf={formula_info['Nf']:.12g})"
        )
        if abs(diff) > warn_tol:
            print("Warning: " + msg, flush=True)
        else:
            print(msg, flush=True)

    if mu_source == "formula":
        if formula_mu is None:
            raise ValueError("Can not use --mu-source formula because formula value is unavailable")
        return formula_mu
    if mu_source == "ef":
        if ef_mu is not None:
            return ef_mu
        if formula_mu is not None:
            print("Warning: EF.dat is missing; using formula value for mu", flush=True)
            return formula_mu
        return 0.0
    raise ValueError(f"unknown mu_source {mu_source!r}")


def interp_complex(x: np.ndarray, xp: np.ndarray, fp: np.ndarray) -> np.ndarray:
    return np.interp(x, xp, fp.real) + 1j * np.interp(x, xp, fp.imag)


def parse_channel(spec: str) -> Channel:
    """Parse NAME=SIGMA or NAME=SIGMA:COLUMN.

    NAME is also the key in Uproject.npz, for example UAA or UAB.
    """
    if "=" not in spec:
        raise ValueError(f"channel must be KEY=SIGMA_PATH[:COLUMN], got {spec!r}")
    key, rhs = spec.split("=", 1)
    if rhs.lower() in {"0", "zero", "none"}:
        return Channel(name=key, projector_key=key, sigma_path="zero", zero=True)
    if ":" in rhs:
        sigma_path, col = rhs.rsplit(":", 1)
        sigma_column = int(col)
    else:
        sigma_path = rhs
        sigma_column = 0
    return Channel(name=key, projector_key=key, sigma_path=sigma_path, sigma_column=sigma_column)


def default_channels(workdir: str) -> list[Channel]:
    channels = []
    p0 = os.path.join(workdir, "imp.0", "Sig.out")
    if os.path.exists(p0):
        channels.append(Channel("UAA", "UAA", p0))
    if not channels:
        raise FileNotFoundError("No default self-energy found: expected imp.0/Sig.out")
    return channels


def load_inputs(args: argparse.Namespace):
    workdir = os.path.abspath(args.workdir)
    ene = np.load(os.path.join(workdir, args.ene))

    vpath = os.path.join(workdir, args.velocity)
    if not os.path.exists(vpath):
        raise FileNotFoundError(
            f"{vpath} does not exist. Generate it by setting par_bs['ComputeVelocity']=True "
            "in tbg_params.py and rerunning tbg_dft.py on the DMFT k mesh."
        )
    vz = np.load(vpath)
    vx = vz["vx"]
    vy = vz["vy"]

    pz = np.load(os.path.join(workdir, args.projectors))
    channels = [parse_channel(c) for c in args.channel] if args.channel else default_channels(workdir)

    sigmas = []
    projectors = []
    for ch in channels:
        if ch.projector_key not in pz.files:
            raise KeyError(f"{args.projectors} has no projector key {ch.projector_key!r}; keys={pz.files}")
        if ch.zero:
            sw = None
            ssig = None
        else:
            spath = ch.sigma_path if os.path.isabs(ch.sigma_path) else os.path.join(workdir, ch.sigma_path)
            sw, ssig = read_sigma(spath, ch.sigma_column)
        sigmas.append((sw, ssig))
        projectors.append(pz[ch.projector_key])

    nk, nb = ene.shape
    if vx.shape != (nk, nb, nb) or vy.shape != (nk, nb, nb):
        raise ValueError(f"velocity shapes must be {(nk, nb, nb)}, got vx={vx.shape}, vy={vy.shape}")
    for ch, up in zip(channels, projectors):
        if up.shape[1:] != (nk, nb):
            raise ValueError(f"projector {ch.projector_key} shape {up.shape} incompatible with ene {ene.shape}")

    return workdir, ene, vx, vy, channels, sigmas, projectors


def optics_unit_prefactor(workdir: str, units: str) -> tuple[float, str]:
    bpath = os.path.join(workdir, "basis.npz")
    if not os.path.exists(bpath):
        raise FileNotFoundError(f"{bpath} is needed for physical optics units")
    basis = np.load(bpath)
    if "Rs" not in basis.files:
        raise KeyError(f"{bpath} does not contain Rs")
    area = abs(float(np.linalg.det(basis["Rs"])))
    if area <= 0:
        raise ValueError(f"Invalid moire cell area from {bpath}: {area}")
    pref_e2h = 2.0 * PI / area
    if units == "e2h":
        return pref_e2h, "e^2/h"
    if units == "si":
        return pref_e2h * E2_OVER_H_SI, "S"
    raise ValueError(f"unknown units {units!r}")


def estimate_sigma_spacing(sigmas: list[tuple[np.ndarray | None, np.ndarray | None]]) -> float:
    spacings = []
    for sw, _ in sigmas:
        if sw is None or len(sw) < 2:
            continue
        dw = np.diff(np.sort(sw))
        dw = dw[dw > 0]
        if len(dw):
            spacings.append(float(np.min(dw)))
    if not spacings:
        raise ValueError("Can not estimate self-energy mesh spacing")
    return min(spacings)


def local_sigma_spacing(sigmas: list[tuple[np.ndarray | None, np.ndarray | None]], wmin: float, wmax: float) -> float:
    """Smallest Sig.out spacing in the frequency interval relevant for one Omega."""
    spacings = []
    for sw, _ in sigmas:
        if sw is None or len(sw) < 2:
            continue
        ws = np.sort(sw)
        lo = max(0, np.searchsorted(ws, wmin, side="left") - 1)
        hi = min(len(ws), np.searchsorted(ws, wmax, side="right") + 1)
        if hi - lo >= 2:
            dw = np.diff(ws[lo:hi])
            dw = dw[dw > 0]
            if len(dw):
                spacings.append(float(np.min(dw)))
    if spacings:
        return min(spacings)
    return estimate_sigma_spacing(sigmas)


def adaptive_frequency_mesh(
    omega: float,
    temperature: float,
    sigmas: list[tuple[np.ndarray | None, np.ndarray | None]],
    thermal_window: float,
    points_per_t: float,
    dc_min_points: int,
    points_per_omega: float,
    omega_mesh_crossover_t: float,
    sigma_spacing_factor: float,
    dw_min: float,
    dw_max: float | None,
) -> np.ndarray:
    pad = thermal_window * temperature
    wmin = -omega - pad
    wmax = pad
    if omega <= 0:
        raise ValueError("omega must be positive")

    candidates = []
    if temperature > 0 and points_per_t > 0:
        candidates.append(temperature / points_per_t)
    use_omega_spacing = points_per_omega > 0
    if temperature > 0 and omega_mesh_crossover_t > 0:
        use_omega_spacing = use_omega_spacing and omega > omega_mesh_crossover_t * temperature
    if use_omega_spacing:
        candidates.append(omega / points_per_omega)

    local_dw = local_sigma_spacing(sigmas, wmin, wmax)
    if sigma_spacing_factor > 0:
        candidates.append(sigma_spacing_factor * local_dw)

    edge_dw = min(candidates) if candidates else local_dw
    edge_dw = max(edge_dw, dw_min)
    if dw_max is not None:
        edge_dw = min(edge_dw, dw_max)

    points = [wmin, wmax, -omega, 0.0]

    def add_linear(a: float, b: float, dw: float) -> None:
        aa = max(wmin, a)
        bb = min(wmax, b)
        if bb <= aa:
            return
        n = max(2, int(np.ceil((bb - aa) / dw)) + 1)
        points.extend(np.linspace(aa, bb, n))

    add_linear(-omega - pad, -omega + pad, edge_dw)
    add_linear(-pad, pad, edge_dw)

    for sw, _ in sigmas:
        if sw is None:
            continue
        ws = np.asarray(sw)
        keep_m = (ws >= wmin) & (ws <= wmax)
        points.extend(ws[keep_m])

        shifted = ws - omega
        keep_p = (shifted >= wmin) & (shifted <= wmax)
        points.extend(shifted[keep_p])

    mesh = np.array(points, dtype=float)
    mesh = mesh[(mesh >= wmin - 1e-14) & (mesh <= wmax + 1e-14)]
    mesh.sort()
    if len(mesh) == 0:
        return np.linspace(wmin, wmax, 3)

    cleaned = [mesh[0]]
    for x in mesh[1:]:
        if x - cleaned[-1] > max(1e-12, 0.1 * dw_min):
            cleaned.append(x)
    return np.array(cleaned)


def dc_frequency_mesh(
    temperature: float,
    sigmas: list[tuple[np.ndarray | None, np.ndarray | None]],
    thermal_window: float,
    points_per_t: float,
    dc_min_points: int,
    sigma_spacing_factor: float,
    dw_min: float,
    dw_max: float | None,
) -> np.ndarray:
    pad = thermal_window * temperature
    wmin = -pad
    wmax = pad
    candidates = []
    if temperature > 0 and points_per_t > 0:
        candidates.append(temperature / points_per_t)
    if dc_min_points > 1:
        candidates.append((wmax - wmin) / (dc_min_points - 1))
    local_dw = local_sigma_spacing(sigmas, wmin, wmax)
    if sigma_spacing_factor > 0:
        candidates.append(sigma_spacing_factor * local_dw)

    edge_dw = min(candidates) if candidates else local_dw
    edge_dw = max(edge_dw, dw_min)
    if dw_max is not None:
        edge_dw = min(edge_dw, dw_max)

    points = [wmin, wmax, 0.0]
    n = max(2, int(np.ceil((wmax - wmin) / edge_dw)) + 1)
    points.extend(np.linspace(wmin, wmax, n))
    for sw, _ in sigmas:
        if sw is None:
            continue
        ws = np.asarray(sw)
        keep = (ws >= wmin) & (ws <= wmax)
        points.extend(ws[keep])

    mesh = np.array(points, dtype=float)
    mesh = mesh[(mesh >= wmin - 1e-14) & (mesh <= wmax + 1e-14)]
    mesh.sort()
    cleaned = [mesh[0]]
    for x in mesh[1:]:
        if x - cleaned[-1] > max(1e-12, 0.1 * dw_min):
            cleaned.append(x)
    return np.array(cleaned)


def compress_projector(U: np.ndarray, svd_tol: float) -> np.ndarray:
    """Return X such that U^H U ~= X X^H, with X shaped (nband, rank)."""
    if svd_tol < 0:
        return U.conj().T.copy()
    _, s, vh = np.linalg.svd(U, full_matrices=False)
    if s.size == 0:
        return np.zeros((U.shape[1], 0), dtype=np.complex128)
    keep = s > svd_tol * s[0]
    if not np.any(keep):
        return np.zeros((U.shape[1], 0), dtype=np.complex128)
    v = vh[keep, :].conj().T
    return v * s[keep][None, :]


def green_woodbury(
    z: complex,
    eps: np.ndarray,
    sigma_values: np.ndarray,
    xcols: np.ndarray,
    xchan: np.ndarray,
) -> np.ndarray:
    nb = eps.size
    dinv = 1.0 / (z - eps)
    g = np.diag(dinv.astype(np.complex128))
    if xcols.shape[1] == 0:
        return g

    coeff = sigma_values[xchan]
    xd = xcols.conj().T * dinv[None, :]
    ksmall = xd @ xcols
    mat = np.eye(xcols.shape[1], dtype=np.complex128) - coeff[:, None] * ksmall
    rhs = coeff[:, None] * xd
    corr = (dinv[:, None] * xcols) @ np.linalg.solve(mat, rhs)
    return g + corr


def green_direct(z: complex, eps: np.ndarray, sigma_values: np.ndarray, Us: list[np.ndarray]) -> np.ndarray:
    h = np.diag(z - eps).astype(np.complex128)
    for sig, U in zip(sigma_values, Us):
        h -= sig * (U.conj().T @ U)
    return np.linalg.inv(h)


def spectral_matrix(g: np.ndarray) -> np.ndarray:
    return 1j * (g - g.conj().T) / (2.0 * PI)


def trace_vertices(va: np.ndarray, ap: np.ndarray, vb: np.ndarray, am: np.ndarray) -> np.ndarray:
    """Vectorized trace of va @ ap[t] @ vb @ am[t]."""
    return np.einsum("ab,tbc,cd,tda->t", va, ap, vb, am, optimize=True)


def build_spectral_for_k(
    eps: np.ndarray,
    mu: float,
    eta: float,
    w: np.ndarray,
    sigmas_on_w: list[np.ndarray],
    Uk: list[np.ndarray],
    mode: str,
    svd_tol: float,
) -> np.ndarray:
    nw = len(w)
    nb = len(eps)
    amat = np.empty((nw, nb, nb), dtype=np.complex128)

    if mode == "woodbury":
        xs = []
        xchan = []
        for ic, U in enumerate(Uk):
            x = compress_projector(U, svd_tol)
            xs.append(x)
            xchan.extend([ic] * x.shape[1])
        xcols = np.concatenate(xs, axis=1) if xs else np.zeros((nb, 0), dtype=np.complex128)
        xchan = np.asarray(xchan, dtype=np.int64)

    for iw, ww in enumerate(w):
        sigv = np.array([sig[iw] for sig in sigmas_on_w], dtype=np.complex128)
        z = ww + mu + 1j * eta
        if mode == "woodbury":
            g = green_woodbury(z, eps, sigv, xcols, xchan)
        elif mode == "direct":
            g = green_direct(z, eps, sigv, Uk)
        else:
            raise ValueError(f"unknown mode {mode!r}")
        amat[iw] = spectral_matrix(g)
    return amat


def conductivity_adaptive(
    ene: np.ndarray,
    vx: np.ndarray,
    vy: np.ndarray,
    projectors: list[np.ndarray],
    sigmas: list[tuple[np.ndarray | None, np.ndarray | None]],
    beta: float,
    mu: float,
    eta: float,
    omega_max: float,
    omega_min: float,
    omega_step: float,
    mode: str,
    svd_tol: float,
    directions: Iterable[str],
    spin: float,
    prefactor: float,
    k_limit: int | None,
    thermal_window: float,
    points_per_t: float,
    dc_min_points: int,
    points_per_omega: float,
    omega_mesh_crossover_t: float,
    sigma_spacing_factor: float,
    dw_min: float,
    dw_max: float | None,
    include_dc: bool,
    dc_only: bool,
) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    if dc_only:
        include_dc = True
        omegas = np.array([], dtype=float)
    else:
        omegas = np.arange(omega_min, omega_max + 0.5 * omega_step, omega_step)
        omegas = omegas[omegas > 0]
    if len(omegas) == 0 and not include_dc:
        raise ValueError("No positive optical frequencies selected")

    dirs = list(directions)
    all_omegas = np.concatenate(([0.0], omegas)) if include_dc else omegas
    accum = {d: np.zeros(len(all_omegas), dtype=np.complex128) for d in dirs}
    nk_total = ene.shape[0] if k_limit is None else min(k_limit, ene.shape[0])

    offset = 0
    if include_dc:
        wdc = dc_frequency_mesh(
            1.0 / beta,
            sigmas,
            thermal_window,
            points_per_t,
            dc_min_points,
            sigma_spacing_factor,
            dw_min,
            dw_max,
        )
        sigdc = []
        for sw, ssig in sigmas:
            if sw is None:
                sigdc.append(np.zeros_like(wdc, dtype=np.complex128))
            else:
                sigdc.append(interp_complex(wdc, sw, ssig))
        weights = np.empty_like(wdc)
        weights[0] = 0.5 * (wdc[1] - wdc[0])
        weights[-1] = 0.5 * (wdc[-1] - wdc[-2])
        weights[1:-1] = 0.5 * (wdc[2:] - wdc[:-2])
        fd_weight = minus_fermi_derivative(wdc, beta)
        print(
            f"Omega 0/DC nw={len(wdc)} dw~{np.mean(np.diff(wdc)):.4g} "
            f"window=[{wdc[0]:.4g},{wdc[-1]:.4g}]",
            flush=True,
        )
        for ik in range(nk_total):
            Uk = [U[:, ik, :] for U in projectors]
            amat = build_spectral_for_k(ene[ik], mu, eta, wdc, sigdc, Uk, mode, svd_tol)
            va = {"xx": vx[ik], "yy": vy[ik], "xy": vx[ik], "yx": vy[ik]}
            vb = {"xx": vx[ik], "yy": vy[ik], "xy": vy[ik], "yx": vx[ik]}
            for d in dirs:
                if d == "trace":
                    vals = trace_vertices(vx[ik], amat, vx[ik], amat)
                    vals += trace_vertices(vy[ik], amat, vy[ik], amat)
                else:
                    vals = trace_vertices(va[d], amat, vb[d], amat)
                accum[d][0] += np.dot(weights * fd_weight, vals)
        offset = 1

    for io, omega in enumerate(omegas):
        wint = adaptive_frequency_mesh(
            omega,
            1.0 / beta,
            sigmas,
            thermal_window,
            points_per_t,
            points_per_omega,
            omega_mesh_crossover_t,
            sigma_spacing_factor,
            dw_min,
            dw_max,
        )
        wplus = wint + omega
        sigm = []
        sigp = []
        for sw, ssig in sigmas:
            if sw is None:
                sigm.append(np.zeros_like(wint, dtype=np.complex128))
                sigp.append(np.zeros_like(wint, dtype=np.complex128))
            else:
                sigm.append(interp_complex(wint, sw, ssig))
                sigp.append(interp_complex(wplus, sw, ssig))

        fdiff = fermi(wint, beta) - fermi(wplus, beta)
        weights = np.empty_like(wint)
        weights[0] = 0.5 * (wint[1] - wint[0])
        weights[-1] = 0.5 * (wint[-1] - wint[-2])
        weights[1:-1] = 0.5 * (wint[2:] - wint[:-2])
        print(
            f"Omega {io+1}/{len(omegas)}={omega:.6g} nw={len(wint)} "
            f"dw~{np.mean(np.diff(wint)):.4g} window=[{wint[0]:.4g},{wint[-1]:.4g}]",
            flush=True,
        )

        for ik in range(nk_total):
            Uk = [U[:, ik, :] for U in projectors]
            am = build_spectral_for_k(ene[ik], mu, eta, wint, sigm, Uk, mode, svd_tol)
            ap = build_spectral_for_k(ene[ik], mu, eta, wplus, sigp, Uk, mode, svd_tol)
            va = {"xx": vx[ik], "yy": vy[ik], "xy": vx[ik], "yx": vy[ik]}
            vb = {"xx": vx[ik], "yy": vy[ik], "xy": vy[ik], "yx": vx[ik]}
            for d in dirs:
                if d == "trace":
                    vals = trace_vertices(vx[ik], ap, vx[ik], am)
                    vals += trace_vertices(vy[ik], ap, vy[ik], am)
                else:
                    vals = trace_vertices(va[d], ap, vb[d], am)
                accum[d][io + offset] += np.dot(weights * fdiff, vals)

    norm = prefactor * spin * PI / nk_total
    for d in dirs:
        if include_dc:
            accum[d][0] *= norm
            accum[d][offset:] = norm * accum[d][offset:] / omegas
        else:
            accum[d] = norm * accum[d] / omegas
    return all_omegas, accum


def validate_woodbury(
    ene: np.ndarray,
    projectors: list[np.ndarray],
    sigmas: list[tuple[np.ndarray | None, np.ndarray | None]],
    mu: float,
    eta: float,
    svd_tol: float,
) -> None:
    nonzero = [(sw, ssig) for sw, ssig in sigmas if sw is not None]
    if not nonzero:
        raise ValueError("Need at least one nonzero self-energy channel for validation")
    wref = nonzero[0][0]
    test_ws = [float(wref[len(wref) // 2]), float(wref[len(wref) // 3])]
    tests = [(0, test_ws[0]), (min(ene.shape[0] - 1, ene.shape[0] // 3), test_ws[1])]
    for ik, ww in tests:
        Uk = [U[:, ik, :] for U in projectors]
        xs = []
        xchan = []
        for ic, U in enumerate(Uk):
            x = compress_projector(U, svd_tol)
            xs.append(x)
            xchan.extend([ic] * x.shape[1])
        xcols = np.concatenate(xs, axis=1) if xs else np.zeros((ene.shape[1], 0), dtype=np.complex128)
        xchan = np.asarray(xchan, dtype=np.int64)
        sigv = np.array([
            0.0j if sw is None else interp_complex(np.array([ww]), sw, ssig)[0]
            for sw, ssig in sigmas
        ], dtype=np.complex128)
        z = ww + mu + 1j * eta
        gd = green_direct(z, ene[ik], sigv, Uk)
        gw = green_woodbury(z, ene[ik], sigv, xcols, xchan)
        err = np.max(np.abs(gd - gw))
        rel = err / max(1.0, np.max(np.abs(gd)))
        print(f"validate ik={ik} w={ww:.8g} max|G_direct-G_woodbury|={err:.6e} rel={rel:.6e}")


def write_output(path: str, omegas: np.ndarray, accum: dict[str, np.ndarray], unit_label: str) -> None:
    keys = list(accum.keys())
    sigma_label = "Re_sigma[e^2/h]" if unit_label == "e^2/h" else "Re_sigma[1/Ohm]"
    header_parts = []
    for k in keys:
        if k == "trace":
            header_parts.append(sigma_label)
        else:
            header_parts.append(f"Re_sigma_{k}[{unit_label}]")
    header = "omega[50meV] " + " ".join(header_parts)
    cols = [omegas]
    for k in keys:
        if k == "trace":
            cols.append(accum[k].real)
        else:
            cols.append(accum[k].real)
    np.savetxt(path, np.column_stack(cols), header=header)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Compute TBG optical conductivity from band energies, velocity.npz, "
            "Uproject.npz, and imp.0/Sig.out. Default output is the isotropic "
            "trace Re(sigma_xx+sigma_yy) in Siemens (1/Ohm), with temperature "
            "read from Sig.out and chemical potential read from EF.dat."
        )
    )
    parser.add_argument("--workdir", default=".", help="TBG run directory (default: current directory)")
    parser.add_argument("--ene", default="ene.npy", help="band energies file (default: ene.npy)")
    parser.add_argument("--velocity", default="velocity.npz", help="velocity file with vx, vy (default: velocity.npz)")
    parser.add_argument("--projectors", default="Uproject.npz", help="projector file (default: Uproject.npz)")
    parser.add_argument(
        "--channel",
        action="append",
        help="projector/self-energy channel KEY=SIGMA_PATH[:COLUMN]; can be repeated. Default: UAA=imp.0/Sig.out. Use UAB=zero for an uncorrelated AB channel.",
    )
    parser.add_argument("--sigma-scale", type=float, default=1.0, help="multiply all nonzero self-energies by this factor; use 0 for a noninteracting test")
    parser.add_argument("--mu", type=float, default=None, help="override chemical potential used in G (default: read EF.dat)")
    parser.add_argument("--mu-source", choices=["ef", "formula"], default="ef", help="chemical potential source: ef reads EF.dat (default); formula uses mu-U0*(Nf-0.5)")
    parser.add_argument("--ef-check-tol", type=float, default=1e-3, help="warning tolerance for EF.dat versus formula check (default: 1e-3)")
    parser.add_argument("--beta", type=float, default=None, help="inverse temperature override (default: read T=... from Sig.out)")
    parser.add_argument("--eta", type=float, default=0.0, help="extra positive broadening in G (default: 0)")
    parser.add_argument("--omega-max", type=float, default=5.0, help="maximum positive optical frequency in units of 50 meV (default: 5)")
    parser.add_argument("--omega-min", type=float, default=None, help="minimum positive optical frequency (default: omega-step)")
    parser.add_argument("--omega-step", type=float, default=None, help="spacing of external optical frequencies (default: min local Sig.out spacing and T/2)")
    parser.add_argument("--include-dc", action="store_true", help="also write Omega=0 conductivity using -df/dw")
    parser.add_argument("--dc-only", action="store_true", help="compute only Omega=0 conductivity")
    parser.add_argument("--thermal-window", type=float, default=6.0, help="thermal integration window L_T (default: 6, window [-Omega-L_T*T, L_T*T])")
    parser.add_argument("--points-per-T", type=float, default=2.0, help="thermal mesh density (default: 2 gives dw=T/2)")
    parser.add_argument("--dc-min-points", type=int, default=50, help="minimum linear mesh points in the DC window [-L_T*T,L_T*T] before adding Sig.out points (default: 50)")
    parser.add_argument("--points-per-omega", type=float, default=40.0, help="large-Omega box resolution (default: 40)")
    parser.add_argument("--omega-mesh-crossover-T", type=float, default=2.0, help="use Omega/points_per_omega spacing only for Omega > this*T (default: 2)")
    parser.add_argument("--sigma-spacing-factor", type=float, default=1.0, help="resolve local Sig.out mesh; 1 means dw no larger than local Sig.out spacing")
    parser.add_argument("--dw-min", type=float, default=1e-5, help="lower bound on adaptive dw")
    parser.add_argument("--dw-max", type=float, default=None, help="upper bound on adaptive dw")
    parser.add_argument("--directions", default=None, help="comma list from xx, yy, xy, yx, trace (default: trace)")
    parser.add_argument("--tensor", action="store_true", help="shortcut for --directions xx,yy,xy,yx")
    parser.add_argument("--mode", choices=["woodbury", "direct"], default="woodbury")
    parser.add_argument("--svd-tol", type=float, default=1e-10, help="relative SVD cutoff; negative disables SVD")
    parser.add_argument("--spin", type=float, default=2.0)
    parser.add_argument("--units", choices=["e2h", "si"], default="si", help="output units: si=Siemens=1/Ohm (default), or e2h")
    parser.add_argument("--k-limit", type=int, default=None, help="debug: use only first k-limit points")
    parser.add_argument("--validate", action="store_true", help="compare direct and Woodbury Green functions")
    parser.add_argument("--output", default="tbg_optics.dat")
    args = parser.parse_args()

    workdir, ene, vx, vy, channels, sigmas, projectors = load_inputs(args)
    if args.sigma_scale != 1.0:
        sigmas = [(sw, None if ssig is None else args.sigma_scale * ssig) for sw, ssig in sigmas]
    beta = read_beta(workdir, channels, args.beta)
    mu = resolve_mu(workdir, channels, args.mu, args.mu_source, args.ef_check_tol)
    unit_prefactor, unit_label = optics_unit_prefactor(workdir, args.units)

    print(f"workdir={workdir}")
    print(f"nk={ene.shape[0]} nb={ene.shape[1]} nw=adaptive beta={beta} T={1.0/beta} mu={mu} mode={args.mode}")
    print(f"units={args.units} prefactor={unit_prefactor:.12g} [{unit_label}]")
    print("channels=" + ", ".join([f"{c.projector_key}:{c.sigma_path}" for c in channels]))

    if args.validate:
        validate_woodbury(ene, projectors, sigmas, mu, args.eta, args.svd_tol)

    if args.tensor and args.directions is not None:
        raise ValueError("Use either --tensor or --directions, not both")
    if args.tensor:
        directions = ["xx", "yy", "xy", "yx"]
    elif args.directions is None:
        directions = ["trace"]
    else:
        directions = [d.strip() for d in args.directions.split(",") if d.strip()]
    for d in directions:
        if d not in {"xx", "yy", "xy", "yx", "trace"}:
            raise ValueError(f"unknown direction {d!r}")

    sigma_dw = estimate_sigma_spacing(sigmas)
    default_step = min(sigma_dw, (1.0 / beta) / args.points_per_T)
    omega_step = args.omega_step if args.omega_step is not None else default_step
    omega_min = args.omega_min if args.omega_min is not None else omega_step
    omegas, accum = conductivity_adaptive(
        ene,
        vx,
        vy,
        projectors,
        sigmas,
        beta,
        mu,
        args.eta,
        args.omega_max,
        omega_min,
        omega_step,
        args.mode,
        args.svd_tol,
        directions,
        args.spin,
        unit_prefactor,
        args.k_limit,
        args.thermal_window,
        args.points_per_T,
        args.dc_min_points,
        args.points_per_omega,
        args.omega_mesh_crossover_T,
        args.sigma_spacing_factor,
        args.dw_min,
        args.dw_max,
        args.include_dc,
        args.dc_only,
    )
    out = args.output if os.path.isabs(args.output) else os.path.join(workdir, args.output)
    write_output(out, omegas, accum, unit_label)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
