#!/usr/bin/env python3
"""Directional density response of the implemented frozen-Sigma lattice map.

No impurity vertex or susceptibility matrix is constructed.  At each k and
kept Matsubara frequency we apply dG = G (dH+dSigma_emb-dmu) G.  The static
resolvent is subtracted before summation and restored with the exact Frechet
derivative of the Fermi matrix function.  Quadrature weights are imported
from dmft2 and applied ONCE, after the two resolvents have been multiplied.

The physical embedding is Pi_W phi X, not the untruncated overlap phi X.
Both dPi_W and dX are included.  The actual Fortran overlap-symmetry map is
exported on matrix units rather than inferred from a diagonal self-energy.
Fixed-index windows cutting a degeneracy are rejected, never regularized.

This initial diagnostic uses the saved LAPW basis and stored complement.
It is an action of the charge-map Jacobian, not a full DMFT susceptibility.
"""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import sys
import json
import time
import numpy as np
from scipy.linalg import eigh
from scipy.special import expit
from types import SimpleNamespace


def herm(a):
    return (a + a.swapaxes(-1, -2).conj()) / 2


@dataclass
class Metadata:
    mu: float
    temperature: float
    omega: np.ndarray
    weights: np.ndarray
    dimensions: np.ndarray
    cix: np.ndarray
    nindo: np.ndarray
    isx: np.ndarray
    active: np.ndarray
    normalization: np.ndarray
    sigma_inf: np.ndarray
    sigma: np.ndarray
    sym_real: np.ndarray
    sym_imag: np.ndarray

    @classmethod
    def read(cls, path):
        with open(path, "rb") as f:
            if f.read(8) != b"EDMRSP01":
                raise ValueError("bad metadata magic")
            version, m, c, m2, o, nw = np.fromfile(f, "<i4", 6)
            if version != 2:
                raise ValueError("unsupported version: re-export metadata")

            def a(dtype, shape):
                v = np.fromfile(f, dtype, int(np.prod(shape)))
                if v.size != np.prod(shape):
                    raise ValueError("truncated metadata")
                return v.reshape(shape, order="F")

            mu, temp = a("<f8", (2,))
            om = a("<f8", (nw,))
            w = a("<f8", (nw,))
            dims = a("<i4", (c,))
            cix = a("<i4", (o,)) - 1
            nindo = a("<i4", (o,))
            isx = a("<i4", (m2, o)) - 1
            sigind = a("<i4", (m, m, c))
            active = np.array([np.diag(sigind[:, :, ic]) > 0 for ic in range(c)]).T
            norm = a("<c16", (m, m, c))
            si = a("<c16", (m, m, c))
            sigma = a("<c16", (m, m, c, nw))
            size = m * m * c
            re = a("<c16", (size, size))
            im = a("<c16", (size, size))
            if f.read():
                raise ValueError("trailing metadata")
        return cls(
            mu, temp, om, w, dims, cix, nindo, isx, active, norm, si, sigma, re, im
        )

    def sym(self, s):
        v = s.ravel(order="F")
        return (self.sym_real @ v.real + self.sym_imag @ v.imag).reshape(
            s.shape, order="F"
        )

    def local_columns(self, raw):
        """Map band/orbital storage into separate impurity blocks."""
        m = self.normalization.shape[0]
        c = len(self.dimensions)
        out = np.zeros((*raw.shape[:2], m, c), complex)
        for orb, ic in enumerate(self.cix):
            if ic < 0:
                continue
            for j in range(self.nindo[orb]):
                out[:, :, self.isx[j, orb], ic] = raw[:, :, j, orb]
        return out

    def compact(self, local):
        return np.concatenate(
            [local[..., :n, c] for c, n in enumerate(self.dimensions)], axis=-1
        )

    def sigmas(self):
        n = int(sum(self.dimensions))
        s = np.zeros((len(self.omega), n, n), complex)
        si = np.zeros((n, n), complex)
        offset = 0
        for c, m in enumerate(self.dimensions):
            sl = slice(offset, offset + m)
            s[:, sl, sl] = self.sigma[:m, :m, c, :].transpose(2, 0, 1)
            si[sl, sl] = self.sigma_inf[:m, :m, c]
            offset += m
        return s, si


def invsqrt(s):
    e, v = eigh(herm(s))
    if e.min() < 1e-12:
        raise ValueError("singular orbital overlap")
    return (v / np.sqrt(e)) @ v.conj().T


def inverse_sqrt_derivative(x, ds):
    # X=S^-1/2.  Use X eigenvalues to anchor to the saved normalization.
    ex, v = eigh(herm(x))
    if ex.min() <= 0:
        raise ValueError("nonpositive normalization")
    r = 1 / ex
    return (
        v
        @ (
            -(v.conj().T @ ds @ v)
            / (r[:, None] * r[None, :] * (r[:, None] + r[None, :]))
        )
        @ v.conj().T
    )


class FrozenSigmaResponse:
    @classmethod
    def from_export(cls, prefix, **context_options):
        """Read a standalone EDMFT_RESPONSE_EXPORT from a central TOT run."""
        prefix = Path(prefix)
        m = Metadata.read(str(prefix) + ".meta")
        records = []
        kw = []
        for path in sorted(prefix.parent.glob(prefix.name + ".k????????.bin")):
            with path.open("rb") as fh:
                if fh.read(8) != b"EDMRPK01":
                    raise ValueError("bad k export magic")
                version, k, nbas, iso, n, first, last, wfirst, wlast, m2, o, nlo = (
                    np.fromfile(fh, "<i4", 12)
                )
                if version != 1 or iso != 1:
                    raise ValueError("initial response requires non-SOC export v1")

                def read(dtype, shape):
                    a = np.fromfile(fh, dtype, int(np.prod(shape)))
                    if a.size != np.prod(shape):
                        raise ValueError("truncated k export")
                    return a.reshape(shape, order="F")

                weight = float(read("<f8", (1,))[0])
                e = read("<f8", (n,))
                c = read("<c16", (nbas, n, iso))
                u = read("<c16", (n, m2, o))
                labels = read("<i4", (3, nbas))
                if fh.read():
                    raise ValueError("trailing k export bytes")
            records.append(
                SimpleNamespace(
                    kpoint=k,
                    nbasis=nbas,
                    nlo=nlo,
                    nstored=n,
                    ncomponents=iso,
                    nemin=first,
                    nemax=last,
                    dm_nemin=wfirst,
                    dm_nemax=wlast,
                    maxdim=m2,
                    norbitals=o,
                    energies=e,
                    coefficients=c,
                    projector=u,
                    labels=labels,
                )
            )
            kw.append(weight / 2)
        if not records:
            raise ValueError("no k exports")
        if [r.kpoint for r in records] != list(range(1, len(records) + 1)):
            raise ValueError("missing or out-of-order k-point exports")

        def signature(r):
            return (r.nstored, r.nemin, r.nemax, r.dm_nemin, r.dm_nemax)

        if any(signature(r) != signature(records[0]) for r in records):
            raise ValueError(
                "response requires a fixed band count and fixed-index window"
            )
        if abs(sum(kw) - 1) > 1e-10:
            raise ValueError(
                "k-point weights do not sum to one: incomplete or magnetic export"
            )
        normalized = np.array([r.projector for r in records])
        local = m.local_columns(normalized)
        raw_local = np.zeros_like(local)
        for c, n in enumerate(m.dimensions):
            raw_local[:, :, :n, c] = local[:, :, :n, c] @ np.linalg.inv(
                m.normalization[:n, :n, c]
            )
        raw = np.zeros_like(normalized)
        for orb, ic in enumerate(m.cix):
            if ic < 0:
                continue
            for j in range(m.nindo[orb]):
                raw[:, :, j, orb] = raw_local[:, :, m.isx[j, orb], ic]
        d = records[0]
        obj = cls(
            m,
            np.array([r.energies for r in records]),
            raw,
            normalized,
            np.array(kw),
            np.arange(d.dm_nemin - d.nemin, d.dm_nemax - d.nemin + 1),
            **context_options,
        )
        return obj, records

    def __init__(self, metadata, energies, raw, normalized, kw, window):
        self.meta = metadata
        self.eps = np.asarray(energies)
        self.kw = np.asarray(kw)
        self.W = np.asarray(window)
        self.Q = np.setdiff1d(np.arange(self.eps.shape[1]), self.W)
        self.raw = metadata.local_columns(raw)
        given = metadata.local_columns(normalized)
        # The exported normalization is authoritative.  Inferring it from
        # a projector pseudoinverse is unnecessary and can lose inactive
        # orbital channels when their global raw overlaps vanish.
        self.X = np.array(metadata.normalization, dtype=complex)
        rebuilt = self.apply_x(self.raw, self.X)
        if np.max(abs(rebuilt - given)) > 1e-10:
            raise ValueError("raw/normalized projector mismatch")
        self.u = np.zeros_like(rebuilt)
        self.u[:, self.W] = rebuilt[:, self.W]
        self.sigma, self.sinf = metadata.sigmas()
        self.uc = metadata.compact(self.u)
        self.s0 = self.overlap(self.raw[:, self.W])
        self.offset = np.zeros_like(self.s0)
        for c, m in enumerate(metadata.dimensions):
            xi = np.linalg.inv(self.X[:m, :m, c])
            self.offset[:m, :m, c] = xi @ xi - self.s0[:m, :m, c]
        self.minimum_gap = (
            float(np.min(abs(self.eps[:, self.W, None] - self.eps[:, None, self.Q])))
            if len(self.Q)
            else float("inf")
        )
        if self.minimum_gap < 1e-8:
            raise ValueError(f"window cuts degeneracy: {self.minimum_gap} Ry")

    def apply_x(self, q, x):
        return np.einsum("knic,ijc->knjc", q, x, optimize=True)

    def overlap(self, q):
        return self.meta.sym(
            np.einsum("k,knic,knjc->ijc", self.kw, q.conj(), q, optimize=True)
        )

    def embedding_response(self, dh):
        dq = np.zeros_like(self.raw)
        W, Q = self.W, self.Q
        for k, e in enumerate(self.eps):
            a = dh[k][np.ix_(Q, W)] / (e[W][None, :] - e[Q][:, None])
            dq[k, Q] = np.einsum("an,nic->aic", a, self.raw[k, W])
            dq[k, W] = np.einsum("na,aic->nic", a.conj().T, self.raw[k, Q])
        ds = self.meta.sym(
            np.einsum(
                "k,knic,knjc->ijc",
                self.kw,
                self.raw[:, W].conj(),
                dq[:, W],
                optimize=True,
            )
            + np.einsum(
                "k,knic,knjc->ijc",
                self.kw,
                dq[:, W].conj(),
                self.raw[:, W],
                optimize=True,
            )
        )
        dx = np.zeros_like(self.X)
        for c, m in enumerate(self.meta.dimensions):
            active = np.flatnonzero(self.meta.active[:m, c])
            ix = np.ix_(active, active)
            dx[:, :, c][ix] = inverse_sqrt_derivative(
                self.X[:, :, c][ix], ds[:, :, c][ix]
            )
        du = self.apply_x(dq, self.X)
        du[:, W] += self.apply_x(self.raw[:, W], dx)
        return self.meta.compact(du), ds, dx

    def embedding_endpoint(self, dh, t):
        q = np.zeros_like(self.raw)
        for k, e in enumerate(self.eps):
            _, v = eigh(np.diag(e) + t * dh[k])
            v = v[:, self.W]
            q[k] = np.einsum("ij,jac->iac", v @ v.conj().T, self.raw[k])
        s = self.overlap(q) + self.offset
        x = np.zeros_like(self.X)
        for c, m in enumerate(self.meta.dimensions):
            x[:m, :m, c] = np.eye(m)
            active = np.flatnonzero(self.meta.active[:m, c])
            ix = np.ix_(active, active)
            x[:, :, c][ix] = invsqrt(s[:, :, c][ix])
        return self.meta.compact(self.apply_x(q, x))

    @staticmethod
    def fermi_matrix_derivative(e, v, d, temp):
        f = expit(-e / temp)
        gap = e[:, None] - e[None, :]
        divided = np.empty_like(gap)
        near = abs(gap) < 1e-8
        np.divide(f[:, None] - f[None, :], gap, out=divided, where=~near)
        mid = expit(-(e[:, None] + e[None, :]) / (2 * temp))
        divided[near] = (-mid * (1 - mid) / temp)[near]
        return v @ ((v.conj().T @ d @ v) * divided) @ v.conj().T

    def density(self, dh=None, t=0.0, u=None):
        """Independent nonlinear endpoint density for directional tests."""
        if u is None:
            u = self.uc
        result = []
        z = 1j * self.meta.omega
        T = self.meta.temperature
        for k, eps in enumerate(self.eps):
            h = np.diag(eps).astype(complex)
            if dh is not None:
                h += t * dh[k]
            uk = u[k]
            hi = h + uk @ self.sinf @ uk.conj().T - self.meta.mu * np.eye(len(eps))
            e, v = eigh(hi)
            rho = (v * expit(-e / T)) @ v.conj().T
            for a in range(0, len(z), 32):
                zz = z[a : a + 32]
                sg = self.sigma[a : a + 32]
                g = np.linalg.inv(
                    zz[:, None, None] * np.eye(len(eps))
                    - (h - self.meta.mu * np.eye(len(eps)))[None]
                    - (uk @ sg) @ uk.conj().T
                )
                gi = (v[None] * (1 / (zz[:, None] - e))[:, None, :]) @ v.conj().T
                rho += (
                    2
                    * T
                    * np.einsum(
                        "w,wij->ij", self.meta.weights[a : a + 32], herm(g - gi)
                    )
                )
            result.append(2 * self.kw[k] * herm(rho))
        return np.array(result)

    def sandwich(self, gw, gq, d):
        """G D G with dense W block and diagonal complement; no Nstored^3 GEMM."""
        W, Q = self.W, self.Q
        n = self.eps.shape[1]
        d = np.broadcast_to(d, (len(gw), n, n))
        out = np.zeros_like(d)
        out[:, W[:, None], W] = gw @ d[:, W[:, None], W] @ gw
        out[:, W[:, None], Q] = (gw @ d[:, W[:, None], Q]) * gq[:, None, :]
        out[:, Q[:, None], W] = gq[:, :, None] * (d[:, Q[:, None], W] @ gw)
        out[:, Q[:, None], Q] = gq[:, :, None] * d[:, Q[:, None], Q] * gq[:, None, :]
        return out

    def response(self, dh, fixed_number=True, embedding=True):
        """Apply the analytic resolvent derivative, retaining signed modes."""
        dh = np.asarray(dh, dtype=complex)
        expected = (*self.eps.shape, self.eps.shape[1])
        if dh.shape != expected or not np.all(np.isfinite(dh)):
            raise ValueError(f"potential must be finite with shape {expected}")
        if np.max(abs(dh - dh.swapaxes(-1, -2).conj())) > 1e-10 * max(
            1.0, float(np.max(abs(dh)))
        ):
            raise ValueError("potential derivative is not Hermitian")
        dh = herm(dh)
        t0 = time.monotonic()
        W, Q = self.W, self.Q
        T = self.meta.temperature
        z = 1j * self.meta.omega
        du, ds, dx = self.embedding_response(dh)
        if not embedding:
            du[:] = 0
        answer = []
        mu_answer = []
        for k, eps in enumerate(self.eps):
            uk = self.uc[k]
            duk = du[k]
            n = len(eps)
            hi = np.diag(eps - self.meta.mu) + uk @ self.sinf @ uk.conj().T
            e, v = eigh(hi)
            ew, vw = eigh(hi[np.ix_(W, W)])
            di = dh[k] + duk @ self.sinf @ uk.conj().T + uk @ self.sinf @ duk.conj().T
            dr = self.fermi_matrix_derivative(e, v, di, T)
            dm = self.fermi_matrix_derivative(e, v, -np.eye(n), T)
            for a in range(0, len(z), 32):
                zz = z[a : a + 32]
                sg = self.sigma[a : a + 32]
                # Only W has a dense central resolvent.  Outside it G is
                # diagonal, so no enlarged Dyson inversion is necessary.
                gw = np.linalg.inv(
                    zz[:, None, None] * np.eye(len(W))
                    - hi[np.ix_(W, W)][None]
                    - (uk[W] @ (sg - self.sinf)) @ uk[W].conj().T
                )
                gq = 1 / (zz[:, None] - (eps[Q] - self.meta.mu))
                # The static resolvent has the same W/complement blocks.
                # Building its full Nstored x Nstored spectral product at
                # every frequency would add an unnecessary cubic cost.
                giw = (vw[None] * (1 / (zz[:, None] - ew))[:, None, :]) @ vw.conj().T
                d = dh[k][None] + (duk @ sg) @ uk.conj().T + (uk @ sg) @ duk.conj().T
                # BLAS matrix products are intentionally explicit here; do
                # not square a quadrature-weighted Green function from dmft2.
                dg = self.sandwich(gw, gq, d)
                dgi = self.sandwich(giw, gq, di)
                ww = self.meta.weights[a : a + 32]
                dr += 2 * T * np.einsum("w,wij->ij", ww, herm(dg - dgi))
                dm[np.ix_(W, W)] += (
                    2 * T * np.einsum("w,wij->ij", ww, herm(-gw @ gw + giw @ giw))
                )
            answer.append(2 * self.kw[k] * herm(dr))
            mu_answer.append(2 * self.kw[k] * herm(dm))
        answer = np.array(answer)
        mu_answer = np.array(mu_answer)
        dn = float(np.trace(answer, axis1=1, axis2=2).real.sum())
        compress = float(np.trace(mu_answer, axis1=1, axis2=2).real.sum())
        mu1 = 0.0
        if fixed_number:
            if abs(compress) < 1e-8:
                raise ValueError(
                    "ill-conditioned fixed-N response; select fixed mu explicitly"
                )
            mu1 = -dn / compress
            answer += mu1 * mu_answer
        info = {
            "seconds": time.monotonic() - t0,
            "dN_fixed_mu": dn,
            "dN_dmu": compress,
            "dmu": mu1,
            "dN": float(np.trace(answer, axis1=1, axis2=2).real.sum()),
            "minimum_window_gap_Ry": self.minimum_gap,
            "normalization_derivative_norm": float(np.linalg.norm(dx)),
        }
        return answer, info

    def write_density(self, prefix, rho, dumps):
        for k, (r, d) in enumerate(zip(rho, dumps), 1):
            p = Path(f"{prefix}.k{k:08d}.bin")
            p.parent.mkdir(parents=True, exist_ok=True)
            with p.open("wb") as f:
                f.write(b"EDMRHO01")
                np.array(
                    [1, k, d.nbasis, d.ncomponents, d.nstored, d.nemin, d.nemax], "<i4"
                ).tofile(f)
                np.asarray(d.energies, "<f8").tofile(f)
                np.asarray(d.coefficients, "<c16").ravel(order="F").tofile(f)
                np.asarray(r, "<c16").ravel(order="F").tofile(f)


def main(argv=None):
    """Standalone directional action; the screening driver calls the same API."""
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--export", required=True, type=Path, help="EDMFT_RESPONSE_EXPORT prefix"
    )
    parser.add_argument(
        "--potential",
        required=True,
        type=Path,
        help="Hermitian dV[k,i,j] .npy array in Ry per strain",
    )
    parser.add_argument(
        "--output", required=True, type=Path, help="signed response-density file prefix"
    )
    parser.add_argument(
        "--fixed-mu",
        action="store_true",
        help="hold mu fixed; default constrains total particle number",
    )
    parser.add_argument(
        "--backend",
        choices=["reference", "cached"],
        default="reference",
        help="cached reuses central Matsubara contractions; same stored-band approximation",
    )
    parser.add_argument("--cache-mb", type=float, default=256.0)
    args = parser.parse_args(argv)
    if args.backend == "cached":
        from frozen_sigma_tangent import CachedFrozenSigmaResponse

        model, records = CachedFrozenSigmaResponse.from_export(
            args.export, cache_mb=args.cache_mb
        )
    else:
        model, records = FrozenSigmaResponse.from_export(args.export)
    perturbation = np.load(args.potential, allow_pickle=False)
    rho, report = model.response(perturbation, fixed_number=not args.fixed_mu)
    model.write_density(args.output, rho, records)
    report.update(
        approximation="fixed central LAPW basis and stored complement",
        potential_units="Ry per unit strain",
        ensemble="fixed mu" if args.fixed_mu else "fixed N",
        backend=args.backend,
    )
    Path(str(args.output) + ".json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        raise SystemExit(2)
