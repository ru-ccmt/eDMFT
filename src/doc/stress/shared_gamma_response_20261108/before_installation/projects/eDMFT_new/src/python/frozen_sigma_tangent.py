#!/usr/bin/env python3
"""Cached forward/adjoint of the frozen-Sigma, stored-band density map.

The nonlinear endpoint implementation in frozen_sigma_response is deliberately
kept independent as an oracle. This module differentiates that SAME stored-band
map. response() works in an orthonormal representation; generalized_response()
also accepts explicit changes of the overlap metric and raw orbital overlaps.
All inputs are in the central KS gauge and the complement contains only
exported bands. This module does not construct radial, core, geometry, or
spatial-potential perturbations; native adapters must supply those separately.

All adjoints below use Re sum(conj(y)*x), with no implicit k, spin, orbital,
or Matsubara weights. The output density already contains the 2*w_k factor.
The forward and reverse operations both include window motion, symmetrized
global overlap normalization, and (optionally) the fixed-N mu constraint.
No numerical Jacobian, numerical transpose, or impurity vertex is used.
"""
from __future__ import annotations

import argparse
import copy
import json
import time
from pathlib import Path

import numpy as np
from scipy.linalg import eigh
from scipy.special import expit

from frozen_sigma_response import FrozenSigmaResponse, herm


def real_inner(a, b):
    """Real Frobenius pairing, including every array index exactly once."""
    return float(np.vdot(a, b).real)


class CachedFrozenSigmaResponse(FrozenSigmaResponse):
    """Reusable central context; callers must treat its central data as fixed.

    WINDOW response contractions and/or resolvents are cached, never full-band
    Green functions. cache_mb bounds their combined storage; overflow chunks
    are recomputed. Static
    eigensystems, divided differences, normalization derivatives, and the
    chemical-potential response are computed once. Construct a new context
    when mu, Sigma, projectors, bands, weights, or temperature change.
    """

    def __init__(
        self, *args, cache_mb=256.0, frequency_chunk=32, precontract=True, **kwargs
    ):
        # Isolate the context from caller-owned arrays before caching. Updating
        # an input self-energy or chemical potential requires a NEW context.
        super().__init__(*copy.deepcopy(args), **copy.deepcopy(kwargs))
        if not np.isfinite(cache_mb) or cache_mb < 0 or frequency_chunk < 1:
            raise ValueError("invalid response cache limits")
        self.cache_limit = int(cache_mb * 1024**2)
        self.frequency_chunk = int(frequency_chunk)
        self._greens = {}
        self._kernels = {}
        self._kernel_bytes = 0
        self.precontract = bool(precontract)
        self._cache_bytes = 0
        self._ready = False
        self._central_scalars = (self.meta.mu, self.meta.temperature)
        self._inverse_count = 0
        self._prepare_seconds = 0.0
        self._norm = []
        for c, n in enumerate(self.meta.dimensions):
            active = np.flatnonzero(self.meta.active[:n, c])
            ix = np.ix_(active, active)
            if not len(active):
                self._norm.append((ix, None, None))
                continue
            ex, v = eigh(herm(self.X[:, :, c][ix]))
            if ex.min() <= 0:
                raise ValueError("nonpositive normalization")
            r = 1.0 / ex
            divided = -1.0 / (r[:, None] * r[None, :] * (r[:, None] + r[None, :]))
            self._norm.append((ix, v, divided))
        self._denominators = [e[self.W][None, :] - e[self.Q][:, None] for e in self.eps]
        for value in list(vars(self).values()) + list(vars(self.meta).values()):
            if isinstance(value, np.ndarray):
                value.flags.writeable = False

    def normalization_action(self, ds):
        """Frechet derivative of S^-1/2, self-adjoint on Hermitian blocks."""
        result = np.zeros_like(self.X)
        for c, (ix, v, divided) in enumerate(self._norm):
            if v is not None:
                result[:, :, c][ix] = (
                    v
                    @ ((v.conj().T @ herm(ds[:, :, c][ix]) @ v) * divided)
                    @ v.conj().T
                )
        return result

    def symmetry_adjoint(self, y):
        """Transpose the actual exported REAL-linear Fortran symmetry map.

        In particular, the map can conjugate its input. Multiplication by a
        complex conjugate transpose alone would therefore be incorrect.
        """
        v = y.ravel(order="F")
        result = (self.meta.sym_real.conj().T @ v).real + 1j * (
            self.meta.sym_imag.conj().T @ v
        ).real
        return result.reshape(y.shape, order="F")

    def _raw_direction(self, derivative):
        """Unnormalized orbital overlaps, in the central coefficient gauge."""
        value = np.asarray(derivative, dtype=complex)
        if value.shape != self.raw.shape or not np.all(np.isfinite(value)):
            raise ValueError("invalid raw-projector derivative")
        return value

    def embedding_response(self, dh, raw_derivative=None):
        """Differentiate Pi_W(H) phi X, including an explicit dphi if supplied.

        dPi_W phi and Pi_W dphi are separate terms. Only the W rows of dphi
        survive the central Pi_W; adding its Q rows would double count window
        motion. The same symmetrized k-sum and dX act on their sum.
        """
        W, Q = self.W, self.Q
        dq = np.zeros_like(self.raw)
        for k, denominator in enumerate(self._denominators):
            a = dh[k][np.ix_(Q, W)] / denominator
            dq[k, Q] = np.einsum("an,nic->aic", a, self.raw[k, W])
            dq[k, W] = np.einsum("na,aic->nic", a.conj().T, self.raw[k, Q])
        if raw_derivative is not None:
            dq[:, W] += self._raw_direction(raw_derivative)[:, W]
        ds = self.meta.sym(
            np.einsum("k,knic,knjc->ijc", self.kw, self.raw[:, W].conj(), dq[:, W])
            + np.einsum("k,knic,knjc->ijc", self.kw, dq[:, W].conj(), self.raw[:, W])
        )
        dx = self.normalization_action(ds)
        du = self.apply_x(dq, self.X)
        du[:, W] += self.apply_x(self.raw[:, W], dx)
        return self.meta.compact(du), ds, dx

    def embedding_pullback(self, y):
        """Return covectors of (dH, dphi) for embedding_response.

        Reverse compacting, U=P X, dX, symmetry, the weighted k sum, and
        finally dPi_W. The 1/2 in the last Hermitian projection is necessary:
        both QW and WQ entries participate in the Frobenius pairing.
        """
        if y.shape != self.uc.shape or not np.all(np.isfinite(y)):
            raise ValueError("invalid projector covector")
        W, Q = self.W, self.Q
        local = np.zeros_like(self.raw)
        start = 0
        for c, n in enumerate(self.meta.dimensions):
            local[:, :, :n, c] = y[:, :, start : start + n]
            start += n
        dq_bar = self.apply_x(local, self.X.swapaxes(0, 1).conj())
        dx_bar = np.einsum("knic,knjc->ijc", self.raw[:, W].conj(), local[:, W])
        ds_bar = self.symmetry_adjoint(self.normalization_action(dx_bar))
        ds_bar = ds_bar + ds_bar.swapaxes(0, 1).conj()
        dq_bar[:, W] += self.kw[:, None, None, None] * self.apply_x(
            self.raw[:, W], ds_bar
        )
        result = np.zeros((*self.eps.shape, self.eps.shape[1]), complex)
        for k, denominator in enumerate(self._denominators):
            a_bar = np.einsum("aic,nic->an", dq_bar[k, Q], self.raw[k, W].conj())
            a_bar += np.einsum("aic,nic->an", self.raw[k, Q], dq_bar[k, W].conj())
            result[k][np.ix_(Q, W)] = a_bar / denominator
        raw_bar = np.zeros_like(self.raw)
        raw_bar[:, W] = dq_bar[:, W]
        return herm(result), raw_bar

    def embedding_adjoint(self, y):
        """Compatibility interface: potential covector with dphi held zero."""
        return self.embedding_pullback(y)[0]

    @staticmethod
    def _divided_difference(e, temperature):
        f = expit(-e / temperature)
        gap = e[:, None] - e[None, :]
        near = abs(gap) < 1e-8
        divided = np.empty_like(gap)
        np.divide(f[:, None] - f[None, :], gap, out=divided, where=~near)
        mid = expit(-(e[:, None] + e[None, :]) / (2 * temperature))
        divided[near] = (-mid * (1 - mid) / temperature)[near]
        return divided

    def prepare(self):
        """Build the reusable central context and fixed-mu charge derivative."""
        if (self.meta.mu, self.meta.temperature) != self._central_scalars:
            raise ValueError(
                "central mu/temperature changed; construct a new response context"
            )
        if self._ready:
            return
        started = time.monotonic()
        self._static = []
        W, Q = self.W, self.Q
        for k, eps in enumerate(self.eps):
            u = self.uc[k]
            hi = np.diag(eps - self.meta.mu) + u @ self.sinf @ u.conj().T
            e, v = eigh(hi)
            ew, vw = eigh(hi[np.ix_(W, W)])
            self._static.append(
                (
                    hi[np.ix_(W, W)],
                    e,
                    v,
                    ew,
                    vw,
                    self._divided_difference(e, self.meta.temperature),
                )
            )
        mu = []
        density = []
        projector_gradient = np.zeros_like(self.uc)
        for k in range(len(self.eps)):
            dm = self._fermi(k, -np.eye(self.eps.shape[1]))
            _, e, v, _, _, _ = self._static[k]
            u = self.uc[k]
            rho_inf = (v * expit(-e / self.meta.temperature)) @ v.conj().T
            rho = rho_inf.copy()
            projector_local = 2 * rho_inf @ u @ self.sinf
            for a, stop, gw, gq, giw in self._chunks(k):
                weights = 2 * self.meta.temperature * self.meta.weights[a:stop]
                rho[np.ix_(W, W)] += np.einsum("w,wij->ij", weights, herm(gw - giw))
                dm[np.ix_(W, W)] += np.einsum(
                    "w,wij->ij", weights, herm(-gw @ gw + giw @ giw)
                )
                sigma = self.sigma[a:stop]
                part = (
                    (gw @ u[W]) @ sigma
                    + (gw.swapaxes(-1, -2).conj() @ u[W])
                    @ sigma.swapaxes(-1, -2).conj()
                    - ((giw + giw.swapaxes(-1, -2).conj()) @ u[W]) @ self.sinf
                )
                projector_local[W] += np.einsum("w,wia->ia", weights, part)
            mu.append(2 * self.kw[k] * herm(dm))
            density.append(2 * self.kw[k] * herm(rho))
            projector_gradient[k] = 2 * self.kw[k] * projector_local
            if self.precontract:
                self._prepare_contraction(k)
        self.mu_response = np.asarray(mu)
        self.mu_response.flags.writeable = False
        self.compressibility = float(
            np.trace(self.mu_response, axis1=1, axis2=2).real.sum()
        )
        self._central_density = np.asarray(density)
        self._central_density.flags.writeable = False
        self._projector_potential_gradient, self._projector_raw_gradient = (
            self.embedding_pullback(projector_gradient)
        )
        self._ready = True
        self._prepare_seconds = time.monotonic() - started

    def _fermi(self, k, d):
        _, _, v, _, _, divided = self._static[k]
        return v @ ((v.conj().T @ d @ v) * divided) @ v.conj().T

    def _chunks(self, k):
        hi, e, v, ew, vw, divided = self._static[k]
        u = self.uc[k, self.W]
        for a in range(0, len(self.meta.omega), self.frequency_chunk):
            stop = min(a + self.frequency_chunk, len(self.meta.omega))
            z = 1j * self.meta.omega[a:stop]
            key = (k, a)
            gw = self._greens.get(key)
            if gw is None:
                gw = np.linalg.inv(
                    z[:, None, None] * np.eye(len(self.W))
                    - hi[None]
                    - (u @ (self.sigma[a:stop] - self.sinf)) @ u.conj().T
                )
                self._inverse_count += len(z)
                if (
                    self._cache_bytes + self._kernel_bytes + gw.nbytes
                    <= self.cache_limit
                ):
                    gw.flags.writeable = False
                    self._greens[key] = gw
                    self._cache_bytes += gw.nbytes
            gq = 1 / (z[:, None] - (self.eps[k, self.Q] - self.meta.mu))
            giw = (vw[None] * (1 / (z[:, None] - ew))[:, None, :]) @ vw.conj().T
            yield a, stop, gw, gq, giw

    def cache_info(self):
        return dict(
            resolvent_bytes=self._cache_bytes,
            resolvent_limit_bytes=self.cache_limit,
            contraction_bytes=self._kernel_bytes,
            total_cache_bytes=self._cache_bytes + self._kernel_bytes,
            precontracted_kpoints=len(self._kernels),
            dense_window_inversions=self._inverse_count,
            preparation_seconds=self._prepare_seconds,
            approximation="fixed central LAPW basis; stored complement",
        )

    @staticmethod
    def _pack(matrices):
        n = matrices.shape[-1]
        i, j = np.triu_indices(n, 1)
        return np.concatenate(
            (
                np.diagonal(matrices, axis1=-2, axis2=-1).real,
                np.sqrt(2) * matrices[..., i, j].real,
                np.sqrt(2) * matrices[..., i, j].imag,
            ),
            axis=-1,
        )

    @staticmethod
    def _unpack(vector, n):
        i, j = np.triu_indices(n, 1)
        nu = len(i)
        out = np.zeros((n, n), complex)
        out[np.arange(n), np.arange(n)] = vector[:n]
        out[i, j] = (vector[n : n + nu] + 1j * vector[n + nu :]) / np.sqrt(2)
        out[j, i] = out[i, j].conj()
        return out

    def _prepare_contraction(self, k):
        """Integrate the frequency dependence ONCE, if the memory cap allows.

        Only the WW block is a fourth-rank tensor. WQ blocks are independent
        for each outside band. QQ dynamic/static terms cancel identically.
        Store real independent Hermitian coordinates for WW, not a redundant
        complex susceptibility of the entire LAPW space. This is a local
        frozen-Sigma bubble/embedding derivative, not an impurity vertex or a
        dense spatial screening Jacobian.
        """
        nw, nq, m = len(self.W), len(self.Q), self.uc.shape[-1]
        needed = 8 * (nw**4 + 2 * nw**3 * m) + 16 * nq * (nw * nw + nw * m)
        reclaim = sum(
            value.nbytes for (kk, a), value in self._greens.items() if kk == k
        )
        if self._cache_bytes + self._kernel_bytes - reclaim + needed > self.cache_limit:
            return
        th = np.zeros((nw, nw, nw, nw), complex)
        tu = np.zeros((nw, nw, nw, m), complex)
        tv = np.zeros_like(tu)
        cross_h = np.zeros((nq, nw, nw), complex)
        cross_u = np.zeros((nq, nw, m), complex)
        u = self.uc[k, self.W]
        for a, stop, gw, gq, giw in self._chunks(k):
            weights = 2 * self.meta.temperature * self.meta.weights[a:stop]
            sigma = self.sigma[a:stop]
            th += np.einsum("w,wia,wbj->ijab", weights, gw, gw, optimize=True)
            th -= np.einsum("w,wia,wbj->ijab", weights, giw, giw, optimize=True)
            right = (sigma @ u.conj().T) @ gw
            left = (gw @ u) @ sigma
            right_i = (self.sinf @ u.conj().T) @ giw
            left_i = (giw @ u) @ self.sinf
            tu += np.einsum("w,wia,wbj->ijab", weights, gw, right, optimize=True)
            tu -= np.einsum("w,wia,wbj->ijab", weights, giw, right_i, optimize=True)
            tv += np.einsum("w,wib,waj->ijab", weights, left, gw, optimize=True)
            tv -= np.einsum("w,wib,waj->ijab", weights, left_i, giw, optimize=True)
            cross_h += np.einsum("w,wq,wij->qij", weights, gq, gw - giw, optimize=True)
            conjugate_part = (gw.swapaxes(-1, -2).conj() @ u) @ sigma.swapaxes(
                -1, -2
            ).conj()
            conjugate_part -= (giw.swapaxes(-1, -2).conj() @ u) @ self.sinf
            cross_u += 0.5 * (
                np.einsum("w,wq,wia->qia", weights, gq, left - left_i, optimize=True)
                + np.einsum(
                    "w,wq,wia->qia", weights, gq.conj(), conjugate_part, optimize=True
                )
            )
        # Transform the input and output WW indices to independent real
        # Hermitian coordinates without a dense change-of-basis multiply.
        i, j = np.triu_indices(nw, 1)
        cols_h = np.concatenate(
            (
                np.diagonal(th, axis1=2, axis2=3),
                (th[:, :, i, j] + th[:, :, j, i]) / np.sqrt(2),
                1j * (th[:, :, i, j] - th[:, :, j, i]) / np.sqrt(2),
            ),
            axis=-1,
        )
        kh = self._pack(herm(cols_h.transpose(2, 0, 1))).T.copy()
        cols_u = np.concatenate(
            ((tu + tv).reshape(nw, nw, -1), (1j * (tu - tv)).reshape(nw, nw, -1)),
            axis=-1,
        )
        ku = self._pack(herm(cols_u.transpose(2, 0, 1))).T.copy()
        cross_h = herm(cross_h)
        kernels = (kh, ku, cross_h, cross_u)
        size = sum(value.nbytes for value in kernels)
        if size != needed:
            raise RuntimeError("internal contraction cache accounting error")
        # Resolvents for this k are no longer used by response/adjoint. They
        # can still be recomputed on demand for independent scalar checks.
        for key in list(self._greens):
            if key[0] == k:
                self._cache_bytes -= self._greens.pop(key).nbytes
        for value in kernels:
            value.flags.writeable = False
        self._kernels[k] = kernels
        self._kernel_bytes += size

    def _contract_forward(self, k, dh, du):
        kh, ku, cross_h, cross_u = self._kernels[k]
        W, Q = self.W, self.Q
        result = np.zeros_like(dh)
        v = du[W].ravel()
        result[np.ix_(W, W)] = self._unpack(
            kh @ self._pack(dh[np.ix_(W, W)]) + ku @ np.r_[v.real, v.imag], len(W)
        )
        result[np.ix_(W, Q)] = np.einsum("qij,jq->iq", cross_h, dh[np.ix_(W, Q)])
        result[np.ix_(W, Q)] += np.einsum("qia,qa->iq", cross_u, du[Q].conj())
        result[np.ix_(Q, W)] = result[np.ix_(W, Q)].conj().T
        return result

    def _contract_adjoint(self, k, y):
        kh, ku, cross_h, cross_u = self._kernels[k]
        W, Q = self.W, self.Q
        hbar, ubar = np.zeros_like(y), np.zeros_like(self.uc[k])
        yy = self._pack(y[np.ix_(W, W)])
        hbar[np.ix_(W, W)] = self._unpack(kh.T @ yy, len(W))
        uv = ku.T @ yy
        nm = len(W) * self.uc.shape[-1]
        ubar[W] = (uv[:nm] + 1j * uv[nm:]).reshape(len(W), -1)
        hbar[np.ix_(W, Q)] = np.einsum("qji,jq->iq", cross_h.conj(), y[np.ix_(W, Q)])
        hbar[np.ix_(Q, W)] = hbar[np.ix_(W, Q)].conj().T
        ubar[Q] = 2 * np.einsum("iq,qia->qa", y[np.ix_(W, Q)].conj(), cross_u)
        return hbar, ubar

    def _matrix(self, a, name):
        a = np.asarray(a, complex)
        shape = (*self.eps.shape, self.eps.shape[1])
        if a.shape != shape or not np.all(np.isfinite(a)):
            raise ValueError(f"{name} must be finite with shape {shape}")
        if np.max(abs(a - a.swapaxes(-1, -2).conj())) > 1e-10 * max(
            1.0, float(np.max(abs(a)))
        ):
            raise ValueError(f"{name} is not Hermitian")
        return herm(a)

    def _check_fixed_number(self):
        if abs(self.compressibility) < 1e-8:
            raise ValueError(
                "ill-conditioned fixed-N response; select fixed mu explicitly"
            )

    def response(self, dh, fixed_number=True, embedding=True, raw_derivative=None):
        """JVP (dH, dphi) -> weighted d(rho), in an orthonormal representation.

        raw_derivative is the EXPLICIT unnormalized orbital-overlap change;
        window motion and normalization are applied here, not by the caller.
        Use generalized_response for a changing overlap metric.
        """
        started = time.monotonic()
        dh = self._matrix(dh, "potential")
        self.prepare()
        du, ds, dx = self.embedding_response(dh, raw_derivative)
        if not embedding:
            du[:] = 0
        result = []
        for k, u in enumerate(self.uc):
            v = du[k]
            di = dh[k] + v @ self.sinf @ u.conj().T + u @ self.sinf @ v.conj().T
            dr = self._fermi(k, di)
            if k in self._kernels:
                dr += self._contract_forward(k, dh[k], v)
                result.append(2 * self.kw[k] * herm(dr))
                continue
            for a, stop, gw, gq, giw in self._chunks(k):
                sigma = self.sigma[a:stop]
                d = dh[k][None] + (v @ sigma) @ u.conj().T + (u @ sigma) @ v.conj().T
                dg = self.sandwich(gw, gq, d) - self.sandwich(giw, gq, di)
                dr += np.einsum(
                    "w,wij->ij",
                    2 * self.meta.temperature * self.meta.weights[a:stop],
                    herm(dg),
                )
            result.append(2 * self.kw[k] * herm(dr))
        result = np.asarray(result)
        dn = float(np.trace(result, axis1=1, axis2=2).real.sum())
        mu1 = 0.0
        if fixed_number:
            self._check_fixed_number()
            mu1 = -dn / self.compressibility
            result += mu1 * self.mu_response
        return result, dict(
            seconds=time.monotonic() - started,
            dN_fixed_mu=dn,
            dN_dmu=self.compressibility,
            dmu=mu1,
            dN=float(np.trace(result, axis1=1, axis2=2).real.sum()),
            minimum_window_gap_Ry=self.minimum_gap,
            normalization_derivative_norm=float(np.linalg.norm(dx)),
            cache=self.cache_info(),
        )

    def adjoint(self, y, fixed_number=True, embedding=True, return_raw=False):
        """VJP of response, including the nonlocal projector-normalization chain.

        With return_raw=True, return (H_bar,phi_bar); otherwise retain the
        original potential-only interface. embedding=False freezes U entirely.

        For fixed N, response = (I - b a^T/(a^T b)) response_mu, where b is
        d(rho)/dmu and a^T extracts the total trace. Its transpose is NOT the
        same projector. Apply it first, before the fixed-mu electronic VJP.
        """
        y = self._matrix(y, "density covector")
        self.prepare()
        if fixed_number:
            self._check_fixed_number()
            y -= (real_inner(self.mu_response, y) / self.compressibility) * np.eye(
                self.eps.shape[1]
            )
        answer = np.zeros_like(y)
        du_bar = np.zeros_like(self.uc)
        for k, u in enumerate(self.uc):
            yy = 2 * self.kw[k] * y[k]
            di_bar = self._fermi(k, yy)
            if k in self._kernels:
                hbar, ubar = self._contract_adjoint(k, yy)
                answer[k] = di_bar + hbar
                du_bar[k] = 2 * di_bar @ u @ self.sinf + ubar
                continue
            for a, stop, gw, gq, giw in self._chunks(k):
                weights = 2 * self.meta.temperature * self.meta.weights[a:stop]
                # Positive-frequency Hermitian projection is self-adjoint.
                # The resolvent sandwich itself is not: use both daggers.
                dg_bar = self.sandwich(gw.swapaxes(-1, -2).conj(), gq.conj(), yy)
                dgi_bar = self.sandwich(giw.swapaxes(-1, -2).conj(), gq.conj(), yy)
                dg_bar *= weights[:, None, None]
                di_bar -= np.einsum("w,wij->ij", weights, dgi_bar)
                answer[k] += dg_bar.sum(axis=0)
                sigma = self.sigma[a:stop]
                du_bar[k] += (
                    (dg_bar @ u) @ sigma.swapaxes(-1, -2).conj()
                    + (dg_bar.swapaxes(-1, -2).conj() @ u) @ sigma
                ).sum(axis=0)
            di_bar = herm(di_bar)
            answer[k] += di_bar
            du_bar[k] += 2 * di_bar @ u @ self.sinf
        answer = herm(answer)
        raw_bar = np.zeros_like(self.raw)
        if embedding:
            h_bar, raw_bar = self.embedding_pullback(du_bar)
            answer += h_bar
        if return_raw:
            return herm(answer), raw_bar
        return herm(answer)

    def central_density(self):
        """Same weighted central density as the independent endpoint oracle."""
        self.prepare()
        return self._central_density.copy()

    def generalized_response(self, dh, do, raw_derivative=None, fixed_number=True):
        """Moving-basis JVP (dH,dO,dphi) -> contravariant density coefficients.

        All inputs are in the SAME central KS gauge, O0=I, H0=diag(eps).
        phi is the raw overlap <basis|local orbital>. At first order the
        symmetric orthonormalization O^-1/2 gives
          h = dH - {dO,H0}/2,  phi1 = dphi - dO phi/2.
        Apply the existing normalized-window response to (h,phi1), then
        undo the moving orthonormalization of the density coefficients:
          rho1_coeff = rho1_orth - {dO,rho0}/2.
        Its trace alone is NOT dN: dN=Tr(rho1_coeff + rho0 dO). A density
        reconstruction adapter must ALSO differentiate the basis functions.
        This API supplies no geometry, radial, energy-search or core source.
        """
        dh = self._matrix(dh, "Hamiltonian derivative")
        do = self._matrix(do, "overlap derivative")
        raw = (
            np.zeros_like(self.raw)
            if raw_derivative is None
            else self._raw_direction(raw_derivative)
        )
        effective_h = dh - 0.5 * (do * self.eps[:, None, :] + self.eps[:, :, None] * do)
        effective_raw = raw - 0.5 * np.einsum("knm,kmic->knic", do, self.raw)
        result, info = self.response(
            effective_h, fixed_number=fixed_number, raw_derivative=effective_raw
        )
        result -= 0.5 * (do @ self._central_density + self._central_density @ do)
        info = dict(
            info,
            coefficient_trace=float(np.trace(result, axis1=1, axis2=2).real.sum()),
            metric_charge=real_inner(self._central_density, do),
        )
        info["dN"] = info["coefficient_trace"] + info["metric_charge"]
        return herm(result), info

    def generalized_adjoint(self, y, fixed_number=True):
        """Exact real-Frobenius transpose of generalized_response.

        Return (H_bar,O_bar,phi_bar). The density metric term uses the original
        y, before the fixed-N projection; projecting it again breaks covariance.
        """
        y = self._matrix(y, "density covector")
        hbar, rawbar = self.adjoint(y, fixed_number=fixed_number, return_raw=True)
        obar = -0.5 * (hbar * self.eps[:, None, :] + self.eps[:, :, None] * hbar)
        obar -= 0.5 * herm(np.einsum("knic,kmic->knm", rawbar, self.raw.conj()))
        obar -= 0.5 * (y @ self._central_density + self._central_density @ y)
        return hbar, herm(obar), rawbar

    def projector_gradients(self):
        """(H_bar,phi_bar) of the extra normalized-projector lattice term.

        Contains no direct Tr(rho dH). Input phi_bar is before window selection
        and normalization, suitable for reversing the native radial adapter.
        The scalar is the regularized lattice scalar, not the total Klein
        functional; no new total-stress term is installed by this method.
        """
        self.prepare()
        return (
            self._projector_potential_gradient.copy(),
            self._projector_raw_gradient.copy(),
        )

    def projector_potential_gradient(self):
        """Common band-potential gradient of the normalized-projector term.

        This is the gradient of the regularized lattice scalar with respect
        to U, followed by embedding_adjoint. It contains ONLY the extra
        dependence H -> Pi_W(H) phi X(H), not the direct Tr(rho dH) term.
        Thus it must not be added to that direct term a second time.

        If a spatial potential operator is V_band=B v and v=K q+v_exp,
        the common charge covector is c=K^T B^T R, R=this result. The same R
        is used for all six strain directions. Radial-projector and moving-
        basis derivatives are separate chains not implemented by this map.
        """
        self.prepare()
        return self._projector_potential_gradient.copy()

    def lattice_scalar(self, dh=None, t=0.0, u=None):
        """Independent nonlinear scalar oracle, NOT the total Klein functional.

        L = 2 sum_k w_k [Omega(H_inf) + 2T sum_w w_w Re(log det
            (iw-H_inf) - log det(iw-H-Sigma_emb))]. Energies are relative
        to the fixed central mu. Its H derivative is density(); its extra
        projector derivative is projector_potential_gradient(). The exact
        same exported quadrature and static subtraction define all three.
        This oracle deliberately uses full matrices, independently of the
        block algorithm, and is not used in a production response action.
        """
        if u is None:
            u = self.uc
        total = 0.0
        T = self.meta.temperature
        for k, eps in enumerate(self.eps):
            h = np.diag(eps - self.meta.mu).astype(complex)
            if dh is not None:
                h += t * dh[k]
            hi = h + u[k] @ self.sinf @ u[k].conj().T
            e = np.linalg.eigvalsh(hi)
            value = -T * np.logaddexp(0, -e / T).sum()
            for a in range(0, len(self.meta.omega), self.frequency_chunk):
                stop = min(a + self.frequency_chunk, len(self.meta.omega))
                z = 1j * self.meta.omega[a:stop]
                static = z[:, None, None] * np.eye(len(eps)) - hi
                dynamic = (
                    z[:, None, None] * np.eye(len(eps))
                    - h
                    - (u[k] @ self.sigma[a:stop]) @ u[k].conj().T
                )
                value += (
                    2
                    * T
                    * np.dot(
                        self.meta.weights[a:stop],
                        np.linalg.slogdet(static)[1] - np.linalg.slogdet(dynamic)[1],
                    )
                )
            total += 2 * self.kw[k] * value
        return float(total)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--export", required=True, type=Path)
    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument(
        "--input",
        type=Path,
        help="Hermitian band matrices (.npy), optionally a leading direction axis",
    )
    inputs.add_argument(
        "--projector-gradient",
        action="store_true",
        help="export the common projector-potential covector",
    )
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--adjoint", action="store_true")
    parser.add_argument("--fixed-mu", action="store_true")
    parser.add_argument(
        "--cache-mb",
        type=float,
        default=256.0,
        help="combined resolvent/contraction cache limit",
    )
    parser.add_argument(
        "--no-precontract",
        action="store_true",
        help="keep frequency-by-frequency action as a control",
    )
    args = parser.parse_args(argv)
    model, _ = CachedFrozenSigmaResponse.from_export(
        args.export, cache_mb=args.cache_mb, precontract=not args.no_precontract
    )
    if args.projector_gradient:
        if args.adjoint:
            parser.error("--adjoint requires --input")
        outgoing = model.projector_potential_gradient()
    else:
        incoming = np.load(args.input, allow_pickle=False)
        if incoming.ndim not in (3, 4) or incoming.shape[0] == 0:
            raise ValueError("input must be one band-matrix array or a nonempty batch")
        batch = incoming.ndim == 4
        directions = incoming if batch else incoming[None]
        results = []
        for direction in directions:
            if args.adjoint:
                result = model.adjoint(direction, fixed_number=not args.fixed_mu)
            else:
                result, _ = model.response(direction, fixed_number=not args.fixed_mu)
            results.append(result)
        outgoing = np.asarray(results) if batch else results[0]
    np.save(args.output, outgoing)
    print(json.dumps(model.cache_info(), indent=2))


if __name__ == "__main__":
    main()
