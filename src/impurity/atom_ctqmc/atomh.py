"""
Python port of atomh.cc focused on direct-basis Bold-CTQMC input generation.

This implementation keeps the direct-product occupation basis, builds atomic
Hamiltonian blocks compatible with the supplied hybridization symmetry,
and writes ``bctqmc.cix``. The original |j,mj> occupation basis is retained.

Spin-orbit is taken to be already included in ``Eimp``. To avoid any
double-counting, this generator intentionally does not expose a separate
``cx`` input.
"""

from __future__ import annotations

import argparse
import ast
import math
from collections import defaultdict, deque
from functools import lru_cache
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple, Union

import numpy as np

try:
    from scipy.optimize import linear_sum_assignment
except ImportError:  # pragma: no cover - optional dependency
    linear_sum_assignment = None


def _is_int_like(x: float, tol: float = 1e-12) -> bool:
    return abs(x - round(x)) < tol


@lru_cache(maxsize=None)
def _lgamma_cached(x: float) -> float:
    return math.lgamma(x)


def factorial(x: Union[int, float]) -> float:
    if isinstance(x, int):
        if x < 0:
            raise ValueError("factorial is defined only for non-negative values")
        return math.factorial(x)

    xx = float(x)
    if xx < 0:
        raise ValueError("factorial is defined only for non-negative values")
    if _is_int_like(xx) and xx <= 170:
        return float(math.factorial(int(round(xx))))
    return math.exp(_lgamma_cached(xx + 1.0))


def _log_factorial(x: float) -> float:
    if x < 0:
        return float("nan")
    if _is_int_like(x) and x <= 170:
        return math.log(float(math.factorial(int(round(x)))))
    return _lgamma_cached(x + 1.0)


def _mone(i: int) -> int:
    return 1 if (i & 1) == 0 else -1


def _safe_sqrt(x: float) -> float:
    if x < 0 and x > -1e-15:
        return 0.0
    return math.sqrt(x)


def _bit_count(value: int) -> int:
    return value.bit_count() if hasattr(value, "bit_count") else bin(value).count("1")


class Function1D:
    __slots__ = ("a",)

    def __init__(self, n: int, dtype=float, fill: float = 0.0):
        self.a = np.full((n,), fill, dtype=dtype)

    def __len__(self) -> int:
        return int(self.a.shape[0])

    def __getitem__(self, i: int):
        return self.a[i]

    def __setitem__(self, i: int, v):
        self.a[i] = v

    def view(self) -> np.ndarray:
        return self.a


class Function2D:
    __slots__ = ("a",)

    def __init__(self, n0: int, n1: int, dtype=float, fill: float = 0.0):
        self.a = np.full((n0, n1), fill, dtype=dtype, order="C")

    @property
    def shape(self) -> Tuple[int, int]:
        return int(self.a.shape[0]), int(self.a.shape[1])

    def __call__(self, i: int, j: int):
        return self.a[i, j]

    def set(self, i: int, j: int, v):
        self.a[i, j] = v

    def view(self) -> np.ndarray:
        return self.a


class ClebschGordan:
    def __init__(self, nn: int = 50, jn: int = 15):
        self.nn = int(nn)
        self.jn = int(jn)
        self._stored: Dict[int, float] = {}

    @staticmethod
    def encode_jm(j: float, m: float) -> int:
        return int(2.0 * j * (j + 1.0) + m + 1e-6)

    @staticmethod
    def decode_jm(i: int) -> Tuple[float, float]:
        kk = 0.25 * (-1.0 + math.sqrt(1.0 + 8.0 * (i + 1.0)))
        j = 0.5 * math.floor(2.0 * (kk - 0.0001))
        m = i - 2.0 * j * (j + 1.0)
        return j, m

    @staticmethod
    def encode_full(j: float, m: float, j1: float, m1: float, j2: float, m2: float) -> int:
        nn = 1000

        def pack(jx: float, mx: float) -> int:
            return int(2.0 * jx * (jx + 1.0) + mx + 1e-6)

        i = pack(j, m)
        i = i * nn + pack(j1, m1)
        i = i * nn + pack(j2, m2)
        return int(i)

    @staticmethod
    def _delta(j1: float, j2: float, j: float) -> float:
        a = j1 + j2 - j
        b = j1 - j2 + j
        c = -j1 + j2 + j
        d = j1 + j2 + j + 1.0
        return _safe_sqrt(math.exp(_log_factorial(a) + _log_factorial(b) + _log_factorial(c) - _log_factorial(d)))

    def cg(self, j: float, m: float, j1: float, m1: float, j2: float, m2: float) -> float:
        if (m1 + m2) != m:
            return 0.0
        if abs(m2) > j2 or abs(m1) > j1 or abs(m) > j:
            return 0.0

        tmin = int(max(0.0, j2 - j - m1, j1 - j + m2) + 1e-14)
        tmax = int(min(j1 + j2 - j, j1 - m1, j2 + m2) + 1e-14)
        if tmax < tmin:
            return 0.0

        log_v1 = 0.5 * (
            math.log(2.0 * j + 1.0)
            + _log_factorial(j1 + m1)
            + _log_factorial(j1 - m1)
            + _log_factorial(j2 + m2)
            + _log_factorial(j2 - m2)
            + _log_factorial(j + m)
            + _log_factorial(j - m)
        )

        sum_ = 0.0
        for t in range(tmin, tmax + 1):
            log_v2 = (
                _log_factorial(t)
                + _log_factorial(j1 + j2 - j - t)
                + _log_factorial(j1 - m1 - t)
                + _log_factorial(j2 + m2 - t)
                + _log_factorial(j - j2 + m1 + t)
                + _log_factorial(j - j1 - m2 + t)
            )
            sum_ += _mone(t) * math.exp(log_v1 - log_v2)

        return sum_ * self._delta(j1, j2, j)

    def cgs(self, j: float, m: float, j1: float, m1: float, j2: float, m2: float) -> float:
        if j > self.jn or j1 > self.jn or j2 > self.jn:
            return self.cg(j, m, j1, m1, j2, m2)

        ii = self.encode_full(j, m, j1, m1, j2, m2)
        cached = self._stored.get(ii)
        if cached is not None:
            return cached
        val = self.cg(j, m, j1, m1, j2, m2)
        self._stored[ii] = val
        return val

    def wigner_3j(self, j1: float, m1: float, j2: float, m2: float, j3: float, m3: float) -> float:
        if (m1 + m2 + m3) != 0.0:
            return 0.0
        if abs(j1 - j2) - 1e-14 > j3 or j3 > (j1 + j2) + 1e-14:
            return 0.0
        if abs(m1) > j1 or abs(m2) > j2 or abs(m3) > j3:
            return 0.0

        tmin = int(max(0.0, j2 - j3 - m1, j1 - j3 + m2) + 1e-14)
        tmax = int(min(j1 + j2 - j3, j1 - m1, j2 + m2) + 1e-14)
        if tmax < tmin:
            return 0.0

        sum_ = 0.0
        for t in range(tmin, tmax + 1):
            log_denom = (
                _log_factorial(t)
                + _log_factorial(j3 - j2 + m1 + t)
                + _log_factorial(j3 - j1 - m2 + t)
                + _log_factorial(j1 + j2 - j3 - t)
                + _log_factorial(j1 - m1 - t)
                + _log_factorial(j2 + m2 - t)
            )
            sum_ += _mone(t) * math.exp(-log_denom)

        log_dn = (
            _log_factorial(j1 + m1)
            + _log_factorial(j1 - m1)
            + _log_factorial(j2 + m2)
            + _log_factorial(j2 - m2)
            + _log_factorial(j3 + m3)
            + _log_factorial(j3 - m3)
        )
        pref = _mone(int(round(j1 - j2 - m3))) * self._delta(j1, j2, j3) * math.exp(0.5 * log_dn)
        return pref * sum_

    def gaunt(self, l1: int, m1: int, l2: int, m2: int, l3: int, m3: int) -> float:
        if l1 < 0 or l2 < 0 or l3 < 0:
            raise ValueError("Quantum number l must be non-negative")
        return (
            _mone(m1)
            * math.sqrt((2 * l1 + 1) * (2 * l2 + 1) * (2 * l3 + 1) / (4 * math.pi))
            * self.wigner_3j(l1, 0, l2, 0, l3, 0)
            * self.wigner_3j(l1, -m1, l2, m2, l3, m3)
        )


def compute_coulomb_ratios(l: int) -> List[float]:
    fkoj = [0.0] * (l + 1)
    if l == 1:
        fkoj[1] = 5.0
    elif l == 2:
        fkoj[1] = 14.0 / (1.0 + 0.625)
        fkoj[2] = 0.625 * fkoj[1]
    elif l == 3:
        fkoj[1] = 6435.0 / (286.0 + 195.0 * 0.668 + 250.0 * 0.494)
        fkoj[2] = 0.668 * fkoj[1]
        fkoj[3] = 0.494 * fkoj[1]
    return fkoj


class GauntCK:
    def __init__(self, l: int, cg: ClebschGordan):
        self.l = l
        self.ck = np.zeros((2 * l + 1, 2 * l + 1, l + 1), dtype=float)
        for m1 in range(-l, l + 1):
            for m2 in range(-l, l + 1):
                for kp, k in enumerate(range(0, 2 * l + 1, 2)):
                    c = cg.gaunt(l, m1, k, m1 - m2, l, m2) * math.sqrt(4.0 * math.pi / (2 * k + 1))
                    if abs(c) < 1e-10:
                        c = 0.0
                    self.ck[m1 + l, m2 + l, kp] = c

    def __call__(self, m1: int, m2: int, k: int) -> float:
        if abs(m1) > self.l or abs(m2) > self.l or k < 0 or k > 2 * self.l or (k % 2):
            raise ValueError("Gaunt coefficient does not exist")
        return float(self.ck[m1 + self.l, m2 + self.l, k // 2])


class AtomicDirectBasis:
    """Many-body helper in the |m_l,m_s> product basis.

    We keep this class as the exact backend for building the atomic Hamiltonian,
    because the Coulomb algebra from the original C++ code is already expressed
    naturally in the |m_l,m_s> basis.
    """
    def __init__(self, l: int):
        self.l = l
        self.baths = 2 * (2 * l + 1)
        self.l2p1 = 2 * l + 1
        self.ml_ms: List[Tuple[int, float]] = []
        for ms in (-0.5, 0.5):
            for ml in range(-l, l + 1):
                self.ml_ms.append((ml, ms))

        self.masks = [1 << i for i in range(self.baths)]
        self.mask_u = [0] * self.l2p1
        self.mask_d = [0] * self.l2p1
        self.lz_u = [0] * self.l2p1
        self.lz_d = [0] * self.l2p1
        self.ind_lz_u = [0] * self.l2p1
        self.ind_lz_d = [0] * self.l2p1
        self.lz = [0] * self.baths
        self.sz = [0] * self.baths
        self.lp_fact = [0.0] * self.l2p1
        self.lm_fact = [0.0] * self.l2p1

        iu = 0
        idn = 0
        ib = 0
        for i, (ml, ms) in enumerate(self.ml_ms):
            if ms > 0:
                self.lz_u[iu] = ml
                self.mask_u[iu] = self.masks[i]
                iu += 1
                self.lz[ib] = ml
                self.sz[ib] = int(round(2 * ms))
                ib += 1
            else:
                self.lz_d[idn] = ml
                self.mask_d[idn] = self.masks[i]
                idn += 1
                self.lz[ib] = ml
                self.sz[ib] = int(round(2 * ms))
                ib += 1

        for i in range(iu):
            self.ind_lz_u[self.lz_u[i] + l] = i
            self.ind_lz_d[self.lz_d[i] + l] = i
        for i in range(self.l2p1):
            m = i - l
            self.lp_fact[i] = math.sqrt(l * (l + 1) - m * (m + 1) + 0.0)
            self.lm_fact[i] = math.sqrt(l * (l + 1) - m * (m - 1) + 0.0)

    def estate(self, state: int, lz: int) -> int:
        lup = bool(state & self.mask_u[self.ind_lz_u[lz + self.l]])
        ldo = bool(state & self.mask_d[self.ind_lz_d[lz + self.l]])
        if lup and ldo:
            return 2
        if lup:
            return 1
        if ldo:
            return -1
        return 0

    def nel(self, state: int) -> int:
        return _bit_count(state)

    def sz_of(self, state: int) -> float:
        nu = 0
        nd = 0
        for ilz in range(self.l2p1):
            if self.mask_u[ilz] & state:
                nu += 1
            if self.mask_d[ilz] & state:
                nd += 1
        return 0.5 * (nu - nd)

    def lz_of(self, state: int) -> int:
        total = 0
        for ilz in range(self.l2p1):
            if self.mask_u[ilz] & state:
                total += self.lz_u[ilz]
            if self.mask_d[ilz] & state:
                total += self.lz_d[ilz]
        return total

    def jz_of(self, state: int) -> float:
        return self.lz_of(state) + self.sz_of(state)

    def sign_between(self, state: int, mask_min: int, mask_max: int) -> int:
        mask = mask_min << 1
        n = 0
        while mask < mask_max:
            if mask & state:
                n += 1
            mask <<= 1
        return 1 - 2 * (n & 1)

    def sign_before(self, state: int, mask: int) -> int:
        n = _bit_count(state & (mask - 1))
        return 1 - 2 * (n & 1)

    def n_el_before(self, state: int, i: int) -> int:
        return _bit_count(state & ((1 << i) - 1))

    def ln_sm(self, ln: int, sm: int, state: int) -> Dict[int, float]:
        sts0: Dict[int, float] = {state: 1.0}
        for _ in range(max(0, ln)):
            sts: Dict[int, float] = {}
            for st, fact in sts0.items():
                for ilz in range(self.l2p1 - 1):
                    if (self.mask_u[ilz] & st) and not (self.mask_u[ilz + 1] & st):
                        newst = st ^ self.mask_u[ilz] ^ self.mask_u[ilz + 1]
                        sts[newst] = sts.get(newst, 0.0) + fact * self.lp_fact[ilz]
                    if (self.mask_d[ilz] & st) and not (self.mask_d[ilz + 1] & st):
                        newst = st ^ self.mask_d[ilz] ^ self.mask_d[ilz + 1]
                        sts[newst] = sts.get(newst, 0.0) + fact * self.lp_fact[ilz]
            sts0 = sts
        for _ in range(max(0, -ln)):
            sts = {}
            for st, fact in sts0.items():
                for ilz in range(1, self.l2p1):
                    if (self.mask_u[ilz] & st) and not (self.mask_u[ilz - 1] & st):
                        newst = st ^ self.mask_u[ilz] ^ self.mask_u[ilz - 1]
                        sts[newst] = sts.get(newst, 0.0) + fact * self.lm_fact[ilz]
                    if (self.mask_d[ilz] & st) and not (self.mask_d[ilz - 1] & st):
                        newst = st ^ self.mask_d[ilz] ^ self.mask_d[ilz - 1]
                        sts[newst] = sts.get(newst, 0.0) + fact * self.lm_fact[ilz]
            sts0 = sts
        for _ in range(max(0, sm)):
            sts = {}
            for st, fact in sts0.items():
                for ilz in range(self.l2p1):
                    if (self.mask_d[ilz] & st) and not (self.mask_u[ilz] & st):
                        newst = st ^ self.mask_d[ilz] ^ self.mask_u[ilz]
                        sign = self.sign_between(st, self.mask_d[ilz], self.mask_u[ilz])
                        sts[newst] = sts.get(newst, 0.0) + fact * sign
            sts0 = sts
        for _ in range(max(0, -sm)):
            sts = {}
            for st, fact in sts0.items():
                for ilz in range(self.l2p1):
                    if (self.mask_u[ilz] & st) and not (self.mask_d[ilz] & st):
                        newst = st ^ self.mask_u[ilz] ^ self.mask_d[ilz]
                        sign = self.sign_between(st, self.mask_d[ilz], self.mask_u[ilz])
                        sts[newst] = sts.get(newst, 0.0) + fact * sign
            sts0 = sts
        return sts0

    def j2(self, state: int) -> Dict[int, float]:
        sts0: Dict[int, float] = {state: 2.0 * self.lz_of(state) * self.sz_of(state)}
        for src in (self.l2(state), self.s2(state), self.ln_sm(1, -1, state), self.ln_sm(-1, 1, state)):
            if isinstance(src, dict):
                items = src.items()
            else:
                items = src
            for st, amp in items:
                sts0[st] = sts0.get(st, 0.0) + amp
        return sts0

    def s2(self, state: int) -> List[Tuple[int, float]]:
        sts: List[Tuple[int, float]] = []
        dd = 0.0
        for ilz in range(self.l2p1):
            up = 1 if (self.mask_u[ilz] & state) else 0
            dn = 1 if (self.mask_d[ilz] & state) else 0
            if up + dn == 1:
                dd += 0.5
        sts.append((state, self.sz_of(state) ** 2 + dd))
        for ilz in range(self.l2p1):
            im1 = self.mask_u[ilz]
            im2 = self.mask_d[ilz]
            ib1 = bool(state & im1)
            ib2 = bool(state & im2)
            if ib1 and not ib2:
                isig = self.sign_between(state, min(im1, im2), max(im1, im2))
                istate = state ^ im1 ^ im2
                for jlz in range(self.l2p1):
                    if ilz == jlz:
                        continue
                    jm1 = self.mask_d[jlz]
                    jm2 = self.mask_u[jlz]
                    jb1 = bool(state & jm1)
                    jb2 = bool(state & jm2)
                    if jb1 and not jb2:
                        jsig = self.sign_between(istate, min(jm1, jm2), max(jm1, jm2))
                        jstate = istate ^ jm1 ^ jm2
                        sts.append((jstate, isig * jsig))
        return sts

    def _l2_off_diag(
        self,
        state: int,
        maski: Sequence[int],
        lzi: Sequence[int],
        maskj: Sequence[int],
        lzj: Sequence[int],
        cont: bool,
    ) -> List[Tuple[int, float]]:
        out: List[Tuple[int, float]] = []
        for ilz in range(1, self.l2p1):
            im1 = maski[ilz]
            im2 = maski[ilz - 1]
            if (state & im1) and not (state & im2):
                m = lzi[ilz]
                ifct = math.sqrt(self.l * (self.l + 1) - m * (m - 1)) * self.sign_between(state, min(im1, im2), max(im1, im2))
                istate = state ^ im1 ^ im2
                for jlz in range(self.l2p1 - 1):
                    jm1 = maskj[jlz]
                    jm2 = maskj[jlz + 1]
                    if cont and (ilz == jlz or im2 == jm2):
                        continue
                    if (state & jm1) and not (state & jm2):
                        m2 = lzj[jlz]
                        jfct = math.sqrt(self.l * (self.l + 1) - m2 * (m2 + 1)) * self.sign_between(istate, min(jm1, jm2), max(jm1, jm2))
                        jstate = istate ^ jm1 ^ jm2
                        out.append((jstate, ifct * jfct))
        return out

    def l2(self, state: int) -> List[Tuple[int, float]]:
        dd = 0.0
        for ilz in range(1, self.l2p1):
            if (self.mask_u[ilz] & state) and not (self.mask_u[ilz - 1] & state):
                dd += 0.5 * (self.l * (self.l + 1) - self.lz_u[ilz] * (self.lz_u[ilz] - 1))
            if (self.mask_d[ilz] & state) and not (self.mask_d[ilz - 1] & state):
                dd += 0.5 * (self.l * (self.l + 1) - self.lz_d[ilz] * (self.lz_d[ilz] - 1))
        for ilz in range(self.l2p1 - 1):
            if (self.mask_u[ilz] & state) and not (self.mask_u[ilz + 1] & state):
                dd += 0.5 * (self.l * (self.l + 1) - self.lz_u[ilz] * (self.lz_u[ilz] + 1))
            if (self.mask_d[ilz] & state) and not (self.mask_d[ilz + 1] & state):
                dd += 0.5 * (self.l * (self.l + 1) - self.lz_d[ilz] * (self.lz_d[ilz] + 1))
        sts = [(state, self.lz_of(state) ** 2 + dd)]
        sts.extend(self._l2_off_diag(state, self.mask_u, self.lz_u, self.mask_u, self.lz_u, True))
        sts.extend(self._l2_off_diag(state, self.mask_d, self.lz_d, self.mask_d, self.lz_d, True))
        sts.extend(self._l2_off_diag(state, self.mask_u, self.lz_u, self.mask_d, self.lz_d, False))
        sts.extend(self._l2_off_diag(state, self.mask_d, self.lz_d, self.mask_u, self.lz_u, False))
        return sts

    def li_dot_si(self, state: int) -> Dict[int, float]:
        sts0: Dict[int, float] = {}
        diag = 0.0
        for i in range(self.baths):
            if self.masks[i] & state:
                diag += self.lz[i] * 0.5 * self.sz[i]
        sts0[state] = diag
        for ilz in range(1, self.l2p1):
            if (self.mask_d[ilz] & state) and not (self.mask_u[ilz - 1] & state):
                newst = state ^ self.mask_d[ilz] ^ self.mask_u[ilz - 1]
                sts0[newst] = sts0.get(newst, 0.0) + 0.5 * self.lm_fact[ilz] * self.sign_between(state, self.mask_d[ilz], self.mask_u[ilz - 1])
        for ilz in range(self.l2p1 - 1):
            if (self.mask_u[ilz] & state) and not (self.mask_d[ilz + 1] & state):
                newst = state ^ self.mask_u[ilz] ^ self.mask_d[ilz + 1]
                sts0[newst] = sts0.get(newst, 0.0) + 0.5 * self.lp_fact[ilz] * self.sign_between(state, self.mask_d[ilz + 1], self.mask_u[ilz])
        return sts0

    def coulomb_u(self, state: int, gck: GauntCK, fkoj: Sequence[float]) -> List[Tuple[int, Tuple[float, float]]]:
        sts: List[Tuple[int, Tuple[float, float]]] = []
        ni = -1
        for i in range(self.baths):
            if not (self.masks[i] & state):
                continue
            ni += 1
            state1 = state ^ self.masks[i]
            m1 = self.lz[i]
            s1 = self.sz[i]
            nj = -1
            for j in range(self.baths):
                if not (self.masks[j] & state1):
                    continue
                nj += 1
                state2 = state1 ^ self.masks[j]
                m2 = self.lz[j]
                s2 = self.sz[j]
                for a in range(self.baths):
                    if (self.masks[a] & state2) or self.sz[a] != s2:
                        continue
                    na = self.n_el_before(state2, a)
                    state3 = state2 ^ self.masks[a]
                    m3 = self.lz[a]
                    for b in range(self.baths):
                        if (self.masks[b] & state3) or self.sz[b] != s1:
                            continue
                        m4 = self.lz[b]
                        if m4 - m1 != m2 - m3:
                            continue
                        nb = self.n_el_before(state3, b)
                        state4 = state3 ^ self.masks[b]
                        sign = 1 - 2 * ((ni + nj + na + nb) & 1)
                        u0 = sign * gck(m4, m1, 0) * gck(m2, m3, 0)
                        us = 0.0
                        for k in range(1, self.l + 1):
                            us += gck(m4, m1, 2 * k) * gck(m2, m3, 2 * k) * fkoj[k]
                        sts.append((state4, (u0, sign * us)))
        return sts

    def fp(self, state: int, j: float, jz: float, cg: ClebschGordan) -> List[Tuple[int, float]]:
        out: List[Tuple[int, float]] = []
        sz = 0.5
        lz = int(round(jz - 0.5))
        if abs(lz) <= self.l and not (self.mask_u[self.l + lz] & state):
            mask = self.mask_u[self.l + lz]
            out.append((state ^ mask, self.sign_before(state, mask) * cg.cgs(j, jz, self.l, lz, 0.5, sz)))
        sz = -0.5
        lz = int(round(jz + 0.5))
        if abs(lz) <= self.l and not (self.mask_d[self.l + lz] & state):
            mask = self.mask_d[self.l + lz]
            out.append((state ^ mask, self.sign_before(state, mask) * cg.cgs(j, jz, self.l, lz, 0.5, sz)))
        return out

    def one_body(self, state: int, h1: np.ndarray) -> Dict[int, float]:
        out: Dict[int, float] = {}
        n_orb = h1.shape[0]
        for b in range(n_orb):
            if not (state & self.masks[b]):
                continue
            state_removed = state ^ self.masks[b]
            sign_remove = self.sign_before(state, self.masks[b])
            for a in range(n_orb):
                hab = float(h1[a, b])
                if abs(hab) < 1e-14:
                    continue
                if a == b:
                    out[state] = out.get(state, 0.0) + hab
                    continue
                if state_removed & self.masks[a]:
                    continue
                sign_add = self.sign_before(state_removed, self.masks[a])
                new_state = state_removed ^ self.masks[a]
                out[new_state] = out.get(new_state, 0.0) + hab * sign_remove * sign_add
        return out

    def printn(self, state: int) -> str:
        chars: List[str] = []
        for lz in range(-self.l, self.l + 1):
            e = self.estate(state, lz)
            chars.append({0: "0", 2: "2", 1: "u", -1: "d"}[e])
        return "".join(chars)


class JmProductBasis:
    """Direct-product basis in the |j,m_j> single-particle basis.

    Bold-CTQMC expects the many-body basis to be a Fock basis built from the
    same one-particle orbitals that appear in the bath definition. In the atomic
    spin-orbit case those orbitals are |j,m_j>, so F^+ must act as a pure
    bit-flip in this basis.
    """

    def __init__(self, one_electron_basis: Sequence[Tuple[float, float]]):
        self.orbitals = list(one_electron_basis)
        self.baths = len(self.orbitals)
        self.masks = [1 << i for i in range(self.baths)]
        self.mj = [mj for _, mj in self.orbitals]
        self.index_of = {(j, mj): i for i, (j, mj) in enumerate(self.orbitals)}

    def nel(self, state: int) -> int:
        return _bit_count(state)

    def jz_of(self, state: int) -> float:
        total = 0.0
        for i in range(self.baths):
            if state & self.masks[i]:
                total += self.mj[i]
        return total

    def sign_before(self, state: int, mask: int) -> int:
        n = _bit_count(state & (mask - 1))
        return 1 - 2 * (n & 1)

    def one_body(self, state: int, h1: np.ndarray) -> Dict[int, float]:
        # This occupation-bit operation is independent of orbital quantum numbers.
        return AtomicDirectBasis.one_body(self, state, h1)

    def fp(self, state: int, j: float, jz: float, _cg: ClebschGordan) -> List[Tuple[int, float]]:
        idx = self.index_of[(j, jz)]
        mask = self.masks[idx]
        if state & mask:
            return []
        return [(state ^ mask, float(self.sign_before(state, mask)))]

    def printn(self, state: int) -> str:
        return "".join("1" if state & self.masks[i] else "0" for i in range(self.baths))


def build_one_electron_basis(l: int) -> List[Tuple[float, float]]:
    out: List[Tuple[float, float]] = []
    for ic in (-1, 1):
        jc = l + 0.5 * ic
        if jc < 0:
            continue
        j2max = int(round(2 * jc))
        for j2z in range(-j2max, j2max + 1, 2):
            out.append((jc, j2z / 2.0))
    return out


def build_jm_to_mlms_transform(
    l: int,
    ml_ms: Sequence[Tuple[int, float]],
    one_electron_basis: Sequence[Tuple[float, float]],
    cg: ClebschGordan,
) -> np.ndarray:
    transform = np.zeros((len(one_electron_basis), len(ml_ms)), dtype=float)
    for ia, (j, jz) in enumerate(one_electron_basis):
        for ib, (ml, ms) in enumerate(ml_ms):
            transform[ia, ib] = cg.cgs(j, jz, l, ml, 0.5, ms)
    return transform


def occupied_orbitals(state: int, baths: int) -> List[int]:
    return [i for i in range(baths) if state & (1 << i)]


def build_sector_transform(
    jm_to_mlms: np.ndarray,
    basis_jm: Sequence[int],
    basis_mlms: Sequence[int],
    baths: int,
) -> np.ndarray:
    """Build the exact sector unitary from |j,m_j> Slater states to |m_l,m_s>.

    If a one-particle basis changes by the unitary matrix U, an N-particle
    Slater determinant transforms by determinants of the occupied submatrices of
    U. We use that sector-by-sector to obtain the exact Hamiltonian in the
    |j,m_j> product basis while reusing the tested |m_l,m_s> backend.
    """

    transform = np.zeros((len(basis_jm), len(basis_mlms)), dtype=float)
    for i, state_jm in enumerate(basis_jm):
        occ_jm = occupied_orbitals(state_jm, baths)
        if not occ_jm:
            transform[i, 0] = 1.0
            continue
        for j, state_mlms in enumerate(basis_mlms):
            occ_mlms = occupied_orbitals(state_mlms, baths)
            if len(occ_mlms) != len(occ_jm):
                continue
            sub = jm_to_mlms[np.ix_(occ_jm, occ_mlms)]
            transform[i, j] = float(np.linalg.det(sub))
    return transform


def default_impqs(one_electron_basis: Sequence[Tuple[float, float]]) -> List[int]:
    impqs: List[int] = []
    cur_j = None
    cur_ind = -1
    for j, _ in one_electron_basis:
        if cur_j is None or j != cur_j:
            cur_j = j
            cur_ind += 1
        impqs.append(cur_ind)
    return impqs


def compute_global_flip(one_electron_basis: Sequence[Tuple[float, float]]) -> List[int]:
    global_flip = [0] * len(one_electron_basis)
    kk = 0
    ii = 0
    start = 0
    while start < len(one_electron_basis):
        j = one_electron_basis[start][0]
        stop = start
        while stop < len(one_electron_basis) and one_electron_basis[stop][0] == j:
            stop += 1
        j2z_vals = [int(round(2 * mj)) for _, mj in one_electron_basis[start:stop]]
        mid = len(j2z_vals) // 2
        local = []
        for _ in range(mid + 1):
            local.append(ii)
            ii += 1
        for _ in range(mid + 1, len(j2z_vals)):
            ii -= 1
            local.append(ii)
        for offset, val in enumerate(local):
            global_flip[start + offset] = val
        ii = max(ii, max(global_flip[:stop])) + 1
        kk = stop
        start = stop
    return global_flip


def parse_list_arg(text: str, cast) -> List:
    stream = text.strip()
    if not stream:
        return []
    if stream[0] in "[(" and stream[-1] in "])":
        stream = stream[1:-1]
    if not stream:
        return []
    return [cast(x.strip()) for x in stream.split(",") if x.strip()]


def construct_reduced_base_njz(op: AtomicDirectBasis, n: int, jz: float) -> List[int]:
    return [state for state in range(1 << op.baths) if op.nel(state) == n and abs(op.jz_of(state) - jz) < 1e-10]


def _angular_operator_matrix(op: AtomicDirectBasis, basis: Sequence[int], name: str) -> np.ndarray:
    """Return a many-electron L², S², or J² matrix in a fixed (N,Jz) sector."""
    index = {state: i for i, state in enumerate(basis)}
    matrix = np.zeros((len(basis), len(basis)), dtype=float)
    apply = getattr(op, name)
    for col, state in enumerate(basis):
        terms = apply(state)
        for dest, amplitude in (terms.items() if isinstance(terms, dict) else terms):
            row = index.get(dest)
            if row is not None:
                matrix[row, col] += amplitude
    return 0.5 * (matrix + matrix.T)


def _twice_angular_momentum(square: float) -> int:
    """Convert the eigenvalue q(q+1) to the integer 2q, checking its accuracy."""
    twice_q = int(round(math.sqrt(max(0.0, 1.0 + 4.0 * square)) - 1.0))
    q = twice_q / 2.0
    if abs(square - q * (q + 1.0)) > 1e-6:
        raise ValueError(f"Angular-momentum operator has nonquantized eigenvalue {square}")
    return twice_q


def build_nlsj_basis(op: AtomicDirectBasis, basis: Sequence[int], n: int) -> Tuple[np.ndarray, List[Tuple[int, int, int, int]]]:
    """Simultaneously diagonalize commuting L², S², J² in the |m_l,m_s> basis.

    Degenerate vectors with the same (N,L,S,J) can be chosen arbitrarily: the
    solver sums their *squared* overlaps, so only the whole labeled subspace
    matters. The tuple stores twice each angular momentum to represent half
    integers exactly in the sidecar file.
    """
    dimension = len(basis)
    operators = {name: _angular_operator_matrix(op, basis, name) for name in ("l2", "s2", "j2")}
    groups: List[Tuple[np.ndarray, Tuple[int, ...]]] = [(np.eye(dimension), (n,))]
    for name in ("l2", "s2", "j2"):
        next_groups: List[Tuple[np.ndarray, Tuple[int, ...]]] = []
        for vectors, labels in groups:
            reduced = vectors.T @ operators[name] @ vectors
            eigenvalues, rotation = np.linalg.eigh(0.5 * (reduced + reduced.T))
            twice_q = [_twice_angular_momentum(float(value)) for value in eigenvalues]
            for quantum_number in sorted(set(twice_q)):
                columns = [i for i, value in enumerate(twice_q) if value == quantum_number]
                next_groups.append((vectors @ rotation[:, columns], labels + (quantum_number,)))
        groups = next_groups
    columns = [vectors[:, i] for vectors, _ in groups for i in range(vectors.shape[1])]
    labels = [label for vectors, label in groups for _ in range(vectors.shape[1])]
    result = np.column_stack(columns)
    if result.shape != (dimension, dimension) or np.max(np.abs(result.T @ result - np.eye(dimension))) > 1e-8:
        raise ValueError("The constructed NLSJ basis is not complete and orthonormal")
    return result, labels


def _fnv1a_file(path: str) -> str:
    """A portable file identifier, checked by bctqmc before using the sidecar."""
    value = 14695981039346656037
    with open(path, "rb") as source:
        while chunk := source.read(1 << 20):
            for byte in chunk:
                value = ((value ^ byte) * 1099511628211) & 0xFFFFFFFFFFFFFFFF
    return f"{value:016x}"


def write_nlsj_sidecar(
    output: str,
    sectors: Sequence[Dict[str, object]],
    natst: int,
    orbitals: Sequence[Tuple[float, float]],
) -> None:
    """List physical energy eigenstates first, then their complete CIX-basis vectors.

    The CIX writes each Hamiltonian entry with ten decimal places. Diagonalize
    that same rounded matrix here, so the reported spectrum matches bctqmc's
    input. The solver receives all six-significant-digit vector coefficients,
    including small ones; no occupation annotations are needed in this format.
    """
    eigensystems = []
    for sector in sectors:
        ham = np.round(np.asarray(sector["ham"], dtype=float), 10)
        energies, vectors = np.linalg.eigh(ham)
        assignments = []
        qn_vectors = sector.get("qn_vectors")
        if qn_vectors is None:
            assignments = [(None, 0.0)] * len(energies)
        else:
            overlaps = np.asarray(qn_vectors, dtype=float).T @ vectors
            for eigen in range(len(energies)):
                weights: Dict[Tuple[int, int, int, int], float] = defaultdict(float)
                for column, qn in enumerate(sector["qn_labels"]):
                    weights[tuple(qn)] += float(abs(overlaps[column, eigen]) ** 2)
                total = sum(weights.values())
                if abs(total - 1.0) > 1e-5:
                    raise ValueError(f"Incomplete NLSJ projection in block {sector['block_id']}")
                # Sorting makes ties deterministic; label_weight exposes cases
                # where the one printed label is not a good quantum number.
                best = max(sorted(weights), key=lambda qn: weights[qn])
                assignments.append((best, min(1.0, weights[best] / total)))
        eigensystems.append((energies, vectors, assignments))

    with open(output + ".nlsj", "w", encoding="ascii") as sidecar:
        print("# NLSJ_EIGENSTATES_V3", file=sidecar)
        print(natst, len(sectors), len(orbitals), _fnv1a_file(output),
              "# nstates nblocks norbitals cix_fingerprint", file=sidecar)
        print(" ".join(str(len(sector["basis"])) for sector in sectors),
              "# block dimensions, in CIX order", file=sidecar)
        jz_heading = "Jz_basis" if any("basis_jz" in sector for sector in sectors) else "Jz"
        print(f"# {'idx':>4} {'blk':>3} {'N':>2} {'L':>3} {'S':>3} {'J':>3} "
              f"{jz_heading:>4} {'E_CIX':>18} {'label_weight':>12}", file=sidecar)
        for sector, (energies, _, assignments) in zip(sectors, eigensystems):
            for local, (label, weight) in enumerate(assignments):
                if label is None:
                    quantum_numbers = (-1, -1, -1)
                else:
                    quantum_numbers = (label[1] / 2, label[2] / 2, label[3] / 2)
                l_value, s_value, j_value = quantum_numbers
                # The V3 reader checks this field against the corresponding
                # CIX product-state row. It is not a conserved eigenstate Jz
                # when the block contains several Jz sectors.
                row_jz = sector["basis_jz"][local] if "basis_jz" in sector else sector["jz"]
                line = (f"  {sector['state_ids'][local]:4d} {sector['block_id']:3d} "
                        f"{sector['n']:2d} {l_value:3g} {s_value:3g} {j_value:3g} "
                        f"{row_jz:4g} {energies[local] + 0.0:18.12g} "
                        f"{weight:12.8g}")
                print(line + (" # unlabeled" if label is None else ""), file=sidecar)

        print("# state eigenvector_components_in_CIX_block_order", file=sidecar)
        for sector, (_, vectors, _) in zip(sectors, eigensystems):
            for column in range(len(sector["basis"])):
                components = vectors[:, column]
                print(f"  {sector['state_ids'][column]:4d} " +
                      " ".join(f"{value:.6g}" for value in components),
                      file=sidecar)


def label_existing_cix(path: str, l: int) -> None:
    """Create only the NLSJ sidecar for an existing atomh.py-style CIX file.

    Its state lines already contain the |j,m_j> occupation bitstrings. Reading
    those makes the sidecar match the original CIX byte-for-byte: no atomic
    energies, interaction parameters, or impurity-solver input is regenerated.
    Blocks that are strict subsets of a (N,Jz) sector are left unlabeled rather
    than falsely projecting angular-momentum eigenvectors into them.
    """
    with open(path, "r", encoding="utf-8") as cix:
        cix.readline()
        cix.readline()
        header = cix.readline().split()
        if len(header) < 6:
            raise ValueError("Not an atomh.py-style bctqmc.cix header")
        natst, npsi, nbath, nblock = map(int, (header[1], header[2], header[3], header[4]))
        if npsi != 2 * (2 * l + 1) or nbath < npsi:
            raise ValueError(f"CIX has {npsi} operators, inconsistent with shell l={l}")
        cix.readline()  # bath heading
        for _ in range(nbath):
            cix.readline()
        cix.readline()  # energy heading
        cix.readline()  # bath energies
        cix.readline()  # state heading
        blocks: List[Dict[str, object]] = [
            {"block_id": ib + 1, "basis": [], "basis_jz": [], "state_ids": [], "n": None, "jz": None}
            for ib in range(nblock)
        ]
        for expected_id in range(1, natst + 1):
            line = cix.readline()
            if "#" not in line:
                raise ValueError(f"Missing occupation bitstring for CIX state {expected_id}")
            left, comment = line.split("#", 1)
            fields = left.split()
            if len(fields) < 5 + npsi or int(fields[0]) != expected_id:
                raise ValueError(f"Malformed CIX state {expected_id}")
            block_id, n, jz = int(fields[1]), int(fields[2]), float(fields[4])
            if not 1 <= block_id <= nblock:
                raise ValueError(f"Invalid block for CIX state {expected_id}")
            bits = comment.split()[0]
            if len(bits) != npsi or any(bit not in "01" for bit in bits):
                raise ValueError(f"Invalid occupation bitstring for CIX state {expected_id}")
            state = sum((1 << index) for index, bit in enumerate(bits) if bit == "1")
            block = blocks[block_id - 1]
            if block["n"] is None:
                block["n"], block["jz"] = n, jz
            elif block["n"] != n:
                raise ValueError(f"Block {block_id} mixes particle numbers")
            block["basis"].append(state)
            block["basis_jz"].append(jz)
            block["state_ids"].append(expected_id)

        if not cix.readline().startswith("# All blocks of local Hamiltonian"):
            raise ValueError("Missing CIX Hamiltonian section")
        for block in blocks:
            dimension = len(block["basis"])
            header = cix.readline().split("#", 1)[0].split()
            if (len(header) < 2 + dimension or
                    int(header[0]) != block["block_id"] or
                    int(header[1]) != dimension or
                    [int(x) for x in header[2:2 + dimension]] != block["state_ids"]):
                raise ValueError(f"Invalid CIX Hamiltonian header for block {block['block_id']}")
            rows = []
            for _ in range(dimension):
                fields = cix.readline().split("#", 1)[0].split()
                if len(fields) != dimension:
                    raise ValueError(f"Incomplete CIX Hamiltonian block {block['block_id']}")
                rows.append([float(x) for x in fields])
            block["ham"] = np.asarray(rows, dtype=float)
            if not np.allclose(block["ham"], block["ham"].T, rtol=0, atol=1e-8):
                raise ValueError(f"Non-Hermitian CIX Hamiltonian block {block['block_id']}")

    op_mlms = AtomicDirectBasis(l)
    op_jm = JmProductBasis(build_one_electron_basis(l))
    cg = ClebschGordan()
    jm_to_mlms = build_jm_to_mlms_transform(l, op_mlms.ml_ms, build_one_electron_basis(l), cg)
    for block in blocks:
        n = int(block["n"])
        basis_jm = block["basis"]
        block["qn_vectors"] = np.zeros((len(basis_jm), len(basis_jm)))
        block["qn_labels"] = []
        for jz in sorted(set(block["basis_jz"])):
            rows = [i for i, value in enumerate(block["basis_jz"]) if value == jz]
            subset = [basis_jm[i] for i in rows]
            full_basis_jm = construct_reduced_base_njz(op_jm, n, jz)
            if set(subset) != set(full_basis_jm):
                block["qn_vectors"] = None
                break
            basis_mlms = construct_reduced_base_njz(op_mlms, n, jz)
            transform = build_sector_transform(jm_to_mlms, subset, basis_mlms, op_jm.baths)
            qn_mlms, labels = build_nlsj_basis(op_mlms, basis_mlms, n)
            column = len(block["qn_labels"])
            block["qn_vectors"][np.ix_(rows, range(column, column + len(rows)))] = transform @ qn_mlms
            block["qn_labels"].extend(labels)
        if len(set(block["basis_jz"])) == 1:
            del block["basis_jz"]
    write_nlsj_sidecar(path, blocks, natst, build_one_electron_basis(l))


def connected_components_from_matrix(mat: np.ndarray, tol: float = 1e-10) -> List[List[int]]:
    n = mat.shape[0]
    seen = [False] * n
    comps: List[List[int]] = []
    for i in range(n):
        if seen[i]:
            continue
        comp: List[int] = []
        q: deque[int] = deque([i])
        seen[i] = True
        while q:
            u = q.popleft()
            comp.append(u)
            row = np.where(np.abs(mat[u, :]) > tol)[0]
            col = np.where(np.abs(mat[:, u]) > tol)[0]
            for v in np.unique(np.concatenate((row, col))):
                if not seen[int(v)]:
                    seen[int(v)] = True
                    q.append(int(v))
        comps.append(comp)
    return comps


def build_hamiltonian_block(
    op: AtomicDirectBasis,
    basis: Sequence[int],
    j_coulomb: float,
    h1: np.ndarray,
    gck: GauntCK,
    fkoj: Sequence[float],
) -> np.ndarray:
    index = {state: i for i, state in enumerate(basis)}
    ham = np.zeros((len(basis), len(basis)), dtype=float)
    for col, state in enumerate(basis):
        for new_state, amp in op.one_body(state, h1).items():
            row = index.get(new_state)
            if row is not None:
                ham[row, col] += amp
        for new_state, (_, u1) in op.coulomb_u(state, gck, fkoj):
            row = index.get(new_state)
            if row is not None:
                ham[row, col] += 0.5 * j_coulomb * u1
    return ham


def build_jz_grid(n: int, l: int) -> List[float]:
    jmax2 = int(round(2 * (l + 0.5) * n))
    parity = n % 2
    return [j2 / 2.0 for j2 in range(-jmax2, jmax2 + 1) if abs(j2) % 2 == parity]


def default_epsk(impqs: Sequence[int], eimp: Sequence[float]) -> List[float]:
    if not eimp:
        n = (max(impqs) + 1) if impqs else 0
        return [0.0] * n
    return list(eimp)


def read_trans_sigind(path: str) -> np.ndarray:
    """Read the one-based Sigind section of the wrapper's Trans.dat."""
    with open(path, encoding="utf-8") as source:
        data = (line.split("#", 1)[0].strip() for line in source)
        data = (line for line in data if line)
        dims = [int(x) for x in next(data).split()]
        if len(dims) != 2 or dims[0] != dims[1] or dims[0] <= 0:
            raise ValueError(f"Invalid Sigind dimensions in {path}")
        return np.asarray([[int(x) for x in next(data).split()]
                           for _ in range(dims[0])], dtype=int)


def prepare_hybridization(
    orbitals: Sequence[Tuple[float, float]],
    impq: Optional[Sequence[int]],
    eimp: Optional[Sequence[float]],
    sigind: Optional[Sequence[Sequence[int]]],
) -> Tuple[List[int], List[float], np.ndarray, List[Tuple[int, int, int]]]:
    """Separate physical operators, directed bath entries, and unique functions.

    Sigind uses the eDMFT convention: zero is absent and positive labels are
    one-based. CIX function labels and legacy Impq entries are zero-based.
    Diagonal bath entries stay first to preserve historical operator ordering.
    """
    npsi = len(orbitals)
    if sigind is None:
        impqs = list(impq) if impq is not None else default_impqs(orbitals)
        if len(impqs) != npsi or any(q < 0 for q in impqs):
            raise ValueError(f"Impq must contain {npsi} nonnegative indices")
        mapping = np.diag(np.asarray(impqs, dtype=int) + 1)
    else:
        raw = np.asarray(sigind)
        if raw.shape != (npsi, npsi) or not np.all(np.isfinite(raw)):
            raise ValueError(f"Sigind must be a finite {npsi} by {npsi} matrix")
        if not np.all(raw == np.asarray(raw, dtype=int)) or np.any(raw < 0):
            raise ValueError("Sigind requires nonnegative integer labels; zero means absent")
        mapping = np.asarray(raw, dtype=int)
        if np.any(np.diag(mapping) <= 0):
            raise ValueError("Every physical orbital needs a positive diagonal Sigind label")
        labels = set(int(q) for q in mapping.flat if q)
        if labels != set(range(1, int(mapping.max()) + 1)):
            raise ValueError("Sigind component labels must be contiguous starting at 1")
        impqs = [int(q) - 1 for q in np.diag(mapping)]
        if impq is not None and list(impq) != impqs:
            raise ValueError("Impq disagrees with the diagonal of Sigind")
    ncomponents = int(mapping.max())
    epsk = default_epsk(impqs, list(eimp) if eimp is not None else [])
    if len(epsk) > ncomponents:
        raise ValueError("Eimp has more entries than the hybridization component mapping")
    epsk.extend([0.0] * (ncomponents - len(epsk)))
    if np.any(np.abs(np.asarray(epsk, dtype=complex).imag) > 1e-14):
        raise ValueError("Complex impurity levels are not supported by the real CIX Hamiltonian")
    epsk = [float(complex(value).real) for value in epsk]
    h1 = np.zeros((npsi, npsi), dtype=float)
    present = mapping > 0
    h1[present] = np.asarray(epsk)[mapping[present] - 1]
    if not np.all(np.isfinite(h1)) or not np.allclose(h1, h1.T, rtol=0, atol=1e-12):
        raise ValueError("Eimp and Sigind must define a real symmetric impurity Hamiltonian")
    baths = [(i, i, int(mapping[i, i]) - 1) for i in range(npsi)]
    baths.extend((i, j, int(mapping[i, j]) - 1)
                 for i in range(npsi) for j in range(npsi)
                 if i != j and mapping[i, j] > 0)
    return impqs, epsk, h1, baths


def merge_hybridization_sectors(
    sectors: Sequence[Dict[str, object]],
    orbitals: Sequence[Tuple[float, float]],
    baths: Sequence[Tuple[int, int, int]],
    h1: np.ndarray,
    op: JmProductBasis,
) -> List[Dict[str, object]]:
    """Retain Jz modulo the transfers allowed by hybridization and Eimp.

    Slater interactions conserve Jz. A bath entry (a,b) allows transfer
    mj[a]-mj[b], even when its initial self-energy or impurity level is zero.
    The gcd of twice those transfers determines the remaining symmetry.
    Full blocks of that symmetry also give each physical operator a unique
    target block, as required by the CIX reader.
    """
    modulus = 0
    for a, b, _ in baths:
        if a != b:
            modulus = math.gcd(modulus, abs(int(round(2 * (orbitals[a][1] - orbitals[b][1])))))
    off = h1 - np.diag(np.diag(h1))
    for a, b in zip(*np.nonzero(np.abs(off) > 1e-14)):
        modulus = math.gcd(modulus, abs(int(round(2 * (orbitals[a][1] - orbitals[b][1])))))
    if not any(a != b for a, b, _ in baths):
        return list(sectors)
    # A zero transfer can mix orbitals with equal mj without merging Jz sectors.
    groups: Dict[Tuple[int, int], List[Dict[str, object]]] = defaultdict(list)
    for sector in sectors:
        twice_jz = int(round(2 * float(sector["jz"])))
        key = twice_jz % modulus if modulus else twice_jz
        groups[(int(sector["n"]), key)].append(sector)
    result = []
    for (n, key), members in sorted(groups.items()):
        basis = [state for member in members for state in member["basis"]]
        ham = np.zeros((len(basis), len(basis)))
        labeled = all("qn_vectors" in member for member in members)
        qn = np.zeros_like(ham) if labeled else None
        labels = []
        start = 0
        for member in members:
            stop = start + len(member["basis"])
            ham[start:stop, start:stop] = member["ham"]
            if labeled:
                qn[start:stop, start:stop] = member["qn_vectors"]
                labels.extend(member["qn_labels"])
            start = stop
        index = {state: i for i, state in enumerate(basis)}
        for col, state in enumerate(basis):
            for dest, amplitude in op.one_body(state, off).items():
                if dest not in index:
                    raise ValueError("Off-diagonal impurity level escapes its symmetry block")
                ham[index[dest], col] += amplitude
        entry = {"n": n, "jz": key / 2, "basis": basis,
                 "basis_jz": [op.jz_of(state) for state in basis],
                 "ham": ham, "exact_components": [(len(result) + 1, list(range(len(basis))))]}
        if labeled:
            entry["qn_vectors"], entry["qn_labels"] = qn, labels
        result.append(entry)
    return result


def build_one_body_matrix(
    l: int,
    ml_ms: Sequence[Tuple[int, float]],
    one_electron_basis: Sequence[Tuple[float, float]],
    impqs: Sequence[int],
    eimp: Sequence[float],
    cg: ClebschGordan,
) -> np.ndarray:
    if not one_electron_basis:
        return np.zeros((0, 0), dtype=float)

    n_groups = (max(impqs) + 1) if impqs else 0
    eps = list(eimp)
    if len(eps) < n_groups:
        eps.extend([0.0] * (n_groups - len(eps)))

    eps_alpha = np.array([eps[impqs[ia]] for ia in range(len(one_electron_basis))], dtype=float)
    transform = build_jm_to_mlms_transform(l, ml_ms, one_electron_basis, cg)
    return transform.T @ np.diag(eps_alpha) @ transform


def build_coulomb_tensor_mlms(
    op: AtomicDirectBasis,
    gck: GauntCK,
    fkoj: Sequence[float],
    j_coulomb: float,
) -> np.ndarray:
    """Return the multiplet-only Coulomb tensor in the |m_l,m_s> basis.

    The Bold solver stores the local four-fermion interaction as a tensor over
    the chosen one-particle basis.  The atomic backend already knows the Slater
    multiplet structure in the |m_l,m_s> basis, so we first assemble the tensor
    there and later rotate it to the |j,m_j> basis used by the solver.

    The returned tensor contains only the J-scaled multiplet part.  The monopole
    F0 piece is intentionally omitted because the impurity solver adds it
    separately for the single-site case, matching the historical ctqmc/bold
    convention.
    """

    baths = op.baths
    uc_mlms = np.zeros((baths, baths, baths, baths), dtype=float)
    if abs(j_coulomb) < 1e-15:
        return uc_mlms

    for i1 in range(baths):
        m1 = op.lz[i1]
        s1 = op.sz[i1]
        for i2 in range(baths):
            m2 = op.lz[i2]
            s2 = op.sz[i2]
            for i3 in range(baths):
                if op.sz[i3] != s2:
                    continue
                m3 = op.lz[i3]
                for i4 in range(baths):
                    if op.sz[i4] != s1:
                        continue
                    m4 = op.lz[i4]
                    if m4 - m1 != m2 - m3:
                        continue
                    val = 0.0
                    for k in range(1, op.l + 1):
                        val += gck(m4, m1, 2 * k) * gck(m2, m3, 2 * k) * fkoj[k]
                    if abs(val) > 1e-14:
                        # The exact backend loops over two occupied orbitals
                        # (i1,i2), removes them, and then creates two empty
                        # orbitals (i3,i4).  The operator ordering represented
                        # by this loop is therefore:
                        #
                        #   c^+_{i4} c^+_{i3} c_{i2} c_{i1}
                        #
                        # rather than the solver convention
                        #
                        #   U[a,b,c,d] c^+_a c^+_b c_c c_d .
                        #
                        # We keep the backend tensor in its natural order here
                        # and permute it once, explicitly, right before writing
                        # the solver-facing UCoulomb section.
                        uc_mlms[i1, i2, i3, i4] = j_coulomb * val
    return uc_mlms


def rotate_coulomb_tensor_to_jm(uc_mlms: np.ndarray, jm_to_mlms: np.ndarray) -> np.ndarray:
    """Rotate U[p,q,r,s] from |m_l,m_s> orbitals to |j,m_j> orbitals.

    If c_a^\dagger creates the solver orbital |j,m_j>_a and d_p^\dagger creates
    the backend orbital |m_l,m_s>_p, then

      c_a^\dagger = sum_p U[a,p] d_p^\dagger

    with U[a,p] = <m_l,m_s|j,m_j>.  The four-index Coulomb tensor therefore
    transforms with one copy of the one-particle unitary on each leg.
    """

    return np.einsum(
        "ap,bq,cr,ds,pqrs->abcd",
        jm_to_mlms,
        jm_to_mlms,
        jm_to_mlms,
        jm_to_mlms,
        uc_mlms,
        optimize=True,
    )


def backend_to_solver_coulomb_order(uc_backend: np.ndarray) -> np.ndarray:
    """Convert backend U to solver convention U[a,b,c,d] c^+_a c^+_b c_c c_d.

    The backend tensor uses the natural order of the state-update loop:

      uc_backend[i1,i2,i3,i4]  <->  c^+_{i4} c^+_{i3} c_{i2} c_{i1}

    while the cix reader expects

      uc_solver[a,b,c,d]       <->  c^+_a c^+_b c_c c_d .

    Therefore the solver-facing tensor is a simple reversal of the four legs.
    """

    return np.transpose(uc_backend, (3, 2, 1, 0))


def match_direct_states_to_exact_eigenvalues(ham: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    eigvals, eigvecs = np.linalg.eigh(ham)
    overlaps = np.abs(eigvecs) ** 2
    n = ham.shape[0]
    direct_to_eig = np.full(n, -1, dtype=int)

    if linear_sum_assignment is not None:
        # Maximize the total direct-basis weight captured by the one-to-one
        # matching while preserving cheap 1x1 approximate sectors.
        rows, cols = linear_sum_assignment(-overlaps)
        direct_to_eig[rows] = cols
    else:
        pairs: List[Tuple[float, int, int]] = []
        for direct_idx in range(n):
            for eig_idx in range(n):
                pairs.append((float(overlaps[direct_idx, eig_idx]), direct_idx, eig_idx))
        pairs.sort(reverse=True)

        assigned_direct = set()
        assigned_eig = set()
        for overlap, direct_idx, eig_idx in pairs:
            if direct_idx in assigned_direct or eig_idx in assigned_eig:
                continue
            direct_to_eig[direct_idx] = eig_idx
            assigned_direct.add(direct_idx)
            assigned_eig.add(eig_idx)
            if len(assigned_direct) == n:
                break

        for direct_idx in range(n):
            if direct_to_eig[direct_idx] >= 0:
                continue
            best_eig = int(np.argmax(overlaps[direct_idx, :]))
            direct_to_eig[direct_idx] = best_eig

    matched_energies = np.array([eigvals[eig_idx] for eig_idx in direct_to_eig], dtype=float)
    matched_overlaps = np.array([overlaps[direct_idx, eig_idx] for direct_idx, eig_idx in enumerate(direct_to_eig)], dtype=float)
    return matched_energies, matched_overlaps


def compute_admissible_approx_components(
    approx_set: Sequence[int],
    sector_data: Sequence[Dict[str, object]],
    one_electron_basis: Sequence[Tuple[float, float]],
    op: AtomicDirectBasis,
    cg: ClebschGordan,
) -> List[List[int]]:
    approx_n = set(int(x) for x in approx_set)
    if not approx_n:
        return []
    if any("basis_jz" in sector for sector in sector_data):
        # Even with an approximate diagonal atomic Hamiltonian, off-diagonal
        # baths generate pseudo self-energy throughout the residual-symmetry
        # block. Do not split those blocks back into fixed-Jz pieces.
        return [list(sector["basis"]) for sector in sector_data if sector["n"] in approx_n]

    state_to_exact_block: Dict[int, int] = {}
    exact_fixed_blocks: List[Tuple[int, List[int]]] = []
    state_to_sector: Dict[int, Tuple[int, float]] = {}
    all_states: List[int] = []

    for sector in sector_data:
        basis = [int(state) for state in sector["basis"]]
        n = int(sector["n"])
        jz = float(sector["jz"])
        for state in basis:
            state_to_sector[state] = (n, jz)
            all_states.append(state)
        for block_id, comp in sector["exact_components"]:
            members = [basis[i] for i in comp]
            if n in approx_n:
                continue
            exact_fixed_blocks.append((int(block_id), members))
            for state in members:
                state_to_exact_block[state] = int(block_id)

    approx_states = [state for state in all_states if state_to_sector[state][0] in approx_n]
    if not approx_states:
        return []

    parent = {state: state for state in approx_states}
    rank = {state: 0 for state in approx_states}

    def find(state: int) -> int:
        while parent[state] != state:
            parent[state] = parent[parent[state]]
            state = parent[state]
        return state

    def union(a: int, b: int) -> bool:
        ra, rb = find(a), find(b)
        if ra == rb:
            return False
        if rank[ra] < rank[rb]:
            ra, rb = rb, ra
        parent[rb] = ra
        if rank[ra] == rank[rb]:
            rank[ra] += 1
        return True

    out_edges: Dict[int, Dict[int, int]] = defaultdict(dict)
    in_edges: Dict[int, Dict[int, List[int]]] = defaultdict(lambda: defaultdict(list))
    for state in all_states:
        for orb, (j, jz) in enumerate(one_electron_basis):
            for dst, amp in op.fp(state, j, jz, cg):
                if abs(amp) < 1e-12:
                    continue
                if dst in state_to_sector:
                    out_edges[state][orb] = int(dst)
                    in_edges[int(dst)][orb].append(int(state))
                break

    changed = True
    while changed:
        changed = False
        current_blocks: List[List[int]] = [members for _, members in exact_fixed_blocks]
        comps = defaultdict(list)
        for state in approx_states:
            comps[find(state)].append(state)
        current_blocks.extend(comps.values())

        for members in current_blocks:
            for orb in range(len(one_electron_basis)):
                dest_states: List[int] = []
                dest_fixed_blocks: set[int] = set()
                for state in members:
                    dst = out_edges[state].get(orb)
                    if dst is None:
                        continue
                    if state_to_sector[dst][0] in approx_n:
                        dest_states.append(dst)
                    else:
                        dest_fixed_blocks.add(state_to_exact_block[dst])
                if len(dest_fixed_blocks) > 1:
                    raise ValueError("Approximate partition is incompatible with fixed exact blocks under F^+ images")
                if dest_states:
                    base = dest_states[0]
                    for dst in dest_states[1:]:
                        if union(base, dst):
                            changed = True

        current_blocks = [members for _, members in exact_fixed_blocks]
        comps = defaultdict(list)
        for state in approx_states:
            comps[find(state)].append(state)
        current_blocks.extend(comps.values())

        for members in current_blocks:
            for orb in range(len(one_electron_basis)):
                src_states: List[int] = []
                src_fixed_blocks: set[int] = set()
                for dst in members:
                    for src in in_edges[dst].get(orb, []):
                        if state_to_sector[src][0] in approx_n:
                            src_states.append(src)
                        else:
                            src_fixed_blocks.add(state_to_exact_block[src])
                if len(src_fixed_blocks) > 1:
                    raise ValueError("Approximate partition is incompatible with fixed exact blocks under F^+ preimages")
                if src_states:
                    base = src_states[0]
                    for src in src_states[1:]:
                        if union(base, src):
                            changed = True

    comps = defaultdict(list)
    for state in approx_states:
        comps[find(state)].append(state)
    return list(comps.values())


def generate_bctqmc(
    l: int = 3,
    ns: int = 0,
    ne: int = 2,
    j_coulomb: float = 0.0,
    impq: Optional[Sequence[int]] = None,
    eimp: Optional[Sequence[float]] = None,
    approx_n: Optional[Sequence[int]] = None,
    approx_mode: str = "matched",
    output: str = "bctqmc.cix",
    sigind: Optional[Sequence[Sequence[int]]] = None,
) -> None:
    op_mlms = AtomicDirectBasis(l)
    cg = ClebschGordan()
    gck = GauntCK(l, cg)
    fkoj = compute_coulomb_ratios(l)
    one_electron_basis = build_one_electron_basis(l)
    op = JmProductBasis(one_electron_basis)
    impqs, epsk, h1_jm, baths = prepare_hybridization(one_electron_basis, impq, eimp, sigind)
    approx_set = set(approx_n or [])
    global_flip = compute_global_flip(one_electron_basis)
    # Assemble the Jz-conserving part using the established atomic backend.
    # Off-diagonal jm levels are added after merging the symmetry blocks.
    h1_mlms = build_one_body_matrix(l, op_mlms.ml_ms, one_electron_basis, impqs, epsk, cg)
    jm_to_mlms = build_jm_to_mlms_transform(l, op_mlms.ml_ms, one_electron_basis, cg)
    uc_mlms_backend = build_coulomb_tensor_mlms(op_mlms, gck, fkoj, j_coulomb)
    uc_jm_backend = rotate_coulomb_tensor_to_jm(uc_mlms_backend, jm_to_mlms)
    uc_jm = backend_to_solver_coulomb_order(uc_jm_backend)

    sector_data: List[Dict[str, object]] = []
    next_exact_block_id = 1
    for n in range(ns, ne + 1):
        for jz in build_jz_grid(n, l):
            basis_jm = construct_reduced_base_njz(op, n, jz)
            if not basis_jm:
                continue
            basis_mlms = construct_reduced_base_njz(op_mlms, n, jz)
            ham_mlms = build_hamiltonian_block(op_mlms, basis_mlms, j_coulomb, h1_mlms, gck, fkoj)
            transform = build_sector_transform(jm_to_mlms, basis_jm, basis_mlms, op.baths)
            ham = transform @ ham_mlms @ transform.T
            ham = 0.5 * (ham + ham.T)
            # In the |j,m_j> product basis the natural bold-CTQMC exact blocks are
            # the full (N,Jz) sectors. This keeps the solver basis aligned with the
            # conserved quantum number used in the bath/orbital labeling.
            exact_components: List[Tuple[int, List[int]]] = [(next_exact_block_id, list(range(len(basis_jm))))]
            next_exact_block_id += 1
            sector_entry: Dict[str, object] = {
                "n": n,
                "jz": jz,
                "basis": list(basis_jm),
                "ham": ham,
                "exact_components": exact_components,
            }
            if n not in approx_set:
                # The CIX Hamiltonian is in the |j,m_j> product basis. Build
                # NLSJ vectors in |m_l,m_s>, then rotate them into precisely
                # the same CIX basis and ordering. Do not use eigenvalue order
                # to associate labels with the solver's Hamiltonian eigenstates.
                qn_mlms, qn_labels = build_nlsj_basis(op_mlms, basis_mlms, n)
                sector_entry["qn_vectors"] = transform @ qn_mlms
                sector_entry["qn_labels"] = qn_labels
            sector_data.append(sector_entry)

    sector_data = merge_hybridization_sectors(sector_data, one_electron_basis, baths, h1_jm, op)
    for sector in sector_data:
        if sector["n"] in approx_set:
            if approx_mode == "matched":
                energies, overlaps = match_direct_states_to_exact_eigenvalues(sector["ham"])
            else:
                energies = np.diag(sector["ham"]).astype(float)
                overlaps = np.ones(len(sector["basis"]), dtype=float)
            sector["matched_energies"], sector["matched_overlaps"] = energies, overlaps

    approx_components = compute_admissible_approx_components(approx_set, sector_data, one_electron_basis, op, cg)
    approx_component_by_state: Dict[int, int] = {}
    for comp_id, members in enumerate(approx_components):
        for state in members:
            approx_component_by_state[int(state)] = comp_id

    sectors: List[Dict[str, object]] = []
    global_state_id = 1
    block_id = 1

    for sector in sector_data:
        n = int(sector["n"])
        jz = float(sector["jz"])
        basis = [int(state) for state in sector["basis"]]
        ham = np.asarray(sector["ham"], dtype=float)
        index = {state: i for i, state in enumerate(basis)}

        if n in approx_set:
            comp_members: Dict[int, List[int]] = defaultdict(list)
            for state in basis:
                comp_members[approx_component_by_state[state]].append(state)
            matched_energies = np.asarray(sector["matched_energies"], dtype=float)
            matched_overlaps = np.asarray(sector["matched_overlaps"], dtype=float)
            for _, members in sorted(comp_members.items(), key=lambda item: min(index[s] for s in item[1])):
                comp_indices = [index[state] for state in members]
                sub_basis = [basis[i] for i in comp_indices]
                diag_values = matched_energies[comp_indices]
                sub_ham = np.diag(diag_values)
                eigvals = np.linalg.eigvalsh(sub_ham)
                state_ids = list(range(global_state_id, global_state_id + len(sub_basis)))
                overlap_map = {sub_basis[i]: float(matched_overlaps[comp_indices[i]]) for i in range(len(sub_basis))}
                sectors.append(
                    {
                        "block_id": block_id,
                        "n": n,
                        "jz": jz,
                        "basis": sub_basis,
                        "ham": sub_ham,
                        "eigvals": eigvals,
                        "state_ids": state_ids,
                        "approx_mode": approx_mode,
                        "match_overlaps": overlap_map,
                        "qn_vectors": None,
                        **({"basis_jz": [op.jz_of(state) for state in sub_basis]}
                           if "basis_jz" in sector else {}),
                    }
                )
                global_state_id += len(sub_basis)
                block_id += 1
            continue

        for _, comp in sector["exact_components"]:
            sub_basis = [basis[i] for i in comp]
            sub_ham = ham[np.ix_(comp, comp)]
            eigvals = np.linalg.eigvalsh(sub_ham)
            state_ids = list(range(global_state_id, global_state_id + len(sub_basis)))
            sectors.append(
                {
                    "block_id": block_id,
                    "n": n,
                    "jz": jz,
                    "basis": sub_basis,
                    "ham": sub_ham,
                    "eigvals": eigvals,
                    "state_ids": state_ids,
                    "approx_mode": None,
                    # A strict sub-block generally is not invariant under
                    # L²/S²/J²; never attach a sliced, incorrectly labeled
                    # subset of the full-sector angular-momentum basis.
                    "qn_vectors": sector["qn_vectors"] if len(comp) == len(basis) else None,
                    "qn_labels": sector["qn_labels"] if len(comp) == len(basis) else None,
                    **({"basis_jz": [op.jz_of(state) for state in sub_basis]}
                       if "basis_jz" in sector else {}),
                }
            )
            global_state_id += len(sub_basis)
            block_id += 1

    state_lookup: Dict[int, Tuple[int, int]] = {}
    for sector in sectors:
        basis = sector["basis"]
        state_ids = sector["state_ids"]
        for local_idx, state in enumerate(basis):
            state_lookup[int(state)] = (int(sector["block_id"]), int(state_ids[local_idx]))

    max_block = max((len(sector["basis"]) for sector in sectors), default=0)
    with open(output, "w", encoding="utf-8") as fh:
        print("# Cix file for atomic Bold-CTQMC", file=fh)
        print(
            "# cluster_size, number of states, number of psi^+ operators, "
            "number of baths, number of blocks, maximum_matrix_size",
            file=fh,
        )
        print(
            1,
            global_state_id - 1,
            len(one_electron_basis),
            len(baths),
            len(sectors),
            max_block,
            file=fh,
        )
        print("# baths, psi^+, psi, Sigind, # spin, moment", file=fh)
        for ib, (a, b, component) in enumerate(baths):
            j, jz = one_electron_basis[a]
            if a == b:
                print(f"{ib} {a} {b}* {component} {jz:4.1f} ({j:3.1f},{jz:4.1f})", file=fh)
            else:
                jb, mjb = one_electron_basis[b]
                print(f"{ib} {a} {b}* {component} 0.0 "
                      f"({j:3.1f},{jz:4.1f})->({jb:3.1f},{mjb:4.1f})", file=fh)
        print("# cluster energies for non-equivalent baths, eps[k]", file=fh)
        print(" ".join(f"{x:.10f}" for x in epsk), file=fh)
        print(f"#  iblck  n k  Jz   {' '.join(f'F^+_{i:<2}' for i in range(len(one_electron_basis)))}", file=fh)

        for sector in sectors:
            basis = sector["basis"]
            state_ids = sector["state_ids"]
            for local_idx, state in enumerate(basis):
                state_id = state_ids[local_idx]
                row_jz = op.jz_of(state) if "basis_jz" in sector else sector["jz"]
                print(f"{state_id:<4d} {sector['block_id']:3d} {sector['n']:2d} {0:1d} {row_jz:5.1f}  ", end="", file=fh)
                for orb, (j, jz) in enumerate(one_electron_basis):
                    dest = 0
                    for new_state, amp in op.fp(state, j, jz, cg):
                        if abs(amp) < 1e-12:
                            continue
                        found = state_lookup.get(new_state)
                        if found is not None:
                            _, dest_id = found
                            sign = 1 if amp > 0 else -1
                            dest = sign * dest_id
                            break
                    print(f"{dest:6d} ", end="", file=fh)
                comment = op.printn(state)
                if sector.get("approx_mode") == "matched":
                    overlap = sector["match_overlaps"][int(state)]
                    comment += f" approx=matched overlap={overlap:.6f}"
                elif sector.get("approx_mode") == "diagonal":
                    comment += " approx=diagonal"
                print(f"# {comment}", file=fh)

        print("# All blocks of local Hamiltonian follow", file=fh)
        for sector in sectors:
            state_ids = sector["state_ids"]
            eigvals = sector["eigvals"]
            ham = sector["ham"]
            print(
                sector["block_id"],
                len(state_ids),
                "   ",
                " ".join(str(x) for x in state_ids),
                "  ",
                [float(x) for x in eigvals],
                "# iblock, dim, indices, eigvals",
                file=fh,
            )
            for row in ham:
                print("   " + " ".join(f"{float(x):14.10f}" for x in row), file=fh)

        # The reader skips exactly one line here and then expects raw
        # UCoulomb tuples.  Keep the whole description on a single line.
        print("# UCoulomb : (j1,mj1) (j2,mj2) (j3,mj3) (j4,mj4)  Uc[a,b,c,d]", file=fh)
        for a, (j1, mj1) in enumerate(one_electron_basis):
            for b, (j2, mj2) in enumerate(one_electron_basis):
                for c, (j3, mj3) in enumerate(one_electron_basis):
                    for d, (j4, mj4) in enumerate(one_electron_basis):
                        val = float(uc_jm[a, b, c, d])
                        if abs(val) <= 1e-10:
                            continue
                        print(
                            f"{a:2d} {b:2d} {c:2d} {d:2d}  {val: .12f}"
                            f"  # ({j1:3.1f},{mj1:4.1f}) ({j2:3.1f},{mj2:4.1f})"
                            f" ({j3:3.1f},{mj3:4.1f}) ({j4:3.1f},{mj4:4.1f})",
                            file=fh,
                        )

    # A separate file keeps the long-established CIX reader and its column
    # layout unchanged. The solver treats this file as optional, but verifies
    # the CIX fingerprint before it attaches any angular-momentum labels.
    write_nlsj_sidecar(output, sectors, global_state_id - 1, one_electron_basis)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate a direct-basis bctqmc.cix for an isolated atom. "
        "Spin-orbit is assumed to be already included in Eimp."
    )
    parser.add_argument("-l", type=int, default=3, help="Orbital angular momentum of the shell")
    parser.add_argument("-ns", type=int, default=0, help="Lowest occupancy")
    parser.add_argument("-ne", type=int, default=2, help="Highest occupancy")
    parser.add_argument("-J", dest="j_coulomb", type=float, default=0.0, help="Coulomb scale J")
    parser.add_argument("-Impq", type=str, default="", help="Equivalent bath indices, e.g. [0,0,1,1,...]")
    parser.add_argument("-Eimp", type=str, default="", help="Bath energies, e.g. [0.0,0.2]")
    parser.add_argument("-Sigind", "--sigind", default="", help="Full one-based hybridization index matrix; zero means absent")
    parser.add_argument("--trans", default=None, metavar="FILE",
                        help="Read Sigind from FILE; by default read Trans.dat when present")
    parser.add_argument("-a", "--approx", type=str, default="", help="Valences to treat approximately, e.g. [4,5]")
    parser.add_argument("--approx-mode", choices=("diagonal", "matched"), default="matched", help="Approximation mode for valences listed in --approx")
    parser.add_argument("-o", "--output", type=str, default="bctqmc.cix", help="Output cix filename")
    parser.add_argument("--label-existing", type=str, default="", metavar="CIX",
                        help="Write only CIX.nlsj for an existing atomh.py-style CIX; leave CIX untouched")
    args = parser.parse_args()

    if args.label_existing:
        label_existing_cix(args.label_existing, args.l)
        return

    impq = parse_list_arg(args.Impq, int) if args.Impq else None
    eimp = parse_list_arg(args.Eimp, float) if args.Eimp else None
    approx_n = parse_list_arg(args.approx, int) if args.approx else None
    if args.sigind and args.trans:
        parser.error("Use either --sigind or --trans")
    sigind = ast.literal_eval(args.sigind) if args.sigind else None
    trans = Path(args.trans) if args.trans else Path("Trans.dat")
    if sigind is None and (args.trans or trans.is_file()):
        sigind = read_trans_sigind(str(trans))
    generate_bctqmc(
        l=args.l,
        ns=args.ns,
        ne=args.ne,
        j_coulomb=args.j_coulomb,
        impq=impq,
        eimp=eimp,
        approx_n=approx_n,
        approx_mode=args.approx_mode,
        output=args.output,
        sigind=sigind,
    )


if __name__ == "__main__":
    main()
