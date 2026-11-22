#!/usr/bin/env python3
"""Contract a directional derivative exported by native LAPWSO.

LAPWSO reads EDMFT_SOC_DISPLACEMENT (layout and one Cartesian velocity
per physical site for rigid translation, or scalar tangents prepared by
soc_scalar_window.py for the cross-window term). Use a diagonal layout for DFT or a matrix layout
with EDMFT_ONEPARTICLE_DUMP from a nonmagnetic spinor DMFT2 force pass.
Files from each pass must use the same k mesh and the same eigenvectors.
The result is one explicit operator contribution in mRy/Bohr per pattern.
"""

import argparse
from pathlib import Path
import numpy as np
from scipy.special import expit


def read(stream, dtype, count):
    a = np.fromfile(stream, dtype, count)
    if a.size != count or not np.all(np.isfinite(a)):
        raise ValueError(f"truncated or nonfinite array in {stream.name}")
    return a


def _gauge(coeff):
    bits = np.ascontiguousarray(coeff).view(np.uint64).reshape(*coeff.shape, 2)
    value = np.zeros(coeff.shape[1:], np.uint64)
    for real, imag in zip(bits[..., 0], bits[..., 1]):
        value = ((value << 1) | (value >> 63)) ^ real
        value = ((value << 7) | (value >> 57)) ^ imag
    return value


def _dmft_term(path, ik, k, energies, matrix, gauge):
    with path.open("rb") as f:
        if f.read(8) != b"EDMFTRH1" or read(f, "i4", 1)[0] != 0x01020304:
            raise ValueError(f"invalid one-particle export {path}")
        version, index, _, size, _, n, iso, first, wfirst, last, _, _ = read(
            f, "i4", 12
        )
        if (
            (version, index, iso) != (1, ik, 2)
            or not (1 <= first <= wfirst <= last <= len(energies))
            or n != last - first + 1
        ):
            raise ValueError(f"incompatible spinor density dimensions in {path}")
        sl = slice(first - 1, last)
        if not np.allclose(k, read(f, "f8", 3), atol=1e-12, rtol=0):
            raise ValueError(f"k coordinates differ in {path}")
        read(f, "i4", 3 * size)
        if not np.allclose(energies[sl], read(f, "f8", n), atol=1e-11, rtol=0):
            raise ValueError(f"band energies differ in {path}")
        coeff = read(f, "c16", size * n * 2).reshape((size, n, 2), order="F")
        if not np.array_equal(gauge[sl], _gauge(coeff)):
            raise ValueError(f"spinor eigenvector gauges differ in {path}")
        rho = read(f, "c16", n * n).reshape((n, n), order="F")
        if not np.allclose(rho, rho.T.conj(), atol=1e-10, rtol=0):
            raise ValueError(f"non-Hermitian spinor density in {path}")
        if len(f.read(n * n * 16)) != n * n * 16 or f.read(1):
            raise ValueError(f"incorrect energy moment size in {path}")
        # At iso=2, nspin1=1, DMFT2 stores 2*w_k*rho and uses vnorm1=1/2
        # in physical traces. The two components are one coherent spinor.
        return 0.5 * np.vdot(rho, matrix[sl, sl]).real


def force(prefix, klist, *, mu=None, temperature=None, moments=None):
    dft = moments is None
    if dft and (
        mu is None
        or temperature is None
        or not np.isfinite(mu)
        or not 0 < temperature < np.inf
    ):
        raise ValueError(
            "DFT contraction requires a finite chemical potential and positive temperature in Ry"
        )
    from lapw_response_kpool import klist_records

    rows, weight_end, _ = klist_records(klist)
    weights = np.array([float(row[weight_end - 5 : weight_end]) for row in rows])
    if not len(weights) or not np.all(np.isfinite(weights)) or min(weights) <= 0:
        raise ValueError("invalid k mesh")
    total_weight = weights.sum()
    terms = []
    for ik, w in enumerate(weights, 1):
        p = Path(str(prefix) + f".k{ik}")
        with p.open("rb") as f:
            version, index, nsv, n, full = read(f, "i4", 5)
            if version != 2 or index != ik or not 0 < n <= nsv or full not in (0, 1):
                raise ValueError(f"invalid LAPWSO derivative header {p}")
            k = read(f, "f8", 3)
            if not np.isclose(w, read(f, "f8", 1)[0], rtol=0, atol=1e-9):
                raise ValueError(f"k weight differs in {p}")
            energies = read(f, "f8", n)
            diagonal = read(f, "f8", n)
            if full:
                matrix = read(f, "c16", n * n).reshape((n, n), order="F")
                gauge = read(f, "u8", n * 2).reshape((n, 2), order="F")
                if not np.allclose(
                    matrix, matrix.T.conj(), atol=1e-11, rtol=0
                ) or not np.allclose(diagonal, matrix.diagonal(), atol=1e-11, rtol=0):
                    raise ValueError(f"inconsistent SOC derivative in {p}")
            if f.read(1):
                raise ValueError(f"trailing data in {p}")
        if dft:
            terms.append(
                w / total_weight * (expit((mu - energies) / temperature) @ diagonal)
            )
        else:
            if not full:
                raise ValueError("eDMFT needs the full matrix in spinor band space")
            terms.append(
                _dmft_term(
                    Path(str(moments) + f".k{ik:08d}.bin"),
                    ik,
                    k,
                    energies,
                    matrix,
                    gauge,
                )
            )
    return -1000 * sum(terms)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("prefix", type=Path)
    parser.add_argument("--klist", type=Path, required=True)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--moments", type=Path)
    mode.add_argument("--mu", type=float, help="DFT Fermi energy, Ry")
    parser.add_argument("--temperature", type=float, help="DFT kBT, Ry")
    args = parser.parse_args()
    result = force(
        args.prefix,
        args.klist,
        mu=args.mu,
        temperature=args.temperature,
        moments=args.moments,
    )
    layout = args.prefix.read_text().splitlines()[0].strip()
    term = (
        "Scalar SOC subspace motion"
        if layout.startswith("subspace_")
        else "Muffin-tin SOC translation at fixed scalar coefficients"
    )
    print(f"{term}: {result:+.9f} mRy/Bohr")
