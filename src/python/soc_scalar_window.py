#!/usr/bin/env python3
"""Compute the part of the scalar LAPW eigenvector derivative outside the SOC input space.

The input vector is the scalar vector passed to LAPWSO. All vectors, matrices and
directional derivatives must refer to the same ordered LAPW coefficient frame.
H is in Ry, O is dimensionless, and the displacement parameter is in Bohr.
Pass either both derivative matrices directly (--derivative prefix, in the
native H/O dump format) or two snapshots in that frame (--plus, --minus, --step).
"""

import argparse
from pathlib import Path
import struct
import numpy as np
from scipy.linalg import cholesky, eigh, qr, solve_triangular

from lapw_dft_vector import _record, _header
from lapw_native_snapshot import read_basis_matrix, vector_k_from_matrix_header
from soc_translation_force import _gauge
from wstruct import Struct


def outside_tangent(H, O, energies, A, dH, dO):
    """Return dA=Q_out dA in the O metric; omit in-window gauge/metric terms."""
    OA = O @ A
    if np.max(np.abs(A.conj().T @ OA - np.eye(len(energies)))) > 1e-7:
        raise ValueError("scalar vector is not orthonormal in the supplied overlap")
    if np.max(np.abs(H @ A - OA * energies)) > 1e-7:
        raise ValueError("scalar vector does not diagonalize the supplied Hamiltonian")
    if A.shape[0] == A.shape[1]:
        return np.zeros_like(A), np.inf
    # Complete the *actual* stored vectors, including when low-lying states
    # were omitted. Internal rotations of either subspace cancel at fixed H,O.
    L = cholesky(O, lower=True)
    complement = qr(L.conj().T @ A, mode="full")[0][:, len(energies) :]
    complement = solve_triangular(L.conj().T, complement, lower=False)
    outside, rotation = eigh(complement.conj().T @ H @ complement)
    B = complement @ rotation
    gap = energies[None, :] - outside[:, None]
    min_gap = np.min(np.abs(gap))
    if min_gap <= 1e-10:
        raise ValueError(f"degenerate scalar SOC cutoff (minimum gap {min_gap:.8g} Ry)")
    coupling = B.conj().T @ (dH @ A - (dO @ A) * energies)
    return B @ (coupling / gap), min_gap


def scalar_vectors(vector, nat):
    with Path(vector).open("rb") as stream:
        for _ in range(2 * nat):
            if _record(stream, skip=True) is None:
                raise ValueError("incomplete scalar radial header")
        while (header := _header(_record(stream))) is not None:
            *k, _, n, bands, _ = header
            labels = np.frombuffer(_record(stream), "<i4").reshape((3, n), order="F")
            energies = np.empty(bands)
            A = np.empty((n, bands), complex, order="F")
            for i in range(bands):
                number, energies[i] = struct.unpack("<id", _record(stream))
                data = _record(stream)
                if number != i + 1 or len(data) not in (8 * n, 16 * n):
                    raise ValueError("invalid scalar vector band record")
                A[:, i] = np.frombuffer(data, "<f8" if len(data) == 8 * n else "<c16")
            yield np.asarray(k), labels, energies, A


def matrices(prefix, ik, k, labels, s, expected_nlo=None):
    pair = []
    for tag, magic in ((".ham", b"EDMFTHV1"), ("", b"EDMFTOV1")):
        m = read_basis_matrix(f"{prefix}{tag}.k{ik:08d}.bin", magic)
        if expected_nlo is None:
            expected_nlo = m.nlo
        if (
            m.kpoint != ik
            or m.nlo != expected_nlo
            or not np.array_equal(m.labels, labels)
            or not np.allclose(
                vector_k_from_matrix_header(m.reciprocal_k, s), k, atol=1e-11, rtol=0
            )
        ):
            raise ValueError(
                f"k, reciprocal labels or local orbitals differ at {prefix}, k={ik}"
            )
        pair.append(m.matrix)
    return pair[0], pair[1], expected_nlo


def prepare(args):
    from math import isfinite

    s = Struct(inAngs=False)
    s.ReadStruct(str(args.structure))
    finite = args.plus is not None
    if finite != (args.minus is not None and args.step is not None) or (
        not finite and (args.minus is not None or args.step is not None)
    ):
        raise ValueError(
            "two signed snapshots and a positive step are required for differences"
        )
    if finite and not (isfinite(args.step) and args.step > 0):
        raise ValueError("invalid displacement step")
    # Fresh directory: an incomplete calculation has no launchable control
    # file, and output cannot overwrite another displacement's derivatives.
    args.output_dir.mkdir()
    prefix = args.output_dir / "subspace"
    count = 0
    for ik, (k, labels, energies, A) in enumerate(
        scalar_vectors(args.vector, s.nat), 1
    ):
        H, O, nlo = matrices(args.central, ik, k, labels, s)
        if finite:
            dH, dO, _ = matrices(args.plus, ik, k, labels, s, nlo)
            hm, om, _ = matrices(args.minus, ik, k, labels, s, nlo)
            dH -= hm
            dO -= om
            dH /= 2 * args.step
            dO /= 2 * args.step
        else:
            dH, dO, _ = matrices(args.derivative, ik, k, labels, s, nlo)
        tangent, gap = outside_tangent(H, O, energies, A, dH, dO)
        if not np.all(np.isfinite(tangent)):
            raise ValueError(f"nonfinite scalar tangent at k={ik}")
        with Path(f"{prefix}.dA.k{ik}").open("xb") as f:
            f.write(b"SOCWIN01")
            np.array([ik, *A.shape], "<i4").tofile(f)
            k.astype("<f8").tofile(f)
            energies.astype("<f8").tofile(f)
            _gauge(A[:, :, None]).astype("<u8").tofile(f)
            tangent.T.astype("<c16").tofile(f)
        print(
            f"k={ik}  retained={A.shape[1]}  excluded={A.shape[0] - A.shape[1]}  min_gap_Ry={gap:.9g}",
            flush=True,
        )
        count = ik
    if not count:
        raise ValueError("empty scalar vector")
    prefix.write_text("subspace_" + args.layout + "\n")
    print(f"EDMFT_SOC_DISPLACEMENT={prefix.resolve()}")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--structure", type=Path, required=True)
    p.add_argument("--vector", type=Path, required=True)
    p.add_argument(
        "--central", type=Path, required=True, help="central native H and O prefix"
    )
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--derivative", type=Path, help="native dH and dO prefix")
    g.add_argument(
        "--plus", type=Path, help="positive displacement native H and O prefix"
    )
    p.add_argument("--minus", type=Path)
    p.add_argument("--step", type=float, help="positive displacement size in Bohr")
    p.add_argument("--layout", choices=("diagonal", "matrix"), default="matrix")
    p.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="new directory for the tangent files",
    )
    prepare(p.parse_args())
