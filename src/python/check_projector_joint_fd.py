#!/usr/bin/env python3
"""Joint finite difference of the complete normalized DMFT projector.

This checker evaluates the projector part of the stationary Klein derivative

    2 Re sum_k Tr[C_k U_{k,strain}]

without splitting ``U_,strain`` into explicit, complementary-window, and
normalization response tensors.  Endpoint band gauges are parallel transported
to the central LAPW eigenvectors with the central generalized metric.  Exactly
degenerate central bands are transported as one unitary block.

Required default-off dmft2 diagnostics:

* ``EDMFT_PROJECTOR_SNAPSHOT=WRITE`` at minus and plus, which writes the
  complete normalized projector;
* ``EDMFT_ONEPARTICLE_DUMP=PREFIX`` at minus, center, and plus, which supplies
  LAPW coefficients and basis labels;
* ``EDMFT_STRESS_PROJECTOR_DIAG=1`` at the center, which writes the conjugate
  projector gradient;
* ``DFTKS_OVERLAP_DUMP=PREFIX`` at the center, which supplies the central
  coefficient-space metric.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
import sys

import numpy as np

from check_oneparticle_stress_fd import (
    OneParticleDump,
    aligned,
    files,
    matrix_path,
    read_basis_matrix,
    read_one_particle,
)


SNAPSHOT_MAGIC = 20261104
SNAPSHOT_VERSION = 2
GRADIENT_MAGIC = 20261105
GRADIENT_VERSION = 1


@dataclass(frozen=True)
class ProjectorSnapshot:
    kpoint: int
    spin: int
    nbands: int
    maxdim: int
    norbitals: int
    nipc: int
    energies: np.ndarray
    projector: np.ndarray


@dataclass(frozen=True)
class ProjectorGradient:
    kpoint: int
    nbands: int
    maxdim: int
    norbitals: int
    energies: np.ndarray
    gradient: np.ndarray


def exact(handle, count: int, dtype: str, what: str) -> np.ndarray:
    data = np.fromfile(handle, dtype=dtype, count=count)
    if data.size != count:
        raise ValueError(f"truncated {what} in {handle.name}")
    return data


def read_snapshot(path: Path) -> ProjectorSnapshot:
    with path.open("rb") as handle:
        header = exact(handle, 8, "=i4", "snapshot header").astype(int)
        magic, version, kpoint, spin, nbands, maxdim, norbitals, nipc = header
        if (magic, version) != (SNAPSHOT_MAGIC, SNAPSHOT_VERSION):
            raise ValueError(f"unsupported projector snapshot {header.tolist()} in {path}")
        energies = exact(handle, nbands, "=f8", "snapshot energies")
        projector = exact(
            handle,
            nbands * maxdim * norbitals * nipc,
            "=c16",
            "normalized projector",
        ).reshape((nbands, maxdim, norbitals, nipc), order="F")
        if handle.read(1):
            raise ValueError(f"trailing data in {path}")
    return ProjectorSnapshot(
        kpoint, spin, nbands, maxdim, norbitals, nipc, energies, projector
    )


def read_gradient(path: Path) -> ProjectorGradient:
    with path.open("rb") as handle:
        header = exact(handle, 6, "=i4", "gradient header").astype(int)
        magic, version, kpoint, nbands, maxdim, norbitals = header
        if (magic, version) != (GRADIENT_MAGIC, GRADIENT_VERSION):
            raise ValueError(f"unsupported projector gradient {header.tolist()} in {path}")
        energies = exact(handle, nbands, "=f8", "gradient energies")
        gradient = exact(
            handle,
            nbands * maxdim * norbitals,
            "=c16",
            "projector gradient",
        ).reshape((nbands, maxdim, norbitals), order="F")
        if handle.read(1):
            raise ValueError(f"trailing data in {path}")
    return ProjectorGradient(kpoint, nbands, maxdim, norbitals, energies, gradient)


def numbered_files(directory: Path, pattern: str, reader) -> dict[int, object]:
    result: dict[int, object] = {}
    for path in sorted(directory.glob(pattern)):
        value = reader(path)
        if value.kpoint in result:
            raise ValueError(f"duplicate k point {value.kpoint} in {directory}")
        result[value.kpoint] = value
    if not result:
        raise FileNotFoundError(f"no {pattern} files in {directory}")
    return result


def reordered_coefficients(
    reference: OneParticleDump, endpoint: OneParticleDump, role: str
) -> np.ndarray:
    if (endpoint.nbasis, endpoint.nnlo, endpoint.ncomponents) != (
        reference.nbasis,
        reference.nnlo,
        reference.ncomponents,
    ):
        raise ValueError(f"coefficient dimensions differ at k={reference.kpoint} for {role}")
    npw = reference.nbasis - reference.nnlo
    if not np.array_equal(endpoint.labels[:, npw:], reference.labels[:, npw:]):
        raise ValueError(f"local-orbital order differs at k={reference.kpoint} for {role}")
    lookup = {tuple(endpoint.labels[:, row]): row for row in range(npw)}
    wanted = [tuple(reference.labels[:, row]) for row in range(npw)]
    if len(lookup) != npw or set(lookup) != set(wanted):
        raise ValueError(f"fixed-G set differs at k={reference.kpoint} for {role}")
    permutation = [lookup[label] for label in wanted]
    permutation.extend(range(npw, reference.nbasis))
    first = endpoint.dm_nemin - endpoint.nemin
    last = endpoint.dm_nemax - endpoint.nemin + 1
    return endpoint.coefficients[np.asarray(permutation), first:last, :]


def window_coefficients(dump: OneParticleDump) -> np.ndarray:
    first = dump.dm_nemin - dump.nemin
    last = dump.dm_nemax - dump.nemin + 1
    return dump.coefficients[:, first:last, :]


def metric_overlap(left: np.ndarray, metric: np.ndarray, right: np.ndarray) -> np.ndarray:
    result = np.zeros((left.shape[1], right.shape[1]), complex)
    for component in range(left.shape[2]):
        result += left[:, :, component].conj().T @ metric @ right[:, :, component]
    return result


def blocks(energies: np.ndarray, tolerance: float) -> list[slice]:
    result: list[slice] = []
    first = 0
    for index in range(len(energies) - 1):
        if abs(energies[index + 1] - energies[index]) > tolerance:
            result.append(slice(first, index + 1))
            first = index + 1
    result.append(slice(first, len(energies)))
    return result


def transport(
    central_coefficients: np.ndarray,
    endpoint_coefficients: np.ndarray,
    central_metric: np.ndarray,
    central_energies: np.ndarray,
    endpoint_projector: np.ndarray,
    tolerance: float,
    central_projector: np.ndarray | None = None,
) -> tuple[np.ndarray, float, int]:
    if endpoint_projector.shape[0] != central_coefficients.shape[1]:
        raise ValueError("projector and coefficient band counts differ")
    transported = np.empty_like(endpoint_projector)
    residual_squared = 0.0
    largest_block = 0
    for block in blocks(central_energies, tolerance):
        central_block = central_coefficients[:, block, :]
        endpoint_block = endpoint_coefficients[:, block, :]
        cross = metric_overlap(central_block, central_metric, endpoint_block)
        left, _singular, right_h = np.linalg.svd(cross, full_matrices=True)
        rotation = right_h.conj().T @ left.conj().T
        aligned_coefficients = np.einsum(
            "bic,ij->bjc", endpoint_block, rotation, optimize=True
        )
        difference = aligned_coefficients - central_block
        residual_squared += float(
            np.real(np.trace(metric_overlap(difference, central_metric, difference)))
        )
        # DMFTU is stored with band index first as the transpose of
        # <chi|psi>.  If band columns transform as A -> A R, the stored
        # projector therefore transforms as DMFTU -> R^T DMFTU.
        if central_projector is None:
            projector_rotation = rotation.T
        else:
            central_flat = central_projector[block, ...].reshape(
                block.stop - block.start, -1
            )
            endpoint_flat = endpoint_projector[block, ...].reshape(
                block.stop - block.start, -1
            )
            projector_cross = central_flat @ endpoint_flat.conj().T
            pleft, _psingular, pright_h = np.linalg.svd(
                projector_cross, full_matrices=True
            )
            projector_rotation = pleft @ pright_h
        transported[block, ...] = np.einsum(
            "ij,j...->i...", projector_rotation, endpoint_projector[block, ...],
            optimize=True,
        )
        largest_block = max(largest_block, block.stop - block.start)
    return transported, np.sqrt(max(residual_squared, 0.0)), largest_block


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--minus-snapshots", required=True, type=Path)
    parser.add_argument("--center-snapshots", required=True, type=Path)
    parser.add_argument("--plus-snapshots", required=True, type=Path)
    parser.add_argument("--gradient", required=True, type=Path)
    parser.add_argument("--minus-coefficients", required=True, type=Path)
    parser.add_argument("--center-coefficients", required=True, type=Path)
    parser.add_argument("--plus-coefficients", required=True, type=Path)
    parser.add_argument("--center-overlap", required=True, type=Path)
    parser.add_argument("--strain", required=True, type=float)
    parser.add_argument("--degeneracy-tol", type=float, default=1.0e-7)
    parser.add_argument(
        "--transport",
        choices=("coefficients", "euclidean", "projector"),
        default="coefficients",
        help="parallel-transport gauge; projector is useful when cross-basis overlaps are unavailable",
    )
    parser.add_argument("--volume", type=float)
    args = parser.parse_args(argv)
    if args.strain <= 0 or not np.isfinite(args.strain):
        parser.error("--strain must be positive")
    if args.degeneracy_tol <= 0 or not np.isfinite(args.degeneracy_tol):
        parser.error("--degeneracy-tol must be positive")
    if args.volume is not None and (args.volume <= 0 or not np.isfinite(args.volume)):
        parser.error("--volume must be positive")
    return args


def main(argv: list[str] | None = None) -> int:
    args = parse_args(sys.argv[1:] if argv is None else argv)
    minus_snapshots = numbered_files(
        args.minus_snapshots, "dmft_projector_snapshot_MAIN_s01_k*.bin", read_snapshot
    )
    center_snapshots = numbered_files(
        args.center_snapshots, "dmft_projector_snapshot_MAIN_s01_k*.bin", read_snapshot
    )
    plus_snapshots = numbered_files(
        args.plus_snapshots, "dmft_projector_snapshot_MAIN_s01_k*.bin", read_snapshot
    )
    gradients = numbered_files(args.gradient, "dmft_projector_gradient_k*.bin", read_gradient)
    minus_coefficients = files(args.minus_coefficients)
    center_coefficients = files(args.center_coefficients)
    plus_coefficients = files(args.plus_coefficients)

    keys = set(gradients)
    for name, values in (
        ("minus snapshots", minus_snapshots),
        ("center snapshots", center_snapshots),
        ("plus snapshots", plus_snapshots),
        ("minus coefficients", minus_coefficients),
        ("center coefficients", center_coefficients),
        ("plus coefficients", plus_coefficients),
    ):
        if set(values) != keys:
            raise ValueError(f"k-point set of {name} differs from the gradient")

    central = 0.0
    minus_side = 0.0
    plus_side = 0.0
    coefficient_residual = {"minus": 0.0, "plus": 0.0}
    projector_change = {"minus": 0.0, "plus": 0.0}
    largest_block = 0

    for kpoint in sorted(keys):
        center_dump = read_one_particle(center_coefficients[kpoint])
        minus_dump = read_one_particle(minus_coefficients[kpoint])
        plus_dump = read_one_particle(plus_coefficients[kpoint])
        overlap_dump = read_basis_matrix(
            matrix_path(args.center_overlap, kpoint, False), b"EDMFTOV1"
        )
        metric = aligned(center_dump, overlap_dump, "central O")
        transport_metric = (
            np.eye(metric.shape[0], dtype=complex)
            if args.transport == "euclidean"
            else metric
        )
        central_coeff = window_coefficients(center_dump)
        minus_coeff = reordered_coefficients(center_dump, minus_dump, "minus")
        plus_coeff = reordered_coefficients(center_dump, plus_dump, "plus")
        gradient = gradients[kpoint]
        minus_snapshot = minus_snapshots[kpoint]
        center_snapshot = center_snapshots[kpoint]
        plus_snapshot = plus_snapshots[kpoint]
        expected = (
            central_coeff.shape[1],
            gradient.maxdim,
            gradient.norbitals,
        )
        if gradient.gradient.shape != expected:
            raise ValueError(f"gradient dimensions differ at k={kpoint}")
        for role, snapshot in (
            ("minus", minus_snapshot),
            ("center", center_snapshot),
            ("plus", plus_snapshot),
        ):
            if snapshot.projector.shape[:3] != expected or snapshot.nipc < 1:
                raise ValueError(f"{role} projector dimensions differ at k={kpoint}")

        minus_u, minus_residual, minus_block = transport(
            central_coeff,
            minus_coeff,
            transport_metric,
            gradient.energies,
            minus_snapshot.projector[:, :, :, 0],
            args.degeneracy_tol,
            center_snapshot.projector[:, :, :, 0]
            if args.transport == "projector"
            else None,
        )
        plus_u, plus_residual, plus_block = transport(
            central_coeff,
            plus_coeff,
            transport_metric,
            gradient.energies,
            plus_snapshot.projector[:, :, :, 0],
            args.degeneracy_tol,
            center_snapshot.projector[:, :, :, 0]
            if args.transport == "projector"
            else None,
        )
        coefficient_residual["minus"] = max(
            coefficient_residual["minus"], minus_residual
        )
        coefficient_residual["plus"] = max(
            coefficient_residual["plus"], plus_residual
        )
        largest_block = max(largest_block, minus_block, plus_block)

        derivative = (plus_u - minus_u) / (2 * args.strain)
        central += 2.0 * float(np.real(np.sum(gradient.gradient * derivative)))
        center_u = center_snapshot.projector[:, :, :, 0]
        minus_side += 2.0 * float(
            np.real(np.sum(gradient.gradient * (center_u - minus_u) / args.strain))
        )
        plus_side += 2.0 * float(
            np.real(np.sum(gradient.gradient * (plus_u - center_u) / args.strain))
        )
        projector_scale = max(float(np.linalg.norm(center_u)), np.finfo(float).tiny)
        projector_change["minus"] = max(
            projector_change["minus"],
            float(np.linalg.norm(minus_u - center_u)) / projector_scale,
        )
        projector_change["plus"] = max(
            projector_change["plus"],
            float(np.linalg.norm(plus_u - center_u)) / projector_scale,
        )

    print(f"k points                                      {len(keys)}")
    print(f"degeneracy tolerance                          {args.degeneracy_tol:.12e} Ry")
    print(f"transport gauge                               {args.transport}")
    print(f"largest transported block                     {largest_block}")
    print(f"max coefficient transport residual, minus     {coefficient_residual['minus']:.12e}")
    print(f"max coefficient transport residual, plus      {coefficient_residual['plus']:.12e}")
    print(f"max relative projector change, minus           {projector_change['minus']:.12e}")
    print(f"max relative projector change, plus            {projector_change['plus']:.12e}")
    print(f"joint normalized-projector derivative          {central:+.12e} Ry")
    print(f"one-sided projector derivative (-h,0)          {minus_side:+.12e} Ry")
    print(f"one-sided projector derivative (0,+h)          {plus_side:+.12e} Ry")
    print(f"one-sided difference                           {plus_side-minus_side:+.12e} Ry")
    if args.volume is not None:
        conversion = 14710.5164057 / (3.0 * args.volume)
        print(f"hydrostatic mean                              {central*conversion:+.9f} GPa")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        raise SystemExit(2)
