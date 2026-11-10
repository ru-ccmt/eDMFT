#!/usr/bin/env python3
"""Finite-difference check of the DMFT ``M_Sigma * dO`` stress channel.

The analytic diagnostic in ``dmft2`` evaluates, before the final volume and
mixer conventions are applied,

    D_q = 2 Re sum_k Tr[M_Sigma(k) O_,q(k)] .

This program obtains the same object without any LAPW stress kernel:

    O_,q^band(k) = A_0(k)^dagger
                     [O_k(+h)-O_k(-h)]/(2h) A_0(k),
    D_q^FD       = 2 Re sum_k Tr[M_Sigma,0(k) O_,q^band(k)].

``O(+h)`` and ``O(-h)`` are produced by the environment-gated dump in
``src/dftKS/lapw1.FPP``.  ``A_0`` and the *un-diagonalized* ``M_Sigma`` must
come from one central dmft2 run.  Keeping the matrix un-diagonalized removes
all pseudo-band ordering and degeneracy ambiguity.

Example
-------

Run LAPW1 in the two strained directories with distinct absolute prefixes::

    DFTKS_OVERLAP_DUMP=$PWD/fd/lapw_overlap ~/dbin/x_dmft lapw1

After a central dmft2 run has written ``dmft_msigma.k########.bin`` files::

    python check_sigma_overlap_fd.py \
        --minus ../strain_m/fd/lapw_overlap \
        --plus  ../strain_p/fd/lapw_overlap \
        --reference ../center/fd/lapw_overlap \
        --central ../center/fd/dmft_msigma \
        --strain 0.001 \
        --analytic -0.123456

``--analytic`` is optional.  Its value is the corresponding raw
``DMFT_STRESS_DIAG_SIGMA_MOMENT_OVERLAP`` component in Ry, *not* GPa and not
divided by the cell volume.

Required dmft2 companion dump (not implemented by this script)
----------------------------------------------------------------

The proposed environment gate is ``EDMFT_MSIGMA_DUMP=PREFIX``.  It writes
``PREFIX.k########.bin`` at each k point immediately after ``sigmaMoment`` is
formed in ``l2main.f90``, before ``Aweight`` or ``sigmaMoment`` is
diagonalized.  The stream format is:

* character*8: ``EDMFTMS1``
* int32: endian marker ``0x01020304``
* int32[10]: version=1, global-k, MPI-rank, nbasis, nnlo, nbands,
  ncomponents, DM_nemin, DM_nemax, flags
* real64[3]: fractional k point
* int32[3,nbasis]: ordered basis labels, first ``k3(:,1:n0-nnlo)`` and then
  ``k3lo(:,1:nnlo)``
* complex128[nbasis,nbands,ncomponents], Fortran order:
  ``As(1:n0,DM_nemin:DM_nemax,1:iso)``
* complex128[nbands,nbands], Fortran order: ``sigmaMoment``

Flag bit 0 must be set to assert that ``sigmaMoment`` already contains the
normalized k weight.  This is true for the current construction because both
``Aweight`` and ``AEweight`` contain ``wkp``.  The helper consequently adds
no second k weight.  Its component sum also matches the current
``cmp_dmft_overlap_stress_diag`` path and does not insert ``vnorm1``.

The intended call, directly after the nested loop that constructs
``sigmaMoment`` in ``l2main.f90``, is conceptually::

    call WriteSigmaMomentOverlapDump(ikp,myrank,n0,nnlo,iso, \
         DM_nemin,DM_nemax,k3,k3lo,bk3,As,sigmaMoment)

The writer should return immediately when ``EDMFT_MSIGMA_DUMP`` is unset.
For ``iso > 1`` the helper requires ``--allow-shared-overlap-components``:
that assertion is valid only when every spinor component uses the same
generalized overlap.  Otherwise separate overlap dumps are required.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
import sys

import numpy as np


ENDIAN_MARKER = 0x01020304


@dataclass(frozen=True)
class OverlapDump:
    path: Path
    kpoint: int
    rank: int
    nbasis: int
    nlo: int
    kfrac: np.ndarray
    labels: np.ndarray
    overlap: np.ndarray


@dataclass(frozen=True)
class SigmaMomentDump:
    path: Path
    kpoint: int
    rank: int
    nbasis: int
    nnlo: int
    nbands: int
    ncomponents: int
    dm_nemin: int
    dm_nemax: int
    flags: int
    kfrac: np.ndarray
    labels: np.ndarray
    coefficients: np.ndarray
    moment: np.ndarray


def _read_exact(handle, count: int, dtype: np.dtype, what: str) -> np.ndarray:
    data = np.fromfile(handle, dtype=dtype, count=count)
    if data.size != count:
        raise ValueError(
            f"truncated {what}: expected {count} values, found {data.size}"
        )
    return data


def _read_preamble(handle, expected_magic: bytes) -> str:
    magic = handle.read(8)
    if magic != expected_magic:
        raise ValueError(f"bad magic {magic!r}; expected {expected_magic!r}")
    raw_marker = handle.read(4)
    if len(raw_marker) != 4:
        raise ValueError("truncated endian marker")
    if int.from_bytes(raw_marker, "little", signed=False) == ENDIAN_MARKER:
        return "<"
    if int.from_bytes(raw_marker, "big", signed=False) == ENDIAN_MARKER:
        return ">"
    raise ValueError(f"invalid endian marker bytes {raw_marker.hex()}")


def _require_eof(handle, path: Path) -> None:
    trailing = handle.read(1)
    if trailing:
        raise ValueError(f"unexpected trailing data in {path}")


def read_overlap(path: Path) -> OverlapDump:
    with path.open("rb") as handle:
        endian = _read_preamble(handle, b"EDMFTOV1")
        i4 = np.dtype(endian + "i4")
        f8 = np.dtype(endian + "f8")
        c16 = np.dtype(endian + "c16")
        header = _read_exact(handle, 8, i4, "overlap header").astype(np.int64)
        version, kpoint, rank, nbasis, nlo, storage, triangle, convention = header
        if version != 1:
            raise ValueError(f"unsupported overlap version {version}")
        if storage != 2 or triangle != 1 or convention != 1:
            raise ValueError(
                "unsupported overlap encoding "
                f"storage={storage}, triangle={triangle}, convention={convention}"
            )
        if nbasis <= 0 or nlo < 0 or nlo > nbasis:
            raise ValueError(f"invalid overlap dimensions nbasis={nbasis}, nlo={nlo}")
        kfrac = _read_exact(handle, 3, f8, "k point").astype(np.float64)
        labels = _read_exact(handle, 3 * nbasis, i4, "basis labels").reshape(
            (3, nbasis), order="F"
        )
        packed_size = nbasis * (nbasis + 1) // 2
        packed = _read_exact(handle, packed_size, c16, "packed overlap")
        _require_eof(handle, path)

    overlap = np.empty((nbasis, nbasis), dtype=np.complex128)
    offset = 0
    for column in range(nbasis):
        count = column + 1
        values = packed[offset : offset + count]
        overlap[:count, column] = values
        overlap[column, :column] = values[:column].conj()
        offset += count
    return OverlapDump(
        path=path,
        kpoint=int(kpoint),
        rank=int(rank),
        nbasis=int(nbasis),
        nlo=int(nlo),
        kfrac=kfrac,
        labels=labels,
        overlap=overlap,
    )


def read_sigma_moment(path: Path) -> SigmaMomentDump:
    with path.open("rb") as handle:
        endian = _read_preamble(handle, b"EDMFTMS1")
        i4 = np.dtype(endian + "i4")
        f8 = np.dtype(endian + "f8")
        c16 = np.dtype(endian + "c16")
        header = _read_exact(handle, 10, i4, "M_Sigma header").astype(np.int64)
        (
            version,
            kpoint,
            rank,
            nbasis,
            nnlo,
            nbands,
            ncomponents,
            dm_nemin,
            dm_nemax,
            flags,
        ) = header
        if version != 1:
            raise ValueError(f"unsupported M_Sigma version {version}")
        if nbasis <= 0 or nnlo < 0 or nnlo > nbasis:
            raise ValueError(f"invalid basis dimensions nbasis={nbasis}, nnlo={nnlo}")
        if nbands <= 0 or dm_nemax - dm_nemin + 1 != nbands:
            raise ValueError(
                f"inconsistent band range [{dm_nemin},{dm_nemax}] for {nbands} bands"
            )
        if ncomponents <= 0:
            raise ValueError(f"invalid component count {ncomponents}")
        if not flags & 1:
            raise ValueError(
                "M_Sigma dump does not assert that the normalized k weight is included"
            )
        kfrac = _read_exact(handle, 3, f8, "k point").astype(np.float64)
        labels = _read_exact(handle, 3 * nbasis, i4, "basis labels").reshape(
            (3, nbasis), order="F"
        )
        coefficients = _read_exact(
            handle,
            nbasis * nbands * ncomponents,
            c16,
            "central vector coefficients",
        ).reshape((nbasis, nbands, ncomponents), order="F")
        moment = _read_exact(
            handle, nbands * nbands, c16, "M_Sigma matrix"
        ).reshape((nbands, nbands), order="F")
        _require_eof(handle, path)

    return SigmaMomentDump(
        path=path,
        kpoint=int(kpoint),
        rank=int(rank),
        nbasis=int(nbasis),
        nnlo=int(nnlo),
        nbands=int(nbands),
        ncomponents=int(ncomponents),
        dm_nemin=int(dm_nemin),
        dm_nemax=int(dm_nemax),
        flags=int(flags),
        kfrac=kfrac,
        labels=labels,
        coefficients=coefficients,
        moment=moment,
    )


def _files_for_prefix(prefix: Path) -> list[Path]:
    parent = prefix.parent if str(prefix.parent) else Path(".")
    files = sorted(parent.glob(prefix.name + ".k????????.bin"))
    if not files:
        raise FileNotFoundError(f"no files found for prefix {prefix}")
    return files


def _index_overlap_prefix(prefix: Path) -> dict[int, Path]:
    indexed: dict[int, Path] = {}
    for path in _files_for_prefix(prefix):
        dump = read_overlap(path)
        if dump.kpoint in indexed:
            raise ValueError(f"duplicate overlap k point {dump.kpoint} for {prefix}")
        indexed[dump.kpoint] = path
    return indexed


def _index_sigma_prefix(prefix: Path) -> dict[int, Path]:
    indexed: dict[int, Path] = {}
    for path in _files_for_prefix(prefix):
        dump = read_sigma_moment(path)
        if dump.kpoint in indexed:
            raise ValueError(f"duplicate M_Sigma k point {dump.kpoint} for {prefix}")
        indexed[dump.kpoint] = path
    return indexed


def _max_hermiticity_error(matrix: np.ndarray) -> float:
    return float(np.max(np.abs(matrix - matrix.conj().T)))


def _overlap_in_reference_basis(
    reference: SigmaMomentDump, overlap: OverlapDump, role: str
) -> np.ndarray:
    """Return ``overlap`` ordered like the central vector coefficients.

    A hydrostatic strain can reorder reciprocal vectors of equal or nearly
    equal length even when the fixed-G *set* is unchanged.  That is a harmless
    basis permutation, not a failed finite-difference control.  Plane-wave
    labels are unique and can therefore be mapped exactly.  Local-orbital
    rows are not uniquely identified by their three integer labels, so their
    order is required to remain identical.
    """
    if overlap.kpoint != reference.kpoint:
        raise ValueError(
            f"k-point mismatch for {role}: {overlap.kpoint} != {reference.kpoint}"
        )
    if overlap.nbasis != reference.nbasis or overlap.nlo != reference.nnlo:
        raise ValueError(
            f"basis-size mismatch at k={reference.kpoint} for {role}: "
            f"({overlap.nbasis},{overlap.nlo}) != "
            f"({reference.nbasis},{reference.nnlo})"
        )
    plane_waves = reference.nbasis - reference.nnlo
    if not np.array_equal(
        overlap.labels[:, plane_waves:], reference.labels[:, plane_waves:]
    ):
        mismatch = np.argwhere(
            overlap.labels[:, plane_waves:] != reference.labels[:, plane_waves:]
        )[0]
        axis, local_row = (int(mismatch[0]), int(mismatch[1]))
        basis = plane_waves + local_row
        raise ValueError(
            f"local-orbital order mismatch at k={reference.kpoint} for {role}, "
            f"component={axis + 1}, row={basis + 1}: "
            f"{overlap.labels[axis, basis]} != {reference.labels[axis, basis]}"
        )

    endpoint_labels = [
        tuple(int(value) for value in overlap.labels[:, row])
        for row in range(plane_waves)
    ]
    reference_labels = [
        tuple(int(value) for value in reference.labels[:, row])
        for row in range(plane_waves)
    ]
    if len(set(endpoint_labels)) != plane_waves:
        raise ValueError(f"non-unique plane-wave labels at k={reference.kpoint} for {role}")
    endpoint_index = {label: row for row, label in enumerate(endpoint_labels)}
    if set(endpoint_index) != set(reference_labels):
        missing = sorted(set(reference_labels) - set(endpoint_index))
        extra = sorted(set(endpoint_index) - set(reference_labels))
        raise ValueError(
            f"fixed-G set mismatch at k={reference.kpoint} for {role}: "
            f"missing={missing[:5]}, extra={extra[:5]}"
        )

    permutation = [endpoint_index[label] for label in reference_labels]
    permutation.extend(range(plane_waves, reference.nbasis))
    return overlap.overlap[np.ix_(permutation, permutation)]


def _band_metric(coefficients: np.ndarray, overlap: np.ndarray) -> np.ndarray:
    nbands = coefficients.shape[1]
    result = np.zeros((nbands, nbands), dtype=np.complex128)
    for component in range(coefficients.shape[2]):
        vectors = coefficients[:, :, component]
        result += vectors.conj().T @ overlap @ vectors
    return result


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Check 2 Re Tr[M_Sigma O_,strain] by fixed-coefficient FD."
    )
    parser.add_argument("--minus", required=True, type=Path, help="-h overlap prefix")
    parser.add_argument("--plus", required=True, type=Path, help="+h overlap prefix")
    parser.add_argument(
        "--reference",
        type=Path,
        help="central overlap prefix; midpoint (-h,+h) is used if omitted",
    )
    parser.add_argument(
        "--central", required=True, type=Path, help="central dmft2 M_Sigma prefix"
    )
    parser.add_argument(
        "--strain", required=True, type=float, help="positive one-sided strain h"
    )
    parser.add_argument(
        "--analytic",
        type=float,
        help="raw analytic SIGMA_MOMENT_OVERLAP component in Ry",
    )
    parser.add_argument(
        "--allow-shared-overlap-components",
        action="store_true",
        help="assert all spinor components use this same generalized overlap",
    )
    parser.add_argument(
        "--quiet-kpoints", action="store_true", help="print only aggregate diagnostics"
    )
    parser.add_argument(
        "--tolerance",
        type=float,
        default=1.0e-8,
        help="warning tolerance for normalization and Hermiticity (default: 1e-8)",
    )
    args = parser.parse_args(argv)
    if not np.isfinite(args.strain) or args.strain <= 0.0:
        parser.error("--strain must be finite and positive")
    if not np.isfinite(args.tolerance) or args.tolerance <= 0.0:
        parser.error("--tolerance must be finite and positive")
    return args


def main(argv: list[str] | None = None) -> int:
    args = parse_args(sys.argv[1:] if argv is None else argv)
    minus_files = _index_overlap_prefix(args.minus)
    plus_files = _index_overlap_prefix(args.plus)
    reference_files = (
        _index_overlap_prefix(args.reference) if args.reference is not None else None
    )
    central_files = _index_sigma_prefix(args.central)

    required = set(central_files)
    for role, available in (("minus", set(minus_files)), ("plus", set(plus_files))):
        if available != required:
            raise ValueError(
                f"{role} k-point set differs from central: "
                f"missing={sorted(required - available)}, extra={sorted(available - required)}"
            )
    if reference_files is not None and set(reference_files) != required:
        available = set(reference_files)
        raise ValueError(
            "reference k-point set differs from central: "
            f"missing={sorted(required - available)}, extra={sorted(available - required)}"
        )

    total = 0.0
    diagonal_total = 0.0
    offdiagonal_total = 0.0
    worst_norm = 0.0
    worst_overlap_herm = 0.0
    worst_moment_herm = 0.0
    warned_kfrac = False

    for kpoint in sorted(required):
        central = read_sigma_moment(central_files[kpoint])
        minus = read_overlap(minus_files[kpoint])
        plus = read_overlap(plus_files[kpoint])
        minus_matrix = _overlap_in_reference_basis(central, minus, "minus")
        plus_matrix = _overlap_in_reference_basis(central, plus, "plus")
        if central.ncomponents != 1 and not args.allow_shared_overlap_components:
            raise ValueError(
                f"k={kpoint} has {central.ncomponents} vector components; pass "
                "--allow-shared-overlap-components only after verifying that they "
                "use the same endpoint overlap"
            )
        if not np.allclose(minus.kfrac, plus.kfrac, rtol=0.0, atol=1.0e-12):
            raise ValueError(
                f"endpoint fractional k points differ at k={kpoint}: "
                f"{minus.kfrac} vs {plus.kfrac}"
            )
        # dmft2 reads the vector-file coordinate, which has a historical CXZ
        # representation change.  KZZ is the decisive basis check; report a
        # k-coordinate mismatch only as a warning.
        if not np.allclose(central.kfrac, plus.kfrac, rtol=0.0, atol=1.0e-12):
            warned_kfrac = True

        overlap_herm = max(
            _max_hermiticity_error(minus_matrix),
            _max_hermiticity_error(plus_matrix),
        )
        moment_herm = _max_hermiticity_error(central.moment)
        derivative = (plus_matrix - minus_matrix) / (2.0 * args.strain)

        if reference_files is None:
            reference_overlap = 0.5 * (plus.overlap + minus.overlap)
        else:
            reference = read_overlap(reference_files[kpoint])
            reference_overlap = _overlap_in_reference_basis(
                central, reference, "reference"
            )
            overlap_herm = max(
                overlap_herm, _max_hermiticity_error(reference_overlap)
            )

        metric = _band_metric(central.coefficients, reference_overlap)
        norm_error = float(
            np.max(np.abs(metric - np.eye(central.nbands, dtype=np.complex128)))
        )
        band_derivative = _band_metric(central.coefficients, derivative)
        value = float(2.0 * np.real(np.trace(central.moment @ band_derivative)))
        diagonal_moment = np.diag(np.diag(central.moment))
        diagonal_value = float(
            2.0 * np.real(np.trace(diagonal_moment @ band_derivative))
        )
        offdiagonal_value = value - diagonal_value
        total += value
        diagonal_total += diagonal_value
        offdiagonal_total += offdiagonal_value
        worst_norm = max(worst_norm, norm_error)
        worst_overlap_herm = max(worst_overlap_herm, overlap_herm)
        worst_moment_herm = max(worst_moment_herm, moment_herm)

        if not args.quiet_kpoints:
            print(
                f"k={kpoint:6d}  nbasis={central.nbasis:5d}  "
                f"nbands={central.nbands:4d}  norm_err={norm_error:.3e}  "
                f"M_herm={moment_herm:.3e}  D_k={value:+.12e} Ry"
            )

    print(f"k points                              {len(required)}")
    print(f"max |A^dagger O A-I|                  {worst_norm:.12e}")
    print(f"max overlap Hermiticity error         {worst_overlap_herm:.12e}")
    print(f"max M_Sigma Hermiticity error         {worst_moment_herm:.12e}")
    print(f"2 Re sum_k Tr[M_Sigma O_,strain]      {total:+.12e} Ry")
    print(f"  diagonal M_Sigma contribution       {diagonal_total:+.12e} Ry")
    print(f"  off-diagonal M_Sigma contribution   {offdiagonal_total:+.12e} Ry")
    if args.analytic is not None:
        difference = total - args.analytic
        scale = max(abs(total), abs(args.analytic), np.finfo(float).tiny)
        print(f"raw analytic channel                  {args.analytic:+.12e} Ry")
        print(f"finite difference - analytic          {difference:+.12e} Ry")
        print(f"relative difference                   {difference / scale:+.12e}")
    if warned_kfrac:
        print(
            "warning: central and LAPW1 fractional k coordinates differ; "
            "ordered KZZ labels matched (this can be the CXZ convention)",
            file=sys.stderr,
        )
    if max(worst_norm, worst_overlap_herm, worst_moment_herm) > args.tolerance:
        print(
            "warning: at least one consistency diagnostic exceeds --tolerance",
            file=sys.stderr,
        )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        raise SystemExit(2)
