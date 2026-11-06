#!/usr/bin/env python3
"""Direct fixed-basis check of the stationary DFT+eDMFT one-particle stress.

For each irreducible k point this program evaluates

    D = sum_k Re Tr[rho_k H_{,strain,k} - rhoE_k O_{,strain,k}]

from central ``rho/rhoE`` matrices and frozen-augmentation LAPW1 matrices at
``-h, 0, +h``.  It bypasses every analytic LAPW stress kernel and therefore
tests the complete one-particle Pulay rearrangement as a single object.

The required files are default-off diagnostics:

* ``DFTKS_OVERLAP_DUMP=PREFIX`` in LAPW1 writes
  ``PREFIX.k########.bin`` and ``PREFIX.ham.k########.bin``;
* ``EDMFT_ONEPARTICLE_DUMP=PREFIX`` in dmft2 writes
  ``PREFIX.k########.bin``.

The endpoint LAPW1 runs must use the same reciprocal-vector set/order, local-
orbital construction, radial potential, linearization energies, and RMT as
the central run.  ``DFTKS_CANONICAL_G_ORDER=1`` makes the coefficient ordering
strain-stable; RKMAX must still be scaled so the G set itself is unchanged.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
import re
import sys

import numpy as np


ENDIAN_MARKER = 0x01020304


@dataclass(frozen=True)
class BasisMatrix:
    kpoint: int
    nbasis: int
    nlo: int
    labels: np.ndarray
    matrix: np.ndarray


@dataclass(frozen=True)
class OneParticleDump:
    kpoint: int
    nbasis: int
    nnlo: int
    nfull: int
    ncomponents: int
    nemin: int
    dm_nemin: int
    dm_nemax: int
    labels: np.ndarray
    energies: np.ndarray
    coefficients: np.ndarray
    rho: np.ndarray
    rhoe: np.ndarray


def _preamble(handle, magic: bytes) -> str:
    found = handle.read(8)
    if found != magic:
        raise ValueError(f"bad magic {found!r} in {handle.name}; expected {magic!r}")
    marker = handle.read(4)
    if int.from_bytes(marker, "little") == ENDIAN_MARKER:
        return "<"
    if int.from_bytes(marker, "big") == ENDIAN_MARKER:
        return ">"
    raise ValueError(f"bad endian marker in {handle.name}")


def _exact(handle, count: int, dtype: str, what: str) -> np.ndarray:
    data = np.fromfile(handle, dtype=dtype, count=count)
    if data.size != count:
        raise ValueError(f"truncated {what} in {handle.name}")
    return data


def read_basis_matrix(path: Path, magic: bytes) -> BasisMatrix:
    with path.open("rb") as handle:
        endian = _preamble(handle, magic)
        header = _exact(handle, 8, endian + "i4", "matrix header").astype(int)
        version, kpoint, _rank, nbasis, nlo, storage, triangle, convention = header
        if (version, storage, triangle, convention) != (1, 2, 1, 1):
            raise ValueError(f"unsupported matrix encoding {header.tolist()} in {path}")
        _exact(handle, 3, endian + "f8", "k point")
        labels = _exact(handle, 3 * nbasis, endian + "i4", "basis labels").reshape(
            (3, nbasis), order="F"
        )
        packed = _exact(
            handle, nbasis * (nbasis + 1) // 2, endian + "c16", "matrix"
        )
        if handle.read(1):
            raise ValueError(f"trailing data in {path}")

    matrix = np.empty((nbasis, nbasis), dtype=np.complex128)
    offset = 0
    for column in range(nbasis):
        count = column + 1
        values = packed[offset : offset + count]
        matrix[:count, column] = values
        matrix[column, :column] = values[:column].conj()
        offset += count
    return BasisMatrix(int(kpoint), int(nbasis), int(nlo), labels, matrix)


def read_one_particle(path: Path) -> OneParticleDump:
    with path.open("rb") as handle:
        endian = _preamble(handle, b"EDMFTRH1")
        header = _exact(handle, 12, endian + "i4", "rho/rhoE header").astype(int)
        (
            version,
            kpoint,
            _rank,
            nbasis,
            nnlo,
            nfull,
            ncomponents,
            nemin,
            _dm_nemin,
            dm_nemax,
            flags,
            _reserved,
        ) = header
        if version != 1 or flags & 1 == 0 or dm_nemax - nemin + 1 != nfull:
            raise ValueError(f"unsupported rho/rhoE encoding {header.tolist()} in {path}")
        _exact(handle, 3, endian + "f8", "k point")
        labels = _exact(handle, 3 * nbasis, endian + "i4", "basis labels").reshape(
            (3, nbasis), order="F"
        )
        energies = _exact(handle, nfull, endian + "f8", "band energies")
        coefficients = _exact(
            handle,
            nbasis * nfull * ncomponents,
            endian + "c16",
            "LAPW coefficients",
        ).reshape((nbasis, nfull, ncomponents), order="F")
        rho = _exact(handle, nfull * nfull, endian + "c16", "rho").reshape(
            (nfull, nfull), order="F"
        )
        rhoe = _exact(handle, nfull * nfull, endian + "c16", "rhoE").reshape(
            (nfull, nfull), order="F"
        )
        if handle.read(1):
            raise ValueError(f"trailing data in {path}")
    return OneParticleDump(
        int(kpoint),
        int(nbasis),
        int(nnlo),
        int(nfull),
        int(ncomponents),
        int(nemin),
        int(_dm_nemin),
        int(dm_nemax),
        labels,
        energies,
        coefficients,
        rho,
        rhoe,
    )


def files(prefix: Path, suffix: str = ".k????????.bin") -> dict[int, Path]:
    found: dict[int, Path] = {}
    for path in sorted(prefix.parent.glob(prefix.name + suffix)):
        if ".ham." in path.name:
            continue
        with path.open("rb") as handle:
            magic = handle.read(8)
        if magic == b"EDMFTRH1":
            kpoint = read_one_particle(path).kpoint
        elif magic == b"EDMFTOV1":
            kpoint = read_basis_matrix(path, magic).kpoint
        else:
            continue
        found[kpoint] = path
    if not found:
        raise FileNotFoundError(f"no diagnostic files for prefix {prefix}")
    return found


def matrix_path(prefix: Path, kpoint: int, hamiltonian: bool) -> Path:
    middle = ".ham.k" if hamiltonian else ".k"
    return prefix.parent / f"{prefix.name}{middle}{kpoint:08d}.bin"


def aligned(reference: OneParticleDump, matrix: BasisMatrix, role: str) -> np.ndarray:
    if (matrix.nbasis, matrix.nlo) != (reference.nbasis, reference.nnlo):
        raise ValueError(
            f"basis-size mismatch at k={reference.kpoint} for {role}: "
            f"({matrix.nbasis},{matrix.nlo}) != ({reference.nbasis},{reference.nnlo})"
        )
    npw = reference.nbasis - reference.nnlo
    if not np.array_equal(matrix.labels[:, npw:], reference.labels[:, npw:]):
        raise ValueError(f"local-orbital order mismatch at k={reference.kpoint} for {role}")
    lookup = {tuple(matrix.labels[:, row]): row for row in range(npw)}
    wanted = [tuple(reference.labels[:, row]) for row in range(npw)]
    if len(lookup) != npw or set(lookup) != set(wanted):
        raise ValueError(f"fixed-G set mismatch at k={reference.kpoint} for {role}")
    permutation = [lookup[label] for label in wanted]
    permutation.extend(range(npw, reference.nbasis))
    return matrix.matrix[np.ix_(permutation, permutation)]


def band_matrix(coefficients: np.ndarray, operator: np.ndarray) -> np.ndarray:
    result = np.zeros((coefficients.shape[1], coefficients.shape[1]), complex)
    for component in range(coefficients.shape[2]):
        vectors = coefficients[:, :, component]
        result += vectors.conj().T @ operator @ vectors
    return result


def read_dft_weights(path: Path) -> dict[int, dict[int, float]]:
    """Read the ordinary WIEN2k ``case.weight`` band occupations."""
    result: dict[int, dict[int, float]] = {}
    for line_number, line in enumerate(path.read_text().splitlines(), 1):
        fields = line.split()
        if len(fields) != 4:
            continue
        try:
            kpoint = int(fields[0])
            band = int(fields[1])
            float(fields[2])  # energy: checked but central dump is authoritative
            weight_text = fields[3]
            try:
                weight = float(weight_text)
            except ValueError:
                # Fortran E format can omit the exponent letter once the
                # exponent needs three digits, e.g. ``7.60-134``.
                match = re.fullmatch(r"(.+?)([+-]\d{3,})", weight_text)
                if match is None:
                    raise
                weight = float(match.group(1) + "E" + match.group(2))
        except ValueError:
            continue
        if kpoint <= 0 or band <= 0:
            raise ValueError(f"invalid band label at {path}:{line_number}")
        if band in result.setdefault(kpoint, {}):
            raise ValueError(f"duplicate k/band at {path}:{line_number}")
        result[kpoint][band] = weight
    if not result:
        raise ValueError(f"no band occupations found in {path}")
    return result


def read_energy_file(path: Path) -> dict[int, dict[int, float]]:
    """Read band energies without relying on whitespace between k coordinates."""
    lines = path.read_text().splitlines()
    result: dict[int, dict[int, float]] = {}
    line_index = 0
    while line_index < len(lines) - 1:
        fields = lines[line_index].split()
        next_fields = lines[line_index + 1].split()
        if len(fields) >= 4 and len(next_fields) == 2 and next_fields[0] == "1":
            try:
                kpoint = int(fields[-4])
                number_bands = int(fields[-2])
            except ValueError:
                line_index += 1
                continue
            bands: dict[int, float] = {}
            for offset in range(1, number_bands + 1):
                if line_index + offset >= len(lines):
                    raise ValueError(f"truncated k point {kpoint} in {path}")
                band_fields = lines[line_index + offset].split()
                if len(band_fields) != 2:
                    raise ValueError(
                        f"bad band record at {path}:{line_index+offset+1}"
                    )
                bands[int(band_fields[0])] = float(band_fields[1])
            result[kpoint] = bands
            line_index += number_bands + 1
        else:
            line_index += 1
    if not result:
        raise ValueError(f"no k-point energy blocks found in {path}")
    return result


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--minus", required=True, type=Path)
    parser.add_argument("--center", required=True, type=Path)
    parser.add_argument("--plus", required=True, type=Path)
    parser.add_argument("--moments", required=True, type=Path)
    parser.add_argument("--strain", required=True, type=float)
    parser.add_argument("--volume", type=float, help="central volume in Bohr^3")
    parser.add_argument(
        "--dft-weight",
        type=Path,
        help="optional central case.weight; also report DFT and DMFT-minus-DFT pairs",
    )
    parser.add_argument("--minus-energy", type=Path)
    parser.add_argument("--plus-energy", type=Path)
    args = parser.parse_args(argv)
    if not np.isfinite(args.strain) or args.strain <= 0:
        parser.error("--strain must be positive")
    if args.volume is not None and (not np.isfinite(args.volume) or args.volume <= 0):
        parser.error("--volume must be positive")
    if (args.minus_energy is None) != (args.plus_energy is None):
        parser.error("--minus-energy and --plus-energy must be supplied together")
    if args.minus_energy is not None and args.dft_weight is None:
        parser.error("endpoint energy comparison also requires --dft-weight")
    return args


def main(argv: list[str] | None = None) -> int:
    args = parse_args(sys.argv[1:] if argv is None else argv)
    moment_files = files(args.moments)
    dft_weights = read_dft_weights(args.dft_weight) if args.dft_weight else None
    endpoint_energies = None
    if args.minus_energy is not None:
        endpoint_energies = (
            read_energy_file(args.minus_energy),
            read_energy_file(args.plus_energy),
        )
    h = args.strain
    totals = {"h": 0.0, "o": 0.0, "pair": 0.0, "minus": 0.0, "plus": 0.0}
    sigma_overlap = 0.0
    sigma_overlap_sides = {"minus": 0.0, "plus": 0.0}
    dft_pair = 0.0
    dft_electrons = 0.0
    eigenvalue_dft_pair = 0.0
    norm_error = 0.0
    eigenpair_error = 0.0
    hermiticity = 0.0

    for kpoint, path in sorted(moment_files.items()):
        reference = read_one_particle(path)
        values: dict[tuple[str, str], np.ndarray] = {}
        for role, prefix in (("minus", args.minus), ("center", args.center), ("plus", args.plus)):
            overlap_dump = read_basis_matrix(
                matrix_path(prefix, kpoint, False), b"EDMFTOV1"
            )
            hamiltonian_dump = read_basis_matrix(
                matrix_path(prefix, kpoint, True), b"EDMFTHV1"
            )
            values[(role, "o")] = aligned(reference, overlap_dump, role + " O")
            values[(role, "h")] = aligned(reference, hamiltonian_dump, role + " H")

        identity = band_matrix(reference.coefficients, values[("center", "o")])
        norm_error = max(norm_error, float(np.max(np.abs(identity - np.eye(reference.nfull)))))
        h_band = band_matrix(reference.coefficients, values[("center", "h")])
        eigenpair_error = max(
            eigenpair_error,
            float(np.max(np.abs(h_band - np.diag(reference.energies)))),
        )
        for value in values.values():
            hermiticity = max(hermiticity, float(np.max(np.abs(value - value.conj().T))))

        do = (values[("plus", "o")] - values[("minus", "o")]) / (2 * h)
        dh = (values[("plus", "h")] - values[("minus", "h")]) / (2 * h)
        do_band = band_matrix(reference.coefficients, do)
        dh_band = band_matrix(reference.coefficients, dh)
        h_part = float(np.real(np.trace(reference.rho @ dh_band)))
        o_part = float(np.real(np.trace(reference.rhoe @ do_band)))
        totals["h"] += h_part
        totals["o"] += o_part
        totals["pair"] += h_part - o_part

        energy_sum = 0.5 * (
            reference.energies[:, None] + reference.energies[None, :]
        ) * reference.rho
        msigma = reference.rhoe - energy_sum
        sigma_overlap += float(-np.real(np.trace(msigma @ do_band)))

        # Keep the two one-sided metric derivatives separate.  The symmetric
        # number alone can look well converged even when a finite strain step
        # is too large and two sizeable curvatures cancel.  These values are
        # consumed by the guarded reference-correction generator; they never
        # alter a dmft2 calculation by themselves.
        for side, left, right in (
            ("minus", "minus", "center"), ("plus", "center", "plus")
        ):
            one_do = (values[(right, "o")] - values[(left, "o")]) / h
            sigma_overlap_sides[side] += float(
                -np.real(
                    np.trace(
                        msigma
                        @ band_matrix(reference.coefficients, one_do)
                    )
                )
            )

        if dft_weights is not None:
            if kpoint not in dft_weights:
                raise ValueError(f"k point {kpoint} is missing from {args.dft_weight}")
            occupations = np.empty(reference.nfull)
            for offset, band in enumerate(
                range(reference.nemin, reference.dm_nemax + 1)
            ):
                try:
                    occupations[offset] = dft_weights[kpoint][band]
                except KeyError as error:
                    raise ValueError(
                        f"band {band} at k point {kpoint} is missing from "
                        f"{args.dft_weight}"
                    ) from error
            dft_electrons += float(np.sum(occupations))
            dft_pair += float(
                np.real(
                    np.sum(occupations * np.diag(dh_band))
                    - np.sum(occupations * reference.energies * np.diag(do_band))
                )
            )
            if endpoint_energies is not None:
                minus_energy, plus_energy = endpoint_energies
                for offset, band in enumerate(
                    range(reference.nemin, reference.dm_nemax + 1)
                ):
                    try:
                        derivative = (
                            plus_energy[kpoint][band]
                            - minus_energy[kpoint][band]
                        ) / (2 * h)
                    except KeyError as error:
                        raise ValueError(
                            f"endpoint energy for band {band}, k point {kpoint} is missing"
                        ) from error
                    eigenvalue_dft_pair += occupations[offset] * derivative

        for side, left, right in (
            ("minus", "minus", "center"), ("plus", "center", "plus")
        ):
            one_do = (values[(right, "o")] - values[(left, "o")]) / h
            one_dh = (values[(right, "h")] - values[(left, "h")]) / h
            totals[side] += float(
                np.real(
                    np.trace(reference.rho @ band_matrix(reference.coefficients, one_dh))
                    - np.trace(reference.rhoe @ band_matrix(reference.coefficients, one_do))
                )
            )

    print(f"k points                                      {len(moment_files)}")
    print(f"max |A^dagger O A-I|                          {norm_error:.12e}")
    print(f"max |A^dagger H A-diag(epsilon)|              {eigenpair_error:.12e} Ry")
    print(f"max H/O Hermiticity error                     {hermiticity:.12e}")
    print(f"Tr[rho H_,strain]                             {totals['h']:+.12e} Ry")
    print(f"Tr[rhoE O_,strain]                            {totals['o']:+.12e} Ry")
    print(f"stationary one-particle pair                  {totals['pair']:+.12e} Ry")
    print(f"-Tr[M_Sigma O_,strain]                        {sigma_overlap:+.12e} Ry")
    print(
        "one-sided M_Sigma overlap (-h,0)             "
        f"{sigma_overlap_sides['minus']:+.12e} Ry"
    )
    print(
        "one-sided M_Sigma overlap (0,+h)             "
        f"{sigma_overlap_sides['plus']:+.12e} Ry"
    )
    print(
        "one-sided M_Sigma difference                  "
        f"{sigma_overlap_sides['plus']-sigma_overlap_sides['minus']:+.12e} Ry"
    )
    if dft_weights is not None:
        print(f"DFT electrons represented by stored bands           {dft_electrons:.12e}")
        print(f"DFT stationary one-particle pair              {dft_pair:+.12e} Ry")
        print(f"DMFT minus DFT one-particle pair               {totals['pair']-dft_pair:+.12e} Ry")
        if endpoint_energies is not None:
            print(f"DFT pair from endpoint eigenvalues             {eigenvalue_dft_pair:+.12e} Ry")
            print(f"matrix minus eigenvalue DFT pair               {dft_pair-eigenvalue_dft_pair:+.12e} Ry")
    print(f"one-sided pair (-h,0)                         {totals['minus']:+.12e} Ry")
    print(f"one-sided pair (0,+h)                         {totals['plus']:+.12e} Ry")
    print(f"one-sided difference                          {totals['plus']-totals['minus']:+.12e} Ry")
    if args.volume is not None:
        conversion = 14710.5164057 / (3.0 * args.volume)
        print(f"hydrostatic mean of pair                      {totals['pair']*conversion:+.9f} GPa")
        print(f"hydrostatic mean of M_Sigma overlap           {sigma_overlap*conversion:+.9f} GPa")
        if dft_weights is not None:
            print(f"hydrostatic mean of DFT pair                  {dft_pair*conversion:+.9f} GPa")
            print(
                "hydrostatic mean of DMFT-minus-DFT pair       "
                f"{(totals['pair']-dft_pair)*conversion:+.9f} GPa"
            )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        raise SystemExit(2)
