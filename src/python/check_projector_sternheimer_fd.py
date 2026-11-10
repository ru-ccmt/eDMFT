#!/usr/bin/env python3
"""Central Sternheimer check of the fixed-index DMFT projector response.

The endpoint-projector finite difference is gauge dependent because the two
strained LAPW eigenproblems may choose unrelated bases inside nearly
degenerate blocks.  This checker never compares endpoint eigenvectors.  It
uses central LAPW eigenvectors/projectors and finite-difference coefficient-
space operators to form

    R[a,n] = <a|H_,q-e_n O_,q|n>/(e_n-e_a),
    U_,q[n] = sum(a outside W) conj(R[a,n]) U[a],

then evaluates the exact fixed-normalization stationary contraction

    2 Re sum(k,n) C[k,n] U_,q[k,n].

``C`` already contains the irreducible-k, spin, and Matsubara weight.  The
result must be compared with dmft2's ``WINDOW_TOTAL`` diagnostic, not with
``NORM_WINDOW_TOTAL``: the latter is the additional response of the global
normalization X=(sum P^+P)^(-1/2).

Required diagnostics:

* ``EDMFT_STRESS_STERNHEIMER_DUMP=PREFIX`` in the central dmft2 run;
* ``DFTKS_OVERLAP_DUMP=PREFIX`` in fixed-potential LAPW1 runs at -h and +h.

With ``--moments`` the checker also evaluates the complete stored-band
identity between the fixed LAPW-coefficient representation and the
instantaneous orthonormal KS-band representation.  The latter needs the
within-window rotations and the diagonal metric connection

    R[n,n] = -O_,q[n,n]/2

in addition to the complementary-band response.  This is the decisive test
for a missing within-window projector term; no endpoint eigenvector gauge is
used.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
import sys

import numpy as np

from check_oneparticle_stress_fd import (
    ENDIAN_MARKER,
    aligned,
    band_matrix,
    matrix_path,
    read_basis_matrix,
    read_one_particle,
    files as oneparticle_files,
)


def aligned_for_sternheimer(reference, matrix, role: str, allow_lo_label_drift: bool):
    """Align a strained coefficient-space matrix to the central basis.

    Plane waves are identified by their integer reciprocal vectors and must
    therefore be permuted explicitly.  LAPW local orbitals, on the other
    hand, have a stable atom/l/m/index order.  The three integers written in
    the diagnostic file for an LO are the auxiliary ``k3lo`` vector used to
    represent that orbital at the particular strained k point; they are not
    an LO identity and may change with strain.  The general one-particle
    checker remains deliberately strict.  This Sternheimer checker permits
    the physically correct index-order alignment only through the explicit
    command-line opt-in below.
    """

    if not allow_lo_label_drift:
        return aligned(reference, matrix, role)
    if (matrix.nbasis, matrix.nlo) != (reference.nbasis, reference.nnlo):
        raise ValueError(
            f"basis-size mismatch at k={reference.kpoint} for {role}: "
            f"({matrix.nbasis},{matrix.nlo}) != "
            f"({reference.nbasis},{reference.nnlo})"
        )
    npw = reference.nbasis - reference.nnlo
    lookup = {tuple(matrix.labels[:, row]): row for row in range(npw)}
    wanted = [tuple(reference.labels[:, row]) for row in range(npw)]
    if len(lookup) != npw or set(lookup) != set(wanted):
        raise ValueError(f"fixed-G set mismatch at k={reference.kpoint} for {role}")
    permutation = [lookup[label] for label in wanted]
    # LO columns follow the stable atom/l/m/index order in both LAPW1 and
    # dmft2.  Do not reorder them by k3lo: those labels are neither unique
    # nor invariant under a strained reciprocal representation.
    permutation.extend(range(npw, reference.nbasis))
    return matrix.matrix[np.ix_(permutation, permutation)]


@dataclass(frozen=True)
class SternheimerDump:
    kpoint: int
    nbasis: int
    nlo: int
    nstored: int
    ncomponents: int
    nemin: int
    nemax: int
    dm_nemin: int
    dm_nemax: int
    maxdim: int
    norbitals: int
    labels: np.ndarray
    energies: np.ndarray
    coefficients: np.ndarray
    projector: np.ndarray
    gradient: np.ndarray


@dataclass(frozen=True)
class NormalizationDump:
    kpoint: int
    nstored: int
    nwindow: int
    maxdim: int
    norbitals: int
    nemin: int
    nemax: int
    dm_nemin: int
    dm_nemax: int
    energies: np.ndarray
    projector: np.ndarray
    gradient: np.ndarray


@dataclass(frozen=True)
class ResponseMatrixDump:
    kpoint: int
    nstored: int
    nemin: int
    nemax: int
    dm_nemin: int
    dm_nemax: int
    density: np.ndarray
    energy: np.ndarray


def exact(handle, count: int, dtype: str, what: str) -> np.ndarray:
    value = np.fromfile(handle, dtype=dtype, count=count)
    if value.size != count:
        raise ValueError(f"truncated {what} in {handle.name}")
    return value


def read_dump(path: Path) -> SternheimerDump:
    with path.open("rb") as handle:
        if handle.read(8) != b"EDMFTST1":
            raise ValueError(f"bad Sternheimer magic in {path}")
        marker = handle.read(4)
        if int.from_bytes(marker, "little") == ENDIAN_MARKER:
            endian = "<"
        elif int.from_bytes(marker, "big") == ENDIAN_MARKER:
            endian = ">"
        else:
            raise ValueError(f"bad endian marker in {path}")
        header = exact(handle, 13, endian + "i4", "header").astype(int)
        (
            version,
            kpoint,
            _rank,
            nbasis,
            nlo,
            nstored,
            ncomponents,
            nemin,
            nemax,
            dm_nemin,
            dm_nemax,
            maxdim,
            norbitals,
        ) = header
        if version != 1 or nstored != nemax - nemin + 1:
            raise ValueError(f"unsupported Sternheimer header {header.tolist()}")
        nwindow = dm_nemax - dm_nemin + 1
        exact(handle, 3, endian + "f8", "k point")
        labels = exact(handle, 3 * nbasis, endian + "i4", "labels").reshape(
            (3, nbasis), order="F"
        )
        energies = exact(handle, nstored, endian + "f8", "energies")
        coefficients = exact(
            handle,
            nbasis * nstored * ncomponents,
            endian + "c16",
            "coefficients",
        ).reshape((nbasis, nstored, ncomponents), order="F")
        projector = exact(
            handle,
            nstored * maxdim * norbitals,
            endian + "c16",
            "projector",
        ).reshape((nstored, maxdim, norbitals), order="F")
        gradient = exact(
            handle,
            nwindow * maxdim * norbitals,
            endian + "c16",
            "gradient",
        ).reshape((nwindow, maxdim, norbitals), order="F")
        if handle.read(1):
            raise ValueError(f"trailing data in {path}")
    return SternheimerDump(
        kpoint,
        nbasis,
        nlo,
        nstored,
        ncomponents,
        nemin,
        nemax,
        dm_nemin,
        dm_nemax,
        maxdim,
        norbitals,
        labels,
        energies,
        coefficients,
        projector,
        gradient,
    )


def dump_files(prefix: Path) -> dict[int, Path]:
    result: dict[int, Path] = {}
    for path in sorted(prefix.parent.glob(prefix.name + ".k????????.bin")):
        value = read_dump(path)
        if value.kpoint in result:
            raise ValueError(f"duplicate k point {value.kpoint}")
        result[value.kpoint] = path
    if not result:
        raise FileNotFoundError(f"no Sternheimer files for prefix {prefix}")
    return result


def read_normalization(path: Path) -> NormalizationDump:
    with path.open("rb") as handle:
        header = exact(handle, 11, "=i4", "normalization header").astype(int)
        (
            magic,
            version,
            kpoint,
            nstored,
            nwindow,
            maxdim,
            norbitals,
            nemin,
            nemax,
            dm_nemin,
            dm_nemax,
        ) = header
        if (magic, version) != (20261107, 1):
            raise ValueError(f"unsupported normalization dump {header.tolist()}")
        energies = exact(handle, nstored, "=f8", "normalization energies")
        projector = exact(
            handle,
            nstored * maxdim * norbitals,
            "=c16",
            "raw projector",
        ).reshape((nstored, maxdim, norbitals), order="F")
        gradient = exact(
            handle,
            nwindow * maxdim * norbitals,
            "=c16",
            "normalization gradient",
        ).reshape((nwindow, maxdim, norbitals), order="F")
        if handle.read(1):
            raise ValueError(f"trailing data in {path}")
    return NormalizationDump(
        kpoint,
        nstored,
        nwindow,
        maxdim,
        norbitals,
        nemin,
        nemax,
        dm_nemin,
        dm_nemax,
        energies,
        projector,
        gradient,
    )


def normalization_files(prefix: Path) -> dict[int, Path]:
    result: dict[int, Path] = {}
    for path in sorted(prefix.parent.glob(prefix.name + ".k????????.bin")):
        value = read_normalization(path)
        if value.kpoint in result:
            raise ValueError(f"duplicate normalization k point {value.kpoint}")
        result[value.kpoint] = path
    if not result:
        raise FileNotFoundError(f"no normalization files for prefix {prefix}")
    return result


def read_response_matrix(path: Path) -> ResponseMatrixDump:
    with path.open("rb") as handle:
        if handle.read(8) != b"EDMFTRS1":
            raise ValueError(f"bad response-matrix magic in {path}")
        marker = handle.read(4)
        if int.from_bytes(marker, "little") == ENDIAN_MARKER:
            endian = "<"
        elif int.from_bytes(marker, "big") == ENDIAN_MARKER:
            endian = ">"
        else:
            raise ValueError(f"bad endian marker in {path}")
        header = exact(handle, 7, endian + "i4", "response header").astype(int)
        version, kpoint, nstored, nemin, nemax, dm_nemin, dm_nemax = header
        if version != 1 or nstored != nemax - nemin + 1:
            raise ValueError(f"unsupported response header {header.tolist()}")
        density = exact(
            handle, nstored * nstored, endian + "c16", "density response"
        ).reshape((nstored, nstored), order="F")
        energy = exact(
            handle, nstored * nstored, endian + "c16", "energy response"
        ).reshape((nstored, nstored), order="F")
        if handle.read(1):
            raise ValueError(f"trailing data in {path}")
    return ResponseMatrixDump(
        kpoint, nstored, nemin, nemax, dm_nemin, dm_nemax, density, energy
    )


def response_matrix_files(prefix: Path) -> dict[int, Path]:
    result: dict[int, Path] = {}
    for path in sorted(prefix.parent.glob(prefix.name + ".k????????.bin")):
        value = read_response_matrix(path)
        if value.kpoint in result:
            raise ValueError(f"duplicate response k point {value.kpoint}")
        result[value.kpoint] = path
    if not result:
        raise FileNotFoundError(f"no response-matrix files for prefix {prefix}")
    return result


def reference_like(dump: SternheimerDump):
    """Minimal object accepted by check_oneparticle_stress_fd.aligned."""

    class Reference:
        pass

    value = Reference()
    value.kpoint = dump.kpoint
    value.nbasis = dump.nbasis
    value.nnlo = dump.nlo
    value.labels = dump.labels
    return value


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--minus", required=True, type=Path)
    parser.add_argument("--plus", required=True, type=Path)
    parser.add_argument("--central", required=True, type=Path)
    parser.add_argument(
        "--moments",
        type=Path,
        help=(
            "optional EDMFT_ONEPARTICLE_DUMP prefix; test the complete "
            "fixed-coefficient versus moving-KS-band identity"
        ),
    )
    parser.add_argument(
        "--normalization",
        type=Path,
        help="optional companion PREFIX_norm; report the X-response branch",
    )
    parser.add_argument(
        "--response-matrices",
        type=Path,
        help=(
            "optional EDMFT_STRESS_RESPONSE_DUMP prefix; compare the exact "
            "Hermitian matrices entering dmft2's LAPW response pass"
        ),
    )
    parser.add_argument("--strain", required=True, type=float)
    parser.add_argument("--volume", type=float)
    parser.add_argument(
        "--denominator-tol",
        type=float,
        default=1.0e-8,
        help="skip cross-boundary denominators below this many Ry",
    )
    parser.add_argument(
        "--report-pairs",
        type=int,
        default=0,
        help=(
            "print this many smallest-denominator and largest-contribution "
            "window/complement pairs (diagnostic only)"
        ),
    )
    parser.add_argument(
        "--allow-lo-label-drift",
        action="store_true",
        help=(
            "align plane waves by G but retain local orbitals by their "
            "stable index order; needed when strained k3lo labels differ"
        ),
    )
    args = parser.parse_args(argv)
    if args.strain <= 0 or not np.isfinite(args.strain):
        parser.error("--strain must be positive")
    if args.denominator_tol <= 0 or not np.isfinite(args.denominator_tol):
        parser.error("--denominator-tol must be positive")
    if args.volume is not None and (args.volume <= 0 or not np.isfinite(args.volume)):
        parser.error("--volume must be positive")
    return args


def main(argv: list[str] | None = None) -> int:
    args = parse_args(sys.argv[1:] if argv is None else argv)
    paths = dump_files(args.central)
    norm_paths = normalization_files(args.normalization) if args.normalization else None
    moment_paths = oneparticle_files(args.moments) if args.moments else None
    response_paths = (
        response_matrix_files(args.response_matrices)
        if args.response_matrices
        else None
    )
    if norm_paths is not None and set(norm_paths) != set(paths):
        raise ValueError("central and normalization k-point sets differ")
    if moment_paths is not None and set(moment_paths) != set(paths):
        raise ValueError("central Sternheimer and one-particle k-point sets differ")
    if response_paths is not None and set(response_paths) != set(paths):
        raise ValueError("central Sternheimer and response-matrix k-point sets differ")
    total = 0.0
    response_hamiltonian_total = 0.0
    response_overlap_total = 0.0
    normalization_total = 0.0
    normalization_hamiltonian_total = 0.0
    normalization_overlap_total = 0.0
    internal_total = 0.0
    diagonal_metric_total = 0.0
    all_band_total = 0.0
    normalization_internal_total = 0.0
    normalization_diagonal_total = 0.0
    normalization_all_total = 0.0
    eigenvalue_total = 0.0
    eigenvalue_below_window = 0.0
    eigenvalue_window = 0.0
    fixed_pair_total = 0.0
    b_rotation_total = 0.0
    sigma_overlap_total = 0.0
    sigma_overlap_diagonal_total = 0.0
    sigma_overlap_offdiagonal_total = 0.0
    b_reconstruction_error = 0.0
    response_norm = 0.0
    maximum_numerator = 0.0
    minimum_denominator = np.inf
    skipped = 0
    boundary_pairs = 0
    pair_diagnostics: list[tuple[float, float, int, int, int, complex, complex]] = []
    normalization_error = 0.0
    eigenpair_error = 0.0
    response_density_error = 0.0
    response_energy_error = 0.0
    response_density_relative = 0.0
    response_energy_relative = 0.0

    for kpoint, path in sorted(paths.items()):
        dump = read_dump(path)
        reference = reference_like(dump)
        matrices: dict[tuple[str, str], np.ndarray] = {}
        for role, prefix in (("minus", args.minus), ("plus", args.plus)):
            overlap = read_basis_matrix(
                matrix_path(prefix, kpoint, False), b"EDMFTOV1"
            )
            hamiltonian = read_basis_matrix(
                matrix_path(prefix, kpoint, True), b"EDMFTHV1"
            )
            matrices[(role, "o")] = aligned_for_sternheimer(
                reference,
                overlap,
                role + " O",
                args.allow_lo_label_drift,
            )
            matrices[(role, "h")] = aligned_for_sternheimer(
                reference,
                hamiltonian,
                role + " H",
                args.allow_lo_label_drift,
            )

        do = (matrices[("plus", "o")] - matrices[("minus", "o")]) / (
            2.0 * args.strain
        )
        dh = (matrices[("plus", "h")] - matrices[("minus", "h")]) / (
            2.0 * args.strain
        )
        do_band = band_matrix(dump.coefficients, do)
        dh_band = band_matrix(dump.coefficients, dh)

        # Endpoint matrices do not include a central O/H file here.  Their
        # symmetric averages are sufficient to expose a mismatched basis or
        # an inconsistent strain step before using the derivative.
        o_mid = 0.5 * (matrices[("plus", "o")] + matrices[("minus", "o")])
        h_mid = 0.5 * (matrices[("plus", "h")] + matrices[("minus", "h")])
        normalization_error = max(
            normalization_error,
            float(
                np.max(
                    np.abs(
                        band_matrix(dump.coefficients, o_mid)
                        - np.eye(dump.nstored)
                    )
                )
            ),
        )
        eigenpair_error = max(
            eigenpair_error,
            float(
                np.max(
                    np.abs(
                        band_matrix(dump.coefficients, h_mid)
                        - np.diag(dump.energies)
                    )
                )
            ),
        )

        first = dump.dm_nemin - dump.nemin
        last = dump.dm_nemax - dump.nemin + 1
        outside = np.r_[0:first, last : dump.nstored]
        window = np.arange(first, last)
        response = np.zeros((window.size,) + dump.projector.shape[1:], complex)
        internal_response = np.zeros_like(response)
        diagonal_response = np.zeros_like(response)
        for local_n, n in enumerate(window):
            denominator = dump.energies[n] - dump.energies[outside]
            numerator = dh_band[outside, n] - dump.energies[n] * do_band[outside, n]
            valid = np.abs(denominator) > args.denominator_tol
            skipped += int(np.count_nonzero(~valid))
            boundary_pairs += denominator.size
            if np.any(valid):
                minimum_denominator = min(
                    minimum_denominator, float(np.min(np.abs(denominator[valid])))
                )
                maximum_numerator = max(
                    maximum_numerator, float(np.max(np.abs(numerator[valid])))
                )
                rotation = numerator[valid] / denominator[valid]
                if args.report_pairs:
                    valid_outside = outside[valid]
                    # This is the individual ordered-pair contribution to
                    # 2 Re sum(C*dU).  Keeping it beside the denominator
                    # reveals when a nominally tiny LAPW-kernel error is
                    # amplified by a band-window boundary degeneracy.
                    k_an = np.einsum(
                        "aij,ij->a",
                        dump.projector[valid_outside],
                        dump.gradient[local_n],
                        optimize=True,
                    )
                    contributions = 2.0 * np.real(k_an * rotation.conj())
                    pair_diagnostics.extend(
                        (
                            float(abs(den)),
                            float(contrib),
                            kpoint,
                            int(n + dump.nemin),
                            int(a + dump.nemin),
                            complex(num),
                            complex(kval),
                        )
                        for den, contrib, a, num, kval in zip(
                            denominator[valid],
                            contributions,
                            valid_outside,
                            numerator[valid],
                            k_an,
                        )
                    )
                response[local_n] = np.einsum(
                    "a,aij->ij", rotation.conj(), dump.projector[outside[valid]], optimize=True
                )
            inside = window[window != n]
            inside_denominator = dump.energies[n] - dump.energies[inside]
            inside_numerator = (
                dh_band[inside, n] - dump.energies[n] * do_band[inside, n]
            )
            inside_valid = np.abs(inside_denominator) > args.denominator_tol
            skipped += int(np.count_nonzero(~inside_valid))
            boundary_pairs += inside_denominator.size
            if np.any(inside_valid):
                minimum_denominator = min(
                    minimum_denominator,
                    float(np.min(np.abs(inside_denominator[inside_valid]))),
                )
                maximum_numerator = max(
                    maximum_numerator,
                    float(np.max(np.abs(inside_numerator[inside_valid]))),
                )
                internal_response[local_n] = np.einsum(
                    "a,aij->ij",
                    inside_numerator[inside_valid].conj()
                    / inside_denominator[inside_valid],
                    dump.projector[inside[inside_valid]],
                    optimize=True,
                )
            # Parallel-transport gauge for a generalized eigenproblem.  This
            # term is compulsory: omitting it differentiates eigenvectors as
            # if the LAPW coefficient metric were strain independent.
            diagonal_response[local_n] = (
                -0.5 * do_band[n, n] * dump.projector[n]
            )
        response_norm = max(response_norm, float(np.linalg.norm(response)))
        total += 2.0 * float(np.real(np.sum(dump.gradient * response)))
        # The same scalar has a second, completely independent form.  Build
        # the Hermitian response density and energy-moment matrices which are
        # passed to dmft2's LAPW stress kernels and contract them directly
        # with the finite-difference coefficient-space H_,q and O_,q.  If the
        # sum below agrees with ``total`` but dmft2's WINDOW_TOTAL does not,
        # the error is in the analytic LAPW kernel used for the response pass,
        # rather than in the Sternheimer/projector algebra.
        rho_response = np.zeros((dump.nstored, dump.nstored), complex)
        energy_response = np.zeros_like(rho_response)
        for local_n, n in enumerate(window):
            denominator = dump.energies[n] - dump.energies[outside]
            valid = np.abs(denominator) > args.denominator_tol
            if not np.any(valid):
                continue
            # K_an/(e_n-e_a) is recovered from the already formed dU by
            # constructing the band response directly from the same H/O
            # numerator.  This avoids relying on a second projector solve.
            numerator = (
                dh_band[outside[valid], n]
                - dump.energies[n] * do_band[outside[valid], n]
            )
            # K_an is the contraction of C_n with the outside projector U_a.
            k_an = np.einsum(
                "aij,nij->an",
                dump.projector[outside[valid]],
                dump.gradient[local_n : local_n + 1],
                optimize=True,
            )[:, 0]
            value = k_an / denominator[valid]
            rows = outside[valid]
            # U differentiates a bra: K[a,n] multiplies conj(R[a,n]).
            # Put K/denominator in (a,n), so Tr[rho H'] uses H'[n,a].
            rho_response[rows, n] += value
            rho_response[n, rows] += np.conj(value)
            energy_response[rows, n] += dump.energies[n] * value
            energy_response[n, rows] += np.conj(dump.energies[n] * value)
        response_hamiltonian_total += float(
            np.real(np.trace(rho_response @ dh_band))
        )
        response_overlap_total += float(
            np.real(np.trace(energy_response @ do_band))
        )
        if response_paths is not None:
            saved = read_response_matrix(response_paths[kpoint])
            expected_header = (
                dump.nstored,
                dump.nemin,
                dump.nemax,
                dump.dm_nemin,
                dump.dm_nemax,
            )
            found_header = (
                saved.nstored,
                saved.nemin,
                saved.nemax,
                saved.dm_nemin,
                saved.dm_nemax,
            )
            if found_header != expected_header:
                raise ValueError(
                    f"response dimensions differ at k={kpoint}: "
                    f"{found_header} != {expected_header}"
                )
            density_delta = saved.density - rho_response
            energy_delta = saved.energy - energy_response
            response_density_error = max(
                response_density_error, float(np.max(np.abs(density_delta)))
            )
            response_energy_error = max(
                response_energy_error, float(np.max(np.abs(energy_delta)))
            )
            response_density_relative = max(
                response_density_relative,
                float(np.linalg.norm(density_delta))
                / max(float(np.linalg.norm(rho_response)), np.finfo(float).tiny),
            )
            response_energy_relative = max(
                response_energy_relative,
                float(np.linalg.norm(energy_delta))
                / max(float(np.linalg.norm(energy_response)), np.finfo(float).tiny),
            )
        internal_total += 2.0 * float(
            np.real(np.sum(dump.gradient * internal_response))
        )
        diagonal_metric_total += 2.0 * float(
            np.real(np.sum(dump.gradient * diagonal_response))
        )
        all_response = response + internal_response + diagonal_response
        all_band_total += 2.0 * float(
            np.real(np.sum(dump.gradient * all_response))
        )

        if norm_paths is not None:
            norm = read_normalization(norm_paths[kpoint])
            expected = (
                dump.nstored,
                dump.dm_nemax - dump.dm_nemin + 1,
                dump.maxdim,
                dump.norbitals,
                dump.nemin,
                dump.nemax,
                dump.dm_nemin,
                dump.dm_nemax,
            )
            found = (
                norm.nstored,
                norm.nwindow,
                norm.maxdim,
                norm.norbitals,
                norm.nemin,
                norm.nemax,
                norm.dm_nemin,
                norm.dm_nemax,
            )
            if found != expected or not np.allclose(norm.energies, dump.energies):
                raise ValueError(f"normalization dimensions differ at k={kpoint}")
            raw_response = np.zeros_like(response)
            raw_internal_response = np.zeros_like(response)
            raw_diagonal_response = np.zeros_like(response)
            for local_n, n in enumerate(window):
                denominator = dump.energies[n] - dump.energies[outside]
                numerator = (
                    dh_band[outside, n]
                    - dump.energies[n] * do_band[outside, n]
                )
                valid = np.abs(denominator) > args.denominator_tol
                if np.any(valid):
                    raw_response[local_n] = np.einsum(
                        "a,aij->ij",
                        numerator[valid].conj() / denominator[valid],
                        norm.projector[outside[valid]],
                        optimize=True,
                    )
                inside = window[window != n]
                inside_denominator = dump.energies[n] - dump.energies[inside]
                inside_numerator = (
                    dh_band[inside, n]
                    - dump.energies[n] * do_band[inside, n]
                )
                inside_valid = np.abs(inside_denominator) > args.denominator_tol
                if np.any(inside_valid):
                    raw_internal_response[local_n] = np.einsum(
                        "a,aij->ij",
                        inside_numerator[inside_valid].conj()
                        / inside_denominator[inside_valid],
                        norm.projector[inside[inside_valid]],
                        optimize=True,
                    )
                raw_diagonal_response[local_n] = (
                    -0.5 * do_band[n, n] * norm.projector[n]
                )
            normalization_total += 2.0 * float(
                np.real(np.sum(norm.gradient * raw_response))
            )
            rho_norm = np.zeros((dump.nstored, dump.nstored), complex)
            energy_norm = np.zeros_like(rho_norm)
            for local_n, n in enumerate(window):
                denominator = dump.energies[n] - dump.energies[outside]
                valid = np.abs(denominator) > args.denominator_tol
                if not np.any(valid):
                    continue
                k_an = np.einsum(
                    "aij,nij->an",
                    norm.projector[outside[valid]],
                    norm.gradient[local_n : local_n + 1],
                    optimize=True,
                )[:, 0]
                value = k_an / denominator[valid]
                rows = outside[valid]
                rho_norm[rows, n] += value
                rho_norm[n, rows] += np.conj(value)
                energy_norm[rows, n] += dump.energies[n] * value
                energy_norm[n, rows] += np.conj(dump.energies[n] * value)
            normalization_hamiltonian_total += float(
                np.real(np.trace(rho_norm @ dh_band))
            )
            normalization_overlap_total += float(
                np.real(np.trace(energy_norm @ do_band))
            )
            normalization_internal_total += 2.0 * float(
                np.real(np.sum(norm.gradient * raw_internal_response))
            )
            normalization_diagonal_total += 2.0 * float(
                np.real(np.sum(norm.gradient * raw_diagonal_response))
            )
            normalization_all_total += 2.0 * float(
                np.real(
                    np.sum(
                        norm.gradient
                        * (
                            raw_response
                            + raw_internal_response
                            + raw_diagonal_response
                        )
                    )
                )
            )

        if moment_paths is not None:
            moment = read_one_particle(moment_paths[kpoint])
            expected = (
                dump.kpoint,
                dump.nbasis,
                dump.nlo,
                dump.ncomponents,
                dump.nemin,
                dump.dm_nemin,
                dump.dm_nemax,
            )
            found = (
                moment.kpoint,
                moment.nbasis,
                moment.nnlo,
                moment.ncomponents,
                moment.nemin,
                moment.dm_nemin,
                moment.dm_nemax,
            )
            if found != expected:
                raise ValueError(
                    f"one-particle dimensions differ at k={kpoint}: "
                    f"{found} != {expected}"
                )
            nfull = moment.nfull
            if not np.allclose(moment.energies, dump.energies[:nfull]):
                raise ValueError(f"one-particle energies differ at k={kpoint}")
            derivative_eigenvalues = np.real(
                np.diag(dh_band[:nfull, :nfull])
                - moment.energies * np.diag(do_band[:nfull, :nfull])
            )
            moment_first = dump.dm_nemin - moment.nemin
            moment_last = dump.dm_nemax - moment.nemin + 1
            eigenvalue_total += float(
                np.real(np.sum(np.diag(moment.rho) * derivative_eigenvalues))
            )
            eigenvalue_below_window += float(
                np.real(
                    np.sum(
                        np.diag(moment.rho)[:moment_first]
                        * derivative_eigenvalues[:moment_first]
                    )
                )
            )
            eigenvalue_window += float(
                np.real(
                    np.sum(
                        np.diag(moment.rho)[moment_first:moment_last]
                        * derivative_eigenvalues[moment_first:moment_last]
                    )
                )
            )
            fixed_pair_total += float(
                np.real(
                    np.trace(moment.rho @ dh_band[:nfull, :nfull])
                    - np.trace(moment.rhoe @ do_band[:nfull, :nfull])
                )
            )
            energy_average = 0.5 * (
                moment.energies[:, None] + moment.energies[None, :]
            )
            sigma_moment = moment.rhoe - energy_average * moment.rho
            sigma_overlap = float(
                np.real(np.trace(sigma_moment @ do_band[:nfull, :nfull]))
            )
            sigma_overlap_diagonal = float(
                np.real(
                    np.sum(
                        np.diag(sigma_moment)
                        * np.diag(do_band[:nfull, :nfull])
                    )
                )
            )
            sigma_overlap_total += sigma_overlap
            sigma_overlap_diagonal_total += sigma_overlap_diagonal
            sigma_overlap_offdiagonal_total += sigma_overlap - sigma_overlap_diagonal
            rho_window = moment.rho[
                moment_first:moment_last, moment_first:moment_last
            ]
            rhoe_window = moment.rhoe[
                moment_first:moment_last, moment_first:moment_last
            ]
            energy_window = moment.energies[moment_first:moment_last]
            b_target = rhoe_window - energy_window[:, None] * rho_window
            u_window = dump.projector[window].reshape((window.size, -1))
            c_window = dump.gradient.reshape((window.size, -1))
            b_reconstructed = u_window @ c_window.T
            b_reconstruction_error = max(
                b_reconstruction_error,
                float(
                    np.linalg.norm(b_target - b_reconstructed)
                    / max(np.linalg.norm(b_target), np.finfo(float).tiny)
                ),
            )
            rotation_window = np.zeros((window.size, window.size), complex)
            for local_n, n in enumerate(window):
                inside = window[window != n]
                local_inside = np.flatnonzero(window != n)
                denominator = dump.energies[n] - dump.energies[inside]
                numerator = (
                    dh_band[inside, n]
                    - dump.energies[n] * do_band[inside, n]
                )
                valid = np.abs(denominator) > args.denominator_tol
                rotation_window[local_inside[valid], local_n] = (
                    numerator[valid] / denominator[valid]
                )
                rotation_window[local_n, local_n] = -0.5 * do_band[n, n]
            b_rotation_total += 2.0 * float(
                np.real(np.sum(b_target * rotation_window.conj()))
            )

    print(f"k points                                      {len(paths)}")
    print(f"max midpoint |A^dag O A-I|                   {normalization_error:.12e}")
    print(f"max midpoint |A^dag H A-e|                   {eigenpair_error:.12e} Ry")
    print(f"cross-window band pairs                       {boundary_pairs}")
    print(f"skipped near-degenerate pairs                 {skipped}")
    print(f"smallest retained |e_n-e_a|                  {minimum_denominator:.12e} Ry")
    print(f"largest |<a|H_q-e_n O_q|n>|                  {maximum_numerator:.12e} Ry")
    print(f"max ||dU_outside|| per k                     {response_norm:.12e}")
    if response_paths is not None:
        print(
            "max dmft2 rho-response matrix error          "
            f"{response_density_error:.12e}  relative {response_density_relative:.12e}"
        )
        print(
            "max dmft2 rhoE-response matrix error         "
            f"{response_energy_error:.12e}  relative {response_energy_relative:.12e}"
        )
    print(f"central fixed-X projector response           {total:+.12e} Ry")
    print(
        "  response Tr[rho H_q]                      "
        f"{response_hamiltonian_total:+.12e} Ry"
    )
    print(
        "  response Tr[rhoE O_q]                     "
        f"{response_overlap_total:+.12e} Ry"
    )
    print(
        "  response H-O pair                         "
        f"{response_hamiltonian_total-response_overlap_total:+.12e} Ry"
    )
    print(f"central within-window response               {internal_total:+.12e} Ry")
    print(f"central diagonal-metric response             {diagonal_metric_total:+.12e} Ry")
    print(f"central all-stored-band fixed-X response     {all_band_total:+.12e} Ry")
    if norm_paths is not None:
        print(
            "central normalization response             "
            f"{normalization_total:+.12e} Ry"
        )
        print(
            "  normalization Tr[rho H_q]                 "
            f"{normalization_hamiltonian_total:+.12e} Ry"
        )
        print(
            "  normalization Tr[rhoE O_q]                "
            f"{normalization_overlap_total:+.12e} Ry"
        )
        print(
            "  normalization H-O pair                    "
            f"{normalization_hamiltonian_total-normalization_overlap_total:+.12e} Ry"
        )
        print(
            "central complete moving-window response    "
            f"{total+normalization_total:+.12e} Ry"
        )
        print(
            "central normalization within-window        "
            f"{normalization_internal_total:+.12e} Ry"
        )
        print(
            "central normalization diagonal-metric      "
            f"{normalization_diagonal_total:+.12e} Ry"
        )
        print(
            "central normalization all-stored-band      "
            f"{normalization_all_total:+.12e} Ry"
        )
        print(
            "central all-stored normalized response     "
            f"{all_band_total+normalization_all_total:+.12e} Ry"
        )
    if moment_paths is not None:
        # The fixed-coefficient H/O pair changes representation through the
        # connection *inside* the stored one-particle space.  The
        # complementary-band response and the response of the projector
        # normalization are separate explicit projector/subspace terms; they
        # must not be folded into this identity.
        band_representation = (
            eigenvalue_total + internal_total + diagonal_metric_total
        )
        print(f"central diagonal eigenvalue derivative       {eigenvalue_total:+.12e} Ry")
        print(f"  occupied bands below DMFT window          {eigenvalue_below_window:+.12e} Ry")
        print(f"  diagonal part inside DMFT window          {eigenvalue_window:+.12e} Ry")
        print(f"central moving-KS-band representation        {band_representation:+.12e} Ry")
        print(f"central fixed-coefficient H/O pair           {fixed_pair_total:+.12e} Ry")
        print(f"central direct B times band connection       {b_rotation_total:+.12e} Ry")
        print(f"central Tr[M_Sigma O_q]                      {sigma_overlap_total:+.12e} Ry")
        print(f"  diagonal M_Sigma overlap                  {sigma_overlap_diagonal_total:+.12e} Ry")
        print(f"  off-diagonal M_Sigma overlap              {sigma_overlap_offdiagonal_total:+.12e} Ry")
        print(f"max relative ||B-U C||                       {b_reconstruction_error:.12e}")
        print(
            "central representation residual            "
            f"{band_representation-fixed_pair_total:+.12e} Ry"
        )
    if args.report_pairs and pair_diagnostics:
        nreport = max(0, args.report_pairs)
        print("smallest cross-window denominators (|de|, contribution, k, n, a, numerator, C.U_a)")
        for item in sorted(pair_diagnostics, key=lambda x: x[0])[:nreport]:
            de, contribution, kpoint, n, a, numerator, kval = item
            print(
                f"  {de:.12e} {contribution:+.12e} {kpoint:5d} {n:5d} {a:5d} "
                f"({numerator.real:+.5e},{numerator.imag:+.5e}) "
                f"({kval.real:+.5e},{kval.imag:+.5e})"
            )
        print("largest cross-window pair contributions")
        for item in sorted(pair_diagnostics, key=lambda x: abs(x[1]), reverse=True)[:nreport]:
            de, contribution, kpoint, n, a, numerator, kval = item
            print(
                f"  {de:.12e} {contribution:+.12e} {kpoint:5d} {n:5d} {a:5d} "
                f"({numerator.real:+.5e},{numerator.imag:+.5e}) "
                f"({kval.real:+.5e},{kval.imag:+.5e})"
            )
    if args.volume is not None:
        conversion = 14710.5164057 / (3.0 * args.volume)
        print(f"hydrostatic mean                              {total*conversion:+.9f} GPa")
        print(
            "hydrostatic within-window mean                 "
            f"{internal_total*conversion:+.9f} GPa"
        )
        print(
            "hydrostatic diagonal-metric mean               "
            f"{diagonal_metric_total*conversion:+.9f} GPa"
        )
        if norm_paths is not None:
            print(
                "hydrostatic normalization mean                "
                f"{normalization_total*conversion:+.9f} GPa"
            )
            print(
                "hydrostatic complete moving-window mean       "
                f"{(total+normalization_total)*conversion:+.9f} GPa"
            )
            print(
                "hydrostatic all-stored normalized mean        "
                f"{(all_band_total+normalization_all_total)*conversion:+.9f} GPa"
            )
        if moment_paths is not None:
            band_representation = (
                eigenvalue_total + internal_total + diagonal_metric_total
            )
            print(
                "hydrostatic moving-KS representation mean     "
                f"{band_representation*conversion:+.9f} GPa"
            )
            print(
                "hydrostatic fixed-coefficient pair mean       "
                f"{fixed_pair_total*conversion:+.9f} GPa"
            )
            print(
                "hydrostatic representation residual mean      "
                f"{(band_representation-fixed_pair_total)*conversion:+.9f} GPa"
            )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        raise SystemExit(2)
