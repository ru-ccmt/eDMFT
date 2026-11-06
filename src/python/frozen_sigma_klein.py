#!/usr/bin/env python3
"""Evaluate the Klein functional for frozen-self-energy strain tests.

The ordinary charge-self-consistent eDMFT free energy is assembled as

    F = ENE - SUM + YSUM - EORB + XEORB,

where WIEN2k's mixer supplies ENE, dmft2 supplies SUM/YSUM/VSUM, and the
impurity solver supplies EORB/XEORB.  At a strained geometry for which the
impurity self-energy and double counting are frozen at their reference values,
XEORB must be replaced by its tangent plane.  Because

    VSUM - YSUM = Tr[(Sigma_imp - Vdc) G_local],

the strain-dependent part of that tangent functional is simply

    K_frozen = ENE - SUM + VSUM - EORB.

This program performs that bookkeeping after the DFT charge has converged.  It
is deliberately read-only: the tangent correction is a scalar needed for the
finite-difference comparison and must not be inserted into mixer.

Typical use, from the unstrained reference directory, is

    frozen_sigma_klein.py
    frozen_sigma_klein.py --minus ../eps_m --plus ../eps_p --strain 0.001

The second command also checks that sig.inp and Eorb_imp.dat are identical in
the reference and strained directories.  The reported derivative uses the
convention +(1/Omega) dGamma/d(strain).

The analogous tangent-functional correction for frozen-Sigma phonon
calculations is related, but is intentionally not implemented by this strain
postprocessor.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
import math
from pathlib import Path
import re
import sys
from typing import Iterable


# 1 Ry/Bohr^3 in GPa, using the 2018 CODATA exact/conventional values used by
# current scientific Python packages.  Keeping the primary result in
# Ry/Bohr^3 avoids making this conversion part of the validation itself.
RY_PER_BOHR3_TO_GPA = 14710.507641078

FLOAT_PATTERN = re.compile(
    r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?"
)


class FunctionalInputError(RuntimeError):
    """Raised when a required WIEN2k/eDMFT output is missing or ambiguous."""


@dataclass(frozen=True)
class FunctionalTerms:
    """Scalar terms, in Ry, entering the implemented Klein functional."""

    directory: str
    case: str
    ene: float
    band_sum: float
    vsum: float
    ysum: float
    eorb: float
    xeorb: float
    volume: float | None
    scf2_files: tuple[str, ...]

    @property
    def sigma_g(self) -> float:
        """Tr[(Sigma_imp-Vdc) G_local] = VSUM-YSUM, in Ry."""

        return self.vsum - self.ysum

    @property
    def on_shell_free_energy(self) -> float:
        """Ordinary fully stationary eDMFT free energy, in Ry."""

        return (
            self.ene
            - self.band_sum
            + self.ysum
            - self.eorb
            + self.xeorb
        )

    @property
    def frozen_reduced(self) -> float:
        """Strain-dependent part of the frozen-Sigma tangent functional."""

        return self.ene - self.band_sum + self.vsum - self.eorb

    def frozen_aligned(self, reference: "FunctionalTerms") -> float:
        """Tangent functional aligned to the physical reference free energy.

        The added constant XEORB_0-Tr[(Sigma_0-Vdc_0)Glocal_0] makes the
        tangent functional equal the ordinary on-shell functional at the
        reference geometry.  This constant cancels from every derivative.
        """

        return self.frozen_reduced + reference.xeorb - reference.sigma_g


def nonempty(path: Path) -> bool:
    return path.is_file() and path.stat().st_size > 0


def detect_case(directory: Path, requested: str | None) -> str:
    """Detect a unique WIEN2k case in *directory*."""

    if requested:
        struct = directory / f"{requested}.struct"
        if not struct.is_file():
            raise FunctionalInputError(f"missing {struct}")
        return requested

    candidates: list[str] = []
    for struct in sorted(directory.glob("*.struct")):
        case = struct.name[: -len(".struct")]
        has_mixer = nonempty(directory / f"{case}.scfm") or nonempty(
            directory / f"{case}.scf"
        )
        has_dmft2 = any(
            nonempty(directory / name)
            for name in (f"{case}.scf2", f"{case}.scf2up")
        )
        if has_mixer and has_dmft2:
            candidates.append(case)

    if not candidates:
        raise FunctionalInputError(
            f"could not detect a case in {directory}: need case.struct, "
            "case.scfm (or case.scf), and case.scf2"
        )
    if len(candidates) > 1:
        raise FunctionalInputError(
            f"ambiguous cases in {directory}: {', '.join(candidates)}; use -f CASE"
        )
    return candidates[0]


def first_float(text: str, context: str) -> float:
    match = FLOAT_PATTERN.search(text)
    if match is None:
        raise FunctionalInputError(f"could not parse a number from {context}: {text!r}")
    value = float(match.group(0).replace("D", "E").replace("d", "e"))
    if not math.isfinite(value):
        raise FunctionalInputError(f"non-finite value in {context}: {value}")
    return value


def last_tag_value(path: Path, tag: str) -> float:
    """Return the numeric value from the last line beginning with *tag*."""

    if not nonempty(path):
        raise FunctionalInputError(f"missing or empty {path}")

    value: float | None = None
    with path.open("r", errors="replace") as stream:
        for line_number, line in enumerate(stream, 1):
            stripped = line.lstrip()
            if not stripped.startswith(tag):
                continue
            # ENE/SUM/VSUM/VOL place the desired scalar after '='.  EORB and
            # XEORB place it immediately after the tag.
            tail = stripped[len(tag) :]
            if "=" in tail:
                tail = tail.rsplit("=", 1)[1]
            value = first_float(tail, f"{path}:{line_number} ({tag})")

    if value is None:
        raise FunctionalInputError(f"tag {tag!r} was not found in {path}")
    return value


def mixer_file(directory: Path, case: str) -> Path:
    scfm = directory / f"{case}.scfm"
    if nonempty(scfm):
        return scfm
    scf = directory / f"{case}.scf"
    if nonempty(scf):
        return scf
    raise FunctionalInputError(f"missing {case}.scfm and {case}.scf in {directory}")


def dmft2_files(directory: Path, case: str) -> tuple[Path, ...]:
    """Return scf2 files using the averaging convention of run_dmft.py."""

    # There are two magnetic naming conventions in existing calculations.
    # A WIEN2k-style spin-polarized run can use scf2up/scf2dn, whereas the
    # antiferromagnetic eDMFT workflow driven by run_dmft.py keeps the first
    # complete-cell mapping in scf2 and the opposite mapping in scf2dn.
    # Both files describe the full electron count, so their extensive scalar
    # terms must be averaged, exactly as scf_merge.py averages SUM/VSUM and
    # the valence-stress tensor.
    up = directory / f"{case}.scf2up"
    default = directory / f"{case}.scf2"
    down = directory / f"{case}.scf2dn"
    if nonempty(up) and nonempty(down):
        return (up, down)

    if nonempty(default) and nonempty(down):
        return (default, down)

    if nonempty(default):
        return (default,)
    raise FunctionalInputError(
        f"missing {case}.scf2, or a magnetic pair {case}.scf2/{case}.scf2dn "
        f"or {case}.scf2up/{case}.scf2dn, in {directory}"
    )


def average_tag(paths: Iterable[Path], tag: str) -> float:
    values = [last_tag_value(path, tag) for path in paths]
    return sum(values) / len(values)


def optional_volume(directory: Path, case: str) -> float | None:
    """Read the last :VOL from case.scf, falling back to case.scfm."""

    for name in (f"{case}.scf", f"{case}.scfm"):
        path = directory / name
        if not nonempty(path):
            continue
        try:
            return last_tag_value(path, ":VOL")
        except FunctionalInputError:
            pass
    return None


def load_terms(directory: Path, requested_case: str | None) -> FunctionalTerms:
    directory = directory.expanduser().resolve()
    if not directory.is_dir():
        raise FunctionalInputError(f"not a directory: {directory}")

    case = detect_case(directory, requested_case)
    scfm = mixer_file(directory, case)
    scf2 = dmft2_files(directory, case)
    eorb_file = directory / "Eorb_imp.dat"

    return FunctionalTerms(
        directory=str(directory),
        case=case,
        ene=last_tag_value(scfm, ":ENE"),
        band_sum=average_tag(scf2, ":SUM"),
        vsum=average_tag(scf2, ":VSUM"),
        ysum=average_tag(scf2, ":YSUM"),
        eorb=last_tag_value(eorb_file, ":EORB"),
        xeorb=last_tag_value(eorb_file, ":XEORB"),
        volume=optional_volume(directory, case),
        scf2_files=tuple(path.name for path in scf2),
    )


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def frozen_file_checks(directories: list[Path]) -> list[dict[str, object]]:
    """Check that the frozen impurity inputs are byte-identical."""

    checks: list[dict[str, object]] = []
    for filename in ("sig.inp", "Eorb_imp.dat"):
        records = []
        for directory in directories:
            path = directory.expanduser().resolve() / filename
            records.append(
                {
                    "directory": str(directory.expanduser().resolve()),
                    "exists": path.is_file(),
                    "sha256": sha256(path) if path.is_file() else None,
                }
            )
        hashes = {record["sha256"] for record in records if record["sha256"]}
        all_exist = all(bool(record["exists"]) for record in records)
        checks.append(
            {
                "file": filename,
                "identical": all_exist and len(hashes) == 1,
                "records": records,
            }
        )
    return checks


def terms_dict(terms: FunctionalTerms, reference: FunctionalTerms) -> dict[str, object]:
    result = asdict(terms)
    result.update(
        {
            "sigma_g": terms.sigma_g,
            "on_shell_free_energy": terms.on_shell_free_energy,
            "frozen_reduced": terms.frozen_reduced,
            "frozen_aligned": terms.frozen_aligned(reference),
        }
    )
    return result


def print_terms(label: str, terms: FunctionalTerms, reference: FunctionalTerms) -> None:
    print(f"{label}: {terms.directory}  (case={terms.case})")
    print(f"  ENE                         {terms.ene: .12f} Ry")
    print(f"  SUM                         {terms.band_sum: .12f} Ry")
    print(f"  VSUM                        {terms.vsum: .12f} Ry")
    print(f"  YSUM                        {terms.ysum: .12f} Ry")
    print(f"  EORB                        {terms.eorb: .12f} Ry")
    print(f"  XEORB                       {terms.xeorb: .12f} Ry")
    print(f"  VSUM-YSUM = Tr[Sigma~ Gloc] {terms.sigma_g: .12f} Ry")
    print(f"  ordinary on-shell F         {terms.on_shell_free_energy: .12f} Ry")
    print(f"  frozen-Sigma reduced K      {terms.frozen_reduced: .12f} Ry")
    print(f"  frozen-Sigma aligned Gamma  {terms.frozen_aligned(reference): .12f} Ry")
    if terms.volume is not None:
        print(f"  volume                      {terms.volume: .12f} Bohr^3")
    print()


def print_checks(checks: list[dict[str, object]]) -> None:
    print("Frozen-input checks:")
    for check in checks:
        state = "identical" if check["identical"] else "WARNING: not identical"
        print(f"  {check['file']}: {state}")
        if not check["identical"]:
            for record in check["records"]:  # type: ignore[index]
                digest = record["sha256"] or "missing"
                print(f"    {record['directory']}: {digest}")
    print()


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Postprocess the stationary or frozen-Sigma DFT+eDMFT Klein functional"
        ),
        epilog=(
            "Example from the reference directory: %(prog)s --minus ../eps_m "
            "--plus ../eps_p --strain 0.001"
        ),
    )
    parser.add_argument(
        "directory",
        nargs="?",
        default=".",
        help="calculation directory; also the default reference directory",
    )
    parser.add_argument("-f", "--case", help="WIEN2k case name (normally detected)")
    parser.add_argument(
        "--reference",
        help="unstrained reference directory (default: DIRECTORY)",
    )
    parser.add_argument("--minus", help="directory for negative strain")
    parser.add_argument("--plus", help="directory for positive strain")
    parser.add_argument(
        "--strain",
        type=float,
        help="positive magnitude of the strain parameter used for +/- directories",
    )
    parser.add_argument(
        "--volume",
        type=float,
        help="reference volume in Bohr^3 (default: last :VOL in reference case.scf)",
    )
    parser.add_argument(
        "--no-check-frozen",
        action="store_true",
        help="do not compare sig.inp and Eorb_imp.dat across directories",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="write machine-readable JSON instead of the explanatory report",
    )
    return parser.parse_args()


def main() -> int:
    options = parse_arguments()
    pair_requested = options.minus is not None or options.plus is not None
    if pair_requested and (options.minus is None or options.plus is None):
        print(
            "frozen_sigma_klein.py: error: --minus and --plus must be used together",
            file=sys.stderr,
        )
        return 2
    if pair_requested and (options.strain is None or options.strain <= 0.0):
        print(
            "frozen_sigma_klein.py: error: positive --strain is required with "
            "--minus/--plus",
            file=sys.stderr,
        )
        return 2

    reference_directory = Path(options.reference or options.directory)
    try:
        reference = load_terms(reference_directory, options.case)
        if pair_requested:
            minus = load_terms(Path(options.minus), options.case)
            plus = load_terms(Path(options.plus), options.case)
        else:
            current = load_terms(Path(options.directory), options.case)
    except FunctionalInputError as error:
        print(f"frozen_sigma_klein.py: error: {error}", file=sys.stderr)
        return 2

    if not pair_requested:
        result = {
            "reference": terms_dict(reference, reference),
            "calculation": terms_dict(current, reference),
            "delta_frozen_aligned": current.frozen_aligned(reference)
            - reference.frozen_aligned(reference),
        }
        if options.json:
            print(json.dumps(result, indent=2, sort_keys=True))
        else:
            print("Frozen-Sigma Klein-functional postprocessor")
            print("  F_on-shell = ENE-SUM+YSUM-EORB+XEORB")
            print("  K_frozen   = ENE-SUM+VSUM-EORB  (strain-dependent part)\n")
            if current.directory != reference.directory:
                print_terms("reference", reference, reference)
            print_terms("calculation", current, reference)
            if current.directory != reference.directory:
                delta = result["delta_frozen_aligned"]
                print(f"Delta frozen-Sigma Gamma = {delta: .12e} Ry")
        return 0

    directories = [reference_directory, Path(options.minus), Path(options.plus)]
    checks = [] if options.no_check_frozen else frozen_file_checks(directories)
    delta = plus.frozen_reduced - minus.frozen_reduced
    derivative = delta / (2.0 * options.strain)
    volume = options.volume if options.volume is not None else reference.volume
    stress_ry = derivative / volume if volume is not None else None
    stress_gpa = stress_ry * RY_PER_BOHR3_TO_GPA if stress_ry is not None else None

    result = {
        "reference": terms_dict(reference, reference),
        "minus": terms_dict(minus, reference),
        "plus": terms_dict(plus, reference),
        "strain_magnitude": options.strain,
        "delta_k_plus_minus": delta,
        "d_k_d_strain": derivative,
        "reference_volume_bohr3": volume,
        "stress_derivative_ry_per_bohr3": stress_ry,
        "stress_derivative_gpa": stress_gpa,
        "stress_sign_convention": "+(1/Omega) dGamma/d(strain)",
        "frozen_input_checks": checks,
    }

    if options.json:
        print(json.dumps(result, indent=2, sort_keys=True))
        return 0

    print("Frozen-Sigma Klein-functional central difference")
    print("  K_frozen = ENE-SUM+VSUM-EORB")
    print("  sign convention: +(1/Omega) dGamma/d(strain)\n")
    print_terms("reference", reference, reference)
    print_terms("minus", minus, reference)
    print_terms("plus", plus, reference)
    if checks:
        print_checks(checks)
    print(f"strain magnitude                  {options.strain: .12e}")
    print(f"K(+)-K(-)                        {delta: .12e} Ry")
    print(f"dK/d(strain)                     {derivative: .12e} Ry")
    if volume is None:
        print("reference volume                  unavailable")
        print("Use --volume BOHR3 to obtain the stress derivative.")
    else:
        print(f"reference volume                  {volume: .12f} Bohr^3")
        print(f"+(1/Omega)dK/d(strain)            {stress_ry: .12e} Ry/Bohr^3")
        print(f"                                    {stress_gpa: .9f} GPa")
    if any(not bool(check["identical"]) for check in checks):
        print(
            "\nWARNING: frozen impurity files differ.  The central difference is not "
            "a controlled frozen-Sigma test.",
            file=sys.stderr,
        )
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
