#!/usr/bin/env python3
"""Merge WIEN2k SCF component files before running mixer.

This is the standalone equivalent of the ``scf()`` and ``scf1()`` steps in
``run_dft.py``.  By default, one newly assembled iteration is appended to
``case.scf`` so that its convergence history is retained.  The current charge
density and potential files are then copied to their ``*_old`` names, because
WIEN2k mixer reads ``case.clmsum_old`` (or the spin-resolved equivalents) and
writes the newly mixed density to ``case.clmsum``.  Use ``--replace`` to create
``case.scf`` from only the current component files.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import shutil
import sys
import tempfile


def nonempty(path: Path) -> bool:
    """Return True only for an existing, nonempty regular file."""
    return path.is_file() and path.stat().st_size > 0


def detect_case(requested: str | None) -> str:
    """Find the unique top-level WIEN2k case having struct and clmsum files."""
    if requested:
        case = requested
        if not Path(f"{case}.struct").is_file():
            raise RuntimeError(f"{case}.struct does not exist")
        if not Path(f"{case}.clmsum").is_file():
            raise RuntimeError(f"{case}.clmsum does not exist")
        return case

    candidates = []
    for struct in sorted(Path.cwd().glob("*.struct")):
        case = struct.name[: -len(".struct")]
        if Path(f"{case}.clmsum").is_file():
            candidates.append(case)

    if not candidates:
        raise RuntimeError(
            "Could not find a case with both case.struct and case.clmsum"
        )
    if len(candidates) > 1:
        names = ", ".join(candidates)
        raise RuntimeError(
            f"Ambiguous WIEN2k case ({names}); select one with -f CASE"
        )
    return candidates[0]


def detect_spin_polarized(case: str) -> bool:
    """Detect the spin-polarized SCF layout from paired component files."""
    return (
        nonempty(Path(f"{case}.scf1up"))
        and nonempty(Path(f"{case}.scf1dn"))
    ) or (
        nonempty(Path(f"{case}.scf2up"))
        and nonempty(Path(f"{case}.scf2dn"))
    )


def component_suffixes(spin_polarized: bool, spin_orbit: bool) -> list[str]:
    """Return the component order used by WIEN2k's run scripts."""
    if not spin_polarized:
        suffixes = ["0", "1"]
        if spin_orbit:
            suffixes.append("so")
        suffixes.extend(["2", "1s", "2s", "c"])
        return suffixes

    suffixes = ["0", "orbup", "orbdn", "orbud", "1up", "1dn"]
    if spin_orbit:
        suffixes.append("so")
    suffixes.extend(["2up", "2dn", "dmup"])
    # This follows runsp_lapw: with SO, dmdn duplicates information in dmup.
    if not spin_orbit:
        suffixes.append("dmdn")
    suffixes.extend(
        ["1sup", "1sdn", "2sup", "2sdn", "cup", "cdn"]
    )
    return suffixes


def required_suffixes(spin_polarized: bool, spin_orbit: bool) -> list[str]:
    if spin_polarized:
        required = ["0", "1up", "1dn", "2up", "2dn"]
    else:
        required = ["0", "1", "2"]
    if spin_orbit:
        required.append("so")
    return required


def selected_components(
    case: str, spin_polarized: bool, spin_orbit: bool
) -> list[Path]:
    missing = [
        Path(f"{case}.scf{suffix}")
        for suffix in required_suffixes(spin_polarized, spin_orbit)
        if not nonempty(Path(f"{case}.scf{suffix}"))
    ]
    if missing:
        raise RuntimeError(
            "Required SCF component files are missing or empty: "
            + ", ".join(str(path) for path in missing)
        )

    components = [
        Path(f"{case}.scf{suffix}")
        for suffix in component_suffixes(spin_polarized, spin_orbit)
        if nonempty(Path(f"{case}.scf{suffix}"))
    ]

    # run_dft.py places the impurity-orbital energy record after scf2 and
    # before the semicore/core records.  Insert it at the same boundary.
    eorb = Path("Eorb_imp.dat")
    if nonempty(eorb):
        first_tail = next(
            (
                i
                for i, path in enumerate(components)
                if path.name.endswith(("scf1s", "scf2s", "scfc", "scf1sup", "scf1sdn", "scf2sup", "scf2sdn", "scfcup", "scfcdn"))
            ),
            len(components),
        )
        components.insert(first_tail, eorb)

    return components


def label(path: Path, case: str) -> str:
    if path.name == "Eorb_imp.dat":
        return "Eorb_imp.dat"
    return path.name if path.name.startswith(case + ".") else str(path)


def merge(case: str, components: list[Path], replace: bool) -> Path:
    """Atomically append or replace case.scf with the selected components."""
    destination = Path(f"{case}.scf")
    old_mode = destination.stat().st_mode if destination.exists() else None
    from stress_response_publish import response_rows, mark_published
    current_iteration = '\n'.join(p.read_text() for p in components)
    correction = response_rows(case, current_iteration)

    fd, temporary_name = tempfile.mkstemp(
        prefix=f".{case}.scf.merge.", dir=Path.cwd()
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(fd, "wb") as output:
            if not replace and destination.exists():
                with destination.open("rb") as previous:
                    shutil.copyfileobj(previous, output)
                if destination.stat().st_size:
                    with destination.open("rb") as previous:
                        previous.seek(-1, os.SEEK_END)
                        if previous.read(1) != b"\n":
                            output.write(b"\n")

            for component in components:
                header = f"\n------- {label(component, case)} --------\n"
                output.write(header.encode("utf-8"))
                with component.open("rb") as source:
                    shutil.copyfileobj(source, output)
                if component.stat().st_size:
                    with component.open("rb") as source:
                        source.seek(-1, os.SEEK_END)
                        if source.read(1) != b"\n":
                            output.write(b"\n")

            output.write(correction.encode('utf-8'))
            output.flush()
            os.fsync(output.fileno())

        if old_mode is not None:
            os.chmod(temporary, old_mode)
        os.replace(temporary, destination)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise

    mark_published(case, correction, destination, current_iteration=current_iteration)
    return destination


def backup_mixer_inputs(case: str) -> list[tuple[Path, Path]]:
    """Advance the charge/potential history expected by WIEN2k mixer.

    This deliberately mirrors ``scf1(case)`` in ``run_dft.py``.  Omitting
    this step leaves ``case.clmsum_old`` frozen at an earlier iteration, so
    repeated manual ``dmft2 -> merge_scf.py -> mixer`` cycles keep mixing
    every new density against that stale reference instead of converging.
    ``copy2`` preserves WIEN2k's text/binary file contents and metadata.
    """
    suffixes = (
        "clmsum",
        "vsp",
        "vns",
        "vrespsum",
        "vspdn",
        "vnsdn",
        "vspup",
        "vnsup",
        "clmup",
        "clmdn",
    )
    copied: list[tuple[Path, Path]] = []
    for suffix in suffixes:
        source = Path(f"{case}.{suffix}")
        if not nonempty(source):
            continue
        destination = Path(f"{source}_old")
        shutil.copy2(source, destination)
        copied.append((source, destination))
    return copied


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Merge current WIEN2k SCF components before mixer"
    )
    parser.add_argument(
        "-f", "--case", help="WIEN2k case name; normally detected automatically"
    )
    parser.add_argument(
        "--replace",
        action="store_true",
        help="replace case.scf instead of appending a new iteration",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="show the detected case and input order without writing case.scf",
    )
    return parser.parse_args()


def main() -> int:
    options = parse_arguments()
    try:
        case = detect_case(options.case)
        spin_polarized = detect_spin_polarized(case)
        spin_orbit = nonempty(Path(f"{case}.inso"))
        components = selected_components(case, spin_polarized, spin_orbit)
    except RuntimeError as error:
        print(f"merge_scf.py: error: {error}", file=sys.stderr)
        return 2

    print(f"case={case}")
    print(f"spin-polarized={'yes' if spin_polarized else 'no'}")
    print(f"spin-orbit={'yes' if spin_orbit else 'no'}")
    print("SCF components:")
    for component in components:
        print(f"  {component}")

    if options.dry_run:
        print("Dry run: no files changed.")
        return 0

    destination = merge(case, components, options.replace)
    backups = backup_mixer_inputs(case)
    action = "Replaced" if options.replace else "Appended one iteration to"
    print(f"{action} {destination}")
    print("Prepared mixer history:")
    for source, backup in backups:
        print(f"  {source} -> {backup}")
    print(f"Now run: x mixer -f {case}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
