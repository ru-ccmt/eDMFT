#!/usr/bin/env python3
"""Select independently controlled stress-response contributions.

This is a report/validation layer, not a new stress formula and not an SCF
file editor.  Its initial scope is the hydrostatic MEAN tensile stress in
GPa; the result must never be used as a general corrected tensor.

Records share a SHA-256 reference_id describing the physical inputs.  A
correction additionally records its component and validation status.
The radial policy is off / diagnostic / on.  Diagnostic displays the radial
candidate and both totals, but does not include it in the selected total.
On requires a correction explicitly validated for production.  Missing or
stale corrections are errors, never silently interpreted as zero.
"""
from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path

CONVENTION = "hydrostatic_mean_tensile_GPa"
MODES = {"off", "diagnostic", "on"}
MATCH_FIELDS = (
    "reference_id",
    "convention",
    "volume_Bohr3",
    "ensemble",
    "radial_energy_policy",
    "potential_coordinates",
)


def number(value, role):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{role} must be a finite number")
    if not math.isfinite(value):
        raise ValueError(f"nonfinite {role}")
    return float(value)


def validate_reference(record):
    if record.get("schema") != 1:
        raise ValueError("unsupported stress-response schema")
    if not re.fullmatch(r"[0-9a-f]{64}", record.get("reference_id", "")):
        raise ValueError("reference_id must be a SHA-256 physics-input fingerprint")
    if record.get("convention") != CONVENTION:
        raise ValueError("only hydrostatic mean tensile stress in GPa is supported")
    if number(record.get("volume_Bohr3"), "volume") <= 0:
        raise ValueError("volume must be positive")
    if record.get("ensemble") not in {"fixed_N", "fixed_mu"}:
        raise ValueError("the thermodynamic ensemble must be explicit")
    for key in ("radial_energy_policy", "potential_coordinates"):
        if not isinstance(record.get(key), str) or not record[key].strip():
            raise ValueError(f"missing {key}")


def correction_value(record, baseline, component, mode):
    if mode not in MODES:
        raise ValueError(f"unknown {component} mode: {mode}")
    if record is None:
        if mode != "off":
            raise ValueError(f"{component} mode {mode} requires a correction record")
        return None
    validate_reference(record)
    if record.get("component") != component:
        raise ValueError(f"expected component {component}")
    for key in MATCH_FIELDS:
        if record[key] != baseline[key]:
            raise ValueError(f"incompatible {component} {key}")
    if component in baseline["includes"]:
        raise ValueError(f"baseline already includes {component}; would double count")
    value = number(record.get("mean_tensile_correction_GPa"), component)
    validation = record.get("validation", {})
    if not isinstance(validation, dict):
        raise ValueError("validation must be an object")
    if mode == "on":
        if validation.get("production_ready") is not True:
            raise ValueError(f"{component} is not validated for production inclusion")
        checks = validation.get("checks")
        if (
            not isinstance(checks, dict)
            or not checks
            or any(v is not True for v in checks.values())
        ):
            raise ValueError(f"{component} has missing or failed validation checks")
    return value


def assemble(
    baseline, projector=None, radial=None, *, projector_mode="off", radial_mode="off"
):
    """Pure, repeatable composition; no in-place changes to input records.

    The fingerprint and conventions are checked even for supplied off-mode
    records, because their values may be shown in diagnostic comparisons.
    Fixed-mu is a reporting convention here, not a claim that the response
    solver currently supports that ensemble.
    """
    validate_reference(baseline)
    included = baseline.get("includes")
    if not isinstance(included, list) or any(not isinstance(x, str) for x in included):
        raise ValueError("baseline must declare its included components")
    if len(included) != len(set(included)):
        raise ValueError("baseline has duplicate included components")
    if "radial" in included or "projector" in included:
        raise ValueError(
            "supply a baseline before response corrections to permit independent toggles"
        )
    base = number(baseline.get("mean_tensile_GPa"), "baseline stress")
    p = correction_value(projector, baseline, "projector", projector_mode)
    r = correction_value(radial, baseline, "radial", radial_mode)
    selected_without_radial = base + (p if projector_mode == "on" else 0.0)
    selected = selected_without_radial + (r if radial_mode == "on" else 0.0)
    # Candidate totals are deliberately separate from the selected total.
    candidate_without_radial = base + (p if p is not None else 0.0)
    candidate_with_radial = None if r is None else candidate_without_radial + r
    return {
        "schema": 1,
        **{key: baseline[key] for key in MATCH_FIELDS},
        "modes": {"projector": projector_mode, "radial": radial_mode},
        "baseline_mean_tensile_GPa": base,
        "projector_correction_GPa": p,
        "radial_correction_GPa": r,
        "candidate_without_radial_GPa": candidate_without_radial,
        "candidate_with_radial_GPa": candidate_with_radial,
        "selected_without_radial_GPa": selected_without_radial,
        "selected_mean_tensile_GPa": selected,
        "selected_pressure_GPa": -selected,
        "applied_components": [
            name
            for name, mode in (("projector", projector_mode), ("radial", radial_mode))
            if mode == "on"
        ],
        "validation": {
            "projector": None if projector is None else projector.get("validation"),
            "radial": None if radial is None else radial.get("validation"),
        },
        "scope": "report_only; no SCF or tensor modification",
    }


def read(path):
    if path is None:
        return None
    with path.open() as stream:
        result = json.load(stream)
    if not isinstance(result, dict):
        raise ValueError(f"expected a JSON object: {path}")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--projector", type=Path)
    parser.add_argument("--radial", type=Path)
    parser.add_argument("--projector-mode", choices=sorted(MODES), default="off")
    parser.add_argument("--radial-mode", choices=sorted(MODES), default="off")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        input_paths = [
            p.resolve() for p in (args.baseline, args.projector, args.radial) if p
        ]
        if args.output.resolve() in input_paths:
            raise ValueError("output must not overwrite an input record")
        result = assemble(
            read(args.baseline),
            read(args.projector),
            read(args.radial),
            projector_mode=args.projector_mode,
            radial_mode=args.radial_mode,
        )
        # Create-only output makes reruns explicit, and prevents using a
        # corrected report as the baseline of a second correction pass.
        with args.output.open("x") as stream:
            json.dump(result, stream, indent=2, allow_nan=False)
            stream.write("\n")
    except (ValueError, OSError) as error:
        parser.error(str(error))
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
