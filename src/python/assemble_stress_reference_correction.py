#!/usr/bin/env python3
"""Assemble a guarded reference correction to the dmft2 valence stress.

The ordinary analytic dmft2 stress already contains an approximate
``M_Sigma*dO`` contribution and three pieces of the normalized-projector
response.  The matrix/projector finite-difference diagnostics provide
reference values for the same two *complete* objects.  This program writes
only their difference,

    delta T = (M_full - M_current) + (U_joint - U_current),

as an extensive 3x3 tensor in Ry.  dmft2 applies that tensor only when the
environment variable ``EDMFT_STRESS_REFERENCE_CORRECTION`` names the output
file.  The production default is therefore unchanged.

At present the diagnostic finite differences are hydrostatic.  Their trace
is distributed isotropically over the three diagonal tensor elements.  This
is sufficient for the cubic/hydrostatic closure tests, but it is deliberately
labelled a reference path and must not be used as a general tensor correction.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import re
import sys

import numpy as np


RY_BOHR3_TO_GPA = 14710.5164057


def tagged_value(text: str, tag: str) -> float:
    match = re.search(rf"(?m)^{re.escape(tag)}\s+([+-]?[0-9.Ee+-]+)\s*$", text)
    if match is None:
        raise ValueError(f"missing required value {tag}")
    return float(match.group(1))


def labelled_value(text: str, label: str) -> float:
    match = re.search(
        rf"(?m)^{re.escape(label)}\s+([+-]?[0-9.Ee+-]+)(?:\s+Ry|\s+GPa)?\s*$",
        text,
    )
    if match is None:
        raise ValueError(f"missing required value {label!r}")
    return float(match.group(1))


def relative_asymmetry(symmetric: float, difference: float) -> float:
    return abs(difference) / max(abs(symmetric), np.finfo(float).tiny)


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outputdmf2", required=True, type=Path)
    parser.add_argument("--metric-report", required=True, type=Path)
    parser.add_argument("--projector-report", type=Path)
    parser.add_argument("--volume", required=True, type=float)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--max-metric-asymmetry",
        type=float,
        default=0.25,
        help="maximum |right-left|/|central| accepted for M_Sigma (default 0.25)",
    )
    parser.add_argument(
        "--max-projector-asymmetry",
        type=float,
        default=0.25,
        help="maximum |right-left|/|central| accepted for U (default 0.25)",
    )
    parser.add_argument(
        "--max-projector-change",
        type=float,
        default=0.25,
        help="maximum endpoint ||delta U||/||U|| accepted (default 0.25)",
    )
    parser.add_argument(
        "--allow-unstable",
        action="store_true",
        help="write the diagnostic tensor even when a one-sided gate fails",
    )
    parser.add_argument("--analytic-total-gpa", type=float)
    parser.add_argument("--numerical-total-gpa", type=float)
    args = parser.parse_args(argv)
    for name in (
        "volume",
        "max_metric_asymmetry",
        "max_projector_asymmetry",
        "max_projector_change",
    ):
        value = getattr(args, name)
        if not np.isfinite(value) or value <= 0:
            parser.error(f"--{name.replace('_', '-')} must be positive")
    return args


def main(argv: list[str] | None = None) -> int:
    args = parse_args(sys.argv[1:] if argv is None else argv)
    output_text = args.outputdmf2.read_text(errors="replace")
    if ":DMFT_STRESS_SIGMA_DIAG_ENABLED T" not in output_text:
        raise ValueError(
            "the central dmft2 run did not print the isolated M_Sigma tensor; "
            "rerun it with EDMFT_STRESS_DETAILED_DIAGNOSTICS=1"
        )

    analytic_m = tagged_value(
        output_text, ":DMFT_STRESS_DIAG_SIGMA_MOMENT_OVERLAP_MEAN_GPA"
    )
    analytic_projector = sum(
        tagged_value(output_text, tag)
        for tag in (
            ":DMFT_STRESS_DIAG_WINDOW_TOTAL_MEAN_GPA",
            ":DMFT_STRESS_DIAG_NORM_WINDOW_TOTAL_MEAN_GPA",
            ":DMFT_STRESS_DIAG_EXPLICIT_PROJECTOR_MEAN_GPA",
        )
    )

    metric_text = args.metric_report.read_text(errors="replace")
    reference_m = labelled_value(
        metric_text, "hydrostatic mean of M_Sigma overlap"
    )
    metric_ry = labelled_value(metric_text, "-Tr[M_Sigma O_,strain]")
    metric_difference_ry = labelled_value(
        metric_text, "one-sided M_Sigma difference"
    )
    metric_asymmetry = relative_asymmetry(metric_ry, metric_difference_ry)

    reference_projector = analytic_projector
    projector_asymmetry = 0.0
    projector_change = 0.0
    if args.projector_report is not None:
        projector_text = args.projector_report.read_text(errors="replace")
        transport_match = re.search(
            r"(?m)^transport gauge\s+(\S+)\s*$", projector_text
        )
        if transport_match is None or transport_match.group(1) != "projector":
            raise ValueError(
                "the joint-projector report must use the explicitly labelled "
                "projector transport gauge"
            )
        reference_projector = labelled_value(projector_text, "hydrostatic mean")
        projector_ry = labelled_value(
            projector_text, "joint normalized-projector derivative"
        )
        projector_difference_ry = labelled_value(
            projector_text, "one-sided difference"
        )
        projector_asymmetry = relative_asymmetry(
            projector_ry, projector_difference_ry
        )
        projector_change = max(
            labelled_value(
                projector_text, "max relative projector change, minus"
            ),
            labelled_value(
                projector_text, "max relative projector change, plus"
            ),
        )

    failures: list[str] = []
    if metric_asymmetry > args.max_metric_asymmetry:
        failures.append(
            f"M_Sigma one-sided asymmetry {metric_asymmetry:.3g} exceeds "
            f"{args.max_metric_asymmetry:.3g}"
        )
    if projector_asymmetry > args.max_projector_asymmetry:
        failures.append(
            f"projector one-sided asymmetry {projector_asymmetry:.3g} exceeds "
            f"{args.max_projector_asymmetry:.3g}"
        )
    if projector_change > args.max_projector_change:
        failures.append(
            f"relative endpoint projector change {projector_change:.3g} exceeds "
            f"{args.max_projector_change:.3g}"
        )
    if failures and not args.allow_unstable:
        raise ValueError("; ".join(failures))

    metric_delta = reference_m - analytic_m
    projector_delta = reference_projector - analytic_projector
    total_delta_gpa = metric_delta + projector_delta
    # A hydrostatic strain derivative is the trace of the extensive tensor.
    # mean_GPa = trace(T)*conversion/(3*volume), hence an isotropic tensor
    # with each diagonal equal to mean_GPa*volume/conversion has that mean.
    diagonal_raw_ry = total_delta_gpa * args.volume / RY_BOHR3_TO_GPA
    correction = np.eye(3) * diagonal_raw_ry

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w") as handle:
        handle.write("EDMFT_STRESS_REFERENCE_CORRECTION_V1\n")
        handle.write("RAW_RY\n")
        for row in correction:
            handle.write(" ".join(f"{value:+.17e}" for value in row) + "\n")

    print(f"analytic M_Sigma overlap                    {analytic_m:+.9f} GPa")
    print(f"reference full-metric M_Sigma              {reference_m:+.9f} GPa")
    print(f"M_Sigma replacement delta                  {metric_delta:+.9f} GPa")
    print(f"M_Sigma one-sided asymmetry                {metric_asymmetry:.6e}")
    print(f"analytic complete projector response       {analytic_projector:+.9f} GPa")
    print(f"reference joint projector response         {reference_projector:+.9f} GPa")
    print(f"projector replacement delta                {projector_delta:+.9f} GPa")
    print(f"projector one-sided asymmetry              {projector_asymmetry:.6e}")
    print(f"maximum relative endpoint projector change {projector_change:.6e}")
    print(f"total additive correction                  {total_delta_gpa:+.9f} GPa")
    print(f"raw diagonal correction                    {diagonal_raw_ry:+.12e} Ry")
    if failures:
        print("WARNING: unstable reference accepted by explicit request:")
        for failure in failures:
            print(f"  {failure}")
    if args.analytic_total_gpa is not None:
        corrected = args.analytic_total_gpa + total_delta_gpa
        print(f"corrected analytic total                   {corrected:+.9f} GPa")
        if args.numerical_total_gpa is not None:
            print(
                "corrected analytic minus numerical          "
                f"{corrected-args.numerical_total_gpa:+.9f} GPa"
            )
    print(f"wrote                                        {args.output}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        raise SystemExit(2)
