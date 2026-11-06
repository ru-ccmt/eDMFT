#!/usr/bin/env python3
"""Numerically compare two extended WIEN2k case.nsh files.

The files are sequential-unformatted Fortran streams.  A bytewise comparison
is unnecessarily strict when LAPW1 implementations use slightly different
radial integrators, while a text conversion can silently lose the integer
angular labels or one of the appended stress blocks.  This utility checks the
record structure, all six labels in every main record, and every complex value
in the V, W, phase-V, phase-W, and local-orbital records.
"""

from __future__ import annotations

import argparse
import math
import struct
import sys
from pathlib import Path


def read_fortran_records(path: Path) -> tuple[str, list[bytes]]:
    data = path.read_bytes()
    for endian in ("<", ">"):
        records: list[bytes] = []
        offset = 0
        try:
            while offset < len(data):
                if offset + 4 > len(data):
                    raise ValueError("truncated leading record marker")
                size = struct.unpack_from(endian + "i", data, offset)[0]
                if size < 0 or offset + 8 + size > len(data):
                    raise ValueError("invalid record size")
                start = offset + 4
                stop = start + size
                tail = struct.unpack_from(endian + "i", data, stop)[0]
                if tail != size:
                    raise ValueError("mismatched record markers")
                records.append(data[start:stop])
                offset = stop + 4
            return endian, records
        except (ValueError, struct.error):
            pass
    raise ValueError(f"{path}: not a 4-byte-marker sequential Fortran file")


def complex_error(endian: str, left: bytes, right: bytes) -> tuple[float, float]:
    if len(left) % 16:
        raise ValueError("complex record length is not divisible by 16")
    values_left = struct.unpack(endian + f"{len(left) // 8}d", left)
    values_right = struct.unpack(endian + f"{len(right) // 8}d", right)
    maximum = 0.0
    scale = 0.0
    for index in range(0, len(values_left), 2):
        zleft = complex(values_left[index], values_left[index + 1])
        zright = complex(values_right[index], values_right[index + 1])
        maximum = max(maximum, abs(zleft - zright))
        scale = max(scale, abs(zleft), abs(zright))
    return maximum, scale


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--atol", type=float, default=1.0e-10)
    parser.add_argument("--rtol", type=float, default=1.0e-8)
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("--max-failures", type=int, default=20)
    args = parser.parse_args()

    endian_ref, reference = read_fortran_records(args.reference)
    endian_new, candidate = read_fortran_records(args.candidate)
    if endian_ref != endian_new:
        print("FAIL: files use different byte ordering")
        return 1
    if len(reference) != len(candidate):
        print(f"FAIL: record counts differ: {len(reference)} != {len(candidate)}")
        return 1

    failures = 0
    numeric_records = 0
    global_error = 0.0
    global_scale = 0.0
    block_names = ("V", "W", "phase-V", "phase-W")
    block_error = [0.0] * 4
    block_scale = [0.0] * 4
    for number, (left, right) in enumerate(zip(reference, candidate), 1):
        if len(left) != len(right):
            print(f"record {number}: size differs: {len(left)} != {len(right)}")
            failures += 1
            continue

        # A main record is six 32-bit integers followed by sixteen complex
        # values.  The all-zero terminator has the identical 280-byte shape.
        if len(left) == 6 * 4 + 16 * 16:
            labels_left = struct.unpack_from(endian_ref + "6i", left)
            labels_right = struct.unpack_from(endian_ref + "6i", right)
            if labels_left != labels_right:
                print(f"record {number}: labels differ: {labels_left} != {labels_right}")
                failures += 1
            error, scale = complex_error(endian_ref, left[24:], right[24:])
            kind = "main"
        elif len(left) > 0 and len(left) % 16 == 0:
            error, scale = complex_error(endian_ref, left, right)
            kind = "LO"
        else:
            # Header strings and the one-integer atom marker are required to
            # match exactly.  No floating data live in either record type.
            if left != right:
                print(f"record {number}: nonnumeric payload differs")
                failures += 1
            continue

        numeric_records += 1
        global_error = max(global_error, error)
        global_scale = max(global_scale, scale)
        # Every numeric record is four equally sized blocks in the extended
        # format, including each of the five optional LO records.
        numeric_left = left[24:] if kind == "main" else left
        numeric_right = right[24:] if kind == "main" else right
        block_bytes = len(numeric_left) // 4
        for block in range(4):
            first = block * block_bytes
            last = first + block_bytes
            error_block, scale_block = complex_error(
                endian_ref,
                numeric_left[first:last],
                numeric_right[first:last],
            )
            block_error[block] = max(block_error[block], error_block)
            block_scale[block] = max(block_scale[block], scale_block)
        tolerance = args.atol + args.rtol * scale
        if error > tolerance:
            if failures < args.max_failures:
                print(
                    f"record {number:6d} {kind:4s}: max |delta|={error:.6e}, "
                    f"scale={scale:.6e}, tolerance={tolerance:.6e}"
                )
            failures += 1
        elif args.verbose and error != 0.0:
            print(
                f"record {number:6d} {kind:4s}: max |delta|={error:.6e}, "
                f"scale={scale:.6e}"
            )

    relative = global_error / global_scale if global_scale else 0.0
    print(f"records checked       {len(reference)}")
    print(f"numeric records       {numeric_records}")
    print(f"maximum |delta|       {global_error:.12e}")
    print(f"maximum relative      {relative:.12e}")
    for name, error, scale in zip(block_names, block_error, block_scale):
        relative_block = error / scale if scale else 0.0
        print(
            f"{name:12s} max delta {error:.12e}  "
            f"relative-to-max {relative_block:.12e}"
        )
    print(f"failed records        {failures}")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
