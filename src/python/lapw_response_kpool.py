#!/usr/bin/env python3
"""Distribute an isolated WIEN LAPW1 export over independent k-point pools.

The unchanged serial WIEN worker assembles/diagonalizes each pool. Radial
searches and potential inputs are identical in all pools. Only the k list is
partitioned; original integration weights, names and global IDs are retained.

Each MPI rank owns a private directory. After successful workers, root merges
vector/energy streams after checking their radial headers and the shared NSH
file. Native H/O files have independent global names; only their integer IDs
are remapped. This avoids a serial eigensolve bottleneck without changing any
LAPW physics or requiring a distributed dense eigensolver for small matrices.
This adapter is for the response workspace, never the live material directory.
"""
from pathlib import Path
import argparse
import hashlib
import os
import shutil
import struct
import subprocess
import sys
import time
import json
from response_parallel import ResponseParallel
from dmft2_step import entries, local_definition
from wstruct import Struct


def klist_records(path):
    lines = Path(path).read_text().splitlines()
    end = next((i for i, line in enumerate(lines) if line[:10].strip() == "END"), None)
    if end is None or end == 0 or any(line.strip() for line in lines[end + 1 :]):
        raise ValueError("incomplete or unsupported k list")
    lines = lines[:end]
    # Match INILPW's new-format probe, then its legacy four-I5 fallback.
    width = None
    for n in (10, 5):
        try:
            [int(lines[0][10 + j * n : 10 + (j + 1) * n]) for j in range(4)]
            float(lines[0][10 + 4 * n : 15 + 4 * n])
            width = n
            break
        except ValueError:
            pass
    if width is None or any(line.startswith(" kp ") for line in lines):
        raise ValueError("unsupported k-list layout")
    weight_end = 15 + 4 * width
    # Energy limits are read from the FIRST point only by INILPW. Preserve
    # these same limits when another original point starts a worker's list.
    suffix = lines[0][weight_end:]
    return lines, weight_end, suffix


def read_record(stream, endian="<"):
    """Bounded header reader for the native four-byte sequential record format."""
    marker = stream.read(4)
    if len(marker) != 4:
        raise ValueError("truncated vector radial header")
    size = struct.unpack(endian + "i", marker)[0]
    if not 0 <= size <= 1024**2:
        raise ValueError("unsupported vector radial record length")
    payload = stream.read(size)
    trailer = stream.read(4)
    if len(payload) != size or trailer != marker:
        raise ValueError("incomplete vector radial record")
    return marker + payload + trailer


def vector_header(stream, nat):
    marker = stream.read(4)
    stream.seek(0)
    if len(marker) != 4:
        raise ValueError("empty vector file")
    endian = "<" if 0 < struct.unpack("<i", marker)[0] <= 1024**2 else ">"
    return b"".join(read_record(stream, endian) for unused in range(2 * nat))


def merge_streams(paths, output, nat, *, binary):
    """Retain identical radial headers once; append ordered k records exactly."""
    header = None
    temporary = Path(str(output) + ".merging")
    try:
        with temporary.open("wb") as out:
            for path in paths:
                with Path(path).open("rb") as stream:
                    current = (
                        vector_header(stream, nat)
                        if binary
                        else b"".join(stream.readline() for unused in range(2 * nat))
                    )
                    if not current or (not binary and current.count(b"\n") != 2 * nat):
                        raise ValueError("missing native energy header")
                    if header is None:
                        header = current
                        out.write(header)
                    elif current != header:
                        raise ValueError("k pools have different radial energies")
                    if not stream.peek(1):
                        raise ValueError("empty native k-point payload")
                    shutil.copyfileobj(stream, out, 4 * 1024**2)
        temporary.replace(output)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def remap_matrix(path, output, global_k, rank):
    """Preserve matrix payload exactly; update the export's integer provenance."""
    with Path(path).open("r+b") as stream:
        magic = stream.read(8)
        marker = stream.read(4)
        if magic not in (b"EDMFTHV1", b"EDMFTOV1"):
            raise ValueError("bad LAPW1 matrix export")
        endian = next(
            (
                sign
                for sign, order in [("<", "little"), (">", "big")]
                if int.from_bytes(marker, order) == 0x01020304
            ),
            None,
        )
        if endian is None:
            raise ValueError("bad matrix byte order")
        header = stream.read(32)
        if len(header) != 32:
            raise ValueError("truncated matrix header")
        version, kpoint, owner, n, nlo, storage, triangle, convention = struct.unpack(
            endian + "8i", header
        )
        if (
            (version, storage, triangle, convention) != (1, 2, 1, 1)
            or not 0 <= nlo <= n
            or n < 1
        ):
            raise ValueError("unsupported matrix header")
        stream.seek(16)
        stream.write(struct.pack(endian + "ii", global_k, rank))
    Path(path).replace(output)


def run(workspace, case, worker, parallel):
    start = time.perf_counter()
    folder = Path(workspace).resolve()
    worker = Path(worker).resolve()
    definition = local_definition((folder / "lapw1.def").read_text())
    defs = entries(definition)
    for unit in (4, 5, 6, 10, 11, 18, 19, 20, 21, 71):
        if unit not in defs:
            raise ValueError(f"missing LAPW1 unit {unit}")
    s = Struct()
    s.ReadStruct(str(folder / defs[20]))
    lines, width, suffix = parallel.root_call(lambda: klist_records(folder / defs[4]))
    if parallel.size > len(lines):
        raise ValueError("LAPW1 MPI ranks must not exceed k points")
    # Contiguous pools make vector concatenation preserve the global k order.
    first = len(lines) * parallel.rank // parallel.size
    last = len(lines) * (parallel.rank + 1) // parallel.size
    poolroot = folder / "lapw1_pools"
    parallel.root_call(lambda: poolroot.mkdir(exist_ok=False))
    local = poolroot / f"rank{parallel.rank:05d}"
    local.mkdir()
    for unit in (5, 18, 19, 20):
        shutil.copy2(folder / defs[unit], local / defs[unit])
    (local / "lapw1.def").write_text(definition)
    subset = list(lines[first:last])
    subset[0] = subset[0][:width].ljust(width) + suffix
    (local / defs[4]).write_text("\n".join(subset) + "\nEND\n")
    env = os.environ | {
        "WIEN_LAPW1_MATRIX_DUMP": str(local / "native"),
        "OMP_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
        "VECLIB_MAXIMUM_THREADS": "1",
    }
    before = time.perf_counter()
    log = local / "worker.log"
    with log.open("w") as stream:
        result = subprocess.run(
            [str(worker), "lapw1.def"],
            cwd=local,
            env=env,
            stdout=stream,
            stderr=subprocess.STDOUT,
        )
    text = log.read_text(errors="replace")
    if (
        result.returncode
        or "LAPW1 END" not in text
        or any(
            word in text
            for word in ["Fortran runtime error", "ERROR STOP", "Error termination"]
        )
    ):
        raise RuntimeError(f"LAPW1 pool {parallel.rank} failed; see {log}")
    nsh = hashlib.sha256((local / defs[71]).read_bytes()).hexdigest()
    for prefix in ("native.k", "native.ham.k"):
        if len(list(local.glob(prefix + "????????.bin"))) != last - first:
            raise ValueError("incomplete LAPW1 pool exports")
        for j, k in enumerate(range(first, last), 1):
            remap_matrix(
                local / f"{prefix}{j:08d}.bin",
                folder / f"{prefix}{k+1:08d}.bin",
                k + 1,
                parallel.rank,
            )
    stats = parallel.gather(
        dict(
            rank=parallel.rank,
            first=first + 1,
            last=last,
            seconds=time.perf_counter() - before,
            nsh_sha256=nsh,
        )
    )

    def merge():
        if len({x["nsh_sha256"] for x in stats}) != 1:
            raise ValueError("k pools have different nonspherical radial matrices")
        pools = [poolroot / f"rank{i:05d}" for i in range(parallel.size)]
        merge_streams(
            [p / defs[10] for p in pools], folder / defs[10], s.nat, binary=True
        )
        merge_streams(
            [p / defs[11] for p in pools], folder / defs[11], s.nat, binary=False
        )
        # Radial searches/NSH are k independent. Full individual output logs
        # remain in the pools; never concatenate duplicate :E radial searches.
        for unit in (6, 21, 71):
            shutil.copy2(pools[0] / defs[unit], folder / defs[unit])
        report = dict(
            ranks=parallel.size,
            kpoints=len(lines),
            pools=stats,
            seconds=time.perf_counter() - start,
            scope="isolated WIEN LAPW1 k pools",
        )
        (folder / "lapw1_pools.json").write_text(json.dumps(report, indent=2) + "\n")
        return report

    report = parallel.root_call(merge)
    if parallel.root:
        print("LAPW1 END: verified and merged MPI k pools", flush=True)
    return report


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--workspace", required=True)
    p.add_argument("--case", required=True)
    p.add_argument("--worker", required=True)
    args = p.parse_args()
    parallel = ResponseParallel.world()
    try:
        run(args.workspace, args.case, args.worker, parallel)
    except BaseException:
        import traceback

        traceback.print_exc()
        parallel.comm.Abort(1)
