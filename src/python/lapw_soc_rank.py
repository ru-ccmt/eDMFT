"""Preserve the parent's scalar and spinor band spaces in a private SOC replay.

The scalar rank passed to LAPWSO and the final spinor output rank can both
vary with k. A buffered LAPW1 input alone does not specify which vectors
were retained by a parent workflow. The native vector headers do.
"""

import hashlib
from pathlib import Path
import struct
from tempfile import TemporaryDirectory
from lapw_dft_vector import _record, _header


def parent_lo_filter(scfso):
    """Files predating the marker used component clipping in LAPWSO."""
    path = Path(scfso)
    policy = "LEGACY_COMPONENT_CUTOFF"
    if path.is_file():
        with path.open(errors="replace") as stream:
            for line in stream:
                if line.startswith(":SOC_LO_MATCHING:"):
                    policy = line.split(":", 2)[2].strip()
    if policy not in ("LINEAR", "LEGACY_COMPONENT_CUTOFF"):
        raise ValueError(f"unknown local-orbital matching policy in {path}: {policy!r}")
    return policy == "LEGACY_COMPONENT_CUTOFF"


def k_key(k):
    return ",".join(str(round(x * 1e9)) for x in k)


def vector_files(name):
    path = Path(name)
    if not path.name.endswith("_x"):
        return [path]
    stem = path.name[:-1]
    paths = sorted(
        (f for f in path.parent.glob(stem + "*") if f.name[len(stem) :].isdigit()),
        key=lambda f: int(f.name[len(stem) :]),
    )
    if not paths:
        raise ValueError(f"missing MPI vector files for {path}")
    return paths


def ranks(name, nat):
    rows = {}
    for path in vector_files(name):
        with path.open("rb") as stream:
            for _ in range(2 * nat):
                if _record(stream, skip=True) is None:
                    raise ValueError(f"truncated radial headers in {path}")
            while (header := _header(_record(stream))) is not None:
                k = k_key(header[:3])
                _, _, _, _, basis, count, weight = header
                labels = _record(stream)
                if k in rows or len(labels) != 12 * basis:
                    raise ValueError(
                        f"duplicate k point or invalid reciprocal labels in {path}"
                    )
                rows[k] = dict(
                    basis=basis,
                    count=count,
                    weight=weight,
                    labels=hashlib.sha256(labels).hexdigest(),
                )
                for _ in range(2 * count):
                    if _record(stream, skip=True) is None:
                        raise ValueError(f"truncated band records in {path}")
    if not rows:
        raise ValueError(f"empty parent mesh at {name}")
    return rows


def parent_ranks(scalar, up, down, nat):
    s, u, d = [ranks(name, nat) for name in (scalar, up, down)]
    if u != d or set(s) != set(u):
        raise ValueError("central scalar and two-component vector meshes differ")
    for k in s:
        if any(s[k][key] != u[k][key] for key in ("basis", "labels", "weight")):
            raise ValueError(f"scalar and spinor LAPW rows differ at k={k}")
    return dict(scalar=s, spinor=u)


def _put(out, data):
    marker = struct.pack("<i", len(data))
    out.write(marker)
    out.write(data)
    out.write(marker)


def enforce_ranks(vector, energy, nat, expected, norm=None, *, check_cut=False):
    """Trim only the private regenerated pool, including energy and norm files."""
    vector, energy = Path(vector), Path(energy)
    actual = ranks(vector, nat)
    if actual.keys() != expected.keys():
        raise ValueError("central and replay k meshes differ")
    for k, record in actual.items():
        target = expected[k]
        if any(record[key] != target[key] for key in ("basis", "labels", "weight")):
            raise ValueError(f"LAPW basis changed during central replay at k={k}")
        if record["count"] < target["count"]:
            raise ValueError(f"regenerated band pool is too small at k={k}")
    if actual == expected:
        return
    with TemporaryDirectory(dir=vector.parent, prefix=".soc-ranks-") as temp:
        tv, te, tn = [Path(temp) / n for n in ("vector", "energy", "norm")]
        counts = []
        with vector.open("rb") as inp, energy.open() as ein, tv.open(
            "wb"
        ) as out, te.open("w") as eout:
            for _ in range(2 * nat):
                _put(out, _record(inp))
                eout.write(ein.readline())
            while (hdr := _header(_record(inp))) is not None:
                k = k_key(hdr[:3])
                keep = expected[k]["count"]
                old = hdr[5]
                line = ein.readline()
                if (
                    k_key([float(line[j : j + 19]) for j in (0, 19, 38)]) != k
                    or int(line[73:79]) != old
                ):
                    raise ValueError(f"energy and vector headers differ at k={k}")
                _put(out, struct.pack("<ddd10siid", *hdr[:5], keep, hdr[6]))
                _put(out, _record(inp))
                eout.write(line[:73] + f"{keep:6d}" + line[79:])
                counts.append((old, keep))
                for j in range(1, old + 1):
                    row = ein.readline()
                    eigen = _record(inp)
                    index, ev = struct.unpack("<id", eigen)
                    fields = row.split()
                    if (
                        index != j
                        or int(fields[0]) != j
                        or abs(float(fields[1]) - ev) > 1e-8
                    ):
                        raise ValueError(f"energy rows inconsistent at k={k}, band {j}")
                    if check_cut:
                        if j == keep:
                            edge = ev
                        elif j == keep + 1 and ev - edge <= 1e-10:
                            raise ValueError(
                                f"signed reference cuts a degenerate band pair at k={k}"
                            )
                    if j <= keep:
                        _put(out, eigen)
                        _put(out, _record(inp))
                        eout.write(row)
                    else:
                        _record(inp, skip=True)
            if ein.read().strip():
                raise ValueError("extra energy records after last k point")
        if norm is not None:
            norm = Path(norm)
            with norm.open() as ins, tn.open("w") as outs:
                for old, keep in counts:
                    for i in range(0, old, 4):
                        line = ins.readline().rstrip("\n")
                        if len(line.rstrip()) != 20 * min(old - i, 4):
                            raise ValueError(f"invalid spinor norm record in {norm}")
                        if i < keep:
                            outs.write(line[: 20 * min(keep - i, 4)] + "\n")
                if ins.read().strip():
                    raise ValueError("extra spinor norm records")
        tv.replace(vector)
        te.replace(energy)
        if norm is not None:
            tn.replace(norm)
