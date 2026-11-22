"""Read native WIEN2k LAPW vectors for projector-free DFT response.

The inventory pass reads only record headers, so MPI ranks may seek directly
to their assigned k points without duplicating all eigenvectors in memory.
"""

from pathlib import Path
from types import SimpleNamespace
import struct

import numpy as np


def _record(stream, *, skip=False):
    head = stream.read(4)
    if not head:
        return None
    if len(head) != 4:
        raise ValueError("truncated WIEN vector record header")
    n = struct.unpack("<i", head)[0]
    if n < 0 or n > 2**31 - 8:
        raise ValueError("invalid WIEN vector record length")
    if skip:
        stream.seek(n, 1)
        payload = b""
    else:
        payload = stream.read(n)
        if len(payload) != n:
            raise ValueError("truncated WIEN vector record")
    tail = stream.read(4)
    if len(tail) != 4 or struct.unpack("<i", tail)[0] != n:
        raise ValueError("mismatched WIEN vector record marker")
    return payload


def _header(data):
    if data is None:
        return None
    if len(data) != struct.calcsize("<ddd10siid"):
        raise ValueError("unsupported WIEN vector k header")
    kx, ky, kz, name, nbasis, nbands, weight = struct.unpack("<ddd10siid", data)
    if nbasis < 1 or nbands < 1 or not np.isfinite(weight) or weight <= 0:
        raise ValueError("invalid WIEN vector k point")
    return kx, ky, kz, name, nbasis, nbands, weight


def load_vectors(path, nat, nlo, parallel):
    path = Path(path)

    def inventory():
        result = []
        with path.open("rb") as stream:
            for _ in range(2 * nat):
                if _record(stream, skip=True) is None:
                    raise ValueError("incomplete WIEN vector radial header")
            while True:
                offset = stream.tell()
                header = _header(_record(stream))
                if header is None:
                    break
                if _record(stream, skip=True) is None:
                    raise ValueError("missing reciprocal labels")
                for _ in range(header[5]):
                    _record(stream, skip=True)
                    _record(stream, skip=True)
                result.append((offset, header[4], header[5], header[6]))
        if not result:
            raise ValueError("empty WIEN vector file")
        if len({r[2] for r in result}) != 1:
            raise ValueError("DFT response requires a fixed stored-band count")
        return result

    all_k = parallel.root_call(inventory)
    ownership = parallel.partition([float(nbasis) ** 3 for _, nbasis, _, _ in all_k])
    total_weight = sum(row[3] for row in all_k)
    records, weights = [], []
    with path.open("rb") as stream:
        for index in ownership:
            offset, nbasis, nbands, weight = all_k[index]
            stream.seek(offset)
            kx, ky, kz, name, n0, n, w = _header(_record(stream))
            if (n0, n, w) != (nbasis, nbands, weight):
                raise ValueError("WIEN vector changed during response inventory")
            labels = np.frombuffer(_record(stream), "<i4").copy()
            if labels.size != 3 * n0:
                raise ValueError("incorrect WIEN vector labels")
            labels = labels.reshape((3, n0), order="F")
            energies = np.empty(n)
            coefficients = np.empty((n0, n, 1), complex)
            for band in range(n):
                number, energy = struct.unpack("<id", _record(stream))
                if number != band + 1:
                    raise ValueError("out-of-order WIEN vector bands")
                data = _record(stream)
                if len(data) == 8 * n0:
                    col = np.frombuffer(data, "<f8")
                elif len(data) == 16 * n0:
                    col = np.frombuffer(data, "<c16")
                else:
                    raise ValueError("unsupported WIEN vector coefficient type")
                energies[band] = energy
                coefficients[:, band, 0] = col
            records.append(SimpleNamespace(
                kpoint=int(index + 1), nbasis=n0, nlo=nlo, nstored=n,
                energies=energies, coefficients=coefficients, labels=labels,
                reciprocal_k=np.array([kx, ky, kz]),
            ))
            weights.append(weight / total_weight)
    if abs(parallel.sum(sum(weights)) - 1) > 1e-10:
        raise ValueError("incomplete DFT k-point partition")
    return records, np.asarray(weights)
