"""Read versioned native LAPW H/O snapshots without diagnostic-script imports.

The existing LAPW1/dftKS exporter stores an endian marker, integer header,
fractional k, native PW/LO labels, and the complex packed upper triangle.
No row reordering or extra conjugation is applied. Matching native/vector
coefficient gauges is checked by the valence preparation layer.
"""

from dataclasses import dataclass
from pathlib import Path
import numpy as np


@dataclass(frozen=True)
class BasisMatrix:
    kpoint: int
    nbasis: int
    nlo: int
    reciprocal_k: np.ndarray
    labels: np.ndarray
    matrix: np.ndarray


def vector_k_from_matrix_header(k, structure):
    """Convert the EDMFTOV1/HV1 k header to the coordinates of its G labels.

    The matrix dump records LAPW1's internal K3. TAPEWF writes the same
    G labels but transforms K3 for a centered monoclinic CXZ cell. DMFT2
    constructs k+G from that vector header. Apply exactly TAPEWF's mapping
    and LATGEN's angular tolerance when reconstructing the basis or projector.
    The version-1 export is already in widespread use.
    """
    if (structure.lattyp.strip() == "CXZ"
            and abs(np.deg2rad(structure.gamma - 90.0)) >= 1.e-4):
        x, y, z = k
        return np.array([x + z, y, z - x])
    return np.asarray(k).copy()


def read_basis_matrix(path, magic):
    """Read EDMFTHV1/EDMFTOV1; reject incomplete, unsupported or nonfinite data."""
    if magic not in (b"EDMFTHV1", b"EDMFTOV1"):
        raise ValueError("unsupported native matrix kind")
    path = Path(path)
    with path.open("rb") as handle:
        if handle.read(8) != magic:
            raise ValueError("wrong native matrix magic: " + str(path))
        marker = handle.read(4)
        endian = next(
            (
                sign
                for sign, order in [("<", "little"), (">", "big")]
                if int.from_bytes(marker, order) == 0x01020304
            ),
            None,
        )
        if endian is None:
            raise ValueError("invalid native matrix byte order")

        def read(dtype, count):
            data = np.fromfile(handle, endian + dtype, count)
            if len(data) != count:
                raise ValueError("truncated native matrix snapshot")
            return data

        header = read("i4", 8)
        version, kpoint, rank, n, nlo, storage, triangle, convention = map(int, header)
        if (version, storage, triangle, convention) != (1, 2, 1, 1):
            raise ValueError("unsupported native matrix encoding")
        if kpoint < 1 or rank < 0 or n < 1 or not 0 <= nlo <= n:
            raise ValueError("invalid native matrix dimensions")
        count = n * (n + 1) // 2
        expected = 8 + 4 + 8 * 4 + 3 * 8 + 3 * n * 4 + count * 16
        if path.stat().st_size != expected:
            raise ValueError("native matrix file length does not match header")
        k = read("f8", 3).astype(float)
        labels = read("i4", 3 * n).reshape((3, n), order="F").astype(int)
        packed = read("c16", count)
        if not np.all(np.isfinite(k)) or not np.all(np.isfinite(packed)):
            raise ValueError("nonfinite native matrix snapshot")
    matrix = np.empty((n, n), complex)
    offset = 0
    for column in range(n):
        values = packed[offset : offset + column + 1]
        matrix[: column + 1, column] = values
        matrix[column, :column] = values[:column].conj()
        offset += column + 1
    return BasisMatrix(kpoint, n, nlo, k, labels, matrix)
