"""Hermitian RR/QR algebra without an unused QQ block.

A column array stores [Y_RR; Y_QR]. Its physical Hermitian completion has
Y_RQ=Y_QR^dagger and Y_QQ=0. Matrix covectors use the real Frobenius pairing
Re Tr(Y_RR^dagger dH_RR)+2 Re Tr(Y_QR^dagger dH_QR). Raw projector arrays do
NOT use this factor of two. These helpers also accept ordinary square matrices.
"""

import numpy as np
from lapw_density_response import checked, herm


def checked_columns(value, size, name):
    value = np.asarray(value)
    if value.ndim != 2 or value.shape[0] != size or not 0 < value.shape[1] <= size:
        raise ValueError(f"invalid {name} columns")
    checked(value, value.shape, name)
    checked(
        value[: value.shape[1]], (value.shape[1],) * 2, name + " RR", hermitian=True
    )
    return value


def hermitian_columns(value):
    """Hermitize only RR; the rectangular QR block has an independent meaning."""
    value = np.array(value, complex, copy=True)
    if value.ndim != 2 or not 0 < value.shape[1] <= value.shape[0]:
        raise ValueError("invalid Hermitian columns")
    n = value.shape[1]
    value[:n] = herm(value[:n])
    return value


def right_multiply_columns(left, columns):
    """left @ HermitianCompletion(columns), without forming that completion."""
    left = np.asarray(left)
    if left.ndim < 1:
        raise ValueError("invalid left factor")
    y = checked_columns(columns, left.shape[-1], "matrix covector")
    n = y.shape[1]
    # Combine all angular/radial leading indices into one BLAS row dimension.
    # A batched (..., nprimitive, nband) matmul dispatches many tiny products;
    # their flattening is one identical GEMM. Also compute the R columns as
    # left @ [Y_RR;Y_QR], avoiding a separate sum of two matrix products.
    flat = left.reshape(-1, left.shape[-1])
    if n == len(y):
        return (flat @ y).reshape(left.shape)
    result = np.concatenate((flat @ y, flat[:, :n] @ y[n:].conj().T), axis=-1)
    return result.reshape(left.shape)


def pair_columns(left, right):
    """Physical pairing, including both QR and RQ legs."""
    a = np.asarray(left)
    if a.ndim != 2:
        raise ValueError("invalid column pairing")
    checked_columns(a, len(a), "left columns")
    b = checked(right, a.shape, "right columns")
    n = a.shape[1]
    return float(np.vdot(a[:n], b[:n]).real + 2 * np.vdot(a[n:], b[n:]).real)
