"""Signed Matsubara contractions for the streamed frozen-self-energy adjoint.

For a Hermitian perturbation the dynamic derivative contains products of G,
G U Σ and Σ U† G at the same frequency. A real orthogonal rotation of the
frequency index preserves every such bilinear sum, including products without
a complex conjugate. Rotate positive and negative quadrature weights separately.
The static subtraction, fixed-N constraint and normalized-projector derivative
are handled by the electronic model outside this local k-point contraction.
"""
from tempfile import TemporaryFile
import errno

import numpy as np
from scipy.linalg import svd
from scipy.linalg.blas import zgemm


class FrequencyModes:
    """One k point in a rank-private mode file, read one signed mode at a time."""

    def __init__(self, stream, offset, shapes, signs, error):
        self.stream = stream
        self.offset = offset
        self.shapes = shapes
        self.ends = np.cumsum([int(np.prod(s)) for s in shapes])
        self.width = int(self.ends[-1])
        self.signs = signs
        self.error = error

    def adjoint(self, model, y):
        w, q = model.W, model.Q
        ww, qw = np.ix_(w, w), np.ix_(q, w)
        h = np.zeros_like(y)
        orbital = np.zeros((len(y), self.shapes[2][1]), dtype=complex)
        self.stream.seek(self.offset)
        for sign in self.signs:
            row = np.fromfile(self.stream, np.complex128, count=self.width)
            if row.size != self.width:
                raise OSError("incomplete SOC frequency mode")
            gw, gq, left, right = (a.reshape(s) for a, s in
                                   zip(np.split(row, self.ends[:-1]), self.shapes))
            gh = gw.conj().T
            h += sign * model.sandwich(gh[None], gq.conj()[None], y)[0]
            # G† Y (Σ U† G)† + G Y (G U Σ). Only its W columns enter.
            orbital[w] += sign * ((gh @ y[ww]) @ right.conj().T + (gw @ y[ww]) @ left)
            orbital[q] += sign * (gq.conj()[:, None] * (y[qw] @ right.conj().T)
                                  + gq[:, None] * (y[qw] @ left))
        return h, orbital


def _signed_basis(samples, weights, max_rank):
    """Use sampled columns to propose modes. Full-factor errors decide retention."""
    vectors, values, signs = [], [], []
    for sign in (1, -1):
        use = sign * weights > 0
        if not use.any():
            continue
        u, s, _ = svd(samples[use], full_matrices=False, check_finite=False)
        keep = min(max_rank, len(s))
        block = np.zeros((len(weights), keep))
        block[use] = u[:, :keep]
        vectors.append(block)
        values.extend(s[:keep])
        signs.extend([sign] * keep)
    order = np.argsort(values)[::-1][:max_rank]
    return np.column_stack(vectors)[:, order], np.asarray(signs)[order]


def _frequency_rows(model, k, u, weights, edges):
    """Generate weighted rows of the four factors in the model's usual chunks."""
    for start, stop, gw, gq, _ in model._chunks(k, need_static=False):
        sigma = model.sigma[start:stop]
        factors = (gw, gq, (gw @ u) @ sigma, (sigma @ u.conj().T) @ gw)
        nr = stop - start
        weight = np.sqrt(abs(weights[start:stop]))[:, None]
        rows = np.empty((nr, edges[-1]), dtype=complex)
        for f, a, b in zip(factors, edges[:-1], edges[1:]):
            np.multiply(f.reshape(nr, b-a), weight, out=rows[:, a:b])
        yield start, stop, rows


def _fit(model, k, tolerance, scratch_limit, stream):
    w, q = model.W, model.Q
    u = model.uc[k, w]
    shapes = ((len(w), len(w)), (len(q),), u.shape, u.T.shape)
    sizes = np.array([np.prod(s, dtype=int) for s in shapes])
    edges = np.r_[0, np.cumsum(sizes)]
    row_bytes = int(edges[-1]) * np.dtype(complex).itemsize
    weights = 2 * model.meta.temperature * model.meta.weights
    nf = len(weights)
    budget = min(model.cache_limit, 96 * 1024**2, scratch_limit)
    rank_limit = min(64, nf // 2, budget // row_bytes)
    if rank_limit < 1:
        return None

    # First pass: use a deterministic sample only to propose the real basis.
    # The full set of columns, including unsampled ones, is checked below.
    rng = np.random.default_rng(351)
    sample_count = min(256, (4 * 1024**2) // (4 * nf * 16))
    if sample_count < 1:
        return None
    columns = [a + rng.choice(size, min(size, sample_count), replace=False)
               for a, size in zip(edges, sizes)]
    sample = [np.empty((nf, len(cols)), complex) for cols in columns]
    for a, b, rows in _frequency_rows(model, k, u, weights, edges):
        for s, cols in zip(sample, columns):
            s[a:b] = rows[:, cols]
    normalized = [s.view(float) / norm for s in sample if (norm := np.linalg.norm(s)) > 0]
    if not normalized:
        return None
    basis, signs = _signed_basis(np.column_stack(normalized), weights, rank_limit)
    del sample, normalized
    rank_limit = len(signs)

    # Second pass: accumulate the coordinates of all columns in the basis.
    # Each column of coeff is one mode. BLAS updates it in place, keeping RAM
    # bounded independently of the number of Matsubara frequencies. Replaying
    # _chunks lets an uncached point trade two Green-function passes for the
    # full nf × width temporary file that would otherwise be required.
    coeff = np.zeros((int(edges[-1]), rank_limit), dtype=complex, order="F")
    for a, b, rows in _frequency_rows(model, k, u, weights, edges):
        coeff = zgemm(1.0, rows.T, basis[a:b], beta=1.0, c=coeff, overwrite_c=1)

    # Third pass: directly measure the residual orthogonal to the largest
    # basis. Subtracting almost equal squared norms here loses precision.
    # For a smaller rank, orthogonality adds the powers of the discarded modes.
    outside = np.zeros(4)
    for a, b, rows in _frequency_rows(model, k, u, weights, edges):
        residual = rows - basis[a:b] @ coeff.T
        for i, (l, r) in enumerate(zip(edges[:-1], edges[1:])):
            part = residual[:, l:r]
            outside[i] += (np.einsum("ij,ij->", part.real, part.real)
                           + np.einsum("ij,ij->", part.imag, part.imag))
    powers = np.zeros((4, rank_limit))
    for i, (l, r) in enumerate(zip(edges[:-1], edges[1:])):
        block = coeff[l:r]
        powers[i] = (np.einsum("ij,ij->j", block.real, block.real)
                     + np.einsum("ij,ij->j", block.imag, block.imag))
    errors = np.tile(outside[:, None], (1, rank_limit))
    errors[:, :-1] += np.cumsum(powers[:, :0:-1], axis=1)[:, ::-1]
    norms = outside + powers.sum(axis=1)
    relative = np.sqrt(np.divide(errors, norms[:, None], out=np.zeros_like(errors),
                                 where=norms[:, None] > 0)).max(axis=0)
    good = np.flatnonzero(relative <= tolerance)
    if not len(good):
        return None
    rank = int(good[0]) + 1

    # One anonymous file per MPI rank, regardless of the k mesh size. Writing
    # only retained rows leaves the uncompressed resolvent grid off disk.
    offset = stream.seek(0, 2)
    for mode in range(rank):
        if stream.write(memoryview(coeff[:, mode]).cast("B")) != row_bytes:
            raise OSError(errno.ENOSPC, "incomplete write of SOC frequency modes")
    return FrequencyModes(stream, offset, shapes, signs[:rank], float(relative[rank-1]))


def prepare_frequency_modes(model, tolerance, scratch_limit, resources):
    """Keep ordinary Matsubara evaluation at points that exceed a budget."""
    tables = {}
    if model.cache_limit <= 0 or scratch_limit <= 0:
        return tables
    space_errors = (errno.ENOSPC, errno.EDQUOT, errno.EFBIG)
    try:
        stream = resources.enter_context(TemporaryFile(buffering=0))
        for k in range(len(model.uc)):
            # Small windows can already have an exact frequency-summed kernel.
            if k in model._kernels:
                continue
            try:
                table = _fit(model, k, tolerance, scratch_limit, stream)
            except np.linalg.LinAlgError:
                # Failed local fit. Its ordinary electronic action is available.
                continue
            if table is not None:
                tables[k] = table
                scratch_limit -= len(table.signs) * table.width * 16
    except OSError as exc:
        if exc.errno not in space_errors:
            raise
        # Earlier records are complete. All subsequent k points use the exact sum.
    return tables
