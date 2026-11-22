"""Read the native spinor export without keeping LAPW vectors for the entire BZ."""

from dataclasses import dataclass
from pathlib import Path
import numpy as np
from scipy.special import expit
from frozen_sigma_response import Metadata
from frozen_sigma_tangent import CachedFrozenSigmaResponse
from response_parallel import SERIAL


def saved_normalization_error(path, meta):
    """Largest change from the parent SOlapm.dat in its native impurity frame."""
    largest = 0.0
    with Path(path).open() as stream:
        for ic, dim in enumerate(meta.dimensions):
            if not stream.readline():
                raise ValueError(f"incomplete orbital normalization in {path}")
            for i in range(dim):
                line = stream.readline().rstrip("\r\n")
                fields = [line[j : j + 20] for j in range(0, 40 * dim, 20)]
                row = np.asarray(
                    [
                        complex(float(re), float(im))
                        for re, im in zip(fields[::2], fields[1::2])
                    ]
                )
                if not np.all(np.isfinite(row)):
                    raise ValueError(f"nonfinite orbital normalization in {path}")
                largest = max(
                    largest, float(np.max(abs(row - meta.normalization[i, :dim, ic])))
                )
    return largest


@dataclass
class SpinRecord:
    path: Path
    index: int
    basis: int
    n: int
    first: int
    last: int
    wfirst: int
    wlast: int
    orbit_size: int
    norb: int
    nlo: int
    weight: float

    @classmethod
    def read(cls, path):
        with path.open("rb") as stream:
            if stream.read(8) != b"EDMRPK01":
                raise ValueError(f"invalid DMFT2 k record {path}")
            v, k, b, iso, n, first, last, wl, wh, m, o, nlo = np.fromfile(
                stream, "<i4", 12
            )
            w = np.fromfile(stream, "<f8", 1)[0] / 2
            if (
                (v, iso) != (1, 2)
                or not first <= wl <= wh <= last
                or last - first + 1 != n
                or not 0 < w <= 1
            ):
                raise ValueError(f"invalid spinor band range or weight at {path}")
            expected = 64 + 8 * n + 16 * n * (2 * b + m * o) + 12 * b
            if path.stat().st_size != expected:
                raise ValueError(f"incomplete spinor record at {path}")
            return cls(path, k, b, n, first, last, wl, wh, m, o, nlo, w)

    def arrays(self):
        # Memory maps release the k vectors at the end of each local k action.
        off = 64
        e = np.memmap(self.path, dtype="<f8", mode="r", offset=off, shape=(self.n,))
        off += e.nbytes
        c = np.memmap(
            self.path,
            dtype="<c16",
            mode="r",
            offset=off,
            shape=(self.basis, self.n, 2),
            order="F",
        )
        off += c.nbytes
        u = np.memmap(
            self.path,
            dtype="<c16",
            mode="r",
            offset=off,
            shape=(self.n, self.orbit_size, self.norb),
            order="F",
        )
        off += u.nbytes
        labels = np.memmap(
            self.path,
            dtype="<i4",
            mode="r",
            offset=off,
            shape=(3, self.basis),
            order="F",
        )
        return e, c, u, labels


def prepare_electronic(
    prefix, *, parallel=SERIAL, cache_mb=128.0, retained=None, require_fresh=True
):
    prefix = Path(prefix)
    meta = parallel.root_call(lambda: Metadata.read(str(prefix) + ".meta"))

    def inventory():
        files = sorted(prefix.parent.glob(prefix.name + ".k????????.bin"))
        records = [SpinRecord.read(f) for f in files]
        if not records or [r.index for r in records] != list(
            range(1, len(records) + 1)
        ):
            raise ValueError("incomplete spinor k-point sequence")
        if (
            len({(r.first, r.wfirst, r.wlast, r.orbit_size, r.norb) for r in records})
            != 1
        ):
            raise ValueError("a common fixed-index DMFT window is required")
        if abs(sum(r.weight for r in records) - 1) > 1e-10:
            raise ValueError("spinor weights do not cover the Brillouin zone")
        return records

    all_r = parallel.root_call(inventory)
    local = [all_r[i] for i in parallel.partition([float(r.basis) ** 3 for r in all_r])]
    nr = (all_r[0].wlast - all_r[0].first + 1) if retained is None else retained
    if not isinstance(nr, (int, np.integer)) or not all(
        r.wlast - r.first < nr <= r.n for r in all_r
    ):
        raise ValueError("retained spinor range does not cover DMFT window")
    energies, raws, norms = [], [], []
    for r in local:
        e, _, u, _ = r.arrays()
        if np.any(np.diff(e) < -1e-9):
            raise ValueError(f"unordered spinor energy spectrum at k={r.index}")
        if np.sum(expit((meta.mu - e[nr:]) / meta.temperature)) > 1e-12 * min(
            1.0, meta.temperature
        ):
            raise ValueError(
                f"more thermally occupied spinors must be retained at k={r.index}"
            )
        norm = np.array(u[:nr])
        # Existing eDMFT normalization acts per impurity, after k-summation.
        # Undo it only in the exported orbital frame; spin rotations are intact.
        loc = meta.local_columns(norm[None])[0]
        raw = np.zeros_like(norm)
        for ic, n in enumerate(meta.dimensions):
            block = loc[:, :n, ic] @ np.linalg.inv(meta.normalization[:n, :n, ic])
            for orb in np.flatnonzero(meta.cix == ic):
                ids = meta.isx[: meta.nindo[orb], orb]
                raw[:, : len(ids), orb] = block[:, ids]
        raws.append(raw)
        norms.append(norm)
        energies.append(e[:nr].copy())
    r = local[0]
    window = np.arange(r.wfirst - r.first, r.wlast - r.first + 1)
    model = CachedFrozenSigmaResponse(
        meta,
        np.asarray(energies),
        np.asarray(raws),
        np.asarray(norms),
        np.array([r.weight for r in local]),
        window,
        parallel=parallel,
        cache_mb=cache_mb,
        frequency_chunk=4,
        precontract=False,
    )
    # A fixed SOlapm in a stale frame would differentiate a different functional.
    for ic, n in enumerate(meta.dimensions):
        idx = np.flatnonzero(meta.active[:n, ic])
        proj = np.ix_(idx, idx)
        x = model.X[:, :, ic][proj]
        if require_fresh and not np.allclose(
            x @ model.s0[:, :, ic][proj] @ x, np.eye(len(idx)), rtol=0, atol=2e-8
        ):
            raise ValueError(
                "exported orbital normalization is stale for the spinor k sum"
            )
    return model, local


def operator_readout(prefix, model, records):
    """Use native rho and C=Sigma U^H G for a direct force, before band rotation.

    The in-memory response already knows the Frechet derivative of the global
    normalization and the window denominators. The compact native moments
    eliminate a repetition of the Matsubara inversion at every k point.
    Historical snapshots without any such records can be checked with the
    Matsubara implementation. A partial export is an error.
    """
    files = [Path(str(prefix) + f".operator.k{r.index:08d}.bin") for r in records]
    count = model.parallel.sum(sum(f.is_file() for f in files))
    if count == 0:
        return None
    if count != model.parallel.sum(len(records)):
        raise ValueError("incomplete native SOC force readout")
    nk, nr = model.eps.shape
    density = np.zeros((nk, nr, nr), dtype=complex)
    compact = np.zeros_like(model.uc)
    for j, (r, path) in enumerate(zip(records, files)):
        n = r.wlast - r.wfirst + 1
        if path.stat().st_size != 44 + 16 * n * (n + r.orbit_size * r.norb):
            raise ValueError(f"invalid SOC force readout size at k={r.index}")
        with path.open("rb") as src:
            if src.read(8) != b"SOCFOR01":
                raise ValueError(f"invalid SOC force readout signature in {path}")
            ids = np.fromfile(src, "<i4", 7)
            weight = np.fromfile(src, "<f8", 1)[0]
            if (
                tuple(ids)
                != (r.index, r.first, r.wfirst, r.wlast, n, r.orbit_size, r.norb)
                or weight != 2 * r.weight
            ):
                raise ValueError(f"incompatible central SOC moment at k={r.index}")
            rho = np.fromfile(src, "<c16", n * n).reshape((n, n), order="F")
            c = np.fromfile(src, "<c16").reshape((n, r.orbit_size, r.norb), order="F")
        if not np.isfinite(c).all() or not np.allclose(
            rho, rho.conj().T, atol=1e-12, rtol=1e-10
        ):
            raise ValueError(
                f"invalid native SOC density or projector gradient at k={r.index}"
            )
        occupied = np.arange(r.wfirst - r.first)
        density[j, occupied, occupied] = weight
        density[j][np.ix_(model.W, model.W)] = rho
        compact[j, model.W] = model.meta.compact(
            model.meta.local_columns(2 * c.conj()[None])
        )[0]
    hbar, raw_bar = model.embedding_pullback(compact)
    return density, 0.5 * (density + hbar), 0.5 * raw_bar
