"""WIEN lattice geometry shared by density, potential and band response.

WIEN's integer G labels do not always refer to a primitive reciprocal basis.
For F/B/C cells BR1 maps stored labels while BR2 spans the primitive reciprocal
lattice. Confusing the two gives incorrect momenta or factors of two/four in
charge and stress. Keep both maps explicit and take their conventions from
the same LATGEN implementation used by the existing input tools.

All matrices below act on COLUMN vectors. Methods accept batches of ROW
vectors for NumPy. Fractional sites/symmetries use the stored-label dual cell.
The integration volume is always the primitive-cell volume.
"""

from dataclasses import dataclass
from functools import lru_cache
from types import SimpleNamespace
import io
import numpy as np
from wlatgen import Latgen


@dataclass(frozen=True)
class NativeLattice:
    reciprocal: np.ndarray
    primitive_reciprocal: np.ndarray
    cell: np.ndarray
    primitive_cell: np.ndarray
    volume: float
    centering: np.ndarray

    @classmethod
    def from_structure(cls, structure):
        s = structure
        return _lattice(s.lattyp.strip(), *map(float, (
            s.a, s.b, s.c, s.alpha, s.beta, s.gamma)))

    def momenta(self, labels):
        return np.asarray(labels) @ self.reciprocal.T

    def cartesian_rotations(self, rotations):
        rotations = np.asarray(rotations)
        result = self.cell @ rotations @ np.linalg.inv(self.cell)
        if not np.allclose(result.swapaxes(-1, -2) @ result, np.eye(3),
                           atol=2e-7, rtol=0):
            raise ValueError("fractional symmetry is not Cartesian orthogonal")
        return result

    def equivalent_displacement(self, fractional, tolerance=2e-7):
        """Periodicity modulo primitive translations, including centering."""
        primitive = np.asarray(fractional) @ (
            np.linalg.solve(self.primitive_cell, self.cell)).T
        return bool(np.max(abs(primitive - np.rint(primitive))) < tolerance)

    def centering_factor(self, labels):
        """Project a conventional FFT box onto allowed reciprocal vectors.

        The full FFT box includes forbidden reflections in centred cells.
        A sphere mask formed from primitive sites alone must be multiplied by
        this factor; otherwise it places spurious holes in the conventional
        FFT cell. Allowed vectors retain the primitive-volume normalization.
        """
        labels = np.asarray(labels)
        phase = np.zeros(labels.shape[:-1], complex)
        for shift in self.centering:
            phase += np.exp(2j * np.pi * (labels @ shift))
        phase /= len(self.centering)
        phase[abs(phase) < 1e-12] = 0
        phase[abs(phase - 1) < 1e-12] = 1
        return phase

    def strained(self, strain, amplitude):
        """Independent geometry endpoint for derivative checks, no LATGEN rerun."""
        f = np.eye(3) + amplitude * np.asarray(strain)
        if f.shape != (3, 3) or not np.isfinite(f).all() or np.linalg.det(f) <= 0:
            raise ValueError("invalid deformation")
        return NativeLattice(np.linalg.solve(f.T, self.reciprocal),
                             np.linalg.solve(f.T, self.primitive_reciprocal),
                             f @ self.cell, f @ self.primitive_cell,
                             self.volume * np.linalg.det(f), self.centering)


@lru_cache(maxsize=64)
def _lattice(kind, a, b, c, alpha, beta, gamma):
    if kind not in {"P", "S", "F", "B", "H", "R", "CXY", "CXZ", "CYZ"}:
        raise ValueError("unsupported WIEN lattice type: " + kind)
    values = np.array([a, b, c, alpha, beta, gamma])
    if not np.isfinite(values).all() or min(a, b, c) <= 0:
        raise ValueError("invalid cell parameters")
    s = SimpleNamespace(lattyp=kind, a=a, b=b, c=c,
                        alpha=alpha, beta=beta, gamma=gamma)
    with np.errstate(invalid="raise", divide="raise"):
        lat = Latgen(s, "", io.StringIO())
    reciprocal = np.array(lat.br1, copy=True)
    primitive = np.array(lat.br2, copy=True)
    cell = 2 * np.pi * np.linalg.inv(reciprocal).T
    pcell = 2 * np.pi * np.linalg.inv(primitive).T
    volume = abs(float(np.linalg.det(pcell)))
    if not np.isclose(volume, lat.Vol, rtol=2e-12, atol=0):
        raise ValueError("LATGEN primitive volume disagrees with determinant")
    shifts = [[0., 0., 0.]]
    if kind == "F":
        shifts += [[0., .5, .5], [.5, 0., .5], [.5, .5, 0.]]
    elif kind == "B":
        shifts += [[.5, .5, .5]]
    elif kind.startswith("C"):
        shift = [0., 0., 0.]
        for axis in kind[1:]:
            shift["XYZ".index(axis)] = .5
        shifts += [shift]
    # R's BR1 already refers to the primitive rhombohedral coordinate system.
    centering = np.array(shifts)
    ratio = abs(np.linalg.det(cell)) / volume
    if not np.isclose(ratio, len(shifts), rtol=2e-12, atol=0):
        raise ValueError("stored-label cell and centering multiplicity disagree")
    for x in [reciprocal, primitive, cell, pcell, centering]:
        x.setflags(write=False)
    return NativeLattice(reciprocal, primitive, cell, pcell, volume, centering)


def native_site_frames(structure):
    """DMFT2's ROTDEF gauge, including its identity-first symmetry swap.

    Several rotations can map a representative site onto the same atom.
    Their charge densities are equivalent, but their fixed local projector
    matrices are not interchangeable. DMFT2 swaps the identity operation to
    slot one BEFORE ROTDEF. Reproduce that ordering, then use the inverse
    Cartesian rotation for momentum and the unreduced real-space position
    for the Bloch phase (including any centering translation).
    """
    s=structure; lattice=NativeLattice.from_structure(s)
    order=list(range(len(s.timat)))
    identities=[i for i,r in enumerate(s.timat) if np.array_equal(r,np.eye(3,dtype=int))]
    if not identities: raise ValueError("missing native identity symmetry")
    i=identities[0];order[0],order[i]=order[i],order[0]
    first=np.r_[0,np.cumsum(s.mult)[:-1]]
    axes=[];locations=[]
    for site,a in enumerate(np.repeat(np.arange(s.nat),s.mult)):
        for i in order:
            rot=s.timat[i];tau=s.tau[i]
            location=rot@s.pos[first[a]]+tau
            if lattice.equivalent_displacement(location-s.pos[site]):
                axes.append(lattice.cartesian_rotations(rot).T)
                locations.append(location)
                break
        else: raise ValueError("missing native equivalent-site rotation")
    return np.asarray(axes),np.asarray(locations)


def native_structure(structure):
    """Return a copy with ROTDEF's finite-precision site repair applied.

    LAPW0/LAPW1/DMFT2 average symmetry images lying within 0.01 Bohr of each
    input site, twice. Eight-decimal STRUCT fractions such as 1/3 otherwise
    disagree with their native phases enough to fail H/O reconstruction.
    This reproduces the native input map; it never edits case.struct.
    """
    import copy
    s=copy.deepcopy(structure)
    s.pos=np.array(s.pos,dtype=float,copy=True)
    lengths=np.array([s.a,s.b,s.c])
    for _ in range(2):
        for a,position in enumerate(s.pos):
            images=(np.einsum('sij,j->si',s.timat,position)+s.tau+1.)%1.
            images=np.where(images>.999999,images-1.,images)
            close=np.all(abs(images-position)*lengths<1e-2,axis=1)
            if np.count_nonzero(close)>1:s.pos[a]=images[close].mean(axis=0)
    return s
