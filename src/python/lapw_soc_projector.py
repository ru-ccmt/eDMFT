"""Spinor projector 5 in the original LAPW rows (single-site impurity blocks).

Each radial projector is carried by its atom, with its local axes held fixed.
The site phase therefore uses k+G even for rows representing local orbitals.
The native rotation is applied before the global orbital normalization X.
"""

from pathlib import Path
import numpy as np

from indmffile import Indmfl
from lapw_lattice import NativeLattice, native_site_frames
from lapw_augmentation_response import PreparedAugmentationChannel, hydrostatic_matching_tangent
from lapw_radial_native import NativeRadialContext
from lapw_projector_response import PreparedProjector5


def spin_frame(path):
    with Path(path).open("rb") as f:
        if f.read(8) != b"RSPIN001":
            raise ValueError("unrecognized DMFT2 spin frame")
        m, n = np.fromfile(f, "<i4", 2)
        r = np.fromfile(f, "<c16", 4 * n).reshape((2, 2, n), order="F")
        cf = np.fromfile(f, "<c16", m * m * n * n).reshape((m, m, n, n), order="F")
        if f.read(1):
            raise ValueError("unexpected bytes after DMFT2 spin frame")
    return r, cf


class SpinorProjectors:
    def __init__(self, mat, directory, case, prefix, metadata, stack, library=None):
        self.mat = mat
        self.meta = metadata
        self.rspin, self.cf = spin_frame(str(prefix) + ".spin")
        inl = Indmfl(str(Path(directory) / case))
        inl.read()
        self.axes, self.locations = native_site_frames(mat.s)
        self.lattice = NativeLattice.from_structure(mat.s)
        atoms = np.repeat(np.arange(mat.s.nat), mat.s.mult)
        with (Path(directory) / "projectorw.dat").open() as f:
            fixed = {}
            for _ in range(int(f.readline().split()[1])):
                n, nr, atom, l = map(int, f.readline().split()[1:])
                if n != nr:
                    raise NotImplementedError(
                        "projector extends outside its muffin tin"
                    )
                fixed[atom - 1, l] = np.loadtxt([next(f) for _ in range(n)])[:, 1:3]
        self.entries = []
        self.contexts = {}
        for ic, orbs in inl.cix.items():
            native = np.flatnonzero(metadata.cix == ic - 1)
            if len(orbs) != 1 or len(native) != 1:
                raise NotImplementedError(
                    "spin response needs a single correlated shell per block"
                )
            site, l, _ = orbs[0]
            site -= 1
            orb = native[0]
            a = int(atoms[site])
            _, xyz, shift, _ = inl.atoms[site + 1]
            if shift is not None:
                raise NotImplementedError("shifted correlated shell")
            xyz = np.asarray(xyz) if len(xyz) else mat.s.rotloc[a]
            if xyz.shape != (3, 3):
                raise NotImplementedError("expected three fixed local axes")
            key = a, l
            if key not in self.contexts:
                e = mat.header[a][l]
                vr = mat.v[a][0, 0]
                s = mat.s
                self.contexts[key] = stack.enter_context(
                    NativeRadialContext(
                        vr,
                        r0=s.r0[a],
                        dx=np.log(s.rmt[a] / s.r0[a]) / (len(vr) - 1),
                        rmt=s.rmt[a],
                        energy=e - 200 if e > 150 else e,
                        lo_energies=mat.lo[a][l, : mat.ilo[a][l]] if l < 4 else [],
                        ell=l,
                        charge=s.Znuc[a],
                        relativistic="RELA" in s.mode,
                        lapw=e < 150,
                        library=library,
                    )
                )
            size = 2 * (2 * l + 1)
            if metadata.nindo[orb] != size or metadata.dimensions[ic - 1] != size:
                raise ValueError("SOC shell and native orbital dimensions disagree")
            # URx in DMFT2 is rotated as a ket and then conjugated.  Native
            # DMFT2 uses (up,down), whereas LAPWSO labels its blocks (down,up).
            t = (
                np.kron(self.rspin[:, :, orb].conj(), np.eye(2 * l + 1))
                @ self.cf[:size, :size, orb, orb].T
            )
            self.entries.append(
                (site, key, ic - 1, orb, xyz @ self.axes[site], fixed[key], t)
            )

    def prepare(self, k, labels):
        frac = labels.T + k
        q = self.lattice.momenta(frac)
        n = len(q)
        unit = np.eye(n)
        npw = n - self.mat.nlo_total
        ops = []
        for site, key, ic, orb, frame, fixed, t in self.entries:
            a, l = key
            ctx = self.contexts[key]
            ch = PreparedAugmentationChannel(
                l,
                ctx.rmt,
                self.lattice.volume,
                q @ frame.T,
                np.linalg.norm(q, axis=1),
                np.exp(2j * np.pi * (frac @ self.locations[site])),
                unit,
                npw,
                [npw + self.mat.lo_indices[a, l, j] for j in range(ctx.nlo)],
            )
            op = PreparedProjector5(
                ch,
                ctx.primitive,
                ctx.lo_coefficients,
                ctx.projector_metric,
                fixed,
                np.eye(2 * l + 1),
                lapw=ctx.lapw,
            )
            from lapw_geometry_response import SpinorMatchingGeometry
            op.strain_geometry = SpinorMatchingGeometry(q, q @ frame.T, frame.T,
                np.exp(2j * np.pi * (frac @ self.locations[site])), self.lattice.volume, npw,
                [npw + self.mat.lo_indices[a, l, j] for j in range(ctx.nlo)])
            ops.append((site, key, ic, op, t))
        return SpinorProjectorAtK(ops, q, self.meta, len(self.locations))


class SpinorProjectorAtK:
    def __init__(self, ops, q, metadata, nsites):
        self.ops, self.q, self.meta = ops, q, metadata
        self.nsites = nsites
        self.phi = {(site, ic): op.value() for site, key, ic, op, t in ops}

    def evaluate(self, psi, *, dpsi=None, radial=None, velocity=None, hydrostatic=False, strain=None):
        """Return raw U or its full tangent, in [band, index, impurity] order."""
        tangent = strain is not None or hydrostatic or any(v is not None for v in (dpsi, radial, velocity))
        m, nc = self.meta.normalization.shape[0], len(self.meta.dimensions)
        out = np.zeros((psi.shape[-1], m, nc), complex)
        for site, key, ic, op, t in self.ops:
            phi = self.phi[site, ic]

            def pair(waves, p):
                return np.column_stack([waves[j].conj().T @ p for j in (1, 0)]) @ t

            if not tangent:
                value = pair(psi, phi)
            else:
                dp = np.zeros_like(phi)
                if radial is not None:
                    dp += op.action(*radial[key])
                if hydrostatic:
                    # Channels here act on the unit matrix of LAPW rows.
                    # Plane-wave columns can be weighted directly; LO columns
                    # of fj vanish, and their matching normalization is separate.
                    q2 = np.einsum("ij,ij->i", self.q, self.q)
                    dot = hydrostatic_matching_tangent(op.channel, op.channel.fj * q2)
                    dz = dot.value(op.phi[-1], op.lo, lapw=op.lapw)
                    dp += np.einsum("min,i,ma->na", dz.conj(), op.q, op.transform)
                if strain is not None:
                    dot = op.strain_geometry.tangent(op.channel, strain)
                    dz = dot.value(op.phi[-1], op.lo, lapw=op.lapw)
                    dp += np.einsum("min,i,ma->na", dz.conj(), op.q, op.transform)
                if velocity is not None:
                    dp -= (1j * (self.q @ velocity[site]))[:, None] * phi
                value = pair(psi, dp)
                if dpsi is not None:
                    value += pair(dpsi, phi)
            out[:, : value.shape[1], ic] = value
        return out

    def pullback(self, psi, y):
        """Adjoint to both spinor legs, radial functions and atomic positions."""
        gw = np.zeros_like(psi)
        gu = np.zeros((self.nsites, 3))
        radial = {}
        for site, key, ic, op, t in self.ops:
            phi = self.phi[site, ic]
            bar = y[:, : t.shape[1], ic] @ t.conj().T
            gphi = np.zeros_like(phi)
            for spin, block in zip((1, 0), np.hsplit(bar, 2)):
                gw[spin] += phi @ block.conj().T
                gphi += psi[spin] @ block
            gu[site] += self.q.T @ (gphi.conj() * (-1j * phi)).real.sum(axis=1)
            pair = op.pullback(gphi)
            if key in radial:
                for target, add in zip(radial[key], pair):
                    target += add
            else:
                radial[key] = pair
        return gw, radial, gu
