"""Ragged-k LAPW potential and spatial response with only RR/QR blocks.

R is the common retained KS space; Q completes the finite LAPW space and may
have a different size at every k point. Electronic response, k weights, global
projector normalization and fixed-N chemical potential remain in ONE shared
EmptyComplementResponse. These adapters neither split that problem into
independent k points nor add weights to its density output.

The caller supplies central native contexts and coefficients in one gauge.
No strain source, screened solve, native-file loader or production switch is
implied by constructing these operators. Only fixed geometry/radial energies
are supported here. QQ density vanishes under the electronic empty-Q checks.
"""

import numpy as np
from lapw_density_response import checked, herm
from lapw_spherical_potential_response import SphericalPotentialBands
from lapw_complement_response import (
    ComplementMatrixCoordinates,
    ComplementPerturbationCoordinates,
)
from stress_response_adjoint import PairedOperator


class CompactSphericalPotentialBands(SphericalPotentialBands):
    """Shared radial work, variable coefficient sizes, no full band H/O arrays.

    Spherical entries use PreparedSphericalBand with a rectangular Wronskian
    kernel (Nr+Nq[k], Nr). Nonspherical matching and raw-projector operators
    may contain R+Q rows, but only the retained raw derivative is needed: the
    central Q overlap supplies the moving-window leg in EmptyComplementResponse.
    No derivative of the central Q overlap is needed at first order.
    """

    def __init__(
        self,
        contexts,
        spherical,
        nonspherical,
        projectors,
        complement_sizes,
        raw_shape,
        *,
        energy_policy,
        projector_contexts=None,
    ):
        if energy_policy != "fixed":
            raise NotImplementedError("only fixed radial energies are supported")
        if len(raw_shape) != 4:
            raise ValueError("invalid raw-projector shape")
        nk, nr = raw_shape[:2]
        matrix = ComplementMatrixCoordinates((nk, nr, nr), complement_sizes)
        self.coordinates = ComplementPerturbationCoordinates(matrix, raw_shape)
        self.contexts = contexts
        # H/O and raw projectors may originate from different radial
        # quadratures. Native WIEN24 LAPW1 uses CHARGE, whereas dmft2's
        # ATPARN uses its historical RINT13 rule. Differentiate each actual
        # map, then sum their covectors on the same physical potential input.
        self.projector_contexts = (
            contexts if projector_contexts is None else projector_contexts
        )
        if set(self.projector_contexts) != set(contexts):
            raise ValueError("projector and matrix atoms differ")
        for a in contexts:
            for l, c in self.projector_contexts[a].items():
                if l not in contexts[a] or not np.array_equal(
                    c.radius, contexts[a][l].radius
                ):
                    raise ValueError("projector and matrix radial grids differ")
        self.spherical, self.nonspherical, self.projectors = map(
            list, (spherical, nonspherical, projectors)
        )
        self.slices = {}
        offset = 0
        for a in sorted(contexts):
            if not contexts[a]:
                raise ValueError("empty radial contexts")
            n = next(iter(contexts[a].values())).n
            if any(c.n != n for c in contexts[a].values()):
                raise ValueError("inconsistent radial grids")
            self.slices[a] = slice(offset, offset + n)
            offset += n
        if not self.slices:
            raise ValueError("no radial contexts")
        self.input_size = offset
        covered = set()
        for a, l, k, op in self.spherical:
            if not 0 <= k < nk or op.triangle.shape != (nr + matrix.shapes[k][0], nr):
                raise ValueError("inconsistent compact spherical layout")
            if not np.array_equal(op.phi, contexts[a][l].primitive):
                raise ValueError("inconsistent spherical radial context")
            covered.add((a, l, k))
        if covered != {
            (a, l, k) for a in contexts for l in contexts[a] for k in range(nk)
        }:
            raise ValueError("missing spherical atom/l/k sectors")
        if {a for a, op, ks in self.nonspherical} != set(contexts):
            raise ValueError("missing nonspherical atom")
        for a, op, ks in self.nonspherical:
            if len(ks) != len(op.entries):
                raise ValueError("missing nonspherical k index")
            for k, entry in zip(ks, op.entries):
                if not 0 <= k < nk or entry[0].shape[1] != nr + matrix.shapes[k][0]:
                    raise ValueError("inconsistent nonspherical band layout")
            for l, p in op.phi.items():
                if not np.array_equal(p, contexts[a][l].primitive):
                    raise ValueError("inconsistent nonspherical radial context")
        for a, l, k, ic, cols, op in self.projectors:
            if (
                not 0 <= k < nk
                or not 0 <= ic < raw_shape[3]
                or op.shape[0] not in (nr, nr + matrix.shapes[k][0])
                or op.shape[1] != len(range(raw_shape[2])[cols])
            ):
                raise ValueError("inconsistent raw-projector layout")
            if not np.array_equal(op.phi, self.projector_contexts[a][l].primitive):
                raise ValueError("inconsistent raw-projector context")
        self.operator = PairedOperator(
            offset,
            self.coordinates.size,
            self.action,
            self.pullback,
            "spherical potential to ragged RR/QR H/O/projector",
        )

    def action(self, vector):
        t = self.tangents(vector)
        pt = t
        if self.projector_contexts is not self.contexts:
            pt = {a: {} for a in self.contexts}
            for a, l in {(a, l) for a, l, *_ in self.projectors}:
                pt[a][l] = self.projector_contexts[a][l].action(vector[self.slices[a]])
        c = self.coordinates
        nk, nr = c.raw_shape[:2]
        h = np.zeros((nk, nr, nr), complex)
        o = np.zeros_like(h)
        hq = [np.zeros(shape, complex) for shape in c.matrix.shapes]
        oq = [np.zeros_like(q) for q in hq]
        raw = np.zeros(c.raw_shape, complex)
        for a, l, k, op in self.spherical:
            dh, do = op.action(*t[a][l])
            h[k] += dh[:nr]
            hq[k] += dh[nr:]
            o[k] += do[:nr]
            oq[k] += do[nr:]
        for a, op, ks in self.nonspherical:
            dp = {l: t[a][l][0] for l in op.phi}
            dm = {l: t[a][l][1] for l in op.phi}
            for k, dh in zip(ks, op.action(dp, dm, retained=nr)):
                h[k] += dh[:nr]
                hq[k] += dh[nr:]
        for a, l, k, ic, cols, op in self.projectors:
            raw[k, :, cols, ic] += op.action(*pt[a][l])[:nr]
        return c.pack(h, o, raw, hq, oq)

    def pullback(self, vector):
        hb, ob, rb, hq, oq = self.coordinates.unpack(vector)
        nr = rb.shape[1]
        # Cache only the needed columns per k. No QQ reconstruction occurs.
        hc = [np.concatenate((h, q)) for h, q in zip(hb, hq)]
        oc = [np.concatenate((o, q)) for o, q in zip(ob, oq)]
        gp = {
            a: {l: np.zeros_like(c.primitive) for l, c in cs.items()}
            for a, cs in self.contexts.items()
        }
        gm = {
            a: {l: np.zeros_like(c.lo_coefficients) for l, c in cs.items()}
            for a, cs in self.contexts.items()
        }
        for a, l, k, op in self.spherical:
            p, m, *_ = op.pullback(hc[k], oc[k])
            gp[a][l] += p
            gm[a][l] += m
        for a, op, ks in self.nonspherical:
            p, m, _ = op.pullback([hc[k] for k in ks])
            for l in op.phi:
                gp[a][l] += p[l]
                gm[a][l] += m[l]
        answer = np.zeros(self.input_size)
        if self.projector_contexts is self.contexts:
            pg, mg = gp, gm
        else:
            pg = {a: {} for a in self.contexts}
            mg = {a: {} for a in self.contexts}
            for a, l in {(a, l) for a, l, *_ in self.projectors}:
                ctx = self.projector_contexts[a][l]
                pg[a][l] = np.zeros_like(ctx.primitive)
                mg[a][l] = np.zeros_like(ctx.lo_coefficients)
        for a, l, k, ic, cols, op in self.projectors:
            y = np.zeros(op.shape, complex)
            y[:nr] = rb[k, :, cols, ic]
            p, m = op.pullback(y)
            pg[a][l] += p
            mg[a][l] += m
        if self.projector_contexts is not self.contexts:
            for a in pg:
                for l in pg[a]:
                    answer[self.slices[a]] += self.projector_contexts[a][l].pullback(
                        pg[a][l], mg[a][l]
                    )[0]
        for a, cs in self.contexts.items():
            for l, c in cs.items():
                answer[self.slices[a]] += c.pullback(gp[a][l], gm[a][l])[0]
        return answer


class CompactSpatialValenceDensity:
    """Add ragged QR/RQ density to an ordinary RETAINED-space spatial map.

    retained_spatial owns the central density and all direct moving-basis
    contributions. mt_complement[atom] supplies one Z_Q per entry, in exactly
    the same k/site order and representative harmonic frame. ir_complement[k]
    supplies PW coefficients C_Q with the same G rows as the retained map.

    Only retained wavefunctions are Fourier transformed. The cross-density
    uses FFT(C_Q d_rho_QR), with Nr columns, not FFT(C_Q) with Nq columns.
    This controls both memory and work when Nq is much larger than Nr.
    """

    def __init__(
        self, retained_spatial, complement_sizes, mt_complement, ir_complement
    ):
        self.retained = retained_spatial
        self.template = retained_spatial.template
        self.coordinates = ComplementMatrixCoordinates(
            retained_spatial.rho.shape, complement_sizes
        )
        nk, nr, _ = retained_spatial.rho.shape
        if (
            set(mt_complement) != set(retained_spatial.atoms)
            or len(ir_complement) != nk
        ):
            raise ValueError("missing complement density sectors")
        self.mt = {}
        for a, atom in retained_spatial.atoms.items():
            qs = mt_complement[a]
            if len(qs) != len(atom.entries):
                raise ValueError("missing complement site entries")
            self.mt[a] = []
            for q, (k, z, w) in zip(qs, atom.entries):
                v = np.array(
                    checked(
                        q,
                        (len(z), self.coordinates.shapes[k][0]),
                        "MT complement matching",
                    ),
                    complex,
                )
                v.setflags(write=False)
                self.mt[a].append(v)
        self.ir = []
        for k, (q, op) in enumerate(zip(ir_complement, retained_spatial.interstitial)):
            v = np.array(
                checked(
                    q,
                    (len(op.c), self.coordinates.shapes[k][0]),
                    "IR complement coefficients",
                ),
                complex,
            )
            v.setflags(write=False)
            self.ir.append(v)
        self.operator = PairedOperator(
            self.coordinates.size,
            len(self.template.vector()),
            self.action,
            self.pullback,
            "ragged RR/QR spatial density",
        )

    def value(self):
        return self.retained.value()

    def action(self, vector):
        rr, qr = self.coordinates.unpack(vector)
        radial = {}
        for a, atom in self.retained.atoms.items():
            dd = np.zeros((atom.kernel.layout.size,) * 2, complex)
            for (k, z, w), q in zip(atom.entries, self.mt[a]):
                cross = (q @ qr[k]) @ z.conj().T
                dd += w * ((z @ rr[k]) @ z.conj().T + cross + cross.conj().T)
            radial[a] = atom.kernel.component_action(herm(dd))
        fourier = np.zeros(len(self.template.pw), complex)
        for k, (op, cq) in enumerate(zip(self.retained.interstitial, self.ir)):
            psi = op.psi if op.psi is not None else op._waves(op.c)
            # The factor two is physical QR+RQ, not an extra spin weight.
            waves = op._waves(op.c @ rr[k] + 2 * cq @ qr[k])
            density = np.einsum("...n,...n->...", waves, psi.conj()).real
            fourier += op.stars.action(np.fft.fftn(density) / (op.ng * op.volume))
        return self.retained._assemble(radial, fourier)

    def pullback(self, vector):
        v = checked(vector, self.template.vector().shape, "density covector", real=True)
        y = self.template.from_vector(v)
        rr = np.zeros_like(self.retained.rho)
        qr = [np.zeros(shape, complex) for shape in self.coordinates.shapes]
        for a, atom in self.retained.atoms.items():
            ym = np.array([y.rad[(a, *h)] for h in atom.kernel.layout.harmonics])
            bar, _ = atom.kernel.component_pullback(ym, include_radial=False)
            for (k, z, w), q in zip(atom.entries, self.mt[a]):
                bz = bar @ z
                rr[k] += w * herm(z.conj().T @ bz)
                # coordinates.pack supplies the QR+RQ pairing, so no factor 2 here.
                qr[k] += w * (q.conj().T @ bz)
        yi = np.array(list(y.pw.values()))
        for k, (op, cq) in enumerate(zip(self.retained.interstitial, self.ir)):
            potential = np.fft.ifftn(op.stars.pullback(yi)).real / op.volume
            psi = op.psi if op.psi is not None else op._waves(op.c)
            vc = op._wave_pullback(potential[..., None] * psi)
            rr[k] += herm(op.c.conj().T @ vc)
            qr[k] += cq.conj().T @ vc
        return self.coordinates.pack(rr, qr)
