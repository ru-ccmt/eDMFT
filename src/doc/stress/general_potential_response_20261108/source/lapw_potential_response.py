"""Native spherical, nonspherical and interstitial potential response inputs.

Spherical inputs are r*delta V_spherical in Ry. Nonspherical inputs are native
real-harmonic delta V_lm in Ry. Interstitial inputs are the ALREADY REANALYZED
case.vns Fourier coefficients of theta_interstitial*delta V, in Ry. Never add
another MT mask or divide by volume/star multiplicity here. Converting an
unmasked physical potential into these native coordinates belongs upstream.

These maps keep geometry, RMT, local axes and radial energies fixed. Only the
spherical sector changes the radial/matching basis; every sector changes the
KS eigenvectors and normalized window through the electronic response.
"""

import numpy as np
from lapw_density_response import checked, herm
from lapw_complement_response import (
    ComplementMatrixCoordinates,
    ComplementPerturbationCoordinates,
)
from stress_response_adjoint import PairedOperator


class NativeInterstitialPotential:
    """Native STERN/SETWAR/upper-triangle PW potential and its transpose.

    coefficients contains ONLY PW rows in native order, including retained and
    complement columns. labels is integer G, not k+G (k cancels in differences).
    Each representative first expands with the conjugated mean translation
    phase, as in SETWAR. TAPEWF writes CONJG(Z) in the complex executable:
    the assembled upper matrix must therefore be CONJ(WARP(G_row-G_col)).
    Conjugating the value, rather than merely reversing G, also reproduces
    the native extension for redundant input coordinates. No star-size factor
    is present. LO rows are absent because LOs have no interstitial leg.

    The linear map reproduces the native Hermitian upper-triangle copy also
    for redundant input coordinates. The physical caller supplies symmetry-
    compatible star coefficients. This map does not invent extra physical
    degrees of freedom or project arbitrary coordinates onto allowed fields.
    """

    def __init__(
        self,
        coefficients,
        labels,
        representatives,
        rotations,
        translations,
        *,
        complex_values
    ):
        c, g, rep, rot, tau = map(
            np.asarray, (coefficients, labels, representatives, rotations, translations)
        )
        if (
            c.ndim != 2
            or min(c.shape) < 1
            or g.shape != (len(c), 3)
            or not np.issubdtype(g.dtype, np.integer)
        ):
            raise ValueError("invalid interstitial PW coefficients/labels")
        checked(c, c.shape, "interstitial coefficients")
        if len(set(map(tuple, g))) != len(g):
            raise ValueError("duplicate PW labels")
        if (
            rep.ndim != 2
            or rep.shape[1] != 3
            or not len(rep)
            or not np.issubdtype(rep.dtype, np.integer)
        ):
            raise ValueError("invalid potential star representatives")
        if (
            rot.ndim != 3
            or rot.shape[1:] != (3, 3)
            or not len(rot)
            or not np.issubdtype(rot.dtype, np.integer)
        ):
            raise ValueError("invalid native reciprocal rotations")
        checked(tau, (len(rot), 3), "native translations", real=True)
        if not isinstance(complex_values, (bool, np.bool_)):
            raise ValueError("complex_values must be explicit")
        self.c = np.array(c, complex, copy=True)
        self.nband = c.shape[1]
        self.nstar = len(rep)
        self.complex_values = bool(complex_values)
        self.representatives = rep.copy()
        members = {}
        for j, k in enumerate(rep):
            star = {}
            for r, t in zip(rot, tau):
                phase = (
                    np.exp(2j * np.pi * (k @ t))
                    if complex_values
                    else complex(np.cos(2 * np.pi * (k @ t)))
                )
                star.setdefault(tuple(r @ k), []).append(phase)
            for key, phases in star.items():
                if key in members:
                    raise ValueError("overlapping potential stars")
                members[key] = (j, np.mean(phases).conjugate())
        self.rows, self.cols = np.triu_indices(len(g))
        differences = g[self.rows] - g[self.cols]
        mapped = [members.get(tuple(d), (-1, 0j)) for d in differences]
        self.indices = np.array([j for j, p in mapped], int)
        self.phases = np.array([p for j, p in mapped], complex)
        self.present = self.indices >= 0
        self.off = self.rows != self.cols
        for array in (
            self.c,
            self.representatives,
            self.rows,
            self.cols,
            self.indices,
            self.phases,
            self.present,
            self.off,
        ):
            array.setflags(write=False)

    def _potential(self, values):
        return checked(
            values,
            (self.nstar,),
            "native interstitial potential",
            real=not self.complex_values,
        )

    def coefficient_matrix(self, values):
        v = self._potential(values)
        upper = np.zeros(len(self.rows), complex)
        # Native HAMilt fills HPANEL(row,col) with WARP(G_row-G_col).
        # Its packed H is conjugated by DumpAssembledMatrices to match the
        # CONJG(Z) coefficients written by TAPEWF. This final conjugation
        # matters for complex potentials even when all star phases are real.
        upper[self.present] = (
            v[self.indices[self.present]] * self.phases[self.present]
        ).conj()
        upper[~self.off] = upper[~self.off].real
        result = np.zeros((len(self.c),) * 2, complex)
        result[self.rows, self.cols] = upper
        result[self.cols[self.off], self.rows[self.off]] = upper[self.off].conj()
        return result

    def action(self, values, retained=None):
        if retained is not None and (
            int(retained) != retained or not 0 < retained < self.nband
        ):
            raise ValueError("invalid retained interstitial band count")
        nr = self.nband if retained is None else int(retained)
        columns = self.c.conj().T @ (self.coefficient_matrix(values) @ self.c[:, :nr])
        return herm(columns) if retained is None else (herm(columns[:nr]), columns[nr:])

    def pullback(self, covector, cross_covector=None):
        y = np.asarray(covector)
        if cross_covector is None:
            checked(y, (self.nband,) * 2, "interstitial band covector", hermitian=True)
            full = (self.c @ y) @ self.c.conj().T
        else:
            if y.ndim != 2 or not 0 < len(y) < self.nband:
                raise ValueError("invalid retained interstitial covector")
            nr = len(y)
            checked(y, (nr, nr), "retained interstitial covector", hermitian=True)
            q = checked(
                cross_covector, (self.nband - nr, nr), "cross interstitial covector"
            )
            cr, cq = self.c[:, :nr], self.c[:, nr:]
            cross = (cq @ q) @ cr.conj().T
            full = (cr @ y) @ cr.conj().T + cross + cross.conj().T
        upper = full[self.rows, self.cols].copy()
        upper[self.off] *= 2
        upper[~self.off] = upper[~self.off].real
        result = np.zeros(self.nstar, complex)
        np.add.at(
            result,
            self.indices[self.present],
            self.phases[self.present].conj() * upper[self.present].conj(),
        )
        return result if self.complex_values else result.real


class NativePotentialCoordinates:
    """Packed native potential sectors; the interstitial part is not density PW.

    The caller must explicitly supply every sector to pack(). Zero arrays may
    be used for an isolated-sector test. The coordinate layout never silently
    substitutes a missing physical operator with a frozen approximation.
    """

    def __init__(
        self, spherical_size, nonspherical_shapes, nstar, *, complex_interstitial
    ):
        self.spherical_size = int(spherical_size)
        if (
            self.spherical_size != spherical_size
            or spherical_size < 1
            or int(nstar) != nstar
            or nstar < 1
        ):
            raise ValueError("invalid potential dimensions")
        if not isinstance(complex_interstitial, (bool, np.bool_)):
            raise ValueError("invalid potential reality flag")
        self.nonspherical_shapes = dict(sorted(nonspherical_shapes.items()))
        self.slices = {}
        pos = self.spherical_size
        for a, shape in self.nonspherical_shapes.items():
            if len(shape) != 2 or any(int(n) != n or n < 1 for n in shape):
                raise ValueError("invalid nonspherical potential shape")
            self.slices[a] = slice(pos, pos + int(np.prod(shape)))
            pos += int(np.prod(shape))
        self.pw_offset = pos
        self.nstar = int(nstar)
        self.complex_interstitial = bool(complex_interstitial)
        self.size = pos + self.nstar * (2 if self.complex_interstitial else 1)

    def pack(self, spherical, nonspherical, interstitial):
        s = checked(spherical, (self.spherical_size,), "spherical potential", real=True)
        if set(nonspherical) != set(self.slices):
            raise ValueError("missing nonspherical atom")
        result = np.empty(self.size)
        result[: self.spherical_size] = s
        for a, sl in self.slices.items():
            result[sl] = checked(
                nonspherical[a],
                self.nonspherical_shapes[a],
                "nonspherical potential",
                real=True,
            ).ravel()
        v = checked(
            interstitial,
            (self.nstar,),
            "native interstitial potential",
            real=not self.complex_interstitial,
        )
        result[self.pw_offset : self.pw_offset + self.nstar] = v.real
        if self.complex_interstitial:
            result[self.pw_offset + self.nstar :] = v.imag
        return result

    def unpack(self, vector):
        v = checked(vector, (self.size,), "native potential coordinates", real=True)
        pw = v[self.pw_offset : self.pw_offset + self.nstar]
        if self.complex_interstitial:
            pw = pw + 1j * v[self.pw_offset + self.nstar :]
        return (
            v[: self.spherical_size],
            {
                a: v[sl].reshape(self.nonspherical_shapes[a])
                for a, sl in self.slices.items()
            },
            pw,
        )


class NativePotentialBands:
    """Joint native input map with complete RR/QR complement output.

    spherical is the existing SphericalPotentialBands. Its nonspherical entries
    supply BOTH their central radial response (inside spherical.action) and
    their explicit potential response here, exactly once each. interstitial
    supplies one NativeInterstitialPotential per k point. Native harmonics
    must exclude (0,0), already represented by the spherical input.

    The nonspherical and IR operations assemble only RR/QR. The spherical
    operation currently supplies a full matrix before extraction; its unused
    QQ assembly remains a separate optimization.
    """

    def __init__(self, spherical, interstitial, *, retained):
        self.spherical = spherical
        self.ir = tuple(interstitial)
        nk, nb, _ = spherical.coordinates.band.shape
        if int(retained) != retained or not 0 < retained < nb:
            raise ValueError("invalid retained count")
        self.nr = int(retained)
        self.nb = nb
        self.nk = nk
        if len(self.ir) != nk or any(op.nband != nb for op in self.ir):
            raise ValueError("inconsistent IR band layout")
        first = self.ir[0]
        if any(
            op.complex_values != first.complex_values
            or not np.array_equal(op.representatives, first.representatives)
            for op in self.ir
        ):
            raise ValueError("inconsistent native potential star layout")
        self.ns = {}
        for a, op, ks in spherical.nonspherical:
            if a in self.ns:
                raise ValueError(
                    "combine each atom's nonspherical entries before preparation"
                )
            if (0, 0) in op.layout.harmonics:
                raise ValueError("spherical harmonic would be counted twice")
            self.ns[a] = (op, ks)
        self.input = NativePotentialCoordinates(
            spherical.input_size,
            {a: op.v.shape for a, (op, ks) in self.ns.items()},
            first.nstar,
            complex_interstitial=first.complex_values,
        )
        matrices = ComplementMatrixCoordinates(
            (nk, self.nr, self.nr), [nb - self.nr] * nk
        )
        self.output = ComplementPerturbationCoordinates(
            matrices, (nk, self.nr, *spherical.coordinates.raw_shape[2:])
        )
        self.operator = PairedOperator(
            self.input.size,
            self.output.size,
            self.action,
            self.pullback,
            "native potential to complete RR/QR H/O/projector",
        )

        def select(v):
            return self.input.unpack(v)[0].copy()

        def select_back(s):
            y = np.zeros(self.input.size)
            y[: spherical.input_size] = checked(
                s, (spherical.input_size,), "spherical covector", real=True
            )
            return y

        self.spherical_selector = PairedOperator(
            self.input.size,
            spherical.input_size,
            select,
            select_back,
            "native spherical input selection",
        )

    def action(self, vector):
        s, ns, pw = self.input.unpack(vector)
        nr = self.nr
        if np.any(s):
            h, o, p = self.spherical.coordinates.unpack(self.spherical.action(s))
            hrr, orr, raw = (
                h[:, :nr, :nr].copy(),
                o[:, :nr, :nr].copy(),
                p[:, :nr].copy(),
            )
            hqr, oqr = list(h[:, nr:, :nr].copy()), list(o[:, nr:, :nr].copy())
        else:
            hrr = np.zeros((self.nk, nr, nr), complex)
            orr = np.zeros_like(hrr)
            hqr = [np.zeros((self.nb - nr, nr), complex) for _ in self.ir]
            oqr = [q.copy() for q in hqr]
            raw = np.zeros(self.output.raw_shape, complex)
        for a, (op, ks) in self.ns.items():
            if np.any(ns[a]):
                for k, (rr, q) in zip(ks, op.potential_action(ns[a], retained=nr)):
                    hrr[k] += rr
                    hqr[k] += q
        if np.any(pw):
            for k, op in enumerate(self.ir):
                rr, q = op.action(pw, retained=nr)
                hrr[k] += rr
                hqr[k] += q
        return self.output.pack(hrr, orr, raw, hqr, oqr)

    def pullback(self, vector):
        hrr, orr, raw, hqr, oqr = self.output.unpack(vector)
        nr = self.nr
        h = np.zeros((self.nk, self.nb, self.nb), complex)
        o = np.zeros_like(h)
        p = np.zeros(self.spherical.coordinates.raw_shape, complex)
        h[:, :nr, :nr] = hrr
        o[:, :nr, :nr] = orr
        p[:, :nr] = raw
        h[:, nr:, :nr] = hqr
        o[:, nr:, :nr] = oqr
        h[:, :nr, nr:] = np.array(hqr).swapaxes(1, 2).conj()
        o[:, :nr, nr:] = np.array(oqr).swapaxes(1, 2).conj()
        sb = self.spherical.pullback(self.spherical.coordinates.pack(h, o, p))
        ns = {
            a: op.potential_pullback([hrr[k] for k in ks], [hqr[k] for k in ks])
            for a, (op, ks) in self.ns.items()
        }
        pw = sum((op.pullback(hrr[k], hqr[k]) for k, op in enumerate(self.ir)))
        return self.input.pack(sb, ns, pw)
