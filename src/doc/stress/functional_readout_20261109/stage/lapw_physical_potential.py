"""Fixed-basis physical-potential matrix used in functional cancellation.

P v changes the physical Hamiltonian only. Radial functions, matching, metric,
raw projectors and geometry are held fixed. In contrast, the installed B v
also differentiates the potential-dependent radial basis and raw projectors.
Both maps use the same native H quadrature and Fourier conventions.

P^T rho is the matrix-defined ordinary potential covector. Identifying it with
the spatial-density term of a discrete Klein implementation still requires a
quadrature/normalization/core audit. This module performs no such cancellation
and publishes no stress correction.
"""

import numpy as np

from lapw_density_response import checked
from stress_response_adjoint import PairedOperator


class FixedBasisPhysicalPotential:
    """Native potential coordinates -> joint RR/QR H-only perturbations.

    Construct from a compact NativePotentialBands. The spherical kernel's
    difference_dot means delta[r(V_physical-V_basis)]; with zero radial
    tangents it supplies the ordinary fixed-basis potential integral. The
    other sectors already have explicit physical-potential actions.

    The reverse selects the matching difference_bar, not the radial covector.
    RR+QR+RQ pairing is inherited from the tested rectangular native kernels.
    There is no overlap/raw derivative in P. Feeding P into the electronic
    window-response map WOULD move the window; only Tr(rho P v) is the ordinary
    physical-potential term to be subtracted from a full lattice derivative.
    """

    def __init__(self, native_bands):
        if not native_bands.compact:
            raise NotImplementedError("requires compact RR/QR native bands")
        self.bands = native_bands
        self.input = native_bands.input
        self.output = native_bands.output
        self.operator = PairedOperator(
            self.input.size,
            self.output.size,
            self.action,
            self.pullback,
            "fixed-basis physical potential to H only",
        )

    def action(self, vector):
        bands = self.bands
        s, ns, pw = self.input.unpack(vector)
        nr = bands.nr
        h = np.zeros((bands.nk, nr, nr), complex)
        hq = [np.zeros((n - nr, nr), complex) for n in bands.sizes]
        for a, l, k, op in bands.spherical.spherical:
            dh, _ = op.action(
                np.zeros_like(op.phi),
                np.zeros_like(op.lo),
                difference_dot=s[bands.spherical.slices[a]],
            )
            h[k] += dh[:nr]
            hq[k] += dh[nr:]
        for a, (op, ks) in bands.ns.items():
            for k, (rr, qr) in zip(ks, op.potential_action(ns[a], retained=nr)):
                h[k] += rr
                hq[k] += qr
        for k, op in enumerate(bands.ir):
            rr, qr = op.action(pw, retained=nr)
            h[k] += rr
            hq[k] += qr
        return self.output.pack(
            h,
            np.zeros_like(h),
            np.zeros(self.output.raw_shape, complex),
            hq,
            [np.zeros_like(q) for q in hq],
        )

    def pullback(self, vector):
        bands = self.bands
        h, _, _, hq, _ = self.output.unpack(vector)
        columns = [np.concatenate((rr, qr)) for rr, qr in zip(h, hq)]
        spherical = np.zeros(bands.spherical.input_size)
        for a, l, k, op in bands.spherical.spherical:
            *_, physical = op.pullback(columns[k], np.zeros_like(columns[k]))
            spherical[bands.spherical.slices[a]] += physical
        nonspherical = {
            a: op.potential_pullback([h[k] for k in ks], [hq[k] for k in ks])
            for a, (op, ks) in bands.ns.items()
        }
        interstitial = sum(op.pullback(h[k], hq[k]) for k, op in enumerate(bands.ir))
        return self.input.pack(spherical, nonspherical, interstitial)

    def density_gradient(self, rho):
        """P^T rho in native potential coordinates; weights are already in rho.

        Compose the result with state.potential_selector.rmatvec to obtain the
        potential/core field-coordinate covector. Central density has only RR
        blocks under the prepared empty-Q policy; no projector gradient enters.
        """
        b = self.bands
        rho = checked(rho, (b.nk, b.nr, b.nr), "weighted density")
        for matrix in rho:
            checked(matrix, (b.nr, b.nr), "weighted density block", hermitian=True)
        qr = [np.zeros((n - b.nr, b.nr), complex) for n in b.sizes]
        covector = self.output.pack(
            rho, np.zeros_like(rho), np.zeros(self.output.raw_shape, complex), qr, qr
        )
        return self.pullback(covector)
