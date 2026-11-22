"""Stationary finite-temperature DFT electronic response in the KS gauge.

The same real-Frobenius conventions as the eDMFT lattice response are used.
There are no local orbitals, projector window, self-energy, or Matsubara sums.
The fixed-N chemical-potential derivative is imposed once for the full k sum.
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy.special import expit

from lapw_complement_response import EmptyComplementResponse
from lapw_density_response import checked


def herm(a):
    return (a + a.swapaxes(-1, -2).conj()) / 2


class DFTFermiModel:
    """Frechet derivative of the Fermi matrix function in a moving LAPW basis.

    The scalar is the Mermin grand potential.  In the orthonormal central KS
    gauge H=diag(eps), O=I and rho=2 w_k f(H-mu).  An overlap derivative is
    applied by the same symmetric orthonormalization used in the DMFT model.
    """
    electronic_model = 'dft'

    def __init__(self, energies, weights, mu, temperature, *, parallel):
        self.eps = np.asarray(energies, float)
        self.kw = np.asarray(weights, float)
        self.parallel = parallel
        if self.eps.ndim != 2 or self.kw.shape != (len(self.eps),):
            raise ValueError("DFT k-point energies and weights disagree")
        if not np.all(np.isfinite(self.eps)) or not np.all(np.isfinite(self.kw)):
            raise ValueError("nonfinite DFT spectrum or weights")
        if abs(parallel.sum(float(self.kw.sum())) - 1) > 1e-10:
            raise ValueError("DFT k-point weights must sum to one globally")
        if not np.isfinite(mu) or not np.isfinite(temperature) or temperature <= 0:
            raise ValueError("DFT response requires finite mu and positive temperature")
        self.meta = SimpleNamespace(
            mu=float(mu), temperature=float(temperature),
            omega=np.empty(0), sigma_inf=np.empty((0, 0, 0), complex),
            sigma=np.empty((0, 0, 0, 0), complex),
            sym=lambda x: np.asarray(x),
        )
        nk, nr = self.eps.shape
        # The shared packed-coordinate type requires a nonempty raw leg.
        # This single zero column is inert; it is never a physical projector.
        self.raw = np.zeros((nk, nr, 1, 1), complex)
        self.W = np.empty(0, int)
        self._f = expit((mu - self.eps) / temperature)
        self._slope = -self._f * (1 - self._f) / temperature
        self._density = 2 * self.kw[:, None, None] * np.array(
            [np.diag(row) for row in self._f], complex
        )
        self._divided = np.empty((nk, nr, nr))
        for k, e in enumerate(self.eps):
            gap = e[:, None] - e[None, :]
            near = abs(gap) < 1e-8
            np.divide(
                self._f[k, :, None] - self._f[k, None, :], gap,
                out=self._divided[k], where=~near,
            )
            mid = expit((mu - (e[:, None] + e[None, :]) / 2) / temperature)
            self._divided[k, near] = (-mid * (1 - mid) / temperature)[near]
        self._charge_slope = float(parallel.sum(np.sum(2 * self.kw[:, None] * self._slope)))

    def central_density(self):
        return self._density.copy()

    def _matrix(self, value, name):
        value = checked(value, self._density.shape, name)
        if np.linalg.norm(value - value.swapaxes(1, 2).conj()) > 1e-11 * max(
            1.0, np.linalg.norm(value)
        ):
            raise ValueError(name + " must be Hermitian at every k point")
        return value

    def generalized_response(self, dh, do, raw_derivative=None, fixed_number=True):
        dh = self._matrix(dh, "DFT dH")
        do = self._matrix(do, "DFT dO")
        if raw_derivative is not None and np.asarray(raw_derivative).shape != self.raw.shape:
            raise ValueError("DFT has no raw projector")
        h = dh - 0.5 * (
            do * self.eps[:, None, :] + self.eps[:, :, None] * do
        )
        chem = 0.0
        if fixed_number:
            if abs(self._charge_slope) < 1e-14:
                raise ValueError("fixed-N DFT response has vanishing compressibility")
            chem = float(self.parallel.sum(np.sum(
                2 * self.kw[:, None] * self._slope * np.diagonal(h, axis1=1, axis2=2).real
            )) / self._charge_slope)
        direction = h - chem * np.eye(h.shape[1])[None]
        density = np.asarray(
            2 * self.kw[:, None, None] * self._divided * direction,
            complex,
        )
        density -= 0.5 * (do @ self._density + self._density @ do)
        density = herm(density)
        number = float(self.parallel.sum(np.trace(
            density + self._density @ do, axis1=1, axis2=2
        ).real.sum()))
        return density, dict(dmu_Ry=chem, dN=number, model="Mermin DFT")

    def generalized_adjoint(self, y, fixed_number=True):
        y = self._matrix(y, "DFT density covector")
        hbar = 2 * self.kw[:, None, None] * self._divided * y
        if fixed_number:
            if abs(self._charge_slope) < 1e-14:
                raise ValueError("fixed-N DFT adjoint has vanishing compressibility")
            slope = np.array([
                np.diag(2 * self.kw[k] * self._slope[k]) for k in range(len(self.kw))
            ])
            contraction = float(self.parallel.sum(np.einsum("kij,kji->", y, slope).real))
            hbar -= slope * (contraction / self._charge_slope)
        obar = -0.5 * (
            hbar * self.eps[:, None, :] + self.eps[:, :, None] * hbar
        ) - 0.5 * (y @ self._density + self._density @ y)
        return herm(hbar), herm(obar), np.zeros_like(self.raw)

    def window_gap_report(self):
        return dict(policy="no DMFT window", omitted_pairs=0)

    def window_gap_warning(self):
        return ""

    def cache_info(self):
        return dict(model="Mermin DFT", retained_bands=self.eps.shape[1])


def native_fermi_parameters(directory, case, *, complex_vectors=None):
    """Read the precise native Fermi energy and the requested Mermin smearing.

    TEMP uses only half the entropy in WIEN's total energy.  Its numerical
    derivative is not the stationary Mermin derivative, so this first DFT
    response path requires TEMPS in every central and strained calculation.
    """
    directory = Path(directory)
    if complex_vectors is None:
        complex_vectors = (directory / (case + ".in1c")).is_file()
    in2 = directory / (case + (".in2c" if complex_vectors else ".in2"))
    fields = in2.read_text().splitlines()[2].split()
    if len(fields) < 2 or fields[0].upper() != "TEMPS":
        raise ValueError("DFT response requires TEMPS in case.in2[c]")
    parameter = float(fields[1])
    temperature = parameter - np.floor(parameter)
    if not 1e-7 < temperature < 0.1:
        raise ValueError("TEMPS requires an explicit temperature between 1e-7 and 0.1 Ry")
    weight = (directory / (case + ".weight")).read_text().splitlines()[0].split()
    if len(weight) < 3 or not weight[2].upper().startswith("TEMP"):
        raise ValueError("native LAPW2 weight file does not use temperature broadening")
    # FERMI5 writes its upper nonzero-weight cutoff, emu+30*T, into the
    # first WEIGHT header field.  The actual chemical potential is emu.
    mu = float(weight[0]) - 30 * temperature
    if not np.isfinite(mu):
        raise ValueError("nonfinite native DFT chemical potential")
    return mu, temperature


class DFTComplementResponse(EmptyComplementResponse):
    """Use the shared finite-LAPW complement with a projector-free gradient."""

    def functional_gradients(self):
        from lapw_complement_response import (
            ComplementMatrixCoordinates, ComplementPerturbationCoordinates,
        )

        coordinates = ComplementPerturbationCoordinates(
            ComplementMatrixCoordinates(
                (len(self.model.kw), self.model.eps.shape[1], self.model.eps.shape[1]),
                [s.nq for s in self.solvers],
            ),
            self.model.raw.shape,
        )
        rho = self.model.central_density()
        eps = self.model.eps
        overlap = -0.5 * (rho * eps[:, None, :] + eps[:, :, None] * rho)
        zero_qr = [np.zeros((s.nq, s.nr), complex) for s in self.solvers]
        direct = coordinates.pack(rho, overlap, np.zeros_like(self.model.raw), zero_qr, zero_qr)
        zero = np.zeros_like(direct)
        return dict(projector=zero, direct=direct, lattice=direct), coordinates
