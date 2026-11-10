#!/usr/bin/env python3
"""Finite-basis diagnostic for the interstitial LAPW overlap derivative.

For a hydrostatic strain ``s`` with fixed muffin-tin radii, the plane-wave
part of the LAPW overlap is the Fourier transform of the interstitial step
function,

    O^I_GG'(s) = Theta_{G-G'}(s).

This program evaluates ``dTheta/ds`` directly from the central WIEN2k
structure and contracts it with the central dmft2 ``rho`` and ``rhoE``
matrices.  It is a decomposition diagnostic only: it neither approximates
the augmentation derivative nor changes an analytic stress calculation.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
import re

import numpy as np

from check_oneparticle_stress_fd import files, read_one_particle


@dataclass(frozen=True)
class Atom:
    positions: np.ndarray
    rmt: float


def read_br1(path: Path) -> np.ndarray:
    lines = path.read_text(errors="replace").splitlines()
    for index, line in enumerate(lines):
        if line.strip() == "BR1=":
            try:
                return np.array(
                    [[float(value) for value in lines[index + row + 1].split()[:3]]
                     for row in range(3)],
                    dtype=float,
                )
            except (IndexError, ValueError) as error:
                raise ValueError(f"malformed BR1 block in {path}") from error
    raise ValueError(f"no BR1 block in {path}")


def _position(line: str) -> np.ndarray:
    match = re.search(
        r"X=\s*([-+0-9.Ee]+)\s+Y=\s*([-+0-9.Ee]+)\s+Z=\s*([-+0-9.Ee]+)",
        line,
    )
    if match is None:
        raise ValueError(f"cannot parse atomic position from {line!r}")
    return np.array([float(match.group(i)) for i in range(1, 4)])


def read_atoms(path: Path) -> list[Atom]:
    lines = path.read_text(errors="replace").splitlines()
    atoms: list[Atom] = []
    index = 4
    while index < len(lines):
        if not lines[index].lstrip().startswith("ATOM"):
            index += 1
            continue
        positions = [_position(lines[index])]
        index += 1
        match = re.search(r"MULT=\s*(\d+)", lines[index])
        if match is None:
            raise ValueError(f"missing MULT after atom in {path}")
        multiplicity = int(match.group(1))
        index += 1
        for _ in range(multiplicity - 1):
            positions.append(_position(lines[index]))
            index += 1
        match = re.search(r"RMT=\s*([-+0-9.Ee]+)", lines[index])
        if match is None:
            raise ValueError(f"missing RMT after atom in {path}")
        atoms.append(Atom(np.array(positions), float(match.group(1))))
        index += 4  # atom record followed by the three local-rotation rows
    if not atoms:
        raise ValueError(f"no atoms found in {path}")
    return atoms


def theta(labels: np.ndarray, br1: np.ndarray, atoms: list[Atom],
          volume: float, scale: float) -> np.ndarray:
    """Return Theta_GG' at hydrostatic scale ``1+scale``.

    Fractional atomic coordinates and therefore the phase factors remain
    invariant.  Reciprocal vectors scale as ``1/(1+scale)`` and the cell
    volume as ``(1+scale)^3``; physical RMT values are held fixed.
    """
    difference = labels[:, :, None] - labels[:, None, :]
    cartesian = np.einsum("ab,bmn->amn", br1 / (1.0 + scale), difference)
    norm = np.sqrt(np.sum(cartesian * cartesian, axis=0))
    strained_volume = volume * (1.0 + scale) ** 3
    result = np.zeros(norm.shape, dtype=np.complex128)
    zero = norm < 1.0e-13
    sphere_fraction = 0.0
    for atom in atoms:
        sphere_volume = 4.0 * np.pi * atom.rmt ** 3 / 3.0
        sphere_fraction += len(atom.positions) * sphere_volume / strained_volume
        q = norm * atom.rmt
        radial = np.empty_like(q)
        radial[zero] = -1.0 / 3.0
        nonzero = ~zero
        radial[nonzero] = (
            q[nonzero] * np.cos(q[nonzero]) - np.sin(q[nonzero])
        ) / q[nonzero] ** 3
        phase = np.zeros(norm.shape, dtype=np.complex128)
        for position in atom.positions:
            argument = 2.0 * np.pi * np.einsum(
                "a,amn->mn", position, difference
            )
            phase += np.exp(-1j * argument)
        result += 3.0 * sphere_volume / strained_volume * radial * phase
    result[zero] = 1.0 - sphere_fraction
    return result


def theta_derivative(labels: np.ndarray, br1: np.ndarray,
                     atoms: list[Atom], volume: float) -> np.ndarray:
    """Analytic ``dTheta/d epsilon_ab`` at fixed spherical RMT.

    Reciprocal vectors obey ``G_,ab=-e_a G_b`` (symmetrized for shear),
    so ``|G|_,ab=-G_a G_b/|G|``.  Fractional atomic positions make the
    phase invariant.  The remaining terms are the derivative of the sphere
    form factor and the explicit ``1/Omega`` normalization.
    """
    difference = labels[:, :, None] - labels[:, None, :]
    cartesian = np.einsum("ab,bmn->amn", br1, difference)
    norm = np.sqrt(np.sum(cartesian * cartesian, axis=0))
    result = np.zeros((3, 3) + norm.shape, dtype=np.complex128)
    zero = norm < 1.0e-13
    sphere_fraction = 0.0
    for atom in atoms:
        sphere_volume = 4.0 * np.pi * atom.rmt ** 3 / 3.0
        sphere_fraction += len(atom.positions) * sphere_volume / volume
        q = norm * atom.rmt
        radial = np.empty_like(q)
        radial_prime = np.empty_like(q)
        radial[zero] = -1.0 / 3.0
        radial_prime[zero] = 0.0
        nonzero = ~zero
        radial[nonzero] = (
            q[nonzero] * np.cos(q[nonzero]) - np.sin(q[nonzero])
        ) / q[nonzero] ** 3
        radial_prime[nonzero] = (
            (3.0 - q[nonzero] ** 2) * np.sin(q[nonzero])
            - 3.0 * q[nonzero] * np.cos(q[nonzero])
        ) / q[nonzero] ** 4
        phase = np.zeros(norm.shape, dtype=np.complex128)
        for position in atom.positions:
            argument = 2.0 * np.pi * np.einsum(
                "a,amn->mn", position, difference
            )
            phase += np.exp(-1j * argument)
        prefactor = 3.0 * sphere_volume / volume
        central = prefactor * radial * phase
        for alpha in range(3):
            for beta in range(3):
                q_derivative = np.zeros_like(norm)
                q_derivative[nonzero] = (
                    -atom.rmt
                    * cartesian[alpha, nonzero]
                    * cartesian[beta, nonzero]
                    / norm[nonzero]
                )
                result[alpha, beta] += (
                    prefactor * radial_prime * q_derivative * phase
                    - (1.0 if alpha == beta else 0.0) * central
                )
    for alpha in range(3):
        for beta in range(3):
            result[alpha, beta, zero] = (
                sphere_fraction if alpha == beta else 0.0
            )
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--moments", required=True, type=Path)
    parser.add_argument("--struct", required=True, type=Path)
    parser.add_argument("--output1", required=True, type=Path)
    parser.add_argument("--strain", required=True, type=float)
    parser.add_argument("--volume", required=True, type=float)
    args = parser.parse_args()
    if args.strain <= 0.0 or args.volume <= 0.0:
        parser.error("--strain and --volume must be positive")

    br1 = read_br1(args.output1)
    atoms = read_atoms(args.struct)
    totals = {"rho": 0.0, "rhoe": 0.0, "energy": 0.0, "msigma": 0.0}
    central_msigma_interstitial = 0.0
    central_msigma_total = 0.0
    hermiticity = 0.0
    analytic_difference = 0.0
    moment_files = files(args.moments)

    for path in moment_files.values():
        reference = read_one_particle(path)
        npw = reference.nbasis - reference.nnlo
        labels = reference.labels[:, :npw]
        derivative = (
            theta(labels, br1, atoms, args.volume, args.strain)
            - theta(labels, br1, atoms, args.volume, -args.strain)
        ) / (2.0 * args.strain)
        analytic_tensor = theta_derivative(labels, br1, atoms, args.volume)
        analytic = sum(analytic_tensor[axis, axis] for axis in range(3))
        analytic_difference = max(
            analytic_difference, float(np.max(np.abs(derivative - analytic)))
        )
        hermiticity = max(
            hermiticity, float(np.max(np.abs(derivative - derivative.conj().T)))
        )
        band_derivative = np.zeros((reference.nfull, reference.nfull), complex)
        for component in range(reference.ncomponents):
            coefficients = reference.coefficients[:npw, :, component]
            band_derivative += coefficients.conj().T @ derivative @ coefficients

        energy_moment = 0.5 * (
            reference.energies[:, None] + reference.energies[None, :]
        ) * reference.rho
        msigma = reference.rhoe - energy_moment
        totals["rho"] += float(np.real(np.trace(reference.rho @ band_derivative)))
        totals["rhoe"] += float(np.real(np.trace(reference.rhoe @ band_derivative)))
        totals["energy"] += float(np.real(np.trace(energy_moment @ band_derivative)))
        totals["msigma"] += float(np.real(np.trace(msigma @ band_derivative)))
        central_theta = theta(labels, br1, atoms, args.volume, 0.0)
        central_theta_band = np.zeros_like(band_derivative)
        for component in range(reference.ncomponents):
            coefficients = reference.coefficients[:npw, :, component]
            central_theta_band += (
                coefficients.conj().T @ central_theta @ coefficients
            )
        central_msigma_interstitial += float(
            np.real(np.trace(msigma @ central_theta_band))
        )
        central_msigma_total += float(np.real(np.trace(msigma)))

    conversion = 14710.5164057 / (3.0 * args.volume)
    print(f"k points                                      {len(moment_files)}")
    print(f"max Hermiticity error                         {hermiticity:.12e}")
    print(f"max hydro analytic-vs-FD error                {analytic_difference:.12e}")
    for name in ("rho", "rhoe", "energy", "msigma"):
        value = totals[name]
        print(f"Tr[{name:7s} dTheta/ds]                        {value:+.12e} Ry")
        print(f"hydro mean {name:7s}                          {value*conversion:+.9f} GPa")
    print(f"closure rhoE-energy-MSigma                    "
          f"{totals['rhoe']-totals['energy']-totals['msigma']:+.12e} Ry")
    print(f"Tr[MSigma Theta] at center                    "
          f"{central_msigma_interstitial:+.12e} Ry")
    print(f"Tr[MSigma O] at center                        "
          f"{central_msigma_total:+.12e} Ry")
    print(f"Tr[MSigma O_aug] at center                    "
          f"{central_msigma_total-central_msigma_interstitial:+.12e} Ry")


if __name__ == "__main__":
    main()
