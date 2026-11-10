#!/usr/bin/env python3
"""Complete, opt-in LAPW tangent engine for frozen-Sigma hydrostatic response.

Read a JSON configuration; copy inputs to a NEW workspace; generate two
unrelaxed affine-strain source evaluations; solve the central-cell linear
screening equation.  No impurity solver or strained SCF loop is called.

Each Jacobian action is [F(n0+t*q)-F(n0)]/t at the CENTRAL cell.  F comprises
LAPW0, LAPW1, dmft2 with fixed Sigma and freshly normalized projector, lcore,
and the ordinary interstitial charge normalization.  This finite-difference
action includes radial augmentation and all unoccupied LAPW states without
storing them.  The algorithm is linear response; halving t and the source
strain is mandatory for a quantitative stress result.

Initial supported scope: nonmagnetic, no SOC, primitive orthogonal cells,
projector 5, fixed band indices, fixed total particle number, fixed RMT.
The exported endpoint potentials are NOT yet a general stress tensor.  A
material-coordinate potential derivative and the existing projector adjoint
are needed for the hydrostatic contraction described in the accompanying note.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import struct
import subprocess
import time
from pathlib import Path

import numpy as np
from scipy.io import FortranFile, FortranEOFError
from scipy.interpolate import CubicSpline
from scipy.special import spherical_jn

from frozen_sigma_screening import solve_response
from wien_response_density import Density
from wstruct import Struct
from dmft2_step import local_definition


def integrated_charge(density, structure, *, quadrature="native"):
    """Total charge with explicit native (default) or legacy spline quadrature.

    Native mode follows MIXER CHARGE/SETN. legacy_spline reproduces archived
    reference runs; it must not be silently mixed with a native fast Jacobian.
    """
    s = Struct()
    s.ReadStruct(str(structure))
    if s.lattyp.strip() != "P" or not np.allclose([s.alpha, s.beta, s.gamma], 90.0):
        raise ValueError("initial response driver supports primitive orthogonal cells")
    if quadrature == "native":
        from lapw_charge_response import native_charge_covectors, native_mixer_inversion

        mt, ir = native_charge_covectors(
            density, s, inversion=native_mixer_inversion(s)
        )
        x = density.vector()
        return float(mt @ x), float(ir @ x)
    if quadrature != "legacy_spline":
        raise ValueError("charge_integration must be native or legacy_spline")
    cell = np.array([s.a, s.b, s.c])
    volume = float(np.prod(cell))
    mt = 0.0
    for (atom, ell, m), values in density.rad.items():
        if ell == m == 0:
            grid = density.grids[atom]
            mt += density.mult[atom] * CubicSpline(
                grid, values * np.exp(grid)
            ).integrate(grid[0], grid[-1])
    interstitial = 0j
    for g, value in density.pw.items():
        g = np.asarray(g)
        magnitude = np.linalg.norm(2 * np.pi * g / cell)
        theta = complex(np.all(g == 0))
        position = 0
        for atom in range(s.nat):
            x = magnitude * s.rmt[atom]
            sphere = 1.0 if x == 0 else 3 * spherical_jn(1, x) / x
            for equivalent in range(s.mult[atom]):
                theta -= (
                    4
                    * np.pi
                    * s.rmt[atom] ** 3
                    / (3 * volume)
                    * sphere
                    * np.exp(-2j * np.pi * np.dot(s.pos[position], g))
                )
                position += 1
        # CLMSUM stores a star sum; the sphere characteristic is constant
        # on the star after including symmetry-equivalent atoms.
        interstitial += value * theta.conjugate() * volume
    return float(mt), float(interstitial.real)


def vector_signature(path, nat, window):
    """Check fixed band/G sets and isolated fixed-index window boundaries."""
    records = []
    with FortranFile(path) as stream:
        for _ in range(2 * nat):
            stream.read_record("u1")
        while True:
            try:
                header = stream.read_record("u1").tobytes()
            except FortranEOFError:
                break
            kx, ky, kz, name, nbasis, nbands, weight = struct.unpack(
                "<ddd10siid", header
            )
            labels = stream.read_record("i4").reshape(3, nbasis, order="F")
            energies = []
            for _ in range(nbands):
                index, energy = struct.unpack("<id", stream.read_record("u1").tobytes())
                energies.append(energy)
                stream.read_record("u1")
            first, last = window
            if first < 1 or last > nbands:
                raise ValueError("DMFT band window extends outside the vector file")
            gaps = []
            if first > 1:
                gaps.append(energies[first - 1] - energies[first - 2])
            if last < nbands:
                gaps.append(energies[last] - energies[last - 1])
            if gaps and min(gaps) < 1e-8:
                raise ValueError("DMFT window cuts a degeneracy")
            # The multiset allows harmless reordering of plane waves.
            records.append((kx, ky, kz, nbands, tuple(sorted(map(tuple, labels.T)))))
    if not records:
        raise ValueError("empty vector file")
    return records


class WienResponse:
    def __init__(self, config):
        self.config = config
        self.case = config["case"]
        self.reference = Path(config["reference"]).resolve()
        self.work = Path(config["workspace"]).resolve()
        self.work.mkdir(parents=True)  # Never overwrite an earlier calculation.
        self.structure = self.reference / (self.case + ".struct")
        self.rho = Density(self.reference / (self.case + ".clmsum"), self.structure)
        self.charge_integration = config.get("charge_integration", "native")
        # Store the resolved choice so the default is visible in saved runs.
        self.config = dict(config, charge_integration=self.charge_integration)
        config = self.config
        mt, ir = integrated_charge(
            self.rho, self.structure, quadrature=self.charge_integration
        )
        self.target_charge = round(mt + ir)
        if abs(mt + ir - self.target_charge) > 1e-3:
            raise ValueError("central total charge is not sufficiently normalized")
        mixer = (
            (self.reference / (self.case + ".inm")).read_text().splitlines()[0].split()
        )
        if len(mixer) < 3 or float(mixer[1]) != 0.0 or mixer[2].upper() != "YES":
            raise ValueError(
                "initial driver requires zero background charge and default mixer normalization"
            )
        s = Struct()
        s.ReadStruct(str(self.structure))
        self.nat = s.nat
        self.cell = np.array([s.a, s.b, s.c])
        header = (
            (self.reference / (self.case + ".indmfl"))
            .read_text()
            .splitlines()[0]
            .split()
        )
        if int(header[2]) != 1 or abs(int(header[3])) != 5:
            raise ValueError("response requires normalized projector 5")
        self.window = tuple(map(int, header[:2]))
        self.wien = Path(config["wien_root"])
        self.dmft2 = str(Path(config["dmft2"]).resolve())
        self.mpi = config.get("mpi", [])
        self.step = float(config.get("probe_step", 1e-5))
        self.strain = float(config.get("strain", 0.00025))
        if not (0 < self.step < 0.01 and 0 < self.strain < 0.01):
            raise ValueError("positive, small probe and strain steps are required")
        self.env = os.environ.copy()
        for key in list(self.env):
            if key.startswith(("EDMFT_", "WIEN_STRESS_", "WIEN_LAPW1_", "DFTKS_")):
                del self.env[key]
        self.env.update(
            WIENROOT=str(self.wien),
            OMP_NUM_THREADS="1",
            OPENBLAS_NUM_THREADS="1",
            VECLIB_MAXIMUM_THREADS="1",
        )
        self.signature = None
        self.timings = []
        (self.work / "config.json").write_text(json.dumps(config, indent=2) + "\n")

    def command(self, folder, label, argv):
        started = time.monotonic()
        logpath = folder / (label + ".log")
        with logpath.open("w") as log:
            result = subprocess.run(
                argv, cwd=folder, env=self.env, stdout=log, stderr=subprocess.STDOUT
            )
        text = logpath.read_text(errors="replace")
        if result.returncode or any(
            token in text
            for token in [
                "Fortran runtime error",
                "Error termination",
                "ERROR STOP",
                "error: command",
            ]
        ):
            raise RuntimeError(f"{label} failed: see {logpath}")
        if label in {"dmft2", "lapw1"} and f"{label.upper()} END" not in text:
            raise RuntimeError(
                f"{label} did not report successful completion: {logpath}"
            )
        for error_file in folder.glob(label + ".error*"):
            if error_file.read_text(errors="replace").strip():
                raise RuntimeError(f"{label} wrote an error: {error_file}")
        self.timings.append(
            {
                "folder": folder.name,
                "program": label,
                "seconds": time.monotonic() - started,
            }
        )

    def stage(self, name, scale=1.0):
        folder = self.work / name
        folder.mkdir()
        suffixes = {
            ".struct",
            ".in0",
            ".in1",
            ".in1c",
            ".in2",
            ".in2c",
            ".inc",
            ".inm",
            ".klist",
            ".kgen",
            ".rotlm",
            ".indmfl",
            ".indmfi",
            ".indmf2",
        }
        for source in self.reference.iterdir():
            if source.is_file() and (
                (source.name.startswith(self.case + ".") and source.suffix in suffixes)
                or source.name in ["EF.dat", "projectorw.dat", "Eorb_imp.dat"]
                or re.fullmatch(r"sig\.inp\d*", source.name)
            ):
                shutil.copy2(source, folder / source.name)
        for name, key in [
            ("lapw1.def", "lapw1_definition"),
            ("dmft2.def", "dmft2_definition"),
        ]:
            text = local_definition(Path(self.config[key]).read_text())
            (folder / name).write_text(text)
        for suffix in ["in0", "in2", "in2c"]:
            path = folder / (self.case + "." + suffix)
            if path.exists():
                lines = path.read_text().splitlines()
                lines[0] = re.sub(r"^\s*\S+", "TOT ", lines[0], count=1)
                if suffix == "in0":
                    lines[1] = lines[1].replace("NR2V", "R2V ")
                path.write_text("\n".join(lines) + "\n")
        path = folder / (self.case + ".indmf2")
        lines = path.read_text().splitlines()
        first = next(
            i
            for i, line in enumerate(lines)
            if line.strip() and not line.startswith("#")
        )
        mu = (folder / "EF.dat").read_text().split()[0]
        lines[first] = f"c {mu} 1 1.0 .T. 0.02"
        path.write_text("\n".join(lines) + "\n")
        path = folder / (self.case + ".struct")
        lines = path.read_text().splitlines()
        lines[3] = "".join(f"{a*scale:10.6f}" for a in self.cell) + lines[3][30:]
        path.write_text("\n".join(lines) + "\n")
        actual_cell = np.array([float(lines[3][i : i + 10]) for i in [0, 10, 20]])
        actual_scale = float(np.mean(actual_cell / self.cell))
        return folder, actual_scale

    def normalize(self, density, structure):
        mt, ir = integrated_charge(
            density, structure, quadrature=self.charge_integration
        )
        if ir <= 0:
            raise ValueError("nonpositive interstitial charge")
        for g in density.pw:
            density.pw[g] *= (self.target_charge - mt) / ir
        return density

    def evaluate(self, name, density, scale=1.0, potential_only=False):
        folder, scale = self.stage(name, scale)
        density.transport(scale).write(folder / (self.case + ".clmsum"))
        self.command(folder, "lapw0", [str(self.wien / "x"), "lapw0", "-f", self.case])
        if potential_only:
            return folder, scale
        # A zero-length in1c can be left by opening an optional unit.  File
        # existence is not a reliable real/complex selector; use the actual
        # LAPW1 unit-5 input named by its definition file.
        definition = (folder / "lapw1.def").read_text()
        input_match = re.search(r"(?m)^\s*5\s*,\s*'([^']+)'", definition)
        if input_match is None:
            raise ValueError("LAPW1 definition has no unit-5 input")
        lapw1 = "lapw1c" if input_match[1].endswith(".in1c") else "lapw1"

        fixed = self.work / "central/lapw1.fixgvec"
        if name == "central":
            (folder / "lapw1.fixgvec.write").touch()
        else:
            if not fixed.exists():
                raise ValueError("LAPW1 lacks the required fixed-G reference support")
            shutil.copy2(fixed, folder / "lapw1.fixgvec")
            (folder / "lapw1.fixgvec.read").touch()
        self.command(folder, "lapw1", [str(self.wien / lapw1), "lapw1.def"])
        signature = vector_signature(
            folder / (self.case + ".vector"), self.nat, self.window
        )
        if self.signature is None:
            self.signature = signature
        elif self.signature != signature:
            raise ValueError(
                "G set, band count, or k points changed; response is not differentiable on this basis"
            )
        self.command(folder, "dmft2", self.mpi + [self.dmft2, "dmft2.def"])
        self.command(folder, "lcore", [str(self.wien / "x"), "lcore", "-f", self.case])
        structure = folder / (self.case + ".struct")
        total = self.rho.zeros().add(
            Density(folder / (self.case + ".clmval"), structure)
        )
        total.add(Density(folder / (self.case + ".clmcor"), structure))
        self.normalize(total, structure)
        total.write(folder / (self.case + ".clmout"))
        pulled = self.normalize(total.transport(scale, True), self.structure)
        pulled.write(folder / (self.case + ".qout"))
        return pulled, scale

    def solve(self):
        base, _ = self.evaluate("central", self.rho)
        if "source_density" in self.config:
            source = Density(self.config["source_density"], self.structure)
        else:
            plus, lp = self.evaluate("source_plus", self.rho, 1 + self.strain)
            minus, lm = self.evaluate("source_minus", self.rho, 1 - self.strain)
            source = plus.add(minus, -1).scale(1 / (lp - lm))
        source.write(self.work / "source.density")
        counter = 0

        def action(vector):
            nonlocal counter
            direction = self.rho.from_vector(vector)
            out, _ = self.evaluate(
                f"action_{counter:03d}", self.rho.copy().add(direction, self.step)
            )
            counter += 1
            return out.add(base, -1).scale(1 / self.step).vector()

        def checkpoint(vector, residual, row):
            self.rho.from_vector(vector).write(self.work / "response.density")
            np.save(self.work / "residual.npy", residual)
            print(json.dumps(row), flush=True)

        initial = None
        if "initial_density" in self.config:
            initial = Density(self.config["initial_density"], self.structure).vector()
        answer, report = solve_response(
            source.vector(),
            action,
            initial=initial,
            absolute_tolerance=float(self.config.get("tolerance", 0.02)),
            max_iterations=int(self.config.get("max_iterations", 25)),
            callback=checkpoint,
        )
        q = self.rho.from_vector(answer)
        q.write(self.work / "response.density")
        report.update(
            case=self.case,
            central_charge_residual_norm=float(
                np.linalg.norm(base.copy().add(self.rho, -1).vector())
            ),
            ensemble="fixed N",
            fixed_G=True,
            kpoints=len(self.signature),
            band_window=self.window,
            stored_bands=sorted({record[3] for record in self.signature}),
            source_strain=self.strain,
            probe_step=self.step,
            potential_convention="physical endpoint potentials; material derivative requires pullback",
        )
        report["timings"] = self.timings
        (self.work / "report.json").write_text(json.dumps(report, indent=2) + "\n")
        if report["converged"]:
            # Export physical endpoint potentials from the predicted material
            # density q0 +/- h*q1.  These require only two LAPW0 evaluations.
            self.evaluate(
                "potential_plus",
                self.rho.copy().add(q, self.strain),
                1 + self.strain,
                True,
            )
            self.evaluate(
                "potential_minus",
                self.rho.copy().add(q, -self.strain),
                1 - self.strain,
                True,
            )
        return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", type=Path)
    args = parser.parse_args()
    report = WienResponse(json.loads(args.config.read_text())).solve()
    return 0 if report["converged"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
