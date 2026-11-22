#!/usr/bin/env python3
"""Screened SOC force reference from directional differences of the native charge map.

For one symmetry-compatible displacement u, M(t,q) performs a single native
LAPW0 / scalar LAPW1 / LAPWSO / frozen-Sigma DMFT2 / LCORE pass from n+t*q
at R+t*u. The response q satisfies q = [M(h,q)-M(-h,q)]/(2*h) to O(h^2).
The same two passes return the corresponding Klein directional derivative.
This is a forward verification method, expensive per mode. The central input
must be a converged fixed-Sigma charge. Self-energy and radial input policy
are held fixed, whereas the normalized projected subspace and KS field follow
each perturbation. The native spinor calculation includes spin flips.
"""

import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

import numpy as np
from dmft2_step import entries, local_definition
from force_response_modes import displaced_struct, internal_modes, check_force_frame
from frozen_sigma_screening import solve_response
from lapw_charge_response import NativeChargeNormalization, native_mixer_inversion
from lapw_dft_vector import _header, _record
from lapw_soc_rank import enforce_ranks, parent_lo_filter
from lapw_lattice import NativeLattice
from wien_response_density import Density
from wstruct import Struct


def _core_energy(folder, case, structure):
    # SCFC prints every (n,kappa) in INC order, once per inequivalent site.
    # INC gives its actual occupation, including the 2j+1 degeneracy.
    lines = iter((folder / (case + ".inc")).read_text().splitlines())
    records = {}
    for line in (folder / (case + ".scfc")).read_text().splitlines():
        match = re.fullmatch(
            r":(\d+[A-Z]+)\s*(\d{3}):\s*\S+\s+([-+\d.Ee]+)\s+Ry\s*", line
        )
        if match:
            records.setdefault(int(match[2]), []).append((match[1], float(match[3])))
    terms = []
    for site, mult in enumerate(structure.mult, 1):
        count = int(next(lines).split()[0])
        states = records.get(site, [])
        if len(states) != count:
            raise ValueError(f"SCFC and INC core-state counts differ at atom {site}")
        for label, energy in states:
            n, kappa, occ = re.split(r"[ ,]+", next(lines).strip())[:3]
            n, kappa = int(n), int(kappa)
            ell = kappa if kappa > 0 else -kappa - 1
            letter = "SPDFGHI"[ell]
            if label != str(n) + letter * (2 if kappa > 0 else 1):
                raise ValueError(f"SCFC and INC order differs at atom {site}, {label}")
            terms.append(mult * float(occ) * energy)
    return math.fsum(terms)


def _scalar(path, tag):
    vals = [
        float(line.split("=")[1].split()[0])
        for line in path.read_text().splitlines()
        if line.startswith(tag)
    ]
    if len(vals) != 1 or not np.isfinite(vals[0]):
        raise ValueError(f"expected one finite {tag} in {path}")
    return vals[0]


def _signature(path, nat, window=None):
    import hashlib
    import struct

    records = []
    with path.open("rb") as f:
        for _ in range(2 * nat):
            _record(f, skip=True)
        while True:
            header = _header(_record(f))
            if header is None:
                break
            kx, ky, kz, _, size, bands, weight = header
            labels = np.frombuffer(_record(f), "<i4").reshape(size, 3)
            labels = labels[np.lexsort(labels.T)]
            edge = {}
            for i in range(1, bands + 1):
                record = _record(f)
                index, energy = struct.unpack("<id", record)
                if index != i or not np.isfinite(energy):
                    raise ValueError(f"invalid band record in {path}")
                if window and i in (window[0] - 1, window[0], window[1], window[1] + 1):
                    edge[i] = energy
                _record(f, skip=True)
            if window:
                first, last = window
                if not 1 <= first <= last <= bands:
                    raise ValueError("projector window exceeds SOC vector")
                for i in (first - 1, last):
                    if i in edge and i + 1 in edge and edge[i + 1] - edge[i] <= 1e-10:
                        raise ValueError(
                            "displaced reference cannot differentiate a split degenerate window"
                        )
            records.append(
                (kx, ky, kz, weight, size, bands, hashlib.sha256(labels).digest())
            )
    if not records:
        raise ValueError(f"no k points in {path}")
    return records


class NativeSOCForceMap:
    def __init__(self, reference, case, root, wien, workspace, native_mpi=()):
        self.reference, self.work = Path(reference).resolve(), Path(workspace).resolve()
        self.root, self.wien = Path(root).resolve(), Path(wien).resolve()
        self.case, self.mpi = case, list(native_mpi)
        self.work.mkdir(parents=True, exist_ok=True)
        rank_file = self.reference / "spinor_ranks.json"
        self.parent_ranks = (
            json.loads(rank_file.read_text()) if rank_file.is_file() else None
        )
        self.legacy_matching = parent_lo_filter(self.reference / (case + ".scfso"))
        self.central_struct = self.reference / (case + ".struct")
        self.structure = Struct(inAngs=False)
        self.structure.ReadStruct(str(self.central_struct))
        self.lattice = NativeLattice.from_structure(self.structure)
        self.density = Density(self.reference / (case + ".clmsum"), self.central_struct)
        self.n = self.density.vector()
        if any(
            (self.reference / (case + ext)).exists()
            and (self.reference / (case + ext)).stat().st_size
            for ext in (".clmup", ".clmdn", ".indmfldn", ".incup", ".incdn")
        ):
            raise ValueError(
                "native SOC force reference requires one nonmagnetic KS potential"
            )
        self.window = tuple(
            map(
                int,
                (self.reference / (case + ".indmfl"))
                .read_text()
                .splitlines()[0]
                .split()[:2],
            )
        )
        self.signature = None
        self.steps = 0
        self.audit = self.work / "native_passes.jsonl"

    def evaluate(self, pattern, t, q):
        # Remove only files created in the fresh private scratch directory on
        # successful extraction. A failed native pass is kept for diagnosis.
        folder = Path(tempfile.mkdtemp(prefix=f"pass_{self.steps:04}_", dir=self.work))
        self.steps += 1
        case = self.case
        suffixes = (
            "in0",
            "in0abp",
            "in1",
            "in1c",
            "in2",
            "in2c",
            "inso",
            "inc",
            "inm",
            "klist",
            "kgen",
            "rotlm",
            "indmfl",
            "indmfi",
            "indmf2",
        )
        for name in [
            *(case + "." + x for x in suffixes),
            "EF.dat",
            "projectorw.dat",
            *(p.name for p in self.reference.glob("sig.inp*")),
        ]:
            source = self.reference / name
            if source.is_file():
                shutil.copy2(source, folder / name)
        for program in ("lapw1", "lapwso", "dmft2"):
            f = program + ".def"
            (folder / f).write_text(local_definition((self.reference / f).read_text()))
        definitions = entries((folder / "lapw1.def").read_text())
        dmft = entries((folder / "dmft2.def").read_text())
        if (Path(dmft[9]).name, Path(dmft[10]).name) != (
            case + ".vectorso",
            case + ".vectorsodn",
        ):
            raise ValueError("expected both native SOC vector components")
        for name in (case + ".in0", dmft[5]):
            p = folder / name
            lines = p.read_text().splitlines()
            lines[0] = "TOT  " + " ".join(lines[0].split()[1:])
            p.write_text("\n".join(lines) + "\n")
        struct = folder / (case + ".struct")
        pos = self.structure.pos + np.linalg.solve(self.lattice.cell, (pattern * t).T).T
        displaced_struct(self.central_struct, struct, pos)
        self.density.from_vector(self.n + t * q).write(folder / (case + ".clmsum"))
        env = {
            k: v
            for k, v in os.environ.items()
            if not k.startswith(("EDMFT_", "WIEN_", "DFTKS_"))
        }
        env.update(
            EDMFT_SOC_LEGACY_LO_FILTER=str(int(self.legacy_matching)),
            WIENROOT=str(self.wien),
            OMP_NUM_THREADS="1",
            OPENBLAS_NUM_THREADS="1",
            VECLIB_MAXIMUM_THREADS="1",
        )

        def run(label, cmd):
            with (folder / (label + ".log")).open("w") as f:
                result = subprocess.run(
                    cmd, cwd=folder, env=env, stdout=f, stderr=subprocess.STDOUT
                )
            log = (folder / (label + ".log")).read_text(errors="replace")
            if result.returncode or any(
                x in log
                for x in (
                    "Error termination",
                    "Fortran runtime error",
                    "ERROR STOP",
                    "DMFT - Error",
                )
            ):
                raise RuntimeError(
                    f"native SOC response failed in {folder}/{label}.log"
                )

        for program in ("lapw0", "lcore"):
            run(program + "_def", [self.wien / "x", program, "-f", case, "-d"])
        run("lapw0", [self.wien / "lapw0", "lapw0.def"])
        run(
            "lapw1",
            self.mpi
            + [
                self.root / ("lapw1c" if definitions[5].endswith(".in1c") else "lapw1"),
                "lapw1.def",
            ],
        )
        # A buffered energy cutoff can select a different number of bands on
        # opposite sides. Follow the parent's fixed rank separately at every k.
        if self.parent_ranks is not None:
            enforce_ranks(
                folder / definitions[10],
                folder / definitions[11],
                self.structure.nat,
                self.parent_ranks["scalar"],
                check_cut=True,
            )
        run("lapwso", self.mpi + [self.root / "lapwso", "lapwso.def"])
        if self.parent_ranks is not None:
            so = entries((folder / "lapwso.def").read_text())
            for vec, energy, norm in ((41, 51, 45), (42, 52, 46)):
                enforce_ranks(
                    folder / so[vec],
                    folder / so[energy],
                    self.structure.nat,
                    self.parent_ranks["spinor"],
                    folder / so[norm],
                    check_cut=True,
                )
        signatures = (
            _signature(folder / definitions[10], self.structure.nat),
            _signature(folder / dmft[9], self.structure.nat, self.window),
        )
        if self.signature is not None and signatures != self.signature:
            raise ValueError(
                f"G labels, band counts, k points or weights changed in {folder}"
            )
        self.signature = signatures
        run("dmft2", self.mpi + [self.root / "dmft2", "dmft2.def"])
        run("lcore", [self.wien / "lcore", "lcore.def"])
        raw = self.density.zeros()
        for ext in (".clmval", ".clmcor"):
            raw.add(Density(folder / (case + ext), struct))
        s = Struct(inAngs=False)
        s.ReadStruct(str(struct))
        inm = (folder / (case + ".inm")).read_text().splitlines()[0].split()
        if float(inm[1]) != 0 or inm[2].upper() not in ("YES", "NO"):
            raise ValueError("invalid or charged INM normalization")
        norm = NativeChargeNormalization(
            raw,
            s,
            np.dot(s.mult, s.Znuc),
            inversion=native_mixer_inversion(s),
            enabled=inm[2].upper() == "YES",
        )
        terms = np.array(
            [
                _scalar(folder / (case + ".scf0"), ":DEN "),
                _core_energy(folder, case, s),
                _scalar(folder / (case + ".scf2"), ":VSUM"),
            ]
        )
        value = norm.central.copy()
        with self.audit.open("a") as audit:
            audit.write(
                json.dumps(
                    dict(
                        pass_index=self.steps,
                        signed_displacement_Bohr=t,
                        density_norm=float(np.linalg.norm(q)),
                        density_integral_Ry=float(terms[0]),
                        core_Ry=float(terms[1]),
                        trace_log_plus_mu_N_Ry=float(terms[2]),
                        normalization=float(norm.alpha),
                        nk=len(signatures[0]),
                    )
                )
                + "\n"
            )
        shutil.rmtree(folder)
        return value, terms

    def solve(self, atom, direction, pattern, h, tol, maxiter, progress, initial=None):
        zero = np.zeros_like(self.n)

        def pair(vector):
            plus, ep = self.evaluate(pattern, h, vector)
            minus, em = self.evaluate(pattern, -h, vector)
            return (plus - minus) / (2 * h), (ep - em) / (2 * h)

        source, explicit = pair(zero)
        energy = explicit

        def action(vector):
            nonlocal energy
            if not np.any(vector):
                energy = explicit
                return zero
            out, energy = pair(vector)
            return out - source

        factor = -1000.0 / self.structure.mult[atom - 1]

        def monitor(x, residual, row):
            row["force_mRy_Bohr"] = factor * math.fsum(energy)
            progress(
                f"atom {atom}, h={h:g} Bohr, iteration {row['iteration']}, "
                f"residual {row['residual_norm']:.7g} (target {tol * np.linalg.norm(source):.7g}), "
                f"force {row['force_mRy_Bohr']:.9f} mRy/Bohr"
            )

        vector, report = solve_response(
            source,
            action,
            initial=zero if initial is None else initial,
            absolute_tolerance=0.0,
            relative_tolerance=tol,
            max_iterations=maxiter,
            callback=monitor,
        )
        return (
            dict(
                atom=atom,
                direction_cartesian=direction.tolist(),
                step_Bohr=h,
                force_component_mRy_Bohr=factor * math.fsum(energy),
                explicit_force_mRy_per_Bohr=factor * math.fsum(explicit),
                screened_force_mRy_per_Bohr=factor * math.fsum(energy - explicit),
                response=report,
            ),
            vector,
        )


def run_reference(
    reference,
    case,
    root,
    wien,
    work,
    *,
    atoms,
    steps,
    tolerance,
    max_iterations,
    mpi,
    progress,
):
    """Two step sizes check that a force can be published at the stated accuracy."""
    model = NativeSOCForceMap(reference, case, root, wien, work, mpi)
    s = model.structure
    selected = atoms or list(range(1, s.nat + 1))
    check_force_frame(s, selected)
    modes = internal_modes(s, selected)
    force = np.zeros((s.nat, 3))
    records = []
    if (
        len(steps) != 2
        or min(steps) <= 0
        or not max(steps) > min(steps)
        or not np.all(np.isfinite(steps))
    ):
        raise ValueError(
            "screened SOC reference needs two distinct positive displacement steps"
        )
    if (
        not 0 < tolerance < 1
        or not isinstance(max_iterations, int)
        or max_iterations < 1
    ):
        raise ValueError("invalid charge-response convergence criterion")
    if modes:
        zeros = np.zeros_like(model.n)
        nout, _ = model.evaluate(modes[0][2], 0.0, zeros)
        # Report the residual in precisely the same native coordinates as the
        # linear solve. A derivative about an unconverged central density is
        # not the force of the central self-consistent state.
        relative = np.linalg.norm(nout - model.n) / max(1.0, np.linalg.norm(model.n))
        progress(f"Central SOC charge-map relative residual {relative:.6g}")
        if relative > 1e-5:
            raise ValueError(
                f"SOC reference central charge is unconverged (relative residual {relative:.6g})"
            )
    else:
        relative = None
    for atom, direction, pattern in modes:
        pair = []
        # All step sizes approximate the same first derivative. The previous
        # solution is an initial guess; reevaluate its residual at the new h.
        initial = None
        for h in sorted(steps, reverse=True):
            row, initial = model.solve(
                atom,
                direction,
                pattern,
                h,
                tolerance,
                max_iterations,
                progress,
                initial,
            )
            if not row["response"]["converged"]:
                raise RuntimeError(
                    f"SOC forward response did not converge for atom {atom}, h={h}"
                )
            pair.append(row)
        big, small = pair
        factor = (big["step_Bohr"] / small["step_Bohr"]) ** 2 - 1.0
        component = dict(atom=atom, direction_cartesian=direction.tolist(), steps=pair)
        for key in (
            "force_component_mRy_Bohr",
            "screened_force_mRy_per_Bohr",
            "explicit_force_mRy_per_Bohr",
        ):
            component[key] = small[key] + (small[key] - big[key]) / factor
        err = abs(
            component["force_component_mRy_Bohr"] - small["force_component_mRy_Bohr"]
        )
        component["extrapolation_shift_mRy_per_Bohr"] = err
        progress(
            f"Atom {atom}, direction {direction}, screened SOC force {component['force_component_mRy_Bohr']:.8f} mRy/Bohr; extrapolation shift {err:.6g}"
        )
        force[atom - 1] += component["force_component_mRy_Bohr"] * direction
        records.append(component)
    return dict(
        observable="screened frozen-Sigma atomic force",
        case=case,
        scope="native nonmagnetic spinor frozen-Sigma charge response",
        source_method="native_soc_forward",
        structure=str(model.central_struct),
        central_relative_charge_residual=relative,
        converged=True,
        selected_atoms=selected,
        force_mRy_Bohr=force.tolist(),
        components=records,
    )
