#!/usr/bin/env python3
"""Isolate two explicit central SOC force contributions on a symmetry orbit.

Run LAPWSO with EDMFT_SOC_RESPONSE_GEOMETRY and DMFT2 with EDMFT_ONEPARTICLE_DUMP
on the same vectors. Supply the central LAPW0/LCORE potential archive. The
output reports the derivative of the SOC matrix at fixed scalar vectors,
split by radial integrals, radial matching, and rigid translation. Optional
dftKS radial-energy searches follow the same native input Fermi reference.
No term here is, in isolation, a correction to the total mixer force.
"""
import argparse
from contextlib import ExitStack
from pathlib import Path
import json

import numpy as np

from force_response_modes import internal_modes
from lapw_potential_core_state import load_potential_core_state
from lapw_soc_response import SOCGeometry, SOCBandHamiltonian
from lapw_soc_subspace import spinor_coefficients
from lapw_strain_field import diagonal_strain_potential
from lapw_valence_state import _radial_fields
from soc_scalar_window import scalar_vectors
from soc_translation_force import read
from response_parallel import SERIAL, ResponseParallel


def spinor_moment(path, ik, labels, vector_k):
    with Path(path).open("rb") as f:
        if f.read(8) != b"EDMFTRH1":
            raise ValueError(f"invalid moment {path}")
        order, v, k, nk, nb, nlo, n, spin, first, dmfirst, last, _, _ = read(
            f, "<i4", 13
        )
        if (order, v, k, nb, spin) != (
            0x01020304,
            1,
            ik,
            labels.shape[1],
            2,
        ) or n != last - first + 1:
            raise ValueError(f"incompatible spinor moment at k={ik}")
        frac = read(f, "<f8", 3)
        g = read(f, "<i4", 3 * nb).reshape((3, nb), order="F")
        if not np.array_equal(labels, g) or not np.allclose(
            vector_k, frac, atol=1e-10, rtol=0
        ):
            raise ValueError(f"inconsistent k point or LAPW row order at k={ik}")
        energies = read(f, "<f8", n)
        coeff = read(f, "<c16", nb * n * 2).reshape((nb, n, 2), order="F")
        rho = read(f, "<c16", n * n).reshape((n, n), order="F")
        if not np.allclose(rho, rho.conj().T, atol=1e-10, rtol=0):
            raise ValueError(f"non-Hermitian density at k={ik}")
        if path.stat().st_size - f.tell() != n * n * 16:
            raise ValueError(f"incorrect energy-moment length at k={ik}")
        return coeff, energies, rho


def components(
    directory,
    case,
    geometry,
    moments,
    potential_core,
    atoms,
    *,
    library=None,
    energy_policy="fixed",
    parallel=SERIAL,
):
    directory = Path(directory)
    pc = load_potential_core_state(potential_core)
    s = pc.structure
    from lapw_potential_core_state import fingerprint_files

    if (
        fingerprint_files({"struct": directory / (case + ".struct")})["struct"]
        != pc.metadata["source_fingerprints"]["struct"]
    ):
        raise ValueError("potential/core archive and the scalar structure differ")
    modes = internal_modes(s, atoms)
    central_v = [a[0, 0] for a in _radial_fields(directory / (case + ".vsp"), s)]
    if energy_policy not in ("fixed", "dftKS"):
        raise ValueError("unrecognized radial-energy policy")
    searches = {}
    if energy_policy == "dftKS":
        from lapw_radial_energy_policy import searched_energy_derivatives

        searches, _ = searched_energy_derivatives(
            directory,
            case,
            complex_vectors=(directory / (case + ".in1c")).exists(),
            library=library,
            scalar_solver="dftKS",
        )
    files = sorted(Path(moments).parent.glob(Path(moments).name + ".k????????.bin"))
    ownership = set(parallel.partition([p.stat().st_size for p in files]))
    with ExitStack() as stack:
        geo = SOCGeometry(geometry, s, central_v, stack, library)
        directions = []
        for atom, direction, pattern in modes:
            v, _ = diagonal_strain_potential(
                pc, np.zeros((3, 3)), displacements=pattern
            )
            dv = {
                a: v[sl]
                for (a, l, m), sl in zip(pc.hartree.keys, pc.hartree.slices)
                if (l, m) == (0, 0)
            }
            dots = {}
            for key, ch in geo.channels.items():
                de = (
                    searches[key] @ dv[key[0]]
                    if key in searches
                    else np.zeros(1 + ch.ctx.nlo)
                )
                dots[key] = ch.action(dv[key[0]], de[0], de[1:])
            directions.append(dots)
        result = np.zeros((len(modes), 3))
        ritz = 0.0
        count = 0
        for ik, (kf, labels, eps, a) in enumerate(
            scalar_vectors(directory / (case + ".vector"), s.nat), 1
        ):
            count = ik
            if ik - 1 not in ownership:
                continue
            expected = Path(str(moments) + f".k{ik:08d}.bin")
            if files[ik - 1] != expected:
                raise ValueError("nonconsecutive spinor moments")
            native, energy, rho = spinor_moment(expected, ik, labels, kf)
            c = spinor_coefficients(a, native)
            op = SOCBandHamiltonian(geo, kf, labels, a)
            h = op.action(central=True)
            h[np.diag_indices_from(h)] += np.tile(eps, 2)
            ritz = max(ritz, float(np.max(abs(h @ c - c * energy))))
            radial, _, gu = op.pullback(c @ rho @ c.conj().T)
            for j, ((atom, direction, pattern), dots) in enumerate(
                zip(modes, directions)
            ):
                # DMFT2 stores 2*w*rho. Sum the two spinor components once
                # and divide by equivalent-site multiplicity to obtain one force.
                scale = -500 / s.mult[atom - 1]
                rr = sum(
                    np.vdot(bar[0], dots[key][0]).real for key, bar in radial.items()
                )
                bb = sum(
                    np.vdot(x, y).real
                    for key, bar in radial.items()
                    for x, y in zip(bar[1:], dots[key][1:])
                )
                result[j] += scale * np.array([rr, bb, np.vdot(gu, pattern).real])
        if count != len(files):
            raise ValueError("scalar mesh and spinor moment count differ")
        ritz = parallel.reduce(ritz, "max")
        if ritz > 1e-7:
            raise ValueError(f"central spinor Ritz residual is {ritz:g} Ry")
        result = parallel.sum(result)
        return dict(
            case=case,
            k_points=count,
            maximum_ritz_residual_Ry=ritz,
            radial_energy_policy=energy_policy,
            units="mRy/Bohr",
            modes=[
                dict(
                    atom=int(a),
                    direction=u.tolist(),
                    radial_integrals=r,
                    radial_matching=b,
                    translation=t,
                    radial_sum=r + b,
                )
                for (a, u, _), (r, b, t) in zip(modes, result)
            ],
        )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for key in ("directory", "geometry", "moments", "potential-core"):
        p.add_argument("--" + key, type=Path, required=True)
    p.add_argument("--case", required=True)
    p.add_argument("--atoms", type=int, nargs="+", required=True)
    p.add_argument("--library", type=Path)
    p.add_argument("--energy-policy", choices=("fixed", "dftKS"), default="fixed")
    p.add_argument("--mpi", action="store_true")
    args = vars(p.parse_args())
    par = ResponseParallel.world() if args.pop("mpi") else SERIAL
    result = components(**args, parallel=par)
    if par.root:
        print(json.dumps(result, indent=2))
