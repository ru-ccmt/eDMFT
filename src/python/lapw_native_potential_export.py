"""Prepare a portable potential/core archive from one native central export.

The isolated LAPW0/LCORE response workers write xc_geometry_*.bin,
xc_native.bin and core_native.bin. These stream formats are native-endian
4-byte integers/8-byte reals, and must come from the same central calculation.
No material names, temporary paths, or external Python scripts are required.
"""

from pathlib import Path
import re
from lapw_lattice import NativeLattice, native_structure
import numpy as np
from wstruct import Struct
from wien_response_density import Density, FLOAT
from lapw_charge_response import native_mixer_inversion
from lapw_potential_core_state import save_potential_core_state


def _read(f, dtype, count):
    values = np.fromfile(f, dtype, count)
    if len(values) != count or not np.all(np.isfinite(values)):
        raise ValueError(f"truncated or nonfinite native export: {f.name}")
    return values


def prepare_native_potential_core(directory, case, output):
    """Validate complete native records before writing the numeric/JSON archive."""
    folder = Path(directory)
    s = Struct()
    s.ReadStruct(str(folder / (case + ".struct")))
    s = native_structure(s)
    density = Density(folder / (case + ".clmsum"), folder / (case + ".struct"))
    in0 = (folder / (case + ".in0")).read_text().splitlines()
    token = in0[0].split()[1].upper()
    functional = {"XC_LDA": "LDA", "XC_PBE": "PBE"}.get(token)
    if functional is None:
        raise NotImplementedError("explicit XC_LDA or XC_PBE required")
    ncon = re.search(r"\bNCON\s+(\d+)", in0[2].upper())
    ncon = int(ncon[1]) if ncon else 9
    with (folder / "xc_native.bin").open("rb") as f:
        nat, nstar, *shape = _read(f, "i4", 5)
        volume = _read(f, "f8", 1)[0]
        if nat != s.nat or min(shape) < 1 or nstar < 1:
            raise ValueError("native XC header differs")
        if not np.isclose(volume, NativeLattice.from_structure(s).volume, rtol=1e-10):
            raise ValueError("native XC volume differs")
        reps = _read(f, "i4", 3 * nstar).reshape(nstar, 3)
        if not np.array_equal(reps, list(density.pw)):
            raise ValueError("native Fourier labels differ")
        for a in range(s.nat):
            nr, nl = _read(f, "i4", 2)
            lm = _read(f, "i4", 2 * nl).reshape(nl, 2)
            if nr != s.jrj[a] or not np.array_equal(
                lm, [k[1:] for k in density.rad if k[0] == a]
            ):
                raise ValueError("native XC radial layout differs")
            _read(f, "f8", nr * nl)
        _read(f, "c16", int(np.prod(shape)))
        if f.read():
            raise ValueError("trailing native XC records")
    angular = []
    for a in range(s.nat):
        with (folder / f"xc_geometry_{a+1:04d}.bin").open("rb") as f:
            nr, nl, ng = _read(f, "i4", 3)
            if nr != s.jrj[a] or nl < 1 or ng < 1:
                raise ValueError("invalid angular export header")
            lm = _read(f, "i4", 2 * nl).reshape(nl, 2)
            r = _read(f, "f8", nr)
            if not np.array_equal(lm, [k[1:] for k in density.rad if k[0] == a]):
                raise ValueError("angular harmonic order differs")
            if not np.allclose(
                r, np.geomspace(s.r0[a], s.rmt[a], s.jrj[a]), rtol=1e-11, atol=1e-14
            ):
                raise ValueError("angular radial grid differs")
            entry = dict(
                angular=_read(f, "f8", nl * ng).reshape(ng, nl),
                weights=_read(f, "f8", ng),
            )
            if functional == "PBE":
                entry.update(
                    sine=_read(f, "f8", ng),
                    cosine=_read(f, "f8", ng),
                    angular_derivatives=_read(f, "f8", 5 * ng * nl).reshape(5, ng, nl),
                )
            if f.read():
                raise ValueError("trailing angular export records")
            angular.append(entry)
    lines = (folder / (case + ".inc")).read_text().splitlines()
    i = 0
    orbitals = []
    for a in range(s.nat):
        words = lines[i].split()
        count = int(words[0])
        shift = float(words[1])
        i += 1
        occ = []
        for j in range(count):
            words = lines[i].replace(",", " ").split()
            i += 1
            occ.append((int(words[0]), int(words[1]), float(words[2])))
        orbitals.append((shift, occ))
    text = (folder / (case + ".vsp")).read_text()
    vsp = []
    for atom in re.finditer(r"ATOMNUMBER\s*=\s*(\d+)", text):
        a = int(atom[1]) - 1
        end = text.find("ATOMNUMBER", atom.end())
        end = len(text) if end < 0 else end
        channel = re.search(r"VLM\(R\) FOR L\s*0\s+M=\s*0", text[atom.end() : end])
        if channel is None or a != len(vsp):
            raise ValueError("missing spherical potential atom")
        values = [
            float(x.group().replace("D", "E"))
            for x in FLOAT.finditer(text[atom.end() + channel.end() : end])
        ][: s.jrj[a]]
        if len(values) != s.jrj[a]:
            raise ValueError("incomplete spherical potential")
        vsp.append(np.array(values))
    if len(vsp) != s.nat:
        raise ValueError("missing spherical potentials")
    cores = []
    with (folder / "core_native.bin").open("rb") as f:
        for a, (shift, occ) in enumerate(orbitals):
            states = []
            for n, k, occupation in occ:
                nn, l, kk, npnt, end, match, nuc = _read(f, "i4", 7)
                energy, z, c, gamma = _read(f, "f8", 4)
                if (
                    nn != n
                    or kk != k
                    or nuc != 0
                    or z != s.Znuc[a]
                    or (c < 1e10) != (s.mode.strip() == "RELA")
                ):
                    raise ValueError("core state/INC policy mismatch")
                if not (1 < match <= end <= npnt) or npnt < s.jrj[a]:
                    raise ValueError("invalid core match grid")
                _read(f, "f8", 4 * npnt)
                states.append(
                    dict(
                        energy_Ha=float(energy),
                        kappa=int(kk),
                        occupation=occupation,
                        match=int(match - 1),
                        end=int(end - 1),
                    )
                )
            cores.append(dict(vsp=vsp[a], shift_Ry=shift, states=states))
        if f.read():
            raise ValueError("trailing core records: remove stale append-only export")
    files = {
        name: folder / (case + "." + name)
        for name in ["struct", "clmsum", "in0", "vsp", "inc"]
    }
    save_potential_core_state(
        output,
        density=density,
        structure=s,
        fft_shape=np.array(shape),
        inversion=native_mixer_inversion(s),
        functional=functional,
        angular_data=angular,
        core_data=cores,
        source_files=files,
        ncon=ncon,
    )
