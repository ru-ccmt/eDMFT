"""Optional native regression for eigenvector phases in LAPWSO local-orbital matching.

Set EDMFT_TEST_SOC_CASE to an unsplit, complex, nonmagnetic LAPWSO case,
EDMFT_TEST_LAPWSO to the built executable, and optionally EDMFT_TEST_SOC_K
to the one-based k index. Only this k point is copied to the temporary case.
"""

import os
from pathlib import Path
import shutil
import struct
import subprocess
from tempfile import TemporaryDirectory
import unittest

import numpy as np

from dmft2_step import entries, local_definition
from lapw_dft_vector import _header, _record
from lapw_soc_rank import _put, parent_lo_filter
from wstruct import Struct


def copy_k_point(vector, energy, destination, names, nat, selected):
    """Keep native radial and k headers and save the two versions of the bands."""
    with vector.open("rb") as v, energy.open() as e:
        radial = [_record(v) for _ in range(2 * nat)]
        energy_header = [e.readline() for _ in range(2 * nat)]
        for index in range(1, selected + 1):
            header = _record(v)
            if header is None:
                raise ValueError(f"k point {selected} is absent from {vector}")
            basis, count = _header(header)[4:6]
            label = _record(v)
            energy_line = e.readline()
            if int(energy_line[73:79]) != count:
                raise ValueError("scalar energy and vector ranks differ")
            energies = [e.readline() for _ in range(count)]
            if index != selected:
                for _ in range(2 * count):
                    _record(v, skip=True)
                continue
            (destination / names[54]).write_text("".join(energy_header + [energy_line] + energies))
            rng = np.random.default_rng(2026)
            with (destination / "original.vector").open("wb") as orig, (
                destination / "rotated.vector"
            ).open("wb") as rotated:
                for stream in (orig, rotated):
                    for record in [*radial, header, label]:
                        _put(stream, record)
                for _ in range(count):
                    band, data = _record(v), _record(v)
                    if len(data) != 16 * basis:
                        raise ValueError("the phase regression needs complex scalar vectors")
                    a = np.frombuffer(data, "<c16")
                    phase = np.exp(1j * rng.uniform(-np.pi, np.pi))
                    _put(orig, band)
                    _put(orig, data)
                    _put(rotated, band)
                    _put(rotated, (a * phase).tobytes())


@unittest.skipUnless(
    os.environ.get("EDMFT_TEST_SOC_CASE") and os.environ.get("EDMFT_TEST_LAPWSO"),
    "set EDMFT_TEST_SOC_CASE and EDMFT_TEST_LAPWSO for a native SOC case",
)
class NativePhaseRegression(unittest.TestCase):
    def test_scalar_phases_leave_native_spinor_spectrum_invariant(self):
        case = Path(os.environ["EDMFT_TEST_SOC_CASE"]).expanduser().resolve()
        binary = Path(os.environ["EDMFT_TEST_LAPWSO"]).expanduser().resolve()
        selected = int(os.environ.get("EDMFT_TEST_SOC_K", "1"))
        self.assertGreaterEqual(selected, 1)
        definition = (case / "lapwso.def").read_text()
        original = entries(definition)
        isolated = local_definition(definition)
        names = entries(isolated)
        with TemporaryDirectory() as temp:
            d = Path(temp)
            (d / "lapwso.def").write_text(isolated)
            for unit in (4, 5, 18, 20, 22):
                shutil.copyfile(case / original[unit], d / names[unit])
            structure = Struct()
            structure.ReadStruct(str(d / names[20]))
            copy_k_point(case / original[9], case / original[54], d, names, structure.nat, selected)
            # Test the production default independently of any diagnostic
            # replay controls or exports in the calling process's environment.
            env = {k: v for k, v in os.environ.items() if not k.startswith("EDMFT_SOC_")}
            spectra = []
            for label in ("original", "rotated"):
                shutil.copyfile(d / (label + ".vector"), d / names[9])
                run = subprocess.run(
                    [binary, "lapwso.def"], cwd=d, env=env, text=True,
                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False,
                )
                self.assertEqual(run.returncode, 0, run.stdout)
                self.assertIn("LAPWSO END", run.stdout)
                self.assertFalse(parent_lo_filter(d / names[8]))
                with (d / names[42]).open("rb") as spinor:
                    for _ in range(2 * structure.nat):
                        _record(spinor, skip=True)
                    count = _header(_record(spinor))[5]
                    self.assertGreater(count, 0)
                    _record(spinor, skip=True)
                    values = []
                    for _ in range(count):
                        values.append(struct.unpack("<id", _record(spinor))[1])
                        _record(spinor, skip=True)
                spectra.append(values)
            np.testing.assert_allclose(spectra[1], spectra[0], atol=1e-11, rtol=0)
