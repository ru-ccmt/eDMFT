"""Cell-potential geometry versus freshly constructed nonlinear operators."""

from copy import deepcopy
from pathlib import Path
import tempfile, unittest
import numpy as np
from lapw_strain_field import harmonic_average_direction, diagonal_strain_potential
from lapw_hartree_periodic import NativeStarTable
from lapw_hydrostatic_potential import hydrostatic_explicit_potential
from lapw_potential_core_state import save_potential_core_state
from lapw_xc_response import InterstitialXCResponse
from test_lapw_potential_core_state import inputs


class StrainFieldTests(unittest.TestCase):
    def test_star_harmonics_with_nontrivial_frame(self):
        reps = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 2], [0, 0, 2]])
        rot = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
        stars = NativeStarTable(
            reps,
            np.array([np.linalg.matrix_power(rot, k) for k in range(4)]),
            np.zeros((4, 3)),
            inversion=False,
        )
        b = np.diag([0.8, 0.8, 0.6])
        e = np.diag([0.3, 0.3, -0.4])
        pos = np.array([0.1, 0.2, 0.4])
        frame = np.array([[0.8, 0, 0.6], [0, 1, 0], [-0.6, 0, 0.8]])
        lm = [(l, m) for l in range(5) for m in range(-l, l + 1)]
        actual = harmonic_average_direction(stars, b, pos, frame, lm, e)
        h = 2e-6
        ends = [
            stars.harmonic_average(np.linalg.inv(np.eye(3) + t * e) @ b, pos, frame, lm)
            for t in [-h, h]
        ]
        np.testing.assert_allclose(
            actual, (ends[1] - ends[0]) / (2 * h), atol=5e-10, rtol=2e-7
        )

    def test_diagonal_coulomb_xc_against_reprepared_operators(self):
        args = inputs()
        d = args["density"]
        for g, z in [
            ((1, 0, 0), 0.002 + 0.001j),
            ((-1, 0, 0), 0.002 - 0.001j),
            ((1, 1, 1), 0.0003 + 0.0002j),
            ((-1, -1, -1), 0.0003 - 0.0002j),
        ]:
            d.pw[g] = z
        e = np.diag([0.2, -0.3, 0.7])
        h = 2e-6
        with tempfile.TemporaryDirectory() as tmp:

            def prepare(t, functional):
                aa = deepcopy(args)
                for name, rate in zip(["a", "b", "c"], np.diag(e)):
                    setattr(
                        aa["structure"],
                        name,
                        getattr(args["structure"], name) * (1 + t * rate),
                    )
                pc = save_potential_core_state(Path(tmp) / "state.npz", **aa)
                # MT is fixed here, so either central MT XC convention gives
                # zero explicit geometry. Replace only IR by a full PBE fit.
                if functional == "PBE":
                    pc.xc.ir = InterstitialXCResponse(
                        pc.hartree.reanalysis,
                        np.array(list(d.pw.values())),
                        [pc.structure.a, pc.structure.b, pc.structure.c],
                        functional="PBE",
                    )
                return pc

            for functional in ["LDA", "PBE"]:
                pc = prepare(0, functional)
                actual, _ = diagonal_strain_potential(pc, e, transport_rate=0)
                ends = [
                    hydrostatic_explicit_potential(prepare(t, functional))[1][
                        "central_potential"
                    ]
                    for t in [-h, h]
                ]
                np.testing.assert_allclose(
                    actual, (ends[1] - ends[0]) / (2 * h), atol=3e-8, rtol=3e-7
                )
                hydro, _ = diagonal_strain_potential(pc, np.eye(3))
                original, _ = hydrostatic_explicit_potential(pc)
                np.testing.assert_allclose(hydro, original, atol=3e-13, rtol=5e-13)


if __name__ == "__main__":
    unittest.main()
