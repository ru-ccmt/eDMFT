"""Screening adapter, explicit source/observable, and stale-state controls."""

from types import SimpleNamespace as NS
from unittest.mock import patch
import unittest
import numpy as np
from lapw_screened_state import NativeScreenedState
from stress_response_adjoint import PairedOperator, ResponseConvergenceError
from wien_response_density import Density


def op(a):
    return PairedOperator(
        a.shape[1],
        a.shape[0],
        lambda x: a @ x,
        lambda y: a.T @ y,
        "dense independent control",
    )


def density(values):
    d = Density.__new__(Density)
    d.rad = {(0, 0, 0): np.asarray(values[:2], float)}
    d.pw = {(0, 0, 0): complex(*values[2:])}
    d.mult = [1]
    d.grids = [np.array([0.0, 1.0])]
    return d


class ScreenedStateTests(unittest.TestCase):
    def test_mixer_normalization_policy_cannot_be_overridden_silently(self):
        for recorded, requested in [(False, True), (True, False)]:
            v = NS(mixer_normalization_enabled=recorded)
            with self.assertRaisesRegex(ValueError, "differs from case.inm"):
                NativeScreenedState(v, None, None, target_charge=1,
                                    fixed_number=False,
                                    normalization_enabled=requested)

    def setup_state(self, fixed=True):
        rng = np.random.default_rng(73)
        a = {k: rng.normal(size=(4, 4)) * 0.08 for k in ["K", "B", "E", "R", "D", "C"]}
        a["select"] = np.diag([1.0, 2.0, 0.5, 1.5])
        a["N"] = np.eye(4)
        a["N"][2, :] = [-0.2, -0.3, 0.0, 0.1]
        base = density([0.8, 0.9, 2.0, 0.1])
        core = density([0.2, 0.3, 0.0, 0.0])
        s = NS(nat=1, lattyp="P", mode="NREL")
        for key in [
            "a",
            "b",
            "c",
            "alpha",
            "beta",
            "gamma",
            "mult",
            "jrj",
            "r0",
            "rmt",
            "pos",
            "timat",
            "tau",
            "rotloc",
            "iatnr",
            "Znuc",
        ]:
            setattr(s, key, np.array([1.0]))
        coords = NS(
            size=4,
            nonspherical_shapes={},
            spherical_size=2,
            nstar=1,
            complex_interstitial=True,
        )
        norm = NS(
            operator=op(a["N"]),
            mt=np.array([1.0, 2.0, 0.0, 0.0]),
            ir=np.array([0.0, 0.0, 1.0, 0.0]),
            alpha=0.9,
            central=base.vector() + core.vector(),
        )
        v = NS(
            structure=s,
            reference=base,
            fingerprints={"case.struct": "s", "case.vsp": "v"},
            bands=NS(
                spherical=NS(slices={0: slice(0, 2)}),
                ns={},
                ir=[NS(representatives=np.array([[0, 0, 0]]))],
                input=coords,
                output=NS(size=4),
                operator=op(a["B"]),
                spherical_selector=op(np.eye(4)),
            ),
            spatial=NS(value=base.vector, operator=op(a["R"]), coordinates=NS(size=4)),
            direct=NS(operator=op(a["D"])),
            complement=object(),
            records=[0, 1],
            model=NS(eps=np.zeros((2, 3)), meta=NS(mu=0.1)),
            matrix_quadrature="wien2k24",
        )
        pc = NS(
            structure=s,
            input_density=base,
            metadata={"source_fingerprints": {"struct": "s", "vsp": "v"}},
            hartree=NS(
                native_potential_output=lambda *args, **kw: (op(a["select"]), coords)
            ),
            core=NS(
                atoms={0: NS(density=core.rad[0, 0, 0])},
                slices={0: slice(0, 2)},
                operator=op(a["C"]),
            ),
            normalization=lambda *args, **kw: norm,
            potential_operator=op(a["K"]),
        )
        electronic = patch(
            "lapw_screened_state.complement_electronic_operator",
            return_value=(op(a["E"]), NS(size=4), NS(size=4)),
        )
        with electronic as mocked:
            state = NativeScreenedState(
                v, pc, core, target_charge=3.0, fixed_number=fixed
            )
        mocked.assert_called_once_with(v.complement, fixed_number=fixed)
        return state, v, pc, core, a, electronic

    def test_native_chain_against_independent_dense_screening(self):
        state, _, _, _, a, _ = self.setup_state()
        j = (
            a["N"]
            @ ((a["R"] @ a["E"] @ a["B"] + a["D"]) @ a["select"] + a["C"])
            @ a["K"]
        )
        b = np.array([0.1, -0.3, 0.2, 0.4])
        actual = state.solve_source(b, rtol=1e-12, atol=1e-14)
        expected = np.linalg.solve(np.eye(4) - j, b)
        np.testing.assert_allclose(actual.density, expected, rtol=1e-12, atol=1e-14)
        np.testing.assert_allclose(
            actual.induced_potential, a["K"] @ expected, rtol=1e-12
        )
        self.assertLess(state.jacobian.check_transpose(), 1e-14)
        self.assertFalse(actual.report["stress_applied"])
        self.assertEqual(actual.report["kpoints"], 2)

    def test_one_observable_adjoint_for_two_sources(self):
        state, _, _, _, a, _ = self.setup_state(False)
        c = np.array([1.0, 0.4, -0.2, 0.1])
        g = np.array([0.2, -0.1, 0.3, 0.7])
        gradient = state.observable_gradient(c, g)
        sources = np.array([[1.0, 2.0], [0.3, -0.4], [-0.1, 0.2], [0.5, 0.8]])
        explicit = np.array([0.17, -0.09])
        values, lam, report = state.screened.contractions(
            gradient, sources, explicit, rtol=1e-12, atol=1e-14
        )
        direct = []
        for b, e in zip(sources.T, explicit):
            q = state.solve_source(b, rtol=1e-12, atol=1e-14)
            direct.append(e + c @ q.density + g @ q.induced_potential)
        np.testing.assert_allclose(values, direct, rtol=1e-12, atol=1e-14)
        self.assertTrue(report["converged"])

    def test_layout_and_stale_snapshot_are_rejected(self):
        _, v, pc, core, a, e = self.setup_state()
        for key in ["struct", "vsp"]:
            old = pc.metadata["source_fingerprints"][key]
            pc.metadata["source_fingerprints"][key] = "stale"
            with self.assertRaisesRegex(ValueError, "fingerprints"):
                NativeScreenedState(v, pc, core, target_charge=3, fixed_number=True)
            pc.metadata["source_fingerprints"][key] = old
        pc.input_density = pc.input_density.copy()
        pc.input_density.grids[0][1] += 0.01
        with self.assertRaisesRegex(ValueError, "layouts"):
            NativeScreenedState(v, pc, core, target_charge=3, fixed_number=True)

    def test_wrong_central_core_and_ensemble_are_rejected(self):
        _, v, pc, core, a, e = self.setup_state()
        with self.assertRaisesRegex(ValueError, "fixed_number"):
            NativeScreenedState(v, pc, core, target_charge=3, fixed_number=None)
        bad = core.copy()
        bad.rad[0, 0, 0][0] += 0.01
        with self.assertRaisesRegex(ValueError, "reconstruction"):
            NativeScreenedState(v, pc, bad, target_charge=3, fixed_number=True)
        bad = core.copy()
        bad.grids[0][1] += 0.01
        with self.assertRaisesRegex(ValueError, "core density layout"):
            NativeScreenedState(v, pc, bad, target_charge=3, fixed_number=True)

    def test_failed_solve_returns_no_density_result(self):
        state, _, _, _, a, _ = self.setup_state()
        from stress_response_adjoint import ScreenedResponse

        state.screened = ScreenedResponse(op(np.eye(4)))
        with self.assertRaises(ResponseConvergenceError):
            state.solve_source(np.ones(4), max_cycles=1)
        with self.assertRaises(ValueError):
            state.observable_gradient(np.zeros(3), np.zeros(4))


if __name__ == "__main__":
    unittest.main()
