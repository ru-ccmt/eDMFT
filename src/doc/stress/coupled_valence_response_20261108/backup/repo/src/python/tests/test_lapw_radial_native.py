"""ABI, units, lifecycle and composed basis-density controls.

Set EDMFT_RADIAL_LIBRARY to a built liblapw_response.so for native tests.
The material validation additionally compares against native dmft2 endpoints;
these small portable tests isolate interface and composition mistakes.
"""

from contextlib import ExitStack
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import os
import unittest
import numpy as np

from lapw_radial_native import NativeRadialContext
from lapw_basis_density_response import BasisDensityResponse
from lapw_augmentation_response import PreparedAugmentationChannel
from lapw_density_response import (
    MuffinTinDensitySum,
    SpatialValenceDensity,
    InterstitialDensity,
    FourierStars,
    herm,
)
from wien_response_density import Density


LIBRARY = Path(
    os.environ.get(
        "EDMFT_RADIAL_LIBRARY",
        str(Path(__file__).resolve().parents[2] / "dmft2/atpar/liblapw_response.so"),
    )
)


def inputs(lapw=True, nlo=1, rel=False, ell=1):
    n = 301
    r0 = 1e-4
    rmt = 2.0
    dx = np.log(rmt / r0) / (n - 1)
    r = r0 * np.exp(dx * np.arange(n))
    v = -4.0 + 0.2 * r * r
    elo = np.array([-0.6, -2.0][:nlo])
    if not lapw and nlo:
        elo[0] = 0.3
    return (
        v,
        dict(
            r0=r0,
            dx=dx,
            rmt=rmt,
            energy=0.3,
            lo_energies=elo,
            ell=ell,
            charge=2.0,
            relativistic=rel,
            lapw=lapw,
            library=LIBRARY,
        ),
        r,
    )


@unittest.skipUnless(LIBRARY.is_file(), "build liblapw_response.so for native tests")
class NativeTests(unittest.TestCase):
    def test_constant_potential_and_energy_shift_cancel(self):
        # r*delta V and delta E must share the same Ry convention. A uniform
        # potential shift with identical shifts of all radial energies leaves
        # the radial differential equation and its discrete construction fixed.
        for rel in [False, True]:
            v, kw, r = inputs(nlo=2, rel=rel)
            with NativeRadialContext(v, **kw) as c:
                dp, dm = c.action(0.11 * r, 0.11, np.full(2, 0.11))
                self.assertLess(np.max(abs(dp)), 2e-12)
                self.assertLess(np.max(abs(dm)), 2e-12)

    def test_joint_potential_energy_and_lo_derivatives(self):
        rng = np.random.default_rng(361)
        for rel in [False, True]:
            for lapw, nlo in [(True, 0), (True, 2), (False, 1), (False, 2)]:
                with self.subTest(rel=rel, lapw=lapw, nlo=nlo), ExitStack() as stack:
                    v, kw, r = inputs(lapw, nlo, rel)
                    c = stack.enter_context(NativeRadialContext(v, **kw))
                    dv = 0.07 * r * np.sin(1.3 * r)
                    de = 0.03
                    delo = 0.02 * np.cos(np.arange(nlo) + 1.0)
                    dp, dm = c.action(dv, de, delo)
                    h = 0.0005
                    ends = [
                        stack.enter_context(
                            NativeRadialContext(
                                v + s * h * dv,
                                **(
                                    kw
                                    | dict(
                                        energy=kw["energy"] + s * h * de,
                                        lo_energies=kw["lo_energies"] + s * h * delo,
                                    )
                                )
                            )
                        )
                        for s in [-1, 1]
                    ]
                    for tangent, p, m in [
                        (dp, ends[1].primitive, ends[0].primitive),
                        (dm, ends[1].lo_coefficients, ends[0].lo_coefficients),
                    ]:
                        fd = (p - m) / (2 * h)
                        error = np.linalg.norm(fd - tangent) / max(
                            np.linalg.norm(tangent), 1e-12
                        )
                        self.assertLess(error, 2e-5)
                    gp = rng.normal(size=dp.shape)
                    gm = rng.normal(size=dm.shape)
                    vb, eb, elb = c.pullback(gp, gm)
                    f = np.vdot(gp, dp).real + np.vdot(gm, dm).real
                    b = vb @ dv + eb * de + elb @ delo
                    self.assertLess(abs(f - b) / max(1.0, abs(f), abs(b)), 2e-11)

    def test_lifetime_array_order_and_thread_serialization(self):
        v, kw, r = inputs()
        dv = np.sin(r)
        with NativeRadialContext(v, **kw) as c:
            expected = c.action(dv)
            with ThreadPoolExecutor(max_workers=3) as pool:
                results = list(pool.map(lambda i: c.action(dv), range(9)))
            for result in results:
                for a, b in zip(result, expected):
                    np.testing.assert_array_equal(a, b)
            self.assertFalse(c.primitive.flags.writeable)
            self.assertFalse(c.lo_coefficients.flags.writeable)
            # C-order/noncontiguous user covectors must be copied safely.
            gp = np.arange(c.primitive.size).reshape(c.primitive.shape)[:, :, ::-1]
            gm = np.ones_like(c.lo_coefficients)
            a = c.pullback(gp, gm)
            b = c.pullback(np.asfortranarray(gp), np.asfortranarray(gm))
            for x, y in zip(a, b):
                np.testing.assert_array_equal(x, y)
        c.close()
        with self.assertRaises(RuntimeError):
            c.action(dv)
        with NativeRadialContext(v, **kw) as second:
            np.testing.assert_array_equal(second.action(dv)[0], expected[0])

    def test_invalid_or_unsupported_input_is_rejected(self):
        v, kw, r = inputs()
        for change in [
            dict(rmt=1.9),
            dict(ell=-1),
            dict(relativistic=1),
            dict(dx=np.nan),
        ]:
            with self.assertRaises(ValueError):
                NativeRadialContext(v, **(kw | change))
        with self.assertRaises(NotImplementedError):
            NativeRadialContext(v, dirac_lo=True, **kw)
        with NativeRadialContext(v, **kw) as c:
            with self.assertRaises(ValueError):
                c.action(v.astype(complex))
            with self.assertRaises(ValueError):
                c.action(v[:-1])
            with self.assertRaises(ValueError):
                c.action(v, energy_dot=np.nan)
            with self.assertRaises(ValueError):
                c.pullback(c.primitive, c.lo_coefficients[:, 0])

    def test_composed_basis_density_against_nonlinear_endpoints(self):
        rng = np.random.default_rng(57)
        v, kw, r = inputs(ell=0)
        nb = 3
        labels = np.array([[0, 0, 0], [1, 0, 0], [-1, 0, 0], [0, 1, 0]])
        q = labels * 0.4
        length = np.linalg.norm(q, axis=1)
        phase = np.exp(0.2j * labels[:, 0])
        coeff = rng.normal(size=(4, nb)) + 1j * rng.normal(size=(4, nb))
        channel = PreparedAugmentationChannel(
            0, kw["rmt"], 30.0, q, length, phase, coeff, 3, [np.array([3])]
        )
        rho = np.array(
            [herm(rng.normal(size=(nb, nb)) + 1j * rng.normal(size=(nb, nb)))]
        )
        stars = FourierStars(
            (9, 9, 9),
            np.array([[0, 0, 0]]),
            np.eye(3, dtype=int)[None],
            np.zeros((1, 3)),
        )
        ir = InterstitialDensity(coeff[:3], labels[:3], 30.0, stars)
        template = Density.__new__(Density)
        template.rad = {(0, 0, 0): np.zeros(len(v))}
        template.pw = {(0, 0, 0): 0j}

        def spatial(ctx):
            z = channel.value(
                ctx.primitive[-1], ctx.lo_coefficients, lapw=ctx.lapw
            ).reshape(-1, nb)
            mt = MuffinTinDensitySum(
                [(0, z, 1.0)], rho, {0: ctx.primitive[:, :2]}, [(0, 0)], 4e-22
            )
            return SpatialValenceDensity(template, {0: mt}, [ir])

        with ExitStack() as stack:
            c = stack.enter_context(NativeRadialContext(v, **kw))
            full = spatial(c)
            op = BasisDensityResponse(
                full, {0: {0: c}}, {0: [{0: channel}]}, energy_policy="fixed"
            )
            dv = 0.07 * r * np.sin(1.3 * r)
            h = 0.0005
            ends = [
                spatial(
                    stack.enter_context(NativeRadialContext(v + s * h * dv, **kw))
                ).value()
                for s in [-1, 1]
            ]
            fd = (ends[1] - ends[0]) / (2 * h)
            actual = op.action(op.pack({0: dv}))
            self.assertLess(np.linalg.norm(fd - actual) / np.linalg.norm(actual), 2e-6)
            self.assertLess(op.operator.check_transpose(), 1e-11)
            with self.assertRaises(NotImplementedError):
                BasisDensityResponse(
                    full, {0: {0: c}}, {0: [{0: channel}]}, energy_policy="automatic"
                )
            with self.assertRaises(ValueError):
                op.pack({})


if __name__ == "__main__":
    unittest.main()
