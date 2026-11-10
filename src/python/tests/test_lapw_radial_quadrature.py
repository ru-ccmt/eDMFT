"""Native H/O versus ATPARN quadrature and storage-policy regression controls."""

from contextlib import ExitStack
import unittest
import numpy as np
from test_lapw_radial_native import inputs
from test_lapw_potential_assembly import fixture
from lapw_radial_native import NativeRadialContext
from lapw_spherical_response import PreparedSphericalBand
from lapw_density_response import herm
from lapw_hartree_periodic import native_charge_weights


class QuadratureTests(unittest.TestCase):
    def test_native_metric_and_legacy_default(self):
        for rel in [False, True]:
            v, kw, r = inputs(nlo=2, rel=rel)
            with NativeRadialContext(v, **kw) as old, NativeRadialContext(
                v, quadrature="atparn", **kw
            ) as legacy, NativeRadialContext(v, quadrature="wien2k24", **kw) as native:
                np.testing.assert_array_equal(old.primitive, legacy.primitive)
                np.testing.assert_array_equal(old.overlap_metric, legacy.overlap_metric)
                np.testing.assert_allclose(
                    native.overlap_metric[:, 0],
                    native_charge_weights(r),
                    atol=3e-15,
                    rtol=2e-14,
                )
                self.assertGreater(
                    np.linalg.norm(native.primitive - old.primitive), 1e-7
                )
                np.testing.assert_allclose(
                    np.sum(native.overlap_metric * native.primitive[:, :, 0] ** 2),
                    1,
                    atol=3e-15,
                )
                np.testing.assert_allclose(
                    np.sum(
                        native.overlap_metric
                        * native.primitive[:, :, 0]
                        * native.primitive[:, :, 1]
                    ),
                    0,
                    atol=3e-15,
                )

    def test_native_quadrature_derivative_reverse_and_gauge(self):
        rng = np.random.default_rng(756)
        for rel in [False, True]:
            for lapw, nlo in [(True, 0), (True, 2), (False, 2)]:
                v, kw, r = inputs(lapw, nlo, rel)
                kw["quadrature"] = "wien2k24"
                with ExitStack() as stack:
                    c = stack.enter_context(NativeRadialContext(v, **kw))
                    dv = 0.07 * r * np.sin(1.3 * r)
                    de = 0.03
                    dlo = 0.02 * np.cos(np.arange(nlo) + 1.0)
                    dp, dm = c.action(dv, de, dlo)
                    h = 0.0005
                    ends = [
                        stack.enter_context(
                            NativeRadialContext(
                                v + t * h * dv,
                                **(
                                    kw
                                    | dict(
                                        energy=kw["energy"] + t * h * de,
                                        lo_energies=kw["lo_energies"] + t * h * dlo,
                                    )
                                )
                            )
                        )
                        for t in [-1, 1]
                    ]
                    for tangent, plus, minus in [
                        (dp, ends[1].primitive, ends[0].primitive),
                        (dm, ends[1].lo_coefficients, ends[0].lo_coefficients),
                    ]:
                        fd = (plus - minus) / (2 * h)
                        self.assertLess(
                            np.linalg.norm(tangent - fd)
                            / max(1, np.linalg.norm(tangent)),
                            2e-7,
                        )
                    p = rng.normal(size=dp.shape)
                    m = rng.normal(size=dm.shape)
                    vb, eb, lb = c.pullback(p, m)
                    lhs = np.vdot(p, dp).real + np.vdot(m, dm).real
                    rhs = vb @ dv + eb * de + lb @ dlo
                    self.assertLess(abs(lhs - rhs) / max(1, abs(lhs), abs(rhs)), 2e-11)
                    p, m = c.action(r, 1.0, np.ones(nlo))
                    self.assertLess(
                        max(np.max(abs(p)), np.max(abs(m), initial=0)), 2e-11
                    )

    def test_low_memory_spherical_value_action_reverse(self):
        rng, ch, p, lo, w, r, e, el, tri, _ = fixture()
        for a in [p, lo, w, r]:
            a.setflags(write=False)
        dp = rng.normal(size=p.shape)
        dm = rng.normal(size=lo.shape)
        y = herm(rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4)))[:, :2]
        for lapw in [False, True]:
            ref = PreparedSphericalBand(ch, p, lo, w, r, e, el, tri[:, :2], lapw=lapw)
            small = PreparedSphericalBand(
                ch,
                p,
                lo,
                w,
                r,
                e,
                el,
                tri[:, :2],
                lapw=lapw,
                cache_products=False,
                share_radial=True,
            )
            self.assertIs(small.phi, p)
            self.assertIs(small.metric, w)
            self.assertIsNone(small._kz)
            self.assertIsNone(small._sz)
            if not lapw:
                self.assertEqual(small.triangle.base.nbytes, 16)
            for a, b in zip(ref.value(), small.value()):
                np.testing.assert_allclose(a, b, rtol=0, atol=2e-13)
            for a, b in zip(ref.action(dp, dm), small.action(dp, dm)):
                np.testing.assert_allclose(a, b, rtol=0, atol=2e-13)
            for a, b in zip(ref.pullback(y, y), small.pullback(y, y)):
                np.testing.assert_allclose(a, b, rtol=0, atol=2e-13)
        with self.assertRaisesRegex(ValueError, "read-only"):
            PreparedSphericalBand(
                ch, p.copy(), lo, w, r, e, el, tri, lapw=True, share_radial=True
            )

    def test_unknown_policy_rejected(self):
        v, kw, r = inputs()
        with self.assertRaises(ValueError):
            NativeRadialContext(v, quadrature="guess", **kw)
