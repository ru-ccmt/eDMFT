"""Complex, off-diagonal product rules and native-triangle assembly controls."""

import unittest
import numpy as np
from scipy.special import spherical_jn
from lapw_augmentation_response import PreparedAugmentationChannel, harmonics
from lapw_response_assembly import NativeBandProjection
from lapw_spherical_response import PreparedSphericalBand, prepare_wronskian_band
from lapw_projector_response import PreparedProjector5
from lapw_nonspherical_response import PreparedNonspherical
from lapw_density_response import herm


def fixture(l=1):
    rng = np.random.default_rng(923 + l)
    r = 2.0
    vol = 37.0
    npw = 5
    nb = 4
    q = rng.normal(size=(7, 3))
    qn = np.linalg.norm(q, axis=1)
    phase = np.exp(1j * q @ np.array([0.2, 0.3, 0.4]))
    c = rng.normal(size=(7, nb)) + 1j * rng.normal(size=(7, nb))
    b = np.array(
        [
            spherical_jn(l, qn[:npw] * r),
            qn[:npw] * spherical_jn(l, qn[:npw] * r, derivative=True),
        ]
    )
    ch = PreparedAugmentationChannel(
        l,
        r,
        vol,
        q,
        qn,
        phase,
        c,
        npw,
        [np.array([5]), np.array([6])],
        matching_values=b,
    )
    p = rng.normal(size=(7, 4, 4))
    lo = rng.normal(size=(4, 2))
    w = np.zeros((7, 4))
    w[:, :2] = rng.uniform(0.1, 0.3, (7, 2))
    radius = np.geomspace(0.01, r, 7)
    e = 0.2
    el = np.array([-0.7, -2.0])
    triangle = prepare_wronskian_band(l, r, vol, q[:npw], phase[:npw], c[:npw], b)
    return rng, ch, p, lo, w, radius, e, el, triangle, (q, phase, c, b, vol, npw)


class PotentialAssemblyTests(unittest.TestCase):
    def test_spherical_value_matches_full_native_triangle(self):
        rng, ch, p, lo, w, radius, e, el, triangle, geo = fixture()
        q, phase, c, b, vol, npw = geo
        difference = rng.normal(size=7)
        for lapw in [False, True]:
            op = PreparedSphericalBand(
                ch, p, lo, w, radius, e, el, triangle, lapw=lapw, difference=difference
            )
            m = np.zeros((4, len(c)))
            m[:, npw:] = ch.rmt**2 * lo
            if lapw:
                m[0, :npw] = ch.rmt * (p[-1, 0, 1] * b[1] - p[-1, 2, 1] * b[0])
                m[1, :npw] = ch.rmt * (p[-1, 2, 0] * b[0] - p[-1, 0, 0] * b[1])
            else:
                m[0, :npw] = ch.rmt * b[0] / p[-1, 0, 0]
            s = np.einsum("rci,rc,rcj->ij", p, w, p)
            s[0, 0] = 1
            s[0, 1] = s[1, 0] = 0
            t = np.diag(np.r_[e, e, el])
            t[0, 1] = 1
            k = (s @ t + t.T @ s) / 2 + np.einsum(
                "rci,rc,rcj,r->ij", p, w, p, difference / radius
            )
            pp = p[-1, 0] @ m / ch.rmt
            pp[npw:] = 0
            dd = p[-1, 2] @ m / ch.rmt
            hr = m.T @ k @ m + ch.rmt**2 / 2 * (np.outer(pp, dd) + np.outer(dd, pp))
            orr = m.T @ s @ m
            wr = p[-1, 2, 0] * p[-1, 0, 1] - p[-1, 0, 0] * p[-1, 2, 1]
            hr[:npw, :npw] += (
                0.5
                * (wr - 1)
                * (np.outer(m[0, :npw], m[1, :npw]) - np.outer(m[1, :npw], m[0, :npw]))
            )
            geometry = (
                4
                * np.pi
                / np.sqrt(vol)
                * 1j**ch.ell
                * harmonics(ch.ell, q).conj()
                * phase
            )
            projection = NativeBandProjection(geometry.conj().T @ geometry, c)
            for actual, expected in zip(
                op.value(), (projection.project(hr), projection.project(orr))
            ):
                np.testing.assert_allclose(actual, expected, rtol=3e-13, atol=2e-11)

    def test_spherical_joint_derivative_and_transpose(self):
        rng, ch, p, lo, w, radius, e, el, triangle, _ = fixture()
        u = rng.normal(size=p.shape) * 0.1
        dm = rng.normal(size=lo.shape) * 0.1
        v = rng.normal(size=7) * 0.1
        dv = rng.normal(size=7) * 0.1
        de = 0.04
        delo = np.array([0.03, -0.02])
        hb = herm(rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4)))
        ob = herm(rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4)))
        for lapw in [False, True]:
            op = PreparedSphericalBand(
                ch, p, lo, w, radius, e, el, triangle, lapw=lapw, difference=v
            )
            h = 1e-5
            ends = [
                PreparedSphericalBand(
                    ch,
                    p + sign * h * u,
                    lo + sign * h * dm,
                    w,
                    radius,
                    e + sign * h * de,
                    el + sign * h * delo,
                    triangle,
                    lapw=lapw,
                    difference=v + sign * h * dv,
                ).value()
                for sign in [-1, 1]
            ]
            dh, do = op.action(u, dm, de, delo, dv)
            for actual, plus, minus in zip((dh, do), ends[1], ends[0]):
                np.testing.assert_allclose(
                    actual, (plus - minus) / (2 * h), rtol=2e-8, atol=5e-8
                )
            gp, gm, ge, gel, gv = op.pullback(hb, ob)
            lhs = np.vdot(hb, dh).real + np.vdot(ob, do).real
            rhs = (
                np.vdot(gp, u).real
                + np.vdot(gm, dm).real
                + ge * de
                + gel @ delo
                + gv @ dv
            )
            self.assertLess(abs(lhs - rhs) / max(1.0, abs(lhs), abs(rhs)), 2e-13)

    def test_projector_radial_matching_and_complex_frame(self):
        rng, ch, p, lo, w, radius, e, el, triangle, _ = fixture()
        fixed = rng.normal(size=(7, 2))
        transform = rng.normal(size=(3, 2)) + 1j * rng.normal(size=(3, 2))
        tail = rng.normal(size=(4, 2)) + 1j * rng.normal(size=(4, 2))
        u = rng.normal(size=p.shape) * 0.1
        dm = rng.normal(size=lo.shape) * 0.1
        op = PreparedProjector5(ch, p, lo, w, fixed, transform, lapw=True, tail=tail)
        y = rng.normal(size=op.shape) + 1j * rng.normal(size=op.shape)
        h = 1e-5
        ends = [
            PreparedProjector5(
                ch,
                p + s * h * u,
                lo + s * h * dm,
                w,
                fixed,
                transform,
                lapw=True,
                tail=tail,
            ).value()
            for s in [-1, 1]
        ]
        actual = op.action(u, dm)
        np.testing.assert_allclose(
            actual, (ends[1] - ends[0]) / (2 * h), rtol=1e-8, atol=2e-8
        )
        gp, gm = op.pullback(y)
        self.assertAlmostEqual(
            np.vdot(y, actual).real,
            np.vdot(gp, u).real + np.vdot(gm, dm).real,
            places=11,
        )

    def test_nonspherical_both_radial_legs_potential_and_matching(self):
        rng, ch, p, lo, w, radius, e, el, triangle, _ = fixture(l=1)
        _, ch0, p0, lo0, *_ = fixture(l=0)
        phi = {0: p0, 1: p}
        los = {0: lo0, 1: lo}
        chs = {0: ch0, 1: ch}
        lapw = {0: False, 1: True}
        hs = [(1, 0), (2, 1), (-2, 1)]
        v = rng.normal(size=(7, 3))
        dp = {l: rng.normal(size=a.shape) * 0.1 for l, a in phi.items()}
        dm = {l: rng.normal(size=a.shape) * 0.1 for l, a in los.items()}
        dv = rng.normal(size=v.shape) * 0.1

        def prepare(ps, ls, vs):
            z = np.concatenate(
                [
                    chs[l].value(ps[l][-1], ls[l], lapw=lapw[l]).reshape(-1, 4)
                    for l in ps
                ]
            )
            return PreparedNonspherical(ps, w, hs, vs, [(z, chs, ls, lapw)])

        op = prepare(phi, los, v)
        h = 1e-5
        ends = [
            prepare(
                {l: a + s * h * dp[l] for l, a in phi.items()},
                {l: a + s * h * dm[l] for l, a in los.items()},
                v + s * h * dv,
            ).value()[0]
            for s in [-1, 1]
        ]
        actual = op.action(dp, dm, dv)[0]
        np.testing.assert_allclose(
            actual, (ends[1] - ends[0]) / (2 * h), rtol=1e-8, atol=2e-8
        )
        y = herm(rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4)))
        gp, gm, gv = op.pullback([y])
        rhs = (
            sum(np.vdot(gp[l], dp[l]).real + np.vdot(gm[l], dm[l]).real for l in phi)
            + np.vdot(gv, dv).real
        )
        self.assertAlmostEqual(np.vdot(y, actual).real, rhs, places=10)


if __name__ == "__main__":
    unittest.main()
