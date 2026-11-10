"""Independent dense/native-kernel controls for ragged RR/QR adapters."""

import copy
import unittest
import numpy as np
from lapw_compact_blocks import right_multiply_columns, pair_columns
from lapw_compact_response import (
    CompactSphericalPotentialBands,
    CompactSpatialValenceDensity,
)
from lapw_spherical_response import PreparedSphericalBand, prepare_wronskian_band
from lapw_nonspherical_response import PreparedNonspherical
from lapw_augmentation_response import PreparedAugmentationChannel
from lapw_projector_response import PreparedProjector5
from lapw_potential_response import NativePotentialBands, NativeInterstitialPotential
from lapw_spherical_potential_response import SphericalPotentialBands
from lapw_density_response import (
    MuffinTinDensitySum,
    SpatialValenceDensity,
    InterstitialDensity,
    FourierStars,
    herm,
)
from wien_response_density import Density
from test_lapw_potential_assembly import fixture
from test_lapw_density_response import fixture as density_fixture
from test_lapw_radial_native import inputs, LIBRARY
from lapw_radial_native import NativeRadialContext


def completion(columns):
    n, nr = columns.shape
    y = np.zeros((n, n), complex)
    y[:, :nr] = columns
    y[:nr, nr:] = columns[nr:].conj().T
    return y


def restrict_channel(ch, nr):
    ch = copy.copy(ch)
    ch.fj = ch.fj[:, :nr].copy()
    ch.dfj = ch.dfj[:, :nr].copy()
    ch.local = ch.local[:, :, :nr].copy()
    ch.nband = nr
    return ch


class CompactTests(unittest.TestCase):
    def test_pairing_and_right_product(self):
        rng = np.random.default_rng(12)
        for n, nr in [(9, 3), (4, 4)]:
            y = herm(rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n)))[:, :nr]
            a = rng.normal(size=(2, 5, n)) + 1j * rng.normal(size=(2, 5, n))
            np.testing.assert_allclose(
                right_multiply_columns(a, y), a @ completion(y), atol=3e-15
            )
            self.assertAlmostEqual(
                pair_columns(y, y),
                np.vdot(completion(y), completion(y)).real,
                places=12,
            )
        with self.assertRaises(ValueError):
            right_multiply_columns(a, np.zeros((3, 2)))

    def test_spherical_rectangular_native_triangle_and_reverse(self):
        rng, ch, p, lo, w, r, e, el, tri, geo = fixture()
        q, phase, c, b, vol, npw = geo
        nr = 2
        compact_tri = prepare_wronskian_band(
            ch.ell, ch.rmt, vol, q[:npw], phase[:npw], c[:npw], b, retained=nr
        )
        np.testing.assert_allclose(compact_tri, tri[:, :nr], atol=2e-14)
        dp = rng.normal(size=p.shape) * 0.1
        dm = rng.normal(size=lo.shape) * 0.1
        de = 0.03
        delocal = np.array([0.02, -0.01])
        dv = rng.normal(size=len(p)) * 0.1
        h = herm(rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4)))[:, :nr]
        o = herm(rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4)))[:, :nr]
        for lapw in [False, True]:
            full = PreparedSphericalBand(
                ch, p, lo, w, r, e, el, tri, lapw=lapw, difference=dv
            )
            small = PreparedSphericalBand(
                ch, p, lo, w, r, e, el, compact_tri, lapw=lapw, difference=dv
            )
            for left, right in zip(small.value(), full.value()):
                np.testing.assert_allclose(left, right[:, :nr], atol=3e-12)
            d = small.action(dp, dm, de, delocal, dv)
            for left, right in zip(d, full.action(dp, dm, de, delocal, dv)):
                np.testing.assert_allclose(left, right[:, :nr], atol=3e-12)
            for left, right in zip(
                small.pullback(h, o), full.pullback(completion(h), completion(o))
            ):
                np.testing.assert_allclose(left, right, atol=5e-12, rtol=2e-13)
            gp, gm, ge, gel, gv = small.pullback(h, o)
            lhs = pair_columns(h, d[0]) + pair_columns(o, d[1])
            rhs = (
                np.vdot(gp, dp).real
                + np.vdot(gm, dm).real
                + ge * de
                + gel @ delocal
                + gv @ dv
            )
            self.assertLess(abs(lhs - rhs) / max(1, abs(lhs), abs(rhs)), 2e-13)

    def test_nonspherical_rectangular_radial_response(self):
        rng, ch, p, lo, w, *_ = fixture()
        z = ch.value(p[-1], lo, lapw=True).reshape(-1, ch.nband)
        op = PreparedNonspherical(
            {1: p},
            w,
            [(2, 0), (-2, 1)],
            rng.normal(size=(len(p), 2)),
            [(z, {1: ch}, {1: lo}, {1: True})],
        )
        dp = {1: rng.normal(size=p.shape)}
        dm = {1: rng.normal(size=lo.shape)}
        dv = rng.normal(size=op.v.shape)
        a = op.action(dp, dm, dv, retained=2)[0]
        np.testing.assert_allclose(a, op.action(dp, dm, dv)[0][:, :2], atol=3e-12)
        y = herm(rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4)))[:, :2]
        small = op.pullback([y])
        large = op.pullback([completion(y)])
        for j in [0, 1]:
            np.testing.assert_allclose(small[j][1], large[j][1], atol=5e-12)
        np.testing.assert_allclose(small[2], large[2], atol=5e-12)

    def test_ragged_spatial_forward_reverse_and_central_density(self):
        rng, phi, z, rho, hs = density_fixture()
        nr = 2
        nqs = [3, 5]
        rr0 = np.array([rho[:nr, :nr], 0.7 * rho[:nr, :nr]])
        zrs = [z[:, :nr], 0.7 * z[:, :nr].conj(), 1.2 * z[:, :nr]]
        entries = [(0, zrs[0], 0.5), (0, zrs[1], 0.5), (1, zrs[2], 1.0)]
        qs = [
            rng.normal(size=(len(z), nqs[k])) + 1j * rng.normal(size=(len(z), nqs[k]))
            for k, _, _ in entries
        ]
        mt = MuffinTinDensitySum(entries, rr0, phi, hs, 0.0001)
        reps = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [-1, 0, 0]])
        stars = FourierStars(
            (9, 9, 9), reps, np.eye(3, dtype=int)[None], np.zeros((1, 3))
        )
        labels = np.array([[0, 0, 0], [1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, 0, 1]])
        cs = [
            rng.normal(size=(5, nr + nq)) + 1j * rng.normal(size=(5, nr + nq))
            for nq in nqs
        ]
        ir = [
            InterstitialDensity(
                c[:, :nr], labels, 123.0, stars, cache_wavefunctions=False
            )
            for c in cs
        ]
        template = Density.__new__(Density)
        template.rad = {(0, *h): np.zeros(7) for h in hs}
        template.pw = {tuple(g): 0j for g in reps}
        retained = SpatialValenceDensity(template, {0: mt}, ir)
        compact = CompactSpatialValenceDensity(
            retained, nqs, {0: qs}, [c[:, nr:] for c in cs]
        )
        np.testing.assert_array_equal(compact.value(), retained.value())
        vector = rng.normal(size=compact.coordinates.size)
        rr, qr = compact.coordinates.unpack(vector)
        # Independent full matrices at each k; never pad unequal complements.
        radial = np.zeros_like(mt.value())
        fourier = np.zeros(len(reps), complex)
        from lapw_density_response import MuffinTinDensity

        for (k, zr, weight), zq in zip(entries, qs):
            d = completion(np.concatenate((rr[k], qr[k])))
            fullz = np.column_stack((zr, zq))
            radial += weight * MuffinTinDensity(fullz, d, phi, hs, 0.0001).value()
        for k, c in enumerate(cs):
            fourier += InterstitialDensity(c, labels, 123.0, stars).value(
                completion(np.concatenate((rr[k], qr[k])))
            )
        expected = retained._assemble({0: radial}, fourier)
        np.testing.assert_allclose(
            compact.action(vector), expected, atol=3e-11, rtol=3e-13
        )
        self.assertLess(
            compact.operator.check_transpose(samples=3, tolerance=3e-12), 3e-12
        )
        # It remains the same physical map if retained FFT caching is enabled.
        for op in ir:
            op.psi = op._waves(op.c)
        np.testing.assert_allclose(
            compact.action(vector), expected, atol=3e-11, rtol=3e-13
        )
        with self.assertRaises(ValueError):
            CompactSpatialValenceDensity(
                retained, [3, 4], {0: qs}, [c[:, nr:] for c in cs]
            )

    @unittest.skipUnless(LIBRARY.is_file(), "native radial library required")
    def test_ragged_native_potential_matches_separate_dense_k_maps(self):
        rng = np.random.default_rng(937)
        nr = 2
        nbs = [5, 7]
        contexts = {0: {}}
        from contextlib import ExitStack

        with ExitStack() as stack:
            for l in [0, 1]:
                v, kw, r = inputs(nlo=1, ell=l)
                contexts[0][l] = stack.enter_context(NativeRadialContext(v, **kw))
            sph = []
            project = []
            entries = []
            full_spheres = []
            full_projects = []
            ir = []
            full_ns = []
            reps = np.array([[0, 0, 0], [1, 0, 0], [-1, 0, 0]])
            for k, nb in enumerate(nbs):
                labels = np.array(
                    [
                        [0, 0, 0],
                        [1, 0, 0],
                        [-1, 0, 0],
                        [0, 1, 0],
                        [0, 0, 1],
                        [1, 1, 0],
                        [1, 0, 1],
                        [-1, 1, 0],
                    ]
                )
                q = labels + rng.normal(size=3) * 0.05
                qn = np.linalg.norm(q, axis=1)
                phase = np.exp(1j * rng.normal(size=8))
                c = rng.normal(size=(8, nb)) + 1j * rng.normal(size=(8, nb))
                channels = {}
                zs = {}
                one = []
                for l, ctx in contexts[0].items():
                    ch = PreparedAugmentationChannel(
                        l, ctx.rmt, 137.0, q, qn, phase, c, 6, [np.array([6 + l])]
                    )
                    channels[l] = ch
                    from scipy.special import spherical_jn

                    match = np.array(
                        [
                            spherical_jn(l, qn[:6] * ctx.rmt),
                            qn[:6] * spherical_jn(l, qn[:6] * ctx.rmt, derivative=True),
                        ]
                    )
                    tri = prepare_wronskian_band(
                        l, ctx.rmt, 137.0, q[:6], phase[:6], c[:6], match
                    )
                    args = (
                        ch,
                        ctx.primitive,
                        ctx.lo_coefficients,
                        ctx.overlap_metric,
                        ctx.radius,
                        ctx.energy,
                        ctx.lo_energies,
                    )
                    one.append((0, l, 0, PreparedSphericalBand(*args, tri, lapw=True)))
                    sph.append(
                        (0, l, k, PreparedSphericalBand(*args, tri[:, :nr], lapw=True))
                    )
                    zs[l] = ch.value(
                        ctx.primitive[-1], ctx.lo_coefficients, lapw=True
                    ).reshape(-1, nb)
                entry = (
                    np.concatenate(list(zs.values())),
                    channels,
                    {l: ctx.lo_coefficients for l, ctx in contexts[0].items()},
                    {0: True, 1: True},
                )
                entries.append(entry)
                full_spheres.append(one)
                ctx = contexts[0][1]
                op = PreparedProjector5(
                    channels[1],
                    ctx.primitive,
                    ctx.lo_coefficients,
                    ctx.projector_metric,
                    ctx.primitive[:, :2, 0],
                    np.eye(3),
                    lapw=True,
                )
                project.append((0, 1, k, 0, slice(0, 3), op))
                full_projects.append([(0, 1, 0, 0, slice(0, 3), op)])
                ir.append(
                    NativeInterstitialPotential(
                        c[:6],
                        labels[:6],
                        reps,
                        np.eye(3, dtype=int)[None],
                        np.zeros((1, 3)),
                        complex_values=True,
                    )
                )
            phi = {l: c.primitive for l, c in contexts[0].items()}
            vns = np.sin(r)[:, None] * 0.01
            ns = PreparedNonspherical(
                phi, contexts[0][0].overlap_metric, [(2, 0)], vns, entries
            )
            b = CompactSphericalPotentialBands(
                contexts,
                sph,
                [(0, ns, [0, 1])],
                project,
                [n - nr for n in nbs],
                (2, nr, 3, 1),
                energy_policy="fixed",
            )
            new = NativePotentialBands(b, ir, retained=nr)
            x = rng.normal(size=new.input.size) * 0.01
            actual = new.output.unpack(new.action(x))
            for k, nb in enumerate(nbs):
                one_ns = PreparedNonspherical(
                    phi, contexts[0][0].overlap_metric, [(2, 0)], vns, [entries[k]]
                )
                full = SphericalPotentialBands(
                    contexts,
                    full_spheres[k],
                    [(0, one_ns, [0])],
                    full_projects[k],
                    (1, nb, nb),
                    (1, nb, 3, 1),
                    energy_policy="fixed",
                )
                old = NativePotentialBands(full, [ir[k]], retained=nr)
                expected = old.output.unpack(old.action(x))
                for left, right in zip(actual, expected):
                    np.testing.assert_allclose(
                        left[k], right[0], atol=4e-12, rtol=2e-12
                    )
            self.assertLess(
                new.operator.check_transpose(samples=2, tolerance=4e-11), 4e-11
            )
            # Zero spherical input exercises compact allocation without radial work.
            x[: b.input_size] = 0
            self.assertTrue(np.all(np.isfinite(new.action(x))))


if __name__ == "__main__":
    unittest.main()
