"""Retained-space and radial SOC derivatives independent of the Ce files."""

from contextlib import ExitStack
from types import SimpleNamespace
import os
from pathlib import Path
import unittest

import numpy as np
from scipy.linalg import eigh

from lapw_radial_native import NativeRadialContext
from lapw_soc_response import SOCRadialChannel, SOCMatching
from lapw_soc_subspace import RetainedEigenspace, hermitian


class SOCPointCacheTests(unittest.TestCase):
    def test_retained_prefix_and_streamed_resource_lifetimes(self):
        from lapw_soc_screened import SOCScreenedForce

        with ExitStack() as lifetime:
            model = SOCScreenedForce.__new__(SOCScreenedForce)
            model.stack = lifetime
            model.records = range(3)
            model.cached_points = {}
            model.cache_bytes = model.cache_scratch = 0
            model.cache_limit = 4 * 1024**2
            model.scratch_limit = 1024
            model.admit_points = True
            opened, closed = [], []

            def prepare(index, record, resources):
                opened.append(index)
                resources.callback(closed.append, index)
                return SimpleNamespace(array=np.ones(2 * 1024**2, dtype=np.uint8),
                                       val=SimpleNamespace(matching_scratch_bytes=512))

            model._prepare_point = prepare
            for _ in range(2):
                self.assertEqual([i for i, _ in model.points()], [0, 1, 2])
            self.assertEqual(opened, [0, 1, 2, 1, 2])
            self.assertEqual(closed, [1, 2, 1, 2])
            self.assertLessEqual(model.cache_bytes, model.cache_limit)
            self.assertEqual(model.cache_scratch, 512)
        self.assertEqual(closed[-1], 0)


class SOCSubspaceTests(unittest.TestCase):
    def test_degenerate_retained_block_and_metric(self):
        rng = np.random.default_rng(191)

        def random_matrix(n, m):
            return rng.normal(size=(n, m)) + 1j * rng.normal(size=(n, m))

        n = 6
        a = np.linalg.qr(random_matrix(n, n))[0]
        energy = np.array([-1.0, -1.0, 0.1, 0.8, 1.5, 2.4])
        o = np.eye(n)
        h = (a * energy) @ a.conj().T
        dh, do = (hermitian(random_matrix(n, n)) for _ in range(2))
        for nr in (2, 4, n):
            space = RetainedEigenspace(h, o, a[:, :nr], energy[:nr])
            band, da = space.action(dh, do)
            ar = a[:, :nr]
            tangent = da @ ar.conj().T + ar @ da.conj().T
            step = 1e-6

            def outer(t):
                _, vectors = eigh(h + t * dh, o + t * do)
                return vectors[:, :nr] @ vectors[:, :nr].conj().T

            np.testing.assert_allclose(
                tangent, (outer(step) - outer(-step)) / (2 * step), atol=3e-8
            )
            cov1, cov2 = hermitian(random_matrix(nr, nr)), random_matrix(n, nr)
            bh, bo = space.pullback(cov1, cov2)
            lhs = np.vdot(cov1, band).real + np.vdot(cov2, da).real
            self.assertAlmostEqual(
                lhs, np.vdot(bh, dh).real + np.vdot(bo, do).real, places=11
            )

    def test_uniform_dilation_in_triclinic_and_centered_cells(self):
        from types import SimpleNamespace
        from lapw_lattice import NativeLattice
        from lapw_soc_response import SOCSite

        rng = np.random.default_rng(731)
        labels = rng.integers(-3, 4, (19, 3)).astype(float)
        labels[0] = 0.0
        c = rng.normal(size=(19, 7)) + 1j * rng.normal(size=(19, 7))
        ctx = SimpleNamespace(rmt=2.35, lapw=True, primitive=rng.normal(size=(1, 4, 3)),
                              lo_coefficients=rng.normal(size=(3, 1)))
        for kind, angles in [("P", (76, 99, 104)), ("CXZ", (90, 90, 96)),
                             ("F", (90, 90, 90))]:
            def build(eta):
                a, b, caxis = np.exp(eta) * np.array([7., 8., 9.])
                s = SimpleNamespace(lattyp=kind, a=a, b=b, c=caxis,
                                    alpha=angles[0], beta=angles[1], gamma=angles[2])
                lattice = NativeLattice.from_structure(s)
                geometry = SimpleNamespace(br1=lattice.reciprocal, nlo=5, volume=lattice.volume,
                    channels={(0, 2): SimpleNamespace(ctx=ctx, slots=[0, 1, 2])},
                    lo_indices={(0, 2): [np.arange(5)]}, legacy_lo_filter=False)
                site = SOCSite(0, lattice.reciprocal, np.array([.13, -.27, .44]), {})
                return SOCMatching(geometry, site, labels, c, 2)
            with self.subTest(kind=kind):
                step = 1e-5
                dz = build(0).action(None, None, np.zeros_like(ctx.primitive[-1]),
                                     np.zeros_like(ctx.lo_coefficients), hydrostatic=True)
                fd = (build(step).z - build(-step).z) / (2 * step)
                np.testing.assert_allclose(dz, fd, atol=1e-8, rtol=1e-7)

    def test_lo_phase_and_differential_at_zero(self):
        # Small components and exact zeros occur even for normalized bands.
        z = np.array([0j, 0.9e-6 + 1.1e-6j, 2e-6 + 0j])
        dz = np.array([0.2 + 0.4j, 1 - 0.1j, 0.7j])
        phase = np.exp(0.78j)
        apply = SOCMatching.component_projection
        np.testing.assert_allclose(apply(phase * z, None), phase * apply(z, None))
        np.testing.assert_array_equal(apply(dz, None), dz)
        masks = (abs(z.real) > 1e-6, abs(z.imag) > 1e-6)
        np.testing.assert_allclose(apply(dz, masks), [0, -0.1j, 0])

    def test_covariant_soc_operator_in_a_moving_scalar_space(self):
        # Ordinary FOR holds the radial operator fixed apart from its explicit
        # geometry derivative. A retained metric leg counted twice here gives
        # a nonzero error even when the scalar space is complete.
        rng = np.random.default_rng(313)

        def random(n):
            return rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))

        n = 7
        y = random(n)
        o = y.conj().T @ y + np.eye(n)
        h, dh, do, w, dw = [hermitian(random(n)) for _ in range(5)]
        e, all_a = eigh(h, o)
        for nr in (3, n):
            a = all_a[:, :nr]
            scalar = RetainedEigenspace(h, o, a, e[:nr])
            b = scalar.outside
            cross = b @ (
                (b.conj().T @ (dh @ a - (do @ a) * scalar.energy))
                * scalar.reciprocal_gaps
            )
            wl = cross.conj().T @ w @ a
            expected = wl + wl.conj().T + a.conj().T @ dw @ a

            def matrix(t):
                ot = o + t * do
                ea, v = eigh(h + t * dh, ot)
                qr = v[:, :nr] @ v[:, :nr].conj().T @ ot
                return (qr @ a).conj().T @ (w + t * dw) @ (qr @ a)

            t = 1e-6
            np.testing.assert_allclose(
                (matrix(t) - matrix(-t)) / (2 * t), expected, atol=3e-9, rtol=1e-8
            )
            if nr == n:
                np.testing.assert_allclose(expected, a.conj().T @ dw @ a)

    def test_only_crossing_pair_is_omitted(self):
        e = np.array([-2.0, 1.0, 1.0, 4.0])
        s = RetainedEigenspace(np.diag(e), None, np.eye(4)[:, :2], e[:2])
        self.assertEqual(s.singular_pairs, 1)
        _, da = s.action(np.ones((4, 4)))
        self.assertEqual(da[2, 1], 0.0)
        self.assertAlmostEqual(da[3, 1], -1 / 3.0)


class SOCRadialTests(unittest.TestCase):
    def test_integrals_matching_and_pullback(self):
        lib = Path(
            os.environ.get(
                "EDMFT_RADIAL_LIBRARY",
                str(
                    Path(__file__).resolve().parents[2]
                    / "dmft2/atpar/liblapw_response.so"
                ),
            )
        )
        if not lib.is_file():
            self.skipTest("build native radial library for radial SOC test")
        for n in (300, 301):
            r = np.geomspace(0.0001, 2.0, n)
            v = -4 + 0.2 * r**2
            dv = 0.07 * r * np.sin(1.3 * r)
            de, dlo = 0.03, np.array([0.02, -0.03])

            def radial(t):
                return NativeRadialContext(
                    v + t * dv,
                    r0=r[0],
                    dx=np.log(r[-1] / r[0]) / (n - 1),
                    rmt=r[-1],
                    energy=0.3 + t * de,
                    lo_energies=np.array([0.3, -2.0]) + t * dlo,
                    ell=1,
                    charge=2.0,
                    relativistic=True,
                    lapw=False,
                    quadrature="lapwso",
                    library=lib,
                )

            step = 0.0001
            with radial(0) as ctx, radial(step) as plus, radial(-step) as minus:
                op = SOCRadialChannel(ctx, v)
                d = op.action(dv, de, dlo)
                sides = [
                    (
                        SOCRadialChannel(c, v + t * dv).matrix,
                        c.primitive[-1],
                        c.lo_coefficients,
                    )
                    for c, t in ((plus, step), (minus, -step))
                ]
                for j in range(3):
                    np.testing.assert_allclose(
                        d[j],
                        (sides[0][j] - sides[1][j]) / (2 * step),
                        rtol=2e-5,
                        atol=1e-7,
                    )
                rng = np.random.default_rng(155)
                cov = [rng.normal(size=x.shape) for x in d]
                gv, ge, glo = op.pullback(*cov)
                self.assertAlmostEqual(
                    sum(np.vdot(x, y).real for x, y in zip(d, cov)),
                    gv @ dv + ge * de + glo @ dlo,
                    places=9,
                )
