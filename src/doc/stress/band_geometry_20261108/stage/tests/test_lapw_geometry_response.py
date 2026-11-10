"""Geometry derivatives against independently rebuilt nonlinear LAPW blocks.

These tests include complex off-diagonal columns, angular motion/shear, ordinary
local orbitals and a deliberately non-unit radial Wronskian. Material/native
checks are separate so none of these tests needs a private calculation folder.
"""

from types import SimpleNamespace
import unittest
import numpy as np
from scipy.special import spherical_jn
from lapw_augmentation_response import PreparedAugmentationChannel, harmonics
from lapw_spherical_response import PreparedSphericalBand, prepare_wronskian_band
from lapw_geometry_response import (
    harmonic_direction,
    matching_direction,
    augmentation_geometry_direction,
    wronskian_geometry_direction,
    spherical_geometry_direction,
    interstitial_geometry,
    native_band_geometry_source,
)
from test_lapw_potential_assembly import fixture


def endpoint(ch, geometry, strain, t):
    """Fresh nonlinear matching at the deformed lattice, fixed C and RMT."""
    q, phase, c, _, volume, npw = geometry
    deformation = np.eye(3) + t * strain
    qp = q @ np.linalg.inv(deformation)
    vp = volume * np.linalg.det(deformation)
    lengths = np.linalg.norm(qp, axis=1)
    groups = [np.array([npw]), np.array([npw + 1])]
    cp = PreparedAugmentationChannel(
        ch.ell, ch.rmt, vp, qp, lengths, phase, c, npw, groups
    )
    b = np.array(
        [
            spherical_jn(ch.ell, lengths[:npw] * ch.rmt),
            lengths[:npw]
            * spherical_jn(ch.ell, lengths[:npw] * ch.rmt, derivative=True),
        ]
    )
    tri = prepare_wronskian_band(
        ch.ell, ch.rmt, vp, qp[:npw], phase[:npw], c[:npw], b, retained=2
    )
    return cp, tri


def tangent(ch, geometry, strain):
    q, phase, c, _, volume, npw = geometry
    dq = -q @ strain
    lengths = np.linalg.norm(q, axis=1)
    dl = np.einsum("ij,ij->i", q, dq) / lengths
    return augmentation_geometry_direction(
        ch,
        q,
        dq,
        lengths,
        dl,
        phase,
        c,
        npw,
        [np.array([npw]), np.array([npw + 1])],
        volume,
        np.trace(strain),
    )


class GeometryTests(unittest.TestCase):
    def test_harmonics_cartesian_and_radial(self):
        rng = np.random.default_rng(903)
        q = rng.normal(size=(11, 3))
        dq = rng.normal(size=q.shape)
        h = 2e-6
        for ell in range(6):
            y, dy = harmonic_direction(ell, q, dq)
            ref = (harmonics(ell, q + h * dq) - harmonics(ell, q - h * dq)) / (2 * h)
            np.testing.assert_allclose(dy, ref, atol=2e-9, rtol=2e-8)
            np.testing.assert_array_equal(y, harmonics(ell, q))
            _, radial = harmonic_direction(ell, q, -q)
            np.testing.assert_allclose(radial, 0, atol=2e-15)

    def test_harmonics_poles_against_cartesian_l1(self):
        q = np.array(
            [[0.0, 0.0, 2.0], [0.0, 0.0, -3.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
        )
        dq = np.array(
            [[0.3, -0.2, 0.1], [-0.2, 0.5, -0.1], [0.2, 0.1, -0.4], [0.0, 0.0, 0.0]]
        )
        _, actual = harmonic_direction(1, q, dq)
        length = np.linalg.norm(q, axis=1)
        n = np.divide(
            q, length[:, None], out=np.zeros_like(q), where=length[:, None] > 0
        )
        dn = np.divide(
            dq - n * np.sum(n * dq, axis=1)[:, None],
            length[:, None],
            out=np.zeros_like(q),
            where=length[:, None] > 0,
        )
        a = np.sqrt(3 / (8 * np.pi))
        expected = np.array(
            [
                a * (dn[:, 0] - 1j * dn[:, 1]),
                np.sqrt(2) * a * dn[:, 2],
                -a * (dn[:, 0] + 1j * dn[:, 1]),
            ]
        )
        np.testing.assert_allclose(actual, expected, atol=2e-16)

    def test_matching_zero_and_bessel_zero(self):
        q = np.array([0.0, np.pi / 2, 0.37, 3.4])
        dq = np.array([0.0, 0.13, -0.22, 0.4])
        h = 2e-6
        for ell in range(5):
            match, actual = matching_direction(ell, 2.0, q, dq)

            def values(t):
                p = q + t * dq
                return np.array(
                    [
                        spherical_jn(ell, 2 * p),
                        p * spherical_jn(ell, 2 * p, derivative=True),
                    ]
                )

            np.testing.assert_allclose(match, values(0), atol=0)
            np.testing.assert_allclose(
                actual, (values(h) - values(-h)) / (2 * h), atol=2e-10, rtol=2e-8
            )

    def test_augmentation_lapw_apw_and_lo(self):
        strain = np.array([[0.3, 0.11, -0.12], [0.11, -0.2, 0.05], [-0.12, 0.05, 0.4]])
        h = 2e-6
        for ell in range(4):
            _, ch, p, lo, _, _, _, _, _, geo = fixture(ell)
            dot, _ = tangent(ch, geo, strain)
            m, _ = endpoint(ch, geo, strain, -h)
            plus, _ = endpoint(ch, geo, strain, h)
            for lapw in [False, True]:
                actual = dot.value(p[-1], lo, lapw=lapw)
                expected = (
                    plus.value(p[-1], lo, lapw=lapw) - m.value(p[-1], lo, lapw=lapw)
                ) / (2 * h)
                np.testing.assert_allclose(actual, expected, atol=3e-8, rtol=2e-7)

    def test_ordered_triangle_not_averaged(self):
        _, ch, _, _, _, _, _, _, _, geo = fixture()
        strain = np.array([[0.4, 0.2, 0], [0.2, -0.3, 0.1], [0, 0.1, 0.2]])
        _, data = tangent(ch, geo, strain)
        c, npw = geo[2], geo[-1]
        actual = wronskian_geometry_direction(data, c[:npw], 2)
        h = 2e-6
        expected = (
            endpoint(ch, geo, strain, h)[1] - endpoint(ch, geo, strain, -h)[1]
        ) / (2 * h)
        np.testing.assert_allclose(actual, expected, atol=3e-9, rtol=2e-8)
        self.assertGreater(np.linalg.norm(actual[2:]), 1e-3)
        # Averaging this coefficient kernel would erase its antisymmetric
        # Wronskian radial factor. The native ordered-copy derivative is nonzero.
        np.testing.assert_allclose(actual[:2], actual[:2].conj().T, atol=2e-15)

    def test_spherical_both_legs_surface_and_nonunit_wronskian(self):
        _, ch, p, lo, w, r, e, el, tri, geo = fixture()
        strain = np.array([[0.2, 0.1, 0], [0.1, -0.3, 0.07], [0, 0.07, 0.4]])
        dot, data = tangent(ch, geo, strain)
        dtri = wronskian_geometry_direction(data, geo[2][: geo[-1]], 2)
        h = 2e-6
        for lapw in [False, True]:
            op = PreparedSphericalBand(ch, p, lo, w, r, e, el, tri[:, :2], lapw=lapw)
            self.assertGreater(abs(op.w - 1), 0.1)
            actual = spherical_geometry_direction(op, dot, dtri)
            ends = []
            for t in [-h, h]:
                cp, tp = endpoint(ch, geo, strain, t)
                ends.append(
                    PreparedSphericalBand(cp, p, lo, w, r, e, el, tp, lapw=lapw).value()
                )
            for value, plus, minus in zip(actual, ends[1], ends[0]):
                np.testing.assert_allclose(
                    value, (plus - minus) / (2 * h), atol=3e-7, rtol=3e-8
                )

    def test_interstitial_mask_and_kinetic_shear(self):
        rng = np.random.default_rng(198)
        labels = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, -1], [-1, 0, 1]])
        q = 2 * np.pi * (labels + np.array([0.1, 0.2, 0.3])) / np.array([7.0, 8.0, 9.0])
        c = rng.normal(size=(5, 4)) + 1j * rng.normal(size=(5, 4))
        pos = np.array([[0.1, 0.2, 0.3], [0.4, 0.5, 0.6]])
        radii = np.array([1.6, 1.9])
        volume = 504.0
        strain = np.array([[0.2, 0.1, 0], [0.1, -0.3, 0.08], [0, 0.08, 0.4]])
        actual = interstitial_geometry(q, labels, c, pos, radii, volume, strain, 2)
        h = 2e-6
        ends = []
        for t in [-h, h]:
            f = np.eye(3) + t * strain
            ends.append(
                interstitial_geometry(
                    q @ np.linalg.inv(f),
                    labels,
                    c,
                    pos,
                    radii,
                    volume * np.linalg.det(f),
                    np.zeros((3, 3)),
                    2,
                )[2:]
            )
        for value, plus, minus in zip(actual[:2], ends[1], ends[0]):
            np.testing.assert_allclose(
                value, (plus - minus) / (2 * h), atol=2e-9, rtol=3e-8
            )
        # Native masked potential coefficients at fixed integer G labels
        # do not change under this geometry derivative.
        before = actual[2:]
        self.assertEqual(before[0].shape, (4, 2))

    def test_unsupported_strain_is_rejected_before_native_loading(self):
        rot = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
        state = SimpleNamespace(
            structure=SimpleNamespace(timat=[rot]),
            bands=SimpleNamespace(spherical=None),
        )
        with self.assertRaisesRegex(ValueError, "symmetric strain"):
            native_band_geometry_source(
                state, "/missing", "case", np.array([[0, 1, 0], [0, 0, 0], [0, 0, 0]])
            )
        with self.assertRaisesRegex(NotImplementedError, "symmetry"):
            native_band_geometry_source(
                state, "/missing", "case", np.diag([1.0, 0.0, 0.0])
            )


if __name__ == "__main__":
    unittest.main()
