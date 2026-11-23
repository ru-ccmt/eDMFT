"""Transpose checks for shared Cartesian force gradients and deferred native exports."""

from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np
from lapw_atomic_geometry import matching_position_gradient, plane_wave_position_gradients
from lapw_augmentation_response import PreparedAugmentationChannel
from lapw_compact_blocks import hermitian_columns, pair_columns
from lapw_complement_response import ComplementMatrixCoordinates, ComplementPerturbationCoordinates
from lapw_geometry_response import (
    atomic_interstitial_columns,
    atomic_wronskian_columns,
    augmentation_geometry_direction,
)
import lapw_potential_core_state as archive
from test_lapw_potential_core_state import inputs


class AtomicGeometryTests(unittest.TestCase):
    def setUp(self):
        self.rng = np.random.default_rng(14290)

    def complex(self, *shape):
        return self.rng.normal(size=shape) + 1j * self.rng.normal(size=shape)

    def test_matching_apw_lapw_with_local_orbitals(self):
        rng = self.rng
        q = rng.normal(size=(20, 3))
        q[0] = 0
        phase = np.exp(1j * rng.normal(size=len(q)))
        lengths = np.linalg.norm(q, axis=1)
        rates = q - [0.2, 0.1, -0.3]
        c = self.complex(20, 6)
        groups = [np.arange(12, 15), np.arange(15, 20)]
        ell, radius, volume = 2, 1.9, 385.0
        channel = PreparedAugmentationChannel(ell, radius, volume, q, lengths, phase, c, 12, groups)
        boundary = rng.normal(size=(4, 4))
        lo = rng.normal(size=(4, 2))
        bar, pw = (self.complex(5, 4, 6) for _ in range(2))
        for lapw in [False, True]:
            gradient = matching_position_gradient(
                ell,
                radius,
                volume,
                q,
                lengths,
                phase,
                rates,
                c,
                12,
                groups,
                boundary,
                lo,
                lapw,
                bar,
                pw,
            )
            for velocity in np.eye(3):
                dz, _ = augmentation_geometry_direction(
                    channel,
                    q,
                    q * 0,
                    lengths,
                    lengths * 0,
                    phase,
                    c,
                    12,
                    groups,
                    volume,
                    0,
                    phase_rate=rates @ velocity,
                )
                reference = (
                    np.vdot(bar, dz.value(boundary, lo, lapw=lapw))
                    + np.vdot(pw, dz.value(boundary, lo * 0, lapw=lapw))
                ).real
                np.testing.assert_allclose(gradient @ velocity, reference, atol=2e-12)

    def test_plane_waves_and_ordered_surface_all_cartesian_directions(self):
        rng = self.rng
        labels = rng.integers(-5, 6, (34, 3))
        q = labels @ np.array([[0.9, 0.1, 0.0], [0.1, 1.0, 0.2], [0.0, 0.3, 0.7]])
        q[0] = 0
        c = self.complex(34, 12)
        nr, volume = 4, 356
        s = SimpleNamespace(nat=2, pos=rng.random((4, 3)), mult=[2, 2], rmt=[1.6, 2.0])
        weights = [{0: -0.1, 1: 0.8, 5: 0.14}, {2: 0.3}]
        active = np.array([True, False, True, True])
        bars = [
            (hermitian_columns(self.complex(12, nr)), hermitian_columns(self.complex(12, nr)), None)
            for _ in range(2)
        ]
        grad = plane_wave_position_gradients(q, labels, c, s, volume, weights, bars, active)
        np.testing.assert_array_equal(grad[:, ~active], 0)
        for site in np.flatnonzero(active):
            for axis in range(3):
                v = np.zeros((4, 3))
                v[site, axis] = 1
                h, o = atomic_interstitial_columns(
                    q, labels, c, s.pos, np.repeat(s.rmt, s.mult), volume, v, nr
                )
                h += atomic_wronskian_columns(q, labels, c, s, v, volume, weights, nr)
                ref = [pair_columns(hbar, h) + pair_columns(obar, o) for hbar, obar, _ in bars]
                np.testing.assert_allclose(grad[:, site, axis], ref, atol=2e-11)

    def test_packed_matrix_stream_variable_complement(self):
        coords = ComplementPerturbationCoordinates(
            ComplementMatrixCoordinates((3, 4, 4), [2, 5, 3]), (3, 4, 2, 3)
        )
        flat = self.rng.normal(size=coords.size)
        h, o, r, hq, oq = coords.unpack(flat)
        for k, (hk, ok, rk) in enumerate(coords.point_covectors(flat)):
            np.testing.assert_array_equal(hk, np.vstack([h[k], hq[k]]))
            np.testing.assert_array_equal(ok, np.vstack([o[k], oq[k]]))
            np.testing.assert_array_equal(rk, r[k])

    def test_deferred_archive_builds_and_validates_on_load(self):
        with TemporaryDirectory() as tmp:
            path = Path(tmp) / "central.npz"
            with patch.object(archive, "_prepare", wraps=archive._prepare) as prepare:
                self.assertIsNone(
                    archive.save_potential_core_state(path, **inputs(), prepare=False)
                )
                prepare.assert_not_called()
                central = archive.load_potential_core_state(path)
                prepare.assert_called_once()
            args = inputs()
            expected = archive.save_potential_core_state(path, **args)
            v = self.rng.normal(size=central.hartree.size)
            np.testing.assert_array_equal(
                central.potential_operator.matvec(v), expected.potential_operator.matvec(v)
            )
            args["angular_data"][0]["weights"][0] = np.nan
            archive.save_potential_core_state(path, **args, prepare=False)
            with self.assertRaises(ValueError):
                archive.load_potential_core_state(path)

    def test_direct_band_product_and_fft_agree_with_nonsymmorphic_stars(self):
        from lapw_density_response import FourierStars, InterstitialDensity

        eye = np.eye(3, dtype=int)
        reps = np.array([[0, 0, 0], [8, 0, 0], [1, 2, 0], [-2, 1, 1], [0, 1, 2]])
        rot = np.array([eye, eye, -eye, -eye])
        tau = np.array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.0], [0.25, 0.0, 0.25], [0.75, 0.5, 0.25]])
        grid = np.indices((5, 5, 5)).reshape(3, -1).T - 2
        labels = self.rng.permutation(grid)[:70]
        c = self.complex(len(labels), 19)
        y = self.complex(len(reps))
        for complex_phases in (False, True):
            stars = FourierStars((16, 20, 18), reps, rot, tau, complex_phases=complex_phases)
            for cache in (False, True):
                op = InterstitialDensity(c, labels, 317.31, stars, cache_wavefunctions=cache)
                self.assertIsNotNone(op.direct_labels)
                candidate = op.potential_times_coefficients(y)
                op.direct_labels = None
                np.testing.assert_allclose(
                    candidate, op.potential_times_coefficients(y), rtol=8e-13, atol=1e-15
                )

    def test_six_strains_and_trace_against_dense_ordered_surface(self):
        from lapw_geometry_response import (
            strain_wronskian_columns,
            wronskian_geometry_direction,
            _matching_factors,
        )

        n, nr, v = 27, 9, 345.0
        g = self.rng.integers(-6, 7, (n, 3))
        g[0], g[1] = [0, 0, 0], [0, 0, 2]
        q = g @ np.array([[1.5, 0.1, 0], [0, 0.7, 0.2], [0, 0, 0.8]])
        s = SimpleNamespace(pos=self.rng.random((3, 3)), rmt=[2.1, 1.8], mult=[2, 1])
        weights = [{0: 0.7, 1: 1.3, 5: 0.03}, {3: -0.6, 8: 0.45}]
        c = self.complex(n, 15)
        eye = np.eye(3)
        strains = [
            (np.outer(a, b) + np.outer(b, a)) / 2 for i, a in enumerate(eye) for b in eye[i:]
        ]
        for e in [*strains, eye]:
            length = np.linalg.norm(q, axis=1)
            dq = -q @ e
            dl = np.divide((q * dq).sum(1), length, out=np.zeros_like(length), where=length != 0)
            reference = np.zeros((c.shape[1], nr), complex)
            for site, a in enumerate(np.repeat([0, 1], s.mult)):
                for ell, weight in weights[a].items():
                    ch = SimpleNamespace(ell=ell, rmt=s.rmt[a])
                    data = _matching_factors(
                        ch,
                        q,
                        dq,
                        length,
                        dl,
                        np.exp(2j * np.pi * (g @ s.pos[site])),
                        n,
                        v,
                        np.trace(e),
                    )
                    reference += weight * wronskian_geometry_direction(data, c, nr)
            np.testing.assert_allclose(
                strain_wronskian_columns(q, g, c, s, e, v, weights, nr),
                reference,
                rtol=2e-12,
                atol=5e-13,
            )

    def test_source_selection_preserves_each_mode_and_both_readouts(self):
        import lapw_atomic_geometry as geometry

        bands = self.rng.normal(size=4), self.rng.normal(size=4)
        charge = self.rng.normal(size=3)
        # One axis per orbit. Test the independent fixed and induced scalars.
        modes = [(i + 1, np.eye(3)[i], np.diag(np.eye(3)[i])) for i in range(3)]
        tangents = self.rng.normal(size=(3, 4))
        spatial = self.rng.normal(size=(3, 3))
        with patch.object(
            geometry, "native_band_geometry_source", side_effect=[(x, None) for x in tangents]
        ), patch.object(geometry, "native_spatial_geometry_source", side_effect=spatial):
            rows = geometry.atomic_geometry_readouts(None, None, None, modes, *bands, charge)
        np.testing.assert_allclose(
            rows, [tangents @ bands[0], tangents @ bands[1] + spatial @ charge]
        )
        # Multiple coordinates of one orbit. Each Cartesian gradient must
        # project onto the corresponding pattern, including equivalent sites.
        modes.append((1, np.eye(3)[1], np.roll(modes[0][2], 1, axis=1)))
        gradients = self.rng.normal(size=(2, 3, 3))
        with patch.object(geometry, "atomic_geometry_gradients", return_value=gradients):
            rows = geometry.atomic_geometry_readouts(None, None, None, modes, *bands, charge)
        ref = [[np.vdot(g, mode[2]).real for mode in modes] for g in gradients]
        np.testing.assert_allclose(rows, ref)
        self.assertEqual(
            geometry.atomic_geometry_readouts(None, None, None, [], *bands, charge).shape, (2, 0)
        )
