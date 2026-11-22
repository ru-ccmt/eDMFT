"""Regression checks for changing the lattice mesh with a frozen impurity."""
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch
import numpy as np

from matsubara_mesh import auto_mesh, cached_auto_mesh, log_indices


class MatsubaraMeshTest(unittest.TestCase):
    def test_fit_rows_are_invariant_and_mc_points_are_preserved(self):
        omega = (2*np.arange(4774) + 1)*np.pi/100
        dynamic = 4/(1j*omega - 1.5)
        expected_fit = log_indices(len(omega), 200, 300)[-24:]
        grids = []
        for tolerance in (1.e-4, 1.e-7):
            rows, fit = auto_mesh(omega, [dynamic, dynamic], .025, 200, tolerance)
            np.testing.assert_array_equal(rows[:200], np.arange(200))
            np.testing.assert_array_equal(rows[fit-1], expected_fit)
            self.assertEqual(rows[-1],len(omega)-1)
            self.assertTrue(np.all(np.diff(rows)>0))
            grids.append(rows)
        self.assertGreater(len(grids[1]), len(grids[0]))

    def test_short_mesh(self):
        for n in (1, 4, 200, 201, 203):
            omega = 2*np.arange(n)+1.
            rows, fit = auto_mesh(omega, [], .02, 200, 1.e-5)
            np.testing.assert_array_equal(rows, np.arange(n))
            self.assertEqual(fit[-1], n)

    def test_invalid_mesh_and_tolerance(self):
        for omega,tol in (([1.,5.],1.e-6),([1.,3.],0.),([1.,3.],np.nan)):
            with self.assertRaises(ValueError):
                auto_mesh(np.array(omega), [], .025, 2, tol)

    def test_cache_tracks_input_and_recovers_from_damage(self):
        omega = 2*np.arange(400)+1.
        sigma = 2/(1j*omega-1)
        with TemporaryDirectory() as work:
            cache = Path(work)/'mesh.npz'
            def select(component):
                return cached_auto_mesh(omega,[component],.02,30,1.e-6,cache)
            rows, fit = select(sigma)
            with patch('matsubara_mesh.auto_mesh', wraps=auto_mesh) as compute:
                np.testing.assert_array_equal(select(sigma.copy())[0], rows)
                compute.assert_not_called()
                select(1.2*sigma)
                compute.assert_called_once()
            cache.write_bytes(b'PK\x03\x04damaged')
            rebuilt = select(sigma)
            np.testing.assert_array_equal(rebuilt[0], rows)
            np.testing.assert_array_equal(rebuilt[1], fit)
