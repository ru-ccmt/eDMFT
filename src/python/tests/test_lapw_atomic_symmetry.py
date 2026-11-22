"""Test internal-coordinate multiplicities and virtual work in a monoclinic cell."""
import unittest
from types import SimpleNamespace

import numpy as np

from force_response_modes import internal_modes, resolve_displaced_modes
from lapw_lattice import NativeLattice


class MonoclinicForceModes(unittest.TestCase):
    def test_rounded_oblique_directions_and_extrapolation(self):
        # Large cancelling partials expose a transverse component too small
        # to see in the coordinates printed in a typical force summary.
        direction = [[1.,0.,0.],[0.,1.,0.]]
        explicit = np.array([73.,-2100.])
        induced = np.array([-70.,2093.])
        exact = np.column_stack([explicit, induced, explicit+induced])
        curvature = np.array([[1.,2.,3.],[80.,-83.,-3.]])
        records = [dict(atom=2,direction_cartesian=d,steps=[]) for d in direction]
        steps = (.01,.005)
        for h, coeff in zip(steps, [4e-5,-2e-5]):
            b = np.array([[1.+1e-6,coeff],[0.,1.-2e-6]])
            obs = b @ (exact + h*h*curvature)
            for row, basis, values in zip(records, b, obs):
                row["steps"].append(dict(sampled_mode_coefficients=basis,
                    explicit_force_mRy_per_Bohr=values[0],
                    screened_force_mRy_per_Bohr=values[1],
                    full_force_mRy_per_Bohr=values[2]))
        force = resolve_displaced_modes(records, steps, 2)
        np.testing.assert_allclose(force, [[0,0,0],[3,-7,0]], atol=2e-12)
        np.testing.assert_allclose([r["explicit_force_mRy_per_Bohr"] for r in records], explicit)
        # For a single step, remove the same directional contamination
        # while retaining the central difference's own truncation error.
        single = resolve_displaced_modes(records, steps[:1], 2)
        np.testing.assert_allclose(single[1,:2],exact[:,2] + steps[0]**2*curvature[:,2])

    def test_screw_related_sites_and_inversion_center(self):
        # Unique b, P2_1/c, with both a four-site general orbit and two sites
        # fixed by inversion. Fractional displacement components are not
        # Cartesian components in this metric.
        identity = np.eye(3, dtype=int)
        twofold = np.diag([-1, 1, -1])
        rotations = np.stack([identity, twofold, -identity, -twofold])
        glide = np.array([0., 0.5, 0.5])
        shifts = np.stack([0*glide, glide, 0*glide, glide])
        general = (rotations @ [.21, .13, .32] + shifts) % 1.
        positions = np.vstack([general, np.zeros(3), glide])
        structure = SimpleNamespace(lattyp="P", a=5., b=7., c=9., alpha=90.,
                                    beta=104., gamma=90., mult=[4, 2], nat=2,
                                    timat=rotations, tau=shifts, pos=positions)
        lattice = NativeLattice.from_structure(structure)
        modes = internal_modes(structure)
        self.assertEqual([atom for atom, _, _ in modes], [1, 1, 1])
        self.assertEqual(internal_modes(structure, [2]), [])

        cartesian = lattice.cartesian_rotations(rotations)
        representative_force = np.array([1.3, -4.1, 2.7])
        forces = cartesian @ representative_force
        reconstructed = np.zeros(3)
        for _, unit, pattern in modes:
            np.testing.assert_allclose(pattern[4:], 0., atol=1e-14)
            # Moving the four symmetry images gives four times the work on
            # the representative. The per-atom force must divide by four.
            work = np.sum(forces * pattern[:4]) / 4.
            reconstructed += unit * work
            fractional = np.linalg.solve(lattice.cell, pattern[:4].T).T
            displaced = general + .001 * fractional
            for j in range(4):
                image = rotations[j] @ displaced[0] + shifts[j]
                self.assertTrue(lattice.equivalent_displacement(image - displaced[j]))
        np.testing.assert_allclose(reconstructed, representative_force, atol=2e-14)
