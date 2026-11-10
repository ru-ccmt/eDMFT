"""Independent metric, centering, and deformation checks for WIEN cells."""
import unittest
from types import SimpleNamespace
import numpy as np
from lapw_lattice import NativeLattice
from lapw_cubic_harmonics import cubic_maps, boundary_factor


def structure(kind, angles=(90.,90.,90.)):
    return SimpleNamespace(lattyp=kind,a=6.,b=6.,c=9.,
                           alpha=angles[0],beta=angles[1],gamma=angles[2])


class LatticeTests(unittest.TestCase):
    def test_volume_and_dual_coordinate_conventions(self):
        for kind, ratio, angles in [('P',1,(90,90,90)),('S',1,(90,90,90)),
          ('F',.25,(90,90,90)),('B',.5,(90,90,90)),
          ('CXY',.5,(90,90,90)),('CXZ',.5,(90,90,90)),('CYZ',.5,(90,90,90)),
          ('H',np.sqrt(3)/2,(90,90,120)),('R',np.sqrt(3)/6,(90,90,120))]:
            with self.subTest(kind=kind):
                g=NativeLattice.from_structure(structure(kind,angles))
                self.assertAlmostEqual(g.volume,6*6*9*ratio,places=11)
                np.testing.assert_allclose(g.reciprocal.T@g.cell,2*np.pi*np.eye(3),atol=2e-15)
                np.testing.assert_allclose(g.primitive_reciprocal.T@g.primitive_cell,2*np.pi*np.eye(3),atol=2e-15)
        s=structure('P',(76,83,97)); g=NativeLattice.from_structure(s)
        cosa,cosb,cosg=np.cos(np.radians([s.alpha,s.beta,s.gamma]))
        expected=s.a*s.b*s.c*np.sqrt(1+2*cosa*cosb*cosg-cosa*cosa-cosb*cosb-cosg*cosg)
        self.assertAlmostEqual(g.volume,expected,places=11)

    def test_centering_and_allowed_reflections(self):
        f=NativeLattice.from_structure(structure('F'))
        np.testing.assert_array_equal(f.centering_factor([[0,0,0],[1,0,0],[1,1,0],[1,1,1],[2,0,0]]),[1,0,0,1,1])
        self.assertTrue(f.equivalent_displacement([.5,.5,0]))
        self.assertFalse(f.equivalent_displacement([.5,0,0]))
        b=NativeLattice.from_structure(structure('B'))
        np.testing.assert_array_equal(b.centering_factor([[1,0,0],[1,1,0],[1,1,1],[2,0,0]]),[0,1,0,1])
        self.assertTrue(b.equivalent_displacement([.5,.5,.5]))

    def test_hexagonal_rotation_and_metric(self):
        g=NativeLattice.from_structure(structure('H',(90,90,120)))
        rotation=np.array([[1,-1,0],[1,0,0],[0,0,1]])
        r=g.cartesian_rotations(rotation)
        np.testing.assert_allclose(r.T@r,np.eye(3),atol=4e-16)
        # a,b are 120 degrees, whereas reciprocal directions are 60 degrees.
        self.assertAlmostEqual(g.cell[:,0]@g.cell[:,1],-18.)
        self.assertAlmostEqual(g.momenta([[1,0,0]])[0]@g.momenta([[0,1,0]])[0],(2*np.pi/6)**2*2/3)

    def test_general_deformation_derivatives(self):
        e=np.array([[.3,.2,-.1],[.2,-.15,.08],[-.1,.08,.4]])
        for kind,angles in [('F',(90,90,90)),('B',(90,90,90)),('H',(90,90,120)),('P',(76,83,97))]:
            g=NativeLattice.from_structure(structure(kind,angles)); h=1e-5
            p,m=g.strained(e,h),g.strained(e,-h)
            np.testing.assert_allclose((p.reciprocal-m.reciprocal)/(2*h),-e.T@g.reciprocal,rtol=2e-9,atol=3e-11)
            self.assertAlmostEqual((p.volume-m.volume)/(2*h)/g.volume,np.trace(e),places=10)

    def test_cubic_coordinate_maps(self):
        c,d=cubic_maps([(0,0),(4,0),(4,4),(6,0),(6,4)])
        np.testing.assert_allclose(c@c.T,np.eye(3),atol=3e-16)
        np.testing.assert_allclose(c@d,np.eye(3),atol=1e-16)
        self.assertAlmostEqual(boundary_factor(4,0),6/7)
        self.assertAlmostEqual(boundary_factor(4,4),6/5)
        with self.assertRaises(ValueError): cubic_maps([(0,0),(4,0)])

class NativeFrameTests(unittest.TestCase):
    def test_identity_is_swapped_first_before_site_mapping(self):
        from lapw_lattice import native_site_frames
        s=structure('F');s.nat=1;s.mult=[1];s.pos=np.array([[0.,0.,0.]])
        s.timat=np.array([np.diag([1,-1,-1]),np.eye(3,dtype=int)])
        s.tau=np.zeros((2,3))
        axes,locations=native_site_frames(s)
        np.testing.assert_allclose(axes[0],np.eye(3),atol=2e-16)
        np.testing.assert_array_equal(locations,s.pos)

    def test_native_rounding_repair_preserves_input(self):
        from lapw_lattice import native_structure
        s=structure('H',(90,90,120));s.pos=np.array([[.33333333,.66666667,.25]])
        s.timat=np.array([np.eye(3,dtype=int),[[-1,1,0],[-1,0,0],[0,0,1]],[[0,-1,0],[1,-1,0],[0,0,1]]])
        s.tau=np.zeros((3,3))
        corrected=native_structure(s)
        np.testing.assert_allclose(corrected.pos,[[1/3,2/3,.25]],atol=2e-16)
        np.testing.assert_array_equal(s.pos,[[.33333333,.66666667,.25]])
        np.testing.assert_allclose(native_structure(corrected).pos,corrected.pos,atol=2e-16)

if __name__=='__main__': unittest.main()
