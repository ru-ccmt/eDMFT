"""Complete tensor bases and nonlinear oblique-cell shear endpoints."""
from copy import deepcopy
from pathlib import Path
from unittest.mock import patch
import tempfile, unittest
import numpy as np
from lapw_total_stress import invariant_strains, stress_from_contractions
from lapw_lattice import NativeLattice
from lapw_strain_field import symmetric_strain_potential
from lapw_hydrostatic_potential import hydrostatic_explicit_potential
from lapw_potential_core_state import save_potential_core_state
from lapw_xc_response import InterstitialXCResponse
from test_lapw_potential_core_state import inputs


class GeneralShearTests(unittest.TestCase):
    def test_triclinic_all_six_and_monoclinic_allowed_shear(self):
        s=inputs()['structure'];s.alpha=81.;s.beta=98.;s.gamma=103.
        s.timat=np.array([np.eye(3,dtype=int),-np.eye(3,dtype=int)])
        e=invariant_strains(s)
        self.assertEqual(len(e),6)
        sigma=np.array([[.1,.21,-.13],[.21,-.2,.32],[-.13,.32,.7]])
        volume=NativeLattice.from_structure(s).volume
        d=volume*np.einsum('ij,aij->a',sigma,e)
        np.testing.assert_allclose(stress_from_contractions(e,d,volume),sigma,atol=1e-14)
        s.alpha=s.gamma=90.;s.beta=98.
        s.timat=np.array([np.eye(3,dtype=int),np.diag([-1,1,-1])])
        e=invariant_strains(s)
        self.assertEqual(len(e),4)
        self.assertTrue(any(abs(x[0,2])>.1 for x in e))
        for x in e:
            for rotation in NativeLattice.from_structure(s).cartesian_rotations(s.timat):
                np.testing.assert_allclose(rotation@x@rotation.T,x,atol=1e-13)

    def test_shear_potential_against_nonlinear_cell_reconstruction(self):
        args=inputs();s=args['structure'];s.alpha=81.;s.beta=98.;s.gamma=103.
        lattice=NativeLattice.from_structure(s)
        d=args['density']
        for g,z in [((1,0,0),.002+.001j),((-1,0,0),.002-.001j),
                    ((1,1,1),.0003+.0002j),((-1,-1,-1),.0003-.0002j)]:d.pw[g]=z
        # Keep the Cartesian orientation exactly F A. Reconstructing WIEN
        # cell lengths/angles would add a rigid rotation to this independent
        # finite-difference oracle; that is a separate coordinate conversion.
        e=np.array([[.2,.3,-.2],[.3,-.1,.4],[-.2,.4,-.1]])
        with tempfile.TemporaryDirectory() as directory:
            def endpoint(t,functional):
                cell=lattice.strained(e,t)
                with patch.object(NativeLattice,'from_structure',return_value=cell):
                    pc=save_potential_core_state(Path(directory)/'state.npz',**deepcopy(args))
                    pc.xc.ir=InterstitialXCResponse(pc.hartree.reanalysis,
                        np.array(list(d.pw.values())),cell.cell,functional=functional)
                    value=hydrostatic_explicit_potential(pc)[1]['central_potential']
                    return pc,value
            for functional in ['LDA','PBE']:
                pc,_=endpoint(0,functional)
                actual,_=symmetric_strain_potential(pc,e,transport_rate=0)
                errors=[]
                for h in [1e-4,5e-5]:
                    minus=endpoint(-h,functional)[1];plus=endpoint(h,functional)[1]
                    errors.append(np.linalg.norm(actual-(plus-minus)/(2*h)))
                self.assertLess(errors[-1],2e-7)
                np.testing.assert_allclose(actual,(plus-minus)/(2*h),atol=4e-8,rtol=3e-6)
