"""The native complement may be transformed once into the tridiagonal solver gauge."""
import unittest
import numpy as np
from scipy.linalg import eigh
from lapw_complement_response import LAPWComplementBasis


class RotatedComplementTests(unittest.TestCase):
    def test_physical_response_is_invariant(self):
        rng=np.random.default_rng(490)
        z=rng.normal(size=(9,9))+1j*rng.normal(size=(9,9))
        o=np.eye(9)+z@z.conj().T/9
        h=z+z.conj().T
        ev,c=eigh(h,o)
        args=h,o,c[:,:3],ev[:3]
        bases=[LAPWComplementBasis.prepare(*args,mu=ev[1],temperature=1e-4,tridiagonal_basis=x) for x in (False,True)]
        native=[]
        for base in bases:
            _, rhs=base.project_matrix(z@z.conj().T)
            t=base.solver.action(rhs)
            self.assertLess(base.solver.residual(t,rhs),5e-14)
            native.append(base.complement@t)
        np.testing.assert_allclose(*native,atol=2e-13)
        self.assertLess(bases[1].solver.info()['stored_bytes'],bases[0].solver.info()['stored_bytes'])
