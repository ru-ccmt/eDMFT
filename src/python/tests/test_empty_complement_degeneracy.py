"""An empty-empty mixing column cannot affect charge or projector readouts."""
import unittest
import numpy as np
from lapw_complement_response import ShiftedComplement

class EmptyDegeneracyTests(unittest.TestCase):
    def test_occupied_columns_remain_exact_at_empty_degeneracy(self):
        h=np.diag([3.,4.,7.]);eps=np.array([-1.,0.,3.]);mask=np.array([True,True,False])
        with self.assertRaises(ValueError):ShiftedComplement(h,eps,mu=.1,temperature=.01)
        op=ShiftedComplement(h,eps,mu=.1,temperature=.01,active_columns=mask)
        rhs=np.random.default_rng(22).normal(size=(3,3));y=rhs.T.copy()
        result=op.action(rhs)
        np.testing.assert_allclose(result[:,:2],rhs[:,:2]/(eps[None,:2]-np.diag(h)[:,None]),atol=1e-16)
        np.testing.assert_array_equal(result[:,2],0)
        self.assertLess(op.residual(result,rhs),1e-16)
        self.assertAlmostEqual(np.vdot(y,result).real,np.vdot(op.action(y),rhs).real,places=14)
    def test_nonempty_column_cannot_be_deactivated(self):
        with self.assertRaisesRegex(ValueError,'not thermally empty'):
            ShiftedComplement(np.diag([3.,4.]),[-1.,0.],mu=.1,temperature=.01,active_columns=np.array([True,False]))

if __name__=='__main__':unittest.main()
