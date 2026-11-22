"""Check the site-resolved ordered surface kernel before band projection."""
from types import SimpleNamespace
import unittest
import numpy as np
from scipy.special import spherical_jn
from lapw_geometry_response import atomic_wronskian_columns
from lapw_spherical_response import prepare_wronskian_band


class AtomicSurfaceTests(unittest.TestCase):
    def test_orbit_derivative_matches_displaced_native_upper_copy(self):
        rng = np.random.default_rng(123)
        n=270
        labels = np.array(np.unravel_index(rng.choice(1000,n,replace=False),(10,10,10))).T - 5
        labels[0] = 0
        reciprocal = np.array([[1.2,.15,0.],[0.,.85,.07],[.03,0.,.7]])
        q = labels @ reciprocal.T
        c = (rng.normal(size=(n,8))+1j*rng.normal(size=(n,8)))/np.sqrt(n)
        pos = np.array([[.16,.15,.21],[.4,.6,.35],[.34,.1,.29]])
        u = np.array([[.5,0,.25],[-.5,.15,-.25],[-.1,.2,.33]])
        structure = SimpleNamespace(nat=2, mult=[2,1],rmt=[1.41,1.63],pos=pos)
        volume=312.
        weights=[{0:.31,1:-.7,4:.15},{1:.4,2:-.2}]
        nr=3
        actual=atomic_wronskian_columns(q,labels,c,structure,u,volume,weights,nr)
        def displaced(t):
            result=np.zeros_like(actual)
            for i, a in enumerate([0,0,1]):
                phase = np.exp(2j*np.pi*(labels @ pos[i])+1j*t*(q @ u[i]))
                length=np.linalg.norm(q,axis=1)
                x=length*structure.rmt[a]
                for ell,w in weights[a].items():
                    matching=np.array([spherical_jn(ell,x),length*spherical_jn(ell,x,derivative=True)])
                    result += w*prepare_wronskian_band(ell,structure.rmt[a],volume,q,phase,c,matching,retained=nr)
            return result
        h=1e-6
        np.testing.assert_allclose(actual,(displaced(h)-displaced(-h))/(2*h),atol=4e-10,rtol=1e-8)
