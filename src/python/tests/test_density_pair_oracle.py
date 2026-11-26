"""Independent finite-grid FFT oracle for direct density and its large-basis fallback."""
import unittest
import numpy as np
from lapw_density_response import FourierStars, InterstitialDensity

class PlaneWaveDensityTests(unittest.TestCase):
    def test_mixed_density_under_screw_symmetry(self):
        rng = np.random.default_rng(1496)
        g = np.array([[-1,0,0],[0,-1,1],[1,1,1],[0,-1,-1],[0,0,0]])
        c, weighted = (rng.normal(size=(len(g),6)) + 1j*rng.normal(size=(len(g),6)) for _ in range(2))
        shape = (6,6,7)
        reps = np.array([[0,0,0],[1,0,1],[-1,0,-1],[3,0,1],[-3,0,-1],[2,1,0]])
        rot = np.array([[0,-1,0],[1,0,0],[0,0,1]])
        rots = np.array([np.linalg.matrix_power(rot,k) for k in range(4)])
        tau = np.column_stack((np.zeros((4,2)), np.arange(4)/4))
        def fft(x):
            grid = np.zeros((*shape,x.shape[1]),complex)
            grid[tuple(g.T)] = x
            return np.fft.ifftn(grid,axes=(0,1,2))*np.prod(shape)
        fourier = np.fft.fftn((fft(weighted)*fft(c).conj()).real.sum(axis=-1))/(np.prod(shape)*91.0)
        both = []
        for complex_phases in (True,False):
            stars = FourierStars(shape,reps,rots,tau,complex_phases=complex_phases)
            reference = stars.action(fourier)
            both.append(reference)
            op = InterstitialDensity(c,g,91.0,stars,cache_wavefunctions=False)
            for direct in (True,False):
                if not direct:
                    op.direct_labels = None
                np.testing.assert_allclose(op.density_from_weighted_coefficients(weighted),reference,atol=3e-15)
        self.assertGreater(np.linalg.norm(both[0]-both[1]),1e-4)
