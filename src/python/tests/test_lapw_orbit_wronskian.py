"""Verify native ordered triangle pooling with complex phases and rectangular band blocks."""
from functools import partial
import unittest
import numpy as np
from lapw_spherical_response import PreparedSphericalBand, prepare_wronskian_band
from lapw_augmentation_response import PreparedAugmentationChannel
from test_lapw_potential_assembly import fixture


class OrbitWronskianTests(unittest.TestCase):
    def test_site_query_forward_and_adjoint(self):
        rng, ch, phi, lo, metric, mesh, energy, lo_energy, _, geometry = fixture()
        q, _, coeff, match, vol, npw = geometry
        positions = rng.normal(size=(3, 3))
        phases = np.exp(1j * positions @ q.T)
        projected = partial(prepare_wronskian_band, ch.ell, ch.rmt, vol, q[:npw],
                            coefficients=coeff[:npw], matching=match, retained=2)
        summed = projected(phases[:, :npw])
        expected = [projected(phase[:npw]) for phase in phases]
        np.testing.assert_allclose(summed, sum(expected), atol=1e-13)
        self.assertGreater(np.linalg.norm(summed[2:]), 1e-2)
        rebuilt = []
        def reader(site):
            rebuilt.append(site)
            return projected(phases[site, :npw])
        original, pooled = [], []
        for i, phase in enumerate(phases):
            channel = PreparedAugmentationChannel(ch.ell, ch.rmt, vol, q, np.linalg.norm(q,axis=1),
                                                  phase, coeff, npw, [[5], [6]], matching_values=match)
            args = channel, phi, lo, metric, mesh, energy, lo_energy
            original.append(PreparedSphericalBand(*args, expected[i], lapw=True))
            pooled.append(PreparedSphericalBand(*args, summed if i == 0 else None,
                                                lapw=True, cache_products=False,
                                                site_triangle=partial(reader,i), retained=2))
        dphi = rng.normal(size=phi.shape)
        dlo = rng.normal(size=lo.shape)
        hbar, obar = original[0].value()
        for method, args in [('value', ()), ('action', (dphi,dlo)), ('pullback', (hbar,obar))]:
            full = [getattr(op, method)(*args) for op in original]
            reduced = [getattr(op, method)(*args,orbit_sum=True) for op in pooled]
            for a,b in zip(zip(*full),zip(*reduced)):
                np.testing.assert_allclose(sum(a),sum(b),rtol=5e-13,atol=2e-10)
        self.assertEqual(rebuilt, [])
        for i,(a,b) in enumerate(zip(original,pooled)):
            np.testing.assert_allclose(b.triangle,expected[i],atol=1e-13)
            for left,right in zip(a.value(),b.value()):
                np.testing.assert_allclose(left,right,atol=1e-12)
        self.assertEqual(rebuilt,[0,0,1,1,2,2])
