"""Immutable scratch data have the same matching derivatives as resident data."""
from contextlib import ExitStack
from copy import copy
import unittest
import numpy as np
from lapw_response_scratch import MatchingScratch
from lapw_spherical_response import PreparedSphericalBand
from test_lapw_potential_assembly import fixture


class ScratchTests(unittest.TestCase):
    def test_site_values_and_transpose(self):
        rng,ch,phi,lo,metric,r,e,el,triangle,_=fixture()
        args=phi,lo,metric,r,e,el,triangle[:,:2]
        op=PreparedSphericalBand(ch,*args,lapw=True,cache_products=False)
        hbar,obar=op.value()
        ref=op.pullback(hbar,obar)
        with ExitStack() as stack:
            scratch=MatchingScratch(stack)
            empty = np.empty((4, 0, 2), dtype=complex)
            np.testing.assert_array_equal(scratch.store(empty).read(), empty)
            stored=PreparedSphericalBand(copy(ch),*args,lapw=True,cache_products=False)
            stored.store_matching(scratch)
            self.assertNotIn("fj",stored.channel.__dict__)
            self.assertFalse(stored.channel.fj.flags.writeable)
            for x,y in zip(ref,stored.pullback(hbar,obar)):
                np.testing.assert_allclose(x,y,rtol=1e-14,atol=1e-13)
            # The moving-basis derivative assigns new F,D to a shallow copy.
            moving=copy(stored.channel)
            moving.fj=np.zeros_like(moving.fj)
            np.testing.assert_array_equal(stored.channel.fj,ch.fj)
        self.assertTrue(scratch.file.closed)
