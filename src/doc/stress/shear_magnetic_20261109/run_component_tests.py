"""Reproduce the staged component suite, including native reference controls."""
from pathlib import Path
import os
import sys
import unittest
root = Path(__file__).resolve().parent
os.environ.setdefault('EDMFT_RADIAL_LIBRARY', '/Users/haule/dbin/liblapw_response.so')
os.environ.setdefault('EDMFT_HARTREE_REFERENCE_LIBRARY', '/private/tmp/edmft_hartree_radial_20261108/native/libhartree_reference.dylib')
os.environ.setdefault('EDMFT_HARTREE_CHARGE_LIBRARY', '/private/tmp/edmft_hartree_periodic_20261108/native/libcharge_reference.dylib')
os.environ.setdefault('EDMFT_MIXER_NORMALIZATION_LIBRARY', '/private/tmp/edmft_charge_closure_20261108/native/libmixer_normalization.dylib')
sys.path[:0] = [str(root/'stage'), str(root/'stage/tests')]
suite = unittest.defaultTestLoader.discover(str(root/'stage/tests'))
result = unittest.TextTestRunner(verbosity=1).run(suite)
raise SystemExit(not result.wasSuccessful())
