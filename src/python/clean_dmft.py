#!/usr/bin/env python
# @Copyright 2007 Kristjan Haule  
import sys
import re, os, glob, shutil
from os.path import getsize

toclean = glob.glob('*.vector*')
toclean += glob.glob('dmft1.error*')+glob.glob('dmft2.error*')+glob.glob('*.outputdmf1.*')+glob.glob('*.outputdmf2.*')+glob.glob('*.broyd*')+glob.glob('dmft')+glob.glob('dmft2')
toclean += glob.glob('dmft1_error')+glob.glob('*.cdos3')+glob.glob('imp.*/status.*')+glob.glob('imp.*/Probability.dat.*')+glob.glob('imp.*/g_hb*.dat')+glob.glob('imp.*/g_qmc.dat')
toclean += glob.glob('imp.*/check.dat.*')+glob.glob('imp.*/Delta.tau.*')+glob.glob('imp.*/ctqmc')+glob.glob('imp.*/s_hb1.dat')+glob.glob('imp.*/new.cix')+glob.glob('imp.*/Sig.out.*')
toclean += glob.glob('dmft1')+glob.glob('lapw1')

for f in toclean: os.remove(f)

for f in glob.glob('*'):
    if getsize(f)==0: os.remove(f)

