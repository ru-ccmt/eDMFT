#!/usr/bin/env python
""" Prepares coherence factors for 2D Fermi surface plot. 
These are orbital characters on bands, but here we store the entire matrix in band space, 
so that we can use it with DMFT. In DFT we would just need a number for each band, while
for DMFT we need matrix in band space.

It extracts information from the DMFT projector, written in the Udmft.0 file, 
and provides cohfactors.dat file with coherence factors. This information is needed for
fat_fsplot.py, for example.

The input is obtained from DMFT projector Udmft.0, which can be obtained by 
  % x_dmft.py dmftu.

The output are coherence factors g_{ij}, where i,j are bands and the relation between 
g and dmft-prjector U is as follows:

  g_{ij} = (U^H @ U)_{ij}

but here the product @ goes only over very limited set of orbitals, which are specified
through the input dictionary. The g_{ij} are thus coherence factors that relate the 
content of each orbital to the bands. 

The input dictionary, which determines the product @ in U^H @ U has to be in the form:
  {icix1:iorb1,icix2:iorb2,....}.
For example, consider two atoms unit cell with cubic harmonics. We want xz orbital from both 
atoms. The first atom has icix=0 and the second icix=1. The xz orbital has index 2, hence the
dictionary would be
  {0:2, 1:2}

"""
# @Copyright 2007 Kristjan Haule
from numpy import *
from numpy import linalg
import struct, sys, re
import rdU, utils
import findEF

if len(sys.argv)<=1:
    print('Give dictionary of coh_orb in the form {icix1:iorb1,icix2:iorb2,....}')
    sys.exit(0)

#coh_orb = {icix=0:xz,icix=1:xz} -> {0:2,1:2}
coh_orb = eval(sys.argv[1])
print('coh_orb=', coh_orb)




fhp = 233
cpuID = 0

dmfe = utils.DmftEnvironment()  # DMFT paths
w2k = utils.W2kEnvironment()    # W2k filenames and paths
case = w2k.case
(EF,NOE) = findEF.FindChemicalPotential(case,'')

# Read Udmft
filename = 'Udmft.'+str(cpuID)
(nkp, nsymop, norbitals) = rdU.fopen(fhp, filename)
nindo = rdU.read1(fhp, norbitals)

print('nsymop=', nsymop)

fc = open('cohfactors.dat', 'w')

print('# Coherence factors for KS states', file=fc)

for ikp in range(nkp):
    (iikp, nbands, tmaxdim2, tnorbitals, nemin) = rdU.read2(fhp)
    print('ikp=', ikp,  'nbands=', nbands, 'maxdim2=', tmaxdim2, 'norbitals=', tnorbitals, 'nemin=', nemin)
    for isym in range(nsymop):
        iisym = rdU.read3(fhp)
        if iisym!=isym+1:
            print('ERROR: isym!=iisym', isym, iisym)
        
        gs = zeros((nbands,nbands),dtype=complex)
        print('norbitals=', norbitals, 'nindo=', nindo)
        for iorb in range(norbitals):
            # Reading DMFT transformation UDMFT
            U = array([rdU.read4(fhp, nbands) for ind in range(nindo[iorb])])
            # print ikp, isym, iorb, shape(U) 
            
            if (isym==0):
                # Prepares coherence factors gs(i,j)=U(orb,i)*U(orb,j).conj()
                if iorb in coh_orb:
                    xorb = coh_orb[iorb]
                    for iband in range(nbands):
                        gs[iband,:] += U[xorb,iband]*U[xorb,:].conj()
                        
        if isym==0:
            # Printing the coherence factors
            print(iikp, nemin, nbands, file=fc)
            for i in range(len(gs)):
                for j in range(len(gs)):
                    print("%15.10f " % abs(gs[j,i]), end=' ', file=fc)
                print(file=fc)

rdU.fclose(fhp)
            
