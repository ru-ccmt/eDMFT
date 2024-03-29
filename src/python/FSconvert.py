#!/usr/bin/env python
# @Copyright 2007 Kristjan Haule
import numpy as np
#from scipy import *
import sys, re
from utils import W2kEnvironment, Ry_in_eV
import argparse


usage = """Reads eigvals.dat and case.outputkgen to construct Fermi surface file for xcrysden.  
Also requires EF.dat, case.struct"""

parser = argparse.ArgumentParser(description=usage)
parser.add_argument('fname', nargs='?', default='eigvals.dat', type=str, help='filename for eigenvalues eigvals.dat')
parser.add_argument('-n', type=str, default='0:50', help='range of bands to print [-n0:50]')
args = parser.parse_args()


print(args)
band0,band1 = list(map(int,args.n.split(':')))
feigvals = args.fname

# Finds what is case
w2k = W2kEnvironment()

#band0=0
#band1=50

f = open(w2k.case+'.outputkgen', 'r')

head0=27 #23
head1=8 # 4


for line in f:
    if re.search('G1\s*G2\s*G3',line) is not None: break

bs=[]
for i in range(3):
    bs.append( list(map(float,next(f).split())) )
bs=np.array(bs).transpose()
b0 = bs[0]
b1 = bs[1]
b2 = bs[2]
print(b0, b1, b2)

for line in f:
    if re.search('NO.\s*OF\s*MESH\s*POINTS\s*IN\s*THE\s*BRILLOUIN\s*ZONE', line) is not None: break


nk = list(map(int, next(f).split()[6:10]))
next(f)

print('nk=', nk)

index = np.zeros((nk[0]+1,nk[1]+1,nk[2]+1),dtype=int)
for k0 in range(nk[0]+1):
    for k1 in range(nk[1]+1):
        for k2 in range(nk[2]+1):
            data = next(f).split()
            #print data, k0, k1, k2
            (ii, i0, i1, i2, tindex) = list(map(int, data[:5]))
            if (k0!=i0): print('ERROR0', k0, i0)
            if (k1!=i1): print('ERROR1', k1, i1)
            if (k2!=i2): print('ERROR2', k2, i2)
            index[k0,k1,k2] = tindex
            #print k0, k1, k2, tindex

next(f)
Nkmax = max(index.flatten())

wind={}
for ik in range(Nkmax):
    line = next(f)
    if line[2:5]=='sum': break
    (r0,r1) = list(map(int,line.split()[:2]))
    if r0-1!=ik: print('ERROR3', r0, ik)
    wind[r1]=r0-1
    #print r0, r1

for k0 in range(nk[0]+1):
    for k1 in range(nk[1]+1):
        for k2 in range(nk[2]+1):
            index[k0,k1,k2] = wind[index[k0,k1,k2]]
            #print "%3d %3d %3d %4d" % (k0, k1, k2, index[k0,k1,k2])

Nkmax = max(index.flatten())+1
print(Nkmax)


g = open(feigvals)
n_emin=0
n_emax=1000
for ik in range(Nkmax):
    (ikp, isym, nbands, nemin, nomega) = list(map(int, next(g).split()[1:6]))
    if ikp-1!=ik: print('ERROR4', ikp-1, ik)
    if isym!=1: print('ERROR isym!=1')
    for im in range(nomega): next(g)
    n_emin = max(n_emin, nemin)
    n_emax = min(n_emax, nemin+nbands)

print('n_emin(eigvals)=', n_emin)
print('n_emax(eigvals)=', n_emax)

g.close()
g = open('eigvals.dat')
eigenvals = np.zeros((Nkmax,n_emax-n_emin),dtype=float)
for ik in range(Nkmax):
    (ikp, isym, nbands, nemin, nomega) = list(map(int, next(g).split()[1:6]))
    if ikp-1!=ik: print('ERROR4', ikp-1, ik)
    if isym!=1: print('ERROR isym!=1')
    dat = list(map(float,next(g).split()))
    for im in range(1,nomega): next(g)
    evals = [dat[1+2*i] for i in range(nbands)]
    eigenvals[ik,:] = evals[n_emin-nemin:n_emax-nemin]


ef = open('EF.dat')
EF = float(next(ef))
ef.close()

#print 'eigvals=', eigenvals
print('EF=', EF)

band0 = max(band0,0)
band1 = min(band1,nbands)


h = open(w2k.case+'.bxsf', 'w')

nbands=n_emax-n_emin
band1 = min(band1,n_emax-n_emin-1)
print(' BEGIN_INFO', file=h)
print('   #', file=h)
print('   # Launch as: xcrysden --bxsf '+w2k.case+'.bxsf', file=h)
print('   #', file=h)
print('   Fermi Energy:', '%7.5f' % (EF/Ry_in_eV), file=h)
print(' END_INFO', file=h)
print(file=h)
print(' BEGIN_BLOCK_BANDGRID_3D', file=h)
print('   bands_energies', file=h)
print('   BEGIN_BANDGRID_3D_BANDS', file=h)
print('    ', band1-band0+1, file=h)
print('    ', nk[0]+1, nk[1]+1, nk[2]+1, file=h)
print('  ', '%10.6f '*3 % (0.0, 0.0, 0.0), file=h)
print('  ', '%10.6f '*3 % tuple(b0), file=h)
print('  ', '%10.6f '*3 % tuple(b1), file=h)
print('  ', '%10.6f '*3 % tuple(b2), file=h)

#print 'band0=', band0, 'band1=', band1
print('band0=', band0, 'band1=', band1)

for iband in range(band0,band1+1):
    print('    BAND:', iband+1, file=h)
    for k0 in range(nk[0]+1):
        for k1 in range(nk[1]+1):
            print('    ', end=' ', file=h)
            for k2 in range(nk[2]+1):
                #print shape(eigenvals), index[k0,k1,k2], iband
                de=0.0
                print("%10.6f " % (eigenvals[index[k0,k1,k2],iband]/Ry_in_eV+de), end=' ', file=h)
            print(file=h)
        if k0<nk[0]: print(file=h)
        
print('   END_BANDGRID_3D', file=h)
print(' END_BLOCK_BANDGRID_3D', file=h)
