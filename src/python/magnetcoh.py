#!/usr/bin/env python
# @Copyright 2007 Kristjan Haule
"""
This code performs the unfolding of the DFT electronic structure by reading eigenvectors from vector file,
creating overlap matrix from eigenvectors (because the eigenvectors satisfy A*O*A^H=1) and than
adding appropriate phase factors to eigenvectors to get unfolded electronic structure. 
The DFT eigenvectors are:
      |psi_{i,k}> = \sum_G A_k[i,G]|chi_{G}>
   where the LAPW overlap comes from non-orthogonal basis functions |chi_{G}>
       O_{GG'} = <chi_{G}|chi_{G'}>
   which immediately requires A*O*A^H = 1 (note O is symmetric).

The unfolded LAPW basis vector has the form 
    |chi_{k+G}> -> \sum_{dR}|chi_G>e^{i*{k+G}*dR} 
where dR are the vectors connecting (approximate) equivalent unit cells, which require unfolding. 
For two equivalent unit cells, we have dR=( [0,0,0], [deltax,deltay,deltaz] )

The eigenvalues are than:

   epsi_{i,k} = <psi_{i,k}|H|psi_{i,k}> = 
              = \sum_{dR1,dR2} e^{i*(G'+k)*dR1-(G+k)*dR2} A*_k[i,G]<chi_{G}|H|chi_{G'}>A_k[i,G']
              = \sum_{dR1} e^{i*(G'+k)*dR1}*A_k[i,G'] O_{G'G} \sum_{dR2} A^H[G,i]e^{-i*(G+k)*dR2}

For two atoms per unit cell, we have:
   A0 == A_k[i,G]
   A1 == A_k[i,G]*e^{(G+k)*dR}
   and 
   epsi_{i,k} = (A0+A1)*O*(A0^H + A1^H) = I + A1*O*(A0^H + A1^H) + A0*O*A1^H
which is implemented below.

Note that because we do not keep all bands in A, we need to compute overlap with SVD not to encounter
singular matrix.

"""
import sys, re, os
from scipy import *
from pylab import *
import rdU, utils, findEF
from wstruct import Struct
########################################
# input that we need to have
delta = array([0.5, -0.5, 0.0])

# input for LDA plot if needed
plotLDA = True
gamma = 0.05
omega = linspace(-1,1,100)

#########################################
mingle_names={'W':'$W$','L':'$L$','LAMBDA':'$\Lambda$','GAMMA':'$\Gamma$','DELTA':'$\Delta$','X':'$X$','Z':'$Z$','W':'$W$','K':'$K$'}

# Vector-file handle
vtape=9

w2k = utils.W2kEnvironment()
print('case=', w2k.case)
EF = findEF.FindChemicalPotential(w2k.case, '')[0]
struct = Struct()
struct.ReadStruct(w2k.case+'.struct')

vectortype=float
vecread = rdU.fvread3
so=''
SCRATCH = '.' if w2k.SCRATCH==None else w2k.SCRATCH
if os.path.isfile(SCRATCH+"/"+w2k.case+".vectorso") and os.path.getsize(SCRATCH+"/"+w2k.case+".vectorso")>0 :
    print('Found '+w2k.case+'.vectorso file, hence assuming so-coupling exists. Switching -so switch!')
    so = 'so'
    vectortype=complex
    vecread = rdU.fvread3c

fc = open('cohfactors.dat', 'w')
print('# Coherence factors for KS states', file=fc)

maxkpoints = 10000
# opens vector file
rdU.fvopen(vtape, w2k.case+'.vector'+so, struct.nat) 
gw=[]
knames={}
for ik in range(maxkpoints):
    # Reads vector file
    (k, kname, wgh, ios, n0, nb) = rdU.fvread1(vtape)
    if ios!=0: break # vector file is finished, no more k-points
    if kname != b'':
        kname = kname.decode('UTF-8')
        print('k=', k, kname)
        if kname in mingle_names:
            kname = mingle_names[kname]
        knames[ik] = kname
    # Reciprocal vectors
    Gs = rdU.fvread2(vtape, n0)

    # Reading KS eigensystem
    As=zeros((nb,n0), dtype=vectortype)
    Ek=zeros(nb, dtype=float)
    for i in range(nb):
        (num, ek, A) = vecread(vtape, n0)
        As[i,:] = A             # KS eigenvector
        Ek[i] = utils.Ry2eV(ek) # KS eigenvalue

    # Building overlap from eigenvectors only.
    # A = u * s * v                        # with SVD
    # A * O * A^H = 1                      # has to be satisified
    # O = v^H * 1/s * u^H * u * 1/s * v    # is than also satisified
    (u,s,v) = linalg.svd( array(As) )
    uu = u.conj().T @ u
    for i in range(len(uu)):
        for j in range(len(uu)):
            uu[i,j] = 1/s[i] * uu[i,j] * 1/s[j]
    N = len(uu)
    vs = v[:N,:] # need to truncate those components that have singular values
    Op = vs.conj().T @ uu @ vs
    
    # Building phase factors for the two sublatices
    ph_iG = delta @ Gs + delta @ k         # phi_iG[iG] = delta[:3]*Gs[:3,iG]
    exph = exp(-1j*2*pi*ph_iG)  # exph[iG] is large vector
        
    # We want to use matrix multiplications only. Need the eigenvector with the phase factor As1.
    As0 = array(As, dtype=complex)
    As1 = zeros(shape(As), dtype=complex) #
    As1 = As * exph                       #  As1[:,ig] = As[:,ig]*exph[ig]

    I = identity(nb)
    gs = I + As1 @ Op @ (As1.conj().T + As0.conj().T) + As0 @ Op @ As1.conj().T

    # Printing the coherence factors
    print(ik+1, 1, len(gs), file=fc)
    for i in range(len(gs)):
        for j in range(len(gs)):
            print("%15.10f " % gs[j,i].real, end=' ', file=fc)
        print(file=fc)
    
    if plotLDA:
        gt = zeros(len(omega), dtype=float)
        for iom,om in enumerate(omega):
            dsum = 0.0
            for i in range(nb):
                dsum -= gs[i,i]*(1./(om+gamma*1j+EF-Ek[i])).imag
                #dsum -= (1./(om+gamma*1j+EF-Ek[i])).imag
            gt[iom] = dsum.real
        gw.append( gt )

nkp=ik
print('nkp=', nkp)
print('knames=', knames)
if plotLDA:
    nkp=ik
    gw = transpose(array(gw))
    
    xmm = [0, nkp]
    ymm = [omega[0], omega[-1]]
    
    imshow(gw, interpolation='bilinear', cmap=cm.hot, origin='lower', extent=[xmm[0],xmm[1],ymm[0],ymm[1]], aspect=(xmm[1]-xmm[0])*0.8/(ymm[1]-ymm[0]) )
    x_min, x_max = xlim()
    y_min, y_max = ylim()
    
    for wi,name in knames.items():
        cs=plot([wi,wi], [y_min,y_max], 'w-')

    xticks( [ k for k in knames], [knames[k] for k in knames], fontsize='x-large' )    
    plot([0,nkp],[0,0], 'w:')
    show()
