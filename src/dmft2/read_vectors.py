import sys, re, os
from numpy import *
#from pylab import *
#from scipy.lib.blas import fblas
import rdVec, utils #, struct1
from wstruct import Struct

if len(sys.argv)<2:
    print('Please give the name of vector file[s]')
    sys.exit(0)
else:
    fnames = sys.argv[1:]

# Vector-file handle
tapes=array(list(range(len(fnames))))+8


w2k = utils.W2kEnvironment()
#struct = struct1.Struct(w2k.case)
struct = Struct(w2k.case)
struct.ReadStruct(w2k.case+'.struct')

vectortype=float
vecread = rdVec.fvread3
so=''
if os.path.isfile(w2k.case+".inso") and os.path.getsize(w2k.case+".inso")>0 :
    print('Found '+w2k.case+'.inso file, hence assuming so-coupling exists. Switching -so switch!')
    so = 'so'
    vectortype=complex
    vecread = rdVec.fvread3c


maxkpoints = 10000
# opens vector file

for fname,tape in zip(fnames,tapes):
    
    Elinear = rdVec.fvopen(tape, fname, struct.nat)
    print('linearization energy=', Elinear)
    gw=[]
    for ik in range(maxkpoints):
        # Reads vector file
        (k, kname, wgh, ios, n0, nb) = rdVec.fvread1(tape)
        if ios!=0: break # vector file is finished, no more k-points
        print('k=', k)
    
        # Reciprocal vectors
        Gs = rdVec.fvread2(tape, n0)
    
        # Reading KS eigensystem
        As=zeros((nb,n0), dtype=vectortype)
        Ek=zeros(nb, dtype=float)
        for i in range(nb):
            (num, ek, A) = vecread(tape, n0)
            As[i,:] = A       # KS eigenvector
            Ek[i] = utils.Ry2eV(ek) # KS eigenvalue
    rdVec.fvclose(tape)
    
