import sys, re, os
from numpy import *
#from pylab import *
import string
#from scipy.lib.blas import fblas
import rdVec, utils#, struct1
from wstruct import Struct

if len(sys.argv)<2:
    print('Please give the name of vector file[s]')
    sys.exit(0)
else:
    fnames = sys.argv[1:]

# Vector-file handle
tapes=array(list(range(len(fnames))))+8


w2k = utils.W2kEnvironment()
struct = Struct(w2k.case)
struct.ReadStruct(w2k.case+'.struct')

maxkpoints = 10000
# opens vector file

Elinears1=[]
Elinears2=[]
heads=[]
Eks=[]
for fname,tape in zip(fnames,tapes):
    
    (Elinear1,Elinear2) = rdVec.feopen(tape, fname, struct.nat)
    #Elinear1 = [(string.join(strng,'')).rstrip() for strng in Elinear1]
    #Elinear2 = [(string.join(strng,'')).rstrip() for strng in Elinear2]
    Elinear1 = [strng.rstrip().decode('UTF-8') for strng in Elinear1]
    Elinear2 = [strng.rstrip().decode('UTF-8') for strng in Elinear2]

    Elinears1.append(Elinear1)
    Elinears2.append(Elinear2)
    
    print('linearization energy 1 =', Elinear1)
    print('linearization energy 2 =', Elinear2)

    for ik in range(maxkpoints):
        # Reads vector file
        head = rdVec.feread1(tape)
        (k, kname, wgh, ios, n0, nb) = head
        if ios!=0: break # vector file is finished, no more k-points
        # Eigenvalues
        Ek = rdVec.feread2(tape, nb)

        heads.append(head)
        Eks.append(Ek)
        print('k=', k)
        
    rdVec.feclose(tape)
    
