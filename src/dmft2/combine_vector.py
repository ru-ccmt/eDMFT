import sys, re, os
from numpy import *
#from pylab import *
#from scipy.lib.blas import fblas
import rdVec, utils, findEF
from wstruct import Struct

if len(sys.argv)<2:
    print('Please give the name of vector file[s]')
    sys.exit(0)
else:
    fnames = sys.argv[1:]

# Vector-file handle
tapes=array(list(range(len(fnames))))+8


w2k = utils.W2kEnvironment()
#print 'case=', w2k.case
#EF = findEF.FindChemicalPotential(w2k.case, '')[0]
#struct = struct1.Struct(w2k.case)
struct = Struct(w2k.case)
struct.ReadStruct(w2k.case+'.struct')

vectortype=float
vecread = rdVec.fvread3
vecwrite = rdVec.fvwrite3
so=''
if os.path.isfile(w2k.case+".inso") and os.path.getsize(w2k.case+".inso")>0 :
    print(('Found '+w2k.case+'.inso file, hence assuming so-coupling exists. Switching -so switch!'))
    so = 'so'
    vectortype=complex
    vecread = rdVec.fvread3c
    vecwrite = rdVec.fvwrite3c


maxkpoints = 10000
# opens vector file

heads=[]
all_Gs=[]
all_As=[]
all_Ek=[]
for fname,tape in zip(fnames,tapes):
    Elinear = rdVec.fvopen(tape, fname, struct.nat)
    print(('linearization energy=', Elinear))
    for ik in range(maxkpoints):
        # Reads vector file
        head = rdVec.fvread1(tape)
        (k, kname, wgh, ios, n0, nb) = head
        if ios!=0: break # vector file is finished, no more k-points
        print(('k=', k, n0))
        heads.append(head)
        
        # Reciprocal vectors
        Gs = rdVec.fvread2(tape, n0)
        all_Gs.append(Gs)
        # Reading KS eigensystem
        As=zeros((nb,n0), dtype=vectortype)
        Ek=zeros(nb, dtype=float)
        for i in range(nb):
            (num, ek, A) = vecread(tape, n0)
            As[i,:] = A       # KS eigenvector
            Ek[i] = ek # KS eigenvalue
        all_As.append(As)
        all_Ek.append(Ek)
        
    rdVec.fvclose(tape)

fh_final=9
rdVec.fwopen(fh_final, w2k.case+".vector_out", struct.nat, Elinear)

for i in range(len(heads)):
    print((heads[i]))
    #print len(all_Gs[i][0]),len(all_Gs[i][1]),len(all_Gs[i][2])
    (k, kname, wgh, ios, n0, nb) = heads[i]

    ios = rdVec.fvwrite1(fh_final, k, kname, wgh, n0, nb)
    ios = rdVec.fvwrite2(fh_final, all_Gs[i])
    
    As = all_As[i]
    Ek = all_Ek[i]
    for j in range(nb):
        vecwrite(fh_final, j+1, Ek[j], As[j,:])
    
rdVec.fwclose(fh_final)
