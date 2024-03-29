#!/usr/bin/env python
# @Copyright 2007 Kristjan Haule

from scipy import *
from scipy import linalg
from pylab import *
import sys
import scipy

from distutils.version import StrictVersion
if StrictVersion(scipy.__version__) > StrictVersion('0.19.0'):
    import weave
else:
    import scipy.weave as weave

def ReadKlist(fklist, ReadBS=False):
    fk = open(fklist,'r')
    data = fk.readlines()
    nkp = [line[:3]=='END' for line in data].index(True)
    if data[nkp][:3]!='END': 
        print 'wrong klist ', fklist
    kp=[]
    for i in range(nkp):
        kp.append( map(int, [data[i][10:15], data[i][15:20], data[i][20:25], data[i][25:30]]) )

    if (ReadBS):
        BS = [map(float,line.split()) for line in data[nkp+1:nkp+4]]
        BSI = matrix(array(BS).T).I
        return (array(kp), array(BS), array(BSI))
    else:
        return array(kp)

class K_index:
    def __init__(self, BSI, kps):
        self.BSI = BSI
        self.SCALE = kps[0][3]
        self.ind1={}
        for ik,k in enumerate(kps):
            wik = tuple(map(int, dot(BSI,k[:3])))
            self.ind1[wik] = ik
    def __call__(self, ik):
        wik = tuple(map(int, dot(self.BSI,ik[:3])%self.SCALE))
        return self.ind1[wik]





def CheckPPHermisity():
    for ik in range(nkp):
        CA=zeros((norb*norb,norb*norb),dtype=float)
        for i1 in range(norb):
            for i2 in range(norb):
                for i3 in range(norb):
                    for i4 in range(norb):
                        CA[findex2(i1,i2,norb,norb),findex2(i3,i4,norb,norb)] = Chi0PP[ik,i1,i2,i3,i4]
        if sum(abs(CA-transpose(CA)))>1e-3: print 'ERROR'
        ei,ev=linalg.eigh(CA)
        print ik, ei.tolist()
    
def CheckTimeReversal():
    for ik in range(nkp):
        for i1 in range(norb):
            for i2 in range(norb):
                for i3 in range(norb):
                    for i4 in range(norb):
                        imk = k_index(-kps[ik])
                        diff = GammaPH[ik,i1,i2,i3,i4]-GammaPH[imk,i3,i4,i1,i2]
                        if abs(diff)>1e-3:
                            print 'DIFF-1=', ik, i1, i2, i3, i4, GammaPH[ik,i1,i2,i3,i4], GammaPH[imk,i3,i4,i1,i2]
                        diff = GammaPH[ik,i1,i2,i3,i4]-GammaPH[ik,i2,i1,i4,i3]
                        if abs(diff)>1e-3:
                            print 'DIFF-2=', ik, i1, i2, i3, i4, GammaPH[ik,i1,i2,i3,i4], GammaPH[ik,i2,i1,i4,i3]
                            
 
def findex3(i1,i2,i3,n1,n2,n3):
    return (i1*n2+i2)*n3+i3
def findex2(i1,i2,n1,n2):
    return i1*n2+i2
    
if __name__ == '__main__':

    if len(sys.argv)<2:
        print 'ERROR : need input filename'
        print 'The input file should contain: '
        print  'case.klist     # filename with k-list'
        print  'Qlist.dat      # filename with Qlist'
        print  'rmesh.dat      # real axis mesh'
        print  'G_k1r_         # file with real axis k-dependent Grens function'
        print  'G_local1r_     # file with real axis local Grens function'
        print  'chi0_real.     # name of the Bubble on real axis'
        print  'G_k1i_         # imaginary axis k-dependent Greens function'
        print  'G_local1i_     # imaginary axis local Greens function'
        print  'tvertex.dat    # ctqmc local vertex function'
        print  '100            # inverse temperature for bose function in Sq(omega)'
        sys.exit(1)
    fin = open(sys.argv[1], 'r')
    
    fin.next() # case.klist
    fQlist      = fin.next().split()[0] # case.qlist
    #fin.next() # rmesh.dat
    #fin.next() # G_k1r_
    #fin.next() # G_local1r_
    #fin.next() # chi0_real.
    fin.next() # G_k1i_
    fin.next() # G_local1i_
    fvertex     = fin.next().split()[0] # tvertex.dat
    fin.close()
    
    fi=open(fvertex)
    fi.next()  # comment # beta, Nvfl, nomv, nOm nom
    beta = float(fi.next().split()[0])
    fi.close()
    
    print 'beta=', beta
    print 'fQlist=', fQlist
    
    fileC0 = 'Chi0pp.dat'
    fileGpm = 'Gpm.dat'
    fileGmm = 'Gmm.dat'
    
    (kps, BS, BSI) = ReadKlist(fQlist,True)
    k_index = K_index(BSI,kps)
    nkp = len(kps)
    
    GammaPM = loadtxt(fileGpm)  # format is (NQ, Norb**4)
    GammaMM = loadtxt(fileGmm)  # format is (NQ, Norb**4)
    Chi0PP = loadtxt(fileC0)    # format is (NQ, Norb**4)

    if shape(GammaPM)[0]!=nkp:
        print 'len('+fileGpm+') should be nkp, but is not compatible with '+fQlist
    if shape(GammaMM)[0]!=nkp:
        print 'len('+fileGmm+') should be nkp, but is not compatible with '+fQlist
    if shape(Chi0PP)[0]!=nkp:
        print 'len('+fileC0+') should be nkp, but is not compatible with '+fQlist
    
    n4 = shape(GammaPM)[1]
    norb = int(sqrt(sqrt(n4)))
    print 'norb=', norb
    
    GammaPM = GammaPM.reshape((nkp,norb,norb,norb,norb))
    GammaMM = GammaMM.reshape((nkp,norb,norb,norb,norb))
    Chi0PP = Chi0PP.reshape((nkp,norb,norb,norb,norb))
    

    
    print 'shape(GammaPM)=', shape(GammaPM)
    print 'shape(GammaMM)=', shape(GammaMM)
    print 'shape(Chi0PP)=',  shape(Chi0PP)
    BCS=zeros((nkp*norb*norb,nkp*norb*norb),dtype=float)
    chi0=zeros((norb*norb, norb*norb), dtype=float)
    Gamma=zeros((norb*norb, norb*norb), dtype=float)
    for ik1 in range(nkp):
        print 'ik=', ik1
        for ik2 in range(nkp):
            k1 = kps[ik1][:3]
            k2 = kps[ik2][:3]
            ik2mk1 = k_index(k2-k1)
            ik1pk2 = k_index(k1+k2)
            imk2mk1 = k_index(-k1-k2)

            support_code="""
            #line 78 "BCS.py"
            int findex3(int i1, int i2, int i3, int n1, int n2, int n3){
                return (i1*n2+i2)*n3+i3;
            }
            int findex2(int i1, int i2, int n1, int n2){
                return i1*n2+i2;
            }
            """
            code="""
            #line 162 "BCS.py"
            for (int i1=0; i1<norb; i1++){
                for (int i2=0; i2<norb; i2++){
                    for (int i3=0; i3<norb; i3++){
                        for (int i4=0; i4<norb; i4++){
                            int i1i2 = findex2(i1,i2,norb,norb);
                            int i3i4 = findex2(i3,i4,norb,norb);
                            chi0(i1i2,i3i4) = Chi0PP(ik2,i1,i2,i3,i4);
                            //Gamma(i1i2,i3i4) = 0.5*(GammaPM(ik1pk2,i3,i1,i2,i4)+GammaMM(ik2mk1,i4,i1,i2,i3));
                            Gamma(i1i2,i3i4) = 0.5*(GammaPM(ik2mk1,i3,i1,i2,i4)+GammaMM(imk2mk1,i4,i1,i2,i3));
                        }
                    }
                }
            }
            """
            weave.inline(code, ['chi0','Gamma','norb','GammaPM','GammaMM','ik2','ik2mk1','ik1pk2','imk2mk1','Chi0PP'],support_code=support_code,type_converters=weave.converters.blitz, compiler='gcc')
            GammaChi0 = dot(Gamma, chi0)
            code="""
            #line 182 "BCS.py"
            for (int i1=0; i1<norb; i1++){
                for (int i2=0; i2<norb; i2++){
                    for (int i3=0; i3<norb; i3++){
                        for (int i4=0; i4<norb; i4++){
                            int index1 = findex3(ik1,i1,i2, nkp,norb,norb);
                            int index2 = findex3(ik2,i3,i4, nkp,norb,norb);
                            int i1i2 = findex2(i1,i2,norb,norb);
                            int i3i4 = findex2(i3,i4,norb,norb);
                            BCS(index1,index2) = -GammaChi0(i1i2,i3i4)/(nkp);
                        }
                    }
                }
            }
            """
            weave.inline(code, ['BCS','GammaChi0','norb','nkp','ik1','ik2'],support_code=support_code,type_converters=weave.converters.blitz, compiler='gcc')
            
    #print 'Diff=', sum(abs(transpose(BCS)-BCS),axis=None)
    print 'Now diagonalizing matrix of size ', shape(BCS)
    
    evalues,vector = linalg.eig(BCS)

    aevals = real(evalues.real)
    ind = range(len(aevals))
    ind = sorted(ind, key=lambda i: aevals[i])
    
    for i in range(len(ind)):
        print i, evalues[ind[i]], vector[:,ind[i]]
        
    for i in range(-1,-6,-1):
        gs=zeros((nkp,norb*norb),dtype=complex)
        for ik in range(nkp):
            for i1 in range(norb):
                for i2 in range(norb):
                    gs[ik,findex2(i1,i2,norb,norb)]=vector[findex3(ik,i1,i2, nkp,norb,norb),ind[i]]
        savetxt('gs_symmetryr.'+str(abs(i)), real(gs))
        savetxt('gs_symmetryi.'+str(abs(i)), imag(gs))
        
