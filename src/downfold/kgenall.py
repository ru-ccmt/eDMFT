#!/usr/bin/env python
import sys, os, re
from scipy import *
from scipy import linalg
from utils import W2kEnvironment
import latgen


def read_klist(kfile):
    fh = open(kfile,'r')
    qpts=[]
    for line in fh:
        if line[:3]=='END': break
        qpts.append( map(int,line.split()[1:5]) )
    return qpts

def print_klist(fh, kpts, NN):
    for ik in range(len(kpts)):
        if ik==0:
            print >> fh, "%10d%10d%10d%10d%10d%5.1f -7.0  1.5         0 k, div: ( %2d %2d %2d)" % tuple([ik+1]+kpts[ik]+NN)
        else:
            print >> fh, "%10d%10d%10d%10d%10d%5.1f" % tuple([ik+1]+kpts[ik])
    print >> fh, 'END'

def DirectVectors(br2,br1):
    v0 = cross(br2[1],br2[2])
    v1 = cross(br2[2],br2[0])
    v2 = cross(br2[0],br2[1])
    vol = dot(v0,br2[0])
    return dot(br1, array([v0/vol,v1/vol,v2/vol]))



class kGen(object):
    def __init__(self, case, largeBZ=False):
        (nat,nsym,ndif,lattic,aa,bb,cc,alpha) = latgen.struct_sizes(case+'.struct')
        (self.br1,self.br2,vol,ortho) = latgen.latgen(lattic,aa,bb,cc,alpha)
        self.BS = dot(linalg.inv(self.br1), self.br2)
        #  Brillouin zone for 1/2 atoms per unit cell
        if largeBZ:
            AA = linalg.inv(self.br2).transpose()
            DV = dot(self.br1, AA)
            #DV = DirectVectors(self.br2,self.br1)
            print 'Direct vectors for original Brillouin zone are'
            print '[[',("%6.3f,"*3) % tuple(DV[0]),'],'
            print ' [',("%6.3f,"*3) % tuple(DV[1]),'],'
            print ' [',("%6.3f,"*3) % tuple(DV[2]),']]'
            print 'Please enter the direct vectors for large Brillouin zone'
            DVn = zeros(1)
            while shape(DVn) != (3,3):
                DVn = array(eval(sys.stdin.readline()))
                print 'DVn=', DVn
                # DV = br1 * br2^{-1}^T
                # bs = br1^{-1} * br2
                #   hence
                # bs = br1^{-1} * br1^T  * DV^{-1}^T
            CC = linalg.inv(DVn).transpose()
            self.BS = dot( dot(linalg.inv(self.br1), self.br1.transpose()) , CC )
        print 'BS=', self.BS
        
    def fLCD(self, N0, N1, L=30):
        NNN = sort([N0,N1])
        L0 = arange(1,L+1)*NNN[0]
        L1 = arange(1,L+1)*NNN[1]
        for i1 in range(L):
            for i0 in range(L):
                if L0[i0] == L1[i1]:
                    return L0[i0]
                
    def KmQ(self, qp,N0,N1,N2):
        # Finds least common denominator
        LCD = self.fLCD(self.fLCD(self.fLCD(N0,N1),N2),qp[3])
        #print 'LCD=', LCD
        q = array(qp[:3])/float(qp[3])
        #print 'q=', q
        # creates list for k-q
        kpts=[]
        for i0 in range(N0):
            for i1 in range(N1):
                for i2 in range(N2):
                    k = array(self.BS[0]*i0/float(N0) + self.BS[1]*i1/float(N1) + self.BS[2]*i2/float(N2))
                    k -= q
                    kp = map(round, k*LCD )
                    kpts.append( kp+[LCD,1.0] )
        return kpts




if __name__ == '__main__':

    if len(sys.argv)<2:
        print 'ERROR:'
        print 'Give the list of k-point divitions, i.e., [Nx,Ny,Nz]'
        sys.exit(1)
    else:
        (N0,N1,N2) = eval(sys.argv[1])

    largeBZ=False
    if len(sys.argv)>2:
        largeBZ = True
        print 'Will produce large Brillouin zone!'
    
    # Finds what is case
    w2k = W2kEnvironment()
    # Class for k-list creation
    kgn = kGen(w2k.case,largeBZ)
    # reads all points from klist and uses them for q-mesh
    
    qpts = kgn.KmQ([0,0,0,1],N0,N1,N2)

    #print_klist(sys.stdout, qpts, [N0,N1,N2])
    #sys.exit(1)
    
    fklst = open(w2k.case+'.klist', 'w')
    print_klist(fklst, qpts, [N0,N1,N2])
    fklst.close()

