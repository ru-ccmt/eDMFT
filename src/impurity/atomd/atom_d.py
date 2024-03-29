#!/usr/bin/env python
# @Copyright 2007 Kristjan Haule
# 
from scipy import *
from numpy import *
import sys, re, os, time
import copy
import getopt
import pickle
import glob
import math
from numpy import linalg
from cubic_harmonics import Spheric2Cubic
import dpybind
import gaunt                

def compress(_groups_):
    """ We are given a two dimensional list of numbers which stands for indices of groups of rows or columns of a matrix.
    We need to return similar list of groups, but simplified, which enumerate all coupled rows and columns of a matrix.
    For example, if we have incoming groups: 
        [[0,1,2],[3,4],[0,5]] 
    we want to return
        [[0,1,2,5],[3,4]]
    """
    #print('in  groups=', _groups_)
    groups = [set(g) for g in _groups_] # these groups can be sets, which have optimized routines for checking overlaps.
    loopsDone = True
    while (loopsDone):
        loopsDone = False
        for i in range(len(groups)):
            if loopsDone: break
            for j in range(i+1,len(groups)):
                if loopsDone: break
                if not groups[i].isdisjoint(groups[j]):    # checking of groups[i] and groups[j] have any overlap
                    groups[i] = groups[i].union(groups[j]) # if yes, we take the union of the two
                    del groups[j]                          # an delete the second group, because we have combined group now.
                    loopsDone = True                       # we can than jump out
    groups=[sorted(list(g)) for g in groups]               # at the end we sort each group for convenience
    groups.sort(key=lambda x: x[0])                        # and than we also sort groups by the first element
    #print('out groups=', groups)                          # so that the output is unique
    return groups

def analizeGroups(A, small = 1e-4):
    """Given matrix A, it finds which columns and rows are coupled by nonzero matrix elements of A.
    grp0 contains indices to columns which are coupled
    grp1 contains indices to rows which are coupled
    """
    nonzeros = where(abs(A[:,:])>small) # gives indices to rows&columns to all nonzero values of A.
    grp0=[]
    for i in range(shape(A)[0]):
        u = nonzeros[1][nonzeros[0]==i] # which columns are coupled through row i
        if len(u):
            grp0.append(set(u))         # these columns are coupled through one of the rows
    grp1=[]
    for j in range(shape(A)[1]):
        u = nonzeros[0][nonzeros[1]==j] # which rows are coupled through column j
        if len(u):
            grp1.append(set(u))         # these rows are coupled through one of the columns

    grp0 = compress(grp0)  # simplifying groups so that all columns which are coupled by matrix A are in the same group.
    grp1 = compress(grp1)  # simplifying groups so that all rows which are coupled by matrix A are in the same group.
    return (grp1, grp0)

def coupled(A, groups0, groups1, small = 1e-4):
    """ given two sets of indeces of a matrix, determines 
        which groups are coupled due to nonzero matrix elements in A
        If a group from groups0 is connected with a group from groups1, 
        fpair gives index of the group in groups1, i.e., 
        return : fpair[<a group from groups0>] = <index of a group from groups1>
        A slower but more readable version of the code is:

           fpair = -ones(len(groups0),dtype=int)
           for i,ig0 in enumerate(groups0):
               for j,ig1 in enumerate(groups1):
                   if any(abs(A[ig0,:][:,ig1])>small):
                       fpair[i]=j
    """
    fpair = [-1]*len(groups0) # -1 means group from groups0 is disjoint with all groups in group1
    for ii,ig0 in enumerate(groups0):
        # fancy indexing A[ig0,:] returns 2D array with a selected set of of rows of matrix A.
        # We checked if any row from a group ig0 has nonzero value, and if yes, with which column.
        # We want to know which columns of matrix A are coupled with one of the rows in ig0.
        # The non-zero columns are stored in nonz. Note that numpy function where 
        # returns two arrays, index in ig0, and index in second dimension. We only need the latter.
        nonz = set(where( abs(A[ig0,:])>small )[1])
        # Now we replace rows in nonz with groups in groups1, and fpair will contain index of the
        # group in groups1 which has overlap with row A.
        for jj,jg1 in enumerate(groups1):
            if not nonz.isdisjoint(jg1):
                fpair[ii] = jj
    
    return fpair

def cprint2(fh, U):
    for i in range(shape(U)[0]):
        for j in range(shape(U)[1]):
            f = U[i,j]
            if abs(f)<1e-10: f = 0j
            print("%12.8f %12.8f*i " % (f.real, f.imag), end=' ', file=fh)
        print(file=fh)
    print(file=fh)
def cprint(fh, U):
    for i in range(shape(U)[0]):
        for j in range(shape(U)[1]):
            f = U[i,j]
            if abs(f)<1e-10: f = 0j
            print("%7.4f %7.4f*i " % (f.real, f.imag), end=' ', file=fh)
        print(file=fh)
    print(file=fh)
def mprint(fh, U):
    for i in range(shape(U)[0]):
        for j in range(shape(U)[1]):
            f = U[i,j]
            if abs(f)<1e-10: f = 0
            print("%7.4f " % f, end=' ', file=fh)
        print(file=fh)
    print(file=fh)


def CoulUsC2_slow(l, T2C):
    # Precomputes gaunt coefficients for speed
    # shape(gck)=(4,7,7,4). It contains gck(l, m4, m1, k)
    gck = gaunt.cmp_all_gaunt()
    #save('gck.npy', gck)
    #save('T2C.npy', T2C)
    
    #print 'Gaunt coefficients precomputed - shape(gck)', shape(gck)
    mw = 2*l+1
    if len(T2C) == mw:
        nw = mw
        ns = 1
    elif len(T2C) == 2*mw:
        nw = 2*(2*l+1)
        ns = 2
    else:
        print("ERROR in atom_d.py: T2C has wrong shape")
        sys.exit(0)
    
    T2Cp = conj(T2C.transpose())
    
    UC = zeros((l+1, nw, nw, nw, nw), dtype=complex)
    shft = 3-l
    shft2 = 2*l+1
    Sum1 = zeros((nw,nw,shft2*2),dtype=complex)
    Sum2 = zeros((nw,nw,shft2*2),dtype=complex)
    
    for k in range(l+1):
        Sum1[:,:,:]=0
        for i4 in range(nw):
            for i1 in range(nw):
                for m4 in range(mw):
                    for m1 in range(mw):
                        for s in range(ns):
                            Sum1[i4,i1,m1-m4+shft2] += T2Cp[i4,m4+s*mw]*gck[l,shft+m4,shft+m1,k]*T2C[m1+s*mw,i1]
        Sum2[:,:,:]=0
        for i3 in range(nw):
            for i2 in range(nw):
                for m3 in range(mw):
                    for m2 in range(mw):
                        for s in range(ns):
                            Sum2[i3,i2,m3-m2+shft2] += T2Cp[i3,m3+s*mw]*gck[l,shft+m2,shft+m3,k]*T2C[m2+s*mw,i2]
        for i4 in range(nw):
            for i3 in range(nw):
                for i2 in range(nw):
                    for i1 in range(nw):
                        csum=0
                        for dm in range(shft2*2):
                            csum += Sum1[i4,i1,dm]*Sum2[i3,i2,dm]
                        UC[k,i4,i3,i2,i1] = csum
    return UC

def CoulUsC2(l, T2C):
    # Precomputes gaunt coefficients for speed
    # shape(gck)=(4,7,7,4). It contains gck(l, m4, m1, k)
    gck = gaunt.cmp_all_gaunt()
    gck = array(gck, order='C')   # This is essential for pybind11 code to work with fortran gck
    T2C = array(T2C, order='C')
    nw = len(T2C)
    UC = zeros((l+1, nw, nw, nw, nw), dtype=complex)
    dpybind.FromSlaterToMatrixU(UC, gck, l, T2C)
    return UC


class operateLS(object):
    def __init__ (self, baths, T2C, Q3d):
        self.baths = baths
        self.Nband = int(self.baths/2)
        self.N = self.baths
        
        self.mask=[]
        for i in range(self.N): self.mask.append(1<<i);

        self.T2C = T2C
        self.Q3d = Q3d

        if not self.Q3d:
        #############################################
        # Here for 5d's where spin-orbit is kept    #
        #############################################
            self.Q3d=False
            M2=[]
            l=(self.Nband-1)/2
            #print 'L is here ', l
            for s in [0.5,-0.5]:
                for m in range(-l,l+1):
                    M2.append( (m+2*s)/2.)
            #print 'M2=',M2
            self.M2a=zeros((len(M2),len(M2)),dtype=float)
            for a in range(len(M2)):
                for b in range(len(M2)):
                    for ms in range(len(M2)):
                        self.M2a[a,b] += real(conj(T2C[ms,a])*T2C[ms,b]*M2[ms])
            #print 'M2a=', self.M2a        
        else:
        ####################################################
        # Here only for 3d's where spin-orbit is neglected #
        ####################################################
            self.Q3d=True
            self.bi=[] # band index
            self.sz=[] # sz
            for i in range(self.Nband):
                self.sz.append(1);
                self.bi.append(i);
            for i in range(self.Nband):
                self.bi.append(i)
                self.sz.append(-1)
                
            self.mask_u = []
            self.mask_d = []
            for i in range(self.Nband):
                self.mask_u.append(self.mask[i])
            for i in range(self.Nband):
                self.mask_d.append(self.mask[self.Nband+i])
        
    def Nel(self, state):
        n=0
        for k in self.mask:
            if (k&state): n+=1
        return n
    def occup(self, state):
        """ gives a list of occupancies per band [n_{band1},n_{band2},...]
        """
        oc=[]
        for i in range(self.N):
            if state & self.mask[i]: oc.append(1)
            else: oc.append(0)
        return oc
        
    def sign(self, state, mask_min, mask_max):
        """ Sign when electron hops from mask_min to mask_max
        Counts number of electrons between the two spaces
        """
        # mask will run between mask_min to mask_max
        mask = mask_min<<1 
        n=0           # number of electrons between mask_min and mask_max
        while (mask<mask_max):    # loop to mask_max
            if (mask&state): n+=1 # found electron between the two places
            mask = mask<<1        # increment the mask
        return 1-2*(n%2)          # (-1)^n
    
    def sign_(self, state, mask_max):
        """ Sign when electron is added to the state (from the left)
        """
        # mask will run between mask_min to mask_max
        mask = 1
        n=0           # number of electrons between mask_min and mask_max
        while (mask<mask_max):    # loop to mask_max
            if (mask&state): n+=1 # found electron between the two places
            mask = mask<<1        # increment the mask
        return 1-2*(n%2)          # (-1)^n
    
    def N_el_before(self, state, i):
        n=0
        for q in range(i):
            if self.mask[q]&state: n+=1
        return n

    def DM(self, state):
        "Density matrix"
        DenM = [[[] for i in range(self.N)] for j in range(self.N)]
        for j in range(self.N):
            if not self.mask[j]&state: continue # c_j operator
            jsig = self.sign_(state, self.mask[j])
            nst = state^self.mask[j]
            for i in range(self.N):
                if self.mask[i]&nst: continue   # c_i^\dagger operator
                nstate = nst^self.mask[i]
                isig = self.sign_(nst, self.mask[i])
                DenM[i][j].append( (nstate, jsig*isig) )
        return DenM
    
    def Fp(self, state, ib):
        """ This implements psi^dagger_{ib} operator acting on state
        indexes are:
          ib - band+spin index
        """
        if state&self.mask[ib]: return (0,1)  # This state is already occupied
        newstate = state^self.mask[ib]
        sig = self.sign_(state, self.mask[ib])
        return (newstate, sig)

        
    def CoulombU(self, state, UC, FkoJ, Ising=False):
        sts=[]
        ni=-1
        maxk=l+1
        if (self.Q3d):
            ### will evaluate  again <sts| U(b,a,j,i) psi^+_b psi^+_a psi_j psi_i | state>, but
            ### this time U is defined in (m4,m3,m2,m1) basis only, and is missing the spin component.
            ### Need to make sure that s_1==s_4 and s_2==s_3
            ### hence  <sts| U(m4,m3,m2,m1) psi^+_{m4,s} psi^+_{m3,s'} psi_{m2,s'} psi_{m1,s} | state>
            for i in range(self.N):
                if not(self.mask[i]&state) : continue # (i,m1) does not exists
                ni+=1
                state1 = state^self.mask[i]
                m1 = self.bi[i]
                s1 = self.sz[i]
                nj=-1
                for j in range(self.N):
                    if not(self.mask[j]&state1) : continue  # (j,m2) does not exists
                    nj+=1
                    # here we have: mask[i]&state && mask[j]&state
                    state2 = state1^self.mask[j]
                    m2 = self.bi[j]
                    s2 = self.sz[j]
                    for a in range(self.N): # (a,m3) exists
                        if self.mask[a]&state2: continue
                        if self.sz[a]!=s2 : continue # s3 == s2
                        na = self.N_el_before(state2,a)
                        state3 = state2^self.mask[a]
                        m3 = self.bi[a]
                        s3 = self.sz[a]
                        for b in range(self.N): # (b,m4) exists
                            if self.mask[b]&state3 : continue
                            if self.sz[b]!=s1: continue # s4 == s1
                            nb = self.N_el_before(state3,b)
                            state4 = state3^self.mask[b]
                            m4 = self.bi[b]
                            s4 = self.sz[b]
            
                            if Ising and state4!=state: continue
                            
                            sign = 1-2*((ni+nj+na+nb)%2)
                            U0 = sign*UC[0,m4,m3,m2,m1]*FkoJ[0]
                            
                            dsum=0
                            for k in range(1,maxk):
                                dsum += UC[k,m4,m3,m2,m1]*FkoJ[k]                            
                            U1 = sign*dsum
            
                            if (abs(U0)>1e-6 or abs(U1)>1e-6): sts.append([state4, [U0,U1]])
            
        else :  # This is used for 5d, but not for 3d
            ### will evaluate  <sts| U(b,a,j,i) psi^+_b psi^+_a psi_j psi_i | state>
            for i in range(self.N):
                if not(self.mask[i]&state) : continue # i should exist, otherwise continue
                ni+=1
                state1 = state^self.mask[i]
                nj=-1
                for j in range(self.N):
                    if not(self.mask[j]&state1) : continue  # j should exist, otherwise continue
                    nj+=1
                    # here we have: mask[i]&state && mask[j]&state
                    state2 = state1^self.mask[j]
                    for a in range(self.N): 
                        if self.mask[a]&state2: continue # a should not exist exist
                        na = self.N_el_before(state2,a)
                        state3 = state2^self.mask[a]
                        for b in range(self.N): 
                            if self.mask[b]&state3: continue # b should not exist
                            nb = self.N_el_before(state3,b)
                            state4 = state3^self.mask[b]
                            
                            if Ising and state4!=state: continue
                            sign = 1-2*((ni+nj+na+nb)%2)
                            U0 = sign*UC[0,b,a,j,i]*FkoJ[0]
                            
                            dsum=0
                            for k in range(1,maxk):
                                dsum += UC[k,b,a,j,i]*FkoJ[k]
                            U1 = sign*dsum
            
                            if (abs(U0)>1e-6 or abs(U1)>1e-6): sts.append([state4, [U0,U1]])
        return sts
    def CoulombUIsing(self, state, UC, FkoJ):
        sts=[]
        ni=-1
        maxk=l+1
        if (self.Q3d):
            ### will evaluate  again <sts| U(b,a,j,i) psi^+_b psi^+_a psi_j psi_i | state>, but
            ### this time U is defined in (m4,m3,m2,m1) basis only, and is missing the spin component.
            ### Need to make sure that s_1==s_4 and s_2==s_3
            ### hence  <sts| U(m4,m3,m2,m1) psi^+_{m4,s} psi^+_{m3,s'} psi_{m2,s'} psi_{m1,s} | state>
            for i in range(self.N):
                if not(self.mask[i]&state) : continue # (i,m1) does not exists
                ni+=1
                state1 = state^self.mask[i]
                m1 = self.bi[i]
                s1 = self.sz[i]
                nj=-1
                for j in range(self.N):
                    if not(self.mask[j]&state1) : continue  # (j,m2) does not exists
                    nj+=1
                    # here we have: mask[i]&state && mask[j]&state
                    state2 = state1^self.mask[j]
                    m2 = self.bi[j]
                    s2 = self.sz[j]
                    for a in [i,j]: # (a,m3) exists
                        if self.mask[a]&state2: continue
                        if self.sz[a]!=s2 : continue # s3 == s2
                        na = self.N_el_before(state2,a)
                        state3 = state2^self.mask[a]
                        m3 = self.bi[a]
                        s3 = self.sz[a]
                        for b in [i,j]: # (b,m4) exists
                            if self.mask[b]&state3 : continue
                            if self.sz[b]!=s1: continue # s4 == s1
                            nb = self.N_el_before(state3,b)
                            state4 = state3^self.mask[b]
                            m4 = self.bi[b]
                            s4 = self.sz[b]
            
                            if state4!=state: continue
                            
                            sign = 1-2*((ni+nj+na+nb)%2)
                            U0 = sign*UC[0,m4,m3,m2,m1]*FkoJ[0]
                            
                            dsum=0
                            for k in range(1,maxk):
                                dsum += UC[k,m4,m3,m2,m1]*FkoJ[k]                            
                            U1 = sign*dsum
            
                            if (abs(U0)>1e-6 or abs(U1)>1e-6): sts.append([state4, [U0,U1]])
                            
                            #if (state==1023): 
                            #    print 'Last:', state4, i, j, a, b, U1 
            
        else :  # This is used for 5d, but not for 3d
            ### will evaluate  <sts| U(b,a,j,i) psi^+_b psi^+_a psi_j psi_i | state>
            for i in range(self.N):
                if not(self.mask[i]&state) : continue # i should exist, otherwise continue
                ni+=1
                state1 = state^self.mask[i]
                nj=-1
                for j in range(self.N):
                    if not(self.mask[j]&state1) : continue  # j should exist, otherwise continue
                    nj+=1
                    # here we have: mask[i]&state && mask[j]&state
                    state2 = state1^self.mask[j]
                    for a in [i,j]: 
                        if self.mask[a]&state2: continue # a should not exist exist
                        na = self.N_el_before(state2,a)
                        state3 = state2^self.mask[a]
                        for b in [i,j]: 
                            if self.mask[b]&state3: continue # b should not exist
                            nb = self.N_el_before(state3,b)
                            state4 = state3^self.mask[b]
                            
                            if state4!=state: continue
                            
                            sign = 1-2*((ni+nj+na+nb)%2)
                            U0 = sign*UC[0,b,a,j,i]*FkoJ[0]
                            
                            dsum=0
                            for k in range(1,maxk):
                                dsum += UC[k,b,a,j,i]*FkoJ[k]
                            U1 = sign*dsum
            
                            if (abs(U0)>1e-6 or abs(U1)>1e-6): sts.append([state4, [U0,U1]])
        return sts

    def printn(self, state):
        sstate=''
        if self.Q3d:
            for i in range(self.Nband):
                if (state & self.mask_u[i]) and (state & self.mask_d[i]) : sstate += '2'
                elif (state & self.mask_u[i]) : sstate += 'u'
                elif (state & self.mask_d[i]) : sstate += 'd'
                else : sstate += '0'
        else:
            for i in range(self.N):
                if state & self.mask[i]: sstate += '1'
                else : sstate += '0'
        return sstate
    
    def Mz(self,state):
        m2=0.0
        for i in range(self.N):
            if state&self.mask[i]: m2 += self.M2a[i,i]
        return m2
        
    #########################
    
    def Sz(self, state):
        if not self.Q3d: return 0
        nu = 0
        nd = 0
        for i in range(self.Nband):
            if state&self.mask_u[i] : nu += 1
            if state&self.mask_d[i] : nd += 1
        return nu-nd
    
    def S2(self, state):
        l2p1 = self.Nband
        sts=[]
        # calculating \sum_{a,b} S^z_a S^z_b + 1/2*(S^+_a S^-_b + S^-_a S^+_b)
        # diagonal part
        dd=0;
        for ilz in range(l2p1):
            up=0; dn=0  # this comes from 1/2(S^+_a S^-_a + S^-_a S^+_a)
            if self.mask_u[ilz] & state: up = 1
            if self.mask_d[ilz] & state: dn = 1
            # if only up or only down in certain lz
            if up+dn==1: dd += 0.5
        # Sz^2 + diagonal part of S^+ S^-
        fct = (0.5*self.Sz(state))**2 + dd
        # store diagonal
        sts.append([state,fct])
        # off diagonal
        for ilz in range(l2p1):
            im1 = self.mask_u[ilz]
            im2 = self.mask_d[ilz]
            ib1 = bool(state & im1)
            ib2 = bool(state & im2)
            if ib1 and not ib2: # S^-_i gives nonzero
                isig = self.sign(state, min(im1,im2), max(im1,im2))
                istate = state^im1^im2
                for jlz in range(l2p1):
                    if (ilz==jlz): continue
                    jm1 = self.mask_d[jlz]
                    jm2 = self.mask_u[jlz]
                    jb1 = bool(state & jm1)
                    jb2 = bool(state & jm2)
                    if jb1 and not jb2: # S^+_j gives nonzero
                        jsig = self.sign(istate, min(jm1,jm2), max(jm1,jm2))
                        jstate = istate^jm1^jm2
                        sts.append([jstate, isig*jsig])
        return sts

    def S_minus(self, state):
        l2p1 = self.Nband
        sts=[]
        # off diagonal
        for ilz in range(l2p1):
            im1 = self.mask_u[ilz]
            im2 = self.mask_d[ilz]
            ib1 = bool(state & im1)
            ib2 = bool(state & im2)
            if ib1 and not ib2: # S^-_i gives nonzero
                isig = self.sign(state, min(im1,im2), max(im1,im2))
                istate = state^im1^im2
                sts.append([istate, isig])
        return sts
    
    def PairHop(self, state):
        """ Computes the pair-hopping term:  D_a^\dagger D_a , where D_a creates or
        anhilates a double occupied site. There is no minus sign in this term!
        """
        doubles=[]
        empty=[]
        for i in range(self.Nband):
            if state & self.mask_u[i] and state & self.mask_d[i]:
                doubles.append(i)
            elif not(state & self.mask_u[i]) and not(state & self.mask_d[i]):
                empty.append(i)

        rst=[]
        for id in doubles:
            nst1 = state^self.mask_u[id]^self.mask_d[id]
            for ie in empty:
                nst2 = nst1^self.mask_u[ie]^self.mask_d[ie]
                rst.append(nst2)
        return rst
                
    def NDouble(self, state):
        ne = 0
        for i in range(self.Nband):
            if state & self.mask_u[i] and state & self.mask_d[i]: ne += 1
        return ne
            
    def OneBodyNab(self, state, Sab):
        """ computing the term Sab[a,i] f^+_a f_i
        returns all matrix elements generated by the above one-body term
        when acting on state
        """
        sts=[]
        ni=-1
        for i in range(self.baths):
            if not(self.mask[i]&state) : continue
            ni+=1
            state1 = state^self.mask[i]
            m1 = self.bi[i]
            s1 = self.sz[i]
            # here we have: mask[i]&state
            for a in range(self.baths):
                if self.mask[a]&state1 : continue # sz_a == sz_j
                # here we have: state&mask[i] and not(state1&mask[a])
                na = self.N_el_before(state1,a)
                state2 = state1^self.mask[a]
                m2 = self.bi[a]
                s2 = self.sz[a]
                        
                sign = 1-2*((ni+na)%2)
                
                nab = sign*Sab[a,i]
                        
                if (abs(nab)>1e-6): sts.append([state2, nab])
        return sts


def baseN(Nrange, prop,Q3d):
    Ntot = len(prop)
    wstates=[]

    if Q3d:
        for n1 in Nrange: # range(Nband*2+1):
            for sz1 in range(-n1,n1+1,2):
                states=[]
                for i in range(Ntot):
                    if prop[i][0]==n1 and prop[i][1]==sz1:
                        states.append(i)
                        
                if (len(states)>0): wstates.append([n1, sz1, states])
    else:
        for n1 in Nrange:  # range(Nband*2+1):
            states=[]
            for i in range(Ntot):
                if prop[i][0]==n1:
                    states.append(i)
            if (len(states)>0): wstates.append([n1, 0, states])
        
    return wstates

def list_to_string(x):
    return str(array(x).flatten().tolist())

def comp(x, y):
    if x[2]!=y[2]: return int(x[2]-y[2])
    else:
        if abs(x[3]-y[3])<1e-5: return 0
        elif (x[3]<y[3]): return -1
        else: return 1

def SpinOrbitM(l,T2C):
    # one electron |l,m,s> base
    ms_base=[]
    for s in [1/2.,-1/2.]:
        for m in range(-l,l+1):
            ms_base.append([m,s])
    #print 'ms_base=', ms_base
    
    # one electron |j,mj> base
    pj = [l-1/2.,l+1/2.]
    if l==0 : pj = [0.5]
    jj_base=[]
    for j in pj:
        for mj in arange(-j,j+1):
            jj_base.append([j, mj])
    #print 'jj_base=', jj_base
    
    # transforms between |lms> and |jmj> base
    Tjls = zeros((len(ms_base),len(jj_base)))
    for ijj,jj in enumerate(jj_base):
        for ims,ms in enumerate(ms_base):
            Tjls[ijj,ims] = gaunt.clebschg(jj[0], jj[1], l, ms[0], 1/2., ms[1])

    # the one-body operator l*s in matrix form
    # in the j-j base
    jSO = zeros((len(jj_base),len(jj_base)))
    for ijj,jj in enumerate(jj_base):
        jSO[ijj,ijj] = 0.5*(jj[0]*(jj[0]+1) - l*(l+1) - 3/4.)
    
    # changing to lms base
    mSO = matrix(Tjls.transpose())*matrix(jSO)*matrix(Tjls)

    # creating large T2C base
    if ( len(T2C) < len(jj_base) ):
        T2Cl = zeros(tuple(array(shape(T2C))*2),dtype=complex)
        T2Cl[:len(T2C),:len(T2C)] = T2C
        T2Cl[len(T2C):,len(T2C):] = T2C
    else:
        T2Cl = T2C
    
    # changing to cubic harmonics base
    cSO = matrix(conj(T2Cl.transpose())) * mSO * matrix(T2Cl)

    print('spin-orbit=', file=fh_info)
    mprint(fh_info,real(cSO))
    
    return cSO


def Diagonalize(Ham, small=1e-4, fh=sys.stdout):
    """ Diagonalization is done in blocks. This is not because of efficiency but because
    the resulting eigenvectors must not mix states of direct base if not absolutely necessary.
    If brute force diagonalization is used in large scale problems, eigenvectors can be seriously
    mix direct states with different symmetry.
    """
    diff = sum(Ham-Ham.conj().T) #transpose(conj(Ham)))
    if abs(diff)>1e-6:
        print('H NOT HERMITIAN!')
    
    # Check block structure of Hamiltonian
    # States which are mixed in Ham will have the same blck[i]
    ndim = len(Ham)
    blck=list(range(ndim))
    for i in range(ndim):
        for j in range(i+1,ndim):
            if (abs(Ham[i][j])>small):
                commonb = min(blck[i],blck[j])
                for k in range(ndim):
                    if blck[k] in [blck[i],blck[j]]: blck[k]=commonb
        #print ('%2d'%i), 'current blck=', '%2d,'*len(blck) % tuple(blck)
        
    # Having blck[i] a new array block[:][:] is created, which contains indexes to all blocks
    # for example [[1,2,3],[4,5,6]] for Hamiltonian containing two blocks
    block=[]
    for i in range(ndim):
        bb=[]
        for j in range(ndim):
            if blck[j]==i: bb.append(j)
        if len(bb)>0:
            block.append(bb)
    #print 'block=', block
    
    # Here we go over all blocks and diagonalize each one.
    eigv=[] # contains all eigenvalues
    eigx=[] # contains full eigenvectors
    for ibl,bl in enumerate(block):
        hs = zeros((len(bl),len(bl)), dtype=complex)
        for i,ib in enumerate(bl):
            for j,jb in enumerate(bl):
                hs[i,j] = Ham[ib,jb]

        eigy = linalg.eigh(hs)

        print('Eigenvalues[',bl,']=',eigy[0], file=fh)
        
        # Checking if eigenvectors are complex!
        for l in range(len(eigy[1])):
            imax=0
            #print 'shape(eigy[1])', shape(eigy[1])
            for iu in range(len(eigy[0])):
                #print iu, imax
                if abs(eigy[1][iu,l])>abs(eigy[1][imax,l]): imax=iu
            z=eigy[1][imax,l]
            phi=math.atan2(z.imag,z.real)
            eigy[1][:,l] *= exp(-phi*1j)
            
            ime = sum([abs(x.imag) for x in eigy[1][:,l]])
            if (abs(ime))<1e-10: ime=0
            print('im=%2d %2d %f' % (ibl, l, ime), file=fh)
            #ime = sum([abs(eigy[1][u,l].imag) for u in range(len(eigy[1]))])
        #    if ime>1e-7: print 'TROUBLES!!! Complex eigenvector! You sould improve that!'
        
        # Creating a big eigenvector with all components
        for l in range(len(eigy[1])):
            large_eig=zeros(ndim, dtype=complex)
            small_eig = eigy[1][:,l]
            for m,mb in enumerate(bl):  large_eig[mb] = small_eig[m]
            eigx.append(large_eig)
        eigv += eigy[0].tolist()

    # Now we need to sort eigenvectors and eigenvalues
    # index is created for sorting
    #indx=list(range(ndim))    
    #indx.sort(key=lambda a: eigv[a]) #lambda a,b: cmp(eigv[a],eigv[b]))
    indx = sorted(range(ndim), key=lambda a: eigv[a])
    # and actual sorting is performed
    seigv = [eigv[i] for i in indx]
    seigx = [eigx[i] for i in indx]
    #seigv=[]
    #seigx=[]
    #for i in range(ndim):
    #    seigv.append(eigv[indx[i]])
    #    seigx.append(eigx[indx[i]])

    # Eigenvectors should be in the form Ham*v[:,i] = w[i]*v[:,i]
    # which means that we need to transpose the list of eigenvectors
    seigx = array(seigx).T
    seigv = array(seigv)
    
    # We also do a brute force diagonalization just to check if something goes wrong with block diagonalization
    # Note that the two resulting eigensystems are not necessary the same due to freedom in choosing eigenvectors
    #eig = linalg.eigh(Ham)
    eig = linalg.eigvalsh(Ham)

    #print('eig=', eig, 'Ham=', Ham)
    # If eigenvalues from block diagonalization and full diagonalization are different, something is wrong
    if sum(abs(eig-seigv))>small:
        print('!!!!!TEZAVE!')
        print('The right eigenvalues are:', eig)
    
    return [seigv,seigx]


def thesame(mx,my,small=1e-3):
    if list(mx.keys()) != list(my.keys()): return False
    for k in list(mx.keys()):
        if abs(mx[k]-my[k])>small: return False
    return True
    
def VEquivalentStates(mps,ind):
    """ Finds which states have the same bubbles """
    wx = [(i,mps[i]) for i in range(len(mps))]
    iequiv=[]
    while len(wx)>0:
        mx = wx[0][1]
        j=0
        rr=[]
        while j < len(wx):
            if thesame(mx,wx[j][1]):
                rr.append(wx[j][0])
                del wx[j]
            else: j+=1
        iequiv.append(rr)
    
    for ik in range(len(iequiv)):
        for ij in range(len(iequiv[ik])):
            iequiv[ik][ij] = ind[iequiv[ik][ij]]
            
    return iequiv

def AverageBubbles(tmps):
    """ Compute average over 'almost' equivalent states """
    trmp=[]
    for mps in tmps:
        all_keys=[]
        for mp in mps:
            all_keys = union(all_keys,list(mp.keys()))
        
        rmp={}
        for k in all_keys:
            sm=0.0
            for mp in mps:
                if k in mp:
                    sm += mp[k]
            #sm/=len(mps)
            rmp[k]=sm
        trmp.append(rmp)    
    return trmp


def EquivalentStates(ipE, ipN):
    iequiv=[]
    equiv = list(range(len(ipE)))
    leq=0
    ju=0
    Nmax = ipN[-1]
    for Ni in range(Nmax+1):
        # all states of the same N are in the interval [ju,je]
        je=ju
        while je<len(ipN) and ipN[je]==Ni: je+=1

        ind = list(range(ju,je))
        ind.sort(lambda x,y: cmp(ipE[x],ipE[y]) )

        #print Ni
        i0=0
        while (i0<len(ind)):
            Ec = ipE[ind[i0]]
            ieq=[]
            #print 'Ec=', Ec
            while i0<len(ind) and abs(ipE[ind[i0]]-Ec)<1e-10:
                #print ind[i0], ipE[ind[i0]], leq
                equiv[ind[i0]] = leq
                ieq.append(ind[i0])
                i0+=1
            leq += 1
            iequiv.append(ieq)
            #print
        #print
        ju=je
    return (equiv, iequiv)

def RenumberStates(pseudostates, Enes, wstates, S2ws):
    # renumbers states such that each of 1024 states has unique index
    # also remembers energy and N for each state
    ij=0
    puniq={}
    ipuniq=[]
    ipE=[]
    ipN=[]
    ipS=[]
    for ii,iwp in enumerate(pseudostates):
        wdim = len(Enes[ii])
        for j in range(wdim):
            puniq[(ii,j)]=ij
            ipuniq.append((ii,j))
            ipE.append(Enes[ii][j])
            ipS.append(S2ws[ii][j])
            wN = sum(wstates[iwp[0]][0])
            ipN.append(wN)
            ij+=1
    return (puniq, ipE, ipN, ipS)


def CreateEmpty3D_Dict(n0,n1,n2):
    return [[[{} for i2 in range(n2)] for i1 in range(n1)] for i0 in range(n0)]
def CreateEmpty2D_Dict(n0,n1):
    return [[{} for i1 in range(n1)] for i0 in range(n0)]

def ReadTrans(filename, fh_info):
    """Read the self-energy index file Sigind and the local transformation matrix CF from a file"""
    with open(filename, 'r') as fh:
        data = fh.readlines()
    
    (n1,n2) = list(map(int, data[0].split()[:2]))

    Sigind=[]
    for i in range(n1):
        Sigind.append( list(map(int, data[i+2].split()[:n2])) )
    Sigind = array(Sigind)

    print('len(data)', len(data), file=fh_info)
    print('n1=', n1, file=fh_info)
    if len(data) >= n1+n1+3:
        n2 = n1
        CF=[]
        for i in range(n2):
            cl = array(list(map(float, data[n1+3+i].split())))
            CF.append( cl[0::2]+cl[1::2]*1j )
        CF = array(CF)
    elif len(data)>=int(n1+n1/2+3):
        n2 = int(n1/2)
        CF=[]
        for i in range(n2):
            cl = array(list(map(float, data[n1+3+i].split())))
            CF.append( cl[0::2]+cl[1::2]*1j )
        CF = array(CF)
        CFN = zeros((2*n2,2*n2), dtype=complex)
        CFN[:n2,:n2] = CF
        CFN[n2:,n2:] = CF
        CF = CFN
    else:
        CF = identify(n1)
    print('shape(CF)=', shape(CF), 'shape(Sigind)=', shape(Sigind), file=fh_info)
    
    return (Sigind, CF)

def SlaterF(U, J, l):
    Fk = zeros((4,4), dtype=float)
    if l==0:
        # F0 for s-electrons
        Fk[0,0] = U
    elif l==1:
        # F2 for p-electrons
        Fk[0,1] = U
        if type(J) is list:
            Fk[1,1] = 5*J[0]
        else:
            Fk[1,1] = 5*J
    elif l==2:
        # F2 and F4 for d-electrons
        Fk[0,2] = U
        if type(J) is list:
            Fk[1,2] = 14./1.625 * J[0]
            Fk[2,2] = 14.*0.625/1.625 * J[1]
        else:
            Fk[1,2] = 14./1.625 * J
            Fk[2,2] = 14.*0.625/1.625 * J
    elif l==3:
        # F2, F4 and F6 for f-electrons
        Fk[0,3] = U
        if type(J) is list:
            Fk[1,3] = 6435./(286+195*0.668+250*0.494) * J[0]
            Fk[2,3] = 0.668*6435./539.76 * J[1]
            Fk[3,3] = 0.494*6435./539.76 * J[2]
        else:
            Fk[1,3] = 6435./(286+195*0.668+250*0.494) * J
            Fk[2,3] = 0.668*6435./539.76 * J
            Fk[3,3] = 0.494*6435./539.76 * J
    return Fk

def Check_T2C_Real(T2C, l, fh_info, small):
    """
     Here we added a routine which checks that cubic harmonics are real.
     Only in this case operators F^+ and F used in ctqmc will be real.
     Otherwise these matrix elements might be complex.
     The condition for cubic harmonics to be real is:
    
      Imag( T2C[m,i] + (-1)**m * T2C[-m,i] ) == 0
        and
      Real( T2C[m,i] - (-1)**m * T2C[-m,i] ) ==0
    
     which follows from the requirement: \sum_m T2C[m,i]*exp(i*m*phi) is real for any phi

     We can also derive this from the definition
        T2C[m,i]==<Y_{lm}|Phi_i>
     and the fact that localized orbital |Phi> is real, hence
       T2C[m,i]^* = <Y_{lm}|Phi_i>^* = <Y_{lm}^*|Phi_i> = (-1)^m <Y_{l,-m}|Phi_i> = (-1)^m T2C[l,-m]
       
     We are free to add any phase to cubic harmonics, hence T2C[m,i] -> T2C[m,i]*exp(i*phi_i)
       with phi_i being arbitrary
    
     This leads to the following 2x2 system of equations:
      ( Rp[m,i], Qp[m,i] ) ( sin(phi_i) )
      ( Qm[m,i], Rm[m,i] ) ( cos(phi_i) ) = 0
    
     where
          Qp[m,i] = Imag( T2C[m,i] + (-1)**m * T2C[-m,i] )
          Rm[m,i] = Real( T2C[m,i] - (-1)**m * T2C[-m,i] )
          Rp[m,i] = Real( T2C[m,i] + (-1)**m * T2C[-m,i] )
          Qm[m,i] = Imag(-T2C[m,i] + (-1)**m * T2C[-m,i] )
    """
    
    for i in range(2*l+1):
        ctg=None
        indexm = list(range(l+1))
        indexm = sorted(indexm, key=lambda m: abs(T2C[m+l,i])+abs(T2C[-m+l,i]), reverse=True)

        print('sorted by magnitude', file=fh_info)
        for m in indexm:
            print(m, T2C[m+l,i],T2C[-m+l,i], file=fh_info)
            
        for m in indexm[:1]:  # We will do it such that only the largest component is OK
            Qp = T2C[m+l,i].imag + (-1)**m * T2C[-m+l,i].imag
            Rm = T2C[m+l,i].real - (-1)**m * T2C[-m+l,i].real
            Rp = T2C[m+l,i].real + (-1)**m * T2C[-m+l,i].real
            Qm =-T2C[m+l,i].imag + (-1)**m * T2C[-m+l,i].imag
            #print 'Qp=', Qp, 'Rm=', Rm
            if abs(Qp) > small or abs(Rm) > small:
                if abs(Qp) > small :
                    ctg = -Rp/Qp
                    xb = -Rp
                    xa = Qp
                if abs(Rm) > small :
                    ctg = -Qm/Rm
                    xb = -Qm
                    xa = Rm
                #break
            
        if ctg is not None:
            phi = arctan2(xa, xb)                        
            print('Correcting T2C because original cubic harmonics were not real. Vector=', i, 'phase=', phi, file=fh_info)
            print('T2C^T before correction:', file=fh_info)
            cprint2(fh_info, T2C.transpose())
            T2C[:,i] = T2C[:,i] * exp(phi*1j)
            print('T2C^T after correction:', file=fh_info)
            cprint2(fh_info, T2C.transpose())
            
            for m in range(0,l+1):
                Rm = T2C[m+l,i].real - (-1)**m * T2C[-m+l,i].real
                Qp = T2C[m+l,i].imag + (-1)**m * T2C[-m+l,i].imag
                #Rp = T2C[m+l,i].real + (-1)**m * T2C[-m+l,i].real
                #Qm =-T2C[m+l,i].imag + (-1)**m * T2C[-m+l,i].imag
                #print 'm=', m, 'Rm=', Rm, 'Qp=', Qp
                
                if abs(Rm)>10*small or abs(Qp)>10*small:
                    print('ERROR: Could not find an angle to make all cubic harmonics real for vector=', i, 'and m=', m, 'D[Real]=', Rm, 'D[Imag]=', Qp)
                    print('ERROR: Could not find an angle to make all cubic harmonics real for vector=', i, 'and m=', m, 'D[Real]=', Rm, 'D[Imag]=', Qp, file=fh_info)
                    
            
def FindCoupledBaths(Sigind):
    # Finds which baths are coupled together through off-diagonal terms
    wcoupled=[]
    left_in = list(range(len(Sigind)))
    while left_in:
        ic=left_in[0]
        if Sigind[ic,ic]==0:
            left_in.remove(ic)
            continue
        coupled=[ic]
        for j in left_in:
            if j!=ic and Sigind[ic,j]!=0:
                coupled.append(j)
        for j in coupled: left_in.remove(j)
    
        wcoupled.append(coupled)
    return wcoupled


def FastIsing(Nrange,UC):
    def inside(i,lim):
        if i<=0 or i>lim:
            return 0
        else:
            return i

    print('Br:', 'Stage00: Fast Ising Computations')
    print('Stage00: Fast Ising Computation', file=fh_info)
    Eimpc = zeros( 2*(2*l+1) )
    for ic in range(len(Sigind)):
        Eimpc[ic] = Eimp[Sigind[ic,ic]-1]
    print('Eimpc=', Eimpc, file=fh_info)

    bath=[]
    for ic in range(len(Sigind)):
        if Eimpc[ic]<1000:
            if Q3d:
                bath.append(op.bi[ic])
            else:
                bath.append(ic)
    bath=array(bath,dtype=int)
    
    print('bath=', bath, file=fh_info)
    
    if l==3:
        global_flip=[0,1,2,2,1,0,3,4,5,6,6,5,4,3]
    else:
        if Q3d:
            #global_flip = range(2*l+1) + range(2*l+1)
            global_flip = bath
        else:
            #global_flip = [i/2 for i in range(2*(2*l+1))]
            global_flip = [i/2 for i in bath]
        
    
    wNrange = Nrange[:]
    if Nrange[0]!=0:
        wNrange = [Nrange[0]-1]+wNrange
    if Nrange[-1]!=(2*l+1)*2:
        wNrange = wNrange+[Nrange[-1]+1]
    print('START: wNrange=', wNrange, file=fh_info)

    temp_states=[]
    temp_occ=[]
    for i in range(Ntot):
        occ = array(op.occup(i))
        if sum(occ) in wNrange and sum(occ*Eimpc)<1000:
            #print i, occ, sum(occ*Eimpc)
            temp_states.append(i)
            temp_occ.append( occ )
    sstates=[]
    Nstates=[]
    for n in wNrange:
        for ind in range(len(temp_states)):
            if sum(temp_occ[ind])==n:
                sstates.append(temp_states[ind])
                Nstates.append(temp_occ[ind])
    Nlimits={}
    iis=iie=0
    for n in wNrange:
        while(iie<len(Nstates) and sum(Nstates[iie])<=n): iie+=1
        Nlimits[n] = (iis,iie)
        iis=iie

    print('Br:', 'Stage0: Exact diagonalization of the atom') 
    print('Stage0: Exact diagonalization of the atom', file=fh_info)

    FORTRAN=True
    maxk=l+1
    FkoJ = array(Fk[:maxk,l])
    Ene=zeros(len(sstates)) # Energy
    for js,st in enumerate(sstates):
        # on-site Energies in base of cubic harmonics, contain crystal-field splittings
        occ = array(Nstates[js])
        EUterms=zeros(2)
        if Q3d:
            if FORTRAN:
                EUterms = dpybind.FastIsing3d(Eimpc, occ, UC, FkoJ, bath)
            else:
                EUterms[0] += sum(Eimpc * occ)
                for i,ii in enumerate(bath):
                    si = int((i*2)/len(bath))
                    for j,jj in enumerate(bath):
                        sj = int((j*2)/len(bath))
                        dsum = sum(UC[:,ii,jj,jj,ii].real*FkoJ[:])
                        if (si==sj): dsum -= sum(UC[:,ii,jj,ii,jj].real*FkoJ[:])
                        EUterms[1] += 0.5*dsum * occ[i]*occ[j]
        else:
            if FORTRAN:
                EUterms = dpybind.FastIsing5d(Eimpc, occ, UC, FkoJ)
            else:
                EUterms[0] += sum(Eimpc * occ)
                for i in range(len(bath)):
                    for j in range(len(bath)):
                        dsum = sum((UC[:,i,j,j,i]-UC[:,i,j,i,j]).real*FkoJ[:])
                        EUterms[1] += 0.5*dsum * occ[i]*occ[j]
            
        Ene[js] = EUterms[0] + EUterms[1]
        print('E['+str(js+1)+']=%f ' % Ene[js], file=fh_info)

    print('Br:', 'Stage1: Computing F^ in direct base')
    print('Stage1: Computing F^ in direct base', file=fh_info)

    Fbp=[]
    for js in range(len(sstates)):
        Fbp.append( [[-1,0.0] for ib in range(len(bath))] )
    
    for js,st in enumerate(sstates):
        Nfinal = sum(Nstates[js])+1
        if Nfinal not in Nlimits: continue
        (ist_start, ist_end) = Nlimits[Nfinal]
        for ib in range(len(bath)):
            (newst, sig) = op.Fp(st, ib)
            if newst>0:
                ii = sstates.index(newst)
                Fbp[js][ib] = (ii,float(sig))
    
    #print 'Fbp=', Fbp
    i0=0
    while (i0<len(sstates) and (sum(Nstates[i0]) not in Nrange)): i0+=1
    i1=len(sstates)-1
    while (i1>=0 and (sum(Nstates[i1]) not in Nrange)): i1-=1;

    fcix = open('actqmc.cix', 'w')
    # ---------------------------------------------------------------------------------------
    # -------------- Below is printing for ctqmc  solver ------------------------------------
    # ---------------------------------------------------------------------------------------
    print('# CIX file for ctqmc! ', file=fcix)
    print('# cluster_size, number of states, number of baths, maximum_matrix_size', file=fcix)
    print(1, i1-i0+1, len(bath), 1, file=fcix)
    print('# baths, dimension, symmetry', file=fcix)
    for ib in range(len(bath)):
        print(ib, '  ', 1, Sigind[ib,ib]-1, '  ', global_flip[ib], file=fcix)
    print('# cluster energies for non-equivalent baths, eps[k]', file=fcix)
    #for E in Eimp: print(E, end=' ', file=fcix)
    print(('{:13.3f} '*len(Eimp)).format(*Eimp), file=fcix)
    print(file=fcix)
    print('#     N  K  Sz    size', file=fcix)
    for i in range(i0,i1+1):
        gs = sstates[i]
        if Q3d:
            Mz = op.Sz(gs)/2.
        else:
            Mz = op.Mz(gs)
        print("%3d  %2d %2d %6.3f %2d " % (i-i0+1, sum(Nstates[i]), 0, 2*Mz, 1), end=' ', file=fcix)
        for ib in range(len(bath)):
            ifinal,m = Fbp[i][ib]
            Nfinal = sum(Nstates[ifinal])
            if Nfinal not in Nrange: ifinal=-1
            if ifinal>=0:
                print("%3d" % (ifinal-i0+1), end=' ', file=fcix)
            else:
                print("%3d" % (0), end=' ', file=fcix)
        S2=0.0
        print("  %12.8f  %3.1f  # %s" % (Ene[i],S2, op.printn(gs)), file=fcix)
        
    print('# matrix elements', file=fcix)
    for i in range(i0,i1+1):
        for ib in range(len(bath)):
            ifinal,m = Fbp[i][ib]
            Nfinal = sum(Nstates[ifinal])
            if Nfinal not in Nrange: ifinal=-1
            if ifinal>=0: 
                print("%3d %3d  %2d %2d %f" % (i-i0+1, ifinal-i0+1,1,1,m), file=fcix)
            else:
                print("%3d %3d  %2d %2d" % (i-i0+1, 0, 0, 0), file=fcix)
    
    if HB2 : print('HB2', file=fcix)
    else: print('HB1', file=fcix)


    if (HB2 and Q3d):
        ii=0
        iind={}
        for i1,bs1 in enumerate(baths):
            m1 = bs1[0]
            s1 = bs1[1]
            if m1 not in bkeep: continue
            iind[i1]=ii
            ii+=1
        print("# UCoulomb : (m1,s1) (m2,s2) (m3,s2) (m4,s1)  Uc[m1,m2,m3,m4]", file=fcix)
        for i1,bs1 in enumerate(baths):
            m1 = bs1[0]
            s1 = bs1[1]
            if m1 not in bkeep: continue
            for i2,bs2 in enumerate(baths):
                m2 = bs2[0]
                s2 = bs2[1]
                if m2 not in bkeep: continue
                for i3,bs3 in enumerate(baths):
                    m3 = bs3[0]
                    s3 = bs3[1]
                    if (s2!=s3): continue
                    if m3 not in bkeep: continue
                    for i4,bs4 in enumerate(baths):
                        m4 = bs4[0]
                        s4 = bs4[1]
                        if (s4!=s1): continue
                        if m4 not in bkeep: continue
                        Uc = 0.0
                        for k in range(0,l+1):
                            Uc += real(UC[k,m1,m2,m3,m4])*Fk[k,l]
                        if abs(Uc)>1e-6:
                            print("%2d %2d %2d %2d  %12.8f" % (iind[i1],iind[i2],iind[i3],iind[i4],Uc), file=fcix)

    if (HB2 and not Q3d):
        print("# UCoulomb : (m1,s1) (m2,s2) (m3,s2) (m4,s1)  Uc[m1,m2,m3,m4]", file=fcix)
        for bs1 in bkeep:
            for bs2 in bkeep:
                for bs3 in bkeep:
                    for bs4 in bkeep:
                        Uc = 0.0
                        for k in range(0,l+1):
                            UC += real(UC[k,bs1,bs2,bs3,bs4])*Fk[k,l]
                        if abs(Uc)>1e-6:
                            print("%2d %2d %2d %2d  %12.8f" % (bs1,bs2,bs3,bs4,Uc), file=fcix)
    
                
    print('# number of operators needed', file=fcix)
    print('1', file=fcix)
    print('# Occupancy ', file=fcix)
    
    bathi=[Sigind[b,b]-1 for b in bath]
    for i in range(i0,i1+1):
        nb=0.0
        nb = zeros(len(set(bathi)))
        for ib in range(len(bath)):
            nb[bathi[ib]] += Nstates[i][ib]
        for ib in range(len(nb)):
            print("%3d  %2d %2d %f" % (i-i0+1, 1, 1, nb[ib]), file=fcix)
    
    print('# Data for HB1', file=fcix)
    
    print(1, len(sstates), len(bath), 1, file=fcix)
    print('#      ind N  K  Jz     size', file=fcix)
    for i in range(len(sstates)):
        gs = sstates[i]
        if Q3d:
            Mz = op.Sz(gs)/2.
        else:
            Mz = op.Mz(gs)
        print("%3d  %3d  %2d %2d %6.3f %2d " % (i+1, inside(i-i0+1,i1+1), sum(Nstates[i]), 0, 2*Mz, 1), end=' ', file=fcix)
        for ib in range(len(bath)):
            ifinal,m = Fbp[i][ib]
            print("%3d" % (ifinal+1), end=' ', file=fcix)
        print("  ", end=' ', file=fcix)
        print("%12.8f  %3.1f  # %s" % (Ene[i], 0, op.printn(gs)), file=fcix)
    print('# matrix elements', file=fcix)
    for i in range(len(sstates)):
        for ib in range(len(bath)):
            ifinal,m = Fbp[i][ib]
            print("%3d %3d " % (i+1, ifinal+1), end=' ', file=fcix) 
            if ifinal>=0:
                print("%2d %2d" % (1,1),  m, file=fcix)
            else:
                print("%2d %2d" % (0, 0), file=fcix)
            
if __name__ == '__main__':
    """ Help here"""
    #n=[1,2,3]  # occupanices used for OCA
    l=1        # angular momentum
    J = 0.3    # Hunds coupling
    #qOCA=1     # OCA diagrams are computed 
    Eoca=10.   # Energy window for OCA diagrams
    #mOCA=1e-3  # matrix element for OCA should be greater than that
    #Ncentral=[5] # OCA diagrams are selected such that central occupancy is in Ncentral
    Ewindow = [-1000,1000]
    max_M_size=500
    add_occupancy=True
    CoulombF = 'Full' #'Bulla_Jarrel' #'Full' # 'Bulla_Jarrel' # 'Full' 'Oles'
    #OCA_G=False
    PrintReal=True
    HB2 = True
    Eimp = [0.0]
    Nmax = 1024
    small = 1e-6
    small_t2c=1e-5
    Nrange=[]
    PrintSminus=False
    ORB=[] #ORB=[0,0,1,-1,0,0,0,1,-1,0] # [z^2,x^2-y^2,xz,yz,xy,...]
    OLD_CTQMC=False

    # l,J, Ewindow=[-1000,1000], max_M_size=500, ORB=[]
    # removed: qOCA, mOCA, n, Eoca, OCA_G, 
    args = sys.argv[1:]
    if ('-h' in args) or ('--help' in args):
        print("""Code for generating impurity cix file for a d-type of material
                 The output cix-file can be used for OCA or CTQMC solvers
                 The outputs:
                    out.cix       -- the oca cix file
                    impurity.cix  -- the ctqmc cix file
                 The input is:
                    filename      -- another python script, which contains some definitions of parameters
                    Eimp          -- list of impurity levels (i.e., [0,0,0...] )
                    Sigind        -- The symmetry of the self-energy and impurity levels given as a 2D-list
                    CF            -- The local rotation matrix given as a 2D list or array
                    l             -- angular momentum (l=0, l=1, l=2 supported)
                    J             -- Hund's coupling
                    Ewindow       -- Energy window for the states kept (used in ctqmc only)
                    max_M_size    -- maximum matrix size kept for ctqmc
                    ORB           -- index for orbtal susceptibility
                 """)
        sys.exit(0)
    Q3d=None
    for arg in args:
        if os.path.isfile(arg):
            exec(compile(open(arg, "rb").read(), arg, 'exec'))
            print('Executed file', arg)
        else:
            exec(arg)

    if CoulombF=='Georges': CoulombF='Ising'
            
    fh_info = open('info_atom_d.dat','w')
    print(' '.join(sys.argv), file=fh_info)
    print('Eimp=', '[','%f, '*len(Eimp) % tuple(Eimp),']', file=fh_info)
    #print('n=', n, file=fh_info)
    print('l=', l, file=fh_info)
    print('J=', J, file=fh_info)
    #print('qOCA=', qOCA, file=fh_info)
    #print('Eoca=', Eoca, file=fh_info)
    #print('mOCA=', mOCA, file=fh_info)
    print('Ewindow=', Ewindow, file=fh_info)
    print('max_M_size=', max_M_size, file=fh_info)
    #print('para=', para, file=fh_info)
    print('ORB=', ORB, file=fh_info)
    
    ftrans='Trans.dat'
    #if (len(glob.glob(ftrans))!=0):
    if os.path.exists(ftrans):
        print('Reading file', ftrans, file=fh_info)
        (Sigind, CF) = ReadTrans(ftrans, fh_info)
        if len(Sigind)==(2*l+1):
            dn = 2*l+1
            SigindN = zeros((2*dn, 2*dn), dtype=int)
            SigindN[:dn,:dn] = Sigind
            SigindN[dn:,dn:] = Sigind
            Sigind=SigindN
        
        if (len(Sigind)!= 2*(2*l+1)):
            print('ERROR: Sigind from file', ftrans, 'does not have correct dimension!')
            sys.exit(1)

        if len(CF)==2*l+1:
            if Q3d==None: Q3d=True # this must be 3d orbital
        elif len(CF)==2*(2*l+1):
            dn=2*l+1
            off_diag = CF[dn:,:dn]
            off = sum(sum(abs(off_diag)))
            print(off)
            if Q3d==None:
                if off>1e-5:
                    Q3d=False
                else:
                    Q3d=True
        else:
            print('ERROR: Transformation CF=T2C does not have correct dimension')
            sys.exit(1)

        # If any diagonal entry in Sigind[] is equal to zero, we want to project it out.
        # This is achieved by making impurity levels Eimp very large for this orbital
        Sigind_orig=copy.copy(Sigind)
        Simax = max(list(map(max,Sigind)))
        for i in range(len(Sigind)):
            if Sigind[i,i]==0:
                Sigind[i,i]=Simax+1
                if len(Eimp)<=Simax : Eimp += [2000.]
                
        
        # Matrix contained in Trans.dat should be T(i,m):
        #   - i runs over real harmonics (z^2, x^2-y^2, yz, xz, xy)
        #   - m runs over complex harmonics (-2, -1, 0, 1, 2)
        # The matrix in Trans.dat, read into CF, is the transpose of
        # what we want in T2C, hence the T2C = transpose(CF) statement
        if Q3d:
            ##### change 2013
            T2C = transpose(CF[:(2*l+1),:(2*l+1)])
            Check_T2C_Real(T2C, l, fh_info, small_t2c)
        else:
            #CFN = zeros(shape(CF),dtype=complex)
            #CFN[:,:(2*l+1)] = CF[:,(2*l+1):]
            #CFN[:,(2*l+1):] = CF[:,:(2*l+1)]
            CFN = CF
            T2C = transpose(CFN)
        ##### change 2013, you will need to generalize this
    else:
        """
        Cubic harmonics have the order compatible with the Wien2K package:
        z^2, x^2-y^2, yz, zx, xy
        """
        print(ftrans, 'file not found; generating Sigind & complex-to-real spherical harmonics transformation inside atom_d.')
        T2C = transpose(Spheric2Cubic(l))
        if Q3d==None: Q3d=True

        Sigind = zeros((2*(2*l+1),2*(2*l+1)), dtype=int)
        if Q3d:
            for i in range(2*l+1):
                Sigind[i,i] = i+1
                Sigind[i+(2*l+1),i+(2*l+1)] = i+1
        else:
            zr = zeros(2*l+1)
            CF = transpose(T2C)
            CFn=[]
            for i in range(2*l+1):
                CFn.append( CF[i].tolist()+zr.tolist()  )
                CFn.append( zr.tolist()+CF[i].tolist()  )
                Sigind[2*i,2*i] = i+1
                Sigind[2*i+1,2*i+1] = i+1
            CFn = array(CFn)
            T2C = transpose(CFn)
            
    if Q3d:
        global_flip = list(range(2*l+1)) + list(range(2*l+1))
    else:
        global_flip=[]
        for i in range(2*l+1):
            global_flip += [i,i]
            
    print('Sigind=', Sigind, file=fh_info)
    
    print('T2C=', file=fh_info)
    for i in range(len(T2C)):
        for j in range(len(T2C)):
            print("%6.3f %6.3f   " % (T2C[i,j].real, T2C[i,j].imag), end=' ', file=fh_info)
        print(file=fh_info)
    
    print('global_flip=', global_flip, file=fh_info)
    
    print('Impurity level structure Sigind is:', file=fh_info)
    print(Sigind, file=fh_info)

    print('T2C^T follows:', file=fh_info)
    print('\n'.join('   '.join('%10f %10f' % (x.real, x.imag) for x in row) for row in T2C.transpose()), file=fh_info)
    print('shape(T2C)=', shape(T2C), file=fh_info)
    print('T2C is Unitary=', sum(abs(matrix(T2C) * matrix(T2C).H - identity(len(T2C)))), file=fh_info)
    
    Nitt=1  # To figure out the symmetry, we itterate Nitt times
    
    Jc = J
    cx = 0. # No spin orbit at the moment

    # Ratio between F2,F4,F6 and J! At the end of the day, we want to use U and J only!
    #Fk = gaunt.slaterf(1., Jc)
    U0=1.
    Fk = SlaterF(U0, Jc, l)
    print('Slater integrals F^k are ', Fk[:,l], file=fh_info)

    # one electron base
    baths=[]
    for s in [1,-1]:
        for b in range(2*l+1):
            baths.append([b,s])
    bathi=[Sigind[b,b]-1 for b in range(len(Sigind))]
    
    dkbth={}
    for i in range(len(bathi)):
        if bathi[i] in dkbth:
            dkbth[bathi[i]].append(i)
        else:
            dkbth[bathi[i]]=[i]
    kbth=[]
    for k in sort(list(dkbth.keys())):
        kbth.append(dkbth[k])
        
    kbth0=[]
    for i in range(len(baths)): kbth0.append([i])
    
    #print max(bathi)
    if max(bathi)>=len(Eimp):
        print('You specified wrong dimension for Eimp! Boiling out.... Need at least '+str(max(bathi)+1)+' components')
        sys.exit(1)
    
    bkeep=[]
    for b in range(len(bathi)):
        if Eimp[bathi[b]]<1000:  bkeep.append(b)

    tEimp = [x for x in Eimp if x<1000]
    tkbth=[]
    for k in kbth:
        if k[0] in bkeep: tkbth.append(k)
    
    print('Some other info in ED:', file=fh_info)
    print('bathi=', bathi, file=fh_info)
    print('kbth=', kbth, file=fh_info)
    print('tkbth=', tkbth, file=fh_info)
    print('Eimp=', Eimp, file=fh_info)
    print('tEimp=', tEimp, file=fh_info)
    print('bkeep=', bkeep, file=fh_info)
    
    Ntot = 2**(len(baths)) # size of the direct base
    
    op = operateLS(2*(2*l+1), T2C, Q3d) # operators class
    
    if op.Q3d:
        print('baths bi=', op.bi, file=fh_info)
        print('spin Sz=', op.sz, file=fh_info)
        print('mask-down=', op.mask_d, file=fh_info)
        print('mask-up  =', op.mask_u, file=fh_info)
    
    SO = SpinOrbitM(l,T2C) # Spin orbit matrix


    UC = CoulUsC2(l,T2C)
    #_UC_ = CoulUsC2_slow(l,T2C)
    #print('are close=', allclose(UC,_UC_))
    
    # These are possible values for CoulombF:
    #   'Ising'  : density-density interaction of Slater type
    #   'Full'   : Fully rotationally invariant interaction of Slater type, but with SU(N) symmetry [like in Kanamori]
    #   'FullS'  : Fully rotationally invariant interaction of Slater type, but with SU(3) symmetry [usual Slater interaction]
    #   'FullK'  : Kanamori type interaction, but U and J are computed so that within t2g block the interaction has the same form as Slater interaction
    if CoulombF in ['Full','Ising','FullS', 'IsingS']:      # Slater parametrization
        UC = CoulUsC2(l,T2C)              # Coulomb repulsion matrix
        Ucn = zeros(shape(UC),dtype=complex)
        if CoulombF in ['Full', 'IsingS']:  # We keep matrix elements which give no sign problem
            for i in range(2*l+1):
                for j in range(2*l+1):
                    Ucn[0,i,j,j,i] = 1.0
                    Ucn[1,i,j,j,i] = UC[1,i,j,j,i]
                    Ucn[2,i,j,j,i] = UC[2,i,j,j,i]
                    if i!=j:
                        Ucn[1,i,j,i,j] = UC[1,i,j,i,j]
                        Ucn[1,i,i,j,j] = UC[1,i,i,j,j]
                        Ucn[2,i,j,i,j] = UC[2,i,j,i,j]
                        Ucn[2,i,i,j,j] = UC[2,i,i,j,j]
            UC = Ucn
    elif CoulombF in ['FullK','IsingK']:  # Kanamori parametrization
        UC = zeros((l+1, 2*l+1, 2*l+1, 2*l+1, 2*l+1), dtype=complex)
        for i in range(5):
            for j in range(5):
                if i==j: 
                    UC[0,i,j,j,i]=1.
                    UC[1,i,j,j,i]=4./49.
                    UC[2,i,j,j,i]=4./49.
                else:
                    UC[0,i,j,j,i]=1.      # F0 
                    UC[1,i,j,j,i]=-2./49. # F2 exact for t2g's
                    UC[2,i,j,j,i]=-4./441.# F4 exact for t2g's
                    UC[1,i,j,i,j]= 3./49. # F2 exact for t2g's
                    UC[2,i,j,i,j]=20./441.# F4 exact for t2g's
                    UC[1,i,i,j,j]= 3./49. # F2 exact for t2g's
                    UC[2,i,i,j,j]=20./441.# F4 exact for t2g's
    else:
        print('ERROR: CoulombF form ', CoulombF, ' not yet implemented!')
        sys.exit(0)
    
    if os.path.isfile('../Uc.dat') and os.path.getsize('../Uc.dat')>0:
        Uc = loadtxt('../Uc.dat')
        for m1 in range(5):
            for m2 in range(5):
                UC[0,m1,m2,m2,m1] = Uc[m1,m2]
                #if abs(UC[0,m1,m2,m2,m1])>1e-3:
                #    print "%2d %2d %2d %2d  %5.2f " % (m1, m2, m2, m1, UC[0,m1,m2,m2,m1])
    else:
        UC[0,:,:,:,:]=0.0
        
    if Q3d:
        nw=2*l+1
    else:
        nw=2*(2*l+1)
    
    if Nrange:
        FastIsing(Nrange)
        sys.exit(0)
    else:
        Nrange = list(range((2*l+1)*2+1))
        
    # some properties of integers which will serve a direct base - partial occupancy and Sz
    prop=[]
    for i in range(Ntot):
        #### 2013 ### all these are wrong for 5d
        occ = op.occup(i)
        prop.append([sum(occ), op.Sz(i),occ])
    # creates direct base from integers having correct properties
    # wstates contains [N, Sz, [all direct states with this N and Sz]]

    if not Nrange: Nrange = list(range((2*l+1)*2+1))
    wstates = baseN(Nrange,prop,op.Q3d)
    
    indx={}
    for ni,ns in enumerate(wstates):
        indx[(ns[0],ns[1])] = ni  # index according to N and Sz of the state

    print('indx=', file=fh_info)
    print(indx, file=fh_info)
    kindx = list(indx.keys())
    
    print('Stage0: Exact diagonalization of the atom', file=fh_info)

    mxs = max(list(map(max,Sigind)))
    if len(Eimp)<mxs:
        print('ERROR: The dimension of the Eimp should be equal to the maximum index of Sigind->',mxs)
        sys.exit(1)
        
    Eimpc = zeros((2*(2*l+1), 2*(2*l+1)), dtype=complex)
    for ic in range(len(Sigind)):
        Eimpc[ic,ic] = Eimp[Sigind[ic,ic]-1]

    print('impurity levels Eimpc0=', file=fh_info)
    mprint(fh_info, real(Eimpc))
    #Eimpc = matrix(T2C) * matrix(Eimpc) * matrix(T2C).H
    #print 'Eimpc1='
    #mprint(Eimpc)

    
    Ene=[] # Energy
    Te=[]  # eigenvectors
    S2w=[] # Spin
    #### NEW
    if PrintSminus:
        Si_minus={}
        Sm_minus={}
    #### NEW
    for ni,ns in enumerate(wstates):

        #print 'Br:', 'n=', ns[0], 'sz=', ns[1]/2.
        
        print('----------------------------------------------------', file=fh_info)
        print('n=', ns[0], 'sz=', ns[1]/2., file=fh_info)
        states = ns[2]
        # printing all states in this sector
        print('states=', end=' ', file=fh_info)
        for ist,st in enumerate(states): print(('%d:'%ist),op.printn(st), end=' ', file=fh_info)
        print(file=fh_info)
        
        S2 = zeros((len(states),len(states)),dtype=complex)
        if CoulombF[:4] == 'Full' and op.Q3d:
            # Computes matrix of S^2
            for js,st in enumerate(states):
                #print st, op.printn(st), "    ",
                stn = op.S2(st)
                #print stn
                for ps in stn:
                    ii = ps[0]
                    iu = states.index(ii)
                    S2[js,iu] += ps[1]
        ##### NEW
        if PrintSminus:
            j_minus=-1
            for j in range(len(wstates)):
                if (wstates[j][0]==ns[0] and wstates[j][1]==ns[1]-2):
                    j_minus = j
                    break
            if j_minus>=0:
                n_new = len(wstates[j_minus][2])
                S_minus = zeros((len(states),n_new),dtype=complex)
                #print 'n=', ns[0], 'sz=', ns[1]/2., 'j_minus=', j_minus, 'n=', wstates[j_minus][0], 'sz=', wstates[j_minus][1]/2., 'states=', wstates[j_minus][2]
                for js,st in enumerate(states):
                    stn = op.S_minus(st)
                    #print 'stn=', stn
                    for ps in stn:
                        ii = ps[0]
                        states_plus = wstates[j_minus][2]
                        iu = states_plus.index(ii)
                        S_minus[js,iu] += ps[1]
        ###### NEW

        
        Ham = zeros((len(states),len(states)),dtype=complex)

        for js,st in enumerate(states):
            # on-site Energies in base of cubic harmonics
            # contain crystal-field splittings
            DM = op.DM(st)
            for i in range(len(Eimpc)):
                for j in range(len(Eimpc)):
                    if abs(Eimpc[i,j])<1e-5:continue
                    for p in DM[i][j]:
                        iu = states.index(p[0])
                        Ham[js,iu] += p[1]*Eimpc[i,j]
            
            if CoulombF[:4] == 'Full':
                ## Coulomb interaction including F2 and F4
                cst = op.CoulombU(st, UC, Fk[:,l])
                for cs in cst:
                    ii = cs[0]
                    U0 = cs[1][0]
                    U1 = cs[1][1]
                    iu = states.index(ii)
                    Ham[js,iu] +=  0.5*U1 # adding only F2,F4,... but not F0
                    Ham[js,iu] +=  0.5*U0 # adding only F2,F4,... but not F0
                    
            elif CoulombF[:5] == 'Ising':
                ## Coulomb interaction including F2 and F4
                cst = op.CoulombUIsing(st, UC, Fk[:,l])
                for cs in cst:
                    ii = cs[0]
                    U0 = cs[1][0]
                    U1 = cs[1][1]
                    #iu = states.index(ii)
                    #if (iu!=js): print 'ERROR', iu, js, 'iu==js', iu==js
                    
                    Ham[js,js] +=  0.5*U1 # adding only F2,F4,... but not F0
                    Ham[js,js] +=  0.5*U0 # adding only F2,F4,... but not F0
                    
            else:
                print('Not yet implemented 2!')
                sys.exit(1)
                
            # Spin-orbit interaction
            if cx>1e-5:
                cst = op.OneBodyNab(st, SO)
                for cs in cst:
                    iu = states.index(cs[0])
                    #if (js>=len(states) or iu>=len(states)): print 'Tezave!'
                    Ham[js,iu] += cs[1]*cx
            
        #if (cx>1e-5): cprint(Ham)
        #else:
        #print >> fh_info, 'H='
        #mprint(fh_info, Ham)
        
        if CoulombF[:4] == 'Full':
            eig = Diagonalize(Ham, small, fh_info)  # Block diagonalization is better!!
            Ex = eig[0]
            Tx = eig[1]
            Ene.append( Ex )
            Te.append( Tx )
        else:
            Ex = [real(Ham[i,i]) for i in range(len(Ham))]
            Tx = eye(len(Ham),len(Ham))
            Ene.append( Ex )
            Te.append( Tx )


        if CoulombF[:4]=='Full' and op.Q3d:
            # Here we compute matrix of S^2 in eigenbase. Should be diagonal if no spin-orbit coupling
            S2e = matrix(conj(Tx.transpose())) * S2 * matrix(Tx)
            
            printS=False
            trouble=[]
            for i0 in range(shape(S2e)[0]):
                for i1 in range(shape(S2e)[1]):
                    if i0!=i1 and abs(S2e[i0,i1])>1e-6 :
                        print('WARNING: Troubles->', i0, i1, S2e[i0,i1])
                        printS=True
                        trouble.append(i0)
            
            printS=False # BRISISSS
            if printS:
                print('S2=', file=fh_info)
                mprint(fh_info, S2e)
                print('H=', file=fh_info)
                cprint(fh_info, Ham)
                for it,t in enumerate(trouble):
                    print('A[%d]=' % t, file=fh_info)
                    print(Tx[t], file=fh_info)
            
            S2w.append([0.5*int(round(-1+sqrt(1+4*S2e[i,i].real))) for i in range(len(S2e))]) # Spin is computed using formula s(s+1)
        else:
            S2w.append( [0 for i in range(len(S2))] )
        


        ##### NEW
        if PrintSminus:
            if (j_minus>=0):
                #print 'len(Te)=', len(Te), 'j_minus=', j_minus
                _S_minus_ = matrix(conj(Tx.transpose())) * S_minus * matrix(Te[j_minus])
                Si_minus[ni]=j_minus
                Sm_minus[ni]= real(_S_minus_)
        ##### NEW
            
        #print >> fh_info, 'E=', '%f '*len(Ex) % tuple(Ex)
        print('state, E, |largest_coeff|', file=fh_info)
        for ie in range(len(Ex)):
            ilargest = argmax(abs(Tx[:,ie]))
            print(op.printn(states[ilargest]),  Ex[ie], abs(Tx[ilargest,ie]), file=fh_info)


    #print 'kindx=', kindx
    #print 'wstates=', wstates

    print('Br:', 'kindx=', kindx)

    # Here we create index for psi^dagger
    iFi = zeros((len(wstates),len(baths)),dtype=int)
    for ni,ns in enumerate(wstates):
        for ib,be in enumerate(baths):
            if op.Q3d:
                st = (ns[0]+1, ns[1] + be[1])  # (n+1,sz+s)
                if st in kindx:
                    iFi[ni,ib] = indx[st]
            else:
                st = (ns[0]+1, 0)  # (n+1,sz+s)
                if st in kindx:
                    iFi[ni,ib] = indx[st]
    
    wgr = [[] for iw in range(len(wstates))]
    print('Br:', 'Stage1: Computing F^ in direct base')
    print('Stage1: Computing F^ in direct base', file=fh_info)
    
    # Below we compute matrix elements of F^ (FKP)
    kindx = list(indx.keys())
    FKP = []
    for ni,ns in enumerate(wstates):
        states = ns[2]
        bFp=[]

        if CoulombF[:5] == 'Ising':
            wgr[ni] += [[ist] for ist in range(len(states))]
            
        for ib,wib in enumerate(baths):
            inew = iFi[ni,ib]
            if inew==0:
                bFp.append([])
                continue
            
            newstates = wstates[inew][2]
            
            Fp = zeros((len(states), len(newstates)), dtype=complex)
            
            for js,st in enumerate(states):
                (newst, sig) = op.Fp(st, ib)
                if newst>0:
                    ii = newstates.index(newst)
                    Fp[js,ii] += sig
                #print 'state=', st, newst, ii

            if CoulombF[:5] == 'Ising':
                bFp.append(Fp)
            else:
                Fpn = Te[ni].conj().T @ Fp @ Te[inew]

                # Set to zero small values
                for i0 in range(shape(Fpn)[0]):
                    Fpn[i0, abs(Fpn[i0,:])<small ] = 0.0
                # slower bu more readable way
                #for i0 in range(shape(Fpn)[0]):
                #    for i1 in range(shape(Fpn)[1]):
                #        if abs(Fpn[i0,i1])<small: Fpn[i0,i1]=0.0
                            
                    
                gr0,gr1 = analizeGroups(Fpn, small)
                
                # |inew> = F^+ |ni>
                wgr[ni]   += gr0 # which states are coupled by F^ in |ni>
                wgr[inew] += gr1 # which states are coupled by F^ in |inew>
                bFp.append(Fpn)
            
        FKP.append(bFp)

    #FKP created!
    
    print('Br:', 'Stage2: Compressing F^+ according to its block diagonal form')
    print('Stage2: Compressing F^+ according to its block diagonal form', file=fh_info)
    
    for i in range(len(wstates)):
        wgr[i] = compress(wgr[i])
        print(i+1, wgr[i], file=fh_info)
    
    print('Stage3: Renumbers states -- creates superstates for ctqmc', file=fh_info)
    print('Br:', 'Stage3: Renumbers states -- creates superstates for ctqmc')

    # Here we store ground state and N for each superstates to be used later for sorting
    tstates=[]
    for iw in range(len(wstates)):
        Nt = sum(wstates[iw][0])
        for ip in range(len(wgr[iw])):
            Eg = Ene[iw][wgr[iw][ip][0]]
            tstates.append([iw,ip,Nt,Eg])

    # Uncomment this if you want to sort according to energy
    #tstates.sort(comp)
    # Here we sort only by particle number
    #tstates.sort(lambda a,b: cmp(a[2],b[2]))

    # tstates contains [index-to-wstates, index-to-state-inside-wstates, N, E]
    
    # superstates == pseudostates are defined
    pseudostates=[]
    indpseudo={}
    jj=0
    for st in tstates:
        iw = st[0]
        ip = st[1]
        pseudostates.append([iw,ip])
        indpseudo[(iw,ip)] = jj
        jj+=1

    iFi_inside=[]
    for iw in range(len(wstates)):
        biFi=[]
        for ib in range(len(baths)):
            ifi = iFi[iw,ib]
            if ifi>0:
                fpair = coupled(FKP[iw][ib], wgr[iw], wgr[ifi], small)
                biFi.append(fpair)
            else:
                biFi.append([])
        iFi_inside.append(biFi)
        
    #print 'iFi_inside=', iFi_inside
    
    # creates arrays containing Energy, occupancy and index table for all superstates
    iFinal = zeros((len(pseudostates),len(baths)),dtype=int)
    Enes = []
    S2ws = []
    Occ=[]
    rS_minus={}
    iS_minus={}
    for ii,iwp in enumerate(pseudostates):
        iw = iwp[0]
        ip = iwp[1]
        wstate = wstates[iw]
        group = wgr[iw][ip]
    
        for ib in range(len(baths)):
            ifi = iFi[iw,ib]
            ifinal = -1
            if (ifi>0):
                ifi_ins = iFi_inside[iw][ib][ip]
                if ifi_ins>=0:
                    ifinal = indpseudo[(ifi,ifi_ins)]
            iFinal[ii,ib] = ifinal

        Ens=[]
        occ=[]
        S2s=[]
        for iq,q in enumerate(group):
            Ens.append(Ene[iw][q])
            occ.append(wstate[0])
            S2s.append(S2w[iw][q])

            
        Enes.append(Ens)
        Occ.append(occ)
        S2ws.append(S2s)
        ### NEW
        if PrintSminus:
            if iw in Si_minus:
                _N_ = wstate[0]
                _Sz_ = wstate[1]
                for jj,jwp in enumerate(pseudostates):
                    jw = jwp[0]
                    jp = jwp[1]
                    if (wstates[jw][0]==_N_ and wstates[jw][1]==_Sz_-2):
                        if jw!=Si_minus[iw]:
                            print('ERROR in S_minus. The two indices are different jw=', jw, 'Si=', Si_minus[iw])
                        group_final = wgr[jw][jp]
                        SM_minus = zeros((len(group),len(group_final)))
                        for iig,ig in enumerate(group):
                            for jjg,jg in enumerate(group_final):
                                SM_minus[iig,jjg] = Sm_minus[iw][ig,jg]
                        if sum(abs(SM_minus))>1e-7:
                            #print 'possible S- between ii+=', ii+1, 'and jj+=', jj+1, 'where jw=', jw, 'Si=', Si_minus[iw], 'wgr=', wgr[jw][jp], 'shape(S_minus)=', shape(Sm_minus[iw])
                            #mprint(sys.stdout, SM_minus)
                            rS_minus[ii] = SM_minus
                            iS_minus[ii] = jj
        ### NEW
    
    #print 'pseudostates=', pseudostates
    #print 'Enes='
    #for ii in range(len(Enes)):
    #    print 'ii=', Enes[ii][0]
    #print 'Occ=', Occ
    #print 'S2=', S2ws

    print('Stage4: F^dagger matrices between superstates evaluated', file=fh_info)
    print('Br:', 'Stage4: F^dagger matrices between superstates evaluated')
    
    # creates small F^dagger matrices between superstates
    maxs = 0
    rFKP = []
    rNn=[]  # This is the matrix F*F^ == 1-N
    
    for ii,iwp in enumerate(pseudostates):
        iw = iwp[0]
        ip = iwp[1]
        wstate = wstates[iw]
        group = wgr[iw][ip]
    
        bNn=[]
        bFKP=[]
        for ib in range(len(baths)):
            ifi = iFi[iw,ib]
            ifinal = iFinal[ii,ib]
            if (ifi>0): ifi_ins = iFi_inside[iw][ib][ip]
    
            if ifinal>=0:
                M = zeros((len(group),len(wgr[ifi][ifi_ins])),dtype=complex)
                for ii0,i0 in enumerate(group):
                    for jj0,j0 in enumerate(wgr[ifi][ifi_ins]):
                        M[ii0,jj0] = FKP[iw][ib][i0,j0]
                        
                if max(shape(M)) > maxs : maxs = max(shape(M))
    
                Nn = zeros((len(group),len(group)))
                for ii0,i0 in enumerate(group):
                    for ii1,i1 in enumerate(group):
                        Nn[ii0,ii1]=sum(M[ii0]*M[ii1]).real
                        
            else:
                M = array([])
                Nn = array([])
            bFKP.append(M)
            bNn.append(Nn)
        rFKP.append(bFKP)
        rNn.append(bNn)


    # Extract low energy states
    lowE=[]
    low_maxsize=0
    for ii in range(len(Enes)):
        size=0
        plowE=[]
        for iq in range(min(len(Enes[ii]),max_M_size)):
            if Enes[ii][iq]>=Ewindow[0] and Enes[ii][iq]<Ewindow[1] and iq<Nmax:
                plowE.append(iq)
                size += 1
        if size>low_maxsize: low_maxsize = size
        if len(plowE)>0:
            lowE.append((ii,plowE))

    # Creates index array between all states and low energy ones
    inv_lowE1={-1:-1}
    for i in range(len(pseudostates)): inv_lowE1[i]=-1
    for i in range(len(lowE)):
        ii = lowE[i][0]
        inv_lowE1[ii]=i

    #print 'lowE=', lowE
    #print 'inv_lowE1=', inv_lowE1

    #print 'Enes='
    #for ii in range(len(Enes)):
    #    print 'ii=', Enes[ii][0]


    nbaths = len(bkeep)
    
    wcoupled = FindCoupledBaths(Sigind_orig)
    print('wcoupled=', wcoupled, file=fh_info)
    DiagonalOnly=True
    if sum(list(map(len,wcoupled)) - ones(len(wcoupled)))!=0:
        DiagonalOnly=False
        nbaths = len(wcoupled)
    print('Diagonal=', DiagonalOnly, file=fh_info)
    
    fcix = open('actqmc.cix', 'w')
    # ---------------------------------------------------------------------------------------
    # -------------- Below is printing for ctqmc  solver ------------------------------------
    # ---------------------------------------------------------------------------------------
    print('# CIX file for ctqmc! ', file=fcix)
    print('# cluster_size, number of states, number of baths, maximum_matrix_size', file=fcix)
    print(1, len(lowE), nbaths, low_maxsize, file=fcix)
    print('# baths, dimension, symmetry', file=fcix)

    if DiagonalOnly:
        for ib in range(nbaths):
            print(ib, '  ', 1, Sigind[bkeep[ib],bkeep[ib]]-1, '  ', global_flip[bkeep[ib]], file=fcix)
    else:
        for ib,w in enumerate(wcoupled):
            print(ib, len(w), end=' ', file=fcix)
            for i in w:
                for j in w:
                    print(Sigind_orig[i,j]-1, end=' ', file=fcix)
            if Q3d:
                print(ib % (len(wcoupled)/2), file=fcix)
            else:
                print(ib/2, file=fcix)
        print('# changing order of psi operators', file=fcix)
        print('FL_FROM_IFL', file=fcix)
        for ii,w in enumerate(wcoupled):
            for i in w:
                print(Sigind_orig[i,i]-1, end=' ', file=fcix)
            print(file=fcix)
    
    print('# cluster energies for non-equivalent baths, eps[k]', file=fcix)
    #for E in tEimp: print(E, end=' ', file=fcix)
    print(('{:.10f} '*len(tEimp)).format(*tEimp), file=fcix)
    #print(file=fcix)
    print('#   N   K   Sz size', file=fcix)

 
    for i in range(len(lowE)):
        ii = lowE[i][0]
        iwp = pseudostates[ii]
        iw = iwp[0]
        ip = iwp[1]
        wstate = wstates[iw]

        if CoulombF[:5] == 'Ising':
            gs=wstate[2][ip]
            nrm = 1.0
        else:
            ilargest = argmax(abs(Te[iw][:,ip]))
            gs = wstate[2][ilargest]
            nrm = abs(Te[iw][ilargest,ip])
        if op.Q3d:
            Mz = sum(wstate[1])/2.
        else:
            Mz = op.Mz(gs)
        
        print("%3d  %2d %2d %6.3f %2d " % (i+1, sum(wstate[0]), 0, Mz, len(lowE[i][1])), end=' ', file=fcix)
        for ib in bkeep:
            ifinal = iFinal[ii,ib]
            print("%3d" % (inv_lowE1[ifinal]+1), end=' ', file=fcix)
        print("  ", end=' ', file=fcix)
        for iq in lowE[i][1]:
            print("%12.8f" % (Enes[ii][iq],), end=' ', file=fcix)
        print("  ", end=' ', file=fcix)
        for iq in lowE[i][1]:
            print(S2ws[ii][iq], end=' ', file=fcix)
        #if CoulombF[:5] == 'Ising':
        #    print >> fcix, "  # ", op.printn(gs),
        print("  # ", op.printn(gs), nrm, end=' ', file=fcix)
        print(file=fcix)
        
    print('# matrix elements', file=fcix)

    for i in range(len(lowE)):
        ii = lowE[i][0]
        for ib in bkeep:
            ifinal = iFinal[ii,ib]
            low_ifinal = inv_lowE1[ifinal]
            print("%3d %3d " % (i+1, low_ifinal+1), end=' ', file=fcix)
            if low_ifinal>=0: 
                ind0 = lowE[i][1]
                ind1 = lowE[low_ifinal][1]
                print("%2d %2d" % (len(ind0), len(ind1)), end=' ', file=fcix) 
                for i0 in ind0:
                    for j0 in ind1:
                        x = rFKP[ii][ib][i0,j0]
                        if abs(x.imag)<1e-4 or PrintReal:
                            print("%20.16f" % (x.real,), end=' ', file=fcix)
                        else:
                            print(x, end=' ', file=fcix)
            else:
                print("%2d %2d" % (0, 0), end=' ', file=fcix)
            print(file=fcix)

    if HB2 : print('HB2', file=fcix)
    else: print('HB1', file=fcix)

    if (HB2 and Q3d):
        ii=0
        iind={}
        for i1,bs1 in enumerate(baths):
            m1 = bs1[0]
            s1 = bs1[1]
            if m1 not in bkeep: continue
            iind[i1]=ii
            ii+=1
        print("# UCoulomb : (m1,s1) (m2,s2) (m3,s2) (m4,s1)  Uc[m1,m2,m3,m4]", file=fcix)
        for i1,bs1 in enumerate(baths):
            m1 = bs1[0]
            s1 = bs1[1]
            if m1 not in bkeep: continue
            for i2,bs2 in enumerate(baths):
                m2 = bs2[0]
                s2 = bs2[1]
                if m2 not in bkeep: continue
                for i3,bs3 in enumerate(baths):
                    m3 = bs3[0]
                    s3 = bs3[1]
                    if (s2!=s3): continue
                    if m3 not in bkeep: continue
                    for i4,bs4 in enumerate(baths):
                        m4 = bs4[0]
                        s4 = bs4[1]
                        if (s4!=s1): continue
                        if m4 not in bkeep: continue
                        Uc = 0.0
                        if CoulombF[:5] == 'Ising' and not ((i1==i4 and i2==i3) or (i1==i3 and i2==i4)):
                            continue
                        for k in range(0,l+1):
                            Uc += real(UC[k,m1,m2,m3,m4])*Fk[k,l]
                        if abs(Uc)>1e-6:
                            print("%2d %2d %2d %2d  %12.8f" % (iind[i1],iind[i2],iind[i3],iind[i4],Uc), file=fcix)

    if (HB2 and not Q3d):
        print("# UCoulomb : (m1,s1) (m2,s2) (m3,s2) (m4,s1)  Uc[m1,m2,m3,m4]", file=fcix)
        for bs1 in bkeep:
            for bs2 in bkeep:
                for bs3 in bkeep:
                    for bs4 in bkeep:
                        Uc = 0.0
                        for k in range(0,l+1):
                            UC += real(UC[k,bs1,bs2,bs3,bs4])*Fk[k,l]
                        if abs(Uc)>1e-6:
                            print("%2d %2d %2d %2d  %12.8f" % (bs1,bs2,bs3,bs4,Uc), file=fcix)

                
    print('# number of operators needed', file=fcix)
    
    if not add_occupancy:
        print('0', file=fcix)
    else:
        print('1', end=' ', file=fcix)
        if ORB: print('1', end=' ', file=fcix)
        elif PrintSminus: print('2', end=' ', file=fcix)
        print(file=fcix)
        print('# Occupancy ', file=fcix)
        
        for i in range(len(lowE)):
            ii = lowE[i][0]
            ind0 = lowE[i][1]

            for ikb,bt in enumerate(tkbth):
                Oub = zeros((len(ind0),len(ind0)),dtype=float)
                for ib in bt:
                    Nm = zeros((len(ind0),len(ind0)),dtype=float)
                    if len(rNn[ii][ib])>0:  
                        for j in range(len(ind0)):
                            for k in range(len(ind0)):
                                Nm[j,k] = rNn[ii][ib][ind0[j],ind0[k]]
                                
                    Oub += identity(len(ind0))-Nm
                print(("%3d " % (i+1)), end=' ', file=fcix)
                print("%2d %2d" % (len(ind0), len(ind0)), end=' ', file=fcix) 
                for iw,i0 in enumerate(ind0):
                    for iz,j0 in enumerate(ind0):
                        ff = Oub[iw,iz]
                        if abs(ff)<small: ff=0.0
                        print(ff, end=' ', file=fcix)
                print(file=fcix)


        #### NEW
        if ORB:
            print('# Orbital susc with ', ORB, file=fcix)
            for i in range(len(lowE)):
                ii = lowE[i][0]
                ind0 = lowE[i][1]
                DM = zeros((len(bkeep),len(ind0),len(ind0)),dtype=float)
                for iib,ib in enumerate(bkeep):
                    N1m = zeros((len(ind0),len(ind0)),dtype=float)
                    if len(rNn[ii][ib])>0:  
                        for j in range(len(ind0)):
                            for k in range(len(ind0)):
                                N1m[j,k] = rNn[ii][ib][ind0[j],ind0[k]]
                                
                    DM[iib,:,:] = identity(len(ind0))-N1m
                
                Dorb = zeros((len(ind0),len(ind0)),dtype=float)
                for ib in range(len(bkeep)):
                    Dorb += DM[ib,:,:]*ORB[ib]
                    
                print("%3d %3d " % (i+1,i+1), end=' ', file=fcix)
                print("%2d %2d" % (len(ind0), len(ind0)), end=' ', file=fcix) 
                for iw,i0 in enumerate(ind0):
                    for iz,j0 in enumerate(ind0):
                        ff = Dorb[iw,iz]
                        if abs(ff)<small: ff=0.0
                        print(ff, end=' ', file=fcix)
                print(file=fcix)
            
        if PrintSminus:
            Sp_i={}
            Sp_m={}
            print('# S_minus followed by S_plus', file=fcix)
            for i in range(len(lowE)):
                ii = lowE[i][0]
                ind0 = lowE[i][1]
                print("%3d " % (i+1,), end=' ', file=fcix)
                if ii in iS_minus:
                    ifinal = iS_minus[ii]
                    low_ifinal = inv_lowE1[ifinal]
                    print("%3d " % (low_ifinal+1,), end=' ', file=fcix)
                    if low_ifinal>=0: 
                        ind0 = lowE[i][1]
                        ind1 = lowE[low_ifinal][1]
                        print("%2d %2d" % (len(ind0), len(ind1)), end=' ', file=fcix) 
                        Sp_i[low_ifinal+1] = i+1
                        Sp_m[low_ifinal+1] = zeros((len(ind1),len(ind0)))
                        for l0,i0 in enumerate(ind0):
                            for l1,j0 in enumerate(ind1):
                                x = rS_minus[ii][i0,j0]
                                if abs(x)<1e-10: x=0.0
                                Sp_m[low_ifinal+1][l1,l0] = x
                                print(x, end=' ', file=fcix)
                    else:
                        print("%2d %2d" % (0, 0), end=' ', file=fcix)
                else:
                    print("%3d %2d %2d" % (0, 0, 0), end=' ', file=fcix)
                print(file=fcix)
                
            #print >> fcix, '# S_plus'
            for i in range(len(lowE)):
                print("%3d " % (i+1,), end=' ', file=fcix)
                if i+1 in Sp_i:
                    ifinal = Sp_i[i+1]
                    M = Sp_m[i+1]
                    print("%3d %2d %2d " % (ifinal,shape(M)[0], shape(M)[1]), end=' ', file=fcix) 
                    for i0 in range(shape(M)[0]):
                        for j0 in range(shape(M)[1]):
                            print(M[i0,j0], end=' ', file=fcix)
                    #print >> fcix
                else:
                    print("%3d %2d %2d" % (0, 0, 0), end=' ', file=fcix)
                print(file=fcix)
        #### NEW
    if (OLD_CTQMC):
        print('# Data for HB1', file=fcix)
        
        PrintAll=False
        if PrintAll:
            print(1, len(pseudostates), len(bkeep), maxs, file=fcix)
            print('# ind   N   K   Jz size', file=fcix)
        
            for ii,iwp in enumerate(pseudostates):
                iw = iwp[0]
                ip = iwp[1]
                wstate = wstates[iw]
                print("%3d  %3d  %2d %2d %4.1f %2d " % (ii+1, inv_lowE1[ii]+1, sum(wstate[0]), 0, sum(wstate[1])/2., len(Enes[ii])), end=' ', file=fcix)
                for ib in bkeep:
                    print("%3d" % (iFinal[ii,ib]+1), end=' ', file=fcix)
                print("  ", end=' ', file=fcix)
                for iq in range(len(Enes[ii])):
                    print(Enes[ii][iq], end=' ', file=fcix)
                print("  ", end=' ', file=fcix)
                for iq in range(len(Enes[ii])):
                    print(0, end=' ', file=fcix)
                print("  ", end=' ', file=fcix)
                print(file=fcix)
                
            print('# matrix elements', file=fcix)
            
            for ii in range(len(pseudostates)):
                for ib in bkeep:
                        print("%3d %3d " % (ii+1, iFinal[ii,ib]+1), end=' ', file=fcix) 
                        ffp = zeros(len(Enes[ii]),dtype=float)
                        if iFinal[ii,ib]>=0:
                            (dim0, dim1) = shape(rFKP[ii][ib])
                            print("%2d %2d" % (dim0,dim1), end=' ', file=fcix) 
                            for i0 in range(dim0):
                                for j0 in range(dim1):
                                    x = rFKP[ii][ib][i0,j0]
                                    if abs(x.imag)<1e-4 or PrintReal:
                                        print(x.real, end=' ', file=fcix)
                                    else:
                                        print(x, end=' ', file=fcix)
                            for i0 in range(dim0):
                                dsum=0
                                for j0 in range(dim1):
                                    dsum += abs(rFKP[ii][ib][i0][j0])**2
                                ffp[i0] += dsum
                        else:
                            print("%2d %2d" % (0, 0), end=' ', file=fcix)
                        print(file=fcix)
        else:
            print(1, len(lowE), nbaths, low_maxsize, file=fcix)
            print('# ind   N   K   Jz size', file=fcix)
            
            for i in range(len(lowE)):
                ii = lowE[i][0]
                iwp = pseudostates[ii]
                iw = iwp[0]
                ip = iwp[1]
                wstate = wstates[iw]
                if op.Q3d:
                    Mz = sum(wstate[1])/2.
                    gs=wstate[2][ip]
                else:
                    gs=wstate[2][ip]
                    Mz = op.Mz(gs)
                print("%3d %3d  %2d %2d %6.3f %2d " % (i+1, i+1, sum(wstate[0]), 0, Mz, len(lowE[i][1])), end=' ', file=fcix)
                for ib in bkeep:
                    ifinal = iFinal[ii,ib]
                    print("%3d" % (inv_lowE1[ifinal]+1), end=' ', file=fcix)
                print("  ", end=' ', file=fcix)
                for iq in lowE[i][1]:
                    print("%10.6f" % (Enes[ii][iq],), end=' ', file=fcix)
                print("  ", end=' ', file=fcix)
                for iq in lowE[i][1]:
                    print(S2ws[ii][iq], end=' ', file=fcix)
                if CoulombF[:5] == 'Ising':
                    print("  # ", op.printn(gs), end=' ', file=fcix)
                print(file=fcix)
                
            print('# matrix elements', file=fcix)
            
            for i in range(len(lowE)):
                ii = lowE[i][0]
                for ib in bkeep:
                    ifinal = iFinal[ii,ib]
                    low_ifinal = inv_lowE1[ifinal]
                    print("%3d %3d " % (i+1, low_ifinal+1), end=' ', file=fcix)
                    if low_ifinal>=0: 
                        ind0 = lowE[i][1]
                        ind1 = lowE[low_ifinal][1]
                        print("%2d %2d" % (len(ind0), len(ind1)), end=' ', file=fcix) 
                        for i0 in ind0:
                            for j0 in ind1:
                                x = rFKP[ii][ib][i0,j0]
                                if abs(x.imag)<1e-4 or PrintReal:
                                    print(x.real, end=' ', file=fcix)
                                else:
                                    print(x, end=' ', file=fcix)
                    else:
                        print("%2d %2d" % (0, 0), end=' ', file=fcix)
                    print(file=fcix)

