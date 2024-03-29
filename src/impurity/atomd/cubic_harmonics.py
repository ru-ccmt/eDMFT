#!/usr/bin/eny python
# @Copyright 2007 Kristjan Haule
from math import *
import numpy as np

def Spheric2Cubic(L):
    """Transformation matrix from spherical to cubic harmonics:
    From |l,m> (l=0,1,2,3) basis to
    basis of triplets, doublets and singlets real functions 
    (see appendices 2.1 and 4.2 of H. Watanabe 'Operator Methods in Ligand Field Theory'
    Prentice Hall, 1966.)
    """
    s2 = np.sqrt(2)
    f=1/s2
    g=1j/s2
    if L == 0:
        return np.array([[1]])
    elif L == 1:
        T = [[ f, 0,-f],  # x
             [ g, 0, g],  # y
             [ 0, 1, 0]]  # z
        return np.array(T)
    elif L == 2:
        T = [[ 0, 0, 1, 0, 0],  # z^2
             [ f, 0, 0, 0, f],  # x^2-y^2
             [ 0, f, 0,-f, 0],  # xz
             [ 0, g, 0, g, 0],  # yz
             [ g, 0, 0, 0,-g]]  # xy
        return np.array(T)
    elif L == 3:
        s3 = np.sqrt(3)
        s5 = np.sqrt(5)
        s8 = np.sqrt(8)
        a=s5/4.
        b=s3/4.
        c=1j*a
        d=1j*b
        T=[[ 0, g, 0, 0, 0,-g, 0], # singlet
           [ a, 0,-b, 0, b, 0,-a], # x
           [-c, 0,-d, 0,-d, 0,-c], # y
           [ 0, 0, 0, 1, 0, 0, 0], # z
           [-b, 0,-a, 0, a, 0, b], # ksi
           [-d, 0, c, 0, c, 0,-d], # eta
           [ 0, f, 0, 0, 0, f, 0]] # zeta
        return np.array(T)
    
def Spheric2jj(L):
    """Transformation matrix from spherical harmonics to jj-coupling scheme: 
     From |L,1/2,mL,mS> (L=0,1,2,3, mS=-1/2,1/2) basis to
     basis |J,L,S,mJ>, J=L-1/2, L+1/2 
     H. Watanabe 'Operator Methods in Ligand Field Theory'
     Prentice Hall, 1966, Table 1.8-1.
     ordering because of the convention used in WIEN is:
                        mS=1/2        mS=-1/2
                      -L .......L  -L ...... L     (2*(2L+1) columns)
             -(L-1/2)
                .
     J=L-1/2    .
                .
              (L-1/2)
              -L-1/2 
                .
     J=L+1/2    .
                .
               L+1/2 
    
    """
    if L == 0:
        return np.array([[0,1],
                         [1,0]])
    else:
        cf = np.zeros((2*(2*L+1),2*(2*L+1)),dtype=complex)
        ms_ml = [(-1,L) for L in range(-L,L+1)] + [(1,L) for L in range(-L,L+1)]
        for k1,(ms,ml) in enumerate(ms_ml):
            ams = -ms/2
            for k2,mj in enumerate(range(-2*L+1, 2*L, 2)):  # L-1/2 states
                amj = mj/2
                d = amj - ml - ams
                if abs(d) < 0.0001:
                    if ms == 1:
                        cf[k2, k1] = -np.sqrt((L+0.5+amj)/(2*L+1))
                    else:
                        cf[k2, k1] = np.sqrt((L+0.5-amj)/(2*L+1))
            for k2,mj in enumerate(range(-2*L-1, 2*L+2, 2)):  # L+1/2 states
                amj = mj/2
                d = amj - ml - ams
                if abs(d) < 0.0001:
                    if ms == 1:
                        cf[k2+2*L, k1] = np.sqrt((L+0.5-amj)/(2*L+1))
                    else:
                        cf[k2+2*L, k1] = np.sqrt((L+0.5+amj)/(2*L+1))
    return cf

def Spheric2EffHalf():
    """ Transforms Spheric harmonics for L=2 to j_eff=1/2 basis for 5d orbitals like Ir
    """
    s23 = sqrt(2./3.)
    s16 = sqrt(1./6.)
    s29 = sqrt(2./9.)
    s12 = sqrt(1./12.)
    s34 = sqrt(3./4.)
    s2  = sqrt(1./2.)
    alp = 1/2.
    T=[
        #  -2,dn,  -1,dn,  0,dn,   1,dn,  2,dn, -2,up, -1,up,  0,up,  1,up,  2,up,                  l.s
        [    0j,    s23,    0j,     0j,    0j,   s16,    0j,    0j,    0j,  -s16],  # sef=-1/2      1.0
        [  -s16,     0j,    0j,     0j,   s16,    0j,    0j,    0j,   s23,    0j],  # sef= 1/2      1.0
        [    0j,    alp,    0j,    alp,    0j,  -alp,    0j,    0j,    0j,   alp],  # sef=-3/2     -0.5
        [   alp,     0j,    0j,     0j,  -alp,    0j,   alp,    0j,   alp,    0j],  # sef= 3/2     -0.5
        [  -s12,     0j,    0j,     0j,   s12,    0j,   s34,    0j,  -s12,    0j],  # sef=-1/2     -0.5
        [    0j,    s12,    0j,   -s34,    0j,  -s12,    0j,    0j,    0j,   s12],  # sef= 1/2     -0.5 
        [    0j,     0j,    1.,     0j,    0j,    0j,    0j,    0j,    0j,    0j],  # z^2,up        0.0
        [    0j,     0j,    0j,     0j,    0j,    0j,    0j,    1.,    0j,    0j],  # z^2,dn        0.0
        [    s2,     0j,    0j,     0j,    s2,    0j,    0j,    0j,    0j,    0j],  # x^2-y^2,up    0.0
        [    0j,     0j,    0j,     0j,    0j,    s2,    0j,    0j,    0j,    s2]]  # x^2-y^2,dn    0.0
    return np.array(T)

#def PrintM(T):
#    for i in range(np.shape(T)[0]):
#        for j in range(np.shape(T)[1]):
#            print('{:8.5f}{:8.5f}i '.format(T[i,j].real,T[i,j].imag), end='')
#        print()
            
if __name__ == '__main__':
    for L in range(4):
        T = Spheric2Cubic(L)
        #T2 = trafo.trafo(L,False)
        check = np.allclose(T @ (T.conj()).T, np.identity(2*L+1))
        #check2 = np.allclose(T, T2)
        print(L, check)
    for L in range(4):
        T = Spheric2jj(L)
        #T2 = trafoso.trafoso(L)
        check = np.allclose(T @ (T.conj()).T, np.identity(2*(2*L+1)))
        #check2 = np.allclose(T,T2)
        print(L, check)
