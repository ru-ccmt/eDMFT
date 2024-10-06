#!/usr/bin/env python

""" Given a list of special points, it constructs w2k list case.klist_band.
"""
# @Copyright 2024 Kristjan Haule
from math import gcd
import re, sys, os
import optparse
from numpy import *
import numpy as np
import numpy.linalg as linalg
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D

from fractions import Fraction
from utils import W2kEnvironment

from wlatgen import Latgen
from wstruct import Struct

def common_denom(i1,i2):
    if ( 2*abs(i1-i2)/(i1+i2) < 0.05 ):  # this is approximation. If i1 and i2 are very close, we do not want to take the product, but approximate
        return max(i1,i2)
    else:
        return int(i1*i2/gcd(i1,i2))

def toStr(ki):
    return '{:10.6f},{:10.6f},{:10.6f}'.format(*ki)

def to1BZ(ki):
    for i in range(3):
        if ki[i]>0.5+1e-16: ki[i] -= 1
    return ki

if __name__ == '__main__':
    # W2k requires
    # Z=[0,0,2/3]
    # X=[1/2,1/2,0]
    # P=[1/2,1/2,1/2]
    # M=[1,0,0]
    # 
    CheckIrreducible = True
    Nxy = 90
    xi = [-1.,1.]
    yi = [-1.,1.]
    # Input in cartesian or conventional-BZ
    #kpth = '[x*2/3.,y*2/3.,0.5*y]'
    #kpth = '[x*2/3.,0,0.5*y]'
    #kpth = '[0,x*2/3.,0.5*y]'
    #kpth = '[0.288675134594813*x,0.5*x,0.5*y]'
    #kpth = '[0.333333333333333*x,0.577350269189626*x,0.5*y]'
    #kpth = '[0.5*x,0.5*y,0.5]'
    #kpth = '[0.5*x,0,0.5*y]'

    #kpth = '[0.5*x,0.666667*y,0]'
    #kpth = '[x,0,0.666666666666667*y]'
    #kpth = '[0.5*(x+y),0.5*(x-y),0]'

    #a=0.5327283770
    #c=0.6933442306
    #x0=sqrt(3)/2.*a
    #y0=0.5*a
    #x0,y0=sqrt(3)/2*a,0.5*a
    
    #kpth = '['+str(x0)+'*x,'+str(a)+'*y,1.5]'
    #kpth = '['+str(x0)+'+'+str(c-x0)+'*0.5*(1+x),'+str(y0)+'*y,1.5-0.5*(1+x)]'
    
    exec(compile(open('2D_params.py', 'rb').read(), '2D_params.py', 'exec'))
    print('kpth=', kpth)
    
    xn = xi[0] + (xi[1]-xi[0])*arange(0,Nxy+1)/Nxy
    yn = yi[0] + (yi[1]-yi[0])*arange(0,Nxy+1)/Nxy

    log = sys.stdout
    w2k = W2kEnvironment()
    strc = Struct()
    strc.ReadStruct(w2k.case+'.struct', log)
    latgen = Latgen(strc, w2k.case, log)

    # pia = [2*pi/a, 2*pi/b, 2pi/c]
    k2icartes = (diag(1/latgen.pia) @ latgen.br2)
    c2f = latgen.k2icartes @ linalg.inv(k2icartes)

    #c2f = latgen.k2icartes @ linag.inv(latgen.br2) @ diag(latgen.pia)
    
    print('k2icartes=', k2icartes)
    print('c2f=', c2f)
    
    imat = zeros((len(strc.timat),3,3))
    for isym in range(len(strc.timat)):
        imat[isym,:,:] = strc.timat[isym,:,:].T
    iind=0
    irred={}
    rxy = zeros((len(xn),len(yn)), dtype=int)
    for j in range(len(yn)):
        for i in range(len(xn)):
            x,y = xn[i],yn[j]
            kc = array(eval(kpth))
            #kf = c2f @ array([xn[i],yn[j],0])
            kf = c2f @ kc

            for l in range(3):
                if abs(kf[l])<1e-5: kf[l]=0
                
            ki = kf % 1.0
            kis = toStr(ki)
            #print(i, j, kis, kf.tolist())
            Found = False
            if kis in irred:
                Found = True
                rxy[i,j] = irred[kis][0]
                print('{:3d} {:3d} [{:10.6f},{:10.6f},{:10.6f}]'.format(i,j, *kf), 'is Reduc #{:-4d}  {:s}'.format(rxy[i,j], kis),
                          'operation', 1, 'transforms to', '[{:10.6f},{:10.6f},{:10.6f}]'.format(*ki))
                
            if not Found and CheckIrreducible:
                ksym = imat @ kf
                ksymi = ksym % 1.0
                for isym in range(len(imat)):
                    ksymii = toStr(ksymi[isym])
                    if ksymii in irred:
                        Found = True
                        rxy[i,j] = irred[ksymii][0]
                        print('{:3d} {:3d} [{:10.6f},{:10.6f},{:10.6f}]'.format(i,j, *kf), 'is Reduc #{:-4d}  {:s}'.format(rxy[i,j], ksymii),
                                  'operation', isym, 'transforms to', '[{:10.6f},{:10.6f},{:10.6f}]'.format(*ksym[isym]))
                        break
                    
            if not Found:
                irred[kis] = (iind,to1BZ(ki))
                rxy[i,j] = iind
                print('{:3d} {:3d} [{:10.6f},{:10.6f},{:10.6f}]'.format(i,j, *kf), 'is Irred #{:-4d}  {:s}'.format(iind, kis))
                iind += 1
    Niid = iind
    print('Number of all irreducible points=', Niid, 'versus all=', len(xn)*len(yn))
    

    # We expect iind to be range(), but just checking
    iind = [irred[kis][0] for kis in irred]
    if sum(abs( arange(0,len(irred)) - array(iind))) > 0.1:
        print('ERROR: index need to be sorted', iind)
        sys.exit(0)

    ind2kis = [kis for kis in irred]

    #max_denom=2000
    max_denom=500
    with open(w2k.case+'.klist_band', 'w') as fo:        
        for kis in irred:
            kf = irred[kis][1]
            frk = [Fraction(str(kf[j])).limit_denominator(max_denom) for j in range(3)]
            f1 = common_denom(common_denom(frk[0].denominator,frk[1].denominator),frk[2].denominator) # common-gcd
            ki = array(round_(kf*f1), dtype=int)
            print(kis, '[{:16.8f},{:16.8f},{:16.8f}]/{:d}'.format(*(kf*f1),f1), ki)
            print('{:10s}{:10d}{:10d}{:10d}{:10d}'.format( str(irred[kis][0]), *ki, f1), file=fo)
        print('END', file=fo)
    with open('kindex.dat','w') as fx:
        print('# kpth=', kpth, 'Nxy=', Nxy, file=fx)
        print('# xi=', xi, file=fx)
        print('# yi=', yi, file=fx)
        for i in range(shape(rxy)[0]):
            for j in range(shape(rxy)[1]):
                print('{:3d} {:3d} {:6d}'.format(i,j,rxy[i,j]), '  # ', ind2kis[rxy[i,j]], file=fx)
                
    print('Number of all irreducible points=', Niid, 'versus all=', len(xn)*len(yn))

