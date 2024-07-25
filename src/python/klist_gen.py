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
#from pymatgen.io.cif import CifParser
#from pymatgen.symmetry.bandstructure import HighSymmKpath
#from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from fractions import Fraction
#from localaxes import *
#from functools import cmp_to_key
from utils import W2kEnvironment

from wlatgen import Latgen
from wstruct import Struct
    
def W2k_klist_band(fname, Nt, labels, Ks, k2icartes, k2cartes, log=sys.stdout):
    """
    fname     -- filename, should be case.klist_band
    Nt        -- total number of k-points (approximately, because we find integer denominated path)
    kpath     -- dictionary of {'label': array(3)} given in primitive BZ
    k2icartes -- transformation from primitive to conventional BZ. It is computed during ci2struct and is in pickle
    k2cartes  -- transformation from primitive BZ to cartesian coordinates. Allows one to compute distance in momentum space
    """
    def common_denom(i1,i2):
        if ( 2*abs(i1-i2)/(i1+i2) < 0.05 ):  # this is approximation. If i1 and i2 are very close, we do not want to take the product, but approximate
            return max(i1,i2)
        else:
            return int(i1*i2/gcd(i1,i2))

    print('Kpoints in the mesh:', file=log)
    dst = np.zeros(len(Ks))             # first we compute distance between momentum points, so that 
    for i in range(len(labels)):           # k-points will be uniformly distributed in cartesian coordinates
        kc = k2cartes @ k2icartes @ Ks[i]  # and two successive k-points will not have the same number of points
        if i>0:                            # ks is momentum point in cartesian coordinates
            dst[i] = dst[i-1] + linalg.norm(kc-kc_p) # cumulative distance from the first k-point
        print('{:10}'.format(labels[i]), 'k-PBZ=[{:6.4g},{:6.4g},{:6.4g}]'.format(*Ks[i]),
                  'k-conBZ=[{:6.4g},{:6.4g},{:6.4g}]'.format(*(k2icartes@Ks[i])), file=log) # primitive and conventional BZ
        kc_p = kc[:]
    # Nkp is number of k-points in each interval
    Nkp = [round(Nt*(dst[ii]-dst[ii-1])/dst[-1]) for ii in range(1,len(dst))]
    # 
    print('suggested and actual number of momentum points:', Nt, sum(Nkp), Nkp, file=log)

    with open(fname, 'w') if fname!=None else sys.stdout as fk:
        ii, dst = 0, 0.0
        kc_p = k2cartes @ k2icartes @ Ks[0]
        for i in range(len(labels)-1):
            if Nkp[i]==0: continue
            k1 = k2icartes @ Ks[i]       # Ks[i] is given in primitive cell, while k1 in conventional
            k2 = k2icartes @ Ks[i+1]     # k2 is next point in conventional, as required by w2k
            frk1 = [Fraction(str(k1[j])).limit_denominator(10) for j in range(3)] # fractions for in representation
            r2 = (k2-k1)/Nkp[i]          # we will need k1 + r2*i
            frk2 = [Fraction(str(r2[j])).limit_denominator(Nkp[i]*10) for j in range(3)] # fraction representation
            f1 = common_denom(common_denom(frk1[0].denominator,frk1[1].denominator),frk1[2].denominator) # common-gcd
            f2 = common_denom(common_denom(frk2[0].denominator,frk2[1].denominator),frk2[2].denominator) # common-gcd
            Dk = common_denom(f1,f2) # finally common denominator of all these fractions
            #print(k1,k2, Dk)
            for ik in range(Nkp[i]):
                ki = k1 + (k2-k1)*ik/Nkp[i]   # k-point in conventional BZ
                kc = k2cartes @ ki            # k-point in cartesian coordinates
                dst += linalg.norm(kc-kc_p)   # distance from first point
                k_int = array(np.round(ki*Dk), dtype=int) # integer representation in conventional BZ
                if ik==0:
                    print('{:10s}{:10d}{:10d}{:10d}{:10d}{:20s}{:6f}'.format(labels[i],*k_int,Dk,'',dst), file=fk)
                else:
                    #print(ii+1, k_int, Dk)
                    print('{:10s}{:10d}{:10d}{:10d}{:10d}{:20s}{:6f}'.format('',*k_int,Dk,'',dst), file=fk)
                ii += 1
                kc_p = kc
        # still need to add the last point
        ki = k2icartes @ Ks[-1]
        kc = k2cartes @ ki
        dst += linalg.norm(kc-kc_p)
        k_int = array(np.round(ki * Dk), dtype=int)
        print('{:10s}{:10d}{:10d}{:10d}{:10d}{:20s}{:6f}'.format(labels[-1],*k_int,Dk,'',dst), file=fk)
        print('END', file=fk)


if __name__ == '__main__':
    import readline
    
    usage = """usage: %klist_gen.py [ options ]
    Generates w2k klist along some path
       No options would start interactive mode.
       Options would prevent interactive mode:
        -n 200         -- number of points along the path
        -p  ([0,0,0],[0.5,0,0],[0.5,0.5,0],[0,0,0])   -- special points in the path expressed in primitive BZ!
        -m  ["Gamma","X","M","Gamma"]
        -o  output fname
    """
    inpt={}  # input from command line
    if len(sys.argv)>1:
        if sys.argv[1] in ['-h', '--help']:
            print(usage)
            sys.exit(0)
        else:
            for i in range(int(len(sys.argv)/2)):
                inpt[sys.argv[2*i+1][1:]] = sys.argv[2*i+2]
            print('given options are', inpt)
            Nkp=200
            path_points=[]
            path_names=[]
            fname = 'case.klist_band'
            if 'n' in inpt:
                Nkp=int(inpt['n'])
            if 'p' in inpt:
                V = eval(inpt['p'])
                path_points=([array(x) for x in V])
            if 'm' in inpt:
                path_names=eval(inpt['m'])
            if 'o' in inpt:
                fname = inpt['o']
    else:
        while True:
            try:
                userin = input("""Give number of k-points along the path (ex: 200): 
  %""")
                Nkp = int(userin)
            except:
                print("---> specified number does not convert to integer {:s}.".format(str(userin)))
                prompt = "Do you want to retry? (y/n): "
                userin = input(prompt).lower().strip()
                if userin.lower() == 'n':
                    exit(1)
                else:
                    continue
            break
        while True:
            try:
                userin = input("""Give special points along the path as two dimensional python list (ex: [[0,0,0],[0.0.5,0],...]): 
  %""")
                path_points = eval(userin)
                path_points = array(path_points)
                #print(path_points)
            except:
                print("---> specified list does not convert to python array {:s}.".format(userin))
                prompt = "Do you want to retry? (y/n): "
                userin = input(prompt).lower().strip()
                if userin.lower() == 'n':
                    exit(1)
                else:
                    continue
            break
        while True:
            try:
                userin = input("""Give name to this special points along the path as python list (ex: ['Gamma','X',...]): 
  %""")
                path_names = eval(userin)
                if len(path_names)!=len(path_points):
                    raise Exception('not compatible lengths')
                    
            except Exception as err:
                print(err)
                print("---> specified list does not seems to work {:s}.".format(userin))
                prompt = "Do you want to retry? (y/n): "
                userin = input(prompt).lower().strip()
                if userin.lower() == 'n':
                    exit(1)
                else:
                    continue
            break
        while True:
            try:
                userin = input("""Give name of the resulting filename (ex: case.klist): 
  %""")
                fname = userin
                f = open(fname,'a')
                f.close()
            except:
                print(err)
                print("---> specified name does not work {:s}.".format(userin))
                prompt = "Do you want to retry? (y/n): "
                userin = input(prompt).lower().strip()
                if userin.lower() == 'n':
                    exit(1)
                else:
                    continue
            break
        
    print('Got the following input:')
    print('Nkp=', Nkp)
    for i in range(len(path_points)):
        print(' {:10s}'.format(path_names[i]), path_points[i].tolist())
    print('writing to', fname)

            
    w2k = W2kEnvironment()
    with open('klist_gen.log','w') as log:
        strc = Struct()
        strc.ReadStruct(w2k.case+'.struct', log)
        latgen = Latgen(strc, w2k.case, log)
        W2k_klist_band(fname, Nkp, path_names, path_points, latgen.k2icartes, latgen.k2cartes, log)
