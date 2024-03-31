#!/usr/bin/env python
"""
sigen.py -- Sigma Index Generator
Also computes legends and crystal field transformations.
Takes as input case.indmf file
Writes case.indmfl file
"""
# @Copyright 2007 Kristjan Haule
from numpy import asarray, zeros, diag, bmat, cumsum, where, array, vstack, hstack
import numpy as np
import optparse, os, traceback
from cubic_harmonics import Spheric2Cubic, Spheric2jj, Spheric2EffHalf
import indmffile
from utils import eV2Ry, W2kEnvironment
from math import *
import sys

# functions to handle nspins in generating sigind from raw lists given in qsplit_table
#
def no_op(sigind, nspins, leg=None):
    if leg==None:
        return sigind
    else:
        return leg

def dup(sigind, nspins, leg=None):
    if leg==None:
        # This is used to produce sigind, because third argument is None
        return [x for x in sigind if x>0]*nspins + [x for x in sigind if x<=0]*nspins
    else:
        # This is used to produce corresponding legends.
        return [leg[i] for i,x in enumerate(sigind)]*nspins + [leg[i] for i,x in enumerate(sigind)]*nspins
        
def dup_shift(sigind, mult, leg=None):  # duplicate and shift
    if leg==None:
        ret = []
        for i in range(mult):
            ret += [x+i*max(sigind) if x!=0 else 0 for x in sigind]
        return ret
    else:
        ret = []
        for i in range(mult):
            ret += [leg[j] if x!=0 else 0 for j,x in enumerate(sigind)]
        return ret

qsplit_table = {
    #  require  action     transtype        L = 0  L = 1          L = 2                  L = 3
    #  2 spins
    #    ------  ------     ---------        -----  -----          -----                  -----
    0 : (False, dup,       "none",         ([1],   [1,1,1],       [1,1,1,1,1],           [1,1,1,1,1,1,1]              )),  # averaged
    1 : (True,  no_op,     "relativistic", ([1,1], [1,1,2,3,3,2], [1,2,2,1,3,4,5,5,4,3], [1,2,3,3,2,1,4,5,6,7,7,6,5,4])),  # |j,mj>
   -1 : (True,  no_op,     "relativistic", ([1,2], list(range(1,7)),    list(range(1,11)),           list(range(1,15)))),  # |j,mj> no time-reversal symmetry
    2 : (False, dup,       "real",         ([1],   [1,2,3],       [1,2,3,4,5],           [1,2,3,4,5,6,7]              )),  # no symmetry, real
   -2 : (False, dup_shift, "real",         ([1],   [1,2,3],       [1,2,3,4,5],           [1,2,3,4,5,6,7]              )),  # no symmetry, real spin-polarized
    3 : (False, dup,       "t2g",          ([1],   [1,2,3],       [1,2,3,0,0],            [1,2,3,4,5,6,7]) ),              # t2gs
   -3 : (False, dup,       "eg",           ([1],   [1,2,3],       [1,2,0,0,0],            [1,2,3,4,5,6,7]) ),              # egs
    4 : (True,  no_op,     "relativistic", ([1,1], [1,1,2,2,2,2], [1,1,1,1,2,2,2,2,2,2], [1,1,1,1,1,1,2,2,2,2,2,2,2,2])),  # |j,mj> jm's equivalent
    5 : (False, dup,       "real",         ([1],   [1,1,2],       [1,2,3,3,4],           [1,2,2,3,4,4,5]              )),  # axial    
    6 : (False, dup,       "real",         ([1],   [1,1,2],       [1,2,3,3,2],           [1,2,3,4,2,3,1]              )),  # hexagonal
    7 : (False, dup,       "real",         ([1],   [1,1,1],       [1,1,2,2,2],           [1,2,2,2,3,3,3]              )),  # cubic
    8 : (True,  dup_shift, "real",         ([1],   [1,1,2],       [1,2,3,3,4],           [1,2,2,3,4,4,5]              )),  # axial spin-polarized
    9 : (True,  dup_shift, "real",         ([1],   [1,1,2],       [1,2,3,3,2],           [1,2,3,4,2,3,1]              )),  # hexagonal spin-polarized
   10 : (True,  dup_shift, "real",         ([1],   [1,1,1],       [1,1,2,2,2],           [1,2,2,2,3,3,3]              )),  # cubic spin-polarized
   11 : (True,  no_op,     "relativistic", ([1,1], [1,1,2,3,3,2], [1,2,2,1,3,4,5,5,4,3], [1,2,3,3,2,1,4,5,6,7,7,6,5,4])),  # |j,mj> + off-diagonal
   12 : (False, dup_shift, "real",         ([1],   [1,2,3],       [1,2,3,4,5],           [1,2,3,4,5,6,7]              )),  # real + off-diagonal
   13 : (True,  no_op,     "5d-jeff",      ([1,1], [1,1,2,3,3,2], [1,1,2,2,3,3,4,4,5,5], [1,2,3,3,2,1,4,5,6,7,7,6,5,4])),  # |J_eff=1/2> 
   14 : (True,  no_op,     "5d-jeff",      ([1,2], [1,2,3,4,5,6], [1,2,3,4,5,6,7,8,9,10], [1,2,3,4,5,6,7,8,9,10,11,12,13,14])),  # |J_eff=1/2> 
    }

def make_legend_table():
    def _j(jj):  return ' '.join(["(%d/2,%d/2)" % (2*jj+1, x) for x in range(2*jj+1, 0, -2)])
    def _jm(jj): return ' '.join([str(x)+'/2' for x in range(-2*jj-1, 2*jj+2, 2)])

    j_labels  = (_j(0),  _j(0)+' '+_j(1),  _j(1)+' '+_j(2), _j(2)+' '+_j(3))
    real      = ('r', 'x y z', 'z^2 x^2-y^2 xz yz xy', 'xyz x^3 y^3 z^3 x(y^2-z^2) y(x^2-z^2) z(x^2-y^2)')
    axial     = ('r', 'x+y z', 'z^2 x^2-y^2 xz+yz xy', 'xyz x^3+y^3 z^3 x(y^2-z^2)+y(x^2-z^2) z(x^2-y^2)')
    hexagonal = ('r', 'x+y z', 'z^2 x^2-y^2+xy yz+xz', 'xyz+zeta(T2) x(T1)+ksi(T2) y(T1)+eta(T2) z(T1)')
    cubic     = ('r', 'x+y+z', 'eg t2g', 'xyz T1 T2')
    jeff_1_2  = (_j(0),  _j(0)+' '+_j(1), '(1/2,1/2) (3/2,3/2) (3/2,1/2) z^2 x^2-y^2', _j(2)+' '+_j(3))
    jeff_1_2_ns  = (_jm(0), _jm(0)+' '+_jm(1), '(1/2,1/2) (1/2,-1/2) (3/2,3/2) (3/2,-3/2) (3/2,1/2) (3/2,-1/2) z^2,up z^2,dn x^2-y^2,up x^2-y^2,dn', _jm(2)+' '+_jm(3))
    t2gs      = ('r', 'x y z', 'xz yz xy', 'xyz x^3 y^3 z^3 x(y^2-z^2) y(x^2-z^2) z(x^2-y^2)')
    egs       = ('r', 'x y z', 'z^2 x^2-y^2', 'xyz x^3 y^3 z^3 x(y^2-z^2) y(x^2-z^2) z(x^2-y^2)')
    
    table = {
        0  : ('s', 'p', 'd', 'f'),
        -1 : (_jm(0), _jm(0)+' '+_jm(1), _jm(1)+' '+_jm(2), _jm(2)+' '+_jm(3)),
        1  : j_labels,
        2  : real,
        -2 : real,
        3  : t2gs,
       -3  : egs,
        4  : ('1/2', '1/2 3/2', '3/2 5/2', '5/2 7/2'),
        5  : axial,
        6  : hexagonal,
        7  : cubic,
        8  : axial,
        9  : hexagonal,
        10 : cubic,
        11 : j_labels,
        12 : real,
        13 : jeff_1_2,
        14 : jeff_1_2_ns,
        }
    return table

legend_table = make_legend_table()



def parse_cmdline_args():
    '''Instantiates command line argument parser, and parses the passed arguments.'''

    # First, instantiate the parser by defining usage of program,
    # and all possible options
    usage = """usage: %prog [options]
    Generates sigma index matrices and transformation matrices implementing desired qsplit.

    Input:  case.indmf

    Output: case.indmfl

    The possible values of qsplit are as follows:
    """ + '\n' + indmffile.qsplit_doc

    parser = optparse.OptionParser(usage)

    parser.add_option("--so", "--spin-orbit", action="store_true", default=False, help="perform calculation with spin-orbit")
    parser.add_option("--sig", "--self-energy", default="sig.inp", help="filename for self-energy: -sig SIGNAME")

    # Next, parse the arguments
    (options, args) = parser.parse_args()

    # There should be no arguments passed (only options)
    if len(args) > 0:
        parser.error("Unexpected argument(s): " + " ".join(args))

    return options


def add_offdiag_sigind(sigind):
    """ Changes Sigind such that the off-diagonal matrix elements are present."""
    icounter = max(diag(sigind)) + 1
    for i in range(len(sigind)):
        for j in range(i+1,len(sigind)):
            sigind[i,j] = icounter
            sigind[j,i] = icounter
            icounter += 1

def add_offdiag_legend(legend):
    n = len(legend)
    icounter = n+1
    for i in range(n):
        for j in range(i+1,n):
            legend.append('('+str(i+1)+','+str(j+1)+')')
            icounter += 1

def cmp_sigind_legend(qsplit, L, nsites, nspins):
    """Computes indices to self-energy, which determines non-zero matrix elements and their symmetry."""
    func = qsplit_table[qsplit][1]                # dup
    sigind_base = qsplit_table[qsplit][3][L]      # [1,1,2,2,2]
    sigind_base_unique = list(set(sigind_base)-{0})
    legend_base = legend_table[qsplit][L].split() # ['eg','t2g']
    sigind_one  = func(sigind_base, nspins)       #
    sigind_one_unique = list(set(sigind_one)-{0})
    #print('sigind_base=', sigind_base, 'sigind_base_unique=', sigind_base_unique, 'legend_base=', legend_base, 'nspins=', nspins)
    #print('sigind_one=', sigind_one, 'sigind_one_unique=', sigind_one_unique)
    legend_one  = func(sigind_base_unique, nspins, legend_base)
    #print('legend_one=', legend_one)
    sigind = diag(dup_shift(sigind_one, nsites))
    #print('sigind_one=',sigind_one) 
    legend = dup_shift(sigind_one_unique, nsites, legend_one)
    if qsplit in [11, 12]:
        add_offdiag_sigind(sigind)
        add_offdiag_legend(legend)
    return (sigind,legend)


def cmp_cftrans(qsplit, L, nsites, nspins):
    """Computes transformation matrix from spherical harmonics to basis specified by qsplit."""
    transtype = qsplit_table[qsplit][2]
    if transtype == 'real':
        T0 = Spheric2Cubic(L)
        if nspins == 1:
            T = T0
        else:
            Z = zeros((2*L+1,2*L+1), dtype=complex)
            T = bmat([[T0, Z], [Z, T0]])
    elif transtype == 'relativistic':
        T = Spheric2jj(L)
    elif transtype == 'none':
        T = np.identity(nspins*(2*L+1), dtype=complex)
    elif transtype == '5d-jeff':
        if (L==2):
            T = Spheric2EffHalf()
        else:
            T = Spheric2jj(L)
    elif transtype == 't2g':
        T0 = Spheric2Cubic(L)
        if L == 2:
            if nspins==1:
                T = vstack((T0[2:,:],T0[:2,:]))
            else:
                z1=zeros(len(T0))
                T1=hstack( (T0[2:,:],vstack((z1,z1,z1))) )  # t2g,up
                T2=hstack( (vstack((z1,z1,z1)),T0[2:,:]) )  # t2g,dn
                T3=hstack( (T0[:2,:],vstack((z1,z1))) )     # eg,up
                T4=hstack( (vstack((z1,z1)),T0[:2,:]) )     # eg,dn
                T = vstack((T1,T2,T3,T4))
        else:
            if nspins==1:
                T = T0
            else:
                Z = zeros((len(T0), len(T0)), dtype=complex)
                T = bmat([[T0, Z], [Z, T0]])
    elif transtype == 'eg':
        T0 = Spheric2Cubic(L)
        if nspins==2:
            if L==2:
                z1=zeros(len(T0))
                T1=hstack( (T0[:2,:],vstack((z1,z1))) )     # eg,up
                T2=hstack( (vstack((z1,z1)),T0[:2,:]) )     # eg,dn
                T3=hstack( (T0[2:,:],vstack((z1,z1,z1))) )  # t2g,up
                T4=hstack( (vstack((z1,z1,z1)),T0[2:,:]) )  # t2g,dn
                T = vstack((T1,T2,T3,T4))
            else:
                Z = zeros((len(T0), len(T0)), dtype=complex)
                T = bmat([[T0, Z], [Z, T0]])
        else:
            T = T0
    else:
        raise Exception('ERROR: Transformation `'+transtype+'` not yet implemented.')
    
    # make block matrix, one block for each site of (potentially) cluster problem
    Z = zeros((len(T), len(T)), dtype=complex)
    S = [[Z]*nsites for i in range(nsites)]
    for i in range(nsites):
        S[i][i] = T
    return asarray(bmat(S))

def check_nspins(qsplit, nspins):
    '''Make sure user-specified nspins is compatible with qsplit.'''
    require_twospins = qsplit_table[qsplit][0]
    if require_twospins and nspins != 2:
        raise Exception('ERROR: qsplit = %d requires spin-polarized/spin-orbit calculation.' % qsplit)

def offset_sigind(sigind, offset):
    sigind += where(sigind == 0, 0, offset)

def cmp_sigind_cftrans(indmf, nspins):
    """Loops over all correlated problems and constructs sigind, cftrans and legend.
    nspins can be 1 or 2
    """
    siginds = {}
    cftrans = {}
    legends = {}

    for icix,cix in indmf.cix.items():
        acix = array(cix)
        iatoms = acix[:,0]
        Ls = set(acix[:,1])
        qsplits = set(acix[:,2])
        #iatoms = [iatom for iatom,L,qsplit in cix]
        #Ls = set([L for iatom,L,qsplit in cix])
        #qsplits = set([qsplit for iatom,L,qsplit in cix])
        # currently only support cix where all Ls and qsplits are the same
        if len(Ls) > 1:
            raise Exception('ERROR: each correlated problem must only have single L.')
        if len(qsplits) > 1:
            raise Exception('ERROR: each correlated problem must only have single qsplit.')
        L = Ls.pop()
        qsplit = qsplits.pop()
        nsites = len(iatoms)
        
        check_nspins(qsplit, nspins)
        siginds[icix], legends[icix] = cmp_sigind_legend(qsplit, L, nsites, nspins)
        cftrans[icix] = cmp_cftrans(qsplit, L, nsites, nspins)
    # offset the siginds so that the range spanned by each nonequivalent
    # correlated problem does not overlap with the range of any other ucp
    #
    maxind = {}   # mapping iucp -> max(sigind)
    for igr,cix in indmf.cixgrp.items():
        maxind[igr] = max(siginds[cix[0]].flat)

    iucps = sorted(maxind.keys())
    maxinds = [maxind[iucp] for iucp in iucps]
    offsets = [0] + list(cumsum(maxinds)[:-1])
    offsets_dict = {}  # mapping iucp -> offset
    for iucp,offset in zip(iucps, offsets):
        offsets_dict[iucp] = offset

    for icp,cp in indmf.cix.items():
        offset_sigind(siginds[icp], offsets_dict[indmf.iucps[icp]])
    
    return siginds, cftrans, legends

def CreateIndmflFromIndmf(indmf, options_so, cmpProjector=True, log=sys.stdout):
    inl = indmffile.Indmfl(indmf.case)
    inl.copy_construct(indmf)
    nspins = 2 if options_so else 1
    # for each correlated problem, compute sigind, cftrans and create labels
    inl.siginds, inl.cftrans, inl.legends = cmp_sigind_cftrans(indmf, nspins)
    # projector in case projector is 5
    if cmpProjector:
        inl.CmpProjector(log)
        inl.CmpBandRange(log)
    # write input file for x_dmft.py dmft1
    return inl

    
if __name__ == '__main__':
    # parse command line arguments
    options = parse_cmdline_args()
    w2kenv = W2kEnvironment()

    # automatically detect so/complex runs
    if not options.so and os.path.isfile(w2kenv.case+".inso") and os.path.getsize(w2kenv.case+".inso")>0:
        print('Found '+w2kenv.case+'.inso file, hence assuming so-coupling exists. Applying -so switch!')
        options.so = True

    # parse main input file case.indmf
    try:
        indmf = indmffile.Indmf(w2kenv.case)
        indmf.read()
    except Exception as e:
        print("ERROR: reading file `%s.indmf` failed." % w2kenv.case) 
        print(e)
        exit(1)

    inl = CreateIndmflFromIndmf(indmf, options.so)
    inl.write()
    
    ## instantiate case.indmfl object and initialize with values read from case.indmf
    #inl = indmffile.Indmfl(w2kenv.case)
    #inl.copy_construct(indmf)
    #
    ## for each correlated problem, compute sigind, cftrans and create labels
    #nspins = 2 if options.so else 1
    #
    #inl.siginds, inl.cftrans, inl.legends = cmp_sigind_cftrans(indmf, nspins)
    #
    ##inl.CmpProjector()
    #
    ## write input file for x_dmft.py dmft0
    #inl.write()



