#!/usr/bin/env python
# @Copyright 2007 Kristjan Haule
import utils,indmffile,sys,re,os
import optparse
from scipy import *
import numpy as np
from math import *
from numpy import zeros,array,shape,loadtxt,savetxt,vstack

def create_log_mesh(sigdata, nom, ntail_):
    """Creates logarithmic mesh on Matsubara axis
       Takes first istart points from mesh om and
       the rest of om mesh is replaced by ntail poinst
       redistribued logarithmically.
       Input:
           om      -- original long mesh
           istart  -- first istart points unchanged
           ntail   -- tail replaced by ntail points only
       Output:
           som             -- smaller mesh created from big om mesh
           sSig[nom,nc]    -- Sig on small mesh
       Also computed but not returned:
           ind_om  -- index array which conatins index to
                      kept Matsubara points
    """
    om = sigdata[0]
    
    istart = min(nom, len(om))
    ntail = min(ntail_, len(om)-istart)
        
    istart = min(nom,len(om))
    ntail = min(ntail, len(om)-istart)

    ind_om=[]
    alpha = log((len(om)-1.)/istart)/(ntail-1.)
    for i in range(istart):
        ind_om.append(i)
    for i in range(ntail):
        t = int(istart*exp(alpha*i)+0.5)
        if (t != ind_om[-1]):
            ind_om.append(t)

    ssigdata = zeros( (shape(sigdata)[0], len(ind_om)), dtype=float )
    for i in range(len(ind_om)):
        ssigdata[:,i] = sigdata[:,ind_om[i]]

    return ssigdata

def union(data):
    " Takes a union of array or list"
    c = []
    for d in data:
        if d not in c:
            c.append(d)
    return c

def Find_s_oo_and_Edc(insig):
    # Searching for s_oo and Edc
    with open(insig, 'r') as fh_sig:
        line = next(fh_sig)
        m = re.search(r'(s_oo\s*=\s*\[.*\])', line)
        if m is not None:
            exec(m.group(1), globals())
        else:
            print('Could not find s_oo in line:', line)
        line = next(fh_sig)
        m = re.search(r'(Edc\s*=\s*\[.*\])', line)
        if m is not None:
            exec(m.group(1), globals())
        else:
            print('Could not find Edc in line:', line)
    return (s_oo, Edc)


if __name__=='__main__':
    """ Takes the self-energy file, which contains all self-energies for all impurity problems.
    It splits the self-energies to insert them into the LDA Hamiltonian, necessary for the dmft0, dmft1 and dmft2 step.
    On the imaginary axis, it also creates a logarithmic mesh in Matsubara frequencies,
    to avoid too many points in the tail.
    """
    usage = """usage: %prog [ options ]

    The script takes the self-energy file, which contains all self-energies
    for all impurity problems. It splits them in a way to insert into
    LDA Hamiltonian, necessary for the dmft0, dmft1 and dmft2 step.
    On the imaginary axis, it also creates a logarithmic mesh in Matsubara
    frequencies, to avoid too many points in the tail.

    In this version, one can also insert some fixed self-energies, which are not computed by the impurity, but
    just split/shift certain orbitals. This is replacement for more common open core calculations.
    Provided the users prepares files sfx.[icix] (icix is a number for correlated block index) they will get 
    properly split into sig.inp[icix] and inserted into the hamiltonian.
    """

    parser = optparse.OptionParser(usage)
    parser.add_option("-n", "--nom",  dest="nom",    type="int", default=100, help="Number of matsubara points in log mesh before the tail")
    parser.add_option("-t", "--tail", dest="ntail",  type="int", default=30,  help="Number of matsubara points in the tail of the log mesh")
    parser.add_option("-i", "--sinp", dest="insig",  default='sig.inp', help="filename of the input file (sig.inp)", metavar="FILE")
    parser.add_option("-o", "--sout", dest="outsig", default='sig.inp', help="filename-part of the output file (sig.inp)")
    parser.add_option("-d", "--outDC", dest="outDC", default='Edc.dat', help="filename of the output DC-file (Edc.dat)")
    parser.add_option("-l", "--lext", dest="m_extn", default='', help="For magnetic calculation, it can be 'dn'.")
    # Next, parse the arguments
    (options, args) = parser.parse_args()

    env = utils.W2kEnvironment()
    case = env.case

    print('case=%s, nom=%s, ntail=%s, insig=%s, outsig=%s' %  (case, options.nom, options.ntail, options.insig, options.outsig+'[1-n]'))

    inl = indmffile.Indmfl(case, 'indmfl'+options.m_extn)
    inl.read()

    # Searching for s_oo and Edc
    (s_oo, Edc) = Find_s_oo_and_Edc(options.insig)
    
    #fh_sig = open(options.insig, 'r')
    #m = re.search(r'(s_oo\s*=\s*\[.*\])', fh_sig.next())
    #if m is not None:
    #    exec(m.group(1))
    #m = re.search(r'(Edc\s*=\s*\[.*\])', fh_sig.next())
    #if m is not None:
    #    exec(m.group(1))
    #fh_sig.close()
    
    sigdata = loadtxt(options.insig).transpose()  # self-energy from 'sig.inp' on long mesh
    # fix for the case of a single frequency in the input
    if len(shape(sigdata))==1 : 
        sigdata = array([sigdata]).T
    print('shape(sigdata)=', shape(sigdata))
    print('s_oo=', s_oo)
    print('Edc=', Edc)
    
    savetxt(options.outDC, array(Edc))
    

    if inl.matsubara:
        print('..... Creating logarithmic mesh')
        sigdata = create_log_mesh(sigdata, options.nom, options.ntail)
    
    
    # next few lines only find the columns which are not written in sig.inp
    missing_columns={}
    for icix in list(inl.siginds.keys()):    # over all imp.problems, even those that are equivalent
        
        Sigind = inl.siginds[icix]
        cols_all = union(array(Sigind).flatten())
        cols_all = sorted(cols_all,key=lambda x: abs(x))
        colsp=[x for x in cols_all if x>0]
        colsm=[x for x in cols_all if x<0]
        if colsm:   # missing columns do not appear in self-energy
            for i in colsm:
                missing_columns[abs(i)]=1

    print('..... Going over all correlated blocks')
    for icix in list(inl.siginds.keys()):    # over all imp.problems, even those that are equivalent
        
        Sigind = inl.siginds[icix]
        cols_all = union(array(Sigind).flatten())
        cols_all = sorted(cols_all,key=lambda x: abs(x))
        
        colsp=[x for x in cols_all if x>0]
        colsm=[x for x in cols_all if x<0]
        
        # if we want to treat an entire orbital as open-core, user has to provide files sfx.icix
        # which will be used for such orbitals
        c_filename = 'sfx.'+str(icix)
        open_core = os.path.isfile(c_filename) and len(colsp)==0 and len(colsm)>0
        
        print('icix=', icix, 'colsp=', colsp, 'colsm=', colsm, end=' ')
        if open_core: print(c_filename, end=' ')
        print()
            

        nom=len(sigdata[0])
        wdata=[]
        wdata.append(sigdata[0])
        sinf = [5000.]
        

        if not open_core: 
            # most of the time, we have at least one self-energy entry computed by the impurity problem
            for col in colsp:
                imiss = len([x for x in missing_columns if x<=col])  # imis== how many columns before this column do not appear in sig.inp
                icol = col-imiss # actual column number in sig.inp
                real_part = sigdata[2*icol-1] + s_oo[col-1]-Edc[col-1] # adding  (s_oo-Edc)
                imag_part = sigdata[2*icol]
                wdata.append(real_part)   # corresponds to real-part of this imp. problem
                wdata.append(imag_part)   # corresponds to imag-part of this imp. problem
                sinf.append(s_oo[col-1]-Edc[col-1]); sinf.append(0) # s_oo also needed in k-sum for Eimp
            for col in colsm:
                sooE = s_oo[abs(col)-1]-Edc[abs(col)-1]
                real_part = sooE*ones(nom)
                imag_part = zeros(nom)
                wdata.append(real_part)
                wdata.append(imag_part)
                sinf.append(sooE); sinf.append(0)
        else:
            # if we are treating an entire atomic orbital as an open core
            (c_s_oo,c_Edc) = Find_s_oo_and_Edc(c_filename)
            c_sigdata = loadtxt(c_filename).transpose()            
            if inl.matsubara:
                c_sigdata = create_log_mesh(c_sigdata, options.nom, options.ntail)

            if len(c_sigdata[0])!=nom:
                print('ERROR: The core file '+c_filename+' should have the same mesh as sig.inp. Length of core-file mesh is', len(c_sigdata[0]), 'and the sig mesh is', nom)
            
            cols = abs(array(colsm))
            for col in cols:
                icol = col-min(cols)+1
                real_part = c_sigdata[2*icol-1] + c_s_oo[icol-1]-c_Edc[icol-1] # adding  (s_oo-Edc)
                imag_part = c_sigdata[2*icol]
                wdata.append(real_part)
                wdata.append(imag_part)
                sinf.append(c_s_oo[icol-1]-c_Edc[icol-1]); sinf.append(0)
            
        wdata = array(wdata).transpose()        
        wdata = vstack((wdata, sinf))                   # adding s_oo as the last line
        
        savetxt(options.outsig+str(icix)+options.m_extn, wdata)
