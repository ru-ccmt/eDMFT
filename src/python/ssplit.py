#!/usr/bin/env python
# @Copyright 2007 Kristjan Haule
import utils,indmffile,sys,re,os
import optparse
#from scipy import *
import numpy as np
from math import *
from numpy import zeros,ones,array,shape,loadtxt,savetxt,vstack

from pathlib import Path
from matsubara_mesh import DEFAULT_TOL, HEADER, cached_auto_mesh, impurity_nom, log_indices

def create_log_mesh(sigdata, nom, ntail_):
    return sigdata[:,log_indices(sigdata.shape[1],nom,int(ntail_))]

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
    parser.add_option("-n", "--nom",  dest="nom",    type="int", default=None, help="Number of matsubara points in log mesh before the tail")
    parser.add_option("-t", "--tail", dest="ntail",  default="30",  help="Tail points on a logarithmic mesh, or auto for adaptive cubic quadrature")
    parser.add_option("--tail-tol", type="float", default=DEFAULT_TOL, help="Relative integrated error of the scalar probes in auto mode (default %default)")
    parser.add_option("-i", "--sinp", dest="insig",  default='sig.inp', help="filename of the input file (sig.inp)", metavar="FILE")
    parser.add_option("-o", "--sout", dest="outsig", default='sig.inp', help="filename-part of the output file (sig.inp)")
    parser.add_option("-d", "--outDC", dest="outDC", default='Edc.dat', help="filename of the output DC-file (Edc.dat)")
    parser.add_option("-l", "--lext", dest="m_extn", default='', help="For magnetic calculation, it can be 'dn'.")
    # Next, parse the arguments
    (options, args) = parser.parse_args()

    env = utils.W2kEnvironment()
    case = env.case

    print('case=%s, nom=%s, ntail=%s, insig=%s, outsig=%s, extn=%s' %  (case, options.nom, options.ntail, options.insig, options.outsig+'[1-n]', options.m_extn))

    inl = indmffile.Indmfl(case, 'indmfl'+options.m_extn)
    inl.read()

    #print('inl=', inl)
    
    # Searching for s_oo and Edc
    (_s_oo_, _Edc_) = Find_s_oo_and_Edc(options.insig)
    
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
    print('s_oo=', _s_oo_)
    print('Edc=', _Edc_)
    
    

    automatic = inl.matsubara and options.ntail == 'auto'
    nom_input = options.nom if options.nom is not None else 100
    if automatic:
        nom_input = impurity_nom(options.nom)
        full_omega = sigdata[0]
        components = list(sigdata[1::2] + 1j*sigdata[2::2])
        # Include open-core curves for both spin blocks so up and down get
        # identical frequency grids even if their local labels differ.
        broadening = inl.broadc
        mappings = [inl]
        for suffix in ('', 'dn'):
            if suffix != options.m_extn and Path(case+'.indmfl'+suffix).is_file():
                other = indmffile.Indmfl(case,'indmfl'+suffix)
                other.read()
                broadening = max(broadening,other.broadc)
                mappings.append(other)
        cores = {ix for mapping in mappings for ix, labels in mapping.siginds.items()
                 if not np.any(labels>0) and Path('sfx.'+str(ix)).is_file()}
        for ix in sorted(cores):
            core = np.loadtxt('sfx.'+str(ix), ndmin=2).T
            if not np.array_equal(core[0],full_omega):
                raise ValueError('open-core and impurity frequency grids differ')
            components.extend(core[1::2]+1j*core[2::2])
        rows,fit = cached_auto_mesh(full_omega,components,broadening,nom_input,options.tail_tol,
                                    Path('.matsubara_mesh.npz'))
        nlow = min(nom_input,len(full_omega))
        marker = ' '.join(map(str,[HEADER,nlow,len(fit),*fit]))
        print(f'Automatic Matsubara mesh: {nlow} consecutive + {len(rows)-nlow} tail of {len(full_omega)-nlow}; '
              f'probe tolerance {options.tail_tol:g}; {len(fit)} fixed moment-fit samples')
        sigdata = sigdata[:,rows]
    elif inl.matsubara:
        sigdata = create_log_mesh(sigdata,nom_input,int(options.ntail))

    savetxt(options.outDC, array(_Edc_))

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
                real_part = sigdata[2*icol-1] + _s_oo_[col-1]-_Edc_[col-1] # adding  (s_oo-Edc)
                imag_part = sigdata[2*icol]
                wdata.append(real_part)   # corresponds to real-part of this imp. problem
                wdata.append(imag_part)   # corresponds to imag-part of this imp. problem
                sinf.append(_s_oo_[col-1]-_Edc_[col-1]); sinf.append(0) # s_oo also needed in k-sum for Eimp
            for col in colsm:
                sooE = _s_oo_[abs(col)-1]-_Edc_[abs(col)-1]
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
                c_sigdata = c_sigdata[:,rows] if automatic else create_log_mesh(c_sigdata, nom_input, options.ntail)

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
        
        filename = options.outsig+str(icix)+options.m_extn
        if automatic:
            # A trailing comment on the static record travels with each
            # self-energy. Standard numeric readers see the usual columns.
            with open(filename,'w') as out:
                savetxt(out,wdata[:-1])
                out.write(' '.join(f'{v:.18e}' for v in wdata[-1])+' '+marker+'\n')
        else:
            savetxt(filename,wdata)
