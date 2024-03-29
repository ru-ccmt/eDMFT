#!/usr/bin/env python
import os, sys
import shutil
from numpy import *
import utils
import indmffile
import find3dRotation as f3dr
import find5dRotation as f5dr
import optparse

if __name__ == '__main__':
    usage = """usage: %findRot.py [ options ]
    
    The script takes the output of dmft1 step, finds hybridization and impurity levels in case.outputdmft1, 
    and diagonalizes it, and modifies case.indmfl file such that the impurity levels with hybridization are diagonal.
    We are diagonalizing E_imp + Delta(omega), which is a frequency dependent object. Delta(omega) vanishes 
    in infinity by definition. 
    
    There are two possible modes:
      We can diagonalize E_imp, which could make Delta(omega==0) somewhat off-diagonal.
      We can force E_imp+Delta(omega==0) to be diagonal, which might make E_imp slighly off diagonal, but hopefully not too bad.
      When E_imp has degeneracies, it is crucail to make E_imp+Delta(omega==0) diagonal, rather than E_imp alone.
    """

    parser = optparse.OptionParser(usage)
    parser.add_option("-i", "--infty", dest="infty", action="store_true", default=False, help="makes E_imp diagonal rather than E_imp+Delta(omega=0), Default=False")
    parser.add_option("-k", "--keep", dest="keep", action="store_true", default=False, help="keep the order of orbitals as close as possible, Default=False")
    
    # Next, parse the arguments
    (options, args) = parser.parse_args()
    
    
    w2k = utils.W2kEnvironment()     # W2k filenames and paths


    outdmft1 = w2k.case+'.outputdmf1'
    if not os.path.isfile(outdmft1):
        print('Can not find file ', outdmft1, 'therefore can not correct local orbitals. Please run "dmft1" first, preferably on the real axis.')
        sys.exit(1)

    # Reads case.indmfl file
    inl = indmffile.Indmfl(w2k.case) # case.indmfl file
    inl.read()                       # case.indmfl read

    f1 = open(outdmft1, 'r')
    Found=False
    for line in f1:
        if line.strip()=='Full matrix of impurity levels follows':
            #print 'Found:', line
            Found=True
            break
    if not Found:
        print('Could not find impurity levels in', outdmft1, '. You should consider to rerun "dmft1", preferably on the real axis.')
        sys.exit(1)

    Hc={}
    for icix in inl.cix:
        line = next(f1).split()
        if int(line[1])!=icix:
            print('It seems that '+w2k.case+'.indmfl and '+outdmft1+' are incompatible...')
            sys.exit(1)
        Sigind = inl.siginds[icix]
        cols = [x for x in Sigind.diagonal() if x!=0]
        dim = len(cols)
        
        Hc[icix] = zeros( (dim,dim), dtype=complex )
        for idim in range(dim):
            line = next(f1)
            if options.infty:
                dat = array(list(map(float,line.split())))
                Hc[icix][idim,:] = dat[::2] + dat[1::2]*1j
            #print 'x: ', line,
        line = next(f1)
        line = next(f1)
        for idim in range(dim):
            line = next(f1)
            if not options.infty:
                dat = array(list(map(float,line.split())))
                Hc[icix][idim,:] = dat[::2] + dat[1::2]*1j
            #print 'y: ', line,
        line = next(f1)
        
        #print Hc[icix]
        #print inl.cftrans[icix]
        l = inl.cix[icix][0][1]
        
        if len(Sigind) == (2*l+1):
            so=False
        elif len(Sigind)==(2*l+1)*2:
            so=True
        else:
            print(':ERROR : Sigind has dimenion', len(Sigind), 'while l=', l)
            sys.exit(1)
        
        hc = 0.5*(Hc[icix]+transpose(conjugate(Hc[icix]))) # Should be Hermitian, so we enforce Hermicity
        T2C = inl.cftrans[icix][:len(hc),:]
        T2Crest = inl.cftrans[icix][len(hc):,:]

        if so:
            T2Cnew = f5dr.GiveNewT2C(hc, T2C, options.keep)
            inext=1
            for i in range(len(Sigind)):
                if Sigind[i,i]>0:
                    Sigind[i,i] = inext;
                    inext += 1
        else:
            T2Cnew = f3dr.GiveNewT2C(hc, T2C)
        
        inl.cftrans[icix][:len(hc),:] = T2Cnew[:,:]
        
    shutil.copy(w2k.case+'.indmfl', w2k.case+'.indmfl_findRot')
    inl.write(w2k.case+'.indmfl')# only_write_stored_data=True)
