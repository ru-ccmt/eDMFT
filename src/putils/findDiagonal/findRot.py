#!/usr/bin/env python
import os, sys
import shutil
from numpy import *
import utils
import indmffile
import find3dRotation as f3dr
import find5dRotation as f5dr
import optparse
#from imp2lattc import SimplifySiginds

def Print(S1, log=sys.stdout):
    for i in range(len(S1)):
        for j in range(len(S1[i])):
            if S1.dtype==complex:
                print('{:10.6f} {:10.6f}'.format(S1[i,j].real,S1[i,j].imag), end='  ', file=log)
            else:
                print('{:10.6f}'.format(S1[i,j]), end='  ', file=log)
        print(file=log)


        
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

    cixs = list(inl.cix.keys())
    imax = max(cixs)
    
    inldn=None
    if os.path.isfile(outdmft1+'dn') and os.path.getsize(outdmft1+'dn')>0:
        print('Have magnetic calculation with', outdmft1+'dn')
        inldn = indmffile.Indmfl(w2k.case, 'indmfl'+'dn')
        inldn.read()
        f2 = open(outdmft1+'dn', 'r')
        Found=False
        for line in f2:
            if line.strip()=='Full matrix of impurity levels follows':
                #print 'Found:', line
                Found=True
                break
        if not Found:
            print('Could not find impurity levels in', outdmft1+'dn', '. You should consider to rerun "dmft1 -l dn", preferably on the real axis.')
            sys.exit(1)

        cixs += [ic+imax for ic in inldn.cix.keys() ]
    print(cixs)
    
    with open('findRot.info', 'w') as log:
        print('cixs=', cixs, 'imax=', imax, file=log)
        Hc={}
        T2C={}
        Cols={}
        for ic in cixs:
            if ic<=imax:
                icix = ic
                inl_ = inl
                fi = f1
            else:
                icix = ic-imax
                inl_ = inldn
                fi = f2
            print('ic=', ic, 'icix=', icix, file=log)
            line = next(fi).split()
            if int(line[1])!=icix:
                print('It seems that '+w2k.case+'.indmfl and '+outdmft1+' are incompatible...', file=log)
                print('It seems that '+w2k.case+'.indmfl and '+outdmft1+' are incompatible...', )
                sys.exit(1)
            Sigind = inl_.siginds[icix]
            print('Sigind=', Sigind, file=log)
            Cols[ic] = [x for x in set(Sigind.ravel()) if x!=0]
            dim = len([x for x in Sigind.diagonal() if x!=0])
            
            Hc[ic] = zeros( (dim,dim), dtype=complex )
            for idim in range(dim):
                line = next(fi)
                if options.infty:
                    dat = array(list(map(float,line.split())))
                    Hc[ic][idim,:] = dat[::2] + dat[1::2]*1j
                #print 'x: ', line,
            line = next(fi)
            line = next(fi)
            for idim in range(dim):
                line = next(fi)
                if not options.infty:
                    dat = array(list(map(float,line.split())))
                    Hc[ic][idim,:] = dat[::2] + dat[1::2]*1j
                #print 'y: ', line,
            line = next(fi)

            #print('Hc_read['+str(ic)+']=', file=log)
            #Print(Hc[ic], log)
            #print inl_.cftrans[icix]
            l = inl_.cix[icix][0][1]
            
            if len(Sigind) == (2*l+1):
                so=False
            elif len(Sigind)==(2*l+1)*2:
                so=True
            else:
                print(':ERROR : Sigind has dimenion', len(Sigind), 'while l=', l)
                sys.exit(1)
            
            hc = 0.5*(Hc[ic]+transpose(conjugate(Hc[ic]))) # Should be Hermitian, so we enforce Hermicity
            Hc[ic] = hc
            T2C[ic] = inl_.cftrans[icix][:len(hc),:]
            T2Crest = inl_.cftrans[icix][len(hc):,:]
            
            print('Hc['+str(ic)+']=', file=log)
            Print(Hc[ic], log)
            print('T2C['+str(ic)+']=', file=log)
            Print(T2C[ic], log)

            #if so:
            #    T2Cnew = f5dr.GiveNewT2C(Hc[icix], T2C[icix], options.keep)
            #    inext=1
            #    for i in range(len(Sigind)):
            #        if Sigind[i,i]>0:
            #            Sigind[i,i] = inext;
            #            inext += 1
            #else:
            #    T2Cnew = f3dr.GiveNewT2C(Hc[icix], T2C[icix])
            #
            #inl_.cftrans[icix][:len(hc),:] = T2Cnew[:,:]
            #print('T2Cnew['+str(icix)+']=', file=log)
            #Print(T2Cnew, file=log)
        groups=[]
        igrp={}
        for i,ic in enumerate(cixs):
            if ic not in igrp:
                groups.append([ic])
                igrp[ic]=len(groups)-1
                #print('groups=', groups)
                #print('igrp=', igrp)
            for j in range(i+1,len(cixs)):
                jc = cixs[j]
                if Cols[ic]==Cols[jc]:
                    if jc not in igrp:
                        groups[igrp[ic]].append(jc)
                        igrp[jc]=igrp[ic]
                    #print(ic,jc,'equal', 'groups=', groups, 'igrp=', igrp)
                    
        print('groups of icix=', groups, file=log)
        print(file=log)
        for ie,ics in enumerate(groups):
            hc= zeros(shape(Hc[ics[0]]), dtype=complex)
            for ic in ics:
                hc += Hc[ic]
            hc *= 1/len(ics)

            
            t2c = T2C[ics[0]]
            average=False
            for ic in ics[1:]:
                diff = sum(abs(T2C[ic]-t2c))
                if diff>1e-6:
                    average=True
                    print('It seems T2C['+str(ic)+']!=T2C['+str(ics[0])+'] dif=', diff, 'we will average', file=log)
            if average:
                for ic in ics[1:]:
                    t2c += T2C[ic]
                t2c *= 1/len(ics)
                U,s,V = linalg.svd(t2c)
                t2c = U @ V
            print('---- Diagonalizing ics=', ics, '--------', file=log)
            if so:
                T2Cnew = f5dr.GiveNewT2C(hc, t2c, options.keep)
            else:
                T2Cnew = f3dr.GiveNewT2C(hc, t2c, log)
            
            print('T2Cnew['+str(ics)+']=', file=log)
            Print(T2Cnew, log)
            
            for ic in ics:
                if ic<=imax:
                    icix = ic
                    inl_ = inl
                else:
                    icix = ic-imax
                    inl_ = inldn
                inl_.cftrans[icix][:len(hc),:] = T2Cnew[:,:]
        if True:        
            shutil.copy(w2k.case+'.indmfl', w2k.case+'.indmfl_findRot')
            inl.write(w2k.case+'.indmfl')# only_write_stored_data=True)
            if inldn is not None:
                shutil.copy(w2k.case+'.indmfldn', w2k.case+'.indmfldn_findRot')
                inldn.write(w2k.case+'.indmfldn')# only_write_stored_data=True)
            
