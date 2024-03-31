#!/usr/bin/env python
from numpy import *
import os, sys
import optparse
import utils
import indmffile
from wstruct import Struct
from itertools import chain
import re
from imp2lattc import ImpurityLatticeConnection, SimplifySiginds

if __name__ == '__main__':
    usage = """usage: %init_proj.py [ options ]
    Finishes initialization of DMFT files:
      - creates projectorw.dat file if projector>=4 and projectorw.dat does not yet exists.
      - computes band range in case.indmfl file in band range is set to 0
      - updates nf0 in params file according to DFT density if flag is -a
    """

    parser = optparse.OptionParser(usage)
    parser.add_option("-a", dest="all", action='store_true', default=False, help="corrects nf0 in params file using DFT density matrix")
    # Next, parse the arguments
    (options, args) = parser.parse_args()
    
    w2k = utils.W2kEnvironment()    # W2k filenames and paths
    wopt={'m_extn':''}
    fh_info = sys.stdout
    # Reading 'case.indmfl' file
    inl = indmffile.Indmfl(w2k.case) # case.indmfl file
    inl.read()                       # case.indmfl read
    # If FM or AFM withouth SO, we read also case.indmfldn file
    inldn = None
    if os.path.exists(w2k.case+'.indmfl'+'dn'): wopt['m_extn'] = 'dn'
    if wopt['m_extn']:
        print('INFO: case.indmfldn present => magnetic calculation with two dmft2 steps', file=fh_info)
        inldn = indmffile.Indmfl(w2k.case, 'indmfl'+wopt['m_extn'])
        inldn.read()
    
    ## produce projectorw.dat if it does not yet exist.
    if inl.projector>=4 and not os.path.exists('projectorw.dat'):
        inl.CmpProjector(log=fh_info)
    ## if the band range is not yet specified, do it now
    if inl.emin==inl.emax==0: # we did not yet determine band ranges
        inl.CmpBandRange(log=fh_info)
        inl.write(w2k.case+'.indmfl')
        if wopt['m_extn']: 
            inldn.emin,inldn.emax = inl.emin, inl.emax
            inldn.write(w2k.case+'.indmfl'+wopt['m_extn'])
    if options.all:
        strc = Struct()
        strc.ReadStruct(w2k.case+'.struct')
        isort=[0]
        for iat in range(1,len(strc.mult)+1):
            isort += [iat]*strc.mult[iat-1]
        #print('isort=', isort)

        # reads case.scf2 file and find occupancy
        with open(w2k.case+'.scf2', 'r') as fsc:
            lines = fsc.readlines()

        # for every correlated atom we find corresponding dft occupancy
        ni_cix={}
        for icix in inl.cix:
            (iatom,l,qsplit) = inl.cix[icix][0]
            isrt = isort[iatom]
            if isrt<10:
                cstr='00'+str(isrt)
            elif isrt<100:
                cstr='0'+str(isrt)
            else:
                cstr = str(isrt)
            for line in lines:
                if line[:8]==':QTL'+cstr+':':
                    #ni = float(line.split()[1+l]) # correlated occupancy
                    nn = [float(line[8+7*i:8+7*(i+1)]) for i in range(4)]
                    ni = nn[l]
                    #nf0 = round(ni)
                    break
            ni_cix[icix]=ni
            print('Found dft-occupancy[icix='+str(icix)+']=', ni)
        
        ## Next find impurity-lattice connection
        iSigind = indmffile.ParsIndmfi(w2k.case)
        cols = SimplifySiginds(inl.siginds)
        icols = SimplifySiginds(iSigind)
        imp2latt = ImpurityLatticeConnection( cols, icols, fh_info)
        print('Impurity-lattice connection: imp2latt=', imp2latt, file=fh_info)
        # setting up nominal occupancy nf0[imp]
        nf0={}
        for imp,icix in imp2latt.items():
            ni = [ni_cix[ic] for ic in icix]
            nf0[imp] = round(sum(ni)/len(ni))
            #print(imp,icix, ni, nf0)
        print('Finally nf0=', nf0)
        
        # Now changing params.dat file.
        # First executing it to see if nf0 is already a parameter
        exec(compile(open('params.dat', "rb").read(), 'params.dat', 'exec'))
        # now reading the entire file
        with open('params.dat','r') as fpar:
            lines = fpar.readlines()
        # going over all impurity problems and changing them
        for imp in nf0:
            ipr = 'iparams'+str(imp) # impurity parameters name
            for iline,line in enumerate(lines):
                if re.search(ipr, line) is not None:
                    break    # Found where in params file iparamsxx starts
            if eval( "'nf0' in iparams"+str(imp) ):  # is nf0 already set to some value in nf0?
                for i in range(iline,iline+100):     # then just change the value in this line
                    if re.search('nf0', lines[i]) is not None:
                        print('orig=',lines[i],end='')
                        lines[i] = re.sub('\[\s*(-?\d\.?\d*)','['+str(nf0[imp]),lines[i])
                        print('repl=',lines[i],end='')
                        break
            else:  # nf0 was not yet set, hence adding a line at the beginning
                line = '     "nf0"               : ['+str(nf0[imp])+'              , "# Nominal occupancy nd for double-counting"],\n'
                lines.insert(iline+1,line)
        # Finally writing params.dat back to file
        with open('params.dat','w') as fpar:
            for line in lines:
                fpar.write(line)
