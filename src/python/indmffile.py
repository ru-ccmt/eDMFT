#!/usr/bin/env python
# @Copyright 2007 Kristjan Haule
"""
Classes to handle reading/writing of case.indmf* files.

note: self.cps  -> self.cix
      self.ucps -> self.cixgrp
"""
import operator, os, re, sys
from copy import deepcopy
from os.path import isfile
from numpy import array, log, sign, mod, zeros
from functools import reduce
from utils import L2str, L2num

qsplit_doc = """    Qsplit  Description
------  ------------------------------------------------------------
     0  average GF, non-correlated
     1  |j,mj> basis, no symmetry, except time reversal (-jz=jz)
    -1  |j,mj> basis, no symmetry, not even time reversal (-jz=jz)
     2  real harmonics basis, no symmetry, except spin (up=dn)
    -2  real harmonics basis, no symmetry, not even spin (up=dn)
     3  t2g orbitals 
    -3  eg orbitals
     4  |j,mj>, only l-1/2 and l+1/2
     5  axial symmetry in real harmonics
     6  hexagonal symmetry in real harmonics
     7  cubic symmetry in real harmonics
     8  axial symmetry in real harmonics, up different than down
     9  hexagonal symmetry in real harmonics, up different than down
    10  cubic symmetry in real harmonics, up different then down
    11  |j,mj> basis, non-zero off diagonal elements
    12  real harmonics, non-zero off diagonal elements
    13  J_eff=1/2 basis for 5d ions, non-magnetic with symmetry
    14  J_eff=1/2 basis for 5d ions, no symmetry
------  ------------------------------------------------------------"""
projector_doc="""  Projector  Description
------  ------------------------------------------------------------
     1  projection to the solution of Dirac equation (to the head)
     2  projection to the Dirac solution, its energy derivative, 
          LO orbital, as described by P2 in PRB 81, 195107 (2010)
     4  similar to projector-2, but takes fixed number of bands in
          some energy range, even when chemical potential and 
          MT-zero moves (folows band with certain index)
     5  fixed projector, which is written to projectorw.dat. You can
        generate projectorw.dat with the tool wavef.py
     6  similar to projector 5, except projector is here momentum dependent, 
        i.e., similarly to Wannier orbital construction
------  ------------------------------------------------------------
"""

def expand_intlist(input):
    '''Expand out any ranges in user input of integer lists.
    Example: input  = "1,2,4-6"
             output = [1, 2, 4, 5, 6]'''
    def parse1(x):
        y = x.split('-')
        return [int(x)] if len(y) == 1 else list(range(int(y[0]), int(y[1])+1))
    return reduce(operator.add, [parse1(x) for x in input.split(',')])

def divmodulo(x,n):
    "We want to take modulo and divide in fortran way, so that it is compatible with fortran code"
    return ( sign(x)* int(abs(x)/n) , sign(x)*mod(abs(x),n))

class IndmfBase:
    '''Conventions used in naming data structures stored in this class:
    i  = index
    u  = unique
    cp = correlated problem (either single-site or cluster)

    The data structures are dictionaries:
    self.atoms[iatom] = (locrot_shift, new_xyz, shift_vec, [Rmt2])
    self.locrot[iatoms] = (locrot,shift)
    self.cix[icp] = [(iatom_1, L_1, qsplit_1), (iatom_2, L_2, qsplit_2), ...]
    self.cixgrp[iucp] = [icp_1, icp_2, ...]

    Derived classes are responsible for filling in self.cixgrp
    '''
    def __init__(self, case):
        self.case = case
        self.extn = 'indmf'      # derived classes should override this
        self.initvars()
        self.__create_inverses()
        
    def member_vars(self):
        # list of tuples (varname, default value)
        # these are deepcopied when we copy_constryct()
        return [
            ('hybr_emin', -10.0 ),  # (eV) range of hybridization to pass on to impurity problem
            ('hybr_emax',  10.0 ),  # (eV)
            ('emin',       0    ),  # starting band index
            ('emax',       0    ),  # ending band index
            ('Qrenorm',    1    ),  # whether or not to renormalize projector. Should be normally 1.
            ('projector',  5    ),  # type of projection onto correlated space (0,1,2,3,4,5,6)
            ('matsubara',  0    ),  # 0 = real axis, 1 = imaginary axis
            ('broadc',     0.025),  # (eV) broadening for correlated orbitals
            ('broadnc',    0.025),  # (eV) broadening for all orbitals, i.e., noncorrelated
            ('om_npts',    200  ),  # number and range of default omega mesh (if no sig.inp file given)
            ('om_emin',    -3.0 ),  # (eV) for plotting spectral function
            ('om_emax',    1.0  ),  # (eV) for plotting spectral function
            ('broken_sym', 0    ),  # FM, AFM or ferrimagnetic run
            ('atoms',      {}   ),
            ('locrot',     {}   ),
            ('cix',        {}   ), 
            ('cixgrp',       {}   ),
            ('symclasses', {}   ),  # group cps forming each ucp into symmetry classes (e.g. spin-up vs. spin-down)
            ('Lsa',        []   ),  # for wave function calculation: the list of L's for which we need projector
            ('icpsa',      []   ),  # for wave function calculation: the list of cix for which we need projector
            ]

    def initvars(self):
        for attr,val in self.member_vars():
            setattr(self, attr, val)

    def copy_construct(self, c):
        myattr = [attr for attr,val in self.member_vars()]
        for attr in dir(c):
            if attr in myattr:
                setattr(self, attr, deepcopy(getattr(c, attr)))
    
    def __create_inverses(self):
        class Iucps:
            def __getitem__(s,icp):
                return [iucp for iucp,icps in self.cixgrp.items() if icp in icps][0]
        self.iucps = Iucps()

        class Icps:
            def __getitem__(s,iatom):
                return [icp for icp,cp in self.cix.items() if iatom in [iat for iat,L,qsplit in cp]]
        self.icps = Icps()

    def filename(self):
        return self.case + '.' + self.extn

    def file_exists(self):
        return os.path.isfile(self.filename())

    def readlines(self, filename = None):
        fname = filename if filename else self.filename()
        findmf = open(fname, 'r')
        lines = [line.split('#')[0].strip() for line in findmf.readlines()] # strip comments
        findmf.close()
        return (line for line in lines if line)  # strip blank lines & create generator expression

    def parse_head(self, lines):
        self.hybr_emin, self.hybr_emax, self.Qrenorm, self.projector = [float(x) for x in next(lines).split()]
        self.matsubara, self.broadc, self.broadnc, self.om_npts, self.om_emin, self.om_emax = [float(e) for e in next(lines).split()]
        self.matsubara = int(self.matsubara)  # recast these to integers
        self.om_npts   = int(self.om_npts) 
        self.Qrenorm   = int(self.Qrenorm)
        self.projector = int(self.projector)
        if self.projector>=5:
            self.hybr_emin = int(self.hybr_emin)
            self.hybr_emax = int(self.hybr_emax)

    def parse_atomlist(self, lines):
        #self.Lsa=[]
        #self.icpsa=[]
        natom = int(next(lines))
        for i in range(natom):
            dat=next(lines).split()
            iatom, nL, locrot_shift = [int(x) for x in dat[:3]]
            Rmt2=0
            if len(dat)>3:
                Rmt2 = float(dat[3])
            (shift,locrot) = divmodulo(locrot_shift,3)
            if locrot<0:
                if locrot==-2: 
                    locrot=3*nL
                else:
                    locrot=3
            (Ls, qsplits, icps) = (zeros(nL,dtype=int), zeros(nL,dtype=int), zeros(nL,dtype=int))
            for il in range(nL):
                (Ls[il], qsplits[il], icps[il]) = list(map(int, next(lines).split()[:3]))

            self.Lsa.append( Ls )
            self.icpsa.append( icps )
            new_xyz = [[float(x) for x in next(lines).split()] for loro in range(abs(locrot))]
            shift_vec = [float(x) for x in next(lines).split()] if shift else None
            self.locrot[iatom] = (locrot, shift)
            self.atoms[iatom] = (locrot_shift, new_xyz, shift_vec, Rmt2)
            for icp, L, qsplit in zip(icps, Ls, qsplits):
                if icp in self.cix:
                    self.cix[icp] += [(iatom, L, qsplit)]
                else:
                    self.cix[icp] = [(iatom, L, qsplit)]

    def write_head(self, lines):
        lines += [
            ("%f %f %d %d" % (self.hybr_emin, self.hybr_emax, self.Qrenorm, self.projector), "hybridization Emin and Emax, measured from FS, renormalize for interstitials, projection type"),
            ("%1d %g %g %d %f %f" % (self.matsubara, self.broadc, self.broadnc, self.om_npts, self.om_emin, self.om_emax),
             "matsubara, broadening-corr, broadening-noncorr, nomega, omega_min, omega_max (in eV)")
            ]
            
    def write_atomlist(self, lines):
        # create flat list of correlated orbitals (tricky because may have cluster problems)
        corbs = [[(icix,iatom,L,qsplit) for iatom,L,qsplit in v] for icix,v in self.cix.items()]
        corbs = reduce(operator.add, corbs)

        # list of atom-indices of correlated atoms
        icatoms = list(set(iatom for icix,iatom,L,qsplit in corbs))
        icatoms.sort()
        
        lines.append((str(len(icatoms)), "number of correlated atoms"))

        for iatom in icatoms:
            locrot_shift, new_xyz, shift_vec, Rmt2 = self.atoms[iatom]
            orbs = [(icix,L,qsplit) for icix,iat,L,qsplit in corbs if iat==iatom]
            
            #locrot_shift = 3+locrot if shift_vec else locrot
            if Rmt2>0:
                atom_header = ("%-3d %3d %3d  %f" % (iatom, len(orbs), locrot_shift, Rmt2), "iatom, nL, locrot Rmt2")
            else:
                atom_header = ("%-3d %3d %3d" % (iatom, len(orbs), locrot_shift), "iatom, nL, locrot")
            lines.append(atom_header)
            
            for icix,L,qsplit in orbs:
                orbstring = ("%3d %3d %3d" % (L, qsplit, icix), "L, qsplit, cix")
                lines.append(orbstring)

            if locrot_shift:
                locrot_labels = ["new x-axis", "new y-axis", "new z-axis"][:len(new_xyz)]
                for vec,label in zip(new_xyz, locrot_labels):
                    lines.append( ("%11.8f %11.8f %11.8f" % tuple(vec), label) )

            if shift_vec:
                lines.append( ("%4s %4s %4s" % tuple(shift_vec), "real-space shift in atom position") )

    def format(self, lines):
        # merge comments with values
        comment_column = max([len(entry) for entry,comment in lines])
        format = '%-' + str(comment_column) + 's  # %s\n'
        return [format % line for line in lines]

    def writelines(self, text, filename = None):
        fname = filename if filename else self.filename()
        f = open(fname, 'w')
        f.writelines(text)
        f.close()
        return fname


class Indmf(IndmfBase):
    """ It is simplified indmfl file, named case.indmf
    It contains header with two lines containing
       [self.hybr_emin, self.hybr_emax, self.Qrenorm, self.projector]
       [self.matsubara, self.broadc, self.broadnc, self.om_npts, self.om_emin, self.om_emax]
    number of correlated groups
       Next we have one line for each group:
       ind1 cixgrp[ind1]
       ind2 cixgrp[ind2]
    next we have the atom list as in indmfl, which contains
      number of correlated atoms
      iatom len(orbs) locrot_shift  [Rmt2]
         L,qsplit,icix
         new_xyz
         shift

    These data are obtained from the following dictionaries:
      self.cix[icix] = [(iatom_1, L_1, qsplit_1), (iatom_2, L_2, qsplit_2), ...] 
      self.atoms[iatom] = (locrot_shift, new_xyz, shift_vec, [Rmt2])
      self.locrot[iatoms] = (locrot,shift)
      self.cixgrp[iucp] = [icp_1, icp_2, ...]
    """
    def __init__(self, case):
        IndmfBase.__init__(self, case)
        self.extn = 'indmf'
    def read(self):
        lines = self.readlines()
        self.parse_head(lines)
        # read ucp = {cp1, cp2, ...} arrays
        nucp = int(next(lines))
        for i in range(nucp):
            line = [int(x) for x in next(lines).split()]
            iucp = line[0]
            self.cixgrp[iucp] = line[1:]
        self.parse_atomlist(lines)
    
    def write(self):
        lines = []
        self.write_head(lines)
        # write ucp = {cp1, cp2, ...} arrays
        lines.append((str(len(self.cixgrp)), "number of nonequivalent correlated problems"))
        for iucp,icps in self.cixgrp.items():
            entry = ("%3d   " % (iucp,) + ' '.join([str(icp) for icp in icps]), "iucp, cix's")
            lines.append(entry)
        self.write_atomlist(lines)
        text = self.format(lines)
        self.writelines(text)

    def Initialize_From_Dict(self, cix, atoms, cixgrp, head_vars):
      """ Instead of reading the file, the class can be initialized from dictionay.
      It requires the following variables:
         cix[icix]=[(iatom_1, L_1, qsplit_1), (iatom_2, L_2, qsplit_2), ...] 
         atoms[iatom]=(locrot_shift, new_xyz, shift_vec, [Rmt2])
         cixgrp[iucp] = [icp_1, icp_2, ...]
         heav_vars={'hybr_emin':x, 'hybr_emax':x, 'projector':5, 'matsubara':0, 'broadc':0.025, 'broadnc':0.025, 
                    'om_npts':200, 'om_emin':-3, 'om_emax':1}
      """
      self.cix = cix            # self.cix[icix]=[(iatom_1, L_1, qsplit_1), (iatom_2, L_2, qsplit_2), ...] 
      self.atoms = atoms        # self.atoms[iatom]=(locrot_shift, new_xyz, shift_vec, [Rmt2])
      self.cixgrp = cixgrp      # self.cixgrp[iucp] = [icp_1, icp_2, ...]
      # setting self.locrot
      for iatom in self.atoms:
        locrot_shift = self.atoms[iatom][0]
        (shift,locrot) = divmodulo(locrot_shift,3)
        if locrot<0:
          if locrot==-2: 
              locrot=3*nL
          else:
              locrot=3
        self.locrot[iatom] = (locrot,shift)

      myattr = [attr for attr,val in self.member_vars()]
      for attr in head_vars:
        if attr in myattr:
          setattr(self, attr, head_vars[attr])
    
    def user_continue(self, prompt = "Do you want to continue; or edit again? (c/e): "):
        while True:
            userin = input(prompt).strip().lower()
            if userin in ['c', '', 'e']:
                break
            else:
                print('Invalid input.')
        return userin == 'c' or userin == ''

    def orb_strings(self, cp, anames):
        '''Given list of orbitals, creates string with atomname, iatom and L for each orbital.'''
        orbstrings = []
        for iatom,L,qsplit in cp:
            orbstrings.append("%s%d %s" % (anames[iatom], iatom, L2str(L)))
        return orbstrings

    def user_input(self, inpt={}):
        '''Conventions used in this function:
        n = nonequivalent        cp  = correlated problem
        c = correlated           orb = orbital
        i = index (into list)

        The intermediate (temporary) data structures:
        catoms[icatom] = iatom
        corbs[icorb] = (iatom, L)
        qsplits[icorb] = qsplit

        Internal indices run from 0, user input indices run from 1.
        '''
        #import wienfile
        from wstruct import Struct
        self.initvars()  # clear old data (if any)
        #w = wienfile.Struct(self.case)     # parse WIEN2k struct file
        strc = Struct(self.case)
        strc.ReadStruct(self.case+'.struct') # parse WIEN2k struct file
        anames = [None] + strc.flat(strc.aname)  # pad to index from 1; flat list of atom names
        if inpt: print("(ca) ", end=' ')
        print("There are %d atoms in the unit cell:" % sum(strc.mult))
        for i,name in enumerate(anames[1:]):
            print("%3d %s" % (i+1, name))

        while True:
            if inpt:
                if 'ca' in inpt:
                    userin = inpt['ca']
                else:
                    userin = '1'
            else:
                userin = input("Specify correlated atoms (ex: 1-4,7,8): ")

            catoms = expand_intlist(userin)

            print("You have chosen the following atoms to be correlated:")
            for iatom in catoms:
                print("%3d %s" % (iatom, anames[iatom]))

            if inpt: break
            if self.user_continue():
                break

        # currently there's no user interface to input local rotations
        for iatom in catoms:
            locrot_shift = 0
            new_xyz = []
            shift_vec = []
            Rmt2=0
            self.atoms[iatom] = (locrot_shift, new_xyz, shift_vec, Rmt2)
            self.locrot[iatom] = (0,0)
        print()
        while True:
            if inpt: print('(ot) ', end=' ')
            print('For each atom, specify correlated orbital(s) (ex: d,f):')
            corbs = []
            
            if inpt:
                if 'ot' in inpt:
                    user_dat = inpt['ot'].split(',') # should be given as d,d,d for three atoms
                else:
                    user_dat = 'd,'*(len(catoms)-1)+'d'
                    
                if len(user_dat) < len(catoms) :
                    print('ERROR in input : There are '+catoms+' correlated atoms and require the same number of orbital-types. Given input=', user_dat)
                for orb in user_dat:
                    if orb not in ['s','p','d','f']:
                        print('ERROR in input : Correlated orbital type '+orb+' is not allowed. Must be one of s,p,d,f') 
                
            for ii,iatom in enumerate(catoms):
                prompt = "%3d %s: " % (iatom, anames[iatom])
                if inpt:
                    userin = user_dat[ii]
                else:
                    userin = input(prompt)
                for orb in userin.split(','):
                    entry = (iatom, L2num(orb.strip()))
                    corbs.append(entry)

            print("You have chosen to apply correlations to the following orbitals:")
            for icorb, (iatom, L) in enumerate(corbs):
                print("%3d  %s-%d %s" % (icorb+1, anames[iatom], iatom, L2str(L)))

            if inpt : break
            if self.user_continue(): break

        print()
        
        while True:
            if inpt: print('(qs) ', end=' ')
            print("Specify qsplit for each correlated orbital (default = 0):")
            print(qsplit_doc)

            if inpt:
                if 'qs' in inpt:
                    user_dat = inpt['qs'].split(',') # should be given as 2,2,2 for three atoms
                else:
                    user_dat = ['0']*len(catoms)
                if len(user_dat) < len(catoms) :
                    print('ERROR in input : There are '+catoms+' correlated atoms and require the same number of Qsplit entries. Given input=', user_dat)
                    
            qsplits = []
            for icorb, (iatom, L) in enumerate(corbs):
                prompt = "%3d  %s-%d %s: " % (icorb+1, anames[iatom], iatom, L2str(L))
                if inpt:
                    userin = user_dat[icorb]
                else:
                    userin = input(prompt).strip()
                
                qsplit = 0 if userin == '' else int(userin)
                qsplits.append(qsplit)

            print("You have chosen the following qsplits:")
            for icorb, (iatom, L) in enumerate(corbs):
                print("%3d  %s-%d %s: %d" % (icorb+1, anames[iatom], iatom, L2str(L), qsplits[icorb]))
            
            if inpt : break
            if self.user_continue(): break
        
        print()
        while True:
            if inpt: print('(p) ', end=' ')
            print("Specify projector type (default = 5):")
            print(projector_doc, end=' ')
            
            if inpt:
                if 'p' in inpt:
                    userin = inpt['p']
                else:
                    userin = '5'
                self.projector = 5 if userin == '' else int(userin)
                print('> ', self.projector)
            else:
                userin = input("> ").strip()
                self.projector = 5 if userin == '' else int(userin)
                print(self.projector)
            
            if self.projector > 4:
                import glob
                strfile = self.case+'.struct'
                enefiles = glob.glob(self.case+'.energyso')+glob.glob(self.case+'.energyso_'+'*')+glob.glob(self.case+'.energy') + glob.glob(self.case+'.energy_'+'*')
                enefiles = [fil for fil in enefiles if os.path.getsize(fil)>0] # Remove empty files
                if len(enefiles)==0:
                    print('WARNING: Energy files are not present in this directory. Please generate/copy case.energy files here when using projector 5.')
                    print()
            
            if inpt: break
            if self.user_continue(): break
        
        print()

        if (len(corbs)>1): # cluster only if more than one atom correlated
            if inpt:  # non-interactive mode
                if 'cl' in inpt:
                    userin = 'y'
                else:  # default is no cluster-dmft
                    userin = 'n'
                print("(cl) Do you want to group any of these orbitals into cluster-DMFT problems? (y/n): ", userin)
                
            else:     # interactive mode
                userin = input("Do you want to group any of these orbitals into cluster-DMFT problems? (y/n): ").strip().lower()
        else:
            userin = 'n'
            
        if userin == 'y':
            while True:
                #if inpt: '(cl) ',
                print("Enter the orbitals forming each cluster-DMFT problem, separated by spaces")
                #userin = inpt['cl']
                if inpt:  # non-interactive mode
                    if 'cl' in inpt:
                        userin = inpt['cl']
                    else: # default = 1 2 3 4 ...
                        userin = ' '.join([str(icorb+1) for icorb,(iatom,L) in enumerate(corbs)])
                    print(userin)
                else:
                    userin = input("(ex: 1,2 3,4 5-8): ")
                
                expanded = [expand_intlist(group) for group in userin.split()]
                expandedflat = reduce(operator.add, expanded)

                # add orbitals not in CDMFT problems
                icp = 1
                for icorb,(iatom,L) in enumerate(corbs):
                    if icorb+1 not in expandedflat:
                        self.cix[icp] = [(iatom, L, qsplits[icorb])]
                        icp += 1

                # then add orbitals that are part of CDMFT problems
                for group in expanded:
                    self.cix[icp] = []
                    for icorb in group:
                        iatom, L = corbs[icorb-1]
                        self.cix[icp] += [(iatom, L, qsplits[icorb-1])]
                    icp += 1

                print("Your choices give the following correlated problems:")
                for icp,cp in self.cix.items():
                    orbstrings = self.orb_strings(cp, anames)
                    print("%2d  (%s)" % (icp, ', '.join(orbstrings)))

                if inpt: break
                if self.user_continue(): break
        else:
            for icorb,(iatom,L) in enumerate(corbs):
                icp = icorb+1
                self.cix[icp] = [(iatom, L, qsplits[icorb])]
            print()

        if (len(corbs)>1):
            while True:
                if inpt: print('(us) ', end=' ')
                print("Enter the correlated problems forming each unique correlated")
                if inpt:  # non-interactive mode
                    if 'us' in inpt: # non-default value
                        userin = inpt['us']
                    else: # if names of two atoms are the same, default = 1,2 otherwise default = 1 2
                        #   self.cixgrp = { 1: [1], 2: [2], 3: [3],... }
                        userin=''
                        atom_names = [anames[iatom] for icorb,(iatom,L) in enumerate(corbs)]
                        #print 'atom_names=', atom_names
                        # If two atoms have the same name, we choose them to be equivalent.
                        # This might not be the case in general, but the the user should give input
                        userin='1'
                        for i in range(1,len(atom_names)):
                            if atom_names[i] == atom_names[i-1]:
                                userin += ','+str(i+1)
                            else:
                                userin += ' '+str(i+1)
                    print(userin)
                else:
                    userin = input("problem, separated by spaces (ex: 1,3 2,4 5-8): ")
                
                for i,group in enumerate(userin.split()):
                    self.cixgrp[i+1] = expand_intlist(group)
                print()
                print("Each set of equivalent correlated problems are listed below:")
                for iucp,ucp in self.cixgrp.items():
                    cpstrings = ['(%s)' % ', '.join(self.orb_strings(self.cix[icp], anames)) for icp in ucp]
                    print("%3d   %s are equivalent." % (iucp, ' '.join(cpstrings)))

                if inpt: break
                if self.user_continue(): break
                    
                self.cixgrp = {}  # reset
            
            print()
        else:
            self.cixgrp = {1: [1]}

        if inpt: print('(hr) ', end=' ')
        print("Range (in eV) of hybridization taken into account in impurity")

        if inpt:  # non-interactive mode
            print("problems; default %.1f, %.1f: " % (self.hybr_emin, self.hybr_emax))
            if 'hr' in inpt: # non-default value
                userin = inpt['hr']
            else: 
                userin = str(self.hybr_emin) +','+str(self.hybr_emax)
            print(userin)
        else:
            userin = input("problems; default %.1f, %.1f: " % (self.hybr_emin, self.hybr_emax))
            
        if userin.strip():
            self.hybr_emin, self.hybr_emax = [float(e) for e in userin.split(',')]
        else:
            print(self.hybr_emin, self.hybr_emax)
        
        print()

        if inpt:  # non-interactive mode
            print("(a) Perform calculation on real; or imaginary axis? (i/r): (default=i)")
            if 'a' in inpt: # non-default value
                userin = inpt['a']
            else:
                userin = 'i'
            print(userin)
        else:
            userin = input("Perform calculation on real; or imaginary axis? (i/r): (default=i)").strip().lower()
        
        if userin=='r':
            self.matsubara = 0
            if not inpt: print('r')
        else:
            self.matsubara = 1
            if not inpt: print('i')
            
        print()
        #self.matsubara = 1 if userin == 'i' else 0

    def __repr__(self):
        string=self.extn+': case={:s} hybr_emin={:f} hybr_emin={:f} Qrenorm={:d} projector={:d} matsubara={:d}\n'.format(self.case,self.hybr_emin,self.hybr_emax,self.Qrenorm,self.projector,self.matsubara)
        string+='broadc={:f} broadnc={:f} om_npts={:d} om_emin={:f} om_emax={:f} broken_sym={:d}\n'.format(self.broadc,self.broadnc,self.om_npts,self.om_emin,self.om_emax,self.broken_sym)
        string+='atoms: \n'
        for iatom,atm in self.atoms.items():
            string+='iatom={:d} locrot_shift={:d} new_xyz='.format(iatom, atm[0]) + str(atm[1])+' shft_vec='+str(atm[1])+'\n'
        string+='correlatex blocks:\n'
        for icix,cix in self.cix.items():
            string+='cix={:d}: '.format(icix)
            for cc in cix:
                string+='(iatom={:d} l={:d} qsplit={:d}) '.format(cc[0],cc[1],cc[2])
            string+='\n'
        for igrp,grp in self.cixgrp.items():
            string+='igrp={:d} contains cix='.format(igrp)+str(grp)+'\n'
        return string
        

class Indmfl(IndmfBase):
    '''Class for case.indmfl file.
    Additional member variables/data structures:
    self.siginds[icp] = sigind
    self.cftrans[icp] = cftrans
    self.legends[icp] = legends
    EF = fermi level in eV
    '''
    def __init__(self, case, extn='indmfl'):
        IndmfBase.__init__(self, case)
        self.extn = extn #'indmfl'

        # Finding the chemical potential
        EF_exists = os.path.isfile('EF.dat')
        scf2_exists = os.path.isfile(case+".scf2")
        scf2up_exists = os.path.isfile(case+".scf2up")
        scf_exists = os.path.isfile(case+".scf")
        self.EF = None
        
        if EF_exists:
            # The previous DMFT chemical potential
            self.EF = float( open('EF.dat','r').readline() )
            
        if self.EF is None and (scf2_exists or scf2up_exists):
            fname = case+".scf2" if scf2_exists else case+".scf2up"
            fscf = open(fname, 'r')
            lines = fscf.readlines()
            for line in lines:
                if re.match(r':FER', line) is not None:
                    Ry2eV = 13.60569193
                    self.EF = float(line[38:])*Ry2eV
                    break
                
        if self.EF is None and scf_exists:
            fname = case+".scf"
            fscf = open(fname, 'r')
            lines = fscf.readlines()
            for line in lines:
                if re.match(r':FER', line) is not None:
                    Ry2eV = 13.60569193
                    self.EF = float(line[38:])*Ry2eV
                    
        if self.EF is None: self.EF = 0
        
    def member_vars(self):
        myvars = [
            ('siginds', {} ),
            ('cftrans', {} ),
            ('legends', {} ),
            ]
        return IndmfBase.member_vars(self) + myvars

    def read(self, filename = None):
        if filename==None: filename=self.case+'.indmfl'
        lines = self.readlines(filename)
        self.parse_head(lines)
        self.emin, self.emax = self.hybr_emin, self.hybr_emax
        # correcting self.hybr_emin with values in case.indmf, if exists.
        if self.hybr_emin==0 and self.hybr_emax==0:
            if os.path.exists(filename[:-1]):
                dat = open(filename[:-1],'r').readline().split()
                self.hybr_emin, self.hybr_emax = float(dat[0]), float(dat[1])
            else:
                self.hybr_emin, self.hybr_emax = -10,10
        self.parse_atomlist(lines)

        # read the big block of siginds and cftrans
        ncp, maxdim, maxsize = [int(e) for e in next(lines).split()]
        for i in range(ncp):
            icp, dim, size = [int(e) for e in next(lines).split()]
            self.legends[icp] = next(lines).split("'")[1::2]
            self.siginds[icp] = array([[int(e) for e in next(lines).split()] for row in range(dim)])
            raw_cftrans = array([[float(e) for e in next(lines).split()] for row in range(dim)])
            self.cftrans[icp] = raw_cftrans[:,0::2] + raw_cftrans[:,1::2]*1j

    def write_head(self, lines):
        if abs(self.projector)<4:
            styp="%f "
            sdoc = "hybridization Emin and Emax, measured from FS, renormalize for interstitials, projection type"
        else:
            styp="%d "
            sdoc = "hybridization band index nemin and nemax, renormalize for interstitials, projection type"
        lines += [( (styp+styp+"%d %d") % (self.emin, self.emax, self.Qrenorm, self.projector), sdoc),
            ("%1d %g %g %d %f %f" % (self.matsubara, self.broadc, self.broadnc, self.om_npts, self.om_emin, self.om_emax),
             "matsubara, broadening-corr, broadening-noncorr, nomega, omega_min, omega_max (in eV)")]


    def write(self, filename = None): #, only_write_stored_data=False):
        # generate text in two chunks, stored in text and text2
        #   text contains basic information about correlated problems
        #   text2 contains all the siginds, legends and crystal-field transformation matrices
        #self.only_write_stored_data = only_write_stored_data
        lines = []
        
        self.write_head(lines)
        self.write_atomlist(lines)
        text = self.format(lines)

        maxdim = max(len(s) for s in list(self.siginds.values())) # dimension of largest sigind matrix
        sizes = {}
        for icp,sigind in self.siginds.items():
            sizes[icp] = len([x for x in set(sigind.flat) if x != 0])
        maxsize = max(sizes.values())  # number of columns in largest

        text2 = [
            '#================ # Siginds and crystal-field transformations for correlated orbitals ================',
            '%-3d %3d %3d       # Number of independent kcix blocks, max dimension, max num-independent-components' % (len(self.cix), maxdim, maxsize)
            ]

        for icp,sigind in self.siginds.items():
            legend = self.legends[icp]
            cftrans = self.cftrans[icp]

            text2 += [
                "%-3d %3d %3d       # %s" % (icp, len(sigind), sizes[icp], 'cix-num, dimension, num-independent-components'),
                '#---------------- # Independent components are --------------',
                "'%s' "*len(legend) % tuple(legend),
                ]

            text2.append('#---------------- # Sigind follows --------------------------')
            max_sigfig = 1 + int(log(max(sigind.flat))/log(10))
            format = '%' + str(max_sigfig) + 'd'
            for row in sigind:
                text2.append(' '.join([format % elem for elem in row]))

            # print local transform matrix (real & imag)
            text2.append('#---------------- # Transformation matrix follows -----------')
            for row in self.cftrans[icp]:
                text2.append(' '.join(["%11.8f %11.8f  " % (elem.real, elem.imag) for elem in row]))

        # join with first half; add \n to each line in text2
        text += [line+'\n' for line in text2]
        self.writelines(text, filename)
        
    def CmpBandRange(self, log=sys.stdout):
        if abs(self.projector)<4: # This is the old scheme, where hybridization is cut-off by energy
            self.emin = self.hybr_emin+self.EF
            self.emax = self.hybr_emax+self.EF
        else:    # In the new scheme, we cut-off at certain band index
            if not(isfile(self.case+'.energy') or isfile(self.case+'.energyso') or isfile(self.case+'.energyso_1') or isfile(self.case+'.energy_1')):
                print('WARNING: energy files not present. Not calculating projector and band range.', file=log)
                return
            
            import findNbands
            import glob
            
            print('Going over all case.energy files to find which bands are used to construct DMFT projector', file=log)
            strfile = self.case+'.struct'
            enefiles = glob.glob(self.case+'.energyso')+glob.glob(self.case+'.energyso_'+'*')
            if not enefiles:  # Not spin-orbit run
                enefiles = glob.glob(self.case+'.energy') + glob.glob(self.case+'.energy_'+'*')
            enefiles = [fil for fil in enefiles if os.path.getsize(fil)>0] # Remove empty files
            
            if len(enefiles)==0:
                print('all enefiles=', enefiles, file=log)
                print("ERROR : The case.energy* files should be present in this directory when using projector 5 or 6. Exiting....", file=log)
                sys.exit(1)
            
            self.emin,self.emax = findNbands.findNbands(self.hybr_emin+self.EF,self.hybr_emax+self.EF,enefiles,strfile,log=log)
        return (self.emin,self.emax)
    
    def CmpProjector(self, log=sys.stdout):
        if abs(self.projector)>=4: # This is the old scheme, where hybridization is cut-off by energy
            print('Computing DMFT real space projector, which is written in projectorw.dat', file=log)
            import wavef
            atms = list(self.atoms.keys())
            Rmt2 = [atm_dat[3] for iatom,atm_dat in self.atoms.items()] # inl.atoms[iatom] contains (locrot_shift, new_xyz, shift_vec, Rmt2)
            projw = wavef.main(self.case, atms, self.Lsa, self.icpsa, Rmt2, 1.0, 'EF', logf=log)
            #Rm2=[self.atoms[iatom][3] for iatom in list(self.atoms.keys())] 
            #print('Rm2=', Rm2, 'atms=', list(self.atoms.keys()), 'Lsa=', self.Lsa, file=log)
            #wavef.main(self.case, list(self.atoms.keys()), self.Lsa, self.icpsa, Rm2)


    def __repr__(self):
        string=self.extn+': case={:s} hybr_emin={:f} hybr_emin={:f} Qrenorm={:d} projector={:d} matsubara={:d}\n'.format(self.case,self.hybr_emin,self.hybr_emax,self.Qrenorm,self.projector,self.matsubara)
        string+='broadc={:f} broadnc={:f} om_npts={:d} om_emin={:f} om_emax={:f} broken_sym={:d}\n'.format(self.broadc,self.broadnc,self.om_npts,self.om_emin,self.om_emax,self.broken_sym)
        string+='atoms: \n'

        for iatom,atm in self.atoms.items():
            string+='iatom={:d} locrot_shift={:d} '.format(iatom, atm[0])
            if atm[0]!=0:
                string+='new_xyz='+ str(atm[1])
            if atm[2]:
                string+=' shft_vec='+str(atm[2])
            if atm[3]>0:
                string+=' Rmt2='+str(atm[3])
            string += '\n'
        string+='correlatex blocks:\n'
        for icix,cix in self.cix.items():
            string+='cix={:d}: '.format(icix)
            for cc in cix:
                string+='(iatom={:d} l={:d} qsplit={:d}) '.format(cc[0],cc[1],cc[2])
            string+='\n'
        for igrp,grp in self.cixgrp.items():
            string+='igrp={:d} contains cix='.format(igrp)+str(grp)+'\n'
            
        for icix,cix in self.cix.items():
            string+='siginds['+str(icix)+']='+str(self.siginds[icix])+'\n'
            string+='cftrans['+str(icix)+']='+str(self.cftrans[icix])+'\n'
            string+='legends['+str(icix)+']='+str(self.legends[icix])+'\n'
        string+='emin='+str(self.emin)+' emax='+str(self.emax)+'\n'
        if self.emin==0 and self.emax==0:
            string += 'WARNING Projector is not yet set!\n'
        return string

class Indmfi:
    def __init__(self, indmfl):
        self.case = indmfl.case
        self.sigind={}
        #print('cixgrp=', indmfl.cixgrp)
        for icix,cixs in indmfl.cixgrp.items():
            icx = cixs[0]  # just take sigind of first correlated problem  (TODO: modify appropriately for broken symmetry runs)
            self.sigind[icx] = indmfl.siginds[icx]
    def write(self):
        with open(self.case + '.indmfi', 'w') as f:
            print(len(self.sigind), '  # number of sigind blocks', file=f)
            #print('sigind=',self.sigind)
            for icix,sigind in self.sigind.items():
                dim = sigind.shape[0]
                print(sigind.shape[0], '  # dimension of this sigind block', file=f)
                max_sigfig = 1 + int(log(max(sigind.flat))/log(10))
                format = '{:'+str(max_sigfig)+'d} '
                for row in sigind:
                    print((format*dim).format(*row), file=f)
            f.close()
    
def ParsIndmfi(case):
    "Parses the input file which relates the impurity problems with the columns in delta files"
    fi = open(case+'.indmfi')
    nuimpurities = int(next(fi).split()[0])  # number of independent impurity problems
    iSiginds={}
    for icix in range(nuimpurities):
        dim = int(next(fi).split()[0])  # size of the impurity problem
        imp=[]
        for i in range(dim):
            imp.append(list(map(int,next(fi).split()[:dim])))
        iSiginds[icix] = array(imp,dtype=int)
    return iSiginds

