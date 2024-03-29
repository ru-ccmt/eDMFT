#!/usr/bin/env python
""" Module which prepares all necessary files for eDMFT execution. 
  From cif file it produces case.struct file, case.indmf and case.indmfl as well as case.indmfi file.
  It also generates params.dat file.
"""

# @Copyright 2024 Kristjan Haule
import re, sys, os
import optparse
from numpy import *
import numpy as np
import numpy.linalg as linalg
from localaxes import FindCageBasis
from cif2struct import WStruct, Latgen, Cif2Struct
from indmffile import Indmf, Indmfl, Indmfi
from sigen import CreateIndmflFromIndmf
from pymatgen.core.periodic_table import Element


sym_elements={'H': ['D','T'], 'D': ['H','T'], 'T': ['H','D'],
              'Li': ['Na','K'], 'B' : ['Al','Ga'], 'C' : ['Si', 'Ge','Al'],
              'N' : ['P'], 'O':  ['S','Se'], 'F':  ['Cl','Br','I'], 'Ne': ['Ar'],
              'Na': ['K'], 'Mg': ['Ca','Sr'], 'Al': ['B','Ga','Ge','Si','C'], 'Si': ['Si','Ge','Al'],
              'P' : ['N'], 'S' : ['O','Se'], 'Cl': ['F','Br','I'], 'Ar': ['Ne'],
              'K' : ['Na'], 'Ca': ['Mg','Sr'], 'Ti': ['Zr'], 'Co':['Ni'], 'Ni':['Co'],
              'Ga': ['B','Al'],
              'Ge': ['C','Si','Al'], 'As': ['Se','Te'], 'Se': ['As','Te'],
              'Br': ['F','Cl','I'], 'Sr': ['Mg','Ca'], 'Y': ['Sc','La'], 'Zr':['Ti'],
              'Sn':['Pb'], 'Te':['Se','As'], 'I':['Cl','Br','F'], 'La':['Ce'], 'Ce':['La','Pr','Nd'],
              'Pr':['La','Pr','Nd'], 'Nd':['La','Pr','Ce']}

cor_elements={3:["V","Cr","Mn","Fe","Co","Ni","Cu"],
              4:["Nb","Mo","Tc","Ru","Rh","Pd"],
              5:["W","Re","Os","Ir","Pt"],
              6:["Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu"],
              7:["Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr"]}
Us={3: 10., 4: 8.0, 5: 7.0, 6: 6.0, 7: 5.5}
Js={3: 0.8, 4: 0.8, 5: 0.8, 6: 0.7, 7: 0.7}

def Element_name(name):
    if len(name)>1 and name[1].islower():
        return name[:2]
    else:
        return name[:1]

def SelectNeighGetPoly(jatom, strc, first_atom, matrix_w2k, log):
    """This function decides how many available neighbors in strc.neighbrs[jatom]
       will be taken to try and construct a polyhedra out of vertices.
       The key idea if to first try to take the group of atoms of the same element 
       which are uninterapted with any other element in-between. We also allow the group
       to conist of very similar elements, like Ga and Al (listed in sym_elements).
       If this group has 3 or more atoms, we just try to find the polyhedra.
       If the group has only two atoms, we add the next shell of atoms which are at the same distance 
       but further than the first two atoms, and try to construct  polyhedra out of all those atoms.
       If the group has a single atom which is electro-positive, we drop that atom and take 
       the next set of atoms that correspond to certain element until uninterapted.
       The idea is that two TM sometimes find each other very close, but they barely hybridize, hence the
       strong hybridization is determined by the other type of elements.
       If the group has a single atom which is electron-negative, we just add to it the next shell
       of atoms with the same distance.
    """
    name0 = Element_name(strc.neighbrs[jatom][0][1]) # which element?
    possible=[name0]                                 # definitely all atoms of this element
    if name0 in sym_elements: possible.extend(sym_elements[name0])
    n=0
    for ngh in strc.neighbrs[jatom]:
        if Element_name(ngh[1]) not in possible:
            break
        n+=1
    if n>1:
        neighbrs = strc.neighbrs[jatom][:n]
        if n>2: # with three or more vertices I can calculate polyhedron
            neighbrs = strc.neighbrs[jatom][:n]
        else: # we need to take next shell to calculate a polyhedron
            dst0 = strc.neighbrs[jatom][n-1][0]
            dst1 = strc.neighbrs[jatom][n][0]
            for ngh in strc.neighbrs[jatom]:
                if ngh[0]-dst1>1e-3:
                    break
                n+=1
            neighbrs = strc.neighbrs[jatom][:n]
    else: # n==1
        #print('electronegativity=',strc.neighbrs[jatom][0][3])
        if strc.neighbrs[jatom][0][3]>0:
            print('WARNING: At analizing', strc.aname[jatom], 'at', '[{:5.3f},{:5.3f},{:5.3f}]'.format(*strc.pos[first_atom[jatom]]),
                      'we disregard the first neighbor {:s} at [{:5.3f},{:5.3f},{:5.3f}]'.format(strc.neighbrs[jatom][0][1],*strc.neighbrs[jatom][0][2]), file=log)
            name0 = Element_name(strc.neighbrs[jatom][1][1])
            possible=[name0]
            if name0 in sym_elements: possible.extend(sym_elements[name0])
            n=0
            for ngh in strc.neighbrs[jatom][1:]:
                if Element_name(ngh[1]) not in possible:
                    break
                n+=1
            neighbrs = strc.neighbrs[jatom][1:n+1] # here we drop the first atom and take only the next shell.
        else:
            print('WARN: Only one atom appears to be valid neighbor', file=log)
            
            dst0 = strc.neighbrs[jatom][n-1][0]
            dst1 = strc.neighbrs[jatom][n][0]
            for ngh in strc.neighbrs[jatom]:
                if ngh[0]-dst1>1e-3:
                    break
                n+=1
            neighbrs = strc.neighbrs[jatom][:n]
            
    print('Analizing', strc.aname[jatom], 'at', strc.pos[first_atom[jatom]], 'with N=', n, file=log)
    return FindCageBasis(neighbrs, matrix_w2k, log)

def Cif2Indmf(fcif, Nkp=300, writek=True, logfile='cif2struct.log'):
    case = os.path.splitext(fcif)[0] # get case
    
    (strc, lat) = Cif2Struct(fcif, Nkp, writek, logfile, cmp_neighbors=True)

    log = open(logfile, 'a')

    matrix_w2k = lat.rbas_conventional
    matrix_vesta = strc.Matrix(vesta=True)
    
    indx=0
    for jatom in range(len(strc.neighbrs)):
        print('ATOM: {:3d} {:3s} AT {:9.5f} {:9.5f} {:9.5f}'.format(jatom+1,strc.aname[jatom],*strc.pos[indx]),file=log)
        for j,ngh in enumerate(strc.neighbrs[jatom]):
            dRj = ngh[2] @ matrix_w2k
            #dRp = ngh[2] @ matrix_vesta
            dst = linalg.norm(dRj)
            #dst2 = linalg.norm(dRp)
            jname = ngh[1]
            if len(jname)>5: jname=jname[:5]
            print('  ATOM: {:5s} AT {:9.5f} {:9.5f} {:9.5f} = {:9.5f} {:9.5f} {:9.5f} IS AWAY {:13.6f} {:13.6f} ANG'.format(jname,*dRj,*ngh[2],ngh[0],dst/strc.tobohr), file=log)
        indx += strc.mult[jatom]

    
    acor=[]
    cor_type={}
    for k in cor_elements:
        acor.extend(cor_elements[k])
        for el in cor_elements[k]:
            cor_type[el]=k
    
    first_atom = zeros(len(strc.mult),dtype=int)
    for jatom in range(len(strc.mult)-1):
        first_atom[jatom+1] = first_atom[jatom] + strc.mult[jatom]
    
    to_frac = linalg.inv(matrix_w2k)
    vesta_vs_w2k = to_frac @ matrix_vesta
    icix=1
    cixgrp={}
    cix={}
    atoms={}
    for jatom in range(strc.nat):
        el = Element_name(strc.aname[jatom])
        if el not in acor: continue
        #if cor_type[el] != 3: continue
        
        R = SelectNeighGetPoly(jatom, strc, first_atom, matrix_w2k, log)

        locrot_shift = 0 if (R is None or allclose(R, identity(3))) else -1
        
        L = 2 if cor_type[el]<6 else 3
        qsplit=2 # for now hard-code this qsplit withouth SO
        iatom0 = first_atom[jatom]
        grp=[]
        for i in range(strc.mult[jatom]):
            iatom = iatom0 + i + 1     # corrected Mar.16.2024 iatom->itaom+1
            cix[icix] = [(iatom,L,qsplit)]
            atoms[iatom] = (locrot_shift, R, [], 0)
            grp.append(icix)
            icix += 1
        cixgrp[jatom] = grp
        
        if R is not None:
            Rv = R @ vesta_vs_w2k
            Rf = R @ to_frac
        
            print('Rotation to input into case.indmfl by locrot=-1 : ', file=log)
            print(file=log)
            for i in range(3): print( ('{:12.8f} '*3).format(*R[i,:]), file=log)
            print(file=log)
            print('Rotation in fractional coords : ', file=log)
            print(file=log)
            for i in range(3): print( ('{:12.8f} '*3).format(*Rf[i,:]*strc.tobohr), file=log)
            print(file=log)
            print('Rotation for VESTA coords : ', file=log)
            print(file=log)
            for i in range(3): print( ('{:12.8f} '*3).format(*Rv[i,:]), file=log)
            print(file=log)

    head_vars={'hybr_emin':-10., 'hybr_emax':10., 'projector':5, 'matsubara':1, 'broadc':0.025, 'broadnc':0.025, 
                'om_npts':200, 'om_emin':-3, 'om_emax':1}

    indmf = Indmf(case)
    indmf.Initialize_From_Dict(cix, atoms, cixgrp, head_vars)
    indmf.write()
    indmfl = CreateIndmflFromIndmf(indmf, options_so=False, cmpProjector=False, log=log)
    #print(indmfl)
    indmfl.write()

    indmfi = Indmfi(indmfl)
    indmfi.write()

    params={}
    params['solver']= ('CTQMC', '# impurity solver')
    params['DCs']= ('exacty', '# double counting scheme')
    params['max_dmft_iterations'] = (1,   '# number of iteration of the dmft-loop only')
    params['max_lda_iterations']  = (100, '# number of iteration of the LDA-loop only')
    params['finish']              = (50,  '# number of iterations of full charge loop (1 = no charge self-consistency')
    params['ntail']               = (300, '# on imaginary axis, number of points in the tail of the logarithmic mesh')
    params['cc']                  = (5e-5,'# the charge density precision to stop the LDA+DMFT run')
    params['ec']                  = (5e-5,'# the energy precision to stop the LDA+DMFT run')
    params['recomputeEF']         = (1,   '# Recompute EF in dmft2 step. If recomputeEF = 2, it tries to find an insulating gap.')
    with open('params.dat','w') as f:
        for k,val in params.items():
            if type(val[0])==str:
                print('{:20s}={:15s}'.format(k,"'"+str(val[0])+"'"), val[1], file=f)
            else:
                print('{:20s}={:15s}'.format(k,str(val[0])), val[1], file=f)
        
        iparams={}
        for i,jatom in enumerate(cixgrp):
            el = Element_name(strc.aname[jatom])
            typ = cor_type[el]
            # determining starting nf0 from the oxidation state and atomic configuration
            nf0 = 0
            elem = Element(el)
            Ls = 'd' if typ<6 else 'f'
            for r in elem.full_electronic_structure[::-1]:
                if r[1]==Ls:
                    nf_conf = r[2]
                    nf_oxi = strc.oxi_state[jatom]
                    nf0 = r[2] - strc.oxi_state[jatom]
                    break
            print('nf0['+el+']=',nf0,' atomic configuration=',nf_conf,' oxidation=', nf_oxi, file=log)
            iparams={"exe"                : ["ctqmc"            , "# Name of the executable"],
                     "U"                  : [Us[typ]              , "# Coulomb repulsion (F0) for "+el],
                     "J"                  : [Js[typ]              , "# Hunds coupling  (J) for "+el],
                     "CoulombF"           : ["'Ising'"            , "# Can be set to Full"],
                     "beta"               : [50                   , "# Inverse temperature"],
                     "svd_lmax"           : [30                   , "# We will use SVD basis to expand G, with this cutoff"],
                     "M"                  : [5e6                  , "# Total number of Monte Carlo steps"],
                     "mode"               : ["SH"                 , "# We will use self-energy sampling, and Hubbard I tail"],
                     "nom"                : [300                  , "# Number of Matsubara frequency points sampled"],
                     "tsample"            : [200                  , "# How often to record measurements"],
                     "GlobalFlip"         : [1000000              , "# How often to try a global flip"],
                     "warmup"             : [3e5                  , "# Warmup number of QMC steps"],
                     "nf0"                : [nf0                  , "# Nominal occupancy nd for double-counting"],
                     }
            print(file=f)
            print('# Impurity problem number '+str(i),file=f)
            print('iparams'+str(i)+'={',file=f)
            for k,val in iparams.items():
                if type(val[0])==str:
                    print('     {:20s}: [{:15s},'.format('"'+k+'"','"'+str(val[0])+'"'), '"'+val[1]+'"],', file=f)
                else:
                    print('     {:20s}: [{:15s},'.format('"'+k+'"',str(val[0])), '"'+val[1]+'"],', file=f)
            print('}',file=f)
    
if __name__ == '__main__':
    usage = """usage: %cif2indmfl.py [ options ] filename.cif

    Converts cif file (Crystallographic Information File) to input files for eDMFT calculation, 
    case.struct, case.indmfl and case.indmfi. It can also produce high symmetry 
    path in momentum space in case.klist_band.
    """
    parser = optparse.OptionParser(usage)
    parser.add_option('-w', '--writek',  dest='wkp', action='store_true', default=False, help="weather to write case.klist_band high symmetry k-path (default False)")
    parser.add_option('-N', '--Nkp',     dest='Nkp', type='int', default=300, help="number of k-points along the high symmetry path")
    parser.add_option('-l', '--log',     dest='log', type='str', default='cif2struct.log', help="info file")
    # Next, parse the arguments
    (options, args) = parser.parse_args()
    if len(args)!=1:
        print('Need exactly one argument: the name of cif file')
        sys.exit(1)
    fcif = args[0]

    Cif2Indmf(fcif, Nkp=options.Nkp, writek=options.wkp, logfile=options.log)
    

