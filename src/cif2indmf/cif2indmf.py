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
              5:["Ta","W","Re","Os","Ir","Pt"],
              6:["Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu"],
              7:["Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr"]}
Us={3: 8., 4: 6.0, 5: 5, 6: 6.0, 7: 5.5}
Js={3: 0.8, 4: 0.8, 5: 0.8, 6: 0.7, 7: 0.7}

def Element_name(name):
    if len(name)>1 and name[1].islower():
        return name[:2]
    else:
        return name[:1]

def Compare_neigbors(dRj,dRs,log):
    N = len(dRj)
    if len(dRs)!=N:
        return False
    for i in range(N):
        diff = sum(abs(dRs-dRj[i]),axis=1)
        iwhich = argmin(diff)
        print('atom['+str(iwhich)+'] -> atom['+str(i)+'] with accuracy=', diff[iwhich], file=log)
        if diff[iwhich]>1e-7:
            return False
    return True
        
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
    name_me = Element_name(strc.aname[jatom])
    name0 = Element_name(strc.neighbrs[jatom][0][1]) # which element?
    possible=[name0]                                 # definitely all atoms of this element
    if name0 in sym_elements: possible.extend(sym_elements[name0])

    print('name0=', name0, 'possible=', possible, 'name_me=', name_me, file=log)
    n=0
    for ngh in strc.neighbrs[jatom]:
        #print('Element_name(ngh[1])=', Element_name(ngh[1]), file=log)
        if Element_name(ngh[1]) not in possible:
            break
        n+=1
    n0=n
    print('n0=', n0, file=log)
    
    n_str=0
    if n>1:
        neighbrs = strc.neighbrs[jatom][:n]
        if n>2: # with three or more vertices I can calculate polyhedron
            neighbrs = strc.neighbrs[jatom][:n]
        else: # we need to take next shell to calculate a polyhedron
            #dst0 = strc.neighbrs[jatom][n-1][0]
            dst1 = strc.neighbrs[jatom][n][0]
            for ngh in strc.neighbrs[jatom]:
                if ngh[0]-dst1>1e-3:
                    break
                n+=1
            if n>6 and name0==name_me:
                neighbrs = strc.neighbrs[jatom][n0:n] # only the next shell
                n_str=n0
            else:
                neighbrs = strc.neighbrs[jatom][:n]
    else: # n==1
        #print('electronegativity=',strc.neighbrs[jatom][0][3])
        if strc.neighbrs[jatom][0][3]>=0:
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
            n_str=1
        else:
            print('WARN: Only one atom appears to be valid neighbor', file=log)
            #dst0 = strc.neighbrs[jatom][n-1][0]
            dst1 = strc.neighbrs[jatom][n][0]
            for ngh in strc.neighbrs[jatom]:
                if ngh[0]-dst1>1e-3:
                    break
                n+=1
            neighbrs = strc.neighbrs[jatom][:n]
    
    print('Analizing', strc.aname[jatom], 'at', strc.pos[first_atom[jatom]], 'with N=', n, file=log)
    print('  with neigbrs:', file=log)
    for i in range(len(neighbrs)):
        print('  n['+str(i+n_str)+']=', neighbrs[i], file=log)
        
    R, N = FindCageBasis(neighbrs, matrix_w2k, log)
    if R is not None:
        return R, (n_str,n_str+N)
    else:
        if N is None:
            return identity(3), (n_str,n_str+1)
        else:
            return identity(3), (n_str,n_str+N)
    

def Cif2Indmf(fcif, input_cor=[3,4,5,6,7], so=False, Qmagnetic=False, DC='exacty', onreal=False,
              fixmu=False,beta=50, QMC_M=5e6, CoulombU=None, CoulombJ=None,
              qsplit=2, Nkp=300, writek=True, logfile='cif2struct.log', keepH=False, ndecimals=3,nradius=4.7,min_Molap=0.7):
    case = os.path.splitext(fcif)[0] # get case

    (strc, lat) = Cif2Struct(fcif, Nkp, writek, Qmagnetic, logfile, cmp_neighbors=True, convertH2R=not keepH, ndecimals=ndecimals,nradius=nradius,min_Molap=min_Molap)

    log = open(logfile, 'a')

    matrix_w2k = lat.rbas_conventional
    #if strc.H2Rtrans:
    #    possible_matrix_w2k = strc.hex2rho @ lat.rbas
    #    
    #    print(('hexagonal setting M=\n'+('{:9.5f} '*3+'\n')*3).format(*ravel(matrix_w2k)), file=log)
    #    print('V_hex=', linalg.det(matrix_w2k), file=log)
    #    print(('rhombohedral settin M=\n'+('{:9.5f} '*3+'\n')*3).format(*ravel(possible_matrix_w2k)), file=log)
    #    print('V_rom=', linalg.det(possible_matrix_w2k), file=log)
    #    matrix_w2k = possible_matrix_w2k
        
    matrix_vesta = strc.Matrix(vesta=True)

    #print('aname=', strc.aname)
    #print('len(strc.neighbrs)=', len(strc.neighbrs))
    
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
        if k not in input_cor: continue
        acor.extend(cor_elements[k])
        for el in cor_elements[k]:
            cor_type[el]=k
    
    first_atom = zeros(len(strc.mult),dtype=int)
    for jatom in range(len(strc.mult)-1):
        first_atom[jatom+1] = first_atom[jatom] + strc.mult[jatom]
    type_jatom=[]
    for jatom in range(strc.nat):
        type_jatom += [jatom]*strc.mult[jatom]
    #type_atom = {}
    #for i in range(len(first_atom)):
    #    type_atom[first_atom[i]]=i
        
    to_frac = linalg.inv(matrix_w2k)
    to_frac_I = matrix_w2k
    vesta_vs_w2k = to_frac @ matrix_vesta
    icix=1
    cixgrp={}
    cix={}
    atoms={}
    Rj = {}
    Nj = {}
    for jatom in range(strc.nat):
        el = Element_name(strc.aname[jatom])
        if el not in acor: continue
        if Qmagnetic:
            Found_equivalent=False
            # check if this atom is related by moment-flipped from another atom
            if first_atom[jatom] in strc.flipped:
                similar = type_jatom[strc.flipped[first_atom[jatom]]]
                if similar in Rj:
                    nnn = Nj[similar] # how many neigbors should we look at?
                    dRj = array([strc.neighbrs[jatom][i][2] for i in range(nnn[0],nnn[1])])   # this atom neigbors
                    dRs = array([strc.neighbrs[similar][i][2] for i in range(nnn[0],nnn[1])]) # similar atom's neighbors, which was already processed
                    print('We are trying to use transformed local axis of',strc.aname[similar], 'for atom', strc.aname[jatom], file=log)
                    print('positions of neigbors around', strc.aname[similar], 'which we are checking:', file=log)
                    for j in range(len(dRs)):
                        print('    [{:10.6f},{:10.6f},{:10.6f}]'.format(*dRs[j,:]), file=log)
                    
                    Found_Rotation=False
                    aRx = strc.rotation_flipped[first_atom[jatom]]
                    aTx = strc.timerevr_flipped[first_atom[jatom]]
                    for ir,Rx in enumerate(aRx):
                        #if aTx[ir]<0:
                        #    print('Changing Rx from', Rx, 'to', -Rx )
                        #    Rx *= -1
                        dRjn = dRj @ Rx.T
                        print('using rotation:',file=log)
                        for ij in range(3):
                            print('    [{:2d},{:2d},{:2d}]'.format(*Rx[ij,:]), file=log)
                        print('transformed positions of neigbors around', strc.aname[jatom], 'which we are checking:', file=log)
                        for j in range(len(dRjn)):
                            print('    [{:10.6f},{:10.6f},{:10.6f}]'.format(*dRjn[j,:]), file=log)
                        
                        Found_Rotation = Compare_neigbors(dRjn,dRs,log)
                        print('Works=', Found_Rotation, file=log)
                        if Found_Rotation:
                            break
                    if Found_Rotation:
                        Rt = linalg.inv(Rx)
                        Rtt = to_frac @ Rt.T @ to_frac_I
                        R = Rj[similar] @ Rtt
                        #R = Rj[similar] @ Rt.T
                        Found_equivalent=True
                        print('Accepting above transformation to transform local axis of '+strc.aname[similar]+' to axis of',
                                  strc.aname[jatom]+' at ',strc.pos[first_atom[jatom]],file=log)
                        print('Using transformation:', file=log)
                        for j in range(3):
                            print('    [{:10.6f},{:10.6f},{:10.6f}]'.format(*Rtt[j,:]), file=log)
            if not Found_equivalent:
                R, nnn = SelectNeighGetPoly(jatom, strc, first_atom, matrix_w2k, log)
        else:
            R, nnn = SelectNeighGetPoly(jatom, strc, first_atom, matrix_w2k, log)
        
        Rj[jatom] = R
        Nj[jatom] = nnn
        locrot_shift = 0 if (R is None or allclose(R, identity(3))) else -1
        
        L = 2 if cor_type[el]<6 else 3
        #qsplit=2 # for now hard-code this qsplit withouth SO
        if qsplit==2 and strc.sgnum>=195:
            # this is cubic system
            qsplit=7
            print('# detected cubic system hence set qplit to',qsplit, file=log)
        
        iatom0 = first_atom[jatom]
        grp=[]
        for i in range(strc.mult[jatom]):
            iatom = iatom0 + i + 1     # corrected Mar.16.2024 iatom->iatom+1
            cix[icix] = [(iatom,L,qsplit)]
            atoms[iatom] = (locrot_shift, R, [], 0)
            grp.append(icix)
            icix += 1
        cixgrp[jatom] = grp
        #cixgrp[iatom0 + 1] = grp
        
        if R is not None:
            Rv = R @ vesta_vs_w2k
            Rf = R @ to_frac
            print('# neigbors used to determine poly:', nnn, file=log)
            print('Rotation to input into case.indmfl by locrot=-1 : ', file=log)
            print(file=log)
            for i in range(3): print( ('{:12.8f} '*3).format(*R[i,:]), file=log)
            print(file=log)
            print('Rotation in fractional coords : ', file=log)
            print(file=log)
            for i in range(3): print( ('{:12.8f} '*3).format(*Rf[i,:]*strc.tobohr), file=log)
            #if strc.H2Rtrans:
            #    Rff = Rf @ strc.hex2rho
            #    print(file=log)
            #    print('Rotation in fractional coords and rhombohedral settings: ', file=log)
            #    print(file=log)
            #    for i in range(3): print( ('{:12.8f} '*3).format(*Rff[i,:]*strc.tobohr), file=log)
            print(file=log)
            print('Rotation for VESTA coords : ', file=log)
            print(file=log)
            for i in range(3): print( ('{:12.8f} '*3).format(*Rv[i,:]), file=log)
            print(file=log)

    matsubara = 0 if onreal else 1
    recomputeEF = 0 if fixmu else 1
    
    head_vars={'hybr_emin':-10., 'hybr_emax':10., 'projector':5, 'matsubara':matsubara, 'broadc':0.025, 'broadnc':0.025, 
                'om_npts':200, 'om_emin':-3, 'om_emax':1}

    indmf = Indmf(case)
    indmf.Initialize_From_Dict(cix, atoms, cixgrp, head_vars)
    indmf.write()
    indmfl = CreateIndmflFromIndmf(indmf, options_so=so, cmpProjector=False, log=log)
    #print(indmfl)
    indmfl.write()

    if Qmagnetic:
        # When moments are canted, it can happen that m1 is mostly opposite to m2 and m2 is mostly opposite to m3,
        # but m1 and m3 appear as quite far away from parallel. In this case we correct the flipping
        # so that m1 and m3 are treated as parallel.
        print('before checking for canted moments type_jatom=', type_jatom, file=log)
        for i1 in strc.flipped:
            for i2 in strc.flipped:
                    j1 = strc.flipped[i1]
                    j2 = strc.flipped[i2]
                    if j1==j2 and type_jatom[i1]!=type_jatom[i2]:
                        print('i1=', i1, 'i2=', i2, 'j1=', j1, 'j2=', j2, 'type_jatom[i1]=', type_jatom[i1], 'type_jatom[i2]=', type_jatom[i2], file=log)
                        type_jatom[i1]=type_jatom[i2] = min(type_jatom[i1],type_jatom[i2])
                    if j1 in strc.flipped:
                        k1 = strc.flipped[j1]
                        if type_jatom[i1]!=type_jatom[k1]:
                            type_jatom[i1]=type_jatom[k1]=min(type_jatom[i1],type_jatom[k1])
                    if j2 in strc.flipped:    
                            k2 = strc.flipped[j2]
                            if type_jatom[i2]!=type_jatom[k2]:
                                type_jatom[i2]=type_jatom[k2]=min(type_jatom[i2],type_jatom[k2])

        print('after checking for canted moments type_jatom=', type_jatom, file=log)
        # With flipped we connected the first atom in supergroup paramegntic structure with all other atoms split into
        # ups and downs in magnetic structure. However, here we want to connect one on one atoms with up and down spins
        # For example, we would have up=[1,2] and down=[3,4], which might show flipped[3]=1 and flipped[4]=1
        # But we would like to have flipped[3]=1 and flipped[4]=2
        remain=[]        
        cix2atom={}
        atom2cix={}
        for icix in indmfl.cix:
            (iatom,L,qsplit) = indmfl.cix[icix][0]
            cix2atom[icix]=iatom
            atom2cix[iatom]=icix
            remain.append(iatom-1)
            
        print('cix2atom=', cix2atom, file=log)
        print('remain=', remain, file=log)

        if qsplit in [-1,1,13,14]:
            # we have spin-orbit, hence there is just one indmfl file, which requires splitting
            print(' before indmfl.cixgrp=', indmfl.cixgrp, file=log)
            for iatm in strc.flipped:
                if iatm not in remain: continue
                remain.remove(iatm)
                jatm = strc.flipped[iatm]
                #print('start with iatm=', iatm, 'jatm=', jatm)
                found=True
                if jatm not in remain:
                    jatm_type = type_jatom[jatm]
                    found=False
                    for jat in remain:
                        if type_jatom[jat]==jatm_type:
                            found=True
                            break
                    if found:
                        jatm = jat
                #print('continue with iatm=', iatm, 'jatm=', jatm, 'found=', found)
                if found:
                    icix1 = atom2cix[iatm+1]
                    icix2 = atom2cix[jatm+1]
                    remain.remove(jatm)
                    icix1, icix2 = sorted([icix1,icix2])
                    iatm, jatm = indmfl.cix[icix1][0][0], indmfl.cix[icix2][0][0]
                    print('flipping[atom='+str(iatm+1)+',icix='+str(icix1)+']=atom='+str(jatm+1)+',icix='+str(icix2), file=log)
                    (iatom,L,qsplit) = indmfl.cix[icix1][0]
                     
                    if qsplit in [-1,1]:
                        #indmfldn.siginds[icix1], indmfldn.siginds[icix2] = indmfl.siginds[icix2], indmfl.siginds[icix1]
                        reorder=ravel([[2*i+1,2*i] for i in range(2*L+1)])
                        print('reorder=', reorder, file=log)
                        
                        # correcting sigind
                        for i in range(len(reorder)):
                            for j in range(len(reorder)):
                                indmfl.siginds[icix2][i,j] = indmfl.siginds[icix1][reorder[i],reorder[j]]
                    else:
                        indmfl.siginds[icix2] = indmfl.siginds[icix1]

                    # correcting indmfl.cixgrp, which contains dictonary indmfl.cixgrp={atom_sort1: [icix1,icix2,...], atom_sort2: [icix3,...]}
                    igrps = list(indmfl.cixgrp.keys())
                    for ir,r in enumerate(igrps):
                        for q in igrps[ir+1:]:
                            #print('r=', r, 'q=', q, 'iatm=', iatm, 'jatm=', jatm)
                            if icix1 in indmfl.cixgrp[r] and icix2 in indmfl.cixgrp[q]:
                                #print('r=', r, 'q=', q, 'indmfl.cixgrp[r]=', indmfl.cixgrp[r], 'indmfl.cixgrp[q]=', indmfl.cixgrp[q])
                                indmfl.cixgrp[r].append(icix2)
                                indmfl.cixgrp[q].remove(icix2)
                            if len(indmfl.cixgrp[q])==0:
                                del indmfl.cixgrp[q]
                    #print('indmfl.cixgrp=', indmfl.cixgrp)
            print(' after indmfl.cixgrp=', indmfl.cixgrp, file=log)
                            
            indmfl.write()
            indmfi = Indmfi(indmfl)
            indmfi.write()
        else:
            # finally go our indmfl file and flip the needed siginds
            indmfldn = indmfl.copy()
            for iatm in strc.flipped:
                if iatm not in remain: continue
                remain.remove(iatm)
                jatm = strc.flipped[iatm]
                #print('start with iatm=', iatm, 'jatm=', jatm)
                found=True
                if jatm not in remain:
                    jatm_type = type_jatom[jatm]
                    found=False
                    for jat in remain:
                        if type_jatom[jat]==jatm_type:
                            found=True
                            break
                    if found:
                        jatm = jat
                #print('continue with iatm=', iatm, 'jatm=', jatm, 'found=', found)
                if found:
                    icix1 = atom2cix[iatm+1]
                    icix2 = atom2cix[jatm+1]
                    remain.remove(jatm)
                    print('flipping[atom='+str(iatm+1)+',icix='+str(icix1)+']=atom='+str(jatm+1)+',icix='+str(icix2), file=log)
                    indmfldn.siginds[icix1], indmfldn.siginds[icix2] = indmfl.siginds[icix2], indmfl.siginds[icix1]
            
            print('remain at the end=', remain, file=log)
            if remain:
                # for those atoms that are not flipped, we need to treat them as different, i.e., different impurity problem
                last_entry_up = max([max(indmfl.siginds[icix].ravel()) for icix in indmfl.cix])
                min_in_remain = min([list(set(indmfldn.siginds[atom2cix[iatm+1]].ravel()))[1] for iatm in remain])
                index_shift = last_entry_up +1 - min_in_remain
                for iatm in remain:
                    icix = atom2cix[iatm+1]
                    indmfldn.siginds[icix] = copy(indmfl.siginds[icix])
                    n = len(indmfldn.siginds[icix])
                    for i in range(n):
                        for j in range(n):
                            if indmfldn.siginds[icix][i,j]!=0:
                                indmfldn.siginds[icix][i,j] += index_shift
                    #print('indmfl.siginds[icix]=  ', indmfl.siginds[icix], file=log)
                    #print('indmfldn.siginds[icix]=', indmfldn.siginds[icix], file=log)
            indmfldn.write(filename=case+'.indmfldn')
            indmfi = Indmfi(indmfl, indmfldn)
            indmfi.write()
    else:
        indmfi = Indmfi(indmfl)
        indmfi.write()

    #print('cixgrp=', cixgrp)
    #print('cix2atom=', indmfi.cix2atom)
    params={}
    params['solver']= ('CTQMC', '# impurity solver')
    params['DCs']= (DC, '# double counting scheme')
    params['max_dmft_iterations'] = (1,   '# number of iteration of the dmft-loop only')
    params['max_lda_iterations']  = (100, '# number of iteration of the LDA-loop only')
    params['finish']              = (50,  '# number of iterations of full charge loop (1 = no charge self-consistency')
    params['ntail']               = (300, '# on imaginary axis, number of points in the tail of the logarithmic mesh')
    params['cc']                  = (1e-5,'# the charge density precision to stop the LDA+DMFT run')
    params['ec']                  = (1e-5,'# the energy precision to stop the LDA+DMFT run')
    params['recomputeEF']         = (recomputeEF,   '# Recompute EF in dmft2 step. If recomputeEF = 2, it tries to find an insulating gap.')
    with open('params.dat','w') as f:
        for k,val in params.items():
            if type(val[0])==str:
                print('{:20s}={:15s}'.format(k,"'"+str(val[0])+"'"), val[1], file=f)
            else:
                print('{:20s}={:15s}'.format(k,str(val[0])), val[1], file=f)
        
        iparams={}
        i_i=0
        #for i,jatom in enumerate(cixgrp):
        for icx,latom in indmfi.cix2atom.items():
            jatom = type_jatom[latom-1]
            #if Qmagnetic and i%2==1:
            #    continue  # only every second impurity problem is needed.
                
            el = Element_name(strc.aname[jatom])
            typ = cor_type[el]

            UCoulomb = Us[typ] if CoulombU is None else CoulombU
            JCoulomb = Js[typ] if CoulombJ is None else CoulombJ
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
            iparams={"exe"                : ["ctqmc"             , "# Name of the executable"],
                     "U"                  : [UCoulomb             , "# Coulomb repulsion (F0) for "+el],
                     "J"                  : [JCoulomb             , "# Hunds coupling  (J) for "+el],
                     "CoulombF"           : ["'Ising'"            , "# Can be set to Full"],
                     "beta"               : [beta                 , "# Inverse temperature"],
                     "svd_lmax"           : [30                   , "# We will use SVD basis to expand G, with this cutoff"],
                     "M"                  : [int(QMC_M)           , "# Total number of Monte Carlo steps"],
                     "mode"               : ["SH"                 , "# We will use self-energy sampling, and Hubbard I tail"],
                     "nom"                : [int(5*beta)          , "# Number of Matsubara frequency points sampled"],
                     "tsample"            : [400                  , "# How often to record measurements"],
                     "GlobalFlip"         : [1000000              , "# How often to try a global flip"],
                     "warmup"             : [int(3e5)             , "# Warmup number of QMC steps"],
                     "nf0"                : [nf0                  , "# Nominal occupancy nd for double-counting"],
                     }
            print(file=f)
            print('# Impurity problem number '+str(i_i),file=f)
            print('iparams'+str(i_i)+'={',file=f)
            for k,val in iparams.items():
                if type(val[0])==str:
                    print('     {:20s}: [{:15s},'.format('"'+k+'"','"'+str(val[0])+'"'), '"'+val[1]+'"],', file=f)
                else:
                    print('     {:20s}: [{:15s},'.format('"'+k+'"',str(val[0])), '"'+val[1]+'"],', file=f)
            print('}',file=f)
            i_i += 1
    
def expand_intlist(input):
    '''Expand out any ranges in user input of integer lists.
    Example: input  = "1,2,4-6"
             output = [1, 2, 4, 5, 6]'''
    from functools import reduce
    import operator
    def parse1(x):
        y = x.split('-')
        return [int(x)] if len(y) == 1 else list(range(int(y[0]), int(y[1])+1))
    return reduce(operator.add, [parse1(x) for x in input.split(',')])

if __name__ == '__main__':
    usage = """usage: %cif2indmf.py [ options ] filename.cif

    Converts cif file (Crystallographic Information File) to input files for eDMFT calculation, 
    case.struct, case.indmfl and case.indmfi. It can also produce high symmetry 
    path in momentum space in case.klist_band.
    """
    parser = optparse.OptionParser(usage)
    parser.add_option('-w', '--writek',  dest='wkp', action='store_true', default=False, help="weather to write case.klist_band high symmetry k-path (default False)")
    parser.add_option('-s', '--so',      dest='so',  action='store_true', default=False, help="is this with spin-orbit or not (default False)")
    parser.add_option('-r', '--real',    dest='real',action='store_true', default=False, help="real axis sets matsubara to 0 in indmffile (default False)")
    parser.add_option('-f', '--fix',     dest='fixmu',action='store_true', default=False, help="fix the chemical potential for Mott insulators (default False)")
    
    parser.add_option('-N', '--Nkp',     dest='Nkp', type='int', default=300, help="number of k-points along the high symmetry path (default 300)")
    parser.add_option('-l', '--log',     dest='log', type='str', default='cif2struct.log', help="info file (default cif2struct.log)")
    parser.add_option('-c', '--cor',     dest='cor', type='str', default='all', help="which atom types are correlated. all or 3,4 or 3-5 or 6-7 (default all)")
    parser.add_option('-d', '--dc',      dest='DC',  type='str', default='exacty', help="The type of double-counting: exact, exacty, exactd, nominal (default exacty)")
    parser.add_option('-U', '--U',       dest='CoulombU', type=float, default=None, help="Coulomb U if we want to set it through input (default None, set inside for each type)")
    parser.add_option('-T', '--T',       dest='Temp', type=float, default=None, help="Temperature in Kelvins (default eV/50)")
    parser.add_option('-J', '--J',       dest='CoulombJ', type=float, default=None, help="Coulomb J if we want to set it through input (default None, set inside for each type)")
    parser.add_option('-b', '--beta',    dest='beta',  type=float, default=50, help="inverse temperature (default 50/eV->232K)")
    parser.add_option('-M', '--QMC_M',   dest='QMC_M', type=float, default=5e6, help="Number of MC steps per core (default = 5M)")
    parser.add_option('-q', '--qsplit',  dest='qsplit', type=int, default=2, help="qsplit (default = 2)")
    parser.add_option('-m', '--magnet',  dest='Qmagnetic', action='store_true', default=False, help="should produce magnetic dft struct that allows broken symmetry from mcif")
    parser.add_option('-H', '--hexagonal',dest='keepH', action='store_true', default=False, help="switches off hexagonal to rhombohedral conversion")
    parser.add_option('-p', '--ndecimals',  dest='ndecimals', type='int', default=3, help="precision when determining the equivalency of atoms we take nn distance with precision ndecimals. Compatible with w2k nn method.")
    parser.add_option('-R', '--nradius',  dest='nradius', type='float', default=4.7, help="How far in AA should nearest neighbors be checked. Compatible with w2k nn method.")
    parser.add_option('-a', '--Molap',   dest='Molap', type='float', default=0.7, help="When checking if spin is flipped we require dot(G*S,S)<-Molap, where G*S is group operation on spin.")
    
    # Next, parse the arguments
    (options, args) = parser.parse_args()
    if len(args)!=1:
        print('Need exactly one argument: the name of cif file')
        sys.exit(1)
    fcif = args[0]
    if options.cor=='all':
        cor=[3,4,5,6,7]
    else:
        cor = expand_intlist(options.cor)
    beta=options.beta
    if options.Temp is not None:
        beta = 11604.5/options.Temp
    Cif2Indmf(fcif, cor, options.so, options.Qmagnetic, options.DC, options.real, options.fixmu, beta, options.QMC_M, options.CoulombU, options.CoulombJ,
                  options.qsplit, Nkp=options.Nkp, writek=options.wkp, logfile=options.log, keepH=options.keepH,
                  ndecimals=options.ndecimals,nradius=options.nradius,min_Molap=options.Molap)

