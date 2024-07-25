#!/usr/bin/env python
""" Converts cif file (Crystallographic Information File) to struct file, 
    as needed for wien2k calculation. It can also produce high symmetry 
    path in momentum space in case.klist_band per wien2k needs.
"""
# @Copyright 2007 Kristjan Haule
from math import gcd
import re, sys, os
import optparse
from copy import deepcopy
from numpy import *
import numpy as np
import numpy.linalg as linalg
from pymatgen.io.cif import CifParser
from pymatgen.symmetry.bandstructure import HighSymmKpath
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from fractions import Fraction
from localaxes import *
from functools import cmp_to_key

__author__ = "Kristan Haule"
__copyright__ = "Copyright 2024, eDMFT"
__version__ = "1.0"
__status__ = "Production"
__date__ = "March 6, 2024"

class WStruct:
    def __init__(self, inAngs=True):
        self.tobohr = 1/0.5291772083
        hex2ort = array([[0,1,0],[sqrt(3.)/2.,-0.5,0],[0,0,1]])
        ort2rho = array([[1/sqrt(3),1/sqrt(3),-2/sqrt(3)],[-1,1,0],[1,1,1]])
        self.hex2rho = hex2ort @ ort2rho
        self.rho2hex = linalg.inv(self.hex2rho)
        self.HRtransf = False
        self.inAngs = inAngs
        
    def Matrix(self, vesta=False):
        " Gives matrix of [a,b,c] in pymatgen or vesta convention (not w2k)."
        sin_alpha,cos_alpha= 1.0, 0.0
        if abs(self.alpha-90)>1e-4:
            _alpha_ = self.alpha*pi/180
            sin_alpha, cos_alpha = sin(_alpha_), cos(_alpha_)
        sin_beta, cos_beta = 1.0, 0.0
        if abs(self.beta-90)>1e-4:
            _beta_ = self.beta*pi/180
            sin_beta, cos_beta = sin(_beta_), cos(_beta_)
        if vesta:
            sin_gamma, cos_gamma = 1.0, 0.0
            if abs(self.gamma-90)>1e-4:
                _gamma_ = self.gamma*pi/180
                sin_gamma, cos_gamma = sin(_gamma_), cos(_gamma_)
            
            c1 = self.c * cos_beta
            c2 = (self.c * (cos_alpha - (cos_beta * cos_gamma)))/sin_gamma
            # matrix of [a,b,c].T
            return array([[self.a, 0.0, 0.0],
                          [self.b*cos_gamma, self.b*sin_gamma, 0],
                          [c1, c2, sqrt(self.c**2 - c1**2 - c2**2)]])
        else:
            tgs=1
            if abs(self.gamma-90)>1e-4 or (abs(self.alpha-90)>1e-4 and abs(self.beta-90)>1e-4):
                _alpha_ = self.alpha*pi/180
                _beta_ = self.beta*pi/180
                _gamma_ = self.gamma*pi/180
                tgs = arccos((cos(_alpha_)*cos(_beta_)-cos(_gamma_))/(sin(_alpha_)*sin(_beta_)))/_gamma_
            sin_gamma, cos_gamma = 1.0, 0.0
            if abs(self.gamma-90)>1e-4 or abs(tgs-1)>1e-5:
                _gamma_ = tgs*self.gamma*pi/180
                sin_gamma, cos_gamma = sin(_gamma_), cos(_gamma_)
            # matrix of [a,b,c].T
            return array([[ self.a*sin_beta,           0,                         self.a*cos_beta],
                          [-self.b*sin_alpha*cos_gamma,self.b*sin_alpha*sin_gamma,self.b*cos_alpha],
                          [0, 0, self.c]])
        
    def ScanCif(self, filein, log=sys.stdout, writek=False, Qmagnetic=True, cmp_neighbors=False, convertH2R=True, ndecimals=3, nradius=4.7):
        """ Here is the majority of the algorithm to convert cif2struct.
        We start by allocating CifParser_W2k, which reads cif using paymatgen, and 
        extracts necessary information.
        Below we than arrange this information in the way W2k needs it in struct file.
        
        """
        cif = CifParser_W2k(filein, log, Qmagnetic, cmp_neighbors, ndecimals, nradius)
        
        nsym = len(cif.parser.symmetry_operations)
        self.timat = zeros((nsym,3,3),dtype=int)
        self.tau = zeros((nsym,3))
        glide_plane = zeros(3, dtype=int)
        print('Symmetry operations available', file=log)
        for isym,op in enumerate(cif.parser.symmetry_operations):
            self.timat[isym,:,:] = array(np.round(op.affine_matrix[:3,:3]),dtype=int)
            self.tau[isym,:] = op.affine_matrix[:3,3]
            print(isym,':', file=log)
            for i in range(3):
                print('{:2d}{:2d}{:2d} {:10.8f}'.format(*self.timat[isym,i,:], self.tau[isym,i]), file=log)
            if array_equal(self.timat[isym],identity(3,dtype=int)):
                if sum(abs(self.tau[isym,:]))!=0:
                    gl = array(np.round(self.tau[isym]*2),dtype=int)
                    print('Found glide plane=', gl, file=log)
                    if max(gl)<2:
                        glide_plane = gl
        print('Final glide_plane for CXY/CXZ/CYZ lattice=', glide_plane, file=log)
        
        self.sgnum = int(cif.sgnum)
        self.sgname = re.sub(r'\s+', '', cif.sgname, flags=re.UNICODE)
        self.lattyp = getlattype(self.sgname, self.sgnum, glide_plane)
        self.title = cif.structure.composition.reduced_formula
        self.a, self.b, self.c, self.alpha, self.beta, self.gamma = cif.ca, cif.cb, cif.cc, cif.calpha, cif.cbeta, cif.cgamma
        self.Znuc   = [cif.Z_element[spec] for spec in cif.w2k_coords]

        if self.sgnum in [5,8,9,12,15,20,21,35,36,37,38,39,40,41,63,64,65,66,67,68] and sum(glide_plane)<2:
            print('ERROR glide plane for group', self.sgnum,'should exist, but could not find it', file=log)
            print('All group operations are', file=log)
            for isym,op in enumerate(cif.parser.symmetry_operations):
                timat = array(np.round(op.affine_matrix[:3,:3]),dtype=int)
                tau = op.affine_matrix[:3,3]
                print(isym,':', file=log)
                for i in range(3):
                    print('{:2d}{:2d}{:2d} {:10.8f}'.format(*timat[i,:], tau[i]), file=log)
        
        self.kpath={}
        if writek:
            sga = SpacegroupAnalyzer(cif.structure)
            primitive = sga.get_primitive_standard_structure()
            kpath = HighSymmKpath(primitive)
            #
            #strc1 = cif.structure.copy()
            #print('strc1=', strc1)
            #strc1.remove_site_property("magmom")
            #print('strc2=', strc1)
            #primitive = strc1.to_primitive()
            #print('strc1=', primitive)
            #kpath = HighSymmKpath(primitive)
            #
            self.kpath = kpath._kpath['kpoints']
            
        for i,spec in enumerate(cif.w2k_coords):
            print(i, spec, cif.w2k_coords[spec], file=log)
            
        print('a=', self.a, 'b=', self.b, 'c=', self.c, 'alpha=', self.alpha, 'beta=', self.beta, 'gamma=', self.gamma, file=log)
        print('sgname=', self.sgname, 'sgnum=', self.sgnum, 'lattty=', self.lattyp, 'Nsym=', nsym, file=log)
        print('Znuc=', ['Z['+spec+']='+str(self.Znuc[i]) for i,spec in enumerate(cif.w2k_coords)], file=log)
        if self.kpath:
            print('high symmetry kpath=', file=log)
            for kname in self.kpath:
                print('{:7s}'.format(kname), self.kpath[kname].tolist(), file=log)
        print('positions found in cif file:', file=log)
        for i in range(len(cif.cname)):
            print('name='+cif.cname[i], '[', cif.cx[i], ',', cif.cy[i], ',', cif.cz[i],']', file=log)

        latt_matrix = cif.structure.lattice.matrix
        if self.sgnum in [5,8,9,12,15] and self.lattyp=='B':
            print('WARNING: Transforming B centered monoclinic -> CXY centered monoclinic', file=log)
            print('Transforms from a,b,c=', [self.a,self.b,self.c], 'alpha,beta,gamma=', [self.alpha,self.beta,self.gamma], file=log)
            if abs(self.alpha-90)<1e-4 and abs(self.gamma-90)<1e-4: # I121
                # beta is special
                P = array([[1,0,-1],[0,1,0],[1,0,0]]) # transformation from B to C centering
            elif abs(self.alpha-90)<1e-4 and abs(self.beta-90)<1e-4: # like I112
                # gamma is special
                P = array([[-1,0,0],[-1,0,1],[0,1,0]]) # transformation from B to C centering
            elif abs(self.beta-90)<1e-4 and abs(self.gamma-90)<1e-4:
                # alpha is special
                P = array([[0,1,0],[-1,0,0],[-1,0,1]]) # transformation from B to C centering
            else:
                print('ERROR it seems at least two angles!=90. Dont know how to convert B->C centered monoclinic')
                sys.exit(0)
            Q = linalg.inv(P)
            abcm = self.Matrix() # it should be the same as latt_matrix, hence we would not need to recompute it.
            tabs = P.T @ abcm  # this is transformation of the lattice box, namely, vectors a,b,c
            abc = [linalg.norm(tabs[i,:]) for i in range(3)]  # from matrix we can compute new a,b,c
            self.gamma = arccos(dot(tabs[0,:],tabs[1,:])/(abc[0]*abc[1]))*180/pi # and alpha
            self.alpha = arccos(dot(tabs[1,:],tabs[2,:])/(abc[1]*abc[2]))*180/pi # beta
            self.beta  = arccos(dot(tabs[0,:],tabs[2,:])/(abc[0]*abc[2]))*180/pi # gamma
            self.a, self.b, self.c = abc[0], abc[1], abc[2]
            # now convert all coordinates
            for spec in cif.w2k_coords:
                for i,pos in enumerate(cif.w2k_coords[spec]):
                    cif.w2k_coords[spec][i] = (Q @ pos) % 1.0
            self.lattyp='CXY'
            for kname in self.kpath: # WARNING: This is something I am not sure about. Should we bring them back in 1BZ? Probably not.
                self.kpath[kname] = Q @ self.kpath[kname]
                
            ss='\n    '.join(['{:9.5f} '*3 for i in range(3)])
            print('Mmy='+ss.format(*abcm.ravel()), file=log)
            print('Mpy='+ss.format(*latt_matrix.ravel()), file=log)
            print('trM='+ss.format(*tabs.ravel()), file=log)
            print('abc=', abc, file=log)
            latt_matrix = self.Matrix()
            print('Transforms to a,b,c=', abc, 'alpha,beta,gamma=', [self.alpha,self.beta,self.gamma],'\n', file=log)

            
        print('conventional unit cell (cif.w2k_coords) for finsym: https://stokes.byu.edu/iso/findsym.php:', file=log)
        #for spec in cif.w2k_coords:
        #    for ipos,pos in enumerate(cif.w2k_coords[spec]):
        #        print('name='+spec, 'Z='+str(cif.Z_element[spec]), pos, file=log)
        print('cartesian coordinates of the basis vectors:',file=log)
        print((('{:9.5f} '*3+'\n')*3).format(*ravel(latt_matrix)),file=log)
        print('Number of atoms in the unit cell: ', np.sum([len(cif.w2k_coords[spec]) for spec in cif.w2k_coords]), file=log)
        print('Type of each atom in the unit cell:',' '.join([str(len(cif.w2k_coords[spec]))+'*'+spec+' ' for spec in cif.w2k_coords]),file=log)
        for spec in cif.w2k_coords:
            for ipos,pos in enumerate(cif.w2k_coords[spec]):
                print('{:10.6f} {:10.6f} {:10.6f}'.format(*pos), file=log)
        
        # Correcting for W2k settings, found in cif2struct
        # Note that this Hexagonal versus Rhombohedral settings is still somewhat misterious in w2k, and this algorithm probably does not work in general
        self.H2Rtrans = False
        if convertH2R and (self.sgnum in [146,148,155,160,161,166,167]): # trigonal groups with R : [R3,R-3,R32,R3m,R3c,R-3m,R-3c]
            if abs(self.alpha-90)<1e-5 and abs(self.beta-90)<1e-5 and abs(self.gamma-120)<1e-5: # trigonal, but w2k wants hexagonal settings
                case1 = array([1,2,2])  # (1/3,2/3,2/3)
                case2 = array([1,1,2])  # (1/3,1/3,2/3)
                hexagonal_copies=0
                for isym in range(len(self.tau)):
                    if array_equal(self.timat[isym],identity(3,dtype=int)):
                        if sum(abs(self.tau[isym,:]))!=0:
                            dR = sort(self.tau[isym,:]*3)
                            if sum(abs(dR-case1))<1e-5 or sum(abs(dR-case2))<1e-5:
                                hexagonal_copies += 1
                print('hexagonal_copies=',hexagonal_copies, file=log)
                if hexagonal_copies==2:
                    print('Hexagonal copies of atoms exist, hence eliminating (1/3,2/3,2/3)&(2/3,1/3,1/3) shifts ', file=log)
                    self.H2Rtrans = True
                    self.lattyp = 'R'
        
        RHtransf = False
        RHexeption = False
        if self.sgnum in [146,148,155,160,161,166,167]: # trigonal groups with R : [R3,R-3,R32,R3m,R3c,R-3m,R-3c]
            small = 1e-6
            if abs(self.alpha-90)<small and abs(self.beta-90)<small and abs(self.gamma-120) < small: # trigonal, but w2k wants hexagonal settings
                if not self.H2Rtrans:
                    self.lattyp = 'H'
                    RHexeption = True
            else:
                #if self.lattyp!='R':
                #    print('ERROR: cif.lattyp should be R but it is ', self.lattyp, file=log)
                self.lattyp = 'R'
                RHtransf = True
                aa = self.a*2.0*sin(self.alpha/2*(pi/180))
                bb = aa
                cc = 3*sqrt(self.a**2 - aa**2/3)
                self.a, self.b, self.c = aa, bb, cc
                self.alpha, self.beta, self.gamma = 90., 90., 120.
                print('WARNING: We are changing trigonal group '+str(self.sgnum)+' to equivalent hexagonal settings. However, we choose alpha,beta,gamma=90,90,120, which might not be appropriate.', file=log)
                print('If errors, please correct this algorithm in cif2struct.py', file=log)
                #for spec in cif.w2k_coords:
                #    for ipos,pos in enumerate(cif.w2k_coords[spec]):
                #        pos1 = (pos @ self.rho2hex) % 1.0
                #        print('Changed pos {:10.6f} {:10.6f} {:10.6f} to {:10.6f} {:10.6f} {:10.6f}'.format(*pos,*pos1), file=log)
                #        pos1 = pos
                
                
        if self.sgnum in [143,144,145,147,149,150,151,152,153,154,156,157,158,159,162,163,164,165,168,169,170,171,172,
                              173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194]:
            sseting = 'hexagonal'
            small = 4e-5
            for spec in cif.w2k_coords:
                for ipos,pos in enumerate(cif.w2k_coords[spec]):
                    for fract in [1/6, 2/6, 3/6, 4/6, 5/6]:
                        for i in range(3):
                            if 1e-14 < abs(pos[i]-fract) < small:
                                print(':WARNING: trigonal or hexagonal SG: position for '+spec+':pos['+str(ipos)+','+str(i)+']='+str(pos[i])+'->'+str(fract), file=log)
                                pos[i]=fract
                            if 1e-14 < abs(pos[i]+fract) < small:
                                print(':WARNING: trigonal or hexagonal SG: position for '+spec+':pos['+str(ipos)+','+str(i)+']='+str(pos[i])+'->'+str(-fract), file=log)
                                pos[i]=-fract
                        if pos[0]!=0.0 and pos[1]!=0 and 1e-15<abs(2.0*pos[0]-pos[1])<small:
                            print(':WARNING: trigonal or hexagonal SG: position for '+spec+':pos['+str(ipos)+',:2]='+str(pos[:2])+'->',pos[1]/2.0,pos[1], file=log)
                            pos[0] = pos[1]/2.0
                        if pos[0]!=0.0 and pos[1]!=0 and 1e-15<abs(2.0*pos[1]-pos[0])<small:
                            print(':WARNING: trigonal or hexagonal SG: position for '+spec+':pos['+str(ipos)+',:2]='+str(pos[:2])+'->',pos[0],pos[0]/2., file=log)
                            pos[1] = pos[0]/2.0
                        if pos[0]!=0.0 and pos[1]!=0 and 1e-15<abs(pos[0]+pos[1])<small:
                            print(':WARNING: trigonal or hexagonal SG: position for '+spec+':pos['+str(ipos)+',:2]='+str(pos[:2])+'->',pos[0],-pos[0], file=log)
                            pos[1] = -pos[0]
        
        # Remove atoms that don't belong into primitive unit cell, as per w2k requirement
        vol_ratio = {'B':2, 'F':4, 'CXY':2, 'CXZ':2, 'CYZ':2, 'R':3 }
        if self.lattyp in ['B', 'F', 'CXY', 'CXZ', 'CYZ'] or self.H2Rtrans:
            ratios=[]
            i_w2k = 0
            To_correct={}
            To_correct_indx={}
            indx_kept=[]
            indx_to_remove=[]
            the_same = TheSame(self.lattyp)
            for ispc,spec in enumerate(cif.w2k_coords):
                all_coords = cif.w2k_coords[spec]
                #print('   all_coords['+spec+']=', len(all_coords), file=log)
                #for i in range(len(all_coords)):
                #    print('      ', all_coords[i].tolist(), file=log)
                kept = [0]
                indx_kept.append(i_w2k)
                i_w2k_start = i_w2k
                i_w2k += 1
                for i in range(1,len(all_coords)):
                    Qunique = True
                    for kt in kept:
                        if the_same(all_coords[i], all_coords[kt]):
                            Qunique = False
                            break
                    if Qunique:
                        kept.append(i)
                        indx_kept.append(i_w2k)
                    i_w2k += 1
                to_remove=set(range(len(all_coords)))-set(kept)
                ratio = len(cif.w2k_coords[spec])/len(kept)
                print(spec, ': ', 'out of ', len(cif.w2k_coords[spec]), 'atoms we keep=', len(kept),':', kept, file=log)
                if (fabs(ratio-vol_ratio[self.lattyp])>1e-3):
                    print('ERROR: conventional/primitive volume ratio is', vol_ratio[self.lattyp], 'while # atoms kept ratio=', ratio, file=log)
                    To_correct[spec] = copy(cif.w2k_coords[spec])
                    To_correct_indx[spec]= range(i_w2k_start,i_w2k_start+len(cif.w2k_coords[spec]))
                    _indx_kept_ = array(indx_kept)
                    indx_to_remove += _indx_kept_[ logical_and(_indx_kept_ >= i_w2k_start, _indx_kept_ < i_w2k) ].tolist()
                    print('indx_to_remove=', indx_to_remove, file=log)
                ratios.append(ratio)    
                for i in sorted(to_remove, reverse=True):
                    del cif.w2k_coords[spec][i]
                #print(spec, ':', 'end up with', len(cif.w2k_coords[spec]), file=log)
            #for i in range(len(indx_kept)):
            #    print('indx_kept['+str(i)+']=', indx_kept[i], file=log)
            if sum(abs(array(ratios)-ratios[0]))==0:
                print('All ratios are equal to ', ratios, 'hence structure given in primitive cell. Ignore ERRORS above', file=log)
                To_correct={}
                
            if len(To_correct)>0:
                #------------------------------------------------------------------------------------------------------------
                # The algorithm below is correcting special rare cases in which not enough group operations are removed during
                # construction of subgroups for magnetic up/dn calculation. Namely, when magnetic *.mcif is used, we 
                # produce structure file for up spin and down spin separately. But many group operations of the magetic space
                # group must be removed. In previous algorithm we already removed most such group operations, because they don't
                # preserve spin orientation. However, some group operations might preserve spin orientations, but would still not
                # be part of the space group for separate up/dn calculation. This issue is reflected in the removal of atoms that 
                # are not part of the primitive unit cell. We might keep too many atoms in the primitive cell, because their grouping
                # into equivalent groups was not correctly performed above. This is because too many group operations were left in 
                # cif.parser.symmetry_operations.
                # We attempt to correct for this problem by the following algorithm:
                #   if number of removed atoms going from conventional to primitive unit cell is equal to the volume reduction,
                #       nothing needs to be done
                #   If reduced number of atoms going from conventional to primitive unit cell does not match volume reduction,
                #       we check which atoms in such subgroup are missplaced. This is achieved by checking if all atoms considered
                #       equivalent can be generated by applying all symmetry operations. If an atom can not generate all equivalent
                #       atoms, we try to place it in anothe subgroup of the same element. If alternative placement now allows one to 
                #       generate all equivalent atoms in the alternative subgroup, we place such atom in the alternative subgroups. 
                #       After such reshufling we eliminate again all atoms not in primitive cell, and check if reduced number of atoms 
                #       matches volume reduction of the primitive cell. If yes, we successfuly regrouped atoms. If no, we stop, 
                #       and algorithm should be improved.
                #-------------------------------------------------------------------------------------------------------------
                Exchange=[]
                for spec in To_correct:
                    all_coords = To_correct[spec]
                    print('To_correct['+str(spec)+']=', len(To_correct[spec]), file=log)
                    for i in range(len(all_coords)):
                        print('      {:2d}'.format(i), all_coords[i].tolist(), 'ind=', To_correct_indx[spec][i], file=log)
                    
                    to_correct=[]
                    for k in range(len(all_coords)):
                        equivs=[k]+[i for i in range(len(all_coords)) if i!=k and the_same(all_coords[i], all_coords[k])]
                        print('equivs['+str(k)+']=', equivs, file=log)
                        if len(equivs)!=vol_ratio[self.lattyp]:
                            to_correct.append(k)
                    print('to_correct=', to_correct, file=log)
                    spcs = [spc for spc in To_correct if spc!=spec and Element_name(spc)==Element_name(spec)]
                    print('spcs=', spcs, file=log)
                    Nsym = len(cif.parser.symmetry_operations)
                    for k in to_correct:
                        coord = all_coords[k]
                        is_equivalent={}
                        ncoords = zeros((Nsym,3))
                        for ig,op in enumerate(cif.parser.symmetry_operations): # loop over all symmetry operations
                            ncoord = op.operate(coord) % 1.0       # apply symmetry operation to this site, and produce ncoord
                            ncoords[ig] = ncoord
                            dx = sum(abs(all_coords - ncoord),axis=1)
                            #print('ncoord=', ncoord, 'dx=', dx, file=log)
                            if any(dx<1e-5):
                                ie = argmin(dx)
                                if ie not in is_equivalent:
                                    is_equivalent[ie] = ig
                                #print('is_equivalent=', is_equivalent, file=log)
                                #print('ncoord=', ncoord, '=coord['+str(ie)+']=', all_coords[ie], file=log)
                        print('Checking', spec+'['+str(k)+'] with r=', coord, 'and len(all_coords)=', len(all_coords),
                                  'is_equivalent to how_many=', len(is_equivalent), file=log)
                        #print('  is_equivalent['+str(k)+']=', len(is_equivalent), is_equivalent, file=log)
                        if len(is_equivalent)<len(all_coords):
                            for spc in spcs:
                                is_equivalent2={}
                                other_coords = To_correct[spc]
                                print('  Trying to put '+spec+'['+str(k)+'] with '+spc, 'which has # atoms=', len(other_coords), file=log)
                                for ig in range(Nsym):
                                    dx = sum(abs(other_coords - ncoords[ig]),axis=1)
                                    if any(dx<1e-5):
                                        ie = argmin(dx)
                                        if ie not in is_equivalent2:
                                            is_equivalent2[ie] = ig
                                print('  in '+spc+' is_equivalent to how many=', len(is_equivalent2)+1, file=log)
                                if len(is_equivalent2)==len(other_coords)-1:
                                    print('  It seems '+spec+'['+str(k)+'] fits better with '+spc+' hence moving it', file=log)
                                    Exchange.append( (spec,k,spc,   To_correct_indx[spec][k]) )
                                    break
                print('Exchange=', Exchange, file=log)
                Exchange = sorted(Exchange, key=lambda x: -x[1])
                print('Exchange=', Exchange, file=log)
                
                To_correct2={}
                for spec in To_correct:
                    all_coords = To_correct[spec]
                    To_correct2[spec]=[[] for i in range(len(all_coords))]
                    for i in range(len(all_coords)):
                        To_correct2[spec][i] = (To_correct[spec][i], To_correct_indx[spec][i])
                
                for spc1,i,spc2,ii in Exchange:
                    To_correct2[spc2].append(To_correct2[spc1][i])
                    del To_correct2[spc1][i]
                
                for spec in To_correct2:
                    print('To_correct['+str(spec)+']=', len(To_correct2[spec]), file=log)
                    for i in range(len(To_correct2[spec])):
                        print('      {:2d}'.format(i), To_correct2[spec][i][0].tolist(), To_correct2[spec][i][1], file=log)

                for spec in To_correct2:
                    all_coords = To_correct2[spec]
                    kept = [0]
                    for i in range(1,len(all_coords)):
                        Qunique = True
                        for kt in kept:
                            if the_same(all_coords[i][0], all_coords[kt][0]):
                                Qunique = False
                                break
                        if Qunique:
                            kept.append(i)
                    to_remove=set(range(len(all_coords)))-set(kept)
                    ratio = len(To_correct2[spec])/len(kept)
                    print(spec, ': ', 'out of ', len(To_correct2[spec]), 'atoms we keep=', len(kept),':', kept, file=log)
                    if (fabs(ratio-vol_ratio[self.lattyp])>1e-3):
                        print('ERROR: conventional/primitive volume ratio is', vol_ratio[self.lattyp], 'while # atoms kept ratio=', ratio, file=log)
                        sys.exit(1)
                    for i in sorted(to_remove, reverse=True):
                        del To_correct2[spec][i]

                print('before modify: indx_kept=', indx_kept, file=log)
                indx_where_to_insert = min([indx_kept.index(i) for i in indx_to_remove])
                print('indx_where_to_insert=', indx_where_to_insert, file=log)
                indx_kept = [i for i in indx_kept if i not in indx_to_remove]
                indx_kept_extention=[]
                for spec in To_correct2:
                    print('To_correct['+str(spec)+']=', len(To_correct2[spec]), file=log)
                    for i in range(len(To_correct2[spec])):
                        print('      {:2d}'.format(i), To_correct2[spec][i][0].tolist(), To_correct2[spec][i][1], file=log)
                    cif.w2k_coords[spec] = array([To_correct2[spec][i][0] for i in range(len(To_correct2[spec]))])
                    indx_kept_extention += [To_correct2[spec][i][1] for i in range(len(To_correct2[spec]))]
                indx_kept = indx_kept[:indx_where_to_insert] + indx_kept_extention + indx_kept[indx_where_to_insert:]
                print('indx_kept=', indx_kept, file=log)
                #-----------------------------------------------------------------------------------------------------------------------------
                        
            in_primitive={}
            for i in range(len(indx_kept)):
                in_primitive[indx_kept[i]]=i
            for i in in_primitive:
                print('in_primitive['+str(i)+']=', in_primitive[i], file=log)
        else:
            Nall = sum([len(cif.w2k_coords[spec]) for spec in cif.w2k_coords])
            in_primitive = range(Nall)

        print('After removing atoms not in primitive cell, we have the following cif.w2k_coords list:', file=log)
        for spec in cif.w2k_coords:
            for ipos,pos in enumerate(cif.w2k_coords[spec]):
                print('name='+spec, 'Z='+str(cif.Z_element[spec]), pos, file=log)
                
        if self.H2Rtrans:
            for spec in cif.w2k_coords:
                for i,pos in enumerate(cif.w2k_coords[spec]):
                    cif.w2k_coords[spec][i] = dot(pos, self.hex2rho) % 1.0
            # Converting symmetry operations from hexagonal to rhombohedral setting
            for isym in range(len(self.tau)):
                #print('operation before trans=', self.timat[isym,:], self.tau[isym,:], file=log)
                self.timat[isym,:,:] = np.round( self.hex2rho.T @ self.timat[isym,:,:] @ self.rho2hex.T )
                self.tau[isym,:] = (self.tau[isym,:] @ self.hex2rho) % 1.0
                #print('operation after  trans=', self.timat[isym,:], self.tau[isym,:], file=log)
            # Take only unique operations
            timat_new = [self.timat[0,:,:]]
            tau_new = [self.tau[0,:]]
            for isym in range(1,len(self.tau)):
                Qunique=True
                for j in range(len(tau_new)):
                    if sum(abs(self.timat[isym,:,:]-timat_new[j]))<1e-7 and sum(abs(self.tau[isym,:]-tau_new[j]))<1e-7:
                        Qunique = False
                if Qunique:
                    timat_new.append(self.timat[isym,:,:])
                    tau_new.append(self.tau[isym,:])
            if len(tau_new)<len(self.tau):
                self.timat = array(timat_new)
                self.tau = array(tau_new)
            
            print('After converting to rombohedral setting, using', file=log)
            for i in range(3):
                print(('{:6.2f}'*3).format(*self.hex2rho[i,:]), file=log)
            print('  we have the following atoms in primitive cell are:', file=log)
            for spec in cif.w2k_coords:
                for ipos,pos in enumerate(cif.w2k_coords[spec]):
                    print('name='+spec, 'Z='+str(cif.Z_element[spec]), pos, file=log)

            print('Unique:', file=log)
            for isym in range(len(tau_new)):
                print('isym=', isym, file=log)
                for i in range(3):
                    print(('{:4.0f}'*3).format(*timat_new[isym][i]), '  {:10.5g}'.format(tau_new[isym][i]), file=log)
                
        self.nat = len(cif.w2k_coords)
        self.mult = [len(cif.w2k_coords[spec]) for spec in cif.w2k_coords]
        self.aname = [spec for spec in cif.w2k_coords]
        self.aname = correct_element_names(self.aname)
        self.pos=[]
        for spec in cif.w2k_coords:
            self.pos.extend(cif.w2k_coords[spec])
        self.pos = array(self.pos)
        
        self.iatnr  = [-(i+1) for i,spec in enumerate(cif.w2k_coords)]
        self.isplit = [15  for spec in cif.w2k_coords]
        self.mode = 'RELA'
        self.jrj    = [781 for spec in cif.w2k_coords]
        self.r0     = [Z2r0(Z) for Z in self.Znuc]
        self.rotloc = [identity(3) for spec in cif.w2k_coords]

        if cmp_neighbors:
            self.rmt = cif.Rmt*self.tobohr
            self.oxi_state = cif.oxi_state
        else:
            self.rmt    = [2.0 for spec in cif.w2k_coords]
            
            
        print(('\npymatgen_conventional=\n'+('\n'.join(['{:9.5f} '*3+' ==a'+str(i+1) for i in range(3)]))).format(
            *ravel(latt_matrix*self.tobohr)), file=log)
        print('paymatgen_Unit cell volume=', linalg.det(latt_matrix)*self.tobohr**3, file=log)
        print('\n', file=log)
        #print('nat=', self.nat, file=log)
        #print('mult=', self.mult, file=log)
        #print('aname=', self.aname, file=log)
        #print('Znuc=', self.Znuc, file=log)
        #print('iatnr=', self.iatnr, file=log)
        #print('isplit=', self.isplit, file=log)
        self.neighbrs = cif.neighbrs

        #if self.lattyp in ['CXY','CXZ','CYZ','F','B']:
        #    # for these Bravais lattices the glide planes should be eliminated from w2k symmetry operations
        #    # because w2k adds them in by hand. Hence it is better to generate symmetry operations by w2k initialization
        #    # then remove half of the symmetry operations here.
        #    self.tau=[]
        if self.lattyp in ['B', 'F', 'CXY', 'CXZ', 'CYZ']:
            #print('Here we have lattyp=', self.lattyp)
            if self.lattyp=='B':
                Ltrans = array([[0.5,0.5,0.5]])
            elif self.lattyp=='F':
                Ltrans= array([[0,0.5,0.5],[0.5,0,0.5],[0.5,0.5,0]])
            elif self.lattyp=='CXY':
                Ltrans=array([[0.5,0.5,0]])
            elif self.lattyp=='CXZ':
                Ltrans=array([[0.5,0,0.5]])
            elif self.lattyp=='CYZ':
                Ltrans=array([[0,0.5,0.5]])
            else:
                print('ERROR: Not implemented Ltrans')
            timat=[]
            tau=[]
            for isym in range(len(self.tau)):
                timat.append(self.timat[isym,:,:])
                tau.append(self.tau[isym,:])
            timat_keep=[]
            tau_keep=[]
            while len(tau)>0:
                imin = argmin([sum(abs(tau[isym])) for isym in range(len(tau))])
                T0=timat[imin]
                t0=tau[imin]
                #print('t0=', t0, 'T0=', T0, file=log)
                remove=[]
                for isym in range(len(tau)):
                    if allclose(timat[isym],T0):
                        dt = (tau[isym]-t0)%1.0
                        d0 = sum(abs(dt))
                        diff = sum(abs(Ltrans-dt),axis=1)
                        if d0<1e-7:
                            timat_keep.append(timat[isym])
                            tau_keep.append(tau[isym])
                            remove.append(isym)
                        elif any(diff<1e-7):
                            remove.append(isym)
                        
                #print('Symm operations remove=', remove, file=log)
                for i in sorted(remove,reverse=True):
                    del(tau[i])
                    del(timat[i])
            print('Symmetry operations keept after modifying for this Bravais lattice:', file=log)
            for i in range(len(tau_keep)):
                for j in range(3):
                    print('{:2d}{:2d}{:2d} {:10.8f}'.format(*timat_keep[i][j,:], tau_keep[i][j]), file=log)
                print('      {:2d}'.format(i+1), file=log)
            self.timat=array(timat_keep)
            self.tau=array(tau_keep)

        self.flipped=None
        if Qmagnetic:
            if self.lattyp in ['B', 'F', 'CXY', 'CXZ', 'CYZ'] or self.H2Rtrans:
                self.flipped = {}
                self.rotation_flipped = {}
                self.timerevr_flipped = {}
                for ii in cif.flipped:
                    if ii in in_primitive and cif.flipped[ii] in in_primitive:
                        jj = in_primitive[ii]
                        self.flipped[jj] = in_primitive[cif.flipped[ii]]
                        self.rotation_flipped[jj] = cif.rotation_flipped[ii]
                        self.timerevr_flipped[jj] = cif.timerevr_flipped[ii]
                print('str.flipped:', file=log)
                for ii in self.flipped:
                    print(str(ii)+'->'+str(self.flipped[ii]), file=log)
                    for i in range(len(self.timerevr_flipped[ii])):
                        print('    T='+str(self.timerevr_flipped[ii][i]), file=log)
                        print('    R={:2d},{:2d},{:2d}'.format(*self.rotation_flipped[ii][i][0,:]), file=log)
                        for j in range(1,3):
                            print('      {:2d},{:2d},{:2d}'.format(*self.rotation_flipped[ii][i][j,:]), file=log)
            else:
                self.flipped = cif.flipped
                self.rotation_flipped = cif.rotation_flipped
                self.timerevr_flipped = cif.timerevr_flipped
                
        
    def ReadStruct(self, fname, log=sys.stdout):
        f = open(fname, 'r')
        self.title = next(f).strip()
        line = next(f)
        self.lattyp = line[0:4].strip()
        self.nat = int(line[27:30])
        self.sgnum = int(line[30:33])
        self.sgname = line[33:].strip()
        self.mode = next(f)[13:17]
        line = next(f)
        self.a, self.b, self.c, self.alpha, self.beta, self.gamma = [float(line[i:i+10]) for i in range(0,60,10)]
        self.iatnr=[]
        self.isplit=[]
        self.mult=[]
        self.pos=[]
        self.aname=[]
        self.jrj=[]
        self.r0=[]
        self.rmt=[]
        self.Znuc=[]
        self.rotloc=[]
        for iat in range(self.nat):
            line = next(f)
            self.iatnr.append( int(line[4:8]) )
            pos = [float(line[col:col+10]) for col in 12+13*arange(3)]
            self.pos.append(pos)
            
            line = next(f)
            _mult_ = int(line[15:17])
            self.mult.append(_mult_)
            self.isplit.append( int(line[34:37]) )
        
            for mu in range(_mult_-1):
                line = next(f)
                pos = [float(line[col:col+10]) for col in 12+13*arange(3)]
                self.pos.append(pos)
                
            #self.pos.append(pos)
                
            line = next(f)
            self.aname.append( line[0:10].strip() )
            self.jrj.append( int(line[15:20]) )
            self.r0.append( float(line[25:35]) )
            self.rmt.append( float(line[40:50]) )
            self.Znuc.append( int(float(line[55:65])) )
            
            rt=[]
            for i in range(3):
                line = next(f)
                rt.append( [float(line[col:col+10]) for col in 20+10*arange(3)] )
            self.rotloc.append(array(rt).T)

        self.pos = array(self.pos)
        
        self.timat=[]
        self.tau=[]
        
        line = next(f)
        dt = line.split()
        if len(line)<2 or dt[1][:6]!='NUMBER':
            f.close()
            print('No symmetry operations yet!', file=log)
            return
        nsym = int(dt[0])
        self.timat  = zeros((nsym,3,3),dtype=int) # it is transpose compared to w2k fortran convention
        self.tau    = zeros((nsym,3))
        for isym in range(nsym):
            for j in range(3):
                line = next(f)
                self.timat[isym,j,:] = [int(line[0:2]),int(line[2:4]),int(line[4:6])]
                self.tau[isym,j]     = float(line[6:16])
            ii = int(next(f))
            if (ii != isym+1):
                print('WARNING : issues with reading symmetry operations in struct file at isym=', isym+1, 'ii=', ii, file=log)
                f.close()
                return
        f.close()
        
    def WriteStruct(self, fname, log=sys.stdout):
        if self.inAngs:
            self.a, self.b, self.c = self.a*self.tobohr, self.b*self.tobohr, self.c*self.tobohr
        
        if self.HRtransf: # Currently it is never done. Wondering when is this needed?
            i = 0
            for jatom in range(self.nat):
                for m in range(self.mult[jatom]):
                    self.pos[i,:] = (self.pos[i,:] @ self.hex2rho) % 1.0
                    i += 1
        #### Warning: W2k can handle only CXZ & gamma!=90 monoclinic case, hence
        # we convert CXY and CYZ monoclinic to CXZ. 
        if self.lattyp[:1]=='C' and (abs(self.alpha-90)>1e-4 or abs(self.beta-90)>1e-4 or abs(self.gamma-90)>1e-4):
            if self.lattyp == 'CXY':
                if abs(self.alpha-90)<1e-4 and abs(self.gamma-90)<1e-4: # beta != 90
                    print('WARN: Conversion of Y  <--> Z for monoclinic CXY lattice', file=log)
                    #self.convert_cxy2cxz()  # W2k does not allow CXY in trigonal, monoclinic
                    self.convert_by_exchange(1,2)
                    self.lattyp = 'CXZ'
                elif abs(self.beta-90)<1e-4 and abs(self.gamma-90)<1e-4: # alpha!=90
                    print('WARN: Conversion of (X,Y,Z) --> (Z,X,Y) for monoclinic CXY lattice', file=log)
                    self.convert_by_cyclic(-1)
                    self.lattyp = 'CXZ'
                else:
                    print('ERROR: Don\'t know how to convert CXY with gamma!=90 to CXZ with gamma!=90.', file=log)
            elif self.lattyp == 'CYZ':
                if abs(self.alpha-90)<1e-4 and abs(self.gamma-90)<1e-4: # beta!=90
                    print('WARN: Conversion of (X,Y,Z) --> (Y,Z,X) from CYZ monoclinic to CXZ lattice', file=log)
                    #self.convert_cyz2cxz()
                    self.convert_by_cyclic(1)
                    self.lattyp = 'CXZ'
                elif abs(self.alpha-90)<1e-4 and abs(self.beta-90)<1e-4: # gamma!=90
                    print('WARN: Conversion of (X,Y) --> (Y,X) from CYZ monoclinic to CXZ lattice', file=log)
                    #self.convert_cyz2cxz_2()
                    self.convert_by_exchange(0,1)
                    self.lattyp = 'CXZ'
                else:
                    print('ERROR: Don\'t know how to convert CYZ with alpha!=90 to CXZ with gamma!=90.', file=log)
            elif self.lattyp == 'CXZ':
                if abs(self.alpha-90)<1e-4 and abs(self.beta-90)<1e-4:
                    pass
                elif abs(self.beta-90)<1e-4 and abs(self.gamma-90)<1e-4:
                    print('WARN: Conversion of (X,Z) --> (Z,X) from CXZ alpha!=90 to gamma!=90 lattice', file=log)
                    self.convert_by_exchange(0,2)
                else:
                    print('ERROR: Don\'t know how to convert CXZ with beta!=90 to CXZ with gamma!=90.', file=log)
            
                
        with open(fname, 'w') as f:
            print(self.title, file=f)
            print('{:3s} LATTICE,NONEQUIV.ATOMS:{:3d} {:3d} {:8s}'.format(self.lattyp,self.nat,self.sgnum,self.sgname), file=f)
            print('MODE OF CALC=RELA unit=bohr', file=f)
            print('{:10.6f}{:10.6f}{:10.6f}{:10.6f}{:10.6f}{:10.6f}'.format(self.a,self.b,self.c,self.alpha,self.beta,self.gamma), file=f)
            index = 0
            for jatom in range(self.nat):
                print('ATOM{:4d}: X={:10.8f} Y={:10.8f} Z={:10.8f}'.format(self.iatnr[jatom],*self.pos[index,:]), file=f)
                print('          MULT={:2d}          ISPLIT={:2d}'.format(self.mult[jatom],self.isplit[jatom]), file=f)
                index += 1
                for m in range(1,self.mult[jatom]):
                    print('    {:4d}: X={:10.8f} Y={:10.8f} Z={:10.8f}'.format(self.iatnr[jatom],*self.pos[index,:]),file=f)
                    index += 1
                
                # fix for label not in position3, PB June 2016
                if len(self.aname[jatom])>=3 and self.aname[jatom][2]==' ':
                    rest_name = self.aname[jatom][3:].strip()
                    self.aname[jatom] = self.aname[jatom][:2]+rest_name
                if len(self.aname[jatom])>10:  # w2k does not allow more than 10 spaces name
                    self.aname[jatom] = self.aname[jatom][:10]
                R0s='{:10.9f}'.format(self.r0[jatom])
                print('{:10s} NPT={:5d}  R0={:10s} RMT={:10.5f}   Z:{:10.5f}'.format(self.aname[jatom], self.jrj[jatom], R0s[1:], self.rmt[jatom], self.Znuc[jatom]), file=f)
                print('LOCAL ROT MATRIX:   {:10.7f}{:10.7f}{:10.7f}'.format(*self.rotloc[jatom][:,0]), file=f)
                print('                    {:10.7f}{:10.7f}{:10.7f}'.format(*self.rotloc[jatom][:,1]), file=f)
                print('                    {:10.7f}{:10.7f}{:10.7f}'.format(*self.rotloc[jatom][:,2]), file=f)
            nsym = shape(self.tau)[0]
            #print('{:4d}      NUMBER OF SYMMETRY OPERATIONS'.format(nsym), file=f)
            print('{:4d}      NUMBER OF SYMMETRY OPERATIONS'.format(0), file=f)  # because w2k should determine the symmetry operations, which sometimes don't agree
            for iord in range(nsym):
                for j in range(3):
                    print('{:2d}{:2d}{:2d} {:10.8f}'.format(*self.timat[iord,j,:], self.tau[iord,j]), file=f)
                print('      {:2d}'.format(iord+1), file=f)

    def convert_by_cyclic(self, direction):
        """ This rouine cyclically exchanges xyz components in positive or negative direction.
            direction==1 : (x,y,z) --> (y,z,x)
            direction==-1: (x,y,z) --> (z,x,y)
            W2k can not handle monoclinic C lattice, except CXZ with gamma!=90.
            We thus need to convert every C lattice to CXZ with alpha=90, beta=90, gamma!=90.
        """
        if direction==1:
            self.sgname = self.sgname + ' con'
            self.b, self.c, self.a = self.a, self.b, self.c
            self.beta, self.gamma, self.alpha = self.alpha, self.beta, self.gamma
            index = 0
            for jatom in range(self.nat):
                for m in range(self.mult[jatom]):
                    self.pos[index,1], self.pos[index,2],self.pos[index,0] = self.pos[index,0], self.pos[index,1], self.pos[index,2]
                    index += 1
            
            for kname in self.kpath:
                self.kpath[kname][1],self.kpath[kname][2],self.kpath[kname][0] = self.kpath[kname][0],self.kpath[kname][1],self.kpath[kname][2]
                
            if self.neighbrs:
                for jatom in range(len(self.neighbrs)):
                    for i in range(len(self.neighbrs[jatom])):
                        R = self.neighbrs[jatom][i][2]
                        R[1],R[2],R[0]=R[0],R[1],R[2]

            # converts also symmetry operations
            Tp  = array([[0,1,0], [0,0,1], [1,0,0]], dtype=int)
            TpI = array([[0,0,1], [1,0,0], [0,1,0]], dtype=int)
            for isym in range(len(self.tau)):
                self.tau[isym,1],self.tau[isym,2],self.tau[isym,0] = self.tau[isym,0],self.tau[isym,1],self.tau[isym,2]
                self.timat[isym,:,:] = TpI @ self.timat[isym,:,:] @ Tp
            if self.flipped is not None:
                for ii in self.flipped:
                    for i in range(len(self.rotation_flipped[ii])):
                        self.rotation_flipped[ii][i][:,:] = TpI @ self.rotation_flipped[ii][i][:,:] @ Tp
                        
        elif direction==-1:
            self.sgname = self.sgname + ' con'
            self.c, self.a, self.b = self.a, self.b, self.c
            self.gamma, self.alpha, self.beta = self.alpha, self.beta, self.gamma
            index = 0
            for jatom in range(self.nat):
                for m in range(self.mult[jatom]):
                    self.pos[index,2], self.pos[index,0],self.pos[index,1] = self.pos[index,0], self.pos[index,1], self.pos[index,2]
                    index += 1
            
            for kname in self.kpath:
                self.kpath[kname][2],self.kpath[kname][0],self.kpath[kname][1] = self.kpath[kname][0],self.kpath[kname][1],self.kpath[kname][2]
                
            if self.neighbrs:
                for jatom in range(len(self.neighbrs)):
                    for i in range(len(self.neighbrs[jatom])):
                        R = self.neighbrs[jatom][i][2]
                        R[2],R[0],R[1]=R[0],R[1],R[2]
                    
            # converts also symmetry operations
            TpI = array([[0,1,0], [0,0,1], [1,0,0]], dtype=int)
            Tp  = array([[0,0,1], [1,0,0], [0,1,0]], dtype=int)
            for isym in range(len(self.tau)):
                self.tau[isym,2],self.tau[isym,0],self.tau[isym,1] = self.tau[isym,0],self.tau[isym,1],self.tau[isym,2]
                self.timat[isym,:,:] = TpI @ self.timat[isym,:,:] @ Tp
            if self.flipped is not None:
                for ii in self.flipped:
                    for i in range(len(self.rotation_flipped[ii])):
                        self.rotation_flipped[ii][i][:,:] = TpI @ self.rotation_flipped[ii][i][:,:] @ Tp
        else:
            print('ERROR: Not defined')
            
    def convert_by_exchange(self, i1, i2):
        """ This routine exchanges i1 and i2 components in all positions, and abc.
            W2k can not handle monoclinic C lattice, except CXZ with gamma!=90.
            We thus need to convert every C lattice to CXZ with alpha=90, beta=90, gamma!=90.
        """
        self.sgname = self.sgname + ' con'
        abc = [self.a, self.b, self.c]
        abg = [self.alpha, self.beta, self.gamma]
        abc[i1],abc[i2] = abc[i2], abc[i1]
        abg[i1],abg[i2] = abg[i2], abg[i1]
        self.a, self.b, self.c = abc
        self.alpha, self.beta, self.gamma = abg
        index = 0
        for jatom in range(self.nat):
            for m in range(self.mult[jatom]):
                self.pos[index,i1], self.pos[index,i2] = self.pos[index,i2], self.pos[index,i1]
                index += 1
        for kname in self.kpath:
            self.kpath[kname][i1], self.kpath[kname][i2] = self.kpath[kname][i2],self.kpath[kname][i1]
        
        if self.neighbrs:
            for jatom in range(len(self.neighbrs)):
                for i in range(len(self.neighbrs[jatom])):
                    R = self.neighbrs[jatom][i][2]
                    R[i1], R[i2] = R[i2], R[i1]
        
        # converts also symmetry operations
        Ts = array(identity(3),dtype=int)
        Ts[i1,:]=0
        Ts[i2,:]=0
        Ts[i1,i2]=Ts[i2,i1]=1
        for isym in range(len(self.tau)):
            self.tau[isym,i1],self.tau[isym,i2] = self.tau[isym,i2],self.tau[isym,i1]
            self.timat[isym,:,:] = Ts @ self.timat[isym,:,:] @ Ts

        if self.flipped is not None:
            for ii in self.flipped:
                for i in range(len(self.rotation_flipped[ii])):
                    self.rotation_flipped[ii][i][:,:] = Ts @ self.rotation_flipped[ii][i][:,:] @ Ts

def Z2r0(Z):
    "This is the choice for real space integration mesh (point closest to nucleous) in w2k."
    if (Z>71): return 5e-6
    if (36<Z<=71): return 1e-5
    if (18<Z<=36): return 5e-5
    if (Z<=18): return 1e-4
    return 5e-6

def Element_name(name):
    if len(name)>1 and name[1].islower():
        return name[:2]
    else:
        return name[:1]

def getlattype(sgname, sgnum, glide_plane):
    """ In W2k we need Bravais lattice type, which is either P,B,F,H,R,CXY=C,CXZ=B,CYZ=A
    It is uniquely given by space group number and glide plane of the form [1,1,0] or [1,0,1] 
    from symmetry operations.
    """
    if sgnum in [1,2,3,4,6,7,10,11,13,14,16,17,18,19,25,26,27,28,29,30,31,32,33,34,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,
                     75,76,77,78,81,83,84,85,86,89,90,91,92,93,94,95,96,99,100,101,102,103,104,105,106,111,112,113,114,115,116,117,
                     118,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,195,198,200,201,205,207,208,212,213,215,218,
                     221,222,223,224]:
        return 'P'
    if sgnum in [5,8,9,12,15,20,21,35,36,37,38,39,40,41,63,64,65,66,67,68]:
        if sum(glide_plane)==2:
            if array_equal(glide_plane,[1,1,0]):
                return 'CXY'
            elif array_equal(glide_plane,[1,0,1]):
                return 'CXZ'
            elif array_equal(glide_plane,[0,1,1]):
                return 'CYZ'
            else:
                print('ERROR wrong glide_plane=', glide_plane, 'and the algorithm to determine Bravais lattices failed!')
                return ''
        elif sum(glide_plane)==3: # I121
            return 'B' 
        else:
            print('ERROR glide_plane=', glide_plane, 'and the algorithm to determine Bravais lattices failed!')
            return ''
    if sgnum in [22,42,43,69,70,196,202,203,209,210,216,219,225,226,227,228]:
        return 'F'
    if sgnum in [23,24,44,45,46,71,72,73,74,79,80,82,87,88,97,98,107,108,109,110,119,
                     120,121,122,139,140,141,142,197,199,204,206,211,214,217,220,229,230]:
        return 'B'
    if sgnum in [143,144,145,147,149,150,151,152,153,154,156,157,158,159,162,163,164,165,168,169,170,171,172,
                     173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194]:
        return 'H'
    if sgnum in [146,148,155,160,161,166,167]:
        return 'R'
    print('ERROR the algorithm to determine Bravais lattices failed!')
    return ''

def RMT_ratios(Z1, Z2, reduction=0.0):
    reduc = (1 - reduction/100)
    _rkmt_fact_ = [60,60,45,50,50,55,60,65,70,60,65,65,65,50,55,60,65,60,65,65,
                   75,75,75,75,80,80,80,80,80,80,75,75,75,75,75,60,65,65,75,75,
                   75,75,60,80,80,80,80,80,80,80,80,80,80,60,65,65,80,80,85,85,
                   85,85,85,85,85,85,85,85,85,85,85,80,80,80,80,85,85,85,85,85,
                   85,85,85,85,85,60,85,85,85,85,85,85,85,85,85,85,85,85,85,85,
                   85,85,85,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60]
    fac12 = [0.497*reduc, 0.497*reduc]
    if (Z1<2 and  Z2>1):
        fac12 = [fac12[0]*0.7, fac12[1]*1.3]
    elif (Z2<2 and Z1>1):
        fac12 = [fac12[0]*1.3, fac12[1]*0.7]
    else:
        f1 = _rkmt_fact_[Z1-1]
        f2 = _rkmt_fact_[Z2-1]
        fac12 = [fac12[0]*(1.0 +(f1-f2)*0.005), fac12[1]*(1.0 -(f1-f2)*0.005)]
    return fac12

class TheSame:
    """ In W2k struct file the Bravais lattice contains only atoms 
    from primitive cell, but the coordinates are written in conventional cell.
    Here we check if an atom r1 and r2 (expressed in conventional cell) are equivalent 
    atoms for primitive cell.
    """
    def __init__(self, lattyp):
        self.tv = []
        if lattyp=='B':
            for i in [-1,1]:
                for j in [-1,1]:
                    for k in [-1,1]:
                        self.tv.append([0.5*i, 0.5*j, 0.5*k])
        elif lattyp[0] == 'F':
            for i in [-1,1]:
                for j in [-1,1]:
                    self.tv.append([0.5*i, 0.5*j, 0.0])
                    self.tv.append([0.5*i, 0.0, 0.5*j])
                    self.tv.append([0.0, 0.5*i, 0.5*j])
        elif lattyp == 'CXY':
            for i in [-1,1]:
                for j in [-1,1]:
                    self.tv.append([0.5*i, 0.5*j, 0.0])
        elif lattyp == 'CXZ':
            for i in [-1,1]:
                for j in [-1,1]:
                    self.tv.append([0.5*i, 0.0, 0.5*j])
        elif lattyp == 'CYZ':
            for i in [-1,1]:
                for j in [-1,1]:
                    self.tv.append([0.0, 0.5*i, 0.5*j])
        elif lattyp == 'R':
            for i in [-1,0,1]:
                for j in [-1,0,1]:
                    for k in [-1,0,1]:
                        self.tv.append([1/3.+i, 2/3.+j, 2/3.+k])
                        self.tv.append([2/3.+i, 1/3.+j, 1/3.+k])
        self.tv = array(self.tv)
    def __call__(self, r1, r2):
        small = 1e-5
        dx = sum( abs (self.tv + (r2-r1)), axis=1)
        return any(dx<small)
        
        # old-fashioned implementation
        Qsame = False
        for tv in self.tv:
            dx = (r1 - r2 - tv) #% 1.0
            if sum(abs(dx))<small:
                Qsame = True
                break
        return Qsame


def correct_element_names(anames):
    """In W2k structure file element is read as first two characters of aname.
    If the first two characters do not represent a name of element, W2k issues an error.
    For example Cr2+ is valid aname, however, O2- is not, and we need to use "O 2-" instead.
    """
    for i in range(len(anames)):
        name = anames[i]
        if len(name)>1 and not name[1].isalpha():
            anames[i] = name[:1]+' '+name[1:]
    return anames
    
def str2float(text):
    """Remove uncertainty brackets from strings and return the float."""
    try:
        # Note that the ending ) is sometimes missing. That is why the code has
        # been modified to treat it as optional. Same logic applies to lists.
        return float(re.sub(r"\(.+\)*", "", text))
    except TypeError:
        if isinstance(text, list) and len(text) == 1:
            return float(re.sub(r"\(.+\)*", "", text[0]))
    except ValueError as exc:
        if text.strip() == ".":
            return 0
        raise exc
    raise ValueError(f"{text} cannot be converted to float")

def is_number(text):
    return text.lstrip('-').replace('.','',1).replace('(','').replace(')','')

def get_matching_coord(coord, parser, coord_to_species):
    keys = list(coord_to_species)
    coords = array(keys)
    for op in parser.symmetry_operations:
        frac_coord = op.operate(coord)
        indices = find_in_coord_list_pbc(coords, frac_coord, atol=parser._site_tolerance)
        if len(indices) > 0:
            return keys[indices[0]]
    return False

    
class CifParser_W2k:
    def __init__(self, fname, log=sys.stdout, Qmagnetic=True, cmp_neighbors=False, ndecimals=3, nradius=4.7):
        self.log = log
        self.parser = CifParser(fname,occupancy_tolerance=3)
        cif_as_dict = self.parser.as_dict()       # get dictionary of cif information
        for k,data in cif_as_dict.items():   # we could have several structures (or items) in cif file
            if '_atom_site_fract_x' in data: # only if _atom_site_fract_x exists, it is a structure specification.
                self.cx = [str2float(a) for a in data['_atom_site_fract_x'] if is_number(a)]
                if len(self.cx)==0: continue       # it was not a valid _atom_site_fract_x, hence skip this data
                self.cy = [str2float(a) for a in data['_atom_site_fract_y'] if is_number(a)]
                self.cz = [str2float(a) for a in data['_atom_site_fract_z'] if is_number(a)]
                if '_atom_site_type_symbol' in data and len(data['_atom_site_type_symbol'])==len(self.cx):
                    self.cname = data['_atom_site_type_symbol']
                elif '_atom_site_label' in data and len(data['_atom_site_label'])==len(self.cx):
                    self.cname = data['_atom_site_label']
                else:
                    print('ERROR: Could not find labels for atoms in cif')
                self.ca = str2float(data['_cell_length_a'])
                self.cb = str2float(data['_cell_length_b'])
                self.cc = str2float(data['_cell_length_c'])
                
                self.calpha = str2float(data['_cell_angle_alpha'])
                self.cbeta = str2float(data['_cell_angle_beta'])
                self.cgamma = str2float(data['_cell_angle_gamma'])
                
                self.ssetting, self.sgname, self.sgnum = None, None, None
                if '_symmetry_cell_setting' in data:
                    self.ssetting = data['_symmetry_cell_setting']
                if '_symmetry_space_group_name_H-M' in data:
                    self.sgname = data['_symmetry_space_group_name_H-M']
                elif '_space_group_name_H-M_alt' in data:
                    self.sgname = data['_space_group_name_H-M_alt']
                elif '_parent_space_group.name_H-M' in data:
                    self.sgname = data['_parent_space_group.name_H-M']
                elif '_parent_space_group.name_H-M_alt' in data:
                    self.sgname = data['_parent_space_group.name_H-M_alt']
                if '_symmetry_Int_Tables_number' in data:
                    self.sgnum = data['_symmetry_Int_Tables_number']
                elif '_symmetry_space_group_IT_number' in data:
                    self.sgnum = data['_symmetry_space_group_IT_number']
                elif '_space_group_IT_number' in data:
                    self.sgnum = data['_space_group_IT_number']
                elif '_parent_space_group.IT_number' in data:
                    self.sgnum = data['_parent_space_group.IT_number']
                    
                print('Information directly extracted from cif file:', file=log)
                for k2, v2 in data.items():
                    print('key=', k+':'+k2,' value=', v2, file=log)
                print(file=log)
                print('structure directly from given input cif file', file=log)
                print('---------------------------------------------', file=log)
                print('ca=', self.ca, 'cb=', self.cb, 'cc=', self.cc, 'calpha=', self.calpha, 'cbeta=', self.cbeta, 'cgamma=', self.cgamma, file=log)
                print('sgname=', self.sgname, 'sgnum=', self.sgnum, 'ssetting=', self.ssetting, file=log)
                for i in range(len(self.cname)):
                    print('{:3d} {:8s} {:f}, {:f}, {:f}'.format(i+1,self.cname[i],self.cx[i], self.cy[i], self.cz[i]), file=log)
                print(file=log)
                
                break # We only allow one structure to be read. If multiple structures in cif, modify cif
                
        # Now we use pymatgen to create structure from cif information.
        # This is easier than figuring out group operations ourselves.
        # We could called "parser.parse_structures", but this can destroy "parser.symmetry_operations"
        # in case there is more than one structure in the cif file. We want to exit the loop as soon as we find
        # the first viable structure, as above. This will preserve symmetry operations
        #
        #self.structure = self.parser.parse_structures(primitive=False,symmetrized=True)[0]  # taking only the first viable structure, just like above with cif
        #
        for idx, dct in enumerate(self.parser._cif.data.values()):
            try:
                #print('dct=', dct)
                self.structure = self.parser._get_structure(dct,primitive=False,symmetrized=False, check_occu=True)
                if self.structure:
                    break
            except (KeyError, ValueError) as exc:
                msg = f"No structure parsed for section {idx + 1} in CIF.\n{exc}"
                #if on_error == "raise":
                #    raise ValueError(msg) from exc
                #if on_error == "warn":
                #    warnings.warn(msg)
                self.parser.warnings.append(msg)
        print('Pymatgen structure information for conventional unit cell:', file=log)
        print('---------------------------------------------', file=log)
        print(self.structure, file=log)
        print('---------------', file=log)
        print('lattice.abc=', self.structure.lattice.abc, file=log)
        print('lattice.angles=', self.structure.lattice.angles, file=log)
        print('lattice.matrix=', self.structure.lattice.matrix, file=log)
        print('composition.formula=', self.structure.composition.formula, file=log)
        print('composition.reduced_formula=', self.structure.composition.reduced_formula, file=log)
        #print('lattice.pbc(periodic boundary conditions)=', self.structure.lattice.pbc, file=log)
        print('len(self.structure)=', len(self.structure), file=log)
        print('structure[i].site.species_string, structure[i].site.frac_coords:', file=log)
        #for site in self.structure:
        #    print(site.species_string, site.frac_coords, file=log)

        Znuc=[]
        for site in self.structure:
            Znuc.append([st.Z for st in site.species][0])
                
        # We now have all sites in a conventional unit cell in self.structure
        # But w2k requires them to be grouped by the symmetry, i.e., all equivalent grouped into list
        # In addition, we will need to eliminate some sites from the list when we have lattice with type F, B or C. But this will be done later.
        
        # First we create distionary of sites, which have the same element. Unfortunately pymatgen looses the information about
        # sites in conventional structure and atoms in the cif file.
        # Sites with different element (site.species_string) are definitely inequivalent
        # We decided to only keep an index on self.structure.sites, because in this way we can access also
        # nearest neigbor and other properties in pymatgen.structure
        indx_by_element={}
        self.Z_element={}
        for ii,site in enumerate(self.structure):
            if site.species_string in indx_by_element: #coord_by_element:
                #coord_by_element[site.species_string].append(site.frac_coords)
                indx_by_element[site.species_string].append(ii)
            else:
                #coord_by_element[site.species_string] = [site.frac_coords]
                indx_by_element[site.species_string] = [ii]
            self.Z_element[site.species_string]=Znuc[ii]

        magnetic=False
        if hasattr(self.parser.symmetry_operations[0], 'time_reversal'):
            magnetic=True
        Qmagnetic = magnetic and Qmagnetic
        
        #print('indx_by_element=', indx_by_element)
        print('Number of symmetry operations available=', len(self.parser.symmetry_operations), file=log)
        print('magnetic cif=', magnetic, file=log)

        for isym,op in enumerate(self.parser.symmetry_operations):
            timat = array(np.round(op.affine_matrix[:3,:3]),dtype=int)
            tau = op.affine_matrix[:3,3]
            if np.allclose(timat,np.identity(3,dtype=int)):
                if sum(abs(tau))==0:
                    symmetry_operations_type='I'
                else:
                    symmetry_operations_type='T'
            elif np.allclose(timat,-np.identity(3,dtype=int)) and sum(abs(tau))==0:
                symmetry_operations_type='P'
            else:
                symmetry_operations_type='R'
            if magnetic:
                print(isym, ':', 'Time='+str(op.time_reversal), symmetry_operations_type, file=log)
            else:
                print(isym,':', symmetry_operations_type, file=log)
            for i in range(3):
                print('{:2d}{:2d}{:2d} {:10.8f}'.format(*timat[i,:], tau[i]), file=log)


        acoords = [site.frac_coords%1.0 for site in self.structure.sites] # all coordinates but always inside home unit cell
        if True:
            # This is for nonmagnetic system, where we don't need to introduce any new splitting
            pgroups=[] # Will contain groups of equivalent sites
            #print('indx_by_element=', indx_by_element, file=log)
            for spec in indx_by_element:  # Next we loop through element in this structure
                indices = indx_by_element[spec] # indices on sites of the same element
                grp = [[indices[0]]]  # we have just one group with the first entry in coords. All entries in coords will have index in grp
                site = self.structure.sites[indices[0]]
                #mfinite=False
                for ii in indices[1:]: # loop over possibly equivalent sites (of the same element)
                    coord = acoords[ii]    # fractional coordinate of that site self.structure.sites[ii]
                    site = self.structure.sites[ii]
                    Qequivalent=False  # for now we thing this might not be equivalent to any other site in grp
                    for ig,op in enumerate(self.parser.symmetry_operations): # loop over all symmetry operations
                        ncoord = op.operate(coord) % 1.0       # apply symmetry operation to this site, and produce ncoord
                        for ty in grp:                         # if coord is equivalent to any existing group in grp, say grp[i]==ty, then one of ncoord should be equal to first entry in grp[i]
                            coord0 = acoords[ty[0]]            # the first entry in the group grp[i][0]
                            if sum(abs(ncoord-coord0))<1e-5:   # is it equal to current coord when a symmetry operation is applied?
                                print('Group['+str(ig)+'] applied to r['+str(ii)+']=',coord,
                                          'gives r['+str(ty[0])+']=', coord0, file=log)
                                Qequivalent=True               # we found which group it corresponds to
                                ty.append(ii)                  # grp[i]=ty is extended with coord[i]
                                break                          # loop over groups and also symmetry operations can be finished
                        if Qequivalent:
                            break
                    if not Qequivalent:                        # if this site (coord) is not equivalent to any existing group in grp, than we create a new group with a single entry [i]
                        grp.append([ii])
                        #print('Starting anew with ii='+str(ii)+' grp=', grp, file=log)
                    #print(spec+' grp=', grp)
                pgroups.append(grp)  # all groups for this element
            print('pgroups=', pgroups, file=log)
            self.pgroups=pgroups
            
                
        amagmoms = [0 for site in self.structure.sites]
        min_Molap = 0.7  # We will treat magnetic moments (anti)parallel if their overlap is more than that
        if Qmagnetic:
            for i,site in enumerate(self.structure.sites):
                if 'magmom' in site.properties:
                    amagmoms[i] = site.properties["magmom"].moment
            #print('amagmoms=', amagmoms, file=log)
            operations_to_remove=set()
            self.flipped={}
            self.rotation_flipped={}
            self.timerevr_flipped={}
            mgroups=[]
            for spec in indx_by_element:  # Next we loop through element in this structure
                indices = indx_by_element[spec] # indices on sites of the same element
                site = self.structure.sites[indices[0]]
                Mi = site.properties["magmom"].moment
                amoments = [sum(abs(amagmoms[i])) for i in indices]
                if sum(amoments)<1e-5:
                    continue
                print(site.species_string+'['+str(indices[0])+']', site.frac_coords, 'M['+str(indices[0])+']=', Mi, file=log)
                grp = [[indices[0]]]  # we have just one group with the first entry in coords. All entries in coords will have index in grp
                grp_compare=[0]
                for ii in indices[1:]: # loop over possibly equivalent sites (of the same element)
                    coord = acoords[ii]    # fractional coordinate of that site self.structure.sites[ii]
                    site = self.structure.sites[ii]
                    Qequivalent=False  # for now we thing this might not be equivalent to any other site in grp
                    Mj = site.properties["magmom"].moment
                    mj = linalg.norm(Mj)
                    print(site.species_string+'['+str(ii)+']', coord, 'M['+str(ii)+']=', site.properties["magmom"].moment, file=log)
                    for ig,op in enumerate(self.parser.symmetry_operations): # loop over all symmetry operations
                        ncoord = op.operate(coord) % 1.0       # apply symmetry operation to this site, and produce ncoord
                        Mjp = op.operate_magmom(Mj).moment
                        for itg,ty in enumerate(grp):                         # if coord is equivalent to any existing group in grp, say grp[i]==ty, then one of ncoord should be equal to first entry in grp[i]
                            for i_grp_compare in range(len(ty)): #@@
                                #@@ i_grp_compare=grp_compare[itg]
                                #print('grp_compare=', grp_compare, 'igrp_compare=', i_grp_compare, 'ic=', ty[i_grp_compare], file=log)
                                coord0 = acoords[ty[i_grp_compare]]            # the first entry in the group grp[i][0]
                                M0 = amagmoms[ty[i_grp_compare]]
                                if sum(abs(ncoord-coord0))<1e-5:   # is it equal to current coord when a symmetry operation is applied?
                                    #print('Group['+str(ig)+'] applied to r['+str(ii)+']=',coord,
                                    #              ' gives r['+str(ty[i_grp_compare])+']=',coord0, file=log)
                                    #print('Group['+str(ig)+'] applied to r['+str(ii)+']=',coord,'and M['+str(ii)+']=[{:.3f},{:.3f},{:.3f}]'.format(*Mj),
                                    #              ' gives r['+str(ty[i_grp_compare])+']=',coord0,
                                    #              ' and M['+str(ty[i_grp_compare])+']=[{:.3f},{:.3f},{:.3f}]'.format(*Mjp),
                                    #              'M0=[{:.3f},{:.3f},{:.3f}]'.format(*M0), file=log)
                                    if sum(abs(M0-Mjp))<1e-5:
                                        print('Group['+str(ig)+'] applied to r['+str(ii)+']=',coord,'and M['+str(ii)+']=', Mj.tolist(),
                                                  ' gives r['+str(ty[i_grp_compare])+']=',coord0,
                                                  ' and M['+str(ty[i_grp_compare])+']=', M0.tolist(), file=log)
                                        m0 = linalg.norm(M0)
                                        Molap = dot(Mj,M0)/(mj*m0) if (mj!=0 and m0!=0) else 0
                                        if (mj==0 and m0==0): Molap = 1.0
                                        print('Molap=', Molap, file=log)
                                        moments_equal = sum(abs(Mj-M0))<1e-5 or (Molap>min_Molap and abs(m0-mj)<1e-5)
                                        if not moments_equal: #sum(abs(Mj-M0))>1e-5 and Molap<min_Molap:
                                            # If we treat up/down together, we should skip that
                                            # But if up is treated separately from down, and spin is flipped, we have to remove such symmetry operation
                                            # and declare that the two atoms are different
                                            moments_opposite = sum(abs(Mj+M0))<1e-5 or (Molap<-min_Molap and abs(m0-mj)<1e-5)
                                            if moments_opposite: #sum(abs(Mj+M0))<1e-5 or Molap < -min_Molap:
                                                self.flipped[ii]=ty[i_grp_compare]
                                                if ii not in self.rotation_flipped:
                                                    self.rotation_flipped[ii]=[ array(np.round(op.affine_matrix[:3,:3]),dtype=int) ]
                                                    self.timerevr_flipped[ii]=[ op.time_reversal ]
                                                else:
                                                    self.rotation_flipped[ii].append( array(np.round(op.affine_matrix[:3,:3]),dtype=int) )
                                                    self.timerevr_flipped[ii].append( op.time_reversal )
                                                if Molap > -0.99999:
                                                    print('WARNING Moments are only approximately antiparallel Molap=', Molap,
                                                              'and we are treating them as antiparallel', file=log)
                                                #if i_grp_compare<len(ty)-1:
                                                #    i_grp_compare += 1
                                                #    print('i_grp_compare increased to', i_grp_compare, 'because ty=', ty, file=log)
                                            #else:
                                            #    non_flipped[ii]=ty[i_grp_compare]
                                            #    if ii not in rotation_non_flipped:
                                            #        rotation_non_flipped[ii] = array(np.round(op.affine_matrix[:3,:3]),dtype=int)
                                            #        timerevr_non_flipped[ii] = op.time_reversal
                                            operations_to_remove.add(ig)
                                            print('Since moments are different ig=', ig, 'should be removed flipped=', self.flipped, file=log)
                                        else:
                                            # Only if spin direction is the same we can say that these are truly
                                            # equivalent
                                            Qequivalent=True               # we found which group it corresponds to
                                            if ii not in ty:
                                                ty.append(ii)                  # grp[i]=ty is extended with coord[i]
                                            if (Molap<0.99999):
                                                print('WARNING: Moments are only approximately equal Molap='+str(Molap),
                                                          'and we are treating them as equal', file=log)
                                                print('Since moments are approx. equal ig=', ig, ' is kept and '+str(ii)+' is in group with ',ty, file=log)
                                            else:
                                                print('Since moments are equal ig=', ig, ' is kept and '+str(ii)+' is in group with ',ty, file=log)
                                            #print('After adding ii='+str(ii)+' grp=', grp, file=log)
                                            #@@break                          # loop over groups and also symmetry operations can be finished
                    if not Qequivalent:                        # if this site (coord) is not equivalent to any existing group in grp, than we create a new group with a single entry [i]
                        grp.append([ii])
                        print('Starting anew with ii='+str(ii)+' grp=', grp, file=log)
                        grp_compare.append(0)
                        
                    if ii in self.flipped:
                        jj = self.flipped[ii]
                        for itg,ty in enumerate(grp):
                            if jj in ty:
                                if len(ty)>grp_compare[itg]+1:
                                    grp_compare[itg]+=1
                mgroups.append(grp)  # all groups for this element
            print('groups for magnetic atoms=', mgroups, file=log)
            
            # Sometimes the symmetry operations for time-reversal symmetry are not given, and moments are antiparallel, but above
            # algorithm does not recognize that, because it is not given as one of the symmetry operations.
            # We here check is we have such a case
            if len(self.flipped)==0:
                for grp in mgroups:
                    used=[]
                    for i in range(len(grp)):
                        grp1 = grp[i]
                        for j in range(len(grp)):
                            grp2 = grp[j]
                            #if i==j or len(grp2)!=len(grp1): continue
                            if i==j: continue
                            for i1 in grp1:
                                site1 = self.structure.sites[i1]
                                Mi1 = site1.properties["magmom"].moment
                                m1 = linalg.norm(Mi1)
                                for i2 in grp2:
                                    site2 = self.structure.sites[i2]
                                    Mi2 = site2.properties["magmom"].moment
                                    m2 = linalg.norm(Mi2)
                                    Molap = dot(Mi1,Mi2)/(m1*m2)
                                    if abs(m1-m2)>1e-5:
                                        continue
                                    print('i1={:2d} i2={:2d} m1={:8.4f} m2={:8.4f}'.format(i1,i2,m1,m2),
                                              'M1=', Mi1, 'M2=', Mi2, 'Molap={:5.3f}'.format(Molap), 'used=', used, file=log)
                                    if Molap<-min_Molap:
                                        if i1 not in self.flipped and i2 not in used:
                                            self.flipped[i1]=i2
                                            used.append(i2)
                                            self.timerevr_flipped[i1]=[]
                                            self.rotation_flipped[i1]=[]
                                        
            print('indx_by_element', indx_by_element, file=log)
            ##-----------------------------------------------------------------------------------------------------------------------------
            # Next we just make atoms of the same element equivalent, and later check the environment around
            # each atom to determine splitting. This is the same as w2k nn routine, and it seems to be most efficient
            # way to determine grouping of atoms. Namely, symmetry operations are sometimes not correctly given, 
            # or we get umbiguous results with symmetry operations.
            for spec in indx_by_element:  
                indices = indx_by_element[spec] # indices on sites of the same element
                #site = self.structure.sites[indices[0]]
                #Mi = site.properties["magmom"].moment
                amoments = [sum(abs(amagmoms[i])) for i in indices]
                if sum(amoments)<1e-5:
                    mgroups.append([indices])
            mgroups = sorted(mgroups, key=lambda a: a[0][0])
            print('starting with groups for magnetic atoms=', mgroups, file=log)
            
            print('atoms which are flipped of each other:', file=log)
            for ii in self.flipped:
                print(str(ii)+'->'+str(self.flipped[ii]), file=log)
                for i in range(len(self.timerevr_flipped[ii])):
                    print('    T='+str(self.timerevr_flipped[ii][i]), file=log)
                    print('    R={:2d},{:2d},{:2d}'.format(*self.rotation_flipped[ii][i][0,:]), file=log)
                    for j in range(1,3):
                        print('      {:2d},{:2d},{:2d}'.format(*self.rotation_flipped[ii][i][j,:]), file=log)


            self.original_symmetry_operations = [deepcopy(self.parser.symmetry_operations[iop])
                                                     for iop in range(len(self.parser.symmetry_operations))]
                        
            print('symmetry operations that need to be removed=', operations_to_remove, file=log)
            operations_to_remove = sorted(list(operations_to_remove))
            for iop in operations_to_remove[::-1]:
                del self.parser.symmetry_operations[iop]
            print('Symmetry operations left:', file=log)
            for isym,op in enumerate(self.parser.symmetry_operations):
                timat = array(np.round(op.affine_matrix[:3,:3]),dtype=int)
                tau = op.affine_matrix[:3,3]
                if np.allclose(timat,np.identity(3,dtype=int)):
                    if sum(abs(tau))==0:
                        symmetry_operations_type='I'
                    else:
                        symmetry_operations_type='T'
                elif np.allclose(timat,-np.identity(3,dtype=int)) and sum(abs(tau))==0:
                    symmetry_operations_type='P'
                else:
                    symmetry_operations_type='R'
                if magnetic:
                    print(isym, ':', 'Time='+str(op.time_reversal), symmetry_operations_type, file=log)
                else:
                    print(isym,':', symmetry_operations_type, file=log)
                for i in range(3):
                    print('{:2d}{:2d}{:2d} {:10.8f}'.format(*timat[i,:], tau[i]), file=log)
            
            ##-----------------------------------------------------------------------------------------------------------------------------
            # When we reduce the symmetry due to magnetism, some atoms might be missplaced.
            # Here we are using algorithm from w2k nn program, which looks at all neighbors of atoms, and
            # if several shells of atoms have equivalent neigbors, they are equivalent. Otherwise they are not.
            print('all_neigbrs:', file=log)
            all_sites = [self.structure.sites[i] for i in range(len(self.structure.sites))]
            center_indices, points_indices, offsets, distances = self.structure.get_neighbor_list(r=nradius, sites=all_sites, numerical_tol=1e-4)
            all_neighbrs=[]
            for i,isite in enumerate(all_sites):
                coord = isite.frac_coords @ self.structure.lattice.matrix
                #print('ATOM: {:3d} {:3s} AT {:9.5f} {:9.5f} {:9.5f}={:9.5f} {:9.5f} {:9.5f}'.
                #          format(i+1,self.structure.sites[i].species_string,*isite.coords,*isite.frac_coords),file=log)
                mask = center_indices==i
                points_indices_ = points_indices[mask]
                offsets_ = offsets[mask]
                distances_ = distances[mask]
                names_ = [self.structure.sites[p_i].species_string for p_i in points_indices_]
                cmp_to_sort = Cmp2Sort(distances_,names_,self.structure.sites[i].species_string)
                idx=sorted(range(len(distances_)), key=cmp_to_key(cmp_to_sort))
                _neighbrs_={}
                for j in idx:
                    if distances_[j]>1e-5:
                        nn = points_indices_[j]
                        dd = distances_[j]
                        r = self.structure.sites[nn]
                        Rj_fract = r.frac_coords + offsets_[j]
                        #print('  ATOM:{:3d} {:5s} AT {:9.5f} {:9.5f} {:9.5f} IS AWAY {:13.6f} ANG'.format(nn,r.species_string,*Rj_fract, distances_[j]), file=log)
                        sdst = str(np.round(distances_[j],decimals=ndecimals))
                        if sdst in _neighbrs_:
                            _neighbrs_[sdst].append(nn)
                        else:
                            _neighbrs_[sdst] = [nn]
                all_neighbrs.append(_neighbrs_)
            print('Going over all neighbors to determine equivalency:', file=log)
            #print('mgroups=', mgroups, file=log)
            
            element=[]
            for i in range(len(mgroups)):
                els = [el for el in indx_by_element if mgroups[i][0][0] in indx_by_element[el]]
                element.append(els[0])
            #element = [el for el in indx_by_element]
            #print('element=', element)
            igroups=[[] for i in range(len(self.structure))]
            ele2grp={}
            for i in range(len(mgroups)):
                for j in range(len(mgroups[i])):
                    ele = element[i]+str(j+1)
                    ele2grp[ele] = (i,j)
                    for k in range(len(mgroups[i][j])):
                        igroups[mgroups[i][j][k]] = ele
            
            print('igroups=', igroups, file=log)
            print('ele2grp=', ele2grp, file=log)
            print('mgroups=', mgroups, file=log)
            
            split_finish=True
            for itt in range(10):
                print('Iteration=', itt+1, file=log)
                fingerprint=[]
                for ii in range(len(all_sites)):
                    fingerp=[]
                    for ds in all_neighbrs[ii]:
                        fgprnt = ds+' '+' '.join([str(g) for g in sorted([igroups[j] for j in all_neighbrs[ii][ds]])])
                        fingerp.append(fgprnt)
                    print('atoms[{:2d}]:'.format(ii), fingerp, file=log)
                    fingerprint.append(fingerp)
                grp_fingerprint={}
                for i in range(len(mgroups)):
                    for j in range(len(mgroups[i])):
                        ele = element[i]+str(j+1)
                        ii = mgroups[i][j][0]
                        grp_fingerprint[ele] = fingerprint[ii]
                #print('grp_fingerprint=', grp_fingerprint)
                print('neigbors of inequivalent atom types after iteration {}'.format(itt+1), file=log)
                for ele in grp_fingerprint:
                    print('  ', ele, grp_fingerprint[ele], file=log)
                    
                split_finish=True
                for ii in range(len(all_sites)):
                    my_ele=igroups[ii]
                    expect = grp_fingerprint[my_ele]
                    have = fingerprint[ii]
                    if have != expect:
                        split_finish=False
                        print('atom[{:2d}] WARNING have != expect'.format(ii), file=log)
                        print('atom[{:2d}] with name={:s} neigbr='.format(ii,my_ele), have, file=log)
                        print('atom[{:2d}] with name={:s} expect='.format(ii,my_ele), expect, file=log)
                        (i_m,j_m) = ele2grp[my_ele]
                        FoundGroup=False
                        for ele,fgp in grp_fingerprint.items():
                            if have==fgp:
                                FoundGroup=True
                                (i,j) = ele2grp[ele]
                                mgroups[i][j].append(ii)
                                mgroups[i_m][j_m].remove(ii)
                                igroups[ii] = ele
                                print('Found', my_ele, 'goes into', ele, 'mgroups=', mgroups, file=log)
                                break
                        if not FoundGroup:
                            j_new = len(mgroups[i_m])
                            mgroups[i_m].append([ii])
                            mgroups[i_m][j_m].remove(ii)
                            my_ele_new = element[i_m]+str(j_new+1)
                            grp_fingerprint[my_ele_new] = have
                            ele2grp[my_ele_new] = (i_m, j_new)
                            igroups[ii] = my_ele_new
                            print('Not found', my_ele, 'changed to', my_ele_new, 'mgroups=', mgroups, file=log)
                if not split_finish:
                    print('igroups=', igroups, file=log)
                    print('ele2grp=', ele2grp, file=log)
                    print('grp_fingerprint=', grp_fingerprint, file=log)
                
                if split_finish:
                    break
            
            print('mgroups=', mgroups, file=log)
            groups = mgroups
            ##-----------------------------------------------------------------------------------------------------------------------------
        else:
            groups = self.pgroups
            if False:
                groups=[] # Will contain groups of equivalent sites
                for spec in indx_by_element:  # Next we loop through element in this structure
                    indices = indx_by_element[spec] # indices on sites of the same element
                    grp = [[indices[0]]]  # we have just one group with the first entry in coords. All entries in coords will have index in grp
                    site = self.structure.sites[indices[0]]
                    if Qmagnetic:
                        Mi = site.properties["magmom"].moment
                        if sum(abs(Mi))>0:
                            mfinite=True
                            print(site.species_string, site.frac_coords, 'Mi=', site.properties["magmom"].moment, file=log)
                        else:
                            mfinite=False
                    for ii in indices[1:]: # loop over possibly equivalent sites (of the same element)
                        coord = acoords[ii]    # fractional coordinate of that site self.structure.sites[ii]
                        site = self.structure.sites[ii]
                        Qequivalent=False  # for now we thing this might not be equivalent to any other site in grp
                        if Qmagnetic and mfinite:
                            Mj = site.properties["magmom"].moment
                            print(site.species_string, coord, 'Mj=', site.properties["magmom"].moment, file=log)
                        for ig,op in enumerate(self.parser.symmetry_operations): # loop over all symmetry operations
                            ncoord = op.operate(coord) % 1.0       # apply symmetry operation to this site, and produce ncoord
                            #if Qmagnetic and mfinite:
                            #    Mjp = op.operate_magmom(Mj).moment
                            for ty in grp:                         # if coord is equivalent to any existing group in grp, say grp[i]==ty, then one of ncoord should be equal to first entry in grp[i]
                                coord0 = acoords[ty[0]]            # the first entry in the group grp[i][0]
                                M0 = amagmoms[ty[0]]
                                #print('coord0=', ty[0], coord0)
                                if sum(abs(ncoord-coord0))<1e-5:   # is it equal to current coord when a symmetry operation is applied?
                                    if Qmagnetic and mfinite:
                                        Molap = dot(Mj,M0)/(linalg.norm(Mj)*linalg.norm(M0))
                                        if Molap > min_Molap: #sum(abs(Mj-M0))<=1e-5:
                                            # Only if spin direction is the same we can say that these are truly
                                            #print('coordinates equal M0=',M0, 'Mj=', Mj, 'ig=', ig, file=log)????
                                            print('Group['+str(ig)+'] applied to r=',coord,'and M=', Mj,' gives r=', coord0, ' and M0=', M0, file=log)
                                            Qequivalent=True               # we found which group it corresponds to
                                            ty.append(ii)                  # grp[i]=ty is extended with coord[i]
                                            #print('Since moments are equal ig=', ig, 'survives magnetic calc', file=log)
                                            break                          # loop over groups and also symmetry operations can be finished
                                    else:# we don't have moments, hence it is equivalent
                                        #print(' ncoord=', ncoord, 'is equivalent to typ', ty) # yes, coord is in group ty
                                        print('Group['+str(ig)+'] applied to r['+str(ii)+']=',coord,
                                                  'gives r['+str(ty[0])+']=', coord0, file=log)
                                        Qequivalent=True               # we found which group it corresponds to
                                        ty.append(ii)                  # grp[i]=ty is extended with coord[i]
                                        break                          # loop over groups and also symmetry operations can be finished
                            if Qequivalent:
                                break
                        if not Qequivalent:                        # if this site (coord) is not equivalent to any existing group in grp, than we create a new group with a single entry [i]
                            grp.append([ii])
                            #print('Starting anew with ii='+str(ii)+' grp=', grp, file=log)
                        #print(spec+' grp=', grp)
                    groups.append(grp)  # all groups for this element
            print('groups=', groups, file=log)
            
        
        # Since groups is only index array to coord_by_element, we will rather create a dictionary self.w2k_coords={},
        # which is easier to use. The keys will be name of the site, which might need to be modified when
        # multiple inequivalent sites have the same name.
        # The values in the dictionary are all equivalent sites,i.e., their fractional coordinates
        # For magnetic calculation, we add additional splitting of atoms in groups (due to orientation of momemnts) and hence we lost connection to cif file.
        # It would need more elaborate sorting, which is not done yet.
        ResortToCif= not Qmagnetic
        if ResortToCif:
            # It is quite inconvenient that pymatgen structure does not sort atoms in the same order as they are sorted in cif file.
            # We here resort all atoms so that they appear in the same order as in cif file.
            # First we create coordinates of atoms in the cif file in the order in which they appear in the cif file.
            cif_coords=zeros((len(self.cx),3))
            for i in range(len(self.cx)):
                cif_coords[i,:] = self.cx[i],self.cy[i],self.cz[i]
            # Next we go over all groups of atoms, which are equivalent, and resort them so that they appear in equal order to cif file
            groups_resort=[[] for i in range(len(cif_coords))] # groups in the right order
            cname=[[] for i in range(len(cif_coords))]  # name for the atom in final struct file
            for idx,spec in enumerate(indx_by_element):

                grp = groups[idx] # group of atoms of the same element. But not necessary equivalent. For example, all oxygen atoms in the structure.
                #print('idx=', idx, 'spec=', spec, 'grp=', grp)
                
                # Instead of comparing with all coordinates from the cif file, we only compare with those
                # that start with the same letter, i.e., likely correspond to the same element.
                which_icif = [i for i in range(len(self.cname)) if self.cname[i][0]==spec[0]] # those are probably for the same element

                #print('which_icif=', which_icif)
                
                for ig in range(len(grp)): # over all groups of oxygen atoms
                    gr = grp[ig]           # here we have a single group of oxygen atoms, which are all equivalent
                    # all coordinates of equivalent atoms and combined together into array
                    coords = array([self.structure.sites[j].frac_coords for j in gr])
                    # Now checking in which order they appear in cif file.
                    icif=-1
                    imin_cif=0
                    for i in which_icif:
                        # This is 2D array in which we compare all coords of this group with one entry in cif file
                        adiff = sum(abs((coords-cif_coords[i]) % 1.0), axis=1) 
                        imin = argmin(adiff)  # which has the smallest distance?
                        if adiff[imin]<1e-10: # self.structure.sites[grp[ig][imin]].frac_coord == cif_coords[i]
                            icif=i
                            imin_cif = imin
                            break
                    #print('imin_cif=', imin_cif)
                    if imin_cif>0: # We want the first atom in struct file to correspond exactly to the atom in cif file. Hence swapping them.
                        #print('imin_cif=', imin_cif, 'hence swaping two entries')
                        gr[0],gr[imin_cif] = gr[imin_cif],gr[0]
                        
                    if icif>=0: # found the equivalent atom from cif
                        groups_resort[icif] = gr
                    else:
                        ResortToCif=False # Can not resort since I failed to find equivalence
                        break
                    # Now creat a unique name for this set of equivalent atoms.
                    _spec_ = spec  # if only one type of atoms for this element
                    if len(grp)>1: # more than one group of atoms, i.e., inequivalent elements
                        if spec[-1] in ['+','-']:        # want to add in index to the name
                            _spec_ = spec[:-1]+str(ig+1)+spec[-1:]  # if using oxidation state, we add integer before oxidation state
                        else:
                            _spec_ = spec+str(ig+1)      # just append integer to previous name
                        self.Z_element[_spec_] = self.Z_element[spec]
                    cname[icif] = _spec_
                    #print(spec, 'grp=', gr, 'icif=', icif)

            # This is needed because cif file can have site-occupancies=0.5, in which case structure has smaller number of sites that cif-file
            #groups_resort = [x for x in groups_resort if x]
            #print('groups_resort=', groups_resort)
            
            first_sites_index=[]
            self.w2k_coords={}
            for ii,group in enumerate(groups_resort):
                if len(group)==0: continue
                psite = self.structure.sites[group[0]]
                self.w2k_coords[cname[ii]] = [self.structure.sites[j].frac_coords for j in group]
                first_sites_index.append( group[0] )
                #print('ii=', ii, 'group=', group)
                #all_sites_order.extend( group )
                #print(cname[ii], coords, psite.label, self.cname[ii], psite.species, psite.species_string)
                #self.Z_element[cname[ii]] = psite.specie.Z
            #print('w2k_coord=', self.w2k_coords)
            #print('all_sites_order=', all_sites_order)
            
        if not ResortToCif:
            self.w2k_coords={}
            for idx,spec in enumerate(indx_by_element): # over all coordinates from pymatgen.structure
                grp = groups[idx]                        # equivalent groups that we just determine above
                if len(grp)==1:                          # if all sites corresponding to this element are equivalent, we have a single group
                    # Just remember all equivalent coordinates, because we need them in strucure file
                    self.w2k_coords[spec] = [self.structure.sites[j].frac_coords for j in grp[0]]       # just store all sites into this group
                else:
                    for ig in range(len(grp)):           # over all inequivalent groups
                        # Let's create a unique name for each group
                        if spec[-1] in ['+','-']:        # want to add in index to the name
                            _spec_ = spec[:-1]+str(ig+1)+spec[-1:]  # if using oxidation state, we add integer before oxidation state
                        else:
                            _spec_ = spec+str(ig+1)      # just append integer to previous name
                        self.w2k_coords[_spec_] = [self.structure.sites[j].frac_coords for j in grp[ig]]  # and save the coordinates that were found equivalent above
                        self.Z_element[_spec_] = self.Z_element[spec]

            first_sites_index=[]
            all_sites_order=[]
            for grp in groups:
                for gr in grp:
                    jatom = len(first_sites_index)
                    first_sites_index.append(gr[0])
                    all_sites_order.extend(gr)
            all_sites_index=zeros(len(all_sites_order),dtype=int)
            for i in range(len(all_sites_order)):
                all_sites_index[all_sites_order[i]]=i
            print('all_sites_order=', all_sites_order,file=log)
            #print('all_sites_index=', all_sites_index.tolist(), file=log)

        if Qmagnetic:
            # We need to resort the index used in flipped to current order in self.w2k_coords
            _flipped_={}
            _timerevr_flipped_={}
            _rotation_flipped_={}
            for ii in self.flipped:
                jj = all_sites_index[ii]
                _flipped_[jj] = all_sites_index[self.flipped[ii]]
                if self.timerevr_flipped[ii]:
                    _timerevr_flipped_[jj] = self.timerevr_flipped[ii]
                    _rotation_flipped_[jj] = self.rotation_flipped[ii]
                else:
                    _timerevr_flipped_[jj] = []
                    _rotation_flipped_[jj] = []
            self.flipped = _flipped_
            self.timerevr_flipped = _timerevr_flipped_
            self.rotation_flipped = _rotation_flipped_
            print('now resorted atoms which are flipped of each other:', file=log)
            for ii in self.flipped:
                print(str(ii)+'->'+str(self.flipped[ii]), file=log)
                for i in range(len(self.timerevr_flipped[ii])):
                    print('    T='+str(self.timerevr_flipped[ii][i]), file=log)
                    print('    R={:2d},{:2d},{:2d}'.format(*self.rotation_flipped[ii][i][0,:]), file=log)
                    for j in range(1,3):
                        print('      {:2d},{:2d},{:2d}'.format(*self.rotation_flipped[ii][i][j,:]), file=log)


        Add_oxi=False
        self.neighbrs=[]
        if cmp_neighbors:
            first_sites = [self.structure.sites[i] for i in first_sites_index]
            if Add_oxi: self.structure.add_oxidation_state_by_guess()
            self.oxi_state={}
            #for site in self.structure:
            #    for specie, amt in site.species.items():
            #        print(specie, amt, (getattr(specie, "oxi_state", 0) or 0) )
            #print('first_sites_index=', first_sites_index)
            #print('first_sites=', first_sites)
            center_indices, points_indices, offsets, distances = self.structure.get_neighbor_list(r=4.0, sites=first_sites, numerical_tol=1e-4)
            self.Rmt = zeros(len(first_sites_index))
            print(center_indices, file=log)
            print(points_indices, file=log)
            for jatom,i in enumerate(first_sites_index):
                t=self.structure.sites[i]
                coord = t.frac_coords @ self.structure.lattice.matrix
                #print('coords=', coord)
                if Add_oxi:
                    try:
                        oxi_st = int(getattr(t.specie, "oxi_state", 0) or 0)
                    except:
                        oxi_st=0
                    self.oxi_state[jatom]=oxi_st
                else:
                    self.oxi_state[jatom]=0
                print('ATOM: {:3d} {:3s} AT {:9.5f} {:9.5f} {:9.5f}={:9.5f} {:9.5f} {:9.5f}'.format(jatom+1,self.structure.sites[i].species_string,*t.coords,*t.frac_coords),file=log)
                mask = center_indices==jatom
                points_indices_ = points_indices[mask]
                offsets_ = offsets[mask]
                distances_ = distances[mask]
                names_ = [self.structure.sites[p_i].species_string for p_i in points_indices_]
                cmp_to_sort = Cmp2Sort(distances_,names_,self.structure.sites[i].species_string)
                idx=sorted(range(len(distances_)), key=cmp_to_key(cmp_to_sort))# lambda i:(distances_[i],names_[i]))
                #
                #print('names_=', names_,file=log)
                #print('distances=', distances_,file=log)
                #print('points_indices=',points_indices_,file=log)
                #print('offsets_=',offsets_,file=log)
                #print('idx=', idx, file=log)
                #
                ni = idx[0]
                if distances_[ni]<1e-6:
                    ni = idx[1]
                nn = points_indices_[ni]
                #print('ni=', ni, 'nn=', nn)
                nn_distance = distances_[ni]
                Z1, Z2 = Znuc[i], Znuc[nn]
                rmt_r = RMT_ratios(Z1, Z2)
                self.Rmt[jatom] = rmt_r[0]*nn_distance
                print('nn_distance=', nn_distance, 'Rmt['+str(jatom)+']=', self.Rmt[jatom], file=log)
                _neighbrs_=[]
                for j in idx:
                    if distances_[j]>1e-5:
                        nn = points_indices_[j]
                        r = self.structure.sites[nn]
                        if Add_oxi:
                            try:
                                oxidation = int(getattr(r.specie, "oxi_state", 0) or 0)
                            except:
                                oxidation=0
                        else:
                            oxidation=0
                        #print('r=', r, 'oxydation=', oxidation)
                        Rj_fract = r.frac_coords + offsets_[j]
                        dRj_fract = Rj_fract - t.frac_coords
                        #Rj = Rj_fract @ self.structure.lattice.matrix
                        #dR = Rj-t.coords
                        dR = dRj_fract @ self.structure.lattice.matrix
                        dst = linalg.norm(dR)
                        #r.specie.symbol,r.species_string
                        #print('    d={:8.5f} d={:8.5f}'.format(distances_[j],dst), r.species, 'off=',offsets_[j], 'dR=', dR)
                        #print('pi=', nn, 'r=', r, 'r+o=', Rj_fract, file=log)
                        print('  ATOM:{:3d} {:5s} AT {:9.5f} {:9.5f} {:9.5f}={:9.5f} {:9.5f} {:9.5f} IS AWAY {:13.6f} {:13.6f} ANG'.format(nn,r.species_string,*dR,*dRj_fract, distances_[j],dst), file=log)
                        _neighbrs_.append([distances_[j],r.species_string,dRj_fract,oxidation])
                self.neighbrs.append(_neighbrs_)

        
    
class Latgen:
    """ This repeats the algorithm of building unit cell inside W2k. This will expose any possible problems or inconsistencies
    in notation between W2k and pymatgen. We will compare the two, and warn the user when inconsistencies occur.
    """
    def __init__(self, strc, case, fout):
        self.pia   = array([2.0*pi/strc.a, 2.0*pi/strc.b, 2.0*pi/strc.c])
        self.alpha = [strc.alpha*pi/180., strc.beta*pi/180., strc.gamma*pi/180.]
        self.ortho = False
        self.br1 = zeros((3,3))   # reciprocal vectors for primitive unit cell
        self.br2 = zeros((3,3))   # reciprocal vectors for conventional unit cell
        if strc.lattyp[:1]=='H': # hexagonal
          print('hexagonal lattice', file=fout)
          self.br2[0,0] = 2.0/sqrt(3.0)
          self.br2[0,1] = 1.0/sqrt(3.0)
          self.br2[1,1] = 1.0
          self.br2[2,2] = 1.0
          self.rvfac = 2.0/sqrt(3.0)
          self.ortho = False
          for j in range(3):
            self.br2[j,:] *= self.pia[j]
          self.br1[:,:] = self.br2[:,:]
        elif strc.lattyp[:1] in ['S','P']: # primitive or simple
          print('primitive or simple lattice', file=fout)
          self.ortho = True
          for i in range(3):
            if( abs(self.alpha[i]-pi/2.) > 0.0001):
              print('alpha['+str(i)+'] not equal 90', file=fout)
              self.ortho = False
              # includes triclinic, monoclinic, simple orthorhombic, tetragonal and cubic
          sinbc = sin(self.alpha[0])
          cosab = cos(self.alpha[2])
          cosac = cos(self.alpha[1])
          cosbc = cos(self.alpha[0])
          wurzel=sqrt(sinbc**2-cosac**2-cosab**2+2*cosbc*cosac*cosab)
          self.br2[0,0]= sinbc/wurzel*self.pia[0]
          self.br2[0,1]= (-cosab+cosbc*cosac)/(sinbc*wurzel)*self.pia[1]
          self.br2[0,2]= ( cosbc*cosab-cosac)/(sinbc*wurzel)*self.pia[2]
          self.br2[1,1]=  self.pia[1]/sinbc
          self.br2[1,2]= -self.pia[2]*cosbc/sinbc
          self.br2[2,2]=  self.pia[2]
          self.rvfac = 1.0/wurzel
          self.br1[:,:] = self.br2[:,:]
        elif strc.lattyp[:1] == 'F': # face centered
          print('face centered lattice', file=fout)
          self.br2[0,0] = -1.0
          self.br2[1,0] =  1.0
          self.br2[2,0] =  1.0
          self.br2[0,1] =  1.0
          self.br2[1,1] = -1.0
          self.br2[2,1] =  1.0
          self.br2[0,2] =  1.0
          self.br2[1,2] =  1.0
          self.br2[2,2] = -1.0
          self.rvfac = 4.0
          self.ortho = True
          for j in range(3):
            self.br2[j,:] *= self.pia[j]
          self.br1[0,0] = self.pia[0]
          self.br1[1,1] = self.pia[1]
          self.br1[2,2] = self.pia[2]
        elif strc.lattyp[:1] == 'B': # body centered
          print('body centered lattice', file=fout)
          self.br2[0,0] = 0.0
          self.br2[1,0] = 1.0
          self.br2[2,0] = 1.0
          self.br2[0,1] = 1.0
          self.br2[1,1] = 0.0
          self.br2[2,1] = 1.0
          self.br2[0,2] = 1.0
          self.br2[1,2] = 1.0
          self.br2[2,2] = 0.0
          self.rvfac = 2.0
          self.ortho = True
          for j in range(3):
            self.br2[j,:] *= self.pia[j]
          self.br1[0,0] = self.pia[0]
          self.br1[1,1] = self.pia[1]
          self.br1[2,2] = self.pia[2]
        elif strc.lattyp[:1] == 'R': # rhombohedral
          print('rhombohedral lattice', file=fout)
          self.br2[0,0] =  1.0/sqrt(3.0)
          self.br2[0,1] =  1.0/sqrt(3.0)
          self.br2[0,2] = -2.0/sqrt(3.0)
          self.br2[1,0] = -1.0
          self.br2[1,1] =  1.0
          self.br2[2,0] =  1.0
          self.br2[2,1] =  1.0
          self.br2[2,2] =  1.0
          #
          self.rvfac = 6.0/sqrt(3.0)
          self.ortho = False
          for j in range(3):
            self.br2[j,:] *= self.pia[j]
          if strc.H2Rtrans:
            # hexagonal conventional cell, but using w2k rhombohedral choice of vectors
            self.br1 = self.br2 @ linalg.inv(strc.hex2rho)
            #self.br1[0,0] = 2.0/sqrt(3.0)
            #self.br1[0,1] = 1.0/sqrt(3.0)
            #self.br1[1,1] = 1.0
            #self.br1[2,2] = 1.0
            #for j in range(3):
            #  self.br1[j,:] *= self.pia[j]
          else:
            self.br1[:,:] = self.br2[:,:]
        elif strc.lattyp[:1] == 'C': # base centered
          print('base centered lattice type', strc.lattyp, file=fout)
          if strc.lattyp[1:3] == 'XZ':   # gamma might not be 90
            ix,iy,iz=0,1,2
          elif strc.lattyp[1:3] == 'YZ': # gamma might not be 90 
            ix,iy,iz=1,0,2  # should be beta not gamma
          elif strc.lattyp[1:3] == 'XY': # beta might not be 90
            ix,iy,iz=0,2,1
          if( abs(self.alpha[iz]-pi/2.0) < 0.0001 ):
            #   orthorombic case
            self.br2[ix,ix] =  1.0
            self.br2[ix,iz] =  1.0
            self.br2[iy,iy] =  1.0
            self.br2[iz,ix] = -1.0
            self.br2[iz,iz] =  1.0
            self.rvfac = 2.0
            self.ortho = True
            for j in range(3):
              self.br2[j,:] *= self.pia[j]
            self.br1[0,0] = self.pia[0]
            self.br1[1,1] = self.pia[1]
            self.br1[2,2] = self.pia[2]
          else:
            #  monoclinic case
            print('alpha['+str(iz)+'] not equal 90 degrees', file=fout)
            sinab = sin(self.alpha[iz])
            cosab = cos(self.alpha[iz])
            self.br2[ix,ix] =  self.pia[ix]/sinab
            self.br2[ix,iy] = -self.pia[iy]*cosab/sinab
            self.br2[ix,iz] =  self.pia[ix]/sinab
            self.br2[iy,iy] =  self.pia[iy]
            self.br2[iz,ix] = -self.pia[iz]
            self.br2[iz,iz] =  self.pia[iz]
            self.rvfac = 2.0/sinab
            self.ortho= False
            #
            self.br1[ix,ix] =  self.pia[ix]/sinab
            self.br1[ix,iy] = -self.pia[iy]*cosab/sinab
            self.br1[iy,iy] =  self.pia[iy]
            self.br1[iz,iz] =  self.pia[iz]
        else:
          print('ERROR wrong lattice=', strc.lattyp, file=fout)
          sys.exit(1)
        
        #  define inverse of cellvolume
        vi = self.rvfac/ (strc.a * strc.b * strc.c)
        self.Vol = 1./vi
        # Calculate the basis vectors of the real lattice
        self.gbas = self.br2[:,:] / (2.0*pi)
        self.rbas = linalg.inv(self.gbas)
        self.rbas_conventional = linalg.inv(self.br1/(2*pi))
        for i in range(3):
            for j in range(3):
                if abs(self.rbas[i,j])<1e-14:
                    self.rbas[i,j]=0
                if abs(self.rbas_conventional[i,j])<1e-14:
                    self.rbas_conventional[i,j]=0
        
        if self.ortho or strc.lattyp=='CXZ':
            self.k2icartes = array((diag(1/self.pia) @ self.br2).round(), dtype=int)
            self.k2cartes = diag(self.pia)
        else:
            # This is a wien2k choice: monoclinic system with CXY and CYZ (but not CXZ) 
            # should not use conventional BZ in defining momentum, but will keep
            # primitive BZ. All other systems use conventional BZ.
            # cif2struct converts all monoclinic CXY to CXZ, hence the only
            # problematic case is CYZ monoclinic structure. Need to check carefuly that it works.
            if strc.H2Rtrans:
                # hexagonal conventional cell, but using w2k rhombohedral choice of vectors
                #R=array([[0, 1, 0], [0, 0, 1], [1, 0, 0]])
                #self.k2icartes = R
                #self.k2cartes = self.br2*linalg.inv(R)
                #self.k2icartes = array((diag(1/self.pia) @ self.br2).round(), dtype=int)
                #self.k2cartes = diag(self.pia)
                self.k2icartes = identity(3,dtype=int)
                self.k2cartes  = self.br2
            else:
                self.k2icartes = identity(3,dtype=int)
                self.k2cartes  = self.br2

        #if writek:
        #    with open(case+'.pkl', 'wb') as f:
        #        pickle.dump(self.k2icartes, f)
        #        pickle.dump(self.k2cartes, f)
            
        #save('k2icartes.npy', self.k2icartes)

        #print('self.br2=', self.br2)
        
        print(('\nw2k_conventional=\n'+('\n'.join(['{:9.5f} '*3+' ==a'+str(i+1) for i in range(3)]))).format(
            *ravel(self.rbas_conventional)), file=fout)
        print('Unit cell volume=', self.Vol, file=fout)
        print('Ortho=', self.ortho, file=fout)
        print(('BR2=\n'+('{:9.5f} '*3+'\n')*3).format(*ravel(self.br2)), file=fout)
        print(('gbas=\n'+('{:9.5f} '*3+'\n')*3).format(*ravel(self.gbas)), file=fout)
        #print(('w2k_primitive(rbas)=\n'+('{:9.5f} '*3+'\n')*3).format(*ravel(self.rbas)), file=fout)
        print(('w2k_primitive(rbas)=\n'+('\n'.join(['{:9.5f} '*3+' ==a'+str(i+1) for i in range(3)]))).format(
            *ravel(self.rbas)), file=fout)
        if strc.H2Rtrans:
            print(('k2icartes=\n'+('{:12.7f} '*3+'\n')*3).format(*ravel(self.k2icartes)), file=fout)
        else:
            print(('k2icartes=\n'+('{:3d} '*3+'\n')*3).format(*ravel(self.k2icartes)), file=fout)
    
def W2k_klist_band(fname, Nt, kpath, k2icartes, k2cartes, H2Rtrans=False, log=sys.stdout):
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
    if H2Rtrans:
        kpath={'\Gamma':[0,0,0],'T':[0.5,0.5,0.5],'M4':[0.37127,0.75747,0.37127],'L':[0,0.5,0],
               'H0':[0.24253,0.5,-0.24253],'M0':[0.37127,0.37127,-0.24253],'FB':[0.5,0.5,0],
               'Q':[0.37127,0.62874,0],'B':[0.5,0.75747,0.24253],"M4'":[0.37127,0.75747,0.37127],
               "B'":[0.24253,0.75747,0.5],"T'":[1/2,1/2,1/2]}
    
    labels = [label for label in kpath]    # it is easier to work here with lists
    Ks = [kpath[label] for label in kpath] # so store labels and K-points
    print('Kpoints in the mesh:', file=log)
    dst = np.zeros(len(kpath))             # first we compute distance between momentum points, so that 
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

    with open(fname, 'w') as fk:
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


def Cif2Struct(fcif, Nkp=300, writek=True, Qmagnetic=True, logfile='cif2struct.log', cmp_neighbors=False, convertH2R=True, ndecimals=3, nradius=4.7):
    log = open(logfile, 'w')
    strc = WStruct()           # w2k structure, for now empty
    strc.ScanCif(fcif, log, writek, Qmagnetic, cmp_neighbors, convertH2R, ndecimals, nradius)    # here is most of the algorithm
    case = os.path.splitext(fcif)[0] # get case
    strc.WriteStruct(case+'.struct', log) # writting out struct file
    lat = Latgen(strc, case, log)  # checking Bravais lattice in wien2k and its consistency.
            
    if writek:
        W2k_klist_band(case+'.klist_band', Nkp, strc.kpath, lat.k2icartes, lat.k2cartes, strc.H2Rtrans, log)
    return (strc, lat)


class Cmp2Sort:
    def __init__(self, distances_,names_,center_name):
        self.distance=distances_
        self.names=names_
        self.center_name = center_name
    def __call__(self, i, j):
        # sorts by distance if distance is different.
        # but if distances are nearly equal, we want to first have neigbors of different element
        # not the same element
        if abs(self.distance[i]-self.distance[j])>1e-5:
            return self.distance[i]-self.distance[j]
        else:
            if self.names[i]==self.center_name and self.names[j]!=self.center_name:
                return 1
            elif self.names[j]==self.center_name and self.names[i]!=self.center_name:
                return -1
            else:
                return (self.names[i]>self.names[j]) - (self.names[i]<self.names[j])
        
if __name__ == '__main__':
    usage = """usage: %cif2struct.py [ options ] filename.cif

    Converts cif file (Crystallographic Information File) to struct file, 
    as needed for wien2k calculation. It can also produce high symmetry 
    path in momentum space in case.klist_band per wien2k needs.
    """
    parser = optparse.OptionParser(usage)
    parser.add_option('-w', '--writek',  dest='wkp', action='store_true', default=False, help="weather to write case.klist_band high symmetry k-path (default False)")
    parser.add_option('-N', '--Nkp',     dest='Nkp', type='int', default=300, help="number of k-points along the high symmetry path")
    parser.add_option('-l', '--log',     dest='log', type='str', default='cif2struct.log', help="info file")
    parser.add_option('-n', '--neigh',   dest='neigh', action='store_false', default=True, help="compute neighbors and evaluate Rmt")
    parser.add_option('-m', '--magnet',  dest='Qmagnetic', action='store_true', default=False, help="should produce magnetic dft struct that allows broken symmetry from mcif")
    parser.add_option('-H', '--hexagonal',  dest='keepH', action='store_true', default=False, help="switches off hexagonal to rhombohedral conversion")
    parser.add_option('-p', '--ndecimals',  dest='ndecimals', type='int', default=3, help="precision when determining the equivalency of atoms we take nn distance with precision ndecimals. Compatible with w2k nn method.")
    parser.add_option('-R', '--nradius',  dest='nradius', type='float', default=4.7, help="How far in AA should nearest neighbors be checked. Compatible with w2k nn method.")
    
    # Next, parse the arguments
    (options, args) = parser.parse_args()
    if len(args)!=1:
        print('Need exactly one argument: the name of cif file')
        sys.exit(1)
    fcif = args[0]

    #print('fcif=', fcif)
    #print('options=', options.Nkp, options.wkp)
    Cif2Struct(fcif, Nkp=options.Nkp, writek=options.wkp, Qmagnetic=options.Qmagnetic, logfile=options.log,
                   cmp_neighbors=options.neigh, convertH2R=not options.keepH, ndecimals=options.ndecimals,nradius=options.nradius)
