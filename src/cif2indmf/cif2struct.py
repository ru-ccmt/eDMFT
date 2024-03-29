#!/usr/bin/env python
""" Converts cif file (Crystallographic Information File) to struct file, 
    as needed for wien2k calculation. It can also produce high symmetry 
    path in momentum space in case.klist_band per wien2k needs.
"""
# @Copyright 2007 Kristjan Haule
from math import gcd
import re, sys, os
import optparse
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
        hex2rho = hex2ort @ ort2rho
        self.rho2hex = linalg.inv(hex2rho)
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
        
    def ScanCif(self, filein, log=sys.stdout, writek=False, cmp_neighbors=False):
        """ Here is the majority of the algorithm to convert cif2struct.
        We start by allocating CifParser_W2k, which reads cif using paymatgen, and 
        extracts necessary information.
        Below we than arrange this information in the way W2k needs it in struct file.
        
        """
        cif = CifParser_W2k(filein, log, cmp_neighbors)
        
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
            self.kpath = kpath._kpath['kpoints']
            
        for i,spec in enumerate(cif.w2k_coords):
            print(i, spec, cif.w2k_coords[spec], file=log)
            
        print('a=', self.a, 'b=', self.b, 'c=', self.c, 'alpha=', self.alpha, 'beta=', self.beta, 'gamma=', self.gamma, file=log)
        print('sgname=', self.sgname, 'sgnum=', self.sgnum, 'lattty=', self.lattyp, 'Nsym=', nsym, file=log)
        print('Znuc=', ['Z['+spec+']='+str(self.Znuc[i]) for i,spec in enumerate(cif.w2k_coords)], file=log)
        print('high symmetry kpath=', self.kpath, file=log)
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
        RHtransf = False
        RHexeption = False
        if self.sgnum in [146,148,155,160,161,166,167]: # trigonal groups with R : [R3,R-3,R32,R3m,R3c,R-3m,R-3c]
            small = 1e-6
            if abs(self.alpha-90)<small and abs(self.beta-90)<small and abs(self.gamma-120) < small: # trigonal, but w2k pretends it is hexagonal
                self.lattyp = 'H'
                RHexeption = True
            else:
                # self.lattyp should be 'R'
                RHtransf = True
                aa = self.a*2.0*sin(self.alpha/180*pi/2)
                bb = aa
                cc = 3*sqrt(self.a**2 - aa**2/3)
                self.a, self.b, self.c = aa, bb, cc
                self.alpha, self.beta, self.gamma = 90., 90., 120.
                print('WARNING: We are changing trigonal group '+str(self.sgnum)+' to equivalent hexagonal settings. However, we choose alpha,beta,gamma=90,90,120, which might not be appropriate.', file=log)
                print('If errors, please correct this algorithm in cif2struct.py', file=log)
                
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
        if self.lattyp in ['B', 'F', 'CXY', 'CXZ', 'CYZ']:
            ll = 0
            the_same = TheSame(self.lattyp)
            for spec in cif.w2k_coords:
                all_coords = cif.w2k_coords[spec]
                #print(spec,'all_coords=', all_coords, file=log)
                #print('shifts=')
                #for i in range(len(the_same.tv)):
                #    print(i, the_same.tv[i])
                kept = [0]
                for i in range(1,len(all_coords)):
                    Qunique = True
                    for kt in kept:
                        if the_same(all_coords[i], all_coords[kt]):
                            Qunique = False
                            break
                    if Qunique:
                        kept.append(i)
                to_remove=set(range(len(all_coords)))-set(kept)
                print(spec, ': ', 'out of ', len(cif.w2k_coords[spec]), 'atoms we keep=', kept, file=log)
                for i in sorted(to_remove, reverse=True):
                    del cif.w2k_coords[spec][i]
                print(spec, ':', 'end up with', len(cif.w2k_coords[spec]), file=log)

        print('After removing atoms not in primitive cell, we have the following cif.w2k_coords list:', file=log)
        for spec in cif.w2k_coords:
            for ipos,pos in enumerate(cif.w2k_coords[spec]):
                print('name='+spec, 'Z='+str(cif.Z_element[spec]), pos, file=log)

        if False: # This is something cif2struct does when it can not find symmetry operations in cif file.
            if RHexeption:
                self.HRtransf = True
            if RHtransf:
                self.HRtransf = True
                
                for spec in cif.w2k_coords:
                    for pos in cif.w2k_coords[spec]:
                        pos = dot(pos, self.rho2hex) % 1.0
            
                # stores something in temp file
            # a lot more code in spacegroup.f
        #usespg = False
            
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
            
        if self.lattyp in ['CXY','CXZ','CYZ','F','B']:
            # for these Bravais lattices the glide planes should be eliminated from w2k symmetry operations
            # because w2k adds them in by hand. Hence it is better to generate symmetry operations by w2k initialization
            # then remove half of the symmetry operations here.
            self.tau=[]
            
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
                    self.pos[i,:] = np.dot(self.pos[i,:], self.hex2rho)
                    self.pos[i,:] = np.where(self.pos[i,:] < 0.0, self.pos[i,:] + 1.0, self.pos[i,:])
                    self.pos[i,:] = np.where(self.pos[i,:] > 1.0, self.pos[i,:] - 1.0, self.pos[i,:])
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
            print('{:4d}      NUMBER OF SYMMETRY OPERATIONS'.format(nsym), file=f)
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
                    
        
def Z2r0(Z):
    "This is the choice for real space integration mesh (point closest to nucleous) in w2k."
    if (Z>71): return 5e-6
    if (36<Z<=71): return 1e-5
    if (18<Z<=36): return 5e-5
    if (Z<=18): return 1e-4
    return 5e-6

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
    def __init__(self, fname, log=sys.stdout, cmp_neighbors=False):
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
                if '_symmetry_Int_Tables_number' in data:
                    self.sgnum = data['_symmetry_Int_Tables_number']
                elif '_symmetry_space_group_IT_number' in data:
                    self.sgnum = data['_symmetry_space_group_IT_number']
                elif '_space_group_IT_number' in data:
                    self.sgnum = data['_space_group_IT_number']
                
                print('Information directly extracted from cif file:', file=log)
                for k2, v2 in data.items():
                    print('key=', k+':'+k2,' value=', v2, file=log)
                print(file=log)
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

        #print('indx_by_element=', indx_by_element)
        print('Number of symmetry operations available=', len(self.parser.symmetry_operations), file=log)
        groups=[] # Will contain groups of equivalent sites
        acoords = [site.frac_coords%1.0 for site in self.structure.sites] # all coordinates but always inside home unit cell
        for spec in indx_by_element:  # Next we loop through element in this structure
            indices = indx_by_element[spec] # indices on sites of the same element
            grp = [[indices[0]]]  # we have just one group with the first entry in coords. All entries in coords will have index in grp
            for ii in indices[1:]: # loop over possibly equivalent sites (of the same element)
                coord = acoords[ii]    # fractional coordinate of that site self.structure.sites[ii]
                Qequivalent=False  # for now we thing this might not be equivalent to any other site in grp
                for op in self.parser.symmetry_operations: # loop over all symmetry operations
                    ncoord = op.operate(coord) % 1.0       # apply symmetry operation to this site, and produce ncoord
                    for ty in grp:                         # if coord is equivalent to any existing group in grp, say grp[i]==ty, then one of ncoord should be equal to first entry in grp[i]
                        coord0 = acoords[ty[0]]            # the first entry in the group grp[i][0]
                        #print('coord0=', ty[0], coord0)
                        if sum(abs(ncoord-coord0))<1e-5:   # is it equal to current coord when a symmetry operation is applied?
                            #print(' ncoord=', ncoord, 'is equivalent to typ', ty) # yes, coord is in group ty
                            Qequivalent=True               # we found which group it corresponds to
                            ty.append(ii)                  # grp[i]=ty is extended with coord[i]
                            break                          # loop over groups and also symmetry operations can be finished
                    if Qequivalent:                        # once we find which group this site corresponds to, we can finish the loop over symmetry operations
                        break
                if not Qequivalent:                        # if this site (coord) is not equivalent to any existing group in grp, than we create a new group with a single entry [i]
                    grp.append([ii])
                #print(spec+' grp=', grp)
            groups.append(grp)  # all groups for this element
        # Since groups is only index array to coord_by_element, we will rather create a dictionary self.w2k_coords={},
        # which is easier to use. The keys will be name of the site, which might need to be modified when
        # multiple inequivalent sites have the same name.
        # The values in the dictionary are all equivalent sites,i.e., their fractional coordinates
        ResortToCif=True
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
                # Instead of comparing with all coordinates from the cif file, we only compare with those
                # that start with the same letter, i.e., likely correspond to the same element.
                which_icif = [i for i in range(len(self.cname)) if self.cname[i][0]==spec[0]] # those are probably for the same element
                grp = groups[idx] # group of atoms of the same element. But not necessary equivalent. For example, all oxygen atoms in the structure.
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
            #print(groups_resort)
            first_sites_index=[]
            self.w2k_coords={}
            for ii,group in enumerate(groups_resort):
                psite = self.structure.sites[group[0]]
                self.w2k_coords[cname[ii]] = [self.structure.sites[j].frac_coords for j in group]
                first_sites_index.append( group[0] )
                #print(cname[ii], coords, psite.label, self.cname[ii], psite.species, psite.species_string)
                #self.Z_element[cname[ii]] = psite.specie.Z

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
            for grp in groups:
                for gr in grp:
                    jatom = len(first_sites_index)
                    first_sites_index.append(gr[0])

        self.neighbrs=[]
        if cmp_neighbors:
            first_sites = [self.structure.sites[i] for i in first_sites_index]
            self.structure.add_oxidation_state_by_guess()
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
                oxi_st = int(getattr(t.specie, "oxi_state", 0) or 0)
                self.oxi_state[jatom]=oxi_st
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
                        oxidation = int(getattr(r.specie, "oxi_state", 0) or 0)
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
          self.rvfac = 6.0/sqrt(3.0)
          self.ortho = False
          for j in range(3):
            self.br2[j,:] *= self.pia[j]
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
            self.k2icartes = identity(3,dtype=int)
            self.k2cartes  = self.br2

        #if writek:
        #    with open(case+'.pkl', 'wb') as f:
        #        pickle.dump(self.k2icartes, f)
        #        pickle.dump(self.k2cartes, f)
            
        #save('k2icartes.npy', self.k2icartes)
        print(('\nw2k_conventional=\n'+('\n'.join(['{:9.5f} '*3+' ==a'+str(i+1) for i in range(3)]))).format(
            *ravel(self.rbas_conventional)), file=fout)
        print('Unit cell volume=', self.Vol, file=fout)
        print('Ortho=', self.ortho, file=fout)
        print(('BR2=\n'+('{:9.5f} '*3+'\n')*3).format(*ravel(self.br2)), file=fout)
        print(('gbas=\n'+('{:9.5f} '*3+'\n')*3).format(*ravel(self.gbas)), file=fout)
        print(('w2k_primitive(rbas)=\n'+('{:9.5f} '*3+'\n')*3).format(*ravel(self.rbas)), file=fout)
        print(('k2icartes=\n'+('{:3d} '*3+'\n')*3).format(*ravel(self.k2icartes)), file=fout)
    
def W2k_klist_band(fname, Nt, kpath, k2icartes, k2cartes, log=sys.stdout):
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


def Cif2Struct(fcif, Nkp=300, writek=True, logfile='cif2struct.log', cmp_neighbors=False):
    log = open(logfile, 'w')
    strc = WStruct()           # w2k structure, for now empty
    strc.ScanCif(fcif, log, writek, cmp_neighbors)    # here is most of the algorithm
    case = os.path.splitext(fcif)[0] # get case
    strc.WriteStruct(case+'.struct', log) # writting out struct file
    lat = Latgen(strc, case, log)  # checking Bravais lattice in wien2k and its consistency.

    #corr = range(strc.nat)
    #
    #first_atom = zeros(len(strc.mult),dtype=int)
    #for jatom in range(len(strc.mult)-1):
    #    first_atom[jatom+1] = first_atom[jatom] + strc.mult[jatom]
    #
    #if cmp_neighbors:
    #    matrix_w2k = lat.rbas_conventional
    #    matrix_vesta = strc.Matrix(vesta=True)
    #
    #    to_frac = linalg.inv(matrix_w2k)
    #    vesta_vs_w2k = to_frac @ matrix_vesta
    #    #print('vesta_vs_w2k=', vesta_vs_w2k)
    #
    #    #print('VESTA=')
    #    #for i in range(3): print( ('{:12.8f} '*3).format(*matrix_vesta[i,:]/strc.tobohr))
    #        
    #    to_fract = linalg.inv(matrix_w2k)
    #    
    #    indx=0
    #    for jatom in range(len(strc.neighbrs)):
    #        print('ATOM: {:3d} {:3s} AT {:9.5f} {:9.5f} {:9.5f}'.format(jatom+1,strc.aname[jatom],*strc.pos[indx]),file=log)
    #        for j,ngh in enumerate(strc.neighbrs[jatom]):
    #            dRj = ngh[2] @ matrix_w2k
    #            dRp = ngh[2] @ matrix_vesta
    #            dst = linalg.norm(dRj)
    #            dst2 = linalg.norm(dRp)
    #            jname = ngh[1]
    #            if len(jname)>3: jname=jname[:3]
    #            print('  ATOM: {:3s} AT {:9.5f} {:9.5f} {:9.5f} IS AWAY {:13.6f} {:13.6f} {:13.6} ANG'.format(jname,*dRj,ngh[0]*strc.tobohr,dst,dst2), file=log)
    #        indx += strc.mult[jatom]
    #    for jatom in corr:
    #        name0 = Element_name(strc.neighbrs[jatom][0][1])
    #        n=0
    #        for ngh in strc.neighbrs[jatom]:
    #            #print('n=',n,'ngh=',ngh,'el=',Element_name(ngh[1]),'name0=',name0)
    #            if Element_name(ngh[1])!=name0:
    #                break
    #            n+=1
    #        #print('n=',n)
    #        if n>1:
    #            neighbrs = strc.neighbrs[jatom][:n]
    #        else: # n==1
    #            if name0==Element_name(strc.aname[jatom]): # skip the first atom and take the rest
    #                print('WARNING: At analizing', strc.aname[jatom], 'at', '[{:5.3f},{:5.3f},{:5.3f}]'.format(*strc.pos[first_atom[jatom]]),
    #                          'we disregard the first neighbor {:s} at [{:5.3f},{:5.3f},{:5.3f}]'.format(strc.neighbrs[jatom][0][1],*strc.neighbrs[jatom][0][2]), file=log)
    #                name0 = Element_name(strc.neighbrs[jatom][1][1])
    #                n=0
    #                for ngh in strc.neighbrs[jatom][1:]:
    #                    if Element_name(ngh[1])!=name0:
    #                        break
    #                    n+=1
    #                neighbrs = strc.neighbrs[jatom][1:n+1]
    #            
    #        print('Analizing', strc.aname[jatom], 'at', strc.pos[first_atom[jatom]], 'with N=', n, file=log)
    #        R = FindCageBasis(neighbrs, matrix_w2k, log)
    #        if R is not None:
    #            Rv = R @ vesta_vs_w2k
    #            Rf = R @ to_frac
    #
    #            print('Rotation to input into case.indmfl by locrot=-1 : ', file=log)
    #            print(file=log)
    #            for i in range(3): print( ('{:12.8f} '*3).format(*R[i,:]), file=log)
    #            print(file=log)
    #            print('Rotation in fractional coords : ', file=log)
    #            print(file=log)
    #            for i in range(3): print( ('{:12.8f} '*3).format(*Rf[i,:]*strc.tobohr), file=log)
    #            print(file=log)
    #            print('Rotation for VESTA coords : ', file=log)
    #            print(file=log)
    #            for i in range(3): print( ('{:12.8f} '*3).format(*Rv[i,:]), file=log)
    #            print(file=log)
    #        
    #        
            
    if writek:
        W2k_klist_band(case+'.klist_band', Nkp, strc.kpath, lat.k2icartes, lat.k2cartes, log)
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
    # Next, parse the arguments
    (options, args) = parser.parse_args()
    if len(args)!=1:
        print('Need exactly one argument: the name of cif file')
        sys.exit(1)
    fcif = args[0]

    #print('fcif=', fcif)
    #print('options=', options.Nkp, options.wkp)
    Cif2Struct(fcif, Nkp=options.Nkp, writek=options.wkp, logfile=options.log, cmp_neighbors=options.neigh)
