#!/usr/bin/env python
""" Given w2k structure file it constructs w2k lattice matrix. 
    For pymatgen or VESTA the matrix is very simple to construct, 
    and can be obtained from wstruct.Struct.Matrix()
    For w2k we need Latgen coded below
"""
from numpy import *

# @Copyright 2007 Kristjan Haule
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

if __name__ == '__main__':
    import sys, os
    from wstruct import Struct
    if len(sys.argv)<2:
        print('Give input structure name case.struct')
        sys.exit(0)
    fname = sys.argv[1]
    strc = Struct()
    strc.ReadStruct(fname)
    case = os.path.splitext(fname)[0] # get case
    latgen = Latgen(strc, case, sys.stdout)
     
