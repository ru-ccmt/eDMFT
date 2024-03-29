#!/usr/bin/env python
""" Reads w2k structure file and extract information about the w2k lattice.
    Can also write back structure file.
"""
# @Copyright 2007 Kristjan Haule
from numpy import *
from numpy import linalg
import sys

class Struct:
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
        
    def ReadStruct(self, fname, log=sys.stdout):
        f = open(fname, 'r')
        self.title = next(f).strip()
        line = next(f)
        self.lattyp = line[0:4].strip()
        self.nat = int(line[27:30])
        for i in range(4): # it seems w2k does not consistently write space group number and name into the same location
            iend=35-i
            if line[30:iend].strip().isdigit():
                self.sgnum = int(line[30:iend])
                break
        #self.sgnum = int(line[30:35]) # had to change from 30:33 to 30:35 for Ta2NiSe5M.struct
        self.sgname = line[iend:].strip()
        #print('nat=', self.nat, 'sgnum=', self.sgnum, 'sgname=', self.sgname)
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

        self.aZnuc=[]
        for iat in range(self.nat):
            self.aZnuc += [self.Znuc[iat]]*self.mult[iat]

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
        
    def WriteStruct(self, fname=None, log=sys.stdout):
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
            
        with open(fname, 'w') if fname!=None else sys.stdout as f:
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

    def flat(self, notflat):
        '''Return a flat view of given data as a list.
        Example: if w.mult = [2,4] and w.aname = ['V', 'O']
        w.flatten(w.aname) -> ['V', 'V', 'O', 'O', 'O', 'O']'''
        from functools import reduce
        import operator
        if notflat is self.pos:
            listoflists = self.pos
        else:
            listoflists = [[elem]*mult for elem,mult in zip(notflat, self.mult)]
        return reduce(operator.add, listoflists)

    def __repr__(self):
        lines = self.title+'\n'
        lines += '{:3s} LATTICE,NONEQUIV.ATOMS:{:3d} {:3d} {:8s}\n'.format(self.lattyp,self.nat,self.sgnum,self.sgname)
        lines += 'MODE OF CALC=RELA unit=bohr\n'
        lines += '{:10.6f}{:10.6f}{:10.6f}{:10.6f}{:10.6f}{:10.6f}\n'.format(self.a,self.b,self.c,self.alpha,self.beta,self.gamma)
        index = 0
        for jatom in range(self.nat):
            lines += 'ATOM{:4d}: X={:10.8f} Y={:10.8f} Z={:10.8f}\n'.format(self.iatnr[jatom],*self.pos[index,:])
            lines += '          MULT={:2d}          ISPLIT={:2d}\n'.format(self.mult[jatom],self.isplit[jatom])
            index += 1
            for m in range(1,self.mult[jatom]):
                lines += '    {:4d}: X={:10.8f} Y={:10.8f} Z={:10.8f}\n'.format(self.iatnr[jatom],*self.pos[index,:])
                index += 1
            R0s='{:10.9f}'.format(self.r0[jatom])
            lines += '{:10s} NPT={:5d}  R0={:10s} RMT={:10.5f}   Z:{:10.5f}\n'.format(self.aname[jatom], self.jrj[jatom], R0s[1:], self.rmt[jatom], self.Znuc[jatom])
            lines += 'LOCAL ROT MATRIX:   {:10.7f}{:10.7f}{:10.7f}\n'.format(*self.rotloc[jatom][:,0])
            lines += '                    {:10.7f}{:10.7f}{:10.7f}\n'.format(*self.rotloc[jatom][:,1])
            lines += '                    {:10.7f}{:10.7f}{:10.7f}\n'.format(*self.rotloc[jatom][:,2])
        nsym = shape(self.tau)[0]
        lines += '{:4d}      NUMBER OF SYMMETRY OPERATIONS\n'.format(nsym)
        for iord in range(nsym):
            for j in range(3):
                lines += '{:2d}{:2d}{:2d} {:10.8f}\n'.format(*self.timat[iord,j,:], self.tau[iord,j])
            lines += '      {:2d}\n'.format(iord+1)
        return lines
        
if __name__ == '__main__':
    if len(sys.argv)<2:
        print('Give input structure name case.struct')
        sys.exit(0)
    fname = sys.argv[1]
    strc = Struct()
    strc.ReadStruct(fname)
    strc.WriteStruct()
    
