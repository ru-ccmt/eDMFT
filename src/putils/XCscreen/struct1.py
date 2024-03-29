# @Copyright 2007 Kristjan Haule
import sys
from scipy import *
from functools import reduce

class Struct:
    def __init__(self, case):
        self.case = case
        self.extn = 'struct'
        self.parse()

    def parse(self):
        '''The file is structured for reading by fortran code,
        so data is positioned by line and by column.'''
        f = open(self.case + '.' + self.extn, 'r')
        self.title = f.next().strip()

        line = next(f)
        self.lattice = line[0:4].strip()
        self.nat = int(line[27:30])
        self.latticename = line[30:].strip()

        self.mode = f.next()[13:17]

        line = next(f)
        self.a, self.b, self.c, self.alpha, self.beta, self.gamma = [float(line[i:i+10]) for i in range(0,60,10)]

        self.iatom  = []
        self.mult   = []
        self.isplit = []
        self.aname  = []
        self.npt    = []
        self.r0     = []
        self.rmt    = []
        self.Znuc   = []
        self.pos    = []
        self.rotloc = []

        for iat in range(self.nat):
            line = next(f)
            self.iatom.append( int(line[4:8]) )
            pos = [[float(line[col:col+10]) for col in 12+13*arange(3)]]
            
            line = next(f)
            mult = int(line[15:17])
            self.mult.append(mult)

            self.isplit.append( int(line[34:37]) )

            for mu in range(mult-1):
                line = next(f)
                pos.append( [float(line[col:col+10]) for col in 12+13*arange(3)] )

            self.pos.append(pos)
                
            line = next(f)
            self.aname.append( line[0:10].strip() )
            self.npt.append( int(line[15:20]) )
            self.r0.append( float(line[25:35]) )
            self.rmt.append( float(line[40:50]) )
            self.Znuc.append( int(float(line[55:65])) )

            rt=[]
            for i in range(3):
                line = next(f)
                rt.append( [float(line[col:col+10]) for col in 20+10*arange(3)] )
            self.rotloc.append(rt)

        self.Ng = int(f.next().split()[0])
        self.iz=[]
        self.tau=[]
        for i in range(self.Ng):
            tiz=zeros((3,3), dtype=float)
            ttau=zeros(3, dtype=float)
            for j in range(3):
                line = next(f)
                tiz[j,:] = array([int(line[0:2]), int(line[2:4]), int(line[4:6])])
                ttau[j] = float(line[6:16])
            self.iz.append(tiz)
            self.tau.append(ttau)
            next(f)
        #print self.iz
        
        f.close()

    def debugprint(self, f=sys.stdout):
        print('***** Structure File Contents *****', file=f)
        print('title =', self.title, file=f)
        print('Number of sorts, nat =', self.nat, file=f)
        print('unit cell parameters (a,b,c) =', self.a, self.b, self.c, file=f)
        print('angles (alpha, beta, gamma) =', self.alpha, self.beta, self.gamma, file=f)


        printlist = [
            (self.aname,  'Atom name, aname ='),
            (self.iatom,  'Atom index, iatom ='),
            (self.mult,   'Number of atoms of this type, mult ='),
            (self.isplit, 'Symmetry of the atom, isplit ='),
            (self.npt,    'Number of radial points, npt ='),
            (self.r0,     'First point in radial mesh, r0 ='),
            (self.rmt,    'Muffin tin radius, rmt ='),
            (self.Znuc,   'Atomic number, Znuc ='),
            ]

        for i in range(self.nat):
            print('---------- atom type', i, '------------', file=f)

            for var,docstr in printlist:
                print(docstr, var[i], file=f)

            print('Position(s) inside the unit cell:', file=f)
            for m in range(self.mult[i]):
                print('    pos =', self.pos[i][m], file=f)

            print('Local rotation matrix:', file=f)
            for row in self.rotloc[i]:
                print(row, file=f)

        print('***********************************', file=f)

    def flat(self, notflat):
        '''Return a flat view of given data as a list.
        Example: if w.mult = [2,4] and w.aname = ['V', 'O']
        w.flatten(w.aname) -> ['V', 'V', 'O', 'O', 'O', 'O']'''
        if notflat is self.pos:
            listoflists = self.pos
        else:
            listoflists = [[elem]*mult for elem,mult in zip(notflat, self.mult)]
        return reduce(operator.add, listoflists)


def DirectVectors(BR,RN):
    v0 = cross(BR[1],BR[2])
    v1 = cross(BR[2],BR[0])
    v2 = cross(BR[0],BR[1])
    vol = dot(v0,BR[0])
    return dot(RN, array([v0/vol,v1/vol,v2/vol]))
    
