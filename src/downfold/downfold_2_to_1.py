#!/usr/bin/env python
from scipy import *
import struct, sys, re
from scipy import linalg
import optparse
import rdU
from utils import W2kEnvironment, Ry_in_eV
import indmffile
import copy

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

        line = f.next()
        self.lattice = line[0:4].strip()
        self.nat = int(line[27:30])
        self.latticename = line[30:].strip()

        self.mode = f.next()[13:17]

        line = f.next()
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
            line = f.next()
            self.iatom.append( int(line[4:8]) )
            pos = [[float(line[col:col+10]) for col in 12+13*arange(3)]]
            
            line = f.next()
            mult = int(line[15:17])
            self.mult.append(mult)

            self.isplit.append( int(line[34:37]) )

            for mu in range(mult-1):
                line = f.next()
                pos.append( [float(line[col:col+10]) for col in 12+13*arange(3)] )

            self.pos.append(pos)
                
            line = f.next()
            self.aname.append( line[0:10].strip() )
            self.npt.append( int(line[15:20]) )
            self.r0.append( float(line[25:35]) )
            self.rmt.append( float(line[40:50]) )
            self.Znuc.append( int(float(line[55:65])) )

            rt=[]
            for i in range(3):
                line = f.next()
                rt.append( [float(line[col:col+10]) for col in 20+10*arange(3)] )
            self.rotloc.append(rt)

        self.Ng = int(f.next().split()[0])
        self.iz=[]
        self.tau=[]
        for i in range(self.Ng):
            tiz=zeros((3,3), dtype=float)
            ttau=zeros(3, dtype=float)
            for j in range(3):
                line = f.next()
                tiz[j,:] = array([int(line[0:2]), int(line[2:4]), int(line[4:6])])
                ttau[j] = float(line[6:16])
            self.iz.append(tiz)
            self.tau.append(ttau)
            f.next()
        #print self.iz
        
        f.close()

    def debugprint(self, f=sys.stdout):
        print >> f, '***** Structure File Contents *****'
        print >> f, 'title =', self.title
        print >> f, 'Number of sorts, nat =', self.nat
        print >> f, 'unit cell parameters (a,b,c) =', self.a, self.b, self.c
        print >> f, 'angles (alpha, beta, gamma) =', self.alpha, self.beta, self.gamma


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
            print >> f, '---------- atom type', i, '------------'

            for var,docstr in printlist:
                print >> f, docstr, var[i]

            print >> f, 'Position(s) inside the unit cell:'
            for m in range(self.mult[i]):
                print >> f, '    pos =', self.pos[i][m]

            print >> f, 'Local rotation matrix:'
            for row in self.rotloc[i]:
                print >> f, row

        print >> f, '***********************************'

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
    #v0 = cross(BR[1],BR[2])
    #v1 = cross(BR[2],BR[0])
    #v2 = cross(BR[0],BR[1])
    #vol = dot(v0,BR[0])
    #array([v0/vol,v1/vol,v2/vol])
    return dot(RN, linalg.inv(BR).transpose() )
    
def Print(fh, R, tHr):
    R0 = R # dot(DVinv,R)
    Rs=[]
    for i in range(3):
        if abs(round(R0[i])-R0[i])<1e-3:
            Rs.append( "%d" % int(round(R0[i])) )
        else:
            Rs.append("%.1f" % R0[i])

    print >> fh, 'Hopping[(%s,%s,%s)]=[' % tuple(Rs)
    for l1 in range(shape(tHr)[0]):
        print >> fh, '[',
        for l2 in range(shape(tHr)[1]):
            print >> fh, "%8.4f+%8.4fj," % (tHr[l1,l2].real, tHr[l1,l2].imag),
        print >> fh, '],'
    print >> fh, ']'


def ReadBRs(case):
    "reads BR1 and BR2 from case.rotlm"
    fr = open(case+'.rotlm') 
    BR = zeros((2,3,3), dtype=float)
    for ii in range(2):
        fr.next()
        for i in range(3):
            BR[ii,:,i] = map(float,fr.next().split())
    
    BRX = array(linalg.inv(BR[1]) * matrix(BR[0]))
    return (BR,BRX)


def FindMaxMinBand(Ek, Emin, Emax):
    "Finds minimum and maximum band to be used in downfolding"

    istart=None
    iend=None
    for i in range(len(Ek)):
        if Ek[i]>Emin:
            istart = i
            break
        
    for i in range(istart,len(Ek)):
        if Ek[i]>Emax:
            iend = i
            break
    return array([istart, iend])
            
def TryAddBands(Nse_global, Ek, Emin, Emax, SmearN, SmearX):
    
    Nse = Nse_global[:]
    
    Em = Ek[Nse[0]-1] # closest energy below the cutoff
    Ep = Ek[Nse[1]]   # closest energy above the cutoff
    
    print Em, Ep

    EmCanAdd, EpCanAdd = False, False
    
    if Emin-Em < SmearN : EmCanAdd = True
    if Ep-Emax < SmearX : EpCanAdd = True

    Succ = True
    
    if EmCanAdd and EpCanAdd:
        if Emin-Em < Ep-Emax:
            Nse[0]-=1
            print 'Both could be added. Adding below.'
        else:
            Nse[1]+=1
            print 'Both could be added. Adding above.'                            
    elif EmCanAdd :
        Nse[0]-=1
        print 'Can add below'
    elif EpCanAdd :
        Nse[1]+=1
        print 'Can add above'
    else:
        print 'Can not find enough bands in your interval. Increase Emin, Emax or smearings'
        Succ = False
        
    print Nse
    
    return (Nse, Succ)


# RS = (i0,i1,i2)*DV[:,:]
# RS*DV^-1 = (i0,i1,i2)

def GiveRs(div, DV):
    "Gives real space vectors"
    Rs=[]
    for i0 in range(div[0]/2,-div[0]/2,-1):
        for i1 in range(div[1]/2,-div[1]/2,-1):
            for i2 in range(div[2]/2,-div[2]/2,-1):
                Rs.append(DV[0]*i0 + DV[1]*i1 + DV[2]*i2)
    return Rs


    

if __name__ == '__main__':
    """ Takes the DMFT projector Udmft.0 and LDA eigenvalues (case.energy) to produces downfolding hopping parameters.
    It also requires two parameters specifying the bands used in downfolding. The lower cutoff is determined by energy
    '--Emin' (counted from EF in eV), while the upped cutoff is determined by the number of bands '--nbnd'.
    """

    usage = """usage: %prog [ options ]
    
    The script takes the DMFT projector Udmft.0 and LDA eigenvalues (case.energy)
    and produces downfolded hopping parameters. It requires parameters
    specifying bands used in downfolding. The lower/upper cutoff is determined
    by energy '--Emin' or '-N' /'--Emax' or '-X' (counted from EF in eV).
    Additional bands are added if needed in the interval
    '--SmearN' or '-n' / '--SmerX' or '-x'

    --------------------   Emax+SmearX
      added if needed
    ---------------------  Emax

      bands always used in downfolding

    ---------------------  Emin
      added if needed
    ---------------------  Emin-SmearN
    """
    
    # Finds what is case
    w2k = W2kEnvironment()

    
    parser = optparse.OptionParser(usage)
    parser.add_option("-N", "--Emin",  dest="Emin", type="float", default=None, help="Bands above Emin (counted from EF in eV) will be consider in downfolding.")
    parser.add_option("-X", "--Emax",  dest="Emax", type="float", default=None, help="Bands below Emax (counted from EF in eV) will be consider in downfolding.")
    parser.add_option("-n", "--SmearN",  dest="SmearN", type="float", default=1., help="When determining Emin cutoff, we use smearing in eV if needed")
    parser.add_option("-x", "--SmearX",  dest="SmearX", type="float", default=1., help="When determining Emax cutoff, we use smearing in eV if needed")
    parser.add_option("-s", "--smalls",  dest="smalls", type="float", default=0.1, help="The smallest acceptable singular value")
    # Next, parse the arguments
    (options, args) = parser.parse_args()
    
    if options.Emin is None:
        print 'You have to specify minimum band energy! Option -E or --Emin'
        sys.exit(1)
    if options.Emax is None:
        print 'You have to specify maximum band energy! Option -M or --Emax'
        sys.exit(1)
        
    print 'case=%s Emin=%s, Emax=%s SmearN=%s SmearX=%s' %  (w2k.case, options.Emin, options.Emax, options.SmearN, options.SmearX)

    # File handlers to exchange files with Fortran
    # These numbers are not relevant, but need to be set
    fhp = 233
    fhe = 60
    Ntry=10    
    
    # Reads indmf-file
    inl = indmffile.Indmfl(w2k.case)
    inl.read()
    EF = inl.EF
    print 'Fermi energy is ', EF
    # atoms in the projection for downfolding

    atms = inl.atoms.keys()
    print 'atoms in the projection for downfolding=', atms
    # Finds actual positions of these atoms
    pos = [] # zeros((max(atms),3),dtype=float)
    fi = open(w2k.case+'.outputdmfu','r')
    for line in fi:
        m = re.search('Actual position of atom *(\d+)',line)
        if m is not None:
            if int(m.group(1)) in atms:
                atom = int(m.group(1))

                print 'ATOM', atom, ' in atms=', atms
                
                m = re.search('Actual position of atom *(\d+) *is: (.*)',line)
                #pos[atom-1,:] = map(float,m.group(2).split())
                pos.append( map(float,m.group(2).split()) )
    pos=array(pos)
    print 'Actual position of atoms for downfolding=',  pos.tolist()
    
    strc = Struct(w2k.case)
    strc.parse()
    
    # reads BR1 and BR2 from case.rotlm
    (BR,BRX) = ReadBRs(w2k.case)
    # Computes direct vectors, knowing the reciprocal basis
    DV = DirectVectors(BR[1],BR[0])
    
    print 'Direct vectors:'
    for i in range(3):
        print ("%7.3f "*3) % tuple(DV[i])

    print 'Please enter the direct vectors for large Brillouin zone'
    DVn = zeros(1)
    while shape(DVn) != (3,3):
        DVn = array(eval(sys.stdin.readline()))
    #DVn = array([[0.5,0.5,0],[0.5,-0.5,0],[0,0,1]])
    print 'But switching to direct vectors of the one atom unit cell:'
    for i in range(3):
        print ("%7.3f "*3) % tuple(DVn[i])
    
    # Reads div from case.klist
    firstline = open(w2k.case+'.klist').next()
    m = re.search('div: *\((\s*\d+\s*\d+\s*\d+)\)', firstline)
    if m is not None:
        div = map(int,m.groups()[0].split())

    Rs = GiveRs(div, DVn)
    print 'Real space vectors used in the Fourier transform to tigh-binding:'
    for i,R in enumerate(Rs):
        print "%3d"%i, "%7.3f %7.3f %7.3f" % tuple(R)

    # Read Udmft
    filename = 'Udmft.0'
    (nkp, nsymop, norbitals) = rdU.fopen(fhp, filename)
    
    nindo = rdU.read1(fhp, norbitals)

    print 'nkp=', nkp, 'nsymop=', nsymop, 'norbitals=', norbitals, 'nindo=', nindo


    delta_A_B = pos[1]-pos[0]
    print 'delta_A_B=', delta_A_B
    
    rdU.feopen(fhe, w2k.case+'.energy', strc.nat)
    
    Hk=[]
    kps=[]
    wgh=[]
    for ikp in range(nkp):
    
        (iikp, nbands, tmaxdim2, tnorbitals, nemin) = rdU.read2(fhp)
        (k, kname, wg, ios, nen) = rdU.feread1(fhe)

        print 'delta*k=', delta_A_B*k
                
        if ios!=0:
            print 'Number of k-points in case.energy is smaller then in Udfmt.0'
            break
        
        # Read LDA energies from case.energy
        Ek = rdU.feread2(fhe, nen) * Ry_in_eV - EF
        
        # finds minimum and maximum band for downfolding
        #(bnd_s, bnd_e) = FindMaxMinBand(Ek, options.Emin, nemin, options.nbnd, options.dnbnd)
        #e = Ek[nemin+bnd_s-1:nemin+bnd_e-1]
        #print 'bnd_s=', bnd_s, 'bnd_e=', bnd_e, 'Emm=', e[0], e[-1]
        #print 'k=', k
        
        for isym in range(nsymop):
            
            kn = dot(strc.iz[isym],k)
            k_nat = [sum(BRX[i,:]*kn) for i in range(3)]
            
            kps.append(kn)
            wgh.append(wg)
            
            iisym = rdU.read3(fhp)
            if iisym!=isym+1:
                print 'ERROR: isym!=iisym', isym, iisym
    
            for iorb in range(norbitals):
    
                U = []
                for ind in range(nindo[iorb]):
                    Udat = rdU.read4(fhp, nbands)
                    U.append( Udat )
                U = array(U)
                
                if iorb==0:
                    Utot = U
                else:
                    Utot = vstack((Utot,U))

            # finds minimum and maximum band for downfolding
            Nse = FindMaxMinBand(Ek, options.Emin, options.Emax)

            print 'len(Utot)=', len(Utot)

            while (Nse[1]-Nse[0] < len(Utot) ):
                print 'Warning: Number of bands in the interval is ', Nse[1]-Nse[0], 'while number of orbitals is ', len(Utot)
                (Nse, Succ) = TryAddBands(Nse, Ek, options.Emin, options.Emax, options.SmearN, options.SmearX)
                if not Succ: sys.exit(1)

            NsingularValues=0
            while (NsingularValues<len(Utot)):
                (bnd_s,bnd_e) = Nse-(nemin-1)
                e = Ek[nemin-1+bnd_s:nemin-1+bnd_e]
                print 'bnd_s=', bnd_s, 'bnd_e=', bnd_e, 'n=', bnd_e-bnd_s, 'Emm=', e[0], e[-1]            
                # Performing SVD on projector
                (u, s, v) = linalg.svd(Utot[:,bnd_s:bnd_e])
                print 'singular-values:', ikp, s.tolist()
                NsingularValues = len(filter(lambda x: x>options.smalls, s))
                print 'NsingularValues=',  NsingularValues

                if (NsingularValues<len(Utot)):
                    (Nse, Succ) = TryAddBands(Nse, Ek, options.Emin, options.Emax, options.SmearN, options.SmearX)
                    if not Succ: sys.exit(1)
                    
                    
            # downfolding transformation
            Rt = dot(u[:,:len(s)],v[:len(s),:])
            # downfolded hamiltonian
            tHk = dot(conjugate(Rt) * e, transpose(Rt))

            N0 = len(Utot)/2
            print 'N0=', N0
            
            zHk = copy.deepcopy(tHk[:N0,:N0])
            zHk += tHk[:N0,:N0]
            zHk += tHk[:N0,N0:]*exp(2*pi*1j*sum(delta_A_B*kn) ) # 1,2
            zHk += tHk[N0:,:N0]*exp(-2*pi*1j*sum(delta_A_B*kn) )# 2,1
            zHk *= 0.5
            Hk.append(zHk)

    
    rdU.fclose(fhp)
    rdU.feclose(fhe)
                
    Hk = array(Hk)     
    
    nrm = sum(wgh)
    wgh = wgh/nrm
    
    # Fourier transform to produce hopping parameters
    Hr=[]
    for R in Rs:
        tHr=zeros(shape(Hk[0]),dtype=complex)
        dsum=0j
        for ik in range(len(kps)):
            expk = exp(1j*2*pi*dot(kps[ik],R)) * wgh[ik]
            tHr += expk * Hk[ik]
            dsum += expk
        print 'R=[', "%4.1f, "*3 % tuple(R), '] exp=', "%10.5f" % abs(dsum)
        Hr.append(tHr)
    

    # Printing of the Hamiltonian
    fh = open('hamiltonian.dat', 'w')
    
    print >> fh, '# Tighbinding Hamiltonian with entries: ', inl.legends
    print >> fh
    #print >> fh, "lattice_unit=[", tuple(DVn[0]), ',', tuple(DVn[1]), ',', tuple(DVn[2]), "]"
    print >> fh, "lattice_unit=", identity(3).tolist()

    print >> fh, "orbital_positions=",zeros((3,3)).tolist()
    print >> fh
    #print '['
    #for atm in inl.siginds.keys():
    #    for i in range(len(inl.siginds[atm])):
    #        print >> fh, pos[atm-1].tolist(),','
    #print 
    #print >> fh, ']\n'


    # sorting such that first come larger hoppings
    for irm,R in enumerate(Rs):
        if sum(abs(Rs[irm]))<1e-10: break
    
    ir_index = range(len(Rs))
    ir_index.sort( lambda a,b: cmp(sum(abs(Hr[b])), sum(abs(Hr[a]))) )
    ir0 = ir_index.index(irm)
    ir_index = [irm] + ir_index[:ir0] + ir_index[ir0+1:]
    
    #print ir_index
    
    DVn_inv = linalg.inv(DVn)
    print >> fh, 'Hopping={}'
    for ir in ir_index :
        # RS*DV^-1 = (i0,i1,i2)
        Rn = dot(Rs[ir] ,DVn_inv)
        Print(fh, Rn, Hr[ir])
        print >> fh
        
    fh.close()
