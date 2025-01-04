import glob
import os, sys
from itertools import chain
import numpy as np

def get_case():
    # directory name is not case (happens when submitting to cluster)
    case = None
    files = glob.glob('*.struct')
    if len(files) < 1:
        raise Exception('No struct file present.')
    elif len(files) > 1:
        # heuristic algorithm to determine case:
        # need in0 and struct present with the same case.
        candidates = [os.path.splitext(f)[0] for f in files]
        for cand in candidates:
            if os.path.isfile(cand+'.in0'):
                case = cand
                break
    else: # just one candidate exists, hence it must be it
        case, ext = os.path.splitext(os.path.basename(files[0]))
    return case

class Indmfl:
    '''Class for case.indmfl file. Simplified/stripped down implementation for plotting.
    '''
    def __init__(self, case, extn='indmfl'):
        self.case = case
        self.extn = extn
        self.broken_sym = 0 # it seems we implemented this, but is not used
        # Finding the chemical potential
        scf2_exists = os.path.isfile(case+".scf2")
        scf2up_exists = os.path.isfile(case+".scf2up")
        scf_exists = os.path.isfile(case+".scf")
        self.EF = None
        
        if  os.path.isfile('EF.dat'):
            # The previous DMFT chemical potential
            self.EF = float( open('EF.dat','r').readline() )
            
        if self.EF is None:
            scf2_exists = os.path.isfile(case+".scf2")
            scf2up_exists = os.path.isfile(case+".scf2up")
            if scf2_exists or scf2up_exists:
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
            for line in lines[::-1]:
                if re.match(r':FER', line) is not None:
                    Ry2eV = 13.60569193
                    self.EF = float(line[38:])*Ry2eV
                    break
            
        if self.EF is None: self.EF = 0
        self.cixgrp={}

    def divmodulo(self, x,n):
        "We want to take modulo and divide in fortran way, so that it is compatible with fortran code"
        return ( np.sign(x)* int(abs(x)/n) , np.sign(x)*np.mod(abs(x),n))

    def read(self, filename = None):
        if filename==None:
            filename = self.case+'.'+self.extn
        with open(filename, 'r') as findmf:
            lines = [line.split('#')[0].strip() for line in findmf.readlines()] # strip comments
        lines_strings = [line for line in lines if line]  # strip blank lines & create generator expression
        lines = iter(lines_strings)
        #
        #self.parse_head(lines)
        self.hybr_emin, self.hybr_emax, self.Qrenorm, self.projector = [float(x) for x in next(lines).split()]
        self.matsubara, self.broadc, self.broadnc, self.om_npts, self.om_emin, self.om_emax = [float(e) for e in next(lines).split()]
        self.matsubara = int(self.matsubara)  # recast these to integers
        self.om_npts   = int(self.om_npts) 
        self.Qrenorm   = int(self.Qrenorm)
        self.projector = int(self.projector)
        if self.projector>=5:
            self.hybr_emin = int(self.hybr_emin)
            self.hybr_emax = int(self.hybr_emax)
        #
        self.emin, self.emax = self.hybr_emin, self.hybr_emax
        # correcting self.hybr_emin with values in case.indmf, if exists.
        if self.hybr_emin==0 and self.hybr_emax==0:
            if os.path.exists(filename[:-1]):
                dat = open(filename[:-1],'r').readline().split()
                self.hybr_emin, self.hybr_emax = float(dat[0]), float(dat[1])
            else:
                self.hybr_emin, self.hybr_emax = -10,10
        #self.parse_atomlist(lines)
        self.Lsa=[]
        self.icpsa=[]
        self.locrot={}
        self.atoms={}
        self.cix={}
        natom = int(next(lines))
        for i in range(natom):
            dat=next(lines).split()
            iatom, nL, locrot_shift = [int(x) for x in dat[:3]]
            Rmt2=0
            if len(dat)>3:
                Rmt2 = float(dat[3])
            (shift,locrot) = self.divmodulo(locrot_shift,3)
            if locrot<0:
                if locrot==-2: 
                    locrot=3*nL
                else:
                    locrot=3
            (Ls, qsplits, icps) = (np.zeros(nL,dtype=int), np.zeros(nL,dtype=int), np.zeros(nL,dtype=int))
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
        
        self.legends={}
        self.siginds={}
        self.cftrans={}
        # read the big block of siginds and cftrans
        ncp, maxdim, maxsize = [int(e) for e in next(lines).split()]
        for i in range(ncp):
            icp, dim, size = [int(e) for e in next(lines).split()]
            self.legends[icp] = next(lines).split("'")[1::2]
            self.siginds[icp] = np.array([[int(e) for e in next(lines).split()] for row in range(dim)])
            raw_cftrans = np.array([[float(e) for e in next(lines).split()] for row in range(dim)])
            self.cftrans[icp] = raw_cftrans[:,0::2] + raw_cftrans[:,1::2]*1j

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
        iSiginds[icix] = np.array(imp,dtype=int)
    return iSiginds

class Struct:
    """Similar to wstruct.py implementation, but stripped down and simplified for plotting.
    """
    def __init__(self, inAngs=True):
        self.tobohr = 1/0.5291772083
        hex2ort = np.array([[0,1,0],[np.sqrt(3.)/2.,-0.5,0],[0,0,1]])
        ort2rho = np.array([[1/np.sqrt(3),1/np.sqrt(3),-2/np.sqrt(3)],[-1,1,0],[1,1,1]])
        hex2rho = hex2ort @ ort2rho
        self.rho2hex = np.linalg.inv(hex2rho)
        self.HRtransf = False
        self.inAngs = inAngs
        
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
            pos = [float(line[col:col+10]) for col in 12+13*np.arange(3)]
            self.pos.append(pos)
            
            line = next(f)
            _mult_ = int(line[15:17])
            self.mult.append(_mult_)
            self.isplit.append( int(line[34:37]) )
        
            for mu in range(_mult_-1):
                line = next(f)
                pos = [float(line[col:col+10]) for col in 12+13*np.arange(3)]
                self.pos.append(pos)
                
            line = next(f)
            self.aname.append( line[0:10].strip() )
            self.jrj.append( int(line[15:20]) )
            self.r0.append( float(line[25:35]) )
            self.rmt.append( float(line[40:50]) )
            self.Znuc.append( int(float(line[55:65])) )
            
            rt=[]
            for i in range(3):
                line = next(f)
                rt.append( [float(line[col:col+10]) for col in 20+10*np.arange(3)] )
            self.rotloc.append(np.array(rt).T)

        self.pos = np.array(self.pos)

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
        self.timat  = np.zeros((nsym,3,3),dtype=int) # it is transpose compared to w2k fortran convention
        self.tau    = np.zeros((nsym,3))
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

def SimplifySiginds(siginds):
    """Takes dictionary of Sigind's and creates list of non-zero values
    that appear in Siginds, cols[icix]=[1,2,3,...]
    """
    cols={}
    for icix in siginds:
        # should remove all negative indices, because they are treated as static self-energy!!!
        cols[icix] = [x for x in set(siginds[icix].ravel()) if x>0] 
    return cols

def ImpurityLatticeConnection( cols, icols, log ):
    """We have correlated index stored in Sigind[icix] and
    impurity index stored in iSigind[iimp].
    The translation between the two is necessary when
    the impurity would like to know in which reference frame
    is impurity problem defined. Local rotation should hence
    be passed to impurity, and for that we need the translator.
    """
    imp2latt={}
    for i in icols: # all impurity problems
        imp2latt[i]=[]
        for j in cols: # all correlated atoms
            imp_cols = set(icols[i])
            atm_cols = set(cols[j])
            dif = imp_cols-atm_cols # checks if the same indices
            imp_contains_atm = (len(imp_cols)==len(atm_cols) and len(dif)==0) or (len(imp_cols)>len(atm_cols) and len(dif)==len(imp_cols)-len(atm_cols))
            #print('i=', i, 'j=', j, 'dif=', dif, 'imp_contains_atm=', imp_contains_atm)
            if imp_contains_atm and len(icols[i])>0:
                imp2latt[i].append(j)
        if len(imp2latt)==0:
            print('WARNING: in ImpurityLatticeConnection impurity['+str(i)+'] has no connection with correlated atom',file=log)
    all_atoms_in_impurities = sorted(chain(*imp2latt.values()))
    all_atoms = sorted([j for j in cols])
    if all_atoms_in_impurities != all_atoms:
        print('ERROR: in ImpurityLatticeConnection all_atoms=',str(all_atoms), file=log)
        print(' while atoms found in impurities', all_atoms_in_impurities, 'are different', file=log)
        sys.exit(1)
    return imp2latt

def FindEF():
    EF = None
    if  os.path.isfile('EF.dat'):
        # The previous DMFT chemical potential, if not fixEF
        EF = float( open('EF.dat','r').readline() )
    
    if EF is None:
        scf2_exists = os.path.isfile(case+".scf2")
        scf2up_exists = os.path.isfile(case+".scf2up")
        if scf2_exists or scf2up_exists:
            fname = case+".scf2" if scf2_exists else case+".scf2up"
            fscf = open(fname, 'r')
            lines = fscf.readlines()
            for line in lines:
                if re.match(r':FER', line) is not None:
                    Ry2eV = 13.60569193
                    EF = float(line[38:])*Ry2eV
                    break
    
    if EF is None and scf_exists:
        fname = case+".scf"
        fscf = open(fname, 'r')
        lines = fscf.readlines()
        for line in lines[::-1]:
            if re.match(r':FER', line) is not None:
                Ry2eV = 13.60569193
                EF = float(line[38:])*Ry2eV
                break
    return EF

class InteractiveLegend:
    def __init__(self,lines,axs,fig,fontsize='medium'):
        self.map_legend_to_ax = {}  # Will map legend lines to original lines.
        pickradius = 5  # Points (Pt). How close the click needs to be to trigger an event.
        self.leg = axs.legend(loc='best',fontsize=fontsize)
        self.leg.set_draggable(True)
        for legend_line, ax_line in zip(self.leg.get_lines(), lines):
            legend_line.set_picker(pickradius)  # Enable picking on the legend line.
            self.map_legend_to_ax[legend_line] = ax_line
        # Works even if the legend is draggable. This is independent from picking legend lines.
        self.fig = fig
        
    def __call__(self, event):
        # On the pick event, find the original line corresponding to the legend
        # proxy line, and toggle its visibility.
        legend_line = event.artist
    
        # Do nothing if the source of the event is not a legend line.
        if legend_line not in self.map_legend_to_ax:
            return
        
        ax_line = self.map_legend_to_ax[legend_line]
        visible = not ax_line.get_visible()
        ax_line.set_visible(visible)
        # Change the alpha on the line in the legend, so we can see what lines
        # have been toggled.
        legend_line.set_alpha(1.0 if visible else 0.2)
        self.fig.canvas.draw()

if __name__ == '__main__':
    
    case = get_case()
    inl = Indmfl(case)
    inl.read()
    mu = FindEF()
    print(inl)
    print('mu=', mu)

    iSiginds = ParsIndmfi(case)
    _icols = SimplifySiginds(iSiginds)
    _lcols = SimplifySiginds(inl.siginds)
    imp2latt = ImpurityLatticeConnection(_lcols, _icols, sys.stdout)
    # these are columns we can fill in with our impurity problems
    print('icols=', _icols)
    print('lcols=', _lcols)
    icols={}
    for ii in _icols:
        icols[ii+1] = list(range(len(_icols[ii])))
    lcols={}
    for icix in _lcols:
        lcols[icix] = (np.array(_lcols[icix])-min(_lcols[icix])).tolist()
    print('imp2lattice connection=', imp2latt)
    for ii in imp2latt:
        print(' imp.'+str(ii)+'/Gf.out icols['+str(ii+1)+']=', icols[ii+1])
        print(' connected to ')
        for icix in imp2latt[ii]:
            print('  gc'+str(icix)+'   lcols['+str(icix)+']=',np.array(lcols[icix]))
    print('All available columns:')
    print(' icols=', icols)
    print(' lcols=', lcols)
    
