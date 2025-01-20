#!/usr/bin/env python
# @Copyright 2007 Kristjan Haule
# 
import sys, re, os, glob, shutil, types
from numpy import *
from scipy import interpolate
from utils import DmftEnvironment
from cmpEimp2 import GiveImpFunctional
import brod  # should be replaced by system call
import subprocess
import dmft_DC
import builtins

Default_CoulombF = 'Ising'
mom=0       # We need those variables to be global in order to use python exect
TrSigmaG=0  # We need those variables to be global in order to use python exect
            
class IMP_CTQMC:
    """Wrapper for ctqmc
    """
    def __init__(self, Znuc, dire, params, Sigind, CF):
        """
         Prepares several variables, which will be passed to ED for impurity solver, and eventualy to ctqmc impurity solver.

         self.sparams --> will be given to ctqmc impurity solver
         self.spr     --> will be given to ED of the atom, which will produce actqmc.cix, and will be eventually used by ctqmc solver.

         atom_d.py (which is ED for l<3) requires:
            Eimp=[], l,J, HB2==(True==mode=S), Ewindow=[-1000,1000], max_M_size=500, CoulombF
                         ORB=[] -- index for orbital susceptibility
                         add_occupancy -- add occupancy to cix file
                         HB2   <-> self.spr['mode']==S
                         CF is given through imp.x/Trans.dat
                         
         atom_d.py does not need anymore:
            qOCA, mOCA, n, Eoca, OCA_G, 

         atom_h (which is ED for l==3) requires ('note that we assume CF is j-j basis'):
            Eimp=[], l,J, HB2==(mode=S), cx=SOC
            ns,ne     <-> self.spr['n']=[ns,ne] -- starting and final occupancy to include
            Ekeep=100                           --  Energy window in eV treated exactly
            Ekeepc=100                          --  Energy window in eV treated exactly in core
            Nmax=1000                           --  Maximum number of states kept in superstate
            Nmaxc=100                           --  Maximum number of states kept in superstate for core
            Impq=[]                             --  List index for equivalent orbitals
            pm                                  --  Prints moment instead of Sz
            ImpTot                              -- Take full Eimp and not just splitting of Eimp - do not subtract Eimp[0]
            HB2       <-> self.spr['mode']==S   -- If we want to compute self-energy from Swinger-like equation and two particle response
        
         The variables for ED will be extracted from params. But not all are relevant for ED, hence
         we will do the following:
           if the variable is in vars and is defined in params, it takes its value from params
           even if the variable is not in vars, but is defined in params (user-defined) as atom_xxx, it is added to spr[xxx]
           This is to allow user to change any varibale for ED through atom_xxx.
        """
        # Environment
        self.env = DmftEnvironment()
        self.mpi_prefix = self.env.MPI
        
        # impurity directory
        self.dir = dire+'/'
        
        # create impurity directory if it does not yet exist
        #if len(glob.glob(dire))==0 : os.mkdir( dire )
        if not os.path.isdir(dire):
            os.mkdir( dire )

        self.Sigind = array(Sigind)
        self.CF = CF
            
        # log-file is opened
        self.fh_info = open(self.dir + 'info.out', 'w')

        # Note that run_dmft.py adds params['icase']<-indmfl.icase
        # and params['UlamJ'] (the latter only for exact DC)
        # to the params[iparams{xx}] from input file
        # 
        # CoulombF is needed later on, so we add one if not existant
        if 'CoulombF' not in params:
            params['CoulombF']=["'"+Default_CoulombF+"'", "# Full Coulomb repulsion"]
        # now all parameters in params are being remembered
        for p in params:
            if type(params[p][0]) is str:
                setattr(self, p, params[p][0])
            else:
                setattr(self, p, params[p][0])
            print('setting ctqmc.self.'+p+'=', params[p][0], file=self.fh_info)
        
        
        #if 'CoulombF' in params:
        #    CoulombF = params['CoulombF'][0].strip("'").strip('"')
        #else:
        #    params['CoulombF']=["'"+Default_CoulombF+"'", "# Full Coulomb repulsion"]
        #    CoulombF = Default_CoulombF
        #
        #if 'icase' in params:
        #    self.icase = params['icase']
        #if 'UlamJ' in params:
        #    self.UlamJ = params['UlamJ']
        
        
        # find smallest sigind to shift siginds so they start by 1
        minsigind = min([s for s in self.Sigind.flat if s != 0])
        self.Sigind_shifted = where(self.Sigind != 0, self.Sigind-minsigind+1, 0)
        
        # Transformation between spheric harmonics and local coordinates needs to be
        # used in exact diagonalization of the atom.
        with open(self.dir + 'Trans.dat', 'w') as fh_T:
            print(len(Sigind), len(Sigind[0]), '#  size of Sigind and CF', file=fh_T)
            print('#---- Sigind follows', file=fh_T)
            max_sigfig = 1 + int(log(max(self.Sigind_shifted.flat))/log(10))
            fmt = '{:'+str(max_sigfig)+'d} '
            for row in self.Sigind_shifted:
                print((fmt*len(row)).format(*row), file=fh_T)
            print('#---- CF follows', file=fh_T)
            for row in CF:
                for elem in row:
                    print('{:12.8f} {:12.8f} '.format(elem.real,elem.imag), end=' ', file=fh_T)
                print(file=fh_T)
        
        # If old distribution of kinks exists, it is copied to imp directory
        status = glob.glob("status.*")
        if self.dir != './':
            print('Copying status files from to', self.dir, file=self.fh_info)
            for fil in status:
                shutil.copy2(fil, self.dir+'/'+fil)
        
        self.PARAMS = 'PARAMS'
        
        # parameters, which will be given to ctqmc and will be written to PARAMS for ctqmc execution
        # This is just a default set of paramerets, which will be modified before the execution
        self.sparams={'exe':'ctqmc',
                      'Sig':'Sig.out',
                      'Delta':'Delta.inp',
                      'cix':'actqmc.cix',
                      'nom': 100 if 'beta' not in params else int(3*params['beta'][0]),
                      'mode':'GH',
                      'aom':1,
                      'Ncout':1000000,
                      'GlobalFlip':1000000,
                      'Gf':'Gf.out',
                      'OffDiagonal':'real'}
        print('At ctqmc.__init__ self.sparams=', self.sparams, file=self.fh_info)
        
        vars_for_ED = ['l','J','mode',                                              # relevant for both EDs
                       'Ewindow', 'max_M_size', 'CoulombF', 'add_occupancy',        # relevant for atom_d.py (l=0,1,2)
                       'cx', 'n', 'Ekeep', 'Ekeepc', 'Nmaxc',                       # relevant for atomh (l=3)
                       'para', 'qOCA', 'kOCA', 'mOCA', 'Ex', 'Ep', 'Eoca', 'qatom'] # probably not relevant anymore, it could be removed probably
        
        self.spr={} # These variables will be given to ED.
        # We check for certain names in vars_for_ED above, and if variable appears in list, it is extracted from params
        # and copied to self.spr.
        # Alternatively user can give any variable to ED algorithm by prepanding name 'atom_xxxx', and variable 'xxxx' will
        # be copied to self.spr.
        for var in vars_for_ED:
            if var in params:
                self.spr[var] = params[var][0]
        for var in params:
            if var[:5]=='atom_':
                self.spr[var[5:]] = params[var][0]
        
        # Some variables have different name in database and input file
        # spr['n'] -- params['nc']
        if 'nc' in params:
            self.spr['n'] = params['nc'][0]
        
        # stores some of just defined variables
        if 'l' in self.spr:
            self.l = self.spr['l']
        else:
            if (Znuc>=57 and Znuc<=71) or (Znuc>=89 and Znuc<=103): self.l=3  # lanthanides & actinied
            elif Znuc in [1,2,3,4,11,12,19,20,37,38,55,56,87,88]: self.l=0
            elif Znuc in [5,6,7,8,9,13,14,15,16,17,31,32,33,34,35,49,50,51,52,53,81,82,83,84,85]: self.l=1
            else: self.l=2
            self.spr['l'] = self.l
            
        # Correcting Hund's coupling if exact double-counting for consistency.
        # UlamJ= {0: (U, lmbda, epsilon, [J2, J4])}
        if 'UlamJ' in params and self.UlamJ[1]>1e-5:  # DCs='exact'
            # If we use Yukawa screening, we need to change Jhunds to be compatible with Yukawa form of screening
            # For dielectric screening we do not change Jhunds from user specified values
            Jh = self.UlamJ[3]
            if self.spr['l']==3: Jh = sum(Jh)/len(Jh)  # Currently does not handle an array of J's. You should correct that.
            self.spr['J']=Jh
            self.J = Jh

        if self.l==3 and ('n' not in self.spr):
            self.spr['n']= list(range(builtins.max(0,params['nf0'][0]-2),builtins.min(params['nf0'][0]+3,2*(2*self.l+1))))
        
        # extra off-diagonal spin-orbit
        self.cx = 0 if 'cx' not in self.spr else self.spr['cx']
        print('At ctqmc.__init__ self.spr=', self.spr, file=self.fh_info)
        print('At ctqmc.__init__ self.l=', self.l, 'self.cx=', self.cx, 'self.J=', self.J, file=self.fh_info)
    

    def _exe(self, params, DCs, extn, UpdateAtom, gbroad=0.0, kbroad=0.0, maxAc=2000.):
        """ Executes the CTQMC impurity solver
        """
        REMOVE_STATUS=False
        Qaverage = False
        # Reading impurity levels
        (Eimp,Olap,Edc) = loadtxt(self.dir+'/Eimp.inp',ndmin=2)   # impurity levels should be produced in run_dmft.py at imp.x/Eimp.inp
        # Reading Delta
        Delta_dat = loadtxt(self.dir+'/Delta.imp').T              # Delta should be at imp.x/Delta.imp
        om = Delta_dat[0]
        Delta = Delta_dat[1::2] + Delta_dat[2::2]*1j
        # subtracting Edc, because the ksum used (Hk+s_oo-Edc)
        self.Eimps = Eimp-Edc                                     # impurity levels with double-counting included, as needed for ED & ctqmc

        if DCs=='fixn0':  # There is a rarely used option (mostly for fs) that do not allow small crystal field splitting in impurity levels.
            self.Eimps = ones(len(self.Eimps))*sum(self.Eimps)/len(self.Eimps)  # All impurity levels should be the same
            print('Eimps(fixn0)=', self.Eimps, file=self.fh_info)
        
        # Adding impurity levels to dictionary params, which is input to this routine
        params['Ed'] = [self.Eimps.tolist(), "# Impurity levels"]
        
        
        # Our choice for mu_QMC is somewhat arbitrary, because we will subtract it from Eimps below before use.
        # The choice for l=3 is:
        #   impurity level for 5/2 is E_{5/2}=E0-2*cx and for 7/2 is E_{7/2}=E0+3/2*cx
        #   Here we correct E_{5/2} for 2*cx, which is already counted in cix file
        #E0 = params['Ed'][0][0]+2*self.cx
        E0 = self.Eimps[0]
        if (self.l == 3):
            if (Qaverage):
                # 5/2 energy levels:
                Eimps_5_2 = [self.Eimps[self.Sigind[i,i]-1] for i in range(2*self.l)]
                aEimps_5_2 = sum(Eimps_5_2)/(2*self.l)
                print('<E_5/2>=', aEimps_5_2)
                #E0 = Eimps[0]+2*self.cx
                E0 = aEimps_5_2 + 2*self.cx
            else:
                E0 = self.Eimps[0]+2*self.cx
        mu_QMC = -E0
        
        # Should be 1 for all impurity levels used. Some are not used, hence are just removed.
        Ident = zeros(len(self.Eimps))
        for i in range(len(self.Sigind)):
            if self.Sigind_shifted[i,i]>0: Ident[self.Sigind_shifted[i,i]-1]=1.0
        print('Ident=', Ident, file=self.fh_info)
        
        # Below exact diagonalization of the atom is run to produce 'actqmc.cix' file
        if UpdateAtom or not os.path.isfile(self.dir+'/'+self.sparams['cix']) : # no cix file yet -> create it
            Eimps = (array(self.Eimps)-E0*Ident).tolist()
            print('Eimps_new=', Eimps, ' muQMC=', mu_QMC, file=self.fh_info)
            Impq = [self.Sigind_shifted[i,i]-1 for i in range(len(self.Sigind_shifted))]
            print('Impq=', Impq, file=self.fh_info)
            
            # Variables which are updated from the current paramateres
            for var in ['CoulombF','mode']:
                if var in params:
                    self.spr[var] = params[var][0]

            if (self.l == 3) :
                # executable for ED of the atom
                self.atom_exe = self.env.ROOT + '/atomh'
                # argument list is created to later run the executable for exact diagonalization of the atom
                atom_arg = '-ns %d -ne %d ' % (self.spr['n'][0], self.spr['n'][-1])
                if (not Qaverage):
                    atom_arg += ' -Eimp \"'+str(Eimps)+'\" '
                    atom_arg += ' -Impq \"'+str(Impq)+'\" '
                    
                for k,p in self.spr.items():
                    if k=='mode':
                        if p[0]=='S':
                            atom_arg += '-HB2 1 '
                        continue
                    if type(p)==list:
                        atom_arg += '-'+k+' "'+str(p)+'" '
                    else:
                        atom_arg += '-'+k+' '+str(p)+' '
                # print to info file
                print('atomh arguments=', atom_arg, file=self.fh_info)
                # runs ED for the atom
                cmd = 'cd '+ self.dir + '; ' + self.atom_exe +' ' +atom_arg+' > nohup.out 2>&1'
                ret = subprocess.call(cmd,shell=True,stdout=self.fh_info,stderr=self.fh_info)
                print('ret=', ret, file=self.fh_info)
            elif (self.l < 3) :
                # executable for ED of the atom
                self.atom_exe = self.env.ROOT + '/atom_d.py'

                if 'U' in params and 'Ud' in params and params['U'][0]==0.0 and params['U'][0]!=params['Ud'][0]:
                    # Must be cluster-DMFT case. We will assume it is 2-site cluster case (link)
                    self.FromSingleImpurity_to_Link(Eimps, params)
                else:
                    # Old good single impurity calculation
                    atom_arg = ''
                    for k,p in list(self.spr.items()):
                        if k=='mode':
                            if p[0]=='S':
                                atom_arg += '\"HB2=True\" '
                        else:
                            atom_arg += '\"'+k+'='+str(p)+'\" '
                    atom_arg += ' \"Eimp='+str(Eimps)+'\" '
                    
                    print('atom_d.py arguments=', atom_arg, file=self.fh_info)
                    # runs ED for the atom
                    cmd = 'cd '+ self.dir + '; ' + self.atom_exe +' ' +atom_arg+' > nohup.out 2>&1'
                    print('.... running ..', cmd, file=self.fh_info)
                    ret = subprocess.call(cmd,shell=True,stdout=self.fh_info,stderr=self.fh_info)
                    print('ret=', ret, file=self.fh_info)
                # If new cix file is created, the eigenvectors are reordered and hence
                # previous ctqmc kinks could be incompatible. One would need to reorder
                # them. This should be relevant for crystal-field switching.
                if (REMOVE_STATUS):
                    # I should work on this, because eigenvectors are similar but not identical.
                    # checking \psi^+_alpha|0> and \psi^+_alpha \psi^+_\beta ...|0> should be obvious
                    for filename in glob.glob(self.dir+'/status.*') :
                        os.remove( filename )
            else:
                print('l=', self.l, ' not yet implemented!', file=self.fh_info)
        
        
        self.fh_info.flush()
        
        if self.l==3:
            params['exe']=['ctqmcf', '#']

        params['mu'] = [mu_QMC, '# QMC chemical potential']
        # executable
        self.q_exe = self.env.ROOT+'/'+params['exe'][0]
        
        # Preparing params.dat file
        with open(self.dir+self.PARAMS, 'w') as fp:
            for p in self.sparams:
                if len(p)>5 and p[:5]=='atom_': continue # no atomic parameters needed here
                if p not in params:  # if this parameter is overwritten by the input parameters, skip it for now
                    print('{:s}  {:s}'.format(p, str(self.sparams[p])), file=fp)
            for p in params:
                if len(p)>5 and p[:5]=='atom_': continue
                print('{:10s}  {:10s}  {:s}'.format(p, str(params[p][0]), params[p][1]), file=fp)

        Sigind = array(self.Sigind_shifted)
        Diagonal={}
        for i in range(len(Sigind)):
            for j in range(len(Sigind)):
                if Sigind[i,j]>0: Diagonal[Sigind[i,j]-1]= (i==j)

        print('Diagonal=', Diagonal, file=self.fh_info)
        

        # Preparing Delta.inp file
        with open(self.dir+self.sparams['Delta'], 'w') as fp:
            # Checking that Delta is causal and not regular
            for ib in range(shape(Delta)[0]):
                if Diagonal[ib]:
                    for im in range(shape(Delta)[1]):
                        if abs(Delta[ib,im])>maxAc*pi : Delta[ib,im] = -1e-10j           # just wrong
                        if Delta[ib,im].imag>0 : Delta[ib,im] = Delta[ib,im].real-1e-10j # non-causal
            
            # creating big mesh
            T=1./params['beta'][0]
            m_max = int(om[-1]/(pi*T)+1e-6) # n of the last mastubara point
            lom = arange(1,m_max+1,2)*pi*T  # all Matsubata points included
            # interpolate on big mesh
            lDelta=[]
            
            print('shape(Delta)=', shape(Delta), file=self.fh_info)
            print('shape(om)=', shape(om), file=self.fh_info)
            print('shape(lom)=', shape(lom), file=self.fh_info)
            print('ss=', shape(Delta[0].real), file=self.fh_info)
            print('om[0]=', om[0], 'lom[0]=', lom[0], 'om[-1]=', om[-1], 'lom[-1]=', lom[-1], file=self.fh_info)
            
            for ib in range(len(Delta)):
                fDr = interpolate.UnivariateSpline(om, Delta[ib].real, s=0)
                fDi = interpolate.UnivariateSpline(om, Delta[ib].imag, s=0)
                lDelta.append( fDr(lom) + fDi(lom)*1j )
                
            for im,ome in enumerate(lom):
                print(ome, end=' ', file=fp)
                for b in range(len(lDelta)):
                    print("{:19.12f}  {:19.12f} ".format(lDelta[b][im].real, lDelta[b][im].imag), end=' ', file=fp)
                print(file=fp)
        
        shutil.copy2(self.dir+self.sparams['Delta'], self.dir+self.sparams['Delta']+'.'+extn)
        
        # saving three previous steps of status files, in case we need to restart from previous steps.
        if os.path.exists(self.dir+'status_2.tgz'):
            shutil.move(self.dir+'status_2.tgz', self.dir+'status_3.tgz')
        if os.path.exists(self.dir+'status_1.tgz'):
            shutil.move(self.dir+'status_1.tgz', self.dir+'status_2.tgz')
        if os.path.exists(self.dir+'status.000'):
            cmd = 'cd '+self.dir+'; tar czvf status_1.tgz status.*'
            ret = subprocess.call(cmd,shell=True,stderr=self.fh_info,stdout=self.fh_info)
        

        # Below we execute ctqmc
        cmd = 'cd '+self.dir+'; '+self.mpi_prefix+' '+self.q_exe+' '+self.PARAMS+' > nohup_imp.out 2>&1 '
        print(cmd, file=self.fh_info)
        print('Running ---- ctqmc -----', file=self.fh_info)
        subprocess.call(cmd,shell=True,stdout=self.fh_info,stderr=self.fh_info)

        
    def HighFrequency(self, params, DCs, nl_imp, extn, wbroad=0.0, kbroad=0.0):
        """ Reads impurity output and prepares for further execution
          Output:
               Sig[b][iom]  -- DMFT dynamic self-energy which vanishes at infinity
               sinf[b]      -- DMFT self-energy at infinity
               Edc[b]       -- double counting using DMFT occupancies
        """
        if kbroad>0.0 or wbroad>0.0:
            bexe = self.env.ROOT + '/broad'
            shutil.move(self.dir+self.sparams['Sig'], self.dir+self.sparams['Sig']+'b')
            broad_cmd = 'cd '+self.dir+'; '+bexe+' -w '+str(wbroad)+' -k '+str(kbroad)+' '+self.sparams['Sig']+'b  > '+self.sparams['Sig']
            print(broad_cmd, file=self.fh_info)
            ret = subprocess.call(broad_cmd,shell=True,stdout=self.fh_info,stderr=self.fh_info)
        OffDiagonalExist = False
        Sigind = array(self.Sigind_shifted)
        self.Diagonal={}
        for i in range(len(Sigind)):
            for j in range(len(Sigind)):
                if Sigind[i,j]>0:
                    self.Diagonal[Sigind[i,j]-1]= (i==j)
                    if i!=j: OffDiagonalExist = True

        # Reading Self-energy of the impurity
        Sg_data = loadtxt(self.dir + self.sparams['Sig']).T
        om = Sg_data[0]
        Sigo = Sg_data[1::2,:]+Sg_data[2::2,:]*1j
        
        # Computing self-energy in alternative way to keep off-diagonal components of self-energy forcing G_{loc}=G_{imp} exactly
        if (self.sparams['OffDiagonal'] and OffDiagonalExist):
            Sign = self.SigmaOffDiagonal(om,Sigo,wbroad,kbroad)
            if Sign is not None: Sigo = Sign

        # WARNING: Gf/Sig.out contains nf, moment, ntot,... We changed recently from Sig->Gf
        with open(self.dir + self.sparams['Gf'],'r') as fS:
            adat = fS.readline().strip().split()
            
        for par in adat: # setting nf from ctqmc output
            m = re.search('nf=', par) or re.search('TrSigmaG=', par) or re.search('mom=', par)
            if m is not None : exec(par, globals())
        ntot = nf
        print('From impurity Sigma: nf=', nf, 'TrSigmaG=', TrSigmaG, 'mom=', mom, file=self.fh_info)
        
        for b in range(shape(Sigo)[0]):   #print b, 'This is diagonal bath, correcting it'
            if self.Diagonal[b]:
                for iom in range(shape(Sigo)[1]):
                    if Sigo[b,iom].imag>0: # violates causality, correct it.
                        print('Disregarding non-causal self-energy at Sig['+str(b)+','+str(iom)+']=',Sigo[b,iom], file=self.fh_info)
                        Sigo[b,iom] = Sigo[b,iom].real
        
        sinf = copy(Sigo[:,-1]).real
        
        print('sinf=', sinf, file=self.fh_info)
        
        # Reading old Edc and impurity levels
        (Eimp_old,Olap_old,Edc_old) = loadtxt(self.dir+'/Eimp.inp',ndmin=2)
        
        print('Diagonal=', self.Diagonal, file=self.fh_info)
        print('len(Edc_old)=', len(Edc_old), 'Eimp_old=', Eimp_old, file=self.fh_info)

        #if 'CoulombF' in params:
        #    CoulombF = params['CoulombF'][0].strip("'").strip('"')
        #else:
        #    CoulombF = Default_CoulombF
        
        # Here we compute self-energy(infinity) and 
        # double-counting using impurity occupancy
        print('DCs=', DCs, file=self.fh_info)
        Ud = params['U'][0]
        if type(self.J) is list:
            Jd = sum(self.J)/len(self.J)
        else:
            Jd = self.J
        # If we want to have U different than the U used in double counting
        # (useful for cluster-DMFT when U can not be simply added in ctqmc)
        if 'Ud' in params: Ud = params['Ud'][0]
        if 'Jd' in params: Jd  = params['Jd'][0]
        
        DC_ones = array([int(self.Diagonal[i]) for i in range(len(Edc_old))])

        if DCs[:5]=='exact':  
            # We allow DCs = ['exact', 'exact_l', 'exactd', 'exactdl']
            # exact and exactd use impurity occupation. exact uses yukawa&dielectric, while exactd uses dielectric screening only
            # exact_l and exactdl use lattice occupations. exact_l uses yukawa&dielectric, while exactdl uses dielectric screening only
            lmbda = self.UlamJ[1]
            epsilon = self.UlamJ[2]
            nf_ = mom               # 'exact': Density determined from impurity occupation. It is way more stable.
            print('Here we start n_imp=', nf_, 'isscalar(nlattice)=', isscalar(nl_imp), file=self.fh_info)
            if len(DCs)>6 and DCs[6]=='l' and not isscalar(nl_imp): 
                nf_ = nl_imp
            print('Here we have nl_imp=', nf_, file=self.fh_info)
            print('lmbda=', lmbda, 'epsilon=', epsilon, file=self.fh_info)
            (Edc0,Phidc0)=dmft_DC.ComputeExactDC(nf_, lmbda, epsilon, self.icase, self.fh_info, projector='projectorw.dat',trans=self.dir+'Trans.dat')
            print('The Exact Double-counting is: Vdc=', Edc0, 'PhiDC=', Phidc0, file=self.fh_info)

        elif DCs=='FLL' or DCs=='default' or DCs=='FLL_l':
            nf_ = nf
            if DCs=='FLL_l' and nf_: 
                nf_ = sum(nl_imp)
            # This is the fully localized limit double-counting by Anisimov
            Edc0 = DC_ones*( Ud*(nf_-0.5) - Jd*(nf_/2.-0.5) )
            Phidc0 = Ud*0.5*nf_*(nf_-1.) - Jd*0.25*nf_*(nf_-2.)
        elif DCs=='AMF':
            nf_ = nf
            if DCs=='AMF_l' and nf_: 
                nf_ = sum(nl_imp)
            Ueff = Ud*(1-1./(2.*(2.*self.l+1.))) - Jd*self.l/(2*self.l+1.)
            Edc0 = DC_ones*Ueff*nf_
            Phidc0 = Ueff*nf_**2/2.
        elif DCs=='fixn' or DCs=='nominal':         # If the scheme if fixn, we keep double-counting potential fixed
            nf0 = params['nf0'][0]
            Vdc = Ud*(nf0-0.5) - Jd*(nf0/2.-0.5)
            Edc0 = DC_ones * Vdc
            Phidc0 = ntot*Vdc
        else:
            print('ERROR: Unkown double-counting', DCs)
            sys.exit(1)
            
        if 'dDC' in params:
            Edc0 += array(params['dDC'][0])
            Phidc0 += sum([params['dDC'][0][i]*nf_[i] for i in range(len(nf_))])
            print('The Double-counting corrected to: Vdc=', Edc0, 'PhiDC=', Phidc0, file=self.fh_info)

        Edc = Edc_old # we want to keep Edc<-1000 for orbitals which are projected out!
        for i in range(len(Edc)):
            if (Edc_old[i]>-1000.): # this is true for all orbitals which are kept!
                Edc[i] = Edc0[i]
            else:                   # these orbitals are projected out -> we want to keep sinf=Edc for orbitals projected out!
                sinf[i] = Edc[i]
                

        # Copying data at each iteration to follow the evolution
        shutil.copy2(self.dir+self.sparams['Sig'], self.dir+self.sparams['Sig']+'.'+extn)
        shutil.copy2(self.dir+self.sparams['Gf'], self.dir+self.sparams['Gf']+'.'+extn)
        shutil.copy2(self.dir+'ctqmc.log', self.dir+'ctqmc.log'+'.'+extn)
        shutil.copy2(self.dir+'Probability.dat', self.dir+'Probability.dat'+'.'+extn)
        
        print('nimp, TrSigmaG=', ntot, TrSigmaG, file=self.fh_info)
        
        (Phi_DMFT, lnGimp, Fimp, Vdc_nd, zeorb) = GiveImpFunctional(self.dir,self.PARAMS,Edc,self.fh_info)
        
        self.fh_info.flush()
        
        Ry2eV = 13.60569193
        Epotential = TrSigmaG
        
        with open(self.dir+'/Eorb.dat', 'w') as fE:
            print('#  Phidc=', Phidc0, file=fE)
            print('#  Tr(Sigma*G)/2=', Epotential, file=fE)
            print('#  Tr(Sigma*G)/2-Phidc=', (Epotential-Phidc0), file=fE)
            print('#  Phi_DMFT=', Phi_DMFT, file=fE)
            print('#  Phi_DMFT-Phi_DC=', Phi_DMFT-Phidc0, file=fE)
            print('#  Tr(log(-Gimp))=', lnGimp, file=fE)
            print('#  Fimpurity+TS=', Fimp, file=fE)
            print('#  E_pot+E_kin-Tr(w*dD/dw)=', zeorb, file=fE)
            print(':IEORB ', (Epotential-Phidc0)/Ry2eV, file=fE)
            print(':XEORB ', (Phi_DMFT-Phidc0)/Ry2eV, file=fE)
            print(':ZEORB ', (zeorb-Phidc0)/Ry2eV, file=fE)

        with open(self.dir+'/sig.out', 'w') as fE:
            print('# s_oo=', sinf.tolist(), file=fE)
            print('# Edc=', Edc.tolist(), file=fE)
            for iom in range(len(om)):
                print(("%20.15f " % om[iom]), end=' ', file=fE)
                for b in range(len(Sigo)):
                    print(("%20.15f %20.15f  " % (Sigo[b,iom].real-Edc[b], Sigo[b,iom].imag)), end=' ', file=fE)
                print(file=fE)
        
        return ntot

    def FromSingleImpurity_to_Link(self, Eimps, params):
        # First create local cix for single-impurity .
        Eimps_local=[]
        for i in range(int(len(Eimps)/2)):
            Eimps_local.append( 0.5*(Eimps[2*i]+Eimps[2*i+1]) )
        
        print('Eimps_local=', Eimps_local, file=self.fh_info)

        # Creats local Counterparts for Trans.dat
        Nt = len(self.CF)
        Sigind_local=zeros((int(Nt/2),int(Nt/2)),dtype=int)
        for i in range(int(Nt/2)):
            for j in range(int(Nt/2)):
                Sigind_local[i,j] = self.Sigind[2*i+1,2*j+1]/2

        floc = open(self.dir+'/Trans.dat_local','w')
        # Creates local CF
        Cp = zeros((int(Nt/2),Nt),dtype=complex)
        Cm = zeros((int(Nt/2),Nt),dtype=complex)
        p_is_local=True
        m_is_local=True
        for i in range(0,int(Nt/2)):
            Cp[i,:] = (self.CF[2*i,:]+self.CF[2*i+1,:])/sqrt(2)
            Cm[i,:] = (self.CF[2*i,:]-self.CF[2*i+1,:])/sqrt(2)
            if sum(abs(Cp[i,int(Nt/2):Nt]))>1e-5: p_is_local=False
            if sum(abs(Cm[i,int(Nt/2):Nt]))>1e-5: m_is_local=False
        if p_is_local:
            CF_local = Cp[:,:int(Nt/2)]
        elif m_is_local:
            CF_local = Cm[:,:int(Nt/2)]
        else:
            print('ERROR: Something wrong with transformation T2C given in case.indmfl file')
            print('ERROR: Something wrong with transformation T2C given in case.indmfl file', file=fh_info)

        print(int(Nt/2), int(Nt/2), '# size of Sigind and CF', file=floc)
        print('#---- Sigind follows', file=floc)
        for i in range(len(Sigind_local)):
            print("%3d "*(len(Sigind_local)) % tuple(Sigind_local[i,:]), file=floc)
        print('#---- CF follows', file=floc)
        for i in range(len(CF_local)):
            for j in range(len(CF_local[i])):
                print("%12.8f "*2 % (CF_local[i,j].real, CF_local[i,j].imag), end=' ', file=floc)
            print(file=floc)
        floc.close()

        shutil.move(self.dir+'/Trans.dat', self.dir+'/Trans.dat_cluster')
        shutil.copy(self.dir+'/Trans.dat_local', self.dir+'/Trans.dat')

        atom_arg = ''
        for k,p in list(self.spr.items()):
            atom_arg += '\"'+k+'='+str(p)+'\" '
        atom_arg += ' \"Eimp='+str(Eimps_local)+'\" '
        
        print('atom_arg=', atom_arg, file=self.fh_info)
        print(atom_arg)
        # runs ED for single site DMFT
        cmd = 'cd '+ self.dir + '; ' + self.atom_exe +' ' +atom_arg+' >> nohup.out 2>&1'
        print('.... running ..', cmd, file=self.fh_info)
        print('.... running ..', cmd)
        ret = subprocess.call(cmd,shell=True,stdout=self.fh_info,stderr=self.fh_info)
        # runs link.py to compute cix file for 2-site cluster
        link_arg=' -U '+str(params['Ud'][0])+' -E "'+str(Eimps)+'"'
        cmd = 'cd '+ self.dir + '; ' + self.env.ROOT + '/link.py' +' ' +link_arg+' >> nohup.out 2>&1'
        print('.... running ..', cmd, file=self.fh_info)
        print('.... running ..', cmd)
        ret = subprocess.call(cmd,shell=True,stdout=self.fh_info,stderr=self.fh_info)
        
        shutil.move(self.dir+'/Trans.dat_cluster', self.dir+'/Trans.dat')
        shutil.move(self.dir+'/link_actqmc.cix', self.dir+'/actqmc.cix')
        
    def SigmaOffDiagonal(self,om,Sigo,wbroad,kbroad,fname='Glatt.imp'):
        glatfile = self.dir + '/'+fname
        if not os.path.isfile(glatfile) or os.path.getsize(glatfile)==0: return None
        
        # Reading Gimp, Glatt, and Delta
        Gimp_dat = loadtxt(self.dir + '/'+self.sparams['Gf'] ).T
        Glat_dat = loadtxt(glatfile).T
        Delta_dat = loadtxt(self.dir+'/Delta.imp').T
        # Two frequency meshes, logarithmic and linear
        om0 = Gimp_dat[0]
        Gimp0 = Gimp_dat[1::2] + Gimp_dat[2::2]*1j
        om1 = Delta_dat[0]
        Glat = Glat_dat[1::2] + Glat_dat[2::2]*1j
        Delta = Delta_dat[1::2] + Delta_dat[2::2]*1j
        # Find correspondence between linear and logarithmic mesh
        iom=0
        Gimp=[] # Gimp on logarithmic mesh
        for i in range(len(om1)):
            while (iom<len(om0) and om0[iom]-1e-6<om1[i] ): iom+=1
            if abs(om1[i]-om0[iom-1])>1e-4:
                print('Problems! Can not find the correspondence between logarithmic and linear meshin Delta.imp and Gf.out!')
            Gimp.append( Gimp0[:,iom-1] )
        Gimp = array(Gimp).T
            
        # Find which columns and rows are correlated?
        Sigind = array(self.Sigind_shifted)
        ind = [i for i in range(len(Sigind)) if Sigind[i,i]>0]
        #ind = []
        #for i in range(len(Sigind)):
        #    if Sigind[i,i]>0:
        #        ind.append(i)
        dim=len(ind)
        # How many times each correlated orbital appears
        noccur = zeros(len(Delta),dtype=int)
        for i in range(len(Sigind)):
            for j in range(len(Sigind)):
                if Sigind[i,j]>0:
                    noccur[Sigind[i,j]-1]+=1
        # Impurity levels in matrix form
        (Eimp,Olap,Edc) = loadtxt(self.dir+'/Eimp.inp',ndmin=2)   # impurity levels should be produced in run_dmft.py at imp.x/Eimp.inp
        # Reading Delta
        Eimps = Eimp-Edc                                     # impurity levels with double-counting included, as needed for ED & ctqmc
        Eimp = matrix(zeros( (dim,dim), dtype=complex))
        for i,ii in enumerate(ind):
            for j,jj in enumerate(ind):
                if (Sigind[ii,jj]>0):
                    Eimp[i,j] = Eimps[Sigind[ii,jj]-1]
        
        #print 'self.Eimps=', self.Eimps
        
        Gc=matrix(zeros((dim,dim),dtype=complex))
        Dlt=matrix(zeros((dim,dim),dtype=complex))
        Id = matrix(identity(dim,dtype=complex))
        Sigma = zeros( (len(Delta), len(om1) ), dtype=complex )

        #print 'len(Delta)=', len(Delta), 'shape(Sigma)=', shape(Sigma), 'len(noccur)=', len(noccur), 'len(ind)=', len(ind), 'ind=', ind
        for iom in range(len(om1)):
            for i,ii in enumerate(ind):
                for j,jj in enumerate(ind):
                    if (Sigind[ii,jj]>0):
                        Dlt[i,j] = Delta[Sigind[ii,jj]-1,iom]
                        Gc[i,j]  = Glat[Sigind[ii,jj]-1,iom]
                    Gc[i,i] = Gimp[Sigind[ii,ii]-1,iom]
                    
            Sigm = Id*(om1[iom]*1j)-Dlt-Eimp-Gc.I
            
            for i,ii in enumerate(ind):
                for j,jj in enumerate(ind):
                    if (Sigind[ii,jj]>0):
                        if (ii!=jj and self.sparams['OffDiagonal']=='real'):
                            Sigma[Sigind[ii,jj]-1, iom] += real(Sigm[i,j])
                        else:
                            Sigma[Sigind[ii,jj]-1, iom] += Sigm[i,j]
                            
            for i in range(len(Sigma)): Sigma[i,iom] *= 1./noccur[i]

        savetxt(self.dir+'S1.dat', vstack((om,  Sigo.real,  Sigo.imag )).T )
        savetxt(self.dir+'S2.dat', vstack((om1, Sigma.real, Sigma.imag)).T )

        #for ib in range(len(Sigma)):
        #    Sigma[ib] = brod.Broad(wbroad, kbroad, om1, Sigma[ib])
        #savetxt(self.dir+'S3.dat', vstack( (om1, Sigma.real, Sigma.imag ) ).T )
        if kbroad>0.0 or wbroad>0.0:
            bexe = self.env.ROOT + '/broad'
            broad_cmd = 'cd '+self.dir+'; '+bexe+' -w '+str(wbroad)+' -k '+str(kbroad)+' S2.dat  > S3.dat'
            print(broad_cmd, file=self.fh_info)
            ret = subprocess.call(broad_cmd,shell=True,stdout=self.fh_info,stderr=self.fh_info)
            dat = loadtxt(self.dir+'S3.dat').T
            N = len(Sigma)
            Sigma = dat[1:N+1] + dat[N+1:2*N+1]*1j

        nom = self.sparams['nom']
        # Interpolate om big mesh
        Sign=[]
        for ib in range(len(Sigma)):
            fSr = interpolate.UnivariateSpline(om1, Sigma[ib].real, s=0)
            fSi = interpolate.UnivariateSpline(om1, Sigma[ib].imag, s=0)
            Sg = fSr(om) + fSi(om)*1j
            if not self.Diagonal[ib]:
                Sg[nom:] *= exp(-(om[nom:]-om[nom]))
            Sign.append( Sg )
            
        Sign = array(Sign)

        #(Eimp,Olap,Edc) = loadtxt(self.dir+'/Eimp.inp',ndmin=2)
        with open(self.dir+'/test.out', 'w') as fE:
            #print(fE, '# s_oo='+str(sinf.tolist()), file=fE)
            #print(fE, '# Edc='+str(Edc.tolist()), file=fE)
            for iom in range(len(om)):
                print("{:20.15f} ".format(om[iom]), end=' ', file=fE)
                for b in range(len(Sign)):
                    print("{:20.15f} {:20.15f}  ".format(Sign[b,iom].real-Edc[b], Sign[b,iom].imag), end=' ', file=fE)
                print(file=fE)
                
        return Sign

            
                
def ferm(x):
    """Fermi function for T=1"""
    if (x>100): return 0
    if (x<-100): return 1
    return 1/(exp(x)+1)


if __name__ == '__main__':
    from readTrans import ReadTrans
    iparams0={"exe"     : ["ctqmc"         , "# Name of the executable"],
          "U"           : [5.0             , "# Coulomb repulsion (F0)"],
          "Ud"          : [5.0             , "# Coulomb repulsion (F0)"],
          "J"           : [0.8             , "# Coulomb repulsion (F0)"],
          "nf0"         : [1               , "# Double counting parameter"],
          "beta"        : [35.0            , "# Inverse temperature"],
          "Nmax"        : [2000            , "# Maximum perturbation order allowed"],
          "M"           : [1e8             , "# Total number of Monte Carlo steps"],
          "Ncout"       : [1000000         , "# How often to print out info"],
          "Naver"       : [1000000000      , "# How often to print out debug info"],
          "nom"         : [50              , "# Number of Matsubara frequency points sampled"],
          "aom"         : [4               , "# Number of frequency points used to determin the value of sigma at nom"],
          "sderiv"      : [0.02            , "# Maximum derivative mismatch accepted for tail concatenation"],
          "Ntau"        : [1000            , "# Number of imaginary time points (only for debugging)"],
          "SampleGtau"  : [1000            , "# How often to update G(tau)"],
          "GlobalFlip"  : [500000          , "# How often to try a global flip"],
          "tsample"     : [10              , "# How often to record measurements"],
          "warmup"      : [100000          , "# Warmup number of QMC steps"],
          "CleanUpdate" : [100000          , "# How often to make clean update"],
          "minM"        : [1e-10           , "# The smallest allowed value for the atomic trace"],
          "minD"        : [1e-10           , "# The smallest allowed value for the determinant"],
          "PChangeOrder": [0.9             , "# Ratio between trial steps: add-remove-a-kink / move-a-kink"],
          "CoulombF"    : ["'Full'"        , "# Georges = rough run"],
          "atom_OCA_G"  : [False           , "# Don't compute OCA diagrams for speed"],
          "UlamJ"       : [[5.0, 1.360508, 1.4459, [0.7506725, 0.8493274]], "#"],
          "icase"       : [0               , "# added by run_dmft.py"],
    }

    
    dire = 'imp.0'
    (Sigind, CF) = ReadTrans(dire+'/Trans.dat', sys.stdout)

    imp = IMP_CTQMC(23, dire, iparams0, Sigind, CF)
    imp._exe(iparams0, 'exact', '_bris', 1)
    
    ntot=imp.HighFrequency(iparams0, 'exact', 0.0, '_bris', 0.05, 0.1)
    print('ntot=', ntot)
