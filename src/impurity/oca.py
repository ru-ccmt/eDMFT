#!/usr/bin/env python
# @Copyright 2007 Kristjan Haule
# 
import sys, re, os, glob, shutil, types
from numpy import *
from utils import DmftEnvironment

#ROOT = os.environ.get('WIEN_DMFT_ROOT')
#if ROOT is not None:
#    sys.path.append( ROOT )
#else:
#    print >> sys.stderr, "Environment variable WIEN_DMFT_ROOT must be set!"
#    print "Environment variable WIEN_DMFT_ROOT must be set!"
#    sys.exit(1)
import ldau
import krams
import strns as sts
import linalg as la


class IMP_OCA(object):
    """One crossing approximation Impurity solver
    
    """
    def __init__(self, Znuc, dire, params, Sigind, CF): #Impq=None):
        self.env = DmftEnvironment()
        self.mpi_prefix = self.env.MPI
        
        # impurity directory
        self.dir = dire+'/'

        # create impurity directory if it does not yet exist
        if len(glob.glob(dire))==0 : os.mkdir( dire )

        # log-file is opened
        self.fh_info = open(self.dir + 'info.out', 'w')


        self.Sigind = array(Sigind)
        self.CF = CF
        # find smallest sigind to shift siginds so they start from 1
        minsigind = min([s for s in self.Sigind.flat if s != 0])
        self.Sigind_shifted = where(self.Sigind != 0, self.Sigind-minsigind+1, 0)
        
        # Transformation between spheric harmonics and local coordinates needs to be
        # used in exact diagonalization of the atom.
        fh_T = open(self.dir + 'Trans.dat', 'w')
        print(len(Sigind), len(Sigind[0]), '#  size of Sigind and CF', file=fh_T)
        print('#---- Sigind follows', file=fh_T)
        for row in self.Sigind_shifted:
            for elem in row:
                print(("%3s " % elem), end=' ', file=fh_T)
            print(file=fh_T)
            
        print('#---- CF follows', file=fh_T)
        for i in range(len(CF)):
            for j in range(len(CF[0])):
                print("%12.8f %12.8f " % (CF[i][j].real, CF[i][j].imag), end=' ', file=fh_T)
            print(file=fh_T)
        
        
        # parameters, which will be given to oca and will be written to PARAMS.oca
        self.sparams={'Sig':'Sigma.000', 'Ac':'Ac.inp', 'cix':'out.cix', 'AlocOut':'Aloc.imp', 'SigOut':'Sigma.000', 'sig': 'sig.inp', 'gloc':'gloc.out', 'pcore':0, 'acore':0}

        # copying Sigma.000 to imp.x subdirectory
        if len(glob.glob(self.dir+'/'+self.sparams['Sig']))==0: 
            shutil.copy2('Sigma.000', self.dir+self.sparams['Sig'])

        
        # file from database which contains atomic information for this particular atom
        self.atom_arg = self.env.ROOT + '/database/' + 'atom.'+str(Znuc)+'.py'
        # executable which broadens
        self.broad = self.env.ROOT + '/broad'
        # input file for oca
        self.PARAMS = 'PARAMS.oca'

        # information from database is read
        exec(compile(open(self.atom_arg, "rb").read(), self.atom_arg, 'exec'))
        database = locals() # it seems python does not know about variables obtained by execfile
        

        # The following variables will be set by the following procedure:
        #   if the variable is defined in params, it takes its value from params otherwise it checks its definition in database
        vars = ['l', 'para', 'qOCA', 'kOCA', 'mOCA', 'Ex', 'Ep', 'J', 'cx', 'Eoca', 'qatom', 'n', 'CoulombF']
        
        self.spr={}
        for var in vars:
            if var in params:
                self.spr[var] = params[var][0]
            elif var in database:
                self.spr[var] = database[var]

        # Some variables have different name in database and input file
        # spr['n'] -- params['nc']  -- database['n']
        # BE CAREFUL - RECENT CHANGE Dec.2009
        if 'nc' in params:
            self.spr['n'] = params['nc'][0]
        else:
            self.spr['n'] = database['nc']
        
        if 'Ncentral' in params:
            self.spr['Ncentral'] = params['Ncentral'][0]
        else:
            self.spr['Ncentral'] = self.spr['n'][1:len(self.spr['n'])-1]
        
        print('self.spr=', self.spr, file=self.fh_info)
        
        #if Impq is not None:
        #    if type(Impq)!= list :  Impq = Impq.tolist()
        #    self.Impq = Impq
        #    print >> self.fh_info, 'Impq=', Impq
        
        # stores some of just defined variables
        self.l = self.spr['l']
        self.J = self.spr['J']
        self.cx = self.spr['cx']

        
        

    def _exe(self, params, DCs, extn, UpdateAtom, gbroad, kbroad=0.0, maxAc=200.):
        """ Calls exact diagonalization of the atom and executes the oca/nca impurity solver
        """

        # Reading impurity levels
        (Eimp,Olap,Edc) = loadtxt(self.dir+'/Eimp.inp')
        # Reading Delta
        Delta_dat = loadtxt(self.dir+'/Delta.imp').transpose()
        om = Delta_dat[0]
        Delta = (Delta_dat[1::2] + Delta_dat[2::2]*1j).transpose()
        
        # subtracting Edc, because the ksum used (Hk+s_oo-Edc)
        Eimps = Eimp-Edc 

        if (DCs=='fixEimp'):
            dE = params['Eimp0'][0]-Eimps[0]
            Eimps += dE
            Edc -= dE
            fE = open(self.dir+'/Eimp.inp', 'w')
            print("%15.10f "*len(Eimp) % tuple(Eimp), ' # Eimp', file=fE)
            print("%15.10f "*len(Olap) % tuple(Olap), ' # Overlap', file=fE)
            print("%15.10f "*len(Edc) % tuple(Edc), ' # Edc', file=fE)
            fE.close()
            

        params['Ed'] = [Eimps.tolist(), "# Impurity levels"]


        #print 'HERE', UpdateAtom, len(glob.glob(self.dir+'/'+self.sparams['cix']))==0

        # Below exact diagonalization of the atom is run to produce 'cix.out' file
        # ED is performed if 'cix'-file does not yet exist, or, UpdateAtom=True
        if UpdateAtom or len(glob.glob(self.dir+'/'+self.sparams['cix']))==0 : # no cix file yet -> create it
        
            if (self.l == 3) :
                # executable for ED of the atom
                self.atom_exe = self.env.ROOT + '/atom'
                # argument list is created to later run the executable for exact diagonalization of the atom

                print('exe=', self.atom_exe)
                
                atom_arg=''
                for k,p in list(self.spr.items()):
                    if type(p)==list:
                        atom_arg += k+'="'+str(p)+'" '
                    else:
                        atom_arg += k+'='+str(p)+' '

                print('atom_arg=', atom_arg, file=self.fh_info)

                cmd = 'cd '+ self.dir + '; ' + self.atom_exe +' ' +atom_arg+' > nohup.out 2>&1'
                print('command: ', cmd, file=self.fh_info)
                # runs ED for the atom
                stdin, stdout, stderr = os.popen3( cmd )
                print(stdout.read(), file=self.fh_info)
                print(stderr.read(), file=self.fh_info)
            elif (self.l==2) :
                # executable for ED of the atom
                self.atom_exe = self.env.ROOT + '/atom_d.py'

                if type(Eimps)!= list :  Eimps = Eimps.tolist()
                print('Eimps=', Eimps)
                print('Ed=', params['Ed'])
                
                atom_arg = self.atom_arg + ' \"Eimp='+str(Eimps)+'\" \"n='+str(self.spr['n'])+'\" \"Ncentral='+str(self.spr['Ncentral'])+'\" \"qOCA='+str(self.spr['qOCA'])+'\"'+' \"para='+str(self.spr['para'])+'\"'

                print('atom_arg=', atom_arg, file=self.fh_info)
                cmd = 'cd '+ self.dir + '; ' + self.atom_exe +' ' +atom_arg+' > nohup.out 2>&1'
                print('command: ', cmd, file=self.fh_info)
                
                # runs ED for the atom
                stdin, stdout, stderr = os.popen3( cmd )
                print(stdout.read(), file=self.fh_info)
                print(stderr.read(), file=self.fh_info)
            else:
                print('l=', self.l, ' not yet implemented!', file=self.fh_info)
        
            
        # OCA executable
        shutil.copy2(self.env.ROOT + '/' + params['exe'][0], self.dir+params['exe'][0])
        self.oca_exe = './'+params['exe'][0]

        if 'qatom' in self.spr and self.spr['qatom']:
            # impurity level for 5/2 is E_{5/2}=E0-2*cx and for 7/2 is E_{7/2}=E0+3/2*cx
            # Here we correct E_{5/2} for 2*cx, which is already counted in cix file
            params['Ed'][0][0] = (params['Ed'][0][0]+2*self.cx)
            for j in range(1,len(params['Ed'][0])): params['Ed'][0][j]=params['Ed'][0][0]
                                    
        
        # Preparing PARAMS.oca file
        fp = open(self.dir+self.PARAMS, 'w')
        for p in list(self.sparams.keys()):
            fp.write(p + '=' + str(self.sparams[p]) + '\n')
        
        for p in list(params.keys()):
            print("%s=%s " % (p, params[p][0]), '\t', params[p][1], file=fp)
        fp.close()

        # Preparing Ac.inp file        
        fp = open(self.dir+self.sparams['Ac']+'b', 'w')
        QcutAc = False
        if 'cutAc' in list(params.keys()):
            QcutAc = True
            cutAc = params['cutAc'][0]
            
        for im in range(len(om)):
            print(om[im], end=' ', file=fp)
            for b in range(len(Delta[im])):
                Ac = -Delta[im][b].imag/pi
                if (Ac<0.0): Ac=0
                if (Ac>maxAc): Ac=0
                if (QcutAc): Ac *= ferm((cutAc[0]-om[im])/params['Th'][0])*ferm((om[im]-cutAc[1])/params['Th'][0])
                print(Ac, end=' ', file=fp)
            print(file=fp)
        fp.close()

        if gbroad>0:
            broad_cmd = 'cd '+self.dir+'; '+self.broad+' -w '+str(gbroad)+' -k '+str(kbroad)+' '+self.sparams['Ac']+'b'+' > '+self.sparams['Ac']
            print(broad_cmd)
            stdin, stdout, stderr = os.popen3(broad_cmd)
            print(stderr.read())
        else:
            shutil.move(self.dir+self.sparams['Ac']+'b', self.dir+self.sparams['Ac'])
            
        shutil.copy2(self.dir+self.sparams['Ac'], self.dir+self.sparams['Ac']+'.'+extn)

        if 'external_exe' in params:
            print('Running external python code ', params['external_exe'][0])
            exec(compile(open(params['external_exe'][0], "rb").read(), params['external_exe'][0], 'exec'))
            
            
        
        # Below we execute oca
        mpi_cmd = self.mpi_prefix if params['exe'][0] == 'oca' else ''  # only oca is parallelized
        oca_cmd = 'cd '+self.dir+'; '+mpi_cmd+' '+self.oca_exe+' '+self.PARAMS+' > nohup_imp.out 2>&1 '
        print(oca_cmd)
        print('Running ---- OCA -----')
        stdin, stdout, stderr = os.popen3(oca_cmd)
        print(stdout.read(), stderr.read(), file=self.fh_info)
        
        
    def HighFrequency(self, params, DCs, nl_imp, extn, wbroad=0.0, kbroad=0.0, Q_ETOT=False):
        """ Corrects high-frequency of OCA spetral functions
        such that it gives correct nf, normalization and sigma_oo
        
        This is achieved by adding two Lorentzians at high frequency. One below EF in the interval (par['epsilon'][0][0],par[epsilon][0][1])
        and one above EF in the interval (par['epsilon'][1][0],par[epsilon][1][1])
        
           Input:
              params       --  parameters must contain
                               T
                               Ed
                               U
                               epsilon    - prefered positions of additional lorentzians
                               Th         - when high-frequency spectral function is cut, it is cut by fermi function with temperature Th
                               Gh         - the width of the two lorentzians added
           Output:
               Sig[b][iom]  -- DMFT dynamic self-energy which vanishes at infinity
               sinf[b]      -- DMFT self-energy at infinity
               Edc[b]       -- double counting using DMFT occupancies

           The 'corrected' Green's function will be of the form:
                      G(om) = int(A(x)/(om-x)) + a1/(om-eps1+i*Gh)  + a2/(om-eps2+i*Gh)
           
           The equations to be satisfied are:
              normalization:    m0 + a1 + a2 = 1
              density:          nc + a1 = nc_{exact}    where nc=int(A(x)f(x)) and nc_{exact} is computed from pseudo spectral functions
              Sigma_oo:         m1 + a1*eps1 + a2*eps2 = Eimp+Sigma_oo == dsinf

           The equations can be brought to the form
              x*u + y*v = w
           with u and v positive numbers and x and y unknown numbers in the interval [0,1]
           
           In terms of above quantities, we have
              u = a1*(eps[0][1]-eps[0][0])
              v = a2*(eps[1][1]-eps[1][0])
              w = Eimp+Sigma_oo - m1 - a1*eps[0][0] - a2*eps[1][0]
              x = (eps1-eps[0][0])/(eps[0][1]-eps[0][0])
              y = (eps2-eps[1][0])/(eps[1][1]-eps[1][0])

           The solution exists if 0<w<u+v
           In this case x is in the interval x=[max((w-v)/u,0), min(w/u,1)] and y is in the interval y=[max((w-u)/v,0), min(w/v,1)]
           The solution choosen is:
                 if (v>u):  x=0.5*(x_min+x_max)
                 else    :  y=0.5*(y_min+y_max)
        """
        ################################################
        # Reading of parameters from impurity cix file #
        ################################################
        # cix file is used to find the following variables: J, l, Ns, baths
        fc = open(self.dir + self.sparams['cix'])
        first_line = next(fc)
        # next line of cix contains Ns and baths
        cixdat = next(fc).split()
        if cixdat[0]=='OFF-DIAGONAL':
            Nl = int(next(fc).split()[0])
            for i in range(Nl):
                next(fc)
            cixdat = next(fc).split()
        
        baths = int(cixdat[0])
        Ns = list(map(int, cixdat[1:1+baths]))

        #print 'Ns=', Ns
        #print 'baths=', baths
        # Reading Aloc.imp which was produced by OCA
        fA = open(self.dir + self.sparams['AlocOut'])

        # Aloc.imp contains nf, moment, ntot,...
        first_line = next(fA).strip()
        adat = first_line.split()
        for par in adat:
            m = re.search('[ntot|nf|moment|dFimpG|Epot]=', par)
            if m is not None : exec(par)
        
        om=[]
        Af=[]
        for p in fA:
            dat = p.split()
            om.append(float(dat[0]))
            Af.append(list(map(float,dat[1:1+baths])))
        Af=array(Af)
        om=array(om)
        
        # reading of bath Weiss field
        fAc = open(self.dir + self.sparams['Ac'])
        Acd=[]
        ii=0
        for p in fAc:
            dat = p.split()
            omega = float(dat[0])
            if (abs(omega-om[ii])>1e-5):
                print('Seems that %s and %s are not compatible. Exiting!'%(self.sparams['AlocOut'],self.sparams['Ac']))
                sys.exit(1)
            ii += 1
            Acd.append(list(map(float,dat[1:1+baths])))
        Acd=array(Acd)

        # define some common variables in this functions
        T = params['T'][0]
        Ed = params['Ed'][0]
        epsilon = params['epsilon'][0]
        # Here we construct few functions for faster computation
        # of moments and occupancies
        #
        # dh    - mesh weights according to trapezoid rule
        # omdh  - omega*dh  for calculation of first moment
        # fedh  - f(omega)*dh for calculation of occupancy
        #
        dh=[0.5*(om[1]-om[0])]
        omdh=[om[0]*dh[-1]]
        fedh=[ferm(om[0]/T)*dh[-1]]
        for im in range(1,len(om)-1):
            dh.append(0.5*(om[im+1]-om[im-1]))
            omdh.append(om[im]*dh[-1])
            fedh.append(ferm(om[im]/T)*dh[-1])
        dh.append(0.5*(om[-1]-om[-2]))
        omdh.append(om[-1]*dh[-1])
        fedh.append(ferm(om[-1]/T)*dh[-1])
        dh=array(dh)
        omdh=array(omdh)
        fedh=array(fedh)
        

        #da = 0
        #if DCs['scheme']=='default' : da = DCs['a']
        
        # Here we compute self-energy(infinity) and 
        # double-counting using impurity occupancy
        (sinf, Edc) = Sinftyv(self.Sigind_shifted, self.CF, params['U'][0], self.J, self.l, baths, Ns, nf, 0, self.fh_info)
        
        Sigind = array(self.Sigind_shifted)
        Diagonal={}
        for i in range(len(Sigind)):
            for j in range(len(Sigind)):
                if Sigind[i,j]>0: Diagonal[Sigind[i,j]-1]= (i==j)


        Ud = params['U'][0]
        Jd = self.J
        DC_ones = array([int(Diagonal[i]) for i in range(len(Edc))])
        # If user fixes double counting by certain predescribed nf0
        # (fixn=True), we need to normalized actual nf[:] such that
        # sum(nf[:]) is equal to predefined nf0
        # Here we do not change Sigma_{infinity} but only Edc.
        # Sigma_{infinity} is determined by actual nf, while Edc is determined by modified nf.
        if DCs=='fixn' or DCs=='nominal':
            nf0 = params['nf0'][0]
            Vdc = Ud*(nf0-0.5) - Jd*(nf0/2.-0.5)
            Edc = DC_ones * Vdc
            Phidc = ntot*Vdc
            print('# Edc=', Edc, file=self.fh_info)
            print('# PhiDc=', Phidc, file=self.fh_info)
        elif DCs=='FLL':
            nf0 = sum(nf)
            Vdc = Ud*(nf0-0.5) - Jd*(nf0/2.-0.5)
            Edc = DC_ones * Vdc
            Phidc = Ud*0.5*nf0*(nf0-1.) - Jd*0.25*nf0*(nf0-2.)
            print('# Edc=', Edc, file=self.fh_info)
            print('# PhiDc=', Phidc, file=self.fh_info)
        elif DCs=='AMF':
            nf0 = sum(nf)
            Ueff = Ud*(1-1./(2.*(2.*self.l+1.))) - Jd*self.l/(2*self.l+1.)
            Edc = DC_ones*Ueff*nf0
            Phidc = Ueff*nf0**2/2.
            print('# Edc=', Edc, file=self.fh_info)
            print('# PhiDc=', Phidc, file=self.fh_info)
        elif DCs=='fixEimp':
            (Eimp,Olap,Edc) = loadtxt(self.dir+'/Eimp.inp')
            Phidc = Edc[0]*ntot
            
        # Some info up to now
        print('# l=%d T=%f U=%f J=%f' % (self.l, T, params['U'][0], self.J), file=self.fh_info)
        print('# Eimp=', Ed, file=self.fh_info)
        print('# baths=%d'%baths, 'Ns=', Ns, file=self.fh_info)
        print('# ntot=%f'%ntot, file=self.fh_info)
        print('# nf=', nf, file=self.fh_info)
        print('# moment=', moment, file=self.fh_info)
        print('# sinfty=', sinf.tolist(), file=self.fh_info)
        self.fh_info.flush()

        _Af=[]
        _Gf=[]
        _Sig=[]
        _Delta=[]
        _eps1=[]
        _eps2=[]
        _a1=[]
        _a2=[]
        for b in range(baths): # over all components of spectral functions (over baths)
        
            nf_exact = nf[b]/Ns[b]
            dsinf = sinf[b]+Ed[b]

            # the limit of vanishing spectral weight, i.e., when we have only two poles
            # this formula should always be obeyed
            #if (dsinf<0 and epsilon[0][0]*nf_exact+epsilon[1][0]*(1-nf_exact)>dsinf):
            #    print >> self.fh_info, 'epsilon was not choosen correctly and the solution can not be found!'
            #    epsilon[1][0] = (dsinf - epsilon[0][0]*nf_exact)/(1-nf_exact) - 0.5
            #    print >> self.fh_info, 'setting epsilon[1][0] to ', epsilon[1][0]
            #if (dsinf>0 and epsilon[0][1]*nf_exact+epsilon[1][1]*(1-nf_exact)<dsinf):
            #    print >> self.fh_info, 'epsilon was not choosen correctly and the solution can not be found!'
            #    epsilon[0][1] = (dsinf - epsilon[1][1]*(1-nf_exact))/nf_exact + 0.5
            #    print >> self.fh_info, 'setting epsilon[0][1] to ', epsilon[0][1]
               
            Afo = array(Af[:,b])
            
            m0 = dot(Afo,dh)
            m1 = dot(Afo,omdh)
            nc = dot(Afo,fedh)
            print('#start [%d]: m0=%f m1=%f nc=%f nf_exact=%f' % (b, m0, m1, nc, nf_exact), file=self.fh_info)

            # By default we correct N, and normalization, but not s_oo!
            Correct_N = True
            Correct_soo=False
            if 'correct' in params and params['correct'][0]=='s_oo':
                Correct_soo=True
            if 'correct' in params and params['correct'][0]==None:
                Correct_N=False
                Correct_soo=False

            print('Correct_N=', Correct_N)
            print('Correct_soo=', Correct_soo)
            
            if Correct_N:
                
                small = 1e-3
                a1 = nf_exact - nc
                if (a1<-small):
                    # weight below Ef is to large -> need to cut it a bit
                    for im in range(len(om)):
                        (a1n, a2n, m1n) = trycutAf(om[im], om, Afo, dh, omdh, fedh, nf_exact, params['Th'][0], self.fh_info)
                        print('#ca1 [%d]: cat L=%f a1=%f a2=%f m1=%f' % (b, om[im], a1n, a2n, m1n), file=self.fh_info)
                        if a1n>0:
                            L1 = om[im]
                            break
                    (a1, a2, m0, m1, nc) = cutAf(L1, om, Afo, dh, omdh, fedh, nf_exact, params['Th'][0], self.fh_info)
                    
                a2 = 1-m0-a1
                if (a2<-small):
                    # lorentzian weight above Ef is to large -> need to cut it a bit
                    for im in range(len(om)-1,0,-1):
                        (a1n, a2n, m1n) = trycutAf(om[im], om, Afo, dh, omdh, fedh, nf_exact, params['Th'][0], self.fh_info)
                        print('#ca2 [%d]: cat L=%f a1=%f a2=%f m1=%f' % (b, om[im], a1n, a2n, m1n), file=self.fh_info)
                        if a2n>0:
                            L2 = om[im]
                            break
                    (a1, a2, m0, m1, nc) = cutAf(L2, om, Afo, dh, omdh, fedh, nf_exact, params['Th'][0], self.fh_info)
                
                # Find u,v,w for this set of parameters
                (a1, a2, u, v, w) = uvw(m0, m1, nf_exact, nc, dsinf, epsilon)
                
                print('# [%d]: miss-nf=%f  miss-weight=%f  s_oo=%f a1=%f a2=%f u=%f v=%f w=%f' % (b, nf_exact-nc, 1-m0, sinf[b], a1, a2, u, v, w), file=self.fh_info)
                self.fh_info.flush()
                
            if Correct_soo:
	            # Here we compute the positions of the two lorentzian
                	# which will allow the exact density, normalization and sigma_oo
	            (success, x, y) = Soluvw(u, v, w)
            else:
                success = True
                m0 = dot(Afo,dh)
                m1 = dot(Afo,omdh)
                nc = dot(Afo,fedh)				    
                a1 = nf_exact - nc
                a2 = 1-m0-a1
                x=0
                y=1
                
                # If is not possible to get exact density, normalization and sigma_oo by adding
                # two lorentzians
                # Will try to cut high-frequency spectral function
                Lb=None
                if (not success):
                    # Here we cut the spectral function to make it more symmetric
                    
                    # start cutting at large frequency
                    # cuts until the condition is met
                    L0 = om[-1] # cut at positive frequency
                    (success, a1, a2, x, y) = ww(L0, params['Th'][0], om, Afo, dh, omdh, fedh, nf_exact, epsilon, dsinf, self.fh_info)
                    L0 /= 1.05
                    for j in range(100):
                        (success, a1, a2, x, y) = ww(L0, params['Th'][0], om, Afo, dh, omdh, fedh, nf_exact, epsilon, dsinf, self.fh_info)
                        if success: break
                        L0 /= 1.05
                    (successn, a1n, a2n, xn, yn) = ww_cut(L0, params['Th'][0], om, Afo, dh, omdh, fedh, nf_exact, epsilon, dsinf, self.fh_info)
                    
                if (not success):
                    print("Can't determin a way to get exact nf, norm and sigma_oo. You have to figure out the way to do that!")
                    sys.exit(1)
                    
                print("# [%d] a1=%f a2=%f x=%f y=%f" %(b, a1, a2, x, y), file=self.fh_info)
                
                eps1_ = epsilon[0][0] + x*(epsilon[0][1]-epsilon[0][0])
                eps2_ = epsilon[1][0] + y*(epsilon[1][1]-epsilon[1][0])
                
                print('# [%d]: a1=%f a2=%f eps1=%f  eps2=%f' % (b, a1, a2, eps1_, eps2_), file=self.fh_info)
                print(file=self.fh_info)
                
                
                # Actual cutting of the functions and adding the Lorentzians
                #Afn = Afo                     # use cutted function
                if ('FUN' in params and params['FUN']=='LOR'):
                    # Adding the Lorentzians
                    if a1>small:
                        for i in range(len(om)): 
                            Afo[i] += -(a1/pi/(om[i]-eps1_+params['Gh'][0]*1j)).imag() 
                    if a2>small:
                        for i in range(len(om)): 
                            Afo[i] += -(a2/pi/(om[i]-eps2_+params['Gh'][0]*1j)).imag()
                else:
                    if a1>small:
                        for i in range(len(om)):      # Adding the Lorentzians
                            Afo[i] += a1 * exp(-((om[i]-eps1_)/params['Gh'][0])**2)/sqrt(pi*params['Gh'][0]**2) 
                    if a2>small:
                        for i in range(len(om)):
                            Afo[i] += a2 * exp(-((om[i]-eps2_)/params['Gh'][0])**2)/sqrt(pi*params['Gh'][0]**2)
                
                _eps1.append(eps1_)
                _eps2.append(eps2_)
                _a1.append(a1)
                _a2.append(a2)
            
            Ac = array(Acd[:,b])
            gf=[]
            sig=[]
            Delta=[]
            for i in range(len(om)):
                gr = -pi*krams.kramarskronig(Afo, om, i)
                gc = gr - pi*Afo[i]*1j
                deltar = -pi*krams.kramarskronig(Ac, om, i)
                delta = deltar - pi*Ac[i]*1j
                
                gf.append(gc)
                Delta.append(delta)
                sigm = om[i] - Ed[b] - 1/gc - delta - sinf[b]
                if sigm.imag > 0:
                    sigm = sigm.real -1e-12j
                sig.append( sigm )
                
            _Af.append(Afo)
            _Gf.append(gf)
            _Sig.append(sig)
            _Delta.append(Delta)
            
        


        if Correct_N:
            print(file=self.fh_info)
            for b in range(baths):
                print('## [%d]: eps1=%f  eps2=%f  a1=%f  a2=%f' % (b, _eps1[b], _eps2[b], _a1[b], _a2[b]), file=self.fh_info)
            self.fh_info.flush()
        
        shutil.move(self.dir+self.sparams['gloc'], self.dir+self.sparams['gloc']+'_o')
        shutil.move(self.dir+self.sparams['sig'], self.dir+self.sparams['sig']+'_o')
        shutil.move(self.dir+self.sparams['AlocOut'], self.dir+self.sparams['AlocOut']+'_o')

        
        fg = open(self.dir+self.sparams['gloc'], 'w')
        fs = open(self.dir+self.sparams['sig'], 'w')
        fa = open(self.dir+self.sparams['AlocOut'], 'w')
        fd = open(self.dir+'Delta.out', 'w')
        print('# s_oo=', sinf.tolist(), file=fs)
        print('# Edc=', Edc.tolist(), file=fs)
        print(first_line, file=fa)
        for i in range(len(om)):
            print(om[i], end=' ', file=fg)
            for b in range(baths): print(_Gf[b][i].real, _Gf[b][i].imag, end=' ', file=fg)
            print(file=fg)
            print(om[i], end=' ', file=fs)
            for b in range(baths): print(_Sig[b][i].real, _Sig[b][i].imag, end=' ', file=fs)
            print(file=fs)
            print(om[i], end=' ', file=fa)
            for b in range(baths): print(_Af[b][i], end=' ', file=fa) 
            print(file=fa)
            print(om[i], end=' ', file=fd)
            for b in range(baths): print(_Delta[b][i].real, _Delta[b][i].imag, end=' ', file=fd) 
            print(file=fd)
            
        fg.close()
        fs.close()
        fa.close()
        fd.close()

        shutil.copy2(self.dir+self.sparams['AlocOut'], self.dir+self.sparams['AlocOut']+'.'+extn)
        shutil.copy2(self.dir+self.sparams['AlocOut']+'_o', self.dir+self.sparams['AlocOut']+'_o.'+extn)

        #shutil.copy2(self.dir+self.sparams['sig'], self.dir+self.sparams['sig']+'.'+extn)
        if os.path.exists(self.dir+'oca_log.000'): 
            shutil.copy2(self.dir+'oca_log.000', self.dir+'oca_log.000'+'.'+extn)
        else:
            shutil.copy2(self.dir+'nohup_imp.out', self.dir+'nohup_imp.out'+'.'+extn)
            
        Ry2eV = 13.60569193
        fE = open(self.dir+'/Eorb.dat', 'w')
        print('  Tr(Sigma*G)/2=', Epot, file=fE)
        print('  Phidc=', Phidc, file=fE)
        print('  Fimp-TrlogGimp=', dFimpG, file=fE)
        print('  Tr(Sigma*G)/2-Phidc=', Epot-Phidc, file=fE)
        if (Q_ETOT):
            print(':EORB ', (Epot-Phidc)/Ry2eV, file=fE)
        else:
            print(':EORB ', dFimpG/Ry2eV, file=fE)
        fE.close()
        

        fE = open(self.dir+'/sig.out', 'w')
        print('# s_oo=', sinf.tolist(), file=fE)
        print('# Edc=', Edc.tolist(), file=fE)

        for iom in range(len(om)):
            print(("%20.15f " % om[iom]), end=' ', file=fE)
            for b in range(len(_Sig)):
                print(("%20.15f %20.15f  " % (_Sig[b][iom].real+sinf[b]-Edc[b], _Sig[b][iom].imag)), end=' ', file=fE)
            print(file=fE)
        
        #return (om, _Sig, sinf, Edc, ntot)
        return ntot
                
    

def printm(a):
    """ Prints matrix in readable form"""
    for i in range(shape(a)[0]):
        for j in range(shape(a)[1]):
            print("%8.4f " % a[i,j].real, end=' ')
        print()
        
def ferm(x):
    """Fermi function for T=1"""
    if (x>100): return 0
    if (x<-100): return 1
    return 1/(exp(x)+1)

def Sinftyv(Sigind, CF, U, J, l, baths, Ns, nf, da_, fh_info):
    """ Computes Sigma_oo using OCA occupancies
        (density matrix given as a vector of numbers)
        Input:
            l             -- quantum number l
            baths         -- number of nonequivalent baths for OCA solver
            Ns[baths]     -- degeneracy of each bath
            nf[baths]     -- occupancy of each bath
            U, J          -- Coulomb interaction
        Output:
            sig_oo[baths] -- hartree-fock value of self-energy using impurity occupancies
            Edc[baths]    -- double counting using impurity occupancies
    """
    # occupation vector from OCA solver
    # OCA gives only few numbers and we need to create vector - matrix with this
    size = 2*(2*l+1)
    if len(Sigind)<size:
        SigindN = zeros((size,size), dtype=int)
        SigindN[:size/2,:size/2] = Sigind
        SigindN[size/2:,size/2:] = Sigind
        Sigind = SigindN
        CFN = zeros((size,size), dtype=complex)
        CFN[:size/2,:size/2] = CF
        CFN[size/2:,size/2:] = CF
        CF = CFN

    T2C = transpose(CF)

    mx = max(list(map(max,Sigind)))
    #print 'mx=', mx
    deg = zeros(mx, dtype=float)
    for i in range(len(Sigind)):
        for j in range(len(Sigind)):
            if Sigind[i,j]>0:
                deg[Sigind[i,j]-1]+=1

    occj = zeros((size,size), dtype=complex, order='F')
    for i in range(len(Sigind)):
        for j in range(len(Sigind)):
            if Sigind[i,j]>0:
                occj[i,j] = nf[Sigind[i,j]-1]/deg[Sigind[i,j]-1]

    #print 'nf=', nf
    #print 'occj='
    #print real(occj)
    
    # occupation - From relativistic -> spheric
    la.ztransform(occj, T2C,'C') 

    #print 'Sigind=', Sigind
    #print 'CF=', CF
    #print 'deg=', deg
    #print 'occjp=', occj

    # correlated index - all are correlated
    corind=ones(2*(2*l+1),order='F',dtype=int)

    # index (atom,l,m,s) in spheric
    bndind=[]
    for s in range(1,3):
        for m in range(-l,l+1):
            bndind.append([1,l,m,s])
    bndind = array(bndind,order='F',dtype=int)
    
    # The code for Sigma_oo works in spheric harmonics
    Uc = ones(2*(2*l+1))*U
    Jc = ones(2*(2*l+1))*J
    (sinfty, Edc) = ldau.hartree(Uc,Jc,occj,bndind,corind,da=da_, subtractdc=0)

    #print 'Uc=', Uc, 'Jc=', Jc
    
    la.ztransform(sinfty, T2C,'N') # From spheric -> relativistic

    #print 'sinfty='
    #for i in range(shape(T2C)[0]):
    #    for j in range(shape(T2C)[1]):
    #        print "%6.3f " % sinfty[i,j].real,
    #    print
    #
    #print 'Edc=', Edc
    
    
    # matrix of Sigma(infinity) is collapsed into few numbers
    # just like OCA occupation was given
    sinftyv=zeros(len(deg), dtype=float)
    Edcv=zeros(len(deg), dtype=float)
    for b1 in range(len(Sigind)):
        for b2 in range(len(Sigind)):
            if Sigind[b1,b2]>0:
                sinftyv[Sigind[b1,b2]-1] += sinfty[b1,b2].real
        if Sigind[b1,b1]>0:
            Edcv[Sigind[b1,b1]-1] += Edc[b1]
        
    for b in range(len(sinftyv)):
        sinftyv[b] /= deg[b]
        Edcv[b] /= deg[b]
        
    print('sinftyv=', sinftyv, file=fh_info)
    print('Edcv=', Edcv, file=fh_info)
    return (array(sinftyv), array(Edcv))


    
def trycutAf(L, om, Af, dh, omdh, fedh, nf_exact, Th, fh_info):
    m0 = 0
    m1 = 0
    nc = 0
    for i in range(len(om)):
        if (L<0): wcut = ferm(-(om[i]-L)/Th)
        else : wcut = ferm((om[i]-L)/Th)
        
        m0 += Af[i]*dh[i]*wcut
        m1 += Af[i]*omdh[i]*wcut
        nc += Af[i]*fedh[i]*wcut
    miss_dp = nf_exact-nc
    a1 = miss_dp
    a2 = 1-m0-a1
    return (a1, a2, m1)

def cutAf(L, om, Af, dh, omdh, fedh, nf_exact, Th, fh_info):
    m0 = 0
    m1 = 0
    nc = 0
    for i in range(len(om)):
        if (L<0): wcut = ferm(-(om[i]-L)/Th)
        else : wcut = ferm((om[i]-L)/Th)
        Af[i] *= wcut
        m0 += Af[i]*dh[i]
        m1 += Af[i]*omdh[i]
        nc += Af[i]*fedh[i]
    miss_dp = nf_exact-nc
    a1 = miss_dp
    a2 = 1-m0-a1
    return (a1, a2, m0, m1, nc)


def Soluvw(u, v, w):
    """ Given positive numbers u and v, gives one solution to the equation x*u + y*v = w
    where x and y are numbers in the interval [0,1].
    If there is no solution, gives False, if solution exists, gives the 'average' solution
    x = (x_min+x_max)/2 or y = (y_min+y_max)/2.
    """
    if (w<0 or w>u+v): return (False, 0.0, 0.0)  # Solution does not exists because u,v,x,y are postive
    if (w==0): return (True, 0.0, 0.0)
    
    x_max = (w<u and w/u or 1.0)  # min(w/u, 1)
    y_max = (w<v and w/v or 1.0)  # min(w/v, 1)
    
    x_min = (w-v)>0 and (w-v)<u  and (w-v)/u or 0.0 # max(0,(w-v)/u)
    y_min = (w-u)>0 and (w-u)<v  and (w-u)/v or 0.0 # max(0,(w-u)/v)
    
    if v>u :
        x = 0.5*(x_max+x_min)
        y = (w-x*u)/v
    else:
        y = 0.5*(y_max+y_min)
        x = (w-y*v)/u
    return (True, x, y)

def uvw(m0, m1, nf_exact, nc, dsinf, epsilon):
    a1 = nf_exact-nc
    if (a1<0): a1=0
    a2 = 1-m0-a1
    if (a2<0): a2=0
    
    w = dsinf-m1-a1*epsilon[0][0]-a2*epsilon[1][0]
    u = a1*(epsilon[0][1]-epsilon[0][0])
    v = a2*(epsilon[1][1]-epsilon[1][0])
    return (a1, a2, u, v, w)

   
def ww(L, Th, om, Af, dh, omdh, fedh, nf_exact, epsilon, dsinf, fh_info):
    m0 = 0
    m1 = 0
    nc = 0
    for i in range(len(om)):
        wcut = ferm(-(L+om[i])/Th)*ferm((om[i]-L)/Th)
        m0 += Af[i]*dh[i]*wcut
        m1 += Af[i]*omdh[i]*wcut
        nc += Af[i]*fedh[i]*wcut

    (a1, a2, u, v, w) = uvw(m0, m1, nf_exact, nc, dsinf, epsilon)
    (succ, x, y) = Soluvw(u, v, w)
        
    print('# L=%9.6f u=%9.6f v=%9.6f w=%9.6f x=%9.6f  y=%9.6f m0=%9.6f m1=%9.6f a1=%9.6f a2=%9.6f dsinf=%9.6f' % (L, u, v, w, x, y, m0, m1, a1, a2, dsinf), file=fh_info)
    return (succ, a1, a2, x, y)

def ww_cut(L, Th, om, Af, dh, omdh, fedh, nf_exact, epsilon, dsinf, fh_info):
    m0 = 0
    m1 = 0
    nc = 0
    for i in range(len(om)):
        wcut = ferm(-(L+om[i])/Th)*ferm((om[i]-L)/Th)
        Af[i] *= wcut        
        m0 += Af[i]*dh[i]
        m1 += Af[i]*omdh[i]
        nc += Af[i]*fedh[i]

    (a1, a2, u, v, w) = uvw(m0, m1, nf_exact, nc, dsinf, epsilon)
    (succ, x, y) = Soluvw(u, v, w)
        
    print('# L=%9.6f u=%9.6f v=%9.6f w=%9.6f x=%9.6f  y=%9.6f m0=%9.6f m1=%9.6f a1=%9.6f a2=%9.6f dsinf=%9.6f' % (L, u, v, w, x, y, m0, m1, a1, a2, dsinf), file=fh_info)
    return (succ, a1, a2, x, y)







if __name__ == '__main__':

    Znuc=58
    dire='.'

    iparams0={"exe"            : ["nca"               , "# Name of the executable"],
          "U"                  : [5.5                 , "# Coulomb repulsion (F0)"],
          "J"                  : [0.60                , "# J"],
          "T"                  : [0.019               , "# Temperature"],
          "nf0"                : [1.0                 , "# Double counting parameter"],
          "nc"                 : [[0, 1, 2]           , "# Impurity occupancies"],
          "Ncentral"           : [[1]                 , "# Central occupancies for OCA diagrams evaluation"],
          "alpha"              : [0.5                 , "# Mixing for bath spectral function"],
          "max_steps"          : [20                  , "# Maximum number of impurity steps"],
          "max_diff"           : [0.001               , "# Maximum difference between steps"],
          "followPeak"         : [-1                  , "# A mode to determin lambda0"],
          "Q"                  : [8.0                 , "# A parameter to determin lambda0"],
          "StartLambda"        : [-20.0               , "# Where to start looking for zero to determin lambda0"],
          "dLambda"            : [0.02                 , "# Step in searching for the lambda"],
          "EndLambda"          : [1.0                 , "# Where to stop searching for the lambda0"],
          "cutAc"              : [[-10.0, 12.0]       , "# Only window [La,Lb] of baths spectra is taken into account"],
          "Gh"                 : [0.5                 , "# Parameter to improve the high frequency self-energy"],
          "epsilon"            : [[[-10.0, -3.], [3., 10.0]], "# Parameter to improve the high frequency self-energy"],
          "Th"                 : [0.5                 , "# Parameter to improve the high frequency self-energy"],
          "lorentz"            : [1                   , "# Weather to subtract lorentz from diverging spectral functions and treat it analytically"],
          "SearchLorentz"      : [2.5                 , "# How far from zero to search fro Lorentz"],
          "LorentzMaxRatio"    : [1.0                 , "# How far from zero to search for Lorentz"],
          "FirstLorentz"       : [0                   , "# First pseudoparticle which could be augmented with lorentz"],
          "LastLorentz"        : [10000               , "# Last pseudoparticle which could be augmented with lorentz"],
          "CmpDiff"            : [-1                  , "# When calculating Difference, only first CmpDiff particles should be taken into account (-1,all)"],
          "correct"            : ['n'                 , "# Do not correct anything at infinity"],
          "Eimp0"              : [-2.0                , "# First impurity level"],
    }

    s_oo= [4.3410663935612748, 5.1905227739208444]
    Edc= [2.749999985125807, 2.749999985125807]
    #Ed=[-1.9206795651000002, -1.5932864251000001]  	# Impurity levels

    Ed=[-1.87103536, -1.5444438300000001]    
    #Ed=[-1.87103536, -1.5444438300000001]
    nf0 = iparams0['nf0'][0]
    Jc = iparams0['J'][0]
    
    iparams0['Ed']=[Ed, '# ']



    def ReadTrans(filename):
        """Read the self-energy index file Sigind and the local transformation matrix CF from a file"""
        fh = open(filename, 'r')
        (n1,n2) = list(map(int, next(fh).split()[:2]))
        next(fh) # comment
        
        Sigind=[]
        for i in range(n1):
            Sigind.append( list(map(int, next(fh).split()[:n2])) )
        Sigind = array(Sigind)
        
        next(fh) # comment
        
        CF=[]
        for i in range(n1):
            cl = array(list(map(double, next(fh).split()[:2*n2])))
            CF.append( cl[0::2]+cl[1::2]*1j )
        CF = array(CF)
        
        return (Sigind, CF)

    (Sigind, CF) = ReadTrans('imp.0/Trans.dat')
    
    solver = IMP_OCA(Znuc, 'imp.0', iparams0, Sigind, CF)
    
    UpdateAtom=False
    gbroad=0.01
    
    DCs = "fixEimp"

    #solver._exe(iparams0, DCs, 'xx', UpdateAtom, gbroad)
    da=0
    extn='xx'
    solver.HighFrequency(iparams0, DCs, da, extn)
    
