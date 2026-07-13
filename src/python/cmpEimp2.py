#!/usr/bin/env python
# @Copyright 2007 Kristjan Haule
from numpy import *
import sys,os,re
import ast
import json
from scipy import interpolate
import builtins

# All variables set through exec call need to be global for Python3 to set them. Hence we are forced to work with global variables.
cix='actqmc.cix'
Delta='Delta.inp'
Gf='Gf.out'
Sig='Sig.out'
#beta,mu=0,0
#nf,mu,Ekin,Epot,TrSigmaG,mom=0,0,0,0,0,0

def ReadPARAMS2(fparams, fout):
    with open(fparams, 'r') as fp:
        lines = fp.readlines()

    stripped = []
    for line in lines:
        n = line.find('#')
        if n>=0:
            line = line[:n]
        stripped.append(line)
    content = ''.join(stripped).strip()

    if content.startswith('{'):
        params = json.loads(content)
        mapping = {'fSig': 'Sig', 'fGf': 'Gf', 'Delta': 'Delta', 'cix': 'cix', 'beta': 'beta', 'mu': 'mu'}
        for key_in, key_out in mapping.items():
            if key_in in params:
                globals()[key_out] = params[key_in]
                print('executed:', key_out+'='+repr(params[key_in]), file=fout)
    else:
        for line in stripped:
            dat = line.split()
            if len(dat)>=2:
                m = re.search('cix|Delta|Gf|Sig', dat[0])
                if m is not None :
                    cmn = dat[0]+'="'+dat[1]+'"'
                    exec(cmn, globals())
                    print('executed:', cmn, file=fout)
                m = re.search('beta|mu',dat[0])
                if m is not None :
                    cmn = dat[0]+'='+dat[1]
                    exec(cmn, globals())
                    print('executed:', cmn, file=fout)
    print('beta=', beta, 'mu=', mu, 'cix=', cix, 'Delta=', Delta, file=fout)
    return (cix,Delta,Gf,Sig,beta,mu)

def ReadCix(cix, fout):
    fcx = open(cix,'r')
    next(fcx);next(fcx)
    header = next(fcx).split()
    Deg={}

    if len(header) == 4:
        (Nclust,Nm,Ns,maxm) = list(map(int, header))
        next(fcx)
        for i in range(Ns):
            dat = list(map(int,next(fcx).split()))
            dim = dat[1]
            if (dim>1): 
                print('WARNING: The impurity functional cmpEimp2.py works for now only when baths are diagonal. You should not trust total free energy. Needs some development...!')
            for j in range(dim):
                for k in range(dim):
                    if k==j:
                        ii = dat[2+k+j*dim]
                        if ii in Deg:
                            Deg[ii] += 1
                        else:
                            Deg[ii]=1
        next(fcx)
        ss = next(fcx)
    else:
        # New Bold-CTQMC cix format:
        # cluster_size, number of states, number of psi^+ operators,
        # number of baths, number of blocks, maximum_matrix_size
        Ns = int(header[3])
        next(fcx)
        for i in range(Ns):
            dat = next(fcx).split()
            ii = int(dat[3])
            if ii in Deg:
                Deg[ii] += 1
            else:
                Deg[ii] = 1
        next(fcx)
        ss = next(fcx)

    print('Deg=', Deg, file=fout)
    if ss[:11]=='FL_FROM_IFL':
        for i in range(Ns):
            next(fcx)
        next(fcx)
        ss = next(fcx)
    Eks = list(map(float,ss.split()))
    print('Eks=', Eks, file=fout)
    if len(Deg)>0:
        deg = zeros(builtins.max(list(Deg.keys())) + 1, dtype=int)
        for ii,val in Deg.items():
            deg[ii] = val
    else:
        deg = zeros(0, dtype=int)
    return (deg, array(Eks))


def ReadGandSigandDelta(Gf, Sig, Delta, fout):
    Delta_data = loadtxt(Delta).transpose()
    omD = Delta_data[0]
    Dlta = Delta_data[1::2]+Delta_data[2::2]*1j
    
    Gf_data = loadtxt(Gf).transpose()
    om = Gf_data[0]
    Gfc = Gf_data[1::2]+Gf_data[2::2]*1j

    Sig_data = loadtxt(Sig).transpose()
    om3 = Sig_data[0]
    Sigc = Sig_data[1::2]+Sig_data[2::2]*1j
    
    print('shapes=', shape(om), shape(omD), shape(om3), file=fout)
    #print 'Diff=', sum(abs(om-omD))
    #if sum(abs(om-omD))>1. or sum(abs(om-om3))>1.:
    #    print >> fout, 'frequency mesh of ',Gf,' and ',Sig, 'and', Delta,' are not equal', om-om2, om-om3
    
    first_line = open(Gf, 'r').readline()
    dat = first_line.split()
    for d in dat:
        m = re.search('nf|mu|T|TrSigmaG|Ekin|Epot|mom', d)
        if m is not None:
            print('executed:', d, file=fout)
            exec(d, globals())
            
    print('nf=', nf, 'mu=', mu, 'Ekin=', Ekin, 'Epot=', Epot, file=fout)
    return (om, omD, Gfc, Sigc, Dlta, nf, mu, Ekin, Epot, TrSigmaG, mom) 


def ferm(x):
    if x>300:
        return 0
    elif x<-300:
        return 1
    else:
        return 1/(exp(x)+1)

def GetHighFrequency(CC,om):
    " Approximates CC ~  A/(i*om-C) "
    A = 1./imag(1/(CC[-1]*om[-1]))
    C = -A*real(1./CC[-1])
    return (A, C)

def ReturnHighFrequency(A,C,beta,E):
    " Returns the value of Tr(A/(iom-C) 1/(iom-E) )"
    return A*(ferm(E*beta)-ferm(C*beta))/(E-C)

def LogGimp(omega,Gfc,EimpS,beta,Deg,fout):
    " Tr(log(-G_imp))"
    def NonIntF0(beta,Ene):
        if beta*Ene>200: return 0.0
        if beta*Ene<-200: return Ene
        return -log(1+exp(-Ene*beta))/beta
    
    lnGimp=0.
    for i in range(len(Deg)): 
        print('EimpS=', EimpS[i], ', real(iom-1/G)', real(-1/Gfc[i][-1]), file=fout)
        if EimpS[i] < 1000 :
            Gd = Gfc[i,:]
            eimps = EimpS[i]
            
            CC = omega*1j-eimps-1/Gd
            A,C = GetHighFrequency(CC,omega)

            ff = log(-Gd)-log(-1./(omega*1j-eimps)) - 1./(omega*1j-eimps)*A/(omega*1j-C)
            lnGimp_i = 2*sum(ff.real)/beta+NonIntF0(beta,eimps)+ReturnHighFrequency(A,C,beta,eimps)
            lnGimp_i *= Deg[i]
            
        else:
            sum_imp = 0.0
            Fnint_imp = 0.0
            lnGimp_i = 0.0
        lnGimp += lnGimp_i
        print('lnGimp: Deg[i]=', Deg[i], 'Tr(log(-Gimp))=', lnGimp_i/Deg[i], 'Tr(log(-Gimp))*deg=', lnGimp_i, file=fout)
    return lnGimp


def LogGimpSig(omega, Gfc, Sigc, Vdc, EimpS, beta, Deg, fout):
    " Tr(log(1+(Sig-Vdc)G))"
    def NonIntF0(beta,Ene):
        if beta*Ene>200: return 0.0
        if beta*Ene<-200: return Ene
        return -log(1+exp(-Ene*beta))/beta
    
    lnGSimp=0.
    for i in range(len(Deg)):#len(Gfc)):
        Gd = Gfc[i,:]
        Sg = Sigc[i,:]
        eimps = EimpS[i]
        s_oo = Sigc[i,-1].real
        vdc = Vdc[i]
        
        #print 'shape(Gd)=', shape(Gd)
        ff0 = log((Sg-vdc)*Gd+1.)
        ff1 = (s_oo-vdc)/(omega*1j-eimps)
        
        lnGimp_i = 2*sum((ff0-ff1).real)/beta + (s_oo-vdc)*ferm(eimps*beta)
        lnGimp_i *= Deg[i]

        lnGSimp += lnGimp_i
        print('lnGSimp: Deg[i]=', Deg[i], 'Tr(log(1+Sig*Gimp))=', lnGimp_i/Deg[i], 'Tr(log(1+Sig*Gimp))*deg=', lnGimp_i, file=fout)

        #plot(omega, real(ff0-ff1), label='diff-real')
        #plot(omega, imag(ff0-ff1), label='diff-imag')
        #legend(loc='best')
        #show()
    return lnGSimp
    
def fTrDeltaG(om, omD, Dlta, Gfc, Deg, EimpS, beta, fout):
    def GetInterpolation(omD, Delta):
        "Returns derivative  d Delta/d om_i"
        Dr = interpolate.UnivariateSpline(omD, Delta.real, s=0,k=3)
        Di = interpolate.UnivariateSpline(omD, Delta.imag, s=0,k=3)
        return (Dr, Di)
    lngh = builtins.min(len(Dlta),len(Gfc),len(Deg))
    GDf=0.0
    for i in range(lngh):
        if len(omD)!=len(om):
            # interpolating Delta
            Dr, Di = GetInterpolation(omD, Dlta[i])
            Delta = Dr(om)+Di(om)*1j
        else:
            Delta = Dlta[i]
        (A,C) = GetHighFrequency(Delta,om)
        print('A_d=', A, 'C_d=', C, file=fout)
        eimps = EimpS[i]
        GD1 = Gfc[i]*Delta - 1./(om*1j-eimps)*A/(om*1j-C)
        dGD = Deg[i]*(2*sum(GD1.real)/beta + ReturnHighFrequency(A,C,beta,eimps))
        GDf += dGD
        #print >> fout, 'dTr(Delta*G)=', dGD/Deg[i]
    return GDf

def fTrSigmaG(om, Gfc, Sigc, Deg, EimpS, beta, fout):
    SDf=0.0
    for i in range(len(Deg)):
        if EimpS[i] < 1000 :
            Gd = Gfc[i,:]
            eimps = EimpS[i]
            Sg = Sigc[i,:]
            s_oo = Sigc[i,-1].real
            Sg0 = Sg[:] - s_oo   # Sigma-s_oo
            (A,C) = GetHighFrequency(Sg0,om)
            #print '::: A,C, for Sigma=', A, C
            C=1.
            ff = Gd[:]*Sg0[:] - 1./(om*1j-eimps)*A/(om*1j-C)
            SDf_i = 2*sum(ff.real)/beta + ReturnHighFrequency(A,C,beta,eimps)
            SDf_i *= Deg[i]
        else:
            SDf_i = 0.0
        SDf += SDf_i
        #print >> fout, 'Tr((Sigma-s_oo)*G): Deg[i]=', Deg[i], 'Tr((Sigma-s_oo)*G)=', SDf_i/Deg[i], 'Tr((Sigma-s_oo)*G)*deg=', SDf_i
    return SDf
            
def CmpImpurityCorrection2(om, omD, Dlt, Gfc, Deg, EimpS, beta, fout):
    def GetInterpolation(omD, Delta):
        "Returns derivative  d Delta/d om_i"
        Dr = interpolate.UnivariateSpline(omD, Delta.real, s=0,k=3)
        Di = interpolate.UnivariateSpline(omD, Delta.imag, s=0,k=3)
        return (Dr, Di)

    Delta_inf=[]
    for b in range(len(Dlt)):
        Delta_inf.append( [Dlt[b,-1].real*omD[-1]**2, Dlt[b,-1].imag*omD[-1] ] )
    
    cxx = 0.0
    lngh = builtins.min(len(Dlt),len(Gfc),len(Deg))
    
    oms = om
    Gfs = Gfc
    if oms[0]<omD[0] and abs(oms[0]-omD[0])<1e-6: oms[0]=omD[0]
    
    print('oms[-1]=', oms[-1], 'omD[-1]=', omD[-1], file=fout)
    
    for i in range(lngh):
        #print('shape(omD)=', shape(omD), 'shape(Dlt[i])=', shape(Dlt[i]), 'len(oms)=', len(oms))
        Dr, Di = GetInterpolation(omD, Dlt[i])
        dDlt=zeros(len(oms),dtype=complex)
        for ix,x in enumerate(oms):
            if omD[0]<=x<=omD[-1]:
                dDlt[ix] = Dr.derivatives(x)[1]+Di.derivatives(x)[1]*1j
            else:
                dDlt[ix] = -2*Delta_inf[i][0]/x**3 - Delta_inf[i][1]/x**2 * 1j
        
        CC = oms*dDlt
        #savetxt( 'brisx', vstack((oms, real(CC), imag(CC))).transpose() )
        #print 'shape(dDlt)=', shape(dDlt), 'shape(CC)=', shape(CC), 'shape(Gfs)=', shape(Gfs)
        (A,C) = GetHighFrequency(CC,oms)
        print('A=', A, 'C=', C, file=fout)
        eimps = EimpS[i]
        GD1 = Gfs[i]*CC - 1./(oms*1j-eimps)*A/(oms*1j-C)
        dcxx = Deg[i]*(2*sum(GD1.real)/beta + ReturnHighFrequency(A,C,beta,eimps))
        cxx += dcxx
        print('dcxx=', dcxx/Deg[i], file=fout)
        
    print('cxx=', cxx, file=fout)
    return cxx

def GiveImpFunctional(dire, fparams, Vdc, fout=sys.stdout):

    if not (os.path.exists(dire+'/'+fparams) and os.path.getsize(dire+'/'+fparams)>0):
        print('Give input PARAMS file', file=fout)
        sys.exit(0)
    
    (cix,Delta,Gf,Sig,beta,mu_qmc) = ReadPARAMS2(dire+'/'+fparams, fout)
    (Deg, Eks) = ReadCix(dire+'/'+cix, fout)
    
    
    (om, omD, Gfc, Sigc, Dlta, nf, mu, Ekin, Epot, TrSigmaG, mom) = ReadGandSigandDelta(dire+'/'+Gf, dire+'/'+Sig, dire+'/'+Delta, fout)
    s_oo = Sigc[:,-1].real
    print('s_oo=', s_oo, file=fout)
    
    lengh = builtins.min(len(Eks),len(s_oo))
    EimpS = Eks[:lengh]-mu_qmc+s_oo[:lengh]
    # Tr(log(-Gimp))
    lnGimp = LogGimp(om,Gfc,EimpS,beta,Deg,fout) 

    #lnGS = LogGimpSig(om, Gfc, Sigc, Vdc, EimpS, beta, Deg, fout)

    trsigg =  fTrSigmaG(om, Gfc, Sigc, Deg, EimpS, beta, fout) # numericall Tr((Sigma-s_oo)*G)
    trsigg_oo = sum(mom[:]*s_oo[:])

    # Tr(Delta*G)
    TrDeltaG = fTrDeltaG(om, omD, Dlta, Gfc, Deg, EimpS, beta, fout)
    print('TrDeltaG=', TrDeltaG, '<k>/T=', Ekin, file=fout)
    # Tr(om*dDelta/dom)
    cxx = CmpImpurityCorrection2(om, omD, Dlta, Gfc, Deg, EimpS, beta, fout)

    # Eimp = Tr(Delta*G)+1/2*Tr(Sigma*G)+eimp*nimp-Tr(om*dDelta/dom)
    #Eimp = Ekin+Epot-cxx-mu*nf                # Impurity Energy
    Eimp = TrDeltaG+Epot-cxx-mu*nf                # Impurity Energy

    # Potthof2[G] == Phi[G]-1/2*Tr(Sigma*G)+T*S_imp = Eimp-Tr(log(-G_imp))+1/2*Tr(Sigma*G)
    #Potthof2 = Eimp-lnGimp+TrSigmaG            # Phi[G]-1/2*Tr(Sigma*G)+T*S_imp

    #nd_numeric=zeros(len(Deg))
    #for i in range(len(Deg)):
    #    Gd = Gfc[i,:]
    #    eimps = EimpS[i]
    #    ff = Gd[:]-1./(om*1j-eimps)
    #    nd_numeric[i] = (2*sum(ff.real)/beta + ferm(eimps*beta))*Deg[i]
    #print 'nd_numeric=', nd_numeric
    
    #Eimp1 = lnGS - lnGimp + Eimp  # Tr(log(1+Sigma*G))-Tr(log(-Gloc))+E_imp

    logZatom=0
    filename=dire+'/ctqmc.log'
    if os.path.isfile(filename):
        fw=open(filename,'r')
        lines=fw.readlines()
        for line in lines:
            m = re.search('log\(Zatom\)=(\d+\.\d+)',line)
            if m is not None:
                logZatom = float(m.groups()[0])
        fw.close()
        #print 'logZatom=', logZatom
    P0=0
    filename=dire+'/histogram.dat'
    if os.path.isfile(filename):
        fw=open(filename,'r')
        next(fw)# comment
        P0 = float(next(fw).split()[1])
        #print 'P0=', P0

    Eimp2 = TrSigmaG   # 1/2*Tr(Sigma*G)
    print('lnGimp=', lnGimp, file=fout)
    print("%-12s "*6 % ('Eimp', 'Ekin', 'Epot', 'cxx', 'mu*nf', 'lnGimp'), file=fout)
    print("%-12.6f "*6 % (Eimp, TrDeltaG, Epot, cxx, mu*nf, lnGimp), file=fout)
    print('#',"%-12.6f "*6 % (Eimp+Ekin-TrDeltaG, Ekin, Epot, cxx, mu*nf, lnGimp), file=fout)
    print('Eimp=', Eimp, file=fout)
    if logZatom and P0>1e-7:
        logZ = logZatom-log(P0)
        Fimp = -logZ/beta
        Entropy = (Eimp-Fimp)*beta
        WARN=''
        if (P0<1e-5): WARN='   ** WARNING P0 is very small, hence Entropy is not accurate'
        print('Entropy=', Entropy, 'Fimp=', Fimp, 'logZatom=', logZatom, 'P0=', P0, WARN, file=fout)
    #print >> fout, 'P[Sigma]=', Potthof2

    print('1/2*Tr(Sigma*G)=', 0.5*(trsigg+trsigg_oo), 'deimp*nf=', sum(Eks*mom), file=fout)

    Vdc_nd = sum(array(mom[:])*array(Vdc[:]))
    zeorb = TrDeltaG+Epot-cxx - sum(Eks*mom) # note that Epot ~ 1/2*Tr(Sigma*G)+Eks*mom
    print('Vdc*nf=', Vdc_nd, file=fout)
    Phi_DMFT = -lnGimp+Eimp+(trsigg+trsigg_oo)
    print('Phi_DMFT=', Phi_DMFT, file=fout)
    print('Phi_DMFT - Vdc*nf=', Phi_DMFT-Vdc_nd, file=fout)
    Ry2eV = 13.60569193
    print('Phi_DMFT - Vdc*nf=', (Phi_DMFT-Vdc_nd)/Ry2eV, 'Ry', file=fout)
    print('Epot - Vdc*nf=', (TrSigmaG-Vdc_nd)/Ry2eV, 'Ry', file=fout)
    print('zeorb[eV]=', zeorb, file=fout)
    return (Phi_DMFT, lnGimp, Eimp, Vdc_nd, zeorb)


if __name__ == '__main__':
    import re
    
    if len(sys.argv)<2:
        dire = '.'
    else:
        dire = sys.argv[1]
    Ry2eV = 13.60569193

    fi = open(dire+'/Eimp.inp','r')
    line = next(fi)  # Eimp
    line = next(fi)  # Overlap
    line = next(fi)  # Vdc
    Vdc = list(map(float,re.sub(r'#.*',' ',line).split()))

    (Phi_DMFT, lnGimp, Fimp, Vdc_nd, zeorb) = GiveImpFunctional(dire, 'PARAMS', Vdc)
    print('ev: Phi_DMFT-Vdc*nd=', Phi_DMFT-Vdc_nd, 'Fimp+TS=', Fimp, 'logGimp=', lnGimp, 'zeorb=', zeorb)
    print('Ry: Phi_DMFT-Vdc*nd=', (Phi_DMFT-Vdc_nd)/Ry2eV, 'FImp+TS=', Fimp/Ry2eV, 'logGimp=', lnGimp/Ry2eV, 'zeorb=', zeorb/Ry2eV)
    


# ===== BHQMC (bold-CTQMC) free-energy helpers ported from fork (uses ctqmc.log) =====
def obtain_last_value(valuename, log_path):
    """
    Reads the last occurrence of `valuename` from ctqmc.log.

    Examples
    --------
    obtain_last_value('Epot =', log_path) -> 80.550044989283
    obtain_last_value('bath_density_matrix=', log_path) -> [0.513..., ...]
    """
    with open(log_path, 'r') as fh:
        lines = fh.readlines()

    for line in reversed(lines):
        pos = line.find(valuename)
        if pos < 0:
            continue

        tail = line[pos + len(valuename):].strip()
        if not tail:
            return None

        # Array-like values, e.g. [1.0,2.0,...]
        if tail[0] == '[':
            end = tail.find(']')
            if end < 0:
                raise ValueError(f"Malformed list for {valuename} in {log_path}")
            return ast.literal_eval(tail[:end+1])

        # Scalar numeric value
        m = re.search(r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?', tail)
        if m is None:
            raise ValueError(f"Could not parse value for {valuename} from line:\n{line}")
        return float(m.group(0))

    raise ValueError(f"Could not find {valuename} in {log_path}")

def ReadCix_bhqmc(cix, fout):
    with open(cix, 'r') as fcx:
        next(fcx)  # # Cix file for cluster DMFT with Bold-CTQMC
        next(fcx)  # # cluster_size, number of states, ...

        (Nclust, Nstates, Npsis, Ns, Nblocks, Mblocksize) = list(map(int, next(fcx).split()))
        print('header=', Nclust, Nstates, Npsis, Ns, Nblocks, Mblocksize, file=fout)

        next(fcx)  # # c^+_i,c_i, Sigind  spin   moment

        Deg = {}
        for i in range(Ns):
            dat = next(fcx).split()
            if len(dat) < 4:
                continue

            ii = int(dat[3])

            if ii in Deg:
                Deg[ii] += 1
            else:
                Deg[ii] = 1

        print('Deg=', Deg, file=fout)

        next(fcx)  # # cluster energies for non-equivalent baths, eps[k]
        ss = next(fcx).strip()

        # for FL_FROM_IFL: not done yet. currently we do not need this.
        # if ss[:11] == 'FL_FROM_IFL':
        #     # TODO:
        #     # for i in range(...):
        #     #     next(fcx)
        #     # next(fcx)
        #     # ss = next(fcx).strip()
        #     pass

        Eks = list(map(float, ss.split()))
        print('Eks=', Eks, file=fout)

        return (Deg, array(Eks))

def ReadGandSigandDelta_bhqmc(Gf, Sig, Delta,boldlog, Eks,Sigind,fout):
    Delta_data = loadtxt(Delta).transpose()
    omD = Delta_data[0]
    Dlta = Delta_data[1::2]+Delta_data[2::2]*1j
    
    Gf_data = loadtxt(Gf).transpose()
    om = Gf_data[0]
    Gfc = Gf_data[1::2]+Gf_data[2::2]*1j

    Sig_data = loadtxt(Sig).transpose()
    om3 = Sig_data[0]
    Sigc = Sig_data[1::2]+Sig_data[2::2]*1j
    
    print('shapes=', shape(om), shape(omD), shape(om3), file=fout)
    #print 'Diff=', sum(abs(om-omD))
    #if sum(abs(om-omD))>1. or sum(abs(om-om3))>1.:
    #    print >> fout, 'frequency mesh of ',Gf,' and ',Sig, 'and', Delta,' are not equal', om-om2, om-om3
    
    # in bhqmc they can be obtained from ctqmc.log
    Epot=obtain_last_value('Epot = ',boldlog)
    mu=obtain_last_value('"mu": ',boldlog)


    # DenMat=obtain_last_value('bath_density_matrix=',boldlog)# density matrix
    # nf= sum(DenMat)
    # # mom is the total occupation on each nonquivalent orbitals.
    # mom=zeros(max(diag(Sigind)))
    # for iorb in arange(len(Sigind[0])):
    #     mom[Sigind[iorb][iorb]-1]+=DenMat[iorb]#up
    #     mom[Sigind[iorb][iorb]-1]+=DenMat[iorb+len(Sigind[0])]#down
    # TrSigmaG=Epot-sum(mom*Eks)


    # nf and the per-component occupations `mom` from the Gf.out Print_Physical2
    # header.  The current bctqmc binary already folds the cluster occupations into
    # mom = [Gamma, X=Y, M] (one entry per independent Sigind component, sum(mom)=nf).
    # The old May-boldc binary instead wrote the raw per-spin bath_density_matrix
    # which this function then refolded; the Jun binary writes the reduced header,
    # so we read mom/nf straight from it (same convention as IMP_BCTQMC.HighFrequency
    # and the ctqmc parent).
    hdr = {}
    with open(Gf, 'r') as _fh:
        for par in _fh.readline().strip().split():
            if '=' in par:
                _k, _v = par.split('=', 1)
                if _k in ('nf', 'mom', 'TrSigmaG'):
                    try:
                        hdr[_k] = eval(_v)
                    except Exception:
                        pass
    nf = hdr['nf']
    mom = array(hdr['mom'])

    TrSigmaG = Epot - sum(mom * Eks)
        
    print('nf=', nf, 'mu=', mu,  'Epot=', Epot, file=fout)
    return (om, omD, Gfc, Sigc, Dlta, nf, mu,  Epot, TrSigmaG, mom)

def GiveImpFunctional_bctqmc(dire, fparams, Vdc,Sigind, fout=sys.stdout):

    if not (os.path.exists(dire+'/'+fparams) and os.path.getsize(dire+'/'+fparams)>0):
        print('Give input PARAMS file', file=fout)
        sys.exit(0)
    
    (cix,Delta,Gf,Sig,beta,mu_qmc) = ReadPARAMS2(dire+'/'+fparams, fout)

    (Deg, Eks) = ReadCix_bhqmc(dire+'/'+cix, fout)
    # note: Ekin in ctqmc is estimated from average perturbation order... there is no such thing in bhqmc i think.
    (om, omD, Gfc, Sigc, Dlta, nf, mu, Epot, TrSigmaG, mom) = ReadGandSigandDelta_bhqmc(dire+'/'+Gf, dire+'/'+Sig, dire+'/'+Delta,dire+'/'+'ctqmc.log', Eks,Sigind,fout)

    s_oo = Sigc[:,-1].real
    print('s_oo=', s_oo, file=fout)
    
    lengh = builtins.min(len(Eks),len(s_oo))
    EimpS = Eks[:lengh]-mu_qmc+s_oo[:lengh]
    # Tr(log(-Gimp))
    lnGimp = LogGimp(om,Gfc,EimpS,beta,Deg,fout) 

    #lnGS = LogGimpSig(om, Gfc, Sigc, Vdc, EimpS, beta, Deg, fout)

    trsigg =  fTrSigmaG(om, Gfc, Sigc, Deg, EimpS, beta, fout) # numericall Tr((Sigma-s_oo)*G)
    trsigg_oo = sum(mom[:]*s_oo[:])

    # Tr(Delta*G)
    TrDeltaG = fTrDeltaG(om, omD, Dlta, Gfc, Deg, EimpS, beta, fout)

    Ekin=TrDeltaG
    print('This is BHQMC. Ekin is set to be the same as TrDeltaG', file=fout)
    print('TrDeltaG=', TrDeltaG, '<k>/T=', Ekin, file=fout)
    # Tr(om*dDelta/dom)
    cxx = CmpImpurityCorrection2(om, omD, Dlta, Gfc, Deg, EimpS, beta, fout)

    # Eimp = Tr(Delta*G)+1/2*Tr(Sigma*G)+eimp*nimp-Tr(om*dDelta/dom)
    #Eimp = Ekin+Epot-cxx-mu*nf                # Impurity Energy
    Eimp = TrDeltaG+Epot-cxx-mu*nf                # Impurity Energy

    # Potthof2[G] == Phi[G]-1/2*Tr(Sigma*G)+T*S_imp = Eimp-Tr(log(-G_imp))+1/2*Tr(Sigma*G)
    #Potthof2 = Eimp-lnGimp+TrSigmaG            # Phi[G]-1/2*Tr(Sigma*G)+T*S_imp

    #nd_numeric=zeros(len(Deg))
    #for i in range(len(Deg)):
    #    Gd = Gfc[i,:]
    #    eimps = EimpS[i]
    #    ff = Gd[:]-1./(om*1j-eimps)
    #    nd_numeric[i] = (2*sum(ff.real)/beta + ferm(eimps*beta))*Deg[i]
    #print 'nd_numeric=', nd_numeric
    
    #Eimp1 = lnGS - lnGimp + Eimp  # Tr(log(1+Sigma*G))-Tr(log(-Gloc))+E_imp

    # in BHQMC, Fimp = lambda0 - T log Q. basically read from projection Q.
    Fimp=obtain_last_value('Fimp = ',dire +'ctqmc.log')
    Entropy = (Eimp-Fimp)*beta
    print('Entropy=', Entropy, 'Fimp=', Fimp, file=fout)
    Eimp2 = TrSigmaG   # 1/2*Tr(Sigma*G)
    print('lnGimp=', lnGimp, file=fout)
    print("%-12s "*6 % ('Eimp', 'Ekin', 'Epot', 'cxx', 'mu*nf', 'lnGimp'), file=fout)
    print("%-12.6f "*6 % (Eimp, TrDeltaG, Epot, cxx, mu*nf, lnGimp), file=fout)
    print('#',"%-12.6f "*6 % (Eimp+Ekin-TrDeltaG, Ekin, Epot, cxx, mu*nf, lnGimp), file=fout)
    print('Eimp=', Eimp, file=fout)


    print('1/2*Tr(Sigma*G)=', 0.5*(trsigg+trsigg_oo), 'deimp*nf=', sum(Eks*mom), file=fout)

    Vdc_nd = sum(array(mom[:])*array(Vdc[:]))
    zeorb = TrDeltaG+Epot-cxx - sum(Eks*mom) # note that Epot ~ 1/2*Tr(Sigma*G)+Eks*mom
    print('Vdc*nf=', Vdc_nd, file=fout)
    Phi_DMFT = -lnGimp+Eimp+(trsigg+trsigg_oo)
    print('Phi_DMFT=', Phi_DMFT, file=fout)
    print('Phi_DMFT - Vdc*nf=', Phi_DMFT-Vdc_nd, file=fout)
    Ry2eV = 13.60569193
    print('Phi_DMFT - Vdc*nf=', (Phi_DMFT-Vdc_nd)/Ry2eV, 'Ry', file=fout)
    print('Epot - Vdc*nf=', (TrSigmaG-Vdc_nd)/Ry2eV, 'Ry', file=fout)
    print('zeorb[eV]=', zeorb, file=fout)
    return (Phi_DMFT, lnGimp, Eimp, Vdc_nd, zeorb,TrSigmaG)
