#!/usr/bin/env python
# @Copyright 2007 Kristjan Haule
from numpy import *
import sys,os,re
from scipy import interpolate

cix='actqmc.cix'
Delta='Delta.inp'
Gf='Gf.out'
Sig='Sig.out'
#beta,mu=0,0
#nf,mu,Ekin,Epot,TrSigmaG,mom=0,0,0,0,0,0

def ReadPARAMS2(fparams, fout):
    fp = open(fparams, 'r')
    for line in fp:
        n = line.find('#')
        if n>=0:
            line = line[:n]
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
    return (cix,Delta,Gf,Sig,beta,mu)

def ReadCix(cix):
    fcx = open(cix,'r')
    next(fcx);next(fcx)
    (Nclust,Nm,Ns,maxm) = list(map(int,next(fcx).split()))
    next(fcx)
    Deg={}
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
    Eks = list(map(float,next(fcx).split()))
    return (Deg, array(Eks))


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
    #fg = open(Sig, 'r')
    #first_line = next(fg)
    #fg.close()
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
    for i in range(len(Gfc)):
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
    for i in range(len(Gfc)):
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
    

def CmpImpurityCorrection2(om, omD, Dlt, Gfc, Deg, EimpS, beta, fout):
    
    def GiveDerivative(om, omD, Delta):
        "Returns derivative  d Delta/d om_i"
        Dr = interpolate.UnivariateSpline(omD, Delta.real, s=0,k=3)
        Di = interpolate.UnivariateSpline(omD, Delta.imag, s=0,k=3)
        dDeltaw = array([Dr.derivatives(x)[1]+Di.derivatives(x)[1]*1j for x in om])
        return dDeltaw

    cxx = 0.0
    lngh = min(len(Dlt),len(Gfc))
    
    #print 'sum(abs(omD-om))=', sum(abs(omD-om))/len(omD)
    if len(omD)!=len(om) or sum(abs(omD-om))>1e-5*len(omD):
        for M in range(len(om)+1):
            if M<len(om) and (om[M]-omD[-1]>2e-5): break
        #print 'om[M]-omD[-1]=', om[M]-omD[-1]>2e-5
        oms = om[:M]
        Gfs = Gfc[:,:M]
    else:
        oms = om
        Gfs = Gfc
    if oms[-1]>omD[-1]: oms[-1]=omD[-1]
    if oms[0]<omD[0]: oms[0]=omD[0]
    
    print('oms[-1]=', oms[-1], 'omD[-1]=', omD[-1], file=fout)
    #print 'oms[0]=', oms[0], 'oms[-1]=', oms[-1], 'om[-1]=', om[-1]
    
    for i in range(lngh):
        dDlt = GiveDerivative(oms, omD, Dlt[i])
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
    (Deg, Eks) = ReadCix(dire+'/'+cix)
    
    
    (om, omD, Gfc, Sigc, Dlta, nf, mu, Ekin, Epot, TrSigmaG, mom) = ReadGandSigandDelta(dire+'/'+Gf, dire+'/'+Sig, dire+'/'+Delta, fout)
    s_oo = Sigc[:,-1].real
    print('s_oo=', s_oo, file=fout)
    
    lengh = min(len(Eks),len(s_oo))
    EimpS = Eks[:lengh]-mu_qmc+s_oo[:lengh]
    # Tr(log(-Gimp))
    lnGimp = LogGimp(om,Gfc,EimpS,beta,Deg,fout) 

    lnGS = LogGimpSig(om, Gfc, Sigc, Vdc, EimpS, beta, Deg, fout)

    # Tr(om*dDelta/dom)
    cxx = CmpImpurityCorrection2(om, omD, Dlta, Gfc, Deg, EimpS, beta, fout)

    # Eimp = Tr(Delta*G)+1/2*Tr(Sigma*G)+eimp*nimp-Tr(om*dDelta/dom)
    Eimp = Ekin+Epot-cxx-mu*nf                # Impurity Energy

    # Potthof2[G] == Phi[G]-1/2*Tr(Sigma*G)+T*S_imp = Eimp-Tr(log(-G_imp))+1/2*Tr(Sigma*G)
    Potthof2 = Eimp-lnGimp+TrSigmaG            # Phi[G]-1/2*Tr(Sigma*G)+T*S_imp
    
    Eimp1 = lnGS - lnGimp + Eimp  # Tr(log(1+Sigma*G))-Tr(log(-Gloc))+E_imp

    Eimp2 = TrSigmaG   # 1/2*Tr(Sigma*G)
    print('lnGimp=', lnGimp, file=fout)
    print("%-12s "*6 % ('Eimp', 'Ekin', 'Epot', 'cxx', 'mu*nf', 'lnGimp'), file=fout)
    print("%-12.6f "*6 % (Eimp, Ekin, Epot, cxx, mu*nf, lnGimp), file=fout)
    print('Eimp=', Eimp, file=fout)
    print('P[Sigma]=', Potthof2, file=fout)
    print('Tr(ln(1+S*G))-Tr(ln(-G))+E_imp=', Eimp1, file=fout)
    #return (Potthof2, lnGimp, Eimp, Eimp2)
    return (Eimp1, lnGimp, Eimp, Eimp2)


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

    (Eimp1, lnGimp, Eimp, Eimp2) = GiveImpFunctional(dire, 'PARAMS', Vdc)
    print('ev: Eimp=', Eimp, 'ln(1+s*g)-ln(g)+E_imp=', Eimp1, 'logGimp=', lnGimp)
    print('Ry: Eimp=', Eimp/Ry2eV, 'ln(1+s*g)-ln(g)+E_imp=', Eimp1/Ry2eV, 'logGimp=', lnGimp/Ry2eV)
    
