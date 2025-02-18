#!/usr/bin/env python
# @Copyright 2007 Kristjan Haule
from numpy import *
import numpy as np
from scipy import integrate
from scipy import interpolate
from scipy import special
import yw_excor                       # LDA exchange-correlation 
from gpoint import tfpoint            # generates points on a unit sphere
from CmpSlaterInt import CmpSlaterInt # Slater integrals
from CmpSlaterInt import  SlaterF2J   # translation from Slater parameters to Hunds coupling
from CmpSlaterInt import GetWaveF     # DMFT projector real space wave function
from CoulUs import CoulUsC2_diagonal  # local DMFT Coulomb U matrix, not four but two index matrix
from readTrans import ReadTrans       # T2C transformation from spherical harmonics to local DMFT basis

Ry2eV = 13.60569193

def Sigind_deg(Sigind):
    csize = max([Sigind[i,i] for i in range(len(Sigind))])
    deg=zeros(csize,dtype=int)
    for i in range(len(Sigind)):
        ii = Sigind[i,i]-1
        if ii>=0:
            deg[ii]+=1
    return deg

def N_occ_large(nocc,Sigind,deg):
    ### creates large array n_occ of size 2*(2*l+1)
    n_occ=zeros(len(Sigind))
    csize = len(deg) 
    for i in range(len(Sigind)): 
        ii = Sigind[i,i]-1
        if ii>=0:
            n_occ[i] = nocc[ii]

    for i in range(len(Sigind)):
        ii = Sigind[i,i]-1
        if ii>=0:
            n_occ[i] *= 1./deg[ii]
    return n_occ

def Anisimov(nocc,Fk,l, fh_info):
    nd = sum(nocc)
    UJ = SlaterF2J(Fk,l)
    #print('U=', UJ[0], 'Js=', UJ[1:], 'n_imp=', nd)
    print('U=', UJ[0], 'Js=', UJ[1:], 'n_imp=', nd, file=fh_info)
    n0 = int(nd+0.1)
    return (UJ[0]*(nd-0.5)-UJ[1]/2.*(nd-1), UJ[0]*(n0-0.5)-UJ[1]/2.*(n0-1) )

def AverageHartreeFock(nocc,Fk,l):
    J_HF = zeros((4,4))
    J_HF[0,:1]=[1.]                  # s-electrons
    J_HF[1,:2]=[1.,2./5.]            # p-electrons
    J_HF[2,:3]=[1.,32./65.,20./117.] # d-electrons
    J_HF[3,:4]=[1.,110./173.,167./1038.,475./4498.] # f-electrons
    
    nd = sum(nocc)
    nd1 = nd/(2*(2*l+1.))
    UJ = array(SlaterF2J(Fk,l))
    Vu = UJ[0]*(nd-nd1)
    Vj = nd1* sum( [UJ[i]*J_HF[l,i] for i in range(1,len(UJ))])
    #print 'Vu=', Vu, 'Vj=', Vj
    return Vu-Vj

def ComputeExactDC(nocc,lmbda,epsilon,icase,fh_info,projector='../projectorw.dat',trans='Trans.dat',o_Nk=8,uniform=False):
    print('ComputeExactDC: icase=', icase, 'of projector "'+projector+'"', file=fh_info)
    print('ComputeExactDC: lambda=', lmbda, file=fh_info)
    print('ComputeExactDC: epsilon=', epsilon, file=fh_info)
    print('ComputeExactDC: nocc=', nocc, file=fh_info)
    print('ComputeExactDC: k radial (2**k+1)=', o_Nk, file=fh_info)
    print('ComputeExactDC: trans=', trans, file=fh_info)
    
    lmax_integration=13 # The same as in wien2k: The number of theta points is lmax_integration+1, and of phi is 2*lmax_integration+1
    
    # Reads the radial wave function for the DMFT projector
    Nk=2**o_Nk+1
    (Rx_, Ag_, Bg_, ls_) = GetWaveF(projector)

    Rx = Rx_[icase]
    Ag = Ag_[icase]
    Bg = Bg_[icase]
    l = ls_[icase]

    # Computes Slater integrals
    print('ComputeExactDC: lambda=', lmbda, file=fh_info)
    print('ComputeExactDC: Computing Slater integrals:', end=' ', file=fh_info)
    Fk = CmpSlaterInt(lmbda,Rx,Ag,Bg,l,Nk)
    Fk = array(Fk)/epsilon
    print('ComputeExactDC: Fk=', Fk.tolist(), file=fh_info)

    if not uniform:
        # Reading Trans.dat file, which contains T2C transformation and Sigind index
        print('ComputeExactDC: Reading file', trans, file=fh_info)
        (Sigind, CF) = ReadTrans(trans, fh_info)
        
        # Matrix contained in Trans.dat should be T(i,m):
        #   - i runs over real harmonics (z^2, x^2-y^2, yz, xz, xy)
        #   - m runs over complex harmonics (-2, -1, 0, 1, 2)
        # The matrix in Trans.dat, read into CF, is the transpose of
        # what we want in T2C, hence the T2C = transpose(CF) statement
        T2C = transpose(CF)
        
        if len(CF)==2*l+1:
            nspin=1
        elif len(CF)==2*(2*l+1):
            nspin=2
        else:
            print('ComputeExactDC: ERROR: Transformation CF=T2C does not have correct dimension', file=fh_info)
            sys.exit(1)
        
        deg = Sigind_deg(Sigind)
        
        print('ComputeExactDC: Sigind is:', file=fh_info)
        print(Sigind, file=fh_info)
        
        print('ComputeExactDC: T2C follows:', file=fh_info)
        print('\n'.join('   '.join('%10f %10f' % (x.real, x.imag) for x in row) for row in T2C), file=fh_info)
        print('ComputeExactDC: shape(T2C)=', shape(T2C), file=fh_info)
        print('ComputeExactDC: T2C is Unitary=', sum(abs(matrix(T2C) * matrix(T2C).H - identity(len(T2C)))), file=fh_info)
        
        n_occ = N_occ_large(nocc,Sigind,deg)
        print('ComputeExactDC: n_occ=', n_occ, file=fh_info)
        
        UC = CoulUsC2_diagonal(l, T2C)
        
        nw=len(n_occ)
        VHartree=zeros(nw)
        EHartree=0
        for k in range(len(Fk)):
            for i in range(nw):
                U_n=0
                for j in range(nw):
                    U_n += UC[k,i,j].real*n_occ[j]
                VHartree[i] += Fk[k]*U_n
                EHartree += 0.5*Fk[k]*U_n*n_occ[i]
        print('ComputeExactDC: VHartree=', VHartree, 'EHartree=', EHartree, file=fh_info)
        
        Vhartree = zeros(len(deg))
        for i in range(nw):
            ii = Sigind[i,i]-1
            if ii>=0: Vhartree[ii] += VHartree[i]
        Vhartree[:] = Vhartree[:]/deg[:]
        print('ComputeExactDC: Vhartree=', Vhartree, 'EHartree=', EHartree, file=fh_info)
    else:
        nf=sum(nocc)
        deg=array([1])
        Vhartree=array([Fk[0]*nf])
        EHartree = 0.5*Fk[0]*nf**2
        print('ComputeExactDC: VHartree=', Vhartree, 'EHartree=', EHartree, file=fh_info)
        
    # We here construct a regular mesh of theta x phi points
    #n_the = lmax_integration+1
    #n_phi = (2*lmax_integration+1)
    (thetas,phis,tweights) = tfpoint(lmax_integration)

    ylm = zeros((2*l+1,len(thetas),len(phis)),dtype=complex)
    for m in range(2*l+1):
        for i in range(len(thetas)):
            ylm[m,i,:] = special.sph_harm(m-l, l, phis, thetas[i])
    
    DEGENERATE=False
    #if abs(min(Vhartree)-max(Vhartree))<1e-4: DEGENERATE=True
    
    CIN = 1/137.0359895**2
    fA = interpolate.UnivariateSpline(Rx,Ag,s=0)
    fB = interpolate.UnivariateSpline(Rx,Bg,s=0)
    r = linspace(1e-12,Rx[-1],Nk)
    An = fA(r)
    Bn = fB(r)
    ul2 = An**2 + CIN*Bn**2
    
    Vdc_Anisimov,Vdc_fixn = Anisimov(nocc,Fk,l,fh_info)
    Vahf = AverageHartreeFock(nocc,Fk,l)
    
    if DEGENERATE or uniform:
        ntot = sum(nocc)
        
        rho = ul2/(4*pi*r**2)*ntot
        rho = clip(rho,1e-10,1e10) # if rho becomes singular it should be removed   
        rs_1 = ((4*pi*rho)/3)**(1./3.)
        
        (epx,Vxr) = yw_excor.exchangelda(rs_1,lmbda) #
        #(epc,Vcr) = yw_excor.corrlda(1./rs_1, lmbda) # here was a bug. We should not use corrlda_2 with epsilon, like in non-DEGENERATE
        (epc,Vcr) = yw_excor.corrlda_2(1./rs_1, lmbda,epsilon) # here was a bug. We should not use corrlda_2 with epsilon, like in non-DEGENERATE
        
        Vx = integrate.romb(ul2*Vxr,dx=r[1]-r[0]) / epsilon
        Vc = integrate.romb(ul2*Vcr,dx=r[1]-r[0]) #/ epsilon   # I think corrlda_2 accounts for epsilon, hence should not divide
        ex = integrate.romb(ul2*epx,dx=r[1]-r[0])*ntot / epsilon
        ec = integrate.romb(ul2*epc,dx=r[1]-r[0])*ntot #/epsilon # I think corrlda_2 accounts for epsilon, hence should not divide

        Vhartree = sum(Vhartree[:]*deg[:])/sum(deg)
        Vdc = Vhartree+(Vx+Vc)*Ry2eV
        PhiDC = EHartree + (ex+ec)*Ry2eV
        #print('ComputeExactDC: PhiDC=', PhiDC, 'Vdc=', Vdc, 'Vdc_Anisimov=', Vdc_Anisimov, 'Vdc_fixn=', Vdc_fixn, 'Vh=', Vhartree, 'Vx=', Vx*Ry2eV, 'Vc=', Vc*Ry2eV, 'Eh=', EHartree, 'Ex=', ex*Ry2eV, 'Ec=', ec*Ry2eV, '<VHF>=', Vahf)
        print('ComputeExactDC: PhiDC=', PhiDC, 'Vdc=', Vdc, 'Vdc_Anisimov=', Vdc_Anisimov, 'Vdc_fixn=', Vdc_fixn, 'Vh=', Vhartree, 'Vx=', Vx*Ry2eV, 'Vc=', Vc*Ry2eV, 'Eh=', EHartree, 'Ex=', ex*Ry2eV, 'Ec=', ec*Ry2eV, '<VHF>=', Vahf, file=fh_info)
        return (Vdc*ones(len(deg)), PhiDC)
    else:
        
        Y_i = zeros((nspin,nw,len(thetas),len(phis)),dtype=complex)
        for i in range(nw):
            for t in range(len(thetas)):
                for f in range(len(phis)):
                    for ispin in range(nspin):  # Y_{i} = sum_m  Y_{lm}*T2C[m,i]
                        Y_i[ispin,i,t,f] = sum( ylm[:,t,f]*T2C[(2*l+1)*ispin:(2*l+1)*(ispin+1),i])  # for botsh spins
        
        Project = zeros((len(deg),len(thetas),len(phis)))
        for i in range(nw):
            ii = Sigind[i,i]-1
            if ii>=0:
                for ispin in range(nspin):
                    Project[ii,:,:] += (Y_i[ispin,i,:,:]*conj(Y_i[ispin,i,:,:])).real
        for ii in range(len(deg)):
            Project[ii,:,:] *= 1./deg[ii]

        
        rho_angle2 = zeros((len(thetas),len(phis)),dtype=float)
        for ii in range(len(deg)):
            rho_angle2 += Project[ii,:,:]*nocc[ii]
        
        VxI = zeros((len(deg),len(r)))
        VcI = zeros((len(deg),len(r)))
        ExcI = zeros(len(r))
        for ir,rx in enumerate(r):
            rho = ul2[ir]/r[ir]**2 * rho_angle2
            rho = clip(rho,1e-10,1e10) # if rho becomes singular it should be removed
            rs_1 = ((4*pi*rho)/3)**(1./3.)
            rs_1 = rs_1.flatten()
            
            (epx,Vx) = yw_excor.exchangelda(rs_1,lmbda)#*renorm_lambda)
            epx *= 1./epsilon
            Vx *= 1./epsilon
            
            #rsi = array([1/rsx if rsx>0 else 1e30 for rsx in rs_1])
            rsi = 1./rs_1
            #rsi = clip(rsi,1e-14,40)  # for rs=40 the value of Vc~1e-300 and sometimes even nan, hence need to cut.
            #(epc,Vc) = yw_excor.corrlda(rsi, lmbda*renorm_lambda)
            #epc *= 1./epsilon
            #Vc *= 1./epsilon
            (epc,Vc) = yw_excor.corrlda_2(rsi,lmbda,epsilon)
            
            #w_nan = np.isnan(Vc)
            #w_rsi = rsi[w_nan]
            #if len(w_rsi)>0:
            #    print('w_rsi=', w_rsi.tolist())
            #    print('lmbda=', lmbda, 'epsilon=', epsilon)
            #    print('w_Vc=', Vc[w_nan])
            #    print('w_epx=', epc[w_nan])
            
            Vxr = Vx.reshape(shape(rho_angle2))
            Vcr = Vc.reshape(shape(rho_angle2))
            epxc = (epx+epc).reshape(shape(rho_angle2))
            #print 'Vxcr=', Vxcr
            
            for ii in range(len(deg)):
                ToIntgx = Project[ii,:,:]*Vxr[:,:]
                ToIntgc = Project[ii,:,:]*Vcr[:,:]
                ToIntg = Project[ii,:,:]*epxc[:,:]
                tresx=0.0
                tresc=0.0
                tres = 0.0
                for t in range(len(thetas)):
                    tresx += sum(ToIntgx[t,:])*(2*pi)/len(phis)*tweights[t]
                    tresc += sum(ToIntgc[t,:])*(2*pi)/len(phis)*tweights[t]
                    tres  += sum(ToIntg[t,:])*(2*pi)/len(phis)*tweights[t]
                VxI[ii,ir]=tresx
                VcI[ii,ir]=tresc
                ExcI[ir] += tres*nocc[ii]
                
            #print 'VxcI=', VxcI[:,ir]

        Exc = integrate.romb(ul2[:]*ExcI[:],dx=r[1]-r[0])*Ry2eV
        PhiDC = EHartree+Exc
        Vedc = zeros(len(Vhartree))
        for ii in range(len(deg)):
            Vx = integrate.romb(ul2[:]*VxI[ii,:],dx=r[1]-r[0])*Ry2eV
            Vc = integrate.romb(ul2[:]*VcI[ii,:],dx=r[1]-r[0])*Ry2eV
            Vxc = Vx+Vc
            Vdc = Vhartree[ii]+Vxc
            Vedc[ii]=Vdc
            #print ii, 'Vdc=', Vhartree[ii]+Vxc, 'Vhartree=', Vhartree[ii], 'Vxc=', Vxc
            #print('ComputeExactDC: ',ii, 'PhiDC=', PhiDC, 'Vdc=', Vdc, 'Vdc_Anisimov=', Vdc_Anisimov, 'Vdc_fixn=', Vdc_fixn, 'Vh=', VHartree[ii], 'Vx=', Vx, 'Vc=', Vc)
            print('ComputeExactDC: ',ii, 'PhiDC=', PhiDC, 'Vdc=', Vdc, 'Vdc_Anisimov=', Vdc_Anisimov, 'Vdc_fixn=', Vdc_fixn, 'Vh=', VHartree[ii], 'Vx=', Vx, 'Vc=', Vc, file=fh_info)
            
        return (Vedc, PhiDC)
        
        
if __name__ == '__main__':
    import argparse
    usage = """
    The script calculates DMFT double-counting using projectorw.dat
    and screening parameter for yukawa potential lambda.

    It evaluates Phi_{xc}(rho_local,Vc_screened)
    """
    parser = argparse.ArgumentParser(description=usage)
    parser.add_argument('nocc', action='store')
    parser.add_argument('-l', '--lmbda', action='store', type=float, default=1.0, help='lambda for screening')
    parser.add_argument('-e', '--epsilon', action='store', type=float, default=1.0, help='epsilon for screening')
    parser.add_argument('-i', '--icase', action='store', type=int, default=0, help='which unique projector shoule be used')
    parser.add_argument('-k', '--Nk', action='store', type=int, default=8, help='number of radial points is 2**k+1')
    parser.add_argument('-p', '--projector', action='store', default='projectorw.dat', help='filename of projector')
    parser.add_argument('-t', '--trans', action='store', default='imp.0/Trans.dat', help='filename which contains Sigind and T2C')
    parser.add_argument('-o', '--info', action='store', type=argparse.FileType('w'), default='DC.info', help='log filename')
    parser.add_argument('-u', '--uniform', action='store', type=bool, default=False, help='should one just take the average over all orbitals')
    
    # Next, parse the arguments
    options = parser.parse_args()
    
    nocc=eval(options.nocc)
    
    #icase=options.icase
    #lmbda=options.lmbda
    #fh_info = options.info 
    #projector = options.projector
    #renorm_lambda=options.renorml
    #trans = options.trans
    #uniform = options.uniform
    #o_Nk = options.Nk
    
    (Vdc,PhiDC)=ComputeExactDC(nocc,options.lmbda,options.epsilon,options.icase,options.info,options.projector,options.trans,options.Nk,options.uniform)#,options.renorml)
    print('Finished with', Vdc, PhiDC)
