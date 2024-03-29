#!/usr/bin/env python
# @Copyright 2007 Kristjan Haule
import sys
from scipy import *
from scipy import special
from scipy import integrate
import readCharge
import ylm
from gpoint import tfpoint
from CmpSlaterInt import GetWaveF
import struct1
import yw_excor
from readIndmfl import ReadIndmfl
import copy

Ry2eV = 13.60569193

def GetProjector(file_projector, file_indmfl, icases, y_lm, fh_info):
    def Sigind_deg(Sigind):
        csize = max([Sigind[i,i] for i in range(len(Sigind))])
        deg=zeros(csize,dtype=int)
        for i in range(len(Sigind)):
            ii = Sigind[i,i]-1
            if ii>=0:
                deg[ii]+=1
        return deg

    (Rx_, Ag_, Bg_, ls_) = GetWaveF(file_projector)

    print('Reading file', file_indmfl, file=fh_info)

    (Siginds,CFs,cps) = ReadIndmfl(file_indmfl,fh_info)
    
    Projects={}
    iatoms={}
    degs={}
    Ag={}
    Bg={}
    Rx={}
    for icix in list(Siginds.keys()):
        Sigind=Siginds[icix]
        CF = CFs[icix]
        icase = icases[icix]

        (iatom,L,qsplit) = cps[icix][0]
        
        Rx[icix] = Rx_[icase]
        Ag[icix] = array(Ag_[icase])
        Bg[icix] = array(Bg_[icase])
        l = ls_[icase]
        
        if (l!=L):
            print('Not equal l. You have wrong icases!')
            sys.exit(1)
        
        T2C = transpose(CF)
        if len(CF)==2*l+1:
            nspin=1
        elif len(CF)==2*(2*l+1):
            nspin=2
        else:
            print('ERROR: Transformation CF=T2C does not have correct dimension', file=fh_info)
            sys.exit(1)
        deg = Sigind_deg(Sigind)
        
        ylm = y_lm[l**2:(l+1)**2,:,:]
        
        nw=len(Sigind)
        Y_i = zeros((nspin,nw,shape(ylm)[1],shape(ylm)[2]),dtype=complex)
        for i in range(nw):
            for t in range(shape(ylm)[2]):
                for f in range(shape(ylm)[1]):
                    for ispin in range(nspin):  # Y_{i} = sum_m  Y_{lm}*T2C[m,i]
                        Y_i[ispin,i,f,t] = sum( ylm[:,f,t]*T2C[(2*l+1)*ispin:(2*l+1)*(ispin+1),i])  # for botsh spins
        
        Project = zeros((len(deg),shape(ylm)[1],shape(ylm)[2]),order='F')
        for i in range(nw):
            ii = Sigind[i,i]-1
            if ii>=0:
                for ispin in range(nspin):
                    Project[ii,:,:] += (Y_i[ispin,i,:,:]*conj(Y_i[ispin,i,:,:])).real
        for ii in range(len(deg)):
            Project[ii,:,:] *= 1./deg[ii]
        
        Projects[icix]=Project
        iatoms[icix]=iatom
        degs[icix]=deg

        # Y_i = sum_m  Y_{lm}*T2C[m,i]
        # Project[i] = Y_i * Y_i^*
    return (Projects, Ag, Bg, Rx, Siginds, iatoms, degs)

if __name__ == '__main__':

    dir='../../'
    case='alpha-Ce'
    
    options_projector = dir+'/projectorw.dat' 
    icases={1:0} # {icix-in-indmfl:icase-in-projectorw}
    
    file_clmsum=dir+'/'+case+'.clmsum'
    options_indmfl=dir+'/'+case+'.indmfl'
    
    
    struct = struct1.Struct(dir+'/'+case)
    
    
    R0 = struct.r0
    Rmt = struct.rmt
    nat = struct.nat
    jri = struct.npt
    iatnr = struct.iatom
    
    iatom2jatom=[]
    for i in range(len(struct.mult)):
        iatom2jatom += [i]*struct.mult[i]
    
    print('iatom2jatom=', iatom2jatom)
    
    print('iatnr=', iatnr)
    print('jri=', jri)
    print('R0=', R0)
    print('Rmt=', Rmt)
    
    
    lmax_integration=13
    
    lmmaxx,lmax2,ncom,nrad,lmMax,lmMax1 = readCharge.checkchargesize(file_clmsum,jri)
    
    print('lmmaxx=', lmmaxx, 'lmax2=', lmax2, 'ncom=', ncom, 'nrad=', nrad)
    print('lmax=', lmMax)
    print('lmax1=', lmMax1)
    
    clm,lm = readCharge.readcharge(file_clmsum,lmmaxx,lmax2,ncom,nrad,lmMax,lmMax1,jri)
    
    print('shape(lm)=', shape(lm))
    print('lm=', transpose(lm))
    print('clm=', shape(clm))
    
    c_kub = readCharge.cmp_c_kub()
    
    lmax = max(lm[0,:,:].ravel())
    
    lmax_integration=max(lmax_integration,lmax)
    
    #n_the = lmax_integration+1
    #n_phi = (2*lmax_integration+1)
    (thetas,phis,tweights) = tfpoint(lmax_integration)
    
    
    # Precomputes spherical harmonics for all l,m -s used later
    # This is because the same Y_{lm} can be used for all atoms and cix.
    y_lm = zeros( ((lmax+1)**2,len(phis),len(thetas)), order='F', dtype=complex)
    for it in range(len(thetas)):
        for ip in range(len(phis)):
            if False:
                for l in range(0,lmax+1):
                    y_lm[l**2:(l+1)**2,ip,it] = special.sph_harm(arange(-l,l+1), l, phis[ip], thetas[it])
            else:
                v=array([sin(thetas[it])*cos(phis[ip]),sin(thetas[it])*sin(phis[ip]),cos(thetas[it])])
                y_lm[:,ip,it] = ylm.ylm(v,lmax)
    
    
    (Projects, Ag, Bg, Rx, Siginds,iatoms,degs) = GetProjector(options_projector, options_indmfl, icases, y_lm, sys.stdout)
    print(shape(Projects[1]))
    print('iatoms=', iatoms)
    CIN = 1/137.0359895**2
    
    for icix in list(Projects.keys()):
        # First compute a combination of spheric harmonics, which are used in the expansion of charge density
        # The correct linear combination of Y_{lm}'s is called YY[lm,phi,theta] where phi=[0,2pi] and theta=[0,pi].
        jatom = iatom2jatom[iatoms[icix]-1]
        cubic = iatnr[jatom]>0
        lmMax_ = lmMax[jatom]
        YY = zeros( (lmMax_,len(phis),len(thetas)), order='F',dtype=float)
        for it in range(len(thetas)):
            for ip in range(len(phis)):
                (YY[:,ip,it],llmm) = readCharge.cmp_rho_spherical(y_lm[:,ip,it],lm[:,:lmMax_,jatom],c_kub,cubic)
    
        YY = reshape(YY, (lmMax_,len(phis)*len(thetas)), order='F')  # For convernience, we make 2D array, such that YY[lm,generalized-angle]
        
        #YY = YY.reshape((lmMax_,len(phis),len(thetas)),order='F')
        #print 'Diff=', sum(abs(YY_backup-YY))
        
        nr=jri[jatom]  # This is the number of radial points
        rmesh = array(Rx[icix]) # Radial mesh should also be equal to
                                # rmesh = R0[jatom]*exp(arange(jri[jatom])/(jri[jatom]-1.)*log(Rmt[jatom]/R0[jatom]))

        # Using the radial charge distribution clm[r,lm,jatom] and proper combination of spherical harmonics YY[lm,generalized-angle], we obtain
        # charge density rho[r,generalized-angle]
        rho = readCharge.cmp_rho_radial(clm[:nr,:lmMax_,jatom],YY,rmesh,lm[:,:lmMax_,jatom],c_kub,cubic,llmm)
        
        #for ir in range(len(rmesh)):
        #    print "%13.8f " % rmesh[ir], ("%13.8f "*shape(rho)[1]) % tuple(rho[ir,:]*rmesh[ir]**2)
    
        #print 'YY='
        #for i in range(shape(YY)[1]):
        #    #print i, YY[0,i], YY[1,i], YY[2,i]
        #    print i, rho[-1,i], clm[nr-1,0,0], clm[nr-1,1,0]
            
        #for i in range(shape(rho)[1]):
        #    print i, rho[-1,i]/2.
    
        
        rho = rho.reshape(nr*len(phis)*len(thetas),order='F')
        # Here we construct one dimensional array of density. We start from rho[r,generalized-angle] to rho[generalized-point]
        # where generalized point is array of [r,phi,theta]
        # To get back charge density in any radial point, one should simply reshape the array:
        #                                rho = reshape(rho, (nr,len(phis),len(thetas)), order='F')
        #

        rs_1 = ((4*pi*rho)/3)**(1./3.)   # rs^{-1} computed from rho for each point is space (inside muffin-thin)
        rsi = array([1/rsx if rsx>1e-30 else 1e30 for rsx in rs_1])  # This should be rs in each point in space
    
        small=0.01
        rsp = rsi+small*rsi                           # rsp = r+epsilon  --- will be needed for taking derivative numerically
        (epx,Vxp) = yw_excor.exchangelda(1/rsp,0.0)   # Here we compute V_x(rs+epsilon) which will be used later to compute derivative 
        (epc,Vcp) = yw_excor.corr_lda(rsp)            # Here we compute V_c(rs+epsilon)
        rsm = rsi-small*rsi                           # rsp = r-epsilon  --- will be needed for taking derivative numerically
        (epx,Vxm) = yw_excor.exchangelda(1/rsm,0.0)   # V_x(rs-epsilon)
        (epc,Vcm) = yw_excor.corr_lda(rsm)            # V_c(rs-epsilin)
        dVxc_drs = (Vxp+Vcp-Vxm-Vcm)/(2*small*rsi)    # dV_{xc}/drs
        dVxc_drho = -(4.*pi/9.)*(rsi**4)*dVxc_drs     # dV_{xc}/drho computed from dV_{xc}/drs
        
        (epx,Vx) = yw_excor.exchangelda(1/rsi,0.0)    # V_{x}(rs)
        (epc,Vc) = yw_excor.corr_lda(rsi)             # V_{c}(rs)
        Vxc = Vx+Vc                                   # V_{xc}(rs)
        
        Vxcr = Vxc.reshape((nr,len(phis),len(thetas)),order='F')             # Finally, V_xc for each point in space
        dVxc_drhor = dVxc_drho.reshape((nr,len(phis),len(thetas)),order='F') # and dVxc/drho for each point in space
        
        Project=Projects[icix]                       # This is the angle part of the projector to correlated orbital, which is constructed from spheric harmonics Y_{lm}
                                                     # and T2C transformation as Y_i = sum_m  Y_{lm}*T2C[m,i]
                                                     # and Project[i] = Y_i * Y_i^*
        deg=degs[icix]                               
        ul2 = Ag[icix]**2 + CIN*Bg[icix]**2          # The radial part of the projector including relativistic correction
        Vxc=[]                                       # Here we compute Project*V_{xc} = <phi|V_{xc}|phi> with |phi> the correlated orbital
        for ii in range(len(deg)):
            VxcI = zeros((len(rmesh)))
            for ir,rx in enumerate(rmesh):
                ToIntgxc = Project[ii,:,:]*Vxcr[ir,:,:]         #  Project*Vxc = Y_i^*(theta,hi)*V_xc(theta,phi)Y_i(theta,phi)
                tresxc=0.0
                for t in range(len(thetas)):                    # For integration over theta and phi
                    tresxc += sum(ToIntgxc[:,t])*(2*pi)/len(phis)*tweights[t]
                VxcI[ir]=tresxc                                 # only integration over radial points is left, (Project*Vxc)(r)
            Vxc.append( integrate.simps(ul2[:]*VxcI[:],x=rmesh)*Ry2eV )
            
        print('Vxc=', Vxc)                                       # This is projection of V_xc to correlated orbital
        
        dVxc_drho=[]                                            # Now we do the same for dV_{xc}/drho, which should be projected to correlated orbital as follows:
        for ii in range(len(deg)):                              # <phi,phi|dV_{xc}/drho|phi,phi> becasue it is a two-particle quantity, which should be added to 
            for jj in range(len(deg)):                          # Coulomb repulsion 1/|r-r'|
                dVxcI = zeros((len(rmesh)))
                for ir,rx in enumerate(rmesh):
                    ToIntg   = Project[ii,:,:]*Project[jj,:,:]*dVxc_drhor[ir,:,:]  # For two particle quantity we have Project*Project*dV_{xc}/rho
                    tres=0.0
                    for t in range(len(thetas)):                # integration over angle
                        tres   += sum(ToIntg  [:,t])*(2*pi)/len(phis)*tweights[t]
                    dVxcI[ir]=tres
                dVxc_drho.append( integrate.simps((ul2[:]**2/rmesh[:]**2)*dVxcI[:],x=rmesh)*Ry2eV )  # integration over r
    
        
        print('dVxc/drho=', dVxc_drho)                          # Final result of <phi,phi|dVxc/drho|phi,phi>
        
        rho = reshape(rho, (nr,len(phis),len(thetas)), order='F')  # Not really needed here
    
