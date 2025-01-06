#!/usr/bin/env python
Fortran = True
import sys, os
from numpy import *
from numpy import linalg
import numpy as np
from scipy import interpolate
from scipy import optimize
from itertools import chain, product
from timeit import default_timer as timer
#from tanmesh import *
from timeit import default_timer as timer
import builtins
if Fortran:
    import tbg_scc as scc

class EmbeddOrbitals:
    def __init__(self, axis, beta, par, plot_bands=False):
        band = ''
        if plot_bands:
            band = '_band'
        self.Eks     = load('ene'+band+'.npy')
        npzfile = load('Uproject'+band+'.npz')
        self.UAA, self.UAB = npzfile['UAA'], npzfile['UAB']
        self.Niik = len(self.Eks)
        self.Nfunc, Nik, self.Nband = shape(self.UAA)
        if (self.Niik != Nik):
            print('ERROR : number of bands in UAA and Eks is inconsistent')
            sys.exit(1)
            
        self.wk = ones(Nik) # k-point weights. For now, we did not yet use symmetry, hence all point have unity weight
        self.beta = beta
        if axis=='real':
            if plot_bands: # need equdistant mesh for imshow
                w_min, w_max, Nw = par['plt_w']
                self.w = linspace(w_min, w_max, Nw)
            else:
                Nmax, x0, wmax = par['Nmax'], par['x0'], par['wmax']
                w_half = GiveTanMesh(x0,wmax,int(Nmax/2))
                self.w = hstack( (-w_half[::-1], w_half[:]) )
            self.nom = 100  # just so that it does not fail
        else:
            wmax, nmax_nom, ntail = par['wmax'], par['nmax_nom'], par['ntail']
            self.nom = int(round((wmax*beta/pi-1)/2.))
            nmax = nmax_nom*self.nom
            #self.w = (2*arange(Nmax)+1)*pi/beta
            self.w , self.ind_om = create_log_mesh(self.nom, nmax, ntail, beta)
    def ComputeGloc0(self, axis, mu):
        if axis=='real':
            broad_values = np.minimum(wbroad*np.abs(self.w), 4.0/unit)
            om = self.w + mu + broad_values*1j
            #wn = self.w + mu + vbroad(self.w)*1j
            #_w_ = tensordot(ones(self.Nband), wn, axes=0)
        else:
            om = self.w*1j + mu
            #_w_ = tensordot(ones(self.Nband), wn, axes=0)
        ones_w = ones(len(self.w))
        UU2 = zeros((self.Nfunc,self.Nfunc,self.Nband), dtype=complex)
        DOS = zeros((3,len(self.w)),dtype=float)
        Gc  = zeros((self.Nfunc,self.Nfunc,len(self.w),2), dtype=complex)
        Eimp = zeros((2,self.Nfunc,self.Nfunc), dtype=complex)
        WHATAA = zeros((self.Nfunc,self.Nfunc), dtype=complex)
        WHATAB = zeros((self.Nfunc,self.Nfunc), dtype=complex)
        for iik in range(self.Niik):
          #_ek_ = tensordot(self.Eks[iik,:], ones_w, axes=0)
          #denom= 1./(_w_ - _ek_)
          #
          ek = self.Eks[iik,:]
          denom = 1.0/(om[np.newaxis,:] - ek[:,np.newaxis]) # denom[iband,om] = 1/(om-ek)
          DOS[0,:] += -1/pi * sum(denom,axis=0).imag * self.wk[iik]
          for ibnd in range(self.Nband):
              UU2[:,:,ibnd] = tensordot(self.UAA[:,iik,ibnd], self.UAA[:,iik,ibnd].conj(), axes=0)
          Gc[:,:,:,0] += dot(UU2, denom) * self.wk[iik]
          Eimp[0,:,:] += dot(UU2, self.Eks[iik,:]-mu) * self.wk[iik]
          WHATAA += sum(UU2[:,:,:],axis=2) * self.wk[iik]
          for ibnd in range(self.Nband):
              UU2[:,:,ibnd] = tensordot(self.UAB[:,iik,ibnd], self.UAB[:,iik,ibnd].conj(), axes=0)
          Gc[:,:,:,1] += dot(UU2, denom) * self.wk[iik]
          Eimp[1,:,:] += dot(UU2, self.Eks[iik,:]-mu) * self.wk[iik]
          WHATAB += sum(UU2[:,:,:],axis=2) * self.wk[iik]
        Nwk = sum(self.wk)
        WHATAA *= 1./Nwk
        WHATAB *= 1./Nwk
        print('WHATAA=')
        PrintM(WHATAA)
        print('WHATAB=')
        PrintM(WHATAB)
        Gc *= 1./Nwk
        DOS  *= 1./Nwk
        Eimp *= 1./Nwk
        pDOS = zeros((2,self.Nfunc,len(self.w)),dtype=float)
        for ifc in range(self.Nfunc):
          pDOS[0,ifc,:] += -1./pi * Gc[ifc,ifc,:,0].imag
          pDOS[1,ifc,:] += -1./pi * Gc[ifc,ifc,:,1].imag
        DOS[1,:]  = sum(pDOS[0,:,:],axis=0) # AA region
        DOS[2,:]  = sum(pDOS[1,:,:],axis=0) # AA region
        return (Gc, DOS, pDOS, Eimp)
        
    def ComputeGloc_internal(self, axis, mu, Sigmas):
        if axis=='real':
            broad_values = np.minimum(wbroad*np.abs(self.w), 4.0/unit)
            wn = self.w + mu + broad_values*1j
        else:
            wn = self.w*1j + mu
        Nd0 = 0.0
        DOS = zeros((3,len(self.w)),dtype=float)
        Gc  = zeros((self.Nfunc,self.Nfunc,len(self.w),2), dtype=complex)
        Eimp = zeros((2,self.Nfunc,self.Nfunc), dtype=complex)
        dos = zeros(len(self.w),dtype=float)
        for iik in range(self.Niik):
            eks = self.Eks[iik,:]
            _UAA_ = self.UAA[:,iik,:]  # _UAA_[ifb,ibnd] = self.UAA[ifb,iik,ibnd]
            _UAB_ = self.UAB[:,iik,:]  # _UAB_[ifb,ibnd] = self.UAB[ifb,iik,ibnd]
            _UAA_H = _UAA_.T.conj()    # _UAA_H[ibnd,ifb]
            _UAB_H = _UAB_.T.conj()    # _UAB_H[ibnd,ifb]
            Eimp[0,:,:] += _UAA_*(eks-mu) @ _UAA_H
            Eimp[1,:,:] += _UAB_*(eks-mu) @ _UAB_H
            for iw in range(len(self.w)):
                sigA_iw = Sigmas[:, iw, 0]
                sigB_iw = Sigmas[:, iw, 1]
                tmpA = sigA_iw[:,None] * _UAA_ # tmpA[ifc,ibnd] = Sigmas[ifc,iw,iAA]*_UAA_[ifc,ibnd]
                tmpB = sigB_iw[:,None] * _UAB_ # tmpB[ifc,ibnd] = Sigmas[ifc,iw,iAB]*_UAB_[ifc,ibnd]
                USigU = _UAA_H @ tmpA + _UAB_H @ tmpB
                gk = linalg.inv( diag(wn[iw]-eks[:]) - USigU )
                dos[iw] = np.trace(gk).real
                Gc[:,:,iw,0] += _UAA_ @ gk @ _UAA_H
                Gc[:,:,iw,1] += _UAB_ @ gk @ _UAB_H
                
            eps = eks + diag(USigU).real # the large omega static limit, eps[iband]
            tdos = np.sum( 1.0/(wn[:,None]-eps[None,:]),axis=1).real              
            dos -= tdos
            Nd0 += np.sum(1.0/(1.0 + np.exp(beta*(eps-mu))))
            DOS[0,:] += dos
        Nwk = self.Niik
        Gc *= 1./Nwk
        DOS  *= 1./Nwk
        Eimp *= 1./Nwk
        Nd0 *= 1/Nwk
        pDOS = zeros((2,self.Nfunc,len(self.w)),dtype=float)
        for ifc in range(self.Nfunc):
            pDOS[0,ifc,:] += -1./pi * Gc[ifc,ifc,:,0].imag
            pDOS[1,ifc,:] += -1./pi * Gc[ifc,ifc,:,1].imag
        DOS[1,:]  = sum(pDOS[0,:,:],axis=0) # AA region
        DOS[2,:]  = sum(pDOS[1,:,:],axis=0) # AA region
        
        if axis=='imag':
            n_max = round((self.w[-1]*self.beta/pi-1.)/2.)
            w_all = (2*arange(n_max+1)+1)*pi/self.beta
            fr = interpolate.UnivariateSpline(self.w, DOS[0,:], s=0)
            nbnd = 2./self.beta * sum(fr(w_all))
            Ntot = 2*(Nd0+nbnd)
            print('Ntot=', Ntot, 'Nd0=', 2*Nd0, '2*nbnd=', 2*nbnd)
        else:
            Ntot=0
        return (Gc, DOS, pDOS, Eimp, Ntot)
    
    def ComputeGloc(self, axis, mu, Sigmas):
        #if sum(abs(Sigmas))<1e-5:
        #    return self.ComputeGloc0(axis, mu)
        if axis=='real':
            wn = self.w + mu + vbroad(self.w)*1j
        else:
            wn = self.w*1j + mu
        Nd0 = 0.0
        DOS = zeros((3,len(self.w)),dtype=float)
        Gc  = zeros((self.Nfunc,self.Nfunc,len(self.w),2), dtype=complex, order='F')
        Eimp = zeros((2,self.Nfunc,self.Nfunc), dtype=complex)
        for iik in range(self.Niik):
          eks = self.Eks[iik,:]
          Uk = array( [self.UAA[:,iik,:], self.UAB[:,iik,:]], order='F').transpose(1,2,0) # Uk[Nfunc,Nband,2]
          Eimp[0,:,:] += dot( Uk[:,:,0]*(eks-mu), Uk[:,:,0].conj().T ) * self.wk[iik]
          Eimp[1,:,:] += dot( Uk[:,:,1]*(eks-mu), Uk[:,:,1].conj().T ) * self.wk[iik]
          gc,dos,_Nd0_ = scc.cmp_local_g(Uk, Sigmas, wn, eks, self.wk[iik], self.beta, axis=='imag')
          Nd0 += _Nd0_
          Gc[:,:,:,:] += gc
          DOS[0,:] += dos
        Nwk = sum(self.wk)
        Gc *= 1./Nwk
        DOS  *= 1./Nwk
        Eimp *= 1./Nwk
        Nd0 *= 1/Nwk
        pDOS = zeros((2,self.Nfunc,len(self.w)),dtype=float)
        for ifc in range(self.Nfunc):
          pDOS[0,ifc,:] += -1./pi * Gc[ifc,ifc,:,0].imag
          pDOS[1,ifc,:] += -1./pi * Gc[ifc,ifc,:,1].imag
        DOS[1,:]  = sum(pDOS[0,:,:],axis=0) # AA region
        DOS[2,:]  = sum(pDOS[1,:,:],axis=0) # AA region

        if axis=='imag':
            n_max = round((self.w[-1]*self.beta/pi-1.)/2.)
            w_all = (2*arange(n_max+1)+1)*pi/self.beta
            fr = interpolate.UnivariateSpline(self.w, DOS[0,:], s=0)
            nbnd = 2./self.beta * sum(fr(w_all))
            Ntot = 2*(Nd0+nbnd)
            print('Ntot=', Ntot, 'Nd0=', 2*Nd0, '2*nbnd=', 2*nbnd)
        else:
            Ntot=0
        return (Gc, DOS, pDOS, Eimp, Ntot)

    def JustComputeDensity(self, Sigmas, mu):
        # First compute frequency dependent eigenvalues
        zek = zeros( (self.Nband,len(self.w),self.Niik), dtype=complex, order='F' )
        for iik in range(self.Niik):
          eks = self.Eks[iik,:]
          Uk = array( [self.UAA[:,iik,:], self.UAB[:,iik,:]], order='F').transpose(1,2,0) # Uk[Nfunc,Nband,2]
          zek[:,:,iik] = scc.cmp_eigenvals(Uk, Sigmas, eks) # calling fast fortran routine

        n_mu = N_mu(0.0, self.w, zek, self.Eks, self.wk, self.beta)
        return n_mu(mu)

    def JustComputeEigenvalues(self, Sigmas, mu):
        # First compute frequency dependent eigenvalues
        zek = zeros( (self.Nband,len(self.w),self.Niik), dtype=complex, order='F' )
        for iik in range(self.Niik):
          eks = self.Eks[iik,:]
          Uk = array( [self.UAA[:,iik,:], self.UAB[:,iik,:]], order='F').transpose(1,2,0) # Uk[Nfunc,Nband,2]
          zek[:,:,iik] = scc.cmp_eigenvals(Uk, Sigmas, eks) # calling fast fortran routine
        return zek

    def FermiLevel_internal(self, n_desired, Sigmas, mu0=0, dmu=0.2):
        zek = zeros( (self.Niik, len(self.w), self.Nband), dtype=complex)
        for iik in range(self.Niik):
            _UAA_ = self.UAA[:,iik,:]
            _UAB_ = self.UAB[:,iik,:]
            _UAA_H = _UAA_.T.conj()
            _UAB_H = _UAB_.T.conj()
            ham0 = np.diag(self.Eks[iik,:])
            for iw in range(len(self.w)):
                sigA_iw = Sigmas[:, iw, 0]
                sigB_iw = Sigmas[:, iw, 1]
                tmpA = sigA_iw[:,None] * _UAA_ # tmpA[ifc,ibnd] = Sigmas[ifc,iw,iAA]*_UAA_[ifc,ibnd]
                tmpB = sigB_iw[:,None] * _UAB_ # tmpB[ifc,ibnd] = Sigmas[ifc,iw,iAB]*_UAB_[ifc,ibnd]
                USigU = _UAA_H @ tmpA + _UAB_H @ tmpB
                zek[iik,iw,:] = linalg.eigvals( ham0+USigU )
        
        mu = mu0
        n_mu = N_mu_internal(n_desired, self.w, zek, beta)
        n1 = n_mu(mu)
        if n1 > 0 : dmu *= -1
        for i in range(100):
            mu_old = mu
            mu += dmu
            n2 = n_mu(mu)
            if n1*n2 < 0: # zero bracketed
                mu = optimize.brentq(n_mu, builtins.min(mu_old,mu), builtins.max(mu_old,mu), xtol=1e-12, rtol=1e-12)
                return mu
        
    def FermiLevel(self, n_desired, Sigmas, mu0=0, dmu=0.2):
        # First compute frequency dependent eigenvalues
        zek = zeros( (self.Nband,len(self.w),self.Niik), dtype=complex, order='F' )
        #_zek_ = zeros( (self.Niik, len(self.w), self.Nband), dtype=complex)
        #tms = zeros(2)
        for iik in range(self.Niik):
          #_t1_ = timer()
          eks = self.Eks[iik,:]
          Uk = array( [self.UAA[:,iik,:], self.UAB[:,iik,:]], order='F').transpose(1,2,0) # Uk[Nfunc,Nband,2]
          zek[:,:,iik] = scc.cmp_eigenvals(Uk, Sigmas, eks) # calling fast fortran routine
          #_t2_ = timer()
          
          #_UAA_ = self.UAA[:,iik,:]
          #_UAB_ = self.UAB[:,iik,:]
          #_UAA_H = _UAA_.T.conj()
          #_UAB_H = _UAB_.T.conj()
          #ham0 = np.diag(eks)
          #for iw in range(len(self.w)):
          #    sigA_iw = Sigmas[:, iw, 0]
          #    sigB_iw = Sigmas[:, iw, 1]
          #    tmpA = sigA_iw[:,None] * _UAA_ # tmpA[ifc,ibnd] = Sigmas[ifc,iw,iAA]*_UAA_[ifc,ibnd]
          #    tmpB = sigB_iw[:,None] * _UAB_ # tmpB[ifc,ibnd] = Sigmas[ifc,iw,iAB]*_UAB_[ifc,ibnd]
          #    USigU = _UAA_H @ tmpA + _UAB_H @ tmpB
          #    _zek_[iik,iw,:] = linalg.eigvals( ham0+USigU )
          #_t3_ = timer()
          #tms[0] += _t2_-_t1_
          #tms[1] += _t3_-_t2_
        #print('times=', tms)
            
        mu = mu0
        n_mu = N_mu(n_desired, self.w, zek, self.Eks, self.wk, beta)

        #n_max = round((self.w[-1]*beta/pi-1.)/2.)
        #w_all = (2*arange(n_max+1)+1)*pi/beta
        #print('n_max=', n_max, 'n=', len(self.w))
        #zek0 = _zek_[:,-1,:].real  # shape (Niik,Nband)
        #Nd0 = np.sum(1.0/(1.0 + np.exp(beta*(zek0-mu))))/self.Niik
        
        #iom = self.w*1j + mu
        #zsum = np.sum( 1.0/(iom[None,:,None]-_zek_), axis=(0, 2) )
        #wsum = np.sum( 1.0/(iom[None,:,None]-zek0[:,None,:]),axis=(0, 2))
        #Nd = (zsum - wsum).real/self.Niik
            
        #fN = interpolate.UnivariateSpline(self.w, Nd, s=0)
        #dn = 2./beta * sum(fN(w_all)) # sum over matsubara points (positive+negative matsubara give 2)
        #Ntot = 2*(Nd0 + dn) # 2 comes from spin degeneracy
        #print('mu=', mu, 'Nd=', Ntot, '2*N0=', 2*Nd0, '2*dn=', 2*dn)
        #_n1_ =  Ntot-n_desired
        
        n1 = n_mu(mu)
        if n1 > 0 : dmu *= -1
        for i in range(100):
            mu_old = mu
            mu += dmu
            n2 = n_mu(mu)
            if n1*n2 < 0: # zero bracketed
                mu = optimize.brentq(n_mu, builtins.min(mu_old,mu), builtins.max(mu_old,mu), xtol=1e-12, rtol=1e-12)
                return mu
        
    def PrintpDOS(self, filename, DOS, pDOS):
        fout = open(filename,'w')
        print('# total  AA+AB  AA   AB  ', end=' ', file=fout)
        for i in range(self.Nfunc) :
              print('funAA['+str(i)+']  ', end=' ', file=fout)
        for i in range(self.Nfunc) :
              print('funAB['+str(i)+']  ', end=' ', file=fout)
        print(file=fout)
        for iw in range(len(self.w)):
              print(self.w[iw], end=' ', file=fout)
              print(DOS[0,iw], DOS[1,iw]+DOS[2,iw], DOS[1,iw], DOS[2,iw], end=' ', file=fout)
              for i in range(len(pDOS)):
                for ifc in range(self.Nfunc):
                  print(pDOS[i,ifc,iw], end=' ', file=fout)
              print(file=fout)
        fout.close()

    def CmpDelta(self, axis, Gc, Sigma, Eimp):
        DltE = zeros(shape(Gc), dtype=complex, order='F')
        Id = eye(self.Nfunc)
        
        for iw in range(len(self.w)):
            if axis=='real':
                _w_ = self.w[iw]*Id
            else:
                _w_ = self.w[iw]*1j*Id
            for iAA in range(2):
                G_inv = linalg.inv(Gc[:,:,iw,iAA])
                Sig = diagflat(Sigma[:,iw,iAA])
                DltE[:,:,iw,iAA] = _w_ - Sig[:,:] - G_inv[:,:]
        
        for iAA in range(2):
            for i in range(self.Nfunc):
                DltE[i,i,:,iAA] -= Eimp[iAA,i,i].real
        return DltE
        
    def PrintGandDelta(self, axis, _dir_, Sigind, fname_G, fname_G_all, Gc, fname_D, fname_D_all, DeltaE):
        fg = open(_dir_+'/'+fname_G_all,'w')
        gd = open(_dir_+'/'+fname_D_all, 'w')
        da = open(_dir_+'/'+fname_G, 'w')
        for iw in range(len(self.w)):
            print(self.w[iw], end=' ', file=fg)
            print(self.w[iw], end=' ', file=da)
            for i in range(self.Nfunc):
                print(Gc[i,i,iw].real, Gc[i,i,iw].imag, '  ', end=' ', file=da)
                for j in range(self.Nfunc):
                    print(Gc[i,j,iw].real, Gc[i,j,iw].imag, '  ', end=' ', file=fg)
            print(file=fg)
            print(file=da)
            
            print(self.w[iw], end=' ', file=gd)
            for i in range(self.Nfunc):
                for j in range(self.Nfunc):
                    print(DeltaE[i,j,iw].real, DeltaE[i,j,iw].imag, '  ', end=' ', file=gd)
            print(file=gd)
        fg.close()
        gd.close()
        da.close()

        if axis=='real':
            sigind = diagonal(Sigind)
            sigind_max = max(sigind)
            aD=[]
            for i in range(sigind_max):
                aDelta = zeros(len(self.w),dtype=complex)
                inum=0
                for j in range(self.Nfunc):
                    if sigind[j]==i+1:
                        aDelta += DeltaE[j,j,:]
                        inum += 1
                aDelta /= inum
                aD.append(aDelta)
            
            fd = open(_dir_+'/'+fname_D, 'w')
            for iw,w in enumerate(self.w):
                print(w, end=' ', file=fd)
                for i in range(sigind_max):
                    print(aD[i][iw].real, aD[i][iw].imag, '  ', end=' ', file=fd)
                print(file=fd)
            fd.close()
        else:
            # forming mesh of all Matsubara points
            n_max = round((self.w[-1]*beta/pi-1.)/2.)
            w_all = (2*arange(n_max+1)+1)*pi/beta
            # interpolating on all Matsubara points
            sigind = diagonal(Sigind)
            sigind_max = max(sigind)
            fdr=[]
            fdi=[]
            aD=[]
            for i in range(sigind_max):
                aDelta = zeros(len(self.w),dtype=complex)
                inum=0
                for j in range(self.Nfunc):
                    if sigind[j]==i+1:
                        aDelta += DeltaE[j,j,:]
                        inum += 1
                aDelta /= inum
                aD.append(aDelta)
                fdr.append( interpolate.UnivariateSpline(self.w[self.nom-1:], aDelta[self.nom-1:].real, s=0) )
                fdi.append( interpolate.UnivariateSpline(self.w[self.nom-1:], aDelta[self.nom-1:].imag, s=0) )
            
            fd = open(_dir_+'/'+fname_D, 'w')
            for iw in range(self.nom):
                w = w_all[iw]
                print(w, end=' ', file=fd)
                for i in range(sigind_max):
                    print(aD[i][iw].real, aD[i][iw].imag, '  ', end=' ', file=fd)
                print(file=fd)
            for iw in range(self.nom,len(w_all)):
                w = w_all[iw]
                print(w, end=' ', file=fd)
                for i in range(sigind_max):
                    print(fdr[i](w), fdi[i](w), '  ', end=' ', file=fd)
                print(file=fd)
            fd.close()

    def CmpEimp(self, which, mu):
        if which==0:
            UXX = self.UAA
        else:
            UXX = self.UAB
        Eimp = zeros((self.Nfunc,self.Nfunc), dtype=complex)
        for iik in range(self.Niik):
          ek = self.Eks[iik,:]-mu
          #             UXX[ifc,iik,ib]*ek[ib] @ UXX.H[ib,iik,ifd]
          Eimp[:,:] += (UXX[:,iik,:]*ek) @ UXX[:,iik,:].conj().T * self.wk[iik]
        Eimp *= 1./self.Niik
        return Eimp
    
    def RotateBasis(self, which, Eimp):
        if which==0:
            UXX = self.UAA
        else:
            UXX = self.UAB
        # diagonalizing Eimp
        (ene, vlp0) = Diagonalize(Eimp) # Note: vlp.H * Eimp * vlp = ene
        # ene = vlp.H * UAA * ek * UAA.H * vlp
        for ik in range(self.Niik):
            UXX[:,ik,:] = dot(vlp0.T.conj(), UXX[:,ik,:])
        return (vlp0, ene)

class N_mu_internal:
    def __init__(self, Ndesired, w, zek, beta):
        self.Ndesired = Ndesired
        self.w = w
        self.zek = zek
        self.zek0 = zek[:,-1,:].real  # shape (Niik,Nband), last frequency point for convergence
        self.Niik = shape(zek)[0]
        self.beta = beta
        # forming mesh of all Matsubara points
        n_max = round((w[-1]*beta/pi-1.)/2.)
        self.w_all = (2*arange(n_max+1)+1)*pi/beta
    def __call__(self, mu):
        Nd0 = np.sum(1.0/(1.0 + np.exp(beta*(self.zek0-mu))))/self.Niik
        iom = self.w*1j + mu
        zsum = np.sum( 1.0/(iom[None,:,None]-self.zek), axis=(0, 2) )
        wsum = np.sum( 1.0/(iom[None,:,None]-self.zek0[:,None,:]),axis=(0, 2))
        Nd = (zsum - wsum).real/self.Niik
        #
        fN = interpolate.UnivariateSpline(self.w, Nd, s=0)
        dn = 2./self.beta * sum(fN(self.w_all)) # sum over matsubara points (positive+negative matsubara give 2)
        Ntot = 2*(Nd0 + dn) # 2 comes from spin degeneracy
        print('mu=', mu, 'Nd=', Ntot, '2*N0=', 2*Nd0, '2*dn=', 2*dn)
        return Ntot-self.Ndesired

class N_mu:
    def __init__(self, Ndesired, w, zek, Eks, wk, beta):
        self.Ndesired = Ndesired
        self.w = w
        #self.wn = w*1j
        self.zek = zek
        self.eks = array(Eks.T, order='F')
        self.wk = wk
        self.beta = beta
        # forming mesh of all Matsubara points
        n_max = round((w[-1]*beta/pi-1.)/2.)
        self.w_all = (2*arange(n_max+1)+1)*pi/beta
    def __call__(self, mu):
        Nd,Nd0 = scc.cmp_dens(mu, self.w, self.zek, self.eks, self.wk, self.beta)
        fN = interpolate.UnivariateSpline(self.w, Nd, s=0)
        #savetxt('dbg.dat', array([self.w,Nd]).T)
        dn = 2./self.beta * sum(fN(self.w_all)) # sum over matsubara points (positive+negative matsubara give 2)
        Ntot = 2*(Nd0 + dn) # 2 comes from spin degeneracy
        print('mu=', mu, 'Nd=', Ntot, '2*N0=', 2*Nd0, '2*dn=', 2*dn)
        return Ntot-self.Ndesired

def Diagonalize(Himp):
    ene,vlp = linalg.eigh(Himp)  # Note: vlp.H * Himp * vlp = ene
    #print 'ene=', ene
    #print 'vlp='
    #PrintM(vlp)
    #print
    #print [abs(vlp[:,j]) for j in range(len(ene))]
    #print

    # sorting eigensystem so that the transformation is close to identity
    av = abs(vlp)
    il=[]
    for i in range(len(ene)):
        avl = av[i,:]
        for j in il: avl[j] = -10000. # set to -infinity, so that we exclude this components, which are already used
        il.append( argmax(avl) )
    #print 'il=', il
    vln = zeros(shape(vlp),dtype=complex)
    enn = zeros(len(ene),dtype=float)
    for j in range(len(il)):
        vln[:,j] = vlp[:,il[j]]
        enn[j] = ene[il[j]]
    return (enn, vln)

        
def broad(w):
      if wbroad*abs(w) < 4.0/unit:
            return wbroad*abs(w)
      else:
            return 4.0/unit
vbroad = vectorize(broad)

def create_log_mesh(nom, nmax, ntail_, beta):
    """Creates logarithmic mesh on Matsubara axis
       Takes first istart points from mesh om and
       the rest of om mesh is replaced by ntail poinst
       redistribued logarithmically.
       Input:
           nom          -- first nom points are all kept
           ntail        -- tail replaced by ntail points only
           higest point -- nmax
       Output:
           som          -- smaller mesh created from big om mesh
           ind_om       -- index array which conatins index to kept Matsubara points
    """
    istart = builtins.min(nom, nmax)
    ntail = builtins.min(ntail_, nmax-istart)
    
    ind_om=[]
    alpha = log((nmax-1.)/istart)/(ntail-1.)
    for i in range(istart):
        ind_om.append(i)
    for i in range(ntail):
        t = int(istart*exp(alpha*i)+0.5)
        if (t != ind_om[-1]):
            ind_om.append(t)

    #ssigdata = zeros( (shape(sigdata)[0], len(ind_om)), dtype=float )
    som = zeros(len(ind_om))
    for i in range(len(ind_om)):
        som[i] = (2*ind_om[i]+1)*pi/beta
    return (som, ind_om)

def PrintM(A):
  
    for i in range(A.shape[0]):
        for j in range(A.shape[1]):
            print("%12.8f %12.8f   " % (A[i,j].real, A[i,j].imag), end=' ')
        print()

def FindSigind(ene,small=1e-4):
    Sigind = zeros((len(ene),len(ene)), dtype=int)
    Sigind[0,0]=1
    ilast=1
    for i in range(len(ene)):
        if (Sigind[i,i]==0):
            ilast+=1
            Sigind[i,i] = ilast
        for j in range(i+1,len(ene)):
            if abs(ene[i]-ene[j])<small:
                Sigind[j,j,]=Sigind[i,i]
    Ene = zeros(ilast)
    for i in range(1,ilast+1):
        inum=0
        for j in range(len(ene)):
            if Sigind[j,j]==i:
                Ene[i-1] += ene[j].real
                inum +=1
        Ene[i-1] /= inum
    return (Sigind, Ene)

def GetEneFromSigind(ene,Sigind,small=1e-4):
    ilast = Sigind.max()
    #print 'ilast=', ilast
    Ene = zeros(ilast)
    for i in range(1,ilast+1):
        inum=0
        for j in range(len(ene)):
            if Sigind[j,j]==i:
                #print 'i=', i, 'j=', j
                Ene[i-1] += ene[j,j].real
                inum +=1
        Ene[i-1] /= inum
    return Ene

#def PrintM(E):
#    for i in range(shape(E)[0]):
        
def GiveTanMesh(x0,L,Nw):
    def fun(x,x0,L,Nw):
        "x[0]=d, x[1]=w"
        d=x[0]
        w=x[1]
        #print 'd=', d, 'w=', w
        return array([L-w/tan(d), x0-w*tan(pi/(2*Nw)-d/Nw) ])
    
    xi=x0/L
    d0 = Nw/2.*(tan(pi/(2*Nw))-sqrt(tan(pi/(2*Nw))**2 - 4*xi/Nw))
    w0 = L*d0

    sol=optimize.root(fun, [d0,w0], args=(x0,L,Nw) )
    (d,w) = sol.x
    om = w*tan(linspace(0.0,1,Nw)*(pi/2-d))
    om[0] = 1e-16
    return om
    
if __name__ == '__main__':
    
    par={'minSigind': 1e-4}
    exec(open("tbg_params.py").read())
    unit = par_bs['unit']
    plot_bands = par_bs['PlotBands']
    wbroad = par['wbroad']
    axis = par['axis']
    beta = par['beta']
    dirs = par['dirs']
    
    # create imp.0 and imp.1 directories
    for d in dirs:
        if not os.path.isdir(d):
            os.mkdir(d)


    Uc = zeros((4,4,4,4))
    for i in range(Uc.shape[0]):
        for j in range(Uc.shape[1]):
            Uc[i,j,j,i] = par['U0']
    
    Uc[0,2,0,2] = par['U0']*par['U_exchange']
    Uc[2,0,2,0] = par['U0']*par['U_exchange']
    Uc[1,3,1,3] = par['U0']*par['U_exchange']
    Uc[3,1,3,1] = par['U0']*par['U_exchange']
    
    em = EmbeddOrbitals(axis, beta, par, plot_bands)

    Siginds=[]
    Eimp=[]
    Sigmas = zeros( (em.Nfunc, len(em.w),2), dtype=complex, order='F')
    for iAA,d in enumerate(dirs):
        eimp = em.CmpEimp(iAA, par['mu'])
        Eimp.append( eimp )
        #print ('eimp=')
        #print(np.array2string(eimp, precision=7, suppress_small=True,max_line_width=200,
        #                      formatter={'complexfloat': lambda x: f"{x.real:+.7f}{x.imag:+.7f}j"}))
        UAA = Uc
        if par['rotateb'][iAA]:
            (vlp, ene) = em.RotateBasis(iAA, eimp)
            if par['rotateU']:
                Unc1 = tensordot(Uc, vlp,axes=(3,0))
                Unc2 = tensordot(conj(vlp),Unc1,axes=(0,0))
                Unc1 = tensordot(Unc2, vlp, axes=(2,0)).transpose((0,1,3,2))
                Unc2 = tensordot(conj(vlp), Unc1, axes=(0,1)).transpose((1,0,2,3))
                UAA = Unc2
        else:
            ene = diagonal(eimp)
        
        print(np.array2string(ene, precision=7, suppress_small=True,max_line_width=200,
                              formatter={'complexfloat': lambda x: f"{x.real:+.7f}{x.imag:+.7f}j"}))
        with open(dirs[iAA]+'/Uc.dat', 'w') as fU1:
            for i0 in range(4):
                for i1 in range(4):
                    for i2 in range(4):
                        for i3 in range(4):
                            if abs(UAA[i0,i1,i2,i3])>1e-4:
                                print("%3d %3d %3d %3d  %17.12f  %12.12f" % (i0, i1, i2, i3, UAA[i0,i1,i2,i3].real, UAA[i0,i1,i2,i3].imag), file=fU1)
                             
        with open(d+'/Trans.dat', 'w') as ftr:
            print(em.Nfunc, em.Nfunc, '# size of Sigind and CF', file=ftr)
            print('#---- Sigind follows', file=ftr)
            if iAA in Sigind0:
                Sigind = diagflat(Sigind0[iAA])
            else:
                Sigind, Ene = FindSigind(ene,par['minSigind'])
            uns = set(Sigind.flatten())
            uns.remove(0)
            #print 'Sigind=', Sigind, uns
            
            Siginds.append( Sigind )
            for i in range(len(Sigind)):
                print('%4d '*len(Sigind) % tuple(Sigind[i,:]), file=ftr)
            print('#---- CF follows', file=ftr)
            Id = eye(len(Sigind))
            for i in range(len(Sigind)):
                for j in range(len(Sigind)):
                    print('%11.8f %11.8f  ' % (Id[i,j], 0.0), end=' ', file=ftr)
                print(file=ftr)
                
        
        sname = d+'/'+par['Sigma']
        sigma_read = False
        if os.path.isfile(sname):
            sdat = loadtxt(sname).T
            if len(sdat)-1 == 2*len(uns): # correct number of columns -- consistent with Sigind
                sigma_read = True
        
            om = sdat[0]                        # first column is frequency mesh
            if (axis=='real' and om[0]>0):
                print('WARNING : Found self-energy on imaginary axis, but have axis=', axis)
                sigma_read = False
            if (axis=='imag' and om[0]<0):
                print('WARNING : Found self-energy on real axis, but have axis=', axis)
                sigma_read = False
        
        sig = zeros( (len(uns),len(em.w)), dtype=complex ) # this is the self-energy we need
        if sigma_read:    # managed to read sigma
            sg = sdat[1::2,:] + sdat[2::2,:]*1j # the rest is put together to complex function
            #print 'Came here', sigma_read
            if sum(abs(om[:em.nom]-em.w[:em.nom]))<1e-3: # Matsubara mesh from the file, and the new generated mesh is the same
                for iw in range(len(em.ind_om)):         # only a subset of points is used here
                    sig[:,iw] = sg[:,em.ind_om[iw]]
            else: # Matsubara mesh has changed, need to interpolate
                print('Now interpolating')
                for icol in range(len(uns)):
                    fsr = interpolate.UnivariateSpline(om, sg[icol].real, s=0)
                    fsi = interpolate.UnivariateSpline(om, sg[icol].imag, s=0)
                    sig[icol,:] = fsr(em.w) + fsi(em.w)*1j
                dd = vstack( (em.w, sig.real, sig.imag) )
                savetxt('debug_Sig', dd.T)
            #print 'sig=', sig
        else:
            print('WARNING : did not manage to read sigma from '+sname+'. Creating empty self-energy')
        # Now taking care of degeneracies.
        for i in range(em.Nfunc):
            Sigmas[i,:,iAA] = sig[Sigind[i,i]-1,:]
        #print 'Sigma=', Sigma

    if axis=='imag':
        if par['recomputeEF']:
            if Fortran:
                mu = em.FermiLevel(par['Nc'], Sigmas, par['mu'], 10./unit)
            else:
                mu = em.FermiLevel_internal(par['Nc'], Sigmas, par['mu'], 10./unit)
            
            open('EF.dat','w').write(str(mu))
        else:
            mu = float(open('EF.dat','r').readline())
    else:
        if os.path.isfile('EF.dat'):
            mu = float(open('EF.dat','r').readline())
        else:
            mu = par['mu']
        print('mu=', mu)

    if Fortran:        
        (Gc, DOS, pDOS, Eimp, Ntot) = em.ComputeGloc(axis, mu, Sigmas)
    else:
        (Gc, DOS, pDOS, Eimp, Ntot) = em.ComputeGloc_internal(axis, mu, Sigmas)
    
    if axis=='imag':
        print('::EF ', mu ,' N= ', Ntot)

    if axis=='real' and par_bs['PlotBands']:
        zek = em.JustComputeEigenvalues(Sigmas, mu)
        savez('zek', zek=zek, w=em.w, mu=mu)
        sys.exit(0)
        
    for iAA,d in enumerate(dirs):
        Ene = GetEneFromSigind(Eimp[iAA],Siginds[iAA],par['minSigind'])
        fe = open(d+'/Eimp','w')
        print(Ene.tolist(), file=fe)
        fe.close()
    
    DltE = em.CmpDelta(axis, Gc, Sigmas, Eimp)
        
    for iAA,d in enumerate(dirs):
        em.PrintGandDelta(axis, dirs[iAA], Siginds[iAA], 'Gloc.inp', 'Gloc.all', Gc[:,:,:,iAA], 'Delta.inp', 'Delta.all', DltE[:,:,:,iAA])
        
    print('Eimp_AA=')
    PrintM(Eimp[0])
    print('Eimp_AB=')
    PrintM(Eimp[1])
    
    if axis=='real':
        em.PrintpDOS('pDOS.out', DOS, pDOS)
    
