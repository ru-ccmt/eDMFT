#!/usr/bin/env python
import sys, os
import numpy as np
from scipy import interpolate
import scipy.linalg as la
import struct
import argparse
import plt_auxiliary as au
from mpi4py import MPI
import time
comm = MPI.COMM_WORLD

import pylab as plt

Ry2eV = 13.60569193

def PrintM(A):
    for i in range(np.shape(A)[0]):
        for j in range(np.shape(A)[1]):
            print(f"({A[i,j].real:12.7f},{A[i,j].imag:12.7f})", end=' ')
        print()
        
def main(args):
    """ Python equivalent of the Fortran program dmftgk """
    mpisize, myrank = comm.size, comm.rank
    #print(f"my rank={myrank} of total {mpisize}")

    Dhybrid = True
    store_Ak = args.store_Ak
    mode = args.mode
    case = au.get_case()
    # read case.indmfl and set imatsubara, gamma, gammac & energy-cutoff from case.indmfl file.
    inl = au.Indmfl(case)
    inl.read()
    imatsubara=inl.matsubara
    gamma, gammac = inl.broadc, inl.broadnc
    emin = inl.om_emin - args.cutoff
    emax = inl.om_emax + args.cutoff
    # is this so run?
    so = ''
    if os.path.isfile(case+".inso") and os.path.getsize(case+".inso")>0 :
        so = 'so'
        if myrank==0:
            print('Found '+case+'.inso file, hence assuming so-coupling exists. Switching -so switch!')
    # what are name of the input files?
    fenergy=case+'.energy'+so
    frotlm=case+'.rotlm'
    fUdmft='Udmft.0'
    # names of the output files
    feigenvals = 'eigenvalues.dat'
    fUR = 'UR.dat'
    fUL = 'UL.dat'
    
    
    (nat,iso,norbitals,nindo,ncix,cixdim,natom,nkpt,nl,ll_arr,iorbital,cix,nind,
     csize,iSx,Sigind,cix_orb,EF,VOL,posc,nmat,nume,Qcomplex,lmax2,maxdim2,maxdim,maxsize) = ReadBasicArrays('BasicArrays.dat',myrank)
    
    # going over possible klists and determine which one is compatible with nkp from projector
    klists = [case+'.klist', case+'.klist_band']
    fklist=None
    for ifi,fname in enumerate(klists):
        if os.path.isfile(fname):
            with open(fname, 'r') as fi:
                lines = fi.readlines()
                for nn,line in enumerate(lines):
                    if line[:3]=='END':
                        break
                if nn==nkpt:
                    fklist = fname
                    break
    band=''
    if ifi==1: band='_band'
    fsigname = ['sig.inp'+str(icix+1)+args.m_ext+band for icix in range(ncix)]
    if band:
        store_Ak = True
        if myrank==0:
            print('Band calculation hence switching on store_Ak.')
    if myrank==0:
        print("mode=", mode)
        print("imatsubara=", imatsubara)
        print("fenergy=", fenergy)
        print("fklist=", fklist)
        print("frotlm=", frotlm)
        print("fUdmft=", fUdmft)
        print("gamma=", gamma)
        print("gammac=", gammac)
        print('fsigname=', fsigname)
        print("feigenvals", feigenvals)
        print("fUR=", fUR)
        print("fUL=", fUL)
        print("emin,emax=", emin, emax)

        
    # Reading k-points from case.klist or case.klist_band
    kvec, wgh = ReadKlist(fklist)
    tweight = np.sum(wgh)
    if myrank==0:
        print('tweight=', tweight)
    
    # output filenames for mode='g'
    fgkout    = ['Gk.'+str(icix+1) for icix in range(ncix)]
    fglocout  = ['Glocal.'+str(icix+1) for icix in range(ncix)]
    fgloc     = ['g_local.'+str(icix+1) for icix in range(ncix)]
    
    # here we decided to create combined index icx_it so that all needed local indices are just one simple vector.
    icx_it=[[[] for it in range(csize[icix])] for icix in range(ncix)]
    nii = sum([csize[icix] for icix in range(ncix)])
    ii=0
    for icix in range(ncix):
        for it in range(csize[icix]):
            icx_it[icix][it] = ii
            ii+=1
    
    norbs = sum(cixdim)
    if myrank==0:
        print('nii=', nii, 'icx_it=', icx_it)
        print('nind=', nind, 'norbs=', norbs)
        print('cixdim=', cixdim)


    small_indx = [{} for icix in range(ncix)]
    for icix in range(ncix):
        for ind1 in range(cixdim[icix]):
            for ind2 in range(cixdim[icix]):
                it = abs(Sigind[icix][ind1,ind2])
                if it!=0:
                    ii = icx_it[icix][it-1]
                    small_indx[icix][(ind1,ind2)]=ii

    #for icix in range(ncix):
    #    for (ind1,ind2) in small_indx[icix]:
    #        ii = small_indx[icix][(ind1,ind2)]
    #        print('icix=', icix, '('+str(ind1)+','+str(ind2)+')','ii=', ii)
    
    # Reading the self-energy for all orbitals
    sig0 = np.loadtxt(fsigname[0])
    (nom,norb_2) = np.shape(sig0)
    sigma = np.zeros((nom,nii),dtype=np.complex128)
    # Read each correlated sub-block from fsigname(icix):
    for icix in range(ncix):
        dat = np.loadtxt(fsigname[icix]).T
        omega = dat[0]
        sig = dat[1::2]+dat[2::2]*1j
        for it in range(csize[icix]):
            sigma[:,icx_it[icix][it]] = sig[it,:]
    
    if args.interp:
        # Interpolates sigmas on equidistant mesh!
        eq_om = np.linspace(inl.om_emin, inl.om_emax, inl.om_npts)  # Equidistant mesh
        sig_band=np.zeros((len(eq_om),nii), dtype=np.complex128)
        for ii in range(nii):
            fr = interpolate.UnivariateSpline(omega, sigma[:,ii].real,s=0)
            fi = interpolate.UnivariateSpline(omega, sigma[:,ii].imag,s=0)
            sig_band[:,ii] = fr(eq_om)+fi(eq_om)*1j
        sigma = sig_band
        omega = eq_om
        nom = len(eq_om)
        
    iaib,degs =  Find_iaib(Sigind, icx_it, cixdim, nii)
    if myrank==0:
        print('iaib=', iaib)
        print('degs=', degs)
    
    # Reading projector from Udmft.myrank
    fh_p = open(fUdmft[:-1]+str(myrank), 'rb')
    data_bytes = FortranReadRecord(fh_p, fUdmft)
    numk, nsymop, tnorbitals = struct.unpack("3i", data_bytes)
    data_bytes = FortranReadRecord(fh_p, fUdmft)
    #print(len(data_bytes))
    tnindo = struct.unpack(str(norbitals)+"i", data_bytes)
    if not np.allclose(tnindo,nindo):
        print('WARNING: when reading', fUdmft, 'in binary mode tnindo=', tnindo, 'while indo=', nindo)
    
    pr_proc  = int(nkpt/mpisize+0.999)
    print(f"rank={myrank} pr_procr={numk} pr_proc={pr_proc} nkpt={nkpt} nsymop={nsymop} tnorbitals={tnorbitals}")

    
    # reading case.energy file on all cores
    fh_ene = open(fenergy, 'r')
    (Eks, kps, wgs) = ReadEnergyFile(fenergy, nat, nkpt)
    #print('len(Eks)=', len(Eks), 'nkp=', nkpt)
    if not np.allclose(wgh,wgs):
         print("ERROR: mismatch in k-point weight", wgs, wgh)
    
    tw_new=0
    if mode=='g':
        # local greens function
        gmloc = np.zeros((nom,nii), dtype=np.complex128)
        if store_Ak:
            Akom = np.zeros((numk,nom,nii+1), dtype=float)
            gtot = np.zeros(nom, dtype=np.complex128)
        Eimp  = np.zeros(nii, dtype=np.complex128)
        for iikp in range(numk): # loop over k-points
            ikp = iikp + myrank*pr_proc # real index of k-point
            Ek = Eks[ikp]*Ry2eV         # Convert Ek from Ry to eV
            # reading projector outside
            data_bytes = FortranReadRecord(fh_p, fUdmft)
            jjkp, nbands, tmaxdim2, tnorbitals, nemin = struct.unpack("5i", data_bytes)
            nemax = nemin + nbands - 1
            #print(f'ikp={ikp} jjkp={jjkp} nbands={nbands} nemin={nemin}')
            e_ks = Ek[nemin-1:nemin-1+nbands]-EF
            UDMFT = np.zeros((norbs,nbands),dtype=np.complex128)
            gk    = np.zeros((nbands, nbands), dtype=np.complex128)
            for isym in range(nsymop):  # For each symmetry op
                data_bytes = FortranReadRecord(fh_p, fUdmft)
                isym2 = struct.unpack("i", data_bytes)
                kweight = wgh[ikp]/(nsymop*tweight)
                DMFTU = ReadDMFTUi(fh_p, fUdmft, nindo, nbands)
                # In dmft1 we still write out projection in non-optimal way, projection to atom,l,orbitals_on_atom_l
                # For cluster calculation, it is better to combined projection that corresponds to the same icix, i.e., cluster
                # We therefore here create UDMFT, which combines blocks for the same cluster into a unit of dimension cixdim[icix]
                # We still take into account interference between different clusters in the unit cell because
                # G_{icix=1,icix=2} could be computed in principle, and is nonzero. In practice we don't compute it due to effieincy reasons.
                # But we do not neglect such terms, namely, H11 below contains all clusters (icix=1,icix=2,...) and H22 contains degrees of freedom that
                # have absolutely no overlap with any cluster (through SVD decomposition). The inversion is than G11 = 1/(om-H11-Sigma_11-H12*1/(om-H22)*H21)
                # If we neglected interference between clusters, we could invert projection to each icix separately, which would be an approximation.
                for iorb in range(len(nindo)):
                    icix = cix_orb[iorb]
                    pre_cix=sum(cixdim[:icix])
                    iind = [iSx[iorb][ind]-1+pre_cix for ind in range(nindo[iorb])]
                    for ind in range(nindo[iorb]):
                        UDMFT[pre_cix+iSx[iorb][ind]-1,:] = DMFTU[iorb][ind][:]
                del DMFTU # done with old projector, from now on use the new UDMFT
                
                t0 = time.time()
                if not args.full_inversion:
                    # Create projector in matrix form for SVD
                    UF = UDMFT.T  # (nbands,nii)
                    nlow = np.shape(UF)[1]  # number of finite singular values is equal to the number of orbitals < nbands
                    # SVD to integrate out the unnecesary matrix part
                    U, S, Vh = np.linalg.svd(UF, full_matrices=True)
                    Unew = S[:,np.newaxis]*Vh  # new projector in rotated basis, taking out unitary matrix U
                    # efficient projection 
                    STrans = np.zeros((nii,nlow,nlow),dtype=np.complex128) # projector for self-energy using Unew
                    for ii in range(len(iaib)):
                        for ia,ib in iaib[ii]:
                            STrans[ii,:,:] += np.outer(Unew[:,ia],Unew[:,ib].conj())
                    # projection to the Green's function when we don't need complete matrix form, but only components which are allowed in Sigind matrix.
                    GTrans = np.reshape(STrans.conj(),(nii,nlow*nlow))
                    # Hc = U^+ @ e_ks @ U is the Kohn-Sham matrix in rotated basis in which self-energy needs to be added only to 11 block,
                    # and is completely absent in 12,21, and 22 block.
                    Hc = U.conj().T @ (e_ks[:,np.newaxis] * U)
                    ekl, v = np.linalg.eigh(Hc[nlow:,nlow:])  # eigensystem of H22 to compute 1/(om-H22) = A*1/(om-ek)*A^+
                    g12 = -Hc[:nlow,nlow:] #
                    gi12 = g12@v           # H12*A, because we need H12*1/(om-H22)*H21=H12*A*1/(om-ek)*(H12*A)^H
                    Hlow = Hc[:nlow,:nlow] # H11
                    if store_Ak:
                        H22 = Hc[nlow:,nlow:]
                        H12 = Hc[:nlow,nlow:]
                        H21 = Hc[nlow:,:nlow]
                    Eimp += (GTrans @ Hlow.reshape((nlow*nlow)))*kweight # impurity levels
                    # loop over frequency
                    for iom in range(nom):
                        xomega = omega[iom]*1j+gamma*1j if imatsubara else omega[iom]+gamma*1j
                        Sigma = np.einsum('i,ijk->jk', sigma[iom, :], STrans) # einstein summation : Sigma[j,k] = \sum_i sigma[iom,i]*STrans[i,j,k]
                        g11 = xomega*np.eye(nlow)-Sigma-Hlow                  # om-Sigma_11-H_11
                        gi11 = np.linalg.inv( g11 - (gi12 * (1/(xomega-ekl))) @ gi12.T.conj() ) # 1/(om-Sigma_11-H_11-H12*1/(om-H22)*H21)
                        gmk = GTrans @ gi11.reshape((nlow*nlow))
                        gmk *= 1./degs
                        gmloc[iom,:] += gmk*kweight
                        if store_Ak: # This is used for plotting with plt_sakplot.py
                            Akom[iikp,iom,1:] -= gmk.imag/(np.pi*nsymop)
                            # We store spectral function projected to orbitals in row 1...nii,
                            # while in row 0 we store the total spectral function Tr(1/(om+mu-ek-Sigma)).
                            # The latter can be computed as Tr(G_{11})+Tr(G_{22}). Note that computing G_22 is more
                            # expensive that computing G_11, hence we want to avoid this step if not absolutely necessary.
                            g22 = xomega*np.eye(nbands-nlow) - H22 - H21@np.linalg.inv(g11)@H12
                            gi22 = np.linalg.inv(g22)
                            gtc = np.trace(gi11)+np.trace(gi22)+np.sum(1/(xomega-Ek[:nemin-1]+EF+gamma*4j))+np.sum(1/(xomega-Ek[nemin-1+nbands:]+EF+gamma*4j))
                            Akom[iikp,iom,0] -= gtc.imag/(np.pi*nsymop)
                            gtot[iom] += gtc*kweight
                else: # Here we don't use SVD trick, hence inversion of the matrix of nbands x nbands size. Conceptually simpler, but more expensive.
                    STrans = np.zeros((nii,nbands,nbands),dtype=np.complex128)
                    for ii in range(len(iaib)):
                        for ia,ib in iaib[ii]:
                            STrans[ii,:,:] += np.outer(UDMFT[ia,:],UDMFT[ib,:].conj()) 
                    GTrans = np.reshape(STrans.conj(),(nii,nbands*nbands))
                    Eimp += (GTrans @ e_ks.reshape((nbands*nbands)))*kweight
                    # loop over frequency
                    for iom in range(nom):
                        xomega = omega[iom]*1j+gamma*1j if imatsubara else omega[iom]+gamma*1j
                        Sigma = np.einsum('i,ijk->jk', sigma[iom, :], STrans)  # einstein summation : Sigma[j,k] = \sum_i sigma[iom,i]*STrans[i,j,k]
                        gk = np.diag(xomega-e_ks) - Sigma
                        gk = np.linalg.inv(gk)
                        gmk = GTrans @ gk.reshape( (nbands*nbands) )
                        gmk *= 1./degs
                        gmloc[iom,:] += gmk*kweight
                        if store_Ak: # This is used for plotting with plt_sakplot.py
                            Akom[iikp,iom,1:] -= gmk.imag/(np.pi*nsymop)
                            # We store spectral function projected to orbitals in row 1...nii,
                            # while in row 0 we store the total spectral function Tr(1/(om+mu-ek-Sigma)).
                            gtc = np.trace(gk)+np.sum(1/(xomega-Ek[:nemin-1]+EF+gamma*4j))+np.sum(1/(xomega-Ek[nemin-1+nbands:]+EF+gamma*4j))
                            # Note that we should add bands which were removed from consideration in DMFT
                            Akom[iikp,iom,0] -= gtc.imag/(np.pi*nsymop)
                            gtot[iom] += gtc*kweight
                            
                t2 = time.time()
                tw_new += t2-t0

        if myrank==0:
            print('Time ', tw_new, 's')
            gmsum = np.empty_like(gmloc)
            Esum = np.empty_like(Eimp)
            if store_Ak: gsm = np.empty_like(gtot)
        else:
            gmsum = None
            Esum = None
            gsm = None
        comm.Reduce(gmloc, gmsum, op=MPI.SUM, root=0)
        comm.Reduce(Eimp, Esum, op=MPI.SUM, root=0)
        if store_Ak: comm.Reduce(gtot, gsm, op=MPI.SUM, root=0)
        if myrank==0:
            Eimp = Esum
            gmloc = gmsum
            if store_Ak: gtot = gsm
            #print('Eimp=', Esum.real)
            for icix in range(ncix):
                with open(fglocout[icix],'w') as fo:
                    for iom in range(len(omega)):
                        print(omega[iom], end=' ', file=fo)
                        #if store_Ak:
                        #    print(gtot[iom].real, gtot[iom].imag, end=' ', file=fo)
                        for it in range(csize[icix]):
                            gc = gmloc[iom,icx_it[icix][it]]
                            print(gc.real, gc.imag, end=' ', file=fo)
                        print(file=fo)
                        
            if Dhybrid:
                delta = np.zeros_like(gmloc)
                for icix in range(ncix):
                    glc = np.zeros((cixdim[icix],cixdim[icix]),dtype=complex)
                    sig = np.zeros((cixdim[icix],cixdim[icix]),dtype=complex)
                    eim = np.zeros((cixdim[icix],cixdim[icix]),dtype=complex)
                    dlt = np.zeros((cixdim[icix],cixdim[icix]),dtype=complex)
                    for iom in range(len(omega)):
                        xomega = omega[iom]*1j+gamma*1j if imatsubara else omega[iom]+gamma*1j
                        for (ind1,ind2) in small_indx[icix]:
                            ii = small_indx[icix][(ind1,ind2)]
                            glc[ind1,ind2] = gmloc[iom,ii]
                            sig[ind1,ind2] = sigma[iom,ii]
                            eim[ind1,ind2] = Esum[ii]
                        ginv = np.linalg.inv(glc)
                        dlt = xomega*np.eye(cixdim[icix]) - eim - sig - ginv
                        for (ind1,ind2) in small_indx[icix]:
                            ii = small_indx[icix][(ind1,ind2)]
                            delta[iom,ii] += dlt[ind1,ind2]/degs[ii]
                            
                for icix in range(ncix):
                    with open('Dlt.'+str(icix+1),'w') as fo:
                        for iom in range(len(omega)):
                            print(omega[iom], end=' ', file=fo)
                            for it in range(csize[icix]):
                                dl = delta[iom,icx_it[icix][it]]
                                print(dl.real, dl.imag, end=' ', file=fo)
                            print(file=fo)
            
        if store_Ak:
            if myrank==0:
                recv_counts = np.array([min(max(0,nkpt-pr_proc*(iproc)),pr_proc) for iproc in range(mpisize)],dtype=int)
                displacements = np.cumsum(recv_counts) - recv_counts # Starting position for each process' data
                recv_counts*=nom*(nii+1)
                displacements*=nom*(nii+1)
                A_kom = np.zeros((nkpt,nom,nii+1),dtype=float)
            else:
                recv_counts = None
                displacements = None
                A_kom = np.empty(1)
            
            comm.Gatherv(sendbuf=Akom.ravel(), recvbuf=(A_kom.ravel(), recv_counts, displacements, MPI.DOUBLE), root=0)
            if myrank==0:
                np.save('A_kom.npy',A_kom)
                np.save('omega.npy',omega)
                #for ik in range(nkpt):
                #    with open('Ak.'+str(ik),'w') as fo:
                #        for iom in range(len(omega)):
                #            print(omega[iom], end=' ', file=fo)
                #            print(('{:12.7f} '*(nii+1)).format(*(A_kom[ik,iom,:])), file=fo)
                        
        
    elif mode=='e':
        zek    = [[] for iikp in range(numk)]
        nemins = np.zeros(numk, dtype=int)
        for iikp in range(numk): # loop over k-points
            ikp = iikp + myrank*pr_proc # real index of k-point
            Ek = Eks[ikp]*Ry2eV         # Convert Ek from Ry to eV
            # reading projector outside
            data_bytes = FortranReadRecord(fh_p, fUdmft)
            jjkp, nbands, tmaxdim2, tnorbitals, nemin = struct.unpack("5i", data_bytes)
            nemax = nemin + nbands - 1
            #print(f'ikp={ikp} jjkp={jjkp} nbands={nbands} nemin={nemin}')
            e_ks = Ek[nemin-1:nemin-1+nbands]-EF
            UDMFT = np.zeros((norbs,nbands),dtype=np.complex128)
            nemins[iikp] = nemin
            zeks = np.zeros((nom,nbands),dtype=np.complex128)
            for isym in range(nsymop): # For each symmetry op
                data_bytes = FortranReadRecord(fh_p, fUdmft)
                isym2 = struct.unpack("i", data_bytes)
                kweight = wgh[ikp]/(nsymop*tweight)
                DMFTU = ReadDMFTUi(fh_p, fUdmft, nindo, nbands)
                for iorb in range(len(nindo)):
                    icix = cix_orb[iorb]
                    pre_cix=sum(cixdim[:icix])
                    iind = [iSx[iorb][ind]-1+pre_cix for ind in range(nindo[iorb])]
                    for ind in range(nindo[iorb]):
                        UDMFT[pre_cix+iSx[iorb][ind]-1,:] = DMFTU[iorb][ind][:]
                del DMFTU
                
                t0 = time.time()
                STrans = np.zeros((nii,nbands,nbands),dtype=np.complex128)
                for ii in range(len(iaib)):
                    for ia,ib in iaib[ii]:
                        STrans[ii,:,:] += np.outer(UDMFT[ia,:],UDMFT[ib,:].conj()) 
                # loop over frequency
                for iom in range(nom):
                    #xomega = omega[iom]*1j+gamma*1j if imatsubara else omega[iom]+gamma*1j
                    Sigma = np.einsum('i,ijk->jk', sigma[iom, :], STrans)
                    hs = np.diag(e_ks) + Sigma
                    _zek_, Al, Ar = la.eig(hs, left=True, right=True)
                    indx = sorted(range(nbands), key=lambda k: _zek_[k].real)
                    _zek_ = _zek_[indx]
                    Al = Al[:,indx]
                    Ar = Ar[:,indx]
                    #print('diagr works=', np.allclose(hs @ Ar,  Ar @ np.diag(_zek_) ) )
                    #print('diagl works=', np.allclose(Al.T.conj() @ hs,  np.diag(_zek_) @ Al.T.conj() ) )
                    #print('diag works=', np.allclose(Al.conj().T @ hs @ Ar , np.diag(_zek_), atol=1e-3, rtol=1e-3) )
                    # Uk^+ 1/(om+mu-ek-Sig) Uk = Uk^+ * Ar * 1/(om-ek) * Al^+ Uk
                    UAr = UDMFT.conj() @ Ar
                    UAl = Al.T.conj() @ UDMFT.T
                    zeks[iom,:] += _zek_ * 1./nsymop
                    
                t2 = time.time()
                tw_new += t2-t0
            if myrank==0:
                print('ikp=', ikp, 'done')
            min_zek = np.min(zeks.real,axis=0)
            max_zek = np.max(zeks.real,axis=0)
            istart = np.where(max_zek>emin)[0][0]
            iend = np.where(min_zek<emax)[0][-1]+1
            zek[iikp] = zeks[:,istart:iend]
            nemins[iikp] = nemin+istart
            
        all_zek = comm.gather(zek, root=0)
        all_nemins = comm.gather(nemins, root=0)
        if myrank==0:
            with open('eigenvalues.dat', 'w') as fe:
                for proc, proc_zek in enumerate(all_zek):
                    #print('proc=', proc, 'len(zek)=', len(proc_zek))
                    for iikp,zeks in enumerate(proc_zek):
                        ikp = iikp + pr_proc*myrank
                        nemin = all_nemins[proc][iikp]
                        #print(ikp,'shape(zek)=', np.shape(zeks))
                        kweight = wgh[ikp]/(nsymop*tweight)
                        nom,nbands = np.shape(zeks)
                        #nemin = nemins[ikp]
                        isym=1
                        nxmin,nxmax=1,nbands
                        print(f'! {ikp} {isym} {nbands} {nemin} {nom} {kweight} : ikp, isym, nbands nemin, nomega, kweight', file=fe)
                        for iom in range(len(omega)):
                            print(f'{omega[iom]:19.14f}  {nxmax-nxmin+1} {nxmin+nemin-1}', end=' ', file=fe)
                            for iband in range(nbands):
                                print(f'{zeks[iom,iband].real+EF:24.15f} {zeks[iom,iband].imag:24.15f}', end=' ', file=fe)
                            print(file=fe)
                    
                    # '(A,5I5,2x,F23.20,2x,A)') '!', ikp, isym, nbands, nemin, nom, wgh(ikp)/(nsymop*tweight), ': ikp, isym, nbands nemin, nomega, kweight'
                    #WRITE(fh_ee+0,'(F19.14,2x,2I5,1x)',advance='no') omega(iom), nxmax-nxmin+1, nxmin+nemin-1
                    #WRITE(fh_ee+0,'(2E24.16,1x)',advance='no') zek(ibnd)

    if myrank==0:
        print("STOP dmftgk")

    
def ReadKlist(fklist):
    """
    Reads k-list of wien2k
    Parameters:
      fklist : str
          The filename of the k-point list.
    Returns:
      kvec : np.ndarray of shape (4, nkpt), dtype=int
          Array of k-vector components.
      wgh : np.ndarray of shape (nkpt,), dtype=float
          Array of k-point weights.
    """
    # Initialize arrays to hold the results.
    kvec=[]
    wgh=[]
    with open(fklist, 'r') as f:
        line = f.readline()
        if line[14:15]==' ':
            # newform
            w = [10,20,30,40,50,55]
        else:
            w = [10,15,20,25,30,34]
        kvec.append( [int(line[w[i]:w[i+1]]) for i in range(4)] )
        wgh.append( float(line[w[4]:w[5]]) )
        for ik in range(1000000):
            line = f.readline()
            if not line:
                print("ERROR1 in reading", fklist)
                break
            KNAME = line[:w[0]]
            if KNAME.strip() == "END":
                break
            kvec.append( [int(line[w[i]:w[i+1]]) for i in range(4)] )
            wgh.append( float(line[w[4]:w[5]]) )
            
    kvec = np.array(kvec).T
    wgh = np.array(wgh)
    return kvec, wgh

def FortranReadRecord(fh_p, fname):
    # Fortran reads or writes "numk, nsymop, tnorbitals, (tnindo(iorb), iorb=1,norbitals)"
    rec_length_bytes = fh_p.read(4)
    rec_length = struct.unpack("i", rec_length_bytes)[0]  # 'i' is a 4-byte integer
    # Read the data
    data_bytes = fh_p.read(rec_length)
    trailing_rec_length = struct.unpack("i", fh_p.read(4))[0]
    if rec_length != trailing_rec_length:
        raise ValueError("Record length mismatch in reading binary file "+fname)
    return data_bytes

def ReadDMFTUi(fh_p, fUdmft, nindo, nbands):
    DMFTU=[]
    for iorb in range(len(nindo)):
        _DMFTU_=[]
        for ind in range(nindo[iorb]):
            data_bytes = FortranReadRecord(fh_p, fUdmft)
            dmftu = np.array( struct.unpack(str(2*nbands)+"d", data_bytes) )
            _DMFTU_.append(dmftu[::2]+dmftu[1::2]*1j)
        DMFTU.append(_DMFTU_)
    return DMFTU
    
 
def ReadBasicArrays(fbasicArrays,myrank):
    # --- Read the header of BasicArrays.dat to get dimensions ---
    with open(fbasicArrays, 'r') as fhb:
        lines_hb = [line.strip() for line in fhb]
    nat, iso, norbitals, ncix, natom = (int(x) for x in lines_hb[1].split()[:5])
    nkpt, nmat, nume = (int(x) for x in lines_hb[3].split()[:3])
    Qcomplex_str = (lines_hb[5].split()[0]).lower()
    Qcomplex = True if Qcomplex_str in ['t', 'true', '.true.'] else False
    lmax2, maxdim2, maxdim, maxsize = (int(x) for x in lines_hb[7].split()[:4])

    # --- Read header "nindo" and then nindo values ---
    tokens = lines_hb[9].split()
    if len(tokens) < norbitals:
        raise ValueError("Not enough values for nindo")
    nindo = np.array(tokens[:norbitals], dtype=int)
    # --- Read header "cidim" (cixdim) ---
    tokens = lines_hb[11].split()
    if len(tokens) < ncix:
        raise ValueError("Not enough values for cixdim")
    cixdim = np.array(tokens[:ncix], dtype=int)
    # --- Read header "nl" ---
    tokens = lines_hb[13].split()
    if len(tokens) < natom:
        raise ValueError("Not enough values for nl")
    nl = np.array(tokens[:natom], dtype=int)
    # Allocate arrays for ll, iorbital, cix, nind.
    ll_arr      = [[] for i in range(natom)]
    iorbital    = [[] for i in range(natom)]
    cix         = [[] for i in range(natom)]
    nind        = [[] for i in range(natom)]
    # Loop over atoms (icase) and over lcase = 1 .. nl(icase)
    ipos = 15
    for icase in range(natom):
        for lcase in range(nl[icase]):
            # Each value is on a separate line.
            ll_arr  [icase].append( int( lines_hb[ipos] ) ); ipos+=1
            iorbital[icase].append( int( lines_hb[ipos] ) ); ipos+=1
            cix     [icase].append( int( lines_hb[ipos] ) ); ipos+=1
            nind    [icase].append( int( lines_hb[ipos] ) ); ipos+=1
    ipos = 15+4*natom+1
    # --- Read header "csize" and then csize values ---
    tokens = lines_hb[ipos].split(); ipos+=1
    if len(tokens) < ncix:
        raise ValueError("Not enough values for csize")
    csize = np.array(tokens[:ncix], dtype=int)
    ipos+=1
    # --- Read header "iSx" and then iSx values ---
    iSx = []
    for iorb in range(norbitals):
        # For each orbital, read nindo[iorb] integers.
        _iSx_=[]
        for ip1 in range(nindo[iorb]):
            _iSx_.append( int(lines_hb[ipos]) )
            ipos+=1
        iSx.append(_iSx_)

    ipos+=1
    # --- Read header "Sigind" and then Sigind values ---
    # Sigind has shape [ncix][dim, dim]. For each icix=0..ncix-1, we read a square block.
    Sigind = [] #np.zeros((maxdim, maxdim, ncix), dtype=int)
    for icix in range(ncix):
        _Sigind_ = np.zeros((cixdim[icix],cixdim[icix]),dtype=int)
        for ip in range(cixdim[icix]):
            for iq in range(cixdim[icix]):
                _Sigind_[ip,iq] = int(lines_hb[ipos])
                ipos+=1
        Sigind.append(_Sigind_)

    ipos+=1
    # --- Read header "EF" and then EF and VOL ---
    c = [float(x) for x in lines_hb[ipos].split()]
    EF, VOL = c[0], c[1]
    EF *= Ry2eV
    ipos+=1
    # --- Read posc for each orbital ---
    posc = [] #np.zeros((norbitals,3), dtype=float)
    for iorb in range(norbitals):
        posc.append( np.array([float(x) for x in lines_hb[ipos].split()]) )

    cix_orb = np.zeros(norbitals,dtype=int)
    for icase in range(natom):
        for lcase in range(nl[icase]):
            # Each value is on a separate line.
            iorb=iorbital[icase][lcase]
            icix=     cix[icase][lcase]
            #          nind[icase][lcase]
            cix_orb[iorb-1]=icix-1
        
    QPrint=False
    if QPrint and myrank==0:
        print("nindo =", nindo)
        print("cixdim =", cixdim)
        print("nl =", nl)
        print("ll =", ll_arr)
        print("iorbital =", iorbital)
        print("cix =", cix)
        print("nind =", nind)
        print("csize =", csize)
        print("iSx =", iSx)
        #print("Sigind =", Sigind)
        print("EF =", EF)
        print("VOL =", VOL)
        print("posc =", posc)
        print('cix_orb=', cix_orb)
    # Up to here in standalone function ReadBasicArray
    #------------------------------------------------------------------------------------------------
    return (nat,iso,norbitals,nindo,ncix,cixdim,natom,nkpt,nl,ll_arr,iorbital,cix,nind,
            csize,iSx,Sigind,cix_orb,EF,VOL,posc,nmat,nume,Qcomplex,lmax2,maxdim2,maxdim,maxsize)

def ReadEnergyFile(fenergy, nat, nkpt):
    """
    Reads band energies from an energy file.
    The energy file is expected to have a header line with fixed-width fields:
      - Three floating-point numbers in E19.12 format (each occupying 19 characters),
      - A 10-character string,
      - Two integers (each in a 6-character field),
      - One floating-point number in E19.12 format.
      
    After the header, there are NE lines, each containing:
      - NUM (an integer band index, assumed to be 1-indexed),
      - E1   (a floating-point energy).

    Parameters:
      filename : str
          Path to the energy file.
      nume : int
          Total number of energy levels (size of the Ek array).

    Returns:
      kp : numpy.ndarray of shape (3,)
          The three k-point coordinates.
      wg : float
          The weight (WG) read from the header.
      NE : int
          The number of energy lines read.
      Ek : numpy.ndarray of shape (nume,)
          The energy array with energies placed at positions (NUM-1).  
          (All values are initially zero except for the ones read from the file.)
    """
    with open(fenergy, 'r') as fh_ene:
        data = fh_ene.readlines()
        #data = fh_ene.read().splitlines()
    
    # In Fortran, read linearization energies for i=1..nat, skipping them
    ipos = 2*nat
    Eks=[]
    kps=[]
    wgs=[]
    #print('nat=', nat, 'nkpt=', nkpt)
    for ikp in range(nkpt):
        # Read the header line.
        if ipos>=len(data):
            print('ERROR: Energy file '+fenergy+' is too short. The #-k points in projector is',nkpt,' while energy file has ',ikp )
            sys.exit(1)
        header = data[ipos]
        ipos+=1
        w = [0,19,38,57,67,73,79,98]
        #if len(header) < w[-1]:
        #    raise ValueError("Header line is shorter than expected.")
        # Extract the three k-point coordinates.
        kp = np.array([float(header[w[i]:w[i+1]].strip()) for i in range(3)])
        KNAME = header[w[3]:w[4]].strip()
        # Extract N and Ne and weight
        N, Ne, wg = int(header[w[4]:w[5]].strip()), int(header[w[5]:w[6]].strip()), float(header[w[6]:w[7]])
        kps.append(kp)
        wgs.append(wg)
        # Loop over Ne lines, each containing NUM and E1.
        Ek = np.zeros(Ne)
        for ie in range(Ne):
            line = data[ipos]; ipos+=1
            parts = line.split()
            num,Ek[ie] = int(parts[0]), float(parts[1])
            if (num-1!=ie):
                raise IndexError(f"Band index {num} out of bounds (1..{Ne}).")
        Eks.append(Ek)
        #print('k=', kp, 'Ne=', Ne, 'Ek=', Ek)
    return Eks, kps, wgs


def Find_iaib(Sigind, icx_it, cixdim, nii):
    iaib = [[] for ic in range(nii)]
    degs = np.zeros(nii,dtype=int)
    for icix in range(len(cixdim)):
        pre_cix = sum(cixdim[:icix])
        for ind1 in range(cixdim[icix]):
            for ind2 in range(cixdim[icix]):
                it = abs(Sigind[icix][ind1,ind2])
                if it!=0:
                    ii = icx_it[icix][it-1]
                    ia = pre_cix+ind1
                    ib = pre_cix+ind2
                    #print('icix=', icix, 'ind1=', ind1, 'ind2=', ind2, 'it=', it, 'ia,ib=', ia,ib)
                    iaib[ii].append( (ia,ib) )
                    degs[ii] += 1
    return (iaib,degs)


if __name__ == "__main__":
    usage = """Recomputes Greens function in mode g or complex eigenvaues/eigenvectors in mode e

 We are computing local Green's function using the numerical trick to speed up the inversion of the matrix.
 The matrix to invert is of the size norbs x norbs rather than nbans x nbands.
 To compute the local Green's function using projector Uk, we normally do the following:
 
    G =  Uk^+ * (om + mu - ek - Uk*Sigma*Uk^+ )^{-1} * Uk
 
 here Uk is projector/embeddor and Sigma is defined in orbital basis, while ek is a 
 Kohn-Sham digonal matrix of band eigenvalues. Here Uk(:nbands,:norbs).
 The idea is to perform SVD decomposition of the projector Uk = U*s*Vt, and 
 U(:nbands,:nbands) is a unitary matrix in large bands space, and Vt(:norbs,:norbs) 
 is unitary in the smaller orbital space. We have only s(:norbs) finite singular values.
 It is straighforward to rewrite the above equation into
 
   G =  Vt^+ * s * (om + mu - U^+*ek*U - s*Vt*Sigma*Vt^+*s)^{-1} * s*Vt
 
 We now recognize that we only need to compute the inverse in the block where singular values 
 s are finite. The rest of the matrix does not contribute to the projected Green's function. 
 It does contribute to the total DOS, or, partial DOS for orbitals which are not correlated 
 (in the projector). Hence, when any icix=0, we need to invert the entire matrix (full_inversion=.True.) 
 But for self-consistent cDMFT, we can invert only a submatrix in the following way:
 
  H11 = U^+[:norbs,:]*(ek-mu)*U[:,:norbs]
  H12 = U^+[:norbs,:]*(ek-mu)*U[:,norbs:]
  H22 = U^+[norbs:,:]*(ek-mu)*U[:,norbs:]
 
  G =         (  om - H11 - s*Vt*Sigma*Vt^+*s,  -H12 )^{-1}
       Vt^+*s*(  -H21                        , om-H22)      * s*Vt
 
  G = Vt^+*s* ( om - H11 - s*Vt*Sigma*Vt^+*s - H12 * 1/(om-H22) * H21 )^{-1} * s*Vt
 
  We can also diagonalize Hermitian matrix H22, i.e., H22 = A*ekp*A^+, and than we can write:
 
  G = Vt^+*s* ( om - H11 - s*Vt*Sigma*Vt^+*s - H12*A * 1/(om-ekp) * (H12*A)^+ )^{-1} * s*Vt
 
 The crucial observation is that the dimension of H11 is norbs x norbs, as opposed to original G^{-1} 
 size of nbands x nbands. We thus need to invert only norbs x norbs matrix to get the local Green's function.
 Note that SVD problem is frequency independent problem, and diagonalization of H22 is also done 
 outside frequency loop, so that for each frequency om, we just need to embbed Sigma with 
 new embeddor U_new==s*Vt and add improper part of the self-energy, which is H12*A * 1/(om-ekp) * (H12*A)^+. 
 This comes from the presence of all other bands in the system, which are not correlated. 
 If we defined hi12== H12*A, we have
 
  G = U_new^+ * (om - H11 - U_new*Sigma*U_new^+ - hi12*1/(om-ekp)*hi12^+ )^{-1} * U_new
 
 and at each frequency we need to effectively only invert the small matrix norbs x norbs, 
 and compute matrix product of hi12*1/(om-ekp)*hi12^+, which is amounts to one matrix-matrix 
 multiplication of norbs x norbs.
 We can also compute total DOS, which is Tr(G), but in this case we need to add back 
 the second part of the matrix. We have:
  
  Tr(G) = ( om - H11 - U_new*Sigma*U_new^+ - hi12*1/(om-ekp)*hi12^+ )^{-1},  ......                                                          )
          (..............................................................., (om - H22 - H21*(om - H11 - U_new*Sigma*U_new^+)^{-1} *H12)^{-1} )
 
 which translates into

  Tr(G) = Tr( ( om - H11 - U_new*Sigma*U_new^+ - hi12*1/(om-ekp)*hi12^+ )^{-1} ) + Tr( (om - H22 - H21*(om - H11 - U_new*Sigma*U_new^+)^{-1} *H12 )^{-1} )

 and requires a bit more work.
    """
    parser = argparse.ArgumentParser(usage=usage)
    parser.add_argument('-m', dest='mode', default='g', type=str, help='mode can be g or e')
    parser.add_argument('-i', dest='interp', default=False, action='store_true', help='interpolate sigma on equidistant mesh from case,indmfl')
    parser.add_argument('-s', dest='store_Ak', default=False, action='store_true', help='store A_kom.npy for plotting with plt_sakplot.py')
    parser.add_argument('-l', dest='m_ext', type=str, default='', help='reading sig.inpxdn/sig.inpx, when set to dn/left-blank')
    parser.add_argument('-f', dest='full_inversion', default=False, action='store_true', help='perform full inversion, which is slower but safer')
    parser.add_argument('-c', dest='cutoff', type=float, default=1.0, help='cutoff energy is taken from case.indmfl and increased by cutoff (1eV by default)') 
    args = parser.parse_args()
    main(args)
            
