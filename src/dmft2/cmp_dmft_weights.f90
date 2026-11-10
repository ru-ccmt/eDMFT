! @Copyright 2007 Kristjan Haule
! 

SUBROUTINE cmp_dmft_weights2(Aweight, gloc, logG, dens, Gdloc, Gdloc0, &
     DMFTU, dNds, Edimp0, AEweight, nedim, STrans, sigma, s_oo, wkp, &
     omega, vnorm1, DM_EF, Temperature, wgamma, gamma, csize, nbands, &
     DM_nemin, DM_nemax, maxsize, ncix, nomega, npomega, nume, &
     matsubara, lgTrans, lg_deg, ntcix, nipc, ikp, maxdim2, iorbital, &
     iSx, norbitals, projectorGradientDirect, cix_orb, cixdim, nindo, &
     Sigind)
  USE defs,  ONLY: pi, IMAG
  USE xa,    ONLY: E
  USE muzero,ONLY: nomq, iomq, jomq, womq, abom, n0_om, matsum_weight
  USE param, ONLY: CHUNK
  use defs,  only: Ry2eV
  USE dmfts, only: maxdim, lmaxp, natom
  IMPLICIT NONE
  COMPLEX*16, intent(out)  :: Aweight(nbands,nbands), AEweight(nedim,nedim)
  REAL*8, intent(inout)    :: dens
  REAL*8, intent(out)      :: logG, dNds(nipc,ntcix)
  COMPLEX*16, intent(inout) :: Gdloc(ntcix,nomega,nipc), gloc(nomega), Gdloc0(maxdim,maxdim,ncix,nomega)
  COMPLEX*16, intent(inout) :: Edimp0(maxdim,maxdim,ncix)
  COMPLEX*16, intent(in) :: STrans(maxsize,ncix,nbands,nbands)
  COMPLEX*16, intent(in) :: sigma(maxsize,ncix,nomega), s_oo(maxsize,ncix), wkp, lgTrans(nbands,nbands,ntcix,nipc)
  REAL*8, intent(in)     :: omega(nomega), vnorm1, DM_EF, Temperature, wgamma, gamma
  INTEGER, intent(in)    :: csize(ncix), nbands, DM_nemin, DM_nemax, maxsize, ncix, nomega, npomega, nume, lg_deg(ntcix)
  INTEGER, intent(in)    :: iorbital(natom,lmaxp+1), iSx(maxdim2, norbitals)
  INTEGER, intent(in)    :: ntcix, nipc, nedim, ikp, maxdim2, norbitals
  INTEGER, intent(in)    :: cix_orb(norbitals), cixdim(ncix), &
       nindo(norbitals), Sigind(maxdim,maxdim,ncix)
  LOGICAL, intent(in)    :: matsubara
  COMPLEX*16, intent(in) :: DMFTU(nbands,maxdim2,norbitals,nipc)
  COMPLEX*16, intent(out) :: &
       projectorGradientDirect(nbands,maxdim2,norbitals)
  ! functions
  Interface
     real*8 FUNCTION ferm(x)
       real*8, intent(in) :: x
     end FUNCTION ferm
     SUBROUTINE GetFreeEnergy(tlogG, tlogG0, tlogGD, tdens, Temperature, DM_EF, iom, omega, nomega, nbands, E, zek, DM_nemin, nume, gamma)
       REAL*8, intent(out) :: tlogG, tlogG0, tlogGD, tdens
       INTEGER, intent(in) :: iom, nomega, nbands, DM_nemin, nume
       REAL*8, intent(in)  :: Temperature, omega(nomega), E(nume), gamma, DM_EF
       COMPLEX*16, intent(in):: zek(nbands)
     end SUBROUTINE GetFreeEnergy
     real*8 FUNCTION FreeE0(Energy, Temp)
       real*8, intent(in) :: Energy, Temp
     end FUNCTION FreeE0
     SUBROUTINE AddSigma_optimized2(gij, sigma, STrans, csize, sign, nbands, ncix, maxsize)
       COMPLEX*16, intent(inout):: gij(nbands,nbands)
       COMPLEX*16, intent(in)   :: sigma(maxsize,ncix)
       COMPLEX*16, intent(in)   :: STrans(maxsize,ncix,nbands,nbands)
       INTEGER, intent(in)      :: csize(ncix), sign, nbands, ncix, maxsize
     end SUBROUTINE AddSigma_optimized2
     SUBROUTINE eigsys(ham, zek, evl, evr, ndim)
       COMPLEX*16, INTENT(inout):: ham(ndim,ndim)  ! Hamiltonian / overwritten
       COMPLEX*16, INTENT(out) :: zek(ndim)        ! Vector of eigenvalues 
       COMPLEX*16, INTENT(out) :: evl(ndim,ndim)   ! left eigenvectors
       COMPLEX*16, INTENT(out) :: evr(ndim,ndim)   ! right eigenvectors
       INTEGER, INTENT(in)     :: ndim             ! Dimension of hamiltonian
     end SUBROUTINE eigsys
     REAL*8 Function ProductTrace(Sigmaij,Gij,nbands)
       COMPLEX*16, intent(in) :: Sigmaij(nbands,nbands), Gij(nbands,nbands)
     end Function ProductTrace
     COMPLEX*16 Function ProductTrace2(A, B, ndim)
       ! returns  sum_{ij} A(i,j)*B(j,i)                                                                                                                                                                                                      
       COMPLEX*16, intent(in) :: A(ndim,ndim), B(ndim,ndim)
       INTEGER, intent(in)    :: ndim
     end Function ProductTrace2
  end interface
  ! locals
  COMPLEX*16 :: sigma_winf(maxsize,ncix)
  COMPLEX*16, allocatable :: Wsum(:,:), WEsum(:,:), w0inf(:), Identity(:,:)
  COMPLEX*16, allocatable :: WSigPair(:,:), SigDirect(:,:), SigScratch(:,:), &
       SigmaTail1(:,:)
  COMPLEX*16, allocatable :: CPair(:,:,:), CStatic(:,:,:), &
       Uactive(:,:,:), SigmaLocal(:,:,:), SigmaInfLocal(:,:,:), &
       SigmaTail1Local(:,:,:), CdiagUdagG(:,:), &
       CdiagUdagGminus(:,:), CdiagProduct(:,:), CdiagProductMinus(:,:)
  COMPLEX*16, allocatable :: Eij(:,:), zek(:), zar(:,:), zal(:,:), tmp(:,:), Sigmaij(:,:), Sigma_inf(:,:), Gij(:,:), tmp2(:,:)
  REAL*8,     allocatable :: zek_inf(:)
  COMPLEX*16, allocatable :: zw1(:), w1(:,:), zar_inf(:,:), zal_inf(:,:), w2(:,:), w3(:,:), Dni(:)
  COMPLEX*16, ALLOCATABLE :: Winf(:,:), Wferm(:,:), Bii(:)
  REAL*8     :: logGD, logG0
  COMPLEX*16 :: cc, gtc, cek, omn, w0, wo, bg, sginf, sginf_!, cc0, cc1, cc2
  REAL*8  :: fer, tlogGD, tlogG, tlogG0, tdens, ek_inf, beta, csm, wkpt, &
       C0, sginf2, second_order_correction, fermi, matsuWeight, tailMoment2!, diff
  INTEGER :: i, iom, num, num0, it, nomq_max, icix, j, ip, iq, iorb, &
       ind, nloc, start_nipc
  LOGICAL :: Qforce, stressMomentDiagnostic, computeDirectC
  CHARACTER(len=64) :: diagnosticEnvironment
  INTEGER :: diagnosticStatus
  REAL*8 :: directNorm, directResidual, directScale
  complex*16 :: czero, cone
  cone  = cmplx(1.d0,0.d0, 8)
  czero = cmplx(0.d0,0.d0, 8)
  
  start_nipc=2
  
  Qforce=.false.
  if (nedim.eq.nbands) Qforce=.true.
  projectorGradientDirect=czero
  AEweight=czero

  diagnosticEnvironment=''
  CALL GET_ENVIRONMENT_VARIABLE('EDMFT_STRESS_PROJECTOR_DIAG', &
       diagnosticEnvironment,STATUS=diagnosticStatus)
  stressMomentDiagnostic=Qforce .and. diagnosticStatus.eq.0 .and. &
       len_trim(diagnosticEnvironment).gt.0 .and. &
       trim(diagnosticEnvironment).ne.'0'

  ! Accumulate C=Sum Sigma_local U^dagger G for every derivative run.
  ! U*C=rhoE-epsilon*rho does not determine ker(U) components of C;
  ! reconstructing C with a pseudoinverse can therefore lose a real force.
  ! Density-only runs still skip these rectangular matrix products.
  computeDirectC=Qforce .and. matsubara

  if (Qforce) then
     allocate( WEsum(nedim,nedim) ) ! Unfortunately we need this for $OMP command to work.
  else
     allocate( WEsum(1,1) )
  endif
  allocate(Wsum(nbands,nbands), w0inf(nbands))
  if (stressMomentDiagnostic) then
     allocate(WSigPair(nbands,nbands),SigDirect(nbands,nbands), &
          SigScratch(nbands,nbands))
     WSigPair=0.d0
  endif

  beta = 1/Temperature
  wkpt = dble(wkp)*vnorm1  
  if (matsubara) then
     C0=1.d0
     allocate( Sigmaij(nbands,nbands) )
     ! Determine Sigma_1 in
     !
     !      Sigma(iw)-Sigma_inf = Sigma_1/(iw)+Sigma_2/(iw)^2+...
     !
     ! from an adaptive two-moment fit.  Smooth tails use the fitted Sigma_1;
     ! noisy tails use a short plateau average; clearly bad/increasing tails
     ! retain the historical endpoint estimate -w_N Im Sigma(iw_N).  The
     ! classification is printed only once (at the first k point).
     CALL DetermineSelfEnergyTailMoment(sigma_winf,sigma,s_oo,omega,csize, &
          nomega,maxsize,ncix,ikp.eq.1)
     allocate(SigmaTail1(nbands,nbands))
     SigmaTail1=0.d0
     CALL AddSigma_optimized2(SigmaTail1,sigma_winf,STrans,csize,1, &
          nbands,ncix,maxsize)
     ! Roundoff, a shared compact Sigind label, or a noisy off-diagonal fit can
     ! leave a tiny anti-Hermitian component.  Only the Hermitian moment enters
     ! the paired Matsubara sum.
     SigmaTail1=0.5d0*(SigmaTail1+conjg(transpose(SigmaTail1)))
     ! The older log(G) completion below calls this matrix Sigmaij.  Keep the
     ! same fitted moment in both the density/stress and free-energy tails;
     ! otherwise two formally equivalent Matsubara completions would use
     ! different estimates of Sigma_1.
     Sigmaij=SigmaTail1

     if (computeDirectC) then
        allocate(CPair(maxdim,nbands,ncix), &
             CStatic(maxdim,nbands,ncix), &
             Uactive(nbands,maxdim,ncix), &
             SigmaLocal(maxdim,maxdim,ncix), &
             SigmaInfLocal(maxdim,maxdim,ncix), &
             SigmaTail1Local(maxdim,maxdim,ncix), &
             CdiagUdagG(maxdim,nbands), &
             CdiagUdagGminus(maxdim,nbands), &
             CdiagProduct(maxdim,nbands), &
             CdiagProductMinus(maxdim,nbands))
        CPair=czero
        Uactive=czero

        ! Assemble the same aggregate local-to-band projector which is
        ! implicit in CompressSigmaTransformation1.  Several orbital records
        ! may belong to one impurity block, hence the addition rather than an
        ! assignment.  This also makes Uactive*Sigma*Uactive^+ reproduce the
        ! embedded self-energy used to construct G.
        do iorb=1,norbitals
           icix=cix_orb(iorb)
           if (icix.le.0) cycle
           do ind=1,nindo(iorb)
              ip=iSx(ind,iorb)
              Uactive(:,ip,icix)=Uactive(:,ip,icix)+ &
                   DMFTU(:,ind,iorb,1)
           enddo
        enddo
        CALL BuildLocalSigmaMatrix(SigmaInfLocal,s_oo,Sigind,cixdim, &
             maxdim,maxsize,ncix)
        CALL BuildLocalSigmaMatrix(SigmaTail1Local,sigma_winf,Sigind, &
             cixdim,maxdim,maxsize,ncix)
        do icix=1,ncix
           nloc=cixdim(icix)
           SigmaTail1Local(1:nloc,1:nloc,icix)=0.5d0*( &
                SigmaTail1Local(1:nloc,1:nloc,icix)+ &
                conjg(transpose( &
                SigmaTail1Local(1:nloc,1:nloc,icix))))
        enddo
     endif

     ! Exact positive-frequency sum from the beginning of the sparse tail:
     !   sum_{m=j0}^inf 1/w_m^2,  w_m=(2m-1)pi/beta.
     ! The closed full sum beta^2/8 avoids any special-function dependency;
     ! subtracting fewer than a few hundred positive terms is stable in double
     ! precision.  This scalar is used to restore the Sigma_1/(iw)^2 term
     ! after it has been removed from the cubic-spline quadrature.
     tailMoment2=beta*beta/8.d0
     do i=1,jomq(n0_om)-1
        tailMoment2=tailMoment2- &
             (beta/(dble(2*i-1)*pi))**2
     enddo

     allocate( zek(nbands), Eij(nbands,nbands) )
     ALLOCATE( zek_inf(nbands) )
     ALLOCATE( zar_inf(nbands,nbands), zal_inf(nbands,nbands) )
     ALLOCATE( Sigma_inf(nbands,nbands) )
     ALLOCATE( Bii(nbands) )
     
     !----------- Computes the infinite frequency limit
     Eij=0
     DO i=1,nbands
        Eij(i,i) = E(i+DM_nemin-1)  !---   preparing hamiltonian: epsk+sigma
     ENDDO
     ! Adding self-energy
     Sigma_inf=0
     CALL AddSigma_optimized2(Sigma_inf, s_oo, STrans, csize, 1, nbands, ncix, maxsize)
     Eij = Eij + Sigma_inf
     ! Here we might want to use Hermitian eigsys!!!!????
     CALL eigsys(Eij, zek, zal_inf, zar_inf, nbands) ! eigenvalues and eigenvectors
     zek_inf(:) = dble(zek(:))-DM_EF
     
     allocate(tmp(nbands,nbands))
     ! Bii = zal_inf * (Sigmaij*iw) * zar_inf                                                                                                                                                                                                      
     CALL zgemm('N','N',nbands,nbands,nbands,cone,zal_inf,nbands,Sigmaij,nbands,czero,tmp,nbands)
     CALL ProductTrace3(Bii, tmp, zar_inf, nbands)
     
     deallocate( Sigmaij )
     
     deallocate( zek, Eij )
     
     w0inf=0.d0
     Wsum=0.d0
     WEsum=0.d0
     Aweight(:,:)=0.d0
     AEweight(:,:)=0.d0
     logG0=0.d0
     logGD=0.d0

     allocate( Identity(nbands,nbands) )
     Identity(:,:) = 0.d0                         ! Initialize the array.
     forall(j=1:nbands) Identity(j,j) = 1.d0     ! Set the diagonal.
     
     sginf=0 ! This will be used for the second-order correction to Tr(Sigma*G) ~ sum_{iw} Bii(i)/(iw-C0) * 1/(iw-zek_inf(i)) = Bii(i)*(f(C0)-f(zek_inf(i)))/(C0-zek_inf(i))   

     ! !$OMP  PARALLEL DO SHARED(gloc,Sigma_inf,Gdloc,ntcix,wkpt,Gdloc0)&
     ! !$OMP& PRIVATE(Eij,Sigmaij,zek,tlogG0,tlogGD,tdens,gtc,num,num0,ek_inf,cek,omn,Gij,it,ip,bg,cc)&
     ! !$OMP& SCHEDULE(STATIC,CHUNK) REDUCTION(+:w0inf,Wsum,WEsum,dens,logG0,logGD,sginf)
     do iom=1,n0_om-1
        ! Necessary allocation due to openMP
        allocate( Eij(nbands,nbands), zek(nbands) )  
        allocate( Sigmaij(nbands,nbands), Gij(nbands,nbands) )
        ! Total hamiltonian (H+Sigma) in Kohn-Sham basis
        Eij=0
        DO i=1,nbands
           Eij(i,i) = E(i+DM_nemin-1)  !---   preparing hamiltonian: epsk+sigma
        ENDDO
        ! Adding self-energy
        Sigmaij=0
        CALL AddSigma_optimized2(Sigmaij, sigma(:,:,iom), STrans, csize, 1, nbands, ncix, maxsize)
        Eij = Eij + Sigmaij

        omn = omega(iom)*IMAG
        ! Computing the lattice Green's function                                                                                                                                                                                              
        Gij(:,:)=-Eij(:,:) ! off diagonal components are just G^{-1} = -(E+Sigma)
        do i=1,nbands
           Gij(i,i)=omn+DM_EF+Gij(i,i) ! diagonal components corrected
        enddo
        CALL zinv(Gij, nbands)
        ! Density from the Green's function                                                                                                                                                                                                   
        Wsum = Wsum + Gij
        
        if (Qforce) then
           WEsum = WEsum + Gij(:,:)*(omn+DM_EF) - Identity(:,:) ! This is summing A^R  eps_{k,w}/(iom+mu-eps_{k,w}) A^L = (iom+mu)*G - I          
        endif
        if (stressMomentDiagnostic) then
           ! Accumulate both positive- and negative-frequency orientations.
           ! The negative partner of Sigma*G is Sigma^+*G^+, not
           ! (Sigma*G)^+=G^+*Sigma^+.  Retaining this order is essential for
           ! comparing the generally non-Hermitian matrix rhoE-epsilon*rho.
           SigDirect=matmul(Sigmaij,Gij)
           SigScratch=matmul(conjg(transpose(Sigmaij)), &
                conjg(transpose(Gij)))
           WSigPair=WSigPair+SigDirect+SigScratch
        endif
        if (computeDirectC) then
           CALL BuildLocalSigmaMatrix(SigmaLocal,sigma(:,:,iom), &
                Sigind,cixdim,maxdim,maxsize,ncix)
           CALL AccumulateDirectProjectorPair(CPair,Gij,Uactive, &
                SigmaLocal,cixdim,nbands,maxdim,ncix,CdiagUdagG, &
                CdiagUdagGminus,CdiagProduct,CdiagProductMinus)
        endif
        
        CALL eigvals(Eij, zek, nbands) ! eigenvalues only

        tlogGD=0.0
        tlogG0=0.0
        tdens=0.0
        gtc=0
        bg=0
        !$OMP  PARALLEL DO SHARED(zek_inf,zek)&
        !$OMP& PRIVATE(num,num0,ek_inf,cek,cc)&
        !$OMP& SCHEDULE(STATIC,CHUNK) REDUCTION(+:w0inf,bg,tlogG0,tlogGD,tdens,gtc)   
        DO num=DM_nemin,DM_nemax
           ! These frequencies are not yet in logarithmic mesh (Where all consecutive matsubara frequencies are present)
           num0 = num-DM_nemin+1
           ek_inf = zek_inf(num0)
           cek = zek(num0)-DM_EF           ! eps_k-mu
           cc = 1.d0/(omn-ek_inf)
           w0inf(num0) = w0inf(num0) + cc  ! w0inf = sum( 1/(i*om+mu-ek_inf) , om)
           bg = bg + Bii(num0)*cc
           tlogGD = tlogGD - dble(log(-omn+cek))
           tlogG0 = tlogG0 - dble(log(-omn+ek_inf))
           tdens  = tdens  + dble( 1/(omn-cek)-1/(omn-ek_inf) )
           gtc = gtc + 1/(omn-cek)
        ENDDO
        !$OMP END PARALLEL DO
        tlogGD = 2/beta*tlogGD
        tlogG0 = 2/beta*tlogG0
        tdens  = 2/beta*tdens
        
        sginf =  sginf + bg/(omn-C0)

        gloc(iom) = gloc(iom) + gtc * (0.5*wkp)

        dens  = dens  + tdens * wkpt
        logG0 = logG0 + tlogG0* wkpt
        logGD = logGD + tlogGD* wkpt

        DO it=1,ntcix  ! sum_{ij} Gij(i,j) * STrans(j,i)
           if (lg_deg(it).eq.0) cycle
           do ip=start_nipc,nipc
              Gdloc(it,iom,ip) = Gdloc(it,iom,ip) + ProductTrace2(Gij, lgTrans(:,:,it,ip), nbands)*(wkp*0.5)  ! 0.5 because wkp=2*wk/sumw
           enddo
        ENDDO
        Call ProjectToLocal(Gdloc0(:,:,:,iom), Gij, DMFTU(:,:,:,1), 0.5*dble(wkp), iorbital, iSx, nbands, maxdim2, norbitals)
        deallocate( Sigmaij, Gij )
        deallocate( Eij, zek )
     enddo   ! iom loop
     ! !$OMP END PARALLEL DO

     deallocate( Identity )
     
     nomq_max=1
     do iom=n0_om,npomega
        nomq_max = max(nomq_max,nomq(iom))
     enddo

     allocate( Eij(nbands,nbands), zek(nbands), zar(nbands,nbands), zal(nbands,nbands) )
     allocate( tmp2(nbands,nbands) )
     allocate( Sigmaij(nbands,nbands), Gij(nbands,nbands) )
     allocate( Dni(nbands) )
     
     do iom=n0_om,npomega
        
        if ( abs((2*jomq(iom)-1)*pi/beta - omega(iom)).gt.1e-6 ) then
           print *, 'ERROR: ', iom, (2*jomq(iom)-1)*pi/beta, omega(iom)
        endif
        
        ! Total hamiltonian (H+Sigma) in Kohn-Sham basis
        Eij=0
        DO i=1,nbands
           Eij(i,i) = E(i+DM_nemin-1)  !---   preparing hamiltonian: epsk+sigma
        ENDDO
        ! Adding self-energy
        Sigmaij=0
        CALL AddSigma_optimized2(Sigmaij, sigma(:,:,iom), STrans, csize, 1, nbands, ncix, maxsize)
        Eij = Eij + Sigmaij
        CALL eigsys(Eij, zek, zal, zar, nbands) ! eigenvalues and eigenvectors                                                                                                                                                                

        tlogGD=0.d0
        tlogG0=0.d0
        tdens=0.d0
        gtc=0.d0
        sginf_=0.d0
        !$OMP PARALLEL DO SHARED(w0inf,tmp,iom,Dni,C0)&
        !$OMP PRIVATE(num0,ek_inf,cek,omn,w0,i,wo,cc,bg)&
        !$OMP& SCHEDULE(STATIC,CHUNK) REDUCTION(+:tlogGD,tlogG0,tdens,sginf_)
        DO num=DM_nemin,DM_nemax
           num0 = num-DM_nemin+1
           ek_inf = zek_inf(num0)           
           cek = zEk(num0)-DM_EF
           w0=0.d0
           bg=0.d0
           wo=0.d0
           ! Cubic-spline quadrature of the final integrand.  matsum_weight
           ! is the sum of the cardinal cubic spline over every omitted
           ! Matsubara point.  In contrast to the historical scheme, no pole
           ! or eigenvalue is interpolated before the nonlinear resolvent is
           ! formed.
           omn = IMAG*omega(iom)
           matsuWeight=matsum_weight(iom)
           w0 = matsuWeight/(omn-cek)
           cc = matsuWeight/(omn-ek_inf)
           wo = cc
           bg = cc/(omn-C0)
           tlogGD = tlogGD-matsuWeight*dble(log(-omn+cek))
           tlogG0 = tlogG0-matsuWeight*dble(log(-omn+ek_inf))
           tdens  = tdens+matsuWeight*dble(1/(omn-cek)-1/(omn-ek_inf))
           w0inf(num0) = w0inf(num0) + wo
           tmp(:,num0) =  zar(:,num0)*w0
           sginf_ = sginf_ + bg*Bii(num0)

           Dni(num0) = 1./(IMAG*omega(iom)-cek)
        ENDDO
        !$OMP END PARALLEL DO

        ! omn is PRIVATE in the OpenMP band loop above.  Set the shared scalar
        ! explicitly before using it in the matrix-valued tail subtraction;
        ! otherwise it retains a stale frequency from the preceding dense
        ! loop and turns the intended 1/omega^2 correction into a nearly
        ! constant, strongly mesh-dependent shift.
        omn=IMAG*omega(iom)

        gtc = sum(Dni)
        sginf = sginf + sginf_
        ! Wsum = Ar * 1/(i*om+mu-ek) * Al
        CALL zgemm('N','N',nbands,nbands,nbands,cone,tmp,nbands,zal,nbands,czero,Gij,nbands)
        Wsum = Wsum + Gij

        if (Qforce) then
           ! WEsum = Ar * 1/(i*om+mu-ek) * ek * Al
           do num0=1,nbands
              tmp(:,num0) = tmp(:,num0)*zek(num0)
           enddo
           CALL zgemm('N','N',nbands,nbands,nbands,cone,tmp,nbands,zal,nbands,cone,WEsum,nbands)
           ! Remove the leading dynamic high-frequency term from the spline
           ! sum.  Since (iw)^2=-w^2 this subtraction numerically adds a
           ! positive Sigma_1/w^2 matrix here.  The exact discrete sum of the
           ! removed term, from this same n0_om boundary to infinity, is
           ! restored below.  Subtracting and restoring over identical limits
           ! prevents double counting and also makes the splined residual
           ! decay two powers faster.
           WEsum=WEsum-matsum_weight(iom)*SigmaTail1/(omn*omn)
        endif
        if (stressMomentDiagnostic) then
           ! Gij here is the quadrature-interpolated positive-frequency
           ! Green function for this logarithmic-mesh node.  Sigmaij is held
           ! at the same node, exactly as in the WEsum construction above.
           SigDirect=matmul(Sigmaij,Gij)- &
                matsum_weight(iom)*SigmaTail1/(omn*omn)
           SigScratch=matmul(conjg(transpose(Sigmaij)), &
                conjg(transpose(Gij)))- &
                matsum_weight(iom)*SigmaTail1/(omn*omn)
           WSigPair=WSigPair+SigDirect+SigScratch
        endif
        if (computeDirectC) then
           CALL BuildLocalSigmaMatrix(SigmaLocal,sigma(:,:,iom), &
                Sigind,cixdim,maxdim,maxsize,ncix)
           CALL AccumulateDirectProjectorPair(CPair,Gij,Uactive, &
                SigmaLocal,cixdim,nbands,maxdim,ncix,CdiagUdagG, &
                CdiagUdagGminus,CdiagProduct,CdiagProductMinus)

           ! Gij already contains matsum_weight at a sparse node.  Remove
           ! the paired Sigma_1 U^+/(iw)^2 asymptote from that same weighted
           ! integrand; its exact discrete tail is restored below.
           do icix=1,ncix
              nloc=cixdim(icix)
              CALL zgemm('N','C',nloc,nbands,nloc,cone, &
                   SigmaTail1Local(1,1,icix),maxdim, &
                   Uactive(1,1,icix),nbands,czero,CdiagProduct,maxdim)
              CPair(1:nloc,:,icix)=CPair(1:nloc,:,icix)- &
                   2.d0*matsum_weight(iom)* &
                   CdiagProduct(1:nloc,:)/(omn*omn)
           enddo
        endif
        
        ! Evaluating G_local
        do j=1,nbands
           tmp2(:,j) = zar(:,j)*Dni(j)
        enddo
        CALL zgemm('N','N',nbands,nbands,nbands,cone,tmp2,nbands,zal,nbands,czero,Gij,nbands)
        DO it=1,ntcix
           if (lg_deg(it).eq.0) cycle
           !CALL zgemm('N','T',nbands,nbands,nbands,cone,zal,nbands,lgTrans(:,:,it),nbands,czero,tmp,nbands)
           !
           ! Here we want to compute   Gdloc = \sum_{ij} (A^R * g * A^L)_{ij} * lgTrans_{ij} = Tr( A^R * g * A^L * lgTrans^T )
           !   We first compute   tmp = lgTrans * {A^L}^T
           !          and next    tmp2_{ij} = A^R_{ij}*g_j
           !          finally     G_loc = \sum_{ij} tmp2_{ji} tmp_{ji} = A^R_{ji} g_i A^L_{il} * lgTrans^T_{lj} 
           do ip=start_nipc,nipc
              !CALL zgemm('N','T',nbands,nbands,nbands,cone,lgTrans(:,:,it,ip),nbands,zal,nbands,czero,tmp,nbands)
              !Gdloc(it,iom,ip) = Gdloc(it,iom,ip) +  ProductTrace2(tmp2,tmp,nbands)* 0.5*wkp
              !
              Gdloc(it,iom,ip) = Gdloc(it,iom,ip) + ProductTrace2(Gij, lgTrans(:,:,it,ip), nbands)*(wkp*0.5)  ! 0.5 because wkp=2*wk/sumw
           enddo
        ENDDO
        Call ProjectToLocal(Gdloc0(:,:,:,iom), Gij, DMFTU(:,:,:,1), 0.5*dble(wkp), iorbital, iSx, nbands, maxdim2, norbitals)
        tlogGD = 2/beta*tlogGD
        tlogG0 = 2/beta*tlogG0
        tdens  = 2/beta*tdens
        
        gloc(iom) = gloc(iom) + gtc * (0.5*wkp)
        dens  = dens  + tdens * wkpt
        logG0 = logG0 + tlogG0* wkpt
        logGD = logGD + tlogGD* wkpt
        
     enddo   ! iom loop
     
     deallocate( Sigmaij )
     deallocate( Eij, zek, zar, zal )
     deallocate( Gij )
     deallocate( Dni )
     
     ALLOCATE( Winf(nbands,nbands), Wferm(nbands,nbands) )
     do num0=1,nbands
        tmp(:,num0) =   zar_inf(:,num0)*w0inf(num0)
     enddo
     Winf = matmul(tmp,zal_inf)

     tlogG=0
     !$OMP PARALLEL DO SHARED(tmp) PRIVATE(fermi)&
     !$OMP& SCHEDULE(STATIC,CHUNK) REDUCTION(+:dens,tlogG)
     do num0=1,nbands
        fermi = ferm(zek_inf(num0)*beta)
        tmp  (:,num0) =   zar_inf(:,num0)*fermi
        dens = dens + fermi*wkpt
        tlogG = tlogG + (FreeE0(zek_inf(num0),Temperature)-FreeE0(E(num0+DM_nemin-1)-DM_EF,Temperature))*wkpt
     enddo
     !$OMP END PARALLEL DO
     Wferm = matmul(tmp,zal_inf)

     if (stressMomentDiagnostic) then
        ! Apply exactly the same static-tail subtraction and Fermi completion
        ! as Wsum/WEsum, but accumulate Sigma*G independently and preserve
        ! the order of the negative-frequency product.
        SigDirect=matmul(Sigma_inf,Winf)
        SigScratch=matmul(conjg(transpose(Sigma_inf)), &
             conjg(transpose(Winf)))
        WSigPair=WSigPair-SigDirect-SigScratch
        WSigPair=WSigPair-2.d0*tailMoment2*SigmaTail1
        WSigPair=WSigPair*(2.d0/beta)
        SigDirect=matmul(Sigma_inf,Wferm)
        SigScratch=matmul(conjg(transpose(Sigma_inf)), &
             conjg(transpose(Wferm)))
        SigDirect=(0.5d0*wkp)*(WSigPair+SigDirect+SigScratch)
     endif

     if (computeDirectC) then
        ! Complete precisely the same Matsubara construction as the
        ! band-space Sigma*G diagnostic, but before multiplying by U.  The
        ! result is the rectangular C conjugate to an explicit dU.
        CStatic=czero
        CALL AccumulateDirectProjectorPair(CStatic,Winf,Uactive, &
             SigmaInfLocal,cixdim,nbands,maxdim,ncix,CdiagUdagG, &
             CdiagUdagGminus,CdiagProduct,CdiagProductMinus)
        CPair=CPair-CStatic
        do icix=1,ncix
           nloc=cixdim(icix)
           CALL zgemm('N','C',nloc,nbands,nloc,cone, &
                SigmaTail1Local(1,1,icix),maxdim, &
                Uactive(1,1,icix),nbands,czero,CdiagProduct,maxdim)
           CPair(1:nloc,:,icix)=CPair(1:nloc,:,icix)- &
                2.d0*tailMoment2*CdiagProduct(1:nloc,:)
        enddo
        CPair=CPair*(2.d0/beta)
        CStatic=czero
        CALL AccumulateDirectProjectorPair(CStatic,Wferm,Uactive, &
             SigmaInfLocal,cixdim,nbands,maxdim,ncix,CdiagUdagG, &
             CdiagUdagGminus,CdiagProduct,CdiagProductMinus)
        CPair=(0.5d0*wkp)*(CPair+CStatic)

        ! Restore the external storage convention used by the existing
        ! Moore--Penrose C: band index first, orbital row second.
        do iorb=1,norbitals
           icix=cix_orb(iorb)
           if (icix.le.0) cycle
           do ind=1,nindo(iorb)
              ip=iSx(ind,iorb)
              projectorGradientDirect(:,ind,iorb)=CPair(ip,:,icix)
           enddo
        enddo
     endif
     
     tmp(:,:) = (Wsum(:,:) - Winf(:,:))*(2/beta) + Wferm(:,:)
     Aweight(:,:) = (0.5*wkp)*(tmp(:,:) + transpose(conjg(tmp(:,:))))

     if (Qforce) then
        !$OMP PARALLEL DO SHARED(tmp,tmp2) SCHEDULE(STATIC)
        do num0=1,nbands
           tmp (:,num0) = zar_inf(:,num0)*(     w0inf(num0)          * (zek_inf(num0)+DM_EF) )
           tmp2(:,num0) = zar_inf(:,num0)*( ferm(zek_inf(num0)*beta) * (zek_inf(num0)+DM_EF) )
        enddo
        !$OMP END PARALLEL DO
        !    Winf = zar_inf * w0inf * eps_inf * zal_inf
        Winf = matmul(tmp, zal_inf)
        !   Wferm = zar_inf * ferm * eps_inf * zal_inf
        Wferm = matmul(tmp2, zal_inf)
        tmp = (WEsum(:,:) - Winf(:,:))*(2/beta) + Wferm(:,:)
        tmp = tmp-(2.d0/beta)*tailMoment2*SigmaTail1
        AEweight(:,:) = (0.5*wkp)*(tmp(:,:) + transpose(conjg(tmp(:,:))))
     endif

     if (stressMomentDiagnostic) then
        ! B=rhoE-epsilon*rho is the source used by the complementary-band
        ! response.  Compare it with a separately accumulated Sigma_emb*G,
        ! including the identical Matsubara tail construction.
        SigScratch=AEweight
        do i=1,nbands
           SigScratch(i,:)=SigScratch(i,:)- &
                E(i+DM_nemin-1)*Aweight(i,:)
        enddo
        directNorm=sqrt(sum(abs(SigScratch)**2))
        directResidual=sqrt(sum(abs(SigScratch-SigDirect)**2))
        directScale=max(directNorm,sqrt(sum(abs(SigDirect)**2)),tiny(1.d0))
        WRITE(6,'(A,I6,5(1X,ES18.9))') ':DMFT_SIGMAG_DIRECT:',ikp, &
             directNorm,directResidual/directScale, &
             dble(sum([(SigScratch(i,i),i=1,nbands)])), &
             dble(sum([(SigDirect(i,i),i=1,nbands)])), &
             aimag(sum([(SigDirect(i,i),i=1,nbands)]))
     endif

     csm=0
     do i=1,nbands
        csm = csm + real(Aweight(i,i))
     enddo

     sginf = dble(sginf)*2/beta
     sginf2=0
     do num0=1,nbands
        ek_inf = zek_inf(num0)
        sginf2 = sginf2 + dble(Bii(num0)*(ferm(C0*beta)-ferm(ek_inf*beta))/(C0-ek_inf))
        !print *, num0, dble(ek_inf), (ferm(C0*beta)-ferm(ek_inf*beta)), (C0-ek_inf), dble(Bii(num0)*(ferm(C0*beta)-ferm(ek_inf*beta))/(C0-ek_inf))                                                                                           
     enddo
     second_order_correction = sginf2 - dble(sginf)
     
     !print *, 'sginf=', dble(sginf), 'sginf2=', dble(sginf2)                                                                                                                                                                                 
     
     logG = logGD - logG0 + tlogG + second_order_correction*wkpt
     
     !WRITE(*,'(A,F12.6,1x,A,F12.6,1x,A,F12.6,1x,A,F12.6)') 'logGD-logG0=', logGD-logG0, 'tlogG=', tlogG, 'logG=', logG
     !WRITE(*, '(A,F12.6,1x,A,F12.6,1x)') 'csm=', csm 
     
     ! Computing impurity levels in the old way
     !DO it=1,ntcix  ! sum_{ij} Gij(i,j) * STrans(j,i)                                                                                                                                                                                         
     !   if (lg_deg(it).eq.0) cycle
     !   csm=0
     !   DO i=1,nbands
     !      csm = csm + E(i+DM_nemin-1)*lgTrans(i,i,it,1)
     !   ENDDO
     !   Edimp(it) = Edimp(it) + (dble(csm)-DM_EF)*(wkp*0.5)  ! 0.5 because wkp=2*wk/sumw                                                                                                                                                      
     !ENDDO
     allocate( Gij(nbands,nbands) )
     Gij(:,:)=0.d0
     DO i=1,nbands
        Gij(i,i) = E(i+DM_nemin-1)-DM_EF
     ENDDO
     call ProjectToLocal(Edimp0, Gij, DMFTU(:,:,:,1), 0.5*dble(wkp), iorbital, iSx, nbands, maxdim2, norbitals)
     deallocate( Gij )
     
     deallocate( Bii )
     DEALLOCATE( Winf, Wferm )
     DEALLOCATE( zek_inf )
     DEALLOCATE( zal_inf, zar_inf )
     deallocate( tmp )
     deallocate( tmp2 )
     deallocate( Sigma_inf )
     deallocate( SigmaTail1 )
     if (computeDirectC) deallocate(CPair,CStatic,Uactive,SigmaLocal, &
          SigmaInfLocal,SigmaTail1Local,CdiagUdagG,CdiagUdagGminus, &
          CdiagProduct,CdiagProductMinus)

     CALL GetLocalNds(dNds, Aweight, lgTrans, vnorm1, lg_deg, nbands, ntcix, nipc)     
  endif

  if (.not. matsubara) then
     Gdloc=0.0
     if (Qforce) then
        ! Real-axis counterpart: integrate the retarded/advanced difference
        ! of Sigma U^dagger G with the same log-bin weights used for rho.
        ! No projector pseudoinverse is needed on this branch either.
        allocate(Uactive(nbands,maxdim,ncix))
        Uactive=czero
        do iorb=1,norbitals
           icix=cix_orb(iorb)
           if(icix.le.0) cycle
           do ind=1,nindo(iorb)
              Uactive(:,iSx(ind,iorb),icix)= &
                   Uactive(:,iSx(ind,iorb),icix)+DMFTU(:,ind,iorb,1)
           enddo
        enddo
        AEweight=czero
     endif
     Aweight(:,:)=0
     logG=0
     logGD=0
     logG0=0
     !$OMP PARALLEL DO PRIVATE(Eij,zek,zal,zar,zw1,tmp,w1,w2,w3, i,fer,num,cc,gtc,tlogG,tlogG0,tlogGD,tdens, &
     !$OMP& SigmaLocal,CPair,CdiagUdagG,CdiagUdagGminus,CdiagProduct,CdiagProductMinus,iorb,icix,ind,ip)&
     !$OMP SHARED(gloc) SCHEDULE(STATIC,CHUNK) REDUCTION(+:Aweight,AEweight,projectorGradientDirect,dens,logG,logG0,logGD)
     do iom=1,npomega
        allocate(Eij(nbands,nbands),zek(nbands),zar(nbands,nbands),zal(nbands,nbands),tmp(nbands,nbands),zw1(nbands),w1(nbands,nbands),w2(nbands,nbands),w3(nbands,nbands))
        ! Total hamiltonian (H+Sigma) in Kohn-Sham basis
        Eij=0
        DO i=1,nbands
           Eij(i,i) = E(i+DM_nemin-1)  !---   preparing hamiltonian: epsk+sigma
        ENDDO
        ! Adding self-energy
        CALL AddSigma_optimized2(Eij, sigma(:,:,iom), STrans, csize, 1, nbands, ncix, maxsize)
        CALL eigsys(Eij, zek, zal, zar, nbands) ! eigenvalues and eigenvectors
        
        zw1=0.0
        !-------- Real axis calculation
        fer = ferm(omega(iom)/Temperature)

        DO num=DM_nemin,DM_nemax
           cc = log(dcmplx(abom(2,iom),wgamma)+DM_EF-zek(num-DM_nemin+1))-log(dcmplx(abom(1,iom),wgamma)+DM_EF-zek(num-DM_nemin+1))
           zw1(num-DM_nemin+1) = wkp*cc*fer
        ENDDO
        ! Aweight = ( Ar.zw1.Al - (Ar.zw1.Al)^+ )/(2*pi*i)
        tmp=0.0
        do num=1,nbands
           tmp(:,num) = zar(:,num)*zw1(num)
        enddo
        w1 = matmul(tmp,zal)
        w2 = transpose(w1)
        w3 = conjg(w2)
        w2 = -(w1-w3)/(2*pi*IMAG)
        Aweight(:,:) = Aweight(:,:) + w2(:,:)
        if (Qforce) then
           allocate(SigmaLocal(maxdim,maxdim,ncix),CPair(maxdim,nbands,ncix), &
                CdiagUdagG(maxdim,nbands),CdiagUdagGminus(maxdim,nbands), &
                CdiagProduct(maxdim,nbands),CdiagProductMinus(maxdim,nbands))
           CALL BuildLocalSigmaMatrix(SigmaLocal,sigma(:,:,iom),Sigind, &
                cixdim,maxdim,maxsize,ncix)
           CPair=czero
           CALL AccumulateDirectProjectorSpectral(CPair,w1,Uactive,SigmaLocal, &
                cixdim,nbands,maxdim,ncix,CdiagUdagG,CdiagUdagGminus, &
                CdiagProduct,CdiagProductMinus)
           do iorb=1,norbitals
              icix=cix_orb(iorb)
              if(icix.le.0) cycle
              do ind=1,nindo(iorb)
                 ip=iSx(ind,iorb)
                 projectorGradientDirect(:,ind,iorb)= &
                      projectorGradientDirect(:,ind,iorb)+CPair(ip,:,icix)
              enddo
           enddo
           ! zG-I = (epsilon+Sigma)G.  Apply the same piecewise-constant
           ! Hamiltonian as the log-bin integral, then its advanced partner.
           do num=1,nbands
              tmp(:,num)=zar(:,num)*zw1(num)*zek(num)
           enddo
           w2=matmul(tmp,zal)
           AEweight=AEweight-(w2-conjg(transpose(w2)))/(2*pi*IMAG)
           deallocate(SigmaLocal,CPair,CdiagUdagG,CdiagUdagGminus,CdiagProduct,CdiagProductMinus)
        endif

        
        gtc=0
        do i=1,nbands
           gtc = gtc + 1/(omega(iom)+DM_EF-zek(i)+gamma*IMAG)
        enddo
        
        CALL GetFreeEnergy(tlogG, tlogG0, tlogGD, tdens, Temperature, DM_EF, iom, omega, nomega, nbands, E, zek, DM_nemin, nume, wgamma)
        
        gloc(iom) = gloc(iom) + gtc * (0.5*wkp)
     
        dens  = dens  + tdens * wkpt
        logG  = logG  + tlogG * wkpt
        logG0 = logG0 + tlogG0* wkpt
        logGD = logGD + tlogGD* wkpt
        deallocate( Eij, zek, zar, zal, tmp, zw1, w1, w2, w3)
     enddo
     !$OMP END PARALLEL DO
     if(Qforce) deallocate(Uactive)
  endif
  if (stressMomentDiagnostic) deallocate(WSigPair,SigDirect,SigScratch)
  deallocate(Wsum, WEsum, w0inf) 

END SUBROUTINE cmp_dmft_weights2


! Expand the compact impurity self-energy labels into one dense local matrix
! for every correlated block.  Keeping this operation in local space is much
! cheaper than reconstructing the band-space embedded self-energy and is the
! purpose of accumulating the direct functional gradient C.
SUBROUTINE BuildLocalSigmaMatrix(localSigma,compactSigma,Sigind,cixdim, &
     maxdim,maxsize,ncix)
  IMPLICIT NONE
  INTEGER,intent(in) :: maxdim,maxsize,ncix,cixdim(ncix)
  INTEGER,intent(in) :: Sigind(maxdim,maxdim,ncix)
  COMPLEX*16,intent(in) :: compactSigma(maxsize,ncix)
  COMPLEX*16,intent(out) :: localSigma(maxdim,maxdim,ncix)
  INTEGER :: icix,ip,iq,it

  localSigma=0.d0
  do icix=1,ncix
     do iq=1,cixdim(icix)
        do ip=1,cixdim(icix)
           it=Sigind(ip,iq,icix)
           if (it.gt.0) localSigma(ip,iq,icix)=compactSigma(it,icix)
        enddo
     enddo
  enddo
END SUBROUTINE BuildLocalSigmaMatrix


! Add the positive and negative Matsubara partners of
!
!   C(iw) = Sigma_local(iw) U^+ G(iw)
!
! to CPair.  Notice that the negative partner is
! Sigma_local^+ U^+ G^+, not C(iw)^+: C is rectangular, and reversing the
! product order would not reproduce U*C=Sigma_embedded*G.
SUBROUTINE AccumulateDirectProjectorPair(CPair,G,Uactive,localSigma, &
     cixdim,nbands,maxdim,ncix,UdagG,UdagGminus,product,productMinus)
  IMPLICIT NONE
  INTEGER,intent(in) :: nbands,maxdim,ncix,cixdim(ncix)
  COMPLEX*16,intent(in) :: G(nbands,nbands), &
       Uactive(nbands,maxdim,ncix),localSigma(maxdim,maxdim,ncix)
  COMPLEX*16,intent(inout) :: CPair(maxdim,nbands,ncix)
  COMPLEX*16,intent(inout) :: UdagG(maxdim,nbands), &
       UdagGminus(maxdim,nbands),product(maxdim,nbands), &
       productMinus(maxdim,nbands)
  INTEGER :: icix,nloc
  COMPLEX*16,parameter :: cone=(1.d0,0.d0),czero=(0.d0,0.d0)

  do icix=1,ncix
     nloc=cixdim(icix)
     if (nloc.le.0) cycle
     CALL zgemm('C','N',nloc,nbands,nbands,cone, &
          Uactive(1,1,icix),nbands,G,nbands,czero,UdagG,maxdim)
     CALL zgemm('N','N',nloc,nbands,nloc,cone, &
          localSigma(1,1,icix),maxdim,UdagG,maxdim,czero,product,maxdim)

     CALL zgemm('C','C',nloc,nbands,nbands,cone, &
          Uactive(1,1,icix),nbands,G,nbands,czero,UdagGminus,maxdim)
     CALL zgemm('C','N',nloc,nbands,nloc,cone, &
          localSigma(1,1,icix),maxdim,UdagGminus,maxdim,czero, &
          productMinus,maxdim)
     CPair(1:nloc,:,icix)=CPair(1:nloc,:,icix)+ &
          product(1:nloc,:)+productMinus(1:nloc,:)
  enddo
END SUBROUTINE AccumulateDirectProjectorPair

! Retarded/advanced contour pair, using the already integrated weighted G.
SUBROUTINE AccumulateDirectProjectorSpectral(CPair,G,Uactive,localSigma, &
     cixdim,nbands,maxdim,ncix,UdagG,UdagGminus,product,productMinus)
  IMPLICIT NONE
  INTEGER,intent(in) :: nbands,maxdim,ncix,cixdim(ncix)
  COMPLEX*16,intent(in) :: G(nbands,nbands), &
       Uactive(nbands,maxdim,ncix),localSigma(maxdim,maxdim,ncix)
  COMPLEX*16,intent(inout) :: CPair(maxdim,nbands,ncix)
  COMPLEX*16,intent(inout) :: UdagG(maxdim,nbands), &
       UdagGminus(maxdim,nbands),product(maxdim,nbands), &
       productMinus(maxdim,nbands)
  INTEGER :: icix,nloc
  COMPLEX*16,parameter :: cone=(1.d0,0.d0),czero=(0.d0,0.d0)

  do icix=1,ncix
     nloc=cixdim(icix)
     if (nloc.le.0) cycle
     CALL zgemm('C','N',nloc,nbands,nbands,cone, &
          Uactive(1,1,icix),nbands,G,nbands,czero,UdagG,maxdim)
     CALL zgemm('N','N',nloc,nbands,nloc,cone, &
          localSigma(1,1,icix),maxdim,UdagG,maxdim,czero,product,maxdim)

     CALL zgemm('C','C',nloc,nbands,nbands,cone, &
          Uactive(1,1,icix),nbands,G,nbands,czero,UdagGminus,maxdim)
     CALL zgemm('C','N',nloc,nbands,nloc,cone, &
          localSigma(1,1,icix),maxdim,UdagGminus,maxdim,czero, &
          productMinus,maxdim)
     CPair(1:nloc,:,icix)=CPair(1:nloc,:,icix)+ &
          -(product(1:nloc,:)-productMinus(1:nloc,:))/cmplx(0.d0,2.d0*acos(-1.d0),8)
  enddo
END SUBROUTINE AccumulateDirectProjectorSpectral


! Determine the leading dynamic self-energy moment used by the Matsubara
! stress and free-energy tails.
!
! For a causal, sufficiently smooth impurity self-energy,
!
!   Sigma(z)-Sigma_infinity = Sigma_1/z + Sigma_2/z^2 + O(z^-3), z=i*omega.
!
! Multiplying by z gives a well-conditioned two-parameter fit,
!
!   z [Sigma(z)-Sigma_infinity] = Sigma_1 + Sigma_2/z.
!
! Real impurity data are not always in their asymptotic regime: a Monte Carlo
! tail can be noisy, or a failed iteration can even grow at high frequency.
! We therefore do not blindly trust the two-moment fit.  Each independent
! self-energy component is classified from its relative residual, stability
! under shortening the fit interval, agreement with the last-point estimate,
! and the fraction of points whose dynamic magnitude grows.  GOOD data use
! the fit; MARGINAL data use a short plateau average of z*(Sigma-Sigma_inf);
! BAD data retain the conservative historical last-point estimate.  The
! classification is diagnostic only and is printed once, at the first k point.
SUBROUTINE DetermineSelfEnergyTailMoment(moment1,sigma,sigmaInfinity, &
     omega,csize,nomega,maxsize,ncix,printDiagnostic)
  IMPLICIT NONE
  INTEGER,intent(in) :: nomega,maxsize,ncix,csize(ncix)
  REAL*8,intent(in) :: omega(nomega)
  COMPLEX*16,intent(in) :: sigma(maxsize,ncix,nomega)
  COMPLEX*16,intent(in) :: sigmaInfinity(maxsize,ncix)
  LOGICAL,intent(in) :: printDiagnostic
  COMPLEX*16,intent(out) :: moment1(maxsize,ncix)

  INTEGER :: icix,it,firstFit,nfit,firstShort,nshort,firstPlateau
  INTEGER :: iom,nIncrease,nComparison,quality
  REAL*8 :: residual,shortResidual,stability,endpointDifference
  REAL*8 :: increaseFraction,scale,tinyScale
  COMPLEX*16 :: fittedMoment,shortMoment,secondMoment
  COMPLEX*16 :: shortSecondMoment,endpointMoment,plateauMoment
  COMPLEX*16 :: dynamicSigma,z
  CHARACTER(len=8) :: qualityName

  tinyScale=100.d0*tiny(1.d0)
  moment1=0.d0

  do icix=1,ncix
     do it=1,csize(icix)
        ! Use at most 24 retained points, all from the last quarter of the
        ! available frequency range.  With a logarithmic mesh this keeps the
        ! fit local to the true high-frequency tail without relying on a
        ! hard-coded energy expressed in eV or Ry.
        firstFit=max(1,nomega-23)
        do while (firstFit.lt.nomega-5 .and. &
             omega(firstFit).lt.0.25d0*omega(nomega))
           firstFit=firstFit+1
        enddo
        nfit=nomega-firstFit+1

        CALL FitTwoSelfEnergyMoments(fittedMoment,secondMoment,residual, &
             sigma(it,icix,:),sigmaInfinity(it,icix),omega,firstFit, &
             nomega,nomega)

        ! A second, nested fit is a cheap sensitivity test.  A formally small
        ! residual is not sufficient when a slowly drifting, non-asymptotic
        ! tail can be fitted over one particular interval.
        nshort=max(6,(nfit+1)/2)
        firstShort=nomega-nshort+1
        CALL FitTwoSelfEnergyMoments(shortMoment,shortSecondMoment, &
             shortResidual,sigma(it,icix,:),sigmaInfinity(it,icix), &
             omega,firstShort,nomega,nomega)

        dynamicSigma=sigma(it,icix,nomega)-sigmaInfinity(it,icix)
        ! This is exactly the old safe estimate for a normal diagonal
        ! self-energy.  Keeping it real also prevents a noisy Re Sigma tail
        ! from manufacturing an anti-Hermitian diagonal moment.
        endpointMoment=dcmplx(-aimag(dynamicSigma)*omega(nomega),0.d0)

        firstPlateau=max(1,nomega-4)
        plateauMoment=0.d0
        do iom=firstPlateau,nomega
           z=dcmplx(0.d0,omega(iom))
           plateauMoment=plateauMoment+z* &
                (sigma(it,icix,iom)-sigmaInfinity(it,icix))
        enddo
        plateauMoment=plateauMoment/dble(nomega-firstPlateau+1)

        scale=max(abs(fittedMoment),abs(shortMoment),abs(endpointMoment), &
             tinyScale)
        stability=abs(fittedMoment-shortMoment)/scale
        endpointDifference=abs(fittedMoment-endpointMoment)/scale

        nIncrease=0
        nComparison=0
        do iom=firstFit+1,nomega
           nComparison=nComparison+1
           if (abs(sigma(it,icix,iom)-sigmaInfinity(it,icix)).gt. &
                1.10d0*abs(sigma(it,icix,iom-1)- &
                sigmaInfinity(it,icix))) nIncrease=nIncrease+1
        enddo
        increaseFraction=dble(nIncrease)/dble(max(1,nComparison))

        ! quality=2: GOOD, quality=1: MARGINAL, quality=0: BAD.
        if (residual.lt.0.05d0 .and. stability.lt.0.10d0 .and. &
             endpointDifference.lt.0.25d0 .and. &
             increaseFraction.lt.0.35d0) then
           quality=2
           qualityName='GOOD'
           moment1(it,icix)=fittedMoment
        else if (residual.lt.0.20d0 .and. stability.lt.0.35d0 .and. &
             increaseFraction.lt.0.50d0) then
           quality=1
           qualityName='MARGINAL'
           moment1(it,icix)=plateauMoment
        else
           quality=0
           qualityName='BAD'
           moment1(it,icix)=endpointMoment
        endif

        ! Guard the production result against NaN/overflow without requiring
        ! compiler-specific IEEE modules.  The comparisons are false for NaN.
        if (.not.(abs(moment1(it,icix)).le.huge(1.d0))) then
           moment1(it,icix)=endpointMoment
           quality=0
           qualityName='BAD'
        endif
        if (.not.(abs(moment1(it,icix)).le.huge(1.d0))) &
             moment1(it,icix)=0.d0

        if (printDiagnostic) then
           WRITE(6,'(A,1X,A,3I6,6(1X,ES12.4),4(1X,ES18.9))') &
                ':SELFENERGY_TAIL:',trim(qualityName),icix,it,nfit, &
                omega(firstFit),omega(nomega),residual,stability, &
                endpointDifference,increaseFraction,dble(moment1(it,icix)), &
                aimag(moment1(it,icix)),dble(secondMoment), &
                aimag(secondMoment)
        endif
     enddo
  enddo

CONTAINS

  ! Complex least-squares fit of y=a+b*x with
  ! y=i*omega*(Sigma-Sigma_inf), x=1/(i*omega).  The normal matrix is only
  ! 2x2, so an explicit Hermitian solution is both cheaper and clearer than
  ! calling a general LAPACK solver for every self-energy component.
  SUBROUTINE FitTwoSelfEnergyMoments(firstMoment,secondMoment,fitResidual, &
       component,staticPart,frequency,first,last,nall)
    IMPLICIT NONE
    INTEGER,intent(in) :: first,last,nall
    REAL*8,intent(in) :: frequency(nall)
    COMPLEX*16,intent(in) :: component(nall),staticPart
    COMPLEX*16,intent(out) :: firstMoment,secondMoment
    REAL*8,intent(out) :: fitResidual
    INTEGER :: i,npoint
    REAL*8 :: sumXX,determinant,normY,normResidual,fitScale
    COMPLEX*16 :: sumX,sumY,sumXY,x,y,yfit,z

    npoint=last-first+1
    sumX=0.d0
    sumXX=0.d0
    sumY=0.d0
    sumXY=0.d0
    do i=first,last
       z=dcmplx(0.d0,frequency(i))
       x=1.d0/z
       y=z*(component(i)-staticPart)
       sumX=sumX+x
       sumXX=sumXX+abs(x)**2
       sumY=sumY+y
       sumXY=sumXY+conjg(x)*y
    enddo
    determinant=dble(npoint)*sumXX-abs(sumX)**2
    if (npoint.lt.2 .or. determinant.le. &
         100.d0*tiny(1.d0)*max(1.d0,sumXX)) then
       firstMoment=0.d0
       secondMoment=0.d0
       fitResidual=huge(1.d0)
       return
    endif
    firstMoment=(sumY*sumXX-sumX*sumXY)/determinant
    secondMoment=(dble(npoint)*sumXY-conjg(sumX)*sumY)/determinant

    normY=0.d0
    normResidual=0.d0
    do i=first,last
       z=dcmplx(0.d0,frequency(i))
       x=1.d0/z
       y=z*(component(i)-staticPart)
       yfit=firstMoment+secondMoment*x
       normY=normY+abs(y)**2
       normResidual=normResidual+abs(y-yfit)**2
    enddo
    fitScale=max(sqrt(normY),sqrt(dble(npoint))*abs(firstMoment), &
         100.d0*tiny(1.d0))
    fitResidual=sqrt(normResidual)/fitScale
  END SUBROUTINE FitTwoSelfEnergyMoments

END SUBROUTINE DetermineSelfEnergyTailMoment


SUBROUTINE cmp_dmft_weights(Aweight, gloc, logG, logGD, logG0, dens, STrans, sigma, wkp, omega, vnorm1, DM_EF, Temperature, wgamma, gamma, csize, nbands, DM_nemin, DM_nemax, maxsize, ncix, nomega, npomega, nume, matsubara)
  USE defs,  ONLY: pi, IMAG
  USE xa,    ONLY: E
  USE muzero,ONLY: nomq, iomq, jomq, womq, abom, n0_om
  USE param, ONLY: CHUNK
  IMPLICIT NONE
  COMPLEX*16, intent(out):: Aweight(nbands,nbands), gloc(nomega)
  REAL*8, intent(out)    :: logG, logGD, logG0, dens
  COMPLEX*16, intent(in) :: STrans(maxsize,ncix,nbands,nbands)
  COMPLEX*16, intent(in) :: sigma(maxsize,ncix,nomega), wkp
  REAL*8, intent(in)     :: omega(nomega), vnorm1, DM_EF, Temperature, wgamma, gamma
  INTEGER, intent(in)    :: csize(ncix), nbands, DM_nemin, DM_nemax, maxsize, ncix, nomega, npomega, nume
  LOGICAL, intent(in)    :: matsubara
  ! functions
  Interface
     real*8 FUNCTION ferm(x)
       real*8, intent(in) :: x
     end FUNCTION ferm
     SUBROUTINE GetFreeEnergy(tlogG, tlogG0, tlogGD, tdens, Temperature, DM_EF, iom, omega, nomega, nbands, E, zek, DM_nemin, nume, gamma)
       REAL*8, intent(out) :: tlogG, tlogG0, tlogGD, tdens
       INTEGER, intent(in) :: iom, nomega, nbands, DM_nemin, nume
       REAL*8, intent(in)  :: Temperature, omega(nomega), E(nume), gamma, DM_EF
       COMPLEX*16, intent(in):: zek(nbands)
     end SUBROUTINE GetFreeEnergy
     real*8 FUNCTION FreeE0(Energy, Temp)
       real*8, intent(in) :: Energy, Temp
     end FUNCTION FreeE0
     SUBROUTINE AddSigma_optimized2(gij, sigma, STrans, csize, sign, nbands, ncix, maxsize)
       COMPLEX*16, intent(inout):: gij(nbands,nbands)
       COMPLEX*16, intent(in)   :: sigma(maxsize,ncix)
       COMPLEX*16, intent(in)   :: STrans(maxsize,ncix,nbands,nbands)
       INTEGER, intent(in)      :: csize(ncix), sign, nbands, ncix, maxsize
     end SUBROUTINE AddSigma_optimized2
     SUBROUTINE eigsys(ham, zek, evl, evr, ndim)
       COMPLEX*16, INTENT(inout):: ham(ndim,ndim)  ! Hamiltonian / overwritten
       COMPLEX*16, INTENT(out) :: zek(ndim)        ! Vector of eigenvalues 
       COMPLEX*16, INTENT(out) :: evl(ndim,ndim)   ! left eigenvectors
       COMPLEX*16, INTENT(out) :: evr(ndim,ndim)   ! right eigenvectors
       INTEGER, INTENT(in)     :: ndim             ! Dimension of hamiltonian
     end SUBROUTINE eigsys
  end interface
  ! locals
  COMPLEX*16 :: wsum(nbands,nbands), w0inf(nbands)
  COMPLEX*16, allocatable :: Eij(:,:), zek(:), zar(:,:), zal(:,:), tmp(:,:), p_tmp(:,:)
  REAL*8,     allocatable :: zek_inf(:)
  COMPLEX*16, allocatable :: zw1(:), w1(:,:), p_zek(:), p_zar(:,:), p_zal(:,:), zar_inf(:,:), zal_inf(:,:), w2(:,:), w3(:,:)
  COMPLEX*16, ALLOCATABLE :: s_inf(:,:), Winf(:,:), Wferm(:,:)
  COMPLEX*16 :: cc, gtc, cek, omn, w0, p_w0, e0, e1
  REAL*8  :: fer, tlogGD, tlogG, tlogG0, tdens, ek_inf, beta
  INTEGER :: i, iom, num, num0, j0, j1, iomt
  complex*16 :: czero, cone
  cone  = cmplx(1.d0,0.d0, 8)
  czero = cmplx(0.d0,0.d0, 8)
  
  beta = 1/Temperature
  
  if (matsubara) then
     allocate( zek(nbands), Eij(nbands,nbands) )
     ALLOCATE( zek_inf(nbands) )
     ALLOCATE( zar_inf(nbands,nbands), zal_inf(nbands,nbands) )

     !----------- Computes the infinite frequency limit
     Eij=0
     DO i=1,nbands
        Eij(i,i) = E(i+DM_nemin-1)  !---   preparing hamiltonian: epsk+sigma
     ENDDO
     ! Adding self-energy
     allocate( s_inf(maxsize,ncix) )
     s_inf(:,:) = dble(sigma(:,:,nomega))
     CALL AddSigma_optimized2(Eij, s_inf, STrans, csize, 1, nbands, ncix, maxsize)
     deallocate(s_inf)
     ! Here we might want to use Hermitian eigsys!!!!????
     CALL eigsys(Eij, zek, zal_inf, zar_inf, nbands) ! eigenvalues and eigenvectors
     zek_inf(:) = dble(zek(:))-DM_EF
     
     deallocate( zek, Eij )
     
     w0inf=0
     Wsum=0
     Aweight(:,:)=0
     logG0=0
     logGD=0
     
     !$OMP  PARALLEL DO SHARED(gloc)&
     !$OMP& PRIVATE(Eij,zek,zal,zar,tlogG0,tlogGD,tlogG,tdens,gtc,num,num0,ek_inf,cek,omn,w0,tmp)&
     !$OMP& SCHEDULE(STATIC,CHUNK) REDUCTION(+:w0inf,Wsum, dens,logG,logG0,logGD)
     do iom=1,n0_om-1
        ! Necessary allocation due to openMP
        allocate( Eij(nbands,nbands), zek(nbands), zar(nbands,nbands), zal(nbands,nbands), tmp(nbands,nbands) )  
        
        ! Total hamiltonian (H+Sigma) in Kohn-Sham basis
        Eij=0
        DO i=1,nbands
           Eij(i,i) = E(i+DM_nemin-1)  !---   preparing hamiltonian: epsk+sigma
        ENDDO
        ! Adding self-energy
        CALL AddSigma_optimized2(Eij, sigma(:,:,iom), STrans, csize, 1, nbands, ncix, maxsize)
        CALL eigsys(Eij, zek, zal, zar, nbands) ! eigenvalues and eigenvectors
        
        tlogG0=0.0
        tlogGD=0.0
        tdens=0.0
        gtc=0
        DO num=DM_nemin,DM_nemax
           num0 = num-DM_nemin+1
           ek_inf = zek_inf(num0)

           ! These frequencies are not yet in logarithmic mesh (Where all consecutive matsubara frequencies are present)
           cek = zek(num0)-DM_EF           ! eps_k-mu
           omn = dcmplx(0, omega(iom) )    
           w0inf(num0) = w0inf(num0) + 1./(omn-ek_inf)  ! w0inf = sum( 1/(i*om+mu-ek_inf) , om)
           w0          = 1./(omn-cek)
           tmp(:,num0) = zar(:,num0)*w0                 ! tmp = Ar * 1/(i*om+mu-ek)

           tlogGD = tlogGD - dble(log(-omn+cek))
           tlogG0 = tlogG0 - dble(log(-omn+ek_inf))
           tdens  = tdens  + dble( 1/(omn-cek)-1/(omn-ek_inf) )

           gtc = gtc + 1/(omega(iom)*IMAG+DM_EF-zek(num0))
        ENDDO
        
        ! Wsum = Ar * 1/(i*om+mu-ek) * Al
        ! Wsum(:,:) = Wsum(:,:) + matmul(tmp,zal)
        CALL zgemm('N','N',nbands,nbands,nbands,cone,tmp,nbands,zal,nbands,cone,Wsum,nbands)


        tlogGD = 2/beta*tlogGD
        tlogG0 = 2/beta*tlogG0
        tdens  = 2/beta*tdens
        !tlogG  = tlogGD - tlogG0


        gloc(iom) = gloc(iom) + gtc * (0.5*wkp)

        dens  = dens  + tdens * dble(wkp)*vnorm1
        !logG  = logG  + tlogG * dble(wkp)*vnorm1
        logG0 = logG0 + tlogG0* dble(wkp)*vnorm1
        logGD = logGD + tlogGD* dble(wkp)*vnorm1
        deallocate( Eij, zek, zar, zal, tmp )  
     enddo   ! iom loop
     !$OMP END PARALLEL DO


     allocate( Eij(nbands,nbands), zek(nbands), zar(nbands,nbands), zal(nbands,nbands), tmp(nbands,nbands) )
     
     ALLOCATE( p_tmp(nbands,nbands), p_zek(nbands), p_zar(nbands,nbands), p_zal(nbands,nbands) )
     p_zek(:) = 0.0   !zek(:)
     p_zal(:,:) = 0.0 !zal(:,:)
     p_zar(:,:) = 0.0 !zar(:,:)
     p_tmp(:,:) = 0.0
     
     do iom=n0_om,npomega
        ! Total hamiltonian (H+Sigma) in Kohn-Sham basis
        Eij=0
        DO i=1,nbands
           Eij(i,i) = E(i+DM_nemin-1)  !---   preparing hamiltonian: epsk+sigma
        ENDDO
        ! Adding self-energy
        CALL AddSigma_optimized2(Eij, sigma(:,:,iom), STrans, csize, 1, nbands, ncix, maxsize)
        CALL eigsys(Eij, zek, zal, zar, nbands) ! eigenvalues and eigenvectors

        tlogG0=0.0
        tlogGD=0.0
        tdens=0.0
        gtc=0
        !$OMP PARALLEL DO SHARED(w0inf,tmp,p_tmp)&
        !$OMP PRIVATE(num0,ek_inf,cek,omn,w0,j0,j1,e0,e1,iomt,i,p_w0)&
        !$OMP& SCHEDULE(STATIC,CHUNK) REDUCTION(+:tlogGD,tlogG0,tdens,gtc)
        DO num=DM_nemin,DM_nemax
           num0 = num-DM_nemin+1
           ek_inf = zek_inf(num0)
           
           ! Frequencies in logarithmic mesh
           j0 = jomq(iom-1)
           j1 = jomq(iom)
           e1 = zEk(num0)-DM_EF
           e0 = e1
           
           if (iom.GT.n0_om) then
              e0 = p_zEk(num0)-DM_EF
              iomt = iom-1
              p_w0=0
              ! Here we interpolate forward. Current and next frequency are needed.
              ! Since we do not know yet quantities in the next frequency, we 
              ! do the calculation for the previous frequency (iom-1).
              do i=1,nomq(iomt)
                 if (iomq(iomt,i).GT.jomq(iomt)) then
                    omn = dcmplx(0, (2*iomq(iomt,i)-1)*pi/beta)
                    cek = e0 + (e1-e0)*(iomq(iomt,i)-j0)/(j1-j0)
                    p_w0 = p_w0 + womq(iomt,i)*1./(omn-cek)
                    w0inf(num0) = w0inf(num0) + womq(iomt,i)*1./(omn-ek_inf)
                 
                    tlogGD = tlogGD - womq(iomt,i)*dble(log(-omn+cek))
                    tlogG0 = tlogG0 - womq(iomt,i)*dble(log(-omn+ek_inf))
                    tdens  = tdens  + womq(iomt,i)*dble(1/(omn-cek)-1/(omn-ek_inf))
                 endif
              enddo
              p_tmp(:,num0) = p_zar(:,num0)*p_w0
           endif

           
           iomt = iom
           w0=0
           ! Here we interpolate backward. Current and previous frequency are needed
           do i=1,nomq(iomt)
              if (iomq(iomt,i).LE.jomq(iomt)) then
                 omn = dcmplx(0, (2*iomq(iomt,i)-1)*pi/beta)
                 cek = e0 + (e1-e0)*(iomq(iomt,i)-j0)/(j1-j0)
                 w0 = w0 + womq(iomt,i)*1./(omn-cek)
                 w0inf(num0) = w0inf(num0) + womq(iomt,i)*1./(omn-ek_inf)
                 
                 tlogGD = tlogGD - womq(iomt,i)*dble(log(-omn+cek))
                 tlogG0 = tlogG0 - womq(iomt,i)*dble(log(-omn+ek_inf))
                 tdens  = tdens  + womq(iomt,i)*dble(1/(omn-cek)-1/(omn-ek_inf))
              endif
           enddo
           tmp(:,num0) =  zar(:,num0)*w0
           
           gtc = gtc + 1/(omega(iom)*IMAG+DM_EF-zek(num0))
        ENDDO
        !$OMP END PARALLEL DO

        ! Wsum = Ar * 1/(i*om+mu-ek) * Al
        ! Wsum(:,:) = Wsum(:,:) + matmul(tmp,zal)
        CALL zgemm('N','N',nbands,nbands,nbands,cone,tmp,nbands,zal,nbands,cone,Wsum,nbands)

        ! iom>n0_om :  Wsum = tmp * Al + p_tmp*p_zal
        if (iom.GT.n0_om) then
           ! Wsum(:,:) = Wsum(:,:) + matmul(p_tmp,p_zal)
           CALL zgemm('N','N',nbands,nbands,nbands,cone,p_tmp,nbands,p_zal,nbands,cone,Wsum,nbands)
        endif

        p_zek(:) = zek(:)
        p_zal(:,:) = zal(:,:)
        p_zar(:,:) = zar(:,:)

        tlogGD = 2/beta*tlogGD
        tlogG0 = 2/beta*tlogG0
        tdens  = 2/beta*tdens
        !tlogG  = tlogGD - tlogG0


        gloc(iom) = gloc(iom) + gtc * (0.5*wkp)

        dens  = dens  + tdens * dble(wkp)*vnorm1
        !logG  = logG  + tlogG * dble(wkp)*vnorm1
        logG0 = logG0 + tlogG0* dble(wkp)*vnorm1
        logGD = logGD + tlogGD* dble(wkp)*vnorm1
     enddo   ! iom loop
     
     DEALLOCATE( p_tmp, p_zek, p_zar, p_zal )
     
     ALLOCATE( Winf(nbands,nbands), Wferm(nbands,nbands) )
     
     do num0=1,nbands
        tmp(:,num0) =   zar_inf(:,num0)*w0inf(num0)              
     enddo
     Winf = matmul(tmp,zal_inf)


     tlogG=0
     !$OMP PARALLEL DO SHARED(tmp)&
     !$OMP& SCHEDULE(STATIC,CHUNK) REDUCTION(+:dens,tlogG)
     do num0=1,nbands
        tmp  (:,num0) =   zar_inf(:,num0)*ferm(zek_inf(num0)*beta)
        dens = dens + ferm(zek_inf(num0)*beta)*dble(wkp)*vnorm1
        tlogG = tlogG + (FreeE0(zek_inf(num0),Temperature)-FreeE0(E(num0+DM_nemin-1)-DM_EF,Temperature))*dble(wkp)*vnorm1
     enddo
     !$OMP END PARALLEL DO
     Wferm = matmul(tmp,zal_inf)
     
     tmp(:,:) = (Wsum(:,:) - Winf(:,:))*(2/beta) + Wferm(:,:)
     Aweight(:,:) = (0.5*wkp)*(tmp(:,:) + transpose(conjg(tmp(:,:))))
     
     logG = logGD - logG0 + tlogG
     WRITE(*,'(A,F12.6,1x,A,F12.6,1x,A,F12.6,1x,A,F12.6)') 'logGD-logG0=', logGD-logG0, 'tlogG=', tlogG, 'logG=', logG

     DEALLOCATE( Winf, Wferm )
     DEALLOCATE( zek_inf )
     DEALLOCATE( zal_inf, zar_inf )
     deallocate( Eij, zek, zar, zal, tmp )
  endif

  if (.not. matsubara) then
     Aweight(:,:)=0
     logG=0
     logGD=0
     logG0=0
     !$OMP PARALLEL DO PRIVATE(Eij,zek,zal,zar,zw1,tmp,w1,w2,w3, i,fer,num,cc,gtc,tlogG,tlogG0,tlogGD,tdens)&
     !$OMP SHARED(gloc) SCHEDULE(STATIC,CHUNK) REDUCTION(+:Aweight,dens,logG,logG0,logGD)
     do iom=1,npomega
        allocate(Eij(nbands,nbands),zek(nbands),zar(nbands,nbands),zal(nbands,nbands),tmp(nbands,nbands),zw1(nbands),w1(nbands,nbands),w2(nbands,nbands),w3(nbands,nbands))
        ! Total hamiltonian (H+Sigma) in Kohn-Sham basis
        Eij=0
        DO i=1,nbands
           Eij(i,i) = E(i+DM_nemin-1)  !---   preparing hamiltonian: epsk+sigma
        ENDDO
        ! Adding self-energy
        CALL AddSigma_optimized2(Eij, sigma(:,:,iom), STrans, csize, 1, nbands, ncix, maxsize)
        CALL eigsys(Eij, zek, zal, zar, nbands) ! eigenvalues and eigenvectors
        
        zw1=0.0
        !-------- Real axis calculation
        fer = ferm(omega(iom)/Temperature)

        DO num=DM_nemin,DM_nemax
           cc = log(dcmplx(abom(2,iom),wgamma)+DM_EF-zek(num-DM_nemin+1))-log(dcmplx(abom(1,iom),wgamma)+DM_EF-zek(num-DM_nemin+1))
           zw1(num-DM_nemin+1) = wkp*cc*fer
        ENDDO
        ! Aweight = ( Ar.zw1.Al - (Ar.zw1.Al)^+ )/(2*pi*i)
        tmp=0.0
        do num=1,nbands
           tmp(:,num) = zar(:,num)*zw1(num)
        enddo
        w1 = matmul(tmp,zal)
        w2 = transpose(w1)
        w3 = conjg(w2)
        w2 = -(w1-w3)/(2*pi*IMAG)
        Aweight(:,:) = Aweight(:,:) + w2(:,:)
        
        gtc=0
        do i=1,nbands
           gtc = gtc + 1/(omega(iom)+DM_EF-zek(i)+gamma*IMAG)
        enddo
        
        CALL GetFreeEnergy(tlogG, tlogG0, tlogGD, tdens, Temperature, DM_EF, iom, omega, nomega, nbands, E, zek, DM_nemin, nume, wgamma)
        
        gloc(iom) = gloc(iom) + gtc * (0.5*wkp)
     
        dens  = dens  + tdens * dble(wkp)*vnorm1
        logG  = logG  + tlogG * dble(wkp)*vnorm1
        logG0 = logG0 + tlogG0* dble(wkp)*vnorm1
        logGD = logGD + tlogGD* dble(wkp)*vnorm1
        deallocate( Eij, zek, zar, zal, tmp, zw1, w1, w2, w3)
     enddo
     !$OMP END PARALLEL DO
  endif
  
END SUBROUTINE cmp_dmft_weights

SUBROUTINE GetLocalNds(Nds, Aweight, lgTrans, vnorm1, lg_deg, nbands, ntcix, nip)
  IMPLICIT NONE
  COMPLEX*16, intent(in) :: lgTrans(nbands,nbands,ntcix,nip), Aweight(nbands,nbands)
  REAL*8, intent(in)     :: vnorm1
  INTEGER, intent(in)    :: lg_deg(ntcix), nbands, ntcix, nip
  REAL*8, intent(out)    :: Nds(nip,ntcix)
  ! locals
  COMPLEX*16 :: ProductTrace2
  INTEGER :: ip, it2
  Nds(:,:)=0.d0
  DO it2=1,ntcix  ! sum_{ij} Gij(i,j) * STrans(j,i)                                                                                                                                                                                      
     if (lg_deg(it2).eq.0) cycle
     do ip=1,nip
        Nds(ip,it2) = Nds(ip,it2) + dble(ProductTrace2(Aweight, lgTrans(:,:,it2,ip), nbands))*lg_deg(it2)*vnorm1
     enddo
  ENDDO
END SUBROUTINE GetLocalNds

REAL*8 Function ProductTrace(Sigmaij,Gij,nbands)
  IMPLICIT NONE
  COMPLEX*16, intent(in) :: Sigmaij(nbands,nbands), Gij(nbands,nbands)
  INTEGER, intent(in)    :: nbands
  !                                                                                                                                                                                                                                                                                           
  COMPLEX*16 :: smSG(nbands,nbands), Gji(nbands,nbands)
  REAL*8     :: res
  INTEGER    :: i, j
  complex*16 :: czero, cone
  cone  = cmplx(1.d0,0.d0, 8)
  czero = cmplx(0.d0,0.d0, 8)

  if (.False.) then
     CALL zgemm('N','N',nbands,nbands,nbands,cone,Sigmaij,nbands,Gij,nbands,cone,smSG,nbands)
     res=0
     do i=1,nbands
        res = res + real(smSg(i,i))
     enddo
  else
     Gji = transpose(Gij)
     smSG = Sigmaij*Gji
     res=0
     do i=1,nbands
        do j=1,nbands
           res = res + real(smSG(i,j))
        enddo
     enddo
  endif
  ProductTrace=res
  RETURN
end Function ProductTrace


COMPLEX*16 Function ProductTrace2(A, B, ndim)
  ! returns  sum_{ij} A(i,j)*B(i,j)                                                                                                                                                                                                                                                           
  IMPLICIT NONE
  COMPLEX*16, intent(in) :: A(ndim,ndim), B(ndim,ndim)
  INTEGER, intent(in)    :: ndim
  ! locals                                                                                                                                                                                                                                                                                    
  INTEGER :: i, j
  COMPLEX*16 :: res
  res=0
  do j=1,ndim
     do i=1,ndim
        res = res + A(i,j)*B(i,j)
     enddo
  enddo
  ProductTrace2 = res
  return
END Function ProductTrace2

SUBROUTINE ProductTrace3_(C, A,B,ndim)
  ! C(i) = sum_j A(j,i)*B(j,i)
  IMPLICIT NONE
  COMPLEX*16, intent(out):: C(ndim)
  COMPLEX*16, intent(in) :: A(ndim,ndim), B(ndim,ndim)
  INTEGER, intent(in)    :: ndim
  ! locals                                                                                                                                                                                                                                                                                    
  INTEGER :: i, j
  COMPLEX*16 :: res
  C=0
  do i=1,ndim
     res=0
     do j=1,ndim
        res = res + A(j,i)*B(j,i)
     enddo
     C(i) = res
  enddo
END SUBROUTINE ProductTrace3_

SUBROUTINE ProductTrace3(C, A,B,ndim)
  ! C(i) = sum_j A(i,j)*B(j,i)                                                                                                                                                                                                                                                                
  IMPLICIT NONE
  COMPLEX*16, intent(out):: C(ndim)
  COMPLEX*16, intent(in) :: A(ndim,ndim), B(ndim,ndim)
  INTEGER, intent(in)    :: ndim
  ! locals                                                                                                                                                                                                                                                                                    
  INTEGER :: i, j
  COMPLEX*16 :: res
  C=0
  do i=1,ndim
     res=0
     do j=1,ndim
        res = res + A(i,j)*B(j,i)
     enddo
     C(i) = res
  enddo
END SUBROUTINE ProductTrace3



  
SUBROUTINE ProjectToLocal(Gloc, Gij, DMFTU, wk, iorbital, iSx, nbands, maxdim2, norbitals)
  USE dmfts, ONLY: natom, nl, cix, ll, iso, maxdim, ncix, lmaxp
  IMPLICIT NONE
  COMPLEX*16, intent(inout):: Gloc(maxdim,maxdim,ncix)
  COMPLEX*16, intent(in)   :: Gij(nbands,nbands)
  COMPLEX*16, intent(in)   :: DMFTU(nbands,maxdim2,norbitals)
  INTEGER, intent(in)      :: iorbital(natom,lmaxp+1), iSx(maxdim2, norbitals), nbands, maxdim2, norbitals
  REAL*8, intent(in)       :: wk
  !                                                                                                                                                                                                                                                                                           
  INTEGER :: icase, jcase, l1case, l2case, icix, l1, l2, nind1, nind2, iorb1, iorb2, ind1, ind2, i1, i2
  COMPLEX*16, allocatable :: tmp(:,:), gd(:,:)
  complex*16 :: czero, cone
  cone  = cmplx(1.d0,0.d0, 8)
  czero = cmplx(0.d0,0.d0, 8)
  
  allocate( tmp(maxdim2,nbands), gd(maxdim2,maxdim2) )
  DO icase=1,natom
     do l1case=1,nl(icase)
        icix = cix(icase,l1case)
        if ( icix.EQ.0 ) CYCLE
        l1 = ll(icase,l1case)
        nind1 = (2*l1+1)*iso
        iorb1 = iorbital(icase,l1case)
        DO jcase=1,natom
           do l2case=1,nl(jcase)
              if ( cix(jcase,l2case).NE.icix ) CYCLE
              l2 = ll(jcase,l2case)
              nind2 = (2*l2+1)*iso
              iorb2 = iorbital(jcase,l2case)
              call zgemm('C','N', nind1, nbands, nbands, cone, DMFTU(:,:,iorb1), nbands, Gij(:,:), nbands, czero, tmp(:,:),maxdim2)
              call zgemm('N','N', nind1, nind2, nbands,  cone, tmp, maxdim2, DMFTU(:,:,iorb2),nbands, czero, gd, maxdim2)
              do ind1=1,nind1
                 do ind2=1,nind2
                    i1 = iSx(ind1,iorb1)
                    i2 = iSx(ind2,iorb2)
                    Gloc(i1,i2,icix) = Gloc(i1,i2,icix) + gd(ind1,ind2)*wk
                 enddo
              enddo
           enddo
        ENDDO
     enddo
  ENDDO
  deallocate( tmp, gd )
END SUBROUTINE ProjectToLocal
