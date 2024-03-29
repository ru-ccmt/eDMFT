! @Copyright 2007 Kristjan Haule
! 

SUBROUTINE cmp_dmft_weights2(Aweight, gloc, logG, dens,  Gdloc, Gdloc0, DMFTU, dNds, Edimp0, AEweight, nedim, STrans, sigma, s_oo, wkp, omega, vnorm1, DM_EF, Temperature, wgamma, gamma, csize, nbands, DM_nemin, DM_nemax, maxsize, ncix, nomega, npomega, nume, matsubara, lgTrans, lg_deg, ntcix, nipc, ikp, maxdim2, iorbital, iSx, norbitals)
  USE defs,  ONLY: pi, IMAG
  USE xa,    ONLY: E
  USE muzero,ONLY: nomq, iomq, jomq, womq, abom, n0_om
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
  LOGICAL, intent(in)    :: matsubara
  COMPLEX*16, intent(in) :: DMFTU(nbands,maxdim2,norbitals,nipc)
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
  COMPLEX*16, allocatable :: Eij(:,:), zek(:), zar(:,:), zal(:,:), tmp(:,:), Sigmaij(:,:), Sigma_inf(:,:), Gij(:,:), tmp2(:,:)
  REAL*8,     allocatable :: zek_inf(:)
  COMPLEX*16, allocatable :: zw1(:), w1(:,:), zar_inf(:,:), zal_inf(:,:), w2(:,:), w3(:,:), Dni(:)
  COMPLEX*16, ALLOCATABLE :: Winf(:,:), Wferm(:,:), Bii(:)
  REAL*8     :: logGD, logG0
  COMPLEX*16 :: cc, gtc, cek, omn, w0, wo, bg, sginf, sginf_!, cc0, cc1, cc2
  REAL*8  :: fer, tlogGD, tlogG, tlogG0, tdens, ek_inf, beta, csm, wkpt, C0, sginf2, second_order_correction, fermi!, diff
  INTEGER :: i, iom, num, num0, it, nomq_max, icix, j, ip, start_nipc
  LOGICAL :: Qforce
  complex*16 :: czero, cone
  cone  = cmplx(1.d0,0.d0, 8)
  czero = cmplx(0.d0,0.d0, 8)
  
  start_nipc=2
  
  Qforce=.false.
  if (nedim.eq.nbands) Qforce=.true.

  if (Qforce) then
     allocate( WEsum(nedim,nedim) ) ! Unfortunately we need this for $OMP command to work.
  else
     allocate( WEsum(1,1) )
  endif
  allocate(Wsum(nbands,nbands), w0inf(nbands))

  beta = 1/Temperature
  wkpt = dble(wkp)*vnorm1  
  if (matsubara) then
     C0=1.d0
     allocate( Sigmaij(nbands,nbands) )
     !  This will calculate the first moment of self-energy, i.e., Sigma= sigma_winf/iw + ...
     sigma_winf=0
     DO icix=1,ncix
        do it=1,csize(icix)
           sigma_winf(it,icix) = -aimag(sigma(it,icix,nomega))*omega(nomega)
        enddo
     ENDDO
     Sigmaij=0
     CALL AddSigma_optimized2(Sigmaij, sigma_winf, STrans, csize, 1, nbands, ncix, maxsize)

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
           ! Frequencies in logarithmic mesh
           ! Here we interpolate backward. Current and previous frequency are needed
           do i=1,nomq(iom)
              omn = IMAG*(2*iomq(iom,i)-1)*pi/beta
              w0 = w0 + womq(iom,i)*1./(omn-cek)
              cc = womq(iom,i)*1./(omn-ek_inf)
              wo = wo + cc
              bg = bg + cc/(omn-C0)
              tlogGD = tlogGD - womq(iom,i)*dble(log(-omn+cek))
              tlogG0 = tlogG0 - womq(iom,i)*dble(log(-omn+ek_inf))
              tdens  = tdens  + womq(iom,i)*dble(1/(omn-cek)-1/(omn-ek_inf))
           enddo
           w0inf(num0) = w0inf(num0) + wo
           tmp(:,num0) =  zar(:,num0)*w0
           sginf_ = sginf_ + bg*Bii(num0)

           Dni(num0) = 1./(IMAG*omega(iom)-cek)
        ENDDO
        !$OMP END PARALLEL DO

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
     
     deallocate( Sigma_inf )
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
        AEweight(:,:) = (0.5*wkp)*(tmp(:,:) + transpose(conjg(tmp(:,:))))
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

     CALL GetLocalNds(dNds, Aweight, lgTrans, vnorm1, lg_deg, nbands, ntcix, nipc)     
  endif

  if (.not. matsubara) then
     Gdloc=0.0
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
     
        dens  = dens  + tdens * wkpt
        logG  = logG  + tlogG * wkpt
        logG0 = logG0 + tlogG0* wkpt
        logGD = logGD + tlogGD* wkpt
        deallocate( Eij, zek, zar, zal, tmp, zw1, w1, w2, w3)
     enddo
     !$OMP END PARALLEL DO
  endif
  deallocate(Wsum, WEsum, w0inf) 

END SUBROUTINE cmp_dmft_weights2


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
