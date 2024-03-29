SUBROUTINE Read_Vec_Spin(E, wgh, As, As_lo, k3, k3lo, bk3, bk3lo, nemin, nemax, DM_nemin, DM_nemax, more_kpoints, n0, emin, emax, DM_Emin, DM_Emax, iso, nmat, nume, nnlo)
  use dmfts, ONLY: Qcomplex, projector
  use param, ONLY: nemin0, nemax0
  IMPLICIT NONE
  !INTEGER, intent(in)  :: ikp
  REAL*8,  intent(out) :: E(nume), wgh !, weight(nume)
  COMPLEX*16,intent(out) :: As(nmat,nume,iso)
  COMPLEX*16,intent(out) :: As_lo(nnlo,nume,iso)
  INTEGER, intent(out) :: k3(3,nmat), k3lo(3,nnlo)
  REAL*8,  intent(out) :: bk3(3,nmat), bk3lo(3,nnlo)
  INTEGER, intent(out) :: nemin, nemax, DM_nemin, DM_nemax
  LOGICAL, intent(out) :: more_kpoints
  INTEGER, intent(out) :: n0
  REAL*8,  intent(in)  :: emin, emax, DM_Emin, DM_Emax
  INTEGER, intent(in)  :: iso, nmat, nume, nnlo
  ! locals
  REAL*8  :: As_tmp(nmat)
  INTEGER :: i, is, itape, ios, ne, num
  REAL*8  :: s,t,z
  CHARACTER*10 :: bname
  
  DO is=1,iso
     itape=8+is
     READ(itape,IOSTAT=ios) s,t,z,bname,n0,ne,wgh
     more_kpoints=.FALSE.
     IF (ios /= 0) RETURN
     more_kpoints=.TRUE.
     
     READ(itape,IOSTAT=ios) (k3(1,i),k3(2,i),k3(3,i),i=1,n0-nnlo),(k3lo(1,i),k3lo(2,i),k3lo(3,i),i=1,nnlo)
     
     DO i=1,n0-nnlo
        bk3(1,i)=(s+k3(1,i))
        bk3(2,i)=(t+k3(2,i))
        bk3(3,i)=(z+k3(3,i))
     ENDDO
     DO i=1,nnlo
        bk3lo(1,i)=(s+k3lo(1,i))
        bk3lo(2,i)=(t+k3lo(2,i))
        bk3lo(3,i)=(z+k3lo(3,i))
     ENDDO
     
     nemin=1
     nemax=0
     DM_nemin=1
     DM_nemax=0
     num=0
     DO WHILE(num.NE.NE)
        READ(itape,IOSTAT=ios) num,E(num)
        if (Qcomplex) then
           READ(itape) (As(i,num,is),i=1,n0)
        else
           READ(itape) (As_tmp(i),i=1,n0)
           As(:n0,num,is) = As_tmp(:n0)
        endif
        As_lo(:nnlo,num,is) = As(n0-nnlo+1:n0,num,is)
        IF(e(num).LT.emin) nemin=nemin+1                           
        IF(e(num).LE.emax) nemax=nemax+1
        if (abs(projector).LT.4) then
           IF(e(num).LT.DM_Emin) DM_nemin=DM_nemin+1                           
           IF(e(num).LE.DM_Emax) DM_nemax=DM_nemax+1     
        !else
        !   DM_nemin=nemin0
        !   DM_nemax=nemax0
        endif
     ENDDO
     if (abs(projector).GE.4) then
        DM_nemin=max(nemin0,nemin)
        DM_nemax=min(nemax0,nemax)
     endif
     
  ENDDO
  RETURN
END SUBROUTINE Read_Vec_Spin


SUBROUTINE Read_Vec_Spin_DontStore(nemin, DM_nemin, DM_nemax, more_kpoints, n0, emin, emax, nnlo, kmax)
  use dmfts, ONLY: DM_Emin, DM_Emax, Qcomplex, projector!, iso
  use param, ONLY: nemin0, nemax0, nmat
  IMPLICIT NONE
  INTEGER, intent(out) :: nemin, DM_nemin, DM_nemax, n0
  LOGICAL, intent(out) :: more_kpoints
  INTEGER, intent(inout):: kmax(3)
  REAL*8, intent(in)   :: emin, emax
  INTEGER, intent(in)  :: nnlo
  REAL*8     :: wAsr
  COMPLEX*16 :: wAsc
  REAL*8   :: wE, wgh
  INTEGER :: nemax!, wkx, wky, wkz, 
  ! locals
  INTEGER :: i, itape, ios, ne, num
  REAL*8  :: s,t,z
  CHARACTER*10 :: bname
  INTEGER, allocatable :: Ks(:,:)
  
  allocate(Ks(3,nmat))
  
  itape=9
  READ(itape,IOSTAT=ios) s,t,z,bname,n0,ne,wgh
  more_kpoints=.FALSE.
  IF (ios /= 0) RETURN
  more_kpoints=.TRUE.
  
  !READ(itape) (wkx,wky,wkz,i=1,n0)
  READ(itape) (Ks(1,i),Ks(2,i),Ks(3,i),i=1,n0)
  
  nemin=1
  nemax=0
  DM_nemin=1
  DM_nemax=0
  num=0
  DO WHILE(num.NE.NE)
     READ(itape) num,wE
     if (Qcomplex) then
        READ(itape) (wAsc,i=1,n0)
     else
        READ(itape) (wAsr,i=1,n0)
     endif
     IF(wE.LT.emin) nemin=nemin+1                           
     IF(wE.LE.emax) nemax=nemax+1
     if (abs(projector).LT.4) then
        IF(wE.LT.DM_Emin) DM_nemin=DM_nemin+1                           
        IF(wE.LE.DM_Emax) DM_nemax=DM_nemax+1     
     endif
  ENDDO
  if (abs(projector).GE.4) then
     DM_nemin=max(nemin0,nemin)
     DM_nemax=min(nemax0,nemax)
  endif

  DO i=n0-nnlo,1,-1
     kmax(1)=MAX(Ks(1,i),kmax(1))
     kmax(2)=MAX(Ks(2,i),kmax(2))
     kmax(3)=MAX(Ks(3,i),kmax(3))
  ENDDO
  
  deallocate(Ks)
  RETURN
END SUBROUTINE Read_Vec_Spin_DontStore

SUBROUTINE CompressSigmaTransformation1(STrans, DMFTU, Sigind, iSx, cix, iorbital, ll, nl, natom, iso, ncix, maxdim, maxdim2, lmaxp, norbitals, nbands, maxsize)
  IMPLICIT NONE
  COMPLEX*16, intent(out) :: STrans(maxsize,ncix,nbands,nbands)
  COMPLEX*16, intent(in)  :: DMFTU(nbands,maxdim2,norbitals)
  INTEGER, intent(in)     :: iSx(maxdim2, norbitals)
  INTEGER, intent(in)     :: Sigind(maxdim,maxdim,ncix), cix(natom,4)!, csize(ncix)
  INTEGER, intent(in)     :: iorbital(natom,lmaxp+1), ll(natom,4), nl(natom)
  INTEGER, intent(in)     :: natom, iso, ncix, maxdim, maxdim2, lmaxp, norbitals, nbands, maxsize
  !----- locals
  INTEGER    :: i, j, icase, jcase, l1case, l2case, iorb1, iorb2, it, ind1, nind1, ind2, nind2, icix, l1, l2
  STrans=0
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
              do ind1=1,nind1
                 do ind2=1,nind2
                    it = Sigind( iSx(ind1,iorb1),iSx(ind2,iorb2), icix )
                    if (it.gt.0) then
                       DO i=1,nbands                     ! over bands-1
                          do j=1,nbands                  ! over bands-2
                             STrans(it,icix,i,j) = STrans(it,icix,i,j) + conjg(DMFTU(j,ind2,iorb2))*DMFTU(i,ind1,iorb1)
                          enddo
                       ENDDO
                    endif
                 enddo
              enddo
           enddo
        ENDDO
     enddo
  ENDDO
END SUBROUTINE CompressSigmaTransformation1

SUBROUTINE CompressSigmaTransformation2(STrans, lgTrans, lg_deg, DMFTU, Sigind, Sigind_orig, iSx, cix, iorbital, ll, nl, natom, iso, ncix, maxdim, maxdim2, lmaxp, norbitals, nbands, maxsize, ntcix, nipc)
  IMPLICIT NONE
  COMPLEX*16, intent(out) :: STrans(maxsize,ncix,nbands,nbands), lgTrans(nbands,nbands,ntcix,nipc)
  INTEGER, intent(out)    :: lg_deg(ntcix)
  COMPLEX*16, intent(in)  :: DMFTU(nbands,maxdim2,norbitals,nipc)
  INTEGER, intent(in)     :: iSx(maxdim2, norbitals)
  INTEGER, intent(in)     :: Sigind(maxdim,maxdim,ncix), Sigind_orig(maxdim,maxdim,ncix), cix(natom,4)!, csize(ncix)
  INTEGER, intent(in)     :: iorbital(natom,lmaxp+1), ll(natom,4), nl(natom)
  INTEGER, intent(in)     :: natom, iso, ncix, maxdim, maxdim2, lmaxp, norbitals, nbands, maxsize, ntcix, nipc
  !----- locals                                                                                                                                                                                                                               
  INTEGER    :: i, j, icase, jcase, l1case, l2case, iorb1, iorb2, it, ind1, nind1, ind2, nind2, icix, l1, l2, it2
  COMPLEX*16 :: cc, dd(3), imag
  ! lgTrans(nbands,nbands,ntcix,nipc)
  imag = cmplx(0.d0,1.d0,8)
  STrans=0
  lgTrans=0
  lg_deg=0
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
              do ind1=1,nind1
                 do ind2=1,nind2
                    it = Sigind( iSx(ind1,iorb1),iSx(ind2,iorb2), icix )
                    it2 = Sigind_orig( iSx(ind1,iorb1),iSx(ind2,iorb2), icix )
                    if (it.gt.0) then
                       lg_deg(it2)=lg_deg(it2)+1
                       DO i=1,nbands                     ! over bands-1                                                                                                                                                                       
                          do j=1,nbands                  ! over bands-2                                                                                                                                                                       
                             cc = conjg(DMFTU(j,ind2,iorb2,1))*DMFTU(i,ind1,iorb1,1)
                             STrans(it,icix,i,j) = STrans(it,icix,i,j) + cc
                             !!! Even more condensed index for computing Tr(log(Gloc))                                                                                                                                 
                             lgTrans(i,j,it2,1) = lgTrans(i,j,it2,1) + conjg(cc)
                             if (nipc.eq.4) then
                                !!! Condensed index for computing force
                                dd(1:3) = (conjg(DMFTU(i,ind1,iorb1,2:4))*DMFTU(j,ind2,iorb2,1) - conjg(DMFTU(i,ind1,iorb1,1))*DMFTU(j,ind2,iorb2,2:4))*imag
                                lgTrans(i,j,it2,2:4) = lgTrans(i,j,it2,2:4) + dd(1:3)
                             endif
                          enddo
                       ENDDO
                       ! Strans  <= DMFTU * Sigma * DMFTU^+
                       ! lgTrans <= DMFTU^+  G * DMFTU
                    endif
                 enddo
              enddo
           enddo
        ENDDO
     enddo
  ENDDO

  if (.False.) then
     do it2=1,ntcix
        dd=0
        do i=1,nbands
           dd(:) = dd(:) + lgTrans(i,i,it2,2:4)
        enddo
        write(6,'(A,I3,1x,A,I2,1x,A,3F10.4)') 'ic=', it2, 'lg_deg=', lg_deg(it2), 'sum=', dble(dd(:))
     enddo
     !DO i=1,nbands
     !   DO it2=1,ntcix
     !      write(6,'(A,I3,1x,A,I4,1x,A,2F10.4,1x,A,2F10.4,1x,A,2F10.4,1x,A,2F10.4)') 'ic=', it2, 'ib=', i, 'c=', lgTrans(i,i,it2,1), 'x=', lgTrans(i,i,it2,2), 'y=', lgTrans(i,i,it2,3), 'z=', lgTrans(i,i,it2,4)
     !   ENDDO
     !ENDDO
  endif
  
  do it=1,ntcix
     if (lg_deg(it).NE.0) lgTrans(:,:,it,:) = lgTrans(:,:,it,:)/lg_deg(it)
  enddo
END SUBROUTINE CompressSigmaTransformation2

SUBROUTINE AddSigma_optimized2(gij, sigma, STrans, csize, sign, nbands, ncix, maxsize)
  !-- Add's self-energy to inverse of the Green's function in band representation
  !-- 50% of all time is spend in this subroutine, the rest 50 % in CmpGk
  IMPLICIT NONE
  COMPLEX*16, intent(inout):: gij(nbands,nbands)
  COMPLEX*16, intent(in)   :: sigma(maxsize,ncix)
  COMPLEX*16, intent(in)  :: STrans(maxsize,ncix,nbands,nbands)
  INTEGER, intent(in)      :: csize(ncix), sign, nbands, ncix, maxsize
  !----- locals
  COMPLEX*16 :: Sigmaij
  INTEGER    :: i, j, icix, it
  !$OMP PARALLEL DO SHARED(gij) PRIVATE(j,Sigmaij,icix,it) SCHEDULE(STATIC)
  DO i=1,nbands                     ! over bands-1
     do j=1,nbands                  ! over bands-2
        Sigmaij=0
        DO icix=1,ncix
           do it=1,csize(icix)
              Sigmaij = Sigmaij + STrans(it,icix,i,j)*sigma(it,icix)
           enddo
        ENDDO
        gij(i,j) = gij(i,j) + sign*Sigmaij
     enddo
  ENDDO
  !$OMP END PARALLEL DO
END SUBROUTINE AddSigma_optimized2


SUBROUTINE Get_A_Dimensions(norbitals, maxdim2, lmaxp, iso, natom, nl, ll)
  IMPLICIT NONE
  INTEGER, intent(out) :: norbitals, maxdim2, lmaxp
  INTEGER, intent(in)  :: iso, natom, nl(natom), ll(natom,4) 
  ! locals
  INTEGER::  icase, lcase, l1, nind
  norbitals=0 !-- number of atom/l cases, called orbitals -------------------------------------!
  maxdim2=0
  lmaxp=0
  do icase=1,natom
     do lcase=1,nl(icase)
        norbitals = norbitals+1
        l1 = ll(icase,lcase)
        nind=(2*l1+1)*iso
        maxdim2 = max(maxdim2, nind )
        lmaxp = max(l1,lmaxp)
     enddo
  enddo
END SUBROUTINE Get_A_Dimensions



SUBROUTINE Set_A_Arrays(iorbital, nindo, cix_orb, iSx, noccur, cixdim, cfX, CF, norbitals, iso, natom, nl, ll, cix, ncix, maxdim, maxdim2, lmaxp, Sigind, csize, maxsize)
  USE com_mpi,ONLY: Qprint!, myrank, master
  IMPLICIT NONE
  INTEGER, intent(out)   :: iorbital(natom,lmaxp+1), nindo(norbitals), cix_orb(norbitals), iSx(maxdim2,norbitals), noccur(maxsize,ncix), cixdim(ncix)
  COMPLEX*16, intent(out):: cfX(maxdim2,maxdim2,norbitals,norbitals)
  COMPLEX*16, intent(in) :: CF(maxdim,maxdim,ncix)
  INTEGER, intent(in)    :: norbitals, iso, natom, nl(natom), ll(natom,4), cix(natom,4), ncix, maxdim, maxdim2, lmaxp, Sigind(maxdim,maxdim,ncix), csize(ncix), maxsize
  ! locals
  INTEGER :: iorb, icase, lcase, l1, nind, icix
  INTEGER :: iorb1, iorb2, nind1, nind2, ip1, ind1, ind2, ip, iq, it

  iorbital=0
  cixdim=0
  iorb=0
  do icase=1,natom
     do lcase=1,nl(icase)
        iorb = iorb+1
        iorbital(icase,lcase)=iorb  !-- index to the orbital number ----!
        l1 = ll(icase,lcase)
        nind=(2*l1+1)*iso
        nindo(iorb) = nind
        icix = cix(icase,lcase)
        cix_orb(iorb) = icix
        do ip1=1,nind
           cixdim(icix) = cixdim(icix) + 1
           iSx(ip1,iorb) = cixdim(icix)
           if (Qprint) WRITE(6,'(A,7I5)') 'icase,lcase,icix,iorb,nind1,ip1,iSx=', icase, lcase, icix, iorb, nind, ip1,cixdim(icix)
        enddo
     enddo
  enddo

  noccur=0
  DO icix=1,ncix
     DO ip=1,cixdim(icix)
        DO iq=1,cixdim(icix)
           it = Sigind(ip,iq,icix)
           if (it.gt.0)  noccur(it,icix) = noccur(it,icix) + 1
        ENDDO
     ENDDO
  ENDDO
  
  if (Qprint) then
     do icix=1,ncix
        WRITE(6,'(A,I2,A,I2)') 'cixdim(', icix, ')=', cixdim(icix)
        DO it=1,csize(icix)
           WRITE(6,'(A,I2,A,I2,A,I4)') 'noccur(',it,',',icix,')=', noccur(it,icix)
        ENDDO
     enddo

     do iorb=1,norbitals
        do ip1=1,nindo(iorb)
           WRITE(6,'(A,I2,A,I2,A,I3)') 'iSx(', ip1, ',', iorb, ')=', iSx(ip1,iorb)
        enddo
     enddo
     do iorb=1,norbitals
        WRITE(6,'(A,I2,A,I3)') 'nindo(', iorb, ')=', nindo(iorb)
     enddo
  end if

  cfX=0.d0
  DO iorb1=1,norbitals
     nind1 = nindo(iorb1)
     if ( cix_orb(iorb1).EQ.0 ) CYCLE
     DO iorb2=1,norbitals
        nind2 = nindo(iorb2)
        if ( cix_orb(iorb1).NE.cix_orb(iorb2) ) CYCLE
        icix = cix_orb(iorb1)
        do ind1=1,nind1
           do ind2=1,nind2
              cfX(ind1,ind2,iorb1,iorb2) = CF( iSx(ind1,iorb1),iSx(ind2,iorb2), icix )
           enddo
        enddo
     ENDDO
  ENDDO


  if (Qprint) then
     do iorb1=1,norbitals
        nind1 = nindo(iorb1)
        do iorb2=1,norbitals
           nind2 = nindo(iorb2)
           if ( cix_orb(iorb1).EQ.0 ) CYCLE
           if ( cix_orb(iorb1).NE.cix_orb(iorb2) ) CYCLE
           WRITE(6,'(A,I2,A,I2)') 'CF for iorb1=', iorb1, ' iorb2=', iorb2
           do ind1=1,nind1
              do ind2=1,nind2
                 WRITE(6,'(f14.8,1x,f14.8,3x)', advance='no') dble(cfX(ind1,ind2,iorb1,iorb2)), aimag(cfX(ind1,ind2,iorb1,iorb2))
              enddo
              WRITE(6,*)
           enddo
        enddo
     enddo
  endif
  
END SUBROUTINE Set_A_Arrays

SUBROUTINE Build_DMFT_Projector(DMFTU, cfX, Rspin, iorbital, norbitals, nipc, n0, nnlo, nbands, cix_orb, nindo, DM_nemin, DM_nemax, maxdim2, lmaxp)
  USE w_atpar
  USE dmfts, ONLY: iso, projector, natom, nl, ll, cix, iatom, shft, isort, crotloc!, lmaxp, Sigind, Qcomplex, ncix, maxdim, maxsize, DM_Emin, DM_Emax, csize, CF
  USE defs, ONLY: PI
  USE param, ONLY: lomax, nloat, iblock, nmat, nume
  USE lo, ONLY: nlo, loor, ilo, lapw!, nlov, nlon, alo, blo, clo, pi12lo, pe12lo, pr12lo!, elo_store, elo, rlo
  USE atspdt,ONLY: P, DP, PE, DPE!, PEI!, e_store, el
  USE xa, ONLY: fj, dfj
  USE xa3, ONLY: bk3, As, k3
  USE structure, ONLY: rotij, tauij, BR1, POS, RMT, VOL, rotloc !, mult
  USE p_project, ONLY: phi_jl, al_ucase, dri, rix_mat!, rix, max_lcase
  USE mod_lo, ONLY: lomain
  !USE w_atpar, ONLY: ri_mat
  USE Forces, ONLY: forcea!, Qforce
  !use com, ONLY: nat
  IMPLICIT NONE
  complex*16, intent(out)  :: DMFTU(nbands,maxdim2,norbitals,nipc)
  complex*16, intent(in)   :: cfX(maxdim2,maxdim2,norbitals,norbitals), Rspin(2,2,norbitals)
  INTEGER, intent(in)      :: nipc, n0, nnlo, DM_nemin, DM_nemax, iorbital(natom,lmaxp+1), maxdim2, nbands, cix_orb(norbitals), nindo(norbitals), norbitals, lmaxp
  ! locals
  COMPLEX*16, allocatable :: alm(:,:),blm(:,:),clm(:,:,:)
  complex*16 :: h_yl((lmaxp+1)*(lmaxp+1),iblock), h_alyl(2*lmaxp+1,iblock), h_blyl(2*lmaxp+1,iblock)
  INTEGER   :: icase, lcase, latom, jatom, lfirst, is, ii, i3, l, M, index, ibb, lda, ldb, ldc, jlo, jlop, iorb, ind, i, nind, icix, num1, nm1, m1, iorb1, nind1, nind2, iorb2, i_h_k
  REAL*8    :: BKROT2(3), BKROT3(3), BKRLOC(3), sm1, sm2
  REAL*8    :: ARG1, ARG2, ARGT2, rmt2, h_al(iblock), h_bl(iblock), dff, TWOPI
  COMPLEX*16:: FAC, PHSHEL, CFAC, cc, dd(3), sm2c
  COMPLEX*16:: YL((lmaxp+1)*(lmaxp+1))
  COMPLEX*16, ALLOCATABLE :: URx(:,:,:,:), tmp(:,:), alm0(:,:)
  COMPLEX*16, PARAMETER:: IMAG=   (0.0D0,1.0D0)
  COMPLEX*16 :: tuu, tud, tdu, tdd
  COMPLEX*16, allocatable :: Uu(:,:), Ud(:,:)
  INTEGER :: no1, isize_, iind, lmmax, ip, lmaxplmaxp, lomaxlomax
  LOGICAL :: nonzero_interst
  ! For p_interstitial
  COMPLEX*16, allocatable :: al_interstitial(:,:,:,:), h_interstitial(:,:)
  ! For Force
  REAL*8, allocatable     :: h_k(:,:)
  complex*16, allocatable :: h_ablyl_hk(:,:), aalm(:,:,:), bblm(:,:,:), cclm(:,:,:,:)
  LOGICAL :: nonzero_shft, Qforce_j
  REAL*8  :: rotloc_x_BR1(3,3), rotloc_x_BR1_x_rotij(3,3), crotloc_x_BR1(3,3)
  complex*16 :: czero, cone
  
  TWOPI=2.0*PI
  lmaxplmaxp = (lmaxp+1)*(lmaxp+1)
  lomaxlomax = (lomax+1)*(lomax+1)
  cone  = cmplx(1.d0,0.d0, 8)
  czero = cmplx(0.d0,0.d0, 8)

  allocate( alm(lmaxplmaxp,nume), blm(lmaxplmaxp,nume), clm(lomaxlomax, nume, nloat) )

  if (nipc.eq.4) then
     allocate( aalm(lmaxplmaxp,nbands,3), bblm(lmaxplmaxp,nbands,3), cclm(lomaxlomax,nbands,nloat,3) )
  else
     allocate( aalm(1,1,1), bblm(1,1,1), cclm(1,1,1,1) ) ! the arrays need to be always allocated
  endif
  allocate( URx(DM_nemax-DM_nemin+1,maxdim2,norbitals,nipc) )
  
  URx = 0.d0
  !-------------------------------------------------!
  !--------- Building DMFT Transformation ----------!
  !-------------------------------------------------!
  do icase=1,natom
     latom = iatom(icase)   ! The succesive number of atom (all atoms counted)
     jatom = isort(latom)   ! The sort of the atom  ( == w_jatom(iucase) )

     Qforce_j = (nipc.eq.4) .AND.forcea(0,jatom)
     
     CALL retrieve_w_atpar(jatom,lfirst,lmmax)
     
     !!! This part is different for projector 5
     nonzero_interst=.false.
     if (abs(projector).ge.5) then
        do lcase=1,nl(icase) !----------- loop over L(jatom) requested in the ionput ---------------!
           iind=al_ucase(icase,lcase)
           if (abs(dri(iind)).gt.1e-10) nonzero_interst=.true.
        enddo
        if (nonzero_interst) then
           allocate( h_interstitial(2*lmaxp+1,iblock) )
           if (Qforce_j) then
              allocate( al_interstitial(2*lmaxp+1,nbands,nl(icase),4) )
           else
              allocate( al_interstitial(2*lmaxp+1,nbands,nl(icase),1) )
           endif
        endif
     endif

     rotloc_x_BR1 = matmul(rotloc(:,:,jatom), BR1)
     rotloc_x_BR1_x_rotij = matmul( rotloc_x_BR1, rotij(:,:,latom) )
     crotloc_x_BR1 = matmul(crotloc(:,:,icase), BR1)
     nonzero_shft = sum(abs(shft(latom,:))) .GT. 1e-10
     DO is=1,iso  !--- over both spins
        alm(:,:) = 0.0         !------  ALM(m,band,nrf,is) will hold product of eigenvectors and a/b expansion coefficients --!
        blm(:,:) = 0.0
        clm(:,:,:) = 0.0
        ALLOCATE( alm0(2*lmaxp+1,nbands) )
        if (abs(projector).ge.5 .and. nonzero_interst) al_interstitial(:,:,:,:)=0
        if (Qforce_j) then
           allocate( h_k(iblock,3) )
           allocate( h_ablyl_hk(2*lmaxp+1,iblock) )
           h_k(:,:) = 0.d0
           aalm(:,:,:)=0.d0
           bblm(:,:,:)=0.d0
           cclm(:,:,:,:)=0.d0
        endif
        !--------- blocks are for efficiency. Matrix is multiplied in block form. This must be important in the past, while modern BLAS should do that better. I think it is obsolete.
        DO ii=1,n0-nnlo,iblock !------ iblock is 128 for 32-bit system -------!
           !-------- nlo-number of local orbitals -----!
           h_yl=0

           !$OMP PARALLEL DO PRIVATE(bkrot2,bkrot3,bkrloc,yl,arg1,arg2,argt2,phshel,lcase,l,m,index)&
           !$OMP& SHARED(h_yl,h_k)&
           !$OMP& SCHEDULE(STATIC)
           do i=ii,min(ii+iblock-1,n0-nnlo)  
              !---------  rotates ylm(k+K) to ylm(k'+K) where k' is in star of irreducible k. ------------!
              ! bk3(:,i) !-----  reciprocal vector and irreducible vector: G=K+k ----!
              ! BKROT2 = R_n.R_a.(k+K), transformation from the first atom to an equivalent atom 
              BKROT2 = matmul(rotij(:,:,latom), bk3(:,I))  ! apply one of the symmetry operations to BKROT=R.(k+K)
              ! !!---- BR1 transforms integer reciprocal lattice vectors, as given in the VECTORLIST of LAPW1, into cartesian system ----!
              ! !! BKROT3 = R_n.R_a.(k+K), but in cartesian coordinate system
              ! !BKROT3 = matmul(BR1, BKROT2)
              ! !!---- BKRLOC = Rotloc.R_g.(k+K),  is rotates Rotloc.R.(k+K) . Rotation Rotloc entered by user
              ! !! Here we use crotloc and not rotloc, hence we use the transformation case.indmfl file, input by the user.
              ! !BKRLOC = matmul(crotloc(:,:,icase),BKROT3)
              BKRLOC = matmul(crotloc_x_BR1, BKROT2)
              !---- YLM = Y_{L}(Rotloc.R_g.(k+K))
              CALL YLM (BKRLOC,lmaxp,YL,lmaxplmaxp)  ! 
              ! ARG1 + ARG2 + ARG3 = (R_g.(k+K)) *  R(lfirst) * 2pi
              ARG1 = dot_product(BKROT2, POS(:,lfirst))*TWOPI 
              ! ARGT2 = (R_a.(k+K)).tau_n * 2pi
              ARGT2= dot_product(BK3(:,I), tauij(:3,latom))*twopi
              ARG2 = 0.d0
              if (nonzero_shft) ARG2 = dot_product(BK3(:,I),shft(latom,:))*TWOPI ! Before user rotation, but already on
              ! PHSEHL = e^{I*2pi*( (R_g.(k+K)) *  R(iatom) + (K+k)*tau(isym))}
              PHSHEL=EXP(IMAG*(ARG1+ARG2+ARGT2))
              do lcase=1,nl(icase) !----------- loop over L(jatom) requested in the ionput ---------------!
                 l=ll(icase,lcase) !------ current L --!
                 h_yl(l*l+1:(l+1)**2, i-ii+1) = dconjg(yl(l*l+1:(l+1)**2))*phshel
              enddo
              IF(Qforce_j) THEN
                 h_k(i-ii+1,:) = matmul(rotloc_x_BR1_x_rotij, k3(:,i))
              ENDIF
           enddo
           !$OMP END PARALLEL DO          

           rmt2=1.D0/(rmt(jatom)**2)
           do lcase=1,nl(icase) !----------- loop over L(jatom) requested in the ionput ---------------!
              l=ll(icase,lcase) !------ current L --!
              h_alyl=0
              h_blyl=0
              !$OMP PARALLEL DO PRIVATE(i3) SHARED(h_al,h_bl) SCHEDULE(STATIC)
              do i=ii,min(ii+iblock-1,n0-nnlo)
                 i3=i-ii+1
                 if (lapw(l)) then
                    h_al(i3)=dfj(l,i,jatom)*PE(l)-fj(l,i,jatom)*DPE(l) 
                    h_bl(i3)=fj(l,i,jatom)*DP(l)-dfj(l,i,jatom)*P(l)
                 else
                    h_al(i3) = rmt2*fj(l,i,jatom)/P(l)
                    h_bl(i3) = 0.0
                 end if
              enddo
              !$OMP END PARALLEL DO
              
              !$OMP PARALLEL DO PRIVATE(i3,m,index) SHARED(h_alyl,h_blyl) SCHEDULE(STATIC)
              do i=ii,min(ii+iblock-1,n0-nnlo)
                 i3 = i-ii+1
                 DO m=1,2*l+1
                    index = l*l+m
                    h_alyl(m,i3)=h_al(i3)*h_yl(index,i3)
                    h_blyl(m,i3)=h_bl(i3)*h_yl(index,i3)
                 enddo
              ENDDO
              !$OMP END PARALLEL DO
              
              if (abs(projector).ge.5 .and. nonzero_interst) then
                 h_interstitial(:,:)=0
                 do i=ii,min(ii+iblock-1,n0-nnlo)
                    !----    h_interstitial =Y*_{lm}(R(k+K))*exp(i*(k+K)*r_latom)  <P_phi|j_l>/Rmt^2
                    h_interstitial(1:2*l+1,i-ii+1) = h_yl(l*l+1:(l+1)**2,i-ii+1)*phi_jl(i,iind)*rmt2
                 enddo
              endif
              !!!  End projector 5
              
              !---- h_alyl(2*lmax+1,iblock)  contains rotated Apw's, such that chi(r) = (Apw*u(r) + Bpw*udot(r))*Ylm(r)
              !---- h_blyl(2*lmax+1,iblock)  contains rotated Bpw's
              !---- A(:N,iband,is)              contains eigenvectors, also named C(k+G,inad,is)
              !---- alm[:,iband][irf=1,is] += sum_{iK\in block} h_alyl[:,iK][is]*A[iK,iband][is]
              !---- alm[:,iband][irf=2,is] += sum_{iK\in block} h_blyl[:,iK][is]*A[iK,iband][is]
              !---- 
              !---- The results is:
              !---- alm[lm,iband][1,is] = sum_G Apw(lm,is,K+G) * C(k+G,iband,is)
              !---- alm[lm,iband][2,is] = sum_G Bpw(lm,is,K+G) * C(k+G,iband,is)
              !---- Where C(k+G,iband,is) are eigenvectors, and Apw and Bpw are expansion coefficients defined in Shick et.al., PRB 60, 10763 (1999).
              ibb=min(iblock,n0-nnlo-ii+1)
              ldb=nmat
              lda=2*lmaxp+1
              ldc=2*lmaxp+1          
              isize_ = (2*l+1)
              CALL zgemm('N','N',isize_,nbands,ibb,cone,h_alyl,lda,As(ii,DM_nemin,is),ldb,czero,alm0,ldc) ! Here was bug, 12/28/2013!
              !$OMP PARALLEL DO PRIVATE(m) SHARED(alm,alm0) SCHEDULE(STATIC)
              DO m=1,2*l+1
                 alm(l*l+m,DM_nemin:DM_nemax) = alm(l*l+m,DM_nemin:DM_nemax) + alm0(m,:nbands)
              ENDDO
              !$OMP END PARALLEL DO
              CALL zgemm('N','N',isize_,nbands,ibb,cone,h_blyl,lda,As(ii,DM_nemin,is),ldb,czero,alm0,ldc) ! Here was bug, 12/28/2013!
              !$OMP PARALLEL DO PRIVATE(m) SHARED(blm,alm0) SCHEDULE(STATIC)
              DO m=1,2*l+1
                 blm(l*l+m,DM_nemin:DM_nemax)=blm(l*l+m,DM_nemin:DM_nemax)+alm0(m,:nbands) ! Here was bug, 12/28/2013!
              ENDDO
              !$OMP END PARALLEL DO
              !! The new
              if (abs(projector).ge.5 .and. nonzero_interst) then 
                 call zgemm('N','N',isize_,nbands,ibb,cone,h_interstitial,lda,As(ii,DM_nemin,is),ldb,cone, al_interstitial(:,:,lcase,1),ldc)
              endif
              !! The new
              
              IF (Qforce_j) THEN 
                 do i_h_k=1,3
                    !$OMP PARALLEL DO PRIVATE(m) SHARED(h_ablyl_hk,h_alyl,h_k) SCHEDULE(STATIC)
                    DO m=1,2*l+1
                       h_ablyl_hk(m,:ibb) = h_alyl(m,:ibb)*h_k(:ibb,i_h_k)    ! h_ablyl <-  alm(lm,K)*K
                    enddo
                    !$OMP END PARALLEL DO
                    
                    CALL zgemm('N','N',isize_,nbands,ibb,cone,h_ablyl_hk,lda,As(ii,DM_nemin,is),ldb,czero,alm0,ldc)
                    !$OMP PARALLEL DO PRIVATE(m) SHARED(aalm,alm0) SCHEDULE(STATIC)
                    DO m=1,2*l+1
                       aalm(l*l+m,:nbands,i_h_k) = aalm(l*l+m,:nbands,i_h_k) + alm0(m,:nbands)
                    ENDDO
                    !$OMP END PARALLEL DO
                    
                    !$OMP PARALLEL DO PRIVATE(m) SHARED(h_ablyl_hk,h_blyl,h_k) SCHEDULE(STATIC)
                    DO m=1,2*l+1
                       h_ablyl_hk(m,:ibb) = h_blyl(m,:ibb)*h_k(:ibb,i_h_k)     ! h_ablyl <-  blm(lm,K)*K
                    enddo
                    !$OMP END PARALLEL DO

                    CALL zgemm('N','N',isize_,nbands,ibb,cone,h_ablyl_hk,lda,As(ii,DM_nemin,is),ldb,czero,alm0,ldc)
                    !$OMP PARALLEL DO PRIVATE(m) SHARED(bblm,alm0) SCHEDULE(STATIC)
                    DO m=1,2*l+1
                       bblm(l*l+m,:nbands,i_h_k) = bblm(l*l+m,:nbands,i_h_k) + alm0(m,:nbands)
                    ENDDO
                    !$OMP END PARALLEL DO

                    !! The new
                    if (abs(projector).ge.5 .and. nonzero_interst) then
                       !$OMP PARALLEL DO PRIVATE(m) SHARED(h_ablyl_hk,h_interstitial,h_k) SCHEDULE(STATIC)
                       DO m=1,2*l+1
                          h_ablyl_hk(m,:ibb) = h_interstitial(m,:ibb)*h_k(:ibb,i_h_k)    ! h_ablyl <-  alm(lm,K)*K
                       enddo
                       !$OMP END PARALLEL DO
                       call zgemm('N','N',isize_,nbands,ibb,cone,h_ablyl_hk,lda,As(ii,DM_nemin,is),ldb,cone, al_interstitial(:,:,lcase,i_h_k+1),ldc)
                    endif
                    !! The new
                 enddo
              ENDIF
           enddo !----------- over lcase --------------!
        ENDDO    !----------- over iblock -------------!
        
        DEALLOCATE( alm0 )
        if (Qforce_j) then
           deallocate( h_k )
           deallocate( h_ablyl_hk )
        endif
        !-------------- Adds localized orbitals to alm. -------------------!
        if (nlo.ne.0) call lomain(crotloc(:,:,icase),is,DM_nemin,DM_nemax,lfirst,latom,n0,jatom,alm,blm,clm,Qforce_j,aalm,bblm,cclm,lmaxp)

        !-------------- Multiplying with a constant all coefficients ------!
        FAC=4.0D0*PI*RMT(jatom)**2/SQRT(VOL)
        do lcase=1,nl(icase) !----------- loop over L(jatom) requested in the ionput ---------------!
           l=ll(icase,lcase) !------ current L --!
           cfac=fac*imag**l
           alm(l*l+1:(l+1)**2,DM_nemin:DM_nemax)=alm(l*l+1:(l+1)**2,DM_nemin:DM_nemax)*cfac
           blm(l*l+1:(l+1)**2,DM_nemin:DM_nemax)=blm(l*l+1:(l+1)**2,DM_nemin:DM_nemax)*cfac
           if (l.le.lomax) clm(l*l+1:(l+1)**2,DM_nemin:DM_nemax,1:ilo(l))=clm(l*l+1:(l+1)**2,DM_nemin:DM_nemax,1:ilo(l))*cfac
           if (abs(projector).ge.5 .and. nonzero_interst) al_interstitial(:,:,lcase,1)=al_interstitial(:,:,lcase,1)*cfac
        enddo
        IF (Qforce_j) THEN
           do lcase=1,nl(icase) !----------- loop over L(jatom) requested in the ionput ---------------!
              l=ll(icase,lcase) !------ current L --!
              cfac=fac*imag**l
              aalm(l*l+1:(l+1)**2,1:nbands,1:3)          = aalm(l*l+1:(l+1)**2,1:nbands,1:3)*cfac
              bblm(l*l+1:(l+1)**2,1:nbands,1:3)          = bblm(l*l+1:(l+1)**2,1:nbands,1:3)*cfac
              if (l.le.lomax) cclm(l*l+1:(l+1)**2,1:nbands,1:ilo(l),1:3) = cclm(l*l+1:(l+1)**2,1:nbands,1:ilo(l),1:3)*cfac
              if (abs(projector).ge.5 .and. nonzero_interst) al_interstitial(:,:,lcase,2:4)=al_interstitial(:,:,lcase,2:4)*cfac
           ENDDO
        END IF
        !---------- Computing the DMFT transformation --------------!
        do lcase=1,nl(icase)
           ! connection between ri_mat and quantities in lapw2
           !
           !       ri_mat(1,1,l)=1.0    ! <u|u>
           !       ri_mat(1,2,l)=0.0    ! <udot|u>
           !       ri_mat(2,1,l)=0.0    ! <u|udot>
           !       ri_mat(2,2,l)=pei(l) ! <udot|udot>
           !
           !       jlo=1,ilo(l)
           !          ri_mat(1,2+jlo,l) = pi12lo(jlo,l)  ! <u | u_lo>
           !          ri_mat(2+jlo,1,l) = pi12lo(jlo,l)  ! <u_lo | u>
           !          ri_mat(2,2+jlo,l) = pe12lo(jlo,l)  ! <udot | u_lo>
           !          ri_mat(2+jlo,2,l) = pe12lo(jlo,l)  ! <u_lo | udot>
           !
           !       jlo=1,ilo(l)
           !          jlop=1,ilo(l)
           !             ri_mat(2+jlo,2+jlop,l) = pr12lo(jlo,jlop,l)  ! <u_lo | u_lo >
           !
           l=ll(icase,lcase)
           nind=(2*l+1)*iso
           icix = cix(icase,lcase)
           iorb = iorbital(icase,lcase)
           do num1=DM_nemin,DM_nemax
              nm1 = num1-DM_nemin+1
              do m1=-l,l           ! over m-1
                 index=l*(l+1)+m1+1
                 ind = l+1+m1+(2*l+1)*(is-1)
                 dff = 1.0
                 if (abs(projector).eq.2 .or. abs(projector).eq.4) then
                    sm1 = abs(alm(index,num1))**2 * ri_mat(1,1,l,jatom) + abs(blm(index,num1))**2 * ri_mat(2,2,l,jatom)  !! 1,1  + 2,2
                    DO jlo=1,ilo(l)
                       IF (loor(jlo,l)) THEN
                          sm1 = sm1 + 2*dble(alm(index,num1)*conjg(clm(index,num1,jlo))) * ri_mat(1,2+jlo,l,jatom) !! 1,3 and 3,1
                          sm1 = sm1 + 2*dble(blm(index,num1)*conjg(clm(index,num1,jlo))) * ri_mat(2,2+jlo,l,jatom) !! 2,3 and 3,2
                          DO jlop=1,ilo(l)
                             IF (loor(jlop,l)) THEN
                                sm1 = sm1 + dble(clm(index,num1,jlo)*conjg(clm(index,num1,jlop))) * ri_mat(2+jlo,2+jlop,l,jatom) !! 3,3
                             ENDIF
                          ENDDO
                       ENDIF
                    ENDDO
                    sm2c = alm(index,num1)
                    DO jlo=1,ilo(l)
                       if (loor(jlo,l)) sm2c = sm2c + clm(index,num1,jlo)*ri_mat(1,2+jlo,l,jatom)
                    ENDDO
                    sm2 = abs(sm2c)**2
                    if (abs(sm2)>0.0) dff = sqrt(abs(sm1/sm2))
                 endif
                 if (abs(projector).le.4) then
                    ! Common part for projection-1, Projection-2, Projection-3 and Projection-4
                    cc = alm(index,num1)                                           !! 1,1
                    DO jlo=1,ilo(l)
                       IF(loor(jlo,l)) THEN
                          cc = cc + clm(index,num1,jlo)*ri_mat(1,2+jlo,l,jatom)    !!  1,3
                       ENDIF
                    ENDDO
                    URx(nm1,ind,iorb,1) = cc * dff

                    if (Qforce_j) then
                       dd(1:3) = aalm(index,nm1,1:3)                                                !!  1,1
                       DO jlo=1,ilo(l)
                          IF(loor(jlo,l)) THEN
                             dd(1:3) = dd(1:3) + cclm(index,nm1,jlo,1:3)*ri_mat(1,2+jlo,l,jatom)    !!  1,3
                          ENDIF
                       ENDDO
                       URx(nm1,ind,iorb,2:4) = dd(1:3) * dff
                    endif
                    
                 else if (abs(projector).eq.5 .or.  abs(projector).eq.6) then
                    iind=al_ucase(icase,lcase)
                    cc = alm(index,num1)*rix_mat(1,iind) + blm(index,num1)*rix_mat(2,iind)  !! 1,1 + 1,2
                    DO jlo=1,ilo(l)
                       IF(loor(jlo,l)) THEN
                          cc = cc + clm(index,num1,jlo)*rix_mat(2+jlo,iind)                 !! 1,3
                       ENDIF
                    ENDDO
                    if (nonzero_interst) cc = cc + al_interstitial(l+m1+1, nm1, lcase, 1)
                    URx(nm1,ind,iorb,1) = cc
                    
                    if (Qforce_j) then
                       dd(1:3) = aalm(index,nm1,1:3)*rix_mat(1,iind) + bblm(index,nm1,1:3)*rix_mat(2,iind)  !! 1,1 + 1,2
                       DO jlo=1,ilo(l)
                          IF(loor(jlo,l)) THEN
                             dd(1:3) = dd(1:3) + cclm(index,nm1,jlo,1:3)*rix_mat(2+jlo,iind)                 !! 1,3
                          ENDIF
                       ENDDO
                       if (nonzero_interst) dd(1:3) = dd(1:3) + al_interstitial(l+m1+1, nm1, lcase, 2:4)
                       URx(nm1,ind,iorb,2:4) = dd(1:3)
                    endif
                    
                 else
                    print *, 'Only projector=[1,2,3,4,5,6,-1,-2,-3,-4,-5,-6] is allowed!'
                    stop
                 endif
              enddo
           enddo
        enddo
     enddo !----------- over both spins ---------!
     if (abs(projector).ge.5 .and. nonzero_interst) then
        deallocate( h_interstitial, al_interstitial )
     endif
  ENDDO    !----------- over all atoms  ---------!


  deallocate( alm, blm, clm )
  deallocate( aalm, bblm, cclm ) 
  

  allocate( tmp(DM_nemax-DM_nemin+1,maxdim2) )

  !nbands = DM_nemax-DM_nemin+1
  !allocate( DMFTU(nbands,maxdim2,norbitals) )
  
  DMFTU = 0
  DO iorb1=1,norbitals
     icix = cix_orb(iorb1)
     if (icix.EQ.0) CYCLE
     nind1 = nindo(iorb1)

     ! Spin rotation
     if (iso.eq.2) then
        no1 = nind1/iso
        ALLOCATE( Uu(DM_nemax-DM_nemin+1,no1), Ud(DM_nemax-DM_nemin+1,no1))
        tuu = Rspin(1,1,iorb1)
        tud = Rspin(1,2,iorb1)
        tdu = Rspin(2,1,iorb1)
        tdd = Rspin(2,2,iorb1)
        do ip=1,nipc
           Uu(:,:) = URx(:,1:no1,iorb1,ip)
           Ud(:,:) = URx(:,no1+1:2*no1,iorb1,ip)
           URx(:,1:no1,iorb1,ip)       = Uu*tuu + Ud*tdu
           URx(:,no1+1:2*no1,iorb1,ip) = Uu*tud + Ud*tdd
        enddo
        DEALLOCATE( Uu, Ud )
     endif

     ! Rotation to local coordinate system specified by the user in case.indmfl
     DO iorb2=1,norbitals
        if ( cix_orb(iorb2).NE. icix ) CYCLE
        nind2 = nindo(iorb2)
        do ip=1,nipc
           call zgemm('N','C', nbands, nind2, nind1, cone, URx(:,:,iorb1,ip), nbands, cfX(:,:,iorb2,iorb1),maxdim2,czero, tmp, nbands)
           DMFTU(:,:nind2,iorb2,ip) = DMFTU(:,:nind2,iorb2,ip) + conjg(tmp(:,:nind2))
        enddo
     ENDDO
  ENDDO


  
  DEALLOCATE( URx, tmp )

END SUBROUTINE Build_DMFT_Projector


SUBROUTINE Diagonalize_DMFT_WEIGHTS(zw2, Aweight, nbands, DM_nemin, DM_nemaxx, determine_DM_nemaxx)
  IMPLICIT NONE
  REAL*8,     intent(out)  :: zw2(nbands)
  COMPLEX*16, intent(inout):: Aweight(nbands,nbands)
  INTEGER,    intent(in)   :: nbands, DM_nemin
  INTEGER,    intent(inout):: DM_nemaxx
  LOGICAL,    intent(in)   :: determine_DM_nemaxx
  ! locals
  INTEGER    :: i, j, n, nbandsx, num
  INTEGER, allocatable    :: imx(:)
  COMPLEX*16, allocatable :: zweight(:,:), znew(:), Anew(:,:)
  COMPLEX*16 :: cc
  INTEGER    ::  lwork, lrwork, liwork
  complex*16, allocatable :: work(:)
  real*8,     allocatable :: rwork(:)
  integer,    allocatable :: iwork(:)
  integer :: info
  logical :: debug, negative
  !
  lwork = 2*nbands+nbands*nbands
  lrwork =  5*nbands + 2*nbands*nbands + 1
  liwork = 3+5*nbands
  ALLOCATE ( work(lwork), rwork(lrwork), iwork(liwork) )

  debug = .false.
  
  if (debug) then
     ! Remember original Aweight
     allocate( zweight(nbands,nbands) )
     zweight(:,:) = Aweight(:,:)
  endif
  
  ! Aweight = Aw . zw2 . Aw^+
  Aweight(:,:) = -Aweight(:,:)

  CALL ZHEEVD('V','U', nbands, Aweight, nbands, zw2, work, lwork, rwork, lrwork, iwork, liwork, info )
  if (info .ne. 0) then
     print *, 'Diagonalization of weights failed. Info-zheevd=', info
  endif
  ! We defines zw2 with minus before to sort eigenvalues from largest to smallest
  ! To make them positive, change sign here
  zw2 = -zw2

  allocate( imx(nbands) )

  negative=.false.
  do num=1,nbands
     if (zw2(num).lt.0) negative=.true.
  enddo
  if (negative) then
     CALL eig_order_abs_val(zw2, imx, nbands)
     CALL permute_eigensystem1(imx, zw2, Aweight, nbands)
  endif
  
  if (determine_DM_nemaxx) then
     ! How many bands have finite weigh and thus need to be considered?
     ! Fron now on, we will only use bands between DM_nemin and DM_nemaxx, rather than DM_nemax.
     nbandsx=nbands
     do nbandsx=nbands,1,-1
        !print *, 'nbandsx=', nbandsx, zw2(nbandsx)
        if (abs(zw2(nbandsx)).gt.1e-9) exit
     enddo
     DM_nemaxx = DM_nemin + nbandsx - 1
  else
     nbandsx = DM_nemaxx - DM_nemin + 1
     do num=nbandsx+1,nbands
        if (abs(zw2(num)).gt.1e-6) nbandsx=num
     enddo
     DM_nemaxx = DM_nemin + nbandsx - 1
  endif
  !print *, 'DM_nemaxx=', DM_nemaxx, 'nbandsx=', nbandsx
  
  ! It turns out that most of occupied bands have almost equal occupancy zw2
  ! and hence large degeneracy tends to reorder eigenvectors arbitrary. While
  ! this is not crucial for the final result, it is useful to reorder eigenvectors
  ! so that the transformation is as close to identiy as possible.
  ! 
  ! Hence we take the eigenvectors with finite weights and find index of
  ! transformation imx(:) so that Aweight is close to identity.
  call sort_to_diagonal(Aweight(:nbandsx,:nbandsx), imx(:nbandsx), nbandsx)

  allocate( znew(nbandsx), Anew(nbands,nbandsx) )
  
  ! Now permuting eigenvalues and eigenvectors.
  do i=1,nbandsx
     Anew(:,imx(i)) = Aweight(:,i)
     znew(imx(i)) = zw2(i)
  enddo
  Aweight(:,:nbandsx) = Anew(:,:nbandsx)
  zw2(:nbandsx) = znew(:nbandsx)
  
  

  if (debug) then

     WRITE(6,*) 'Aweight after'
     do i=1,nbands
        do j=1,nbands
           WRITE(6,'(2f10.5,2x)', advance='no') Aweight(i,j)
        enddo
        WRITE(6,*)
     enddo
     
     ! checking exact diagonalization
     do i=1,nbands
        do j=1,nbands
           cc=0
           do n=1,nbands
              cc = cc + Aweight(i,n)*zw2(n)*conjg(Aweight(j,n))
           enddo
           if (abs(zweight(i,j)-cc).GT.1e-8) print *, 'd',i,j,abs(zweight(i,j)-cc)
        enddo
     enddo
  endif
  DEALLOCATE ( work, rwork, iwork )
  
  deallocate( imx )
  deallocate( znew, Anew )
  if (debug) deallocate( zweight)
END SUBROUTINE Diagonalize_DMFT_WEIGHTS


!---------------------------------------------
! To compute overlap needed by renormalization
!---------------------------------------------
SUBROUTINE cmp_overlap(Olapm, SOlapm, max_nbands, nnlo, norbitals, natom, maxdim2, emin, emax, iorbital, cix_orb, iSx, nindo, noccur, cixdim, cfX, Rspin, sumw, lmaxp, time_dmf0, time_dmf0w, SIMPLE)
  USE xa,    ONLY: fj, dfj
  USE xa3,   ONLY: bk3, bk3lo, aK, k3, k3lo, As, As_lo
  USE param, ONLY: nume, nmat, nloat, iblock, lmax2
  USE structure, ONLY: rmt, BR1
  USE com,   ONLY: nat
  USE defs,  ONLY: PI, IMAG
  USE dmfts,   ONLY: DM_Emin, DM_Emax, maxsize, iso, ncix, Sigind, csize, maxdim, projector!, lmaxp, nl, cix, DM_EF, ll, 
  USE com_mpi,ONLY: myrank, master, FilenameMPI, Reduce_MPI, FindMax_MPI, reduce_MPI0,pr_proc,vector_para,vectors,nvector,fvectors,Qprint
  USE p_project,ONLY: phi_jl, rix, P_rfi, kNri, n_al_ucase, dri, l_al_ucase, j_al_ucase
  !USE sym2,  ONLY: iord
  !USE MProjector, ONLY: Build_DMFT_Projector2
  IMPLICIT NONE
  COMPLEX*16, intent(out):: Olapm(maxdim,maxdim,ncix), SOlapm(maxdim,maxdim,ncix)
  INTEGER,    intent(out):: max_nbands
  REAL*8,     intent(out):: time_dmf0, time_dmf0w
  INTEGER,    intent(in) :: nnlo, norbitals, natom, maxdim2, lmaxp
  REAL*8,     intent(in) :: emin, emax, sumw
  INTEGER,    intent(in) :: iorbital(natom,lmaxp+1), cix_orb(norbitals), iSx(maxdim2, norbitals), nindo(norbitals), noccur(maxsize,ncix), cixdim(ncix)
  complex*16, intent(in) :: cfX(maxdim2,maxdim2,norbitals,norbitals), Rspin(2,2,norbitals)
  LOGICAL, intent(in)    :: SIMPLE
  ! interfaces
  interface
     REAL*8 Function Romb(y,N,dx)
       IMPLICIT NONE
       REAL*8, intent(in) :: y(N)
       REAL*8, intent(in) :: dx
       INTEGER, intent(in):: N
     end Function Romb
  end interface
  ! locals
  complex*16, ALLOCATABLE :: DMFTU(:,:,:,:), Olapmk(:,:,:), olp(:,:)
  INTEGER    :: ikp, nemin, nemax, DM_nemin, DM_nemax, n0
  REAL*8     :: E(nume), wgh0, emist, Olapc(maxsize,ncix), Kn(3), rx
  LOGICAL    :: more_kpoints, Tcompute
  INTEGER    :: jatom, isize, nbands, itape, is, i, iorb1, iorb2, icix, nind1, nind2, ind1, ind2, cixdm, ip, iq, it, iikp, iks, ivector, nkp
  REAL*8     :: t1c, t2c, t1w, t2w, wg
  !
  COMPLEX*16, allocatable :: olocf(:,:), work(:),tmp1(:,:), solap(:,:)
  REAL*8, allocatable     :: ws(:), rwork(:)
  INTEGER, allocatable    :: iwork(:)
  INTEGER :: lwork, lrwork, liwork, info, iind, ir, l, Nri, nipc!, isym
  ! for p_interstitial
  REAL*8, allocatable     :: aKR(:), jlr(:), jlp(:)
  !
  INTEGER, allocatable :: cind(:), cini(:)
  INTEGER :: cixdms
  complex*16 :: czero, cone

  cone  = cmplx(1.d0, 0.d0, 8)
  czero = cmplx(0.d0, 0.d0, 8)
  allocate( Olapmk(maxdim,maxdim,ncix), olp(maxdim2,maxdim2) )
  Olapm=0

  call cputim(t1c)
  call walltim(t1w)


  max_nbands=0
  iikp=0
  DO ivector=1,nvector
     DO is=1,iso    !------ over up/dn ---------------------!
        itape=8+is
        if (vector_para) then
           open(itape,FILE=fvectors(ivector,is),STATUS='old',FORM='unformatted')
        else
           rewind(itape) !--- both vector files: 9 and 10 rewind -------------------------!
        endif
        DO I=1,NAT
           READ(itape) EMIST !---- At the beginninge we have linearization energies --!
           READ(itape) EMIST !---- We just skip them ---------------------------------!
        ENDDO
     END DO
     nkp = vectors(ivector,2)

     DO iks=1,nkp ! kpoint loop begin
        if (vector_para) then
           ikp = vectors(ivector,3)+iks  ! successive index in k-point table from case.klist
           iikp = iikp+1                 ! successive index in k-point table on this processor
        else
           ikp = iks                     ! successive index in k-point table from case.klist
           !--- We need to go over all k-points even though we will compute only some of them on this processor.
           !--- This is because we need to read vector file sequentially.
           iikp = ikp-myrank*pr_proc            ! The index of the point to compute. If negative, do not compute! 
        endif

        Tcompute=.FALSE.                                                                                                                                                                        
        if (iikp.gt.0) Tcompute=.TRUE.       ! If Tcompute is true, the point needs to be computed.                                                                                             
        if (iikp.GT.pr_proc) EXIT            ! Processor finished. The rest of the points will be taken care of by the other processors.                                                        

        CALL Read_Vec_Spin(E, wgh0, As, As_lo, k3, k3lo, bk3, bk3lo, nemin, nemax, DM_nemin, DM_nemax, more_kpoints, n0, emin, emax, DM_Emin, DM_Emax, iso, nmat, nume, nnlo)

        nbands = DM_nemax-DM_nemin+1
        max_nbands = max(nbands, max_nbands)
        IF (.NOT.Tcompute) CYCLE                                                                                                                                                                
        IF (.NOT.more_kpoints) EXIT                                                                                                                                                              

        isize=n0-nnlo
        DO I=1,isize                                !----- Over reciprocal vectors -----------------------!
           Kn(:) = matmul(BR1,BK3(:,I))             !----  Creates K+k in cartesian coordinates ----------!
           aK(I) = sqrt(Kn(1)**2+Kn(2)**2+Kn(3)**2) !----  calculates |K+k| for use in bessel functions --!
        ENDDO
        do jatom=1,nat
           CALL HARMON(isize,aK(:isize),lmax2,fj(:,:isize,jatom),dfj(:,:isize,jatom),rmt(jatom))
        enddo
        ! New
        if (abs(projector).EQ.5) then
           allocate( phi_jl(nmat, n_al_ucase) )
           
           Nri=2**kNri+1    ! Number of radial points in the interstitials
           phi_jl(:,:)=0
           allocate( aKR(Nri), jlr(Nri), jlp(Nri) )
           DO iind=1,n_al_ucase
              if (abs(dri(iind)).gt.1e-10) then
                 l     = l_al_ucase(iind)
                 jatom = j_al_ucase(iind)
                 DO i=1,isize         ! over all reciprocal vectors K
                    rx = RMT(jatom)
                    do ir=1,Nri
                       aKR(ir)=rx*aK(i)    ! |k+K|*r
                       rix(ir)=rx          !  r
                       rx = rx + dri(iind)
                    enddo
                    CALL sphbes2(l,Nri,aKR,jlr)  ! spherical bessel : j_l(|k+K|*r)
                    jlr(:) = jlr(:)*rix(:)       ! r*j_l . We need to do that, because the integral is Int[ (phi(r)/r)*j_l(r)*r^2]=Int[phi(r)*j_l(r)*r]
                    jlp(:) = jlr(:)*P_rfi(:,iind)
                    phi_jl(i,iind) = romb(jlp, Nri, dri(iind))  ! Integration over r on rix mesh
                 ENDDO
              endif
           ENDDO
           deallocate( aKR, jlr, jlp )
        endif
        ! New

        nipc=1
        allocate( DMFTU(nbands,maxdim2,norbitals,nipc) )

        CALL Build_DMFT_Projector(DMFTU, cfX, Rspin, iorbital, norbitals, nipc, n0, nnlo, nbands, cix_orb, nindo, DM_nemin, DM_nemax, maxdim2, lmaxp)
        
        Olapmk=0
        DO iorb1=1,norbitals
           icix = cix_orb(iorb1)
           if ( icix.EQ.0 ) CYCLE
           nind1 = nindo(iorb1)
           DO iorb2=1,norbitals
              if ( icix.NE.cix_orb(iorb2) ) CYCLE
              nind2 = nindo(iorb2)
              olp=0
              call zgemm('C','N', nind1, nind2, nbands, cone, DMFTU(:,:,iorb1,1),nbands, DMFTU(:,:,iorb2,1),nbands, czero, olp,maxdim2)
              do ind1=1,nind1
                 do ind2=1,nind2
                    Olapmk( iSx(ind1,iorb1), iSx(ind2,iorb2), icix ) = olp(ind1,ind2)
                 enddo
              enddo
           enddo
        ENDDO
           
        wg = (wgh0/sumw)
        Olapm(:,:,:) = Olapm(:,:,:) + Olapmk(:,:,:) * wg ! proper weight for the reducible k-point
        deallocate( DMFTU )

        if (abs(projector).EQ.5) deallocate( phi_jl )
        WRITE(*,'(I3,A,1x,I3,1x,I3,1x,A,I4)') myrank, ') Finished k-point number', ikp, iikp, 'with #bands=', nbands

     ENDDO         ! kpoint loop ends
     DO is=1,iso    !------ over up/dn ---------------------!
        itape=8+is
        if (vector_para) then
           close(itape)
        else
           rewind(itape)
        endif
     ENDDO
  ENDDO
  
  CALL Reduce_MPI0(Olapm, maxdim, ncix, max_nbands)
  
  deallocate( olapmk, olp )
  
  if (Qprint) WRITE(*,*) 'Renormalizing Gloc to account for the interstitials'

  if (SIMPLE) then
     Olapc=0
     DO icix=1,ncix
        cixdm = cixdim(icix)
        !---- Olap to vector form, and s_oo to matrix form
        DO ip=1,cixdm
           do iq=1,cixdm
              it = Sigind(ip,iq,icix)
              if (it.gt.0) then
                 Olapc(it, icix) =  Olapc(it, icix) + dble(Olapm(ip,iq,icix))
              endif
           enddo
        ENDDO
        DO it=1,csize(icix)
           Olapc(it, icix) = Olapc(it, icix)/noccur(it,icix)
        ENDDO
     ENDDO
     if (Qprint) WRITE(*,*) 'Z due to missing-interstitials=', Olapc
  else if (.True.) then
     !!! We copute SOlapm = 1/sqrt(Olapm) for later use
     lwork  = 4*maxdim + maxdim*maxdim
     lrwork = 10*maxdim + 2*maxdim*maxdim + 1
     liwork = 3+10*maxdim
     allocate( work(lwork), rwork(lrwork), iwork(liwork) )
     
     SOlapm(:,:,:) = 0.d0
     DO icix=1,ncix
        cixdm = cixdim(icix)

        allocate( cind(cixdm) )
        cixdms=0 ! real dimension, excluding states which must be projected out, because they are treated as non-correlated                                                                                    
        cind=0   ! In most cases, cind(i)=i, however, when some states are projected out, cind points to smaller block                                                                                         
        DO ip=1,cixdm
           if (Sigind(ip,ip,icix) .ne. 0) then
              cixdms = cixdms + 1
              cind(ip) = cixdms
           endif
        ENDDO
        allocate( cini(cixdms))
        do ip=1,cixdm
           if (cind(ip).gt.0) cini(cind(ip))=ip
        enddo
        deallocate( cind )
        
        !allocate( olocf(cixdm,cixdm), ws(cixdm), tmp1(cixdm,cixdm) )
        !olocf(:,:) = Olapm(:cixdm,:cixdm,icix)
        !CALL ZHEEVD('V','U', cixdm, olocf, cixdm, ws, work, lwork, rwork, lrwork, iwork, liwork, info )
        !if (info .ne. 0) then
        !   print *, 'Diagonalization of renormalization factor failed. Info-zheevd=', info
        !endif
        !do ip=1,cixdm
        !   tmp1(:,ip) = olocf(:,ip)*(1./sqrt(abs(ws(ip))))
        !enddo
        !call zgemm('N','C', cixdm, cixdm, cixdm, cone, tmp1, cixdm, olocf, cixdm, czero, SOlapm(:,:,icix), maxdim)

        allocate( olocf(cixdms,cixdms), ws(cixdms), tmp1(cixdms,cixdms), solap(cixdms,cixdms) )
        do ip=1,cixdms
           do iq=1,cixdms
              olocf(ip,iq) = Olapm(cini(ip),cini(iq),icix)
           enddo
        enddo
        CALL ZHEEVD('V','U', cixdms, olocf, cixdms, ws, work, lwork, rwork, lrwork, iwork, liwork, info )
        if (info .ne. 0) then
           print *, 'Diagonalization of renormalization factor failed. Info-zheevd=', info
        endif
        do ip=1,cixdms
           tmp1(:,ip) = olocf(:,ip)*(1./sqrt(abs(ws(ip))))
        enddo
        call zgemm('N','C', cixdms, cixdms, cixdms, cone, tmp1, cixdms, olocf, cixdms, czero, solap, cixdms)


        do ip=1,cixdm
           SOlapm(ip,ip,icix)=1.0d0
        enddo
        do ip=1,cixdms
           do iq=1,cixdms
              SOlapm(cini(ip),cini(iq),icix) = solap(ip,iq)
           enddo
        enddo
        
        if (myrank.eq.master) then
           !WRITE(*,*) 'Z to renormalize=', ws(:cixdm)
           WRITE(*,'(A)',advance='no') 'Z to renormalize='
           do ip=1,cixdms
              WRITE(*,'(F16.10)',advance='no') 1./dble(SOlapm(cini(ip),cini(ip),icix))
           enddo
           WRITE(*,*)
        endif

        deallocate( cini )
        deallocate( olocf, ws, tmp1, solap )
     ENDDO
     deallocate( work, rwork, iwork )
  else
     !!! We copute SOlapm = 1/sqrt(Olapm) for later use
     lwork  = 4*maxdim + maxdim*maxdim
     lrwork = 10*maxdim + 2*maxdim*maxdim + 1
     liwork = 3+10*maxdim
     allocate( work(lwork), rwork(lrwork), iwork(liwork) )
     allocate( olocf(maxdim,maxdim), ws(maxdim), tmp1(maxdim,maxdim) )
     
     DO iorb1=1,norbitals
        icix = cix_orb(iorb1)
        if ( icix.EQ.0 ) CYCLE
        nind1 = nindo(iorb1)
        olocf=0
        do ind1=1,nind1
           do ind2=1,nind1
              olocf(ind1,ind2) = Olapm( iSx(ind1,iorb1), iSx(ind2,iorb1), icix )
           enddo
        enddo

        !if (myrank.eq.master) then
        !   WRITE(*,*) (olocf(ind1,ind1),ind1=1,nind1)
        !endif
        
        CALL ZHEEVD('V','U', nind1, olocf, maxdim, ws, work, lwork, rwork, lrwork, iwork, liwork, info )
        if (info .ne. 0) then
           print *, 'Diagonalization of renormalization factor failed. Info-zheevd=', info
        endif
        
        do ip=1,nind1
           tmp1(:,ip) = olocf(:,ip)*(1./sqrt(abs(ws(ip))))
        enddo
        call zgemm('N','C', nind1, nind1, nind1, cone, tmp1, maxdim, olocf,maxdim, czero, SOlapm(:,:,icix),maxdim)

        if (myrank.eq.master) then
           WRITE(*,*) 'Z due to missing-interstitials=', ws(:nind1)
        endif
        
     ENDDO
     deallocate( work, rwork, iwork )
     deallocate( olocf, ws, tmp1 )
  endif
  call cputim(t2c)
  call walltim(t2w)
  time_dmf0 = t2c-t1c
  time_dmf0w= t2w-t1w
END SUBROUTINE cmp_overlap

SUBROUTINE read_overlap_from_file(info, Olapm, SOlapm, cixdim, maxdim, ncix)
  USE com_mpi,ONLY: myrank, master
  IMPLICIT NONE
  INTEGER,    intent(out) :: info
  COMPLEX*16, intent(out) :: Olapm(maxdim,maxdim,ncix), SOlapm(maxdim,maxdim,ncix)
  INTEGER,    intent(in)  :: cixdim(ncix), maxdim, ncix
  ! locals
  LOGICAL :: there
  INTEGER :: cixdm, ind1, ind2, icix

  INQUIRE( FILE='SOlapm.dat', EXIST=there)
  if (.not.there) then
     info=1
     return
  endif
  
  open(996, FILE='SOlapm.dat', status='old')
  DO icix=1,ncix
     cixdm = cixdim(icix)
     READ(996,*,ERR=991) !icix
     do ind1=1,cixdm
        do ind2=1,cixdm
           READ(996,'(2F20.16)',advance='no',ERR=991)  SOlapm(ind1,ind2,icix)
        enddo
        READ(996,*)
     enddo
  ENDDO
  close(996)
  info=0
  Olapm = SOlapm ! Olapm is currently not needed, since we use orthonormal basis. We thus set overlap to the orhonormal equivalent. Maybe we want in future to read Olapm separately....
  if (myrank.eq.master) then
     WRITE(6,*) 'SOlapm succesfully read from file'
     WRITE(6,*) 'Z will be used:'
     DO icix=1,ncix
        cixdm = cixdim(icix)
        do ind1=1,cixdm
           WRITE(6,'(F16.10)',advance='no')  1./dble(SOlapm(ind1,ind1,icix))
        enddo
        WRITE(6,*)
     ENDDO
     WRITE(6,*)
  endif

  return
991 CONTINUE
  info=2
  close(996)
END SUBROUTINE read_overlap_from_file

!------------------------------------------------------
! To renormalize the DMFT transformation
!------------------------------------------------------
SUBROUTINE RenormalizeTrans(DMFTU, Olapm0, SOlapm, cix_orb, cixdim, nindo, iSx, nbands, nipc, maxdim2, norbitals, maxdim, ncix, SIMPLE)
  IMPLICIT NONE
  COMPLEX*16, intent(inout) :: DMFTU(nbands,maxdim2,norbitals,nipc)
  COMPLEX*16, intent(in)    :: Olapm0(maxdim,maxdim,ncix), SOlapm(maxdim,maxdim,ncix)
  INTEGER, intent(in)       :: cix_orb(norbitals), nindo(norbitals), iSx(maxdim2,norbitals), cixdim(ncix)
  INTEGER, intent(in)       :: nbands, maxdim2, norbitals, maxdim, ncix, nipc
  LOGICAL, intent(in)       :: SIMPLE
  ! locals
  INTEGER :: iorb1, icix, nind1, ind1, ip, cixdm
  REAL*8  :: olocef
  COMPLEX*16, allocatable :: tmp3(:,:), Ucix(:,:), Ucix2(:,:)
  complex*16 :: cone, czero
  cone  = cmplx(1.d0, 0.d0, 8)
  czero = cmplx(0.d0, 0.d0, 8)

  if (SIMPLE) then
     do ip=1,nipc
        DO iorb1=1,norbitals
           icix = cix_orb(iorb1)
           if ( icix.EQ.0 ) CYCLE
           nind1 = nindo(iorb1)
        
           do ind1=1,nind1
              olocef = 1/sqrt(real(Olapm0( iSx(ind1,iorb1), iSx(ind1,iorb1), icix )))
              DMFTU(:,ind1,iorb1,ip) = DMFTU(:,ind1,iorb1,ip) * olocef
           enddo
        ENDDO
     enddo
  else if (.True.) then
     DO icix=1,ncix
        cixdm = cixdim(icix)
        allocate( Ucix(nbands,cixdm), Ucix2(nbands,cixdm) )
        do ip=1,nipc
           Ucix(:,:) = 0.d0
           DO iorb1=1,norbitals
              if ( cix_orb(iorb1).NE.icix ) CYCLE
              nind1 = nindo(iorb1)
              do ind1=1,nind1
                 Ucix(:,iSx(ind1,iorb1)) = DMFTU(:,ind1,iorb1,ip)
              enddo
           ENDDO
           call zgemm('N','N', nbands, cixdm, cixdm, cone, Ucix, nbands, SOlapm(:,:,icix), maxdim, czero, Ucix2, nbands)
           DO iorb1=1,norbitals
              if ( cix_orb(iorb1).NE.icix ) CYCLE
              nind1 = nindo(iorb1)
              do ind1=1,nind1
                 DMFTU(:,ind1,iorb1,ip) = Ucix2(:,iSx(ind1,iorb1))
              enddo
           ENDDO
        enddo
        deallocate( Ucix, Ucix2 )
     ENDDO
  else
     allocate( tmp3(nbands,maxdim2) )
     do ip=1,nipc
        DO iorb1=1,norbitals
           icix = cix_orb(iorb1)
           if ( icix.EQ.0 ) CYCLE
           nind1 = nindo(iorb1)
           call zgemm('N','N', nbands, nind1, nind1, cone, DMFTU(:,:,iorb1,ip),nbands, SOlapm(:,:,icix),maxdim, czero, tmp3, nbands)
           DMFTU(:,:nind1,iorb1,ip) = tmp3(:,:nind1)
        ENDDO
     enddo
     deallocate( tmp3 )
  endif
END SUBROUTINE RenormalizeTrans

SUBROUTINE RenormalizeTransK(DMFTU, cix_orb, cixdim, nindo, iSx, Sigind, projector, nbands, nipc, maxdim2, norbitals, maxdim, ncix)
  IMPLICIT NONE
  COMPLEX*16, intent(inout) :: DMFTU(nbands,maxdim2,norbitals,nipc)
  INTEGER, intent(in)       :: cix_orb(norbitals), nindo(norbitals), iSx(maxdim2,norbitals), cixdim(ncix), Sigind(maxdim,maxdim,ncix)
  INTEGER, intent(in)       :: nbands, maxdim2, norbitals, maxdim, ncix, projector, nipc
  ! locals
  INTEGER :: iorb1, icix, nind1, ind1, cixdm
  REAL*8  :: olocef
  !
  INTEGER :: cixdms, cixdms_m, info, ip, lwork, lrwork, i, ipc
  INTEGER,    allocatable :: cind(:), cini(:), iwork(:)
  COMPLEX*16, allocatable :: Ucix(:,:), Ucix2(:,:), UU(:,:), work(:), Uw(:,:), Vw(:,:), Vw2(:,:)
  REAL*8,     allocatable :: rwork(:), ws(:)
  complex*16 :: cone, czero
  cone  = cmplx(1.d0, 0.d0, 8)
  czero = cmplx(0.d0, 0.d0, 8)
  !
  DO icix=1,ncix
     cixdm = cixdim(icix)
     
     ! If we project out some orbitals, for example eg-orbitals, we might have
     ! Sigind= [[0,0,0,0,0], [0,0,0,0,0], [0,0,1,0,0], [0,0,0,2,0], [0,0,0,0,3]]
     ! then
     !     cind[1:5] = [0,0,1,2,3] and
     !     cini[1:3] = [3,4,5]
     allocate( cind(cixdm) )
     cixdms=0 ! real dimension, excluding states which must be projected out, because they are treated as non-correlated                                                                                    
     cind=0   ! In most cases, cind(i)=i, however, when some states are projected out, cind points to smaller block                                                                                         
     DO ip=1,cixdm
        if (Sigind(ip,ip,icix) .ne. 0) then
           cixdms = cixdms + 1
           cind(ip) = cixdms
        endif
     ENDDO
     
     allocate( cini(cixdms))
     do ip=1,cixdm
        if (cind(ip).gt.0) cini(cind(ip))=ip
     enddo
     ! Now cini[1:3] = [3,4,5] contains the small index of non-zero components
     
     ! If we have cluster-DMFT calculations, we need several orbitals combined into cix block
     allocate( Ucix(nbands,cixdms), Ucix2(nbands,cixdms) )
     Ucix(:,:) = 0.d0
     DO iorb1=1,norbitals
        if ( cix_orb(iorb1).NE.icix ) CYCLE
        nind1 = nindo(iorb1)
        do ind1=1,nind1
           ip = iSx(ind1,iorb1)
           if (cind(ip).gt.0) Ucix(:,cind(ip)) = DMFTU(:,ind1,iorb1,1)
        enddo
     ENDDO

     cixdms_m = min(cixdms,nbands) ! should normally be cixdms, as the number of bands should be larger
     allocate( ws(cixdms_m), Uw(nbands,cixdms_m), Vw(cixdms_m,cixdms) )
     
     lwork = 2*cixdms*(cixdms+nbands)
     lrwork = 7*cixdms*(cixdms + 1)
     allocate( work(lwork), rwork(lrwork), iwork(8*cixdms) )

     !do i=1,nbands
     !   WRITE(6,'(A)', advance='no') "["
     !   do ip=1,cixdms
     !      WRITE(6, '(f14.10,"+",f8.5,3x,"*1j, ")',advance='no') real(Ucix(i,ip)), aimag(Ucix(i,ip))
     !   enddo
     !   WRITE(6,*) "],"
     !enddo

     
     call ZGESDD('S', nbands, cixdms, Ucix, nbands, ws, Uw, nbands, Vw, cixdms_m, work, lwork, rwork, iwork, info )
     if (info .ne. 0) then
        print *, 'SVD decomposition of the projector failed. Info-zgesdd=', info
        if (info.lt.0) print *, 'The ', abs(info),' th argument had an illegal value.'
        if (info.gt.0) print *, 'The updating process of DBDSDC did not converge.'
     endif
     call zgemm('N', 'N', nbands, cixdms, cixdms_m, cone, Uw, nbands, Vw, cixdms_m, czero, Ucix2, nbands)


     !WRITE(6,*)
     !do i=1,nbands
     !   do ip=1,cixdms
     !      WRITE(6, '(f14.10,1x,f8.5,3x)',advance='no') real(Ucix(i,ip)), aimag(Ucix(i,ip))
     !   enddo
     !   WRITE(6,*)
     !enddo
     print *, 'Singular values=', ws
     
     deallocate( work, rwork, iwork )

     DO iorb1=1,norbitals
        if ( cix_orb(iorb1).NE.icix ) CYCLE
        nind1 = nindo(iorb1)
        do ind1=1,nind1
           ip = iSx(ind1,iorb1)
           if (cind(ip).gt.0) DMFTU(:,ind1,iorb1,1) = Ucix2(:,cind(ip))
        enddo
     ENDDO

     
     if (nipc.gt.1) then
        ! Here we create UU, which can be used to get from original projector DMFTU the normalized projector by DMFTU*UU
        ! Since renormalized DMFTU was obtained by svd, i.e., DMFTU = U * s * Vt, we have
        ! UU = Vt^H * 1/s * Vt since renormalized DMFTU = U * Vt
        allocate( Vw2(cixdms_m,cixdms) )
        allocate( UU(cixdms,cixdms) )
        Vw2 = 0
        do ip=1,cixdms_m
           if ( abs(ws(ip)).ge.1e-10) then
              Vw2(ip,:) = 1./ws(ip) * Vw(ip,:)
           endif
        enddo
        call zgemm('C', 'N', cixdms, cixdms, cixdms_m, cone, Vw, cixdms_m, Vw2, cixdms_m, czero, UU, cixdms)
        deallocate(Vw2)
        
        do ipc=2,nipc
           Ucix(:,:) = 0.d0
           DO iorb1=1,norbitals
              if ( cix_orb(iorb1).NE.icix ) CYCLE
              nind1 = nindo(iorb1)
              do ind1=1,nind1
                 ip = iSx(ind1,iorb1)
                 if (cind(ip).gt.0) Ucix(:,cind(ip)) = DMFTU(:,ind1,iorb1,ipc)
              enddo
           ENDDO

           ! Ucix2 = dot(Ucix,UU)
           call zgemm('N', 'N', nbands, cixdms, cixdms, cone, Ucix, nbands, UU, cixdms, czero, Ucix2, nbands)
           
           DO iorb1=1,norbitals
              if ( cix_orb(iorb1).NE.icix ) CYCLE
              nind1 = nindo(iorb1)
              do ind1=1,nind1
                 ip = iSx(ind1,iorb1)
                 if (cind(ip).gt.0) DMFTU(:,ind1,iorb1,ipc) = Ucix2(:,cind(ip))
              enddo
           ENDDO
        enddo
        deallocate( UU )
     endif
     deallocate( ws, Uw, Vw )
     deallocate( cind, cini )
     deallocate( Ucix, Ucix2 )
  ENDDO
END SUBROUTINE RenormalizeTransK

REAL*8 FUNCTION ferm(x)
  IMPLICIT NONE
  REAL*8, intent(in) :: x
  if (x.LT.-100.) then
     ferm = 1.0
     RETURN
  endif
  if (x.GT.100.) then
     ferm = 0.0
     RETURN
  endif
  ferm = 1/(exp(x)+1)
  RETURN
END FUNCTION ferm

REAL*8 FUNCTION FreeE0(Energy, Temperature)
  IMPLICIT NONE
  REAL*8, intent(in) :: Energy, Temperature
  if (Energy/Temperature.LT.-200) then
     FreeE0 = Energy
     return
  endif
  if (Energy/Temperature.GT.200) then
     FreeE0 = 0.0
     return
  endif
  FreeE0 = -Temperature*log(1.+exp(-Energy/Temperature))
END FUNCTION FreeE0


REAL*8 FUNCTION Density(EF)
  USE defs, ONLY: PI, IMAG
  USE com,  ONLY: elecn
  USE muzero,ONLY: w_sum, w_gamma, w_norm1, nkp, w_nomega, w_npomega, wgh, abom, zEk, nemm, wprint, wmatsubara, n0_om, w_omega, jomq, nomq, iomq, w_beta, womq, wprint, wprint1, Qcheckbands
  USE dmfts, ONLY: Temperature
  IMPLICIT NONE
  REAL*8, intent(in)    :: EF
  ! functions
  REAL*8 :: ferm
  ! locals
  REAL*8,PARAMETER :: Ry2eV= 13.60569253d0 
  REAL*8    :: cc, zw0, wkp, ek0, Om, drest, ci, fer
  INTEGER   :: ikp, num, num0, iom, nemin, DM_nemin, DM_nemax, j0, j1, j2, i, nbands
  COMPLEX*16:: cE, ccn, cek, e0, e1, e2, ca, cb, omn, csum1, csum
  
  IF (wmatsubara) THEN
     if (wprint1) WRITE(*,'(A2,1x,A2,1x,A4,1x,A4,1x,A4,1x,A14,1x,A40,1x,A20,1x,A20)') 'ik','i','iom','j1','iomq','omn', 'cek', '1./(omn-cek)', 'dsum'
     zw0=0
     
     !$OMP PARALLEL DO PRIVATE(nemin,DM_nemin,DM_nemax,wkp,cc,num0,ek0,Om,drest,ci,iom,j0,j1,j2,e0,e1,e2,i,omn,cek,csum,csum1) SCHEDULE(STATIC) REDUCTION(+:zw0)
     do ikp=1,nkp
        nemin    = nemm(1,ikp)
        DM_nemin = nemm(2,ikp)
        DM_nemax = nemm(3,ikp)
        wkp = (wgh(ikp)/w_sum*2.0)*w_norm1
        if (DM_nemin.GT.nemin) zw0 = zw0 + wkp*(DM_nemin-nemin) ! These bands are fully filled
        
        cc=0
        DO num=DM_nemin,DM_nemax
           num0 = num-DM_nemin+1
           ek0 = dble(zEk(num0,w_nomega,ikp))-EF
           csum1=0
           csum=0
           do iom=1,n0_om-1
              cek = zEk(num0,iom,ikp)-EF
              omn = dcmplx(0, w_omega(iom) )
              csum1 = csum1 + 1./(omn-cek) - 1./(omn-ek0)
              if (wprint1) WRITE(*,'(I2,1x,I2,1x,I4,1x,I4,1x,I4,1x,f14.10,1x,2f20.15,1x,f20.15,1x,f20.15)') ikp,num,iom,iom, iom, aimag(omn), cek, dble(1./(omn-cek)-1./(omn-ek0)), dble(csum1)
           enddo
           
           do iom=n0_om,w_nomega-1
              j0 = jomq(iom-1)
              j1 = jomq(iom)
              j2 = jomq(iom+1)
              e0 = zEk(num0,iom-1,ikp)-EF
              e1 = zEk(num0,iom,  ikp)-EF
              e2 = zEk(num0,iom+1,ikp)-EF
              do i=1,nomq(iom)
                 omn = dcmplx(0, (2*iomq(iom,i)-1)*pi/w_beta)
                 if (iomq(iom,i).LT.j1) then
                    cek = e1 + (e0-e1)*(iomq(iom,i)-j1)/(j0-j1) ! linear interpolation
                 else if (iomq(iom,i).GT.j1) then
                    cek = e1 + (e2-e1)*(iomq(iom,i)-j1)/(j2-j1) ! linear interolation
                 else
                    cek = e1
                 endif
                 csum = csum + womq(iom,i)*(1./(omn-cek) - 1./(omn-ek0))
                 if (wprint1) WRITE(*,'(I2,1x,I2,1x,I4,1x,I4,1x,I4,1x,f14.10,1x,2f20.15,1x,f20.15,1x,f20.15)') ikp,num,iom,j1,iomq(iom,i),aimag(omn), cek, dble(womq(iom,i)*(1./(omn-cek) - 1/(omn-ek0))), dble(csum)
              enddo
           enddo
           iom = w_nomega
           j0 = jomq(iom-1)
           j1 = jomq(iom)
           e0 = zEk(num0,iom-1,ikp)-EF
           e1 = zEk(num0,iom,  ikp)-EF
           do i=1,nomq(iom)
              omn = dcmplx(0, (2*iomq(iom,i)-1)*pi/w_beta)
              cek = e1 + (e0-e1)*(iomq(iom,i)-j1)/(j0-j1)
              csum = csum + womq(iom,i)*(1./(omn-cek) - 1./(omn-ek0))
              if (wprint1) WRITE(*,'(I2,1x,I2,1x,I4,1x,I4,1x,I4,1x,f14.10,1x,2f20.15,1x,f20.15,1x,f20.15)') ikp,num,iom,j1,iomq(iom,i),aimag(omn), cek, dble(1./(omn-cek)), dble(csum)
           enddo
        
           Om = (2*jomq(w_nomega)*pi/w_beta)
           drest = (atan(ek0/Om)-atan(dble(cek)/Om))/pi
           ci = 2*dble(csum+csum1)/w_beta + drest + ferm(ek0*w_beta)
        
           if (Qcheckbands) then
              if (ci.GT.1.0) ci=1.0
              if (ci.LT.0.0) ci=0.0
           endif
           cc = cc + ci
        ENDDO
        zw0 = zw0 + wkp * cc
     enddo
     !$OMP END PARALLEL DO
  ELSE
     zw0=0
     do ikp=1,nkp
        nemin    = nemm(1,ikp)
        DM_nemin = nemm(2,ikp)
        DM_nemax = nemm(3,ikp)
        nbands = DM_nemax-DM_nemin+1
        wkp = (wgh(ikp)/w_sum*2.0)*w_norm1
        if (DM_nemin.GT.nemin) zw0 = zw0 + wkp*(DM_nemin-nemin) ! These bands are fully filled
        do iom=1,w_npomega
           fer = ferm(w_omega(iom)/Temperature)
           ca = abom(1,iom) + w_gamma*IMAG
           cb = abom(2,iom) + w_gamma*IMAG
           DO num=1,nbands
              cE = zEk(num,iom,ikp)
              ccn = log(cE-cb-EF)-log(cE-ca-EF)
              zw0 = zw0 - wkp * fer * aimag(ccn)/pi
              !if (wprint1) WRITE(991,*) ikp, iom, num, cE, ca, cb, EF, aimag(ccn)/pi!, zw0 
           ENDDO
        enddo
     enddo
  ENDIF
  if (wprint) WRITE(*,'(A,f13.8,A,f14.8,A,f10.4,A,f14.8)') 'EF[eV]=', EF*Ry2eV, ' Density=', zw0, ' N=', elecn, ' EF[Ry]=', EF
  Density = zw0-elecn
  RETURN
END FUNCTION Density

module Fermi
  contains
SUBROUTINE cmp_EF(EF, recomputeEF, wmax_nbands, Olapm, SOlapm, omega, sigma, nnlo, norbitals, natom, maxdim2, nomega, npomega, emin, emax, iorbital, cix_orb, cixdim, iSx, nindo, cfX, Rspin, sumw, wgamma, beta, vnorm1, LowE, lmaxp, SIMPLE, kmax, time_dmf0, time_dmf0w)
  USE xa,    ONLY: fj, dfj
  USE xa3,   ONLY: bk3, bk3lo, k3, k3lo, As, As_lo, aK
  USE param, ONLY: lmax2, nume, nmat, nloat, iblock, chunk
  USE structure, ONLY: rmt, BR1
  USE com,   ONLY: nat, elecn
  USE defs,  ONLY: PI, IMAG
  USE dmfts,   ONLY: DM_Emin, DM_Emax, DM_EF, maxsize, iso, ncix, Sigind, cix, csize, ll, nl, maxdim, Qrenormalize, matsubara, mixEF, projector, mode
  USE muzero,ONLY: w_sum, w_gamma, w_norm1, nkp, w_nomega, w_npomega, Ek, wgh, zEk, nemm, wprint, w_beta, w_omega, wmatsubara, wprint1, max_nbands, Qcheckbands, print_muzero
  USE p_project, ONLY: phi_jl, rix, P_rfi, kNri, l_al_ucase, j_al_ucase, phi_jl, n_al_ucase, dri
  USE com_mpi
  IMPLICIT NONE
  !COMPLEX*16, intent(in) :: Olapm(maxdim,maxdim,ncix), SOlapm(maxdim,maxdim,ncix)
  COMPLEX*16, allocatable, intent(in) :: Olapm(:,:,:), SOlapm(:,:,:)
  INTEGER, intent(in)    :: recomputeEF, lmaxp
  INTEGER, intent(inout) :: wmax_nbands
  REAL*8, intent(out)    :: EF, time_dmf0, time_dmf0w
  INTEGER, intent(out)   :: kmax(3)
  REAL*8, intent(in)     :: omega(nomega)
  COMPLEX*16, intent(in) :: sigma(maxsize,ncix,nomega)
  INTEGER, intent(in)    :: nnlo, norbitals, natom, maxdim2, nomega, npomega
  REAL*8, intent(in)     :: emin, emax, sumw, wgamma, beta, vnorm1
  INTEGER, intent(in)    :: iorbital(natom,lmaxp+1), cix_orb(norbitals), iSx(maxdim2, norbitals), nindo(norbitals), cixdim(ncix) !, noccur(maxsize,ncix)
  complex*16, intent(in) :: cfX(maxdim2,maxdim2,norbitals,norbitals), Rspin(2,2,norbitals)
  REAL*8, intent(in)     :: LowE(5,maxsize,ncix)
  LOGICAL, intent(in)    :: SIMPLE
  !
  REAL*8 :: Density
  ! interfaces
  interface
     REAL*8 Function Romb(y,N,dx)
       IMPLICIT NONE
       REAL*8, intent(in) :: y(N)
       REAL*8, intent(in) :: dx
       INTEGER, intent(in):: N
     end Function Romb
  end interface
  ! locals
  CHARACTER*100  :: Efilename
  REAL*8,PARAMETER :: Ry2eV= 13.60569253d0 
  complex*16, ALLOCATABLE :: DMFTU(:,:,:,:), STrans(:,:,:,:), Eij(:,:), zw1(:), tzek(:)
  INTEGER     :: ikp, nemin, nemax, DM_nemin, DM_nemax, n0!, num
  REAL*8      :: E(nume), emist, Kn(3)
  COMPLEX*16  :: gtot(nomega)
  LOGICAL     :: more_kpoints, Tcompute
  INTEGER     :: jatom, isize, nbands, itape, is, i, iorb, icix, ip, it, iom, itmax, fh_EF, iikp, ind, iband
  REAL*8      :: tEF, t1c, t2c, t1w, t2w, dd, dd2, dEF, Ea, Eb, fatol, xatol, xrtol, twgh, dsm, wkp, zw0, EFnew, rx
  LOGICAL     :: Need_EF_update
  INTEGER, ALLOCATABLE :: indHa(:,:)
  INTEGER     :: dNHa, ncandidates
  COMPLEX*16, ALLOCATABLE :: Ha(:,:)
  INTEGER     ::  lwork, lrwork, liwork, lhda
  complex*16, allocatable :: work(:)
  real*8,     allocatable :: rwork(:), Ebnd(:), weib(:), wg_bnd(:), Ebmin(:), Ebmax(:), candidate(:,:)
  integer,    allocatable :: iwork(:)
  integer     :: info, icn, iks, ivector, max_nume, ir, iind, l, Nri, nip
  COMPLEX*16  :: soo(maxsize,ncix)
  REAL*8, PARAMETER :: small_density_mismatch = 1.0
  REAL*8, allocatable     :: aKR(:), jlr(:), jlp(:)

  if (abs(projector).le.5) then ! Need Olapm, SOlapm
     if (.not.allocated(SOlapm)) then
        print *, 'ERROR: SOlapm should have been allocated for this projector. Check cmp_EF in dmft2 step.'
        call stop_MPI
        STOP 'ERROR: SOlapm not allocated'
     endif
     if ( UBOUND(SOlapm, DIM = 1) .ne. maxdim .or. UBOUND(SOlapm, DIM = 2) .ne. maxdim .or. UBOUND(SOlapm, DIM = 3) .ne. ncix ) then
        print *, 'ERROR: SOlapm does not have correct dimensions. Check cmp_EF in dmft2 step: ', UBOUND(SOlapm, DIM = 1), UBOUND(SOlapm, DIM = 2), UBOUND(SOlapm, DIM = 3)
        call stop_MPI
        STOP 'ERROR: SOlapm wrong dimension'
     endif
  endif

  if (recomputeEF.gt.1) then
     soo(:,:) = LowE(5,:,:)  ! self-energy at zero frequency
     ! Will use a quasiparticle approximation for Mott insulators 
     ! to compute the chemical potential for insulators!
     ALLOCATE( indHa(maxdim*ncix,5) )
     dNHa=0 ! Number of poles which will be added
     DO iorb=1,norbitals
        icix = cix_orb(iorb)
        if (icix.EQ.0) CYCLE
        do ind=1,nindo(iorb)
           ip = iSx(ind,iorb)
           it = Sigind(ip,ip,icix)
           if (it.gt.0 .and. abs(LowE(1,it,icix)-1.0).LT.1e-6) then
              dNHa = dNHa + 1
              indHa(dNHa,1) = ip
              indHa(dNHa,2) = ind
              indHa(dNHa,3) = iorb
              indHa(dNHa,4) = icix
              !WRITE(*,'(A,I2,1x,A,I2,1x,A,I2,1x,A,I2,1x,A,I2,1x,A,I2)') 'iorb=', iorb, 'icix=', icix, 'ind=', ind, 'ip=', ip, 'it=', it, 'dNHa=', dNHa
           endif
        ENDDO
     ENDDO
  endif
  
  Qcheckbands=.TRUE.
  def = 0.02
  itmax = 500
  fatol = 1e-10
  xatol = 1e-9
  xrtol = 1e-9
  kmax(:)=0
  if (wmax_nbands.EQ.0) then
     DO ivector=1,nvector
        itape=9
        if (vector_para) then
           open(itape,FILE=fvectors(ivector,1),STATUS='old',FORM='unformatted')
        else
           rewind(itape) !--- both vector files: 9 and 10 rewind -------------------------!
        endif
        DO I=1,NAT
           READ(itape) EMIST !---- At the beginninge we have linearization energies --!
           READ(itape) EMIST !---- We just skip them ---------------------------------!
        ENDDO
        DO iks=1,vectors(ivector,2)   ! kpoint loop begin
           CALL Read_Vec_Spin_DontStore(nemin, DM_nemin, DM_nemax, more_kpoints, n0, emin, emax, nnlo, kmax)
           IF (.NOT.more_kpoints) EXIT
           wmax_nbands = max(DM_nemax-DM_nemin+1, wmax_nbands)
        ENDDO
        if (vector_para) then
           close(itape)
        else
           rewind(itape)
        endif
     END DO
     CALL Reduce_maxbands(wmax_nbands)
     CALL FindMaxK_MPI(kmax)
  endif
  max_nbands=wmax_nbands
  
  !print *, 'MAX_NBANDS=', max_nbands
  CALL Reduce_maxnume(max_nume, nume)
  
  if (recomputeEF.gt.1) then
     lhda = max_nbands+dNHa
     ALLOCATE( Ha(lhda,lhda) )
     ALLOCATE( Ebnd(lhda), weib(max_nume+dNHa), wg_bnd(lhda), Ebmin(max_nume+dNHa), Ebmax(max_nume+dNHa) )
  
     lwork = 2*lhda+lhda*lhda
     lrwork =  5*lhda + 2*lhda*lhda + 1
     liwork = 3+5*lhda
     ALLOCATE ( work(lwork), rwork(lrwork), iwork(liwork) )
  endif

  w_sum = sumw; w_gamma = wgamma; w_norm1 = vnorm1
  w_nomega = nomega; w_npomega = npomega; w_beta=beta
  w_omega(:) = omega(:)

  call cputim(t1c)
  call walltim(t1w)
  
  ! Ek, zEk, nemm
  ALLOCATE( wgh(pr_procr), Ek(max_nume,pr_procr), zEk(max_nbands,nomega,pr_procr), nemm(3,pr_procr) )
  wprint = .TRUE.
  
  Ek(:,:) = 0.0
  zEk(:,:,:) = 0.0
  wgh(:) = 0

  if (recomputeEF.gt.1) then
     weib(:) = 0
     Ebmin(:) = 1e10
     Ebmax(:) = -1e10
     gtot=0
  endif

  nkp=0
  iikp=0
  DO ivector=1,nvector
     DO is=1,iso    !------ over up/dn ---------------------!
        itape=8+is
        if (vector_para) then
           open(itape,FILE=fvectors(ivector,is),STATUS='old',FORM='unformatted')
        else
           rewind(itape) !--- both vector files: 9 and 10 rewind -------------------------!
        endif
        DO I=1,NAT
           READ(itape) EMIST !---- At the beginninge we have linearization energies --!
           READ(itape) EMIST !---- We just skip them ---------------------------------!
        ENDDO
     END DO
     
     DO iks=1,vectors(ivector,2) ! kpoint loop begin
        if (vector_para) then
           ikp = vectors(ivector,3)+iks  ! successive index in k-point table from case.klist
           iikp = iikp+1                 ! successive index in k-point table on this processor
        else
           ikp = iks                     ! successive index in k-point table from case.klist
           !--- We need to go over all k-points even though we will compute only some of them on this processor.
           !--- This is because we need to read vector file sequentially.
           iikp = ikp-myrank*pr_proc            ! The index of the point to compute. If negative, do not compute! 
        endif
        Tcompute=.FALSE.                                                                                                                                                                        
        if (iikp.gt.0) Tcompute=.TRUE.       ! If Tcompute is true, the point needs to be computed.                                                                                             
        if (iikp.GT.pr_proc) EXIT            ! Processor finished. The rest of the points will be taken care of by the other processors.                                                        
        
        CALL Read_Vec_Spin(E, twgh, As, As_lo, k3, k3lo, bk3, bk3lo, nemin, nemax, DM_nemin, DM_nemax, more_kpoints, n0, emin, emax, DM_Emin, DM_Emax, iso, nmat, nume, nnlo)
        IF (.NOT.Tcompute) CYCLE
        IF (.NOT.more_kpoints) EXIT

        wgh(iikp) = twgh
        
        isize=n0-nnlo
        DO I=1,isize                                !----- Over reciprocal vectors -----------------------!
           Kn(:) = matmul(BR1,BK3(:,I))             !----  Creates K+k in cartesian coordinates ----------!
           aK(I) = sqrt(Kn(1)**2+Kn(2)**2+Kn(3)**2) !----  calculates |K+k| for use in bessel functions --!
        ENDDO
        do jatom=1,nat
           CALL HARMON(isize,aK(:isize),lmax2,fj(:,:isize,jatom),dfj(:,:isize,jatom),rmt(jatom))
        enddo
        
        ! New
        if (abs(projector).ge.5) then
           allocate( phi_jl(nmat, n_al_ucase) )
           
           Nri=2**kNri+1    ! Number of radial points in the interstitials
           phi_jl(:,:)=0
           allocate( aKR(Nri), jlr(Nri), jlp(Nri) )
           DO iind=1,n_al_ucase
              if (abs(dri(iind)).gt.1e-10) then
                 l     = l_al_ucase(iind)
                 jatom = j_al_ucase(iind)
                 DO i=1,isize         ! over all reciprocal vectors K
                    rx = RMT(jatom)
                    do ir=1,Nri
                       aKR(ir)=rx*aK(i)    ! |k+K|*r
                       rix(ir)=rx          !  r
                       rx = rx + dri(iind)
                    enddo
                    CALL sphbes2(l,Nri,aKR,jlr)  ! spherical bessel : j_l(|k+K|*r)
                    jlr(:) = jlr(:)*rix(:)       ! r*j_l . We need to do that, because the integral is Int[ (phi(r)/r)*j_l(r)*r^2]=Int[phi(r)*j_l(r)*r]
                    jlp(:) = jlr(:)*P_rfi(:,iind)
                    phi_jl(i,iind) = romb(jlp, Nri, dri(iind))  ! Integration over r on rix mesh
                 ENDDO
              endif
           ENDDO
           deallocate( aKR, jlr, jlp )
        endif
        ! New
        
        nbands = DM_nemax-DM_nemin+1
        nip=1
        allocate( DMFTU(nbands,maxdim2,norbitals,nip) )
        CALL Build_DMFT_Projector(DMFTU, cfX, Rspin, iorbital, norbitals, nip, n0, nnlo, nbands, cix_orb, nindo, DM_nemin, DM_nemax, maxdim2, lmaxp)

        if (abs(projector).ge.5) deallocate( phi_jl )
        
        if (Qrenormalize) then
           if (abs(projector).le.5) then
              CALL RenormalizeTrans(DMFTU, Olapm, SOlapm, cix_orb, cixdim, nindo, iSx, nbands, nip, maxdim2, norbitals, maxdim, ncix, SIMPLE)
           else
              CALL RenormalizeTransK(DMFTU, cix_orb, cixdim, nindo, iSx, Sigind, projector, nbands, nip, maxdim2, norbitals, maxdim, ncix)
           endif
        endif


        
        allocate( STrans(maxsize,ncix,nbands,nbands) )
        CALL CompressSigmaTransformation1(STrans, DMFTU(:,:,:,1), Sigind, iSx, cix, iorbital, ll, nl, natom, iso, ncix, maxdim, maxdim2, lmaxp, norbitals, nbands, maxsize)
        
        allocate( Eij(nbands,nbands) )
        
        if (recomputeEF.gt.1) then
        
           ! Using pole approximation for the self-energy and the trick described in 
           ! PRL 96, 036404 (2006)
           ! to obtain good approximation for the chemical potential. 
           ! 
           Eij=0
           DO iband=1,nbands
              Eij(iband,iband) = E(iband+DM_nemin-1)  !---   preparing hamiltonian: epsk+sigma
           ENDDO
           CALL AddSigma_optimized2(Eij, soo, STrans, csize, 1, nbands, ncix, maxsize)
           Ha=0
           Ha(:nbands,:nbands) = Eij(:,:)
           DO i=1,dNHa
              ip  = indHa(i,1)
              ind = indHa(i,2)
              iorb = indHa(i,3)
              icix = indHa(i,4)
              it = Sigind(ip,ip,icix)
              do iband=1,nbands
                 Ha(iband,nbands+i) = DMFTU(iband,ind,iorb,1)*sqrt(LowE(2,it,icix))
                 Ha(nbands+i,iband) = conjg(DMFTU(iband,ind,iorb,1))*sqrt(LowE(2,it,icix))
              enddo
           ENDDO
           DO i=1,dNHa
              ip  = indHa(i,1)
              ind = indHa(i,2)
              iorb = indHa(i,3)
              icix = indHa(i,4)
              it = Sigind(ip,ip,icix)
              Ha(nbands+i,nbands+i) = LowE(3,it,icix)+DM_EF
           ENDDO
           
           CALL ZHEEVD('V','U', nbands+dNHa, Ha, lhda, Ebnd, work, lwork, rwork, lrwork, iwork, liwork, info )
           if (info .ne. 0) then
              print *, 'Diagonalization of extended Hamiltonian in computing EF failed. Info-zheevd=', info
           endif
        
           wkp = (twgh/w_sum*2.0)*w_norm1
           zw0=0
           do i=1,DM_nemin-1 ! Fully filled bands not treated by DMFT
              zw0 = zw0+1.
              weib(i) = weib(i)+ zw0*wkp
              if (E(i).GT.Ebmax(i)) Ebmax(i) = E(i)
              if (E(i).LT.Ebmin(i)) Ebmin(i) = E(i)
              !WRITE(*,'(I3,1x,F12.7,1x,F10.7,1x,F10.7,1x,F10.7)') i, E(i)*Ry2eV, 1., zw0, weib(i)
           enddo
           DO i=1,nbands+dNHa ! Bands treated by DMFT
              dsm=0.0
              DO iband=1,nbands
                 dsm = dsm + abs(Ha(iband,i))**2 ! DMFT bands carry a fractional charge
              ENDDO
              wg_bnd(i) = dsm  ! The fractional charge of Mott-bands
              zw0 = zw0 + dsm  ! The total weight up to this band
              weib(DM_nemin-1+i) = weib(DM_nemin-1+i) + zw0*wkp  ! storing the total weight (summed over k) up to this band
              if (Ebnd(i).GT.Ebmax(DM_nemin-1+i)) Ebmax(DM_nemin-1+i) = Ebnd(i)  ! Useful to estimate the gap
              if (Ebnd(i).LT.Ebmin(DM_nemin-1+i)) Ebmin(DM_nemin-1+i) = Ebnd(i)  ! Useful to estimate the gap
              !WRITE(*,'(I3,1x,F12.7,1x,F10.7,1x,F10.7,1x,F10.7)') i+DM_nemin-1, Ebnd(i)*Ry2eV, dsm, zw0, weib(DM_nemin-1+i)
           ENDDO
        endif
        
        deallocate( DMFTU )
        allocate( zw1(nbands) )
        deallocate( Eij )
        
        !! OMP
        !! private: Eij, tzek, i,
        !! shared(zEk)

        !$OMP  PARALLEL DO SHARED(zEk) PRIVATE(i,Eij,tzek) SCHEDULE(STATIC,CHUNK)
        do iom=1,npomega
           allocate( Eij(nbands,nbands), tzek(nbands) )
           Eij=0
           DO i=1,nbands
              Eij(i,i) = E(i+DM_nemin-1)  !---   preparing hamiltonian: epsk+sigma
           ENDDO
           ! Adding self-energy
           CALL AddSigma_optimized2(Eij, sigma(:,:,iom), STrans, csize, 1, nbands, ncix, maxsize)
           CALL eigvals(Eij, tzek, nbands)
           zEk(:nbands,iom,iikp) = tzek(:)
           deallocate( Eij, tzek )
        enddo
        !$OMP END PARALLEL DO
        
        DEALLOCATE( zw1 )
        
        nemm(1,iikp) = nemin
        nemm(2,iikp) = DM_nemin
        nemm(3,iikp) = DM_nemax
        Ek(:nume,iikp) = E(:nume)
        
        DEALLOCATE( STrans )
        WRITE(*,'(I3,A,1x,I3,1x,I3,1x,A,I4)') myrank, ') Finished k-point number', ikp, iikp, 'with #bands=', nbands
        nkp=nkp+1
        if (Qprint) WRITE(6,'(I3,A,1x,I3,1x,I3,1x,A,I4)') myrank, ') Finished k-point number', ikp, iikp, 'with #bands=', nbands
     ENDDO         ! kpoint loop ends
     DO is=1,iso    !------ over up/dn ---------------------!
        itape=8+is
        if (vector_para) then
           close(itape)
        else
           rewind(itape)
        endif
     ENDDO
  END DO           ! ivector loop ends

  ! Ebmax <- max(Ebmax)
  ! Ebmin <- min(Ebmin)
  ! weib(i) <- sum(weib)
  ! gtot(iom) <- sum(gtot)
  
  CALL Reduce_MPI1(pr_procr, nomega, max_nume)
  IF (recomputeEF.GT.1) CALL Reduce_MPI1b(nomega, Ebmax, Ebmin, weib, gtot, max_nume+dNHa)

  if (mode.EQ.'e') then
     if (myrank.EQ.master) then
        wmatsubara = matsubara
        wprint=.TRUE.
        wprint1=.FALSE.
        !dd = Density(DM_EF)
        !print *, 'at', DM_EF, 'dens-Z=',  dd
        call print_muzero(15,elecn,DM_EF)
     endif
     return
  endif
     
  !wprint1=.false.
  if (myrank.EQ.master) then
     wmatsubara = matsubara

     Need_EF_update=.True.
     
     if (recomputeEF.gt.1) then
        ! Print g computed by quasiparticle approximation.
        ! Only for debugging purposes.
        !open(599,file='gtot.debug1',status='unknown',form='formatted')
        !DO iom=1,nomega
        !   WRITE(599,*) omega(iom), dble(gtot(iom))/2., aimag(gtot(iom))/2.
        !ENDDO
        !close(599)
        
        wprint=.FALSE.
        ncandidates=0
        ALLOCATE(candidate(DM_nemax+dNHa-1,5))
        candidate=0
        DO i=DM_nemin,DM_nemax+dNHa-1
           IF (Ebmax(i).LT.Ebmin(i+1) .AND. abs(weib(i)-elecn)<5.0) THEN ! There is a gap and the density is not crazy!
              ncandidates = ncandidates + 1       ! One more gap. Is this the best gap?
              candidate(ncandidates,1) = 0.5*(Ebmin(i+1)+Ebmax(i)) ! The middle of the gap.
              candidate(ncandidates,2) = Density(candidate(ncandidates,1))  ! Actual density in the middle of the gap
              candidate(ncandidates,3) = weib(i)-elecn ! The estimated density by quasiparticle approximation
              candidate(ncandidates,4) = (Ebmin(i+1)-Ebmax(i)) ! Quasiparticle estimation for the gap size
              WRITE(*,'(I3,1x,F13.7,1x,F13.7,1x,F13.7,2x,A,F8.3,2x,A,F8.3,2x,A,F8.3)') i, Ebmax(i)*Ry2eV, Ebmin(i+1)*Ry2eV, weib(i), 'gap-size=', (Ebmin(i+1)-Ebmax(i))*Ry2eV, 'middle-point=', 0.5*(Ebmin(i+1)+Ebmax(i))*Ry2eV, 'charge-density=', elecn+candidate(ncandidates,2)
           ELSE
              WRITE(*,'(I3,1x,F13.7,1x,F13.7,1x,F13.7)') i, Ebmax(i)*Ry2eV, Ebmin(i+1)*Ry2eV, weib(i)
           ENDIF
        ENDDO
        if (ncandidates.gt.0) then
           ! Which candidate has best density
           IF (recomputeEF.GT.2) THEN ! For recomputeEF=3 we take take the one with density closer to charge neutrality, but density of electrons should larger or equal than required by charge neutrality
              icn=ncandidates
              !DO ip=1,ncandidates ! Finding the best candidate
              !   IF (candidate(ip,2).GE.0.0 ) THEN
              !      icn=ip
              !      EXIT
              !   ENDIF
              !ENDDO
              !DO ip=1,ncandidates ! Finding the best candidate
              !   IF (abs(candidate(ip,2)).LT.abs(candidate(icn,2)) .AND. candidate(ip,2).GE.0.0 ) icn=ip
              !ENDDO
           ELSE  ! For recomputeEF=2 we always take the one with density closer to charge neutrality
              icn=1 
              DO ip=2,ncandidates ! Finding the best candidate
                 IF (abs(candidate(ip,2)).LT.abs(candidate(icn,2))) icn=ip
              ENDDO
           ENDIF
           IF (abs(candidate(icn,2)).LT.small_density_mismatch) THEN ! Tolerable mismatch in density
              tEF = candidate(icn,1)
              Need_EF_update=.False.
              WRITE(*,*) 'Found a quasiparticle-gap. EF is placed in the middle of this gap!'
              WRITE(*,'(A,F13.7,2x,A,F13.7,2x,A,F13.7,2x)') 'Estimated-gap-size=', candidate(icn,4)*Ry2eV, 'mismatch-density-integration=', candidate(icn,2), 'mismatch-density-quasiparticle=', candidate(icn,3)
           ENDIF
        endif
        DEALLOCATE( candidate )
        IF (Need_EF_update) print *, 'Failed to put EF in the gap. Computing EF in the regular way.'
     endif

     if (Need_EF_update) then
        wprint=.TRUE.
        wprint1=.FALSE.
        
        tEF = DM_EF
        dd = Density(tEF)
        
        def = abs(dd)/elecn*1.
        if (def.LT.0.00001) def=0.00001
        if (def.GT.0.05)  def=0.05
        ! new
        !print *, 'dd=', dd, 'de=', def
        dd2 = Density(tEF-def*sign(1.0d0,dd))
        !print *, 'dd2=', dd2, 'de=', def
        def = 0.5 * abs(dd)*def/abs(dd2-dd)  ! 0.5*|n-n0|*1/(dn/de)
        !print *, 'new-def=', def
        if (def.LT.0.00001) def=0.00001
        ! new

        if (abs(dd).gt.1e-10) then
           if (dd.lt.0) then
              do while (dd.lt.0)
                 tEF = tEF + def
                 dd = Density(tEF)
              end do
              Ea = tEF-def
              Eb = tEF
           else
              do while (dd.gt.0)
                 tEF = tEF - def
                 dd = Density(tEF)
              enddo
              Ea = tEF
              Eb = tEF+def
           endif
           wprint = .FALSE.
           CALL brent ( Density, fatol, itmax, Ea, Eb, xatol, xrtol )
           tEF = Ea
        endif
        
     endif
     
     if (abs(tEF-DM_EF)*Ry2eV .lt. 0.03) then ! If we are close to convergence, we start mixing EF
        EFnew = DM_EF + (tEF-DM_EF)*mixEF !!! mixing
     else ! When EF is changing a lot, one should just take the new EF.
        EFnew = tEF
        WRITE(*,*) 'Difference in EF too large to mix ', (tEF-DM_EF)*Ry2eV, 'cutoff=', 0.03
     endif
     WRITE(*,'(A,f15.10,A,f15.10,A,f15.10)') 'Difference in EF=', (tEF-DM_EF)*Ry2eV, ' EF set to : ', EFnew*Ry2eV, ' new EF (no-mixing): ', tEF*Ry2eV
     WRITE(*,*)
     WRITE(21,'(A,f20.15,1x,f20.15)') ':DEF', (tEF-DM_EF)*Ry2eV, EF*Ry2eV
     EF = EFnew
     
     fh_EF=1000
     Efilename = 'EF.dat'
     OPEN(fh_EF,file=Efilename,status='unknown',form='formatted')
     WRITE(fh_EF,'(f25.16)') EF*Ry2eV
     CLOSE(fh_EF)
     
  endif


  CALL Bcast_MPI(EF)
  
  if (recomputeEF.gt.1) then
     DEALLOCATE( Ha, Ebnd, weib, wg_bnd, Ebmin, Ebmax )
     DEALLOCATE ( work, rwork, iwork )
     DEALLOCATE( indHa )
  endif
  
  DEALLOCATE( wgh, Ek, zEk, nemm )
  call cputim(t2c)
  call walltim(t2w)
  time_dmf0 = t2c-t1c
  time_dmf0w= t2w-t1w

END SUBROUTINE cmp_EF
end module Fermi


SUBROUTINE CountFrequencyArrays(n0_om, qmax, beta, omega, nomega, matsubara)
  USE defs,  ONLY: PI
  USE com_mpi
  IMPLICIT NONE
  INTEGER,intent(out):: n0_om, qmax
  REAL*8, intent(out):: beta
  REAL*8, intent(in) :: omega(nomega)
  INTEGER,intent(in) :: nomega
  LOGICAL,intent(in) :: matsubara
  ! locals
  REAL*8,PARAMETER :: Ry2eV= 13.60569253d0 
  INTEGER :: i,j,ij,k,iom
  INTEGER :: nomq(nomega)
  
  qmax=1
  beta=0
  IF (matsubara) THEN
     beta = PI/omega(1)
     do iom=1,nomega
        if (abs((omega(iom)*beta/PI+1)/2. - iom).GT.1e-1) EXIT
     enddo
     n0_om = iom-1

     nomq(:)=0
     do iom=n0_om,nomega-1
        i = nint((omega(iom)*beta/PI+1)/2.)
        j = nint((omega(iom+1)*beta/PI+1)/2.)
        ij = (i+j)/2
        do k=i,ij
           nomq(iom)=nomq(iom)+1
        enddo
        do k=ij,j-1
           nomq(iom+1)=nomq(iom+1)+1
        enddo
     enddo
     nomq(nomega)=nomq(nomega)+1
     do iom=n0_om,nomega
        qmax = max(qmax,nomq(iom))
     enddo

     if (Qprint) WRITE(6,*) 'beta....=', beta/Ry2eV, 'n0_om...=', n0_om
  ENDIF
END SUBROUTINE CountFrequencyArrays

SUBROUTINE CreateFrequencyArrays(nomq, jomq, iomq, womq, abom, npomega, n0_om, qmax, beta, omega, nomega, matsubara, Temperature)
  USE defs,  ONLY: PI
  IMPLICIT NONE
  INTEGER,intent(in) :: n0_om, qmax, nomega
  INTEGER,intent(out):: nomq(nomega), jomq(nomega), iomq(nomega,qmax)
  REAL*8, intent(out):: womq(nomega,qmax), abom(2,nomega)
  INTEGER,intent(out):: npomega
  REAL*8, intent(in) :: beta, Temperature, omega(nomega)
  LOGICAL,intent(in) :: matsubara
  ! externals
  REAL*8  :: ferm
  ! locals
  REAL*8,PARAMETER :: Ry2eV= 13.60569253d0 
  REAL*8  :: Large, small, fer
  INTEGER :: i,j,ij,k,iom,ic
  ! jomq(iom)   :: iom is index in logarithmic mesh, and jomq(iom) gives the true sequential number in Matsubara mesh
  ! nomq(iom)   :: how many Matsubara points will need this point from logarithmic mesh to compute value
  ! iomq(iom,:) :: all Matsubara points
  ! womq(iom,:) :: weight
  nomq(:)=0; iomq(:,:)=0; jomq(:)=0; womq(:,:)=0
  abom(:,:)=0
  npomega = nomega

  IF (matsubara) THEN
     do iom=1,n0_om-1
        nomq(iom)=1
        iomq(iom,1)=iom
        womq(iom,1)=1.0
        jomq(iom)=iom
     enddo

     do iom=n0_om,nomega-1
        i = nint((omega(iom)*beta/PI+1)/2.)     ! current matsubara frequency from log mesh (in integer representation)
        j = nint((omega(iom+1)*beta/PI+1)/2.)   ! next matsubara frequency from log mesh (in integer representation)
        ij = (i+j)/2                      ! center 
        jomq(iom) = i

        if (MOD(j-i,2).eq.1) then
           do k=i,ij
              nomq(iom)=nomq(iom)+1
              iomq(iom,nomq(iom))=k
              womq(iom,nomq(iom))=1.0
           enddo
           do k=ij+1,j-1
              nomq(iom+1)=nomq(iom+1)+1
              iomq(iom+1,nomq(iom+1))=k
              womq(iom+1,nomq(iom+1))=1.0
           enddo
        else
           do k=i,ij
              nomq(iom)=nomq(iom)+1
              iomq(iom,nomq(iom))=k
              womq(iom,nomq(iom))=1.0
           enddo
           womq(iom,nomq(iom))=0.5
           ic = nomq(iom+1)+1
           do k=ij,j-1
              nomq(iom+1)=nomq(iom+1)+1
              iomq(iom+1,nomq(iom+1))=k
              womq(iom+1,nomq(iom+1))=1.0
           enddo
           womq(iom+1,ic)=0.5
        endif
     enddo
     k = nint((omega(nomega)*beta/PI+1)/2.)
     nomq(nomega)=nomq(nomega)+1
     iomq(nomega,nomq(nomega))=k
     womq(nomega,nomq(nomega))=1.0
     jomq(nomega)=k

     if (.False.) then
        do iom=1,npomega
           WRITE(990,'(A,I4,A,I4,A,I4,1x,A,F15.8)') 'iom= ',iom,' jomq= ',jomq(iom), ' nomq= ',nomq(iom), 'w=', omega(iom)*Ry2eV
           do i=1,nomq(iom)
              WRITE(990,'(A,I4,A,I4,A,F8.4)') '   ', i, ' iomq= ', iomq(iom,i), ' womq= ', womq(iom,i)
           enddo
        enddo
        close(990)
     endif
     
     
  ELSE
     Large = 10000.
     small = 1e-6
     !-- creates integration intervals
     do npomega=1,nomega
        fer = ferm(omega(npomega)/Temperature)
        if (fer.LT.small) exit
     enddo
     
     abom(1,1) = -Large
     abom(2,1) = 0.5*(omega(1)+omega(2))
     do iom=2,npomega-1
        abom(1,iom) = 0.5*(omega(iom)+omega(iom-1)) ! always start at midpoint
        abom(2,iom) = 0.5*(omega(iom)+omega(iom+1)) ! and end at midpoint
     enddo
     abom(1,npomega) = 0.5*(omega(npomega)+omega(npomega-1))
     if (npomega.LT.nomega) then
        abom(2,npomega) = 0.5*(omega(npomega)+omega(npomega+1))
     else
        abom(2,npomega) = omega(npomega)
     endif
  ENDIF
END SUBROUTINE CreateFrequencyArrays

SUBROUTINE CountFrequencyArrays2(n0_om, qmax, beta, omega, nomega, matsubara)
  USE defs,  ONLY: PI
  USE com_mpi
  IMPLICIT NONE
  INTEGER,intent(out):: n0_om, qmax
  REAL*8, intent(out):: beta
  REAL*8, intent(in) :: omega(nomega)
  INTEGER,intent(in) :: nomega
  LOGICAL,intent(in) :: matsubara
  ! locals                                                                                                                                                                                                  
  REAL*8,PARAMETER :: Ry2eV= 13.60569253d0 
  INTEGER :: i,j,ij,k,iom
  INTEGER :: nomq(nomega)

  qmax=1
  beta=0
  IF (matsubara) THEN
     beta = PI/omega(1)
     do iom=1,nomega
        if (abs((omega(iom)*beta/PI+1)/2. - iom).GT.1e-1) EXIT
     enddo
     n0_om = iom-1

     nomq(:)=0
     do iom=n0_om,nomega-1
        i = nint((omega(iom)*beta/PI+1)/2.)
        j = nint((omega(iom+1)*beta/PI+1)/2.)
        ij = (i+j)/2
        do k=i,j-1
           nomq(iom)=nomq(iom)+1
        enddo
        do k=i+1,j-1
           nomq(iom+1)=nomq(iom+1)+1
        enddo
     enddo
     nomq(nomega)=nomq(nomega)+1
     do iom=n0_om,nomega
        qmax = max(qmax,nomq(iom))
     enddo

     if (Qprint) WRITE(6,*) 'beta....=', beta/Ry2eV, 'n0_om...=', n0_om
  ENDIF
END SUBROUTINE CountFrequencyArrays2

SUBROUTINE CreateFrequencyArrays2(nomq, jomq, iomq, womq, abom, npomega, n0_om, qmax, beta, omega, nomega, matsubara, Temperature)
  USE defs,  ONLY: PI
  IMPLICIT NONE
  INTEGER,intent(in) :: n0_om, qmax, nomega
  INTEGER,intent(out):: nomq(nomega), jomq(nomega), iomq(nomega,qmax)
  REAL*8, intent(out):: womq(nomega,qmax), abom(2,nomega)
  INTEGER,intent(out):: npomega
  REAL*8, intent(in) :: beta, Temperature, omega(nomega)
  LOGICAL,intent(in) :: matsubara
  ! externals                                                                                                                                                                                                                                                                                                                                                        
  REAL*8  :: ferm
  ! locals                                                                                                                                                                                                                                                                                                                                                           
  REAL*8  :: Large, small, fer
  INTEGER :: i,j,k,iom!,ij,ic
  ! jomq(iom)   :: iom is index in logarithmic mesh, and jomq(iom) gives the true sequential number in Matsubara mesh                                                                                                                                                                                                                                                
  ! nomq(iom)   :: how many Matsubara points will need this point from logarithmic mesh to compute value                                                                                                                                                                                                                                                             
  ! iomq(iom,:) :: all Matsubara points                                                                                                                                                                                                                                                                                                                              
  ! womq(iom,:) :: weight                                                                                                                                                                                                                                                                                                                                            
  nomq(:)=0; iomq(:,:)=0; jomq(:)=0; womq(:,:)=0
  abom(:,:)=0
  npomega = nomega

  IF (matsubara) THEN
     do iom=1,n0_om-1
        nomq(iom)=1
        iomq(iom,1)=iom
        womq(iom,1)=1.0
        jomq(iom)=iom
     enddo

     do iom=n0_om,nomega-1
        i = nint((omega(iom)*beta/PI+1)/2.)     ! current matsubara frequency from log mesh (in integer representation)                                                                                                                                                                                                                                              
        j = nint((omega(iom+1)*beta/PI+1)/2.)   ! next matsubara frequency from log mesh (in integer representation)                                                                                                                                                                                                                                                 
        jomq(iom) = i
        do k=i,j-1
           nomq(iom)=nomq(iom)+1
           iomq(iom,nomq(iom))=k
           womq(iom,nomq(iom))=dble(j-k)/dble(j-i)
        enddo
        do k=i+1,j-1
           nomq(iom+1)=nomq(iom+1)+1
           iomq(iom+1,nomq(iom+1))=k
           womq(iom+1,nomq(iom+1))=dble(k-i)/dble(j-i)
        enddo
     enddo
     k = nint((omega(nomega)*beta/PI+1)/2.)
     nomq(nomega)=nomq(nomega)+1
     iomq(nomega,nomq(nomega))=k
     womq(nomega,nomq(nomega))=1.0
     jomq(nomega)=k

     if (.False.) then
        do iom=1,npomega
           WRITE(990,'(A,I4,A,I4,A,I4,A,F15.8)') 'iom= ',iom,' jomq= ',jomq(iom), ' nomq= ',nomq(iom), 'w=', omega(iom)
           do i=1,nomq(iom)
              WRITE(990,'(A,I4,A,I4,A,F8.4)') '   ', i, ' iomq= ', iomq(iom,i), ' womq= ', womq(iom,i)
           enddo
        enddo
     endif
  ELSE
     Large = 10000.
     small = 1e-6
     !-- creates integration intervals                                                                                                                                                                                                                                                                                                                               
     do npomega=1,nomega
        fer = ferm(omega(npomega)/Temperature)
        if (fer.LT.small) exit
     enddo
     if (npomega.GT.nomega) npomega=nomega
     
     abom(1,1) = -Large
     abom(2,1) = 0.5*(omega(1)+omega(2))
     do iom=2,npomega-1
        abom(1,iom) = 0.5*(omega(iom)+omega(iom-1)) ! always start at midpoint                                                                                                                                                                                                                                                                                       
        abom(2,iom) = 0.5*(omega(iom)+omega(iom+1)) ! and end at midpoint                                                                                                                                                                                                                                                                                            
     enddo
     abom(1,npomega) = 0.5*(omega(npomega)+omega(npomega-1))
     if (npomega.LT.nomega) then
        abom(2,npomega) = 0.5*(omega(npomega)+omega(npomega+1))
     else
        abom(2,npomega) = omega(npomega)
     endif

     do iom=1,nomega
        jomq(iom)=iom
     enddo

  ENDIF
END SUBROUTINE CreateFrequencyArrays2



SUBROUTINE GiveIntLimits(a,b,omega,iom, nomega)
  IMPLICIT NONE
  REAL*8, intent(out) :: a, b
  INTEGER, intent(in) :: iom, nomega
  REAL*8, intent(in)  :: omega(nomega)
  !print *, 'Before GiveLimits'
  REAL*8 :: Large
  Large = 10000.
  if (iom.EQ.1) then
     a = -Large
     b = 0.5*(omega(iom)+omega(iom+1))
  else if (iom.EQ.nomega) then
     a = 0.5*(omega(iom)+omega(iom-1))
     b = 100.
  else
     a = 0.5*(omega(iom)+omega(iom-1))
     b = 0.5*(omega(iom)+omega(iom+1))
  endif
  !print *, 'After GiveLimits'
end SUBROUTINE GiveIntLimits


SUBROUTINE GetFreeEnergy(tlogG, tlogG0, tlogGD, tdens, Temperature, DM_EF, iom, omega, nomega, nbands, E, zek, DM_nemin, nume, gamma)
  USE defs,  ONLY: PI, ZERO, TWO, ZEROC, IMAG
  use param, ONLY: CHUNK
  IMPLICIT NONE
  INTEGER, intent(in) :: iom, nomega, nbands, DM_nemin, nume
  REAL*8, intent(in)  :: Temperature, omega(nomega), E(nume), gamma, DM_EF
  REAL*8, intent(out) :: tlogG, tlogG0, tlogGD, tdens
  COMPLEX*16, intent(in):: zek(nbands)
  ! external
  REAL*8  :: ferm
  ! locals
  INTEGER :: i
  REAL*8  :: a, b, E_LDA, fer
  COMPLEX*16 :: clgdmft, clglda, ddens
  !
  CALL GiveIntLimits(a,b,omega,iom, nomega)
  
  fer = ferm(omega(iom)/Temperature)
  tlogG  = 0.0
  tlogG0 = 0.0
  tdens = 0.0
  tlogGD = 0.0
  
  !$OMP PARALLEL DO PRIVATE(E_LDA,clgdmft,clglda,ddens) SCHEDULE(STATIC,CHUNK) REDUCTION(+:tlogG,tlogG0,tdens,tlogGD)
  do i=1,nbands
     E_LDA=E(i+DM_nemin-1)
     
     clgdmft = (b-zek(i)+DM_EF)*log(-b+zek(i)-DM_EF-gamma*IMAG)-(a-zek(i)+DM_EF)*log(-a+zek(i)-DM_EF-gamma*IMAG)
     clglda  = (b- E_LDA+DM_EF)*log(-b+ E_LDA-DM_EF-gamma*IMAG)-(a- E_LDA+DM_EF)*log(-a+ E_LDA-DM_EF-gamma*IMAG)
     ddens = log(-b+zek(i)-DM_EF-gamma*IMAG)-log(-a+zek(i)-DM_EF-gamma*IMAG)
     
     tlogG  = tlogG  + fer*dimag(clgdmft-clglda)/PI
     tlogG0 = tlogG0 + fer*dimag(clglda)/PI
     tdens  = tdens  - fer*dimag(ddens)/PI
     tlogGD = tlogGD + fer*dimag(clgdmft)/PI
  enddo
  !$OMP END PARALLEL DO
end SUBROUTINE GetFreeEnergy

SUBROUTINE GetLocalDensityMatrix(DM, DMFTU, Aweight, vnorm1, iorbital, iSx, nbands, maxdim2, norbitals, lmaxp)
  USE dmfts, ONLY: natom, nl, cix, ll, iso, maxdim, ncix!, lmaxp
  USE param, ONLY: 
  IMPLICIT NONE
  COMPLEX*16, intent(out):: DM(maxdim,maxdim,ncix)
  COMPLEX*16, intent(in) :: DMFTU(nbands,maxdim2,norbitals)
  COMPLEX*16, intent(in) :: Aweight(nbands,nbands)
  INTEGER, intent(in)    :: lmaxp, iorbital(natom,lmaxp+1), iSx(maxdim2, norbitals), nbands, maxdim2, norbitals
  REAL*8, intent(in)     :: vnorm1
  !                                                                                                                                                                                                                                                                                           
  INTEGER :: icase, jcase, l1case, l2case, icix, l1, l2, nind1, nind2, iorb1, iorb2, ind1, ind2!, i
  !INTEGER, allocatable :: cind(:), cini(:)
  COMPLEX*16, allocatable :: tmp(:,:), tdens(:,:)
  complex*16 :: cone, czero
  cone  = cmplx(1.d0, 0.d0, 8)
  czero = cmplx(0.d0, 0.d0, 8)
  DM=0
  allocate( tmp(maxdim2,nbands), tdens(maxdim2,maxdim2) )
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
              call zgemm('C','N', nind1, nbands, nbands, cone, DMFTU(:,:,iorb1), nbands, Aweight(:,:),nbands, czero, tmp(:,:),maxdim2)
              call zgemm('N','N', nind1, nind2, nbands, cone, tmp,maxdim2, DMFTU(:,:,iorb2),nbands, czero, tdens,maxdim2)
              do ind1=1,nind1
                 do ind2=1,nind2
                    DM( iSx(ind1,iorb1), iSx(ind2,iorb2), icix ) = tdens(ind1,ind2)*vnorm1
                 enddo
              enddo
           enddo
        ENDDO
     enddo
  ENDDO
  deallocate( tmp, tdens )

END SUBROUTINE GetLocalDensityMatrix




SUBROUTINE PrintLocalDensityMatrix(DM, fh_DM, s_oo, cixdim, PrintEntire)
  USE dmfts, ONLY: maxdim, ncix, Sigind, maxsize!, cix
  !USE dmfts, ONLY: lmaxp, natom
  IMPLICIT NONE
  COMPLEX*16, intent(inout):: DM(maxdim,maxdim,ncix)
  INTEGER, intent(in)     :: fh_DM, cixdim(ncix)
  COMPLEX*16, intent(in)  :: s_oo(maxsize,ncix)
  LOGICAL, intent(in)     :: PrintEntire
  !
  REAL*8,PARAMETER :: Ry2eV= 13.60569253d0 
  !                                                                                                                                                                                                                                                                                           
  INTEGER :: icix, cixdm, cixdms, ip, iq, it
  INTEGER, allocatable :: cind(:), cini(:)
  REAL*8, allocatable :: soo(:,:), ncorr(:,:), tsg(:,:)
  REAL*8 :: nf, sg, sgtot

  !print *, 'cixdim=', cixdim
  !print *, 's_oo=', s_oo
  sgtot=0
  WRITE(fh_DM,*)
  do icix=1,ncix
     cixdm = cixdim(icix)

     allocate( cind(cixdm) )
     ! If Sigind(i,i)=0, we eliminate the i-th column and rown, because such orbital should be treated as non-correlated                                                                                                                                                                      
     cixdms=0 ! real dimension, excluding states which must be projected out, because they are treated as non-correlated                                                                                                                                                                      
     cind=0   ! In most cases, cind(i)=i, however, when some states are projected out, cind points to smaller block                                                                                                                                                                           
     DO ip=1,cixdm
        it = Sigind(ip,ip,icix)
        if (it.gt.0) then
           cixdms = cixdms + 1
           cind(ip) = cixdms
        endif
     ENDDO
     allocate( cini(cixdms))
     do ip=1,cixdm
        if (cind(ip).gt.0) cini(cind(ip))=ip
     enddo

     allocate( soo(cixdms,cixdms), ncorr(cixdms,cixdms), tsg(cixdms,cixdms) )

     ncorr=0  ! Density matrix                                                                                                                                                                                                                                                                
     soo=0    ! Sigma(inf)                                                                                                                                                                                                                                                                    
     nf=0     ! total nd                                                                                                                                                                                                                                                                      
     DO ip=1,cixdms
        DO iq=1,cixdms
           ncorr(ip,iq) = dble(DM(cini(ip),cini(iq),icix))
           it = Sigind(cini(ip),cini(iq),icix)
           if (it.gt.0) soo(ip,iq) = real(s_oo(it,icix))
        ENDDO
        nf = nf + ncorr(ip,ip)
     ENDDO
     tsg=matmul(soo,ncorr)   ! Sigma(int)*G_loc                                                                                                                                                                                                                                               

     sg=0
     do ip=1,cixdms
        sg = sg + tsg(ip,ip)
     enddo
     sgtot = sgtot + sg

     WRITE(fh_DM,'(A,f14.8,I3,2x,I3,3x,f14.8,2x,A)') ':NCOR', nf, icix, cixdms, sg*Ry2eV, '# nf, icix, dimension Tr((s_oo-Edc)*G)'
     DO ip=1,cixdms
        DO iq=1,cixdms
           WRITE(fh_DM,'(f14.8,3x)',advance='no') ncorr(ip,iq)
        ENDDO
        WRITE(fh_DM,*)
     ENDDO
     deallocate( soo, ncorr, tsg )
     deallocate( cind, cini )
  enddo
  WRITE(fh_DM,*)

  if (PrintEntire) then
     do icix=1,ncix
        WRITE(fh_DM,'(A,1x,I3,3x,A)') ':DM ', icix, 'Density matrix'
        DO ip=1,cixdm
           DO iq=1,cixdm
              WRITE(fh_DM,'(2f14.8,3x)',advance='no') DM(ip,iq,icix)
           ENDDO
           WRITE(fh_DM,*)
        ENDDO
     enddo
  endif

END SUBROUTINE PrintLocalDensityMatrix

SUBROUTINE GetAEweight(wEpsw, Aweight, AEweight, zw2, nbands, nbands_dft, nemin, DM_nemin, DM_nemaxx)
  USE xa,  ONLY: E, weight
  IMPLICIT NONE
  COMPLEX*16, intent(out):: wEpsw(nbands_dft,nbands_dft)
  COMPLEX*16, intent(in) :: Aweight(nbands,nbands), AEweight(nbands,nbands)
  REAL*8, intent(in)     :: zw2(nbands)
  INTEGER, intent(in)    :: nbands, nbands_dft, nemin, DM_nemin, DM_nemaxx
  ! locals
  COMPLEX*16 :: tmp(nbands,nbands), tmp2(nbands,nbands)!, tmp3(nbands,nbands)
  !REAL*8  :: sqt(nbands)!, diff
  INTEGER :: i, nbandsx
  complex*16 :: cone, czero
  cone  = cmplx(1.d0, 0.d0, 8)
  czero = cmplx(0.d0, 0.d0, 8)
  !
  ! Here we calculate Eps, which is defined by the density matrix :
  !                rho=\sum_iw A^R 1/(iw + mu - eps_w) A^L ,
  !  and its cousing :
  !               (rho*E)=\sum_iw A^R eps_w/(iw + mu - eps_w) A^L
  !  which is here called AEweight==(rho*E)
  !  by:
  !       rho = Aweight * zw2 * Aweight^C
  !   (rho*E) = Aweight * wEpsw * Aweight^C
  ! Both rho and (rho*E) are Hermitian, hence Eps is Hermitian too
  ! 
  call zgemm('C','N', nbands, nbands, nbands, cone, Aweight, nbands, AEweight, nbands, czero, tmp2, nbands)
  call zgemm('N','N', nbands, nbands, nbands, cone, tmp2, nbands, Aweight, nbands, czero, tmp, nbands)

  nbandsx = DM_nemaxx-DM_nemin+1
  
  ! For the bands at which the density matrix is too small (empty bands), we just used DFT energies and weights
  wEpsw(:,:)=0.d0
  do i=1,nbands_dft
     wEpsw(i,i) = E(i+nemin-1)*weight(i+nemin-1)
  enddo
  ! Now replace the occupied and partially occupied bands with DMFT formula
  wEpsw((1+DM_nemin-nemin):(1+DM_nemaxx-nemin),(1+DM_nemin-nemin):(1+DM_nemaxx-nemin)) = tmp(1:nbandsx,1:nbandsx)
        
END SUBROUTINE GetAEweight
