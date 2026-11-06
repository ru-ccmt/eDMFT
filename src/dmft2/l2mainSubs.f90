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

SUBROUTINE Build_DMFT_Projector(DMFTU, cfX, Rspin, iorbital, &
     norbitals, nipc, n0, nnlo, nbands, cix_orb, nindo, DM_nemin, &
     DM_nemax, maxdim2, lmaxp, CalcStrain)
  USE w_atpar
  USE dmfts, ONLY: iso, projector, natom, nl, ll, cix, iatom, shft, isort, crotloc!, lmaxp, Sigind, Qcomplex, ncix, maxdim, maxsize, DM_Emin, DM_Emax, csize, CF
  USE defs, ONLY: PI
  USE param, ONLY: lomax, nloat, iblock, nmat, nume
  USE lo, ONLY: nlo, loor, ilo, lapw!, nlov, nlon, alo, blo, clo, pi12lo, pe12lo, pr12lo!, elo_store, elo, rlo
  USE atspdt,ONLY: P, DP, PE, DPE!, PEI!, e_store, el
  USE xa, ONLY: fj, dfj
  USE xa3, ONLY: bk3, As, k3
  USE structure, ONLY: rotij, tauij, BR1, POS, RMT, VOL, rotloc !, mult
  USE p_project, ONLY: phi_jl, dphi_jl, al_ucase, dri, rix_mat, &
       drix_mat!, rix, max_lcase
  USE mod_lo, ONLY: lomain,lomain_strain
  !USE w_atpar, ONLY: ri_mat
  USE Forces, ONLY: forcea!, Qforce
  USE dmft_stress, ONLY: projector_strain
  !use com, ONLY: nat
  IMPLICIT NONE
  complex*16, intent(out)  :: DMFTU(nbands,maxdim2,norbitals,nipc)
  complex*16, intent(in)   :: cfX(maxdim2,maxdim2,norbitals,norbitals), Rspin(2,2,norbitals)
  INTEGER, intent(in)      :: nipc, n0, nnlo, DM_nemin, DM_nemax, iorbital(natom,lmaxp+1), maxdim2, nbands, cix_orb(norbitals), nindo(norbitals), norbitals, lmaxp
  LOGICAL, intent(in)      :: CalcStrain
  ! locals
  COMPLEX*16, allocatable :: alm(:,:),blm(:,:),clm(:,:,:)
  complex*16 :: h_yl((lmaxp+1)*(lmaxp+1),iblock), h_alyl(2*lmaxp+1,iblock), h_blyl(2*lmaxp+1,iblock)
  INTEGER   :: icase, lcase, latom, jatom, lfirst, is, ii, i3, l, &
       M, index, index1, ibb, lda, ldb, ldc, jlo, jlop, iorb, ind, &
       i, nind, icix, num1, nm1, m1, iorb1, nind1, nind2, iorb2, &
       i_h_k, iq, ia, ib, igamma, iv, it, nchan, strain0
  REAL*8    :: BKROT2(3), BKROT3(3), BKRLOC(3), sm1, sm2, akg
  REAL*8    :: ARG1, ARG2, ARGT2, rmt2, h_al(iblock), h_bl(iblock), &
       dff, TWOPI, jprime, jsecond, dal, dbl, radial, dq
  COMPLEX*16:: FAC, PHSHEL, CFAC, cc, dd(3), sm2c, dy
  COMPLEX*16:: YL((lmaxp+2)*(lmaxp+2))
  COMPLEX*16, ALLOCATABLE :: URx(:,:,:,:), tmp(:,:), alm0(:,:)
  COMPLEX*16, PARAMETER:: IMAG=   (0.0D0,1.0D0)
  COMPLEX*16 :: tuu, tud, tdu, tdd
  COMPLEX*16, allocatable :: Uu(:,:), Ud(:,:)
  INTEGER :: no1, isize_, iind, lmmax, ip, lmaxplmaxp, lomaxlomax
  INTEGER :: strain_a(6),strain_b(6)
  LOGICAL :: nonzero_interst
  ! For p_interstitial
  COMPLEX*16, allocatable :: al_interstitial(:,:,:,:), &
       sal_interstitial(:,:,:,:), h_interstitial(:,:), &
       sh_interstitial(:,:,:)
  ! For Force
  REAL*8, allocatable     :: h_k(:,:)
  complex*16, allocatable :: h_ablyl_hk(:,:), aalm(:,:,:), bblm(:,:,:), cclm(:,:,:,:)
  complex*16, allocatable :: salm(:,:,:),sblm(:,:,:),sclm(:,:,:,:), &
       dh_yl(:,:,:),dh_alyl(:,:,:),dh_blyl(:,:,:)
  complex*16 :: ca_st_lm(3,-1:1,-1:1)
  LOGICAL :: nonzero_shft, Qforce_j
  REAL*8  :: rotloc_x_BR1(3,3), rotloc_x_BR1_x_rotij(3,3), crotloc_x_BR1(3,3)
  complex*16 :: czero, cone
  
  TWOPI=2.0*PI
  lmaxplmaxp = (lmaxp+1)*(lmaxp+1)
  lomaxlomax = (lomax+1)*(lomax+1)
  strain_a=(/1,2,3,2,1,1/)
  strain_b=(/1,2,3,3,3,2/)
  nchan=nipc
  strain0=0
  if (CalcStrain) then
     strain0=nipc
     nchan=nipc+6
     if (allocated(projector_strain)) deallocate(projector_strain)
     allocate(projector_strain(nbands,maxdim2,norbitals,6))
     projector_strain=0.d0
  endif
  cone  = cmplx(1.d0,0.d0, 8)
  czero = cmplx(0.d0,0.d0, 8)

  allocate( alm(lmaxplmaxp,nume), blm(lmaxplmaxp,nume), clm(lomaxlomax, nume, nloat) )

  if (nipc.eq.4) then
     allocate( aalm(lmaxplmaxp,nbands,3), bblm(lmaxplmaxp,nbands,3), cclm(lomaxlomax,nbands,nloat,3) )
  else
     allocate( aalm(1,1,1), bblm(1,1,1), cclm(1,1,1,1) ) ! the arrays need to be always allocated
  endif
  if (CalcStrain) then
     allocate(salm(lmaxplmaxp,nbands,6),sblm(lmaxplmaxp,nbands,6), &
          sclm(lomaxlomax,nbands,nloat,6))
     allocate(dh_yl(lmaxplmaxp,iblock,6), &
          dh_alyl(2*lmaxp+1,iblock,6), &
          dh_blyl(2*lmaxp+1,iblock,6))
  else
     allocate(salm(1,1,1),sblm(1,1,1),sclm(1,1,1,1))
     allocate(dh_yl(1,1,1),dh_alyl(1,1,1),dh_blyl(1,1,1))
  endif
  allocate( URx(DM_nemax-DM_nemin+1,maxdim2,norbitals,nchan) )
  
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
           if (CalcStrain) then
              allocate(sh_interstitial(2*lmaxp+1,iblock,6))
              allocate(sal_interstitial(2*lmaxp+1,nbands,nl(icase),6))
           endif
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
        if (CalcStrain) then
           salm=0.d0; sblm=0.d0; sclm=0.d0
        endif
        ALLOCATE( alm0(2*lmaxp+1,nbands) )
        if (abs(projector).ge.5 .and. nonzero_interst) then
           al_interstitial(:,:,:,:)=0
           if (CalcStrain) sal_interstitial=0
        endif
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
           if (CalcStrain) dh_yl=0

           !$OMP PARALLEL DO PRIVATE(bkrot2,bkrot3,bkrloc,akg,yl,arg1,arg2,argt2,phshel,lcase,l,m,index,index1,ca_st_lm,dy,iq,ia,ib,igamma,iv,it)&
           !$OMP& SHARED(h_yl,h_k,dh_yl)&
           !$OMP& SCHEDULE(STATIC)
           do i=ii,min(ii+iblock-1,n0-nnlo)  
              !---------  rotates ylm(k+K) to ylm(k'+K) where k' is in star of irreducible k. ------------!
              ! bk3(:,i) !-----  reciprocal vector and irreducible vector: G=K+k ----!
              ! BKROT2 = R_n.R_a.(k+K), transformation from the first atom to an equivalent atom 
              BKROT2 = matmul(rotij(:,:,latom), bk3(:,I))  ! apply one of the symmetry operations to BKROT=R.(k+K)
              BKROT3 = matmul(BR1,BKROT2)
              akg=sqrt(dot_product(BKROT3,BKROT3))
              ! !!---- BR1 transforms integer reciprocal lattice vectors, as given in the VECTORLIST of LAPW1, into cartesian system ----!
              ! !! BKROT3 = R_n.R_a.(k+K), but in cartesian coordinate system
              ! !BKROT3 = matmul(BR1, BKROT2)
              ! !!---- BKRLOC = Rotloc.R_g.(k+K),  is rotates Rotloc.R.(k+K) . Rotation Rotloc entered by user
              ! !! Here we use crotloc and not rotloc, hence we use the transformation case.indmfl file, input by the user.
              ! !BKRLOC = matmul(crotloc(:,:,icase),BKROT3)
              BKRLOC = matmul(crotloc(:,:,icase),BKROT3)
              !---- YLM = Y_{L}(Rotloc.R_g.(k+K))
              ! One extra angular momentum is required by the analytic
              ! derivative of Y_lm.  The ordinary coefficients use only the
              ! first (lmaxp+1)^2 entries, so this does not alter P itself.
              CALL YLM (BKRLOC,lmaxp+1,YL,(lmaxp+2)*(lmaxp+2))
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
                 ! Keep this scalar.  gfortran 15.2 miscompiles the
                 ! vectorized DCONJG array expression at -O2, producing
                 ! incorrect muffin-tin charges (and mixer NEC01 errors).
                 index=l*l
!GCC$ novector
!DIR$ NOVECTOR
                 do m=1,2*l+1
                    index=index+1
                    h_yl(index,i-ii+1)=dconjg(yl(index))*phshel
                    if (CalcStrain) then
                       call c_factor(l,m-l-1,ca_st_lm)
                       do iq=1,6
                          ia=strain_a(iq); ib=strain_b(iq)
                          dy=czero
                          if (akg.gt.1.d-14) then
                             do igamma=1,3
                                do iv=-1,1,2
                                   if (l+iv.lt.0) cycle
                                   do it=-1,1
                                      if (m-l-1+it.lt.-(l+iv) .or. &
                                           m-l-1+it.gt.(l+iv)) cycle
                                      index1=(l+iv)*(l+iv+1)+(m-l-1+it)+1
                                      ! Fixed Cartesian correlated axes.  Under
                                      ! a symmetric homogeneous strain,
                                      ! dK=-eps*K and q=C*K, hence
                                      ! dq_gamma=-1/2(C_gamma,a K_b+
                                      ! C_gamma,b K_a).  c_factor supplies the
                                      ! angular part |K| grad_K Y_lm.  The
                                      ! division by |K| is essential because a
                                      ! spherical harmonic depends on K-hat,
                                      ! not on the magnitude of K.  This is the
                                      ! same 1/sqrt(KPLUSG) factor used by the
                                      ! WIEN2k stress kernels.
                                      dy=dy-0.5d0* &
                                           conjg(ca_st_lm(igamma,iv,it))* &
                                           (crotloc(igamma,ia,icase)* &
                                            BKROT3(ib)+ &
                                            crotloc(igamma,ib,icase)* &
                                            BKROT3(ia))* &
                                           conjg(yl(index1))*phshel/akg
                                   enddo
                                enddo
                             enddo
                          endif
                          dh_yl(index,i-ii+1,iq)=dy
                       enddo
                    endif
                 enddo
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

              if (CalcStrain) then
                 dh_alyl=0.d0; dh_blyl=0.d0
                 do i=ii,min(ii+iblock-1,n0-nnlo)
                    i3=i-ii+1
                    BKROT2=matmul(rotij(:,:,latom),bk3(:,i))
                    BKROT3=matmul(BR1,BKROT2)
                    akg=sqrt(dot_product(BKROT3,BKROT3))
                    if (akg.gt.1.d-14) then
                       jprime=dfj(l,i,jatom)/akg
                       jsecond=((dble(l*(l+1))/(akg*rmt(jatom))**2)-1.d0)* &
                            fj(l,i,jatom)-2.d0*jprime/(akg*rmt(jatom))
                    else
                       jprime=0.d0; jsecond=0.d0
                    endif
                    if (lapw(l)) then
                       ! Minus the radial derivative dh/d|K|.  Multiplying
                       ! by K_a K_b/|K| below gives dh/d epsilon_ab because
                       ! d|K|/d epsilon_ab=-K_a K_b/|K|.
                       dal=rmt(jatom)*DPE(l)*jprime-PE(l)*jprime- &
                            PE(l)*akg*rmt(jatom)*jsecond
                       dbl=-rmt(jatom)*DP(l)*jprime+P(l)*jprime+ &
                            P(l)*akg*rmt(jatom)*jsecond
                    else
                       dal=0.d0; dbl=0.d0
                       if (akg.gt.1.d-14) &
                            dal=-jprime/(rmt(jatom)*P(l))
                    endif
                    do iq=1,6
                       ia=strain_a(iq); ib=strain_b(iq)
                       radial=0.d0
                       if (akg.gt.1.d-14) &
                            radial=BKROT3(ia)*BKROT3(ib)/akg
                       do m=1,2*l+1
                          index=l*l+m
                          dh_alyl(m,i3,iq)=dh_yl(index,i3,iq)*h_al(i3)+ &
                               h_yl(index,i3)*radial*dal
                          dh_blyl(m,i3,iq)=dh_yl(index,i3,iq)*h_bl(i3)+ &
                               h_yl(index,i3)*radial*dbl
                       enddo
                    enddo
                 enddo
              endif
              
              if (abs(projector).ge.5 .and. nonzero_interst) then
                 h_interstitial(:,:)=0
                 if (CalcStrain) sh_interstitial=0
                 do i=ii,min(ii+iblock-1,n0-nnlo)
                    !----    h_interstitial =Y*_{lm}(R(k+K))*exp(i*(k+K)*r_latom)  <P_phi|j_l>/Rmt^2
                    h_interstitial(1:2*l+1,i-ii+1) = h_yl(l*l+1:(l+1)**2,i-ii+1)*phi_jl(i,iind)*rmt2
                    if (CalcStrain) then
                       BKROT2=matmul(rotij(:,:,latom),bk3(:,i))
                       BKROT3=matmul(BR1,BKROT2)
                       akg=sqrt(dot_product(BKROT3,BKROT3))
                       do iq=1,6
                          ia=strain_a(iq); ib=strain_b(iq)
                          radial=0.d0
                          if (akg.gt.1.d-14) radial=-BKROT3(ia)* &
                               BKROT3(ib)/akg*dphi_jl(i,iind)
                          do m=1,2*l+1
                             index=l*l+m
                             sh_interstitial(m,i-ii+1,iq)= &
                                  (dh_yl(index,i-ii+1,iq)*phi_jl(i,iind)+ &
                                   h_yl(index,i-ii+1)*radial)*rmt2
                          enddo
                       enddo
                    endif
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

              if (CalcStrain) then
                 do iq=1,6
                    call zgemm('N','N',isize_,nbands,ibb,cone, &
                         dh_alyl(:,:,iq),lda,As(ii,DM_nemin,is),ldb, &
                         czero,alm0,ldc)
                    do m=1,2*l+1
                       salm(l*l+m,:,iq)=salm(l*l+m,:,iq)+alm0(m,:)
                    enddo
                    call zgemm('N','N',isize_,nbands,ibb,cone, &
                         dh_blyl(:,:,iq),lda,As(ii,DM_nemin,is),ldb, &
                         czero,alm0,ldc)
                    do m=1,2*l+1
                       sblm(l*l+m,:,iq)=sblm(l*l+m,:,iq)+alm0(m,:)
                    enddo
                    if (abs(projector).ge.5 .and. nonzero_interst) &
                         call zgemm('N','N',isize_,nbands,ibb,cone, &
                         sh_interstitial(:,:,iq),lda, &
                         As(ii,DM_nemin,is),ldb,cone, &
                         sal_interstitial(:,:,lcase,iq),ldc)
                 enddo
              endif
              
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
        if (nlo.ne.0) then
           call lomain(crotloc(:,:,icase),is,DM_nemin,DM_nemax,lfirst, &
                latom,n0,jatom,alm,blm,clm,Qforce_j,aalm,bblm,cclm,lmaxp)
           if (CalcStrain) call lomain_strain(crotloc(:,:,icase),is, &
                DM_nemin,DM_nemax,lfirst,latom,jatom,salm,sblm,sclm, &
                lmaxp)
        endif

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
        if (CalcStrain) then
           do lcase=1,nl(icase)
              l=ll(icase,lcase)
              cfac=fac*imag**l
              salm(l*l+1:(l+1)**2,:,1:6)= &
                   salm(l*l+1:(l+1)**2,:,1:6)*cfac
              sblm(l*l+1:(l+1)**2,:,1:6)= &
                   sblm(l*l+1:(l+1)**2,:,1:6)*cfac
              if (l.le.lomax) sclm(l*l+1:(l+1)**2,:,1:ilo(l),1:6)= &
                   sclm(l*l+1:(l+1)**2,:,1:ilo(l),1:6)*cfac
              if (abs(projector).ge.5 .and. nonzero_interst) &
                   sal_interstitial(:,:,lcase,1:6)= &
                   sal_interstitial(:,:,lcase,1:6)*cfac

              ! All matching coefficients carry 1/sqrt(Omega).  At fixed
              ! muffin-tin radii, d Omega/Omega=Tr(epsilon), hence the three
              ! diagonal strain channels receive -P/2.  This term is easy to
              ! miss because it is not contained in the angular/radial
              ! derivative above; it also supplies the entire explicit
              ! derivative of a fixed radial local orbital.
              do iq=1,3
                 salm(l*l+1:(l+1)**2,:,iq)= &
                      salm(l*l+1:(l+1)**2,:,iq)-0.5d0* &
                      alm(l*l+1:(l+1)**2,DM_nemin:DM_nemax)
                 sblm(l*l+1:(l+1)**2,:,iq)= &
                      sblm(l*l+1:(l+1)**2,:,iq)-0.5d0* &
                      blm(l*l+1:(l+1)**2,DM_nemin:DM_nemax)
                 if (l.le.lomax) sclm(l*l+1:(l+1)**2,:,1:ilo(l),iq)= &
                      sclm(l*l+1:(l+1)**2,:,1:ilo(l),iq)-0.5d0* &
                      clm(l*l+1:(l+1)**2,DM_nemin:DM_nemax,1:ilo(l))
                 if (abs(projector).ge.5 .and. nonzero_interst) &
                      sal_interstitial(:,:,lcase,iq)= &
                      sal_interstitial(:,:,lcase,iq)-0.5d0* &
                      al_interstitial(:,:,lcase,1)
              enddo
           enddo
        endif
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
                    if (CalcStrain) then
                       do iq=1,6
                          cc=salm(index,nm1,iq)
                          do jlo=1,ilo(l)
                             if (loor(jlo,l)) cc=cc+ &
                                  sclm(index,nm1,jlo,iq)* &
                                  ri_mat(1,2+jlo,l,jatom)
                          enddo
                          ! For projectors 1 and 3 dff=1 and this is exact.
                          ! The historical per-band scalar normalization of
                          ! projectors 2 and 4 is retained here; global matrix
                          ! normalization (the production path) is
                          ! differentiated separately below.
                          URx(nm1,ind,iorb,strain0+iq)=cc*dff
                       enddo
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
                    if (CalcStrain) then
                       do iq=1,6
                          cc=salm(index,nm1,iq)*rix_mat(1,iind)+ &
                               sblm(index,nm1,iq)*rix_mat(2,iind)
                          do jlo=1,ilo(l)
                             if (loor(jlo,l)) cc=cc+ &
                                  sclm(index,nm1,jlo,iq)* &
                                  rix_mat(2+jlo,iind)
                          enddo
                          ! Optional projector-5 radial response.  The normal
                          ! production path has drix_mat=0.  When a diagnostic
                          ! derivative is supplied, add the missing chain-rule
                          ! branch a*d<u|phi>/deps+b*d<udot|phi>/deps+LO.
                          cc=cc+alm(index,num1)*drix_mat(1,iind,iq)+ &
                               blm(index,num1)*drix_mat(2,iind,iq)
                          do jlo=1,ilo(l)
                             if (loor(jlo,l)) cc=cc+ &
                                  clm(index,num1,jlo)* &
                                  drix_mat(2+jlo,iind,iq)
                          enddo
                          if (nonzero_interst) cc=cc+ &
                               sal_interstitial(l+m1+1,nm1,lcase,iq)
                          URx(nm1,ind,iorb,strain0+iq)=cc
                       enddo
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
        if (CalcStrain) deallocate(sh_interstitial,sal_interstitial)
     endif
  ENDDO    !----------- over all atoms  ---------!


  deallocate( alm, blm, clm )
  deallocate( aalm, bblm, cclm ) 
  deallocate(salm,sblm,sclm,dh_yl,dh_alyl,dh_blyl)
  

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
        do ip=1,nchan
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
        do ip=1,nchan
           call zgemm('N','C', nbands, nind2, nind1, cone, URx(:,:,iorb1,ip), nbands, cfX(:,:,iorb2,iorb1),maxdim2,czero, tmp, nbands)
           if (ip.le.nipc) then
              DMFTU(:,:nind2,iorb2,ip)= &
                   DMFTU(:,:nind2,iorb2,ip)+conjg(tmp(:,:nind2))
           else
              projector_strain(:,:nind2,iorb2,ip-strain0)= &
                   projector_strain(:,:nind2,iorb2,ip-strain0)+ &
                   conjg(tmp(:,:nind2))
           endif
        enddo
     ENDDO
  ENDDO


  
  DEALLOCATE( URx, tmp )

END SUBROUTINE Build_DMFT_Projector


!===============================================================================
! OPTIONAL FINITE-DIFFERENCE CHECK OF THE EXPLICIT PROJECTOR STRAIN DERIVATIVE
!===============================================================================
!
! This is a diagnostic, not part of the production stress evaluation.  It is
! activated only by a file named projector_strain_fd.in in the case directory.
!
! The analytic array projector_strain contains the explicit basis derivative
! P_,q at fixed LAPW eigenvector coefficients.  Comparing projectors obtained
! from independently strained LAPW1 calculations would not test that object:
! the eigenvectors have arbitrary phases and may rotate within degenerate
! subspaces, and they also contain the physical outside-window response.  Here
! we instead keep As fixed, perturb only BR1 and VOL, rebuild the matching and
! projector-tail integrals, and central-difference the resulting raw P.  Thus
!
!          [P(+h)-P(-h)]/(2h)  ->  P_,q  (fixed eigenvectors).
!
! For off-diagonal Voigt channels q=yz,xz,xy, the scalar parameter h denotes
! the symmetric strain with epsilon_ab=epsilon_ba=h/2.  This is the convention
! used by Build_DMFT_Projector and by the stress kernels.
!===============================================================================
SUBROUTINE CheckProjectorStrainFiniteDifference(P0,dP,cfX,Rspin, &
     iorbital,norbitals,n0,nnlo,nbands,cix_orb,nindo,nemin,nemax, &
     maxdim2,lmaxp)
  USE dmfts,ONLY: iso,projector,natom
  USE param,ONLY: nmat,lmax2
  USE com,ONLY: nat
  USE com_mpi,ONLY: myrank
  USE structure,ONLY: BR1,VOL,RMT
  USE xa,ONLY: fj,dfj
  USE xa3,ONLY: bk3,aK
  USE p_project,ONLY: phi_jl,P_rfi,rix,kNri,n_al_ucase,dri, &
       l_al_ucase,j_al_ucase
  IMPLICIT NONE
  INTEGER,intent(in) :: norbitals,n0,nnlo,nbands,nemin,nemax,maxdim2,lmaxp
  INTEGER,intent(in) :: iorbital(natom,lmaxp+1),cix_orb(norbitals), &
       nindo(norbitals)
  COMPLEX*16,intent(in) :: P0(nbands,maxdim2,norbitals), &
       dP(nbands,maxdim2,norbitals,6)
  COMPLEX*16,intent(in) :: cfX(maxdim2,maxdim2,norbitals,norbitals), &
       Rspin(2,2,norbitals)
  INTEGER,parameter :: strain_a(6)=(/1,2,3,2,1,1/), &
       strain_b(6)=(/1,2,3,3,3,2/)
  CHARACTER(len=2),parameter :: qname(6)=(/'xx','yy','zz','yz','xz','xy'/)
  REAL*8,parameter :: h=1.d-6
  LOGICAL :: enabled
  INTEGER :: iq,i,isize,jatom,Nri,iind,l,ir,fh
  REAL*8 :: BR10(3,3),VOL0,Kn(3),rx,nA,nFD,nerr,scale,ncurve
  REAL*8,allocatable :: aKR(:),jlr(:),jlp(:)
  COMPLEX*16,allocatable :: Pplus(:,:,:,:),Pminus(:,:,:,:),fdP(:,:,:)

  INQUIRE(file='projector_strain_fd.in',exist=enabled)
  if (.not.enabled) return

  isize=n0-nnlo
  BR10=BR1
  VOL0=VOL
  allocate(Pplus(nbands,maxdim2,norbitals,1), &
       Pminus(nbands,maxdim2,norbitals,1), &
       fdP(nbands,maxdim2,norbitals))
  Nri=2**kNri+1
  if (abs(projector).ge.5) allocate(aKR(Nri),jlr(Nri),jlp(Nri))

  fh=998
  open(fh,file='projector_strain_fd.dat',status='replace',action='write')
  write(fh,'(A)') '# Explicit raw-projector strain finite difference'
  write(fh,'(A,ES24.15)') '# h = ',h
  write(fh,'(A)') &
       '# q  ||analytic||  ||finite-difference||  ||difference||  relative  curvature'
  write(6,'(/,A)') 'PROJECTOR_STRAIN_FD: fixed-eigenvector raw-projector test'

  do iq=1,6
     call SetGeometry(+h,iq)
     CALL Build_DMFT_Projector(Pplus,cfX,Rspin,iorbital,norbitals,1, &
          n0,nnlo,nbands,cix_orb,nindo,nemin,nemax,maxdim2,lmaxp,.false.)
     call SetGeometry(-h,iq)
     CALL Build_DMFT_Projector(Pminus,cfX,Rspin,iorbital,norbitals,1, &
          n0,nnlo,nbands,cix_orb,nindo,nemin,nemax,maxdim2,lmaxp,.false.)

     fdP=(Pplus(:,:,:,1)-Pminus(:,:,:,1))/(2.d0*h)
     nA=sqrt(sum(abs(dP(:,:,:,iq))**2))
     nFD=sqrt(sum(abs(fdP)**2))
     nerr=sqrt(sum(abs(fdP-dP(:,:,:,iq))**2))
     scale=max(nA,nFD,1.d-30)
     ! This dimensionless second difference warns when h is too large or a
     ! non-smooth branch has entered the diagnostic.
     ncurve=sqrt(sum(abs(Pplus(:,:,:,1)+Pminus(:,:,:,1)-2.d0*P0)**2)) &
          /max(h*scale,1.d-30)
     write(fh,'(A2,1X,5ES24.15)') qname(iq),nA,nFD,nerr,nerr/scale,ncurve
     write(6,'(A,1X,A2,A,ES12.4,A,ES12.4,A,ES12.4)') &
          'PROJECTOR_STRAIN_FD:',qname(iq),'  rel.err=',nerr/scale, &
          '  ||analytic||=',nA,'  ||FD||=',nFD
  enddo

  ! No diagnostic is allowed to alter the subsequent production calculation.
  call SetGeometry(0.d0,1)
  close(fh)
  if (abs(projector).ge.5) deallocate(aKR,jlr,jlp)
  deallocate(Pplus,Pminus,fdP)

CONTAINS

  SUBROUTINE SetGeometry(delta,iq_local)
    REAL*8,intent(in) :: delta
    INTEGER,intent(in) :: iq_local
    INTEGER :: ia,ib
    REAL*8 :: Romb

    ia=strain_a(iq_local)
    ib=strain_b(iq_local)
    BR1=BR10
    if (ia.eq.ib) then
       ! To first order B(epsilon)=(I-epsilon)B(0).
       BR1(ia,:)=(1.d0-delta)*BR10(ia,:)
       VOL=VOL0*(1.d0+delta)
    else
       BR1(ia,:)=BR10(ia,:)-0.5d0*delta*BR10(ib,:)
       BR1(ib,:)=BR10(ib,:)-0.5d0*delta*BR10(ia,:)
       VOL=VOL0
    endif

    do i=1,isize
       Kn=matmul(BR1,bk3(:,i))
       aK(i)=sqrt(dot_product(Kn,Kn))
    enddo
    do jatom=1,nat
       CALL HARMON(isize,aK(:isize),lmax2,fj(:,:isize,jatom), &
            dfj(:,:isize,jatom),RMT(jatom))
    enddo

    ! Projector 5 has a finite interstitial tail.  Its radial Bessel integral
    ! depends on |k+G| and must be rebuilt together with the LAPW matching
    ! coefficients for the numerical perturbation.
    if (abs(projector).ge.5) then
       phi_jl=0.d0
       do iind=1,n_al_ucase
          if (abs(dri(iind)).le.1.d-10) cycle
          l=l_al_ucase(iind)
          jatom=j_al_ucase(iind)
          do i=1,isize
             rx=RMT(jatom)
             do ir=1,Nri
                aKR(ir)=rx*aK(i)
                rx=rx+dri(iind)
             enddo
             CALL sphbes2(l,Nri,aKR,jlr)
             jlp=jlr*rix*P_rfi(:,iind)
             phi_jl(i,iind)=Romb(jlp,Nri,dri(iind))
          enddo
       enddo
    endif
  END SUBROUTINE SetGeometry

END SUBROUTINE CheckProjectorStrainFiniteDifference


SUBROUTINE Diagonalize_DMFT_WEIGHTS(zw2, Aweight, nbands, DM_nemin, &
     DM_nemaxx, determine_DM_nemaxx, diagnostic_label, diagnostic_kpoint)
  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
  USE com_mpi, ONLY: myrank, stop_MPI
  IMPLICIT NONE
  REAL*8,     intent(out)  :: zw2(nbands)
  COMPLEX*16, intent(inout):: Aweight(nbands,nbands)
  INTEGER,    intent(in)   :: nbands, DM_nemin
  INTEGER,    intent(inout):: DM_nemaxx
  LOGICAL,    intent(in)   :: determine_DM_nemaxx
  CHARACTER(len=*), intent(in) :: diagnostic_label
  INTEGER, intent(in) :: diagnostic_kpoint
  ! locals
  INTEGER    :: i, j, n, nbandsx, num, invalid_count
  INTEGER, allocatable    :: imx(:)
  COMPLEX*16, allocatable :: zweight(:,:), znew(:), Anew(:,:), &
       Ainput(:,:), work_fallback(:)
  COMPLEX*16 :: cc
  INTEGER    :: lwork, lrwork, liwork, lwork_fallback
  complex*16, allocatable :: work(:)
  real*8,     allocatable :: rwork(:), rwork_fallback(:)
  integer,    allocatable :: iwork(:)
  integer :: info_zheevd, info_zheev
  logical :: debug, negative
  REAL*8 :: matrix_scale, hermiticity_error, relative_hermiticity_error
  !
  lwork = 2*nbands+nbands*nbands
  lrwork =  5*nbands + 2*nbands*nbands + 1
  liwork = 3+5*nbands
  ALLOCATE ( work(lwork), rwork(lrwork), iwork(liwork) )
  allocate(Ainput(nbands,nbands))

  debug = .false.

  ! ZHEEVD assumes a finite Hermitian matrix and is allowed to overwrite its
  ! input even when it fails.  The stress response matrices are deliberately
  ! low rank and can contain a large exactly degenerate null space.  Preserve
  ! a clean copy so that a divide-and-conquer failure can be retried with the
  ! more conservative QR driver ZHEEV instead of continuing with corrupted
  ! eigenvectors.
  invalid_count=0
  do j=1,nbands
     do i=1,nbands
        if (.not.IEEE_IS_FINITE(dble(Aweight(i,j))) .or. &
             .not.IEEE_IS_FINITE(aimag(Aweight(i,j)))) &
             invalid_count=invalid_count+1
     enddo
  enddo
  if (invalid_count.ne.0) then
     WRITE(*,'(A,1X,A,3(A,I8))') ':DMFT_WEIGHT_DIAG_INVALID:', &
          trim(diagnostic_label),' kpoint=',diagnostic_kpoint, &
          ' rank=',myrank,' nonfinite=',invalid_count
     FLUSH(6)
     CALL stop_MPI
     STOP 'non-finite DMFT weight matrix'
  endif

  matrix_scale=maxval(abs(Aweight))
  hermiticity_error=maxval(abs(Aweight-transpose(conjg(Aweight))))
  relative_hermiticity_error=hermiticity_error/max(matrix_scale,tiny(1.d0))
  if (relative_hermiticity_error.gt.1.d-10) then
     WRITE(*,'(A,1X,A,2(A,I8),3(A,ES14.6))') &
          ':DMFT_WEIGHT_HERMITICITY_REPAIR:',trim(diagnostic_label), &
          ' kpoint=',diagnostic_kpoint,' rank=',myrank, &
          ' scale=',matrix_scale,' error=',hermiticity_error, &
          ' relative=',relative_hermiticity_error
  endif

  ! The historical routine diagonalizes -Aweight so LAPACK's ascending
  ! eigenvalue order becomes descending physical occupations after the sign
  ! is restored below.  Explicit symmetrization also makes the result
  ! independent of whether a caller populated the upper or lower triangle.
  Ainput(:,:)=-0.5d0*(Aweight(:,:)+transpose(conjg(Aweight(:,:))))
  Aweight(:,:)=Ainput(:,:)

  if (debug) then
     ! Remember the repaired physical matrix for the reconstruction check.
     allocate( zweight(nbands,nbands) )
     zweight(:,:) = -Ainput(:,:)
  endif

  CALL ZHEEVD('V','U', nbands, Aweight, nbands, zw2, work, lwork, &
       rwork, lrwork, iwork, liwork, info_zheevd)
  if (info_zheevd.ne.0) then
     WRITE(*,'(A,1X,A,3(A,I8),3(A,ES14.6))') &
          ':DMFT_WEIGHT_ZHEEVD_FAILURE:',trim(diagnostic_label), &
          ' kpoint=',diagnostic_kpoint,' rank=',myrank, &
          ' info=',info_zheevd,' scale=',matrix_scale, &
          ' herm_error=',hermiticity_error, &
          ' herm_relative=',relative_hermiticity_error

     ! ZHEEV uses the QR algorithm.  It is slower asymptotically, but these
     ! matrices contain only the stored LAPW/DMFT bands, so the retry is tiny
     ! compared with constructing the stress kernels.  Most importantly it
     ! is substantially less sensitive to the large degenerate null spaces
     ! of complementary-window response matrices on Accelerate/LAPACK.
     lwork_fallback=max(1,2*nbands-1)
     allocate(work_fallback(lwork_fallback))
     allocate(rwork_fallback(max(1,3*nbands-2)))
     Aweight(:,:)=Ainput(:,:)
     CALL ZHEEV('V','U',nbands,Aweight,nbands,zw2,work_fallback, &
          lwork_fallback,rwork_fallback,info_zheev)
     if (info_zheev.eq.0) then
        WRITE(*,'(A,1X,A,2(A,I8),A,I8)') &
             ':DMFT_WEIGHT_ZHEEV_FALLBACK_OK:',trim(diagnostic_label), &
             ' kpoint=',diagnostic_kpoint,' rank=',myrank, &
             ' zheevd_info=',info_zheevd
     else
        WRITE(*,'(A,1X,A,4(A,I8))') &
             ':DMFT_WEIGHT_DIAGONALIZATION_FATAL:',trim(diagnostic_label), &
             ' kpoint=',diagnostic_kpoint,' rank=',myrank, &
             ' zheevd_info=',info_zheevd,' zheev_info=',info_zheev
        FLUSH(6)
        CALL stop_MPI
        STOP 'DMFT weight diagonalization failed'
     endif
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
  DEALLOCATE ( work, rwork, iwork, Ainput )
  if (allocated(work_fallback)) deallocate(work_fallback)
  if (allocated(rwork_fallback)) deallocate(rwork_fallback)
  
  deallocate( imx )
  deallocate( znew, Anew )
  if (debug) deallocate( zweight)
END SUBROUTINE Diagonalize_DMFT_WEIGHTS


!---------------------------------------------
! To compute overlap needed by renormalization
!---------------------------------------------
SUBROUTINE cmp_overlap(Olapm,SOlapm,max_nbands,nnlo,norbitals, &
     natom,maxdim2,emin,emax,iorbital,cix_orb,iSx,nindo,noccur, &
     cixdim,cfX,Rspin,sumw,lmaxp,time_dmf0,time_dmf0w,SIMPLE, &
     dXrel,dXstress,CalcDerivative,CalcStrain,KeepNormalization)
  USE xa,    ONLY: fj, dfj
  USE xa3,   ONLY: bk3, bk3lo, aK, k3, k3lo, As, As_lo
  USE param, ONLY: nume, nmat, nloat, iblock, lmax2
  USE structure, ONLY: rmt, BR1
  USE com,   ONLY: nat
  USE defs,  ONLY: PI, IMAG
  USE dmfts,   ONLY: DM_Emin, DM_Emax, maxsize, iso, ncix, Sigind, csize, maxdim, projector!, lmaxp, nl, cix, DM_EF, ll, 
  USE com_mpi,ONLY: myrank, master, FilenameMPI, Reduce_MPI, FindMax_MPI, reduce_MPI0,pr_proc,vector_para,vectors,nvector,fvectors,Qprint
  USE p_project,ONLY: phi_jl,dphi_jl,rix,P_rfi,kNri,n_al_ucase, &
       dri,l_al_ucase,j_al_ucase
  USE dmft_stress,ONLY: projector_strain
  !USE sym2,  ONLY: iord
  !USE MProjector, ONLY: Build_DMFT_Projector2
  IMPLICIT NONE
  COMPLEX*16, intent(out):: Olapm(maxdim,maxdim,ncix), SOlapm(maxdim,maxdim,ncix)
  COMPLEX*16, intent(out):: dXrel(maxdim,maxdim,ncix,3)
  COMPLEX*16, intent(out):: dXstress(maxdim,maxdim,ncix,6)
  INTEGER,    intent(out):: max_nbands
  REAL*8,     intent(out):: time_dmf0, time_dmf0w
  INTEGER,    intent(in) :: nnlo, norbitals, natom, maxdim2, lmaxp
  REAL*8,     intent(in) :: emin, emax, sumw
  INTEGER,    intent(in) :: iorbital(natom,lmaxp+1), cix_orb(norbitals), iSx(maxdim2, norbitals), nindo(norbitals), noccur(maxsize,ncix), cixdim(ncix)
  complex*16, intent(in) :: cfX(maxdim2,maxdim2,norbitals,norbitals), Rspin(2,2,norbitals)
  LOGICAL, intent(in)    :: SIMPLE,CalcDerivative,CalcStrain, &
       KeepNormalization
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
  complex*16, ALLOCATABLE :: DMFTU(:,:,:,:), Olapmk(:,:,:), olp(:,:), olp2(:,:)
  complex*16, ALLOCATABLE :: dOlapm(:,:,:,:),dOlapmk(:,:,:,:), &
       dOlapStress(:,:,:,:),dOlapStressK(:,:,:,:)
  complex*16, ALLOCATABLE :: Xsaved(:,:,:)
  INTEGER    :: ikp, nemin, nemax, DM_nemin, DM_nemax, n0, idir
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
  REAL*8, allocatable     :: aKR(:),jlr(:),jlrp(:),jlp(:)
  !
  INTEGER, allocatable :: cind(:), cini(:)
  INTEGER :: cixdms
  complex*16 :: czero, cone

  cone  = cmplx(1.d0, 0.d0, 8)
  czero = cmplx(0.d0, 0.d0, 8)
  allocate( Olapmk(maxdim,maxdim,ncix), olp(maxdim2,maxdim2) )
  Olapm=0
  dXrel=0
  dXstress=0
  if (KeepNormalization) then
     allocate(Xsaved(maxdim,maxdim,ncix))
     Xsaved=SOlapm
  endif
  if (CalcDerivative .or. CalcStrain) then
     allocate(olp2(maxdim2,maxdim2))
  endif
  if (CalcDerivative) then
     allocate(dOlapm(maxdim,maxdim,ncix,3), dOlapmk(maxdim,maxdim,ncix,3))
     dOlapm=0
  endif
  if (CalcStrain) then
     allocate(dOlapStress(maxdim,maxdim,ncix,6), &
          dOlapStressK(maxdim,maxdim,ncix,6))
     dOlapStress=0
  endif

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
        if (abs(projector).ge.5) then
           allocate(phi_jl(nmat,n_al_ucase),dphi_jl(nmat,n_al_ucase))
           
           Nri=2**kNri+1    ! Number of radial points in the interstitials
           phi_jl=0; dphi_jl=0
           allocate(aKR(Nri),jlr(Nri),jlrp(Nri),jlp(Nri))
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
                    CALL sphbes2(l,Nri,aKR,jlr)
                    CALL sphbes2(l+1,Nri,aKR,jlrp)
                    jlp(:)=jlr(:)*rix(:)*P_rfi(:,iind)
                    phi_jl(i,iind) = romb(jlp, Nri, dri(iind))  ! Integration over r on rix mesh
                    if (aK(i).gt.1.d-14) then
                       do ir=1,Nri
                          jlp(ir)=(dble(l)*jlr(ir)/aKR(ir)-jlrp(ir))* &
                               rix(ir)**2*P_rfi(ir,iind)
                       enddo
                       dphi_jl(i,iind)=romb(jlp,Nri,dri(iind))
                    endif
                 ENDDO
              endif
           ENDDO
           deallocate(aKR,jlr,jlrp,jlp)
        endif
        ! New

        ! The ordinary overlap construction needs only projector channel 1.
        ! For a force calculation channels 2:4 contain D_A, with U_A=-iD_A.
        ! They let us build S_A analytically while traversing the same k-points
        ! and using exactly the same raw projector definition as for S.
        nipc=1
        if (CalcDerivative) nipc=4
        allocate( DMFTU(nbands,maxdim2,norbitals,nipc) )

        CALL Build_DMFT_Projector(DMFTU,cfX,Rspin,iorbital,norbitals, &
             nipc,n0,nnlo,nbands,cix_orb,nindo,DM_nemin,DM_nemax, &
             maxdim2,lmaxp,CalcStrain)
        
        Olapmk=0
        if (CalcDerivative) dOlapmk=0
        if (CalcStrain) dOlapStressK=0
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
              if (CalcDerivative) then
                 do idir=1,3
                    ! olp  = D_A^+ U and olp2 = U^+ D_A.  Since U_A=-iD_A,
                    !
                    !   S_A = U_A^+U + U^+U_A
                    !       = i(D_A^+U-U^+D_A).
                    !
                    ! This combination is Hermitian up to roundoff, as a
                    ! derivative of the Hermitian overlap S must be.
                    call zgemm('C','N', nind1, nind2, nbands, cone, DMFTU(:,:,iorb1,idir+1),nbands, DMFTU(:,:,iorb2,1),nbands, czero, olp,maxdim2)
                    call zgemm('C','N', nind1, nind2, nbands, cone, DMFTU(:,:,iorb1,1),nbands, DMFTU(:,:,iorb2,idir+1),nbands, czero, olp2,maxdim2)
                    do ind1=1,nind1
                       do ind2=1,nind2
                          ! Build dS=U_A^+ U+U^+ U_A.  The force channel D
                          ! is defined by U_A=-iD in Build_DMFT_Projector.
                          dOlapmk(iSx(ind1,iorb1),iSx(ind2,iorb2),icix,idir) = IMAG*(olp(ind1,ind2)-olp2(ind1,ind2))
                       enddo
                    enddo
                 enddo
              endif
              if (CalcStrain) then
                 do idir=1,6
                    call zgemm('C','N',nind1,nind2,nbands,cone, &
                         projector_strain(:,:,iorb1,idir),nbands, &
                         DMFTU(:,:,iorb2,1),nbands,czero,olp,maxdim2)
                    call zgemm('C','N',nind1,nind2,nbands,cone, &
                         DMFTU(:,:,iorb1,1),nbands, &
                         projector_strain(:,:,iorb2,idir),nbands, &
                         czero,olp2,maxdim2)
                    do ind1=1,nind1
                       do ind2=1,nind2
                          dOlapStressK(iSx(ind1,iorb1), &
                               iSx(ind2,iorb2),icix,idir)= &
                               olp(ind1,ind2)+olp2(ind1,ind2)
                       enddo
                    enddo
                 enddo
              endif
           enddo
        ENDDO
           
        wg = (wgh0/sumw)
        Olapm(:,:,:) = Olapm(:,:,:) + Olapmk(:,:,:) * wg ! proper weight for the reducible k-point
        if (CalcDerivative) dOlapm(:,:,:,:) = dOlapm(:,:,:,:) + dOlapmk(:,:,:,:) * wg
        if (CalcStrain) dOlapStress=dOlapStress+dOlapStressK*wg
        deallocate( DMFTU )

        if (abs(projector).ge.5) deallocate(phi_jl,dphi_jl)
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
  if (CalcDerivative) then
     do idir=1,3
        CALL Reduce_MPI0(dOlapm(:,:,:,idir), maxdim, ncix, max_nbands)
     enddo
  endif
  if (CalcStrain) then
     do idir=1,6
        CALL Reduce_MPI0(dOlapStress(:,:,:,idir),maxdim,ncix,max_nbands)
     enddo
  endif

  ! wgh0 already contains the multiplicity of each irreducible k-point, but
  ! P_k^+P_k is a matrix in the correlated-orbital basis and is therefore not
  ! a symmetry scalar.  Reconstruct the missing members of every k-star now,
  ! after the MPI-reduced k sum.  For derivatives the Cartesian vector/tensor
  ! indices must be rotated together with the local-orbital indices.  This has
  ! to precede S^(-1/2) and its Frechet derivative: symmetrizing the nonlinear
  ! normalization matrix afterward would not be mathematically equivalent.
  CALL SymmetrizeOverlap(Olapm,cfX,cix_orb,cixdim,iSx,iorbital,nindo, &
       norbitals,maxdim2)
  if (CalcDerivative) then
     CALL SymmetrizeOverlapDerivative(dOlapm,3,cfX,cix_orb,cixdim, &
          iSx,iorbital,nindo,norbitals,maxdim2)
  endif
  if (CalcStrain) then
     CALL SymmetrizeOverlapDerivative(dOlapStress,6,cfX,cix_orb, &
          cixdim,iSx,iorbital,nindo,norbitals,maxdim2)
  endif
  
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
  ! If SOlapm.dat was present, Xsaved is the normalization used by the
  ! converged DMFT calculation.  The overlap traversal above is needed to
  ! obtain S_A, not to redefine X.  Restoring X also provides a useful test:
  ! adding force derivatives must not shift EF in an otherwise identical run.
  if (KeepNormalization) SOlapm=Xsaved
  if (CalcDerivative) then
     CALL ComputeNormalizationDerivative(dXrel,dOlapm,SOlapm,Sigind, &
          cixdim,maxdim,ncix,SIMPLE,3)
     deallocate(dOlapm,dOlapmk)
  endif
  if (CalcStrain) then
     CALL ComputeNormalizationDerivative(dXstress,dOlapStress,SOlapm, &
          Sigind,cixdim,maxdim,ncix,SIMPLE,6)
     deallocate(dOlapStress,dOlapStressK)
  endif
  if (CalcDerivative .or. CalcStrain) deallocate(olp2)
  if (KeepNormalization) deallocate(Xsaved)
  call cputim(t2c)
  call walltim(t2w)
  time_dmf0 = t2c-t1c
  time_dmf0w= t2w-t1w
END SUBROUTINE cmp_overlap


!===============================================================================
! DERIVATIVE OF THE GLOBAL PROJECTOR NORMALIZATION
!===============================================================================
!
! Why this term is needed
! -----------------------
! Let P_k be the raw LAPW-to-local-orbital projector.  In this implementation
! the projector used by DMFT is not P_k itself, but
!
!                         U_k = P_k X,
!
! where the normalization matrix is common to all k-points,
!
!               S = sum_k w_k P_k^+ P_k,       X = S^(-1/2).
!
! The old force code differentiated P_k and subsequently multiplied that
! derivative by X.  This gives the first term in
!
!                 U_{k,A} = P_{k,A} X + P_k X_A,
!
! but misses P_k X_A.  The omitted term is especially important when the
! correlated and uncorrelated atoms move differently: translating the entire
! crystal must leave the normalized projector unchanged, and that cancellation
! involves both terms above.  Ignoring X_A therefore violates the acoustic sum
! rule even when the explicit derivative P_{k,A} is accurate.
!
! Derivative of X=S^(-1/2)
! ------------------------
! For a matrix S, X_A is NOT obtained by differentiating every matrix element
! as if S were diagonal.  Introduce Y=X^(-1)=S^(1/2).  Differentiating
!
!                            X S X = 1
!
! and eliminating S gives the Sylvester equation
!
!                  Y X_A + X_A Y = -X S_A X,                  (1)
!
! with
!
!          S_A = sum_k w_k (P_{k,A}^+ P_k + P_k^+ P_{k,A}).   (2)
!
! The k-point weights are fixed band-integration weights, hence w_{k,A}=0.
! Equation (1) is solved in the eigenbasis of X.  If x_i are its eigenvalues,
!
!          (X_A)'_ij = -x_i x_j (S_A)'_ij/(1/x_i+1/x_j).
!
! It is convenient for the projector update to return the relative derivative
!
!                         B_A = X^(-1) X_A.
!
! In the eigenbasis of X this is
!
!              (B_A)'_ij = -x_j (S_A)'_ij/(1/x_i+1/x_j),      (3)
!
! which is exactly the expression implemented below.
!
! Relation to the force-channel convention
! ----------------------------------------
! Build_DMFT_Projector stores the three force channels D_A with the convention
!
!                            U_A = -i D_A.
!
! Once the ordinary channels have been normalized, Ubar=P X, their missing
! normalization contribution is Ubar B_A.  Therefore the stored D channel must
! receive
!
!                         D_A <- D_A + i Ubar B_A.
!
! This is performed by AddProjectorNormalizationDerivative.
!
! SOlapm.dat subtlety
! -------------------
! SOlapm.dat contains the normalization X used by the converged DMFT run.  We
! recompute S_A from the current LAPW vectors, but deliberately restore and use
! the saved X.  Replacing X by a newly reconstructed value in a force-only step
! changes the physical projector and can move the chemical potential.  A force
! derivative must differentiate the projector that was actually used, not
! silently define a nearby projector.
!===============================================================================
SUBROUTINE ComputeNormalizationDerivative(dXrel,dS,X,Sigind,cixdim, &
     maxdim,ncix,SIMPLE,nq)
  IMPLICIT NONE
  INTEGER, intent(in) :: maxdim,ncix,nq,cixdim(ncix), &
       Sigind(maxdim,maxdim,ncix)
  LOGICAL, intent(in) :: SIMPLE
  COMPLEX*16, intent(in)  :: dS(maxdim,maxdim,ncix,nq), &
       X(maxdim,maxdim,ncix)
  COMPLEX*16, intent(out) :: dXrel(maxdim,maxdim,ncix,nq)
  INTEGER :: icix, cixdm, cixdms, ip, iq, idir, info
  INTEGER :: lwork, lrwork, liwork
  INTEGER, allocatable :: cini(:), iwork(:)
  REAL*8, allocatable :: xeig(:), rwork(:)
  COMPLEX*16, allocatable :: Xvec(:,:), dSt(:,:), Bt(:,:), Bout(:,:), work(:)

  dXrel=0.d0
  if (SIMPLE) then
     do icix=1,ncix
        do ip=1,cixdim(icix)
           if (Sigind(ip,ip,icix).eq.0) cycle
           do idir=1,nq
              ! X=S^{-1/2}: X^{-1} X_A = -1/2 X^2 S_A.
              dXrel(ip,ip,icix,idir) = -0.5d0*X(ip,ip,icix)**2*dS(ip,ip,icix,idir)
           enddo
        enddo
     enddo
     return
  endif

  do icix=1,ncix
     cixdm=cixdim(icix)
     allocate(cini(cixdm))
     cixdms=0
     do ip=1,cixdm
        if (Sigind(ip,ip,icix).ne.0) then
           cixdms=cixdms+1
           cini(cixdms)=ip
        endif
     enddo
     if (cixdms.eq.0) then
        deallocate(cini)
        cycle
     endif

     lwork=4*cixdms+cixdms*cixdms
     lrwork=10*cixdms+2*cixdms*cixdms+1
     liwork=3+10*cixdms
     allocate(Xvec(cixdms,cixdms), dSt(cixdms,cixdms), Bt(cixdms,cixdms), Bout(cixdms,cixdms))
     allocate(xeig(cixdms), work(lwork), rwork(lrwork), iwork(liwork))
     do ip=1,cixdms
        do iq=1,cixdms
           Xvec(ip,iq)=X(cini(ip),cini(iq),icix)
        enddo
     enddo
     CALL ZHEEVD('V','U',cixdms,Xvec,cixdms,xeig,work,lwork,rwork,lrwork,iwork,liwork,info)
     if (info.ne.0 .or. minval(xeig).le.1.d-14) then
        write(6,*) 'ERROR: diagonalization of S^{-1/2} failed in normalization derivative',icix,info
        STOP 'ComputeNormalizationDerivative failed'
     endif

     do idir=1,nq
        do ip=1,cixdms
           do iq=1,cixdms
              Bout(ip,iq)=dS(cini(ip),cini(iq),icix,idir)
           enddo
        enddo
        dSt=matmul(conjg(transpose(Xvec)),matmul(Bout,Xvec))
        do ip=1,cixdms
           do iq=1,cixdms
              ! Solve Y X_A+X_A Y=-X S_A X with Y=X^{-1}, then
              ! return the relative derivative B_A=X^{-1}X_A.  The actual
              ! normalized-projector derivative is Ubar_A=...+Ubar B_A.
              Bt(ip,iq)=-xeig(iq)*dSt(ip,iq)/(1.d0/xeig(ip)+1.d0/xeig(iq))
           enddo
        enddo
        Bout=matmul(Xvec,matmul(Bt,conjg(transpose(Xvec))))
        do ip=1,cixdms
           do iq=1,cixdms
              dXrel(cini(ip),cini(iq),icix,idir)=Bout(ip,iq)
           enddo
        enddo
     enddo
     deallocate(Xvec,dSt,Bt,Bout,xeig,work,rwork,iwork,cini)
  enddo
END SUBROUTINE ComputeNormalizationDerivative


! Dump the analytic derivative of X=S^(-1/2) for comparison with independently
! strained SOlapm.dat files.  The internal dXrel is B=X^(-1)X_,q, hence the
! actual derivative written below is X_,q=X*B.  Nothing is written unless the
! marker file normalization_strain_fd.in exists in the case directory.
SUBROUTINE WriteNormalizationStrainDiagnostic(X,dXrel,cixdim,maxdim,ncix)
  USE com_mpi,ONLY: myrank,master
  IMPLICIT NONE
  INTEGER,intent(in) :: maxdim,ncix,cixdim(ncix)
  COMPLEX*16,intent(in) :: X(maxdim,maxdim,ncix)
  COMPLEX*16,intent(in) :: dXrel(maxdim,maxdim,ncix,6)
  LOGICAL :: enabled
  INTEGER :: icix,cdim,iq,i,j,fh
  COMPLEX*16,allocatable :: dX(:,:)

  INQUIRE(file='normalization_strain_fd.in',exist=enabled)
  if (.not.enabled .or. myrank.ne.master) return
  fh=997
  open(fh,file='normalization_strain_derivative.dat',status='replace', &
       action='write')
  write(fh,'(A)') '# Analytic derivative of global projector normalization'
  write(fh,'(A)') '# Voigt order: xx yy zz yz xz xy; entries are Re Im pairs'
  write(fh,'(A,I8)') 'NCIX ',ncix
  do icix=1,ncix
     cdim=cixdim(icix)
     allocate(dX(cdim,cdim))
     write(fh,'(A,2I8)') 'CIX ',icix,cdim
     write(fh,'(A)') 'X'
     do i=1,cdim
        write(fh,'(*(2ES25.16,1X))') (dble(X(i,j,icix)), &
             aimag(X(i,j,icix)),j=1,cdim)
     enddo
     do iq=1,6
        dX=matmul(X(:cdim,:cdim,icix), &
             dXrel(:cdim,:cdim,icix,iq))
        write(fh,'(A,I8)') 'DX ',iq
        do i=1,cdim
           write(fh,'(*(2ES25.16,1X))') (dble(dX(i,j)),aimag(dX(i,j)), &
                j=1,cdim)
        enddo
     enddo
     deallocate(dX)
  enddo
  close(fh)
  WRITE(6,'(A)') &
       ':DMFT_NORMALIZATION_STRAIN_DUMP: normalization_strain_derivative.dat'
END SUBROUTINE WriteNormalizationStrainDiagnostic


SUBROUTINE AddProjectorNormalizationDerivative(DMFTU, dXrel, cix_orb, cixdim, nindo, iSx, nbands, maxdim2, norbitals, maxdim, ncix)
  IMPLICIT NONE
  INTEGER, intent(in) :: nbands, maxdim2, norbitals, maxdim, ncix
  INTEGER, intent(in) :: cix_orb(norbitals), cixdim(ncix), nindo(norbitals), iSx(maxdim2,norbitals)
  COMPLEX*16, intent(in) :: dXrel(maxdim,maxdim,ncix,3)
  COMPLEX*16, intent(inout) :: DMFTU(nbands,maxdim2,norbitals,4)
  INTEGER :: icix, cixdm, iorb, nind, ind, idir
  COMPLEX*16, allocatable :: Ucix(:,:), corr(:,:)
  COMPLEX*16, parameter :: cone=(1.d0,0.d0), czero=(0.d0,0.d0), IMAG=(0.d0,1.d0)

  do icix=1,ncix
     cixdm=cixdim(icix)
     allocate(Ucix(nbands,cixdm),corr(nbands,cixdm))
     Ucix=0.d0
     do iorb=1,norbitals
        if (cix_orb(iorb).ne.icix) cycle
        nind=nindo(iorb)
        do ind=1,nind
           Ucix(:,iSx(ind,iorb))=DMFTU(:,ind,iorb,1)
        enddo
     enddo
     do idir=1,3
        call zgemm('N','N',nbands,cixdm,cixdm,cone,Ucix,nbands,dXrel(:,:,icix,idir),maxdim,czero,corr,nbands)
        do iorb=1,norbitals
           if (cix_orb(iorb).ne.icix) cycle
           nind=nindo(iorb)
           do ind=1,nind
              ! The force channel D uses U_A=-iD, hence D_norm=i Ubar B_A.
              DMFTU(:,ind,iorb,idir+1)=DMFTU(:,ind,iorb,idir+1)+IMAG*corr(:,iSx(ind,iorb))
           enddo
        enddo
     enddo
     deallocate(Ucix,corr)
  enddo
END SUBROUTINE AddProjectorNormalizationDerivative


SUBROUTINE AddProjectorStressNormalization(dU,U,dXrel,cix_orb, &
     cixdim,nindo,iSx,nbands,maxdim2,norbitals,maxdim,ncix)
  ! Complete dU=P_,q X by adding the derivative of the common matrix
  ! normalization, U B_q with B_q=X^{-1}X_,q.  Unlike the historical force
  ! channels, dU here is the derivative itself and therefore carries no
  ! extra factor of i.
  IMPLICIT NONE
  INTEGER,intent(in) :: nbands,maxdim2,norbitals,maxdim,ncix
  INTEGER,intent(in) :: cix_orb(norbitals),cixdim(ncix), &
       nindo(norbitals),iSx(maxdim2,norbitals)
  COMPLEX*16,intent(in) :: U(nbands,maxdim2,norbitals)
  COMPLEX*16,intent(in) :: dXrel(maxdim,maxdim,ncix,6)
  COMPLEX*16,intent(inout) :: dU(nbands,maxdim2,norbitals,6)
  INTEGER :: icix,cixdm,iorb,nind,ind,iq
  COMPLEX*16,allocatable :: Ucix(:,:),corr(:,:)
  COMPLEX*16,parameter :: cone=(1.d0,0.d0),czero=(0.d0,0.d0)

  do icix=1,ncix
     cixdm=cixdim(icix)
     allocate(Ucix(nbands,cixdm),corr(nbands,cixdm))
     Ucix=czero
     do iorb=1,norbitals
        if (cix_orb(iorb).ne.icix) cycle
        nind=nindo(iorb)
        do ind=1,nind
           Ucix(:,iSx(ind,iorb))=U(:,ind,iorb)
        enddo
     enddo
     do iq=1,6
        call zgemm('N','N',nbands,cixdm,cixdm,cone,Ucix,nbands, &
             dXrel(:,:,icix,iq),maxdim,czero,corr,nbands)
        do iorb=1,norbitals
           if (cix_orb(iorb).ne.icix) cycle
           do ind=1,nindo(iorb)
              dU(:,ind,iorb,iq)=dU(:,ind,iorb,iq)+ &
                   corr(:,iSx(ind,iorb))
           enddo
        enddo
     enddo
     deallocate(Ucix,corr)
  enddo
END SUBROUTINE AddProjectorStressNormalization


!===============================================================================
! RESPONSE OF A PROJECTOR DEFINED IN A FIXED BAND-INDEX WINDOW
!===============================================================================
!
! The missing Sternheimer contribution
! -------------------------------------
! The DMFT window W is a fixed set of band indices, for example bands 12--31 at
! every k-point.  Even though the membership of W does not change under an
! infinitesimal displacement, a state inside W has the usual first-order
! response to every state outside W:
!
!   |psi_{n,A}>_outside = sum_{a not in W} |psi_a>
!       <psi_a|H_A-eps_n O_A|psi_n>/(eps_n-eps_a),       n in W.       (4)
!
! H_A and O_A are the LAPW Hamiltonian and overlap derivatives already derived
! and implemented for the ordinary Pulay force.  The original DMFT force code
! differentiated projector matrix elements whose band index lies inside W, but
! did not include the outside-band resolution of identity in (4).  This is why
! increasing the number of LAPW1 bands improves the force sum: it improves the
! numerical representation of this missing response.
!
! A direct Sternheimer solve for every atom and Cartesian direction would be
! expensive.  We can instead rewrite its contribution as an additional
! band-space density matrix and reuse the existing H_A and O_A contractions.
!
! Quantities already available from the DMFT frequency sum
! --------------------------------------------------------
! Use U_W for the normalized projectors restricted to the DMFT window and
! U_full for the same projectors evaluated for every band stored in case.vector.
! At each frequency the embedded self-energy is
!
!                       Sigma_emb = U_W Sigma U_W^+.
!
! cmp_dmft_weights2 already supplies the two Hermitian matrices
!
!      rho = sum_omega G(omega),
!      EG  = sum_omega [epsilon + Sigma_emb(omega)] G(omega),
!
! including all temperature and k-point weights.  Therefore no new frequency
! loop and no two-particle vertex are required:
!
!       B = EG - epsilon rho
!         = sum_omega Sigma_emb(omega) G(omega)
!         = U_W C,                                                   (5)
!
! where C=sum_omega Sigma U_W^+ G.  Stationarity of the DFT+eDMFT Klein
! functional is what allows the implicit self-energy response (and hence the
! vertex) to cancel; only the explicit one-particle projector response remains.
!
! Reconstructing C and extending it outside W
! -------------------------------------------
! When the active projector columns are linearly independent in W,
!
!                  C = (U_W^+ U_W)^(-1) U_W^+ B.
!
! In practice symmetries or a narrow window can make the metric rank deficient,
! so the code uses an eigenvalue-cutoff Moore-Penrose pseudoinverse:
!
!          K = U_full (U_W^+ U_W)^+ U_W^+ B.                         (6)
!
! K_{an} is the frequency-summed DMFT coefficient multiplying the
! Sternheimer matrix element between an outside state a and a window state n.
!
! Response density and energy-density matrices
! ---------------------------------------------
! Combining (4) and (6), define for a outside W and n inside W
!
!       rhoResponse_{a n} = K_{a n}/(eps_n-eps_a),
!       rhoResponse_{n a} = rhoResponse_{a n}^*,                     (7)
!
!       eResponse_{a n}   = eps_n rhoResponse_{a n},
!       eResponse_{n a}   = eResponse_{a n}^*.                       (8)
!
! The omitted contribution to the stationary-functional derivative is then
!
!       Tr[rhoResponse H_A] - Tr[eResponse O_A].                     (9)
!
! This has exactly the same algebraic form as the Pulay contractions already
! evaluated by Force1_DMFT, Force2, and Force3.  l2main diagonalizes the signed
! Hermitian rhoResponse only as a computational representation, transforms
! eResponse to the same basis, and calls cmp_MT_force_only.
!
! rhoResponse is NOT a correction to the physical charge density
! --------------------------------------------------------------
! Equations (7)--(9) are a device for evaluating a derivative of the projector.
! Adding rhoResponse to w_RHOLM or the interstitial density would change the
! state whose force is being evaluated.  Consequently cmp_MT_force_only reuses
! only the force kernels and intentionally performs no charge accumulation.
!
! Finite stored-band approximation
! --------------------------------
! The exact sum in (4) contains all bands.  Here "outside" means all additional
! states already present in case.vector/case.vectorso.  Convergence can therefore
! be checked, and improved when necessary, by asking LAPW1 to write more bands;
! the DMFT window itself remains unchanged.  Exact cross-window degeneracies
! require degenerate perturbation theory.  Isolated denominators below 1e-10 Ry
! are skipped and reported rather than divided by an unstable number.
!===============================================================================
SUBROUTINE BuildWindowResponseMatrices(rhoResponse,eResponse, &
     projectorGradient,DMFTUfull,Aweight,AEweight,E,nemin,nemax, &
     DM_nemin,DM_nemax,nbands_full,nbands,maxdim2,norbitals,maxdim, &
     ncix,cix_orb,cixdim,nindo,iSx,Sigind)
  USE com_mpi, ONLY: Qprint
  IMPLICIT NONE
  INTEGER, intent(in) :: nemin, nemax, DM_nemin, DM_nemax
  INTEGER, intent(in) :: nbands_full, nbands, maxdim2, norbitals, maxdim, ncix
  INTEGER, intent(in) :: cix_orb(norbitals), cixdim(ncix), nindo(norbitals)
  INTEGER, intent(in) :: iSx(maxdim2,norbitals), Sigind(maxdim,maxdim,ncix)
  REAL*8, intent(in) :: E(*)
  COMPLEX*16, intent(in) :: DMFTUfull(nbands_full,maxdim2,norbitals)
  COMPLEX*16, intent(in) :: Aweight(nbands,nbands), AEweight(nbands,nbands)
  COMPLEX*16, intent(out) :: rhoResponse(nbands_full,nbands_full)
  COMPLEX*16, intent(out) :: eResponse(nbands_full,nbands_full)
  COMPLEX*16, intent(out) :: projectorGradient(nbands,maxdim2,norbitals)
  INTEGER :: icix, ip, iq, iorb, ind, ia, iw, in, iw0, iw1
  INTEGER :: nactive, icolumn, info, lwork, lrwork, liwork, nsmall, ndegenerate
  INTEGER, allocatable :: activeColumn(:,:), iwork(:)
  REAL*8 :: eigmax, cutoff, denominator, epsn, responseNorm
  REAL*8 :: sourceNorm, reconstructionResidual, smallestKept
  REAL*8, allocatable :: seig(:), rwork(:)
  LOGICAL :: projectorDiagnostic
  CHARACTER(LEN=64) :: diagnosticEnvironment
  INTEGER :: diagnosticStatus
  COMPLEX*16 :: response
  COMPLEX*16, allocatable :: Uall(:,:), Uwindow(:,:), Svec(:,:), Sinv(:,:)
  COMPLEX*16, allocatable :: B(:,:), reconstructedB(:,:), localPart(:,:), &
       Kfull(:,:), work(:)
  COMPLEX*16, parameter :: cone=(1.d0,0.d0), czero=(0.d0,0.d0)

  rhoResponse=czero
  eResponse=czero
  projectorGradient=czero
  diagnosticEnvironment=''
  CALL GET_ENVIRONMENT_VARIABLE('EDMFT_STRESS_PROJECTOR_DIAG', &
       diagnosticEnvironment,STATUS=diagnosticStatus)
  projectorDiagnostic=diagnosticStatus.eq.0 .and. &
       len_trim(diagnosticEnvironment).gt.0 .and. &
       trim(diagnosticEnvironment).ne.'0'
  iw0=DM_nemin-nemin+1
  iw1=DM_nemax-nemin+1
  if (iw0.lt.1 .or. iw1.gt.nbands_full .or. iw1-iw0+1.ne.nbands) then
     WRITE(6,*) 'ERROR: inconsistent fixed DMFT window in BuildWindowResponseMatrices', &
          nemin,nemax,DM_nemin,DM_nemax,nbands_full,nbands
     STOP 'BuildWindowResponseMatrices: inconsistent band bounds'
  endif

  ! Keep every local state which occurs in at least one self-energy matrix
  ! element.  Usually testing Sigind(ip,ip) is sufficient, but the row/column
  ! test also covers a purely off-diagonal block.
  allocate(activeColumn(maxdim,ncix))
  activeColumn=0
  nactive=0
  do icix=1,ncix
     do ip=1,cixdim(icix)
        if (any(Sigind(ip,1:cixdim(icix),icix).gt.0) .or. &
            any(Sigind(1:cixdim(icix),ip,icix).gt.0)) then
           nactive=nactive+1
           activeColumn(ip,icix)=nactive
        endif
     enddo
  enddo
  if (nactive.eq.0) then
     deallocate(activeColumn)
     return
  endif

  allocate(Uall(nbands_full,nactive),Uwindow(nbands,nactive))
  Uall=czero
  do iorb=1,norbitals
     icix=cix_orb(iorb)
     if (icix.eq.0) cycle
     do ind=1,nindo(iorb)
        ip=iSx(ind,iorb)
        icolumn=activeColumn(ip,icix)
        if (icolumn.gt.0) Uall(:,icolumn)=DMFTUfull(:,ind,iorb)
     enddo
  enddo
  Uwindow=Uall(iw0:iw1,:)

  allocate(Svec(nactive,nactive),Sinv(nactive,nactive),seig(nactive))
  call zgemm('C','N',nactive,nactive,nbands,cone,Uwindow,nbands,Uwindow,nbands, &
       czero,Svec,nactive)
  lwork=2*nactive+nactive*nactive
  lrwork=5*nactive+2*nactive*nactive+1
  liwork=3+5*nactive
  allocate(work(lwork),rwork(lrwork),iwork(liwork))
  CALL ZHEEVD('V','U',nactive,Svec,nactive,seig,work,lwork,rwork,lrwork,iwork,liwork,info)
  if (info.ne.0) then
     WRITE(6,*) 'ERROR: projector metric diagonalization failed in window response',info
     STOP 'BuildWindowResponseMatrices: ZHEEVD failed'
  endif

  eigmax=maxval(seig)
  cutoff=max(1.d-13,1.d-11*eigmax)
  smallestKept=huge(1.d0)
  Sinv=Svec
  nsmall=0
  do ip=1,nactive
     if (seig(ip).gt.cutoff) then
        Sinv(:,ip)=Sinv(:,ip)/seig(ip)
        smallestKept=min(smallestKept,seig(ip))
     else
        Sinv(:,ip)=czero
        nsmall=nsmall+1
     endif
  enddo
  Sinv=matmul(Sinv,conjg(transpose(Svec)))

  allocate(B(nbands,nbands),localPart(nactive,nbands),Kfull(nbands_full,nbands))
  do iw=1,nbands
     do in=1,nbands
        ! eps multiplies the Green's function from the left.
        B(iw,in)=AEweight(iw,in)-E(DM_nemin+iw-1)*Aweight(iw,in)
     enddo
  enddo
  call zgemm('C','N',nactive,nbands,nbands,cone,Uwindow,nbands,B,nbands, &
       czero,localPart,nactive)
  localPart=matmul(Sinv,localPart)
  if (projectorDiagnostic .and. Qprint) then
     ! B=rhoE-epsilon*rho must lie in the column space of the normalized
     ! projector U_W if it is to be represented as U_W C.  This check is
     ! independent of the subsequent Sternheimer denominator and therefore
     ! cleanly distinguishes an algebra/convention error in C from an
     ! incomplete complementary-band sum.
     allocate(reconstructedB(nbands,nbands))
     reconstructedB=matmul(Uwindow,localPart)
     sourceNorm=sqrt(sum(abs(B)**2))
     reconstructionResidual=sqrt(sum(abs(B-reconstructedB)**2))/ &
          max(sourceNorm,tiny(1.d0))
     WRITE(6,'(A,4ES18.8,I6)') ':DMFT_PROJECTOR_RECONSTRUCTION: ', &
          sourceNorm,reconstructionResidual,eigmax,smallestKept,nsmall
     deallocate(reconstructedB)
  endif
  ! localPart(m,n)=sum_omega [Sigma U^+ G]_{m n}, including the
  ! irreducible-k weight.  It is the one-particle gradient conjugate to the
  ! explicit projector derivative: dGamma_proj=2 Re sum_mn C_mn dU_nm.
  ! This is an energy derivative.  WIEN2k's analytic stress uses the
  ! +dGamma/d(strain)/Omega convention, so l2main.f90 adds this contraction.
  ! (The corresponding atomic force has the opposite sign.)
  do iorb=1,norbitals
     icix=cix_orb(iorb)
     if (icix.eq.0) cycle
     do ind=1,nindo(iorb)
        ip=iSx(ind,iorb)
        icolumn=activeColumn(ip,icix)
        if (icolumn.gt.0) projectorGradient(:,ind,iorb)= &
             localPart(icolumn,:)
     enddo
  enddo
  call zgemm('N','N',nbands_full,nbands,nactive,cone,Uall,nbands_full, &
       localPart,nactive,czero,Kfull,nbands_full)

  ndegenerate=0
  do ia=1,nbands_full
     if (ia.ge.iw0 .and. ia.le.iw1) cycle
     do iw=1,nbands
        in=iw0+iw-1
        epsn=E(DM_nemin+iw-1)
        denominator=epsn-E(nemin+ia-1)
        if (abs(denominator).le.1.d-10) then
           ndegenerate=ndegenerate+1
           cycle
        endif
        response=Kfull(ia,iw)/denominator
        rhoResponse(ia,in)=response
        rhoResponse(in,ia)=conjg(response)
        eResponse(ia,in)=epsn*response
        eResponse(in,ia)=conjg(eResponse(ia,in))
     enddo
  enddo

  responseNorm=sqrt(sum(abs(rhoResponse)**2))
  if (Qprint .and. (nsmall.gt.0 .or. ndegenerate.gt.0)) then
     WRITE(6,'(A,2I6,A,ES12.4)') 'Window-response pseudoinverse: small/degenerate=', &
          nsmall,ndegenerate,'  ||rho_resp||=',responseNorm
  endif

  deallocate(activeColumn,Uall,Uwindow,Svec,Sinv,seig,work,rwork,iwork)
  deallocate(B,localPart,Kfull)
END SUBROUTINE BuildWindowResponseMatrices


! Save the one-particle gradient conjugate to the normalized DMFT projector.
! This default-off diagnostic permits a direct finite-difference check of
!
!       delta Gamma_projector = 2 Re sum(C * delta U)
!
! without routing either side through an LAPW stress kernel.  It is enabled
! by EDMFT_STRESS_PROJECTOR_DIAG and writes one rank-independent file per
! global k point.  The band energies are included so independently
! diagonalized endpoint projectors can be gauge-aligned inside exactly
! degenerate subspaces before forming their central difference.
SUBROUTINE WriteProjectorGradientDiagnostic(gradient,bandEnergy,ikp, &
     nbands,maxdim2,norbitals)
  USE com_mpi,ONLY: myrank,master,stop_MPI
  IMPLICIT NONE
  INTEGER,intent(in) :: ikp,nbands,maxdim2,norbitals
  COMPLEX*16,intent(in) :: gradient(nbands,maxdim2,norbitals)
  REAL*8,intent(in) :: bandEnergy(nbands)
  INTEGER,parameter :: gradientMagic=20261105,gradientVersion=1
  CHARACTER(len=64) :: environment
  CHARACTER(len=256) :: filename
  INTEGER :: status,fh,ios

  environment=''
  CALL GET_ENVIRONMENT_VARIABLE('EDMFT_STRESS_PROJECTOR_DIAG', &
       environment,STATUS=status)
  if (status.ne.0 .or. len_trim(environment).eq.0 .or. &
       trim(environment).eq.'0') return

  WRITE(filename,'(A,I8.8,A)') 'dmft_projector_gradient_k',ikp,'.bin'
  fh=972
  open(fh,file=trim(filename),status='replace',form='unformatted', &
       access='stream',action='write',iostat=ios)
  if (ios.eq.0) write(fh,iostat=ios) gradientMagic,gradientVersion, &
       ikp,nbands,maxdim2,norbitals
  if (ios.eq.0) write(fh,iostat=ios) bandEnergy
  if (ios.eq.0) write(fh,iostat=ios) gradient
  close(fh)
  if (ios.ne.0) then
     WRITE(6,'(A,I6,2A,I8)') 'ERROR: rank ',myrank, &
          ' cannot write projector-gradient diagnostic ',trim(filename),ios
     CALL stop_MPI
     STOP 'projector-gradient diagnostic write failed'
  endif
  if (myrank.eq.master .and. ikp.eq.1) WRITE(6,'(A,1X,ES22.13)') &
       ':DMFT_PROJECTOR_GRADIENT_DUMP:',sqrt(sum(abs(gradient)**2))
END SUBROUTINE WriteProjectorGradientDiagnostic


! Save the central LAPW coefficients and the un-diagonalized embedded
! self-energy moment for a kernel-independent finite difference of
!
!       2 Re sum_k Tr[M_Sigma(k) O_,epsilon(k)].
!
! DFTKS_OVERLAP_DUMP supplies O(+epsilon) and O(-epsilon).  Keeping M_Sigma
! in the original Kohn-Sham band basis avoids all pseudo-band sorting and
! degeneracy choices made later by Diagonalize_DMFT_WEIGHTS.  Ordinary runs
! return before touching the arrays because EDMFT_MSIGMA_DUMP is unset.
SUBROUTINE WriteSigmaMomentOverlapDump(ikp,n0,nnlo,iso,DM_nemin, &
     DM_nemax,nmat,nume,k3,k3lo,bk3,As,sigmaMoment,nbands)
  USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: INT32,REAL64
  USE com_mpi,ONLY: myrank,master,stop_MPI
  IMPLICIT NONE
  INTEGER,intent(in) :: ikp,n0,nnlo,iso,DM_nemin,DM_nemax
  INTEGER,intent(in) :: nmat,nume,nbands
  INTEGER,intent(in) :: k3(3,nmat),k3lo(3,nnlo)
  REAL*8,intent(in) :: bk3(3,nmat)
  COMPLEX*16,intent(in) :: As(nmat,nume,iso)
  COMPLEX*16,intent(in) :: sigmaMoment(nbands,nbands)

  CHARACTER(len=8),parameter :: magic='EDMFTMS1'
  CHARACTER(len=1024) :: environment,prefix
  CHARACTER(len=1200) :: filename
  INTEGER :: status,length,unit,ios,isize
  INTEGER(INT32) :: marker,header(10)
  INTEGER(INT32),allocatable :: labels(:,:)
  REAL(REAL64) :: kfrac(3)

  environment=''
  prefix=''
  length=0
  status=1
  CALL GET_ENVIRONMENT_VARIABLE('EDMFT_MSIGMA_DUMP',environment, &
       LENGTH=length,STATUS=status)
  if (status.ne.0 .or. length.le.0) return
  environment=adjustl(environment)
  select case(trim(environment))
  case('0','OFF','off','FALSE','false','NO','no')
     return
  case('1','ON','on','TRUE','true','YES','yes')
     prefix='dmft_msigma'
  case default
     prefix=trim(environment)
  end select

  isize=n0-nnlo
  if (isize.le.0 .or. nbands.ne.DM_nemax-DM_nemin+1) then
     if (myrank.eq.master) WRITE(6,'(A,5I10)') &
          ':DMFT_MSIGMA_DUMP_BAD_DIMENSIONS:',n0,nnlo,nbands, &
          DM_nemin,DM_nemax
     CALL stop_MPI
     STOP 'invalid M_Sigma dump dimensions'
  endif

  allocate(labels(3,n0))
  labels(:,1:isize)=INT(k3(:,1:isize),INT32)
  if (nnlo.gt.0) labels(:,isize+1:n0)=INT(k3lo(:,1:nnlo),INT32)
  kfrac=REAL(bk3(:,1)-DBLE(k3(:,1)),REAL64)
  marker=INT(Z'01020304',INT32)
  ! Flag bit zero asserts that sigmaMoment already contains the normalized
  ! irreducible-k weight wkp.  The Python helper must not multiply it again.
  header=(/ INT(1,INT32),INT(ikp,INT32),INT(myrank,INT32), &
       INT(n0,INT32),INT(nnlo,INT32),INT(nbands,INT32), &
       INT(iso,INT32),INT(DM_nemin,INT32),INT(DM_nemax,INT32), &
       INT(1,INT32) /)
  WRITE(filename,'(A,".k",I8.8,".bin")') trim(prefix),ikp

  OPEN(NEWUNIT=unit,FILE=trim(filename),ACCESS='STREAM', &
       FORM='UNFORMATTED',STATUS='REPLACE',ACTION='WRITE',IOSTAT=ios)
  if (ios.eq.0) WRITE(unit,IOSTAT=ios) magic
  if (ios.eq.0) WRITE(unit,IOSTAT=ios) marker
  if (ios.eq.0) WRITE(unit,IOSTAT=ios) header
  if (ios.eq.0) WRITE(unit,IOSTAT=ios) kfrac
  if (ios.eq.0) WRITE(unit,IOSTAT=ios) labels
  if (ios.eq.0) WRITE(unit,IOSTAT=ios) &
       As(1:n0,DM_nemin:DM_nemax,1:iso)
  if (ios.eq.0) WRITE(unit,IOSTAT=ios) sigmaMoment
  CLOSE(unit)
  deallocate(labels)

  if (ios.ne.0) then
     WRITE(6,'(A,I6,2A,I8)') 'ERROR: rank ',myrank, &
          ' cannot write M_Sigma diagnostic ',trim(filename),ios
     CALL stop_MPI
     STOP 'M_Sigma diagnostic write failed'
  endif
  if (myrank.eq.master .and. ikp.eq.1) WRITE(6,'(A,1X,A)') &
       ':DMFT_MSIGMA_DUMP:',trim(filename)
END SUBROUTINE WriteSigmaMomentOverlapDump


! Save the two stationary one-particle matrices in a common stored-band
! basis for a direct finite difference of
!
!       Tr[rho H_,q] - Tr[rhoE O_,q].
!
! Bands below the DMFT window are the ordinary fully occupied LAPW bands;
! the DMFT window contains the un-diagonalized Aweight and AEweight blocks.
! Both already contain the normalized irreducible-k and spin weight.  This
! path is entirely inert unless EDMFT_ONEPARTICLE_DUMP is set.
SUBROUTINE WriteOneParticleMomentDump(ikp,n0,nnlo,iso,nemin, &
     DM_nemin,DM_nemax,nmat,nume,k3,k3lo,bk3,As,rhoWindow, &
     rhoEWindow,bandEnergy,nbands,wkp)
  USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: INT32,REAL64
  USE com_mpi,ONLY: myrank,master,stop_MPI
  IMPLICIT NONE
  INTEGER,intent(in) :: ikp,n0,nnlo,iso,nemin,DM_nemin,DM_nemax
  INTEGER,intent(in) :: nmat,nume,nbands
  INTEGER,intent(in) :: k3(3,nmat),k3lo(3,nnlo)
  REAL*8,intent(in) :: bk3(3,nmat),bandEnergy(nume)
  COMPLEX*16,intent(in) :: As(nmat,nume,iso)
  COMPLEX*16,intent(in) :: rhoWindow(nbands,nbands)
  COMPLEX*16,intent(in) :: rhoEWindow(nbands,nbands),wkp

  CHARACTER(len=8),parameter :: magic='EDMFTRH1'
  CHARACTER(len=1024) :: environment,prefix
  CHARACTER(len=1200) :: filename
  INTEGER :: status,length,unit,ios,isize,nfull,nlower,iw0,i
  INTEGER(INT32) :: marker,header(12)
  INTEGER(INT32),allocatable :: labels(:,:)
  REAL(REAL64) :: kfrac(3)
  COMPLEX(REAL64),allocatable :: rho(:,:),rhoE(:,:)

  environment=''; prefix=''; length=0; status=1
  CALL GET_ENVIRONMENT_VARIABLE('EDMFT_ONEPARTICLE_DUMP',environment, &
       LENGTH=length,STATUS=status)
  if (status.ne.0 .or. length.le.0) return
  environment=adjustl(environment)
  select case(trim(environment))
  case('0','OFF','off','FALSE','false','NO','no')
     return
  case('1','ON','on','TRUE','true','YES','yes')
     prefix='dmft_rho_rhoe'
  case default
     prefix=trim(environment)
  end select

  isize=n0-nnlo
  nfull=DM_nemax-nemin+1
  nlower=DM_nemin-nemin
  iw0=nlower+1
  if (isize.le.0 .or. nfull.le.0 .or. nlower.lt.0 .or. &
       nbands.ne.DM_nemax-DM_nemin+1) then
     if (myrank.eq.master) WRITE(6,'(A,6I10)') &
          ':DMFT_ONEPARTICLE_DUMP_BAD_DIMENSIONS:',n0,nnlo,nfull, &
          nbands,DM_nemin,DM_nemax
     CALL stop_MPI
     STOP 'invalid one-particle dump dimensions'
  endif

  allocate(labels(3,n0),rho(nfull,nfull),rhoE(nfull,nfull))
  labels(:,1:isize)=INT(k3(:,1:isize),INT32)
  if (nnlo.gt.0) labels(:,isize+1:n0)=INT(k3lo(:,1:nnlo),INT32)
  rho=(0.0_REAL64,0.0_REAL64)
  rhoE=(0.0_REAL64,0.0_REAL64)
  do i=1,nlower
     rho(i,i)=CMPLX(DBLE(wkp),AIMAG(wkp),KIND=REAL64)
     rhoE(i,i)=CMPLX(DBLE(wkp)*bandEnergy(nemin+i-1), &
          AIMAG(wkp)*bandEnergy(nemin+i-1),KIND=REAL64)
  enddo
  rho(iw0:nfull,iw0:nfull)=rhoWindow
  rhoE(iw0:nfull,iw0:nfull)=rhoEWindow

  kfrac=REAL(bk3(:,1)-DBLE(k3(:,1)),REAL64)
  marker=INT(Z'01020304',INT32)
  header=(/ INT(1,INT32),INT(ikp,INT32),INT(myrank,INT32), &
       INT(n0,INT32),INT(nnlo,INT32),INT(nfull,INT32), &
       INT(iso,INT32),INT(nemin,INT32),INT(DM_nemin,INT32), &
       INT(DM_nemax,INT32),INT(1,INT32),INT(0,INT32) /)
  WRITE(filename,'(A,".k",I8.8,".bin")') trim(prefix),ikp
  OPEN(NEWUNIT=unit,FILE=trim(filename),ACCESS='STREAM', &
       FORM='UNFORMATTED',STATUS='REPLACE',ACTION='WRITE',IOSTAT=ios)
  if (ios.eq.0) WRITE(unit,IOSTAT=ios) magic
  if (ios.eq.0) WRITE(unit,IOSTAT=ios) marker
  if (ios.eq.0) WRITE(unit,IOSTAT=ios) header
  if (ios.eq.0) WRITE(unit,IOSTAT=ios) kfrac
  if (ios.eq.0) WRITE(unit,IOSTAT=ios) labels
  if (ios.eq.0) WRITE(unit,IOSTAT=ios) &
       REAL(bandEnergy(nemin:DM_nemax),REAL64)
  if (ios.eq.0) WRITE(unit,IOSTAT=ios) &
       As(1:n0,nemin:DM_nemax,1:iso)
  if (ios.eq.0) WRITE(unit,IOSTAT=ios) rho
  if (ios.eq.0) WRITE(unit,IOSTAT=ios) rhoE
  CLOSE(unit)
  deallocate(labels,rho,rhoE)
  if (ios.ne.0) then
     WRITE(6,'(A,I6,2A,I8)') 'ERROR: rank ',myrank, &
          ' cannot write one-particle diagnostic ',trim(filename),ios
     CALL stop_MPI
     STOP 'one-particle diagnostic write failed'
  endif
  if (myrank.eq.master .and. ikp.eq.1) WRITE(6,'(A,1X,A)') &
       ':DMFT_ONEPARTICLE_DUMP:',trim(filename)
END SUBROUTINE WriteOneParticleMomentDump


! Save the complete explicit strain derivative of the normalized DMFT
! projector.  This is the object conjugate to projectorGradient in
!
!        dGamma_projector/d epsilon_q
!          = 2 Re sum(projectorGradient * projectorStrain_q).
!
! The six components use the dmft2 stress ordering xx, yy, zz, yz, xz, xy.
! The first three can consequently be summed to form a hydrostatic derivative
! and compared directly with independently strained endpoint projectors.
! Nothing is written unless EDMFT_STRESS_PROJECTOR_DIAG is explicitly set.
SUBROUTINE WriteProjectorStrainDiagnostic(projectorStrain,bandEnergy,ikp, &
     nbands,maxdim2,norbitals)
  USE com_mpi,ONLY: myrank,master,stop_MPI
  IMPLICIT NONE
  INTEGER,intent(in) :: ikp,nbands,maxdim2,norbitals
  COMPLEX*16,intent(in) :: projectorStrain(nbands,maxdim2,norbitals,6)
  REAL*8,intent(in) :: bandEnergy(nbands)
  INTEGER,parameter :: strainMagic=20261106,strainVersion=1,nstrain=6
  CHARACTER(len=64) :: environment
  CHARACTER(len=256) :: filename
  INTEGER :: status,fh,ios

  environment=''
  CALL GET_ENVIRONMENT_VARIABLE('EDMFT_STRESS_PROJECTOR_DIAG', &
       environment,STATUS=status)
  if (status.ne.0 .or. len_trim(environment).eq.0 .or. &
       trim(environment).eq.'0') return

  WRITE(filename,'(A,I8.8,A)') 'dmft_projector_strain_k',ikp,'.bin'
  fh=973
  open(fh,file=trim(filename),status='replace',form='unformatted', &
       access='stream',action='write',iostat=ios)
  if (ios.eq.0) write(fh,iostat=ios) strainMagic,strainVersion,ikp, &
       nbands,maxdim2,norbitals,nstrain
  if (ios.eq.0) write(fh,iostat=ios) bandEnergy
  if (ios.eq.0) write(fh,iostat=ios) projectorStrain
  close(fh)
  if (ios.ne.0) then
     WRITE(6,'(A,I6,2A,I8)') 'ERROR: rank ',myrank, &
          ' cannot write projector-strain diagnostic ',trim(filename),ios
     CALL stop_MPI
     STOP 'projector-strain diagnostic write failed'
  endif
  if (myrank.eq.master .and. ikp.eq.1) WRITE(6,'(A,1X,ES22.13)') &
       ':DMFT_PROJECTOR_STRAIN_DUMP:',sqrt(sum(abs(projectorStrain)**2))
END SUBROUTINE WriteProjectorStrainDiagnostic


!===============================================================================
! NORMALIZATION / FIXED-WINDOW CROSS RESPONSE
!===============================================================================
!
! The normalized projector is U_k=P_k X with
!
!                 X=S^(-1/2),   S=sum_k w_k P_k^+ P_k.
!
! BuildWindowResponseMatrices includes the complementary-band response of P_k
! while holding X fixed.  ComputeNormalizationDerivative includes X_,q caused
! by the explicit LAPW derivative P_,q.  The remaining chain-rule branch is
!
!     wave-function response -> P_k response -> S response -> X response.
!
! A literal implementation would construct S_,q for all six strains.  That
! would require exposing every off-diagonal H_,q-epsilon O_,q matrix element.
! The adjoint implementation below is equivalent and reuses the existing LAPW
! contractions.  If the stationary-functional projector variation is
!
!             delta Phi = 2 Re sum_k Tr(C_k^+ delta U_k),
!
! first accumulate N=sum_k U_k^+ C_k.  The code stores conjugated gradients,
! because its direct contraction is sum(projectorGradient*projector_strain),
! so normalizationMoment below is conjg(N).
!
! For B=X^(-1) delta X, ComputeNormalizationDerivative implements, in the
! eigenbasis X=V diag(x_i) V^+,
!
!       B_ij = -x_j deltaS_ij/(1/x_i+1/x_j).
!
! Applying the adjoint of this map to N and retaining its Hermitian part gives
! Q such that 2 Re Tr(N^+ B)=2 Re Tr(Q^+ delta S).  Since
!
!       delta S = sum_k w_k(delta P_k^+ P_k+P_k^+ delta P_k),
!
! the extra gradient with respect to a raw projector is
!
!                       Cnorm_k = 2 w_k P_k Q.
!
! A second vector-only traversal uses Cnorm_k in exactly the same finite-band
! response construction as BuildWindowResponseMatrices.  It repeats neither
! the Matsubara sum nor the density calculation.
!===============================================================================

SUBROUTINE AccumulateNormalizationMoment(moment,Ufull,gradient,iw0,nbands, &
     nbands_full,maxdim2,norbitals,maxdim,ncix,cix_orb,cixdim,nindo,iSx)
  IMPLICIT NONE
  INTEGER,intent(in) :: iw0,nbands,nbands_full,maxdim2,norbitals,maxdim,ncix
  INTEGER,intent(in) :: cix_orb(norbitals),cixdim(ncix),nindo(norbitals), &
       iSx(maxdim2,norbitals)
  COMPLEX*16,intent(in) :: Ufull(nbands_full,maxdim2,norbitals)
  COMPLEX*16,intent(in) :: gradient(nbands,maxdim2,norbitals)
  COMPLEX*16,intent(inout) :: moment(maxdim,maxdim,ncix)
  INTEGER :: icix,cixdm,iorb,ind
  COMPLEX*16,allocatable :: Ucix(:,:),Gbar(:,:)

  do icix=1,ncix
     cixdm=cixdim(icix)
     allocate(Ucix(nbands,cixdm),Gbar(nbands,cixdm))
     Ucix=0.d0; Gbar=0.d0
     do iorb=1,norbitals
        if (cix_orb(iorb).ne.icix) cycle
        do ind=1,nindo(iorb)
           Ucix(:,iSx(ind,iorb))= &
                Ufull(iw0:iw0+nbands-1,ind,iorb)
           Gbar(:,iSx(ind,iorb))=gradient(:,ind,iorb)
        enddo
     enddo
     ! Gbar=conjg(C).  Therefore transpose(Ucix)*Gbar=conjg(Ucix^+ C).
     moment(1:cixdm,1:cixdm,icix)= &
          moment(1:cixdm,1:cixdm,icix)+matmul(transpose(Ucix),Gbar)
     deallocate(Ucix,Gbar)
  enddo
END SUBROUTINE AccumulateNormalizationMoment


! Write the conjugated normalization gradient used by the adjoint window
! response.  This is a marker-controlled regression aid: together with the
! independently finite-differenced SOlapm.dat it tests the complete scalar
! normalization chain without changing any production stress component.
SUBROUTINE WriteNormalizationMomentDiagnostic(momentBar,cixdim,maxdim,ncix)
  USE com_mpi,ONLY: myrank,master
  IMPLICIT NONE
  INTEGER,intent(in) :: maxdim,ncix,cixdim(ncix)
  COMPLEX*16,intent(in) :: momentBar(maxdim,maxdim,ncix)
  LOGICAL :: enabled
  INTEGER :: icix,cdim,i,j,fh

  INQUIRE(file='normalization_strain_fd.in',exist=enabled)
  if (.not.enabled .or. myrank.ne.master) return
  fh=996
  open(fh,file='normalization_moment.dat',status='replace',action='write')
  write(fh,'(A)') '# conjugated normalization gradient momentBar=conjg(N)'
  write(fh,'(A)') '# entries are Re Im pairs; contraction is 2 Re sum(momentBar*B)'
  write(fh,'(A,I8)') 'NCIX ',ncix
  do icix=1,ncix
     cdim=cixdim(icix)
     write(fh,'(A,2I8)') 'CIX ',icix,cdim
     do i=1,cdim
        write(fh,'(*(2ES25.16,1X))') (dble(momentBar(i,j,icix)), &
             aimag(momentBar(i,j,icix)),j=1,cdim)
     enddo
  enddo
  close(fh)
  WRITE(6,'(A)') ':DMFT_NORMALIZATION_MOMENT_DUMP: normalization_moment.dat'
END SUBROUTINE WriteNormalizationMomentDiagnostic


SUBROUTINE BuildNormalizationAdjoint(Qbar,momentBar,X,Sigind,cixdim, &
     maxdim,ncix)
  IMPLICIT NONE
  INTEGER,intent(in) :: maxdim,ncix,cixdim(ncix), &
       Sigind(maxdim,maxdim,ncix)
  COMPLEX*16,intent(in) :: momentBar(maxdim,maxdim,ncix), &
       X(maxdim,maxdim,ncix)
  ! Qbar=conjg(Q) is returned because all projector-gradient contractions in
  ! l2main use the conjugated-gradient convention.
  COMPLEX*16,intent(out) :: Qbar(maxdim,maxdim,ncix)
  INTEGER :: icix,ip,iq,na,info,lwork,lrwork,liwork
  INTEGER,allocatable :: active(:),iwork(:)
  REAL*8,allocatable :: xeig(:),rwork(:)
  REAL*8 :: lij,directTest,adjointTest,adjointScale
  LOGICAL :: projectorDiagnostic
  CHARACTER(LEN=64) :: diagnosticEnvironment
  INTEGER :: diagnosticStatus
  COMPLEX*16,allocatable :: Xvec(:,:),Nmat(:,:),Np(:,:),Q0p(:,:), &
       Qp(:,:),Qmat(:,:),work(:),dStest(:,:),dStestPrime(:,:), &
       BtestPrime(:,:),Btest(:,:),Mtest(:,:)

  Qbar=0.d0
  diagnosticEnvironment=''
  CALL GET_ENVIRONMENT_VARIABLE('EDMFT_STRESS_PROJECTOR_DIAG', &
       diagnosticEnvironment,STATUS=diagnosticStatus)
  projectorDiagnostic=diagnosticStatus.eq.0 .and. &
       len_trim(diagnosticEnvironment).gt.0 .and. &
       trim(diagnosticEnvironment).ne.'0'
  do icix=1,ncix
     allocate(active(cixdim(icix)))
     na=0
     do ip=1,cixdim(icix)
        if (any(Sigind(ip,1:cixdim(icix),icix).gt.0) .or. &
             any(Sigind(1:cixdim(icix),ip,icix).gt.0)) then
           na=na+1
           active(na)=ip
        endif
     enddo
     if (na.eq.0) then
        deallocate(active)
        cycle
     endif

     allocate(Xvec(na,na),Nmat(na,na),Np(na,na),Q0p(na,na), &
          Qp(na,na),Qmat(na,na),xeig(na))
     do ip=1,na
        do iq=1,na
           Xvec(ip,iq)=X(active(ip),active(iq),icix)
           Nmat(ip,iq)=conjg(momentBar(active(ip),active(iq),icix))
        enddo
     enddo

     lwork=4*na+na*na
     lrwork=10*na+2*na*na+1
     liwork=3+10*na
     allocate(work(lwork),rwork(lrwork),iwork(liwork))
     CALL ZHEEVD('V','U',na,Xvec,na,xeig,work,lwork,rwork,lrwork, &
          iwork,liwork,info)
     if (info.ne.0 .or. minval(xeig).le.1.d-14) then
        WRITE(6,*) 'ERROR: normalization adjoint diagonalization failed', &
             icix,info,minval(xeig)
        STOP 'BuildNormalizationAdjoint failed'
     endif

     Np=matmul(conjg(transpose(Xvec)),matmul(Nmat,Xvec))
     do ip=1,na
        do iq=1,na
           lij=-xeig(iq)/(1.d0/xeig(ip)+1.d0/xeig(iq))
           Q0p(ip,iq)=lij*Np(ip,iq)
        enddo
     enddo
     ! delta S is Hermitian.  Only the Hermitian part of the unconstrained
     ! adjoint can contribute to its real Frobenius contraction.
     Qp=0.5d0*(Q0p+conjg(transpose(Q0p)))
     Qmat=matmul(Xvec,matmul(Qp,conjg(transpose(Xvec))))
     do ip=1,na
        do iq=1,na
           Qbar(active(ip),active(iq),icix)=conjg(Qmat(ip,iq))
        enddo
     enddo

     if (projectorDiagnostic) then
        ! Verify the adjoint without reference to any numerical stress.  For
        ! a deterministic Hermitian trial dS, form B=X^{-1}dX with the same
        ! Sylvester map used by ComputeNormalizationDerivative.  The direct
        ! projector contraction and the adjoint contraction must agree:
        !
        !   2 Re sum(M*B) = 2 Re sum(Qbar*dS).
        !
        ! This catches a misplaced transpose or conjugation in the global
        ! normalization response even when X happens to be nearly diagonal.
        allocate(dStest(na,na),dStestPrime(na,na), &
             BtestPrime(na,na),Btest(na,na),Mtest(na,na))
        do ip=1,na
           do iq=1,na
              dStest(ip,iq)=cmplx(0.1d0*dble(ip+iq), &
                   0.03d0*dble(ip-iq),kind=8)
              Mtest(ip,iq)=momentBar(active(ip),active(iq),icix)
           enddo
        enddo
        dStestPrime=matmul(conjg(transpose(Xvec)), &
             matmul(dStest,Xvec))
        do ip=1,na
           do iq=1,na
              lij=-xeig(iq)/(1.d0/xeig(ip)+1.d0/xeig(iq))
              BtestPrime(ip,iq)=lij*dStestPrime(ip,iq)
           enddo
        enddo
        Btest=matmul(Xvec,matmul(BtestPrime, &
             conjg(transpose(Xvec))))
        directTest=2.d0*dble(sum(Mtest*Btest))
        adjointTest=2.d0*dble(sum(conjg(Qmat)*dStest))
        adjointScale=max(abs(directTest),abs(adjointTest),tiny(1.d0))
        WRITE(6,'(A,I4,3ES22.12)') ':DMFT_NORMALIZATION_ADJOINT: ', &
             icix,directTest,adjointTest, &
             abs(directTest-adjointTest)/adjointScale
        deallocate(dStest,dStestPrime,BtestPrime,Btest,Mtest)
     endif

     deallocate(active,Xvec,Nmat,Np,Q0p,Qp,Qmat,xeig,work,rwork,iwork)
  enddo
END SUBROUTINE BuildNormalizationAdjoint


SUBROUTINE BuildNormalizationWindowResponseMatrices(rhoResponse,eResponse, &
     Pfull,Qbar,wg,E,nemin,nemax,DM_nemin,DM_nemax,nbands_full,nbands, &
     maxdim2,norbitals,maxdim,ncix,cix_orb,cixdim,nindo,iSx)
  USE com_mpi,ONLY: Qprint
  IMPLICIT NONE
  INTEGER,intent(in) :: nemin,nemax,DM_nemin,DM_nemax,nbands_full,nbands
  INTEGER,intent(in) :: maxdim2,norbitals,maxdim,ncix
  INTEGER,intent(in) :: cix_orb(norbitals),cixdim(ncix),nindo(norbitals), &
       iSx(maxdim2,norbitals)
  REAL*8,intent(in) :: wg,E(*)
  COMPLEX*16,intent(in) :: Pfull(nbands_full,maxdim2,norbitals), &
       Qbar(maxdim,maxdim,ncix)
  COMPLEX*16,intent(out) :: rhoResponse(nbands_full,nbands_full), &
       eResponse(nbands_full,nbands_full)
  INTEGER :: iw0,iw1,icix,cixdm,iorb,ind,ia,iw,in,ndegenerate
  REAL*8 :: denominator,epsn,responseNorm
  COMPLEX*16 :: response
  COMPLEX*16,allocatable :: Pall(:,:),Pwin(:,:),Gbar(:,:),Kfull(:,:)

  rhoResponse=0.d0; eResponse=0.d0
  iw0=DM_nemin-nemin+1
  iw1=DM_nemax-nemin+1
  if (iw0.lt.1 .or. iw1.gt.nbands_full .or. iw1-iw0+1.ne.nbands) &
       STOP 'BuildNormalizationWindowResponseMatrices: band bounds'
  allocate(Kfull(nbands_full,nbands))
  Kfull=0.d0

  do icix=1,ncix
     cixdm=cixdim(icix)
     allocate(Pall(nbands_full,cixdm),Pwin(nbands,cixdm), &
          Gbar(nbands,cixdm))
     Pall=0.d0
     do iorb=1,norbitals
        if (cix_orb(iorb).ne.icix) cycle
        do ind=1,nindo(iorb)
           Pall(:,iSx(ind,iorb))=Pfull(:,ind,iorb)
        enddo
     enddo
     Pwin=Pall(iw0:iw1,:)
     ! Qbar=conjg(Q), hence conjg(Pwin)*Qbar=conjg(Pwin*Q).
     ! The factor two follows from differentiating both factors in P^+P;
     ! wg is exactly the normalization weight used in cmp_overlap.
     Gbar=2.d0*wg*matmul(conjg(Pwin),Qbar(1:cixdm,1:cixdm,icix))
     Kfull=Kfull+matmul(Pall,transpose(Gbar))
     deallocate(Pall,Pwin,Gbar)
  enddo

  ndegenerate=0
  do ia=1,nbands_full
     if (ia.ge.iw0 .and. ia.le.iw1) cycle
     do iw=1,nbands
        in=iw0+iw-1
        epsn=E(DM_nemin+iw-1)
        denominator=epsn-E(nemin+ia-1)
        if (abs(denominator).le.1.d-10) then
           ndegenerate=ndegenerate+1
           cycle
        endif
        response=Kfull(ia,iw)/denominator
        rhoResponse(ia,in)=response
        rhoResponse(in,ia)=conjg(response)
        eResponse(ia,in)=epsn*response
        eResponse(in,ia)=conjg(eResponse(ia,in))
     enddo
  enddo
  responseNorm=sqrt(sum(abs(rhoResponse)**2))
  if (Qprint .and. ndegenerate.gt.0) WRITE(6,'(A,I6,A,ES12.4)') &
       'Normalization-window response: degenerate=',ndegenerate, &
       '  ||rho_resp||=',responseNorm
  deallocate(Kfull)
END SUBROUTINE BuildNormalizationWindowResponseMatrices


SUBROUTINE cmp_normalization_window_stress(Qbar,cfX,Rspin,iorbital, &
     cix_orb,cixdim,iSx,nindo,nnlo,norbitals,natom,maxdim2,emin,emax, &
     sumw,lmaxp,time_bl,time_bl_w,time_lo,time_lo_w,max_response_modes, &
     savedWindowRho,savedWindowEnergy,max_nbands,fuse_window_response)
  USE param,ONLY: nume,nmat,lmax2
  USE com,ONLY: nat
  USE dmfts,ONLY: iso,ncix,maxdim,projector,DM_Emin,DM_Emax
  USE xa,ONLY: E,fj,dfj
  USE xa3,ONLY: bk3,bk3lo,aK,k3,k3lo,As,As_lo
  USE structure,ONLY: rmt,BR1
  USE com_mpi,ONLY: myrank,vector_para,vectors,nvector,fvectors,pr_proc,Qprint
  USE p_project,ONLY: phi_jl,dphi_jl,rix,P_rfi,kNri,n_al_ucase, &
       dri,l_al_ucase,j_al_ucase
  USE dmft_stress,ONLY: stress_val_corr,stress_kin_mt,stress_apw, &
       stress_kin_int,stress_norm_window,stress_diag_norm
  IMPLICIT NONE
  INTEGER,intent(in) :: nnlo,norbitals,natom,maxdim2,lmaxp,max_nbands
  INTEGER,intent(out) :: max_response_modes
  LOGICAL,intent(in) :: fuse_window_response
  INTEGER,intent(in) :: iorbital(natom,lmaxp+1),cix_orb(norbitals), &
       cixdim(ncix),iSx(maxdim2,norbitals),nindo(norbitals)
  REAL*8,intent(in) :: emin,emax,sumw
  REAL*8,intent(inout) :: time_bl,time_bl_w,time_lo,time_lo_w
  COMPLEX*16,intent(in) :: Qbar(maxdim,maxdim,ncix), &
       cfX(maxdim2,maxdim2,norbitals,norbitals),Rspin(2,2,norbitals), &
       savedWindowRho(nume,max_nbands,pr_proc), &
       savedWindowEnergy(nume,max_nbands,pr_proc)
  INTERFACE
     REAL*8 FUNCTION Romb(y,N,dx)
       IMPLICIT NONE
       INTEGER,intent(in) :: N
       REAL*8,intent(in) :: y(N),dx
     END FUNCTION Romb
  END INTERFACE
  INTEGER :: ivector,is,itape,iks,ikp,iikp,i,n0,nemin,nemax,DM_nemin, &
       DM_nemax,isize,jatom,Nri,iind,l,ir,nbands_full,nbands, &
       response_nemaxx,rho_response_nemaxx,energy_response_nemaxx,ia,iw,iw0
  REAL*8 :: emist,wgh0,Kn(3),rx,wg
  REAL*8 :: before_val(3,3),before_kinmt(3,3),before_kinint(3,3), &
       before_apw(3,3)
  LOGICAL :: more_kpoints,Tcompute
  REAL*8,allocatable :: aKR(:),jlr(:),jlrp(:),jlp(:),zResponse(:), &
       zEResponse(:),weightResponse(:)
  COMPLEX*16,allocatable :: Pfull(:,:,:,:),rhoResponse(:,:), &
       eResponse(:,:),AResponse(:,:),eResponseDiag(:,:)

  ! Snapshot every linear LAPW kernel separately.  In particular this exposes
  ! whether the signed normalization-response pseudo-bands are accidentally
  ! contracted with physical-band weights in any one kernel.
  before_val=stress_val_corr
  before_kinmt=stress_kin_mt
  before_kinint=stress_kin_int
  before_apw=stress_apw
  stress_norm_window=0.d0
  max_response_modes=0
  iikp=0
  do ivector=1,nvector
     do is=1,iso
        itape=8+is
        if (vector_para) then
           open(itape,FILE=fvectors(ivector,is),STATUS='old',FORM='unformatted')
        else
           rewind(itape)
        endif
        do i=1,nat
           read(itape) emist
           read(itape) emist
        enddo
     enddo

     do iks=1,vectors(ivector,2)
        if (vector_para) then
           ikp=vectors(ivector,3)+iks
           iikp=iikp+1
        else
           ikp=iks
           iikp=ikp-myrank*pr_proc
        endif
        Tcompute=(iikp.gt.0)
        if (iikp.gt.pr_proc) exit
        CALL Read_Vec_Spin(E,wgh0,As,As_lo,k3,k3lo,bk3,bk3lo,nemin, &
             nemax,DM_nemin,DM_nemax,more_kpoints,n0,emin,emax, &
             DM_Emin,DM_Emax,iso,nmat,nume,nnlo)
        if (.not.Tcompute) cycle
        if (.not.more_kpoints) exit

        isize=n0-nnlo
        do i=1,isize
           Kn=matmul(BR1,bk3(:,i))
           aK(i)=sqrt(dot_product(Kn,Kn))
        enddo
        do jatom=1,nat
           CALL HARMON(isize,aK(:isize),lmax2,fj(:,:isize,jatom), &
                dfj(:,:isize,jatom),rmt(jatom))
        enddo

        if (abs(projector).ge.5) then
           allocate(phi_jl(nmat,n_al_ucase),dphi_jl(nmat,n_al_ucase))
           phi_jl=0.d0; dphi_jl=0.d0
           Nri=2**kNri+1
           allocate(aKR(Nri),jlr(Nri),jlrp(Nri),jlp(Nri))
           do iind=1,n_al_ucase
              if (abs(dri(iind)).le.1.d-10) cycle
              l=l_al_ucase(iind); jatom=j_al_ucase(iind)
              do i=1,isize
                 rx=rmt(jatom)
                 do ir=1,Nri
                    aKR(ir)=rx*aK(i); rix(ir)=rx; rx=rx+dri(iind)
                 enddo
                 CALL sphbes2(l,Nri,aKR,jlr)
                 jlp=jlr*rix*P_rfi(:,iind)
                 phi_jl(i,iind)=romb(jlp,Nri,dri(iind))
              enddo
           enddo
           deallocate(aKR,jlr,jlrp,jlp)
        endif

        nbands_full=nemax-nemin+1
        nbands=DM_nemax-DM_nemin+1
        allocate(Pfull(nbands_full,maxdim2,norbitals,1))
        CALL Build_DMFT_Projector(Pfull,cfX,Rspin,iorbital,norbitals,1, &
             n0,nnlo,nbands_full,cix_orb,nindo,nemin,nemax,maxdim2, &
             lmaxp,.false.)
        if (abs(projector).ge.5) deallocate(phi_jl,dphi_jl)

        allocate(rhoResponse(nbands_full,nbands_full), &
             eResponse(nbands_full,nbands_full))
        wg=wgh0/sumw
        CALL BuildNormalizationWindowResponseMatrices(rhoResponse, &
             eResponse,Pfull(:,:,:,1),Qbar,wg,E,nemin,nemax,DM_nemin, &
             DM_nemax,nbands_full,nbands,maxdim2,norbitals,maxdim,ncix, &
             cix_orb,cixdim,nindo,iSx)
        deallocate(Pfull)

        if (fuse_window_response) then
           ! Both complementary-space responses are linear density matrices
           ! with nonzero window--complement blocks.  Add them before the
           ! Hermitian diagonalization so the expensive LAPW stress kernels
           ! see one rank-at-most-2*Nwindow object instead of two.  The saved
           ! object contains only its window columns; restore the conjugate
           ! block explicitly to preserve Hermiticity to roundoff.
           iw0=DM_nemin-nemin+1
           do iw=1,nbands
              do ia=1,nbands_full
                 if (ia.ge.iw0 .and. ia.lt.iw0+nbands) cycle
                 rhoResponse(ia,iw0+iw-1)=rhoResponse(ia,iw0+iw-1)+ &
                      savedWindowRho(ia,iw,iikp)
                 rhoResponse(iw0+iw-1,ia)= &
                      conjg(rhoResponse(ia,iw0+iw-1))
                 eResponse(ia,iw0+iw-1)=eResponse(ia,iw0+iw-1)+ &
                      savedWindowEnergy(ia,iw,iikp)
                 eResponse(iw0+iw-1,ia)= &
                      conjg(eResponse(ia,iw0+iw-1))
              enddo
           enddo
        endif

        if (maxval(abs(rhoResponse)).gt.1.d-14 .or. &
             maxval(abs(eResponse)).gt.1.d-14) then
           allocate(AResponse(nbands_full,nbands_full), &
                eResponseDiag(nbands_full,nbands_full), &
                zResponse(nbands_full),zEResponse(nbands_full), &
                weightResponse(nume))
           AResponse=rhoResponse; eResponseDiag=eResponse
           rho_response_nemaxx=nemax
           CALL Diagonalize_DMFT_WEIGHTS(zResponse,AResponse,nbands_full, &
                nemin,rho_response_nemaxx,.true., &
                'normalization_window_density_response',ikp)
           energy_response_nemaxx=nemax
           CALL Diagonalize_DMFT_WEIGHTS(zEResponse,eResponseDiag, &
                nbands_full,nemin,energy_response_nemaxx,.true., &
                'normalization_window_energy_response',ikp)
           ! The response matrices contain the two conjugate cross-window
           ! blocks explicitly, so their trace has already formed the 2*Re
           ! projector contraction.  CMP_DMFT_STRESS subsequently uses LAPW
           ! quadratic-form kernels with another factor two.  Halve only the
           ! signed pseudo-band weights to avoid double counting.  This is
           ! also correct when FUSE_WINDOW_RESPONSE is active, because both
           ! fixed-X and normalization responses are combined before this
           ! single diagonalization.
           zResponse(:rho_response_nemaxx-nemin+1)=0.5d0* &
                zResponse(:rho_response_nemaxx-nemin+1)
           zEResponse(:energy_response_nemaxx-nemin+1)=0.5d0* &
                zEResponse(:energy_response_nemaxx-nemin+1)
           response_nemaxx=max(rho_response_nemaxx,energy_response_nemaxx)
           if (rho_response_nemaxx.lt.nemax) &
                zResponse(rho_response_nemaxx-nemin+2:)=0.d0
           if (energy_response_nemaxx.lt.nemax) &
                zEResponse(energy_response_nemaxx-nemin+2:)=0.d0
           weightResponse=0.d0
           weightResponse(nemin:rho_response_nemaxx)= &
                zResponse(:rho_response_nemaxx-nemin+1)
           max_response_modes=max(max_response_modes,response_nemaxx-nemin+1)
           CALL cmp_dmft_stress(AResponse,eResponseDiag,weightResponse, &
                zEResponse,nemin,nemax,nemin,nemax,response_nemaxx, &
                nbands_full,nbands_full,n0,nnlo,time_bl,time_bl_w, &
                time_lo,time_lo_w)
           deallocate(AResponse,eResponseDiag,zResponse,zEResponse, &
                weightResponse)
        endif
        deallocate(rhoResponse,eResponse)
        if (Qprint) WRITE(6,'(A,I5)') &
             ' Finished normalization-window stress k-point ',ikp
     enddo

     do is=1,iso
        itape=8+is
        if (vector_para) then
           close(itape)
        else
           rewind(itape)
        endif
     enddo
  enddo

  stress_diag_norm(:,:,1)=stress_val_corr-before_val
  stress_diag_norm(:,:,2)=stress_kin_mt-before_kinmt
  stress_diag_norm(:,:,3)=stress_kin_int-before_kinint
  stress_diag_norm(:,:,4)=stress_apw-before_apw
  stress_norm_window=sum(stress_diag_norm,dim=3)
END SUBROUTINE cmp_normalization_window_stress


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


! Evaluate only the k-points owned by this rank.  The chemical-potential search
! combines these scalar contributions; it never needs the eigenvalue arrays on
! the master rank (except in the explicit mode='e' band-data export).
REAL*8 FUNCTION Density_local(EF)
  USE defs, ONLY: PI, IMAG
  USE muzero,ONLY: w_sum, w_gamma, w_norm1, nkp, w_nomega, w_npomega, &
       wgh, abom, zEk, nemm, wmatsubara, n0_om, w_omega, jomq, &
       w_beta, matsum_weight, wprint1, Qcheckbands
  USE dmfts, ONLY: Temperature
  IMPLICIT NONE
  REAL*8, intent(in)    :: EF
  ! functions
  REAL*8 :: ferm
  ! locals
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
           
           ! The sparse tail is summed with the precomputed discrete cubic-
           ! spline weights.  Interpolate the final resolvent integrand, not
           ! its complex pole/eigenvalue: nonlinear interpolation of the pole
           ! can change the answer under a band rotation and cannot share one
           ! set of weights with the matrix density and stress moments.
           do iom=n0_om,w_nomega
              omn = dcmplx(0.d0,w_omega(iom))
              cek = zEk(num0,iom,ikp)-EF
              csum = csum + matsum_weight(iom)* &
                   (1./(omn-cek)-1./(omn-ek0))
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
  Density_local = zw0
  RETURN
END FUNCTION Density_local

! Brent and the gap-candidate search run on the master.  Each trial EF is
! broadcast to the workers, which evaluate their local bands in
! Serve_Density_requests.  Only one real occupancy per rank is reduced.
REAL*8 FUNCTION Density(EF)
  USE com, ONLY: elecn
  USE muzero, ONLY: wprint
  USE com_mpi, ONLY: Bcast_Density_request, Reduce_Density_sum
  IMPLICIT NONE
  REAL*8, intent(in) :: EF
  REAL*8, PARAMETER :: Ry2eV=13.60569253d0
  REAL*8 :: trial_EF, local_charge, total_charge, Density_local
  LOGICAL :: evaluate

  trial_EF = EF
  evaluate = .TRUE.
  CALL Bcast_Density_request(trial_EF, evaluate)
  local_charge = Density_local(trial_EF)
  CALL Reduce_Density_sum(local_charge, total_charge)
  if (wprint) WRITE(*,'(A,f13.8,A,f14.8,A,f10.4,A,f14.8)') &
       'EF[eV]=', EF*Ry2eV, ' Density=', total_charge, ' N=', elecn, ' EF[Ry]=', EF
  Density = total_charge-elecn
END FUNCTION Density

SUBROUTINE Serve_Density_requests()
  USE com_mpi, ONLY: Bcast_Density_request, Reduce_Density_sum
  IMPLICIT NONE
  REAL*8 :: trial_EF, local_charge, total_charge, Density_local
  LOGICAL :: evaluate

  do
     trial_EF = 0.d0
     evaluate = .FALSE.
     CALL Bcast_Density_request(trial_EF, evaluate)
     if (.NOT.evaluate) exit
     local_charge = Density_local(trial_EF)
     CALL Reduce_Density_sum(local_charge, total_charge)
  enddo
END SUBROUTINE Serve_Density_requests

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
  LOGICAL     :: more_kpoints, Tcompute, density_request
  INTEGER     :: jatom, isize, nbands, itape, is, i, iorb, icix, ip, it, iom, itmax, fh_EF, iikp, ind, iband
  REAL*8      :: tEF, t1c, t2c, t1w, t2w, dd, dd2, dEF, Ea, Eb, fatol, xatol, xrtol, twgh, dsm, wkp, zw0, EFnew, rx, density_request_EF
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
  
  Qcheckbands=.FALSE. ! .TRUE.  ! 2026 We changed this because Matsubara sum can mix different green's functions of bands, and hence clipping of density
  ! of eigenvalues can be problematic. The issue was found by Benjamin Bacq.
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
        CALL Build_DMFT_Projector(DMFTU,cfX,Rspin,iorbital,norbitals, &
             nip,n0,nnlo,nbands,cix_orb,nindo,DM_nemin,DM_nemax, &
             maxdim2,lmaxp,.False.)

        if (abs(projector).ge.5) deallocate( phi_jl )
        
        if (Qrenormalize) then
           if (abs(projector).le.5) then
              CALL RenormalizeTrans(DMFTU, Olapm, SOlapm, cix_orb, cixdim, nindo, iSx, nbands, nip, maxdim2, norbitals, maxdim, ncix, SIMPLE)
           else
              CALL RenormalizeTransK(DMFTU, cix_orb, cixdim, nindo, iSx, Sigind, projector, nbands, nip, maxdim2, norbitals, maxdim, ncix)
           endif
        endif


        ! Keep the chemical-potential pass consistent with the main dmft2
        ! pass when the opt-in frozen-projector finite-difference diagnostic
        ! is active.  Freezing U only after cmp_EF would leave an implicit
        ! projector dependence through the chemical potential and would not
        ! cleanly isolate the projector derivative.
        CALL SnapshotDMFTProjector(DMFTU, &
             E(DM_nemin:DM_nemax),nbands,maxdim2,norbitals,nip, &
             ikp,iso,'EF')


        
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
  
  ! The ordinary EF search sums one local occupancy per rank for each trial.
  ! Only mode='e' writes the full frequency-dependent eigenvalue dataset and
  ! therefore still needs to gather zEk, Ek, weights, and band bounds.
  if (mode.EQ.'e') CALL Reduce_MPI1(pr_procr, nomega, max_nume)
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
     
  ! All ranks must use the same Matsubara/real-axis occupancy formula.
  wmatsubara = matsubara
  !wprint1=.false.
  if (myrank.EQ.master) then
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
     ! Release workers waiting for another Density trial before the EF bcast.
     density_request_EF = 0.d0
     density_request = .FALSE.
     CALL Bcast_Density_request(density_request_EF, density_request)
  else
     CALL Serve_Density_requests()
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


! Save or restore the fully normalized DMFT projector for one k point.
!
! This is a deliberately opt-in regression diagnostic controlled by
!
!   EDMFT_PROJECTOR_SNAPSHOT=WRITE
!   EDMFT_PROJECTOR_SNAPSHOT=READ
!
! WRITE is used once at the reference structure.  READ is then used at the
! plus/minus strained structures after copying the snapshot files from the
! reference directory.  Files are keyed by stage, spin and the global k-point
! index, rather than by MPI rank.  Consequently the diagnostic is independent
! of the MPI distribution and can even be read with a different process
! count.  The EF and MAIN stages are stored separately because they are two
! independent projector constructions and may in general use different array
! extents.
!
! The complete allocated array is stored instead of only active correlated
! columns.  This makes the read operation an exact replacement and ensures
! that STrans, lgTrans, the local observables and the density all see precisely
! the same reference projector.  In the normal (environment variable unset)
! path this routine returns before any file operation or array traversal.
SUBROUTINE SnapshotDMFTProjector(DMFTU,bandEnergy,nbands,maxdim2, &
     norbitals,nipc,ikp,ispin,stage)
  USE com_mpi, ONLY: myrank,master,stop_MPI
  IMPLICIT NONE
  INTEGER,intent(in) :: nbands,maxdim2,norbitals,nipc,ikp,ispin
  CHARACTER(len=*),intent(in) :: stage
  COMPLEX*16,intent(inout) :: DMFTU(nbands,maxdim2,norbitals,nipc)
  REAL*8,intent(in) :: bandEnergy(nbands)
  INTEGER,parameter :: snapshotMagic=20261104,snapshotVersion=2
  REAL*8,parameter :: defaultDegeneracyTolerance=1.d-8
  CHARACTER(len=32) :: mode
  CHARACTER(len=256) :: filename
  INTEGER :: status,letter,i,j,fh,ios,firstBand,lastBand,blockSize
  INTEGER :: ii,jj,ind,iorb,ip,info,lwork
  INTEGER :: fileMagic,fileVersion,fileIkp,fileSpin
  INTEGER :: fileBands,fileDim,fileOrbitals,fileNipc
  LOGICAL :: exists
  REAL*8 :: projectorNorm,gaugeResidual,gaugeScale,overlapScale
  REAL*8 :: degeneracyTolerance
  CHARACTER(len=64) :: toleranceText
  COMPLEX*16 :: phase
  COMPLEX*16,allocatable :: currentU(:,:,:,:),cross(:,:), &
       leftVectors(:,:),rightVectorsH(:,:),unitary(:,:),work(:), &
       alignedBlock(:,:,:,:)
  REAL*8,allocatable :: referenceEnergy(:),singularValues(:),rwork(:)

  mode=''
  CALL GET_ENVIRONMENT_VARIABLE('EDMFT_PROJECTOR_SNAPSHOT',mode, &
       STATUS=status)
  if (status.ne.0 .or. len_trim(mode).eq.0) return
  mode=adjustl(mode)
  do i=1,len_trim(mode)
     letter=iachar(mode(i:i))
     if (letter.ge.iachar('a') .and. letter.le.iachar('z')) &
          mode(i:i)=achar(letter-iachar('a')+iachar('A'))
  enddo
  if (trim(mode).ne.'WRITE' .and. trim(mode).ne.'READ') then
     if (myrank.eq.master) WRITE(6,'(A,A)') &
          'ERROR: EDMFT_PROJECTOR_SNAPSHOT must be WRITE or READ, not ', &
          trim(mode)
     CALL stop_MPI
     STOP 'invalid EDMFT_PROJECTOR_SNAPSHOT mode'
  endif

  ! Independently diagonalized strained structures can rotate bands whose
  ! splittings are small but not exactly zero.  A single hard-coded tolerance
  ! therefore makes a projector finite difference depend on the eigensolver's
  ! arbitrary gauge.  Keep the historical value as the production default,
  ! but permit a diagnostic sweep without recompiling dmft2.
  degeneracyTolerance=defaultDegeneracyTolerance
  toleranceText=''
  CALL GET_ENVIRONMENT_VARIABLE( &
       'EDMFT_PROJECTOR_SNAPSHOT_DEGENERACY_TOL',toleranceText, &
       STATUS=status)
  if (status.eq.0 .and. len_trim(toleranceText).gt.0) then
     read(toleranceText,*,iostat=ios) degeneracyTolerance
     if (ios.ne.0 .or. degeneracyTolerance.le.0.d0) then
        if (myrank.eq.master) WRITE(6,'(2A)') &
             'ERROR: invalid EDMFT_PROJECTOR_SNAPSHOT_DEGENERACY_TOL=', &
             trim(toleranceText)
        CALL stop_MPI
        STOP 'invalid projector snapshot degeneracy tolerance'
     endif
  endif
  if (myrank.eq.master .and. ikp.eq.1 .and. ispin.eq.1) &
       WRITE(6,'(A,1X,ES22.13)') &
       ':DMFT_PROJECTOR_SNAPSHOT_DEGENERACY_TOL:',degeneracyTolerance

  WRITE(filename,'(A,A,A,I2.2,A,I8.8,A)') &
       'dmft_projector_snapshot_',trim(stage),'_s',ispin,'_k',ikp,'.bin'
  fh=971
  allocate(referenceEnergy(nbands))

  if (trim(mode).eq.'WRITE') then
     open(fh,file=trim(filename),status='replace',form='unformatted', &
          access='stream',action='write',iostat=ios)
     if (ios.ne.0) then
        WRITE(6,'(A,I6,2A,I8)') 'ERROR: rank ',myrank, &
             ' cannot create projector snapshot ',trim(filename),ios
        CALL stop_MPI
        STOP 'projector snapshot open failed'
     endif
     write(fh,iostat=ios) snapshotMagic,snapshotVersion,ikp,ispin, &
          nbands,maxdim2,norbitals,nipc
     if (ios.eq.0) write(fh,iostat=ios) bandEnergy
     if (ios.eq.0) write(fh,iostat=ios) DMFTU
     close(fh)
     if (ios.ne.0) then
        WRITE(6,'(A,I6,2A,I8)') 'ERROR: rank ',myrank, &
             ' cannot write projector snapshot ',trim(filename),ios
        CALL stop_MPI
        STOP 'projector snapshot write failed'
     endif
  else
     ! Preserve the endpoint projector before loading the reference.  It is
     ! used only to infer the arbitrary band eigenvector gauge; its amplitudes
     ! are not mixed into the frozen reference projector.
     allocate(currentU(nbands,maxdim2,norbitals,nipc))
     currentU=DMFTU
     inquire(file=trim(filename),exist=exists)
     if (.not.exists) then
        WRITE(6,'(A,I6,2A)') 'ERROR: rank ',myrank, &
             ' cannot find projector snapshot ',trim(filename)
        CALL stop_MPI
        STOP 'projector snapshot is missing'
     endif
     open(fh,file=trim(filename),status='old',form='unformatted', &
          access='stream',action='read',iostat=ios)
     if (ios.eq.0) read(fh,iostat=ios) fileMagic,fileVersion,fileIkp, &
          fileSpin,fileBands,fileDim,fileOrbitals,fileNipc
     if (ios.ne.0) then
        close(fh)
        WRITE(6,'(A,I6,2A,I8)') 'ERROR: rank ',myrank, &
             ' cannot read projector snapshot header ',trim(filename),ios
        CALL stop_MPI
        STOP 'projector snapshot header read failed'
     endif
     if (fileMagic.ne.snapshotMagic .or. &
          fileVersion.ne.snapshotVersion .or. fileIkp.ne.ikp .or. &
          fileSpin.ne.ispin .or. fileBands.ne.nbands .or. &
          fileDim.ne.maxdim2 .or. fileOrbitals.ne.norbitals .or. &
          fileNipc.ne.nipc) then
        close(fh)
        WRITE(6,'(A,I6,2A)') 'ERROR: rank ',myrank, &
             ' found incompatible projector snapshot ',trim(filename)
        WRITE(6,'(A,8I10)') ' expected: ',snapshotMagic, &
             snapshotVersion,ikp,ispin,nbands,maxdim2,norbitals,nipc
        WRITE(6,'(A,8I10)') ' found:    ',fileMagic,fileVersion, &
             fileIkp,fileSpin,fileBands,fileDim,fileOrbitals,fileNipc
        CALL stop_MPI
        STOP 'incompatible projector snapshot'
     endif
     read(fh,iostat=ios) referenceEnergy
     read(fh,iostat=ios) DMFTU
     close(fh)
     if (ios.ne.0) then
        WRITE(6,'(A,I6,2A,I8)') 'ERROR: rank ',myrank, &
             ' cannot read projector snapshot data ',trim(filename),ios
        CALL stop_MPI
        STOP 'projector snapshot data read failed'
     endif

     ! Express the frozen reference projector in the endpoint eigenvector
     ! gauge.  For an isolated band this is the usual phase alignment.  Exact
     ! degeneracies permit an arbitrary unitary rotation, so for each block
     ! which is degenerate in both structures we solve the unitary Procrustes
     ! problem
     !
     !        min_W || W U_reference - U_endpoint ||_F ,  W^+W=1.
     !
     ! If C=U_endpoint U_reference^+=L s R^+, the minimizing rotation is
     ! W=L R^+.  Applying W to every projector channel changes only the band
     ! gauge.  The reference singular values/amplitudes remain frozen.
     firstBand=1
     do while (firstBand.le.nbands)
        lastBand=firstBand
        do while (lastBand.lt.nbands)
           if (abs(bandEnergy(lastBand+1)-bandEnergy(lastBand)).le. &
                degeneracyTolerance .and. &
                abs(referenceEnergy(lastBand+1)- &
                referenceEnergy(lastBand)).le.degeneracyTolerance) then
              lastBand=lastBand+1
           else
              exit
           endif
        enddo
        blockSize=lastBand-firstBand+1

        if (blockSize.eq.1) then
           phase=sum(conjg(DMFTU(firstBand,:,:,1))* &
                currentU(firstBand,:,:,1))
           overlapScale=sqrt(sum(abs(DMFTU(firstBand,:,:,1))**2)* &
                sum(abs(currentU(firstBand,:,:,1))**2))
           if (abs(phase).gt.1.d-12*max(overlapScale,tiny(1.d0))) &
                DMFTU(firstBand,:,:,:)=DMFTU(firstBand,:,:,:)* &
                phase/abs(phase)
        else
           allocate(cross(blockSize,blockSize), &
                leftVectors(blockSize,blockSize), &
                rightVectorsH(blockSize,blockSize), &
                unitary(blockSize,blockSize), &
                singularValues(blockSize),rwork(5*blockSize))
           lwork=max(1,4*blockSize)
           allocate(work(lwork))
           do j=1,blockSize
              jj=firstBand+j-1
              do i=1,blockSize
                 ii=firstBand+i-1
                 cross(i,j)=sum(currentU(ii,:,:,1)* &
                      conjg(DMFTU(jj,:,:,1)))
              enddo
           enddo
           CALL ZGESVD('A','A',blockSize,blockSize,cross,blockSize, &
                singularValues,leftVectors,blockSize,rightVectorsH, &
                blockSize,work,lwork,rwork,info)
           if (info.ne.0) then
              WRITE(6,'(A,4I8)') &
                   'ERROR: projector gauge-alignment SVD failed: ', &
                   ikp,firstBand,lastBand,info
              CALL stop_MPI
              STOP 'projector gauge alignment failed'
           endif
           unitary=matmul(leftVectors,rightVectorsH)
           allocate(alignedBlock(blockSize,maxdim2,norbitals,nipc))
           alignedBlock=0.d0
           do ip=1,nipc
              do iorb=1,norbitals
                 do ind=1,maxdim2
                    do i=1,blockSize
                       do j=1,blockSize
                          alignedBlock(i,ind,iorb,ip)= &
                               alignedBlock(i,ind,iorb,ip)+unitary(i,j)* &
                               DMFTU(firstBand+j-1,ind,iorb,ip)
                       enddo
                    enddo
                 enddo
              enddo
           enddo
           DMFTU(firstBand:lastBand,:,:,:)=alignedBlock
           deallocate(alignedBlock,work,rwork,singularValues,unitary, &
                rightVectorsH,leftVectors,cross)
        endif
        firstBand=lastBand+1
     enddo

     gaugeResidual=sqrt(sum(abs(DMFTU(:,:,:,1)- &
          currentU(:,:,:,1))**2))
     gaugeScale=max(sqrt(sum(abs(currentU(:,:,:,1))**2)),tiny(1.d0))
     if (myrank.eq.master .and. ikp.eq.1 .and. ispin.eq.1) &
          WRITE(6,'(A,1X,A,1X,ES22.13)') &
          ':DMFT_PROJECTOR_SNAPSHOT_GAUGE:',trim(stage), &
          gaugeResidual/gaugeScale
     deallocate(currentU)
  endif

  ! Print only the first k point from the master rank.  The norm provides a
  ! compact confirmation that WRITE and READ used byte-compatible data while
  ! avoiding one line of output for every k point.
  if (myrank.eq.master .and. ikp.eq.1 .and. ispin.eq.1) then
     projectorNorm=sqrt(sum(abs(DMFTU)**2))
     WRITE(6,'(A,1X,A,1X,A,1X,ES22.13)') &
          ':DMFT_PROJECTOR_SNAPSHOT:',trim(mode),trim(stage),projectorNorm
  endif
  deallocate(referenceEnergy)
END SUBROUTINE SnapshotDMFTProjector
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
  USE muzero, ONLY: matsum_weight
  USE splines, ONLY: spline3
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
  INTEGER :: i,j,k,iom, inode, ntail, ntarget!,ij,ic
  REAL*8, ALLOCATABLE :: splineNode(:), targetPoint(:), cardinalValue(:), &
       interpolatedValue(:)
  ! jomq(iom)   :: iom is index in logarithmic mesh, and jomq(iom) gives the true sequential number in Matsubara mesh                                                                                                                                                                                                                                                
  ! nomq(iom)   :: how many Matsubara points will need this point from logarithmic mesh to compute value                                                                                                                                                                                                                                                             
  ! iomq(iom,:) :: all Matsubara points                                                                                                                                                                                                                                                                                                                              
  ! womq(iom,:) :: weight                                                                                                                                                                                                                                                                                                                                            
  nomq(:)=0; iomq(:,:)=0; jomq(:)=0; womq(:,:)=0
  matsum_weight(:)=0.d0
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

     ! The old logarithmic-tail quadrature distributed every omitted
     ! Matsubara point linearly between its two neighbouring retained nodes.
     ! A cubic spline is also a linear function of the retained values, so
     ! the complete discrete sum can instead be written once and for all as
     !
     !       sum_n F(iw_n) = sum_j matsum_weight(j) F(iw_j).
     !
     ! Build these cardinal-spline weights once for the frequency mesh.  The
     ! first n0_om-1 frequencies are consecutive and retain unit weight.  On
     ! the sparse tail we spline in x=1/(2*n-1), not in n.  This choice is
     ! essential: every causal high-frequency expansion is a polynomial in
     ! 1/(i*omega), whereas a cubic in n has no asymptotic justification and
     ! becomes inaccurate when the last logarithmic intervals are wide.  We
     ! evaluate each cardinal cubic at every actual Matsubara point and sum
     ! those values.  The same scalar weights can then be reused for scalar,
     ! vector and matrix integrands.  Quantity-specific high-frequency
     ! moments are subtracted before this quadrature and restored analytically
     ! in the consuming routines.
     if (n0_om.gt.1) matsum_weight(:n0_om-1)=1.d0
     ntail=nomega-n0_om+1
     ntarget=jomq(nomega)-jomq(n0_om)+1
     if (ntail.eq.1) then
        matsum_weight(n0_om)=dble(ntarget)
     else if (ntail.lt.4) then
        ! spline3 uses four points to construct its consistent end
        ! conditions.  Very short artificial/test meshes therefore retain
        ! the already constructed piecewise-linear cardinal weights.
        do iom=n0_om,nomega
           matsum_weight(iom)=sum(womq(iom,1:nomq(iom)))
        enddo
     else
        allocate(splineNode(ntail),targetPoint(ntarget), &
             cardinalValue(ntail),interpolatedValue(ntarget))
        ! spline3 requires an ascending grid.  Since reciprocal frequency
        ! decreases as the Matsubara index increases, store both the retained
        ! and target grids in reverse frequency order.  inode is mapped back
        ! to the original retained-frequency index when the weight is saved.
        do i=1,ntail
           splineNode(i)=1.d0/dble(2*jomq(nomega-i+1)-1)
        enddo
        do i=1,ntarget
           targetPoint(i)=1.d0/ &
                dble(2*(jomq(nomega)-i+1)-1)
        enddo
        do inode=1,ntail
           cardinalValue(:)=0.d0
           cardinalValue(inode)=1.d0
           interpolatedValue=spline3(splineNode,cardinalValue,targetPoint)
           matsum_weight(nomega-inode+1)=sum(interpolatedValue)
        enddo
        deallocate(splineNode,targetPoint,cardinalValue,interpolatedValue)
     endif
     WRITE(6,'(A,2I8,4(1X,ES16.8))') ':MATSUBARA_CUBIC_WEIGHTS:', &
          nomega,jomq(nomega),sum(matsum_weight), &
          dble(jomq(nomega)),minval(matsum_weight(n0_om:)), &
          maxval(matsum_weight(n0_om:))

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
