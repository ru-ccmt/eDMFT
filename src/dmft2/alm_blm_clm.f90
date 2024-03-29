module alm_blm_clm
contains
  !SUBROUTINE cmp_alm_blm_clm(alm,blm,clm,aalm,bblm,cclm,Aweight,latom,lfirst,jatom,nemin,nemax,n0,nnlo,is,DM_nemin,DM_nemax,DM_nemaxx,nbands,nbandsx,Qforce_j,lmax2lmax2_,nbands_dft_,dim_,lomaxlomax_,nloat_,time_bl,time_bl_w,time_lo,time_lo_w)
  SUBROUTINE cmp_alm_blm_clm(alm,blm,clm,aalm,bblm,cclm,Aweight,latom,lfirst,jatom,nemin,nemax,n0,nnlo,is,DM_nemin,DM_nemax,DM_nemaxx,nbands,nbandsx,Qforce_j,time_bl,time_bl_w,time_lo,time_lo_w)
    USE defs,  ONLY: pi, imag
    USE param, ONLY: lmax2, lomax, nume, nloat, iblock, nmat
    USE xa,    ONLY: bkrloc, fj, dfj
    USE xa3,   ONLY: bk3, k3, As
    USE structure, ONLY: rmt, vol, rotij, br1, pos, tauij, rotloc
    USE atspdt,ONLY: p, dp, pe, dpe
    USE lo,    ONLY: nlo, lapw, ilo
    USE dmfts, ONLY: iso, shft
    USE mod_lo, ONLY: lomain
    IMPLICIT NONE
    COMPLEX*16, intent(out) :: alm((lmax2+1)*(lmax2+1),nume),blm((lmax2+1)*(lmax2+1),nume),clm((lomax+1)*(lomax+1),nume,nloat)
    COMPLEX*16, intent(in)  :: Aweight(nbands,nbands)
    INTEGER, intent(in)     :: latom,lfirst,jatom,nemin,nemax,n0,nnlo,is,DM_nemin,DM_nemax,DM_nemaxx,nbands,nbandsx
    LOGICAL, intent(in)     :: Qforce_j
    REAL*8, intent(inout)   :: time_bl, time_bl_w, time_lo, time_lo_w
    ! These data can not be passed as arguments to the function, because they might not be allocated, hence is should be packed inside the module.
    COMPLEX*16, ALLOCATABLE :: aalm(:,:,:), bblm(:,:,:), cclm(:,:,:,:)
    !INTEGER, intent(in)     :: lmax2lmax2_,nbands_dft_,dim_,lomaxlomax_,nloat_
    !COMPLEX*16, intent(out) :: aalm(lmax2lmax2_,nbands_dft_,dim_), bblm(lmax2lmax2_,nbands_dft_,dim_), cclm(lomaxlomax_,nbands_dft_,nloat_,dim_)
    ! locals
    REAL*8  :: fac, twopi, arg1, argt2, arg2, rmt2
    COMPLEX*16 :: phshel, cfac
    REAL*8  :: bkrot2(3)
    REAL*8  :: time1, time1_w, time2, time2_w
    INTEGER :: ii,ibb,i,i3, index, lmx, l, m, maxx, lda, ldb, ldc, jlo, num, lmax2lmax2, lomaxlomax, i_h_k, ik, ilo_max, ind, nbands_dft, ip
    REAL*8 :: h_al(iblock),h_bl(iblock)
    COMPLEX*16:: h_yl
    COMPLEX*16 :: yl((LMAX2+1)*(LMAX2+1)), h_alyl((LMAX2+1)*(LMAX2+1),iblock), h_blyl((LMAX2+1)*(LMAX2+1),iblock)
    COMPLEX*16, allocatable :: xlm_tmp(:,:)
    REAL*8, allocatable     :: h_k(:,:)
    complex*16, allocatable :: h_ablyl_hk(:,:)
    complex*16, allocatable :: ylm_tmp(:,:)
    LOGICAL :: nonzero_shft
    REAL*8  :: rotloc_x_BR1(3,3), rotloc_x_BR1_x_rotij(3,3)
    INTEGER :: j
    complex*16 :: cone, czero
    cone  = cmplx(1.d0, 0.d0, 8)
    czero = cmplx(0.d0, 0.d0, 8)
  
    
    CALL cputim(time1)
    CALL walltim(time1_w)
    
    lmax2lmax2 = (lmax2+1)*(lmax2+1)
    lomaxlomax = (lomax+1)*(lomax+1)
    nbands_dft = nemax-nemin+1
    twopi = pi*2.d0
    fac=4.0d0*pi*rmt(jatom)**2/SQRT(vol)
    lmx=lmax2
    ilo_max=0
    DO l=0,lmax2
       if (ilo(l).gt.ilo_max) ilo_max=ilo(l)
    ENDDO
       
    alm(:,:)=0.d0
    blm(:,:)=0.d0
    clm(:,:,:)=0.d0

    if (Qforce_j) then
       allocate( h_ablyl_hk(lmax2lmax2,iblock) )
       allocate( h_k(iblock,3) )
       
       if (size(aalm,3).NE.3 .or. size(cclm,4).NE.3) print*, 'ERROR alm_blm_clm: Third dimension of aalm should be 3 but is', size(aalm,3)
       if (size(aalm,1).NE.lmax2lmax2 .or. size(cclm,1).NE.lomaxlomax) print*, 'ERROR alm_blm_clm: first dimension of aalm should be', lmax2lmax2, 'but is', size(aalm,1)
       if (size(cclm,3).NE.nloat) print*, 'ERROR alm_blm_clm: Second dimension of cclm should be nloat but is', size(cclm,3)
       if (size(aalm,2).NE.size(cclm,2)) print*, 'ERROR alm_blm_clm: nbands_dft should be dimension in aalm and cclm', size(aalm,3)
       if (size(aalm,2).NE.nbands_dft) print*, 'ERROR alm_blm_clm: size(aalm,2)!=nbands_dft', size(aalm,2), nbands_dft
       
       aalm(1:lmax2lmax2,1:nbands_dft,1:3)=0.d0
       bblm(1:lmax2lmax2,1:nbands_dft,1:3)=0.d0
       cclm(1:lomaxlomax,1:nbands_dft,1:nloat,1:3)=0.d0
    endif
    
    rotloc_x_BR1(:,:) = matmul(rotloc(:,:,jatom), BR1)
    rotloc_x_BR1_x_rotij = matmul( rotloc_x_BR1, rotij(:,:,latom) )
    nonzero_shft = sum(abs(shft(latom,:))) .GT. 1e-10
    
    DO ii=1,n0-nnlo,iblock  !! ii=1,iblock,2*iblock,3*ibloc,... n0-nnlo+d
       ibb=MIN(iblock,n0-nnlo-ii+1)  ! size of this block
       if(ibb.gt.0) then
          if (Qforce_j) h_k(:,:) = 0.d0
          !$OMP PARALLEL DO PRIVATE(bkrot2,bkrloc,yl,arg1,arg2,argt2,phshel)&
          !$OMP& SHARED(h_k,h_alyl)&
          !$OMP& SCHEDULE(STATIC)
          DO i=ii,min(ii+iblock-1,n0-nnlo)
             !---------  rotates ylm(k+K) to ylm(k'+K) where k' is in star of irreducible k. ------------!
             ! bk3(:,i) !-----  reciprocal vector and irreducible vector: G=K+k ----!
             ! BKROT2 = R_n.R_a.(k+K), transformation from the first atom to an equivalent atom 
             bkrot2 = matmul(rotij(:,:,latom), bk3(:,i))
             ! ! !---- BR1 transforms integer reciprocal lattice vectors, as given in the VECTORLIST of LAPW1, into cartesian system ----!
             ! ! ! BKROT3 = R_n.R_a.(k+K), but in cartesian coordinate system
             ! ! bkrot3= matmul(br1,bkrot2)
             ! ! !---- BKRLOC = Rotloc.R_g.(k+K),  is rotates Rotloc.R.(k+K) . Rotation Rotloc entered by user
             ! ! ! Here we use rotloc and not crotloc, hence we rather use the transformation
             ! ! ! from the structure file and not the user input. This is because the wave functions are
             ! ! ! compatible with the rotation in the structure file!
             ! ! bkrloc = matmul(rotloc(:,:,jatom),bkrot3)
             bkrloc = matmul( rotloc_x_BR1,bkrot2)
             !---- YLM = Y_{L}(Rotloc.R_g.(k+K))
             CALL YLM (bkrloc,lmax2,yl,lmax2lmax2)
             ! (R_n.R_a.(k+K)) *  R(first) * 2pi
             arg1 = dot_product(bkrot2, pos(:,lfirst))*twopi
             ! ARGT2 = (R_a.(k+K)).tau_n * 2pi
             argt2= dot_product(bk3(:,i), tauij(:,latom))*twopi
             ! ARG2 = (R_a.(k+K)) *  shft * 2pi
             arg2=0.d0
             if (nonzero_shft) arg2 = dot_product(bk3(:,i), shft(latom,:))*twopi
             ! PHSEHL = e^{I*2pi*( (R_a.(k+K))*tau_n + (K+k)*tau(isym) + (R_n.R_a.(k+K)*R(first)))}
             phshel=EXP( imag*(arg1+arg2+argt2) )
             h_alyl(:(lmx+1)**2,i-ii+1) = dconjg(yl(:(lmx+1)**2))*phshel
             IF(Qforce_j) THEN 
                h_k(i-ii+1,:) = matmul( rotloc_x_BR1_x_rotij, k3(:,i) )
             ENDIF
          ENDDO
          !$OMP END PARALLEL DO          
          
          rmt2=1.0D0/(rmt(jatom)**2)
          DO l=0,lmx
             ! !$OMP PARALLEL DO PRIVATE(i3) SHARED(h_al,h_bl) SCHEDULE(STATIC)
             DO i=ii,min(ii+iblock-1,n0-nnlo)
                i3=i-ii+1
                if(lapw(l)) then
                   h_al(i3)=dfj(l,i,jatom)*pe(l)-fj(l,i,jatom)*dpe(l) 
                   h_bl(i3)=fj(l,i,jatom)*dp(l)-dfj(l,i,jatom)*p(l)
                ELSE
                   h_al(i3) = rmt2*fj(l,i,jatom)/p(l)
                   h_bl(i3) = 0.d0
                ENDIF
             ENDDO
             ! !$OMP END PARALLEL DO
             DO i3=1,ibb ! size of this block                    
                h_blyl(l**2+1:(l+1)**2,i3)=h_bl(i3)*h_alyl(l**2+1:(l+1)**2,i3)
                h_alyl(l**2+1:(l+1)**2,i3)=h_al(i3)*h_alyl(l**2+1:(l+1)**2,i3)
             ENDDO
          ENDDO
          !
          lda=lmax2lmax2
          ldc=lda
          ldb=nmat
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
          CALL zgemm('N','N',lmax2lmax2,nemax-nemin+1,ibb,cone,h_alyl,lda,As(ii,nemin,is),ldb,cone,alm(1,nemin),ldc)
          CALL zgemm('N','N',lmax2lmax2,nemax-nemin+1,ibb,cone,h_blyl,lda,As(ii,nemin,is),ldb,cone,blm(1,nemin),ldc)
          IF (Qforce_j) THEN 
             do i_h_k=1,3
                !$OMP PARALLEL DO SHARED(h_ablyl_hk,h_alyl,h_k) SCHEDULE(STATIC)
                do i3=1,ibb
                   h_ablyl_hk(:,i3) = h_alyl(:,i3)*h_k(i3,i_h_k)   ! h_ablyl <-  alm(lm,K)*K
                enddo
                !$OMP END PARALLEL DO
                CALL zgemm('N','N',lmax2lmax2,nemax-nemin+1,ibb,cone,h_ablyl_hk,lda,As(ii,nemin,is),ldb,cone,aalm(1,1,i_h_k),ldc)
             enddo
             do i_h_k=1,3
                !$OMP PARALLEL DO SHARED(h_ablyl_hk,h_blyl,h_k) SCHEDULE(STATIC)
                do i3=1,ibb
                   h_ablyl_hk(:,i3) = h_blyl(:,i3)*h_k(i3,i_h_k)   ! h_ablyl <-  blm(lm,K)*K
                enddo
                !$OMP END PARALLEL DO
                CALL zgemm('N','N',lmax2lmax2,nemax-nemin+1,ibb,cone,h_ablyl_hk,lda,As(ii,nemin,is),ldb,cone,bblm(1,1,i_h_k),ldc)
             enddo
          ENDIF
       endif
    ENDDO
  
    if (Qforce_j) then
       deallocate( h_k )
       deallocate( h_ablyl_hk )
    endif

    CALL cputim(time2)
    CALL walltim(time2_w)
  
    time_bl=time_bl+time2-time1
    time_bl_w=time_bl_w+time2_w-time1_w
    
    CALL cputim(time1)
    CALL walltim(time1_w)

    !! output: alm, blm, clm, aalm, bblm, cclm
    if (nlo.ne.0) call lomain(rotloc(:,:,jatom),is,nemin,nemax,lfirst,latom,n0,jatom,alm,blm,clm,Qforce_j,aalm,bblm,cclm,lmax2)
    
    !$OMP PARALLEL DO PRIVATE(cfac) SHARED(alm,blm) SCHEDULE(STATIC)
    DO l=0,lmax2
       cfac=fac*imag**l
       alm(l**2+1:(l+1)**2,nemin:nemax)=alm(l**2+1:(l+1)**2,nemin:nemax)*cfac
       blm(l**2+1:(l+1)**2,nemin:nemax)=blm(l**2+1:(l+1)**2,nemin:nemax)*cfac
    ENDDO
    !$OMP END PARALLEL DO
    !$OMP PARALLEL DO PRIVATE(cfac) SHARED(clm) SCHEDULE(STATIC)
    DO l=0,lomax
       cfac=fac*imag**l
       clm(l**2+1:(l+1)**2,nemin:nemax,1:ilo(l))=clm(l**2+1:(l+1)**2,nemin:nemax,1:ilo(l))*cfac
    ENDDO
    !$OMP END PARALLEL DO
    
    IF (Qforce_j) THEN
       !$OMP PARALLEL DO PRIVATE(cfac) SHARED(aalm,bblm) SCHEDULE(STATIC)
       DO l=0,lmax2
          cfac = fac*imag**l
          aalm(l**2+1:(l+1)**2,1:nbands_dft,1:3) = aalm(l**2+1:(l+1)**2,1:nbands_dft,1:3)*cfac
          bblm(l**2+1:(l+1)**2,1:nbands_dft,1:3) = bblm(l**2+1:(l+1)**2,1:nbands_dft,1:3)*cfac
       ENDDO
       !$OMP END PARALLEL DO
       !$OMP PARALLEL DO PRIVATE(cfac) SHARED(cclm) SCHEDULE(STATIC)
       DO l=0,lomax
          cfac = fac*imag**l
          cclm(l**2+1:(l+1)**2,1:nbands_dft,1:ilo(l),1:3) = cclm(l**2+1:(l+1)**2,1:nbands_dft,1:ilo(l),1:3)*cfac
       ENDDO
       !$OMP END PARALLEL DO
    END IF
    
    ! Chenging alm,blm,clm to make them DMFT-like
    !   alm_{lm,i}     <--  alm_{lm,j} * A^{DMFT}_{ji}
    !   blm_{lm,i}     <--  blm_{lm,j} * A^{DMFT}_{ji}
    !   clm_{lm,i,jlo} <--  clm_{lm,j,jlo} * A^{DMFT}_{ji}
    !
    ALLOCATE( xlm_tmp(lmax2lmax2,nbandsx) )
    CALL zgemm('N','N',lmax2lmax2,nbandsx,nbands,cone,alm(1,DM_nemin),lmax2lmax2,Aweight,nbands,czero,xlm_tmp,lmax2lmax2)
    alm(:,DM_nemin:DM_nemaxx)=xlm_tmp(:,:)
    CALL zgemm('N','N',lmax2lmax2,nbandsx,nbands,cone,blm(1,DM_nemin),lmax2lmax2,Aweight,nbands,czero,xlm_tmp,lmax2lmax2)
    blm(:,DM_nemin:DM_nemaxx)=xlm_tmp(:,:)
    IF (Qforce_j) THEN
       ! Chenging aalm,bblm,cclm to make them DMFT-like
       do ik=1,3
          CALL zgemm('N','N',lmax2lmax2,nbandsx,nbands,cone,aalm(1,1+DM_nemin-nemin,ik),lmax2lmax2,Aweight,nbands,czero,xlm_tmp,lmax2lmax2)
          aalm(:,DM_nemin-nemin+1:DM_nemaxx-nemin+1,ik)=xlm_tmp(:,:)
          CALL zgemm('N','N',lmax2lmax2,nbandsx,nbands,cone,bblm(1,1+DM_nemin-nemin,ik),lmax2lmax2,Aweight,nbands,czero,xlm_tmp,lmax2lmax2)
          bblm(:,DM_nemin-nemin+1:DM_nemaxx-nemin+1,ik)=xlm_tmp(:,:)
       enddo
    ENDIF
    DEALLOCATE( xlm_tmp )
    ALLOCATE( xlm_tmp(lomaxlomax,nbandsx) )
    DO jlo=1,ilo_max
       CALL zgemm('N','N',lomaxlomax,nbandsx,nbands,cone,clm(1,DM_nemin,jlo),lomaxlomax,Aweight,nbands,czero,xlm_tmp,lomaxlomax)
       clm(:,DM_nemin:DM_nemaxx,jlo)=xlm_tmp(:,:)
    ENDDO
    IF (Qforce_j) THEN
       do ik=1,3
          DO jlo=1,ilo_max
             CALL zgemm('N','N',lomaxlomax,nbandsx,nbands,cone,cclm(1,1+DM_nemin-nemin,jlo,ik),lomaxlomax,Aweight,nbands,czero,xlm_tmp,lomaxlomax)
             cclm(:,DM_nemin-nemin+1:DM_nemaxx-nemin+1,jlo,ik)=xlm_tmp(:,:)
          ENDDO
       enddo
    ENDIF
    DEALLOCATE( xlm_tmp )
    
    CALL cputim(time2)
    CALL walltim(time2_w)
  
    time_lo=time_lo+time2-time1
    time_lo_w=time_lo_w+time2_w-time1_w

  END SUBROUTINE cmp_alm_blm_clm
END module alm_blm_clm
