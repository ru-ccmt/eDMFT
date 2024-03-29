module w_atpar
  IMPLICIT NONE
  REAL*8 , ALLOCATABLE :: w_alo(:,:,:)
  REAL*8 , ALLOCATABLE :: w_blo(:,:,:)
  REAL*8 , ALLOCATABLE :: w_clo(:,:,:)

  REAL*8,  ALLOCATABLE :: ul_Rmt(:,:,:)
  REAL*8,  ALLOCATABLE :: dul_Rmt(:,:,:)
  REAL*8,  ALLOCATABLE  :: ri_mat(:,:,:,:)
  
  INTEGER, ALLOCATABLE :: w_lm(:,:,:)
  INTEGER, ALLOCATABLE :: w_lmmax(:)
  !REAL*8 , ALLOCATABLE :: w_R(:,:)

  REAL*8 , ALLOCATABLE :: w_xwt1(:,:)
  REAL*8 , ALLOCATABLE :: w_xwt1l(:,:)
  REAL*8 , ALLOCATABLE :: w_xwt1h(:,:)
  REAL*8 , ALLOCATABLE :: w_xwteh(:,:)
  REAL*8 , ALLOCATABLE :: w_xwtel(:,:)

  REAL*8 , ALLOCATABLE :: w_a1lo(:,:,:,:)
  REAL*8 , ALLOCATABLE :: w_b1lo(:,:,:,:)
  REAL*8 , ALLOCATABLE :: w_RRAD1(:,:,:)
  REAL*8 , ALLOCATABLE :: w_RRAD2(:,:,:)
  REAL*8 , ALLOCATABLE :: w_RADE1(:,:,:)
  REAL*8 , ALLOCATABLE :: w_RADE2(:,:,:)

  REAL*8 , ALLOCATABLE :: w_rholm(:,:,:)
  
  COMPLEX*16, allocatable :: Vts(:,:,:,:)
  INTEGER,    allocatable :: lpv(:,:), lv(:,:), mpv(:,:), mv(:,:), Nhmx(:)
  
CONTAINS

  SUBROUTINE w_allocate0(nat)
    USE param
    IMPLICIT NONE
    INTEGER, intent(in) :: nat
    ! 
    ALLOCATE( w_lm(1:2,1:NCOM,1:nat) )
    ALLOCATE( w_lmmax(1:nat) )
  END SUBROUTINE w_allocate0

  SUBROUTINE w_allocate(LM_MAX, nat)
    USE param
    IMPLICIT NONE
    INTEGER, intent(in) :: LM_MAX, nat

    ALLOCATE( w_alo(0:lomax,1:nloat,1:nat) )
    ALLOCATE( w_blo(0:lomax,1:nloat,1:nat) )
    ALLOCATE( w_clo(0:lomax,1:nloat,1:nat) )

    ALLOCATE( ul_Rmt(nloat+2,0:lmax2,1:nat) )
    ALLOCATE( dul_Rmt(nloat+2,0:lmax2,1:nat) )
    ALLOCATE( ri_mat(2+nloat,2+nloat,0:lmax2,1:nat) )
    
    ALLOCATE( w_xwt1(0:21,1:nat) )
    ALLOCATE( w_xwt1l(0:3,1:nat) )
    ALLOCATE( w_xwt1h(0:3,1:nat) )
    ALLOCATE( w_xwteh(0:3,1:nat) )
    ALLOCATE( w_xwtel(0:3,1:nat) )

    ALLOCATE( w_a1lo(1:nrad,1:nloat,0:lomax,1:nat) )
    ALLOCATE( w_b1lo(1:nrad,1:nloat,0:lomax,1:nat) )
    ALLOCATE( w_RRAD1(1:nrad,0:lmax2,1:nat) )
    ALLOCATE( w_RRAD2(1:nrad,0:lmax2,1:nat) )
    ALLOCATE( w_RADE1(1:nrad,0:lmax2,1:nat) )
    ALLOCATE( w_RADE2(1:nrad,0:lmax2,1:nat) )

    ALLOCATE( w_rholm(1:NRAD,1:LM_MAX,1:nat) )

  END SUBROUTINE w_allocate

  SUBROUTINE w_deallocate
    IMPLICIT NONE
    DEALLOCATE( w_alo )
    DEALLOCATE( w_blo )
    DEALLOCATE( w_clo )

    DEALLOCATE( ul_Rmt )
    DEALLOCATE( dul_Rmt )
    DEALLOCATE( ri_mat )

    DEALLOCATE( w_lm )
    DEALLOCATE( w_lmmax )
    !DEALLOCATE( w_R )

    DEALLOCATE( w_xwt1 )
    DEALLOCATE( w_xwt1l )
    DEALLOCATE( w_xwt1h )
    DEALLOCATE( w_xwteh )
    DEALLOCATE( w_xwtel )

    DEALLOCATE( w_a1lo )
    DEALLOCATE( w_b1lo )
    DEALLOCATE( w_RRAD1 )
    DEALLOCATE( w_RRAD2 )
    DEALLOCATE( w_RADE1 )
    DEALLOCATE( w_RADE2 )

    DEALLOCATE( w_rholm )
  END SUBROUTINE w_deallocate

  SUBROUTINE w_deallocate1
    IMPLICIT NONE
    DEALLOCATE( w_alo )
    DEALLOCATE( w_blo )
    DEALLOCATE( w_clo )
    DEALLOCATE( ul_Rmt )
    DEALLOCATE( dul_Rmt )
    DEALLOCATE( ri_mat )
    !DEALLOCATE( w_lm )
    !DEALLOCATE( w_lmmax )
    DEALLOCATE( w_xwt1 )
    DEALLOCATE( w_xwt1l )
    DEALLOCATE( w_xwt1h )
    DEALLOCATE( w_xwteh )
    DEALLOCATE( w_xwtel )
    DEALLOCATE( w_a1lo )
    DEALLOCATE( w_b1lo )
    DEALLOCATE( w_RRAD1 )
    DEALLOCATE( w_RRAD2 )
    DEALLOCATE( w_RADE1 )
    DEALLOCATE( w_RADE2 )
  END SUBROUTINE w_deallocate1
  SUBROUTINE w_deallocate2
    IMPLICIT NONE
    DEALLOCATE( w_lm )
    DEALLOCATE( w_lmmax )
    DEALLOCATE( w_rholm )
  END SUBROUTINE w_deallocate2
  
  SUBROUTINE precmp_w_atpar(cform,zz,nnlo,ISCF)
    USE com,   ONLY: nat, rel
    USE structure, ONLY: R0,  mult, JRI, DX, Rmt
    USE param, ONLY: nrad, lmax2, lomax, NLOAT, filename_V_sph
    USE lo,    ONLY: elo_store, elo, loor, lapw, nlo, nlov, nlon, ilo, rlo
    USE atspdt,ONLY: e_store, el 
    USE readPotential, ONLY: init_V_vsp, read_V_vsp, close_V_vsp
    USE com_mpi
    IMPLICIT NONE
    CHARACTER*4, intent(in)  :: cform
    REAL*8, intent(in)       :: ZZ(*)
    INTEGER, intent(out)     :: nnlo, ISCF
    !locals
    INTEGER :: jatom, lfirst, i, imax, l, jlo, jlop
    INTEGER :: nr
    REAL*8  :: Vr(NRAD)
    LOGICAL :: qprint_
    
    qprint_ = myrank.EQ.master
    
    CALL init_V_vsp(filename_V_sph, 18, ISCF)
!!! BRISI
    !open(1999,FILE='_atpar.dat',status='unknown')
!!! BRISI
    
    do jatom=1,nat
       ! Here we calculate all radial wave functions and other quantities
       ! computed by atpar, which do not depend on k vector.

       ! This computes on the fly the following quantities: el,elo,loor,rlo,lapw,nlov,nlo,nlon,ilo
       CALL get_ilo(e_store,elo_store,jatom,mult,nat,lmax2,lomax,nloat, el,elo,loor,rlo,lapw,nlov,nlo,nlon,ilo)
       ! Here we read the spherically symetric part of Kohn-Sham potential
       Vr(:)=0.d0
       nr = jri(jatom)
       CALL read_V_vsp(Vr(:nr),jatom,nr)
       
       if (qprint_) then
          lfirst=1
          do i=1,jatom-1
             lfirst=lfirst + mult(i)
          enddo
          CALL Print_pre_atpar(jatom,lfirst,el,lmax2)
       endif
       ! calculates radial functions U(r), UE(R), ..., overlaps, coefficients,....
       CALL ATPARN(w_RRAD1(:,:,jatom),w_RADE1(:,:,jatom),w_RRAD2(:,:,jatom),w_RADE2(:,:,jatom),w_a1lo(:,:,:,jatom),w_b1lo(:,:,:,jatom),ul_Rmt(:,:,jatom),dul_Rmt(:,:,jatom),ri_mat(:,:,:,jatom),w_alo(:,:,jatom),w_blo(:,:,jatom),w_clo(:,:,jatom), Vr,nr,el,elo,loor,rlo,lapw,nlov,nlo,nlon,ilo,rel,ZZ(jatom),qprint_,jri(jatom),r0(jatom),dx(jatom),Rmt(jatom),lmax2,lomax,nloat,NRAD)
       
       nnlo=nlo+nlon+nlov

       if (qprint_) then
          DO l=0,lmax2
             WRITE(6,'(A,I2,2x,A,I1,2x,A,F12.8,2x,A,L1,2x,A,F10.5)') 'jatom=', jatom, 'l=', l, 'e(l)=', el(l), 'lapw=', lapw(l), '<dot_u|dot_u>=', ri_mat(2,2,l,jatom)
          ENDDO
          DO l=0,lomax
             if (ilo(l).GT.0) then
                WRITE(6,'(A,I2,2x,A,I1,2x,A,I1,2x,A,3F15.8)') 'jatom=', jatom, 'l=', l, 'ilo=', ilo(l), 'elo(l)=', elo(l,:)
                WRITE(6,'(A)',advance='no') '(rlo,loor): '
                do jlo=1,ilo(l)
                   WRITE(6,'(A,L1,A,L1,A)',advance='no') '(',rlo(jlo,l),',',loor(jlo,l),')  '
                enddo
                WRITE(6,*)
                WRITE(6,'(A12,A12,A12)',advance='no') '<u|u_lo>', '<udot|u_lo>', '<u_lo|u_lo>'
                do jlo=1,ilo(l)
                   WRITE(6,'(F10.5,2x,F10.5,2x,F10.5,2x)',advance='no') ri_mat(2+jlo,1,l,jatom), ri_mat(2+jlo,2,l,jatom), ri_mat(2+jlo,2+jlo,l,jatom)
                enddo
                WRITE(6,*)
             endif
          ENDDO
       endif
       
!!! BRISI
       !WRITE(1999,*) jatom, lmax2, lomax
       !do l=0,lmax2
       !   WRITE(1999,'(I4,F20.16,1x,F20.16)') l, el(l), ri_mat(2,2,l,jatom)
       !enddo
       !do l=0,lomax
       !   DO jlo=1,ilo(l)
       !      WRITE(1999,'(I3,1x,I3,1x,F20.16,1x,F20.16,1x,F20.16,1x)') l, jlo, elo(l,jlo), ri_mat(1,2+jlo,l,jatom), ri_mat(2,2+jlo,l,jatom)
       !   ENDDO
       !enddo
       !do l=0,lomax
       !   DO jlo=1,ilo(l)
       !      DO jlop=1,ilo(l)
       !         WRITE(1999,'(3I3,1x,F20.16)') l, jlo, jlop, ri_mat(2+jlo,2+jlop,l,jatom)
       !      ENDDO
       !   ENDDO
       !enddo
       !WRITE(1999,*)
!!! BRISI
       
       !ul_Rmt(:,:,:)=0.d0
       !dul_Rmt(:,:,:)=0.d0
       !do l=0,lmax2
       !   ul_Rmt(1,l,jatom) = P(l) 
       !   ul_Rmt(2,l,jatom) = PE(l)
       !   dul_Rmt(1,l,jatom) = DP(l)
       !   dul_Rmt(2,l,jatom) = DPE(l)
       !enddo
       !do l=0,lomax
       !   DO jlo=1,ilo(l)
       !      ul_Rmt(2+jlo,l,jatom) = PLO(jlo,l)
       !      dul_Rmt(2+jlo,l,jatom) = DPLO(jlo,l)
       !   ENDDO
       !enddo
       
       !ri_mat(:,:,:,:)=0.0
       !do l=0,lmax2
       !   ri_mat(1,1,l,jatom)=1.0    ! <u|u>
       !   ri_mat(1,2,l,jatom)=0.0    ! <udot|u>
       !   ri_mat(2,1,l,jatom)=0.0    ! <u|udot>
       !   ri_mat(2,2,l,jatom)=pei(l) ! <udot|udot>
       !enddo
       !do l=0,lomax
       !   DO jlo=1,ilo(l)
       !      ri_mat(1,2+jlo,l,jatom) = pi12lo(jlo,l)  ! <u | u_lo>
       !      ri_mat(2+jlo,1,l,jatom) = pi12lo(jlo,l)  ! <u_lo | u>
       !      ri_mat(2,2+jlo,l,jatom) = pe12lo(jlo,l)  ! <udot | u_lo>
       !      ri_mat(2+jlo,2,l,jatom) = pe12lo(jlo,l)  ! <u_lo | udot>
       !   ENDDO
       !   DO jlo=1,ilo(l)
       !      DO jlop=1,ilo(l)
       !         ri_mat(2+jlo,2+jlop,l,jatom) = pr12lo(jlo,jlop,l)  ! <u_lo | u_lo >
       !      ENDDO
       !   ENDDO
       !enddo

!!!BRISI DEBUGGING
       !print *, 'Debugging'
       !DO l=0,lmax2
       !   print *, 'jatom=', jatom, 'l=', l, 'ilo=', ilo(l)
       !   DO jlo=1,ilo(l)
       !      DO jlop=1,ilo(l)
       !         print *, jlo, jlop, ri_mat(2+jlo,2+jlop,l,jatom)
       !      ENDDO
       !   ENDDO
       !   print *, '<u|u_lo>', '<dot u|u_lo>'
       !   DO jlo=1,ilo(l)
       !      print *, jlo, ri_mat(2+jlo,1,l,jatom), ri_mat(2+jlo,2,l,jatom)
       !   ENDDO
       !ENDDO
    enddo
    CALL close_V_vsp()
!!! BRISI
    !close(1999)
!!! BRISI
  END SUBROUTINE precmp_w_atpar

  SUBROUTINE retrieve_w_atpar(jatom,lfirst,lmmax)
    USE xa,    ONLY: LM, R
    USE param, ONLY: NCOM, lmax2, nloat, lomax, nrad
    USE lo,    ONLY: elo_store, elo, rlo, loor, lapw, nlo, nlov, nlon, ilo, a1lo, b1lo, alo, blo, clo, plo, dplo, pi12lo, pe12lo, pr12lo
    USE atspdt,ONLY: e_store, el, P, DP, PE, DPE, PEI
    USE structure, ONLY: mult, r0, jri, dx
    USE com,   ONLY: nat
    IMPLICIT NONE
    INTEGER, intent(in) :: jatom
    INTEGER, intent(out) :: lfirst,lmmax
    COMMON /RADFU/   RRAD1(NRAD,0:LMAX2),RADE1(NRAD,0:LMAX2),RRAD2(NRAD,0:LMAX2),RADE2(NRAD,0:LMAX2)
    REAL*8 :: RRAD1, RRAD2, RADE1, RADE2
    ! locals
    INTEGER :: l, jlo, jlop, i
    ! Quantities computed by "atpar" which do not depend on k-point

    ! This computes on the fly the following quantities: el,elo,loor,rlo,lapw,nlov,nlo,nlon,ilo
    CALL get_ilo(e_store,elo_store,jatom,mult,nat,lmax2,lomax,nloat, el,elo,loor,rlo,lapw,nlov,nlo,nlon,ilo)
    
    lm(1:2,1:NCOM) = w_lm(1:2,1:NCOM,jatom)
    lmmax = w_lmmax(jatom)

    alo(0:lomax,1:nloat) = w_alo(0:lomax,1:nloat,jatom)
    blo(0:lomax,1:nloat) = w_blo(0:lomax,1:nloat,jatom)
    clo(0:lomax,1:nloat) = w_clo(0:lomax,1:nloat,jatom)
    a1lo(1:nrad,1:nloat,0:lomax) = w_a1lo(1:nrad,1:nloat,0:lomax,jatom)
    b1lo(1:nrad,1:nloat,0:lomax) = w_b1lo(1:nrad,1:nloat,0:lomax,jatom)
    RRAD1(1:nrad,0:lmax2) = w_RRAD1(1:nrad,0:lmax2,jatom)
    RRAD2(1:nrad,0:lmax2) = w_RRAD2(1:nrad,0:lmax2,jatom)
    RADE1(1:nrad,0:lmax2) = w_RADE1(1:nrad,0:lmax2,jatom)
    RADE2(1:nrad,0:lmax2) = w_RADE2(1:nrad,0:lmax2,jatom)
    
    do l=0,lmax2
       P(l) = ul_Rmt(1,l,jatom)
       PE(l) = ul_Rmt(2,l,jatom)
       DP(l) = dul_Rmt(1,l,jatom)
       DPE(l) = dul_Rmt(2,l,jatom)
    enddo
    do l=0,lomax
       DO jlo=1,ilo(l)
          PLO(jlo,l) = ul_Rmt(2+jlo,l,jatom)
          DPLO(jlo,l) = dul_Rmt(2+jlo,l,jatom)
       ENDDO
    enddo

    pei(0:lmax2) = ri_mat(2,2,0:lmax2,jatom)       ! <udot|udot>
    do l=0,lomax
       DO jlo=1,ilo(l)
          pi12lo(jlo,l) = ri_mat(1,2+jlo,l,jatom)  ! <u | u_lo>
          pe12lo(jlo,l) = ri_mat(2,2+jlo,l,jatom)  ! <udot | u_lo>
       ENDDO
       DO jlo=1,ilo(l)
          DO jlop=1,ilo(l)
             pr12lo(jlo,jlop,l) = ri_mat(2+jlo,2+jlop,l,jatom) ! <u_lo | u_lo >
          ENDDO
       ENDDO
    enddo
    
    ! which atom is first of that type?
    lfirst=1
    do i=1,jatom-1
       lfirst=lfirst + mult(i)
    enddo
    ! Radial mesh
    R(:)=0.d0
    DO i=1,jri(jatom)
       R(i)=r0(jatom)*Exp((i-1.)*dx(jatom)) ! Radial mesh
    ENDDO
  END SUBROUTINE retrieve_w_atpar
  
  SUBROUTINE Read_nsh()
    USE param, ONLY: filename_V_nsh, lomax, ngau
    USE com, ONLY: nat
    USE readPotential, ONLY: init_V_nsh, read_V_nsh, close_V_nsh
    IMPLICIT NONE
    ! locals
    INTEGER :: jatom, ihmx, ihmx_max
    COMPLEX*16, allocatable :: Vts_(:,:,:)
    INTEGER,    allocatable :: lpv_(:), lv_(:), mpv_(:), mv_(:)
    ! We first read the entire potential to determine maximal number of nonzero components (ngau)
    allocate( Vts_(3,3,ngau) )      ! non-spherical potential matrix elements
    allocate( lpv_(ngau), lv_(ngau), mpv_(ngau), mv_(ngau) )
    ihmx_max=0
    CALL init_V_nsh(filename_V_nsh, 9902) 
    !WRITE(6,*) 'First read_V_nsh'
    do jatom=1,nat
       !WRITE(6,*) 'jatom=', jatom
       CALL read_V_nsh(Vts_,lv_,lpv_,mv_,mpv_,ihmx, jatom,lomax,ngau)
       !WRITE(6,*) 'ihmx=', ihmx
       if (ihmx.GT.ihmx_max) ihmx_max=ihmx
    enddo
    CALL close_V_nsh()
    ngau = ihmx_max+1  ! CAREFUL: NGAU was constant in original code. Here we change it to the smallest acceptable number.
    deallocate( Vts_ )
    deallocate( lpv_, lv_, mpv_, mv_ )

    ! Here we reread non-spherical potential to store it in Vts(:,:,:,:)
    allocate( Vts(3,3,ngau,nat) )      ! non-spherical potential matrix elements
    allocate( lpv(ngau,nat), lv(ngau,nat), mpv(ngau,nat), mv(ngau,nat) )
    allocate( Nhmx(nat) )
    CALL init_V_nsh(filename_V_nsh, 9902) 
    !WRITE(6,*) 'Second read_V_nsh'
    do jatom=1,nat
       !WRITE(6,*) 'jatom=', jatom
       CALL read_V_nsh(Vts(:,:,:,jatom),lv(:,jatom),lpv(:,jatom),mv(:,jatom),mpv(:,jatom),Nhmx(jatom), jatom,lomax,ngau)
       !WRITE(6,*) 'ihmx=', ihmx
    enddo
    CALL close_V_nsh()
  END SUBROUTINE Read_nsh
  
  SUBROUTINE Deallocate_nsh()
    IMPLICIT NONE
    deallocate( Vts )
    deallocate( lpv, lv, mpv, mv )
    deallocate( Nhmx )
  END SUBROUTINE Deallocate_nsh
  
end module w_atpar
