! @Copyright 2007 Kristjan Haule
! 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This module builds transformation matrix between cubis<->spheric harmonics !!
!! and between (l,m,s) base and (l,j,jm) base                                 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE cmp_t2c(T2C, is, ilsb, ltpan, lmtind, ndim, norb, natom, nsort, lmaxu, ilsbmax2)
  IMPLICIT NONE
  !---------- Passed variables ----------
  COMPLEX*16, intent(out) :: T2C(ndim,ndim)
  INTEGER,    intent(in)  :: is(natom)
  INTEGER,    intent(in)  :: ilsb(nsort)
  INTEGER,    intent(in)  :: ltpan(0:lmaxu,natom)
  INTEGER,    intent(in)  :: lmtind(ilsbmax2,natom)
  INTEGER,    intent(in)  :: ndim, norb, natom, nsort, lmaxu, ilsbmax2
!f2py integer intent(hide), depend(is)             :: natom=shape(is,0)
!f2py integer intent(hide), depend(ilsb)           :: nsort=shape(ilsb,0)
!f2py integer intent(hide), depend(lmtind)         :: ilsbmax2=shape(lmtind,0)
!f2py integer intent(hide), depend(ltpan)          :: lmaxu=shape(ltpan,0)-1
  ! Temporary variables
  COMPLEX*16 :: TC(-3:3,7,0:3)             ! sph->cubic convertion matrix
  INTEGER    :: iatom, isort, l, lm, inds, inde, i, j
!!! T2C is matrix to covert the self-energy or Hamiltonian 
!!! matrixes from spherical to cubic: 
!!! SE_c=(T2C)^{+}*SE_s*T2C
!!! SE_s=T2C*SE_c*(T2C)^{+}

!!! Builts transformation matrix to convert cubic harmonics to spherical harmonics
  CALL TRANSC(TC)
  T2C = 0.0 
  DO iatom = 1,natom
     isort = is(iatom)
     DO l = 0,ilsb(isort)
        !if (.NOT.ltpan(l,iatom)) cycle
        if (ltpan(l,iatom).EQ.0) cycle
        lm=l**2+1 ! the first index in the [-l,l] range
        inds = lmtind(lm,iatom)  ! was m
        inde = inds + 2*l        ! was m1
        IF (l .EQ. 0) THEN
           T2C(inds:inde,inds:inde) = (1.0,0.0)
        ELSE IF (l.EQ.1) THEN
           T2C(inds:inde,inds:inde) = TC(-1:1,1:3,1)
        ELSE IF (l.EQ.2) THEN
           T2C(inds:inde,inds:inde) = TC(-2:2,1:5,2)
        ELSE IF (l.EQ.3) THEN
           T2C(inds:inde,inds:inde) = TC(-3:3,1:7,3)
        ELSE 
           WRITE(0,*) "spheric to cubic implemented only for l<=3 !"
        ENDIF
     ENDDO
  ENDDO
  IF (ndim.EQ.2*norb) T2C(norb+1:,norb+1:) = T2C(1:norb,1:norb)
  
END SUBROUTINE cmp_t2c

SUBROUTINE cmp_t2j(T2C, bndindJ, bndind, ndim)
!!! Transformation matrix to state with good quantum number J
!  USE clbsch, ONLY : ClebschG
  IMPLICIT NONE
  REAL*8, intent(in)      :: bndindJ(ndim,4)
  REAL*8, intent(in)      :: bndind(ndim,4)
  COMPLEX*16, intent(out) :: T2C(ndim,ndim)
  INTEGER, intent(in)     :: ndim
!f2py integer intent(hide), depend(bndind)  :: ndim=shape(bndind,0)
  ! External function
  REAL*8 :: ClebschG
! Temporary
  REAL*8     :: j0, mj, s0, l, ml, ms
  INTEGER    :: iat, i, j
  s0=0.5
  T2C=0
  do i=1,ndim
     do j=1,ndim
        j0  = bndindJ(i,3)
        mj  = bndindJ(i,4)
        l   = bndindJ(i,2)
        iat = bndindJ(i,1)
        if (l.NE.bndind(j,2) .OR. iat.NE.bndind(j,1)) CYCLE
        ml  = bndind(j,3)
        ms = 1.5-bndind(j,4)
        T2C(i,j) = ClebschG(j0,mj,l,ml,s0,ms)
     enddo
  enddo
  T2C = transpose(T2C)
END SUBROUTINE cmp_t2j


!  SUBROUTINE cmp_rot(T2C, is, ilsb, ltpan, lmtind, rotateq, rot_vec_angle, & 
!       rot_num, ndim, norb, add_rot, rot_cor, correlated, bcorrelated, &
!       max_nrot, natom, nsort, lmaxu, ilsbmax2, ncorr_all)
!    IMPLICIT NONE
!    !---------- Passed variables ----------
!    COMPLEX*16, intent(inout) :: T2C(ndim,ndim)
!    INTEGER,    intent(in)  :: is(natom)
!    INTEGER,    intent(in)  :: ilsb(nsort)
!    INTEGER,    intent(in)  :: ltpan(0:lmaxu,natom)
!    INTEGER,    intent(in)  :: lmtind(ilsbmax2,natom)
!    LOGICAL,    intent(in)  :: rotateq(natom)
!    REAL*8,     intent(in)  :: rot_vec_angle(4,max_nrot,natom)
!    INTEGER,    intent(in)  :: rot_num(natom)
!    INTEGER,    intent(in)  :: ndim, norb
!    LOGICAL,    intent(in)  :: add_rot
!    INTEGER,    intent(in)  :: ncorr_all
!    REAL*8,     intent(in)  :: rot_cor(ncorr_all,ncorr_all)
!    INTEGER,    intent(in)  :: max_nrot, natom, nsort, lmaxu, ilsbmax2
!    INTEGER,    intent(in)  :: correlated(ndim)
!    LOGICAL,    intent(in)  :: bcorrelated(ndim)
!  !f2py integer intent(hide), depend(is)             :: natom=shape(is,0)
!  !f2py integer intent(hide), depend(ilsb)           :: nsort=shape(ilsb,0)
!  !f2py integer intent(hide), depend(lmtind)         :: ilsbmax2=shape(lmtind,0)
!  !f2py integer intent(hide), depend(ltpan)          :: lmaxu=shape(ltpan,0)-1
!  !f2py integer intent(hide), depend(rot_vec_angle)  :: max_nrot=shape(rot_vec_angle,1)
!  !f2py integer intent(hide), depend(rot_cor)       :: ncorr_all = shape(rot_cor,0)
!    ! Temporary variables
!    REAL*8     :: R0(3,3,natom), Runity(3,3), Rd(5,5)
!    COMPLEX*16 :: Rc3(3,3), Rc5(5,5), Uc_rotate(ndim,ndim)
!    INTEGER    :: iatom, r, isort, l, lm, inds, inde, i, j
!  !!! Here we create Rotation matrix R0(3,3) which rotates each atom in space as specified in the input file
!    IF (ndim.GT.norb) THEN
!       WRITE(0,*) 'Seems to me you are using spin-orbit coupling and want to rotate atoms. This is not yet implemented. Boiling out!'
!       CALL EXIT(1)
!    ENDIF
!    Runity=0;  Runity(1,1)=1.;  Runity(2,2)=1.;  Runity(3,3)=1.
!    DO iatom=1,natom
!       R0(:,:,iatom) = Runity
!       DO r=1,rot_num(iatom)
!          CALL Rotatea(rot_vec_angle(4,r,iatom), rot_vec_angle(1:3,r,iatom), R0(:,:,iatom))
!       ENDDO
!    ENDDO
!  !!! We add rotation to each atom as specified in the input file
!  !!! to have local coordinate system in each correlated atom simple
!    DO iatom = 1,natom
!       IF (.NOT.rotateQ(iatom)) CYCLE
!       isort = is(iatom)
!       DO l = 0,ilsb(isort)
!          IF (.NOT.ltpan(l,iatom)) CYCLE
!          lm=l**2+1               ! the first index in the [-l,l] range
!          inds = lmtind(lm,iatom) ! was m
!          inde = inds + 2*l       ! was m1
!          IF (l .EQ. 1) THEN
!             Rc3 = R0(:,:,iatom)
!             T2C(inds:inde,inds:inde) = matmul(T2C(inds:inde,inds:inde),transpose(conjg(Rc3)))
!          ELSE IF (l .EQ. 2) THEN
!             CALL Rotate_d(Rd, R0(:,:,iatom))
!             Rc5 = Rd
!             T2C(inds:inde,inds:inde) = matmul(T2C(inds:inde,inds:inde),transpose(conjg(Rc5)))
!          ELSE  IF (l.Gt.2) THEN
!             ! Rotation for f-s not yet implemented
!             WRITE(0,*) 'ERROR: Rotation for fs not yet implemented! '
!             CALL EXIT(1)
!          ENDIF
!       ENDDO
!    ENDDO
!  !!! At the end we rotate correlated block according to input cix file.
!    IF (add_rot) THEN
!       Uc_rotate = 0
!       DO i=1,ndim
!          IF (.NOT.bcorrelated(i)) Uc_rotate(i,i)=1
!       ENDDO
!       DO i=1,ncorr_all
!          DO j=1,ncorr_all
!             Uc_rotate(correlated(i),correlated(j)) = rot_cor(i,j)
!          ENDDO
!       ENDDO
!       T2C = matmul(T2C, Uc_rotate)
!    ENDIF
!  END SUBROUTINE cmp_rot


!SUBROUTINE cmp_t2c(T2C, is, ilsb, ltpan, lmtind, rotateq, rot_vec_angle, & 
!     rot_num, ndim, norb, add_rot, rot_cor, correlated, bcorrelated, &
!     max_nrot, natom, nsort, lmaxp1, ilsbmax2, ncorr_all)
!  IMPLICIT NONE
!  !---------- Passed variables ----------
!  COMPLEX*16, intent(out) :: T2C(ndim,ndim)
!  INTEGER,    intent(in)  :: is(natom)
!  INTEGER,    intent(in)  :: ilsb(nsort)
!  INTEGER,    intent(in)  :: ltpan(0:lmaxu,natom)
!  INTEGER,    intent(in)  :: lmtind(ilsbmax2,natom)
!  LOGICAL,    intent(in)  :: rotateq(natom)
!  REAL*8,     intent(in)  :: rot_vec_angle(4,max_nrot,natom)
!  INTEGER,    intent(in)  :: rot_num(natom)
!  INTEGER,    intent(in)  :: ndim, norb
!  LOGICAL,    intent(in)  :: add_rot
!  INTEGER,    intent(in)  :: ncorr_all
!  REAL*8,     intent(in)  :: rot_cor(ncorr_all,ncorr_all)
!  INTEGER,    intent(in)  :: max_nrot, natom, nsort, lmaxp1, ilsbmax2
!  INTEGER,    intent(in)  :: correlated(ndim)
!  LOGICAL,    intent(in)  :: bcorrelated(ndim)
!!f2py integer intent(hide), depend(is)             :: natom=shape(is,0)
!!f2py integer intent(hide), depend(ilsb)           :: nsort=shape(ilsb,0)
!!f2py integer intent(hide), depend(lmtind)         :: ilsbmax2=shape(lmtind,0)
!!f2py integer intent(hide), depend(ltpan)          :: lmaxu=shape(ltpan,0)-1
!!f2py integer intent(hide), depend(rot_vec_angle)  :: max_nrot=shape(rot_vec_angle,1)
!!f2py integer intent(hide), depend(rot_cor)       :: ncorr_all = shape(rot_cor,0)
!  ! Temporary variables
!  COMPLEX*16 :: TC(-3:3,7,0:3)             ! sph->cubic convertion matrix
!  REAL*8     :: R0(3,3,natom), Runity(3,3), Rd(5,5)
!  COMPLEX*16 :: Rc3(3,3), Rc5(5,5), Uc_rotate(ndim,ndim)
!  INTEGER    :: iatom, r, isort, l, m, lm, m1, i, j
!! Here we create Rotation matrix R0(3,3) which rotates each atom in space as specified in the input file
!  Runity=0;  Runity(1,1)=1.;  Runity(2,2)=1.;  Runity(3,3)=1.
!  DO iatom=1,natom
!     R0(:,:,iatom) = Runity
!     DO r=1,rot_num(iatom)
!        CALL Rotatea(rot_vec_angle(4,r,iatom), rot_vec_angle(1:3,r,iatom), R0(:,:,iatom))
!     ENDDO
!  ENDDO
!!!! T2C is matrix to covert the self-energy or Hamiltonian 
!!!! matrixes from spherical to cubic: 
!!!! SE_c=(T2C)^{+}*SE_s*T2C
!!!! SE_s=T2C*SE_c*(T2C)^{+}
!!!! Then there is no need to warry about which block of the SE
!!!! one needs to convert
!
!!!! Builts transformation matrix to convert cubic harmonics to spherical harmonics
!  CALL TRANSC(TC)
!  T2C = 0.0 
!  DO iatom = 1,natom
!     isort = is(iatom)
!     DO l = 0,ilsb(isort)
!        if (.NOT.ltpan(l,iatom)) cycle
!        lm=l**2+1
!        m = lmtind(lm,iatom)
!        m1=m+2*l
!        IF      (l .EQ. 0) THEN
!           T2C(m:m1,m:m1) = (1.0,0.0)
!        ELSE IF (l .EQ. 1) THEN
!           IF (.NOT.rotateQ(iatom)) THEN
!              T2C(m:m1,m:m1) = TC(-1:1,1:3,1)
!           ELSE
!              Rc3 = R0(:,:,iatom)
!              T2C(m:m1,m:m1) = matmul(TC(-1:1,1:3,1),transpose(conjg(Rc3)))
!           ENDIF
!        ELSE IF (l .EQ. 2) THEN
!           IF (.NOT.rotateQ(iatom)) THEN
!              T2C(m:m1,m:m1) = TC(-2:2,1:5,2)
!           ELSE
!              CALL Rotate_d(Rd, R0(:,:,iatom))
!              Rc5 = Rd
!              T2C(m:m1,m:m1) = matmul(TC(-2:2,1:5,2),transpose(conjg(Rc5)))
!           ENDIF
!        ELSE 
!           T2C(m:m1,m:m1) = TC(-3:3,1:7,3)
!           ! Rotation for f-s not yet implemented
!        ENDIF
!     ENDDO
!  ENDDO
!  IF (ndim.gt.norb) T2C(norb+1:,norb+1:) = T2C(1:norb,1:norb)
!!!! At the end we rotate correlated block according to input cix file.
!  IF (add_rot) THEN
!     Uc_rotate = 0
!     DO i=1,ndim
!        IF (.NOT.bcorrelated(i)) Uc_rotate(i,i)=1
!     ENDDO
!     DO i=1,ncorr_all
!        DO j=1,ncorr_all
!           Uc_rotate(correlated(i),correlated(j)) = rot_cor(i,j)
!        ENDDO
!     ENDDO
!     T2C = matmul(T2C, Uc_rotate)
!  ENDIF
!END SUBROUTINE cmp_t2c

SUBROUTINE GET_TCUBC(string)
  IMPLICIT NONE
! Passed variable
  CHARACTER*280, intent(out) :: string
! temporary
  CHARACTER*10 :: tcubc(-3:3,0:3,2)
  INTEGER      :: i, l, m
  CALL INTCUBC(TCUBC)
  DO l=0,3
     DO m=-3,3
        string  = TRIM(string) //" "// TRIM(tcubc(m,l,1))
     ENDDO
  ENDDO
END SUBROUTINE GET_TCUBC

SUBROUTINE INTCUBC(TCUBC)
  IMPLICIT NONE
! Passed variable
  CHARACTER*10, intent(out) :: TCUBC(-3:3,0:3,2)
! Local variables
  INTEGER I, L, M
  !**************************************************************
  !*   Initializes character string.                            *
  !**************************************************************
  DO I=1,2
     DO L=0,3
        DO M=-3,3
           TCUBC( M,L,I)='  error!  '
        ENDDO
     ENDDO
  ENDDO
  ! Cubic harmonics
  ! S-wave
  TCUBC( 0,0,1)='     s    '
  ! P-wave
  TCUBC(-1,1,1)='     x    '
  TCUBC( 0,1,1)='     y    '
  TCUBC(+1,1,1)='     z    '
  ! D-wave
  TCUBC(-2,2,1)='     yz   '
  TCUBC(-1,2,1)='     zx   '
  TCUBC( 0,2,1)='     xy   '
  TCUBC(+1,2,1)='   x2-y2  '
  TCUBC(+2,2,1)='   3z2-1  '
  ! F-wave
  TCUBC(-3,3,1)=' x(5x2-3) '
  TCUBC(-2,3,1)=' y(5y2-3) '
  TCUBC(-1,3,1)=' z(5z2-3) '
  TCUBC( 0,3,1)=' y(x2-z2) '
  TCUBC( 1,3,1)=' z(x2-y2) '
  TCUBC(+2,3,1)=' x(y2-z2) '
  TCUBC(+3,3,1)='    xyz   '
  ! Spherical harmonics:
  ! S-wave
  TCUBC( 0,0,2)='     0    '
  ! P-wave
  TCUBC(-1,1,2)='    -1    '
  TCUBC( 0,1,2)='     0    '
  TCUBC(+1,1,2)='    +1    '
  ! D-wave
  TCUBC(-2,2,2)='    -2    '
  TCUBC(-1,2,2)='    -1    '
  TCUBC( 0,2,2)='     0    '
  TCUBC(+1,2,2)='    +1    '
  TCUBC(+2,2,2)='    +2    '
  ! F-wave
  TCUBC(-3,3,2)='    -3    '
  TCUBC(-2,3,2)='    -2    '
  TCUBC(-1,3,2)='    -1    '
  TCUBC( 0,3,2)='     0    '
  TCUBC(+1,3,2)='    +1    '
  TCUBC(+2,3,2)='    +2    '
  TCUBC(+3,3,2)='    +3    '
  !
END  SUBROUTINE INTCUBC

SUBROUTINE TRANSC(TC)
!**********************************************************
!*  Builts transformation matrix which transforms         *
!*  cubic harmonics to spherical harmonics.               *
!*  up to Lmax = 3 (f - electrons) ONLY!                  *
!*  Transformation matrix TC(mu,i) is given by            *
!*       _   --              _                            *
!*  Y   (r)= >  TC(mu,i)*Y  (r)                           *
!*   1mu     --           1i                              *
!*          i=xyz                                         *
!**********************************************************
  IMPLICIT REAL*8 (A-H,O-Z)
  COMPLEX*16, intent(out) :: TC(-3:3,7,0:3)
  COMPLEX*16 :: T1(-1:1,3),T2(-2:2,5),T3(-3:3,7)
  DO L=0,3
     DO M=-3,3
        DO K=1,7
           TC(M,K,L)=(0.D0,0.D0)
        ENDDO
     ENDDO
  ENDDO
  CALL TRANS1(T1,1)
  CALL TRANS2(T2,1)
  CALL TRANS3(T3,1)
  !  S-WAVE
  TC(0,1,0)=(1.D0,0.D0)
  !  P-WAVE
  DO MY=-1,1
     DO IX=1,3
	TC(MY,IX,1)=T1(MY,IX)
     ENDDO
  ENDDO
  !  D-WAVE
  DO MY=-2,2
     DO IX=1,5
        TC(MY,IX,2)=T2(MY,IX)
     ENDDO
  ENDDO
  !  F-WAVE
  DO MY=-3,3
     DO IX=1,7
        TC(MY,IX,3)=T3(MY,IX)
     ENDDO
  ENDDO
END SUBROUTINE TRANSC

SUBROUTINE TRANS1(A1,KEY)
!**********************************************************
!*  Generates transformation matrix from complex          *
!*  spherical harmonic to cubic harmonic for l=1.         *
!*                                                        *
!*  Complex spherical harmonics are (after Varshalovich): *
!*       _                                                *
!*  Y   (r)=+1/sqrt(2)*(x-iy)/r * a                       *
!*   1-1                                                  *
!*       _                                                *
!*  Y   (r)= z/r                * a                       *
!*   1 0                                                  *
!*       _                                                *
!*  Y   (r)=-1/sqrt(2)*(x+iy)/r * a                       *
!*   1+1                                                  *
!*                                                        *
!*  Cubic harmonics are:                                  *
!*       _                                                *
!*  Y   (r)= x/r * a                                      *
!*   1 1                                                  *
!*       _                                                *
!*  Y   (r)= y/r * a                                      *
!*   1 2                                                  *
!*       _                                                *
!*  Y   (r)= z/r * a                                      *
!*   1 3                                                  *
!*  where a=sqrt(3/pi)/2 is a normalization constant.     *
!*                                                        *
!*  For a=1 it is equivalent to transition between        *
!*  cartesian r=(x,y,z) and cyclic r=(r  ,r  ,r  )        *
!*  coordinates of vector r.           -1   0  +1         *
!*                                                        *
!*  Transformation matrix T(mu,i) is given by             *
!*       _   --              _                            *
!*  Y   (r)= >   T(mu,i)*Y  (r)                           *
!*   1mu     --           1i                              *
!*          i=xyz                                         *
!*                                                        *
!*  Returns T if key=1, returns reciprocal(T) if key=2.   *
!**********************************************************
  IMPLICIT REAL*8 (A-H,O-Z)
  COMPLEX*16 T1(-1:1,3),T1R(3,-1:1),WORK(3),A1(3,3)
  ! mu=-1
  T1(-1,1)=1.D0/SQRT(2.D0)
  T1(-1,2)=1.D0/SQRT(2.D0)*(0.D0,-1.D0)
  T1(-1,3)=(0.D0,0.D0)
  ! mu=0
  T1( 0,1)=(0.D0,0.D0)
  T1( 0,2)=(0.D0,0.D0)
  T1( 0,3)=(1.D0,0.D0)
  ! mu=+1
  T1(+1,1)=-1.D0/SQRT(2.D0)
  T1(+1,2)=-1.D0/SQRT(2.D0)*(0.D0,1.D0)
  T1(+1,3)=(0.D0,0.D0)
  ! SET RECIPROCAL MATRIX
  DO M=-1,1
     DO I=1,3
        T1R(I,M)=DCONJG(T1(M,I))
     ENDDO
  ENDDO
  IF(KEY.EQ.1)THEN
     !  RETURN DIRECT MATRIX
     DO M=-1,1
        DO I=1,3
           A1(M+2,I)=T1(M,I)
        ENDDO
     ENDDO
  ELSEIF(KEY.EQ.2)THEN
     !  RETURN RECIPROCAL MATRIX
     DO M=-1,1
        DO I=1,3
           A1(I,M+2)=T1R(I,M)
        ENDDO
     ENDDO
  ENDIF
  RETURN
END SUBROUTINE TRANS1

SUBROUTINE TRANS2(A2,KEY)
!**********************************************************
!*  Generates transformation matrix from complex          *
!*  spherical harmonic to cubic harmonic for l=2.         *
!*                                                        *
!*  Complex spherical harmonics are (after Varshalovich)  *
!*       _                                                *
!*  Y   (r)=+1/sqrt(2)*(x^2-y^2-2ixy)/r^2 * a             *
!*   2-2                                                  *
!*       _                                                *
!*  Y   (r)=+1/sqrt(2)*(2zx-2izy)/r^2     * a             *
!*   2-1                                                  *
!*       _                                                *
!*  Y   (r)=(3z^2/r^2-1)/sqrt(3)          * a             *
!*   2 0                                                  *
!*       _                                                *
!*  Y   (r)=-1/sqrt(2)*(2zx+2izy)/r^2     * a             *
!*   2 1                                                  *
!*       _                                                *
!*  Y   (r)=+1/sqrt(2)*(x^2-y^2+2ixy)/r^2 * a             *
!*   2 2                                                  *
!*                                                        *
!*  Cubic harmonics are:                                  *
!*       _                                                *
!*  Y   (r)= 2yz/r^2              * a                     *
!*   2 1                                                  *
!*       _                                                *
!*  Y   (r)= 2zx/r^2              * a                     *
!*   2 2                                                  *
!*       _                                                *
!*  Y   (r)= 2xy/r^2              * a                     *
!*   2 3                                                  *
!*       _                                                *
!*  Y   (r)= (x^2-y^2)/r^2        * a                     *
!*   2 4                                                  *
!*       _                                                *
!*  Y   (r)= (3z^2/r^2-1)/sqrt(3) * a                     *
!*   2 5                                                  *
!*                                                        *
!*  where a=sqrt(3*5/pi)/4 is a normalization constant.   *
!*                                                        *
!*  Transformation matrix T(mu,i) is given by             *
!*       _   --              _                            *
!*  Y   (r)= >   T(mu,i)*Y  (r)                           *
!*   2mu     --           2i                              *
!*          i=xyz                                         *
!*                                                        *
!*  Returns T if key=1, returns reciprocal(T) if key=2.   *
!**********************************************************
  IMPLICIT REAL*8 (A-H,O-Z)
  COMPLEX*16 T2(-2:2,5),T2R(5,-2:2),WORK(5),A2(5,5)
  ! mu=-2
  T2(-2,1)=(0.D0,0.D0)
  T2(-2,2)=(0.D0,0.D0)
  T2(-2,3)=+1.D0/SQRT(2.D0)*(0.D0,-1.D0)
  T2(-2,4)=+1.D0/SQRT(2.D0)
  T2(-2,5)=(0.D0,0.D0)
  ! mu=+2
  T2(+2,1)=(0.D0,0.D0)
  T2(+2,2)=(0.D0,0.D0)
  T2(+2,3)=+1.D0/SQRT(2.D0)*(0.D0,+1.D0)
  T2(+2,4)=+1.D0/SQRT(2.D0)
  T2(+2,5)=(0.D0,0.D0)
  ! mu=-1
  T2(-1,1)=+1.D0/SQRT(2.D0)*(0.D0,-1.D0)
  T2(-1,2)=+1.D0/SQRT(2.D0)
  T2(-1,3)=(0.D0,0.D0)
  T2(-1,4)=(0.D0,0.D0)
  T2(-1,5)=(0.D0,0.D0)
  ! mu=+1
  T2(+1,1)=-1.D0/SQRT(2.D0)*(0.D0,+1.D0)
  T2(+1,2)=-1.D0/SQRT(2.D0)
  T2(+1,3)=(0.D0,0.D0)
  T2(+1,4)=(0.D0,0.D0)
  T2(+1,5)=(0.D0,0.D0)
  ! mu= 0
  T2( 0,1)=(0.D0,0.D0)
  T2( 0,2)=(0.D0,0.D0)
  T2( 0,3)=(0.D0,0.D0)
  T2( 0,4)=(0.D0,0.D0)
  T2( 0,5)=(1.D0,0.D0)
  ! SET RECIPROCAL MATRIX
  DO M=-2,2
     DO I=1,5
        T2R(I,M)=DCONJG(T2(M,I))
     ENDDO
  ENDDO
  IF(KEY.EQ.1)THEN
     !  RETURN DIRECT MATRIX
     DO M=-2,2
        DO I=1,5
           A2(M+3,I)=T2(M,I)
        ENDDO
     ENDDO
  ELSEIF(KEY.EQ.2)THEN
     !  RETURN RECIPROCAL MATRIX
     DO M=-2,2
        DO I=1,5
           A2(I,M+3)=T2R(I,M)
        ENDDO
     ENDDO
  ENDIF
  RETURN
END  SUBROUTINE TRANS2


SUBROUTINE TRANS3(A3,KEY)
!**********************************************************
!*  Generates transformation matrix from complex          *
!*  spherical harmonic to cubic harmonic for l=3.         *
!*                                                        *
!*  Complex spherical harmonics are (after Varshalovich)  *
!*       _                                                *
!*  Y   (r)=+sqrt( 5/16)*(x-iy)^3/r^3                * a  *
!*   3-3                                                  *
!*       _                                                *
!*  Y   (r)=+sqrt(15/ 8)*z(x-iy)^2/r^3               * a  *
!*   3-2                                                  *
!*       _                                                *
!*  Y   (r)=+sqrt( 3/16)*(5z^2-r^2)(x-iy)^2/r^3      * a  *
!*   3-1                                                  *
!*       _                                                *
!*  Y   (r)=+sqrt( 1/ 4)*(5z^2-3r^2)z/r^3            * a  *
!*   3 0                                                  *
!*       _                                                *
!*  Y   (r)=-sqrt( 3/16)*(5z^2-r^2)(x+iy)^2/r^3      * a  *
!*   3+1                                                  *
!*       _                                                *
!*  Y   (r)=+sqrt(15/ 8)*z(x+iy)^2/r^3               * a  *
!*   3+2                                                  *
!*       _                                                *
!*  Y   (r)=-sqrt( 5/16)*(x+iy)^3/r^3                * a  *
!*   3+3                                                  *
!*                                                        *
!*  Cubic harmonics are:                                  *
!*       _                                                *
!*  Y   (r)=+sqrt( 1/ 4)*(5x^2-3r^2)x/r^3            * a  *
!*   3 1                                                  *
!*       _                                                *
!*  Y   (r)=+sqrt( 1/ 4)*(5y^2-3r^2)y/r^3            * a  *
!*   3 2                                                  *
!*       _                                                *
!*  Y   (r)=+sqrt( 1/ 4)*(5z^2-3r^2)z/r^3            * a  *
!*   3 3                                                  *
!*       _                                                *
!*  Y   (r)=+sqrt(15/ 4)*(x^2-z^2)y/r^3              * a  *
!*   3 4                                                  *
!*       _                                                *
!*  Y   (r)=+sqrt(15/ 4)*(x^2-y^2)z/r^3              * a  *
!*   3 5                                                  *
!*       _                                                *
!*  Y   (r)=+sqrt(15/ 4)*(y^2-z^2)x/r^3              * a  *
!*   3 6                                                  *
!*       _                                                *
!*  Y   (r)=+sqrt(15   )*xyz/r^3                     * a  *
!*   3 7                                                  *
!*                                                        *
!*  where a=sqrt(7/4/pi) is a normalization constant.     *
!*                                                        *
!*  Transformation matrix T(mu,i) is given by             *
!*       _   --              _                            *
!*  Y   (r)= >   T(mu,i)*Y  (r)                           *
!*   2mu     --           2i                              *
!*          i=xyz                                         *
!*                                                        *
!*  Returns T if key=1, returns reciprocal(T) if key=2.   *
!**********************************************************
  IMPLICIT REAL*8 (A-H,O-Z)
  COMPLEX*16 T3(-3:3,7),T3R(7,-3:3),WORK(7),A3(7,7)
  ! SET RECIPROCAL TRANSFORMATION FIRST
  T3R(1,-3)=+sqrt(5.d0/16.d0)*(1.d0,0.d0)
  T3R(1,-2)=(0.d0,0.d0)
  T3R(1,-1)=-sqrt(3.d0/16.d0)*(1.d0,0.d0)
  T3R(1, 0)=(0.d0,0.d0)
  T3R(1,+1)=+sqrt(3.d0/16.d0)*(1.d0,0.d0)
  T3R(1,+2)=(0.d0,0.d0)
  T3R(1,+3)=-sqrt(5.d0/16.d0)*(1.d0,0.d0)
  !
  T3R(2,-3)=-sqrt(5.d0/16.d0)*(0.d0,1.d0)
  T3R(2,-2)=(0.d0,0.d0)
  T3R(2,-1)=-sqrt(3.d0/16.d0)*(0.d0,1.d0)
  T3R(2, 0)=(0.d0,0.d0)
  T3R(2,+1)=-sqrt(3.d0/16.d0)*(0.d0,1.d0)
  T3R(2,+2)=(0.d0,0.d0)
  T3R(2,+3)=-sqrt(5.d0/16.d0)*(0.d0,1.d0)
  !
  T3R(3,-3)=(0.d0,0.d0)
  T3R(3,-2)=(0.d0,0.d0)
  T3R(3,-1)=(0.d0,0.d0)
  T3R(3, 0)=(1.d0,0.d0)
  T3R(3,+1)=(0.d0,0.d0)
  T3R(3,+2)=(0.d0,0.d0)
  T3R(3,+3)=(0.d0,0.d0)
  !
  T3R(4,-3)=+sqrt(3.d0/16.d0)*(0.d0,1.d0)
  T3R(4,-2)=(0.d0,0.d0)
  T3R(4,-1)=-sqrt(5.d0/16.d0)*(0.d0,1.d0)
  T3R(4, 0)=(0.d0,0.d0)
  T3R(4,+1)=-sqrt(5.d0/16.d0)*(0.d0,1.d0)
  T3R(4,+2)=(0.d0,0.d0)
  T3R(4,+3)=+sqrt(3.d0/16.d0)*(0.d0,1.d0)
  !
  T3R(5,-3)=(0.d0,0.d0)
  T3R(5,-2)=+sqrt(1.d0/2.0d0)*(1.d0,0.d0)
  T3R(5,-1)=(0.d0,0.d0)
  T3R(5, 0)=(0.d0,0.d0)
  T3R(5,+1)=(0.d0,0.d0)
  T3R(5,+2)=+sqrt(1.d0/2.0d0)*(1.d0,0.d0)
  T3R(5,+3)=(0.d0,0.d0)
  !
  T3R(6,-3)=-sqrt(3.d0/16.d0)*(1.d0,0.d0)
  T3R(6,-2)=(0.d0,0.d0)
  T3R(6,-1)=-sqrt(5.d0/16.d0)*(1.d0,0.d0)
  T3R(6, 0)=(0.d0,0.d0)
  T3R(6,+1)=+sqrt(5.d0/16.d0)*(1.d0,0.d0)
  T3R(6,+2)=(0.d0,0.d0)
  T3R(6,+3)=+sqrt(3.d0/16.d0)*(1.d0,0.d0)
  !
  T3R(7,-3)=(0.d0,0.d0)
  T3R(7,-2)=+sqrt(1.d0/2.0d0)*(0.d0,1.d0)
  T3R(7,-1)=(0.d0,0.d0)
  T3R(7, 0)=(0.d0,0.d0)
  T3R(7,+1)=(0.d0,0.d0)
  T3R(7,+2)=-sqrt(1.d0/2.0d0)*(0.d0,1.d0)
  T3R(7,+3)=(0.d0,0.d0)
  !
  DO M=-3,3
     DO I=1,7
        T3(M,I)=DCONJG(T3R(I,M))
     ENDDO
  ENDDO
  IF(KEY.EQ.1)THEN
     !  RETURN DIRECT MATRIX
     DO M=-3,3
        DO I=1,7
           A3(M+4,I)=T3(M,I)
        ENDDO
     ENDDO
  ELSEIF(KEY.EQ.2)THEN
     !  RETURN RECIPROCAL MATRIX
     DO M=-3,3
        DO I=1,7
           A3(I,M+4)=T3R(I,M)
        ENDDO
     ENDDO
  ENDIF
END SUBROUTINE TRANS3

SUBROUTINE Rotate_d(Rd, R)
! Rotation of the d electron cubic harmonic orbitals. Real space 
! rotation matrix is R(3,3). The cubic wave functions are defined 
! as follows:
! b1 = 2yz/r^2 * a
! b2 = 2xz/r^2 * a
! b3 = 2xy/r^2 * a
! b4 = (x^2-y^2)/r^2 * a
! b5 = (3z^2-r^2)/r^2 * a
! Rotation changes them as follows R.b_1 = R(2,i)R(3,j)*2*u_i*u_j/r^2 * a
! which leads to the following equations
! R.b1 = (R(2,2)*R(3,3) + R(2,3)*R(3,2))*b1 + (R(2,1)*R(3,3) + R(2,3)*R(3,1))*b2 
!        + (R(2,1)*R(3,2) + R(2,2)*R(3,1))*b3
!        + (R(2,1)*R(3,1) - R(2,2)*R(3,2))*b4 + sqrt(3.)*R(2,3)*R(3,3)*b5
  IMPLICIT NONE
  REAL*8, intent(in)  :: R(3,3)
  REAL*8, intent(out) :: Rd(5,5)
  !
  Rd(1,1) = R(2,2)*R(3,3) + R(2,3)*R(3,2);
  Rd(1,2) = R(2,1)*R(3,3) + R(2,3)*R(3,1);
  Rd(1,3) = R(2,1)*R(3,2) + R(2,2)*R(3,1);
  Rd(1,4) = R(2,1)*R(3,1) - R(2,2)*R(3,2);
  Rd(1,5) = sqrt(3.)*R(2,3)*R(3,3);
  Rd(2,1) = R(1,2)*R(3,3) + R(1,3)*R(3,2);
  Rd(2,2) = R(1,1)*R(3,3) + R(1,3)*R(3,1);
  Rd(2,3) = R(1,1)*R(3,2) + R(1,2)*R(3,1);
  Rd(2,4) = R(1,1)*R(3,1) - R(1,2)*R(3,2);
  Rd(2,5) = sqrt(3.)*R(1,3)*R(3,3);
  Rd(3,1) = R(1,2)*R(2,3) + R(1,3)*R(2,2);
  Rd(3,2) = R(1,1)*R(2,3) + R(1,3)*R(2,1);
  Rd(3,3) = R(1,1)*R(2,2) + R(1,2)*R(2,1);
  Rd(3,4) = R(1,1)*R(2,1) - R(1,2)*R(2,2);
  Rd(3,5) = sqrt(3.)*R(1,3)*R(2,3);
  Rd(4,1) = R(1,2)*R(1,3) - R(2,2)*R(2,3);
  Rd(4,2) = R(1,1)*R(1,3) - R(2,1)*R(2,3);
  Rd(4,3) = R(1,1)*R(1,2) - R(2,1)*R(2,2);
  Rd(4,4) = 0.5*(R(1,1)**2+R(2,2)**2) -0.5*(R(1,2)**2 + R(2,1)**2);
  Rd(4,5) = 0.5*sqrt(3.)*(R(1,3)**2-R(2,3)**2);
  Rd(5,1) = sqrt(3.)*R(3,2)*R(3,3);
  Rd(5,2) = sqrt(3.)*R(3,1)*R(3,3);
  Rd(5,3) = sqrt(3.)*R(3,1)*R(3,2);
  Rd(5,4) = 0.5*sqrt(3.)*(R(3,1)**2-R(3,2)**2);
  Rd(5,5) = -0.5 + 1.5*R(3,3)**2;
END SUBROUTINE Rotate_d

SUBROUTINE Rotate_nt(t, m, R)
! Rotation matrix R(3,3) around normal vector 
! (m0,m1,m2) and angle t
  IMPLICIT NONE
  REAL*8, intent(in) :: t
  REAL*8, intent(in) :: m(3)
  REAL*8, intent(out):: R(3,3)
  REAL*8 :: norm, c, s, omc
  REAL*8 :: n0(3)
  norm = sqrt(m(1)**2+m(2)**2+m(3)**2)
  n0(1) = m(1)/norm
  n0(2) = m(2)/norm
  n0(3) = m(3)/norm
  c = cos(t)
  s = sin(t) 
  omc=1-c
  R(1,1) = n0(1)*n0(1)*omc + c;
  R(1,2) = n0(1)*n0(2)*omc - n0(3)*s;
  R(1,3) = n0(1)*n0(3)*omc + n0(2)*s;
  R(2,1) = n0(2)*n0(1)*omc + n0(3)*s;
  R(2,2) = n0(2)*n0(2)*omc + c;
  R(2,3) = n0(2)*n0(3)*omc - n0(1)*s;
  R(3,1) = n0(3)*n0(1)*omc - n0(2)*s;
  R(3,2) = n0(3)*n0(2)*omc + n0(1)*s;
  R(3,3) = n0(3)*n0(3)*omc + c;
END SUBROUTINE Rotate_nt

SUBROUTINE Rotatea(t, m, R)
! This subroutine takes "old" rotation R and adds 
! "new" rotation which is specified by vector (m0,m1,m2)
! and angle t.
  IMPLICIT NONE
  REAL*8, intent(in) :: t
  REAL*8, intent(in) :: m(3)
  REAL*8, intent(out):: R(3,3)
  REAL*8 :: Rt(3,3)

  CALL Rotate_nt(t,m,Rt)
  R = matmul(Rt,R)
END SUBROUTINE Rotatea

