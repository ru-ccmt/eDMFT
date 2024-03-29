!**************************************************!
!******* § Copyright by Kristjan Haule, 2002 ******!
!**************************************************!
SUBROUTINE eigsys(ham, zek, evl, evr, ndim)
!!!----------------------------------------------------------------------!!!
!!! This routine solves the generalized right/left eigenvalue problems   !!!
!!!       H.Ev^{R} = E Ev^{R}      and    Ev^{L}H = E Ev^{L}             !!!
!!! ( . denotes matrix multiplication )  where H is a general matrix and.!!!
!!!                                                                      !!!
!!! The problem is solved in several steps:                              !!!
!!!                                                                      !!!
!!!                                                                      !!!
!!! Variable contents on entry:                                          !!!
!!! ===========================                                          !!!
!!!                                                                      !!!
!!! ham        : Hamiltonian                                             !!!
!!! ndim       : Dimension of matrices                                   !!!
!!! fh_stderr  : Standard error file handle (Unit number)                !!!
!!!                                                                      !!!
!!! Variable contents on exit:                                           !!!
!!! ==========================                                           !!!
!!!                                                                      !!!
!!! evr        : Right eigenvectors as columns                           !!!
!!! evl        : Left eigenvectors as rows                               !!!
!!! zek        : Eigenvalues                                             !!!
!!!----------------------------------------------------------------------!!!
  IMPLICIT NONE
  !---------- Passed variables ----------
  COMPLEX*16, INTENT(inout):: ham(ndim,ndim)  ! Hamiltonian / overwritten
  COMPLEX*16, INTENT(out) :: zek(ndim)       ! Vector of eigenvalues 
  COMPLEX*16, INTENT(out) :: evl(ndim,ndim)   ! left eigenvectors
  COMPLEX*16, INTENT(out) :: evr(ndim,ndim)   ! right eigenvectors
  INTEGER, INTENT(in)     :: ndim            ! Dimension of hamiltonian
  !f2py integer intent(hide), depend(ham)  :: ndim=shape(ham,0)
  !---------- Local variables ----------
  REAL*8, PARAMETER :: smalleps = 1e-5
  CHARACTER*1, PARAMETER :: uplo  = "U"  ! Job parameter for zpotrf/ztrtri
  CHARACTER*1, PARAMETER :: diag  = "N"  ! Job parameter for ztrtri
  CHARACTER*1, PARAMETER :: jobvl = "V"  ! Jop parameter for zgeev
  CHARACTER*1, PARAMETER :: jobvr = "V"  ! Job parameter for zgeev
  CHARACTER*1, PARAMETER :: transn  = "N" ! Job parameter for ztrmm
  CHARACTER*1, PARAMETER :: transc  = "C" ! Job parameter for ztrmm
  CHARACTER*1, PARAMETER :: left   = "L" ! Job parameter for ztrmm
  CHARACTER*1, PARAMETER :: right  = "R" ! Job parameter for ztrmm
  INTEGER :: ierr                        ! Error parameter for lapack
  INTEGER :: p,q                         ! Loop index 
  COMPLEX*16 :: ctmp                     ! Temporary variable 
  !COMPLEX*16 :: htmp(ndim,ndim)  ! Transformed Hamiltonian
  COMPLEX*16 :: scaler(ndim)     ! Array of normalization parameters
  COMPLEX*16 :: cworkvec(8*ndim) ! Work array for zgeev
  REAL*8     :: rworkvec(8*ndim) ! Work array for zgeev
  INTEGER    :: idxarr(ndim)

  !htmp = ham
  !========== Step 4, Solve eigenvalue problem, H'v' = Ev' ==========
  CALL zgeev(jobvl,jobvr,ndim,ham,ndim,zek,evl,ndim,evr,ndim,cworkvec,8*ndim,rworkvec,ierr)
  IF(ierr.NE.0) THEN
     WRITE(6,*) 'Error code of zgeev ', ierr
     WRITE(0,*) 'Error in dia_gho! Stopp the code!'
  ENDIF
  ! transpose left eigenvectors
  evl = dconjg(TRANSPOSE(evl))

  !========== Step 5, Make degenerate eigenvectors orthogonal
  DO q=1,ndim
     DO p=1,q-1
        IF (abs(zek(p)-zek(q)).LT.smalleps .AND. abs(scalprod(evl(p,:),evr(:,q),ndim)) .GT.smalleps) THEN
           evr(:,q) = evr(:,q) - scalprod(evl(p,:),evr(:,q),ndim)/scalprod(evl(p,:),evr(:,p),ndim) * evr(:,p)
        ENDIF
     ENDDO
     DO p=1,q-1
        IF (abs(zek(p)-zek(q)).LT.smalleps .AND. abs(scalprod(evl(q,:),evr(:,p),ndim)) .GT.smalleps) THEN
           evl(q,:) = evl(q,:) - scalprod(evl(q,:),evr(:,p),ndim)/scalprod(evl(p,:),evr(:,p),ndim) * evl(p,:)
        ENDIF
     ENDDO
  ENDDO

  !========= Step 6, Normalize eigenvectors
  DO p = 1,ndim
     ctmp = 0.d0
     DO q = 1,ndim
        ctmp = ctmp+evl(p,q)*evr(q,p)
     ENDDO
     scaler(p) = SQRT(ctmp)
  ENDDO
  DO p = 1,ndim
     evl(p,:) = evl(p,:)/scaler(p)
     evr(:,p) = evr(:,p)/scaler(p)
  ENDDO

  !========= Sorting acording to real parts
  CALL eig_order_real_part(zek, idxarr, ndim)
  CALL permute_eigensystem(idxarr, zek, evl, evr, ndim)

  RETURN

CONTAINS
  COMPLEX*16 FUNCTION scalprod(a,b,ndim)
    IMPLICIT NONE
    COMPLEX*16 :: a(:), b(:)
    INTEGER    :: ndim
    INTEGER    :: i
    scalprod = 0.0
    DO i=1,ndim
       scalprod = scalprod + a(i)*b(i)
    ENDDO
  END FUNCTION scalprod

END SUBROUTINE eigsys


SUBROUTINE eigvals(ham, zek, ndim)
!!!----------------------------------------------------------------------!!!
!!! Eigenvalues only, no eigenvectors                                    !!!
!!! This routine solves the generalized right/left eigenvalue problems   !!!
!!!       H.Ev^{R} = E Ev^{R}      and    Ev^{L}H = E Ev^{L}             !!!
!!! ( . denotes matrix multiplication )  where H is a general matrix and.!!!
!!!                                                                      !!!
!!! The problem is solved in several steps:                              !!!
!!!                                                                      !!!
!!!                                                                      !!!
!!! Variable contents on entry:                                          !!!
!!! ===========================                                          !!!
!!!                                                                      !!!
!!! ham        : Hamiltonian                                             !!!
!!! ndim       : Dimension of matrices                                   !!!
!!! fh_stderr  : Standard error file handle (Unit number)                !!!
!!!                                                                      !!!
!!! Variable contents on exit:                                           !!!
!!! ==========================                                           !!!
!!!                                                                      !!!
!!! evr        : Right eigenvectors as columns                           !!!
!!! evl        : Left eigenvectors as rows                               !!!
!!! zek        : Eigenvalues                                             !!!
!!!----------------------------------------------------------------------!!!
  IMPLICIT NONE
  !---------- Passed variables ----------
  COMPLEX*16, INTENT(in)  :: ham(ndim,ndim)  ! Hamiltonian
  COMPLEX*16, INTENT(out) :: zek(ndim)       ! Vector of eigenvalues 
  INTEGER, INTENT(in)     :: ndim            ! Dimension of hamiltonian
  !f2py integer intent(hide), depend(ham)  :: ndim=shape(ham,0)
  !---------- Local variables ----------
  REAL*8, PARAMETER :: smalleps = 1e-5
  CHARACTER*1, PARAMETER :: uplo  = "U"  ! Job parameter for zpotrf/ztrtri
  CHARACTER*1, PARAMETER :: diag  = "N"  ! Job parameter for ztrtri
  CHARACTER*1, PARAMETER :: jobvl = "N"  ! Jop parameter for zgeev
  CHARACTER*1, PARAMETER :: jobvr = "N"  ! Job parameter for zgeev
  CHARACTER*1, PARAMETER :: transn  = "N" ! Job parameter for ztrmm
  CHARACTER*1, PARAMETER :: transc  = "C" ! Job parameter for ztrmm
  CHARACTER*1, PARAMETER :: left   = "L" ! Job parameter for ztrmm
  CHARACTER*1, PARAMETER :: right  = "R" ! Job parameter for ztrmm
  INTEGER :: ierr                        ! Error parameter for lapack
  COMPLEX*16 :: cworkvec(8*ndim) ! Work array for zgeev
  REAL*8     :: rworkvec(8*ndim) ! Work array for zgeev
  COMPLEX*16 :: evl(1,1)   ! left eigenvectors
  COMPLEX*16 :: evr(1,1)   ! right eigenvectors
  INTEGER    :: idxarr(ndim)

  !========== Step 4, Solve eigenvalue problem, H'v' = Ev' ==========
  CALL zgeev(jobvl,jobvr,ndim,ham,ndim,zek,evl,ndim,evr,ndim,cworkvec,8*ndim,rworkvec,ierr)
  IF(ierr.NE.0) THEN
     WRITE(6,*) 'Error code of zgeev ', ierr
     WRITE(0,*) 'Error in dia_gho! Stopp the code!'
  ENDIF

  !========= Sorting acording to real parts
  CALL eig_order_real_part(zek, idxarr, ndim)
  CALL permute_eigenvals(idxarr, zek, ndim)

  RETURN
END SUBROUTINE eigvals

!===========================================================================
SUBROUTINE eig_order_real_part(ev, idxarr, ndim)
  IMPLICIT NONE 
!!!-----------------------------------------------------------------!!!
!!! This routine sorts complex eigenvalues of a matrix according to !!!
!!! its real parts with the smallest in the first slot and reorders !!!
!!! the matrices of left (row) and right (column) eigenvectors in a !!!
!!! corresponding manner.                                           !!!
!!!-----------------------------------------------------------------!!!
  !---------- Passed variables ----------
  COMPLEX*16, intent(in) :: ev(ndim)         ! Array of eigenvalues
  INTEGER, intent(out)   :: idxarr(ndim)     ! Index array which gives proper order
  INTEGER :: ndim                            ! Dimension of matrices 
  !f2py integer intent(hide), depend(ev)  :: ndim=shape(ev,0)
  !---------- Parameters ----------
  REAL*8, PARAMETER :: maxval = 1000.d0
  !---------- Local variables ----------
  LOGICAL, ALLOCATABLE :: sorted(:)
  REAL*8,  ALLOCATABLE :: sortonr(:)
  INTEGER :: p
  INTEGER :: q
  INTEGER :: idx
  REAL*8  :: min
  !---------- Allocate dynamic memory storage ----------
  ALLOCATE(sortonr(ndim), sorted(ndim))
  !---------- Initialize arrays ----------
  idxarr = 0
  sortonr = DBLE(ev)
  sorted = .FALSE.
  !---------- Create index array for real value ----------
  sorted = .FALSE.
  DO p = 1,ndim
     min = maxval
     DO q = 1,ndim
        IF(.NOT.sorted(q).AND.min.GT.sortonr(q)) THEN
           min = sortonr(q)
           idx = q
        ENDIF
     ENDDO
     idxarr(p) = idx
     sorted(idx) = .TRUE.
  ENDDO
  DEALLOCATE(sortonr, sorted)
  RETURN
END SUBROUTINE eig_order_real_part

SUBROUTINE permute_eigensystem(idxarr, ev, evl, evr, ndim)
  IMPLICIT NONE
  !---------- Passed variables ----------
  INTEGER, intent(in)       :: idxarr(ndim)     ! Index array which gives proper order
  COMPLEX*16, intent(inout) :: ev(ndim)         ! Array of eigenvalues
  COMPLEX*16, intent(inout) :: evl(ndim,ndim)   ! Matrix of left eigenvectors  (row)
  COMPLEX*16, intent(inout) :: evr(ndim,ndim)   ! Matrix of right eigenvectors (column)
  INTEGER :: ndim                               ! Dimension of matrices 
  !f2py integer intent(hide), depend(ev)  :: ndim=shape(ev,0)
  !---------- Local variables ------------------
  INTEGER :: p
  COMPLEX*16, ALLOCATABLE :: eval(:)
  COMPLEX*16, ALLOCATABLE :: evec(:,:)
  ALLOCATE(eval(ndim), evec(ndim,ndim)) 
  !---------- Permute the eigenvalues ----------
  DO p = 1,ndim
     eval(p) = ev(idxarr(p))
  ENDDO
  ev = eval
  !---------- Permute the right eigenvectors ----------
  DO p = 1,ndim
     evec(:,p) = evr(:,idxarr(p))
  ENDDO
  evr = evec
  !---------- Permute the left eigenvectors ----------
  DO p = 1,ndim
     evec(p,:) = evl(idxarr(p),:)
  ENDDO
  evl = evec
  !---------- Deallocate dynamic memory storage ----------
  DEALLOCATE(eval, evec) 
  RETURN 
END SUBROUTINE permute_eigensystem

SUBROUTINE permute_eigenvals(idxarr, ev, ndim)
  IMPLICIT NONE
  !---------- Passed variables ----------
  INTEGER, intent(in)       :: idxarr(ndim)     ! Index array which gives proper order
  COMPLEX*16, intent(inout) :: ev(ndim)         ! Array of eigenvalues
  INTEGER :: ndim                               ! Dimension of matrices 
  !f2py integer intent(hide), depend(ev)  :: ndim=shape(ev,0)
  !---------- Local variables ------------------
  INTEGER :: p
  COMPLEX*16, ALLOCATABLE :: eval(:)
  ALLOCATE(eval(ndim))
  !---------- Permute the eigenvalues ----------
  DO p = 1,ndim
     eval(p) = ev(idxarr(p))
  ENDDO
  ev = eval
  !---------- Deallocate dynamic memory storage ----------
  DEALLOCATE(eval)
  RETURN 
END SUBROUTINE permute_eigenvals

SUBROUTINE zheigsys(ham, zek, evl, evr, ndim)
!!!----------------------------------------------------------------------!!!
!!! This routine solves the hermitian eigenvalue problems                !!!
!!!       H.Ev^{R} = E Ev^{R}                                            !!!
!!! but acts like the non-hermitian solver above                         !!!
!!!----------------------------------------------------------------------!!!
  IMPLICIT NONE
  !---------- Passed variables ----------
  COMPLEX*16, INTENT(inout):: ham(ndim,ndim)  ! Hamiltonian / overwritten with right eigenvectors
  COMPLEX*16, INTENT(out)  :: zek(ndim)       ! Vector of eigenvalues
  COMPLEX*16, INTENT(out) :: evl(ndim,ndim)   ! left eigenvectors
  COMPLEX*16, INTENT(out) :: evr(ndim,ndim)   ! right eigenvectors
  INTEGER, INTENT(in)     :: ndim            ! Dimension of hamiltonian
  ! locals
  INTEGER :: lwork, liwork, lrwork, info
  COMPLEX*16, allocatable :: work(:)
  REAL*8,     allocatable :: rwork(:)
  INTEGER,    allocatable :: iwork(:)
  REAL*8  :: ek(ndim)
  
  lwork  = 1+4*ndim+ndim**2
  liwork = 4+5*ndim
  lrwork = 1+6*ndim+2*ndim**2
  allocate( work(lwork), iwork(liwork), rwork(lrwork))
  
  ! ham = evr * ek * evr^C
  !   or
  ! ek = evr^C * ham * evr
  !
  ! zheevd(jobz, uplo, n, a, lda, w, work, lwork, rwork, lrwork, iwork, liwork, info)
  CALL ZHEEVD('V','U', ndim, ham, ndim, ek, work, lwork, rwork, lrwork, iwork, liwork, info )
  evr(:,:) = ham(:,:)
  evl(:,:) = transpose(conjg(ham))
  zek(:) = ek(:)
  deallocate( work, iwork, rwork )
END SUBROUTINE zheigsys

SUBROUTINE dheigsys(ham, ek, evr, ndim)
!!!----------------------------------------------------------------------!!!
!!! This routine solves the hermitian eigenvalue problems                !!!
!!!       H.Ev^{R} = E Ev^{R}                                            !!!
!!! but acts like the non-hermitian solver above                         !!!
!!!----------------------------------------------------------------------!!!
  IMPLICIT NONE
  !---------- Passed variables ----------
  COMPLEX*16, INTENT(inout):: ham(ndim,ndim)  ! Hamiltonian / overwritten with right eigenvectors
  REAL*8, INTENT(out)     :: ek(ndim)       ! Vector of eigenvalues
  COMPLEX*16, INTENT(out) :: evr(ndim,ndim)   ! right eigenvectors
  INTEGER, INTENT(in)     :: ndim            ! Dimension of hamiltonian
  ! locals
  INTEGER :: lwork, liwork, lrwork, info
  COMPLEX*16, allocatable :: work(:)
  REAL*8,     allocatable :: rwork(:)
  INTEGER,    allocatable :: iwork(:)
  lwork  = 1+4*ndim+ndim**2
  liwork = 4+5*ndim
  lrwork = 1+6*ndim+2*ndim**2
  allocate( work(lwork), iwork(liwork), rwork(lrwork))
  
  ! ham = evr * ek * evr^C
  !   or
  ! ek = evr^C * ham * evr
  !
  ! zheevd(jobz, uplo, n, a, lda, w, work, lwork, rwork, lrwork, iwork, liwork, info)
  CALL ZHEEVD('V','U', ndim, ham, ndim, ek, work, lwork, rwork, lrwork, iwork, liwork, info )
  evr(:,:) = ham(:,:)
  deallocate( work, iwork, rwork )
END SUBROUTINE dheigsys

SUBROUTINE zheigvals(ham, zek, ndim)
  IMPLICIT NONE
  !---------- Passed variables ----------
  COMPLEX*16, INTENT(in)  :: ham(ndim,ndim)  ! Hamiltonian
  COMPLEX*16, INTENT(out) :: zek(ndim)       ! Vector of eigenvalues 
  INTEGER, INTENT(in)     :: ndim            ! Dimension of hamiltonian
  !---------- Local variables ----------
  INTEGER :: lwork, liwork, lrwork, info
  COMPLEX*16, allocatable :: work(:)
  REAL*8,     allocatable :: rwork(:)
  INTEGER,    allocatable :: iwork(:)
  REAL*8  :: ek(ndim)
  lwork  = 1+2*ndim
  liwork = 1+ndim
  lrwork = 1+2*ndim
  allocate( work(lwork), iwork(liwork), rwork(lrwork))
  !
  ! zheevd(jobz, uplo, n, a, lda, w, work, lwork, rwork, lrwork, iwork, liwork, info)
  CALL ZHEEVD('N','U', ndim, ham, ndim, ek, work, lwork, rwork, lrwork, iwork, liwork, info )
  zek(:) = ek(:)
  deallocate( work, iwork, rwork )
END SUBROUTINE zheigvals

SUBROUTINE heigsys(ham, ek, ndim)
!!!----------------------------------------------------------------------!!!
!!! This routine solves the hermitian eigenvalue problems                !!!
!!!       H.Ev^{R} = E Ev^{R}                                            !!!
!!!----------------------------------------------------------------------!!!
  IMPLICIT NONE
  !---------- Passed variables ----------
  COMPLEX*16, INTENT(inout):: ham(ndim,ndim)  ! Hamiltonian / overwritten with right eigenvectors
  REAL*8,     INTENT(out)  :: ek(ndim)       ! Vector of eigenvalues 
  INTEGER, INTENT(in)      :: ndim            ! Dimension of hamiltonian
  ! locals
  INTEGER :: lwork, liwork, lrwork, info
  COMPLEX*16, allocatable :: work(:)
  REAL*8,     allocatable :: rwork(:)
  INTEGER,    allocatable :: iwork(:)
  
  lwork  = 1+6*ndim+2*ndim**2
  liwork = 3+5*ndim
  lrwork = 1+5*ndim+2*ndim**2
  allocate( work(lwork), iwork(liwork) )
  allocate( rwork(lrwork))
  
  ! ham = evr * ek * evr^C
  !   or
  ! ek = evr^C * ham * evr
  !
  CALL ZHEEVD('V','U', ndim, ham, ndim, ek, work, lwork, rwork, lrwork, iwork, liwork, info )
  
  deallocate( rwork )
  deallocate( work, iwork )
END SUBROUTINE heigsys
