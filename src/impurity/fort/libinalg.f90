! @Copyright 2007 Kristjan Haule
! 

subroutine ZProduct_MM(C,A,B,transa,transb, na1, na2, nb1, nb2, nc1, nc2)
  IMPLICIT NONE
  COMPLEX*16, intent(out) :: C(nc1,nc2)
  COMPLEX*16, intent(in)  :: A(na1,na2)
  COMPLEX*16, intent(in)  :: B(nb1,nb2)
  CHARACTER,  intent(in)  :: transa
  CHARACTER,  intent(in)  :: transb
  INTEGER,    intent(in)  :: na1, na2, nb1, nb2, nc1, nc2
  !f2py integer intent(hide), depend(A)  :: na1=shape(A,0)
  !f2py integer intent(hide), depend(A)  :: na2=shape(A,1)
  !f2py integer intent(hide), depend(B)  :: nb1=shape(B,0)
  !f2py integer intent(hide), depend(B)  :: nb2=shape(B,1)  
  ! local variables
  COMPLEX*16 :: alpha_, beta_
  INTEGER :: m, n, k, kp
  alpha_ = 1.
  beta_ = 0
  IF (transa.EQ.'N' .AND. transb.EQ.'N') THEN
     m = na1
     k = na2
     n = nb2
     kp= nb1
  ELSEIF (transa.EQ.'C' .AND. transb.EQ.'N') THEN
     m = na2
     k = na1
     n = nb2
     kp= nb1
  ELSEIF (transa.EQ.'N' .AND. transb.EQ.'C') THEN
     m = na1
     k = na2
     n = nb1
     kp= nb2
  ELSEIF (transa.EQ.'C' .AND. transb.EQ.'C') THEN
     m = na2
     k = na1
     n = nb1
     kp= nb2
  ELSE
     print *, 'The call to ZProduct_MM was wrong!'
     CALL EXIT(1)
  ENDIF
  IF (kp.NE.k) print *, 'Error in ProductMM, sizes not correct!'
  IF (nc1.NE.m) print *, 'Error in ProductMM, sizes not correct!'
  IF (nc2.NE.n) print *, 'Error in ProductMM, sizes not correct!'
  CALL zgemm(transa,transb,m,n,k,alpha_,A,na1,B,nb1,beta_,C,nc1)
end subroutine ZProduct_MM


subroutine ZProduct_sum_MM(C,A,B,transa,transb, na1, na2, nb1, nb2, nc1, nc2)
  IMPLICIT NONE
  COMPLEX*16, intent(out) :: C(nc1,nc2)
  COMPLEX*16, intent(in)  :: A(na1,na2)
  COMPLEX*16, intent(in)  :: B(nb1,nb2)
  CHARACTER,  intent(in)  :: transa
  CHARACTER,  intent(in)  :: transb
  INTEGER,    intent(in)  :: na1, na2, nb1, nb2, nc1, nc2
  !f2py integer intent(hide), depend(A)  :: na1=shape(A,0)
  !f2py integer intent(hide), depend(A)  :: na2=shape(A,1)
  !f2py integer intent(hide), depend(B)  :: nb1=shape(B,0)
  !f2py integer intent(hide), depend(B)  :: nb2=shape(B,1)  
  ! local variables
  COMPLEX*16 :: alpha_, beta_
  INTEGER :: m, n, k, kp
  alpha_ = 1.
  beta_ = 1.
  IF (transa.EQ.'N' .AND. transb.EQ.'N') THEN
     m = na1
     k = na2
     n = nb2
     kp= nb1
  ELSEIF (transa.EQ.'C' .AND. transb.EQ.'N') THEN
     m = na2
     k = na1
     n = nb2
     kp= nb1
  ELSEIF (transa.EQ.'N' .AND. transb.EQ.'C') THEN
     m = na1
     k = na2
     n = nb1
     kp= nb2
  ELSEIF (transa.EQ.'C' .AND. transb.EQ.'C') THEN
     m = na2
     k = na1
     n = nb1
     kp= nb2
  ELSE
     print *, 'The call to ZProduct_MM was wrong!'
     CALL EXIT(1)
  ENDIF
  IF (kp.NE.k) print *, 'Error in ProductMM, sizes not correct!'
  IF (nc1.NE.m) print *, 'Error in ProductMM, sizes not correct!'
  IF (nc2.NE.n) print *, 'Error in ProductMM, sizes not correct!'
  CALL zgemm(transa,transb,m,n,k,alpha_,A,na1,B,nb1,beta_,C,nc1)
end subroutine ZProduct_sum_MM

subroutine ZTransform(A,U,trans,isize)
  IMPLICIT NONE
  COMPLEX*16, intent(inout) :: A(isize,isize)
  COMPLEX*16, intent(in)    :: U(isize,isize)
  CHARACTER,  intent(in)    :: trans
  INTEGER,    intent(in)    :: isize
  !f2py integer intent(hide), depend(U)  :: isize=shape(U,0)
  COMPLEX*16 :: temp(isize,isize)  
  IF (trans.EQ.'N') THEN
     CALL Zproduct_MM(temp,A,U,'N','N',isize,isize,isize,isize,isize,isize)
     CALL Zproduct_MM(A,U,temp,'C','N',isize,isize,isize,isize,isize,isize)
  ELSEIF (trans.EQ.'C') THEN
     CALL Zproduct_MM(temp,A,U,'N','C',isize,isize,isize,isize,isize,isize)
     CALL Zproduct_MM(A,U,temp,'N','N',isize,isize,isize,isize,isize,isize)
  ELSE     
     print *, 'The call to ZTransform was wrong!'
     CALL EXIT(1)
  ENDIF
end subroutine ZTransform

subroutine ZTransformN(Anew,A,U,isize1,isize2)
  IMPLICIT NONE
  COMPLEX*16, intent(out)   :: Anew(isize2,isize2)
  COMPLEX*16, intent(in)    :: A(isize1,isize1)
  COMPLEX*16, intent(in)    :: U(isize1,isize2)
  INTEGER,    intent(in)    :: isize1, isize2
  !f2py integer intent(hide), depend(U)  :: isize2=shape(U,1)
  !f2py integer intent(hide), depend(U)  :: isize1=shape(U,0)
  COMPLEX*16 :: temp(isize1,isize2)
  CALL Zproduct_MM(temp,A,U,'N','N',   isize1,isize1,isize1,isize2,isize1,isize2)
  CALL Zproduct_MM(Anew,U,temp,'C','N',isize1,isize2,isize1,isize2,isize2,isize2)
end subroutine ZTransformN

subroutine ZTransformC(Anew, A, U, isize1, isize2)
  IMPLICIT NONE
  COMPLEX*16, intent(out)   :: Anew(isize2,isize2)
  COMPLEX*16, intent(in)    :: A(isize1,isize1)
  COMPLEX*16, intent(in)    :: U(isize2,isize1)
  INTEGER,    intent(in)    :: isize1,isize2
  !f2py integer intent(hide), depend(U)   :: isize1=shape(U,1)
  !f2py integer intent(hide), depend(U)   :: isize2=shape(U,0)
  COMPLEX*16 :: temp(isize1,isize2)
  CALL Zproduct_MM(temp,A,U,   'N','C',isize1,isize1,isize2,isize1,isize1,isize2)
  CALL Zproduct_MM(Anew,U,temp,'N','N',isize2,isize1,isize1,isize2,isize2,isize2)
end subroutine ZTransformC

subroutine ZProduct_NN(C, A, B, na1, na2, nb2)
  IMPLICIT NONE
  COMPLEX*16, intent(in)  :: A(na1,na2)
  COMPLEX*16, intent(in)  :: B(na2,nb2)
  COMPLEX*16, intent(out) :: C(na1,nb2)
  INTEGER,    intent(in)  :: na1, na2, nb2
  !f2py integer intent(hide), depend(A)  :: na1=shape(A,0)
  !f2py integer intent(hide), depend(A)  :: na2=shape(A,1)
  !f2py integer intent(hide), depend(B)  :: nb2=shape(B,1)  
  ! local variables
  COMPLEX*16 :: alpha_, beta_
  alpha_ = 1.
  beta_ = 0
  CALL zgemm('N','N',na1,nb2,na2,alpha_,A,na1,B,na2,beta_,C,na1)
end subroutine ZProduct_NN

subroutine ZProduct_NC(C, A, B, na1, na2, nb1)
  IMPLICIT NONE
  COMPLEX*16, intent(in)  :: A(na1,na2)
  COMPLEX*16, intent(in)  :: B(nb1,na2)
  COMPLEX*16, intent(out) :: C(na1,nb1)
  INTEGER,    intent(in)  :: na1, na2, nb1
  !f2py integer intent(hide), depend(A)  :: na1=shape(A,0)
  !f2py integer intent(hide), depend(A)  :: na2=shape(A,1)
  !f2py integer intent(hide), depend(B)  :: nb1=shape(B,0)  
  ! local variables
  COMPLEX*16 :: alpha_, beta_
  alpha_ = 1.
  beta_ = 0
  CALL zgemm('N','C',na1,nb1,na2,alpha_,A,na1,B,nb1,beta_,C,na1)
end subroutine ZProduct_NC

subroutine ZProduct_CN(C, A, B, na1, na2, nb2)
  IMPLICIT NONE
  COMPLEX*16, intent(in)  :: A(na1,na2)
  COMPLEX*16, intent(in)  :: B(na1,nb2)
  COMPLEX*16, intent(out) :: C(na2,nb2)
  INTEGER,    intent(in)  :: na1, na2, nb2
  !f2py integer intent(hide), depend(A)  :: na1=shape(A,0)
  !f2py integer intent(hide), depend(A)  :: na2=shape(A,1)
  !f2py integer intent(hide), depend(B)  :: nb2=shape(B,1)  
  ! local variables
  COMPLEX*16 :: alpha_, beta_
  alpha_ = 1.
  beta_ = 0
  CALL zgemm('C','N',na2,nb2,na1,alpha_,A,na1,B,na1,beta_,C,na2)
end subroutine ZProduct_CN

subroutine ZProduct_CC(C, A, B, na1, na2, nb1)
  IMPLICIT NONE
  COMPLEX*16, intent(in)  :: A(na1,na2)
  COMPLEX*16, intent(in)  :: B(nb1,na1)
  COMPLEX*16, intent(out) :: C(na2,nb1)
  INTEGER,    intent(in)  :: na1, na2, nb1
  !f2py integer intent(hide), depend(A)  :: na1=shape(A,0)
  !f2py integer intent(hide), depend(A)  :: na2=shape(A,1)
  !f2py integer intent(hide), depend(B)  :: nb1=shape(B,0)  
  ! local variables
  COMPLEX*16 :: alpha_, beta_
  alpha_ = 1.
  beta_ = 0
  CALL zgemm('C','C',na2,nb1,na1,alpha_,A,na1,B,nb1,beta_,C,na2)
end subroutine ZProduct_CC

subroutine seigval(w, vr, A, isize)
  IMPLICIT NONE
  COMPLEX*16, intent(in)  :: A(isize,isize)
  REAL*8, intent(out)     :: w(isize)
  COMPLEX*16, intent(out) :: vr(isize,isize)
  INTEGER, intent(in)     :: isize
  !f2py integer intent(hide), depend(A)  :: isize=shape(A,0)
  ! temporaries
  COMPLEX*16, allocatable :: work(:)
  REAL*8, allocatable     :: rwork(:)
  INTEGER, allocatable    :: iwork(:)
  INTEGER :: lwork, lrwork, liwork, info

  lwork  = 2*isize + isize*isize + 5
  lrwork = 5*isize + 2*isize*isize + 5
  liwork = 3 + 5*isize + 5

  ALLOCATE(work(lwork), rwork(lrwork), iwork(liwork))

  vr = A
  CALL ZHEEVD('V', 'U', isize, vr, isize, w, work, lwork, rwork, lrwork, iwork, liwork, info)
  
  IF (info .NE. 0) THEN
     WRITE(0,*) 'ERROR in ZHEEVD inside seigval'
  ENDIF

  DEALLOCATE(work, rwork, iwork)
end subroutine seigval

subroutine zeigval(zek, A, isize)
  IMPLICIT NONE
  COMPLEX*16, intent(in)  :: A(isize,isize)
  COMPLEX*16, intent(out) :: zek(isize)
  INTEGER, intent(in)     :: isize
  !f2py integer intent(hide), depend(A)  :: isize=shape(A,0)
  ! temporaries
  COMPLEX*16, allocatable :: work(:), At(:,:)
  REAL*8, allocatable     :: rwork(:)
  INTEGER                 :: lwork, lrwork, info
  REAL*8, PARAMETER       :: smalleps = 1e-5
  COMPLEX*16, allocatable :: evl(:,:), evr(:,:)

  lwork  = 3*isize + 1
  lrwork = 2*isize
  
  ALLOCATE(work(lwork), rwork(lrwork), At(isize,isize))

  At = A
  CALL zgeev('N','N',isize,At,isize,zek,evl,isize,evr,isize,work,lwork,rwork,info)

  IF (info .NE. 0) THEN
     WRITE(0,*) 'ERROR in ZGEEV inside seigval'
  ENDIF
  DEALLOCATE(At, work, rwork)
END subroutine zeigval

subroutine zeigsys(zek, evl, evr, A, isize)
  IMPLICIT NONE
  COMPLEX*16, intent(in)  :: A(isize,isize)
  COMPLEX*16, intent(out) :: zek(isize)
  COMPLEX*16, intent(out) :: evl(isize,isize)
  COMPLEX*16, intent(out) :: evr(isize,isize)
  INTEGER, intent(in)     :: isize
  !f2py integer intent(hide), depend(A)  :: isize=shape(A,0)
  ! temporaries
  COMPLEX*16, allocatable :: work(:), At(:,:)
  REAL*8, allocatable     :: rwork(:)
  INTEGER    :: lwork, lrwork, info
  REAL*8, PARAMETER       :: smalleps = 1e-5
  COMPLEX*16, allocatable :: scaler(:)
  INTEGER    :: q, p
  COMPLEX*16 :: ctmp
  lwork  = 3*isize + 1
  lrwork = 2*isize
  
  ALLOCATE(work(lwork), rwork(lrwork))
  ALLOCATE(At(isize,isize))
  ALLOCATE(scaler(isize))

  At = A
  CALL zgeev('V','V',isize,At,isize,zek,evl,isize,evr,isize,work,lwork,rwork,info)

  IF (info .NE. 0) THEN
     WRITE(0,*) 'ERROR in ZGEEV inside seigval'
  ENDIF
  ! transpose left eigenvectors
  evl = dconjg(TRANSPOSE(evl))

  ! Maybe this is not really necessary casue zgeev gives already orthogonal eigensystem
  !========== Step 5, Make degenerate eigenvectors orthogonal
  DO q=1,isize
     DO p=1,q-1
        IF (abs(zek(p)-zek(q)).LT.smalleps .AND. abs(scalprod(evl(p,:),evr(:,q),isize)) .GT.smalleps) THEN
           evr(:,q) = evr(:,q) - scalprod(evl(p,:),evr(:,q),isize)/scalprod(evl(p,:),evr(:,p),isize) * evr(:,p)
        ENDIF
     ENDDO
     DO p=1,q-1
        IF (abs(zek(p)-zek(q)).LT.smalleps .AND. abs(scalprod(evl(q,:),evr(:,p),isize)) .GT.smalleps) THEN
           evl(q,:) = evl(q,:) - scalprod(evl(q,:),evr(:,p),isize)/scalprod(evl(p,:),evr(:,p),isize) * evl(p,:)
        ENDIF
     ENDDO
  ENDDO
  !========= Step 6, Normalize eigenvectors
  DO p = 1,isize
     ctmp = 0.d0
     DO q = 1,isize
        ctmp = ctmp+evl(p,q)*evr(q,p)
     ENDDO
     scaler(p) = SQRT(ctmp)
  ENDDO
  DO p = 1,isize
     evl(p,:) = evl(p,:)/scaler(p)
     evr(:,p) = evr(:,p)/scaler(p)
  ENDDO

  !========== Deallocate dynamic arrays ==========
  DEALLOCATE(scaler)
  DEALLOCATE(At)
  DEALLOCATE(work, rwork)

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

end subroutine zeigsys

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


subroutine ZProduct_MD(C, A, D, na1, na2)
  IMPLICIT NONE
  COMPLEX*16, intent(out) :: C(na1,na2)
  COMPLEX*16, intent(in)  :: A(na1,na2)
  COMPLEX*16, intent(in)  :: D(na2)
  INTEGER,    intent(in)  :: na1, na2
  !f2py integer intent(hide), depend(A)  :: na1=shape(A,0)
  !f2py integer intent(hide), depend(A)  :: na2=shape(A,1)
  ! local variables
  INTEGER :: i, j
  DO i=1,na1
     DO j=1,na2
        C(j,i) = A(j,i)*D(i)
     ENDDO
  ENDDO
end subroutine ZProduct_MD

subroutine check_causality(E, small_positive, isize)
  IMPLICIT NONE
  COMPLEX*16, intent(inout) :: E(isize)
  REAL*8, intent(in)        :: small_positive
  INTEGER, intent(in)       :: isize  
  !f2py integer intent(hide), depend(E)  :: isize=shape(E,0)
  INTEGER :: i
  REAL*8  :: a
  DO i=1,isize
     IF (dimag(E(i)).GT.-small_positive) E(i) = dcmplx(dreal(E(i)),-small_positive)
  ENDDO
end subroutine check_causality

subroutine check_causality2(sig, small_positive, isize)
  IMPLICIT NONE
  COMPLEX*16, intent(inout) :: sig(isize,isize)
  REAL*8, intent(in)        :: small_positive
  INTEGER, intent(in)       :: isize  
  !f2py integer intent(hide), depend(sig)  :: isize=shape(sig,0)
  INTEGER :: i, j
  REAL*8  :: a
  DO i=1,isize
     IF (dimag(sig(i,i)).GT.-small_positive) sig(i,i) = dcmplx(dreal(sig(i,i)),-small_positive)
  ENDDO
end subroutine check_causality2

subroutine ZProduct_ADA(C, A, D, B, na1, na2, nb2)
  IMPLICIT NONE
  COMPLEX*16, intent(in)  :: A(na1,na2)
  COMPLEX*16, intent(in)  :: D(na2)
  COMPLEX*16, intent(in)  :: B(na2,nb2)
  COMPLEX*16, intent(out) :: C(na1,nb2)
  INTEGER,    intent(in)  :: na1, na2, nb2
  !f2py integer intent(hide), depend(A)  :: na1=shape(A,0)
  !f2py integer intent(hide), depend(A)  :: na2=shape(A,1)
  !f2py integer intent(hide), depend(B)  :: nb2=shape(B,1)  
  ! local variables
  INTEGER    :: i, j
  COMPLEX*16 :: tmp(na1,na2)
  
  DO j=1,na2
     DO i=1,na1
        tmp(i,j) = A(i,j)*D(j)
     ENDDO
  ENDDO

  CALL ZProduct_NN(C, tmp, B, na1, na2, nb2)

end subroutine ZProduct_ADA

subroutine ZProduct_ADAt(C, A, D, B, na1, na2, nb2, temp)
  IMPLICIT NONE
  COMPLEX*16, intent(in)  :: A(na1,na2)
  COMPLEX*16, intent(in)  :: D(na2)
  COMPLEX*16, intent(in)  :: B(na2,nb2)
  COMPLEX*16, intent(out) :: C(na1,nb2)
  INTEGER,    intent(in)  :: na1, na2, nb2
  COMPLEX*16, intent(inout) :: temp(na1,na2)
  !f2py integer intent(hide), depend(A)  :: na1=shape(A,0)
  !f2py integer intent(hide), depend(A)  :: na2=shape(A,1)
  !f2py integer intent(hide), depend(B)  :: nb2=shape(B,1)  
  ! local variables
  INTEGER    :: i, j
  COMPLEX*16 :: tmp(na1,na2)
  
  DO j=1,na2
     DO i=1,na1
        tmp(i,j) = A(i,j)*D(j)
     ENDDO
  ENDDO

  CALL ZProduct_NN(C, tmp, B, na1, na2, nb2)

end subroutine ZProduct_ADAt
