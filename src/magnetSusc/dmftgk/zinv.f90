subroutine zinv(A, ndim)
  IMPLICIT NONE
  COMPLEX*16, intent(inout) :: A(ndim,ndim)
  INTEGER, intent(in)       :: ndim
  ! locals
  INTEGER    :: info, lwork, lda
  INTEGER    :: ipiv(ndim)
  COMPLEX*16 :: work(ndim*64)
  lwork = ndim*64
  lda = ndim
  
  CALL ZGETRF( ndim, ndim, A, lda, ipiv, info )

  if (info.ne.0) then 
     print *, 'zgetrf info=', info
  endif

  CALL ZGETRI( ndim, A, lda, ipiv, work, lwork, info )
  
  if (info.ne.0) then
     print *, 'zgetri info=', info
  endif
  
end subroutine zinv

subroutine dinv(A, ndim)
  IMPLICIT NONE
  REAL*8, intent(inout) :: A(ndim,ndim)
  INTEGER, intent(in)   :: ndim
  ! locals
  INTEGER    :: info, lwork, lda
  INTEGER    :: ipiv(ndim)
  REAL*8 :: work(ndim*64)
  lwork = ndim*64
  lda = ndim
  
  CALL DGETRF( ndim, ndim, A, lda, ipiv, info )

  if (info.ne.0) then 
     print *, 'dgetrf info=', info
  endif

  CALL DGETRI( ndim, A, lda, ipiv, work, lwork, info )
  
  if (info.ne.0) then
     print *, 'dgetri info=', info
  endif
  
end subroutine dinv
