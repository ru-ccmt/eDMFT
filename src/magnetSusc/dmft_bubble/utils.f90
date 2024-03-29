! @Copyright 2007 Kristjan Haule
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


SUBROUTINE kramarskronig(gr, fi, om, i0, nom)
  IMPLICIT NONE
  REAL*8, intent(out):: gr
  REAL*8, intent(in) :: fi(nom)
  REAL*8, intent(in) :: om(nom)
  INTEGER, intent(in) :: i0, nom
  !f2py integer intent(hide), depend(om)  :: nom = shape(om,0)
  !locals
  REAL*8, PARAMETER :: pi=3.14159265
  REAL*8 :: sm
  REAL*8 :: dh(nom)
  REAL*8 :: S0, x
  INTEGER :: j, i0f
  i0f = i0
  S0 = fi(i0f)
  x = om(i0f)
  dh(1) = 0.5*(om(2)-om(1))
  do j=2,nom-1
     dh(j) = 0.5*(om(j+1)-om(j-1))
  enddo
  dh(nom) = 0.5*(om(nom)-om(nom-1))
  sm=0
  do j=1,i0f-2
     sm = sm + (fi(j)-S0)*dh(j)/(om(j)-x)
  enddo
  if (i0f>1)   sm = sm + (fi(i0f-1)-S0)*(dh(i0f-1)+0.5*dh(i0f))/(om(i0f-1)-x)
  if (i0f<nom) sm = sm + (fi(i0f+1)-S0)*(dh(i0f+1)+0.5*dh(i0f))/(om(i0f+1)-x)
  do j=i0f+2,nom
     sm = sm + (fi(j)-S0)*dh(j)/(om(j)-x)
  enddo
  if (.NOT.(i0f .EQ. 1 .OR. i0f.EQ.nom)) then
     sm = sm + S0*log(abs((om(nom)-x)/(x-om(1))))
  endif
  gr = sm/pi
END SUBROUTINE kramarskronig

SUBROUTINE complex_kramarskronig(gr, fi, om, i0, nom)
  IMPLICIT NONE
  COMPLEX*16, intent(out):: gr
  COMPLEX*16, intent(in) :: fi(nom)
  REAL*8, intent(in) :: om(nom)
  INTEGER, intent(in) :: i0, nom
  !f2py integer intent(hide), depend(om)  :: nom = shape(om,0)
  !locals
  REAL*8, PARAMETER :: pi=3.14159265
  REAL*8 :: x
  REAL*8 :: dh(nom)
  COMPLEX*16 :: S0, sm
  INTEGER :: j, i0f
  i0f = i0
  S0 = fi(i0f)
  x = om(i0f)
  dh(1) = 0.5*(om(2)-om(1))
  do j=2,nom-1
     dh(j) = 0.5*(om(j+1)-om(j-1))
  enddo
  dh(nom) = 0.5*(om(nom)-om(nom-1))
  sm=0
  do j=1,i0f-2
     sm = sm + (fi(j)-S0)*dh(j)/(om(j)-x)
  enddo
  if (i0f>1)   sm = sm + (fi(i0f-1)-S0)*(dh(i0f-1)+0.5*dh(i0f))/(om(i0f-1)-x)
  if (i0f<nom) sm = sm + (fi(i0f+1)-S0)*(dh(i0f+1)+0.5*dh(i0f))/(om(i0f+1)-x)
  do j=i0f+2,nom
     sm = sm + (fi(j)-S0)*dh(j)/(om(j)-x)
  enddo
  if (.NOT.(i0f .EQ. 1 .OR. i0f.EQ.nom)) then
     sm = sm + S0*log(abs((om(nom)-x)/(x-om(1))))
  endif
  gr = sm/pi
END SUBROUTINE complex_kramarskronig

SUBROUTINE integrate_trapz(res, f, x, n)
  IMPLICIT NONE
  REAL*8,  intent(out) :: res
  REAL*8,  intent(in)  :: f(n), x(n)
  INTEGER, intent(in)  :: n
  !
  INTEGER :: i
  res = f(1)*0.5*(x(2)-x(1)) + f(n)*0.5*(x(n)-x(n-1))
  do i=2,n-1
     res = res + f(i)*0.5*(x(i+1)-x(i-1))
  enddo
END SUBROUTINE integrate_trapz

SUBROUTINE complex_integrate_trapz(res, f, x, n)
  IMPLICIT NONE
  COMPLEX*16,intent(out) :: res
  COMPLEX*16,intent(in)  :: f(n)
  REAL*8, intent(in)     :: x(n)
  INTEGER, intent(in)    :: n
  !
  INTEGER :: i
  res = f(1)*0.5*(x(2)-x(1)) + f(n)*0.5*(x(n)-x(n-1))
  do i=2,n-1
     res = res + f(i)*0.5*(x(i+1)-x(i-1))
  enddo
END SUBROUTINE complex_integrate_trapz
