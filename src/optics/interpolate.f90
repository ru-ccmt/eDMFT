! @Copyright 2007 Kristjan Haule

SUBROUTINE findNext(x, ix, R, n)
  IMPLICIT NONE
  INTEGER, intent(inout) :: ix
  REAL*8, intent(in)  :: R(n)
  REAL*8, intent(in)  :: x
  INTEGER, intent(in) :: n
  !f2py integer intent(hide), depend(R)      :: n = shape(R,0)
  INTEGER :: i

  if (x.le.R(1)) then
     ix=1
     return 
  endif
  if (x.gt.R(n-1)) then
     ix=n-1
     return
  endif
  if (ix.lt.1) then
     ix=1
  endif
  if (ix.gt.n) then
     ix=n
  endif
  if (x.ge.R(ix)) then
     do i=ix,n
        if (R(i).ge.x) then
           ix = i-1
           return
        endif
     enddo
     ix = n-1
     return
  else
     do i=ix,1,-1
        if (x.gt.R(i)) then
           ix = i
           return
        endif
     enddo
     ix=1
     return
  endif
end SUBROUTINE findNext

subroutine interp(res, f, x, ix, R, n)
  IMPLICIT NONE
  COMPLEX*16, intent(out) :: res
  COMPLEX*16, intent(in)  :: f(n)
  REAL*8, intent(in)  :: x
  INTEGER, intent(in) :: ix
  REAL*8, intent(in)  :: R(n)
  INTEGER, intent(in) :: n
  !f2py integer intent(hide), depend(R)      :: n = shape(R,0)
  if (x.lt.R(ix) .or. x.gt.R(ix+1)) then
     print *, 'Tezave:', x, R(ix), R(ix+1), ix, n
  endif
  res = f(ix) + (f(ix+1)-f(ix))*(x-R(ix))/(R(ix+1)-R(ix))
  return
end subroutine interp
