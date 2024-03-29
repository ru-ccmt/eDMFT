! @Copyright 2007 Kristjan Haule
! 

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
  
  i0f = i0+1
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
  
  RETURN

END SUBROUTINE kramarskronig
