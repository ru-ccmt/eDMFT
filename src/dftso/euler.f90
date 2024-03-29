Subroutine Euler(Rot_new,a,b,c)
  IMPLICIT NONE
  REAL*8, intent(in) :: Rot_new(3,3)
  REAL*8, intent(out):: a,b,c
  ! external
  REAL*8 :: ang
  ! locals
  REAL*8 :: Rot_old(3,3), z(3), zz(3)
  REAL*8 :: y(3),yy(3),yyy(3),pom(3)
  REAL*8 :: pi, y_norm
  INTEGER :: i, j
  !
  pi=acos(-1d0)
  do i=1,3
     do j=1,3
        Rot_old(i,j)=0
        if (i.eq.j) Rot_old(i,i)=1
     end do
  end do

  do j=1,3
     y(j)=Rot_old(j,2)
     yyy(j)=Rot_new(j,2)
     z(j)=Rot_old(j,3)
     zz(j)=Rot_new(j,3)
  end do

  call vecprod(z,zz,yy)
  y_norm=dsqrt(dot_product(yy,yy))

  if (y_norm.lt.1d-10) then
     a=ang(y,yyy)
     if (dot_product(z,zz).gt.0.d0) then
        c=0.d0
        b=0.d0
        if (yyy(1).gt.0.d0) a=2*pi-a
     else
        c=a
        a=0.d0
        b=pi
        if (yyy(1).lt.0.d0) c=2*pi-c
     end if
  else
     do j=1,3
        yy(j)=yy(j)/y_norm
     end do
     a=ang(y,yy)
     b=ang(z,zz)
     c=ang(yy,yyy)
     if (yy(1).gt.0.d0) a=2*pi-a
     call vecprod(yy,yyy,pom)
     if (dot_product(pom,zz).lt.0.d0) c=2*pi-c
  end if
END Subroutine Euler


Subroutine vecprod(a,b,c)
  IMPLICIT NONE
  REAL*8, intent(out) :: c(3)
  REAL*8, intent(in)  :: a(3), b(3)
  c(1)=a(2)*b(3)-a(3)*b(2)
  c(2)=a(3)*b(1)-a(1)*b(3)
  c(3)=a(1)*b(2)-a(2)*b(1)
end Subroutine vecprod

REAL*8 function ang(a,b)
  IMPLICIT NONE
  REAL*8, intent(in) :: a(3), b(3)
  ! locals
  REAL*8 :: c(3), pi, aa, bb, e, cc
  pi=dacos(-1.d0)
  aa=dsqrt(dot_product(a,a))
  bb=dsqrt(dot_product(b,b))
  e=dot_product(a,b)/(aa*bb)
  if (abs(e).gt.0.8) then
     call vecprod(a,b,c)
     cc=dsqrt(dot_product(c,c))
     ang=dasin(cc/(aa*bb))
     if (e.lt.0.d0) ang=pi-ang
  else
     ang=dacos(e)
  end if
end function ang
