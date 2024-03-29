subroutine vsincos(x,y,z,n)
  implicit none
  real*8 x(*),y(*),z(*)
  integer n, j
  !$OMP PARALLEL DO 
  do j=1,n
     x(j)=sin(z(j))
     y(j)=cos(z(j))
  enddo
  !$OMP END PARALLEL DO
  return
end subroutine vsincos
subroutine vcos(y,x,n)
  implicit none
  real*8 x(*),y(*)
  integer n,j
  !$OMP PARALLEL DO 
  do j=1,n
     y(j)=cos(x(j))
  enddo
  !$OMP END PARALLEL DO
  return
end subroutine vcos
subroutine vcosisin(y,x,n)
  implicit none
  complex*16 y(*)
  real*8 x(*)
  integer n,j
  !$OMP PARALLEL DO 
  do j=1,n
     y(j)=dcmplx(cos(x(j)),sin(x(j)))
  enddo
  !$OMP END PARALLEL DO
  return
end subroutine vcosisin
