subroutine native_charge(n,r,dx,f,s) bind(c)
use iso_c_binding
implicit none
integer(c_int),value :: n
real(c_double),intent(in) :: r(n),f(n)
real(c_double),value :: dx
real(c_double),intent(out) :: s
call charge(r,dx,f,1,n,s)
end subroutine
