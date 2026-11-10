module radial
implicit none
real(8) :: A(881),B(881),AP(881),BP(881)
end module
module TestCases
implicit none
logical :: StartWin=.false., WinFour=.false., SimpInt=.false., ConOpt=.true.
end module
subroutine native_outwin(n,rel,vr,radius,dx,e,l,z,result) bind(c)
use iso_c_binding
use radial
implicit none
integer(c_int), value :: n,rel,l
real(c_double),intent(in) :: vr(n),radius(n)
real(c_double),value :: dx,e,z
real(c_double),intent(out) :: result(n,4)
real(8) :: vpad(881),rpad(881),val,slo
integer :: nodes
vpad=0;rpad=0;vpad(:n)=vr;rpad(:n)=radius
call OUTWINB(rel/=0,vpad,rpad,dx,n,e,real(l,8),val,slo,nodes,z)
result(:,1)=A(:n);result(:,2)=B(:n);result(:,3)=AP(:n);result(:,4)=BP(:n)
end subroutine
