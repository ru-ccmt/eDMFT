subroutine native_normalization(nr,nat,nkk,r,rad,npts,mult,rmt,pos,cell,kzz,pw,target,inversion,charges) bind(c)
use iso_c_binding
use struct,only:br1,vi
implicit none
integer(c_int),value:: nr,nat,nkk,inversion
integer(c_int):: npts(nat),mult(nat),kzz(3,nkk)
real(c_double):: r(nr,nat),rad(nr,nat),rmt(nat),pos(3,nat*48),cell(3),target,charges(nat+1)
complex(c_double_complex)::pw(nkk)
real(8)::vols(nat),qel(nat,2),qpw(2),zz(nat),delta,pi
complex(8)::mask(nkk)
real(8),external::zzdotk
integer i
logical err
pi=acos(-1d0);vi=product(cell);br1=0
open(21,status='scratch')
do i=1,3
 br1(i,i)=2*pi/cell(i)
enddo
vols=4*pi*rmt**3/3
call setn(nkk,mask,kzz,nat,mult,rmt,vols,pos,inversion/=0)
qpw=0;qel=0;zz=0;zz(1)=target/mult(1)
qpw(1)=zzdotk(nkk,mask,1,pw,1)
do i=1,nat
 delta=log(r(npts(i),i)/r(1,i))/(npts(i)-1)
 call charge(r(:,i),delta,rad(:,i),1,npts(i),qel(i,1))
enddo
charges(1:nat)=qel(:,1);charges(nat+1)=qpw(1)
err=.false.
! IsVresp suppresses leakage warnings in synthetic/one-k controls. NORMAL's
! rescaling itself is unchanged with NORM=YES, as in ordinary output assembly.
call normal(pw,nkk,nkk,1,nat,mult,zz,qpw,qel,'01','YES  ',0d0,.true.,err)
close(21)
end subroutine

subroutine zdscal(n,alpha,x,inc)
implicit none
integer n,inc,i
real(8) alpha
complex(8) x(*)
do i=1,n*inc,inc
 x(i)=x(i)*alpha
enddo
end subroutine
