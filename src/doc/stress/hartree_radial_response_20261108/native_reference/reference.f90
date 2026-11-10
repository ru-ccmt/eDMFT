subroutine radial_dirichlet_ref(n,r,dx,l,c,bdy,out)
  implicit none
  integer,intent(in)::n,l
  real(8),intent(in)::r(n),dx,c(n),bdy
  real(8),intent(out)::out(n)
  real(8)::aa(n),bb(n),fa(n),fb(n),fac,pi
  integer::i
  pi=acos(-1.d0)
  do i=1,n
    aa(i)=r(i)**l
    bb(i)=(1.d0-(r(i)/r(n))**(2*l+1))/r(i)**(l+1)
  enddo
  fa=aa*c
  fb=bb*c
  call CoulInt(r,dx,n,fa)
  call CoulInt(r,dx,n,fb)
  fac=8.d0*pi/(2*l+1)
  if(l==0) fac=fac/(4.d0*pi)
  do i=1,n
    out(i)=fac*(bb(i)*fa(i)+aa(i)*(fb(n)-fb(i)))+(r(i)/r(n))**l*bdy
  enddo
  if(l==0) out=out*r
end subroutine
