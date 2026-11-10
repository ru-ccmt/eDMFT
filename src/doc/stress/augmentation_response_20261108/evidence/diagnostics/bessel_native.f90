module parallel
contains
subroutine abort_parallel
 error stop 'native Bessel error'
end subroutine
end module
subroutine outerr(a,b)
 character(*)::a,b
 error stop 'native Bessel outerr'
end subroutine
subroutine native_bessel(n,lmax,x,j,dj) bind(C)
 use iso_c_binding
 implicit none
 integer(c_int),value::n,lmax
 real(c_double),intent(in)::x(n)
 real(c_double),intent(out)::j(lmax+1,n),dj(lmax+1,n)
 integer::i
 do i=1,n
  call sphbes(lmax,x(i),j(:,i))
  call dvbes1(j(:,i),dj(:,i),x(i),1.0_c_double,lmax+1)
 enddo
end subroutine
