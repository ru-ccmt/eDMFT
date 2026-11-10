program test_lo_response
  use iso_fortran_env, only: real64
  use outwin_response
  implicit none
  integer, parameter :: dp=real64
  real(dp) :: x(10),xp(10),xm(10),c(3),reference(3),plus(3),minus(3),jac(3,10),fd(3,10),h,error
  integer :: mode,jlo,i
  logical :: lapw
  do mode=1,4
     lapw=mode==1.or.mode==4
     jlo=1
     if(mode==3) jlo=2
     x=[.6_dp,.3_dp,-.4_dp,.8_dp,1.4_dp,.5_dp,-.6_dp,.2_dp,.1_dp,2.0_dp]
     if(mode==4) x=[1.0_dp,0.0_dp,0.0_dp,-1.0_dp,1.0_dp,1.0_dp,0.0_dp, &
                   1.0_dp-1e-6_dp,0.0_dp,1.0_dp]
     call local_orbital_coefficients(lapw,jlo,x,c,jac)
     call native(x,lapw,jlo,reference)
     if(maxval(abs(c-reference))>1e-12_dp) error stop 'ABC central mismatch'
     h=1e-7_dp
     do i=1,10
        xp=x;xm=x;xp(i)=xp(i)+h;xm(i)=xm(i)-h
        call native(xp,lapw,jlo,plus)
        call native(xm,lapw,jlo,minus)
        fd(:,i)=(plus-minus)/(2*h)
     enddo
     error=maxval(abs(fd-jac))/max(1.0_dp,maxval(abs(jac)))
     write(*,'(A,I2,A,ES15.7)') 'ABC branch ',mode,' FD Jacobian relative=',error
     if(error>1e-7_dp) error stop 'ABC Jacobian mismatch'
  enddo
  print *, 'PASS ABC tangent and transpose matrix'
contains
  subroutine native(x,lapw,jlo,c)
    real(dp),intent(in) :: x(10)
    logical,intent(in) :: lapw
    integer,intent(in) :: jlo
    real(dp),intent(out) :: c(3)
    call abc0(c(1),c(2),c(3),0,jlo,lapw,x(1),x(2),x(3),x(4),x(5), &
              x(6),x(7),x(8),x(9),1,1,1,1,x(10),.false.)
  end subroutine
end program
