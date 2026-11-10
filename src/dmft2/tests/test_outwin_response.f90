program test_outwin_response
  use iso_fortran_env, only: real64
  use outwin_response
  implicit none
  integer, parameter :: dp=real64, maximum=401
  type(radial_pair_context) :: center, plus, minus
  real(dp) :: potential(maximum),dv(maximum),bar_v(maximum),r(maximum)
  real(dp) :: tangent(maximum,4,2),covector(maximum,4,2),fd(maximum,4,2)
  real(dp) :: values(maximum,4,2),step,dx,r0,de,bar_e,left,right,err,rev,maxfd,maxrev
  real(dp) :: energy,ell,charge,first_error,last_error,zero(maximum,4,2),orthogonality
  integer :: n,mode,l,i,j,k,ih
  logical :: rel
  maxfd=0;maxrev=0
  r0=1.0e-4_dp;energy=0.17_dp;charge=26;de=0.03_dp
  do n=400,401
     dx=log(2.1_dp/r0)/(n-1)
     do i=1,n
        r(i)=r0*exp(dx*(i-1))
        potential(i)=-charge*exp(-r(i))+0.15_dp*r(i)**2
        dv(i)=0.07_dp*r(i)*sin(1.3_dp*r(i))
        do j=1,4
           do k=1,2
              covector(i,j,k)=sin(0.2_dp*i+j+0.4_dp*k)*exp(-r(i))
           enddo
        enddo
     enddo
     do mode=0,1
        rel=mode==1
        do l=0,3
           ell=real(l,dp)
           call prepare_radial_pair(center,rel,potential(:n),r0,dx,energy,ell,charge)
           values(:n,:,:)=radial_pair_value(center)
           orthogonality=abs(sum(center%center%metric*values(:n,:,1)*values(:n,:,2)))
           if (orthogonality>1e-12_dp) error stop 'udot orthogonality failure'
           call radial_pair_action(center,dv(:n),de,tangent(:n,:,:))
           call radial_pair_pullback(center,covector(:n,:,:),bar_v(:n),bar_e)
           left=sum(covector(:n,:,:)*tangent(:n,:,:))
           right=sum(bar_v(:n)*dv(:n))+bar_e*de
           rev=abs(left-right)/max(1.0_dp,abs(left),abs(right))
           maxrev=max(maxrev,rev)
           do ih=1,2
              step=0.002_dp/ih
              call prepare_radial_pair(plus,rel,potential(:n)+step*dv(:n),r0,dx,energy+step*de,ell,charge)
              call prepare_radial_pair(minus,rel,potential(:n)-step*dv(:n),r0,dx,energy-step*de,ell,charge)
              fd(:n,:,:)=(radial_pair_value(plus)-radial_pair_value(minus))/(2*step)
              err=sqrt(sum((fd(:n,:,:)-tangent(:n,:,:))**2)/sum(tangent(:n,:,:)**2))
              if(ih==1) first_error=err
              last_error=err
           enddo
           maxfd=max(maxfd,last_error)
           call radial_pair_action(center,0*dv(:n),0.0_dp,zero(:n,:,:))
           if(maxval(abs(zero(:n,:,:)))>0) error stop 'nonzero zero-direction response'
           write(*,'(A,I4,A,L1,A,I2,A,ES12.4,A,ES12.4,A,F8.3)') &
                'n=',n,' RELA=',rel,' l=',l,' FD=',last_error,' adjoint=',rev,' step ratio=',first_error/last_error
        enddo
     enddo
  enddo
  if(maxfd>1e-5_dp.or.maxrev>1e-10_dp) error stop 'radial tangent validation failed'
  write(*,'(A,ES15.7,A,ES15.7)') 'PASS max FD relative=',maxfd,' max transpose=',maxrev
end program
