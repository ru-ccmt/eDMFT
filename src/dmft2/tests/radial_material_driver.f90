program radial_material_driver
  use iso_fortran_env, only: real64
  use outwin_response
  implicit none
  integer, parameter :: dp=real64
  type(radial_pair_context) :: center,plus,minus
  real(dp), allocatable :: v(:),dv(:),bv(:),direction(:,:,:),covector(:,:,:),fd(:,:,:)
  real(dp) :: r0,dx,energy,ell,z,de,be,r,step,error,left,right,rev,t0,t1,t2
  integer :: n,mode,i,j,k,iteration
  character(1024) :: path
  call get_command_argument(1,path)
  open(10,file=trim(path),status='old')
  read(10,*) n,mode,r0,dx,energy,ell,z
  allocate(v(n),dv(n),bv(n),direction(n,4,2),covector(n,4,2),fd(n,4,2))
  read(10,*) v
  close(10)
  do i=1,n
     r=r0*exp(dx*(i-1))
     dv(i)=0.07_dp*r*sin(1.3_dp*r)
     do j=1,4
        do k=1,2
           covector(i,j,k)=sin(.2_dp*i+j+.4_dp*k)*exp(-r)
        enddo
     enddo
  enddo
  de=.03_dp;step=.002_dp
  call prepare_radial_pair(center,mode==1,v,r0,dx,energy,ell,z)
  call radial_pair_action(center,dv,de,direction)
  call radial_pair_pullback(center,covector,bv,be)
  left=sum(covector*direction);right=sum(bv*dv)+be*de
  rev=abs(left-right)/max(1.0_dp,abs(left),abs(right))
  call prepare_radial_pair(plus,mode==1,v+step*dv,r0,dx,energy+step*de,ell,z)
  call prepare_radial_pair(minus,mode==1,v-step*dv,r0,dx,energy-step*de,ell,z)
  fd=(radial_pair_value(plus)-radial_pair_value(minus))/(2*step)
  error=sqrt(sum((fd-direction)**2)/sum(direction**2))
  call cpu_time(t0)
  do iteration=1,1000
     call radial_pair_action(center,dv,de,direction)
  enddo
  call cpu_time(t1)
  do iteration=1,1000
     call radial_pair_pullback(center,covector,bv,be)
  enddo
  call cpu_time(t2)
  write(*,'(4ES24.15)') error,rev,(t1-t0)/1000,(t2-t1)/1000
  if(error>2e-5_dp.or.rev>2e-10_dp) error stop 'material radial validation failed'
end program
