program test_radial_potential
 use iso_fortran_env, only: real64
 use lapw_radial_response
 use lapw_radial_potential
 implicit none
 integer,parameter::dp=real64,nv=3
 type inputs
  integer::n,mode,ell,lapw,nlo
  real(dp)::r0,dx,rmt,e,z
  real(dp),allocatable::v(:),dv(:),elo(:),delo(:)
 end type
 type(inputs)::x,y
 type(radial_basis_context)::a,b,ap,bp,am,bm
 real(dp),allocatable::v(:,:),vd(:,:),val(:,:,:),ref(:,:,:),vdot(:,:,:),vp(:,:,:),vm(:,:,:),bar(:,:,:)
 real(dp),allocatable::ad(:,:,:),bd(:,:,:),ab(:,:,:),bb(:,:,:),vb(:,:),ca(:,:),cb(:,:),av(:),bv(:),ae(:),be(:)
 real(dp)::left,right,ce,fd,adj,ae0,be0,h,r
 real(dp),external::RINT13
 integer::i,j,k,q,it
 character(1024)::path
 call get_command_argument(1,path);call read_inputs(path,x)
 call get_command_argument(2,path);call read_inputs(path,y)
 call prepare(x,0._dp,a);call prepare(y,0._dp,b)
 allocate(v(x%n,nv),vd(x%n,nv),val(2+x%nlo,2+y%nlo,nv),ref(2+x%nlo,2+y%nlo,nv))
 allocate(vdot(2+x%nlo,2+y%nlo,nv),vp(2+x%nlo,2+y%nlo,nv),vm(2+x%nlo,2+y%nlo,nv),bar(2+x%nlo,2+y%nlo,nv))
 allocate(ad(x%n,4,2+x%nlo),bd(y%n,4,2+y%nlo),ab(x%n,4,2+x%nlo),bb(y%n,4,2+y%nlo),vb(x%n,nv))
 allocate(ca(2+x%nlo,x%nlo),cb(2+y%nlo,y%nlo),av(x%n),bv(y%n),ae(x%nlo),be(y%nlo))
 do q=1,nv
  do i=1,x%n
   r=x%r0*exp(x%dx*(i-1));v(i,q)=sin(q*r)*exp(-.3_dp*r);vd(i,q)=.07_dp*cos(q*r+.2_dp)*r
  enddo
  do j=1,2+y%nlo
   do i=1,2+x%nlo
    bar(i,j,q)=sin(real(i+3*j+7*q,dp))
   enddo
  enddo
 enddo
 call radial_components_action(a,reshape([real(dp)::],[2,0]),x%dv,.03_dp,x%delo,ad,ca)
 call radial_components_action(b,reshape([real(dp)::],[2,0]),y%dv,.03_dp,y%delo,bd,cb)
 call primitive_potential_value(a,b,v,val)
 call oracle(a,b,v,ref)
 ce=sqrt(sum((val-ref)**2)/max(1e-30_dp,sum(ref**2)))
 call primitive_potential_action(a,b,v,ad,bd,vd,vdot)
 call primitive_potential_pullback(a,b,v,bar,ab,bb,vb)
 call radial_components_pullback(a,reshape([real(dp)::],[2,0]),ab,0*ca,av,ae0,ae)
 call radial_components_pullback(b,reshape([real(dp)::],[2,0]),bb,0*cb,bv,be0,be)
 left=sum(bar*vdot)
 right=sum(av*x%dv)+.03_dp*ae0+sum(ae*x%delo)+sum(bv*y%dv)+.03_dp*be0+sum(be*y%delo)+sum(vb*vd)
 adj=abs(left-right)/max(1._dp,abs(left),abs(right));fd=huge(1._dp)
 do it=1,3
  h=.004_dp/2**it
  call prepare(x,h,ap);call prepare(y,h,bp);call prepare(x,-h,am);call prepare(y,-h,bm)
  call oracle(ap,bp,v+h*vd,vp);call oracle(am,bm,v-h*vd,vm)
  fd=min(fd,sqrt(sum(((vp-vm)/(2*h)-vdot)**2)/max(1e-30_dp,sum(vdot**2))))
 enddo
 write(*,'(3ES24.15)')ce,fd,adj
 if(ce>2e-11_dp.or.fd>2e-5_dp.or.adj>2e-9_dp)error stop 'native nonspherical radial response failed'
contains
 subroutine read_inputs(path,d)
  character(*),intent(in)::path
  type(inputs),intent(out)::d
  integer::u,i
  real(dp)::r,boundary(2,3)
  open(newunit=u,file=path,status='old')
  read(u,*)d%n,d%mode,d%ell,d%lapw,d%nlo
  read(u,*)d%r0,d%dx,d%rmt,d%e,d%z
  allocate(d%v(d%n),d%dv(d%n),d%elo(d%nlo),d%delo(d%nlo))
  if(d%nlo>0)read(u,*)d%elo
  read(u,*)boundary;read(u,*)d%v;close(u)
  do i=1,d%n
   r=d%r0*exp(d%dx*(i-1));d%dv(i)=.07_dp*r*sin((1.1_dp+.1_dp*d%ell)*r)
  enddo
  do i=1,d%nlo
   d%delo(i)=.02_dp*cos(real(i,dp))
  enddo
 end subroutine
 subroutine prepare(d,step,ctx)
  type(inputs),intent(in)::d
  real(dp),intent(in)::step
  type(radial_basis_context),intent(out)::ctx
  call prepare_radial_basis(ctx,d%mode==1,d%v+step*d%dv,d%r0,d%dx,d%rmt,d%e+step*.03_dp, &
    d%elo+step*d%delo,d%lapw==1,real(d%ell,dp),d%z)
 end subroutine
 subroutine oracle(left,right,potential,result)
  type(radial_basis_context),intent(in)::left,right
  real(dp),intent(in)::potential(:,:)
  real(dp),intent(out)::result(2+x%nlo,2+y%nlo,nv)
  integer::i,j,q
  do q=1,nv
   do j=1,2+y%nlo
    do i=1,2+x%nlo
     result(i,j,q)=RINT13(x%mode==1,left%primitive(:,1,i)*potential(:,q), &
      left%primitive(:,2,i)*potential(:,q),right%primitive(:,1,j),right%primitive(:,2,j),x%n,x%dx,x%n,x%r0)
    enddo
   enddo
  enddo
 end subroutine
end program
