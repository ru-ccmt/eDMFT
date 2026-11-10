program test_lapw_radial_response
  ! Independent endpoint oracle: ATPARN, ABC, RINT13 and Rint13n. No endpoint
  ! invokes the new assembly module. Test both primitive and complete LO
  ! factors, all bras/kets of O, raw projector factors, and their transpose.
  use iso_fortran_env, only: real64
  use lapw_radial_response
  implicit none
  integer,parameter :: dp=real64,np=3
  type(radial_basis_context) :: ctx
  integer :: n,mode,ell,lapw_flag,nlo,m,i,j,k,it,unit
  real(dp) :: r0,dx,rmt,energy,z,step,r,de,be,left,right,t0,t1,t2
  real(dp) :: central_error,fd_error,form_error,adj_error,chain_error,wronskian,wronskian_dot
  real(dp),allocatable :: v(:),dv(:),elo(:),delo(:),bev(:),belo(:),projector(:,:),boundary(:,:)
  real(dp),allocatable :: value(:,:,:),reference(:,:,:),plus(:,:,:),minus(:,:,:),derivative(:,:,:),bar(:,:,:)
  real(dp),allocatable :: ov(:,:),proj(:),dov(:,:),dproj(:),op(:,:),om(:,:),pp(:),pm(:),ob(:,:),pb(:)
  real(dp),external :: RINT13,Rint13n
  character(1024) :: path
  call get_command_argument(1,path)
  open(newunit=unit,file=trim(path),status='old')
  read(unit,*) n,mode,ell,lapw_flag,nlo
  read(unit,*) r0,dx,rmt,energy,z
  m=np+nlo
  allocate(v(n),dv(n),elo(nlo),delo(nlo),bev(n),belo(nlo),projector(n,2),boundary(2,np))
  allocate(value(n,4,m),reference(n,4,m),plus(n,4,m),minus(n,4,m),derivative(n,4,m),bar(n,4,m))
  allocate(ov(m,m),proj(m),dov(m,m),dproj(m),op(m,m),om(m,m),pp(m),pm(m),ob(m,m),pb(m))
  if(nlo>0) read(unit,*) elo
  read(unit,*) boundary
  read(unit,*) v
  read(unit,*) projector
  close(unit)
  do i=1,n
     r=r0*exp(dx*(i-1))
     dv(i)=.07_dp*r*sin(1.3_dp*r)
     do j=1,4
        do k=1,m
           bar(i,j,k)=sin(.2_dp*i+j+.4_dp*k)*exp(-r)
        enddo
     enddo
  enddo
  de=.03_dp
  do i=1,nlo
     delo(i)=.023_dp*cos(real(i,dp))
  enddo
  do i=1,m
     pb(i)=sin(real(i,dp))
     do j=1,m
        ob(i,j)=cos(2.3_dp*i+.7_dp*j) ! deliberately nonsymmetric covector
     enddo
  enddo
  call prepare_radial_basis(ctx,mode==1,v,r0,dx,rmt,energy,elo,lapw_flag==1,real(ell,dp),z)
  value=radial_basis_value(ctx,boundary)
  call endpoint(v,energy,elo,reference)
  central_error=sqrt(sum((reference-value)**2)/sum(reference**2))
  call radial_basis_forms(ctx,value,projector,ov,proj)
  call reference_forms(reference,op,pp)
  central_error=max(central_error,sqrt(sum((ov-op)**2)/sum(op**2)),sqrt(sum((proj-pp)**2)/sum(pp**2)))
  call radial_basis_action(ctx,boundary,dv,de,delo,derivative)
  call radial_basis_pullback(ctx,boundary,bar,bev,be,belo)
  left=sum(bar*derivative);right=sum(dv*bev)+de*be+sum(delo*belo)
  adj_error=abs(left-right)/max(1.0_dp,abs(left),abs(right))
  call radial_basis_forms_action(ctx,value,derivative,projector,dov,dproj)
  call radial_basis_forms_pullback(ctx,value,projector,ob,pb,bar)
  call radial_basis_pullback(ctx,boundary,bar,bev,be,belo)
  left=sum(ob*dov)+sum(pb*dproj);right=sum(dv*bev)+de*be+sum(delo*belo)
  chain_error=abs(left-right)/max(1.0_dp,abs(left),abs(right))
  ! Probe three steps: energy-stencil cancellation eventually sets a floor.
  fd_error=huge(1.0_dp);form_error=huge(1.0_dp)
  do it=1,3
     step=.004_dp/(2**it)
     call endpoint(v+step*dv,energy+step*de,elo+step*delo,plus)
     call endpoint(v-step*dv,energy-step*de,elo-step*delo,minus)
     call reference_forms(plus,op,pp);call reference_forms(minus,om,pm)
     fd_error=min(fd_error,sqrt(sum(((plus-minus)/(2*step)-derivative)**2)/sum(derivative**2)))
     form_error=min(form_error,max(sqrt(sum(((op-om)/(2*step)-dov)**2)/sum(dov**2)), &
                                  sqrt(sum(((pp-pm)/(2*step)-dproj)**2)/sum(dproj**2))))
  enddo
  ! In LAPW the boundary factor is W=R^2(dp*pe-p*dpe), not an
  ! exact algebraic inverse. The actual derivative need not vanish there.
  wronskian=ctx%primitive(n,3,1)*ctx%primitive(n,1,2)-ctx%primitive(n,1,1)*ctx%primitive(n,3,2)
  wronskian_dot=0
  if(lapw_flag==1) wronskian_dot=derivative(n,1,1)/(rmt*boundary(1,1))
  call cpu_time(t0)
  do it=1,200
     call radial_basis_action(ctx,boundary,dv,de,delo,derivative)
  enddo
  call cpu_time(t1)
  do it=1,200
     call radial_basis_pullback(ctx,boundary,bar,bev,be,belo)
  enddo
  call cpu_time(t2)
  write(*,'(9ES24.15)') central_error,fd_error,form_error,adj_error,chain_error, &
       wronskian,wronskian_dot,(t1-t0)/200,(t2-t1)/200
  if(central_error>2e-10_dp.or.fd_error>2e-5_dp.or.form_error>3e-5_dp) error stop 'radial assembly FD failed'
  if(max(adj_error,chain_error)>2e-9_dp) error stop 'radial assembly adjoint failed'
contains
  subroutine endpoint(vr,e,elocal,result)
    real(dp),intent(in) :: vr(n),e,elocal(nlo)
    real(dp),intent(out) :: result(n,4,m)
    integer :: nl,j,nt
    real(dp) :: al,bl
    real(dp),allocatable :: u(:,:),du(:,:),b(:,:),db(:,:),up(:,:),dup(:,:),bp(:,:),dbp(:,:)
    real(dp),allocatable :: ul(:,:,:),bul(:,:,:),upl(:,:,:),bpl(:,:,:),p(:,:),s(:,:),ri(:,:,:)
    real(dp),allocatable :: ca(:,:),cb(:,:),cc(:,:),es(:),els(:,:)
    integer :: ilo(0:ell)
    logical,allocatable :: loor(:,:),rlo(:,:),lapw(:)
    nl=max(nlo,1)
    allocate(u(n,0:ell),du(n,0:ell),b(n,0:ell),db(n,0:ell),up(n,0:ell),dup(n,0:ell),bp(n,0:ell),dbp(n,0:ell))
    allocate(ul(n,nl,0:ell),bul(n,nl,0:ell),upl(n,nl,0:ell),bpl(n,nl,0:ell))
    allocate(p(nl+2,0:ell),s(nl+2,0:ell),ri(nl+2,nl+2,0:ell),ca(0:ell,nl),cb(0:ell,nl),cc(0:ell,nl))
    allocate(es(0:ell),els(0:ell,nl),loor(nl,0:ell),rlo(nl,0:ell),lapw(0:ell))
    es=e;els=999;els(ell,:nlo)=elocal
    loor=.false.;loor(:nlo,ell)=.true.;rlo=.false.;lapw=.true.;lapw(ell)=lapw_flag==1
    if(.not.lapw(ell)) loor(1,ell)=.false.
    ilo=0;ilo(ell)=nlo
    nt=nlo*(2*ell+1)
    call ATPARN(u,du,b,db,up,dup,bp,dbp,ul,bul,upl,bpl,p,s,ri,ca,cb,cc, &
         vr,n,es,els,loor,rlo,lapw,0,nt,0,ilo,mode==1,z,.false.,n,r0,dx,rmt,ell,ell,nl,n)
    do j=1,np
       if(lapw(ell)) then
          al=rmt*rmt*(boundary(2,j)*p(2,ell)-boundary(1,j)*s(2,ell))
          bl=rmt*rmt*(boundary(1,j)*s(1,ell)-boundary(2,j)*p(1,ell))
       else
          al=boundary(1,j)/p(1,ell);bl=0
       endif
       result(:,1,j)=al*u(:,ell)+bl*du(:,ell);result(:,2,j)=al*b(:,ell)+bl*db(:,ell)
       result(:,3,j)=al*up(:,ell)+bl*dup(:,ell);result(:,4,j)=al*bp(:,ell)+bl*dbp(:,ell)
    enddo
    do j=1,nlo
       result(:,1,np+j)=ca(ell,j)*u(:,ell)+cb(ell,j)*du(:,ell)+cc(ell,j)*ul(:,j,ell)
       result(:,2,np+j)=ca(ell,j)*b(:,ell)+cb(ell,j)*db(:,ell)+cc(ell,j)*bul(:,j,ell)
       result(:,3,np+j)=ca(ell,j)*up(:,ell)+cb(ell,j)*dup(:,ell)+cc(ell,j)*upl(:,j,ell)
       result(:,4,np+j)=ca(ell,j)*bp(:,ell)+cb(ell,j)*dbp(:,ell)+cc(ell,j)*bpl(:,j,ell)
    enddo
  end subroutine
  subroutine reference_forms(f,o,p)
    real(dp),intent(in) :: f(n,4,m)
    real(dp),intent(out) :: o(m,m),p(m)
    integer :: i,j
    do i=1,m
       p(i)=Rint13n(n,f(:,1,i),f(:,2,i),projector(:,1),projector(:,2),r0,dx,mode==1)
       do j=1,m
          o(i,j)=RINT13(mode==1,f(:,1,i),f(:,2,i),f(:,1,j),f(:,2,j),n,dx,n,r0)
       enddo
    enddo
  end subroutine
end program
