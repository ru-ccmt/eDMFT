program test_radial_hamiltonian
  ! Independent endpoint oracle: ATPARN, ABC, RINT13 and literal HAMILTON
  ! radial coefficients. No endpoint invokes the new assembly module.
  ! Test complete LO factors, both bras/kets of H/O, independent physical
  ! and basis potentials, energy parameters, and the exact discrete transpose.
  use iso_fortran_env, only: real64
  use lapw_radial_response
  use lapw_radial_hamiltonian
  implicit none
  integer,parameter :: dp=real64,np=3
  type(radial_basis_context) :: ctx
  integer :: n,mode,ell,lapw_flag,nlo,m,i,j,k,it,unit
  real(dp) :: r0,dx,rmt,energy,z,step,r,de,be,left,right,t0,t1,t2
  real(dp) :: central_error,fd_error,form_error,adj_error,chain_error,wronskian,wronskian_dot
  real(dp),allocatable :: v(:),dv(:),elo(:),delo(:),bev(:),belo(:),projector(:,:),boundary(:,:)
  real(dp),allocatable :: value(:,:,:),reference(:,:,:),plus(:,:,:),minus(:,:,:),derivative(:,:,:),bar(:,:,:)
  real(dp),allocatable :: ov(:,:),proj(:),dov(:,:),dproj(:),op(:,:),om(:,:),pp(:),pm(:),ob(:,:),pb(:)
  real(dp),allocatable :: hd(:),hd1(:),hdb(:),h(:,:),h0(:,:),hp(:,:),hm(:,:),dh(:,:),hb(:,:)
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
  allocate(hd(n),hd1(n),hdb(n),h(m,m),h0(m,m),hp(m,m),hm(m,m),dh(m,m),hb(m,m))
  if(nlo>0) read(unit,*) elo
  read(unit,*) boundary
  read(unit,*) v
  read(unit,*) projector
  close(unit)
  do i=1,n
     r=r0*exp(dx*(i-1))
     dv(i)=.07_dp*r*sin(1.3_dp*r)
     hd(i)=.021_dp*r*cos(r);hd1(i)=-.05_dp*r*cos(2*r)
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
        hb(i,j)=sin(.9_dp*i+1.4_dp*j)
     enddo
  enddo
  call prepare_radial_basis(ctx,mode==1,v,r0,dx,rmt,energy,elo,lapw_flag==1,real(ell,dp),z)
  call radial_hamiltonian_value(ctx,boundary,energy,elo,hd,h,ov)
  call endpoint(v,energy,elo,hd,reference,h0,op)
  central_error=max(sqrt(sum((h-h0)**2)/sum(h0**2)),sqrt(sum((ov-op)**2)/sum(op**2)))
  call radial_hamiltonian_action(ctx,boundary,energy,elo,hd,dv,de,delo,hd1,dh,dov)
  call radial_hamiltonian_pullback(ctx,boundary,energy,elo,hd,hb,ob,bev,be,belo,hdb)
  left=sum(hb*dh)+sum(ob*dov);right=sum(dv*bev)+de*be+sum(delo*belo)+sum(hd1*hdb)
  adj_error=abs(left-right)/max(1.0_dp,abs(left),abs(right))
  fd_error=huge(1.0_dp);form_error=huge(1.0_dp)
  do it=1,3
     step=.004_dp/(2**it)
     call endpoint(v+step*dv,energy+step*de,elo+step*delo,hd+step*hd1,plus,hp,op)
     call endpoint(v-step*dv,energy-step*de,elo-step*delo,hd-step*hd1,minus,hm,om)
     fd_error=min(fd_error,sqrt(sum(((hp-hm)/(2*step)-dh)**2)/sum(dh**2)))
     form_error=min(form_error,sqrt(sum(((op-om)/(2*step)-dov)**2)/sum(dov**2)))
  enddo
  call cpu_time(t0)
  do it=1,100
     call radial_hamiltonian_action(ctx,boundary,energy,elo,hd,dv,de,delo,hd1,dh,dov)
  enddo
  call cpu_time(t1)
  do it=1,100
     call radial_hamiltonian_pullback(ctx,boundary,energy,elo,hd,hb,ob,bev,be,belo,hdb)
  enddo
  call cpu_time(t2)
  write(*,'(6ES24.15)') central_error,fd_error,form_error,adj_error,(t1-t0)/100,(t2-t1)/100
  if(central_error>2e-10_dp.or.fd_error>2e-5_dp.or.form_error>3e-5_dp) error stop 'native radial H/O FD failed'
  if(adj_error>2e-9_dp) error stop 'native radial H/O adjoint failed'
contains
  subroutine endpoint(vr,e,elocal,delta,result,hn,on)
    real(dp),intent(in) :: vr(n),e,elocal(nlo),delta(n)
    real(dp),intent(out) :: result(n,4,m),hn(m,m),on(m,m)
    integer :: nl,j,nt,i,jlo,klo
    real(dp) :: al,bl,aa(np),bb(np),c1,c2,c3,c11,c12,c13,flux,rad(n),wd(n)
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
       aa(j)=al;bb(j)=bl
       result(:,1,j)=al*u(:,ell)+bl*du(:,ell);result(:,2,j)=al*b(:,ell)+bl*db(:,ell)
       result(:,3,j)=al*up(:,ell)+bl*dup(:,ell);result(:,4,j)=al*bp(:,ell)+bl*dbp(:,ell)
    enddo
    do j=1,nlo
       result(:,1,np+j)=ca(ell,j)*u(:,ell)+cb(ell,j)*du(:,ell)+cc(ell,j)*ul(:,j,ell)
       result(:,2,np+j)=ca(ell,j)*b(:,ell)+cb(ell,j)*db(:,ell)+cc(ell,j)*bul(:,j,ell)
       result(:,3,np+j)=ca(ell,j)*up(:,ell)+cb(ell,j)*dup(:,ell)+cc(ell,j)*upl(:,j,ell)
       result(:,4,np+j)=ca(ell,j)*bp(:,ell)+cb(ell,j)*dbp(:,ell)+cc(ell,j)*bpl(:,j,ell)
    enddo
    ! Literal scalar radial factors of native HAMILTON, before upper-triangle
    ! folding. This independent oracle uses ATPARN and the explicit C11/C12/C13
    ! expressions, not the compact S*T + T^T*S formulation under test.
    do i=1,np
       do j=1,np
          on(i,j)=aa(i)*aa(j)+bb(i)*bb(j)*ri(2,2,ell)
          hn(i,j)=aa(i)*((e+p(1,ell)*s(1,ell)*rmt**2)*aa(j)+p(2,ell)*s(1,ell)*rmt**2*bb(j)) &
               +bb(i)*((e*ri(2,2,ell)+p(2,ell)*s(2,ell)*rmt**2)*bb(j) &
                       +(1+p(1,ell)*s(2,ell)*rmt**2)*aa(j))
       enddo
    enddo
    do jlo=1,nlo
       c1=ca(ell,jlo)+ri(1,2+jlo,ell)*cc(ell,jlo)
       c2=cb(ell,jlo)*ri(2,2,ell)+ri(2,2+jlo,ell)*cc(ell,jlo)
       c11=ca(ell,jlo)*e+.5_dp*(cb(ell,jlo)+cc(ell,jlo)*ri(1,2+jlo,ell)*(elocal(jlo)+e))
       c12=cb(ell,jlo)*ri(2,2,ell)*e+.5_dp*(ca(ell,jlo)+cc(ell,jlo)*ri(1,2+jlo,ell) &
            +cc(ell,jlo)*ri(2,2+jlo,ell)*(elocal(jlo)+e))
       flux=.5_dp*rmt**2*(ca(ell,jlo)*s(1,ell)+cb(ell,jlo)*s(2,ell)+cc(ell,jlo)*s(2+jlo,ell))
       do i=1,np
          on(i,np+jlo)=aa(i)*c1+bb(i)*c2;on(np+jlo,i)=on(i,np+jlo)
          hn(i,np+jlo)=aa(i)*(c11+flux*p(1,ell))+bb(i)*(c12+flux*p(2,ell))
          hn(np+jlo,i)=hn(i,np+jlo)
       enddo
       do klo=1,nlo
          c3=cc(ell,jlo)*ri(2+klo,2+jlo,ell)+ri(2,2+klo,ell)*cb(ell,jlo)+ca(ell,jlo)*ri(1,2+klo,ell)
          c13=cc(ell,jlo)*.5_dp*(elocal(jlo)+elocal(klo))*ri(2+klo,2+jlo,ell) &
               +.5_dp*cb(ell,jlo)*ri(1,2+klo,ell) &
               +.5_dp*(elocal(klo)+e)*(ca(ell,jlo)*ri(1,2+klo,ell)+cb(ell,jlo)*ri(2,2+klo,ell))
          on(np+klo,np+jlo)=ca(ell,klo)*c1+cb(ell,klo)*c2+cc(ell,klo)*c3
          hn(np+klo,np+jlo)=ca(ell,klo)*c11+cb(ell,klo)*c12+cc(ell,klo)*c13
       enddo
    enddo
    do i=1,n
       rad(i)=r0*exp(dx*(i-1));wd(i)=delta(i)/rad(i)
    enddo
    do i=1,m
       do j=1,m
          hn(i,j)=hn(i,j)+RINT13(mode==1,wd*result(:,1,i),wd*result(:,2,i), &
                result(:,1,j),result(:,2,j),n,dx,n,r0)
       enddo
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
