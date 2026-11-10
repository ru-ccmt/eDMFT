module lapw_radial_hamiltonian
  ! Differentiate the radial algebra used by LAPW1 HAMILTON, not a continuum
  ! replacement H*udot=E*udot+u applied only to an arbitrary response bra.
  ! Includes ordinary LO energy averages, APW boundary flux, matching/ABC
  ! changes, and Vphysical-Vbasis restoration on a fixed radial mesh.
  !
  ! H is the UNFOLDED radial contribution, before angular factors and the
  ! native upper-triangle-to-Hermitian copy. LAPW's two u/udot surface
  ! orientations differ by W_R-1 for the implemented finite stencil. Do not
  ! silently symmetrize this small difference before coefficient assembly.
  ! A caller must assemble it in native PW/LO order, fold the same triangle,
  ! and only then transform to the central KS basis. This matters for the
  ! exact transpose too. This module does not claim to fix native conventions.
  !
  ! Primitive ordering is (u, udot, u_LO1, ...); M maps these primitives
  ! to complete matched PW radial factors followed by complete LOs.
  ! With S the primitive Gram matrix and T=diag(E,E,E_LO,...)+e1*e2^T,
  !   O = M^T S M,
  !   H_volume = M^T [(S T + T^T S)/2 + DeltaV] M.
  ! The symmetric boundary flux and the ordered Wronskian correction are
  ! added below. DeltaV integrates Vphysical-Vbasis, in Ry. All energy
  ! arguments are Ry; potential arguments store r*V, also in Ry.
  ! For a physical potential change in both operators, d(DeltaV)=0.
  ! For a basis-only test with the physical H fixed, d(DeltaV)=-dVbasis.
  ! These are different probes and must not be confused in validation.
  use iso_fortran_env, only: real64
  use lapw_radial_response
  implicit none
  private
  integer,parameter :: dp=real64
  public :: radial_hamiltonian_value,radial_hamiltonian_action,radial_hamiltonian_pullback
  public :: radial_hamiltonian_component_action,radial_hamiltonian_component_pullback
contains
  subroutine algebra(ctx,boundary,energy,lo_energy,rv_difference,m,s,t,k,p,d,w)
    type(radial_basis_context),intent(in) :: ctx
    real(dp),intent(in) :: boundary(:,:),energy,lo_energy(ctx%nlo),rv_difference(ctx%n)
    real(dp),intent(out) :: m(2+ctx%nlo,size(boundary,2)+ctx%nlo)
    real(dp),intent(out) :: s(2+ctx%nlo,2+ctx%nlo),t(2+ctx%nlo,2+ctx%nlo),k(2+ctx%nlo,2+ctx%nlo)
    real(dp),intent(out) :: p(size(boundary,2)+ctx%nlo),d(size(boundary,2)+ctx%nlo),w
    real(dp) :: weight(ctx%n,4)
    integer :: i,j,np
    np=size(boundary,2)
    m=radial_component_coefficients(ctx,boundary)
    do i=1,size(s,1)
       do j=1,size(s,2)
          s(i,j)=sum(ctx%overlap_metric*ctx%primitive(:,:,i)*ctx%primitive(:,:,j))
       enddo
    enddo
    ! LAPW1 uses normalization/orthogonalization identities algebraically in
    ! its central 2x2 block. Keep those entries constant also in the tangent.
    s(1,1)=1;s(1,2)=0;s(2,1)=0
    t=0;t(1,1)=energy;t(2,2)=energy;t(1,2)=1
    do i=1,ctx%nlo
       t(2+i,2+i)=lo_energy(i)
    enddo
    k=.5_dp*(matmul(s,t)+matmul(transpose(t),s))
    ! Integrate in the SMALL primitive space once. Expanding radial functions
    ! for every plane wave would cost O(Nrad*Npw^2) per response action.
    weight=ctx%overlap_metric*spread(rv_difference/ctx%pair%center%r,2,4)
    do i=1,size(s,1)
       do j=1,size(s,2)
          k(i,j)=k(i,j)+sum(weight*ctx%primitive(:,:,i)*ctx%primitive(:,:,j))
       enddo
    enddo
    p=matmul(ctx%primitive(ctx%n,1,:)/ctx%rmt,m)
    d=matmul(ctx%primitive(ctx%n,3,:)/ctx%rmt,m)
    ! HAMILTON explicitly omits LO boundary values, even if the finite udot
    ! stencil makes the computed ABC combination nonzero at roundoff/stencil
    ! accuracy. Retain its actual definition for an exact algorithm derivative.
    p(np+1:)=0
    w=ctx%primitive(ctx%n,3,1)*ctx%primitive(ctx%n,1,2) &
      -ctx%primitive(ctx%n,1,1)*ctx%primitive(ctx%n,3,2)
  end subroutine

  subroutine radial_hamiltonian_value(ctx,boundary,energy,lo_energy,rv_difference,h,o)
    type(radial_basis_context),intent(in) :: ctx
    real(dp),intent(in) :: boundary(:,:),energy,lo_energy(ctx%nlo),rv_difference(ctx%n)
    real(dp),intent(out) :: h(size(boundary,2)+ctx%nlo,size(boundary,2)+ctx%nlo),o(size(h,1),size(h,2))
    real(dp) :: m(2+ctx%nlo,size(h,1)),s(2+ctx%nlo,2+ctx%nlo),t(2+ctx%nlo,2+ctx%nlo),k(2+ctx%nlo,2+ctx%nlo)
    real(dp) :: p(size(h,1)),d(size(h,1)),w
    integer :: i,j,np
    np=size(boundary,2)
    call algebra(ctx,boundary,energy,lo_energy,rv_difference,m,s,t,k,p,d,w)
    o=matmul(transpose(m),matmul(s,m));h=matmul(transpose(m),matmul(k,m))
    do i=1,size(h,1)
       do j=1,size(h,2)
          h(i,j)=h(i,j)+.5_dp*ctx%rmt**2*(p(i)*d(j)+d(i)*p(j))
       enddo
    enddo
    do i=1,np
       do j=1,np
          h(i,j)=h(i,j)+.5_dp*(w-1)*(m(1,i)*m(2,j)-m(2,i)*m(1,j))
       enddo
    enddo
  end subroutine

  subroutine radial_hamiltonian_action(ctx,boundary,energy,lo_energy,rv_difference, &
       dv,de,delo,delta_rv_difference,dh,doverlap)
    type(radial_basis_context),intent(in) :: ctx
    real(dp),intent(in) :: boundary(:,:),energy,lo_energy(ctx%nlo),rv_difference(ctx%n)
    real(dp),intent(in) :: dv(ctx%n),de,delo(ctx%nlo),delta_rv_difference(ctx%n)
    real(dp),intent(out) :: dh(size(boundary,2)+ctx%nlo,size(boundary,2)+ctx%nlo),doverlap(size(dh,1),size(dh,2))
    real(dp) :: u(ctx%n,4,2+ctx%nlo),dm(2+ctx%nlo,size(dh,1))
    call radial_components_action(ctx,boundary,dv,de,delo,u,dm)
    call radial_hamiltonian_component_action(ctx,boundary,energy,lo_energy,rv_difference, &
         u,dm,de,delo,delta_rv_difference,dh,doverlap)
  end subroutine

  ! Low-level shared component interface: evaluate the radial tangent ONCE
  ! per atom/l and direction, then feed it to H/O, projector and density
  ! assembly. These routines do not call the radial differential equation.
  subroutine radial_hamiltonian_component_action(ctx,boundary,energy,lo_energy,rv_difference, &
       u,dm,de,delo,delta_rv_difference,dh,doverlap)
    type(radial_basis_context),intent(in) :: ctx
    real(dp),intent(in) :: boundary(:,:),energy,lo_energy(ctx%nlo),rv_difference(ctx%n)
    real(dp),intent(in) :: u(ctx%n,4,2+ctx%nlo),dm(2+ctx%nlo,size(boundary,2)+ctx%nlo)
    real(dp),intent(in) :: de,delo(ctx%nlo),delta_rv_difference(ctx%n)
    real(dp),intent(out) :: dh(size(boundary,2)+ctx%nlo,size(boundary,2)+ctx%nlo),doverlap(size(dh,1),size(dh,2))
    real(dp) :: m(2+ctx%nlo,size(dh,1))
    real(dp) :: s(2+ctx%nlo,2+ctx%nlo),t(2+ctx%nlo,2+ctx%nlo),k(2+ctx%nlo,2+ctx%nlo)
    real(dp) :: ds(2+ctx%nlo,2+ctx%nlo),dt(2+ctx%nlo,2+ctx%nlo),dk(2+ctx%nlo,2+ctx%nlo)
    real(dp) :: p(size(dh,1)),d(size(dh,1)),pd(size(dh,1)),dd(size(dh,1)),w,dw
    real(dp) :: weight(ctx%n,4),dwgt(ctx%n,4)
    integer :: i,j,b,np,n
    np=size(boundary,2);n=ctx%n
    call algebra(ctx,boundary,energy,lo_energy,rv_difference,m,s,t,k,p,d,w)
    do i=1,size(s,1)
       do j=1,size(s,2)
          ds(i,j)=sum(ctx%overlap_metric*(u(:,:,i)*ctx%primitive(:,:,j)+ctx%primitive(:,:,i)*u(:,:,j)))
       enddo
    enddo
    ds(1,1)=0;ds(1,2)=0;ds(2,1)=0
    dt=0;dt(1,1)=de;dt(2,2)=de
    do i=1,ctx%nlo
       dt(2+i,2+i)=delo(i)
    enddo
    dk=.5_dp*(matmul(ds,t)+matmul(s,dt)+matmul(transpose(dt),s)+matmul(transpose(t),ds))
    weight=ctx%overlap_metric*spread(rv_difference/ctx%pair%center%r,2,4)
    dwgt=ctx%overlap_metric*spread(delta_rv_difference/ctx%pair%center%r,2,4)
    do i=1,size(s,1)
       do j=1,size(s,2)
          dk(i,j)=dk(i,j)+sum(weight*(u(:,:,i)*ctx%primitive(:,:,j)+ctx%primitive(:,:,i)*u(:,:,j)) &
                    +dwgt*ctx%primitive(:,:,i)*ctx%primitive(:,:,j))
       enddo
    enddo
    doverlap=matmul(transpose(dm),matmul(s,m))+matmul(transpose(m),matmul(ds,m)+matmul(s,dm))
    dh=matmul(transpose(dm),matmul(k,m))+matmul(transpose(m),matmul(dk,m)+matmul(k,dm))
    pd=(matmul(u(n,1,:),m)+matmul(ctx%primitive(n,1,:),dm))/ctx%rmt
    dd=(matmul(u(n,3,:),m)+matmul(ctx%primitive(n,3,:),dm))/ctx%rmt
    pd(np+1:)=0
    dw=u(n,3,1)*ctx%primitive(n,1,2)+ctx%primitive(n,3,1)*u(n,1,2) &
       -u(n,1,1)*ctx%primitive(n,3,2)-ctx%primitive(n,1,1)*u(n,3,2)
    do i=1,size(dh,1)
       do j=1,size(dh,2)
          dh(i,j)=dh(i,j)+.5_dp*ctx%rmt**2*(pd(i)*d(j)+p(i)*dd(j)+dd(i)*p(j)+d(i)*pd(j))
       enddo
    enddo
    do i=1,np
       do j=1,np
          dh(i,j)=dh(i,j)+.5_dp*dw*(m(1,i)*m(2,j)-m(2,i)*m(1,j)) &
               +.5_dp*(w-1)*(dm(1,i)*m(2,j)+m(1,i)*dm(2,j)-dm(2,i)*m(1,j)-m(2,i)*dm(1,j))
       enddo
    enddo
  end subroutine

  subroutine radial_hamiltonian_pullback(ctx,boundary,energy,lo_energy,rv_difference, &
       hbar,obar,dvbar,debar,delobar,rv_difference_bar)
    type(radial_basis_context),intent(in) :: ctx
    real(dp),intent(in) :: boundary(:,:),energy,lo_energy(ctx%nlo),rv_difference(ctx%n)
    real(dp),intent(in) :: hbar(size(boundary,2)+ctx%nlo,size(boundary,2)+ctx%nlo),obar(size(hbar,1),size(hbar,2))
    real(dp),intent(out) :: dvbar(ctx%n),debar,delobar(ctx%nlo),rv_difference_bar(ctx%n)
    real(dp) :: gp(ctx%n,4,2+ctx%nlo),gm(2+ctx%nlo,size(hbar,1)),eb,elb(ctx%nlo)
    call radial_hamiltonian_component_pullback(ctx,boundary,energy,lo_energy,rv_difference, &
         hbar,obar,gp,gm,eb,elb,rv_difference_bar)
    call radial_components_pullback(ctx,boundary,gp,gm,dvbar,debar,delobar)
    debar=debar+eb;delobar=delobar+elb
  end subroutine

  ! Sum gp/gm with the projector and density covectors before ONE call to
  ! radial_components_pullback. eb/elb here contain only the explicit energy
  ! coefficients of HAMILTON; the radial energy response is added separately.
  subroutine radial_hamiltonian_component_pullback(ctx,boundary,energy,lo_energy,rv_difference, &
       hbar,obar,gp,gm,debar,delobar,rv_difference_bar)
    type(radial_basis_context),intent(in) :: ctx
    real(dp),intent(in) :: boundary(:,:),energy,lo_energy(ctx%nlo),rv_difference(ctx%n)
    real(dp),intent(in) :: hbar(size(boundary,2)+ctx%nlo,size(boundary,2)+ctx%nlo),obar(size(hbar,1),size(hbar,2))
    real(dp),intent(out) :: gp(ctx%n,4,2+ctx%nlo),gm(2+ctx%nlo,size(hbar,1))
    real(dp),intent(out) :: debar,delobar(ctx%nlo),rv_difference_bar(ctx%n)
    real(dp) :: m(2+ctx%nlo,size(hbar,1))
    real(dp) :: s(2+ctx%nlo,2+ctx%nlo),t(2+ctx%nlo,2+ctx%nlo),k(2+ctx%nlo,2+ctx%nlo)
    real(dp) :: gs(2+ctx%nlo,2+ctx%nlo),gt(2+ctx%nlo,2+ctx%nlo),gk(2+ctx%nlo,2+ctx%nlo)
    real(dp) :: p(size(hbar,1)),d(size(hbar,1)),pb(size(hbar,1)),db(size(hbar,1)),w,wb,anti
    real(dp) :: weight(ctx%n,4),temp(ctx%n,4)
    integer :: i,j,b,np,n
    np=size(boundary,2);n=ctx%n
    call algebra(ctx,boundary,energy,lo_energy,rv_difference,m,s,t,k,p,d,w)
    gm=matmul(k,matmul(m,hbar+transpose(hbar)))+matmul(s,matmul(m,obar+transpose(obar)))
    gk=matmul(m,matmul(hbar,transpose(m)))
    gk=.5_dp*(gk+transpose(gk))
    gs=matmul(m,matmul(obar,transpose(m)))+.5_dp*(matmul(gk,transpose(t))+matmul(t,gk))
    gt=matmul(s,gk)
    gs(1,1)=0;gs(1,2)=0;gs(2,1)=0
    gp=0
    do i=1,size(s,1)
       do j=1,size(s,2)
          gp(:,:,i)=gp(:,:,i)+(gs(i,j)+gs(j,i))*ctx%overlap_metric*ctx%primitive(:,:,j)
       enddo
    enddo
    pb=.5_dp*ctx%rmt**2*matmul(hbar+transpose(hbar),d)
    db=.5_dp*ctx%rmt**2*matmul(hbar+transpose(hbar),p)
    pb(np+1:)=0
    do b=1,size(m,2)
       gm(:,b)=gm(:,b)+(ctx%primitive(n,1,:)*pb(b)+ctx%primitive(n,3,:)*db(b))/ctx%rmt
    enddo
    gp(n,1,:)=gp(n,1,:)+matmul(m,pb)/ctx%rmt
    gp(n,3,:)=gp(n,3,:)+matmul(m,db)/ctx%rmt
    wb=0
    do i=1,np
       do j=1,np
          anti=.5_dp*hbar(i,j)*(w-1)
          gm(1,i)=gm(1,i)+anti*m(2,j);gm(2,j)=gm(2,j)+anti*m(1,i)
          gm(2,i)=gm(2,i)-anti*m(1,j);gm(1,j)=gm(1,j)-anti*m(2,i)
          wb=wb+.5_dp*hbar(i,j)*(m(1,i)*m(2,j)-m(2,i)*m(1,j))
       enddo
    enddo
    gp(n,3,1)=gp(n,3,1)+wb*ctx%primitive(n,1,2);gp(n,1,2)=gp(n,1,2)+wb*ctx%primitive(n,3,1)
    gp(n,1,1)=gp(n,1,1)-wb*ctx%primitive(n,3,2);gp(n,3,2)=gp(n,3,2)-wb*ctx%primitive(n,1,1)
    rv_difference_bar=0
    weight=ctx%overlap_metric*spread(rv_difference/ctx%pair%center%r,2,4)
    do i=1,size(s,1)
       do j=1,size(s,2)
          gp(:,:,i)=gp(:,:,i)+(gk(i,j)+gk(j,i))*weight*ctx%primitive(:,:,j)
          temp=ctx%overlap_metric*ctx%primitive(:,:,i)*ctx%primitive(:,:,j)
          rv_difference_bar=rv_difference_bar+gk(i,j)*sum(temp,dim=2)/ctx%pair%center%r
       enddo
    enddo
    debar=gt(1,1)+gt(2,2)
    do i=1,ctx%nlo
       delobar(i)=gt(2+i,2+i)
    enddo
  end subroutine
end module lapw_radial_hamiltonian
