module outwin_response
  ! Exact directional and reverse derivatives of the DISCRETE OUTWIN solver.
  !
  ! We reuse OUTWIN for the central solution and differentiate its implicit
  ! Adams--Moulton steps. No continuous (H-E)udot=u identity, potential finite
  ! difference, or dense radial response matrix is substituted for the actual
  ! algorithm. The first propagated point retains its fourth-order startup.
  ! Both NREL and scalar RELA use exactly OUTWIN's small-component convention.
  !
  ! Scope: fixed radial grid, l, Z and REL flag; derivatives with respect to
  ! V=r*V_spherical in Hartree and the supplied energy E in Hartree. Automatic
  ! LAPW1 energy searches, DIRACOUT p1/2 orbitals, local-orbital ABC matching,
  ! and moving-sphere derivatives are separate operations, not included here.
  !
  ! Normalized and energy-derivative functions use ATPARN's RINT13 metric,
  ! +/-0.002 Hartree energy stencil, and final orthogonalization. Optional
  ! quadrature=1 uses default WIEN2k24 CHARGE weights instead, consistently
  ! in every normalization and derivative; ordinary ATPARN stays unchanged.
  ! Returned
  ! channels are (A,B,AP,BP), with A=r*u and AP=r*d(u)/dr. Boundary values
  ! VAL=A(n)/r(n), SLO=AP(n)/r(n) can be folded into the last-point covectors.
  use iso_fortran_env, only: real64
  implicit none
  private
  integer, parameter :: dp=real64
  real(dp), parameter :: clight=137.0359895_dp
  public :: radial_context, radial_pair_context
  public :: prepare_radial, radial_value, radial_action, radial_pullback
  public :: prepare_radial_pair, radial_pair_value
  public :: radial_pair_action, radial_pair_pullback
  public :: local_orbital_coefficients

  type :: radial_context
     integer :: n=0
     real(dp) :: dx=0, c=0, norm=0
     real(dp), allocatable :: r(:), y(:,:), f(:,:), output(:,:), metric(:,:)
     real(dp), allocatable :: matrix(:,:,:), phi_slope(:), phi_v(:), phi_e(:)
  end type

  type :: radial_pair_context
     type(radial_context) :: center, plus, minus
     real(dp), allocatable :: u(:,:), raw_dot(:,:)
     real(dp) :: cross=0
     real(dp) :: stencil=0.25_dp/0.002_dp
  end type

contains

  subroutine prepare_radial(ctx, rel, potential, r0, dx, energy, ell, charge, quadrature)
    type(radial_context), intent(out) :: ctx
    logical, intent(in) :: rel
    real(dp), intent(in) :: potential(:), r0, dx, energy, ell, charge
    real(dp) :: a(size(potential)), b(size(potential)), ap(size(potential)), bp(size(potential))
    real(dp) :: val, slo, phi, u, x, yy, cin, radius, ratio, rfirst
    integer :: n, k, nodes, j, first
    integer, intent(in), optional :: quadrature
    integer :: policy
    n=size(potential)
    policy=0
    if(present(quadrature)) policy=quadrature
    if(policy<0.or.policy>1) error stop "unknown radial quadrature"
    if(policy==1.and.n<8) error stop "WIEN24 radial quadrature needs at least eight points"
    if (n<4.or.r0<=0.or.dx<=0) error stop 'invalid radial response grid'
    ctx%n=n; ctx%dx=dx
    ctx%c=2*clight
    if (.not.rel) ctx%c=1.0e10_dp
    allocate(ctx%r(n), ctx%y(n,2), ctx%f(n,2), ctx%output(n,4), ctx%metric(n,4))
    allocate(ctx%matrix(2,2,n), ctx%phi_slope(n), ctx%phi_v(n), ctx%phi_e(n))
    call outwin(a,b,ap,bp,nodes,val,slo,rel,potential,r0,dx,n,energy,ell,charge,n)
    do k=1,n
       ctx%r(k)=r0*exp(dx*real(k-1,dp))
    enddo
    ctx%output(:,1)=a; ctx%output(:,2)=b
    ctx%output(:,3)=ap; ctx%output(:,4)=bp
    ctx%y(:,1)=a; ctx%y(:,2)=b*2/ctx%c
    ctx%f(:,1)=dx*(ctx%r*ap+a)
    ctx%f(:,2)=dx*(ctx%r*bp+b)*2/ctx%c
    do k=1,n
       radius=ctx%r(k)
       phi=(2*energy-2*potential(k)/radius)*dx*radius/ctx%c
       u=dx*radius*ctx%c+phi
       x=-dx
       yy=-ell*(ell+1)*x*x/u+phi
       ctx%matrix(1,1,k)=-x; ctx%matrix(1,2,k)=u
       ctx%matrix(2,1,k)=-yy; ctx%matrix(2,2,k)=x
       ctx%phi_slope(k)=1+ell*(ell+1)*x*x/(u*u)
       ctx%phi_v(k)=-2*dx/ctx%c
       ctx%phi_e(k)=2*dx*radius/ctx%c
    enddo
    ! Reproduce RINT13's actual end weights, including its origin term and
    ! even-grid startup correction. Its NREL small metric is 1e-22, not zero.
    ctx%metric=0
    cin=1/(clight*clight)
    if (.not.rel) cin=1.0e-22_dp
    j=3-mod(n,2); first=j-1
    ratio=exp(dx); radius=r0*ratio**(j-1); rfirst=radius/ratio
    do
       ctx%metric(j,1)=ctx%metric(j,1)+4*dx*radius/3
       radius=radius*ratio; j=j+1
       if (j>=n) exit
       ctx%metric(j,1)=ctx%metric(j,1)+2*dx*radius/3
       radius=radius*ratio; j=j+1
    enddo
    ctx%metric(j,1)=ctx%metric(j,1)+dx*radius/3
    ctx%metric(first,1)=ctx%metric(first,1)+dx*rfirst/3
    ctx%metric(1,1)=ctx%metric(1,1)+r0/3
    if (first>1) then
       ctx%metric(1,1)=ctx%metric(1,1)+0.5_dp*dx*r0
       ctx%metric(first,1)=ctx%metric(first,1)+0.5_dp*dx*rfirst
    endif
    if(policy==1) then
       ! WIEN2k 24 RINT13 calls CHARGE despite the historical UseSpline name.
       ! Its default is the flat extended-Simpson rule in log(r), with NO
       ! origin correction. This metric must enter normalization, udot
       ! orthogonalization, ABC factors and their forward/reverse derivatives.
       ! dmft2 ATPARN still uses policy=0 above; do not silently change it.
       ctx%metric(:,1)=dx*ctx%r
       ctx%metric(1:4,1)=ctx%metric(1:4,1)*[17.0_dp,59.0_dp,43.0_dp,49.0_dp]/48
       ctx%metric(n-3:n,1)=ctx%metric(n-3:n,1)*[49.0_dp,43.0_dp,59.0_dp,17.0_dp]/48
    endif
    ctx%metric(:,2)=cin*ctx%metric(:,1)
    ctx%norm=1/sqrt(sum(ctx%metric*ctx%output**2))
  end subroutine

  subroutine history(k, scale, weight, first)
    integer, intent(in) :: k
    real(dp), intent(out) :: scale, weight(4)
    integer, intent(out) :: first
    if (k==4) then
       scale=8.0_dp/3
       weight=[1.0_dp,-5.0_dp,19.0_dp,0.0_dp]/9
       first=1
    else
       scale=720.0_dp/251
       weight=[-19.0_dp,106.0_dp,-264.0_dp,646.0_dp]/251
       first=k-4
    endif
  end subroutine

  function solve_step(matrix, scale, rhs, transposed) result(answer)
    real(dp), intent(in) :: matrix(2,2), scale, rhs(2)
    logical, intent(in) :: transposed
    real(dp) :: answer(2), a(2,2), determinant
    a=-matrix
    a(1,1)=a(1,1)+scale; a(2,2)=a(2,2)+scale
    if (transposed) a=transpose(a)
    determinant=a(1,1)*a(2,2)-a(1,2)*a(2,1)
    answer(1)=(a(2,2)*rhs(1)-a(1,2)*rhs(2))/determinant
    answer(2)=(a(1,1)*rhs(2)-a(2,1)*rhs(1))/determinant
  end function

  function radial_value(ctx, normalized) result(value)
    type(radial_context), intent(in) :: ctx
    logical, intent(in) :: normalized
    real(dp) :: value(ctx%n,4)
    value=ctx%output
    if (normalized) value=ctx%norm*value
  end function

  subroutine radial_action(ctx, dv, de, derivative, normalized)
    type(radial_context), intent(in) :: ctx
    real(dp), intent(in) :: dv(ctx%n), de
    real(dp), intent(out) :: derivative(ctx%n,4)
    logical, intent(in) :: normalized
    real(dp) :: dy(ctx%n,2), df(ctx%n,2), dl_y(2), rhs(2), phi, scale, weight(4), dn
    integer :: k,j,first
    ! Origin seeds are independent of V and E in the implemented OUTWIN.
    dy=0; df=0
    do k=4,ctx%n
       call history(k,scale,weight,first)
       rhs=scale*dy(k-1,:)
       do j=first,k-1
          rhs=rhs+weight(j-first+1)*df(j,:)
       enddo
       phi=ctx%phi_v(k)*dv(k)+ctx%phi_e(k)*de
       dl_y=[phi*ctx%y(k,2), -ctx%phi_slope(k)*phi*ctx%y(k,1)]
       dy(k,:)=solve_step(ctx%matrix(:,:,k),scale,rhs+dl_y,.false.)
       df(k,:)=matmul(ctx%matrix(:,:,k),dy(k,:))+dl_y
    enddo
    derivative(:,1)=dy(:,1)
    derivative(:,2)=dy(:,2)*ctx%c/2
    derivative(:,3)=(df(:,1)/ctx%dx-dy(:,1))/ctx%r
    derivative(:,4)=(df(:,2)/ctx%dx-dy(:,2))/ctx%r*ctx%c/2
    if (normalized) then
       dn=-ctx%norm**3*sum(ctx%metric*ctx%output*derivative)
       derivative=ctx%norm*derivative+dn*ctx%output
    endif
  end subroutine

  subroutine radial_pullback(ctx, covector, dv_bar, de_bar, normalized)
    type(radial_context), intent(in) :: ctx
    real(dp), intent(in) :: covector(ctx%n,4)
    real(dp), intent(out) :: dv_bar(ctx%n), de_bar
    logical, intent(in) :: normalized
    real(dp) :: gradient(ctx%n,4), by(ctx%n,2), bf(ctx%n,2)
    real(dp) :: rhs_bar(2), combined(2), scale, weight(4), phi_bar
    integer :: k,j,first
    gradient=covector
    if (normalized) gradient=ctx%norm*gradient-ctx%norm**3 &
         *sum(covector*ctx%output)*ctx%metric*ctx%output
    by(:,1)=gradient(:,1)-gradient(:,3)/ctx%r
    by(:,2)=(gradient(:,2)-gradient(:,4)/ctx%r)*ctx%c/2
    bf(:,1)=gradient(:,3)/(ctx%dx*ctx%r)
    bf(:,2)=gradient(:,4)/(ctx%dx*ctx%r)*ctx%c/2
    dv_bar=0; de_bar=0
    do k=ctx%n,4,-1
       call history(k,scale,weight,first)
       ! f_k=L_k*y_k followed the implicit solve. Reverse it first, then
       ! solve (scale*I-L_k)^T for the implicit-step covector.
       by(k,:)=by(k,:)+matmul(transpose(ctx%matrix(:,:,k)),bf(k,:))
       combined=by(k,:)
       rhs_bar=solve_step(ctx%matrix(:,:,k),scale,combined,.true.)
       combined=bf(k,:)+rhs_bar
       phi_bar=combined(1)*ctx%y(k,2)-ctx%phi_slope(k)*combined(2)*ctx%y(k,1)
       dv_bar(k)=ctx%phi_v(k)*phi_bar
       de_bar=de_bar+ctx%phi_e(k)*phi_bar
       by(k-1,:)=by(k-1,:)+scale*rhs_bar
       do j=first,k-1
          bf(j,:)=bf(j,:)+weight(j-first+1)*rhs_bar
       enddo
    enddo
  end subroutine

  subroutine prepare_radial_pair(ctx,rel,potential,r0,dx,energy,ell,charge,quadrature)
    type(radial_pair_context), intent(out) :: ctx
    logical, intent(in) :: rel
    real(dp), intent(in) :: potential(:),r0,dx,energy,ell,charge
    integer, intent(in), optional :: quadrature
    call prepare_radial(ctx%center,rel,potential,r0,dx,energy,ell,charge,quadrature)
    call prepare_radial(ctx%plus,rel,potential,r0,dx,energy+0.002_dp,ell,charge,quadrature)
    call prepare_radial(ctx%minus,rel,potential,r0,dx,energy-0.002_dp,ell,charge,quadrature)
    ctx%u=radial_value(ctx%center,.true.)
    ctx%raw_dot=ctx%stencil*(radial_value(ctx%plus,.true.)-radial_value(ctx%minus,.true.))
    ctx%cross=sum(ctx%center%metric*ctx%u*ctx%raw_dot)
  end subroutine

  function radial_pair_value(ctx) result(value)
    type(radial_pair_context), intent(in) :: ctx
    real(dp) :: value(ctx%center%n,4,2)
    value(:,:,1)=ctx%u
    value(:,:,2)=ctx%raw_dot-ctx%cross*ctx%u
  end function

  subroutine radial_pair_action(ctx,dv,de,derivative)
    type(radial_pair_context), intent(in) :: ctx
    real(dp), intent(in) :: dv(ctx%center%n),de
    real(dp), intent(out) :: derivative(ctx%center%n,4,2)
    real(dp) :: du(ctx%center%n,4),dz(ctx%center%n,4),temp(ctx%center%n,4),dcross
    call radial_action(ctx%center,dv,de,du,.true.)
    call radial_action(ctx%plus,dv,de,dz,.true.)
    call radial_action(ctx%minus,dv,de,temp,.true.)
    dz=ctx%stencil*(dz-temp)
    dcross=sum(ctx%center%metric*(du*ctx%raw_dot+ctx%u*dz))
    derivative(:,:,1)=du
    derivative(:,:,2)=dz-dcross*ctx%u-ctx%cross*du
  end subroutine

  subroutine radial_pair_pullback(ctx,covector,dv_bar,de_bar)
    type(radial_pair_context), intent(in) :: ctx
    real(dp), intent(in) :: covector(ctx%center%n,4,2)
    real(dp), intent(out) :: dv_bar(ctx%center%n),de_bar
    real(dp) :: gu(ctx%center%n,4),gz(ctx%center%n,4),temp(ctx%center%n),te,s
    s=sum(covector(:,:,2)*ctx%u)
    gu=covector(:,:,1)-ctx%cross*covector(:,:,2)-s*ctx%center%metric*ctx%raw_dot
    gz=covector(:,:,2)-s*ctx%center%metric*ctx%u
    call radial_pullback(ctx%center,gu,dv_bar,de_bar,.true.)
    call radial_pullback(ctx%plus,ctx%stencil*gz,temp,te,.true.)
    dv_bar=dv_bar+temp; de_bar=de_bar+te
    call radial_pullback(ctx%minus,-ctx%stencil*gz,temp,te,.true.)
    dv_bar=dv_bar+temp; de_bar=de_bar+te
  end subroutine

  subroutine local_orbital_coefficients(lapw,jlo,x,coefficient,jacobian)
    ! Differentiate the actual ABC/ABC0 formula, INCLUDING its CLO cap.
    ! x=(p,dp,pe,dpe,pei,plo,dplo,pi12lo,pe12lo,RMT).
    ! Jacobian is tiny (3 by 10), so store it directly. Its transpose is the
    ! exact reverse action. All radial overlaps and boundary values must be
    ! differentiated by the caller with the same radial/matching convention.
    !
    ! In particular, the LAPW branch follows ABC's RMT**2 factors; replacing
    ! these by an exact boundary-matrix inverse would change the underlying
    ! basis algorithm. This routine does not make that unvalidated change.
    logical,intent(in) :: lapw
    integer,intent(in) :: jlo
    real(dp),intent(in) :: x(10)
    real(dp),intent(out) :: coefficient(3),jacobian(3,10)
    real(dp) :: eye(10,10),a,b,c,norm,da(10),db(10),dc(10),dn(10),r
    integer :: i
    eye=0
    do i=1,10
       eye(i,i)=1
    enddo
    coefficient=0;jacobian=0;r=x(10)
    if (lapw) then
       a=(x(6)*x(4)-x(7)*x(3))*r*r
       b=-(x(6)*x(2)-x(7)*x(1))*r*r
       da=(x(4)*eye(6,:)+x(6)*eye(4,:)-x(3)*eye(7,:)-x(7)*eye(3,:))*r*r &
          +2*r*(x(6)*x(4)-x(7)*x(3))*eye(10,:)
       db=-(x(2)*eye(6,:)+x(6)*eye(2,:)-x(1)*eye(7,:)-x(7)*eye(1,:))*r*r &
          -2*r*(x(6)*x(2)-x(7)*x(1))*eye(10,:)
       norm=1+a*(a+2*x(8))+b*(b*x(5)+2*x(9))
       if(norm<=0) error stop 'nonpositive LO normalization'
       dn=(2*a+2*x(8))*da+2*a*eye(8,:) &
          +(2*b*x(5)+2*x(9))*db+b*b*eye(5,:)+2*b*eye(9,:)
       c=1/sqrt(norm)
       if(abs(c-200.0_dp)<1e-8_dp) error stop 'LO derivative at nondifferentiable CLO cap'
       if(c>200) then
          c=200;dc=0
       else
          dc=-0.5_dp*c**3*dn
       endif
       coefficient=[c*a,c*b,c]
       jacobian(1,:)=c*da+a*dc
       jacobian(2,:)=c*db+b*dc
       jacobian(3,:)=dc
    else if(jlo==1) then
       if(abs(x(3))<tiny(1.0_dp)) error stop 'zero udot boundary value in APW+lo'
       b=-x(1)/x(3)
       db=-eye(1,:)/x(3)+x(1)*eye(3,:)/x(3)**2
       norm=1+b*b*x(5)
       if(norm<=0) error stop 'nonpositive APW+lo normalization'
       dn=2*b*x(5)*db+b*b*eye(5,:)
       a=1/sqrt(norm);da=-0.5_dp*a**3*dn
       coefficient=[a,a*b,0.0_dp]
       jacobian(1,:)=da;jacobian(2,:)=a*db+b*da
    else
       if(abs(x(6))<tiny(1.0_dp)) error stop 'zero LO boundary value in APW+LO'
       b=-x(1)/x(6)
       db=-eye(1,:)/x(6)+x(1)*eye(6,:)/x(6)**2
       norm=1+b*b+2*b*x(8)
       if(norm<=0) error stop 'nonpositive APW+LO normalization'
       dn=(2*b+2*x(8))*db+2*b*eye(8,:)
       a=1/sqrt(norm);da=-0.5_dp*a**3*dn
       coefficient=[a,0.0_dp,a*b]
       jacobian(1,:)=da;jacobian(3,:)=a*db+b*da
    endif
  end subroutine
end module outwin_response
