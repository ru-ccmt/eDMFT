module lapw_radial_response
  ! Assemble the discrete radial tangent into the ACTUAL LAPW/APW+lo basis.
  ! One context represents one atom/l channel. Complex angular/phase factors
  ! and KS coefficients are external linear maps; this module handles their
  ! real radial factors and can be pulled back through those maps unchanged.
  !
  ! The LAPW matching amplitudes are those used in ALM_BLM_CLM, including the
  ! RMT**2 prefactor. Do not replace their derivative by a boundary repair
  ! using an exact Wronskian inverse: ATPARN's finite udot stencil and RELA
  ! convention need not obey the continuous Wronskian identity exactly.
  !
  ! Input potentials are r*V_spherical in Ry, energies are Ry. The radial
  ! grid, RMT, ell, Z, and interstitial matching data are fixed. Ordinary LOs
  ! and all three ABC branches are included. DIRACOUT p1/2 LOs are rejected;
  ! their separate tangent has not been implemented. Energy derivatives are
  ! supplied by the caller, not inferred from an automatic energy search.
  use iso_fortran_env, only: real64
  use outwin_response
  implicit none
  private
  integer,parameter :: dp=real64
  public :: radial_basis_context, prepare_radial_basis
  public :: radial_basis_value, radial_basis_action, radial_basis_pullback
  public :: radial_basis_forms, radial_basis_forms_action, radial_basis_forms_pullback
  public :: radial_component_coefficients, radial_components_action, radial_components_pullback

  type :: radial_basis_context
     type(radial_pair_context) :: pair
     type(radial_context),allocatable :: local(:)
     integer :: n=0,nlo=0
     logical :: lapw=.true.,rel=.false.
     real(dp) :: rmt=0
     real(dp),allocatable :: primitive(:,:,:),coefficient(:,:),abc_jac(:,:,:)
     real(dp),allocatable :: overlap_metric(:,:),projector_metric(:,:)
  end type

contains

  subroutine prepare_radial_basis(ctx,rel,potential,r0,dx,rmt,energy,lo_energy,lapw,ell,charge,dirac_lo,quadrature)
    type(radial_basis_context),intent(out) :: ctx
    logical,intent(in) :: rel,lapw
    logical,intent(in),optional :: dirac_lo(:)
    integer,intent(in),optional :: quadrature
    real(dp),intent(in) :: potential(:),r0,dx,rmt,energy,lo_energy(:),ell,charge
    real(dp) :: x(10)
    integer :: j
    if(present(dirac_lo)) then
       if(size(dirac_lo)/=size(lo_energy)) error stop 'DIRACOUT flag dimensions disagree'
       if(any(dirac_lo)) error stop 'DIRACOUT LO radial response is not implemented'
    endif
    if(rmt<=0) error stop 'invalid RMT in radial basis response'
    ctx%n=size(potential);ctx%nlo=size(lo_energy);ctx%rmt=rmt
    ctx%lapw=lapw;ctx%rel=rel
    call prepare_radial_pair(ctx%pair,rel,potential/2,r0,dx,energy/2,ell,charge,quadrature)
    allocate(ctx%local(ctx%nlo),ctx%primitive(ctx%n,4,2+ctx%nlo))
    allocate(ctx%coefficient(3,ctx%nlo),ctx%abc_jac(3,10,ctx%nlo))
    ctx%primitive=0
    ctx%primitive(:,:,1:2)=radial_pair_value(ctx%pair)
    ctx%overlap_metric=ctx%pair%center%metric
    ctx%projector_metric=ctx%overlap_metric
    ! p_cmp_rixmat uses Rint13n: unlike ATPARN's RINT13 it omits the NREL
    ! small component altogether. Preserve that convention in BOTH actions.
    if(.not.rel) ctx%projector_metric(:,2)=0
    do j=1,ctx%nlo
       if(lapw.or.j>1) then
          call prepare_radial(ctx%local(j),rel,potential/2,r0,dx,lo_energy(j)/2,ell,charge,quadrature)
          ctx%primitive(:,:,2+j)=radial_value(ctx%local(j),.true.)
       endif
       call abc_inputs(ctx,j,x)
       call local_orbital_coefficients(lapw,j,x,ctx%coefficient(:,j),ctx%abc_jac(:,:,j))
    enddo
  end subroutine

  subroutine abc_inputs(ctx,j,x)
    type(radial_basis_context),intent(in) :: ctx
    integer,intent(in) :: j
    real(dp),intent(out) :: x(10)
    integer :: n
    n=ctx%n
    x(1:4)=[ctx%primitive(n,1,1),ctx%primitive(n,3,1), &
            ctx%primitive(n,1,2),ctx%primitive(n,3,2)]/ctx%rmt
    x(5)=sum(ctx%overlap_metric*ctx%primitive(:,:,2)**2)
    x(6:7)=[ctx%primitive(n,1,2+j),ctx%primitive(n,3,2+j)]/ctx%rmt
    x(8)=sum(ctx%overlap_metric*ctx%primitive(:,:,1)*ctx%primitive(:,:,2+j))
    x(9)=sum(ctx%overlap_metric*ctx%primitive(:,:,2)*ctx%primitive(:,:,2+j))
    x(10)=ctx%rmt
  end subroutine

  subroutine matching(ctx,boundary,amplitude)
    type(radial_basis_context),intent(in) :: ctx
    real(dp),intent(in) :: boundary(:,:)
    real(dp),intent(out) :: amplitude(2,size(boundary,2))
    real(dp) :: p,slope,pe,slope_e,r
    r=ctx%rmt
    p=ctx%primitive(ctx%n,1,1)/r;slope=ctx%primitive(ctx%n,3,1)/r
    pe=ctx%primitive(ctx%n,1,2)/r;slope_e=ctx%primitive(ctx%n,3,2)/r
    if(size(boundary,1)/=2) error stop 'boundary data must be (value,radial derivative)'
    if(ctx%lapw) then
       amplitude(1,:)=r*r*(boundary(2,:)*pe-boundary(1,:)*slope_e)
       amplitude(2,:)=r*r*(boundary(1,:)*slope-boundary(2,:)*p)
    else
       if(abs(p)<tiny(1.0_dp)) error stop 'zero APW matching value'
       amplitude(1,:)=boundary(1,:)/p;amplitude(2,:)=0
    endif
  end subroutine

  function radial_basis_value(ctx,boundary) result(value)
    type(radial_basis_context),intent(in) :: ctx
    real(dp),intent(in) :: boundary(:,:)
    real(dp) :: value(ctx%n,4,size(boundary,2)+ctx%nlo),amplitude(2,size(boundary,2))
    integer :: j,np
    call matching(ctx,boundary,amplitude)
    np=size(boundary,2)
    do j=1,np
       value(:,:,j)=amplitude(1,j)*ctx%primitive(:,:,1)+amplitude(2,j)*ctx%primitive(:,:,2)
    enddo
    do j=1,ctx%nlo
       value(:,:,np+j)=ctx%coefficient(1,j)*ctx%primitive(:,:,1) &
            +ctx%coefficient(2,j)*ctx%primitive(:,:,2)+ctx%coefficient(3,j)*ctx%primitive(:,:,2+j)
    enddo
  end function

  function radial_component_coefficients(ctx,boundary) result(matrix)
    ! Factor each complete basis function as f_b = sum_i primitive_i M_ib.
    ! Exposing M and its paired derivative lets H, O, density and projectors
    ! use the SAME matching/LO chain, with no inverse or inferred connection.
    type(radial_basis_context),intent(in) :: ctx
    real(dp),intent(in) :: boundary(:,:)
    real(dp) :: matrix(2+ctx%nlo,size(boundary,2)+ctx%nlo)
    integer :: j,np
    np=size(boundary,2);matrix=0
    call matching(ctx,boundary,matrix(1:2,:np))
    do j=1,ctx%nlo
       matrix(1:2,np+j)=ctx%coefficient(1:2,j)
       matrix(2+j,np+j)=ctx%coefficient(3,j)
    enddo
  end function

  subroutine radial_components_action(ctx,boundary,dv,de,delo,primitive_dot,matrix_dot)
    type(radial_basis_context),intent(in) :: ctx
    real(dp),intent(in) :: boundary(:,:),dv(ctx%n),de,delo(ctx%nlo)
    real(dp),intent(out) :: primitive_dot(ctx%n,4,2+ctx%nlo)
    real(dp),intent(out) :: matrix_dot(2+ctx%nlo,size(boundary,2)+ctx%nlo)
    real(dp) :: amplitude(2,size(boundary,2)),da(2,size(boundary,2)),xdot(10),cdot(3),r,p
    integer :: j,n,np
    n=ctx%n;np=size(boundary,2);r=ctx%rmt
    primitive_dot=0
    call radial_pair_action(ctx%pair,dv/2,de/2,primitive_dot(:,:,1:2))
    do j=1,ctx%nlo
       if(ctx%lapw.or.j>1) call radial_action(ctx%local(j),dv/2,delo(j)/2,primitive_dot(:,:,2+j),.true.)
    enddo
    call matching(ctx,boundary,amplitude)
    if(ctx%lapw) then
       da(1,:)=r*(boundary(2,:)*primitive_dot(n,1,2)-boundary(1,:)*primitive_dot(n,3,2))
       da(2,:)=r*(boundary(1,:)*primitive_dot(n,3,1)-boundary(2,:)*primitive_dot(n,1,1))
    else
       p=ctx%primitive(n,1,1)/r
       da(1,:)=-boundary(1,:)*primitive_dot(n,1,1)/(r*p*p);da(2,:)=0
    endif
    matrix_dot=0;matrix_dot(1:2,:np)=da
    do j=1,ctx%nlo
       xdot=0
       xdot(1:4)=[primitive_dot(n,1,1),primitive_dot(n,3,1), &
                  primitive_dot(n,1,2),primitive_dot(n,3,2)]/r
       xdot(5)=2*sum(ctx%overlap_metric*ctx%primitive(:,:,2)*primitive_dot(:,:,2))
       xdot(6:7)=[primitive_dot(n,1,2+j),primitive_dot(n,3,2+j)]/r
       xdot(8)=sum(ctx%overlap_metric*(primitive_dot(:,:,1)*ctx%primitive(:,:,2+j) &
                                     +ctx%primitive(:,:,1)*primitive_dot(:,:,2+j)))
       xdot(9)=sum(ctx%overlap_metric*(primitive_dot(:,:,2)*ctx%primitive(:,:,2+j) &
                                     +ctx%primitive(:,:,2)*primitive_dot(:,:,2+j)))
       cdot=matmul(ctx%abc_jac(:,:,j),xdot)
       matrix_dot(1:2,np+j)=cdot(1:2)
       matrix_dot(2+j,np+j)=cdot(3)
    enddo
  end subroutine

  subroutine radial_components_pullback(ctx,boundary,primitive_bar,matrix_bar,dv_bar,de_bar,delo_bar)
    type(radial_basis_context),intent(in) :: ctx
    real(dp),intent(in) :: boundary(:,:),primitive_bar(ctx%n,4,2+ctx%nlo)
    real(dp),intent(in) :: matrix_bar(2+ctx%nlo,size(boundary,2)+ctx%nlo)
    real(dp),intent(out) :: dv_bar(ctx%n),de_bar,delo_bar(ctx%nlo)
    real(dp) :: gp(ctx%n,4,2+ctx%nlo),ga,gb,gc(3),gx(10),tmp(ctx%n),te,r,p
    integer :: j,n,np
    gp=primitive_bar;n=ctx%n;np=size(boundary,2);r=ctx%rmt
    do j=1,np
       ga=matrix_bar(1,j);gb=matrix_bar(2,j)
       if(ctx%lapw) then
          gp(n,1,2)=gp(n,1,2)+ga*r*boundary(2,j)
          gp(n,3,2)=gp(n,3,2)-ga*r*boundary(1,j)
          gp(n,3,1)=gp(n,3,1)+gb*r*boundary(1,j)
          gp(n,1,1)=gp(n,1,1)-gb*r*boundary(2,j)
       else
          p=ctx%primitive(n,1,1)/r
          gp(n,1,1)=gp(n,1,1)-ga*boundary(1,j)/(r*p*p)
       endif
    enddo
    do j=1,ctx%nlo
       gc=[matrix_bar(1,np+j),matrix_bar(2,np+j),matrix_bar(2+j,np+j)]
       gx=matmul(transpose(ctx%abc_jac(:,:,j)),gc)
       gp(n,1,1)=gp(n,1,1)+gx(1)/r;gp(n,3,1)=gp(n,3,1)+gx(2)/r
       gp(n,1,2)=gp(n,1,2)+gx(3)/r;gp(n,3,2)=gp(n,3,2)+gx(4)/r
       gp(:,:,2)=gp(:,:,2)+2*gx(5)*ctx%overlap_metric*ctx%primitive(:,:,2)
       gp(n,1,2+j)=gp(n,1,2+j)+gx(6)/r;gp(n,3,2+j)=gp(n,3,2+j)+gx(7)/r
       gp(:,:,1)=gp(:,:,1)+gx(8)*ctx%overlap_metric*ctx%primitive(:,:,2+j)
       gp(:,:,2)=gp(:,:,2)+gx(9)*ctx%overlap_metric*ctx%primitive(:,:,2+j)
       gp(:,:,2+j)=gp(:,:,2+j)+ctx%overlap_metric*(gx(8)*ctx%primitive(:,:,1)+gx(9)*ctx%primitive(:,:,2))
    enddo
    call radial_pair_pullback(ctx%pair,gp(:,:,1:2),dv_bar,de_bar)
    delo_bar=0
    do j=1,ctx%nlo
       if(.not.ctx%lapw.and.j==1) cycle
       call radial_pullback(ctx%local(j),gp(:,:,2+j),tmp,te,.true.)
       dv_bar=dv_bar+tmp;delo_bar(j)=te/2
    enddo
    dv_bar=dv_bar/2;de_bar=de_bar/2
  end subroutine


  subroutine radial_basis_action(ctx,boundary,dv,de,delo,derivative)
    type(radial_basis_context),intent(in) :: ctx
    real(dp),intent(in) :: boundary(:,:),dv(ctx%n),de,delo(ctx%nlo)
    real(dp),intent(out) :: derivative(ctx%n,4,size(boundary,2)+ctx%nlo)
    real(dp) :: primitive_dot(ctx%n,4,2+ctx%nlo)
    real(dp) :: matrix(2+ctx%nlo,size(boundary,2)+ctx%nlo),matrix_dot(2+ctx%nlo,size(boundary,2)+ctx%nlo)
    integer :: i,b
    matrix=radial_component_coefficients(ctx,boundary)
    call radial_components_action(ctx,boundary,dv,de,delo,primitive_dot,matrix_dot)
    derivative=0
    do b=1,size(matrix,2)
       do i=1,size(matrix,1)
          derivative(:,:,b)=derivative(:,:,b)+matrix(i,b)*primitive_dot(:,:,i)+matrix_dot(i,b)*ctx%primitive(:,:,i)
       enddo
    enddo
  end subroutine

  subroutine radial_basis_pullback(ctx,boundary,covector,dv_bar,de_bar,delo_bar)
    type(radial_basis_context),intent(in) :: ctx
    real(dp),intent(in) :: boundary(:,:),covector(ctx%n,4,size(boundary,2)+ctx%nlo)
    real(dp),intent(out) :: dv_bar(ctx%n),de_bar,delo_bar(ctx%nlo)
    real(dp) :: gp(ctx%n,4,2+ctx%nlo)
    real(dp) :: matrix(2+ctx%nlo,size(boundary,2)+ctx%nlo),gm(2+ctx%nlo,size(boundary,2)+ctx%nlo)
    integer :: i,b
    matrix=radial_component_coefficients(ctx,boundary);gp=0;gm=0
    do b=1,size(matrix,2)
       do i=1,size(matrix,1)
          gp(:,:,i)=gp(:,:,i)+matrix(i,b)*covector(:,:,b)
          gm(i,b)=sum(ctx%primitive(:,:,i)*covector(:,:,b))
       enddo
    enddo
    call radial_components_pullback(ctx,boundary,gp,gm,dv_bar,de_bar,delo_bar)
  end subroutine

  subroutine radial_basis_forms(ctx,value,projector,overlap,projection)
    ! Radial overlap and raw projector-5 factors. The finite projector may be
    ! nonzero at RMT; no boundary value is discarded. Angular, Bloch, window
    ! selection, k-sum, and normalization are later operations.
    type(radial_basis_context),intent(in) :: ctx
    real(dp),intent(in) :: value(:,:,:),projector(ctx%n,2)
    real(dp),intent(out) :: overlap(size(value,3),size(value,3)),projection(size(value,3))
    integer :: i,j
    do i=1,size(value,3)
       projection(i)=sum(ctx%projector_metric(:,1:2)*value(:,1:2,i)*projector)
       do j=1,size(value,3)
          overlap(i,j)=sum(ctx%overlap_metric*value(:,:,i)*value(:,:,j))
       enddo
    enddo
  end subroutine

  subroutine radial_basis_forms_action(ctx,value,derivative,projector,overlap_dot,projection_dot)
    type(radial_basis_context),intent(in) :: ctx
    real(dp),intent(in) :: value(:,:,:),derivative(:,:,:),projector(ctx%n,2)
    real(dp),intent(out) :: overlap_dot(size(value,3),size(value,3)),projection_dot(size(value,3))
    integer :: i,j
    do i=1,size(value,3)
       projection_dot(i)=sum(ctx%projector_metric(:,1:2)*derivative(:,1:2,i)*projector)
       do j=1,size(value,3)
          overlap_dot(i,j)=sum(ctx%overlap_metric*(derivative(:,:,i)*value(:,:,j)+value(:,:,i)*derivative(:,:,j)))
       enddo
    enddo
  end subroutine

  subroutine radial_basis_forms_pullback(ctx,value,projector,overlap_bar,projection_bar,basis_bar)
    type(radial_basis_context),intent(in) :: ctx
    real(dp),intent(in) :: value(:,:,:),projector(ctx%n,2)
    real(dp),intent(in) :: overlap_bar(size(value,3),size(value,3)),projection_bar(size(value,3))
    real(dp),intent(out) :: basis_bar(size(value,1),size(value,2),size(value,3))
    integer :: i,j
    basis_bar=0
    do i=1,size(value,3)
       basis_bar(:,1:2,i)=projection_bar(i)*ctx%projector_metric(:,1:2)*projector
       do j=1,size(value,3)
          basis_bar(:,:,i)=basis_bar(:,:,i)+(overlap_bar(i,j)+overlap_bar(j,i))*ctx%overlap_metric*value(:,:,j)
       enddo
    enddo
  end subroutine
end module lapw_radial_response
