module lapw_radial_potential
  ! Nonspherical potential integrals in the SMALL primitive radial space.
  !
  ! I(i,j,L) = sum_(r,c) W(r,c) phi_left(r,c,i) V_L(r) phi_right(r,c,j).
  ! The two channels may have different l and different numbers of LOs.
  ! W is the actual discrete RINT13 metric, including the origin term and
  ! RELA small component. V_L is the potential itself in Ry, NOT r*V_L.
  ! Angular Gaunt factors and matching/ABC coefficients are applied later.
  !
  ! The action includes changes on BOTH radial legs and in the potential.
  ! The reverse action returns covectors of these three independent inputs.
  ! Sum the radial covectors with spherical H/O, projector, and density
  ! contributions before ONE radial_components_pullback per atom/l.
  use iso_fortran_env, only: real64
  use lapw_radial_response, only: radial_basis_context
  implicit none
  private
  integer,parameter :: dp=real64
  public :: primitive_potential_value, primitive_potential_action, primitive_potential_pullback
contains
  subroutine compatible(left,right)
    type(radial_basis_context),intent(in) :: left,right
    if(left%n/=right%n.or.(left%rel.neqv.right%rel)) error stop 'incompatible radial potential contexts'
    if(maxval(abs(left%pair%center%r-right%pair%center%r))>1e-13_dp) &
         error stop 'radial potential contexts have different grids'
  end subroutine

  subroutine primitive_potential_value(left,right,potential,value)
    type(radial_basis_context),intent(in) :: left,right
    real(dp),intent(in) :: potential(:,:)
    real(dp),intent(out) :: value(2+left%nlo,2+right%nlo,size(potential,2))
    real(dp) :: product(left%n)
    integer :: i,j,c
    call compatible(left,right)
    if(size(potential,1)/=left%n) error stop 'radial potential length differs'
    do j=1,size(value,2)
       do i=1,size(value,1)
          product=0
          do c=1,2
             product=product+left%overlap_metric(:,c)*left%primitive(:,c,i)*right%primitive(:,c,j)
          enddo
          value(i,j,:)=matmul(product,potential)
       enddo
    enddo
  end subroutine

  subroutine primitive_potential_action(left,right,potential,left_dot,right_dot,potential_dot,value_dot)
    type(radial_basis_context),intent(in) :: left,right
    real(dp),intent(in) :: potential(:,:),left_dot(left%n,4,2+left%nlo),right_dot(right%n,4,2+right%nlo)
    real(dp),intent(in) :: potential_dot(size(potential,1),size(potential,2))
    real(dp),intent(out) :: value_dot(2+left%nlo,2+right%nlo,size(potential,2))
    real(dp) :: product(left%n),product_dot(left%n)
    integer :: i,j,c
    call compatible(left,right)
    if(size(potential,1)/=left%n) error stop 'radial potential length differs'
    do j=1,size(value_dot,2)
       do i=1,size(value_dot,1)
          product=0;product_dot=0
          do c=1,2
             product=product+left%overlap_metric(:,c)*left%primitive(:,c,i)*right%primitive(:,c,j)
             product_dot=product_dot+left%overlap_metric(:,c)*( &
                  left_dot(:,c,i)*right%primitive(:,c,j)+left%primitive(:,c,i)*right_dot(:,c,j))
          enddo
          value_dot(i,j,:)=matmul(product_dot,potential)+matmul(product,potential_dot)
       enddo
    enddo
  end subroutine

  subroutine primitive_potential_pullback(left,right,potential,value_bar,left_bar,right_bar,potential_bar)
    type(radial_basis_context),intent(in) :: left,right
    real(dp),intent(in) :: potential(:,:),value_bar(2+left%nlo,2+right%nlo,size(potential,2))
    real(dp),intent(out) :: left_bar(left%n,4,2+left%nlo),right_bar(right%n,4,2+right%nlo)
    real(dp),intent(out) :: potential_bar(size(potential,1),size(potential,2))
    real(dp) :: weighted(left%n),product(left%n)
    integer :: i,j,c,q
    call compatible(left,right)
    if(size(potential,1)/=left%n) error stop 'radial potential length differs'
    left_bar=0;right_bar=0;potential_bar=0
    do j=1,size(value_bar,2)
       do i=1,size(value_bar,1)
          weighted=matmul(potential,value_bar(i,j,:));product=0
          do c=1,2
             left_bar(:,c,i)=left_bar(:,c,i)+left%overlap_metric(:,c)*weighted*right%primitive(:,c,j)
             right_bar(:,c,j)=right_bar(:,c,j)+left%overlap_metric(:,c)*weighted*left%primitive(:,c,i)
             product=product+left%overlap_metric(:,c)*left%primitive(:,c,i)*right%primitive(:,c,j)
          enddo
          do q=1,size(potential,2)
             potential_bar(:,q)=potential_bar(:,q)+product*value_bar(i,j,q)
          enddo
       enddo
    enddo
  end subroutine
end module lapw_radial_potential
