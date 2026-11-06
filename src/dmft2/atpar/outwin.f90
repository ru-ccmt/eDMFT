subroutine outwin(a, b, ap, bp, nodes, val, slo, rel, v, rnot, dh, jri, eh, fl, z, nrad)
  ! Integrate the scalar-relativistic radial equations on a logarithmic mesh.
  !
  ! The propagated quantities A and B are r times the large and (unscaled)
  ! small radial components.  The independent variable is the radial-mesh
  ! index, so the derivative history below contains
  !
  !       dA/di = (dh*r) dA/dr,       dB/di = (dh*r) dB/dr .
  !
  ! We use an implicit Adams--Moulton scheme.  The first propagated point
  ! (k=4) uses fourth order because only three derivative values are then
  ! available.  Every subsequent point uses the fifth-order formula
  !
  !   y_k = y_{k-1} + (251 f_k + 646 f_{k-1} - 264 f_{k-2}
  !                    + 106 f_{k-3} - 19 f_{k-4}) / 720 .
  !
  ! Solving the coupled first-order radial equations implicitly at point k
  ! gives the compact 2x2 solve implemented below.  This fifth-order update
  ! matches the accuracy used by modern LAPW radial solvers while retaining
  ! dmft2's AP/BP arrays, which are required by the RELA stress calculation.
  !
  ! Units and conventions:
  !   * EH and V are supplied in Hartree; E=2*EH enters the Rydberg equation.
  !   * A and B are normalized later by ATPAR.
  !   * AP and BP are r*d(A/r)/dr and r*d(B/r)/dr, respectively.
  !   * B and BP are converted to the stored small-component convention only
  !     after the integration has finished.

  implicit none

  integer, parameter :: dp = selected_real_kind(15, 300)

  ! Fourth-order Adams--Moulton coefficients, divided by the implicit
  ! coefficient 9/24.  These are used only for the k=4 startup step.
  real(dp), parameter :: am4_scale = 8.0_dp / 3.0_dp
  real(dp), parameter :: am4_w1    = 1.0_dp / 9.0_dp
  real(dp), parameter :: am4_w2    = -5.0_dp / 9.0_dp
  real(dp), parameter :: am4_w3    = 19.0_dp / 9.0_dp

  ! Fifth-order Adams--Moulton coefficients, divided by the implicit
  ! coefficient 251/720.  Derivative history runs from oldest to newest.
  real(dp), parameter :: am5_scale = 720.0_dp / 251.0_dp
  real(dp), parameter :: am5_w0    = -19.0_dp / 251.0_dp
  real(dp), parameter :: am5_w1    = 106.0_dp / 251.0_dp
  real(dp), parameter :: am5_w2    = -264.0_dp / 251.0_dp
  real(dp), parameter :: am5_w3    = 646.0_dp / 251.0_dp

  real(dp), parameter :: twice_speed_of_light = 2.0_dp * 137.0359895_dp
  real(dp), parameter :: nonrelativistic_c     = 1.0e10_dp

  integer, intent(in) :: jri, nrad
  logical, intent(in) :: rel
  real(dp), intent(in) :: v(nrad)
  real(dp), intent(in) :: rnot, dh, eh, fl, z

  real(dp), intent(out) :: a(nrad), b(nrad)
  real(dp), intent(out) :: ap(nrad), bp(nrad)
  integer, intent(out) :: nodes
  real(dp), intent(out) :: val, slo

  integer :: i, k
  real(dp) :: radial_mesh(nrad)
  real(dp) :: d_large(0:3), d_small(0:3)
  real(dp) :: energy_ryd, c, ell_factor
  real(dp) :: exponent_large, exponent_small, small_amplitude
  real(dp) :: alpha_z, radius, dr_di
  real(dp) :: phi, coupling_u, coupling_x, coupling_y
  real(dp) :: step_scale, determinant, rhs_large, rhs_small
  real(dp) :: new_d_large, new_d_small

  if (jri < 4) error stop 'OUTWIN requires at least four radial mesh points'

  energy_ryd = 2.0_dp * eh
  c = twice_speed_of_light
  if (.not. rel) c = nonrelativistic_c

  a  = 0.0_dp
  b  = 0.0_dp
  ap = 0.0_dp
  bp = 0.0_dp
  d_large = 0.0_dp
  d_small = 0.0_dp
  nodes = 0

  do i = 1, jri
     radial_mesh(i) = rnot * exp(dh * real(i - 1, dp))
  end do

  ! Regular near-nucleus solution.  For a finite/weak nuclear charge this
  ! reduces smoothly to the nonrelativistic r**(l+1) behavior.
  ell_factor = fl * (fl + 1.0_dp)
  if (z < 0.9_dp) then
     exponent_large = fl + 1.0_dp
     exponent_small = fl
     small_amplitude = fl / c
  else
     alpha_z = 2.0_dp * z / c
     exponent_large = sqrt(ell_factor + 1.0_dp - alpha_z * alpha_z)
     exponent_small = exponent_large
     small_amplitude = (exponent_large - 1.0_dp) / alpha_z
  end if

  ! Analytic values at the first three mesh points seed the multistep method.
  do k = 1, 3
     radius = radial_mesh(k)
     dr_di = dh * radius

     a(k) = radius**exponent_large
     b(k) = radius**exponent_small * small_amplitude

     d_large(k) = dh * exponent_large * a(k)
     d_small(k) = dh * exponent_small * b(k)

     ap(k) = d_large(k) / dr_di - a(k) / radius
     bp(k) = d_small(k) / dr_di - b(k) / radius
  end do

  do k = 4, jri
     radius = radial_mesh(k)
     dr_di = dh * radius

     ! Coefficients of the coupled first-order Koelling--Harmon equations.
     ! V/r is the spherical potential in the convention used by ATPAR.
     phi = (energy_ryd - 2.0_dp * v(k) / radius) * dr_di / c
     coupling_u = dr_di * c + phi
     coupling_x = -dr_di / radius
     coupling_y = -ell_factor * coupling_x * coupling_x / coupling_u + phi

     if (k == 4) then
        ! A fifth-order multistep formula needs four previous derivatives.
        ! At k=4 only f_1, f_2, and f_3 exist, so start with fourth order.
        step_scale = am4_scale
        rhs_large = step_scale * a(k - 1) + am4_w1 * d_large(1) &
                    + am4_w2 * d_large(2) + am4_w3 * d_large(3)
        rhs_small = step_scale * b(k - 1) + am4_w1 * d_small(1) &
                    + am4_w2 * d_small(2) + am4_w3 * d_small(3)
     else
        step_scale = am5_scale
        rhs_large = step_scale * a(k - 1) + am5_w0 * d_large(0) &
                    + am5_w1 * d_large(1) + am5_w2 * d_large(2) &
                    + am5_w3 * d_large(3)
        rhs_small = step_scale * b(k - 1) + am5_w0 * d_small(0) &
                    + am5_w1 * d_small(1) + am5_w2 * d_small(2) &
                    + am5_w3 * d_small(3)
     end if

     ! Implicit solve for A_k and B_k.  Keeping this solve common to both
     ! integration orders makes the startup and production paths identical
     ! except for their history weights.
     determinant = step_scale * step_scale - coupling_x * coupling_x &
                   + coupling_u * coupling_y
     a(k) = (rhs_large * (step_scale - coupling_x) &
             + rhs_small * coupling_u) / determinant
     b(k) = (rhs_small * (step_scale + coupling_x) &
             - rhs_large * coupling_y) / determinant

     if (a(k) * a(k - 1) < 0.0_dp) nodes = nodes + 1

     new_d_large = coupling_u * b(k) - coupling_x * a(k)
     new_d_small = coupling_x * b(k) - coupling_y * a(k)

     d_large(0:2) = d_large(1:3)
     d_small(0:2) = d_small(1:3)
     d_large(3) = new_d_large
     d_small(3) = new_d_small

     ! These derivatives come directly from the first-order equations.
     ! Numerical differentiation here produces a visible error in RELA stress.
     ap(k) = new_d_large / dr_di - a(k) / radius
     bp(k) = new_d_small / dr_di - b(k) / radius
  end do

  ! Convert the raw small solution to the convention used by RINT13 and the
  ! remainder of dmft2.  Apply the same conversion to its radial derivative.
  b(1:jri)  = b(1:jri) * c / 2.0_dp
  bp(1:jri) = bp(1:jri) * c / 2.0_dp

  val = a(jri) / radial_mesh(jri)
  slo = ap(jri) / radial_mesh(jri)

end subroutine outwin
