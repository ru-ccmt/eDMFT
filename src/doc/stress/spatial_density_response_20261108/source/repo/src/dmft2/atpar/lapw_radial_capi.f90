module lapw_radial_capi
  ! In-memory C interface for the native radial response. A context is prepared
  ! ONCE per atom/l, then reused for many forward and reverse operations.
  ! Potentials are r*V in Ry; energies are Ry. This calls the same discrete
  ! implementation already tested against ATPARN, not a new radial solver.
  !
  ! Handles belong to this process. Calls to this small registry must be
  ! serialized; the Python adapter supplies a lock. MPI processes have their
  ! own independent registries. Ordinary ATPARN/dmft2 does not call this API.
  use iso_c_binding, only: c_int,c_double
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use lapw_radial_response
  implicit none
  private
  integer,parameter :: max_contexts=512
  type slot
     type(radial_basis_context),pointer :: context=>null()
  end type
  type(slot),save :: contexts(max_contexts)
contains
  integer(c_int) function api_version() bind(C,name='edmft_radial_api_version')
    api_version=1
  end function

  logical function valid(handle)
    integer(c_int),intent(in) :: handle
    valid=.false.
    if(handle<1.or.handle>max_contexts)return
    valid=associated(contexts(handle)%context)
  end function

  integer(c_int) function create_context(n,nlo,rel,ell,lapw,r0,dx,rmt,energy,charge,v,elo,handle) &
       bind(C,name='edmft_radial_create') result(status)
    integer(c_int),value :: n,nlo,rel,ell,lapw
    real(c_double),value :: r0,dx,rmt,energy,charge
    real(c_double),intent(in) :: v(n),elo(nlo)
    integer(c_int),intent(out) :: handle
    integer :: i
    status=2;handle=0
    if(n<6.or.nlo<0.or.ell<0.or.(rel/=0.and.rel/=1).or.(lapw/=0.and.lapw/=1))return
    if(.not.all(ieee_is_finite([r0,dx,rmt,energy,charge])))return
    if(.not.all(ieee_is_finite(v)).or..not.all(ieee_is_finite(elo)))return
    if(r0<=0.or.dx<=0.or.rmt<=r0.or.charge<=0)return
    if(abs(log(rmt/r0)-dx*(n-1))>1e-10_c_double)return
    status=3
    do i=1,max_contexts
       if(associated(contexts(i)%context))cycle
       allocate(contexts(i)%context)
       call prepare_radial_basis(contexts(i)%context,rel==1,v,r0,dx,rmt,energy,elo, &
            lapw==1,real(ell,c_double),charge)
       handle=i;status=0
       return
    enddo
  end function

  integer(c_int) function release_context(handle) bind(C,name='edmft_radial_destroy') result(status)
    integer(c_int),value :: handle
    status=1
    if(.not.valid(handle))return
    deallocate(contexts(handle)%context)
    nullify(contexts(handle)%context)
    status=0
  end function

  integer(c_int) function get_value(handle,primitive,matrix) bind(C,name='edmft_radial_value') result(status)
    integer(c_int),value :: handle
    real(c_double),intent(out) :: primitive(*),matrix(*)
    type(radial_basis_context),pointer :: p
    real(c_double) :: boundary(2,0)
    integer :: n,m
    status=1
    if(.not.valid(handle))return
    p=>contexts(handle)%context;n=p%n;m=2+p%nlo
    primitive(1:n*4*m)=reshape(p%primitive,[n*4*m])
    matrix(1:m*p%nlo)=reshape(radial_component_coefficients(p,boundary),[m*p%nlo])
    status=0
  end function

  integer(c_int) function forward(handle,dv,de,delo,primitive_dot,matrix_dot) &
       bind(C,name='edmft_radial_action') result(status)
    integer(c_int),value :: handle
    real(c_double),value :: de
    real(c_double),intent(in) :: dv(*),delo(*)
    real(c_double),intent(out) :: primitive_dot(*),matrix_dot(*)
    type(radial_basis_context),pointer :: p
    real(c_double) :: boundary(2,0)
    real(c_double),allocatable :: dp(:,:,:),dm(:,:)
    integer :: n,m
    status=1
    if(.not.valid(handle))return
    p=>contexts(handle)%context;n=p%n;m=2+p%nlo
    status=2
    if(.not.ieee_is_finite(de).or..not.all(ieee_is_finite(dv(1:n))).or. &
         .not.all(ieee_is_finite(delo(1:p%nlo))))return
    allocate(dp(n,4,m),dm(m,p%nlo))
    call radial_components_action(p,boundary,dv(1:n),de,delo(1:p%nlo),dp,dm)
    primitive_dot(1:n*4*m)=reshape(dp,[n*4*m])
    matrix_dot(1:m*p%nlo)=reshape(dm,[m*p%nlo])
    status=0
  end function

  integer(c_int) function reverse(handle,primitive_bar,matrix_bar,dv_bar,de_bar,delo_bar) &
       bind(C,name='edmft_radial_pullback') result(status)
    integer(c_int),value :: handle
    real(c_double),intent(in) :: primitive_bar(*),matrix_bar(*)
    real(c_double),intent(out) :: dv_bar(*),de_bar,delo_bar(*)
    type(radial_basis_context),pointer :: p
    real(c_double) :: boundary(2,0)
    integer :: n,m
    status=1
    if(.not.valid(handle))return
    p=>contexts(handle)%context;n=p%n;m=2+p%nlo
    status=2
    if(.not.all(ieee_is_finite(primitive_bar(1:n*4*m))).or. &
         .not.all(ieee_is_finite(matrix_bar(1:m*p%nlo))))return
    call radial_components_pullback(p,boundary,reshape(primitive_bar(1:n*4*m),[n,4,m]), &
         reshape(matrix_bar(1:m*p%nlo),[m,p%nlo]),dv_bar(1:n),de_bar,delo_bar(1:p%nlo))
    status=0
  end function
end module lapw_radial_capi
