! Optional metadata written by ssplit on the static (last) self-energy record.
! The fit samples are selected before compression and are identical at every
! quadrature tolerance. The coarser check reuses the fine Green functions.
module adaptive_matsubara
  implicit none
  private
  public :: read_mesh_header, prepare_coarse_check, release_mesh
  public :: low_count, fit_rows, weight_difference, moment, check_sum
  integer :: low_count=0
  integer, allocatable :: fit_rows(:)
  real*8, allocatable :: weight_difference(:)
  complex*16, allocatable :: moment(:,:)
  real*8 :: check_sum(3)=0.d0
contains
  subroutine read_mesh_header(unit,ncorr,nomega,first)
    integer,intent(in) :: unit,ncorr,nomega
    logical,intent(in) :: first
    character(len=32) :: tag
    real*8 :: numeric(1+2*ncorr)
    integer :: ios,nlow,nfit,indices(24)
    logical :: tagged

    backspace(unit)
    tag=''
    read(unit,*,iostat=ios) numeric,tag
    tagged=ios==0 .and. tag=='#EDMFT_AUTO_TAIL_V1'
    if (.not.first .and. (tagged .neqv. allocated(fit_rows))) &
         error stop 'mixed automatic and conventional self-energy inputs'
    if (.not.tagged) return
    backspace(unit)
    read(unit,*,iostat=ios) numeric,tag,nlow,nfit
    if (ios/=0) error stop 'invalid automatic Matsubara header'
    if (nfit<1 .or. nfit>24 .or. nlow<1 .or. nlow>nomega) error stop 'invalid automatic Matsubara dimensions'
    backspace(unit)
    read(unit,*,iostat=ios) numeric,tag,nlow,nfit,indices(1:nfit)
    if (ios/=0) error stop 'invalid Matsubara moment-fit rows'
    if (indices(1)<1 .or. indices(nfit)/=nomega .or. &
         any(indices(2:nfit)<=indices(1:nfit-1))) error stop 'invalid Matsubara moment-fit ordering'
    if (first) then
       low_count=nlow
       fit_rows=indices(:nfit)
    else
       if (low_count/=nlow .or. size(fit_rows)/=nfit) error stop 'inconsistent Matsubara headers'
       if (any(fit_rows/=indices(:nfit))) error stop 'inconsistent Matsubara moment-fit rows'
    endif
  end subroutine

  subroutine prepare_coarse_check(fine,indices)
    use splines, only: spline3
    real*8,intent(in) :: fine(:)
    integer,intent(in) :: indices(:)
    real*8,allocatable :: x(:),targets(:),basis(:)
    integer,allocatable :: nodes(:)
    integer :: n,i,j,start
    if (.not.allocated(fit_rows)) return
    start=low_count+1
    n=size(fine)
    ! Four distinct coarse knots are required by the production spline.
    if (n-start+1<7) return
    nodes=[(i,i=start,n,2)]
    if (nodes(size(nodes))/=n) nodes=[nodes,n]
    nodes=nodes(size(nodes):1:-1)
    x=1.d0/(2.d0*indices(nodes)-1.d0)
    targets=[(1.d0/(2.d0*i-1.d0),i=indices(n),indices(start),-1)]
    allocate(basis(size(x)))
    weight_difference=fine
    weight_difference(:low_count)=0.d0
    do i=1,size(nodes)
       basis=0.d0
       basis(i)=1.d0
       j=nodes(i)
       weight_difference(j)=fine(j)-sum(spline3(x,basis,targets))
    enddo
  end subroutine

  subroutine release_mesh()
    if (allocated(fit_rows)) deallocate(fit_rows)
    if (allocated(moment)) deallocate(moment)
    if (allocated(weight_difference)) deallocate(weight_difference)
    low_count=0
    check_sum=0.d0
  end subroutine
end module
