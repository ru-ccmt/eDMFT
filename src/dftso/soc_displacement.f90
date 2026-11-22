! Explicit derivatives of the spherical spin-orbit matrix in second variation.
! At fixed scalar LAPW coefficients, a displacement of site a by t*u_a
! differentiates both matching legs of B^dagger R B. R is the radial SOC
! integral on the grid moving with a. Its limits remain fixed in that frame.
! At fixed R and site positions the alternative scalar_window mode varies
! both matching legs with Q_out*dA. Neither diagnostic includes radial
! relaxation. Each is an isolated matrix term in the free-energy derivative.
module soc_displacement
  use iso_fortran_env, only: real64, int32, int64
  use ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: init_soc_displacement, velocity, accumulate_soc_translation, write_soc_translation
  public :: scalar_window, read_scalar_tangent, export_soc_geometry
  public :: init_soc_matching, legacy_lo_filter
  logical, protected :: legacy_lo_filter = .false.
  logical :: scalar_window = .false.
  real(real64), allocatable :: velocity(:,:)
  character(1024), private :: prefix
  logical, private :: full_matrix
contains
  subroutine init_soc_matching(print_policy)
    logical, intent(in) :: print_policy
    character(16) :: setting
    integer :: status
    call get_environment_variable('EDMFT_SOC_LEGACY_LO_FILTER',setting,status=status)
    if (status/=0 .and. status/=1) error stop 'invalid EDMFT_SOC_LEGACY_LO_FILTER'
    select case (trim(setting))
    case ('','0')
       legacy_lo_filter=.false.
    case ('1')
       legacy_lo_filter=.true.
    case default
       error stop 'EDMFT_SOC_LEGACY_LO_FILTER must be 0 or 1'
    end select
    if (print_policy) then
       if (legacy_lo_filter) then
          write(8,'(a)') ':SOC_LO_MATCHING: LEGACY_COMPONENT_CUTOFF'
       else
          write(8,'(a)') ':SOC_LO_MATCHING: LINEAR'
       end if
    end if
  end subroutine

  ! Compact central data for the analytic B^H R B derivative. ROTIJ and
  ! the spin coupling must be those chosen by LAPWSO, including the gauge on
  ! equivalent sites. Reading them once avoids reconstructing the spin frame.
  subroutine export_soc_geometry(ri,e,jspin,nrel_lo,iorbpot)
    use param, only: labc,lmax,nato,nloat
    use structure, only: br1,vol,rotij,rotloc,tauij,pos,mult
    use couples, only: couplo
    use lolog, only: lso,mrf,lapw,ilo,count_lo=>nlo,nlov
    use loabc, only: elo
    integer, intent(in) :: jspin,nrel_lo,iorbpot
    real(real64), intent(in) :: ri(nloat,nloat,0:labc,nato,2,2),e(0:lmax,nato,2)
    character(1024) :: path
    integer :: unit,status,site,atom,l,nr,first,j,nlo
    call get_environment_variable('EDMFT_SOC_RESPONSE_GEOMETRY',path,status=status)
    if(status==1 .or. len_trim(path)==0) return
    if(status/=0) error stop 'SOC response geometry path is too long'
    if(jspin/=1.or.nrel_lo/=0.or.iorbpot/=0) error stop 'SOC response requires nonmagnetic scalar input without Dirac LO or orbital potential'
    open(newunit=unit,file=trim(path),status='replace',form='unformatted',access='stream',action='write')
    write(unit) 'SOCGEOM3',int([nato,sum(mult),labc,sum(count_lo),merge(1,0,legacy_lo_filter)],int32),vol,br1
    first=1
    do atom=1,nato
       write(unit) int([mult(atom),merge(1,0,lso(atom)),nlov(atom)],int32)
       do l=0,labc
          nr=mrf(l,atom)
          nlo=0
          if(l<=ubound(ilo,1)) nlo=ilo(l,atom)
          write(unit) int([nr,nlo,merge(1,0,lapw(l,atom))],int32),e(l,atom,1)
          if(nlo>0) write(unit) elo(l,1:nlo,atom,1)
          write(unit) ri(1:nr,1:nr,l,atom,1,1)
       end do
       do site=first,first+mult(atom)-1
          write(unit) matmul(rotloc(:,:,atom),matmul(br1,rotij(:,:,site))), &
            matmul(transpose(rotij(:,:,site)),pos(:,first))+tauij(:,site)
          do l=1,labc
             write(unit) couplo(site,l,-l:l,-l:l,:,:)
          end do
       end do
       first=first+mult(atom)
    end do
    close(unit)
  end subroutine

  subroutine init_soc_displacement(nsites, nrel_lo, jspin, have_vectors)
    integer, intent(in) :: nsites, nrel_lo, jspin
    logical, intent(in) :: have_vectors
    integer :: status, unit, site
    character(64) :: layout
    call get_environment_variable('EDMFT_SOC_DISPLACEMENT',prefix,status=status)
    if (status == 1 .or. len_trim(prefix)==0) return
    if (status/=0) error stop 'EDMFT_SOC_DISPLACEMENT path is too long'
    if (nrel_lo/=0 .or. jspin/=1 .or. .not.have_vectors) &
         error stop 'SOC displacement requires scalar nonmagnetic orbitals, vectors and no p1/2 RLO'
    open(newunit=unit,file=trim(prefix),status='old',action='read')
    read(unit,*) layout
    select case(layout)
    case('subspace_matrix','subspace_diagonal')
       scalar_window=.true.
    case('matrix','diagonal')
       scalar_window=.false.
    case default
       error stop 'unknown SOC derivative layout'
    end select
    full_matrix=layout=='matrix' .or. layout=='subspace_matrix'
    allocate(velocity(3,nsites))
    velocity=0.d0
    if (.not.scalar_window) then
       do site=1,nsites
          read(unit,*) velocity(:,site)
       end do
    endif
    read(unit,*,iostat=status) layout
    if(status>=0) error stop 'unexpected additional SOC displacement rows'
    close(unit)
    if (.not.all(ieee_is_finite(velocity))) error stop 'nonfinite SOC displacement'
  end subroutine

  ! Read Q_out*dA from the generalized scalar eigenproblem. The producer
  ! fixes the retained-subspace gauge by projecting the tangent exclusively
  ! into excluded scalar states. Eigenvalues alone cannot identify the
  ! original eigenvectors inside degenerate subspaces, so check both.
  subroutine read_scalar_tangent(k,frac,energies,a,n,b,tangent)
    integer, intent(in) :: k,n,b
    real(real64), intent(in) :: frac(3),energies(:)
    complex(real64), intent(in) :: a(:,:)
    complex(real64), allocatable, intent(out) :: tangent(:,:)
    real(real64), allocatable :: e(:)
    real(real64) :: kfile(3)
    integer(int64), allocatable :: gauge(:)
    integer(int32) :: header(3)
    character(8) :: magic
    character(16) :: id
    integer :: unit,j,ios
    write(id,'(i0)') k
    open(newunit=unit,file=trim(prefix)//'.dA.k'//trim(id),status='old',access='stream',form='unformatted',action='read')
    read(unit) magic,header,kfile
    if(magic/='SOCWIN01' .or. any(header/=[k,n,b]) .or. maxval(abs(kfile-frac))>1.d-11) &
         error stop 'SOC scalar-window tangent: incompatible k or dimensions'
    allocate(e(b),gauge(b),tangent(size(a,1),size(a,2)))
    read(unit) e,gauge,tangent(1:n,1:b)
    read(unit,iostat=ios) magic(1:1)
    close(unit)
    if(ios>=0 .or. .not.all(ieee_is_finite(e)) .or. maxval(abs(e-energies(1:b)))>1.d-11) &
         error stop 'SOC scalar-window tangent: incompatible record or eigenvalues'
    do j=1,b
       if(gauge(j)/=column_gauge(a(1:n,j))) error stop 'SOC scalar-window tangent: scalar gauge differs'
       if(.not.all(ieee_is_finite(real(tangent(1:n,j),real64))) .or. &
            .not.all(ieee_is_finite(aimag(tangent(1:n,j))))) error stop 'nonfinite SOC scalar tangent'
    end do
  end subroutine

  pure function column_gauge(a) result(hash)
    complex(real64), intent(in) :: a(:)
    integer(int64) :: hash
    integer :: j
    hash=0_int64
    do j=1,size(a)
       hash=ieor(ishftc(hash,1),transfer(real(a(j),real64),0_int64))
       hash=ieor(ishftc(hash,7),transfer(aimag(a(j)),0_int64))
    end do
  end function

  subroutine accumulate_soc_translation(hdot, base, tangent, radial, nband, site, atom)
    use param, only: labc, nloat, nato
    use lolog, only: mrf
    use couples, only: couplo
    integer, intent(in) :: nband,site,atom
    complex(real64), intent(inout) :: hdot(2*nband,2*nband)
    complex(real64), intent(in) :: base(:,:,:), tangent(:,:,:)
    real(real64), intent(in) :: radial(nloat,nloat,0:labc,nato,2,2)
    integer :: l,m,nr,nrow,m1,m2,i,j,s1,s2,offset1,offset2
    complex(real64), allocatable :: b(:,:),db(:,:),r(:,:),rb(:,:),rdb(:,:)
    complex(real64), parameter :: one=(1.d0,0.d0)
    ! Angular spin matrices are already rotated to the native quantization
    ! axis. Each radial block couples only equal l. hu,hd labels have the
    ! same ordering as HMSEC, including for a complex scalar vector.
    do l=1,labc
       nr=mrf(l,atom)
       nrow=nr*(2*l+1)
       allocate(b(nrow,nband),db(nrow,nband),r(nrow,nrow),rb(nrow,nband),rdb(nrow,nband))
       do m=-l,l
          i=(m+l)*nr
          b(i+1:i+nr,:)=base(1:nr,l*l+l+m+1,1:nband)
          db(i+1:i+nr,:)=tangent(1:nr,l*l+l+m+1,1:nband)
       end do
       do s2=1,2
          offset2=(s2-1)*nband
          do s1=1,2
             offset1=(s1-1)*nband
             do m2=-l,l
                j=(m2+l)*nr
                do m1=-l,l
                   i=(m1+l)*nr
                   r(i+1:i+nr,j+1:j+nr)=couplo(site,l,m1,m2,s1,s2)*transpose(radial(1:nr,1:nr,l,atom,1,1))
                end do
             end do
             rb=matmul(r,b)
             rdb=matmul(r,db)
             ! Accumulate directly into each spin block with the leading
             ! dimension of the full matrix; avoid a dense temporary per l.
             call zgemm('C','N',nband,nband,nrow,one,db,nrow,rb,nrow,one, &
                  hdot(offset1+1,offset2+1),size(hdot,1))
             call zgemm('C','N',nband,nband,nrow,one,b,nrow,rdb,nrow,one, &
                  hdot(offset1+1,offset2+1),size(hdot,1))
          end do
       end do
       deallocate(b,db,r,rb,rdb)
    end do
  end subroutine

  subroutine write_soc_translation(k,frac,w,eval,vec,hdot,nbands_so,nbasis,spinors,nlapw)
    integer, intent(in) :: k,nbands_so,nbasis,nlapw
    real(real64), intent(in) :: frac(3),w,eval(:)
    complex(real64), intent(in) :: vec(:,:),hdot(:,:),spinors(:,:,:)
    complex(real64), allocatable :: applied(:,:),projected(:,:)
    real(real64), allocatable :: diagonal(:)
    integer(int64), allocatable :: gauge(:,:)
    integer :: j,spin,n,unit
    character(12) :: id
    write(id,'(i0)') k
    applied=matmul(hdot,vec(1:nbasis,1:nbands_so))
    allocate(diagonal(nbands_so))
    do n=1,nbands_so
       diagonal(n)=real(dot_product(vec(1:nbasis,n),applied(:,n)),real64)
    end do
    open(newunit=unit,file=trim(prefix)//'.k'//trim(id),form='unformatted',access='stream',status='replace')
    write(unit) int([2,k,nbasis,nbands_so,merge(1,0,full_matrix)],int32),frac,w,eval(1:nbands_so),diagonal
    if(full_matrix) then
       projected=matmul(conjg(transpose(vec(1:nbasis,1:nbands_so))),applied)
       ! A degenerate eigenspace can rotate without changing a single
       ! eigenvalue. Pair the gradient with the exact two-component vectors
       ! read by DMFT2. Print the hash in DMFT2 order (up, down); vectors
       ! are assembled here in LAPWSO order (down, up). No p1/2 RLOs in
       ! this export, hence K and LO rows retain their final file order.
       allocate(gauge(nbands_so,2),source=0_int64)
       do spin=1,2
          do n=1,nbands_so
             gauge(n,spin)=column_gauge(spinors(1:nlapw,n,3-spin))
          end do
       end do
       write(unit) projected,gauge
    end if
    close(unit)
  end subroutine
end module
