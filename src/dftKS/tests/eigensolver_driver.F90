! Standalone harness: compile before diagonalize.FPP, without LAPW modules.o.
module comi
  integer :: nemax_it
  logical :: fixed_nbands
end module
program test_solver
  use comi
  implicit none
  integer :: n,ldz,fixed,nlo,ne,unit,labels_dummy,i
  integer, allocatable :: labels(:,:)
  real(8) :: lo,hi,kfrac(3)
  real(8), allocatable :: eig(:)
#ifdef _REAL
  real(8), allocatable :: raw(:,:),hd(:),z(:,:)
#else
  complex(8), allocatable :: raw(:,:),hd(:),z(:,:)
#endif
  character(2048) :: input,output,cache
  call get_command_argument(1,input)
  call get_command_argument(2,output)
  call get_command_argument(3,cache)
  open(newunit=unit,file=trim(input),access='stream',form='unformatted',status='old')
  read(unit) n,ldz,fixed,nemax_it,nlo,lo,hi,kfrac
  fixed_nbands=fixed/=0
  allocate(labels(3,n),raw(n,n),hd(n),z(n,ldz),eig(n))
  read(unit) labels,raw,hd
  close(unit)
  call Diagonalize(lo,hi,ne,raw,hd,z,eig,n,ldz,.true.,1,labels,kfrac,nlo,trim(cache))
  open(newunit=unit,file=trim(output),access='stream',form='unformatted',status='replace')
  write(unit) n,ne,eig(1:ne),z(:,1:ne)
  close(unit)
end program
subroutine OUTERR(where,message)
  character(*) :: where,message
  write(*,'(A,1X,A)') trim(where),trim(message)
end subroutine
