subroutine get_nloat(lomax,nat,nloat)
  IMPLICIT NONE
  integer, intent(in) :: lomax, nat
  integer, intent(out):: nloat
  ! locals
  integer ::   i,ios, ii
  real*8, allocatable :: elo(:,:)
  nloat=1000
  allocate (elo(0:lomax,1:nloat))
  ii=1
  nloat=0
  DO i=1,nat
     read(9)  ! linearization energy for lapw
     do
        read(9,iostat=ios) elo(0:lomax,1:ii)
        if (ios.ne.0) exit ! no more entries
        ii=ii+1
        backspace(9)
     enddo
     nloat=ii
  ENDDO
  deallocate(elo)
  rewind(9)
end subroutine get_nloat
