REAL*8 Function Romb(y,N,dx)
  IMPLICIT NONE
  REAL*8, intent(in) :: y(N)
  REAL*8, intent(in) :: dx
  INTEGER, intent(in):: N
  ! locals
  INTEGER :: Ninterv, ni, k, istart, istop, istep, i, j, nn
  REAL*8, allocatable:: R(:,:)
  REAL*8 :: h, dsum
  !maybe we should use JISHFT instead of ISHFT!
  !
  Ninterv = N-1
  h = Ninterv*dx
  ni = 1
  k = 0
  do while(ni.LT.Ninterv)
     ni = ISHFT(ni,1)
     k = k+1
  enddo
  
  allocate( R(k+1,k+1) )
  
  R(1,1) = (y(1) + y(N))/2.0*h
    
  istart = Ninterv
  istop  = Ninterv
  istep = Ninterv
  do i=2,k
     istart = ISHFT(istart,-1)

     dsum=0.0
     do j=istart,istop-1,istep
        dsum = dsum + y(j+1)
     enddo
     
     istep  = ISHFT(istep,-1)
        
     R(i,1) = 0.5*(R(i-1,1) + h*dsum)

     do j=2,i
        nn = ISHFT(1, (2*(j-1)))
        R(i,j) = R(i,j-1) + (R(i,j-1)-R(i-1,j-1))/(nn-1.)
     enddo
     h = h/2.0
  enddo
  Romb = R(k,k)
  deallocate( R)
  return
END Function Romb
