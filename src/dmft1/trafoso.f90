!-*- f90 -*- - file contains F90 code in free format

subroutine trafoso(L,cf)
  implicit none
  integer, intent(in)  :: L
  real*8,  intent(out) :: cf(2*(2*L+1),2*(2*L+1))
  ! locals
  integer :: m1, m2, k1, k2, ms, ml, mj
  real*8  :: ams, amj, d
! trafoso provides transformation matrices from
! |L,1/2,mL,mS> (L=0,1,2,3, mS=-1/2,1/2) basis to
! basis |J,L,S,mJ>, J=L-1/2, L+1/2 
! H. Watanabe 'Operator Methods in Ligand Field Theory'
! Prentice Hall, 1966, Table 1.8-1.
! ordering because of the convention used in WIEN is:
!                    mS=1/2        mS=-1/2
!                  -L .......L  -L ...... L     (2*(2L+1) columns)
!         -(L-1/2)
!            .
! J=L-1/2    .
!            .
!          (L-1/2)
!          -L-1/2 
!            .
! J=L+1/2    .
!            .
!           L+1/2 
!
  cf=0
  if(L.eq.0)then
     cf(1,2)=1.
     cf(2,1)=1.
  else 
     k1=0
     do ms=-1,1,2
        ams=-ms/2.
        do ml=-L,L
           k1=k1+1
           k2=0
           do mj=-2*L+1,2*L-1,2   ! L-1/2 states
              amj=mj/2.
              k2=k2+1
              d=amj-mL-ams
              if(abs(d).lt.0.0001)then
                 if(ms.eq.1)then
                    cf(k2,k1)=-sqrt((L+0.5+amj)/(2*L+1))
                 else
                    cf(k2,k1)= sqrt((L+0.5-amj)/(2*l+1))
                 endif
              endif                ! L-1/2 states end
           enddo
           do mj=-2*L-1,2*L+1,2   ! L+1/2 states
              amj=mj/2.
              k2=k2+1
              d=amj-mL-ams
              if(abs(d).lt.0.0001)then
                 if(ms.eq.1)then
                    cf(k2,k1)= sqrt((L+0.5-amj)/(2*L+1))
                 else
                    cf(k2,k1)= sqrt((L+0.5+amj)/(2*l+1))
                 endif
              endif                ! L+1/2 states end
           enddo
        enddo
     enddo
  endif
  return
end subroutine trafoso
