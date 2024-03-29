      SUBROUTINE symop (NAT,OUTME)
      use struk
!*
!*     transforms symmetri matrices to cartesian system
!*
!*     for .orhto. systems   opimat == imat
!* 
!*     for not.ortho. systems opimat == BR1 . imat . BR1^-1 
!*
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'param.inc'
      REAL*8   BR1in(3,3),te(3,3)
      CHARACTER*3 OUTME
      LOGICAL         ORTHO
!ad 
      COMMON /ORTH/   ORTHO
!      COMMON /STRUK/  POS(3,NDIF),AA,BB,CC,ALPHA(3),RMT(NATO),V(NATO), &
!                      PIA(3),VOL,ZZ(NATO), &
!                      IATNR(NATO),MULT(NATO),ISPLIT(NATO)
      COMMON /SYM2/   TAU(3,NSYM),IORD,IMAT(3,3,NSYM)
      COMMON /SYMo/   opimat(3,3,NSYM)
      COMMON /GENER/  BR1(3,3),BR2(3,3)
!ad
!.........inverssymdef is a subroutine in sph-UP.frc and
!.........calculates the inverse of an 3*3 matrix.......
!ad
      CALL INVERSSYMDEF(BR1,BR1in)
!ad
!ad   WRITE(6,*) 'BRAVAIS matrix ......'
!ad   WRITE(6,12) ( (br1(J2,J1),J1=1,3), J2=1,3 )
!ad   WRITE(6,*) 'inverse BRAVAIS matrix ......'
!ad   WRITE(6,12) ( (br1in(J2,J1),J1=1,3), J2=1,3 )
         
          do j=1,3
          do k=1,3
          te(j,k)= 0.d0
          do l=1,3
          te(j,k)= te(j,k) + &
               BR1(j,l)  * BR1in(l,k)
            end do
            end do
            end do
! 
!.......define symmetry matrices for kartesian system
!
        do i=1,IORD
          do j=1,3
          do k=1,3
          opimat(j,k,i)=0.d0 
          end do
          end do
        end do
 
        do i=1,IORD
          do j=1,3
          do k=1,3
            do l=1,3
            do m=1,3
          opimat(j,k,i)= opimat(j,k,i) +  &
               BR1(j,l) * IMAT(l,m,i) * BR1in(m,k)
            end do
            end do
          end do
          end do
        end do
! 
!.......define symmetry matrices for kartesian system
!
!ole ##### Begin #####
!ole store optmat(line,col)
!ole
      WRITE (25,100) IORD
      DO k = 1,IORD
         WRITE(25,110) ((opimat(j,i,k),i=1,3),j=1,3)
      END DO

!      WRITE (25,100) IORD
!      DO k = 1,IORD
!         WRITE(25,110) ((DBLE(IMAT(j,i,k)),i=1,3),j=1,3)
!      END DO


 100  FORMAT (I6)
 110  FORMAT (3(3f8.5/))
!ole #####  End  #####
!sas
!sas ________________ WRITE OUT SYMMETRY OPERATIONS __________________
!sas                        needed for NLO
!sas
       if(OUTME.eq.'ON ') then
       write(8,*) iord
       do l=1,iord
       write(8,'(3(3f8.5))') ((opimat(i,ii,l),ii=1,3),i=1,3)
       end do
       endif
!sas
!sas
!sas __________________________________________________________________
       return
 11   FORMAT(3(3I2,F10.5/))
 12   FORMAT(3(3f8.5/))

       end
