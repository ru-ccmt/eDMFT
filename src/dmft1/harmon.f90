SUBROUTINE HARMON(N,aK,LMAX2,F,DF,RI)
  ! Input: N,X,Y,Z,LMAX2,RI
  ! Output: F,DF
  !
  IMPLICIT NONE
  REAL*8, intent(in)  :: aK(N)
  REAL*8, intent(in)  :: RI
  INTEGER, intent(in) :: N, LMAX2
  REAL*8, intent(out) :: F(LMAX2+1,N),DF(LMAX2+1,N)
  !locals
  INTEGER :: i
  REAL*8  :: XA

  DO I=1,N 
     XA=RI*aK(I)
     CALL SPHBES(LMAX2,XA,F(1,I))
     CALL DVBES1(F(1,I),DF(1,I),XA,RI,LMAX2+1)
     DF(:,I) = aK(I)*DF(:,I)
  ENDDO
  
  RETURN
END SUBROUTINE HARMON
