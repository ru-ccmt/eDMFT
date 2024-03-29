SUBROUTINE HARMON(N,aK,LMAX2,F,DF,RI)                           
  use structure, ONLY: BR1
  IMPLICIT REAL*8 (A-H,O-Z)
  REAL*8, intent(out):: F(LMAX2+1,N),DF(LMAX2+1,N)
  INTEGER, intent(in):: N, LMAX2
  REAL*8, intent(in) :: aK(N), RI
  ! locals
  LMX=LMAX2+1                                                        
  DO I=1,N
     XA=RI*aK(I)
     CALL SPHBES(LMAX2,XA,F(1,I))
     CALL DVBES1(F(1,I),DF(1,I),XA,RI,LMX)
     DO J=1,LMX
        DF(J,I)=aK(I)*DF(J,I)
     ENDDO
  ENDDO
  RETURN
END SUBROUTINE HARMON

