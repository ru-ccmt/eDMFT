SUBROUTINE LAGDER(N,X,Y,DER)
  IMPLICIT NONE
  REAL*8, intent(in) :: X(N),Y(N)
  INTEGER, intent(in):: N
  REAL*8, intent(out):: DER
  ! locals
  INTEGER:: i, j, k
  REAL*8 :: ade, anu, de, xu
  XU=X(1)
  DE=0.D0
  DO K=2,N
     DE=DE+1.D0/(XU-X(K))
  ENDDO
  DE=DE*Y(1)
  
  DO I=2,N
     ANU=1.D0
     DO J=2,N
        IF(J.EQ.I) CYCLE
        ANU=ANU*(XU-X(J))
     ENDDO
     ADE=1.D0
     DO J=1,N
        IF(J.EQ.I) CYCLE
        ADE=ADE*(X(I)-X(J))
     ENDDO
     DE=DE+Y(I)*ANU/ADE
  ENDDO
  DER=DE
  RETURN
END SUBROUTINE LAGDER
