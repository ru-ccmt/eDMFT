INTEGER FUNCTION NOTRI(K,L,M)
  !IMPLICIT REAL*8 (A-H,O-Z)
  IMPLICIT NONE
  INTEGER, intent(in) :: K,L,M
  NOTRI=-1
  IF(MOD((K+L+M),2).EQ.1) RETURN
  IF((K+L-M).LT.0) RETURN
  IF((K-L+M).LT.0) RETURN
  IF((M+L-K).LT.0) RETURN
  NOTRI=1
  RETURN 
END FUNCTION NOTRI
