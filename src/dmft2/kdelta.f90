LOGICAL FUNCTION KDELTA(K,G,NST,STG)
  !                                                                       
  !.... TEST, IF K IS IN STAR OF G  (GENERATED IN STERN)
  ! Routine used in recpr
  !                                                                       
  !
  USE param
  IMPLICIT REAL*8 (A-H,O-Z)
  INTEGER        G,NST,STG(3,NSYM)
  DIMENSION      K(3)                                               
!---------------------------------------------------------------------  
  DO 1 I=1,NST                                                      
     DO 2 J=1,3 
        IF(STG(J,I).NE.K(J)) GOTO 1
2    CONTINUE                                                          
     KDELTA=.TRUE.                                                     
     RETURN                                                            
1 CONTINUE                                                          
  KDELTA=.FALSE.                                                    
  RETURN                                                            
END FUNCTION KDELTA

