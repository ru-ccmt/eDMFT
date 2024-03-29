SUBROUTINE ROTATE(VECTOR,ROTMAT,ROTVEC)
  !                                                                       
  !     ROTATE PERFORMS A ROTATION OF THE VECTOR FROM THE GENERAL         
  !     CARTESIAN COORDINATION SYSTEM INTO THE  LOCAL ONE  OF THE         
  !     JATOM-TH SPHERE.                                                  
  !     THIS SUBROUTINE IS ONLY REQUIRED FOR NONSYMMORPHIC CASES.         
  !                                                                       
  IMPLICIT NONE
  REAL*8, intent(in) :: VECTOR(3), ROTMAT(3,3)
  REAL*8, intent(out):: ROTVEC(3)
  !---------------------------------------------------------------------- 
  !
  ROTVEC = matmul(ROTMAT,VECTOR)
  !DO JCOORD=1,3     ! 10                                                
  !   DOTPRO=0.0D0                                                   
  !   DO J=1,3       ! 20
  !      DOTPRO = DOTPRO + ROTMAT(JCOORD,J)*VECTOR(J)
  !   ENDDO
  !   ROTVEC(JCOORD)=DOTPRO
  !ENDDO !10   CONTINUE                                                          
  RETURN                                                            
END SUBROUTINE ROTATE
