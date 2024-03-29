      SUBROUTINE ROTATE (VECTOR,ROTMAT,ROTVEC)                          
!                                                                       
!     ROTATE PERFORMS A ROTATION OF THE VECTOR FROM THE GENERAL         
!     CARTESIAN COORDINATION SYSTEM INTO THE  LOCAL ONE  OF THE         
!     JATOM-TH SPHERE.                                                  
!     THIS SUBROUTINE IS ONLY REQUIRED FOR NONSYMMORPHIC CASES.         
!                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION VECTOR(3),ROTVEC(3),ROTMAT(3,3)                         
!---------------------------------------------------------------------- 
!                                                                       
      DO 10 JCOORD=1,3                                                  
         DOTPRO=0.0D0                                                   
         DO 20 J=1,3                                                    
            DOTPRO=DOTPRO + VECTOR(J)*ROTMAT(JCOORD,J)                  
 20      CONTINUE                                                       
         ROTVEC(JCOORD)=DOTPRO                                          
 10   CONTINUE                                                          
      RETURN                                                            
      END                                                               
