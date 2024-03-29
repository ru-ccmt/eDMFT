SUBROUTINE RINT13(REL,A,B,X,Y,S,JATOM)                            
  use potnlc
  use xrpar
  !
  !LO S = INT { AX + CIN BY } = INT { AX + (1/ALPHA) BY }   
  !LO                                                                    
  INCLUDE 'param.inc'
  IMPLICIT REAL*8 (A-H,O-Z)
  LOGICAL REL                                                       
  DIMENSION A(NRAD),B(NRAD),X(NRAD),Y(NRAD)                         
  !      COMMON /POTNLC/ Z(NRAD),RNOT(NATO),DX(NATO),JRI(NATO)             
  D=EXP(DX(JATOM))                                                  
  if (xmcd.eq.1) then
     CIN=7.29927D-3    ! for core, CIN=ALPHA=1/137    
  else
     CIN=1.331258D-5   ! CIN=1/137^2
     IF(.NOT.REL) CIN=1E-22
     CIN=CIN*4
  endif
  !       CIN=1.0D0                                            
  !      IF(.NOT.REL) CIN=1E-22                                            
  !                                                                       
  !     CONVERT FROM RYDBERG TO HARTREE UNITS                             
  !                                                                       
  !      CIN=CIN*4                                                         
  J=3-MOD(JRI(JATOM),2)                                             
  J1=J-1                                                            
  R=RNOT(JATOM)*(D**(J-1))                                          
  R1=R/D                                                            
  Z4=0                                                              
  Z2=0                                                              
10 Z4=Z4+R*(A(J)*X(J)+CIN*B(J)*Y(J))                                 
  R=R*D                                                             
  J=J+1                                                             
  IF(J.GE.JRI(JATOM)) GOTO 20                                       
  Z2=Z2+R*(A(J)*X(J)+CIN*B(J)*Y(J))                                 
  R=R*D                                                             
  J=J+1                                                             
  GOTO 10                                                           
20 P1=RNOT(JATOM)*(A(1)*X(1)+CIN*B(1)*Y(1))                          
  P2=R1*(A(J1)*X(J1)+CIN*B(J1)*Y(J1))                               
  S=2*Z2+4*Z4+R*(A(J)*X(J)+CIN*B(J)*Y(J))+P2                        
  S=(DX(JATOM)*S+P1)/3.0D0                                          
  IF(J1.GT.1) S=S+0.5D0*DX(JATOM)*(P1+P2)                           
  RETURN                                                            
END SUBROUTINE RINT13
