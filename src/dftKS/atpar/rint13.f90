FUNCTION RINT13(REL,A,B,X,Y,NRAD,DX_,JRI_,R0_) result(S)
  ! Input: A, B, X, Y, jatom
  !
  ! Output: S
  !                                                                       
  !     PERFORM RADIAL INTEGRALS REQUIRED BY BHDK13                       
  !                            D.D.KOELLING                               
  ! We are integrating f(r) = A(r)*X(r)+ B(r)*Y(r)*1/137 
  !  where r-mesh is exponential in the form r=e^{i*dx}
  !  therefore the integral is Int f(r)dr = Int f(e^{i*dx})e^{i*dx}*dx
  !  which can also be written Sum f(r_i)*r_i*dx
  ! The simpson method then collects the even and odd terms and
  ! adds 2/3 and 4/3 to each, and also takes care of the end points.
  !
  IMPLICIT NONE
  !
  REAL*8 :: S
  LOGICAL, intent(in) :: REL   ! relativistic or not
  REAL*8, intent(in)  :: A(NRAD),B(NRAD),X(NRAD),Y(NRAD)
  INTEGER, intent(in) :: NRAD, JRI_
  REAL*8, intent(in)  :: DX_, R0_
  !locals
  REAL*8 :: D, CIN, R, R1, P1, P2, Z2, Z4
  INTEGER :: J, J1
  !
  D=EXP(DX_)                                                  
  CIN=1.d0/137.0359895d0**2                           
  IF(.NOT.REL) CIN=1E-22                                            
  !                                                                       
  J=3-MOD(JRI_,2)                                             
  J1=J-1                                                            
  R=r0_*(D**(J-1))                                          
  R1=R/D                                                            
  Z4=0                                                              
  Z2=0                                                              
  DO
     Z4=Z4+R*(A(J)*X(J)+CIN*B(J)*Y(J)) ! Z4 ~ R*A*X (r)
     R=R*D                                                             
     J=J+1                                                             
     IF(J.GE.JRI_) EXIT
     Z2=Z2+R*(A(J)*X(J)+CIN*B(J)*Y(J)) ! Z2 ~ R*A*X (r+dr)                               
     R=R*D                                                             
     J=J+1                                                             
  ENDDO
  P1=r0_*(A(1)*X(1)+CIN*B(1)*Y(1))
  P2=R1*(A(J1)*X(J1)+CIN*B(J1)*Y(J1))
  S=2*Z2+4*Z4+R*(A(J)*X(J)+CIN*B(J)*Y(J))+P2
  S=(DX_*S+P1)/3.0D0
  IF(J1.GT.1) S=S+0.5D0*DX_*(P1+P2)
  RETURN                                                            
END FUNCTION RINT13

FUNCTION RINT13g(C1,C2,A,B,X,Y,NRAD,DX_,JRI_,R0_) result(S)
  ! C1=1, C2=[1.d0/137.0359895d0**2 or 0]
  IMPLICIT NONE
  !     PERFORM RADIAL INTEGRALS REQUIRED BY BHDK13
  !     Nornally C1=1 and C2=1/speed_of_light^2
  !        but here C1 and C2 are arbitrary
  !                            D.D.KOELLING                               
  REAL*8 :: S
  REAL*8, intent(in) :: C1, C2
  REAL*8, intent(in) :: A(NRAD), B(NRAD), X(NRAD), Y(NRAD)
  INTEGER, intent(in):: NRAD, JRI_
  REAL*8, intent(in) :: DX_, R0_
  ! locals
  INTEGER::            J, J1
  REAL*8 ::   D, P1, P2, R, R1, Z2, Z4
  D = EXP(DX_)
  J = 3-MOD(JRI_,2)
  J1=J-1
  R = r0_*(D**(J-1))
  R1 = R/D
  Z4 = 0
  Z2 = 0
  DO
     Z4=Z4+R*(C1*A(J)*X(J) + C2*B(J)*Y(J))     ! Z4 ~ R*A*X (r)
     R=R*D
     J=J+1
     IF (J .GE. JRI_) EXIT !GOTO 20
     Z2 = Z2 + R*(C1*A(J)*X(J) + C2*B(J)*Y(J))  ! Z2 ~ R*A*X (r+dr)   
     R=R*D
     J=J+1
  enddo
  P1 = R0_*(C1*A(1)*X(1) + C2*B(1)*Y(1))
  P2 = R1*(C1*A(J1)*X(J1) + C2*B(J1)*Y(J1))
  S = 2*Z2 + 4*Z4 + R*(C1*A(J)*X(J) + C2*B(J)*Y(J)) + P2
  S = (DX_*S+P1)/3.0D+0
  IF (J1.GT.1) S=S+0.5D0*DX_*(P1+P2)
  RETURN
END FUNCTION RINT13g

REAL*8 FUNCTION Rint13n(N,A,B,X,Y,r0,dx,rel)
  ! PERFORM RADIAL INTEGRALS SIMILAR TO RINT13 in modern form
  !  D.D.KOELLING                               
  IMPLICIT NONE
  INTEGER, intent(in) :: N
  REAL*8, intent(in)  :: A(N), B(N), X(N), Y(N)
  REAL*8, intent(in)  :: r0, dx
  LOGICAL, intent(in) :: rel
  ! locals
  REAL*8  :: S
  REAL*8  :: CIN, D, Z4, Z2, P1, P2, PN
  INTEGER :: j, i
  REAL*8 :: FM(N), RX(N)
  
  D=exp(dx)
  RX(1)=r0
  DO j=2,N
     RX(j)=RX(j-1)*D
  ENDDO
  
  IF (rel) then  ! relativistic or no
     CIN=1.d0/137.0359895d0**2
     FM(:) = RX(:)*(A(:)*X(:)+CIN*B(:)*Y(:))
  else
     FM(:) = RX(:)*A(:)*X(:)
  endif
  
  j=3-MOD(N,2)
  
  P1=FM(1)
  P2=FM(j-1)
  
  Z4=0
  Z2=0
  DO i=j,N-3,2
     Z4=Z4+FM(i)
     Z2=Z2+FM(i+1)
  ENDDO
  Z4=Z4+FM(N-1)
  PN=FM(N)
  
  S=(dx*(2*Z2+4*Z4+PN+P2)+P1)/3.0D0
  if (j-1.GT.1) S=S+0.5D0*dx*(P1+P2)
  Rint13n=S
  Return
END FUNCTION Rint13n

