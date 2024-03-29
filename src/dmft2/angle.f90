SUBROUTINE ANGLE(XMS,THETA,PHI)
  USE structure, ONLY: LATTIC, alpha, ORTHO, aa, bb, cc
  USE param
  USE com_mpi,ONLY: myrank, master
  IMPLICIT NONE
!*******************************************************************
  CHARACTER*10     ANA
  REAL*8, intent(in) :: XMS(3)
  REAL*8, intent(out):: THETA, PHI
  ! locals                                                                                                                                                                                                  
  REAL*8 :: PI, XA, XB, XC, TT, cosg1, gamma0, XX
  
  !LOGICAL          ORTHO
  !
  !COMMON/ORTH/     ORTHO
  !DIMENSION XMS(3)

  PI=ACOS(-1.D0)
  !---------------------------------------------------------------------
  IF (ORTHO) THEN
     XA=AA*XMS(1)
     XB=BB*XMS(2)
     XC=CC*XMS(3)
     GOTO 200
  END IF

  write(6,*)'LATTICE:',lattic	
  IF(LATTIC(1:1).EQ.'H') THEN
     XA=XMS(1)*AA*SQRT(3.D0)/2.d0
     XB=AA*(XMS(2)-XMS(1)/2.d0)
     XC=CC*XMS(3)

  ELSE IF(LATTIC(1:1).EQ.'R') THEN
     XA=(XMS(1)+XMS(2)-2.d0*XMS(3))*AA/(2.d0*SQRT(3.D0))
     XB=(-XMS(1)+XMS(2))*AA/2.d0
     XC=(XMS(1)+XMS(2)+XMS(3))*CC/3.d0

  ELSE 
     ! alpha1=alpha2=Pi/2 , alpha3!=Pi/2 
     IF ((ABS(ALPHA(3)-PI/2.d0).GT.1.D-4).and.(ABS(ALPHA(2)-PI/2.d0).LT.1.D-4).and.(ABS(ALPHA(1)-PI/2.d0).LT.1.D-4)) THEN
        XA=XMS(1)*AA*SIN(ALPHA(3))
        XB=XMS(1)*AA*COS(ALPHA(3))+BB*XMS(2)
        XC=CC*XMS(3)
     ! alpha1=alpha3=Pi/2,  alpha2!=Pi/2
     ELSE IF ((ABS(ALPHA(2)-PI/2.d0).GT.1.D-4).and.(ABS(ALPHA(1)-PI/2.d0).LT.1.D-4).and.(ABS(ALPHA(3)-PI/2.d0).LT.1.D-4)) THEN
	XA=XMS(1)*AA*SIN(ALPHA(2))
        XB=XMS(2)*BB
	XC=XMS(1)*AA*COS(ALPHA(2))+CC*XMS(3)
     ELSE
        !	WRITE(6,*)'EXCHANGE THE LATTICE VECTORS, ALPHA(1) MUST BE PI/2'
        !        END IF

	TT=ABS((ALPHA(1)-PI/2.d0)*(ALPHA(2)-PI/2.d0)*(ALPHA(3)-PI/2.d0))
        !	IF (TT.GT.1.D-3) STOP 'TRICLINIC NOT IMPLEMENTED'
        write(6,*) ' Triclinic implemented, but never tested'
        cosg1=(cos(ALPHA(3))-cos(alpha(1))*cos(ALPHA(2)))/sin(alpha(1))/sin(alpha(2))
        gamma0=acos(cosg1)
        !      from lapw5
        !  vec{a} = a(sin(gamma0)*sin(beta),cos(gamma0)*sin(beta),cos(beta))
        !  vec{b} = b(0,sin(alpha),cos(alpha))
        !  vec{c} = c(0,0,1)
        !  (xa,xb,xc)= xms(1)*vec{a}+xms(2)*vec{b}+xms(3)*vec{c}
        !      BR2(1,1)=A(1)*sin(gamma0)*sin(beta)
        !      BR2(1,2)=A(1)*cos(gamma0)*sin(beta)
        !      BR2(1,3)=A(1)*cos(beta)    
        !      BR2(2,1)=0.0d0              
        !      BR2(2,2)=A(2)*sin(alpha)
        !      BR2(2,3)=A(2)*cos(alpha)
        !      BR2(3,1)=0.0d0              
        !      BR2(3,2)=0.0d0              
        !      BR2(3,3)=A(3)*1.0d0              
        XA=XMS(1)*sin(gamma0)*sin(alpha(2))*AA
        XB=XMS(1)*cos(gamma0)*sin(alpha(2))*AA + XMS(2)*sin(alpha(1))*BB
        XC=XMS(1)*cos(alpha(2))*AA + XMS(2)*cos(alpha(1))*BB + XMS(3)*CC
     END IF

  ENDIF

200 CONTINUE  
  XX=SQRT(XA**2+XB**2+XC**2)
  THETA=ACOS(XC/XX)
  XX=SQRT(XA**2+XB**2)
  IF (XX.LT.1.D-5) THEN
     PHI=0.D0
  ELSE
     PHI=ACOS(XA/XX)
     IF (ABS(XB).GT.1.D-5) PHI=PHI*XB/ABS(XB)
  END IF

END SUBROUTINE ANGLE
	



