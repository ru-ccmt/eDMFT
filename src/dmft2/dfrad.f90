subroutine dfrad(df,f,r,nrad)
  ! Computes the radial derivative of function f(r)
  IMPLICIT NONE
  REAL*8, intent(out):: df(nrad)
  INTEGER, intent(in):: nrad
  REAL*8, intent(in) :: r(nrad), f(nrad)
  ! locals
  REAL*8 ::  C(nrad), D(nrad)
  call spline(nrad, r, f, df, C, D)
  return
end subroutine dfrad

REAL*8 FUNCTION int_radial(funct,RX,DELTA,nr)
  ! Integrates radial function by Simpson: int funct(rx)*drx
  ! Notice that rx is logarithmic mesh, hence int[func(r)] requires summation of func(r)*r
  IMPLICIT NONE
  INTEGER, intent(in):: nr
  REAL*8, intent(in) :: RX(nr), funct(nr), DELTA
  ! locals
  REAL*8 :: dsum
  INTEGER :: jstart, jfin, j
  dsum=0.0D0
  jstart = 2-MOD(nr,2)
  jfin = nr-2
  ! Integration by Simpson  Integrate[f(r) r dr]
  ! dsum = (1/3.*f_0*r0 + 4/3.*f1*r1 + 2/3.*f2*r2 + 4/3.*f3*r3+....+fn*rn)*dx
  DO j=jstart,jfin,2
     dsum = dsum + funct(j)*RX(j)+4.*funct(j+1)*RX(j+1)+funct(j+2)*RX(j+2)
  ENDDO
  dsum = dsum * DELTA/3.0D0
  dsum = dsum + (1-MOD(nr,2))*(funct(1)+funct(jstart))/2.0D0*(RX(jstart)-RX(1)) ! correction for even number of mesh points
  int_radial = dsum
  RETURN                                                            
end FUNCTION int_radial


SUBROUTINE SPLINE(N,X,Y,B,C,D)
  IMPLICIT NONE
  INTEGER, intent(in) :: N
  REAL*8, intent(in)  :: X(N), Y(N)
  REAL*8, intent(out) :: B(N), C(N), D(N)
  !---------------------------------------------------------
  !     CALCULATE COEFFICIENTS B(I),C(I),D(I),I=1,2,...,N
  !     FOR CUBIC INTERPOLATIONAL SPLINE:
  !     S(X) = Y(I) + B(I) * (X-X(I)) + C(I)*(X-X(I))**2 + D(I)*(X-X(I))**3
  !     FOR X(I).LE.X.LE.X(I+1),
  !     WHERE: N - NUMBER OF POINTS (N.GE.2),
  !     X - ABSCISES (IN STRICTLY INCREASING ORDER),
  !     Y - ORDINATES.
  !---------------------------------------------------------
  !     
  INTEGER :: NM1, IB, I
  REAL*8  :: T
  !     
  NM1=N-1
  IF(N .LT. 2) RETURN
  IF(N .LT. 3) GO TO 50
  !---------------------------------------------------------
  !     CONSTRUCT 3-DIAGONAL SYSTEM WITH B - DIAGONAL,
  !     D - OVERDIAGONAL, C - RIGHT SIDES
  !---------------------------------------------------------
  !     
  D(1)=X(2)-X(1)
  C(2)=(Y(2)-Y(1))/D(1)
  DO I=2,NM1
     D(I)=X(I+1)-X(I)
     B(I)=2.d0*(D(I-1)+D(I))
     C(I+1)=(Y(I+1)-Y(I))/D(I)
     C(I)=C(I+1)-C(I)
  ENDDO
  !----------------------------------------------------------
  !     BOUNDARY CONDITIONS
  !----------------------------------------------------------
  !     
  B(1)= -D(1)
  B(N)= -D(N-1)
  C(1)=0.D00
  C(N)=0.D00
  IF(N .EQ. 3) GO TO 15
  C(1)=C(3)/(X(4)-X(2))-C(2)/(X(3)-X(1))
  C(N)=C(N-1)/(X(N)-X(N-2))-C(N-2)/(X(N-1)-X(N-3))
  C(1)=C(1)*D(1)*D(1)/(X(4)-X(1))
  C(N)= -C(N)*D(N-1)*D(N-1)/(X(N)-X(N-3))
  !----------------------------------------------------------
  !     UPWARD STEP
  !----------------------------------------------------------
  !     
15 CONTINUE
  DO I=2,N
     T=D(I-1)/B(I-1)
     B(I)=B(I)-T*D(I-1)
     C(I)=C(I)-T*C(I-1)
  ENDDO
  !----------------------------------------------------------
  !     BACKWARD SUBSTITUTION
  !----------------------------------------------------------
  !     
  C(N)=C(N)/B(N)
  DO IB =1,NM1
     I=N-IB
     C(I)=(C(I)-D(I)*C(I+1))/B(I)
  ENDDO
  !----------------------------------------------------------
  !     CALCULATE COEFFICIENTS OF THE POLINOMIALS
  !----------------------------------------------------------
  !     
  B(N)=(Y(N)-Y(NM1))/D(NM1)+D(NM1)*(C(NM1)+2.D00*C(N))
  DO I=1,NM1
     B(I)=(Y(I+1)-Y(I))/D(I)-D(I)*(C(I+1)+2.D00*C(I))
     D(I)=(C(I+1)-C(I))/D(I)
     C(I)=3.D00*C(I)
  ENDDO
     
  C(N)=3.D00*C(N)
  D(N)=D(N-1)
  RETURN
!     
50 CONTINUE
  B(1)=(Y(2)-Y(1))/(X(2)-X(1))
  C(1)=0.D00
  D(1)=0.D00
  B(2)=B(1)
  C(2)=0.D00
  D(2)=0.D00
  !     
  RETURN
END SUBROUTINE SPLINE
