! @Copyright 2007 Kristjan Haule
!*****************************************************************!
! This module implements exchange and LDA-correlation
! functional in case of yukawa screened interaction
!    2*exp(-lambda*r)/r
! Here units are implemented in Rydbergs.
! The relevant subroutines are:
!
!   subroutine ExchangeLDA(Ex, Vx, rs_1, lambda, N)
!      needs:
!         real*8 :: Ex(N), Vx(N), rs_1(N), lambda  ! rs_1==1/rs
!
!   subroutine CorrLDA(Ec, Vc, rsi, lambda, N)
!      needs:
!         real*8 :: Ec(N_, Vc(N), rsi(N), lambda  ! rsi=rs
!*****************************************************************!
subroutine ExchangeLDA(Ex, Vx, rs_1, lambda, N)
  IMPLICIT NONE
  REAL*8, intent(in) :: lambda
  REAL*8, intent(in) :: rs_1(N)
  REAL*8, intent(out):: Ex(N), Vx(N)
  INTEGER, intent(in):: N
  ! locals
  REAL*8 :: kf_rs, c0, pi
  REAL*8 :: xs(N), dEx(N)
  !
  pi = 4.*atan(1.)
  !
  kf_rs = (9*pi/4.)**(1./3.)
  c0 = (3./(2.*pi))*(9.*pi/4.)**(1./3.)
  if (lambda.ne.0) then
     xs(:) = rs_1(:)*(kf_rs/lambda)       ! kf/(rs*lambda)
     CALL fexchange(xs,Ex,dEx,N)
  else
     Ex(:)=1.0
     dEx(:)=0.0
     xs(:)=0.0
  endif
  Ex(:) = -(c0*rs_1)*Ex(:)
  Vx(:) = 4./3.*Ex(:) - (c0/3.)*rs_1(:)*xs(:)*dEx(:)
end subroutine ExchangeLDA

SUBROUTINE fexchange(xs,ex,dex,N)
  !* Evaluating a function (and its derivative) for exchange energy, which has a formula
  !*  f(x) = 1-1/(6*x^2)-4/(3*x)*atan(2*x)+(1+1/(12*x^2))/(2*x^2)*log(1+4*x^2)
  !* df/dx = 2/(3*x^3) + 4/(3*x^2)*atan(2*x) - (1+6*x^2)/(6*x^5)*log(1+4*x^2)
  IMPLICIT NONE
  REAL*8, intent(in) :: xs(N)
  REAL*8, intent(out):: ex(N), dex(N)
  INTEGER, intent(in):: N
  ! locals
  INTEGER :: i
  REAL*8 :: x, x2, x3, at2, lg2
  do i=1,N
     x = xs(i)
     x2 = x*x
     x3 = x2*x
     if (x<0.01) then
        ex(i) = 4.*x2/9.*(1-6*x2/5.)
        dex(i) = 8.*x/9.*(1-12*x2/5.)
     else
        at2 = atan(2*x)*4./(3.*x)
        lg2 = log(1+4.*x2)/(2.*x2)
        ex(i) = 1-1/(6.*x2) - at2+ (1.+1/(12.*x2))*lg2
        dex(i) = 2./(3.*x3) + at2/x - (1.+6*x2)/(3.*x3)*lg2
     end if
  end do
END SUBROUTINE fexchange
  
REAL*8 Function f1(coef,x,x2,x4,x6)
  IMPLICIT NONE
  REAL*8, intent(in) :: x, x2, x4, x6, coef(5)
  f1 = (x*coef(1)+x2*coef(2))/(1+x2*coef(3)+x4*coef(4)+x6*coef(5))
  return
END Function f1
REAL*8 Function f2(coef,x2,x3,x4)
  IMPLICIT NONE
  REAL*8, intent(in) :: x2, x3, x4, coef(4)
  ! locals
  f2 = (x2*coef(1)+x3*coef(2))/(1+x2*coef(3)+x4*coef(4))
  return
END Function f2
REAL*8 Function f3(coef,x2,x3,x4)
  IMPLICIT NONE
  REAL*8, intent(in) :: x2, x3, x4, coef(3)
  f3 = (x3*coef(1)+x4*coef(2))/(1+coef(3)*x2)
  return
END Function f3
REAL*8 Function f4(coef,x2,x4)
  IMPLICIT NONE
  REAL*8, intent(in) :: x2, x4, coef(2)
  f4 = x4*(coef(1)+coef(2)*x2)
  return
END Function f4

subroutine EcVc_reduce(rsi,lambda,A,C,N)
  IMPLICIT NONE
  INTEGER, intent(in):: N
  REAL*8, intent(in) :: rsi(N)
  REAL*8, intent(in) :: lambda
  REAL*8, intent(out):: A(N), C(N)
  !
  interface
     REAL*8 Function f1(coef,x,x2,x4,x6)
       REAL*8, intent(in) :: x, x2, x4, x6, coef(5)
     END Function f1
     REAL*8 Function f2(coef,x2,x3,x4)
       REAL*8, intent(in) :: x2, x3, x4, coef(4)
     END Function f2
     REAL*8 Function f3(coef,x2,x3,x4)
       REAL*8, intent(in) :: x2, x3, x4, coef(3)
     END Function f3
     REAL*8 Function f4(coef,x2,x4)
       REAL*8, intent(in) :: x2, x4, coef(2)
     END Function f4
  end interface
  ! locals
  INTEGER:: i
  REAL*8 :: coef1(5), coef2(4), coef3(3), coef4(2)
  REAL*8 :: cfs(4)
  REAL*8 :: x, x2, x3, x4, x6, rs, Br
  coef1 = (/0.12238912,  0.73648662,  0.96044695, -0.07501634,  0.00207808/)
  coef2 = (/0.05839362,  0.11969474,  0.10156124,  0.01594125/)
  coef3 = (/0.00827519,  0.00557133,  0.01725079/)
  coef4 = (/5.29134419e-04, 4.49628225e-06/)

  x=lambda
  x2=x*x
  x3=x*x2
  x4=x2*x2
  x6=x4*x2

  cfs(1)=exp(f1(coef1,x,x2,x4,x6))-1.
  cfs(2)=exp(f2(coef2,x2,x3,x4))-1.
  cfs(3)=exp(f3(coef3,x2,x3,x4))-1.
  cfs(4)=exp(f4(coef4,x2,x4))-1.

  !print *, 'cfs=', cfs
  do i=1,N
     rs=rsi(i)
     A(i)  = 1.0+rs*(cfs(1)+rs*(cfs(2)+rs*(cfs(3)+rs*cfs(4))))
     Br = rs*(cfs(1)+rs*(2*cfs(2)+rs*(3*cfs(3)+rs*4*cfs(4))))
     C(i) = 3*A(i)**2/Br
  enddo
end subroutine EcVc_reduce

SUBROUTINE  Gcor(GG,GGrs,rs, A,A1,B1,B2,B3,B4,P)
  IMPLICIT NONE
  REAL*8, intent(in)  :: rs, A, A1, B1, B2, B3, B4
  INTEGER, intent(in) :: P
  REAL*8, intent(out) :: GG, GGrs
  ! locals
  REAL*8 :: P1, Q0, RS12, RS32, RSP, Q1, Q2, Q3
  P1 = P + 1.
  Q0 = -2.*A*(1.+A1*rs)
  RS12 = sqrt(rs)
  RS32 = RS12**3
  RSP = rs**P
  Q1 = 2.*A*(B1*RS12+B2*rs+B3*RS32+B4*rs*RSP)
  Q2 = log(1.+1./Q1)
  GG = Q0*Q2
  Q3 = A*(B1/RS12+2.*B2+3.*B3*RS12+2.*B4*P1*RSP)
  GGRS = -2.*A*A1*Q2-Q0*Q3/(Q1**2+Q1)
END SUBROUTINE Gcor

subroutine CorLDA(Ec,Vc,rs)
  IMPLICIT NONE
  !  UNITS OF Rydberg
  !  UNIFORM-GAS CORRELATION OF PERDEW AND WANG 1991
  !  INPUT: SEITZ RADIUS (rs)
  !  OUTPUT: CORRELATION ENERGY PER ELECTRON (Ec) and POTENTIALS (Vc)
  !
  REAL*8, intent(out):: Ec, Vc
  REAL*8, intent(in) :: rs
  ! locals
  REAL*8 :: H2Ry, Ecrs
  H2Ry = 2.0    ! Conversion from Hartree to Rydberg
  CALL Gcor(Ec, Ecrs,    rs, 0.0310907D0,  0.21370D0,  7.5957D0, 3.5876D0, 1.6382D0, 0.49294D0, 1)
  Vc = Ec - rs*Ecrs/3.
  Vc = Vc*H2Ry  ! From Hartree to Rydbergs
  Ec = Ec*H2Ry  ! From Hartree to Rydbergs
end subroutine CorLDA

subroutine Corr_LDA(Ec,Vc,rsi,N)
  IMPLICIT NONE
  ! Vector form of CorLDA
  !  UNITS OF Rydberg
  INTEGER, intent(in):: N
  REAL*8, intent(in) :: rsi(N)
  REAL*8, intent(out):: Ec(N), Vc(N)
  ! locals
  INTEGER :: i
  do i=1,N
     CALL CorLDA(Ec(i),Vc(i),rsi(i))
  enddo
end subroutine Corr_LDA

subroutine CorrLDA(Ec, Vc, rsi, lambda, N)
  IMPLICIT NONE
  INTEGER, intent(in):: N
  REAL*8, intent(in) :: rsi(N)
  REAL*8, intent(in) :: lambda
  REAL*8, intent(out):: Ec(N), Vc(N)
  !locals
  INTEGER :: i
  REAL*8 :: A(N), C(N)
  
  do i=1,N
     CALL CorLDA(Ec(i),Vc(i),rsi(i))
  enddo
  
  CALL EcVc_reduce(rsi,lambda,A,C,N)
  Vc(:) = Vc(:)/A(:)+Ec(:)/C(:)
  Ec(:) = Ec(:)/A(:)
  
end subroutine CorrLDA

! program exc
!   IMPLICIT NONE
!   REAL*8 :: rs(4), A(4), C(4), rs_1(4)
!   REAL*8 :: lambda, Ec, Vc
!   INTEGER:: i
!   REAL*8 :: rx, Ec1, Vc1, Ec2, Vc2, dr
!   REAL*8 :: Ecs(4), Vcs(4), Ecs1(4), Vcs1(4), Ecs2(4), Vcs2(4)
!   REAL*8 :: A0(4), C0(4), A1(4), C1(4), A2(4), C2(4), dA(4), dE(4), E0(4), Ex(4), Vx(4)
!   rs=(/1.5,2.5,3.5,4.5/)
!   
!   lambda=1.5
! 
!   dr=0.001
!   
!   CALL EcVc_reduce(rs,lambda,A,C,4)
!   print *, A
!   print *, C
!   
!   do i=1,4
!      rx=rs(i)
!      CALL CorLDA(Ec,Vc,rx)
! 
!      CALL CorLDA(Ec1,Vc1,rx+dr)
!      CALL CorLDA(Ec2,Vc2,rx-dr)
!      
!      Vc1 = Ec - 1./3.*rx*(Ec1-Ec2)/(2*dr)
!      
!      print *, 'rx=', rx, 'Vc=', Vc, 'Ec=', Ec, 'Vapprox=', Vc1
!   enddo
! 
!   CALL  CorrLDA(Ecs1, Vcs1, rs, 0.d0, 4)
!   print *, 'Vc0=', Vcs1
!   print *, 'Ec0=', Ecs1
!   
!   CALL  CorrLDA(Ecs, Vcs, rs, lambda, 4)
!   print *, 'Vc=', Vcs
!   print *, 'Ec=', Ecs
! 
!   print *, 'Vc_ratio=', Vcs/Vcs1
!   
!   !CALL  CorrLDA(Ecs1, Vcs1, rs+dr, lambda, 4)
!   !CALL  CorrLDA(Ecs2, Vcs2, rs-dr, lambda, 4)
!   !Vcs1 = Ecs - 1./3.*rs*(Ecs1-Ecs2)/(2*dr)
!   !print *, 'Vapprox=', Vcs1
! 
! 
! 
!   !CALL EcVc_reduce(rs,lambda,A0,C0,4)
!   !CALL EcVc_reduce(rs+dr,lambda,A1,C1,4)
!   !CALL EcVc_reduce(rs-dr,lambda,A2,C2,4)
!   !
!   !dA = (1/A1(:)-1/A2(:))/(2*dr)
!   !print *, 'dA=', dA
!   !print *, 'dA=', -3./(C0*rs)
! 
!   
!   !CALL  CorrLDA(Ecs1, Vcs1, rs+dr, lambda, 4)
!   !CALL  CorrLDA(Ecs2, Vcs2, rs-dr, lambda, 4)
!   !dE = (Ecs1-Ecs2)/(2*dr)
!   !print *, 'dE/dr=', dE
!   !
!   !do i=1,4
!   !   CALL CorLDA(Ec,Vc,rs(i))
!   !   CALL CorLDA(Ec1,Vc1,rs(i)+dr)
!   !   CALL CorLDA(Ec2,Vc2,rs(i)-dr)
!   !   dE(i) = (Ec1-Ec2)/(2*dr)
!   !   E0(i)=Ec
!   !enddo
!   !CALL EcVc_reduce(rs,lambda,A0,C0,4)
!   !dE(:) = dE(:)/A0(:)-3*E0(:)/(C0*rs)
!   !print *, 'dE/dr=', dE
! 
! 
!   rs_1=1./rs
!   CALL ExchangeLDA(Ex, Vx, rs_1, lambda, 4)
!   print *, 'Ex=', Ex
!   print *, 'Vx=', Vx
!   
! end program exc
