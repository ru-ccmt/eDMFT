! @Copyright 2007 Kristjan Haule
!*****************************************************************!
! This module implements exchange and LDA-correlation
!  functional for Coulomb interaction screened by
!  yukawa form or dielectric form:
!    2*exp(-lambda*r)/(eps*r)
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
!
!  subroutine CorrLDA_2(Ec, Vc, rs, lambda, eps, N)
!      needs:
!          rs(N), lambda, eps
!      returns:
!          Ec(N), Vc(N)
!
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
  if (abs(lambda).gt.1e-10) then
     !print *, 'lambda too small'
     xs(:) = rs_1(:)*(kf_rs/lambda)       ! kf/(rs*lambda)
     CALL fexchange(xs,Ex,dEx,N)
  else
     !print *, 'lambda zero', lambda, abs(lambda).gt.1e-10
     xs(:) = rs_1(:)*kf_rs                ! kf/(rs*lambda)
     Ex(:)=1.0
     dEx(:)=0.0
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
  integer, volatile :: iCopy
  REAL*8 :: x, x2, x3, at2, lg2
  do i=1,N
     x = xs(i)
     x2 = x*x
     x3 = x2*x
     iCopy = i  ! A fix for intel compiler bug, so that it does not optimize out the loop
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
  !if (ISNAN(GG)) then
  !print *, 'rs=', rs, 'RS12=', RS12, 'RS32=', RS32, 'RSP=', RSP, 'Q1=', Q1, 'Q2=', Q2, 'GG=', GG, 'GGRS=', GGRS
  !endif
  !print *, 'rs=', rs, Ec, Vc
end subroutine CorLDA

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



subroutine EcVc_reduce_di_2(rs,eps,fn,qn,N)
  IMPLICIT NONE
  INTEGER, intent(in):: N
  REAL*8, intent(in) :: rs(N)
  REAL*8, intent(in) :: eps
  REAL*8, intent(out):: fn(N), qn(N)
  ! locals
  INTEGER:: i
  REAL*8 :: a, b, d, eps_a, eps_d, eps_abd, b2r_9, r_da_dr, r_db_dr, r_dd_dr
  REAL*8 :: ca(3), cb(3), cd(4)
  ! The coefficients a, b, d depend on rs, and here is their more precise fit:
  ca = (/1.74596971, -0.0892907,   0.00658866/)
  cb = (/ 1.63289109,  1.15291480, 0.149402/)
  cd = (/3.64370598, 0.03636027, -0.03886317, 0.00693599/)

  do i=1,N
     if (rs(i) .gt. 40) then
        qn(i)=0.0
        fn(i)=0.0
        CYCLE
     endif
     a = ca(1) + rs(i) * ca(2) + rs(i)**2 * ca(3)
     b = 0.001*(cb(1)*sqrt(rs(i))+cb(2)*rs(i))/(1+ (cb(3)*rs(i))**9)
     d = cd(1) + rs(i) * cd(2) + rs(i)**2 * cd(3) + rs(i)**3 * cd(4)
     eps_a = eps**a
     eps_d = eps**d
     eps_abd = eps_a+b*eps_d
     b2r_9 = (cb(3)*rs(i))**9
     fn(i) = (1.+b)/eps_abd
     r_da_dr = rs(i) * ca(2) + 2*rs(i)**2 * ca(3)
     r_dd_dr = rs(i) * cd(2) + 2*rs(i)**2 * cd(3) + 3*rs(i)**3 * cd(4)
     r_db_dr = 0.001*(cb(1)*sqrt(rs(i))*(0.5 - 17./2.*b2r_9)+cb(2)*rs(i)*(1-b2r_9*8))/(1+ b2r_9)**2
     qn(i) = -1./3.*(r_db_dr*(eps_a-eps_d)-(eps_a*r_da_dr + b*eps_d*r_dd_dr)*(1+b)*log(eps))/eps_abd**2
  enddo
  ! Ec = Ev * fn
  ! Vc = Vc * fn + Ec * qn
end subroutine EcVc_reduce_di_2

subroutine EcVc_reduce_yw_2(rs,lambda,fn,qn,N)
  IMPLICIT NONE
  INTEGER, intent(in):: N
  REAL*8, intent(in) :: rs(N)
  REAL*8, intent(in) :: lambda
  REAL*8, intent(out):: fn(N), qn(N)
  ! locals
  INTEGER:: i, m
  REAL*8 :: lmrs, dlms, te, das
  REAL*8 :: C(7,6)
  REAL*8 :: an(6)
  C = Reshape((/0.15805009, -0.77391602, 1.23971169, -1.04865383, 0.47809619, -0.11057964, 0.01016968,&
               -0.306851, -0.77296572, 0.8791705, -0.69185034, 0.33779654, -0.08858483, 0.00935635,&
               0.13215843, -0.2776552, 0.45727548, -0.31469164, 0.10787374, -0.01661214, 0.0007591,&
               -0.03086548, 0.0549528, -0.07252823, 0.04177618, -0.01084882, 0.00062192, 0.0001177,&
               0.00273230889, -0.00357007233, 0.00425309814, -0.00198811211, 0.000233761378, 0.000106803015, -2.50612307e-05,&
               -9.28530649e-05, 8.09009085e-05, -9.43747991e-05, 3.89520548e-05, -3.10149723e-07, &
               -4.23041605e-06, 8.02291467e-07/),(/7,6/))
  do m=1,6
     an(m) = lambda * (c(1,m) + lambda * (c(2,m) + lambda * (c(3,m) + lambda * (c(4,m) + lambda * (c(5,m) + lambda * (c(6,m) &
             + lambda*c(7,m)))))))
  enddo
  !print *, an

  do i=1,N
     te = exp(an(1)+rs(i)*(an(2)+rs(i)*(an(3)+rs(i)*(an(4)+rs(i)*(an(5)+rs(i)*an(6))))))
     lmrs = 0.008 - 0.00112 * rs(i)**2
     dlms = -2*0.00112*rs(i)
     fn(i) = te*(1-lmrs) + lmrs
     das = rs(i)*(an(2) + rs(i)*(2*an(3) + rs(i)*(3*an(4) + rs(i)*(4*an(5) + rs(i)*5*an(6)))))
     qn(i) = -1./3.*te*(1-lmrs)*das - 1./3.*rs(i)*(1-te)*dlms
  enddo
  ! Ec = Ev * fn
  ! Vc = Vc * fn + Ec * qn
end subroutine EcVc_reduce_yw_2

subroutine EcVc_reduce_yw_2_bad(rs,lambda,fn,qn,N)
  IMPLICIT NONE
  INTEGER, intent(in):: N
  REAL*8, intent(in) :: rs(N)
  REAL*8, intent(in) :: lambda
  REAL*8, intent(out):: fn(N), qn(N)
  ! locals
  INTEGER:: i, m
  REAL*8 :: lmrs, dlms, te, das
  REAL*8 :: C(6,5)
  REAL*8 :: an(5)
  C = Reshape((/0.09285467, -0.47775779, 0.56979491, -0.3400053, 0.09811778, -0.01084068,&
            -2.71279489e-01, -6.09031829e-01, 4.13668946e-01, -1.37495480e-01, 1.91211856e-02, -5.07505731e-04,&
             0.07712942,-0.19993112,0.32952836,-0.21156433,0.06355575,-0.00726003,&
             -0.01562657,0.04090609,-0.05576869,0.03415515,-0.01021151,0.00117334,&
             7.84146655e-04,-2.02700552e-03,2.70487819e-03,-1.67140617e-03,5.12525527e-04,-6.01123272e-05/),(/6,5/))
  do m=1,5
     an(m) = lambda * (c(1,m) + lambda * (c(2,m) + lambda * (c(3,m) + lambda * (c(4,m) + lambda * (c(5,m) + c(6,m) * lambda)))))
  enddo
  !print *, an

  do i=1,N
     te = exp(an(1)+rs(i)*(an(2)+rs(i)*(an(3)+rs(i)*(an(4)+an(5)*rs(i)))))
     lmrs = 0.008 - 0.00112 * rs(i)**2
     dlms = -2*0.00112*rs(i)
     fn(i) = te*(1-lmrs) + lmrs
     das = rs(i)*(an(2) + rs(i)*(2*an(3) + rs(i)*(3*an(4) + 4*an(5)*rs(i) )))
     qn(i) = -1./3.*te*(1-lmrs)*das - 1./3.*rs(i)*(1-te)*dlms
  enddo
  ! Ec = Ev * fn
  ! Vc = Vc * fn + Ec * qn
end subroutine EcVc_reduce_yw_2_bad

subroutine CorrLDA_2(Ec, Vc, rs, lambda, eps, N)
  IMPLICIT NONE
  INTEGER, intent(in):: N
  REAL*8, intent(in) :: rs(N)
  REAL*8, intent(in) :: lambda, eps
  REAL*8, intent(out):: Ec(N), Vc(N)
  !locals
  INTEGER :: i
  REAL*8 :: fn_l(N), qn_l(N), fn_e(N), qn_e(N), fn(N), qn(N)
  
  do i=1,N
     CALL CorLDA(Ec(i),Vc(i),rs(i))
  enddo
  CALL EcVc_reduce_yw_2(rs,lambda,fn_l,qn_l,N)
  !print *, 'l=', lambda, 'fn_l=', fn_l
  !print *, 'rs=', rs, 'qn_l=', qn_l
  CALL EcVc_reduce_di_2(rs,eps,fn_e,qn_e,N)
  !print *, 'e=', eps, 'fn_e=', fn_e
  !print *, 'rs=', rs, 'qn_e=', qn_e

  fn(:) = fn_l(:)*fn_e(:)
  qn(:) = qn_l(:)*fn_e(:)+fn_l(:)*qn_e(:)
  
  Vc(:) = Vc(:)*fn(:) + Ec(:)*qn(:)
  Ec(:) = Ec(:)*fn(:)
end subroutine CorrLDA_2

!program exc
!  IMPLICIT NONE
!  REAL*8 :: lambda, eps
!  INTEGER:: i
!  REAL*8 :: rx,  dr
!  REAL*8 :: Ecs(4), Vcs(4), Ecs1(4), Vcs1(4), Ecs2(4), Vcs2(4)
!  REAL*8 :: A0(4), C0(4), A1(4), C1(4), A2(4), C2(4), dA(4), dE(4), E0(4), Ex(4), Vx(4)
!  REAL*8 :: Ec(9), Vc(9), Ec1(9), Ec2(9), Vc1(9), Vc2(9)
!  REAL*8 :: fn(9), qn(9)
!  REAL*8 :: rs(9)
!  rs=(/0.5,1.,2.,3.,4.,5.,6.,7.,8./)
!  lambda=0.5
!  eps=1.5
!  
!  dr=0.0001
!  
!  !CALL EcVc_reduce(rs,lambda,A,C,size(rs,1))
!  CALL EcVc_reduce_yw_2(rs,lambda,fn,qn,size(rs,1))
!  print *, fn
!  print *, qn
!  CALL EcVc_reduce_di_2(rs,eps,fn,qn,size(rs,1))
!  
!
!  CALL CorrLDA_2(Ec,Vc,rs,lambda,eps,size(rs,1))
!  CALL CorrLDA_2(Ec1,Vc1,rs+dr,lambda,eps,size(rs,1))
!  CALL CorrLDA_2(Ec2,Vc2,rs-dr,lambda,eps,size(rs,1))
!   
!  do i=1,size(rs,1)
!     rx=rs(i)
!     Vc1(i) = Ec(i) - 1./3.*rx*(Ec1(i)-Ec2(i))/(2*dr)
!     print *, 'rx=', rx, 'Vc=', Vc(i), 'Ec=', Ec(i), 'Vapprox=', Vc1(i)
!  enddo
!
!end program exc
