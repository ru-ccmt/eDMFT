Function Force4(Clm_new,Vlm,R,lmmax,lm,jri,dx,nrad,lmax2,ncom) result(fvdrho)
  ! ************ Calculates the following three dimensional integral ****************
  !     Int[  V_{ks} * (nabla)rho ]
  !
  !     Both V_{ks} and rho are given in terms of spherical harmonics expansion
  !     The integral over angle is performed analytically, while the radial integral
  !     is performed numerically.
  !
  !    For derivation of the formula, see: Computer Physics Communications 94, 31-48 (1996).
  ! *********************************************************************************
  IMPLICIT NONE
  !use param, ONLY: lmax2, ncom, nrad
  REAL*8 :: fvdrho(3)
  REAL*8, intent(in)  :: Clm_new(nrad,0:ncom-1)  ! density rho(r)
  REAL*8, intent(in)  :: Vlm(nrad,0:ncom-1)     ! potential V(r)
  REAL*8, intent(in)  :: R(nrad), dx
  INTEGER, intent(in) :: lm(2,0:NCOM-1), lmax2, nrad, ncom, lmmax, jri
  ! functions
  Interface
     real*8 FUNCTION int_radial(func,RX,DELTA,nr)
       INTEGER, intent(in):: nr
       REAL*8, intent(in) :: RX(nr), funct(nr), DELTA
     end FUNCTION int_radial
  end Interface
  ! locals
  REAL*8 :: intfeld1(0:lmax2,0:lmax2,-1:1,0:lmax2,0:lmax2,-1:1)
  REAL*8 :: intfeld2(0:lmax2,0:lmax2,-1:1,0:lmax2,0:lmax2,-1:1) 
  REAL*8  :: toint(NRAD),derv(NRAD)
  INTEGER :: l1, m1, l2, m2, mz, ik1, ik2, lm1, lm2, llmax, L, M
  REAL*8 :: VFA, VFC, VFF, fakt, faktmo, faktpo, forvrx, forvry, forvrz, sqhalf
  !*******************Statement Functions***************************
  VFA(L,M) = SQRT((L+M+1.)*(L+M+2.)/((2.*L+3.)*(2.*L+1.)))
  !VFC(L,M) = SQRT((L-M+1.)*(L-M+2.)/((2.*L+3.)*(2.*L+1.)))
  VFF(L,M) = SQRT((L+M+1.)*(L-M+1.)/((2.*L+3.)*(2.*L+1.)))
  !**************************************************************
  SQHALF=SQRT(0.5D0) 
  !***************** Integrals of V*nabla*rho ********  
  forvrz=0.D0
  forvry=0.D0
  forvrx=0.D0
  intfeld1(0:lmax2,0:lmax2,-1:1,0:lmax2,0:lmax2,-1:1)=0.D0   ! Int( V_{l1,m1}(r) * drho_{l2,m2}(r)/dr r**2 dr)
  intfeld2(0:lmax2,0:lmax2,-1:1,0:lmax2,0:lmax2,-1:1)=0.D0   ! Int( V_{l1,m1}(r) * rho_{l2,m2}(r) / r r**2 dr)

  llmax=iabs(lm(1,lmmax-1))
  DO lm1=0,lmmax-2 
     l1=abs(lm(1,lm1)) 
     m1=lm(2,lm1)
     ik1 = sign(1, lm(1,lm1) )
     DO lm2=lm1+1,lmmax-1
        l2 = abs(lm(1,lm2))
        m2 = lm(2,lm2)
        ik2 = sign(1, lm(1,lm2) )
        if ( (l2 .eq. (l1+1)) .and. ( (m2.eq.m1) .or. (m2.eq.(m1-1)) .or. (m2.eq.(m1+1)) )) then
           ! Computes intfeld2(l1,m1,l2,m2) = Integrate( V_{lm1}(r)*rho_{lm2}(r)/r * r**2 dr)
           !    Note that rho(r)*r**2 = Clm(r)
           !    intfeld2(l1,m1,l2,m2) = Int( V_{lm}(r) * (r**2*rho(r)_{l+1,m+(0,1,-1)})/r )
           toint(:jri) = Vlm(:jri,lm1) * Clm_new(:jri,lm2) / R(:jri)
           intfeld2(l1,m1,ik1,l2,m2,ik2) = int_radial(toint(:jri),r(:jri),dx,jri)
           !    Here we just replace lm1 <-> lm2 so that the loop is only over l2>l1
           toint(:jri) = Vlm(:jri,lm2) * Clm_new(:jri,lm1) / R(:jri)
           intfeld2(l2,m2,ik2,l1,m1,ik1) = int_radial(toint(:jri),r(:jri),dx,jri)
           !*************************************************************************
           ! Computes intfeld1(l1,m1,l2,m2) = Integrate( V_{lm1}(r)* d/dr rho_{lm2}(r) * r**2 dr)
           !    Note that rho(r)*r**2 = Clm(r) hence  d/dr ( Clm/r**2) = ( dClm/dr - 2 Clm/r )/r**2 and
           !    intfeld1(l1,m1,l2,m2) = Integrate( V_{lm1}(r)* ( dClm/dr - 2 Clm/r ) dr)
           call dfrad(derv(:jri),Clm_new(:jri,lm2),r(:jri),jri)
           toint(:jri) = Vlm(:jri,lm1) * ( derv(:jri) - 2.D0*Clm_new(:jri,lm2)/R(:jri) )
           intfeld1(l1,m1,ik1,l2,m2,ik2) = int_radial(toint(:jri),r(:jri),dx,jri)
           !    Here we just replace lm1 <-> lm2 so that the loop is only over l2>l1 
           call dfrad(derv(:jri),Clm_new(:jri,lm1),r(:jri),jri)
           toint(:jri) = Vlm(:jri,lm2) * ( derv(:jri) - 2.D0*Clm_new(:jri,lm1)/R(:jri) )
           intfeld1(l2,m2,ik2,l1,m1,ik1) = int_radial(toint(:jri),r(:jri),dx,jri)
           !*************************************************************************
        endif
     ENDDO
  ENDDO
  DO l=0,llmax-1
     DO mz=-l,l
        m=abs(mz)
        fakt=1.d0
        if (m.eq.0) fakt=2.d0
        ! Fz = SQRT((L+M+1)*(L-M+1)/((2*L+3)*(2*L+1)))*0.5*(
        !      (2+l)* V_{l,m,+}*rho_{l+1,m,+}/r
        !      (2+l)* V_{l,m,-}*rho_{l+1,m,-}/r
        !        -l * V_{l+1,m,+}*rho_{l,m,+}/r
        !        -l * V_{l+1,m,-}*rho_{l,m,-}/r
        !             V_{l,m,+} * (drho/dr)_{l+1,m,+}
        !             V_{l,m,-} * (drho/dr)_{l+1,m,-}
        !             V_{l+1,m,+} * (drho/dr)_{l,m,+}
        !             V_{l+1,m,-} * (drho/dr)_{l,m,-}
        !      )
        forvrz = forvrz + 0.5d0*fakt*vff(l,mz)*( (2+l)*( intfeld2(l,  m,+1,l+1,m,+1) + intfeld2(l,  m,-1,l+1,m,-1)) &
                                                    -l*( intfeld2(l+1,m,+1,l,  m,+1) + intfeld2(l+1,m,-1,l,  m,-1)) &
                                                      +( intfeld1(l,  m,+1,l+1,m,+1) + intfeld1(l,  m,-1,l+1,m,-1)) &
                                                      +( intfeld1(l+1,m,+1,l,  m,+1) + intfeld1(l+1,m,-1,l,  m,-1)))
        
        if (mz.gt.0)    faktpo = -1.0d0
        if (mz.lt.(-1)) faktpo =  1.0d0
        if (mz.eq.0)    faktpo = -sqrt(2.0d0)
        if (mz.eq.(-1)) faktpo =  sqrt(2.0d0)

        if (mz.gt.1) faktmo = -1.0d0
        if (mz.lt.0) faktmo =  1.0d0
        if (mz.eq.0) faktmo =  sqrt(2.0d0)
        if (mz.eq.1) faktmo = -sqrt(2.0d0)
        !
        ! Fx = SQRT((L+M+1)*(L+M+2)/((2*L+3)*(2*L+1)))*0.5*(
        !      (2+l)* V_{l,|m|,+} * rho_{l+1,|mz+1|,+}/r
        !      (2+l)* V_{l,|m|,-} * rho_{l+1,|mz+1|,-}/r
        !       -l  * V_{l+1,|m+1|,+} * rho_{l,|m|,+}/r
        !       -l  * V_{l+1,|m+1|,-} * rho_{l,|m|,-}/r
        !             V_{l,|m|,+} * (drho/dr)_{l+1,|m+1|,+}
        !             V_{l,|m|,-} * (drho/dr)_{l+1,|m+1|,-}
        !             V_{l+1,|m+1|,+} * (drho/dr)_{l,|m|,+}
        !             V_{l+1,|m+1|,-} * (drho/dr)_{l,|m|,-}
        !      )
        ! Fx += SQRT((L-M+1)*(L-M+2)/((2*L+3)*(2*L+1)))*0.5*(
        !       (2+l)* V_{l,|m|,+} * rho_{l+1,|m-1|,+}
        !       (2+l)* V_{l,|m|,-} * rho_{l+1,|m-1|,-}
        !         -l * V_{l+1,|m-1|,+} * rho_{l,|m|,+}
        !         -l * V_{l+1,|m-1|,-} * rho_{l,|m|,-}
        !              V_{l,|m|,+} * (drho/dr)_{l+1,|m-1|,+}
        !              V_{l,|m|,-} * (drho/dr)_{l+1,|m-1|,-}
        !              V_{l+1,|m-1|,+} * (drho/dr)_{l,|m|,+}
        !              V_{l+1,|m-1|,-} * (drho/dr)_{l,|m|,-}
        !          )
        forvrx = forvrx - vfa(l, mz)*0.25D0*faktpo*(  (2+l)*( intfeld2(l,  abs(mz),  +1,l+1,abs(mz+1),+1) + intfeld2(l,  abs(mz),  -1,l+1,abs(mz+1),-1)) &
                                                         -l*( intfeld2(l+1,abs(mz+1),+1,l,  abs(mz),  +1) + intfeld2(l+1,abs(mz+1),-1,l,  abs(mz),  -1)) &
                                                           +( intfeld1(l,  abs(mz),  +1,l+1,abs(mz+1),+1) + intfeld1(l,  abs(mz),  -1,l+1,abs(mz+1),-1)) &
                                                           +( intfeld1(l+1,abs(mz+1),+1,l,  abs(mz),  +1) + intfeld1(l+1,abs(mz+1),-1,l,  abs(mz),  -1)))
        forvrx = forvrx + vfa(l,-mz)*0.25D0*faktmo*(  (2+l)*( intfeld2(l,  abs(mz),  +1,l+1,abs(mz-1),+1) + intfeld2(l,  abs(mz),  -1,l+1,abs(mz-1),-1)) &
                                                         -l*( intfeld2(l+1,abs(mz-1),+1,l,  abs(mz),  +1) + intfeld2(l+1,abs(mz-1),-1,l,  abs(mz),  -1)) &
                                                           +( intfeld1(l,  abs(mz),  +1,l+1,abs(mz-1),+1) + intfeld1(l,  abs(mz),  -1,l+1,abs(mz-1),-1)) &
                                                           +( intfeld1(l+1,abs(mz-1),+1,l,  abs(mz),  +1) + intfeld1(l+1,abs(mz-1),-1,l,  abs(mz),  -1)))
        faktmo = 1.0D0
        if (mz.eq.1 .or. mz.eq.0) faktmo = sqrt(2.0d0)
        faktpo = 1.0d0
        if (mz.eq.0 .or. mz.eq.(-1)) faktpo = sqrt(2.0d0)
        forvry = forvry + vfa(l, mz)*0.25D0*faktpo*( (2+l)*( intfeld2(l,  abs(mz),  +1,l+1,abs(mz+1),-1) - intfeld2(l,  abs(mz),  -1,l+1,abs(mz+1),+1)) &
                                                        +l*( intfeld2(l+1,abs(mz+1),+1,l,  abs(mz),  -1) - intfeld2(l+1,abs(mz+1),-1,l,  abs(mz),  +1)) &
                                                          +( intfeld1(l,  abs(mz),  +1,l+1,abs(mz+1),-1) - intfeld1(l,  abs(mz),  -1,l+1,abs(mz+1),+1)) &
                                                          -( intfeld1(l+1,abs(mz+1),+1,l,  abs(mz),  -1) - intfeld1(l+1,abs(mz+1),-1,l,  abs(mz),  +1)))
        forvry = forvry + vfa(l,-mz)*0.25D0*faktmo*( (2+l)*( intfeld2(l,  abs(mz),  +1,l+1,abs(mz-1),-1) - intfeld2(l,  abs(mz),  -1,l+1,abs(mz-1),+1)) &
                                                        +l*( intfeld2(l+1,abs(mz-1),+1,l,  abs(mz),  -1) - intfeld2(l+1,abs(mz-1),-1,l,  abs(mz),  +1)) &
                                                          +( intfeld1(l,  abs(mz),  +1,l+1,abs(mz-1),-1) - intfeld1(l,  abs(mz),  -1,l+1,abs(mz-1),+1)) &
                                                          -( intfeld1(l+1,abs(mz-1),+1,l,  abs(mz),  -1) - intfeld1(l+1,abs(mz-1),-1,l,  abs(mz),  +1)) )
     ENDDO
  ENDDO
  fvdrho(1)=forvrx
  fvdrho(2)=forvry
  fvdrho(3)=forvrz
  return
end function Force4



