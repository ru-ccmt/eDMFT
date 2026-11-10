subroutine kin_str_integrals(REL, rrad1,rrad2,rade1,rade2,a1lo,b1lo,jatom,imax,dx,VR)
!! The following integrals are used to program the kinetic stress tensor, equation 6.71 at the moment. The integrals dependent only in the atomic index and corresponding, therefore
! it is execute inside the atomic index loop but outside of k-points loop. This subroutine is being called from atpar.F (at the moment).
!TC1S00  = integral with double derivative           r2*(u8ddu)  for LOs juddu
!TC2S00  = integral with single derivative           r1*(u*du)   for LOs judu
!TC3S00  = integral cont product of two functions.   (u*u)       for LOs juu
   use param
   use struk,     only: r0
   use lo,        only: ilo,loor
!    use radial_functions, only: 
   USE stress_kb, only:  tc1s00,tc1s01,tc1s10,tc1s11,tc2s00,tc2s01,&
                         tc2s10,tc2s11,tc3s00,tc3s01,tc3s10,tc3s11, &   
                        juddu_02,juddu_20,juddu_12,juddu_21,juddu_22, &
                        judu_02, judu_20, judu_12, judu_21, judu_22,  &
                        juu_02,  juu_20,  juu_12,  juu_21,  juu_22
   IMPLICIT NONE
   logical             :: REL
   integer,intent(in)  :: jatom,imax
   real*8 ,intent(in)  :: rrad1(nrad,0:lmax2),rrad2(nrad,0:lmax2), &     !before increase the radial points
                          rade1(nrad,0:lmax2),rade2(nrad,0:lmax2), &
                          a1lo(nrad,nloat,0:lomax),b1lo(nrad,nloat,0:lomax), &
                          VR(nrad),dx

   integer             :: i,l,jlo,lp,jlop
   
   real*8              :: uppl2(nrad,0:lmax2) ,upps2(nrad,0:lmax2), &
                          udppl2(nrad,0:lmax2),udpps2(nrad,0:lmax2)

   real*8              :: dul(nrad),dus(nrad),dudl(nrad),duds(nrad), &
                          dulol(nrad),dulos(nrad), ddulol(nrad),ddulos(nrad)
                          
   real*8              :: UAL(nrad),UAS(nrad),UBL(nrad),UBS(nrad), &    !old
                          UALP(nrad),UASP(nrad),UBLP(nrad),UBSP(nrad), &
                          UCL(nrad),UCS(nrad),UCLP(nrad),UCSP(nrad),   &
                          UAB(nrad),UBA(nrad), &
                          UAA(nrad), r1(nrad)                                                     
   real*8              :: CIN,CINP,g0,TC, VALUE(nrad)                          
   
   
   
    CIN = 1.d0/137.0359895d0**2
    CINP = 1.d0/137.0359895d0
    if(.not.REL) CIN = 1E-22 
    
    DO i = 1, imax
       r1(i) = r0(jatom)*exp(dx*(i-1))
    ENDDO
        
    DO lp = 0, lmax2n
       DO i = 1, imax          
          UAL(i) = rrad1(i,lp)/r1(i) 
          UAS(i) = rrad2(i,lp)/r1(i)
          UBL(i) = rade1(i,lp)/r1(i) 
          UBS(i) = rade2(i,lp)/r1(i)           
       ENDDO       
             
       call g1_radial_derivative(imax,r1,UAL,UAA)   !in new dergl last logic is absent              
       call g1_radial_derivative(imax,r1,UAS,UAB)       
       
       call g1_radial_derivative(imax,r1,UAA,UAL)
       call g1_radial_derivative(imax,r1,UAB,UAS)
       
       DO i = 1, imax
          uppl2(i,lp) = UAL(i)*r1(i)
          upps2(i,lp) = UAS(i)*r1(i)
       ENDDO       
       
       call g1_radial_derivative(imax,r1,UBL,UAA)
       call g1_radial_derivative(imax,r1,UBS,UAB)
       
       call g1_radial_derivative(imax,r1,UAA,UBL)
       call g1_radial_derivative(imax,r1,UAB,UBS)
       
       DO i = 1, imax
          udppl2(i,lp) = UBL(i)*r1(i)
          udpps2(i,lp) = UBS(i)*r1(i)
       ENDDO
              
    ENDDO    
    
    
    DO l = 0, lmax2n        
       DO lp = 0, lmax2n          
          
          CALL RINT13(rel,rrad1(1,l),rrad2(1,l),uppl2(1,lp),upps2(1,lp),  TC1S00(l,lp),jatom) !integrand contains just product of two functions: r*r*u_llmd*ddu_llmdp; dd double derivative
          CALL RINT13(rel,rrad1(1,l),rrad2(1,l),udppl2(1,lp),udpps2(1,lp),TC1S01(l,lp),jatom)
          CALL RINT13(rel,rade1(1,l),rade2(1,l),uppl2(1,lp),upps2(1,lp),  TC1S10(l,lp),jatom)
          CALL RINT13(rel,rade1(1,l),rade2(1,l),udppl2(1,lp),udpps2(1,lp),TC1S11(l,lp),jatom)                                                 
          
          DO I = 1, Imax
             UAL(i) = rrad1(i,lp)/r1(i) 
             UAS(i) = rrad2(i,lp)/r1(i)
             UBL(i) = rade1(i,lp)/r1(i) 
             UBS(i) = rade2(i,lp)/r1(i) 
             
             UALP(i) = rrad1(i,l)/r1(i) 
             UASP(i) = rrad2(i,l)/r1(i)
             UBLP(i) = rade1(i,l)/r1(i) 
             UBSP(i) = rade2(i,l)/r1(i)              
          ENDDO        
          call g1_radial_derivative(imax,r1,UAL,dul)
          call g1_radial_derivative(imax,r1,UAS,dus)
             
          call g1_radial_derivative(imax,r1,UBL,dudl)
          call g1_radial_derivative(imax,r1,UBS,duds)
          
          CALL RINT13(rel,rrad1(1,l),rrad2(1,l),dul(1),dus(1),  TC2S00(l,lp), jatom)   !integrand contains just product of two functions: r*u_llmd*du_llmdp; d derivative
          CALL RINT13(rel,rrad1(1,l),rrad2(1,l),dudl(1),duds(1),TC2S01(l,lp), jatom)
          CALL RINT13(rel,rade1(1,l),rade2(1,l),dul(1),dus(1),  TC2S10(l,lp), jatom)
          CALL RINT13(rel,rade1(1,l),rade2(1,l),dudl(1),duds(1),TC2S11(l,lp), jatom)                     
           
          CALL RINT13(rel,UALP(1),UASP(1),UAL(1),UAS(1),  TC3S00(l,lp), jatom)   !integrand contains just product of two functions: u_llmd*u_llmdp
          CALL RINT13(rel,UALP(1),UASP(1),UBL(1),UBS(1),  TC3S01(l,lp), jatom)
          CALL RINT13(rel,UBLP(1),UBSP(1),UAL(1),UAS(1),  TC3S10(l,lp), jatom)
          CALL RINT13(rel,UBLP(1),UBSP(1),UBL(1),UBS(1),  TC3S11(l,lp), jatom)                    
          
          DO jlop = 1, ilo(lp)
             IF(.NOT. loor(jlop,lp)) CYCLE                          
             DO jlo = 1, ilo(l)
                IF(.NOT. loor(jlo,l)) CYCLE
                DO i = 1, imax
                   UCLP(i) = a1lo(i,jlo,l)/r1(i)                         !large comp of LO fun ul(e2)
                   UCSP(i) = b1lo(i,jlo,l)/r1(i)                         !small comp of LO fun us(e2) 
                   
                   UCL(i) = a1lo(i,jlop,lp)/r1(i)                         !large comp of LO fun ul(e2)
                   UCS(i) = b1lo(i,jlop,lp)/r1(i)                         !small comp of LO fun us(e2)
                ENDDO
                call g1_radial_derivative(imax,r1,UCL,dulol)         !first derivative of  ul(r,e2)
                call g1_radial_derivative(imax,r1,UCS,dulos)         !first derivative of  uS(r,e2)
                
                call g1_radial_derivative(imax,r1,dulol,ddulol)      !2nd derivative of  ul(r,e2)
                call g1_radial_derivative(imax,r1,dulos,ddulos)      !2nd derivative of  us(r,e2)  
             
                DO i = 1, imax
                   ddulol(i) = r1(i)*ddulol(i)
                   ddulos(i) = r1(i)*ddulos(i)
                ENDDO
                
                CALL RINT13(REL,a1lo(1,jlo,l),b1lo(1,jlo,l),ddulol(1),ddulos(1),juddu_22(jlo,jlop,l,lp),jatom)
                CALL RINT13(REL,a1lo(1,jlo,l),b1lo(1,jlo,l),dulol(1),dulos(1),  judu_22(jlo,jlop,l,lp),jatom)    
                
                DO I = 1, Imax                   
                   UAA(I) = (a1lo(i,jlop,lp)*a1lo(i,jlo,l) + CIN*b1lo(i,jlop,lp)*b1lo(i,jlo,l))/r1(i)/r1(i)
                ENDDO           
                CALL CHARGEL2(r1,1.D0,1.D0,UAA(1),dx,imax,TC);juu_22(jlo,jlop,l,lp) = TC
                                   
!                 CALL RINT13(REL,UCL(1),UCS(1),UCLP(1),UCSP(1),juu_22(jlo,jlop,l,lp),jatom)
             ENDDO
          ENDDO
          
       ENDDO
    ENDDO  
    
    DO l = 0, lmax2n
       DO jlo = 1, ilo(l)
          IF(.NOT. loor(jlo,l)) CYCLE
          DO i = 1, imax
             UCL(i) = a1lo(i,jlo,l)/r1(i)                         !large comp of LO fun ul(e2)
             UCS(i) = b1lo(i,jlo,l)/r1(i)                         !small comp of LO fun us(e2)
          ENDDO
          
          DO lp = 0,lmax2n
             DO i = 1, imax                                                        
                UAL(i) = rrad1(i,lp)/r1(i) 
                UAS(i) = rrad2(i,lp)/r1(i)
                UBL(i) = rade1(i,lp)/r1(i) 
                UBS(i) = rade2(i,lp)/r1(i) 
            ENDDO                  
            call g1_radial_derivative(imax,r1,UAL,dul)
            call g1_radial_derivative(imax,r1,UAS,dus)
             
            call g1_radial_derivative(imax,r1,UBL,dudl)
            call g1_radial_derivative(imax,r1,UBS,duds)
            
            CALL RINT13(REL,a1lo(1,jlo,l),b1lo(1,jlo,l),uppl2(1,lp),upps2(1,lp),  juddu_20(jlo,l,lp),jatom)     !integral u times double radia derivative
            CALL RINT13(REL,a1lo(1,jlo,l),b1lo(1,jlo,l),udppl2(1,lp),udpps2(1,lp),juddu_21(jlo,l,lp),jatom)
            
            CALL RINT13(REL,a1lo(1,jlo,l),b1lo(1,jlo,l),dul(1),dus(1),  judu_20(jlo,l,lp),jatom)     
            CALL RINT13(REL,a1lo(1,jlo,l),b1lo(1,jlo,l),dudl(1),duds(1),judu_21(jlo,l,lp),jatom)
             
            DO I = 1, Imax                   
               UAA(I) = (a1lo(i,jlo,l)*rrad1(i,lp) + CIN*b1lo(i,jlo,l)*rrad2(i,lp))/r1(i)/r1(i)
            ENDDO           
            CALL CHARGEL2(r1,1.D0,1.D0,UAA(1),dx,imax,TC); juu_20(jlo,l,lp) = TC
                
                
            DO I = 1, Imax                   
               UAA(I) = (a1lo(i,jlo,l)*rade1(i,lp) + CIN*b1lo(i,jlo,l)*rade2(i,lp))/r1(i)/r1(i)
            ENDDO           
            CALL CHARGEL2(r1,1.D0,1.D0,UAA(1),dx,imax,TC); juu_21(jlo,l,lp) = TC                             
             
!             CALL RINT13(REL,UCL(1),UCS(1),UAL(1),UAS(1),juu_20(jlo,l,lp),jatom)
!             CALL RINT13(REL,UCL(1),UCS(1),UBL(1),UBS(1),juu_21(jlo,l,lp),jatom)
            
          ENDDO          !lp loop
       ENDDO
    ENDDO
    
    
    DO lp = 0, lmax2n
       DO jlop = 1, ilo(lp)
          IF(.NOT. loor(jlop,lp)) CYCLE
          DO i = 1, imax
             UCL(i) = a1lo(i,jlop,lp)/r1(i)                         !large comp of LO fun ul(e2)
             UCS(i) = b1lo(i,jlop,lp)/r1(i)                         !small comp of LO fun us(e2)
          ENDDO
          call g1_radial_derivative(imax,r1,UCL,dulol)         !first derivative of  ul(r,e2)
          call g1_radial_derivative(imax,r1,UCS,dulos)         !first derivative of  uS(r,e2)
                
          call g1_radial_derivative(imax,r1,dulol,ddulol)      !2nd derivative of  ul(r,e2)
          call g1_radial_derivative(imax,r1,dulos,ddulos)      !2nd derivative of  us(r,e2)  
             
          DO i = 1, imax
             ddulol(i) = r1(i)*ddulol(i)
             ddulos(i) = r1(i)*ddulos(i)
          ENDDO
          
          DO l = 0, lmax2n
             DO i = 1, imax
                UALP(i) = rrad1(i,l)/r1(i) 
                UASP(i) = rrad2(i,l)/r1(i)
                UBLP(i) = rade1(i,l)/r1(i) 
                UBSP(i) = rade2(i,l)/r1(i) 
             ENDDO
          
             CALL RINT13(REL,rrad1(1,l),rrad2(1,l),ddulol(1),ddulos(1),juddu_02(jlop,lp,l),jatom)   
             CALL RINT13(REL,rade1(1,l),rade2(1,l),ddulol(1),ddulos(1),juddu_12(jlop,lp,l),jatom) 
             CALL RINT13(REL,rrad1(1,l),rrad2(1,l),dulol(1),dulos(1),judu_02(jlop,lp,l),jatom)
             CALL RINT13(REL,rade1(1,l),rade2(1,l),dulol(1),dulos(1),judu_12(jlop,lp,l),jatom)
             
            DO I = 1, Imax                   
               UAA(I) = (a1lo(i,jlop,lp)*rrad1(i,l) + CIN*b1lo(i,jlop,lp)*rrad2(i,l))/r1(i)/r1(i)
            ENDDO           
            CALL CHARGEL2(r1,1.D0,1.D0,UAA(1),dx,imax,TC); juu_02(jlop,lp,l) = TC
                                
            DO I = 1, Imax                   
               UAA(I) = (a1lo(i,jlop,lp)*rade1(i,l) + CIN*b1lo(i,jlop,lp)*rade2(i,l))/r1(i)/r1(i)
            ENDDO           
            CALL CHARGEL2(r1,1.D0,1.D0,UAA(1),dx,imax,TC); juu_12(jlop,lp,l) = TC
                 
             
!              CALL RINT13(REL,UCL(1),UCS(1),UALP(1),UASP(1),juu_02(jlop,lp,l),jatom)
!              CALL RINT13(REL,UCL(1),UCS(1),UBLP(1),UBSP(1),juu_12(jlop,lp,l),jatom)
                                    
          ENDDO !l loop          
       ENDDO
    ENDDO
    
return
end subroutine kin_str_integrals

subroutine g1_radial_derivative(n,r,f,df)
  ! Stress-only derivative on the physical positive-r mesh.
  ! The old G1 path passed an UNINITIALIZED IsEven to DERGL, which
  ! actually uses the flag to mirror f across r=0. Even a specified
  ! parity is not generally valid for Coulomb-cusp/RELA radial functions,
  ! or for both a radial function and its derivative. No extrapolated
  ! negative-radius samples are needed: SPLINE already supplies a
  ! one-sided endpoint condition from positive-mesh divided differences.
  ! Keep the existing cubic interpolation and twice-differentiated
  ! strategy; this repair does NOT add physical radial response to strain.
  ! Generic DERGL and non-stress density/force routines are untouched.
  implicit none
  integer,intent(in) :: n
  real*8,intent(in) :: r(n),f(n)
  real*8,intent(out) :: df(n)
  real*8 :: c(n),d(n)
  if(n.lt.4) stop 'G1 radial derivative requires at least four points'
  call spline(n,r,f,df,c,d)
end subroutine g1_radial_derivative
