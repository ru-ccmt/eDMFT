! Atom-dependent radial integrals for the analytic LAPW stress.
!
! Adapted from WIEN2k 24 SRC_lapw2.  The helpers are renamed so they do
! not collide with the unrelated DERGL/RINT13 implementation in libatpar.

subroutine val_corr_integrals(REL,rrad1,rrad2,rade1,rade2,a1lo,b1lo,jatom,imax,dx)   
   use param
   use structure,     only: r0
   use lo,        only: ilo,loor
   USE stress_kb, only: jllp00,jllp10,jllp01,jllp11, &
                        Illp00,Illp01,Illp10,Illp11, & 
                        jr3upu_20,jr3upu_21,jr3upu_22 , &
                        jr3upu_02,jr3upu_12,jlaplo1_12, &
                        jlaplo1_02
                                         
   IMPLICIT NONE
   logical             :: REL
   integer,intent(in)  :: jatom,imax
   real*8 ,intent(in)  :: rrad1(nrad,0:lmax2),rrad2(nrad,0:lmax2), &     !before increase the radial points
                          rade1(nrad,0:lmax2),rade2(nrad,0:lmax2), &
                          a1lo(nrad,nloat,0:lomax),b1lo(nrad,nloat,0:lomax), &
                          dx

   integer             :: i,l,jlo,lp,jlop

   real*8              :: dul(nrad),dus(nrad),dudl(nrad),duds(nrad), r1(nrad), &
                          dulol(nrad),dulos(nrad), ddulol(nrad),ddulos(nrad)
                          
   real*8              :: UAL(nrad),UAS(nrad),UBL(nrad),UBS(nrad), &    !old
                          UCL(nrad),UCS(nrad),UAB(nrad),UBA(nrad), &
                          UAA(nrad)                                                     
                          
   real*8              :: CIN,CINP,g0,TC                          
                          
                          
    CIN = 1.d0/137.0359895d0**2
    CINP = 1.d0/137.0359895d0
    if(.not.REL) CIN = 1E-22                       
    
    DO i = 1, imax
       r1(i) = r0(jatom)*exp(dx*(i-1))
    ENDDO    
                          
    DO l = 0, lmax2n
       DO lp = 0, lmax2n
          
          DO I = 1, Imax
             UAL(i) = rrad1(i,l)/r1(i) 
             UAS(i) = rrad2(i,l)/r1(i)
             UBL(i) = rade1(i,l)/r1(i) 
             UBS(i) = rade2(i,l)/r1(i) 
          ENDDO        
          call dergl_stress(imax,r1,UAL,dul,.true.,g0)
          call dergl_stress(imax,r1,UAS,dus,.true.,g0)
             
          call dergl_stress(imax,r1,UBL,dudl,.true.,g0)
          call dergl_stress(imax,r1,UBS,duds,.true.,g0)
          
          DO I = 1, Imax   
             dul(i) = r1(i)*r1(i)*dul(i)
             dus(i) = r1(i)*r1(i)*dus(i)
             
             dudl(i) = r1(i)*r1(i)*dudl(i)
             duds(i) = r1(i)*r1(i)*duds(i)
          ENDDO   
          
          CALL rint13_stress(rel, dul(1), dus(1),rrad1(1,lp),rrad2(1,lp),  Jllp00(l,lp), jatom)
          CALL rint13_stress(rel,dudl(1),duds(1),rrad1(1,lp),rrad2(1,lp),  Jllp10(l,lp), jatom)
          CALL rint13_stress(rel, dul(1), dus(1),rade1(1,lp),rade2(1,lp),  Jllp01(l,lp), jatom)
          CALL rint13_stress(rel,dudl(1),duds(1),rade1(1,lp),rade2(1,lp),  Jllp11(l,lp), jatom)
             
          CALL rint13_stress(rel,rrad1(1,l),rrad2(1,l),rrad1(1,lp),rrad2(1,lp),  Illp00(l,lp), jatom)
          CALL rint13_stress(rel,rrad1(1,l),rrad2(1,l),rade1(1,lp),rade2(1,lp),  Illp01(l,lp), jatom)
          CALL rint13_stress(rel,rade1(1,l),rade2(1,l),rrad1(1,lp),rrad2(1,lp),  Illp10(l,lp), jatom)
          CALL rint13_stress(rel,rade1(1,l),rade2(1,l),rade1(1,lp),rade2(1,lp),  Illp11(l,lp), jatom)                    
          
       ENDDO
    ENDDO
    
    
    DO l = 0, lmax2n
       DO jlo = 1, ilo(l)
          IF(.NOT. loor(jlo,l)) CYCLE  
          
          call rint13_stress(rel,rade1(1,l),rade2(1,l),a1lo(1,jlo,l),b1lo(1,jlo,l),TC,jatom)
            jlaplo1_12(jlo,l)  = TC
          call rint13_stress(rel,rrad1(1,l),rrad2(1,l),a1lo(1,jlo,l),b1lo(1,jlo,l),TC,jatom)  
            jlaplo1_02(jlo,l)  = TC
            
          DO lp = 0, lmax2n          
             DO I = 1, Imax
                UAL(i) = rrad1(i,lp)/r1(i) 
                UAS(i) = rrad2(i,lp)/r1(i)
                UBL(i) = rade1(i,lp)/r1(i) 
                UBS(i) = rade2(i,lp)/r1(i) 
             ENDDO        
             call dergl_stress(imax,r1,UAL,dul,.true.,g0)
             call dergl_stress(imax,r1,UAS,dus,.true.,g0)
             
             call dergl_stress(imax,r1,UBL,dudl,.true.,g0)
             call dergl_stress(imax,r1,UBS,duds,.true.,g0)
          
             DO I = 1, Imax   
                dul(i) = r1(i)*r1(i)*dul(i)
                dus(i) = r1(i)*r1(i)*dus(i)
             
                dudl(i) = r1(i)*r1(i)*dudl(i)
                duds(i) = r1(i)*r1(i)*duds(i)
             ENDDO   
             CALL rint13_stress(REL,a1lo(1,jlo,l),b1lo(1,jlo,l),dul(1),dus(1),TC,jatom)     
                jr3upu_20(jlo,l,lp) = TC
             CALL rint13_stress(REL,a1lo(1,jlo,l),b1lo(1,jlo,l),dudl(1),duds(1),TC,jatom) 
                jr3upu_21(jlo,l,lp) = TC
          ENDDO  !lp loop
          
       ENDDO
    ENDDO
    
    
    DO lp = 0, lmax2n
       DO jlop = 1, ilo(lp)
          IF(.NOT. loor(jlop,lp)) CYCLE
          DO l = 0, lmax2n
             
             DO i = 1, imax
                UCL(i) = a1lo(i,jlop,lp)/r1(i)                         !large comp of LO fun ul(e2)
                UCS(i) = b1lo(i,jlop,lp)/r1(i)                         !small comp of LO fun us(e2)                  
             ENDDO
             call dergl_stress(imax,r1,UCL,dulol,.true.,g0)         !first derivative of  ul(r,e2)
             call dergl_stress(imax,r1,UCS,dulos,.true.,g0)         !first derivative of  uS(r,e2)
             
             DO i = 1, imax
                dulol(i) = dulol(i)*r1(i)*r1(i)
                dulos(i) = dulos(i)*r1(i)*r1(i)
             ENDDO
             
             CALL rint13_stress(REL,rrad1(1,l),rrad2(1,l),dulol(1),dulos(1),TC,jatom)   
               jr3upu_02(jlop,lp,l) = TC
             CALL rint13_stress(REL,rade1(1,l),rade2(1,l),dulol(1),dulos(1),TC,jatom)
               jr3upu_12(jlop,lp,l) = TC 
             
             DO jlo = 1, ilo(l)
                IF(.NOT. loor(jlo,l)) CYCLE
                  CALL rint13_stress(REL,a1lo(1,jlo,l),b1lo(1,jlo,l),dulol(1),dulos(1),TC,jatom)
                  jr3upu_22(jlo,jlop,l,lp) = TC
             ENDDO
             
          ENDDO
       ENDDO
    ENDDO

return
end subroutine val_corr_integrals

subroutine kin_str_integrals(REL, rrad1,rrad2,rade1,rade2,a1lo,b1lo,jatom,imax,dx,VR)
!! The following integrals are used to program the kinetic stress tensor, equation 6.71 at the moment. The integrals dependent only in the atomic index and corresponding, therefore
! it is execute inside the atomic index loop but outside of k-points loop. This subroutine is being called from atpar.F (at the moment).
!TC1S00  = integral with double derivative           r2*(u8ddu)  for LOs juddu
!TC2S00  = integral with single derivative           r1*(u*du)   for LOs judu
!TC3S00  = integral cont product of two functions.   (u*u)       for LOs juu
   use param
   use structure,     only: r0
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
             
       ! u_l(r)/r behaves as r**l at the origin.  DERGL_STRESS extends
       ! the radial mesh to negative r, so its parity flag must follow l.
       ! The distributed WIEN2k source leaves this assignment commented,
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
          
          CALL rint13_stress(rel,rrad1(1,l),rrad2(1,l),uppl2(1,lp),upps2(1,lp),  TC1S00(l,lp),jatom) !integrand contains just product of two functions: r*r*u_llmd*ddu_llmdp; dd double derivative
          CALL rint13_stress(rel,rrad1(1,l),rrad2(1,l),udppl2(1,lp),udpps2(1,lp),TC1S01(l,lp),jatom)
          CALL rint13_stress(rel,rade1(1,l),rade2(1,l),uppl2(1,lp),upps2(1,lp),  TC1S10(l,lp),jatom)
          CALL rint13_stress(rel,rade1(1,l),rade2(1,l),udppl2(1,lp),udpps2(1,lp),TC1S11(l,lp),jatom)                                                 
          
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
          
          CALL rint13_stress(rel,rrad1(1,l),rrad2(1,l),dul(1),dus(1),  TC2S00(l,lp), jatom)   !integrand contains just product of two functions: r*u_llmd*du_llmdp; d derivative
          CALL rint13_stress(rel,rrad1(1,l),rrad2(1,l),dudl(1),duds(1),TC2S01(l,lp), jatom)
          CALL rint13_stress(rel,rade1(1,l),rade2(1,l),dul(1),dus(1),  TC2S10(l,lp), jatom)
          CALL rint13_stress(rel,rade1(1,l),rade2(1,l),dudl(1),duds(1),TC2S11(l,lp), jatom)                     
           
          CALL rint13_stress(rel,UALP(1),UASP(1),UAL(1),UAS(1),  TC3S00(l,lp), jatom)   !integrand contains just product of two functions: u_llmd*u_llmdp
          CALL rint13_stress(rel,UALP(1),UASP(1),UBL(1),UBS(1),  TC3S01(l,lp), jatom)
          CALL rint13_stress(rel,UBLP(1),UBSP(1),UAL(1),UAS(1),  TC3S10(l,lp), jatom)
          CALL rint13_stress(rel,UBLP(1),UBSP(1),UBL(1),UBS(1),  TC3S11(l,lp), jatom)                    
          
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
                
                CALL rint13_stress(REL,a1lo(1,jlo,l),b1lo(1,jlo,l),ddulol(1),ddulos(1),juddu_22(jlo,jlop,l,lp),jatom)
                CALL rint13_stress(REL,a1lo(1,jlo,l),b1lo(1,jlo,l),dulol(1),dulos(1),  judu_22(jlo,jlop,l,lp),jatom)    
                
                DO I = 1, Imax                   
                   UAA(I) = (a1lo(i,jlop,lp)*a1lo(i,jlo,l) + CIN*b1lo(i,jlop,lp)*b1lo(i,jlo,l))/r1(i)/r1(i)
                ENDDO           
                CALL CHARGE(r1,1.D0,1.D0,UAA(1),dx,imax,TC);juu_22(jlo,jlop,l,lp) = TC
                                   
!                 CALL rint13_stress(REL,UCL(1),UCS(1),UCLP(1),UCSP(1),juu_22(jlo,jlop,l,lp),jatom)
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
            
            CALL rint13_stress(REL,a1lo(1,jlo,l),b1lo(1,jlo,l),uppl2(1,lp),upps2(1,lp),  juddu_20(jlo,l,lp),jatom)     !integral u times double radia derivative
            CALL rint13_stress(REL,a1lo(1,jlo,l),b1lo(1,jlo,l),udppl2(1,lp),udpps2(1,lp),juddu_21(jlo,l,lp),jatom)
            
            CALL rint13_stress(REL,a1lo(1,jlo,l),b1lo(1,jlo,l),dul(1),dus(1),  judu_20(jlo,l,lp),jatom)     
            CALL rint13_stress(REL,a1lo(1,jlo,l),b1lo(1,jlo,l),dudl(1),duds(1),judu_21(jlo,l,lp),jatom)
             
            DO I = 1, Imax                   
               UAA(I) = (a1lo(i,jlo,l)*rrad1(i,lp) + CIN*b1lo(i,jlo,l)*rrad2(i,lp))/r1(i)/r1(i)
            ENDDO           
            CALL CHARGE(r1,1.D0,1.D0,UAA(1),dx,imax,TC); juu_20(jlo,l,lp) = TC
                
                
            DO I = 1, Imax                   
               UAA(I) = (a1lo(i,jlo,l)*rade1(i,lp) + CIN*b1lo(i,jlo,l)*rade2(i,lp))/r1(i)/r1(i)
            ENDDO           
            CALL CHARGE(r1,1.D0,1.D0,UAA(1),dx,imax,TC); juu_21(jlo,l,lp) = TC                             
             
!             CALL rint13_stress(REL,UCL(1),UCS(1),UAL(1),UAS(1),juu_20(jlo,l,lp),jatom)
!             CALL rint13_stress(REL,UCL(1),UCS(1),UBL(1),UBS(1),juu_21(jlo,l,lp),jatom)
            
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
          
             CALL rint13_stress(REL,rrad1(1,l),rrad2(1,l),ddulol(1),ddulos(1),juddu_02(jlop,lp,l),jatom)   
             CALL rint13_stress(REL,rade1(1,l),rade2(1,l),ddulol(1),ddulos(1),juddu_12(jlop,lp,l),jatom) 
             CALL rint13_stress(REL,rrad1(1,l),rrad2(1,l),dulol(1),dulos(1),judu_02(jlop,lp,l),jatom)
             CALL rint13_stress(REL,rade1(1,l),rade2(1,l),dulol(1),dulos(1),judu_12(jlop,lp,l),jatom)
             
            DO I = 1, Imax                   
               UAA(I) = (a1lo(i,jlop,lp)*rrad1(i,l) + CIN*b1lo(i,jlop,lp)*rrad2(i,l))/r1(i)/r1(i)
            ENDDO           
            CALL CHARGE(r1,1.D0,1.D0,UAA(1),dx,imax,TC); juu_02(jlop,lp,l) = TC
                                
            DO I = 1, Imax                   
               UAA(I) = (a1lo(i,jlop,lp)*rade1(i,l) + CIN*b1lo(i,jlop,lp)*rade2(i,l))/r1(i)/r1(i)
            ENDDO           
            CALL CHARGE(r1,1.D0,1.D0,UAA(1),dx,imax,TC); juu_12(jlop,lp,l) = TC
                 
             
!              CALL rint13_stress(REL,UCL(1),UCS(1),UALP(1),UASP(1),juu_02(jlop,lp,l),jatom)
!              CALL rint13_stress(REL,UCL(1),UCS(1),UBLP(1),UBSP(1),juu_12(jlop,lp,l),jatom)
                                    
          ENDDO !l loop          
       ENDDO
    ENDDO
    
return
end subroutine kin_str_integrals

!     Calculate the first derivate of f
!
      subroutine dergl_stress(n,r,f,der1,even,g0)
!
!     Spline derivative, can be unstable at end points
!
      implicit real*8 (a-h,o-z)
!
      real*8 r(N), f(N), der1(N)
      real*8,allocatable  :: W1(:),W2(:),W3(:),W4(:),rw(:)
!
!     For origin corrections
      logical even, correct
!
       allocate (w1(-N+1:N),W2(-N+1:N),rw(-N+1:N),w3(-N+1:N),w4(-N+1:N))
       N2=2*N
       N0=-N+1
       rw(1:n)=r(1:n)
       w1(1:n)=f(1:n)
       if(.not.even)then
                do J=1,n
                        rw(-J+1)=-r(j)
                        w1(-J+1)=-f(J)
                enddo
        else
                do J=1,n
                        rw(-J+1)=-r(j)
                        w1(-J+1)= f(J)
                enddo
        endif
        call spline(n2,rw,w1,w2,w3,w4)
        DER1(1:N)=w2(1:N)

        G0=Der1(1)
!
      deallocate (w1,W2,rw,w3,w4)
      return
      end


subroutine rint13_stress(rel,a,b,x,y,s,jatom)
  ! Radial overlap used only by the stress precomputation.  This is the
  ! WIEN2k RINT13 integrand evaluated with dmft2's CHARGE calling sequence.
  use param,     only: nrad
  use structure, only: r0, dx, jri
  implicit none
  logical, intent(in) :: rel
  integer, intent(in) :: jatom
  real*8, intent(in)  :: a(nrad), b(nrad), x(nrad), y(nrad)
  real*8, intent(out) :: s
  real*8 :: rx(nrad), funct(nrad), cin
  integer :: j

  cin=1.d0/137.0359895d0**2
  if (.not.rel) cin=1.d-22
  do j=1,jri(jatom)
     rx(j)=r0(jatom)*exp(dx(jatom)*dble(j-1))
     funct(j)=a(j)*x(j)+cin*b(j)*y(j)
  enddo
  call charge(rx,1.d0,1.d0,funct,dx(jatom),jri(jatom),s)
end subroutine rint13_stress

subroutine rotate_stress(vector,rotmat,rotvec)
  implicit none
  real*8, intent(in)  :: vector(3), rotmat(3,3)
  real*8, intent(out) :: rotvec(3)
  rotvec=matmul(rotmat,vector)
end subroutine rotate_stress


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
