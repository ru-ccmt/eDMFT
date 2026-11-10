subroutine stress_gga(lm,lmmax,jspin,gga_str_int,sigma_str,lmmaxx)
    use densit  ,   only: clmsp, clm    !clmsp(j,lm1,jatom,ispin)
    use parallel,   only: nstart, nstop, myid
    use struct  ,   only: mult,jri,r0,dx,rmt,iatnr,nat
    use gga_str_kb, only: gga_fin_str   ! , sigma_str
    use work      , only: exc,vxc
    implicit none
   include 'param.inc'
#ifdef Parallel
      include 'mpif.h'
#endif
    
    integer, intent(in)    :: lm(2,ncom+3,nat), lmmax(nat), jspin,lmmaxx    
    real*8                 ::sigma_str(1:nrad,1:lmmaxx,nstart(myid):nstop(myid),1:3)

    logical                :: HEXA, CUBIC, even, spin1
    complex*16             :: cabt(3,-1:1), clm2(1:nrad,ncom+3,nat,1:2), clm3(1:nrad,ncom+3), &                            
                              dclm2(1:nrad,ncom+3,nat,1:2), &
                              fc(ncom+3,nat), value_c(1:nrad), gga_str1(nat,1:2),gga_str2(1:9,nat,1:2), &
                              d_exc_darho(1:3,1:nrad,ncom+3),d_exc(1:nrad,ncom+3), &
                              INT1(1:3,ncom+3,ncom+3), gga_str3(1:9,nat,1:2), &
                              gga_str_2nd(1:9), &
                              sigma_buf(1:nrad,ncom+3,nat,1:3),gga_str_int(1:9)    
                              
    complex*16             :: exc_buf(1:nrad,ncom+3,nat,1:2), vxc_buf(1:nrad,ncom+3,nat,1:2)                              
    real*8                 :: exc_stress_LDA
    
    complex*16             :: zeroc,one_cmplx, out_c, exh_int                          
    integer                :: jatom, mu, imax, lm1, ispin, L0, M0, lmmax_2, lmtot1,                                     &
                              lmtot(2,ncom+3), alpha, beta, index1, lm1p,lm0,l1,m1,index_a,                             &
                              l0p,m0p, s1, t,tp,t0,t0p, v1, pone, J, lm1pp,lm1ppp, v2, s2,l0pp,l0ppp,m0pp,m0ppp,v,s,ierr
                              
    real*8                 :: R(1:nrad),R3(1:nrad), sqrt_4pii, pi,out1,out2, value(1:nrad), value2(1:nrad), g0,          &
                              g1(1:nrad), g2(1:nrad), gaunt1, gaunt2, drho_a(1:nrad), gaunt3, gaunt4,buf_str4(1:9),      &
                              gga_str4(1:9,nat,1:2), dx1, on_thrd,for_thrd,fac_22, buf_spin_ud(1:9),gga_fin_str1(1:9)
    integer, external      :: NOTRI
    real*8 , external      :: GAUNT
 
!define all constant  
    one_cmplx = (1.d0,0.d0)
    zeroc     = (0.d0,0.d0)
    pi        = 4.0d+0*ATAN(1.d0)
    sqrt_4pii = 1.d0/sqrt(4.d0*pi)
    on_thrd   = 1.0/3.0
    for_thrd  = 4.D0*on_thrd
    fac_22    = 1.d0/(3.0*pi*pi)**on_thrd
    
!exc_m_vxc calculated in xcpot3, iloop 8
!initialized variables
    gga_str1    = zeroc
    gga_str2    = zeroc
    gga_str3    = zeroc
    gga_str_2nd = zeroc
    gga_fin_str = zeroc        
    buf_str4    = 0.d0
    gga_str4    = 0.d0
    buf_spin_ud = 0.d0
    exc_stress_LDA = 0.d0
    
!cabt is the c-coefficient that requires when a unit vector is expanded in term of the unit spherical harmonics.     
    call c_alpha_m(cabt)                          
    
    IF(jspin .eq. 2) then
       spin1 = .TRUE.
    ELSE
       spin1 = .FALSE.
    ENDIF
    
    ! ROTIJ is indexed over all atoms, not inequivalent atom types.  Account
    ! for atom types assigned to earlier MPI ranks before rotating this
    ! rank's GGA stress tensors.  Otherwise parallel LAPW0 uses rotations
    ! belonging to the wrong equivalent atoms whenever nstart(myid)>1.
    index_a = 0
    DO jatom = 1, nstart(myid)-1
       index_a = index_a + mult(jatom)
    ENDDO
    DO 29 jatom = nstart(myid), nstop(myid)
    DO 29 mu    = 1, mult(jatom)
          index_a = index_a + 1
          imax  = jri(jatom)
          dx1   = dx(jatom)
          HEXA  = .FALSE.
          CUBIC = .FALSE.
          IF(iatnr(jatom) .lt. 0) HEXA  = .TRUE.
          IF(iatnr(jatom) .gt. 0) CUBIC = .TRUE.
          
          DO J = 1, JRI(jatom)
             R(J) = R0(jatom)*exp( dble(j-1)*dx(jatom) )
             R3(j) = R(J)*R(J)*R(J)
          ENDDO
    
     SPIN: DO ispin = 1,2 !jspin              
              DO lm1 = 1,lmmax(jatom)
                 clm2(1:imax,lm1,jatom,ispin) = clmsp(1:imax,lm1,jatom,ispin)   
                 sigma_buf(1:imax,lm1,jatom,ispin) = sigma_str(1:imax,lm1,jatom,ispin)
                 
!                  exc_buf(1:imax,lm1,jatom,ispin)   = exc(1:imax,lm1,jatom,ispin)
!                  vxc_buf(1:imax,lm1,jatom,ispin)   = vxc(1:imax,lm1,jatom,ispin)
                 
                 IF( spin1 ) sigma_buf(1:imax,lm1,jatom,3) = sigma_str(1:imax,lm1,jatom,3)
                 IF(lm1 .eq. 1) clm2(:,lm1,jatom,ispin) = clm2(:,lm1,jatom,ispin)*sqrt_4pii                     
              ENDDO                             
!        write(6,*)'ispin',ispin,sigma_str(imax-1:imax,1,jatom,ispin),clmsp(781,1,jatom,ispin)                    
!        IF( spin1 ) write(6,*)sigma_str(imax-1:imax,1,jatom,3)
              IF(iatnr(jatom) .lt. 0) THEN   !NON-CUBIC CASE
                 lmmax_2 = lmmax(jatom)
                 DO lm1 = 1, lmmax(jatom)
                    M0 = lm(2,lm1,jatom)
                    IF(M0 .NE. 0) THEN
                      lmmax_2 = lmmax_2 + 1
                      IF(lmmax_2 .gt. ncom+3) STOP 'error in "stress_gga.f" for lm combination of NON-CUBIC CASE '
                      clm2(:,lmmax_2,jatom,ispin) = clm2(:,lm1,jatom,ispin)                      
                      sigma_buf(1:imax,lmmax_2,jatom,ispin) = sigma_buf(1:imax,lm1,jatom,ispin)
                      
!                       exc_buf(1:imax,lmmax_2,jatom,ispin) = exc_buf(1:imax,lm1,jatom,ispin)
!                       vxc_buf(1:imax,lmmax_2,jatom,ispin) = vxc_buf(1:imax,lm1,jatom,ispin)
                      
                      IF( spin1 ) sigma_buf(1:imax,lmmax_2,jatom,3) = sigma_buf(1:imax,lm1,jatom,3)
                    ENDIF
                 ENDDO
!               ENDIF
              
                 call multfc(fc,jatom,lmmax_2,lm(1,1,jatom),iatnr(jatom))
              
                 DO lm1 = 1, lmmax_2
                    clm2(:,lm1,jatom,ispin) = clm2(:,lm1,jatom,ispin)*fc(lm1,jatom)                    
                    sigma_buf(1:imax,lm1,jatom,ispin)  = sigma_buf(1:imax,lm1,jatom,ispin)*fc(lm1,jatom)
                    
!                     exc_buf(1:imax,lmmax_2,jatom,ispin) = exc_buf(1:imax,lm1,jatom,ispin)*fc(lm1,jatom)
!                     vxc_buf(1:imax,lmmax_2,jatom,ispin) = vxc_buf(1:imax,lm1,jatom,ispin)*fc(lm1,jatom)
                    
                    IF( spin1 ) sigma_buf(1:imax,lm1,jatom,3)  = sigma_buf(1:imax,lm1,jatom,3)*fc(lm1,jatom)                                                           
                 ENDDO
                 
                 lmtot1 = 0
                 DO 39 lm1 = 1, lmmax_2
                       lmtot1 = lmtot1 + 1
                       lmtot(1,lmtot1) = lm(1,lm1,jatom)
                       lmtot(2,lmtot1) = lm(2,lm1,jatom)
                       
                       clm2(:,lmtot1,jatom,ispin) = clm2(:,lm1,jatom,ispin)                       
                       sigma_buf(1:nrad,lmtot1,jatom,ispin)  = sigma_buf(1:nrad,lm1,jatom,ispin)
                       
                       exc_buf(1:nrad,lmtot1,jatom,ispin)  = exc_buf(1:nrad,lm1,jatom,ispin)
                       vxc_buf(1:nrad,lmtot1,jatom,ispin)  = vxc_buf(1:nrad,lm1,jatom,ispin)
                       
                       IF(spin1) sigma_buf(1:nrad,lmtot1,jatom,3)  = sigma_buf(1:nrad,lm1,jatom,3)
                       IF(lmtot1 .eq. 1) GOTO 39
                       IF( (IABS( lmtot(1,lmtot1) ) .EQ. IABS( lmtot(1,lmtot1-1) ) ) .AND. &
                           (  lmtot(2,lmtot1) .EQ. lmtot(2,lmtot1-1) )     ) THEN
                           clm2(:,lmtot1-1,jatom,ispin) = clm2(:,lmtot1-1,jatom,ispin)+clm2(:,lmtot1,jatom,ispin)                                                                                     
                           
                           sigma_buf(:,lmtot1-1,jatom,ispin) = sigma_buf(:,lmtot1-1,jatom,ispin) + &
                                                               sigma_buf(:,lmtot1,jatom,ispin)
                                                               
!                            exc_buf(:,lmtot1-1,jatom,ispin) = exc_buf(:,lmtot1-1,jatom,ispin) + &
!                                                                exc_buf(:,lmtot1,jatom,ispin)    
!                            vxc_buf(:,lmtot1-1,jatom,ispin) = vxc_buf(:,lmtot1-1,jatom,ispin) + &
!                                                                vxc_buf(:,lmtot1,jatom,ispin)                                        
                                                               
                           IF(spin1) sigma_buf(:,lmtot1-1,jatom,3) = sigma_buf(:,lmtot1-1,jatom,3) + &
                                                               sigma_buf(:,lmtot1,jatom,3)                                          
                           lmtot1 = lmtot1 - 1                                
                       ENDIF    
                        
39   continue                 
                 
              ELSE   !CUBIC CASE 
                   call lm_combine(imax,lmmax(jatom),clm2(1,1,jatom,ispin),lm(1,1,jatom),jatom)  
                   call lm_combine(imax,lmmax(jatom),sigma_buf(1,1,jatom,ispin),lm(1,1,jatom),jatom)
                   
!                    call lm_combine(imax,lmmax(jatom),exc_buf(1,1,jatom,ispin),lm(1,1,jatom),jatom)
!                    call lm_combine(imax,lmmax(jatom),vxc_buf(1,1,jatom,ispin),lm(1,1,jatom),jatom)
                   
                   IF(spin1) call lm_combine(imax,lmmax(jatom),sigma_buf(1,1,jatom,3),lm(1,1,jatom),jatom)
                  
                  lmtot1 = lmmax(jatom)
                  DO lm1 = 1, lmmax(jatom)
                     lmtot(1,lm1) = lm(1,lm1,jatom)
                     lmtot(2,lm1) = lm(2,lm1,jatom)
                     M0 = lm(2,lm1,jatom)
                     IF( M0 .ne. 0) THEN
                        lmtot1 = lmtot1 + 1
                        IF(lmtot1 .gt. ncom+3) STOP 'error in "stress_gga.f" for lm combination for CUBIC CASE '
                        lmtot(1,lmtot1) = lm(1,lm1,jatom)
                        lmtot(2,lmtot1) =-lm(2,lm1,jatom)
                        
                        clm2(:,lmtot1,jatom,ispin) = clm2(:,lm1,jatom,ispin)  
                        sigma_buf(:,lmtot1,jatom,ispin) = sigma_buf(:,lm1,jatom,ispin)
                        
!                         exc_buf(:,lmtot1,jatom,ispin) = exc_buf(:,lm1,jatom,ispin)
!                         vxc_buf(:,lmtot1,jatom,ispin) = vxc_buf(:,lm1,jatom,ispin)
                        
                        IF(spin1) sigma_buf(:,lmtot1,jatom,3) = sigma_buf(:,lm1,jatom,3)                                                                                           
                     ENDIF
                  ENDDO
              ENDIF                                              
              
             
!================================from gradient correction============================
    DO lm1 = 1, lmtot1                                                
       l0 = lmtot(1,lm1)
       m0 = lmtot(2,lm1)                                  
       pone = (-1)**m0
       even = .FALSE.
       IF(2*(l0/2) .eq. l0) even = .TRUE.
       DO J = 1, imax
          VALUE(j) = REAL(clm2(j,lm1,jatom,ispin) )/( R(j)*R(j) )
       ENDDO   
       call dergl(imax,r,VALUE,g1,even,g0)       
       DO J = 1, imax
          VALUE(j) = AIMAG(clm2(j,lm1,jatom,ispin) )/( R(j)*R(j) )
       ENDDO   
       call dergl(imax,r,VALUE,g2,even,g0)                                
       DO J = 1, imax
          out1 = g1(j)
          out2 = g2(j)
          dclm2(j,lm1,jatom,ispin) = dcmplx(out1,out2)          
       ENDDO   
    ENDDO
    
!==========This block calculate the same thing which is in the 'exchange_str.f'=====
!----------To use this block uncomment above exc_buf and vxc_buf, and it is just for dbugging purpose.
!     DO lm1 = 1, lmtot1                                                
!        l0 = lmtot(1,lm1)
!        m0 = lmtot(2,lm1)                                  
!        DO J = 1, imax
!           VALUE(j) =  clm2(j,lm1,jatom,ispin)*(exc_buf(j,lm1,jatom,ispin) - vxc_buf(j,lm1,jatom,ispin) )
!        ENDDO   
!        call charge(r,dx1,VALUE,1,imax,out1)
!        exc_stress_LDA = exc_stress_LDA + out1
!     ENDDO    
!_------____ end similar term as in LDA    
   

        buf_str4 = 0.0d0
        DO lm0 = 1,lmtot1
           l0 = lmtot(1,lm0)
           m0 = lmtot(2,lm0)
           DO t0 = -1,1
              DO t0p = -1,1
                 IF(NOTRI(l0,1,1) .lt. 0) CYCLE
                 IF((m0+t0+t0p) .ne. 0) CYCLE
                 pone = (-1)**m0
                 gaunt1 = GAUNT(l0,1,1,-m0,t0,t0p)
                 DO j = 1,imax
                    VALUE(j) = R(j)*R(j)*REAL(sigma_buf(j,lm0,jatom,ispin))* &
                            REAL(dclm2(j,lm0,jatom,ispin))*REAL(dclm2(j,lm0,jatom,ispin))                       
                 ENDDO
                 call charge(r,dx1,VALUE,1,imax,out1)                 
                 index1 = 0
                 DO alpha = 1,3
                    DO beta = 1,3
                       index1 = index1 + 1
                       buf_str4(index1) = buf_str4(index1) + &
                       pone*cabt(alpha,t0)*cabt(beta,t0p)*gaunt1*out1*sqrt_4pii*sqrt_4pii
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO

       !======================symmetrization using rotloc and rotij=======
        call symmetry_stress_rotij(index_a,jatom,buf_str4) 
        gga_str4(1:9,jatom,ispin) = buf_str4(1:9)                            

                          
          DO index1 = 1, 9
             gga_fin_str(index1) = gga_fin_str(index1)+ gga_str4(index1,jatom,ispin)*2                       
          ENDDO                                                 
          
     ENDDO SPIN           
     
     IF( SPIN1 ) THEN !This block is for mix spin state
       buf_str4 = 0.0
       DO lm0 = 1, lmtot1
          l0 = lmtot(1,lm0)
          m0 = lmtot(2,lm0)
          DO t0 = -1,1
             DO t0p = -1,1
                IF(NOTRI(l0,1,1) .lt. 0) CYCLE
                IF((m0+t0+t0p) .ne. 0) CYCLE
                gaunt1 = GAUNT(l0,1,1,-m0,t0,t0p)
                DO j = 1,imax
                    VALUE(j) = R(j)*R(j)*REAL(sigma_buf(j,lm0,jatom,3))* &
                               2.d0*REAL(dclm2(j,lm0,jatom,1))*REAL(dclm2(j,lm0,jatom,2))
                ENDDO
                out1=0.d0
                call charge(r,dx1,VALUE,1,imax,out1)
                index1 = 0
                DO alpha = 1, 3
                   DO beta = 1, 3
                      index1 = index1 + 1
                      buf_str4(index1) = buf_str4(index1) + &
                      pone*cabt(alpha,t0)*cabt(beta,t0p)*gaunt1*out1*sqrt_4pii*sqrt_4pii
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
!       write(*,*) 'buf_str4',lm0,out1,buf_str4(1),sigma_buf(imax-3:imax,lm0,jatom,3)
       ENDDO
       call symmetry_stress_rotij(index_a,jatom,buf_str4)
       buf_spin_ud = buf_spin_ud + buf_str4                                   
     ENDIF
       
29  continue    

#ifdef Parallel
   gga_fin_str1=gga_fin_str
!
   call mpi_reduce(gga_fin_str1,gga_fin_str,9,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
#endif

     write(6,*) spin1,' gga-total:',2.d0*gga_fin_str(1)+gga_str_int(1)+ 2.d0*buf_spin_ud(1),' gga-int:',gga_str_int(1)
     write(6,*) '  gga-sphere:',2.d0*gga_fin_str(1),' gga_ud:',buf_spin_ud(1)*2.d0
     gga_fin_str = 2.d0*gga_fin_str+gga_str_int  
     IF( SPIN1 ) gga_fin_str = gga_fin_str + buf_spin_ud *2.0d0   ! to make spin-pol = NON-SPINPOL! 
     
return
END
