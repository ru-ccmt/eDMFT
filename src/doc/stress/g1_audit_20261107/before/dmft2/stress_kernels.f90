! Analytic LAPW valence-stress kernels.
!
! These routines are adapted from WIEN2k 24 SRC_lapw2.  Their angular
! algebra is intentionally kept unchanged so that the DFT limit can be
! checked term by term against lapw2.  dmft2-specific contractions with
! rho and rhoE are performed by cmp_dmft_stress, outside these kernels.

  subroutine c_alpha_m(c_matrix)
     implicit none
     complex*16         :: c_matrix(3,-1:1)
     complex*16, parameter :: imag  = (0.D0,1.d0)
     complex*16, parameter :: zeroc = (0.d0,0.d0)
     real*8, parameter     :: pi    = acos(-1.d0)
     real*8, parameter     :: two   = 2.D0
     real*8, parameter     :: three = 3.d0
     real*8, parameter     :: twopi = 2.D0*pi
    
       c_matrix(1,-1) = sqrt(twopi/three)
       c_matrix(1,0)  = zeroc
       c_matrix(1,1)  = -sqrt(twopi/three)

       c_matrix(2,-1) = imag*sqrt(twopi/three)    !our work doesn't have minus, kuppel does
       c_matrix(2,0)  = zeroc
       c_matrix(2,1)  = imag*sqrt(twopi/three)    !our work doesn't have minus, kuppel does

       c_matrix(3,-1) = zeroc
       c_matrix(3,0)  = sqrt(two*twopi/three)
       c_matrix(3,1)  = zeroc
  
  RETURN
  end subroutine c_alpha_m
  
      
  subroutine C_factor(ll,mm, c_ast_lm)
     !=====subroutine to calculate c_ast(l,m) which being used to calculate derivative of the spherical hermonics==========     
     implicit none
     integer,    intent(in)   :: ll,mm     
     complex*16, intent(out)  :: c_ast_lm(3,-1:1,-1:1)
     complex*16, parameter    :: zeroc = (0.d0,0.d0)
     complex*16, parameter    :: imag  = (0.D0,1.D0)
     integer,    parameter    :: one = 1.D0, two = 2.D0, three = 3.D0, f4r = 4.D0      
     integer                  :: index1, ss, tt, index, alpha, beta
     real*8                   :: delta(3,3), C1, C2, FF(6)

                    
     C1  = (ll-mm+one)*(ll-mm+two)
     C2  = f4r*(two*ll+one)*(two*ll+three)
     FF(1) = sqrt(C1/C2)

     C1 = (ll+mm-one)*(ll+mm)
     C2 = f4r*(two*ll-one)*(two*ll+one)
     FF(2) = sqrt(C1/C2)

     C1 = (ll+mm+one)*(ll+mm+two)
     C2 = f4r*(two*ll+one)*(two*ll+three)
     FF(3) = sqrt(C1/C2)

     C1 = (ll-mm-one)*(ll-mm)
     C2 = f4r*(two*ll-one)*(two*ll+one)
     FF(4) = sqrt(C1/C2)

     C1 = (ll+one)*(ll+one) - mm*mm
     C2 = (two*ll+one)*(two*ll+three)
     FF(5) = sqrt(C1/C2)

     C1 = ll*ll - mm*mm
     C2 = (two*ll-one)*(two*ll+one)
     FF(6) = sqrt(C1/C2)
     
                                                  
     delta  = 0
     DO alpha = 1, 3
        DO beta = 1, 3
           if(alpha .eq. beta) then
               delta(alpha,beta) = 1.d0
           else
               delta(alpha,beta) = 0.d0
           endif
        ENDDO
     ENDDO
     
     c_ast_lm(:,0,:) = zeroc
         
     DO alpha = 1, 3
           
        c_ast_lm(alpha,-1,-1) = -(ll+1)*(delta(alpha,1)+imag*delta(alpha,2))*FF(2)
        c_ast_lm(alpha,-1, 0) =  (ll+1)*delta(alpha,3)*FF(6)
        c_ast_lm(alpha,-1, 1) =  (ll+1)*(delta(alpha,1)-imag*delta(alpha,2))*FF(4)
              
        c_ast_lm(alpha, 1,-1) = -ll*(delta(alpha,1)+imag*delta(alpha,2))*FF(1)
        c_ast_lm(alpha, 1, 0) = -ll*delta(alpha,3)*FF(5)
        c_ast_lm(alpha, 1, 1) =  ll*(delta(alpha,1) - imag*delta(alpha,2))*FF(3)
                           
     ENDDO     
  return
  end subroutine C_factor



! SUBROUTINE val_corr_stress(isize,ibpp_max,ibb,nemin,nemax,ALM,BLM,clm,CAPD0LM,CAPD1LM,buf_sph_val_fin, jatom, rmt,latom, &
!                             lfirst)
SUBROUTINE val_corr_stress(isize,ibpp_max,ibb,nemin,nemax,ALM,BLM,clm, &
     buf_sph_val_fin,buf_nsp_match,buf_nsp_radial,buf_nsp_angular, &
     buf_sph_match,buf_sph_radial,buf_sph_angular,do_nsp,jatom,rmt, &
     latom,lfirst)
    use param
    use xa3   ,    only: a,bkx, bky, bkz
                         
    use structure ,    only: vol, rotloc, rotij, br1, pos, tauij
    use defs  ,    only: zeroc, zero, pi, TWO, imag
    use atspdt,    only: pei, el, P, DP, PE, DPE
    use xa    ,    only: e, fj, dfj        
    use stress_pw, only: ddfj, dfj_st,ddfj_st, KPLUSG 
    use stress_kb, only: jllp00,jllp10,jllp01,jllp11,&
                         Illp00,Illp01,Illp10,Illp11, &
                         jr3upu_20,jr3upu_21,jr3upu_22 , &
                         jr3upu_02,jr3upu_12,jlaplo1_12, &
                         jlaplo1_02
    use lo,        only: loor,elo,ilo,plo,dplo,lapw            !plo and dplo are the radial fn and derivative for lo at rmt                         
    implicit none
    REAL*8, parameter       :: hf = 1.0/2.D0
    integer, intent(in)     :: nemin,nemax,isize,ibpp_max,ibb,jatom,latom,lfirst
    logical, intent(in)     :: do_nsp
    complex*16, intent(in)  :: ALM((lmax2+1)*(lmax2+1),nume), &
                               BLM((lmax2+1)*(lmax2+1),nume)
    complex*16, intent(in)  :: clm((lmax2+1)*(lmax2+1),nume,nloat)                                   
    
    complex*16, intent(out) :: buf_sph_val_fin(1:3,1:3,nume), &
                               buf_nsp_match(1:3,1:3,nume), &
                               buf_nsp_radial(1:3,1:3,nume), &
                               buf_nsp_angular(1:3,1:3,nume), &
                               buf_sph_match(1:3,1:3,nume), &
                               buf_sph_radial(1:3,1:3,nume), &
                               buf_sph_angular(1:3,1:3,nume)
    real*8    , intent(in)  :: rmt
    
    complex*16              :: CAPD0LM(3,3,(lmax2+1)*(lmax2+1),nume), &
                               CAPD1LM(3,3,(lmax2+1)*(lmax2+1),nume)
    complex*16 before_lo                           
    integer                 :: delta(0:1,0:1), lmd, lmdp, ii
    integer                 :: sp,t,tp,u,v,v1, l, lp, m, mp , index, index1, &
                               index2, alpha, beta, num, i3, index3, index4, &
                               ig, igp, jlo, jlop, maxx
                               
    real*8                  :: gaunt1, gaunt2, jlpl(0:1,0:1,0:lmax2,0:lmax2), &
                               Ilpl(0:1,0:1,0:lmax2,0:lmax2)                               
    integer, external       :: NOTRI
    real*8,  external       :: GAUNT
    complex*16              :: c_am(3,-1:1), ca_st_lm(3,-1:1,-1:1), &
                               ca_st_lm11(3,-1:1,-1:1), &
                               buf_sph_val1(3,3,nume), buf_sph_val2(3,3,nume), &
                               buf_sph_val0(3,3,nume), buf_sph_val1_buf1(3,3,nume), &
                               buf_sph_val1_buf2(3,3,nume),buf_sph_val1_buf3(3,3,nume), &
                               ABLM(0:1,(lmax2+1)*(lmax2+1),nume)
                               
    integer                 :: maxindex, lda, ldb, ldc, lmx_tmp_st                          
    real*8                  :: twopi,fac,fac1, d_al1(isize), d_bl1(isize), &
                               h_al(isize),h_bl(isize), rmt2p
                                       
    complex*16              :: pcfac,pcfac1,pcfac2
    complex*16              :: h_yl1((lmax2+2)*(lmax2+2),isize), &                                                              
                               dl0ab_buf((lmax2+1)*(lmax2+1), isize), &
                               dl1ab_buf((lmax2+1)*(lmax2+1), isize), &
                               CAPD0LM_BUF((lmax2+1)*(lmax2+1),nume), &
                               CAPD1LM_BUF((lmax2+1)*(lmax2+1),nume)
                                   
!     complex*16, intent(in)  :: h_yl((lmax2+1)*(lmax2+1),isize) !, &
                             !  yl1((lmax2+2)*(lmax2+2)), phshel(isize), &
                             !  h_yll11((lmax2+2)*(lmax2+2),isize)
    COMPLEX*16              :: yl1((lmax2+2)*(lmax2+2)), phshel(isize), &
                               h_yll11((lmax2+2)*(lmax2+2),isize) , &
                               h_yl2(isize)
                               
    REAL*8                  :: bkrloc_p(3,isize), bk(3), bkr(3), bkr1(3)
    REAL*8                  :: tim_tmp1_D,tim_tmp2_D, time_cal_D
                
!   write(6,*)'in 2dd',ddfj_st(0,1) 
    twopi = two*pi
    fac   = 4.D0*pi*rmt**2/sqrt(vol)
    fac1  = 4.D0*pi/sqrt(vol)
    rmt2p = 1.d0/(rmt*rmt)
   
        
    buf_sph_val0 = zeroc
    buf_sph_val1 = zeroc    
    buf_sph_val2 = zeroc
    buf_nsp_match = zeroc
    buf_nsp_radial = zeroc
    buf_nsp_angular = zeroc
    
    buf_sph_val1_buf1 = zeroc
    buf_sph_val1_buf2 = zeroc
    buf_sph_val1_buf3 = zeroc
        
    
    call c_alpha_m(c_am)
    
    index = (lmax2n+1)*(lmax2n+1)
    DO L = 1, index
       DO num = nemin, nemax
          ABLM(0,L,num) = ALM(L,num)
          ABLM(1,L,num) = BLM(L,num)
       ENDDO
    ENDDO  
    
    delta = 0
    DO lmd = 0, 1
       DO lmdp = 0, 1
          IF(lmd .eq. lmdp) delta(lmd,lmdp) = 1
       ENDDO
    ENDDO
        
    DO l = 0, lmax2n
       DO lp = 0, lmax2n                    
          jlpl(0,0,l,lp) = jllp00(l,lp)          
          jlpl(0,1,l,lp) = jllp01(l,lp)          
          jlpl(1,0,l,lp) = jllp10(l,lp)          
          jlpl(1,1,l,lp) = jllp11(l,lp) 
          
          Ilpl(0,0,l,lp) = Illp00(l,lp)
          Ilpl(0,1,l,lp) = Illp01(l,lp)
          Ilpl(1,0,l,lp) = Illp10(l,lp)
          Ilpl(1,1,l,lp) = Illp11(l,lp)
       ENDDO
    ENDDO
    
!===========to calculate the phase factor, and spherical harmonics================
!       call cputim(tim_tmp1_D)
       DO i3 = 1, isize
          bk(1) = bkx(i3)
          bk(2) = bky(i3)
          bk(3) = bkz(i3)
          
          CALL ROTATE_STRESS (bk,rotij(1,1,latom), bkr)
          
          bk(1) = bkr(1)*BR1(1,1) + bkr(2)*BR1(1,2) + bkr(3)*BR1(1,3)
          bk(2) = bkr(1)*BR1(2,1) + bkr(2)*BR1(2,2) + bkr(3)*BR1(2,3)
          bk(3) = bkr(1)*BR1(3,1) + bkr(2)*BR1(3,2) + bkr(3)*BR1(3,3)
          
         CALL ROTATE_STRESS (bk,rotloc(1,1,jatom), bkr1)
         
         CALL YLM(bkr1, lmax2+1, yl1)
         h_yll11(:,i3) = yl1(:)
         
         phshel(i3) = exp(imag*( bkr(1)*pos(1,lfirst)*twopi + &
                                 bkr(2)*pos(2,lfirst)*twopi + &
                                 bkr(3)*pos(3,lfirst)*twopi + &
                               bkx(i3)*tauij(1,latom)*twopi + &
                               bky(i3)*tauij(2,latom)*twopi + &
                               bkz(i3)*tauij(3,latom)*twopi ) )
         
          bkrloc_p(1,i3) = bkr1(1)
          bkrloc_p(2,i3) = bkr1(2)
          bkrloc_p(3,i3) = bkr1(3)   
       ENDDO                                                                      !END i3
        
 !================================================================================
 
 
 !  go to 977  
 !----calculating DLM(a,b) -------------
     maxindex = (lmax2n+1)**2
 !_REAL           lda=2*(LMAX2+1)*(LMAX2+1)
 !_COMPLEX        lda=(LMAX2+1)*(LMAX2+1)             
!              ldb = ibpp_max
!              ldc = lda

              ldb = isize
              ldc = ldb
              
              CAPD0LM = zeroc
              CAPD1LM = zeroc
            
     DO alpha = 1,3
        DO beta = 1,3
              h_yl1 = zeroc                        
              h_yl2 = zeroc                         
                 index = 0
                 DO l = 0, lmax2
                    pcfac = imag**l
                    pcfac = pcfac*fac             
                    DO m = -l,l
                       index = index+1                      
                       call c_factor(l,m,ca_st_lm)
                       
                     DO i3 = 1, isize  
                       DO v = -1,1,2
                          DO t = -1,1
                             IF((l+v) .lt. 0) CYCLE
                             IF((m+t) .lt. -(l+v)) CYCLE
                             IF((m+t) .gt.  (l+v)) CYCLE
                             index1 = (l+v)*(l+v+1)+(t+m) + 1
                             
                             h_yl1(index,i3) = h_yl1(index,i3) - &
                                ( conjg(ca_st_lm(alpha,v,t))*bkrloc_p(beta,i3) + &
                                 conjg(ca_st_lm(beta,v,t) )*bkrloc_p(alpha,i3) )* &
                                 conjg(h_yll11(index1,i3))*phshel(i3)/two                                                                                                                                                                     
                          ENDDO
                       ENDDO
                                             
                       
                     ENDDO                                !i3 END
                     
                     DO i3 = 1, isize
                       
                       d_al1(i3) = rmt*DPE(l)*dfj_st(l,i3) - &
                                   PE(l)*dfj_st(l,i3)-PE(l)*ddfj_st(l,i3)
!                           if(d_al1(i3).gt.10.d0.and.index.eq.1) write(6,'("d_al1",2i4,7f20.10)')  l,i3,ddfj_st(l,i3)
                                   
                       d_bl1(i3) =-rmt*DP(l)*dfj_st(l,i3) + &
                                   P(l)*dfj_st(l,i3)+P(l)*ddfj_st(l,i3)                                                                              
                       
                       h_al(i3)  = dfj(l,i3,jatom)*pe(l)-fj(l,i3,jatom)*dpe(l) 
                       h_bl(i3)  = fj(l,i3,jatom)*dp(l)-dfj(l,i3,jatom)*p(l)                      
                       IF(KPLUSG(i3) .eq. 0.D0) THEN
                         dl0ab_buf(index,i3) = zeroc                            !This is because we have the |k+G| in the denominator
                         dl1ab_buf(index,i3) = zeroc                                                
                       ELSE                                                                    
                               if(abs(kplusg(i3)).lt.1.10-6) print*,'kplusg',i3,kplusg(i3) 
                         dl0ab_buf(index,i3) = ( h_yl1(index,i3)*h_al(i3) + &
                         bkrloc_p(alpha,i3)*bkrloc_p(beta,i3)*conjg(h_yll11(index,i3))*&
                         phshel(i3)*d_al1(i3) )*pcfac/sqrt(KPLUSG(i3))
                         
                         dl1ab_buf(index,i3) = ( h_yl1(index,i3)*h_bl(i3) + &
                         bkrloc_p(alpha,i3)*bkrloc_p(beta,i3)*conjg(h_yll11(index,i3))*&
                         phshel(i3)*d_bl1(i3) )*pcfac/sqrt(KPLUSG(i3))
                          
                       ENDIF 
!                        IF(i3.eq.1.and.index.eq.1.and.abs(dl0ab_buf(index,i3)).gt.1.d0) write(6, '("dl0",6f20.10)') dl0ab_buf(index,i3), h_yll11(index,i3),d_al1(i3)
                     ENDDO                                                                                              !END i3
                     
                    ENDDO                                                                                               !END m
                 ENDDO                                                                                                  !END l 
                      
                      !commented before
                      !!ucomment above initialization also
                      DO index = 1, maxindex
                        DO num = nemin,nemax
                           DO i3 = 1, isize
                              CAPD0LM(alpha,beta,index,num) = CAPD0LM(alpha,beta,index,num) + &
                                             a(i3,num)*dl0ab_buf(index,i3)
                                             
                              CAPD1LM(alpha,beta,index,num) = CAPD1LM(alpha,beta,index,num) + &                  
                                             a(i3,num)*dl1ab_buf(index,i3)
!    if(i3.eq.1.and.index.eq.1.and.num.eq.1.and.abs(CAPD0LM(alpha,beta,index,num)).gt.1.d0.and.abs(CAPD0LM(alpha,beta,index,num)).lt.2000.d0) then
!                                            write(6,*) 'i3',i3,a(i3,num),dl0ab_buf(index,i3),CAPD0LM(alpha,beta,index,num)
!                                    endif
                           ENDDO                  
                        ENDDO  
                      ENDDO  
                                                                                          
        ENDDO                                                                                                           !END beta
     ENDDO                                                                                                              !END alpha
     
!     call cputim(tim_tmp2_D)
!     time_cal_D = time_cal_D + (tim_tmp2_D-tim_tmp1_D)
         
    index = 0
    DO l = 0, lmax2n
       DO m = -l,l
          index = index +1
          DO alpha = 1,3
             DO beta = 1,3
              !  DO lmd = 0,1
              !     DO lmdp = 0,1
                   
                      DO num = nemin, nemax                                                                                                
                         buf_sph_val0(alpha,beta,num) = buf_sph_val0(alpha,beta,num) + &
                        (el(l)-e(num))*( CONJG(CAPD0LM(alpha,beta,index,num))*ABLM(0,index,num) + &
                            CONJG(CAPD1LM(alpha,beta,index,num))*ABLM(1,index,num)*PEI(l) ) + &
                            CONJG(CAPD0LM(alpha,beta,index,num))*ABLM(1,index,num)
                                              
                      ENDDO   
!                       !For the moment this below part is only for C-Diamond test
      before_lo=buf_sph_val0(1,1,1)
                      DO num = nemin,nemax
                         DO jlo = 1, ilo(l)
                            IF(loor(jlo,l)) THEN
                              buf_sph_val0(alpha,beta,num) = buf_sph_val0(alpha,beta,num) + &
                              ( CONJG(CAPD0LM(alpha,beta,index,num))*jlaplo1_02(jlo,l) + &
                                CONJG(CAPD1LM(alpha,beta,index,num))*jlaplo1_12(jlo,l) )* &
                                (elo(l,jlo)-e(num))*clm(index,num,jlo)
                            ENDIF
                         ENDDO
                      ENDDO
                         
              !     ENDDO !lmd
              !  ENDDO   !lmdp
                
             ENDDO
          ENDDO
!          if(index.lt.5.and.abs(buf_sph_val0(1,1,1)).lt.1.d0.and.abs(buf_sph_val0(1,1,1)).gt.0.1d0) then
!                  write(6,*) index,before_lo,buf_sph_val0(1,1,1)
!                  write(6,*) (el(l)-e(1)),CONJG(CAPD0LM(1,1,index,1)),ABLM(0,index,1)
!                  write(6,*)CONJG(CAPD1LM(1,1,index,1)),ABLM(1,index,1),PEI(l)
!          endif       

       ENDDO
    ENDDO
         
    DO sp = 0, 2
       IF(NOTRI(sp,1,1) .lt. 0) CYCLE
       DO tp = -sp, sp
          DO t = -1,1
             DO u = -1,1
                IF((-tp+t+u) .ne. 0) CYCLE
                DO l = 0, lmax2n
                   DO lp = 0, lmax2n
                      IF(NOTRI(lp,l,sp) .lt. 0) CYCLE
                      DO m = -l,l
                         index1 = l*(l+1)+m+1
                         DO mp = -lp,lp
                            index2 = lp*(lp+1) + mp+ 1
                            IF((-mp+m+tp) .ne. 0) CYCLE
                            gaunt1 = GAUNT(sp,1,1,tp,t,u)
                            gaunt2 = GAUNT(lp,l,sp,mp,m,tp)                            
                            DO alpha = 1, 3
                               DO beta = 1,3
                                  DO lmd = 0, 1
                                     DO lmdp = 0,1
                                        
                                        DO num = nemin, nemax                                         
                                         
                                            buf_sph_val1(alpha,beta,num) = buf_sph_val1(alpha,beta,num) + &
                                            conjg(ABLM(lmd,index1,num))*ABLM(lmdp,index2,num)*gaunt1*gaunt2*&
                                            conjg(c_am(alpha,t))*conjg(c_am(beta,u))* &
                                            ((el(lp)-e(num))*jlpl(lmd,lmdp,l,lp) + delta(lmdp,1)*jlpl(lmd,0,l,lp))
                                            
                                        ENDDO   
                                     ENDDO
                                  ENDDO                                  
                                  
                                  DO num = nemin, nemax
                                     DO jlo = 1, ilo(l)
                                        IF(loor(jlo,l)) THEN
                                             buf_sph_val1(alpha,beta,num) = buf_sph_val1(alpha,beta,num) + &
                                           !  buf_sph_val1_buf1(alpha,beta,num) = buf_sph_val1_buf1(alpha,beta,num) + &
                                             gaunt1*gaunt2*conjg(c_am(alpha,t))*conjg(c_am(beta,u))*&
                                          conjg(clm(index1,num,jlo))*(el(lp)-e(num))* &
                                          ( ABLM(0,index2,num)*jr3upu_02(jlo,l,lp) + &
                                                ABLM(1,index2,num)*jr3upu_12(jlo,l,lp) ) + &
                                           gaunt1*gaunt2*conjg(c_am(alpha,t))*conjg(c_am(beta,u))*&     
                                          conjg(clm(index1,num,jlo))*ABLM(1,index2,num)*jr3upu_02(jlo,l,lp)
                                        ENDIF
                                     ENDDO
                                  ENDDO
                                  
                                  DO num = nemin, nemax
                                     DO jlop = 1, ilo(lp)
                                        IF(loor(jlop,lp)) THEN
                                             buf_sph_val1(alpha,beta,num) = buf_sph_val1(alpha,beta,num) + &
                                           ! buf_sph_val1_buf2(alpha,beta,num) = buf_sph_val1_buf2(alpha,beta,num) + & 
                                             gaunt1*gaunt2*conjg(c_am(alpha,t))*conjg(c_am(beta,u))*&
                                           clm(index2,num,jlop)*(elo(lp,jlop)-e(num))* &
                                          (conjg( ABLM(0,index1,num) )*jr3upu_20(jlop,lp,l) + &
                                           conjg( ABLM(1,index1,num) )*jr3upu_21(jlop,lp,l))     
                                        ENDIF
                                     ENDDO
                                  ENDDO
                                  
                                  DO num = nemin, nemax
                                     DO jlo = 1, ilo(l)
                                        IF(loor(jlo,l)) THEN
                                           DO jlop = 1, ilo(lp)
                                              IF(loor(jlop,lp)) THEN
                                                 buf_sph_val1(alpha,beta,num) = buf_sph_val1(alpha,beta,num) + &
                                           !     buf_sph_val1_buf3(alpha,beta,num) = buf_sph_val1_buf3(alpha,beta,num) + & 
                                                 gaunt1*gaunt2*conjg(c_am(alpha,t))*conjg(c_am(beta,u))*&
                                                 conjg(clm(index1,num,jlo))*clm(index2,num,jlop)* &
                                                 ! jr3upu_22 stores the derivative on its second index.
                                                 ! G2 requires a differentiated bra, with the ket energy.
                                                 (elo(lp,jlop)-e(num))*jr3upu_22(jlop,jlo,lp,l)
                                              ENDIF
                                           ENDDO                                           
                                        ENDIF
                                     ENDDO
                                  ENDDO
                                  
                                  
                               ENDDO                                                                                             !beta
                            ENDDO                                                                                                !alpha
                         ENDDO                                                                                                   !mp
                      ENDDO                                                                                                      !END m
                   ENDDO                                                                                                         !END lp
                ENDDO                                                                                                            !END l
             ENDDO                                                                                                               !END u
          ENDDO                                                                                                                  !END t
       ENDDO                                                                                                                     !END tp
    ENDDO                                                                                                                        !END sp 
    
    DO l = 0, lmax2n
       DO v = -1,1,2
          IF((l+v) .lt. 0)   CYCLE
          DO lp = 0, lmax2n
             IF(NOTRI(lp,1,l+v) .lt. 0) CYCLE              
             DO m = -l,l
                index1 = l*(l+1)+m+1
                DO t = -1,1
                   IF((m+t) .lt. -(l+v)) CYCLE
                   IF((m+t) .gt.  (l+v)) CYCLE
                   DO mp = -lp,lp
                      index2 = lp*(lp+1)+mp+1
                      DO u = -1,1
                         IF((-mp+u+m+t) .ne. 0) CYCLE
                         gaunt1 = GAUNT(lp,1,l+v,mp,u,m+t)                         
                         CALL c_factor(l,m,ca_st_lm)
                         DO alpha = 1, 3
                            DO beta = 1,3                               
                               DO lmd = 0,1
                                  DO lmdp = 0,1
                                  
                                     DO num = nemin, nemax                                           
                                    !    buf_sph_val2(alpha,beta,num) = buf_sph_val2(alpha,beta,num) + &
                                    !    (el(lp)-e(num))*conjg(ABLM(lmd,index1,num))*ABLM(lmdp,index2,num)*&
                                    !  ( conjg(ca_st_lm(beta,v,t))*conjg(c_am(alpha,u)) +  &
                                    !    conjg(ca_st_lm(alpha,v,t))*conjg(c_am(beta,u)) )*gaunt1*hf*Ilpl(lmd,lmdp,l,lp)+&
                                        
                                    !    delta(lmdp,1)*conjg(ABLM(lmd,index1,num))*ABLM(lmdp,index2,num)*&
                                    !    (conjg(ca_st_lm(beta,v,t))*conjg(c_am(alpha,u)) +  &
                                    !    conjg(ca_st_lm(alpha,v,t))*conjg(c_am(beta,u)))*gaunt1*hf*Ilpl(lmd,0,l,lp)
                                        
                                        buf_sph_val2(alpha,beta,num) = buf_sph_val2(alpha,beta,num) + &
                                       hf* conjg(ABLM(lmd,index1,num))*ABLM(lmdp,index2,num)*gaunt1* &
                                       (  conjg(ca_st_lm(alpha,v,t))*conjg(c_am(beta,u))  +  & 
                                          conjg(ca_st_lm(beta,v,t))* conjg(c_am(alpha,u))  )* &
                                       ( ( el(lp)-e(num) )*Ilpl(lmd,lmdp,l,lp)+delta(lmdp,1)*Ilpl(lmd,0,l,lp) )   
                                     ENDDO
                                     
                                  ENDDO
                               ENDDO
                            ENDDO
                         ENDDO
                   !     print*,l,v,lp,m,t,mp,u
                   !     print*,buf_sph_val2(1,1,1) 
                      ENDDO
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO  
    
    DO alpha = 1,3
       DO beta = 1,3
          DO num = nemin,nemax
             ! Keep the three lines of the spherical Blaha kernel available
             ! separately.  Their sum is unchanged below.  This split is
             ! needed for the off-shell DMFT audit: in DFT only the on-shell
             ! H-epsilon O combination is observable, whereas M_Sigma probes
             ! the overlap derivative itself and can reveal a term that was
             ! legitimately eliminated in the on-shell derivation.
             buf_sph_match(alpha,beta,num)=buf_sph_val0(alpha,beta,num)
             buf_sph_radial(alpha,beta,num)=buf_sph_val1(alpha,beta,num)
             buf_sph_angular(alpha,beta,num)=buf_sph_val2(alpha,beta,num)
             buf_sph_val_fin(alpha,beta,num) = buf_sph_val0(alpha,beta,num)+buf_sph_val1(alpha,beta,num) +&
                                               buf_sph_val2(alpha,beta,num)                                                                                             
          ENDDO
       ENDDO
    ENDDO

    ! The first three lines of Eq. (G2), above, use only the spherical
    ! augmentation data.  Its last three lines require the nonspherical
    ! matrices written by LAPW1/dftKS to case.nsh.  They are part of
    ! H_{,epsilon}, hence only the rho (Hamiltonian) contraction requests
    ! them; the rhoE pass below needs only O_{,epsilon}.
    IF (do_nsp) CALL val_nsp_stress(nemin,nemax,jatom,ALM,BLM,clm, &
         CAPD0LM,CAPD1LM,buf_nsp_match,buf_nsp_radial,buf_nsp_angular)
!    write(6,*) buf_sph_val_fin(1,1,1),buf_sph_val0(1,1,1)
!  977 continue        
        
  RETURN
END SUBROUTINE val_corr_stress


! Last three (nonspherical-potential) lines of Blaha et al. Eq. (G2).
!
! This routine deliberately consumes the same case.nsh matrices as Force2.
! The radial basis ordering is
!
!                 (u, udot, LO_1, ..., LO_n),
!
! and the first matrix index is the bra radial function.  Vts and Wts
! already contain the Hamiltonian Gaunt coefficient.  VPts and WPts retain
! only the real/cubic-harmonic phase; the new strain angular factors below
! replace the original Gaunt coefficient and therefore must use these bare
! phase-weighted matrices.  Dividing Vts or Wts by the original Gaunt is not
! valid: that Gaunt may vanish for channels made nonzero by the two inserted
! rank-one strain tensors.
SUBROUTINE val_nsp_stress(nemin,nemax,jatom,ALM,BLM,CLM, &
     CAPD0LM,CAPD1LM,buf_match,buf_radial,buf_angular)
  use param, only: lmax2,lmax2n,nume,nloat
  use defs, only: zeroc
  use lo, only: ilo
  use w_atpar, only: Vts,VPts,WPts,lpv,lv,mpv,mv,Nhmx,nsp_geometry
  implicit none
  INTEGER, intent(in) :: nemin,nemax,jatom
  COMPLEX*16, intent(in) :: ALM((lmax2+1)*(lmax2+1),nume), &
       BLM((lmax2+1)*(lmax2+1),nume), &
       CLM((lmax2+1)*(lmax2+1),nume,nloat)
  COMPLEX*16, intent(in) :: CAPD0LM(3,3,(lmax2+1)*(lmax2+1),nume), &
       CAPD1LM(3,3,(lmax2+1)*(lmax2+1),nume)
  COMPLEX*16, intent(out) :: buf_match(3,3,nume), &
       buf_radial(3,3,nume),buf_angular(3,3,nume)

  INTEGER :: ih,l,lp,ly,lpy,index1,index2,iterm
  INTEGER :: alpha,beta,num,ip,iq,np,nq,jlo
  COMPLEX*16 :: cbra(2+nloat),cket(2+nloat)
  COMPLEX*16 :: afac,bfac,wmat,vmat

  buf_match=zeroc
  buf_radial=zeroc
  buf_angular=zeroc

  ! Fourth line: the strain derivative D of the APW matching coefficients
  ! contracted with the ordinary nonspherical Hamiltonian matrix.  D acts
  ! on the ket APW coefficients, hence only columns u and udot are needed;
  ! every APW and LO bra row is retained.
  DO ih=1,Nhmx(jatom)
     lpy=lpv(ih,jatom)*(lpv(ih,jatom)+1)+mpv(ih,jatom)+1
     ly =lv(ih,jatom) *(lv(ih,jatom)+1) +mv(ih,jatom) +1
     np=2+ilo(lpv(ih,jatom))
     DO num=nemin,nemax
        cbra=zeroc
        cbra(1)=ALM(lpy,num)
        cbra(2)=BLM(lpy,num)
        DO jlo=1,ilo(lpv(ih,jatom))
           cbra(2+jlo)=CLM(lpy,num,jlo)
        ENDDO
        afac=zeroc
        bfac=zeroc
        DO ip=1,np
           afac=afac+conjg(cbra(ip))*Vts(ip,1,ih,jatom)
           bfac=bfac+conjg(cbra(ip))*Vts(ip,2,ih,jatom)
        ENDDO
        DO beta=1,3
           DO alpha=1,3
              buf_match(alpha,beta,num)=buf_match(alpha,beta,num)+ &
                   afac*CAPD0LM(alpha,beta,ly,num)+ &
                   bfac*CAPD1LM(alpha,beta,ly,num)
           ENDDO
        ENDDO
     ENDDO
  ENDDO

  ! Fifth and sixth lines: the geometry is independent of k and of the
  ! density matrix.  Build it once per atom, then retain only the nonzero
  ! (record,bra-m,ket-m) triples.  This turns the dominant repeated Gaunt
  ! algebra into a small matrix contraction over the current pseudo-bands.
  IF (.NOT.nsp_geometry(jatom)%ready) CALL build_nsp_geometry(jatom)
  DO iterm=1,nsp_geometry(jatom)%nterm
     ih=nsp_geometry(jatom)%ih(iterm)
     index1=nsp_geometry(jatom)%ibra(iterm)
     index2=nsp_geometry(jatom)%iket(iterm)
     l=lpv(ih,jatom)
     lp=lv(ih,jatom)
     np=2+ilo(l)
     nq=2+ilo(lp)
     DO num=nemin,nemax
        cbra=zeroc
        cket=zeroc
        cbra(1)=ALM(index1,num)
        cbra(2)=BLM(index1,num)
        cket(1)=ALM(index2,num)
        cket(2)=BLM(index2,num)
        DO jlo=1,ilo(l)
           cbra(2+jlo)=CLM(index1,num,jlo)
        ENDDO
        DO jlo=1,ilo(lp)
           cket(2+jlo)=CLM(index2,num,jlo)
        ENDDO
        wmat=zeroc
        vmat=zeroc
        DO iq=1,nq
           DO ip=1,np
              wmat=wmat+conjg(cbra(ip))*WPts(ip,iq,ih,jatom)*cket(iq)
              vmat=vmat+conjg(cbra(ip))*VPts(ip,iq,ih,jatom)*cket(iq)
           ENDDO
        ENDDO
        buf_radial(:,:,num)=buf_radial(:,:,num)+ &
             nsp_geometry(jatom)%radial(:,:,iterm)*wmat
        buf_angular(:,:,num)=buf_angular(:,:,num)+ &
             nsp_geometry(jatom)%angular(:,:,iterm)*vmat
     ENDDO
  ENDDO
END SUBROUTINE val_nsp_stress


SUBROUTINE build_nsp_geometry(jatom)
  ! Precompute the k-independent three-Gaunt tensors in the fifth and sixth
  ! lines of Eq. (G2).  The cache is sparse in (radial channel,m_bra,m_ket)
  ! and is owned by w_atpar together with the case.nsh matrices whose labels
  ! define it.
  use param, only: lmax2n
  use defs, only: zeroc
  use com_mpi, only: myrank,master
  use w_atpar, only: VPts,WPts,lpv,lv,llv,mmv,Nhmx,nsp_geometry
  implicit none
  REAL*8, parameter :: hf=0.5d0
  INTEGER, intent(in) :: jatom
  INTEGER :: ih,l,lp,mbra,mket,index1,index2,maxterm,iterm
  INTEGER :: alpha,beta,sp,tp,t,u,v,v1,lint,nu
  INTEGER, external :: NOTRI
  REAL*8, external :: GAUNT
  REAL*8 :: gaunt1,gaunt_pot,gaunt_couple
  REAL*8 :: channel_diff,channel_scale
  LOGICAL :: duplicate_radial_channel
  COMPLEX*16 :: c_am(3,-1:1),ca_st_lm(3,-1:1,-1:1)
  COMPLEX*16 :: radial_geom(3,3),angular_geom(3,3)
  COMPLEX*16 :: radial_factor,angular_factor

  IF (nsp_geometry(jatom)%ready) RETURN
  call c_alpha_m(c_am)

  ! First determine a strict upper bound.  The same radial channel is
  ! repeated consecutively in case.nsh for the original Hamiltonian-allowed
  ! magnetic pairs.  Its phase-weighted matrices and (l,L,M,l') labels are
  ! identical, so comparing with the preceding record removes the replica.
  maxterm=0
  DO ih=1,Nhmx(jatom)
     l=lpv(ih,jatom)
     lp=lv(ih,jatom)
     duplicate_radial_channel=.FALSE.
     IF (ih.GT.1) THEN
     IF (lpv(ih-1,jatom).EQ.l .AND. &
          lv(ih-1,jatom).EQ.lp .AND. &
          llv(ih-1,jatom).EQ.llv(ih,jatom) .AND. &
          mmv(ih-1,jatom).EQ.mmv(ih,jatom)) THEN
        channel_diff=SUM(ABS(VPts(:,:,ih-1,jatom)-VPts(:,:,ih,jatom)))+ &
             SUM(ABS(WPts(:,:,ih-1,jatom)-WPts(:,:,ih,jatom)))
        channel_scale=SUM(ABS(VPts(:,:,ih,jatom)))+ &
             SUM(ABS(WPts(:,:,ih,jatom)))
        duplicate_radial_channel= &
             channel_diff.LE.1.d-12*(1.d0+channel_scale)
     ENDIF
     ENDIF
     IF (.NOT.duplicate_radial_channel) &
          maxterm=maxterm+(2*l+1)*(2*lp+1)
  ENDDO

  allocate(nsp_geometry(jatom)%ih(maxterm))
  allocate(nsp_geometry(jatom)%ibra(maxterm))
  allocate(nsp_geometry(jatom)%iket(maxterm))
  allocate(nsp_geometry(jatom)%radial(3,3,maxterm))
  allocate(nsp_geometry(jatom)%angular(3,3,maxterm))
  nsp_geometry(jatom)%radial=zeroc
  nsp_geometry(jatom)%angular=zeroc

  iterm=0
  DO ih=1,Nhmx(jatom)
     l=lpv(ih,jatom)
     lp=lv(ih,jatom)
     duplicate_radial_channel=.FALSE.
     IF (ih.GT.1) THEN
     IF (lpv(ih-1,jatom).EQ.l .AND. &
          lv(ih-1,jatom).EQ.lp .AND. &
          llv(ih-1,jatom).EQ.llv(ih,jatom) .AND. &
          mmv(ih-1,jatom).EQ.mmv(ih,jatom)) THEN
        channel_diff=SUM(ABS(VPts(:,:,ih-1,jatom)-VPts(:,:,ih,jatom)))+ &
             SUM(ABS(WPts(:,:,ih-1,jatom)-WPts(:,:,ih,jatom)))
        channel_scale=SUM(ABS(VPts(:,:,ih,jatom)))+ &
             SUM(ABS(WPts(:,:,ih,jatom)))
        duplicate_radial_channel= &
             channel_diff.LE.1.d-12*(1.d0+channel_scale)
     ENDIF
     ENDIF
     IF (duplicate_radial_channel) CYCLE

     DO mbra=-l,l
        index1=l*(l+1)+mbra+1
        CALL c_factor(l,mbra,ca_st_lm)
        DO mket=-lp,lp
           index2=lp*(lp+1)+mket+1
           radial_geom=zeroc
           DO sp=0,2
              IF (NOTRI(sp,1,1).LT.0) CYCLE
              DO tp=-sp,sp
                 DO t=-1,1
                    DO u=-1,1
                       IF (-tp+t+u.NE.0) CYCLE
                       gaunt1=GAUNT(sp,1,1,tp,t,u)
                       IF (ABS(gaunt1).LT.1.d-30) CYCLE
                       DO v1=0,llv(ih,jatom)+lp
                          IF (NOTRI(v1,llv(ih,jatom),lp).LT.0) CYCLE
                          IF (NOTRI(l,sp,v1).LT.0) CYCLE
                          DO v=-v1,v1
                             IF (-v+mmv(ih,jatom)+mket.NE.0) CYCLE
                             gaunt_pot=GAUNT(v1,llv(ih,jatom),lp,v, &
                                  mmv(ih,jatom),mket)
                             IF (ABS(gaunt_pot).LT.1.d-30) CYCLE
                             IF (-mbra+tp+v.NE.0) CYCLE
                             gaunt_couple=GAUNT(l,sp,v1,mbra,tp,v)
                             IF (ABS(gaunt_couple).LT.1.d-30) CYCLE
                             DO beta=1,3
                                DO alpha=1,3
                                   radial_factor=c_am(alpha,t)* &
                                        c_am(beta,u)*gaunt1* &
                                        gaunt_pot*gaunt_couple
                                   radial_geom(alpha,beta)= &
                                        radial_geom(alpha,beta)+radial_factor
                                ENDDO
                             ENDDO
                          ENDDO
                       ENDDO
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO

           angular_geom=zeroc
           DO v=-1,1,2
              ! l+1 is an intermediate harmonic even at the radial cutoff.
              IF (l+v.LT.0) CYCLE
              DO t=-1,1
                 IF (mbra+t.LT.-(l+v) .OR. mbra+t.GT.l+v) CYCLE
                 DO tp=-1,1
                    DO lint=0,llv(ih,jatom)+lp
                       IF (NOTRI(lint,llv(ih,jatom),lp).LT.0) CYCLE
                       IF (NOTRI(lint,1,l+v).LT.0) CYCLE
                       DO nu=-lint,lint
                          IF (-nu+mmv(ih,jatom)+mket.NE.0) CYCLE
                          gaunt_pot=GAUNT(lint,llv(ih,jatom),lp,nu, &
                               mmv(ih,jatom),mket)
                          IF (ABS(gaunt_pot).LT.1.d-30) CYCLE
                          IF (-nu+tp+mbra+t.NE.0) CYCLE
                          gaunt_couple=GAUNT(lint,1,l+v,nu,tp,mbra+t)
                          IF (ABS(gaunt_couple).LT.1.d-30) CYCLE
                          DO beta=1,3
                             DO alpha=1,3
                                angular_factor=hf*( &
                                     conjg(c_am(alpha,tp))* &
                                     conjg(ca_st_lm(beta,v,t))+ &
                                     conjg(c_am(beta,tp))* &
                                     conjg(ca_st_lm(alpha,v,t)))
                                angular_geom(alpha,beta)= &
                                     angular_geom(alpha,beta)+ &
                                     angular_factor*gaunt_pot*gaunt_couple
                             ENDDO
                          ENDDO
                       ENDDO
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO

           IF (MAX(MAXVAL(ABS(radial_geom)), &
                MAXVAL(ABS(angular_geom))).LT.1.d-30) CYCLE
           iterm=iterm+1
           nsp_geometry(jatom)%ih(iterm)=ih
           nsp_geometry(jatom)%ibra(iterm)=index1
           nsp_geometry(jatom)%iket(iterm)=index2
           nsp_geometry(jatom)%radial(:,:,iterm)=radial_geom
           nsp_geometry(jatom)%angular(:,:,iterm)=angular_geom
        ENDDO
     ENDDO
  ENDDO

  nsp_geometry(jatom)%nterm=iterm
  nsp_geometry(jatom)%ready=.TRUE.
  IF (myrank.EQ.master) WRITE(6,'(A,3I8)') &
       ':DMFT_G2_GEOMETRY_CACHE ',jatom,iterm,maxterm
END SUBROUTINE build_nsp_geometry

SUBROUTINE val_kinstrmt(nemin,nemax,alm,blm,clm,F_buf_total,RMT)
  !  use param, only:lmax=>lmax2
    use param
    use atspdt,    only:p,dp,pe,dpe
    use defs  ,    only:zeroc   
    USE stress_kb, only: tc1s00,tc1s01,tc1s10,tc1s11,tc2s00,tc2s01,&
                         tc2s10,tc2s11,tc3s00,tc3s01,tc3s10,tc3s11, &
                        juddu_02,juddu_20,juddu_12,juddu_21,juddu_22, &
                        judu_02, judu_20, judu_12, judu_21, judu_22,  &
                        juu_02,  juu_20,  juu_12,  juu_21,  juu_22 
                
    use lo      ,  only: ilo,loor
    use com     ,  only: REL
!   use gaunt_num
    implicit none
!     integer, parameter         :: lmax = lmax2
    integer                    :: lmax 
    integer, intent(in)        :: nemin,nemax
    COMPLEX*16, intent(in)     :: alm((lmax2+1)*(lmax2+1),nume), &
                                  blm((lmax2+1)*(lmax2+1),nume), &
                                  clm((lmax2+1)*(lmax2+1),nume,nloat)                                  
    REAL*8, intent(in)         :: RMT                              
    REAL*8                     :: FF(0:1,0:lmax2),DFF(0:1,0:lmax2)
    COMPLEX*16, intent(out)    :: F_buf_total(3,3,nume)
    integer                    :: l,ll,lp,llp,m,mp,t,tp,u,s,sp,v1,v2,v, &
                                   jlo,jlop
    integer                    :: alpha,beta, num, lmda,lmdap
    integer                    :: notri1,notri2, index1, index2, index
    integer, external          :: NOTRI
    REAL*8,  external          :: GAUNT
    real*8                     :: gaunt1,gaunt2,gaunt3,RMT2
    real*8                     :: g1_deriv_scale
    integer                    :: g1_legacy_status
    character(len=32)          :: g1_legacy_env
    REAL*8, PARAMETER          :: ohalf = 1.D0/2.D0    
    complex*16                 :: c_am(3,-1:1), ca_st_lm(3,-1:1,-1:1), &
                                  ca_st_lm11(3,-1:1,-1:1)
    complex*16                 :: F1_buf(3,3,nume), F2_buf(3,3,nume), &
                                  F3_buf(3,3,nume), F4_buf(3,3,nume), &
                                  F5_buf(3,3,nume), F6_buf(3,3,nume), &                                  
                                  sum_buf(4,nume),buf_tempo , &
                                  F6_buf1(3,3,nume),F6_buf2(3,3,nume), &
                                  F6_buf3(3,3,nume)

                                  
!=====================    
    lmax = lmax2n        !at least lmax = 7 is needed
!=====================                                      

    ! Equation G1 carries one overall factor 1/2.  In the historical
    ! implementation the first-radial-derivative part duplicated its angular
    ! factor and consequently lost this 1/2, while the accompanying -u/r**2
    ! piece retained it.  Use the corrected factor by default, matching the
    ! patched WIEN2k-24 stress kernel.  The environment switch is diagnostic
    ! only and permits exact reproduction of the legacy result.
    g1_deriv_scale = 0.5d0
    g1_legacy_env = ''
    call get_environment_variable('WIEN_STRESS_G1_LEGACY',g1_legacy_env, &
         status=g1_legacy_status)
    if (g1_legacy_status.eq.0 .and. len_trim(g1_legacy_env).gt.0 .and. &
         trim(g1_legacy_env).ne.'0') g1_deriv_scale = 1.d0
                                  
    DO l = 0, lmax
       FF(0,l)   = P(l)
       FF(1,l)   = PE(l)
       DFF(0,l)  = DP(l)
       DFF(1,l)  = DPE(l)       
    ENDDO


   
    call c_alpha_m(c_am)        
    buf_tempo  = zeroc
    F1_buf  = zeroc
    F2_buf  = zeroc
    F3_buf  = zeroc
    F4_buf  = zeroc
    F5_buf  = zeroc
    F6_buf  = zeroc
    F6_buf1  = zeroc
    F6_buf2  = zeroc
    F6_buf3  = zeroc
    sum_buf = zeroc
    
  !=================================valence kinetice stress inside atomic sphere using wave function not density==============================                
 ! go to 7777
       DO sp = 0, 2
          IF(NOTRI(sp,1,1) .lt. 0) CYCLE
          DO l = 0, lmax
             DO lp = 0, lmax
                IF(NOTRI(l,lp,sp) .lt. 0) CYCLE
                DO tp = -sp,sp
                   DO u = -1,1
                      DO t = -1,1
                         IF((-tp+u+t) .ne. 0) CYCLE
                         DO m = -l,l
                            index1 = l*(l+1)+m+1
                            DO mp = -lp,lp
                               index2 = lp*(lp+1)+mp+1
!                                IF (index1 .eq. 1 .or. index2 .eq. 1) CYCLE
                               IF((-m+mp+tp) .ne. 0) CYCLE
                               gaunt1 = GAUNT(sp,1,1,tp,u,t)
                               gaunt2 = GAUNT(l,lp,sp,m,mp,tp)
                               DO alpha = 1,3
                                  DO beta = 1,3
                                     
                                     DO num = nemin,nemax
                                        !HERE IS THE FORMULA F1
                                        F1_buf(alpha,beta,num) = F1_buf(alpha,beta,num) + &                                        
                                        (c_am(alpha,t)*c_am(beta,u)+c_am(beta,t)*c_am(alpha,u))*ohalf*gaunt1*gaunt2*&
                                        (ALM(index2,num)*(conjg(ALM(index1,num))*TC1S00(l,lp)+&
                                        +conjg(BLM(index1,num))*TC1S10(l,lp))+BLM(index2,num)*&
                                        (conjg(ALM(index1,num))*TC1S01(l,lp)+conjg(BLM(index1,num))*TC1S11(l,lp)) )
                                        
                                       DO jlo = 1, ilo(l) 
                                         ! IF(.NOT. loor(jlo,l)) CYCLE
                                          IF(loor(jlo,l)) THEN
                                          F1_buf(alpha,beta,num) = F1_buf(alpha,beta,num) + & 
                                         ! F6_buf1(alpha,beta,num) = F6_buf1(alpha,beta,num) + &                                           
                                          (c_am(alpha,t)*c_am(beta,u)+c_am(beta,t)*c_am(alpha,u))*ohalf*gaunt1*gaunt2*&
                                          conjg(clm(index1,num,jlo))*( ALM(index2,num)*juddu_20(jlo,l,lp) + &
                                          BLM(index2,num)*juddu_21(jlo,l,lp)  )
                                          ENDIF
                                       ENDDO
                                       
                                       DO jlop = 1, ilo(lp)
                                         ! IF(.NOT.loor(jlop,lp)) CYCLE
                                          IF(loor(jlop,lp)) THEN
                                          F1_buf(alpha,beta,num) = F1_buf(alpha,beta,num) + & 
                                         ! F6_buf2(alpha,beta,num) = F6_buf2(alpha,beta,num) + &
                                          (c_am(alpha,t)*c_am(beta,u)+c_am(beta,t)*c_am(alpha,u))*ohalf*gaunt1*gaunt2*&
                                          clm(index2,num,jlop)*( conjg(ALM(index1,num))*juddu_02(jlop,lp,l) + &
                                          conjg(BLM(index1,num))*juddu_12(jlop,lp,l)  )
                                          ENDIF
                                       ENDDO
                                       
                                       DO jlo = 1, ilo(l)
                                         ! IF(.NOT. loor(jlo,l)) CYCLE
                                         IF(loor(jlo,l)) THEN
                                          DO jlop = 1, ilo(lp)
                                            ! IF(.NOT.loor(jlop,lp)) CYCLE
                                             IF(loor(jlop,lp)) THEN
                                             F1_buf(alpha,beta,num) = F1_buf(alpha,beta,num) + &                                               
                                          !   F6_buf3(alpha,beta,num) = F6_buf3(alpha,beta,num) + &
                                             (c_am(alpha,t)*c_am(beta,u)+c_am(beta,t)*c_am(alpha,u))*ohalf*gaunt1*gaunt2*&
                                             conjg(clm(index1,num,jlo))*clm(index2,num,jlop)*juddu_22(jlo,jlop,l,lp)
                                             ENDIF
                                          ENDDO
                                         ENDIF 
                                       ENDDO                                                                              
                                       
                                     ENDDO                                                                                !END num
                                  ENDDO                                                                                   !END beta
                               ENDDO                                                                                      !END alpha
                            ENDDO                                                                                         !END mp
                         ENDDO                                                                                            !END m
                      ENDDO                                                                                               !END t
                   ENDDO                                                                                                  !END u
                ENDDO                                                                                                     !END tp
             ENDDO                                                                                                        !END lp
          ENDDO                                                                                                           !END l
       ENDDO                                                                                                              !END sp
       
       DO l = 0, lmax
          DO lp = 0, lmax
             DO v = -1,1,2
                IF((1+v) .lt. 0) CYCLE
                IF(NOTRI(l,lp,1+v) .lt. 0) CYCLE
                DO m = -l,l
                   index1 = l*(l+1)+m+1
                   DO mp = -lp,lp
                      index2 = lp*(lp+1)+mp+1
!                       IF (index1 .eq. 1 .or. index2 .eq. 1) CYCLE
                      DO u = -1,1
                         DO t = -1,1                         
                            IF((t+u) .lt. -(1+v)) CYCLE
                            IF((t+u) .gt. (1+v)) CYCLE
                            IF((-m+mp+t+u) .ne. 0) CYCLE
                                                        
                            gaunt1 = GAUNT(l,lp,1+v,m,mp,t+u)                             
                            CALL c_factor(1,t,ca_st_lm)                            
                            DO alpha = 1,3
                               DO beta = 1,3
                                  DO num = nemin,nemax
                                  
                                     !HERE IS THE SUM OF F2
                                     F2_buf(alpha,beta,num) = F2_buf(alpha,beta,num) + &
                                     (c_am(alpha,t)*ca_st_lm(beta,v,u)+c_am(beta,t)*ca_st_lm(alpha,v,u))*ohalf*gaunt1*&
                                     (ALM(index2,num)*(conjg(ALM(index1,num))*TC2S00(l,lp)+&
                                    +conjg(BLM(index1,num))*TC2S10(l,lp))+BLM(index2,num)*&
                                    (conjg(ALM(index1,num))*TC2S01(l,lp)+conjg(BLM(index1,num))*TC2S11(l,lp)))
                                    
                                     DO jlo = 1, ilo(l)
                                        IF(loor(jlo,l)) THEN
                                        F2_buf(alpha,beta,num) = F2_buf(alpha,beta,num) + &
                                       ! F6_buf1(alpha,beta,num) = F6_buf1(alpha,beta,num) + &
                                        (c_am(alpha,t)*ca_st_lm(beta,v,u)+c_am(beta,t)*ca_st_lm(alpha,v,u))*ohalf*gaunt1*&
                                        conjg( clm(index1,num,jlo) )*( ALM(index2,num)*judu_20(jlo,l,lp) + &
                                        BLM(index2,num)*judu_21(jlo,l,lp) )
                                        ENDIF
                                     ENDDO
                                     
                                     DO jlop = 1, ilo(lp)
                                        IF(loor(jlop,lp)) THEN
                                           F2_buf(alpha,beta,num) = F2_buf(alpha,beta,num) + & 
                                         ! F6_buf2(alpha,beta,num) = F6_buf2(alpha,beta,num) + &
                                           (c_am(alpha,t)*ca_st_lm(beta,v,u)+c_am(beta,t)*ca_st_lm(alpha,v,u))*ohalf*gaunt1*&
                                            clm(index2,num,jlop) *(conjg( ALM(index1,num))*judu_02(jlop,lp,l) + &
                                              conjg( BLM(index1,num))*judu_12(jlop,lp,l) )
                                        ENDIF
                                     ENDDO
                                     
                                     DO jlo = 1, ilo(l)
                                        IF(loor(jlo,l)) THEN
                                          DO jlop = 1, ilo(lp)
                                             IF(loor(jlop,lp)) THEN
                                                F2_buf(alpha,beta,num) = F2_buf(alpha,beta,num) + & 
                                             !  F6_buf3(alpha,beta,num) = F6_buf3(alpha,beta,num) + &
                                                (c_am(alpha,t)*ca_st_lm(beta,v,u)+c_am(beta,t)* &
                                                 ca_st_lm(alpha,v,u))*ohalf*gaunt1*&
                                                 conjg( clm(index1,num,jlo) )*clm(index2,num,jlop)*judu_22(jlo,jlop,l,lp)
                                             ENDIF
                                          ENDDO
                                        ENDIF
                                     ENDDO
                                     
                                  ENDDO                                                                                   !END num
                               ENDDO                                                                                      !END beta
                            ENDDO                                                                                         !END alpha                            
                         ENDDO                                                                                            !END t
                      ENDDO                                                                                               !END u
                   ENDDO                                                                                                  !END mp
                ENDDO                                                                                                     !END m
             ENDDO                                                                                                        !END v
          ENDDO                                                                                                           !END lp
       ENDDO                                                                                                              !END l
       
       
       DO lp = 0,lmax
          DO v = -1,1,2
             IF((lp+v) .lt. 0) CYCLE
             DO l = 0, lmax
                IF(NOTRI(l,1,lp+v) .lt. 0) CYCLE
                DO mp = -lp,lp
                   index2 = lp*(lp+1)+mp+1
                   DO u = -1,1                                          
                      IF((mp+u) .lt. -(lp+v)) CYCLE
                      IF((mp+u) .gt. (lp+v) ) CYCLE
                      DO m = -l,l
                         index1 = l*(l+1)+m+1
!                          IF (index1 .eq. 1 .or. index2 .eq. 1) CYCLE
                         DO t = -1,1
                            IF((-m+t+mp+u) .ne. 0) CYCLE
                                                        
                            gaunt1 = GAUNT(l,1,lp+v,m,t,mp+u)                                                         
                            CALL c_factor(lp,mp,ca_st_lm)                                                       
                            
                            DO alpha =1,3
                               DO beta = 1,3
                                  DO num = nemin,nemax
                                     !HERE IS THE FORMULA F3+F4-F5                                     
                                     F3_buf(alpha,beta,num) = F3_buf(alpha,beta,num) + &
                                     (c_am(alpha,t)*ca_st_lm(beta,v,u)+c_am(beta,t)*ca_st_lm(alpha,v,u) +&
                                      c_am(beta,t)*ca_st_lm(alpha,v,u)+c_am(alpha,t)*ca_st_lm(beta,v,u))*&
                                      ohalf*g1_deriv_scale*gaunt1*(ALM(index2,num)*(conjg(ALM(index1,num))*TC2S00(l,lp)+&
                                      +conjg(BLM(index1,num))*TC2S10(l,lp)) + BLM(index2,num)*&                                      
                                    (conjg(ALM(index1,num))*TC2S01(l,lp)+conjg(BLM(index1,num))*TC2S11(l,lp))) - &
                                     (c_am(beta,t)*ca_st_lm(alpha,v,u)+c_am(alpha,t)*ca_st_lm(beta,v,u))*ohalf*gaunt1*&
                                     (ALM(index2,num)*(conjg(ALM(index1,num))*TC3S00(l,lp)+&
                                     +conjg(BLM(index1,num))*TC3S10(l,lp))+BLM(index2,num)*&
                                     (conjg(ALM(index1,num))*TC3S01(l,lp)+conjg(BLM(index1,num))*TC3S11(l,lp)))
                                     
                                 !commented lines below are for testing purpose and F3_buf is zero in LAPW and LO case. 
                               !     F4_buf(alpha,beta,num) = F4_buf(alpha,beta,num) + &
                               !      (c_am(beta,t)*ca_st_lm(alpha,v,u)+c_am(alpha,t)*ca_st_lm(beta,v,u))*ohalf*gaunt1*&
                               !      (ALM(index2,num)*(conjg(ALM(index1,num))*TC2S00(l,lp)+&
                               !     +conjg(BLM(index1,num))*TC2S10(l,lp))+BLM(index2,num)*&
                               !     (conjg(ALM(index1,num))*TC2S01(l,lp)+conjg(BLM(index1,num))*TC2S11(l,lp)))
                                   
                               !     F5_buf(alpha,beta,num) = F5_buf(alpha,beta,num) - &
                               !     (c_am(beta,t)*ca_st_lm(alpha,v,u)+c_am(alpha,t)*ca_st_lm(beta,v,u))*ohalf*gaunt1*&
                               !     (ALM(index2,num)*(conjg(ALM(index1,num))*TC3S00(l,lp)+&
                               !      +conjg(BLM(index1,num))*TC3S10(l,lp))+BLM(index2,num)*&
                               !       (conjg(ALM(index1,num))*TC3S01(l,lp)+conjg(BLM(index1,num))*TC3S11(l,lp)))
                                      
                                 !     F3_buf(alpha,beta,num) = F3_buf(alpha,beta,num) + &
                                 !    (c_am(alpha,t)*ca_st_lm(beta,v,u)+c_am(beta,t)*ca_st_lm(alpha,v,u) )*ohalf*gaunt1* &
                                 !     (ALM(index2,num)*(conjg(ALM(index1,num))*TC2S00(l,lp)+&
                                 !     +conjg(BLM(index1,num))*TC2S10(l,lp)) + BLM(index2,num)*&                                      
                                 !   (conjg(ALM(index1,num))*TC2S01(l,lp)+conjg(BLM(index1,num))*TC2S11(l,lp)))                                     
                                      
                                      DO jlo = 1, ilo(l)
                                         IF(loor(jlo,l)) THEN
                                           F3_buf(alpha,beta,num) = F3_buf(alpha,beta,num) + &
                                           (c_am(alpha,t)*ca_st_lm(beta,v,u)+c_am(beta,t)*ca_st_lm(alpha,v,u) +&
                                           c_am(beta,t)*ca_st_lm(alpha,v,u)+c_am(alpha,t)*ca_st_lm(beta,v,u))*&
                                           ohalf*g1_deriv_scale*gaunt1* &
                                           conjg( clm(index1,num,jlo) )*(ALM(index2,num)*judu_20(jlo,l,lp) + &
                                           BLM(index2,num)*judu_21(jlo,l,lp) ) - &
                                           (c_am(beta,t)*ca_st_lm(alpha,v,u)+c_am(alpha,t)*ca_st_lm(beta,v,u))*ohalf*gaunt1*&
                                           conjg( clm(index1,num,jlo) )*( ALM(index2,num)*juu_20(jlo,l,lp) + &
                                           BLM(index2,num)*juu_21(jlo,l,lp) )                                           
                                         ENDIF
                                      ENDDO
                                      
                                      DO jlop = 1, ilo(lp)
                                         IF(loor(jlop,lp)) THEN
                                            F3_buf(alpha,beta,num) = F3_buf(alpha,beta,num) + &
                                           (c_am(alpha,t)*ca_st_lm(beta,v,u)+c_am(beta,t)*ca_st_lm(alpha,v,u) +&
                                           c_am(beta,t)*ca_st_lm(alpha,v,u)+c_am(alpha,t)*ca_st_lm(beta,v,u))*&
                                           ohalf*g1_deriv_scale*gaunt1* &
                                           clm(index2,num,jlop)*(conjg(ALM(index1,num))*judu_02(jlop,lp,l) + &
                                           conjg(BLM(index1,num))*judu_12(jlop,lp,l) ) - &
                                           (c_am(beta,t)*ca_st_lm(alpha,v,u)+c_am(alpha,t)*ca_st_lm(beta,v,u))*ohalf*gaunt1*&
                                           clm(index2,num,jlop)*( conjg(ALM(index1,num))*juu_02(jlop,lp,l) + &
                                           conjg(BLM(index1,num))*juu_12(jlop,lp,l) )                                                                                                                       
                                         ENDIF
                                      ENDDO
                                      
                                      DO jlo = 1, ilo(l)
                                         IF(loor(jlo,l)) THEN
                                          DO jlop = 1, ilo(lp)
                                            IF(loor(jlop,lp)) THEN
                                               F3_buf(alpha,beta,num) = F3_buf(alpha,beta,num) + &
                                               (c_am(alpha,t)*ca_st_lm(beta,v,u)+c_am(beta,t)*ca_st_lm(alpha,v,u) +&
                                               c_am(beta,t)*ca_st_lm(alpha,v,u)+c_am(alpha,t)*ca_st_lm(beta,v,u))*&
                                               ohalf*g1_deriv_scale*gaunt1*2.d0* &
                                               conjg( clm(index1,num,jlo) )*clm(index2,num,jlop)*&
                                               judu_22(jlo,jlop,l,lp)- &
                                               (c_am(beta,t)*ca_st_lm(alpha,v,u)+c_am(alpha,t)*ca_st_lm(beta,v,u))*ohalf*gaunt1*&
                                               2.d0*conjg( clm(index1,num,jlo) )*clm(index2,num,jlop)*&
                                               juu_22(jlo,jlop,l,lp) 
                                            ENDIF
                                          ENDDO
                                         ENDIF
                                      ENDDO
                                      
                                      
                                  ENDDO
                               ENDDO                                                                                      !END beta
                            ENDDO                                                                                         !END alpha                            
                         ENDDO                                                                                            !END t
                      ENDDO                                                                                               !END m
                   ENDDO                                                                                                  !END u
                ENDDO                                                                                                     !END mp
             ENDDO                                                                                                        !END l
          ENDDO                                                                                                           !END v
       ENDDO                                                                                                              !END lp
       
       
       DO lp = 0, lmax
          DO v = -1,1,2
             IF((lp+v) .lt. 0) CYCLE
             DO l = 0, lmax
                DO v1 = -1,1,2
                   IF(l .ne. (lp+v+v1)) CYCLE
                   DO mp = -lp,lp
                      index2 = lp*(lp+1)+mp+1
                      DO m = -l,l
                         index1 = l*(l+1)+m+1
!                          IF (index1 .eq. 1 .or. index2 .eq. 1) CYCLE
                         DO u = -1,1                           
                            IF((mp+u) .lt. -(lp+v)) CYCLE
                            IF((mp+u) .gt. (lp+v) ) CYCLE
                            DO t = -1,1                               
                               IF(m .ne. (mp+u+t)) CYCLE
                                                                                             
                               CALL c_factor(lp,mp,ca_st_lm)
                               CALL c_factor(lp+v,mp+u,ca_st_lm11)                                                                                                                            
                               
                               DO alpha = 1,3
                                  DO beta = 1,3                                                                          
                                     DO num = nemin,nemax
                                     
                                        !HERE IS THE FORMULA FROM F6                                                                                                                        
                                        F6_buf(alpha,beta,num) = F6_buf(alpha,beta,num) + &
                                        (ca_st_lm(alpha,v,u)*ca_st_lm11(beta,v1,t)+ca_st_lm(beta,v,u)*ca_st_lm11(alpha,v1,t))*&
                                        ohalf*(ALM(index2,num)*(conjg(ALM(index1,num))*TC3S00(l,lp)+&
                                        +conjg(BLM(index1,num))*TC3S10(l,lp))+BLM(index2,num)*&
                                        (conjg(ALM(index1,num))*TC3S01(l,lp)+conjg(BLM(index1,num))*TC3S11(l,lp)))
                                       
                                       DO jlo = 1, ilo(l)
                                          IF(loor(jlo,l)) THEN
                                            F6_buf(alpha,beta,num) = F6_buf(alpha,beta,num) + &
                                            (ca_st_lm(alpha,v,u)*ca_st_lm11(beta,v1,t)+ca_st_lm(beta,v,u)*ca_st_lm11(alpha,v1,t))*&
                                            ohalf*conjg( clm(index1,num,jlo) )*( ALM(index2,num)*juu_20(jlo,l,lp) + &
                                            BLM(index2,num)*juu_21(jlo,l,lp))
                                          ENDIF
                                       ENDDO
                                       
                                       DO jlop = 1, ilo(lp)
                                          IF(loor(jlop,lp)) THEN
                                            F6_buf(alpha,beta,num) = F6_buf(alpha,beta,num) + &
                                            (ca_st_lm(alpha,v,u)*ca_st_lm11(beta,v1,t)+ca_st_lm(beta,v,u)*ca_st_lm11(alpha,v1,t))*&
                                            ohalf*clm(index2,num,jlop)*( conjg(ALM(index1,num))*juu_02(jlop,lp,l) + &
                                            conjg(BLM(index1,num))*juu_12(jlop,lp,l))
                                          ENDIF
                                       ENDDO
                                       
                                       DO jlo = 1, ilo(l)
                                          IF(loor(jlo,l)) THEN
                                            DO jlop = 1, ilo(lp)
                                               IF(loor(jlop,lp)) THEN
                                                  F6_buf(alpha,beta,num) = F6_buf(alpha,beta,num) + &
                                                  (ca_st_lm(alpha,v,u)*ca_st_lm11(beta,v1,t)+ca_st_lm(beta,v,u)*ca_st_lm11(alpha,v1,t))*&
                                                  ohalf*conjg( clm(index1,num,jlo) )*clm(index2,num,jlop)*juu_22(jlo,jlop,l,lp)
                                               ENDIF
                                            ENDDO
                                          ENDIF
                                       ENDDO
                                       
                                       
                                     ENDDO                                                                                !END num
                                  ENDDO                                                                                   !END beta
                               ENDDO                                                                                      !END alpha
                               
                            ENDDO                                                                                         !END t
                         ENDDO                                                                                            !END u
                      ENDDO                                                                                               !END l
                   ENDDO                                                                                                  !END mp
                ENDDO                                                                                                     !END v1
             ENDDO                                                                                                        !END l
          ENDDO                                                                                                           !END v
       ENDDO                                                                                                              !END lp
              
       !THESE FOMRULAS ARE CALCULATED USING SECTION (5.6) BUT DERIVED IN MY WAY.
       DO alpha = 1, 3
          DO beta = 1, 3
             DO num = nemin, nemax                              
               F_buf_total(alpha,beta,num) = F1_buf(alpha,beta,num) +F2_buf(alpha,beta,num)+&
                                              F6_buf(alpha,beta,num)+F3_buf(alpha,beta,num)                
                                              
                                             !For the calculations of Al, F3_buf is zero               
             ENDDO
          ENDDO
       ENDDO
  !  7777 continue   
 !=============================END valence kinetice stress inside atomic sphere==============================        
 
! These below lines are for the surafce term coming from the kinetic energy term    


  RETURN
END SUBROUTINE val_kinstrmt

subroutine stress_APW_surface(nemin,nemax,jatom,rmt,ALM,BLM,CLM,A_surf2,buf_APW)
   
   use param   
   use atspdt, only: P,DP,PE,DPE
   use lo    , only: loor,ilo,plo,dplo
   use xa    , only: weight
   implicit none
   
   real*8 ,    intent(in)   :: rmt
   integer,    intent(in)   :: nemin, nemax,jatom
   complex*16, intent(in)   :: CLM((lmax2+1)*(lmax2+1),nume,nloat),&
                               ALM((lmax2+1)*(lmax2+1),nume),  &
                               BLM((lmax2+1)*(lmax2+1),nume),  &                               
                               A_surf2((lmax2+1)*(lmax2+1),nume)
                               
   complex*16, intent(out)  :: buf_APW(3,3,nume)   
   
   complex*16          :: c_am(3,-1:1), buf_temp1(3,3,nume),       &
                          buf_temp2(3,3,nume),buf_temp3(3,3,nume), &
                          buf_temp4(3,3,nume), zeroc
                          
   integer             :: l0,lp,s0,v0,t0,tp, alpha, beta, m0,mp,num, &
                          index1,index2,jlo,jlop, index, index3,lmda,lmdap 
                           
   real*8              :: gaunt1,gaunt2,rmt3,rmt2,tw_by_r                        
   integer, external   :: NOTRI
   real*8 , external   :: GAUNT
   
   
   zeroc = (0.d0,0.0d0)
   call c_alpha_m(c_am)
   
   buf_temp1 = zeroc
   buf_temp2 = zeroc
   buf_temp3 = zeroc
   buf_temp4 = zeroc
   
   buf_APW   = zeroc        
   
   rmt3 = rmt*rmt*rmt
   rmt2 = rmt*rmt
   tw_by_r = 2.d0/rmt
   
   index1 = 0
   index2 = 0
   DO s0 = 0,2
      IF(NOTRI(s0,1,1) .lt. 0) CYCLE
      DO v0 = -s0,s0
         DO t0 = -1,1
            DO tp = -1,1
               IF(-v0+t0+tp .ne. 0) CYCLE
               index1 = 0
               DO l0 = 0,lmax2
                  DO m0 = -l0,l0
                     index1 = index1 + 1                     
                     index2 = 0
                     DO lp = 0,lmax2
                        DO mp = -lp,lp
                           index2 = index2+1
                           IF(NOTRI(l0,lp,s0) .lt. 0) CYCLE
                           IF(-m0+mp+v0 .ne. 0) CYCLE
                           gaunt1 = GAUNT(s0,1,1,v0,t0,tp)
                           gaunt2 = GAUNT(l0,lp,s0,m0,mp,v0)
                           
                           DO alpha = 1, 3
                              DO beta = 1,3
                                 DO num = nemin,nemax
                                    buf_APW(alpha,beta,num) = buf_APW(alpha,beta,num) + &
                                    c_am(alpha,t0)*c_am(beta,tp)*gaunt1*gaunt2*rmt2*&
                                    (conjg(ALM(index1,num))*P(l0)+conjg(BLM(index1,num))*PE(l0))*&
                                    (ALM(index2,num)*DP(lp)+BLM(index2,num)*DPE(lp)-A_surf2(index2,num))

                                    ! INDEX1 is the bra (l0,m0), whereas
                                    ! INDEX2 is the ket (lp,mp).  The rank-two
                                    ! strain tensor couples l0 to l0 and
                                    ! l0+/-2, so the ket radial derivative
                                    ! must be evaluated in its own lp channel.
                                    ! Using DP(l0)/DPE(l0) here was the old
                                    ! WIEN2k stress bug and mixes the bra
                                    ! radial function into off-diagonal
                                    ! angular blocks.
                                 ENDDO
                                 
                                 DO num = nemin,nemax
                                    DO jlo = 1,ilo(l0)
                                       IF(loor(jlo,l0)) THEN
                                         buf_APW(alpha,beta,num) = buf_APW(alpha,beta,num) + &
                                         c_am(alpha,t0)*c_am(beta,tp)*gaunt1*gaunt2*rmt2*&
                                         conjg(clm(index1,num,jlo))*plo(jlo,l0)*   &
                                         (ALM(index2,num)*DP(lp)+BLM(index2,num)*DPE(lp)-A_surf2(index2,num))                                                                                  
                                       ENDIF
                                    ENDDO
                                 ENDDO
                                 
                                 DO num = nemin,nemax
                                    DO jlop = 1,ilo(lp)
                                       IF(loor(jlop,lp)) THEN
                                          buf_APW(alpha,beta,num) = buf_APW(alpha,beta,num) + &
                                          c_am(alpha,t0)*c_am(beta,tp)*gaunt1*gaunt2*rmt2*&
                                          clm(index2,num,jlop)*dplo(jlop,lp)*&
                                          (conjg(ALM(index1,num))*P(l0)+conjg(BLM(index1,num))*PE(l0))
                                       ENDIF
                                    ENDDO
                                 ENDDO
                                 
                                 DO num = nemin,nemax
                                    DO jlo = 1,ilo(l0)
                                       IF(loor(jlo,l0))THEN
                                          DO jlop = 1,ilo(lp)
                                             IF(loor(jlop,lp)) THEN
                                                buf_APW(alpha,beta,num) = buf_APW(alpha,beta,num) + &
                                                c_am(alpha,t0)*c_am(beta,tp)*gaunt1*gaunt2*rmt2*&
                                                conjg(clm(index1,num,jlo))*plo(jlo,l0)*clm(index2,num,jlop)*dplo(jlop,lp)
                                             ENDIF
                                          ENDDO
                                       ENDIF
                                    ENDDO
                                 ENDDO
                                 
                              ENDDO
                           ENDDO                                                                                 
                           
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO

!    !---from inside the sphere
!    index1 = 0
!    DO l0 = 0,lmax2      
!       DO m0 = -l0,l0
!          index1 = index1 + 1
!          DO num = nemin,nemax
!             ex_surf_buf(num) = ex_surf_buf(num) + &              
!               (conjg(ALM(index1,num))*P(l0)+conjg(BLM(index1,num))*PE(l0) )* &
!                   ( ALM(index1,num)*DP(l0)+BLM(index1,num)*DPE(l0) )-    &                                                                               
!                   (conjg(ALM(index1,num))*P(l0)+conjg(BLM(index1,num))*PE(l0) )*A_surf2(index1,num)                  !             
!          ENDDO
!          
!          DO num = nemin,nemax
!             DO jlo = 1, ilo(l0)
!                IF(loor(jlo,l0)) THEN
!                      ex_surf_buf(num) =  ex_surf_buf(num) + &                   
!                    conjg(CLM(index1,num,jlo))*plo(jlo,l0)*&
!                    ( ALM(index1,num)*DP(l0)+BLM(index1,num)*DPE(l0) )+&
!                     CLM(index1,num,jlo)*dplo(jlo,l0)*&
!                     (conjg(ALM(index1,num))*P(l0)+conjg(BLM(index1,num))*PE(l0) )+&
!                     conjg(CLM(index1,num,jlo))*plo(jlo,l0)*CLM(index1,num,jlo)*dplo(jlo,l0)-&
!                    
!                       conjg(CLM(index1,num,jlo))*plo(jlo,l0)*A_surf2(index1,num)                       
!                ENDIF
!             ENDDO
!          ENDDO
!       ENDDO      
!    ENDDO
   
   
   

   
return
end subroutine


