SUBROUTINE val_kinstrmt(nemin,nemax,alm,blm,clm,F_buf_total,F_buf_parts,RMT)
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
    COMPLEX*16, intent(out)    :: F_buf_total(3,3,nume), &
                                  F_buf_parts(3,3,nume,4)
    integer                    :: l,ll,lp,llp,m,mp,t,tp,u,s,sp,v1,v2,v, &
                                   jlo,jlop
    integer                    :: alpha,beta, num, lmda,lmdap
    integer                    :: notri1,notri2, index1, index2, index
    integer, external          :: NOTRI
    REAL*8,  external          :: GAUNT
    real*8                     :: gaunt1,gaunt2,gaunt3,RMT2
    REAL*8, PARAMETER          :: ohalf = 1.D0/2.D0    
    real*8                     :: g1_deriv_scale
    integer                    :: g1_legacy_status
    character(len=32)          :: g1_legacy_env
    complex*16                 :: c_am(3,-1:1), ca_st_lm(3,-1:1,-1:1), &
                                  ca_st_lm11(3,-1:1,-1:1)
    complex*16                 :: angular_buf(3,3),angular_buf2(3,3), &
                                  radial_buf1,radial_buf2
    complex*16                 :: F1_buf(3,3,nume), F2_buf(3,3,nume), &
                                  F3_buf(3,3,nume), F4_buf(3,3,nume), &
                                  F5_buf(3,3,nume), F6_buf(3,3,nume), &                                  
                                  sum_buf(4,nume),buf_tempo , &
                                  F6_buf1(3,3,nume),F6_buf2(3,3,nume), &
                                  F6_buf3(3,3,nume)

                                  
!=====================    
    lmax = lmax2n        !at least lmax = 7 is needed
!=====================                                      

    ! Eq. (G1) has one overall factor 1/2 in front of all four
    ! muffin-tin terms.  In the historical F3 implementation the
    ! first-radial-derivative piece duplicated its angular factor and
    ! therefore lost this 1/2, whereas the accompanying -u/r**2 piece kept
    ! it.  Restoring the common 1/2 removes the full traceless discrepancy
    ! in the exact fixed-G FeSe finite-difference test.  The legacy switch is
    ! retained only to reproduce and diagnose historical WIEN2k results.
    g1_deriv_scale = 0.5d0
    g1_legacy_env = ''
    call get_environment_variable('WIEN_STRESS_G1_LEGACY',g1_legacy_env, &
                                  status=g1_legacy_status)
    if(g1_legacy_status.eq.0 .and. len_trim(g1_legacy_env).gt.0 .and. &
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
                               ! The angular 3x3 tensor is independent of
                               ! band number and of all radial/LO products.
                               ! Historically the complete, expensive radial
                               ! contraction below was repeated for every
                               ! alpha,beta pair.  Factor it once and form the
                               ! outer tensor product afterwards.  This is an
                               ! exact rearrangement of Eq. (G1), not an
                               ! approximation or a diagnostic shortcut.
                               angular_buf = zeroc
                               DO alpha = 1,3
                                  DO beta = 1,3
                                     angular_buf(alpha,beta) = &
                                          (c_am(alpha,t)*c_am(beta,u)+ &
                                           c_am(beta,t)*c_am(alpha,u))* &
                                          ohalf*gaunt1*gaunt2
                                  ENDDO                                                                                   !END beta
                               ENDDO                                                                                      !END alpha
                               DO num = nemin,nemax
                                  radial_buf1 = &
                                       ALM(index2,num)*(conjg(ALM(index1,num))*TC1S00(l,lp)+ &
                                       conjg(BLM(index1,num))*TC1S10(l,lp))+BLM(index2,num)* &
                                       (conjg(ALM(index1,num))*TC1S01(l,lp)+ &
                                       conjg(BLM(index1,num))*TC1S11(l,lp))
                                  DO jlo = 1,ilo(l)
                                     IF(loor(jlo,l)) radial_buf1 = radial_buf1 + &
                                          conjg(clm(index1,num,jlo))* &
                                          (ALM(index2,num)*juddu_20(jlo,l,lp)+ &
                                           BLM(index2,num)*juddu_21(jlo,l,lp))
                                  ENDDO
                                  DO jlop = 1,ilo(lp)
                                     IF(loor(jlop,lp)) radial_buf1 = radial_buf1 + &
                                          clm(index2,num,jlop)* &
                                          (conjg(ALM(index1,num))*juddu_02(jlop,lp,l)+ &
                                           conjg(BLM(index1,num))*juddu_12(jlop,lp,l))
                                  ENDDO
                                  DO jlo = 1,ilo(l)
                                     IF(.NOT.loor(jlo,l)) CYCLE
                                     DO jlop = 1,ilo(lp)
                                        IF(loor(jlop,lp)) radial_buf1 = radial_buf1 + &
                                             conjg(clm(index1,num,jlo))*clm(index2,num,jlop)* &
                                             juddu_22(jlo,jlop,l,lp)
                                     ENDDO
                                  ENDDO
                                  F1_buf(:,:,num)=F1_buf(:,:,num)+angular_buf(:,:)*radial_buf1
                               ENDDO
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
                            angular_buf = zeroc
                            DO alpha = 1,3
                               DO beta = 1,3
                                  angular_buf(alpha,beta) = &
                                       (c_am(alpha,t)*ca_st_lm(beta,v,u)+ &
                                        c_am(beta,t)*ca_st_lm(alpha,v,u))*ohalf*gaunt1
                               ENDDO                                                                                      !END beta
                            ENDDO                                                                                         !END alpha                            
                            DO num = nemin,nemax
                               radial_buf1 = &
                                    ALM(index2,num)*(conjg(ALM(index1,num))*TC2S00(l,lp)+ &
                                    conjg(BLM(index1,num))*TC2S10(l,lp))+BLM(index2,num)* &
                                    (conjg(ALM(index1,num))*TC2S01(l,lp)+ &
                                    conjg(BLM(index1,num))*TC2S11(l,lp))
                               DO jlo = 1,ilo(l)
                                  IF(loor(jlo,l)) radial_buf1 = radial_buf1 + &
                                       conjg(clm(index1,num,jlo))* &
                                       (ALM(index2,num)*judu_20(jlo,l,lp)+ &
                                        BLM(index2,num)*judu_21(jlo,l,lp))
                               ENDDO
                               DO jlop = 1,ilo(lp)
                                  IF(loor(jlop,lp)) radial_buf1 = radial_buf1 + &
                                       clm(index2,num,jlop)* &
                                       (conjg(ALM(index1,num))*judu_02(jlop,lp,l)+ &
                                        conjg(BLM(index1,num))*judu_12(jlop,lp,l))
                               ENDDO
                               DO jlo = 1,ilo(l)
                                  IF(.NOT.loor(jlo,l)) CYCLE
                                  DO jlop = 1,ilo(lp)
                                     IF(loor(jlop,lp)) radial_buf1 = radial_buf1 + &
                                          conjg(clm(index1,num,jlo))*clm(index2,num,jlop)* &
                                          judu_22(jlo,jlop,l,lp)
                                  ENDDO
                               ENDDO
                               F2_buf(:,:,num)=F2_buf(:,:,num)+angular_buf(:,:)*radial_buf1
                            ENDDO
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
                            angular_buf = zeroc
                            angular_buf2 = zeroc
                            DO alpha =1,3
                               DO beta = 1,3
                                  angular_buf(alpha,beta) = &
                                       (c_am(alpha,t)*ca_st_lm(beta,v,u)+ &
                                        c_am(beta,t)*ca_st_lm(alpha,v,u)+ &
                                        c_am(beta,t)*ca_st_lm(alpha,v,u)+ &
                                        c_am(alpha,t)*ca_st_lm(beta,v,u))* &
                                       ohalf*g1_deriv_scale*gaunt1
                                  angular_buf2(alpha,beta) = &
                                       (c_am(beta,t)*ca_st_lm(alpha,v,u)+ &
                                        c_am(alpha,t)*ca_st_lm(beta,v,u))*ohalf*gaunt1
                               ENDDO                                                                                      !END beta
                            ENDDO                                                                                         !END alpha                            
                            DO num = nemin,nemax
                               radial_buf1 = &
                                    ALM(index2,num)*(conjg(ALM(index1,num))*TC2S00(l,lp)+ &
                                    conjg(BLM(index1,num))*TC2S10(l,lp))+BLM(index2,num)* &
                                    (conjg(ALM(index1,num))*TC2S01(l,lp)+ &
                                    conjg(BLM(index1,num))*TC2S11(l,lp))
                               radial_buf2 = &
                                    ALM(index2,num)*(conjg(ALM(index1,num))*TC3S00(l,lp)+ &
                                    conjg(BLM(index1,num))*TC3S10(l,lp))+BLM(index2,num)* &
                                    (conjg(ALM(index1,num))*TC3S01(l,lp)+ &
                                    conjg(BLM(index1,num))*TC3S11(l,lp))
                               DO jlo = 1,ilo(l)
                                  IF(.NOT.loor(jlo,l)) CYCLE
                                  radial_buf1 = radial_buf1 + conjg(clm(index1,num,jlo))* &
                                       (ALM(index2,num)*judu_20(jlo,l,lp)+ &
                                        BLM(index2,num)*judu_21(jlo,l,lp))
                                  radial_buf2 = radial_buf2 + conjg(clm(index1,num,jlo))* &
                                       (ALM(index2,num)*juu_20(jlo,l,lp)+ &
                                        BLM(index2,num)*juu_21(jlo,l,lp))
                               ENDDO
                               DO jlop = 1,ilo(lp)
                                  IF(.NOT.loor(jlop,lp)) CYCLE
                                  radial_buf1 = radial_buf1 + clm(index2,num,jlop)* &
                                       (conjg(ALM(index1,num))*judu_02(jlop,lp,l)+ &
                                        conjg(BLM(index1,num))*judu_12(jlop,lp,l))
                                  radial_buf2 = radial_buf2 + clm(index2,num,jlop)* &
                                       (conjg(ALM(index1,num))*juu_02(jlop,lp,l)+ &
                                        conjg(BLM(index1,num))*juu_12(jlop,lp,l))
                               ENDDO
                               DO jlo = 1,ilo(l)
                                  IF(.NOT.loor(jlo,l)) CYCLE
                                  DO jlop = 1,ilo(lp)
                                     IF(.NOT.loor(jlop,lp)) CYCLE
                                     radial_buf1 = radial_buf1 + 2.d0* &
                                          conjg(clm(index1,num,jlo))*clm(index2,num,jlop)* &
                                          judu_22(jlo,jlop,l,lp)
                                     radial_buf2 = radial_buf2 + 2.d0* &
                                          conjg(clm(index1,num,jlo))*clm(index2,num,jlop)* &
                                          juu_22(jlo,jlop,l,lp)
                                  ENDDO
                               ENDDO
                               F3_buf(:,:,num)=F3_buf(:,:,num)+ &
                                    angular_buf(:,:)*radial_buf1- &
                                    angular_buf2(:,:)*radial_buf2
                            ENDDO
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
                               angular_buf = zeroc
                               DO alpha = 1,3
                                  DO beta = 1,3                                                                          
                                     angular_buf(alpha,beta) = &
                                          (ca_st_lm(alpha,v,u)*ca_st_lm11(beta,v1,t)+ &
                                           ca_st_lm(beta,v,u)*ca_st_lm11(alpha,v1,t))*ohalf
                                  ENDDO                                                                                   !END beta
                               ENDDO                                                                                      !END alpha
                               DO num = nemin,nemax
                                  radial_buf1 = &
                                       ALM(index2,num)*(conjg(ALM(index1,num))*TC3S00(l,lp)+ &
                                       conjg(BLM(index1,num))*TC3S10(l,lp))+BLM(index2,num)* &
                                       (conjg(ALM(index1,num))*TC3S01(l,lp)+ &
                                       conjg(BLM(index1,num))*TC3S11(l,lp))
                                  DO jlo = 1,ilo(l)
                                     IF(loor(jlo,l)) radial_buf1 = radial_buf1 + &
                                          conjg(clm(index1,num,jlo))* &
                                          (ALM(index2,num)*juu_20(jlo,l,lp)+ &
                                           BLM(index2,num)*juu_21(jlo,l,lp))
                                  ENDDO
                                  DO jlop = 1,ilo(lp)
                                     IF(loor(jlop,lp)) radial_buf1 = radial_buf1 + &
                                          clm(index2,num,jlop)* &
                                          (conjg(ALM(index1,num))*juu_02(jlop,lp,l)+ &
                                           conjg(BLM(index1,num))*juu_12(jlop,lp,l))
                                  ENDDO
                                  DO jlo = 1,ilo(l)
                                     IF(.NOT.loor(jlo,l)) CYCLE
                                     DO jlop = 1,ilo(lp)
                                        IF(loor(jlop,lp)) radial_buf1 = radial_buf1 + &
                                             conjg(clm(index1,num,jlo))*clm(index2,num,jlop)* &
                                             juu_22(jlo,jlop,l,lp)
                                     ENDDO
                                  ENDDO
                                  F6_buf(:,:,num)=F6_buf(:,:,num)+angular_buf(:,:)*radial_buf1
                               ENDDO
                               
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
               ! Keep the four lines of Blaha et al. Eq. (G1) separately for
               ! diagnostics.  The production result below is still their
               ! unchanged sum.  The part ordering is:
               !   1 = F1 (second radial derivative),
               !   2 = F2 (first angular-gradient term),
               !   3 = F3 (combined third line, including F4-F5 notation),
               !   4 = F6 (two angular-gradient operators).
               F_buf_parts(alpha,beta,num,1) = F1_buf(alpha,beta,num)
               F_buf_parts(alpha,beta,num,2) = F2_buf(alpha,beta,num)
               F_buf_parts(alpha,beta,num,3) = F3_buf(alpha,beta,num)
               F_buf_parts(alpha,beta,num,4) = F6_buf(alpha,beta,num)
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
!CKAMAL STRESS/PRESS/NEW SURFACE CORREECTION 9 MAY 2018


! SUBROUTINE rela_stress_corr(nemin,nemax,jatom,alm,blm,clm,buf_rela_corr)
!   use struk  ,only: r0,rmt,dx,jri
!   use param
!   use stress_kb, only : RI_00,RI_01,RI_10,RI_11,     &
!                         RI2_00,RI2_01,RI2_10,RI2_11, &
!                         RI3_00,RI3_01,RI3_10,RI3_11, &   !this comes from the derivative of M
!                         RI_LO_02,RI_LO_12,RI_LO_20 , &
!                         RI_LO_21, RI_LO_22,         &
!                         RI2_LO_02,RI2_LO_12,RI2_LO_20 , &
!                         RI2_LO_21, RI2_LO_22     , &
!                         RI3_LO_02,RI3_LO_12,RI3_LO_20 , &
!                         RI3_LO_21, RI3_LO_22     , &
!                          rela_int_00, rela_int_01, rela_int_02, &
!                          rela_int_10, rela_int_11, rela_int_12, &
!                          rela_int_20, rela_int_21, rela_int_22, &
!                          K_int_00, K_int_01, K_int_10, K_int_11, &
!                          K_int_02, K_int_20, K_int_12, K_int_21, &
!                          K_int_22
!   use defs     , only : zeroc
!   use lo       , only : ilo,loor
!   implicit none 
!   COMPLEX*16, intent(in)     :: alm((lmax2+1)*(lmax2+1),nume), &
!                                 blm((lmax2+1)*(lmax2+1),nume), &
!                                 clm((lmax2+1)*(lmax2+1),nume,nloat)
!   integer   , intent(in)     :: nemin,nemax, jatom                                
!   COMPLEX*16, intent(out)    :: buf_rela_corr(3,3,nume)                                
!   complex*16                 :: c_am(3,-1:1), REL_BUF(3,3,nume), rela_test(3,3,nume)
!   real*8                     :: gaunt1,gaunt2,ohalf  
!   integer                    :: l0,l0p,m0,m0p,s0,v0,t0,t0p,index1,index2, &
!                                 alpha, beta, num, jlo, jlop
!                                 
!   integer, external          :: NOTRI
!   real*8 , external          :: GAUNT
!   
!   
!   REL_BUF = zeroc
!   rela_test = zeroc
!   ohalf  = 0.5d0
!   
!   call c_alpha_m(c_am)    
!   
!   DO s0 = 0,2,2
!      DO v0 = -s0,s0
!         DO t0 = -1,1
!            DO t0p = -1,1
!               IF((-v0+t0+t0p) .ne. 0) CYCLE
!               index1 = 0
!               DO l0 = 0,lmax2
!                  DO m0 = -l0,l0
!                     index1 = index1 + 1
!                     index2 = 0 
!                     DO l0p = 0,lmax2
!                        DO m0p = -l0p,l0p
!                           index2 = index2 + 1
!                           
! !                           IF( index1 .ne. 1 .or. index2 .ne. 1) CYCLE
!                           
!                           IF(NOTRI(l0,l0p,s0) .lt. 0) CYCLE
!                           IF((-m0+m0p+v0) .ne. 0) CYCLE
!                           gaunt1 = GAUNT(s0,1,1,v0,t0,t0p)
!                           gaunt2 = GAUNT(l0,l0p,s0,m0,m0p,v0)
!                           DO alpha = 1, 3
!                              DO beta = 1,3
!                                 DO num = nemin, nemax
!                                    REL_BUF(alpha,beta,num) = REL_BUF(alpha,beta,num) - &
!                                       gaunt1*gaunt2*c_am(alpha,t0)*c_am(beta,t0p)*(&
! !                                     conjg(ALM(index1,num))*ALM(index2,num)*(RI_00(l0,l0p)+ RI2_00(l0,l0p) )+ &                                      
! !                                     conjg(BLM(index1,num))*ALM(index2,num)*(RI_10(l0,l0p)+ RI2_10(l0,l0p) )+ &                                      
! !                                     conjg(ALM(index1,num))*BLM(index2,num)*(RI_01(l0,l0p)+ RI2_01(l0,l0p) )+ &                                      
! !                                     conjg(BLM(index1,num))*BLM(index2,num)*(RI_11(l0,l0p)+ RI2_11(l0,l0p) ) ) ! + & 
! !                                       
!                                       conjg(ALM(index1,num))*ALM(index2,num)*(RI_00(l0,l0p) ) + &
!                                       conjg(BLM(index1,num))*ALM(index2,num)*(RI_10(l0,l0p) )+ &
!                                       conjg(ALM(index1,num))*BLM(index2,num)*(RI_01(l0,l0p) )+ &
!                                       conjg(BLM(index1,num))*BLM(index2,num)*(RI_11(l0,l0p) )      + &
! ! 
!                                       conjg(ALM(index1,num))*ALM(index2,num)*(RI2_00(l0,l0p) ) + &
!                                       conjg(BLM(index1,num))*ALM(index2,num)*(RI2_10(l0,l0p) )+ &
!                                       conjg(ALM(index1,num))*BLM(index2,num)*(RI2_01(l0,l0p) )+ &
!                                       conjg(BLM(index1,num))*BLM(index2,num)*(RI2_11(l0,l0p) )     + &
! 
!                                       conjg(ALM(index1,num))*ALM(index2,num)*(RI3_00(l0,l0p) ) + &
!                                       conjg(BLM(index1,num))*ALM(index2,num)*(RI3_10(l0,l0p) )+ &
!                                       conjg(ALM(index1,num))*BLM(index2,num)*(RI3_01(l0,l0p) )+ &
!                                       conjg(BLM(index1,num))*BLM(index2,num)*(RI3_11(l0,l0p) )     + &
! ! ! 
!                                       conjg(ALM(index1,num))*ALM(index2,num)*(K_int_00(l0,l0p) ) + &
!                                       conjg(BLM(index1,num))*ALM(index2,num)*(K_int_10(l0,l0p) )+ &
!                                       conjg(ALM(index1,num))*BLM(index2,num)*(K_int_01(l0,l0p) )+ &
!                                       conjg(BLM(index1,num))*BLM(index2,num)*(K_int_11(l0,l0p) )   )
!                                       
!                                    DO jlo = 1, ilo(l0)
!                                       IF( loor(jlo,l0) ) THEN
!                                         REL_BUF(alpha,beta,num) = REL_BUF(alpha,beta,num) - &
!                                                 gaunt1*gaunt2*c_am(alpha,t0)*c_am(beta,t0p)*conjg(clm(index1,num,jlo))*( &
! !                                               ALM(index2,num)*( RI_LO_20(jlo,l0,l0p)+RI2_LO_20(jlo,l0,l0p) ) + &
! !                                               BLM(index2,num)*( RI_LO_21(jlo,l0,l0p)+RI2_LO_21(jlo,l0,l0p) ) ) 
!                                                 
!                                                 ALM(index2,num)*( RI_LO_20(jlo,l0,l0p)  ) + &
!                                                 BLM(index2,num)*( RI_LO_21(jlo,l0,l0p)  )      + &                                      
! !                                                 
!                                                 ALM(index2,num)*( RI2_LO_20(jlo,l0,l0p)  ) + &
!                                                 BLM(index2,num)*( RI2_LO_21(jlo,l0,l0p)  )   + &                                      
!                                                 
!                                                 ALM(index2,num)*( RI3_LO_20(jlo,l0,l0p)  ) + &
!                                                 BLM(index2,num)*( RI3_LO_21(jlo,l0,l0p)  )      + &
! ! 
!                                                 ALM(index2,num)*( K_int_20(jlo,l0,l0p)  ) + &
!                                                 BLM(index2,num)*( K_int_21(jlo,l0,l0p)  )  )
!                                       ENDIF
!                                    ENDDO
!                                       
!                                    DO jlop = 1, ilo(l0p)
!                                       IF( loor(jlop,l0p) ) THEN
!                                         REL_BUF(alpha,beta,num) = REL_BUF(alpha,beta,num) - &
!                                                 gaunt1*gaunt2*c_am(alpha,t0)*c_am(beta,t0p)*clm(index2,num,jlop)*(&  
! !                                            conjg(ALM(index1,num))*( RI_LO_02(jlop,l0p,l0)+RI2_LO_02(jlop,l0p,l0) ) + &
! !                                           conjg(BLM(index1,num))*( RI_LO_12(jlop,l0p,l0)+RI2_LO_12(jlop,l0p,l0) )  ) 
!                                                 
!                                                 conjg(ALM(index1,num))*( RI_LO_02(jlop,l0p,l0)  ) + &
!                                                 conjg(BLM(index1,num))*( RI_LO_12(jlop,l0p,l0)  )     + &
! !                                                 
!                                                 conjg(ALM(index1,num))*( RI2_LO_02(jlop,l0p,l0)  ) + &
!                                                 conjg(BLM(index1,num))*( RI2_LO_12(jlop,l0p,l0)  )   + &  
!                                                 
!                                                 conjg(ALM(index1,num))*( RI3_LO_02(jlop,l0p,l0)  ) + &
!                                                 conjg(BLM(index1,num))*( RI3_LO_12(jlop,l0p,l0)  )   + &
! 
!                                                 conjg(ALM(index1,num))*( K_int_02(jlop,l0p,l0)  ) + &
!                                                 conjg(BLM(index1,num))*( K_int_12(jlop,l0p,l0)  )   )
!                                       ENDIF
!                                    ENDDO                                      
! !                                    
!                                    DO jlo = 1, ilo(l0)
!                                       IF( loor(jlo,l0) ) THEN
!                                          DO jlop = 1, ilo(l0p)
!                                             IF( loor(jlop,l0p) ) THEN
!                                                REL_BUF(alpha,beta,num) = REL_BUF(alpha,beta,num) - &
!                                                               gaunt1*gaunt2*c_am(alpha,t0)*c_am(beta,t0p)* (&
! !                                                             conjg(clm(index1,num,jlo))*clm(index2,num,jlop)*&
! !                                                             ( RI_LO_22(jlo,jlop,l0,l0p)+RI2_LO_22(jlo,jlop,l0,l0p) ) ) !+ &
! !                                                             
!                                                               conjg(clm(index1,num,jlo))*clm(index2,num,jlop)*( RI_LO_22(jlo,jlop,l0,l0p) )  + &
!                                                               conjg(clm(index1,num,jlo))*clm(index2,num,jlop)*( RI2_LO_22(jlo,jlop,l0,l0p) )  + &
!                                                               conjg(clm(index1,num,jlo))*clm(index2,num,jlop)*( RI3_LO_22(jlo,jlop,l0,l0p) )  + &
!                                                               conjg(clm(index1,num,jlo))*clm(index2,num,jlop)*( K_int_22(jlo,jlop,l0,l0p) )  )
!                                             ENDIF
!                                          ENDDO
!                                       ENDIF
!                                    ENDDO
!                                       
!                                 ENDDO
!                              ENDDO
!                           ENDDO
!                           
!                        ENDDO
!                     ENDDO
!                  ENDDO
!               ENDDO
!            ENDDO
!         ENDDO
!      ENDDO
!   ENDDO  
!   buf_rela_corr = REL_BUF*ohalf      
!   
! RETURN
! END SUBROUTINE rela_stress_corr
