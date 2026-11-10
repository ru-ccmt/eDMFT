! SUBROUTINE val_corr_stress(isize,ibpp_max,ibb,nemin,nemax,ALM,BLM,clm,CAPD0LM,CAPD1LM,buf_sph_val_fin, jatom, rmt,latom, &
!                             lfirst)
SUBROUTINE val_corr_stress(isize,ibpp_max,ibb,nemin,nemax,ALM,BLM,clm, &
                           buf_sph_val_fin,buf_sph_val_parts, &
                           buf_nsp_match,buf_nsp_radial, &
                           buf_nsp_angular,buf_nsp_radial_parts, &
                           buf_nsp_radial_direct,buf_nsp_radial_phase, &
                           buf_nsp_v_direct,buf_nsp_v_parts, &
                           jatom,rmt,latom,lfirst, &
                           ihmx,tuu,tdd,tud,tdu,wuu,wdd,wud,wdu, &
                           vpuu,vpdd,vpud,vpdu,wpuu,wpdd,wpud,wpdu, &
                           tuu12,tuu21,tuu22,tud21,tdu12, &
                           wuu12,wuu21,wuu22,wud21,wdu12, &
                           vpuu12,vpuu21,vpuu22,vpud21,vpdu12, &
                           wpuu12,wpuu21,wpuu22,wpud21,wpdu12, &
                           lv,lpv,llv,mv,mpv,mmv,g2_diag_enabled)
    use param
    use xa3   ,    only: a,bkx, bky, bkz
                         
    use struk ,    only: vol, rotloc, rotij, br1, pos, tauij
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
    integer, intent(in)     :: nemin,nemax,isize,ibpp_max,ibb,jatom,latom,lfirst,ihmx
    logical, intent(in)     :: g2_diag_enabled
    complex*16, intent(in)  :: ALM((lmax2+1)*(lmax2+1),nume), &
                               BLM((lmax2+1)*(lmax2+1),nume)
    complex*16, intent(in)  :: clm((lmax2+1)*(lmax2+1),nume,nloat)                                   
    
    complex*16, intent(out) :: buf_sph_val_fin(1:3,1:3,nume), &
                               buf_sph_val_parts(1:3,1:3,nume,3), &
                               buf_nsp_match(1:3,1:3,nume), &
                               buf_nsp_radial(1:3,1:3,nume), &
                               buf_nsp_angular(1:3,1:3,nume), &
                               buf_nsp_radial_parts(1:3,1:3,nume,4), &
                               buf_nsp_radial_direct(nume), &
                               buf_nsp_radial_phase(nume), &
                               buf_nsp_v_direct(nume), &
                               buf_nsp_v_parts(nume,4)
    real*8    , intent(in)  :: rmt
    integer   , intent(in)  :: lv(ngau),lpv(ngau),llv(ngau), &
                               mv(ngau),mpv(ngau),mmv(ngau)
    complex*16, intent(in)  :: tuu(ngau),tdd(ngau),tud(ngau),tdu(ngau), &
                               wuu(ngau),wdd(ngau),wud(ngau),wdu(ngau), &
                               vpuu(ngau),vpdd(ngau),vpud(ngau),vpdu(ngau), &
                               wpuu(ngau),wpdd(ngau),wpud(ngau),wpdu(ngau), &
                               tuu12(nloat,ngau),tuu21(nloat,ngau), &
                               tuu22(nloat,nloat,ngau), &
                               tud21(nloat,ngau),tdu12(nloat,ngau), &
                               wuu12(nloat,ngau),wuu21(nloat,ngau), &
                               wuu22(nloat,nloat,ngau), &
                               wud21(nloat,ngau),wdu12(nloat,ngau), &
                               vpuu12(nloat,ngau),vpuu21(nloat,ngau), &
                               vpuu22(nloat,nloat,ngau), &
                               vpud21(nloat,ngau),vpdu12(nloat,ngau), &
                               wpuu12(nloat,ngau),wpuu21(nloat,ngau), &
                               wpuu22(nloat,nloat,ngau), &
                               wpud21(nloat,ngau),wpdu12(nloat,ngau)
    
    complex*16              :: CAPD0LM(3,3,(lmax2+1)*(lmax2+1),nume), &
                               CAPD1LM(3,3,(lmax2+1)*(lmax2+1),nume)
    complex*16              :: buf_rela_pulay_d(3,3,nume)
    complex*16 before_lo,afac_nsp,bfac_nsp,vmat_nsp,wmat_nsp, &
               wmat_parts(4), &
               angular_factor,radial_factor, &
               radial_geom(3,3),angular_geom(3,3)
    integer                 :: delta(0:1,0:1), lmd, lmdp, ii
    integer                 :: sp,t,tp,u,v,v1, l, lp, m, mp , index, index1, &
                               index2, alpha, beta, num, i3, index3, index4, &
                               ig, igp, jlo, jlop, maxx,ih,ly,lpy,klo, &
                               lint,nu,mbra,mket
                               
    real*8                  :: gaunt1, gaunt2,gaunt_orig,gaunt_pot, &
                               gaunt_couple,gaunt_ratio, &
                               channel_diff,channel_scale, &
                               g2_d_scale, &
                               jlpl(0:1,0:1,0:lmax2,0:lmax2), &
                               Ilpl(0:1,0:1,0:lmax2,0:lmax2)                               
    logical                 :: duplicate_radial_channel, d_fd_check_enabled, &
                               d_component_energy,d_component_cross, &
                               d_component_lo
    integer                 :: g2_d_scale_status, d_fd_check_status, &
                               d_component_status
    character(len=32)       :: g2_d_scale_env, d_fd_check_env, &
                               d_component_env,rela_pulay_d_env
    integer                 :: rela_pulay_d_status
    logical                 :: rela_pulay_d_apply
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
                               h_al(isize),h_bl(isize), rmt2p, &
                               fd_step, fd_k0, fd_kplus, fd_kminus, &
                               fd_fjplus(0:lmax2),fd_dfjplus(0:lmax2), &
                               fd_fjminus(0:lmax2),fd_dfjminus(0:lmax2), &
                               fd_an_norm0,fd_an_norm1, &
                               fd_ref_norm0,fd_ref_norm1, &
                               fd_diff_norm0,fd_diff_norm1, &
                               fd_cross0,fd_cross1, &
                               fd_scale0,fd_scale1,fd_relerr0,fd_relerr1
                                       
    complex*16              :: pcfac,pcfac1,pcfac2
    complex*16              :: h_yl1((lmax2+2)*(lmax2+2),isize), &                                                              
                               dl0ab_buf((lmax2+1)*(lmax2+1), isize), &
                               dl1ab_buf((lmax2+1)*(lmax2+1), isize), &
                               dl0ab_trace((lmax2+1)*(lmax2+1), isize), &
                               dl1ab_trace((lmax2+1)*(lmax2+1), isize), &
                               CAPD0LM_BUF((lmax2+1)*(lmax2+1),nume), &
                               CAPD1LM_BUF((lmax2+1)*(lmax2+1),nume)
    complex*16              :: fd_coeff0_plus,fd_coeff0_minus,fd_deriv0, &
                               fd_coeff1_plus,fd_coeff1_minus,fd_deriv1, &
                               fd_analytic0,fd_analytic1
                                   
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
    buf_sph_val_parts = zeroc
    buf_nsp_match = zeroc
    buf_nsp_radial = zeroc
    buf_nsp_angular = zeroc
    buf_nsp_radial_parts = zeroc
    buf_nsp_radial_direct = zeroc
    buf_nsp_radial_phase = zeroc
    buf_nsp_v_direct = zeroc
    buf_nsp_v_parts = zeroc
    
    buf_sph_val1_buf1 = zeroc
    buf_sph_val1_buf2 = zeroc
    buf_sph_val1_buf3 = zeroc
    buf_rela_pulay_d = zeroc

    ! Scalar-relativistic correction to the first (matching-coefficient)
    ! Pulay line.  The historical identity used below assumes that H is
    ! independent of the radial solution energy.  Koelling-Harmon has
    ! H_SR(epsilon)=-div[A(r,epsilon)grad]+V, so that identity is not valid
    ! when a basis radial function generated at E_l is tested against a
    ! band eigenvalue epsilon.  APPLY evaluates the exact strong radial
    ! matrix element and adds exact-minus-historical.  Keep this opt-in until
    ! the independent finite-strain pressure establishes the sign and size.
    rela_pulay_d_apply=.false.
    rela_pulay_d_env=''
    call get_environment_variable('WIEN_STRESS_RELA_PULAY_D', &
                                  rela_pulay_d_env, &
                                  status=rela_pulay_d_status)
    if(rela_pulay_d_status.eq.0 .and. &
       len_trim(rela_pulay_d_env).gt.0) then
       select case(trim(adjustl(rela_pulay_d_env)))
       case('APPLY','apply','1')
          rela_pulay_d_apply=.true.
       case('NONE','none','0')
          rela_pulay_d_apply=.false.
       case default
          write(6,*) 'Invalid WIEN_STRESS_RELA_PULAY_D: ', &
                     trim(rela_pulay_d_env)
          stop 'RELA_PULAY_D must be NONE or APPLY'
       end select
    endif

    ! Diagnostic scale for the first (matching-coefficient D) line of the
    ! spherical part of Eq. (G2).  An early total-trace comparison happened
    ! to improve when this line was multiplied by 1/4, but the direct
    ! coefficient-level finite difference below subsequently showed that D
    ! itself has unit scale (relative errors about 1.e-10 in the FeSe test).
    ! The quarter factor was therefore an accidental cancellation with a
    ! different, still missing augmentation Pulay contribution.  Keep this
    ! switch only for reproducing that diagnostic; the default and the
    ! published formula both use scale 1.
    g2_d_scale = 1.d0
    g2_d_scale_env = ''
    call get_environment_variable('WIEN_STRESS_G2_D_SCALE', &
                                  g2_d_scale_env, &
                                  status=g2_d_scale_status)
    if(g2_d_scale_status.eq.0 .and. &
       len_trim(g2_d_scale_env).gt.0) then
       read(g2_d_scale_env,*,err=998) g2_d_scale
    endif

    ! Optional, non-invasive verification of the matching-coefficient
    ! derivative D in Eq. (G4).  Under an isotropic strain eps=t*I the
    ! reciprocal vector changes as q(t)=(1-t)q, whereas its direction and
    ! the phase at a fixed fractional atomic coordinate do not change.  We
    ! can therefore finite-difference the *same* matching coefficients used
    ! by LAPW without rerunning lapw1.  The comparison below is between
    !
    !   sum_alpha D^{alpha alpha}_{lm,G}
    !
    ! assembled by the analytic code and a central finite difference of the
    ! Bessel matching formula.  It changes no stress contribution.  The
    ! fitted scale printed in :DBG_STR_D_FD is the number by which the
    ! analytic D would have to be multiplied to best match the numerical D;
    ! it must be one if D itself is correct.  Enable only for debugging with
    ! WIEN_STRESS_D_FD_CHECK=1 because it evaluates extra Bessel functions.
    d_fd_check_enabled = .false.
    d_fd_check_env = ''
    call get_environment_variable('WIEN_STRESS_D_FD_CHECK', &
                                  d_fd_check_env, &
                                  status=d_fd_check_status)
    if(d_fd_check_status.eq.0 .and. &
       len_trim(d_fd_check_env).gt.0 .and. &
       trim(d_fd_check_env).ne.'0') d_fd_check_enabled = .true.

    ! Component selector for diagnosing the first line of Eq. (G2).  The
    ! three algebraically distinct pieces test different possible sources
    ! of a normalization error:
    !
    !   ENERGY : (El-epsilon)[D0^* A + D1^* B <udot|udot>]
    !   CROSS  : D0^* B, from <u|(H-epsilon)|udot>=1
    !   LO     : the local-orbital extension of the ENERGY piece.
    !
    ! ALL is the published expression and is always the default.  Selecting
    ! one component is intentionally diagnostic and changes the computed
    ! stress, but only when WIEN_STRESS_G2_D_COMPONENT is explicitly set.
    d_component_energy = .true.
    d_component_cross = .true.
    d_component_lo = .true.
    d_component_env = ''
    call get_environment_variable('WIEN_STRESS_G2_D_COMPONENT', &
                                  d_component_env, &
                                  status=d_component_status)
    if(d_component_status.eq.0 .and. &
       len_trim(d_component_env).gt.0) then
       d_component_energy = .false.
       d_component_cross = .false.
       d_component_lo = .false.
       select case(trim(adjustl(d_component_env)))
       case('ALL','all')
          d_component_energy = .true.
          d_component_cross = .true.
          d_component_lo = .true.
       case('ENERGY','energy')
          d_component_energy = .true.
       case('CROSS','cross')
          d_component_cross = .true.
       case('LO','lo')
          d_component_lo = .true.
       case default
          write(6,*) 'Invalid WIEN_STRESS_G2_D_COMPONENT: ', &
                     trim(d_component_env)
          stop 'G2 D component must be ALL, ENERGY, CROSS, or LO'
       end select
    endif
        
    
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
          
          CALL ROTATE (bk,rotij(1,1,latom), bkr)
          
          bk(1) = bkr(1)*BR1(1,1) + bkr(2)*BR1(1,2) + bkr(3)*BR1(1,3)
          bk(2) = bkr(1)*BR1(2,1) + bkr(2)*BR1(2,2) + bkr(3)*BR1(2,3)
          bk(3) = bkr(1)*BR1(3,1) + bkr(2)*BR1(3,2) + bkr(3)*BR1(3,3)
          
         CALL ROTATE (bk,rotloc(1,1,jatom), bkr1)
         
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
              dl0ab_trace = zeroc
              dl1ab_trace = zeroc
            
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
                       
                       h_al(i3)  = dfj(l,i3)*pe(l)-fj(l,i3)*dpe(l) 
                       h_bl(i3)  = fj(l,i3)*dp(l)-dfj(l,i3)*p(l)                      
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
                       IF(d_fd_check_enabled .and. alpha.eq.beta) THEN
                          dl0ab_trace(index,i3) = dl0ab_trace(index,i3) + &
                                                 dl0ab_buf(index,i3)
                          dl1ab_trace(index,i3) = dl1ab_trace(index,i3) + &
                                                 dl1ab_buf(index,i3)
                       ENDIF
!                        IF(i3.eq.1.and.index.eq.1.and.abs(dl0ab_buf(index,i3)).gt.1.d0) write(6, '("dl0",6f20.10)') dl0ab_buf(index,i3), h_yll11(index,i3),d_al1(i3)
                     ENDDO                                                                                              !END i3
                     
                    ENDDO                                                                                               !END m
                 ENDDO                                                                                                  !END l 
                      
                      ! D_{lm,nu}^{alpha beta} is the dense contraction
                      !
                      !   sum_G dA_{lm,G}^{alpha beta} A_{G,nu} .
                      !
                      ! This used to be evaluated by the three scalar loops
                      ! over lm, band, and G above.  At realistic RKMAX the
                      ! contraction dominates the LAPW2 stress run, although
                      ! it is precisely a matrix multiplication.  Use BLAS
                      ! here, as the older blocked implementation in L2MAIN
                      ! was designed to do.  This changes only the summation
                      ! order; the mathematical expression and all phase
                      ! conventions are unchanged.
                      !
                      ! In the REAL executable DL* and CAPD* remain complex,
                      ! while A is real.  Treating the complex arrays as two
                      ! adjacent real rows lets DGEMM form real and imaginary
                      ! parts together.  The COMPLEX executable uses ZGEMM
                      ! directly.  The leading dimensions must be those of
                      ! the allocated arrays, not MAXINDEX, because LMAX2N
                      ! can be smaller than LMAX2.
!_REAL                lda=2*(LMAX2+1)*(LMAX2+1)
!_REAL                ldc=lda
!_REAL                ldb=ibpp_max
!_REAL                CALL dgemm('N','N',2*maxindex,nemax-nemin+1,isize, &
!_REAL                     1.d0,dl0ab_buf,lda,a(1,nemin),ldb,0.d0, &
!_REAL                     CAPD0LM_BUF(1,nemin),ldc)
!_REAL                CALL dgemm('N','N',2*maxindex,nemax-nemin+1,isize, &
!_REAL                     1.d0,dl1ab_buf,lda,a(1,nemin),ldb,0.d0, &
!_REAL                     CAPD1LM_BUF(1,nemin),ldc)
!_COMPLEX             lda=(LMAX2+1)*(LMAX2+1)
!_COMPLEX             ldc=lda
!_COMPLEX             ldb=ibpp_max
!_COMPLEX             CALL zgemm('N','N',maxindex,nemax-nemin+1,isize, &
!_COMPLEX                  (1.d0,0.d0),dl0ab_buf,lda,a(1,nemin),ldb, &
!_COMPLEX                  (0.d0,0.d0),CAPD0LM_BUF(1,nemin),ldc)
!_COMPLEX             CALL zgemm('N','N',maxindex,nemax-nemin+1,isize, &
!_COMPLEX                  (1.d0,0.d0),dl1ab_buf,lda,a(1,nemin),ldb, &
!_COMPLEX                  (0.d0,0.d0),CAPD1LM_BUF(1,nemin),ldc)
                      CAPD0LM(alpha,beta,1:maxindex,nemin:nemax) = &
                           CAPD0LM_BUF(1:maxindex,nemin:nemax)
                      CAPD1LM(alpha,beta,1:maxindex,nemin:nemax) = &
                           CAPD1LM_BUF(1:maxindex,nemin:nemax)
                                                                                          
        ENDDO                                                                                                           !END beta
     ENDDO                                                                                                              !END alpha

     IF(d_fd_check_enabled) THEN
        fd_step = 1.d-6
        fd_an_norm0 = 0.d0
        fd_an_norm1 = 0.d0
        fd_ref_norm0 = 0.d0
        fd_ref_norm1 = 0.d0
        fd_diff_norm0 = 0.d0
        fd_diff_norm1 = 0.d0
        fd_cross0 = 0.d0
        fd_cross1 = 0.d0
        DO i3 = 1,isize
           IF(KPLUSG(i3).le.0.d0) CYCLE
           fd_k0 = sqrt(KPLUSG(i3))
           ! '+' and '-' refer to the sign of the real-space strain.  A
           ! positive isotropic strain reduces the reciprocal length.
           fd_kplus  = fd_k0*(1.d0-fd_step)
           fd_kminus = fd_k0*(1.d0+fd_step)
           CALL SPHBES(lmax2,rmt*fd_kplus,fd_fjplus)
           CALL DVBES1(fd_fjplus,fd_dfjplus,rmt*fd_kplus,rmt,lmax2+1)
           CALL SPHBES(lmax2,rmt*fd_kminus,fd_fjminus)
           CALL DVBES1(fd_fjminus,fd_dfjminus,rmt*fd_kminus,rmt,lmax2+1)
           index = 0
           DO l = 0,lmax2
              pcfac = imag**l
              pcfac = pcfac*fac
              DO m = -l,l
                 index = index+1
                 fd_coeff0_plus = pcfac*conjg(h_yll11(index,i3))* &
                                  phshel(i3)*(fd_kplus*fd_dfjplus(l)*PE(l) - &
                                  fd_fjplus(l)*DPE(l))
                 fd_coeff0_minus = pcfac*conjg(h_yll11(index,i3))* &
                                   phshel(i3)*(fd_kminus*fd_dfjminus(l)*PE(l) - &
                                   fd_fjminus(l)*DPE(l))
                 fd_coeff1_plus = pcfac*conjg(h_yll11(index,i3))* &
                                  phshel(i3)*(fd_fjplus(l)*DP(l) - &
                                  fd_kplus*fd_dfjplus(l)*P(l))
                 fd_coeff1_minus = pcfac*conjg(h_yll11(index,i3))* &
                                   phshel(i3)*(fd_fjminus(l)*DP(l) - &
                                   fd_kminus*fd_dfjminus(l)*P(l))
                 fd_deriv0 = (fd_coeff0_plus-fd_coeff0_minus)/(2.d0*fd_step)
                 fd_deriv1 = (fd_coeff1_plus-fd_coeff1_minus)/(2.d0*fd_step)
                 fd_analytic0 = dl0ab_trace(index,i3)
                 fd_analytic1 = dl1ab_trace(index,i3)
                 fd_an_norm0 = fd_an_norm0 + abs(fd_analytic0)**2
                 fd_an_norm1 = fd_an_norm1 + abs(fd_analytic1)**2
                 fd_ref_norm0 = fd_ref_norm0 + abs(fd_deriv0)**2
                 fd_ref_norm1 = fd_ref_norm1 + abs(fd_deriv1)**2
                 fd_diff_norm0 = fd_diff_norm0 + &
                                 abs(fd_analytic0-fd_deriv0)**2
                 fd_diff_norm1 = fd_diff_norm1 + &
                                 abs(fd_analytic1-fd_deriv1)**2
                 fd_cross0 = fd_cross0 + &
                             dble(conjg(fd_analytic0)*fd_deriv0)
                 fd_cross1 = fd_cross1 + &
                             dble(conjg(fd_analytic1)*fd_deriv1)
              ENDDO
           ENDDO
        ENDDO
        fd_scale0 = 0.d0
        fd_scale1 = 0.d0
        fd_relerr0 = 0.d0
        fd_relerr1 = 0.d0
        IF(fd_an_norm0.gt.0.d0) fd_scale0 = fd_cross0/fd_an_norm0
        IF(fd_an_norm1.gt.0.d0) fd_scale1 = fd_cross1/fd_an_norm1
        IF(fd_ref_norm0.gt.0.d0) &
           fd_relerr0 = sqrt(fd_diff_norm0/fd_ref_norm0)
        IF(fd_ref_norm1.gt.0.d0) &
           fd_relerr1 = sqrt(fd_diff_norm1/fd_ref_norm1)
        write(6,'(a,3i6,6es22.12)') ':DBG_STR_D_FD ', &
             jatom,latom,lfirst,fd_scale0,fd_relerr0, &
             sqrt(fd_an_norm0),fd_scale1,fd_relerr1,sqrt(fd_an_norm1)
     ENDIF
     
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
                         IF(d_component_energy) THEN
                            buf_sph_val0(alpha,beta,num) = &
                                 buf_sph_val0(alpha,beta,num) + &
                                 (el(l)-e(num))*( &
                                 CONJG(CAPD0LM(alpha,beta,index,num))* &
                                 ABLM(0,index,num) + &
                                 CONJG(CAPD1LM(alpha,beta,index,num))* &
                                 ABLM(1,index,num)*PEI(l) )
                         ENDIF
                         IF(d_component_cross) THEN
                            buf_sph_val0(alpha,beta,num) = &
                                 buf_sph_val0(alpha,beta,num) + &
                                 CONJG(CAPD0LM(alpha,beta,index,num))* &
                                 ABLM(1,index,num)
                         ENDIF
                      ENDDO
!                       !For the moment this below part is only for C-Diamond test
      before_lo=buf_sph_val0(1,1,1)
                      IF(d_component_lo) THEN
                         DO num = nemin,nemax
                            DO jlo = 1, ilo(l)
                               IF(loor(jlo,l)) THEN
                                 buf_sph_val0(alpha,beta,num) = &
                                 buf_sph_val0(alpha,beta,num) + &
                                 (CONJG(CAPD0LM(alpha,beta,index,num))* &
                                  jlaplo1_02(jlo,l) + &
                                  CONJG(CAPD1LM(alpha,beta,index,num))* &
                                  jlaplo1_12(jlo,l))* &
                                  (elo(l,jlo)-e(num))*clm(index,num,jlo)
                               ENDIF
                            ENDDO
                         ENDDO
                      ENDIF
                         
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
                                                 (elo(lp,jlop)-e(num))*jr3upu_22(jlo,jlop,l,lp)
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

    ! Fourth line of Blaha et al. Eq. (G2): contraction of the strain
    ! derivative D of the APW matching coefficients with the nonspherical
    ! effective-potential matrix.  LAPW1 already writes these matrix elements
    ! for forces; reusing them here guarantees the same Gaunt, phase, and
    ! local-orbital conventions as the Hamiltonian.  Return this contribution
    ! separately because the following W and angular-V lines of Eq. (G2) must
    ! still be implemented before this can become a production contribution.
    IF (g2_diag_enabled) THEN
    DO ih = 1,ihmx
       lpy = lpv(ih)*(lpv(ih)+1)+mpv(ih)+1
       ly  = lv(ih) *(lv(ih)+1) +mv(ih) +1
       DO alpha = 1,3
          DO beta = 1,3
             DO num = nemin,nemax
                afac_nsp = conjg(ALM(lpy,num))*tuu(ih) + &
                            conjg(BLM(lpy,num))*tdu(ih)
                bfac_nsp = conjg(ALM(lpy,num))*tud(ih) + &
                            conjg(BLM(lpy,num))*tdd(ih)
                DO jlo = 1,ilo(lpv(ih))
                   afac_nsp = afac_nsp + &
                        conjg(clm(lpy,num,jlo))*tuu21(jlo,ih)
                   bfac_nsp = bfac_nsp + &
                        conjg(clm(lpy,num,jlo))*tud21(jlo,ih)
                ENDDO
                buf_nsp_match(alpha,beta,num) = &
                     buf_nsp_match(alpha,beta,num) + &
                     afac_nsp*CAPD0LM(alpha,beta,ly,num) + &
                     bfac_nsp*CAPD1LM(alpha,beta,ly,num)
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    ! Fifth line of Eq. (G2): radial-coordinate variation of the
    ! nonspherical potential matrix, containing
    ! W = integral r*u_l*V_L*u_l' dr.
    !
    ! A crucial distinction from the force/Hamiltonian contraction is that
    ! strain inserts two rank-one tensors.  The resulting three-Gaunt product
    ! can be nonzero even when the original Hamiltonian Gaunt
    !
    !                 G_{l,L,l'}^{m,M,m'}
    !
    ! vanishes.  The force file, however, contains one record only for the
    ! Hamiltonian-allowed pairs -m+M+m'=0.  The radial integral in that record
    ! is independent of m and m'.  We must therefore take one representative
    ! record for each (l,L,M,l',potential-phase) radial channel and explicitly
    ! sum *all* bra and ket magnetic quantum numbers.  Looping directly over
    ! the force records (the former implementation below) silently omitted
    ! most of Eq. (G2) and generated a large spurious anisotropy.
    DO ih = 1,ihmx
       l  = lpv(ih)
       lp = lv(ih)

       ! LQIND repeats the same radial channel for every Hamiltonian-allowed
       ! (m,m') pair.  Retain its first occurrence.  V/W values include the
       ! real/cubic-harmonic phase, so the comparison keeps distinct potential
       ! harmonics separate even if they have the same L and |M| labels.
       duplicate_radial_channel = .FALSE.
       ! ATPAR/FORFHS writes all (m,m') replicas of one radial channel
       ! consecutively, so comparison with the preceding record is sufficient
       ! and avoids an O(ihmx**2) search at every k point.
       IF (ih.GT.1 .AND. lpv(ih-1).EQ.l .AND. lv(ih-1).EQ.lp .AND. &
           llv(ih-1).EQ.llv(ih) .AND. mmv(ih-1).EQ.mmv(ih)) THEN
          channel_diff = ABS(vpuu(ih-1)-vpuu(ih))+ &
                         ABS(vpdd(ih-1)-vpdd(ih))+ &
                         ABS(vpud(ih-1)-vpud(ih))+ &
                         ABS(vpdu(ih-1)-vpdu(ih))+ &
                         ABS(wpuu(ih-1)-wpuu(ih))+ &
                         ABS(wpdd(ih-1)-wpdd(ih))+ &
                         ABS(wpud(ih-1)-wpud(ih))+ &
                         ABS(wpdu(ih-1)-wpdu(ih))
          channel_scale = ABS(vpuu(ih))+ABS(vpdd(ih))+ABS(vpud(ih))+ &
                          ABS(vpdu(ih))+ABS(wpuu(ih))+ABS(wpdd(ih))+ &
                          ABS(wpud(ih))+ABS(wpdu(ih))
          duplicate_radial_channel = &
               channel_diff.LE.1.d-12*(1.d0+channel_scale)
       ENDIF
       IF (duplicate_radial_channel) CYCLE

       DO mbra = -l,l
          index1 = l*(l+1)+mbra+1
          DO mket = -lp,lp
             index2 = lp*(lp+1)+mket+1
             radial_geom = zeroc
             DO sp = 0,2
                IF (NOTRI(sp,1,1).LT.0) CYCLE
                DO tp = -sp,sp
                   DO t = -1,1
                      DO u = -1,1
                         ! GAUNT does not reliably return zero for a
                         ! magnetically forbidden triplet.  Every historical
                         ! stress loop tests -m1+m2+m3=0 before calling it;
                         ! do the same here.  Omitting this guard fed
                         ! undefined values into the experimental G2 term.
                         IF (-tp+t+u.NE.0) CYCLE
                         gaunt1 = GAUNT(sp,1,1,tp,t,u)
                         IF (ABS(gaunt1).LT.1.d-30) CYCLE
                         ! The recoupled angular momentum is not a basis
                         ! expansion index.  It may exceed lmax2n; the exact
                         ! finite upper bound follows from L'' x l'.
                         DO v1 = 0,llv(ih)+lp
                            IF (NOTRI(v1,llv(ih),lp).LT.0) CYCLE
                            IF (NOTRI(l,sp,v1).LT.0) CYCLE
                            DO v = -v1,v1
                               IF (-v+mmv(ih)+mket.NE.0) CYCLE
                               gaunt_pot = GAUNT(v1,llv(ih),lp,v, &
                                                mmv(ih),mket)
                               IF (ABS(gaunt_pot).LT.1.d-30) CYCLE
                               IF (-mbra+tp+v.NE.0) CYCLE
                               gaunt_couple = GAUNT(l,sp,v1,mbra,tp,v)
                               IF (ABS(gaunt_couple).LT.1.d-30) CYCLE
                               DO alpha = 1,3
                                  DO beta = 1,3
                                     radial_factor = &
                                          conjg(c_am(alpha,t))* &
                                          conjg(c_am(beta,u))*gaunt1* &
                                          gaunt_pot*gaunt_couple
                                     radial_geom(alpha,beta) = &
                                          radial_geom(alpha,beta)+radial_factor
                                  ENDDO
                               ENDDO
                            ENDDO
                         ENDDO
                      ENDDO
                   ENDDO
                ENDDO
             ENDDO
             ! The radial matrix has no dependence on the intermediate
             ! Gaunt indices.  Contract it only once per (m,m',band), after
             ! the complete geometrical tensor has been assembled above.
             DO num = nemin,nemax
                ! Keep the four radial-function sectors separate.  Their
                ! sum is the published fifth line of Eq. (G2), while the
                ! split identifies whether a bad result enters through the
                ! ordinary APW block or through a particular LO channel.
                wmat_parts = zeroc
                wmat_parts(1) = &
                     conjg(ALM(index1,num))*( &
                       wpuu(ih)*ALM(index2,num) + &
                       wpud(ih)*BLM(index2,num) ) + &
                     conjg(BLM(index1,num))*( &
                       wpdu(ih)*ALM(index2,num) + &
                       wpdd(ih)*BLM(index2,num) )
                DO jlo = 1,ilo(l)
                   wmat_parts(2) = wmat_parts(2) + &
                        conjg(clm(index1,num,jlo))*( &
                        wpuu21(jlo,ih)*ALM(index2,num) + &
                        wpud21(jlo,ih)*BLM(index2,num) )
                ENDDO
                DO klo = 1,ilo(lp)
                   wmat_parts(3) = wmat_parts(3) + ( &
                        conjg(ALM(index1,num))*wpuu12(klo,ih) + &
                        conjg(BLM(index1,num))*wpdu12(klo,ih) )* &
                        clm(index2,num,klo)
                   DO jlo = 1,ilo(l)
                      wmat_parts(4) = wmat_parts(4) + &
                           conjg(clm(index1,num,jlo))* &
                           wpuu22(jlo,klo,ih)*clm(index2,num,klo)
                   ENDDO
                ENDDO
                wmat_nsp = SUM(wmat_parts)
                buf_nsp_radial(:,:,num) = buf_nsp_radial(:,:,num) + &
                                           radial_geom(:,:)*wmat_nsp
                DO ii=1,4
                   buf_nsp_radial_parts(:,:,num,ii) = &
                        buf_nsp_radial_parts(:,:,num,ii) + &
                        radial_geom(:,:)*wmat_parts(ii)
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    ! Independent hydrostatic check of the recoupled fifth line of Eq. (G2).
    ! Under an isotropic strain, sum_alpha rhat_alpha*rhat_alpha = 1.  The
    ! three-Gaunt tensor above must therefore collapse to the ordinary Gaunt
    ! coefficient already stored in the historical force matrices WUU,... .
    ! This one-Gaunt contraction uses those matrices directly and shares no
    ! phase conversion or recoupling algebra with BUF_NSP_RADIAL.  LAPW2
    ! prints the fully k/atom-summed comparison under WIEN_STRESS_G2_DIAG.
    DO ih=1,ihmx
       l=lpv(ih)
       lp=lv(ih)
       index1=l*(l+1)+mpv(ih)+1
       index2=lp*(lp+1)+mv(ih)+1
       gaunt_orig=GAUNT(l,llv(ih),lp,mpv(ih),mmv(ih),mv(ih))
       DO num=nemin,nemax
          ! Before testing the radial derivative W, contract the ordinary
          ! nonspherical potential matrix written by LAPW1.  After the
          ! k-point/occupation sum this must equal integral rho*V_ns over
          ! all muffin tins.  This is an independent diagnostic of radial
          ! normalization, real-harmonic phases, LO sectors, and weights.
          wmat_parts=zeroc
          wmat_parts(1)=conjg(ALM(index1,num))*( &
               tuu(ih)*ALM(index2,num)+tud(ih)*BLM(index2,num)) + &
               conjg(BLM(index1,num))*( &
               tdu(ih)*ALM(index2,num)+tdd(ih)*BLM(index2,num))
          DO jlo=1,ilo(l)
             wmat_parts(2)=wmat_parts(2)+conjg(clm(index1,num,jlo))*( &
                  tuu21(jlo,ih)*ALM(index2,num)+ &
                  tud21(jlo,ih)*BLM(index2,num))
          ENDDO
          DO klo=1,ilo(lp)
             wmat_parts(3)=wmat_parts(3)+( &
                  conjg(ALM(index1,num))*tuu12(klo,ih)+ &
                  conjg(BLM(index1,num))*tdu12(klo,ih))* &
                  clm(index2,num,klo)
             DO jlo=1,ilo(l)
                wmat_parts(4)=wmat_parts(4)+ &
                     conjg(clm(index1,num,jlo))*tuu22(jlo,klo,ih)* &
                     clm(index2,num,klo)
             ENDDO
          ENDDO
          buf_nsp_v_direct(num)=buf_nsp_v_direct(num)+SUM(wmat_parts)
          ! Preserve the four basis sectors of the ordinary V_ns matrix
          ! element independently.  Their sum is exactly the scalar above:
          !   1 APW/APW, 2 LO-bra/APW-ket, 3 APW-bra/LO-ket, 4 LO/LO.
          ! This is diagnostic bookkeeping only.  Comparing each occupation-
          ! summed sector against the density integral identifies whether a
          ! failure of <V_ns>=integral rho*V_ns is introduced by local
          ! orbitals before the more delicate radial derivative W is tested.
          DO ii=1,4
             buf_nsp_v_parts(num,ii)=buf_nsp_v_parts(num,ii)+ &
                                      wmat_parts(ii)
          ENDDO

          wmat_parts=zeroc
          wmat_parts(1)=conjg(ALM(index1,num))*( &
               wuu(ih)*ALM(index2,num)+wud(ih)*BLM(index2,num)) + &
               conjg(BLM(index1,num))*( &
               wdu(ih)*ALM(index2,num)+wdd(ih)*BLM(index2,num))
          DO jlo=1,ilo(l)
             wmat_parts(2)=wmat_parts(2)+conjg(clm(index1,num,jlo))*( &
                  wuu21(jlo,ih)*ALM(index2,num)+ &
                  wud21(jlo,ih)*BLM(index2,num))
          ENDDO
          DO klo=1,ilo(lp)
             wmat_parts(3)=wmat_parts(3)+( &
                  conjg(ALM(index1,num))*wuu12(klo,ih)+ &
                  conjg(BLM(index1,num))*wdu12(klo,ih))* &
                  clm(index2,num,klo)
             DO jlo=1,ilo(l)
                wmat_parts(4)=wmat_parts(4)+ &
                     conjg(clm(index1,num,jlo))*wuu22(jlo,klo,ih)* &
                     clm(index2,num,klo)
             ENDDO
          ENDDO
          buf_nsp_radial_direct(num)=buf_nsp_radial_direct(num)+ &
                                      SUM(wmat_parts)

          ! Repeat the same ordinary-Gaunt contraction using the new
          ! phase-weighted bare W matrices.  Comparing this with the legacy
          ! contraction above isolates the real/cubic-harmonic conversion
          ! from the subsequent three-Gaunt recoupling.
          wmat_parts=zeroc
          wmat_parts(1)=conjg(ALM(index1,num))*( &
               wpuu(ih)*ALM(index2,num)+wpud(ih)*BLM(index2,num)) + &
               conjg(BLM(index1,num))*( &
               wpdu(ih)*ALM(index2,num)+wpdd(ih)*BLM(index2,num))
          DO jlo=1,ilo(l)
             wmat_parts(2)=wmat_parts(2)+conjg(clm(index1,num,jlo))*( &
                  wpuu21(jlo,ih)*ALM(index2,num)+ &
                  wpud21(jlo,ih)*BLM(index2,num))
          ENDDO
          DO klo=1,ilo(lp)
             wmat_parts(3)=wmat_parts(3)+( &
                  conjg(ALM(index1,num))*wpuu12(klo,ih)+ &
                  conjg(BLM(index1,num))*wpdu12(klo,ih))* &
                  clm(index2,num,klo)
             DO jlo=1,ilo(l)
                wmat_parts(4)=wmat_parts(4)+ &
                     conjg(clm(index1,num,jlo))*wpuu22(jlo,klo,ih)* &
                     clm(index2,num,klo)
             ENDDO
          ENDDO
          buf_nsp_radial_phase(num)=buf_nsp_radial_phase(num)+ &
                                     gaunt_orig*SUM(wmat_parts)
       ENDDO
    ENDDO

    ! Sixth (last) line of Eq. (G2): change of the angular part of the
    ! nonspherical potential matrix element.  The published equation uses
    ! the symbol s both for the outer shift s=+/-1 and for an independent
    ! intermediate angular momentum.  The latter cannot be negative; here
    ! V is the +/-1 shift and LINT is the independent angular momentum.
    !
    ! LAPW1 supplies the phase-weighted bare V radial matrix, so GAUNT_POT
    ! installs the new intermediate angular coupling without any division by
    ! the original Hamiltonian Gaunt.  APW and LO channels are both present.
    ! Keeping this term separate prevents an unvalidated G2 block from
    ! entering production.
    DO ih = 1,ihmx
       l  = lpv(ih)
       lp = lv(ih)
       duplicate_radial_channel = .FALSE.
       IF (ih.GT.1 .AND. lpv(ih-1).EQ.l .AND. lv(ih-1).EQ.lp .AND. &
           llv(ih-1).EQ.llv(ih) .AND. mmv(ih-1).EQ.mmv(ih)) THEN
          channel_diff = ABS(vpuu(ih-1)-vpuu(ih))+ &
                         ABS(vpdd(ih-1)-vpdd(ih))+ &
                         ABS(vpud(ih-1)-vpud(ih))+ &
                         ABS(vpdu(ih-1)-vpdu(ih))+ &
                         ABS(wpuu(ih-1)-wpuu(ih))+ &
                         ABS(wpdd(ih-1)-wpdd(ih))+ &
                         ABS(wpud(ih-1)-wpud(ih))+ &
                         ABS(wpdu(ih-1)-wpdu(ih))
          channel_scale = ABS(vpuu(ih))+ABS(vpdd(ih))+ABS(vpud(ih))+ &
                          ABS(vpdu(ih))+ABS(wpuu(ih))+ABS(wpdd(ih))+ &
                          ABS(wpud(ih))+ABS(wpdu(ih))
          duplicate_radial_channel = &
               channel_diff.LE.1.d-12*(1.d0+channel_scale)
       ENDIF
       IF (duplicate_radial_channel) CYCLE

       DO mbra = -l,l
          index1 = l*(l+1)+mbra+1
          CALL c_factor(l,mbra,ca_st_lm)
          DO mket = -lp,lp
             index2 = lp*(lp+1)+mket+1
             angular_geom = zeroc
             DO v = -1,1,2
                IF (l+v.LT.0 .OR. l+v.GT.lmax2n) CYCLE
                DO t = -1,1
                   IF (mbra+t.LT.-(l+v) .OR. &
                       mbra+t.GT. (l+v)) CYCLE
                   DO tp = -1,1
                      ! Likewise this intermediate momentum belongs to the
                      ! four-harmonic product, not to the LAPW expansion.
                      ! Truncating it at lmax2n violates the hydrostatic trace
                      ! identity of the angular derivative.
                      DO lint = 0,llv(ih)+lp
                         IF (NOTRI(lint,llv(ih),lp).LT.0) CYCLE
                         IF (NOTRI(lint,1,l+v).LT.0) CYCLE
                         DO nu = -lint,lint
                            IF (-nu+mmv(ih)+mket.NE.0) CYCLE
                            gaunt_pot = GAUNT(lint,llv(ih),lp,nu, &
                                             mmv(ih),mket)
                            IF (ABS(gaunt_pot).LT.1.d-30) CYCLE
                            IF (-nu+tp+mbra+t.NE.0) CYCLE
                            gaunt_couple = GAUNT(lint,1,l+v,nu,tp, &
                                                mbra+t)
                            IF (ABS(gaunt_couple).LT.1.d-30) CYCLE
                            DO alpha = 1,3
                               DO beta = 1,3
                                  angular_factor = hf*( &
                                       conjg(c_am(alpha,tp))* &
                                       conjg(ca_st_lm(beta,v,t)) + &
                                       conjg(c_am(beta,tp))* &
                                       conjg(ca_st_lm(alpha,v,t)) )
                                  angular_geom(alpha,beta) = &
                                       angular_geom(alpha,beta) + &
                                       angular_factor*gaunt_pot*gaunt_couple
                               ENDDO
                            ENDDO
                         ENDDO
                      ENDDO
                   ENDDO
                ENDDO
             ENDDO
             DO num = nemin,nemax
                vmat_nsp = &
                     conjg(ALM(index1,num))*( &
                       vpuu(ih)*ALM(index2,num) + &
                       vpud(ih)*BLM(index2,num) ) + &
                     conjg(BLM(index1,num))*( &
                       vpdu(ih)*ALM(index2,num) + &
                       vpdd(ih)*BLM(index2,num) )
                DO jlo = 1,ilo(l)
                   vmat_nsp = vmat_nsp + &
                        conjg(clm(index1,num,jlo))*( &
                        vpuu21(jlo,ih)*ALM(index2,num) + &
                        vpud21(jlo,ih)*BLM(index2,num) )
                ENDDO
                DO klo = 1,ilo(lp)
                   vmat_nsp = vmat_nsp + ( &
                        conjg(ALM(index1,num))*vpuu12(klo,ih) + &
                        conjg(BLM(index1,num))*vpdu12(klo,ih) )* &
                        clm(index2,num,klo)
                   DO jlo = 1,ilo(l)
                      vmat_nsp = vmat_nsp + &
                           conjg(clm(index1,num,jlo))* &
                           vpuu22(jlo,klo,ih)*clm(index2,num,klo)
                   ENDDO
                ENDDO
                buf_nsp_angular(:,:,num) = buf_nsp_angular(:,:,num) + &
                                            angular_geom(:,:)*vmat_nsp
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    ENDIF
    
    if(rela_pulay_d_apply) then
       call rela_pulay_d_corr(nemin,nemax,jatom,alm,blm,clm, &
                              capd0lm,capd1lm,buf_rela_pulay_d)
       buf_sph_val0=buf_sph_val0+buf_rela_pulay_d
    endif

    DO alpha = 1,3
       DO beta = 1,3
          DO num = nemin,nemax
             ! Keep the three lines of the spherical part of Eq. (G2)
             ! separately for geometry diagnostics.  WIEN2k keeps R_MT
             ! fixed and spherical under strain, whereas the published
             ! derivation differentiates an affinely transported local
             ! coordinate.  A missing shape derivative would cancel some
             ! or all of the radial/angular pieces.  These extra outputs do
             ! not change BUF_SPH_VAL_FIN or the production stress.
             buf_sph_val_parts(alpha,beta,num,1) = &
                  g2_d_scale*buf_sph_val0(alpha,beta,num)
             buf_sph_val_parts(alpha,beta,num,2) = &
                  buf_sph_val1(alpha,beta,num)
             buf_sph_val_parts(alpha,beta,num,3) = &
                  buf_sph_val2(alpha,beta,num)
             buf_sph_val_fin(alpha,beta,num) = &
                  sum(buf_sph_val_parts(alpha,beta,num,:))
          ENDDO
       ENDDO
    ENDDO
!    write(6,*) buf_sph_val_fin(1,1,1),buf_sph_val0(1,1,1)
!  977 continue        
        
  RETURN
998 write(6,*) 'Invalid WIEN_STRESS_G2_D_SCALE: ',trim(g2_d_scale_env)
    stop 'WIEN_STRESS_G2_D_SCALE must be a real number'
END SUBROUTINE val_corr_stress
