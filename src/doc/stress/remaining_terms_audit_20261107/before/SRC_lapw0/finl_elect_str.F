subroutine finl_elect_str(lm,lmmax,coulomb_str_int,E_xc_m_V_xc,exc_vxc)
!subroutine finl_elect_str(clm1,lmmax1,lm,imax,lmmax,r,dx,iatnr,nkk)
!    use pw_mult,  only : cabt
   use densit ,    only : lmmaxx,clm
   use work   ,    only : V, exc, vxc
   use struct ,    only : zz, mult, jri,RMT,dx,iatnr,nat,r0,vol,irel
   use parallel !,   only : nstart,nstop,myid 
   use gga_str_kb, only : gga_fin_str
   implicit none
   include 'param.inc'
#ifdef Parallel
      include 'mpif.h'
#endif

   integer, intent(in)  :: lm(2,ncom+3,nat),lmmax(nat)
!   real*8, intent(in)   :: clm1(1:nrad,1:lmmaxx), dx, r(1:nrad)!,vpot(1:nrad,1:lmmaxx) potential taken from module
   complex*16            :: coulomb_str_int(1:9),ctest
   real*8    ,intent(in) :: E_xc_m_V_xc
   real*8               ::  exc_vxc(1:9)
   real*8               ::  exc_vxc_buf(1:9)
   logical              :: HEXA,CUBIC, stress_debug
   character(len=32)    :: stress_debug_env
   integer              :: stress_debug_status
   integer, external    :: NOTRI
   real*8 , external    :: GAUNT
   integer              :: lmmax_2,lm1,lm1p, m0,l0,l0p,m0p,t,tp,v1,s, lmtot1,        &
                           lmtot(2,ncom+3), alpha, beta, j,        &
                           index2,jatom,mu,imax, lmmax1(nat), index_a
                           
   real*8               :: VALUE(1:nrad), out1,out2,gaunt1,       & 
                           gaunt2, test_col_sum, zero
                                       
                           
                           
   real*8               :: pi, RMT1,RMT_i, sqrt_4pi, sqrt_4pii,fopi_i, dx1
   complex*16           :: zeroc, fc(ncom+3,nat), clm2(1:nrad,ncom+3,nat),vpot(1:nrad,ncom+3,nat), &
                           value_c(1:nrad), out_c                             
   COMPLEX*16           :: cabt(3,-1:1)                           
   real*8               :: r(1:nrad), pone                           
   real*8               :: tot_stress(1:3,1:3)   !for exachange-correlation lda   
   real*8               :: col_buf_pot(1:9), col_pot_1(1:9), madelung_1(1:9),madelung_2(1:9),madelung_3(1:9), &
                           madelung_1_buf(nat,1:9),madelung_2_buf(nat,1:9),madelung_3_buf(nat,1:9), &
                           col_pot_final(1:9),madelung_final(1:9),madelung_final_sum(1:9), sum_density
   
   zeroc     = (0.d0,0.d0)
   zero      = 0.d0
   pi        = 4.0d+0*ATAN(1.d0)
   fopi_i    = 1.d0/( 4.d0*pi )
   sqrt_4pi  = sqrt( 4.d0*pi )
   sqrt_4pii = 1.d0/sqrt( 4.d0*pi )
   !====intialized the variables=======
   madelung_1     = zero
   madelung_2     = zero
   madelung_3     = zero
   madelung_final = zero
   col_pot_1      = zero   
   col_pot_final  = zero

   ! Optional, machine-readable stress diagnostics.  Keep the normal
   ! production output unchanged unless WIEN_STRESS_DEBUG is set to a
   ! nonempty value other than "0".
   stress_debug_env = ''
   call get_environment_variable('WIEN_STRESS_DEBUG',stress_debug_env, &
                                 status=stress_debug_status)
   stress_debug = stress_debug_status.eq.0 .and. &
                  len_trim(stress_debug_env).gt.0 .and. &
                  trim(stress_debug_env).ne.'0'
   
   madelung_1_buf = zero
   madelung_2_buf = zero
   madelung_3_buf = zero         
   
   call c_alpha_m(cabt)                     !this is c matrix of unit vector expansion
   
!  index_a addresses ROTIJ for every atom in the full unit cell, whereas
!  jatom addresses only inequivalent atom types.  In an MPI calculation a
!  rank may start at an inequivalent type later than one; initializing the
!  full-atom index to zero then applies the rotations belonging to earlier
!  atoms.  Advance over the multiplicities owned by preceding ranks before
!  entering this rank's local atom-type loop.  The serial result is
!  unchanged because nstart=1.
   index_a = 0
   DO jatom = 1, nstart(myid)-1
      index_a = index_a + mult(jatom)
   ENDDO
   DO 29 jatom = nstart(myid), nstop(myid)         
   DO 30 mu    = 1, mult(jatom)         
         index_a   = index_a + 1
         RMT1      = RMT(jatom)
         RMT_i     = 1.d0/RMT1   
         dx1       = dx(jatom)
         imax      = jri(jatom)
         HEXA      = .FALSE.
         CUBIC     = .FALSE.
         IF(iatnr(jatom) .lt. 0) HEXA  = .TRUE.
         IF(iatnr(jatom) .gt. 0) CUBIC = .TRUE.
      
         DO J = 1, imax
            r(j) = r0(jatom)*exp( dx(jatom)*dble(j-1) )
         ENDDO
   
   DO lm1 = 1, lmmax(jatom)     
      clm2(1:imax,lm1,jatom) = clm(1:imax,lm1,jatom)
      vpot(1:imax,lm1,jatom) = V(1:imax,lm1,jatom)
      if(lm1 .eq. 1) clm2(1:imax,lm1,jatom) = clm2(1:imax,lm1,jatom)*sqrt_4pii            
   ENDDO
   
   IF(iatnr(jatom) .lt. 0) THEN         
   lmmax_2 = lmmax(jatom)   
     DO lm1 = 1, lmmax(jatom)
        M0 = lm(2,lm1,jatom)
        IF(M0 .NE. 0) THEN
          lmmax_2 = lmmax_2 + 1
          if(lmmax_2 .gt. ncom+3) STOP  'error in "finl_elect_str.f" for lm combination'
          clm2(1:imax,lmmax_2,jatom) = clm2(1:imax,lm1,jatom)
          vpot(1:imax,lmmax_2,jatom) = vpot(1:imax,lm1,jatom)
        ENDIF
     ENDDO
!   ENDIF

    call multfc(fc,jatom,lmmax_2,lm(1,1,jatom),iatnr(jatom))    
    
    DO lm1 = 1, lmmax_2       
       clm2(1:imax,lm1,jatom) = clm2(1:imax,lm1,jatom)*fc(lm1,jatom)      !note this is r2*rho
       vpot(1:imax,lm1,jatom) = vpot(1:imax,lm1,jatom)*fc(lm1,jatom)      !it is only v, but when it print in output file it will multiply by r2. 
    ENDDO
    
    lmtot1 = 0
    DO 39 lm1 = 1, lmmax_2
          lmtot1 = lmtot1 + 1
          lmtot(1,lmtot1) = lm(1,lm1,jatom)
          lmtot(2,lmtot1) = lm(2,lm1,jatom)
          clm2(1:imax,lmtot1,jatom) = clm2(1:imax,lm1,jatom)
          vpot(1:imax,lmtot1,jatom) = vpot(1:imax,lm1,jatom)
          IF(lmtot1 .eq. 1) GOTO 39
          IF( ( IABS( lmtot(1,lmtot1)) .eq. IABS( lmtot(1,lmtot1-1))  ) .AND.   &
              (  lmtot(2,lmtot1) .eq.  lmtot(2,lmtot1-1)  )  ) THEN
              clm2(1:imax,lmtot1-1,jatom) = clm2(1:imax,lmtot1-1,jatom) + clm2(1:imax,lmtot1,jatom)             !previously in last clm2 there was lm1 instead of lmtot1
              vpot(1:imax,lmtot1-1,jatom) = vpot(1:imax,lmtot1-1,jatom) + vpot(1:imax,lmtot1,jatom)
              lmtot1 = lmtot1 - 1
          ENDIF    
39  CONTINUE  
 !   go to 999    
    ELSE     
       !CUBIC CASE       
      ! call rho_combine(imax,lmmax(jatom),clm2(1,1,jatom),vpot(1,1,jatom),lm(1,1,jatom),jatom)              
       call lm_combine(imax,lmmax(jatom),clm2(1,1,jatom),lm(1,1,jatom),jatom)
       call lm_combine(imax,lmmax(jatom),vpot(1,1,jatom),lm(1,1,jatom),jatom)              
       
       lmtot1 = lmmax(jatom)
       DO lm1 = 1, lmmax(jatom)          
          lmtot(1,lm1) =  lm(1,lm1,jatom)
          lmtot(2,lm1) =  lm(2,lm1,jatom)                                        
          IF( lm(2,lm1,jatom) .ne. 0) THEN
            lmtot1 = lmtot1 + 1
            IF( lmtot1 .gt. ncom+3 ) STOP 'error in finl_elect_str '
            lmtot(1,lmtot1) =  lm(1,lm1,jatom)
            lmtot(2,lmtot1) = -lm(2,lm1,jatom)
         
            clm2(1:imax,lmtot1,jatom) = clm2(1:imax,lm1,jatom)           !uncommenting this one leads to the perfect unstymmetry
            vpot(1:imax,lmtot1,jatom) = vpot(1:imax,lm1,jatom)
          ENDIF
       ENDDO
              
!999 CONTINUE                 
    ENDIF    
    
    ! WRITE(21,*) ':COLUMB_RMT==========Coulomb Stress ========================='     
     col_buf_pot = zero
     sum_density = zero
     ctest=zero
     
     DO lm1 = 1, lmtot1             
        l0 = lmtot(1,lm1)
        m0 = lmtot(2,lm1)
        pone = (-1.0)**m0
        
        DO lm1p = 1, lmtot1
           l0p = lmtot(1,lm1p)
           m0p = lmtot(2,lm1p)
           
           DO j = 1, imax
!!              value_c(j) = conjg(clm2(j,lm1,jatom))*vpot(j,lm1p,jatom)                     
!! conjg ?? PB  value_c(j) = conjg(clm2(j,lm1,jatom))*vpot(j,lm1p,jatom)                     
             value_c(j) = (clm2(j,lm1,jatom))*vpot(j,lm1p,jatom)                          
           ENDDO
           VALUE(:) = REAL( value_c(:) )
           call charge(r,dx1,VALUE,1,imax,out1)
           VALUE(:) = AIMAG( value_c(:) )
           call charge(r,dx1,VALUE,1,imax,out2)
           out_c = dcmplx(out1,out2)
!write(6,'("out_c",6i3,2f20.10)') lm1,l0,m0,lm1p,l0p,m0p,out_c           
           DO s = 0,2
              IF( NOTRI(s,1,1) .lt. 0 ) CYCLE
              IF( NOTRI(l0,l0p,s) .lt. 0 ) CYCLE
              DO v1 = -s,s
                 IF( (m0+m0p+v1) .ne. 0 ) CYCLE                 
                 DO t = -1,1
                    DO tp = -1,1
                       IF( (-v1+t+tp) .ne. 0 ) CYCLE
                       gaunt1 = GAUNT(s,1,1,v1,t,tp)
                       gaunt2 = GAUNT(l0,l0p,s,-m0,m0p,v1)                                               
                       
                       index2 = 0
                       DO alpha = 1, 3
                          DO beta = 1, 3
                             index2 = index2 + 1
                             col_buf_pot(index2) = col_buf_pot(index2) - &
                             pone*cabt(alpha,t)*cabt(beta,tp)*gaunt1*gaunt2*out_c              
!                     if(index2.eq.2) then
!                             write(6,'(4i3,6f15.8)') s,v1,t,tp,cabt(alpha,t),cabt(beta,tp),gaunt1,gaunt2                    
!                             ctest=ctest-pone*cabt(alpha,t)*cabt(beta,tp)*gaunt1*gaunt2*out_c
!                     endif
                             ENDDO
                       ENDDO
                    ENDDO
                 ENDDO
              ENDDO
!              write(6,*) 'ctest',ctest
           ENDDO
        ENDDO        
     ENDDO   !lm1        
              write(6,*) 'col,ctest',col_buf_pot(2),ctest
     !symmetriyation among the equivalent atoms using the rotloc and rotij 
     call symmetry_stress_rotij(index_a,jatom,col_buf_pot)
     col_pot_1 = col_pot_1 + col_buf_pot                  
     
    ! WRITE(21,*) ':MADELUNG_STR  ==========Madelung Stress ========================='     
     madelung_1 = zero
     DO lm1 = 1, lmtot1        
        l0 = lmtot(1,lm1)
        m0 = lmtot(2,lm1)
        pone = (-1)**m0
          DO J = 1, imax
!!              value_c(j) = conjg( clm2(j,lm1,jatom) )*(1.d0/r(j)-RMT_i)               
             value_c(j) = clm2(j,lm1,jatom)*(1.d0/r(j)-RMT_i) 
          ENDDO
          VALUE(:) = REAL( value_c(:) )
          call charge(r,dx1,VALUE,1,imax,out1)
          VALUE(:) = AIMAG( value_c(:) )
          call charge(r,dx1,VALUE,1,imax,out2)   !instead of dx1, dx from module was being used which instigated error for the long time.
          out_c = dcmplx(out1,out2)                        
              
        IF(NOTRI(l0,1,1) .lt. 0) CYCLE
        DO t = -1,1
           DO tp = -1,1
              IF( (m0+t+tp) .ne. 0) CYCLE                          
              gaunt1 = GAUNT(l0,1,1,-m0,t,tp)                           
                           
              index2 = 0
              DO alpha = 1,3
                 DO beta = 1,3
                    index2 = index2 + 1                                        
                    
                   madelung_1(index2) = madelung_1(index2) + &
                   pone*cabt(alpha,t)*cabt(beta,tp)*gaunt1*out_c                                      
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
                    
     DO t = -1,1
        DO tp = -1,1
           IF( -t .ne. tp) CYCLE
           index2 = 0
           DO alpha = 1,3
              DO beta = 1,3
                 index2 = index2 + 1
                 madelung_2(index2) = madelung_2(index2) + &
                 (-1)**t*cabt(alpha,t)*cabt(beta,tp)*ZZ(jatom)*fopi_i*RMT_i
                 
                 madelung_2_buf(jatom,index2) = madelung_2_buf(jatom,index2) + &
                 (-1)**t*cabt(alpha,t)*cabt(beta,tp)*ZZ(jatom)*fopi_i*RMT_i
              ENDDO
           ENDDO
        ENDDO
     ENDDO
!====potential involves in the Madelung potential===============    
     madelung_3 = zero
     DO lm1 = 1, lmtot1          
        l0   = lmtot(1,lm1)
        m0   = lmtot(2,lm1)        
        pone = (-1)**m0
        IF( NOTRI(l0,1,1) .lt. 0 ) CYCLE
        DO t = -1,1
           DO tp = -1,1
              IF( (m0+t+tp) .ne. 0 ) CYCLE
              gaunt1 = GAUNT(l0,1,1,-m0,t,tp)
!! PB would be my guess      gaunt1 = GAUNT(l0,1,1,m0,t,tp)
              index2 = 0
              DO alpha = 1,3
                 DO beta = 1,3
                    index2 = index2 + 1
                    madelung_3(index2) = madelung_3(index2) + &
                    pone*cabt(alpha,t)*cabt(beta,tp)*0.5d0*vpot(imax,lm1,jatom)*gaunt1*fopi_i                                        
                 ENDDO
              ENDDO
        write(6,*) lm1,l0,m0,t,tp,madelung_3(1),madelung_3(5),madelung_3(9),vpot(imax,lm1,jatom),gaunt1
           ENDDO
        ENDDO
     ENDDO           
   
    call symmetry_stress_rotij(index_a,jatom,madelung_1)
     madelung_1_buf(jatom,:) = madelung_1_buf(jatom,:) + madelung_1(:)
write(6,'(3e10.3)') madelung_3
    call symmetry_stress_rotij(index_a,jatom,madelung_3)
     madelung_3_buf(jatom,:) = madelung_3_buf(jatom,:) + madelung_3(:)           
write(6,'(3e10.4)') madelung_3
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!=====This block is to calculate individual compoents of Madelung stress and commented since 4 June 2020
    WRITE(6,*)
    WRITE(6,*) ':MADELUNG_STR  ========== jatom',jatom,' mu ',mu   
    DO index2 = 1,9
       WRITE(6,'(I4,3f20.10)') index2, ZZ(jatom)*madelung_1_buf(jatom,index2),        &                                
                                        ZZ(jatom)*madelung_2_buf(jatom,index2),        &
                                        ZZ(jatom)*madelung_3_buf(jatom,index2)
    ENDDO
    WRITE(6,*)    
!====================== individual compoents of Madelung stress==========================================
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

30  CONTINUE
29  CONTINUE

   DO jatom = nstart(myid), nstop(myid)         
      DO index2 = 1, 9                  
         madelung_final(index2) = madelung_final(index2)+ &
          ZZ(jatom)*( madelung_1_buf(jatom,index2) + madelung_2_buf(jatom,index2)+ madelung_3_buf(jatom,index2) )         
         !this line added once we had problem doing NaCl, MgS, GaAs calculation
      ENDDO
   ENDDO 

!   gga_fin_str1=gga_fin_str
#ifdef Parallel
   call mpi_reduce(col_pot_1,col_buf_pot,9,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
   call mpi_reduce(madelung_final,madelung_final_sum,9,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
   call mpi_reduce(exc_vxc,exc_vxc_buf,9,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
!   call mpi_reduce(gga_fin_str1,gga_fin_str1_buf,9,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
   if(myid.eq.0) then
           madelung_final=madelung_final_sum
           col_pot_1=col_buf_pot
           exc_vxc=exc_vxc_buf
!           gga_fin_str1=gga_fin_str1_buf
   endif        
#endif
   
!    exc_vxc = 0.d0
!    call exchange_str(lm,lmmax,exc_vxc)
   
   DO index2 = 1, 9
      IF(index2 .eq. 1 .or. index2 .eq. 5 .or. index2 .eq. 9) THEN
        exc_vxc(index2) = exc_vxc(index2) + E_xc_m_V_xc      
      ENDIF            
   ENDDO   
   
!!***********************************************************************************
!===================individual compoents============================================
!    WRITE(21,*) ':COLUMB_RMT==========Coulomb Stress RMT========================='     
!    DO index2 = 1,9
!       WRITE(21,'(I4,2ES19.12)') index2, 0.5d0*col_pot_1(index2)
!    ENDDO
!    
!    
!    WRITE(21,*) ':COLUMB_INT==========Coulomb Stress Interstitial================='     
!    DO index2 = 1,9
!       WRITE(21,'(I4,2ES19.11)') index2, coulomb_str_int(index2)
!    ENDDO
!    
!    WRITE(21,*)
!    WRITE(21,*) ':MADELUNG_STR  ==========Madelung Stress ========================='     
!    DO index2 = 1,9
!       WRITE(21,'(I4,2f20.10)') index2, madelung_final(index2)        
!    ENDDO         
!    WRITE(21,*)
!    
!    WRITE(21,*)
!    WRITE(21,*) ':EXC-VXC    ==========Exchange correlation stress=============='     
!    DO index2 = 1,9
!       WRITE(21,'(I4,2f20.10)') index2, exc_vxc(index2)        
!    ENDDO         
!    WRITE(21,*)
!    
!    WRITE(21,*)
!    WRITE(21,*) ':EXC-VXC    ==========GGA correction stress=============='     
!    DO index2 = 1,9
!       WRITE(21,79) index2, gga_fin_str(index2)        
!    ENDDO         
!    WRITE(21,*)   
!*************************************************************************************   
   
   write(6,*) '===== Total STRESS DECOMPOSITION ======='
   write(6,'("  i  j     coul_pot_sp      madelung     Coul_pot_int     exc_vxc       GGA_cont")')
   index2 = 0
   DO beta = 1, 3
      DO alpha = 1, 3
         index2 = index2 + 1
         tot_stress(alpha,beta) =  0.5*REAL(col_pot_1(index2))  + &
                                   madelung_final(index2)       + &
                                   coulomb_str_int(index2)      + &
                                   exc_vxc(index2)              - &
                                   gga_fin_str(index2)                                        
         if(stress_debug .and. myid.eq.0) then
            write(21,81) alpha,beta, &
                         0.5d0*REAL(col_pot_1(index2))/vol, &
                         madelung_final(index2)/vol, &
                         REAL(coulomb_str_int(index2))/vol, &
                         exc_vxc(index2)/vol, &
                         -gga_fin_str(index2)/vol, &
                         tot_stress(alpha,beta)/vol
         endif
                         write(6,'(2i3,5f15.6)') alpha,beta, col_pot_1(index2)  , &
                                   madelung_final(index2)       , &
                                   real(coulomb_str_int(index2))      , &
                                   exc_vxc(index2)              ,&
                                   gga_fin_str(index2)
      ENDDO
   ENDDO
   tot_stress = tot_stress /vol
   
   if(myid.eq.0) then
     ! The electrostatic and exchange-correlation stress evaluated in this
     ! routine is a functional of the charge density and potential.  Nothing
     ! in finl_elect_str assumes a nonrelativistic valence Hamiltonian: RELA
     ! merely changes the self-consistent density supplied to LAPW0.  The old
     ! NREL-only output guard therefore discarded a valid contribution and,
     ! worse, allowed mixer to reuse stale :STR_DENS records from an earlier
     ! NREL calculation.  Always write the tensor when the STR path called
     ! this routine.  Relativistic kinetic/Pulay corrections belong to LAPW2,
     ! not to this density-functional contribution.
     WRITE(21,*)
     WRITE(21,80)    
     DO alpha = 1, 3
        WRITE(21,77) alpha, (tot_stress(alpha,beta),beta = 1, 3)
     ENDDO
     WRITE(21,*)
   endif  
!    DO index2 = 1, 9
!       tot_stress(index2) =  0.5*REAL(col_pot_1(index2))   + &
!                             REAL(madelung_final(index2))  + &
!                             REAL(coulomb_str_int(index2)) + &
!                             REAL(exc_vxc(index2))         !- &
! !                             REAL(gga_fin_str(index2))
!                                             
! !       WRITE(21,77) index2, tot_stress(index2)/vol
!    ENDDO         
   
   
   
77 FORMAT(':STR_DENS',i3.3,':',2x,3f20.12)      
78 FORMAT(':ELE_CORE_',i3.3,':',2x,f20.12,2x,f20.12,2x,f20.12)   
79 FORMAT(':GGA_CORR_',i3.3,':',2x,3f20.12)    
80 FORMAT(':STRESS (Ry/Bohr^3)')
81 FORMAT(':DBG_STR_DENS:',2I3,6(1X,ES23.15))
   
RETURN
end subroutine finl_elect_str

! subroutine rho_combine(imax,lmmax,clm2,vpot,lm1,jatom)
subroutine lm_combine(imax,lmmax,clm2,lm1,jatom)
   use struct, only : nat
   implicit none
   include 'param.inc'
   
   integer, intent(in)  :: imax,lmmax,lm1(2,ncom+3), jatom
   complex*16           :: clm2(1:nrad,ncom+3)!,vpot(1:nrad,ncom+3)
   real*8               :: c_kub(0:10,0:10)
   integer              :: i,j
   real*8               :: sq1,sqrt2,c1,c2,c3
!=======coefficient for real spherical harmonics in cubic system========   
  c_kub=0.0d0
  c_kub(0,0)=1.d0
  c_kub(3,2)=1.d0
  c_kub(4,0)=.5d0*SQRT(7.d0/3.d0)
  c_kub(4,4)=.5d0*SQRT(5.d0/3.d0)
  c_kub(6,0)=.5d0*SQRT(.5d0)
  c_kub(6,2)=.25d0*SQRT(11.d0)
  c_kub(6,4)=-.5d0*SQRT(7.d0/2.d0)
  c_kub(6,6)=-.25d0*SQRT(5.d0)
  c_kub(7,2)=.5d0*SQRT(13.d0/6.d0)
  c_kub(7,6)=.5d0*SQRT(11.d0/6.d0)
  c_kub(8,0)=.125d0*SQRT(33.d0)
  c_kub(8,4)=.25d0*SQRT(7.d0/3.d0)
  c_kub(8,8)=.125d0*SQRT(65.d0/3.d0)
  c_kub(9,2)=.25d0*SQRT(3.d0)
  c_kub(9,4)=.5d0*SQRT(17.d0/6.d0)
  c_kub(9,6)=-.25d0*SQRT(13.d0)
  c_kub(9,8)=-.5d0*SQRT(7.d0/6.d0)
  c_kub(10,0)=.125d0*SQRT(65.D0/6.D0)
  c_kub(10,2)=.125d0*SQRT(247.D0/6.D0)
  c_kub(10,4)=-.25d0*SQRT(11.D0/2.D0)
  c_kub(10,6)=0.0625d0*SQRT(19.D0/3.D0)
  c_kub(10,8)=-.125d0*SQRT(187.D0/6.D0)
  c_kub(10,10)=-.0625d0*SQRT(85.d0)
!===========================================================================   
  sqrt2 = sqrt(2.d0)
  
  i = 1
5 CONTINUE
     IF( i .gt. lmmax ) GOTO 6
     IF( ( lm1(1,i) .eq. 0 ) .and. ( lm1(2,i) .eq. 0 ) ) THEN
        i = i + 1
     ELSEIF ( ( lm1(1,i) .eq. -3 ) .and. ( lm1(2,i) .eq. 2 )) THEN
        DO j = 1, imax
           clm2(j,i) = clm2(j,i)/sqrt2
!            vpot(j,i) = vpot(j,i)/sqrt2
        ENDDO
        i = i + 1
     ELSEIF (( lm1(1,i) .eq. 4 ) .OR. ( lm1(1,i) .eq. 6 ) .OR. &
            ( lm1(1,i) .eq. -7 ) .OR. ( lm1(1,i) .eq. -9 )) THEN
            IF( lm1(2,i) .eq. 0 ) THEN
               sq1 = 1.d0
            ELSE
               sq1 = sqrt2
            ENDIF
            c1 = c_kub(IABS(lm1(1,i)),lm1(2,i) )
            c2 = c_kub(IABS(lm1(1,i)),lm1(2,i)+4 )
            DO J = 1, imax
               clm2(j,i)   = clm2(j,i)*c1 + clm2(j,i+1)*c2
               clm2(j,i+1) = clm2(j,i)*c2/sqrt2
               clm2(j,i)   = clm2(j,i)*c1/sq1
               
!                vpot(j,i)   = vpot(j,i)*c1 + vpot(j,i+1)*c2
!                vpot(j,i+1) = vpot(j,i)*c2/sqrt2
!                vpot(j,i)   = vpot(j,i)*c1/sq1
            ENDDO
               i = i + 2
     ELSEIF( ( lm1(1,i) .eq. 8 ) .OR.( lm1(1,i) .eq. 10 ) ) THEN       
           IF( lm1(2,i) .eq. 0 ) THEN
              sq1 = 1.d0
           ELSE
              sq1 = sqrt2
           ENDIF
           c1 = c_kub(IABS(lm1(1,i)),lm1(2,i) )
           c2 = c_kub(IABS(lm1(1,i)),lm1(2,i)+4 )
           ! l=8 and l=10 cubic harmonics contain three ordinary
           ! spherical-harmonic components.  The third coefficient belongs
           ! in c3.  Assigning it to c1 overwrote the first coefficient and
           ! left c3 undefined, corrupting the stress-only conversion.
           c3 = c_kub(IABS(lm1(1,i)),lm1(2,i)+8 )
           DO j = 1, imax
              clm2(j,i)   = clm2(j,i)*c1 + clm2(j,i+1)*c2 + clm2(j,i+2)*c3
              clm2(j,i+1) = clm2(j,i)*c2/sqrt2
              clm2(j,i+2) = clm2(j,i)*c3/sqrt2
              clm2(j,i)   = clm2(j,i)*c1/sq1
              
!               vpot(j,i)   = vpot(j,i)*c1 + vpot(j,i+1)*c2 + vpot(j,i+2)*c3
!               vpot(j,i+1) = vpot(j,i)*c2/sqrt2
!               vpot(j,i+2) = vpot(j,i)*c3/sqrt2
!               vpot(j,i)   = vpot(j,i)*c1/sq1
           ENDDO
     ELSE 
           WRITE(6,*) 'UNCORRECT LM LIST FOR CUBIC STRUCTURE'
           stop '=============incorrect LM in lm_combine========= '
     ENDIF      
  GOTO 5
6 CONTINUE  
  
  
RETURN
end subroutine lm_combine

subroutine exchange_str(lm,lmmax,exc_vxc)
   use densit,   only: clmsp
   use work  ,   only: vxc,exc
   use struct,   only: mult,jri,r0,dx,rmt,iatnr,nat
   use parallel, only: nstart,nstop,myid
   use norm_kub, only: c_kub
   implicit none
   include 'param.inc'   
   
   integer, intent(in)  :: lm(2,ncom+3,nat), lmmax(nat)
   logical              :: HEXA, CUBIC
   integer              :: jatom, lm1, mu, l0,m0,ispin, J, alpha,beta, &
                           index2
   real*8               :: r(nrad), rho,rho1, sqrt_4pi,sqrt_4pii, pi, &
                           VALUE(nrad),out1,exc_vxc(1:9)
   real*8               :: exc_vxc_buf1(ncom,nat), buf_exc_vxc,exc_vxc_buf                            
  
   pi        = 4.0d+0*ATAN(1.d0)
   sqrt_4pi  = sqrt( 4.D0*pi)
   sqrt_4pii = 1.d0/sqrt_4pi
   
   exc_vxc_buf = 0.d0   
   DO 20 jatom = nstart(myid), nstop(myid)
   DO 20 mu = 1, mult(jatom)
         HEXA  = .FALSE.
         CUBIC = .FALSE.
         IF(iatnr(jatom) .lt. 0) HEXA  = .TRUE.
         IF(iatnr(jatom) .gt. 0) CUBIC = .TRUE.
         
         DO J = 1, jri(jatom)
            r(j) = r0(jatom)*exp( dx(jatom)*(j-1) )
         ENDDO
         
         DO 40 lm1 = 1, lmmax(jatom)
            l0 = IABS( lm(1,lm1,jatom) )
            m0 =       lm(2,lm1,jatom)
            
            exc_vxc_buf1(lm1,jatom) = 0.d0
         DO 40 ispin = 1, 2
            
               IF(lm1 .eq. 1) THEN
               
                 DO j = 1, jri(jatom)
                    RHO = clmsp(j,lm1,jatom,ispin)*sqrt_4pii
                    VALUE(j) = rho*( exc(j,lm1,jatom,ispin)-vxc(j,lm1,jatom,ispin) )*0.5d0
                 ENDDO
               ELSEIF(CUBIC .and. (l0 .ne. 9 .and. m0 .eq. 4)) THEN
                  GOTO 40
               ELSEIF(CUBIC.and. (m0.eq.6 .or. m0.eq.8 .or. m0.eq. 10 ) ) THEN
                  GOTO 40
               ELSEIF(CUBIC.and. (l0.eq.4 .or. l0.eq.6 .or. l0.eq.7 .or. l0.eq.8 .or.l0.eq.9 .or.l0.eq.10 ) ) THEN
                 DO j = 1, jri(jatom)
                    rho = clmsp(j,lm1,jatom,ispin)*c_kub(l0,m0) + &
                          clmsp(j,lm1+1,jatom,ispin)*c_kub(l0,m0+4)
                    IF(l0 .eq. 8 .or. l0 .eq. 10) rho = rho+ clmsp(j,lm1+2,jatom,ispin)*c_kub(l0,m0+8)     
                    rho1 = ( Exc(j,lm1,jatom,ispin)-Vxc(j,lm1,jatom,ispin) )*c_kub(l0,m0) + &
                           ( Exc(j,lm1+1,jatom,ispin)-Vxc(j,lm1+1,jatom,ispin) )*c_kub(l0,m0+4)  
                    IF(l0 .eq. 8 .or. l0 .eq. 10) rho1 = rho1 + &
                        ( Exc(j,lm1+2,jatom,ispin)-Vxc(j,lm1+2,jatom,ispin) )*c_kub(l0,m0+8)
                        
                    VALUE(j) = rho*rho1*0.5d0    
                 ENDDO
                 
                 ELSE
                   DO J = 1, jri(jatom)
                      rho = clmsp(j,lm1,jatom,ispin)
                      VALUE(j) = rho*(Exc(j,lm1,jatom,ispin)-Vxc(j,lm1,jatom,ispin))*0.5d0
                   ENDDO
               ENDIF
               
               call charge(r,dx(jatom),VALUE,1,jri(jatom),out1)
               exc_vxc_buf1(lm1,jatom) = exc_vxc_buf1(lm1,jatom) + out1  !sum over all spin states              
    40 CONTINUE                             
    20 CONTINUE                                  
       
       DO 50 jatom = nstart(myid),nstop(myid)
          DO 50 mu = 1, mult(jatom)
             HEXA  = .FALSE.
             CUBIC = .FALSE.
             IF(iatnr(jatom) .lt. 0) HEXA  = .TRUE.
             IF(iatnr(jatom) .gt. 0) CUBIC = .TRUE.
             buf_exc_vxc = 0.d0
             DO 55 lm1 = 1, lmmax(jatom)
                l0 = IABS( lm(1,lm1,jatom) )
                m0 =       lm(2,lm1,jatom)                
                IF(CUBIC) THEN
                   IF( m0.eq.4 .or. m0.eq.6 .or. m0 .eq.8 .or. m0 .eq. 10 ) GOTO 55
                ENDIF
                buf_exc_vxc = buf_exc_vxc + exc_vxc_buf1(lm1,jatom)     !sum over all lm states
       55 CONTINUE     
             exc_vxc_buf = exc_vxc_buf + buf_exc_vxc                    !sum over equivalent and non-equivalent atoms
  50   CONTINUE     
       
     !  WRITE(*,*) '=======exc-uxc stress inside a sphere'
       index2 = 0
       DO alpha = 1, 3
          DO beta = 1, 3
             index2 = index2 + 1
             IF(alpha .eq. beta) exc_vxc(index2) = exc_vxc_buf*2.d0                           
          ENDDO
       ENDDO
RETURN
end subroutine exchange_str


subroutine symmetry_stress_rotij(index_a,jatom,ele_tensor)

use struct, only : rotij, rotloc, lattic, ndif
implicit none
    integer, intent(in)    :: index_a,jatom
    real*8 , intent(inout) :: ele_tensor(1:9) 
    real*8              :: buf_prod1(1:3,1:3),         &
                           buf_prod2(1:3,1:3), A1(3,3),&
                           buf_prod3(1:3,1:3), &
                           detinv,B1(3,3), final_stress(3,3), &
                           rotij_st(3,3,ndif), BR4(3,3)!, B1_rot(3,3), B1_ij(3,3)
    integer             :: alpha, beta, index2, ii
    
    rotij_st = rotij        
    
    IF(LATTIC(1:1) .EQ. 'H') THEN
      BR4(1,1)=SQRT(3.d0)/2.d0                                              
      BR4(1,2)=-.5d0     
      BR4(1,3)=0.0d0                                                      
      BR4(2,1)=0.0d0                                                     
      BR4(2,2)=1.0d0                                                      
      BR4(2,3)=0.0d0                                                      
      BR4(3,1)=0.0d0                                                      
      BR4(3,2)=0.0d0                                                      
      BR4(3,3)=1.d0  

!see PIA3 in lapw2 symmetrization_stress.f      
!       BR4(1,1)=2.D0/SQRT3*PIA3(1)                                        
!       BR4(1,2)=1.D0/SQRT3*PIA3(1)                                        
!       BR4(1,3)=0.0D0                                                    
!       BR4(2,1)=0.0D0                                                    
!       BR4(2,2)=PIA3(2)                                                   
!       BR4(2,3)=0.0D0                                                    
!       BR4(3,1)=0.0D0                                                    
!       BR4(3,2)=0.0D0                                                    
!       BR4(3,3)=PIA3(3)                                                         

    ELSEIF(LATTIC(1:1) .EQ. 'R') THEN
      BR4(1,1)=1/2.d0/sqrt(3.d0)
      BR4(1,2)=-1/2.d0                                                     
      BR4(1,3)=1/3.d0                                                      
      BR4(2,1)=1/2.d0/SQRT(3.d0)                                          
      BR4(2,2)=1*0.5d0                                                
      BR4(2,3)=1/3.d0                                                      
      BR4(3,1)=-1/SQRT(3.d0)                                         
      BR4(3,2)=0.d0                                                
      BR4(3,3)=1/3.d0            
    ENDIF



    IF( (LATTIC(1:1) .eq. 'H' ) .OR. (LATTIC(1:1) .eq. 'R' ) ) THEN   
      A1(:,:) = BR4(:,:)
      detinv = 1/(A1(1,1)*A1(2,2)*A1(3,3) - A1(1,1)*A1(2,3)*A1(3,2)&
              - A1(1,2)*A1(2,1)*A1(3,3) + A1(1,2)*A1(2,3)*A1(3,1)&
              + A1(1,3)*A1(2,1)*A1(3,2) - A1(1,3)*A1(2,2)*A1(3,1))
      IF(detinv .eq. 0.0) STOP 'inverse of sym_struct does not exist'    
      B1(1,1) = +detinv * (A1(2,2)*A1(3,3) - A1(2,3)*A1(3,2))
      B1(2,1) = -detinv * (A1(2,1)*A1(3,3) - A1(2,3)*A1(3,1))
      B1(3,1) = +detinv * (A1(2,1)*A1(3,2) - A1(2,2)*A1(3,1))
      B1(1,2) = -detinv * (A1(1,2)*A1(3,3) - A1(1,3)*A1(3,2))
      B1(2,2) = +detinv * (A1(1,1)*A1(3,3) - A1(1,3)*A1(3,1))
      B1(3,2) = -detinv * (A1(1,1)*A1(3,2) - A1(1,2)*A1(3,1))
      B1(1,3) = +detinv * (A1(1,2)*A1(2,3) - A1(1,3)*A1(2,2))
      B1(2,3) = -detinv * (A1(1,1)*A1(2,3) - A1(1,3)*A1(2,1))
      B1(3,3) = +detinv * (A1(1,1)*A1(2,2) - A1(1,2)*A1(2,1))
   
      buf_prod2(:,:) = matmul(B1(:,:),rotij(:,:,index_a))
      buf_prod1(:,:) = matmul(buf_prod2(:,:),A1(:,:))
      
      rotij_st(:,:,index_a) = buf_prod1(:,:)

    ENDIF
   
    index2 = 0
    DO alpha = 1, 3
       DO beta = 1,3
          index2 = index2 + 1
          buf_prod1(alpha,beta) = ele_tensor(index2)
       ENDDO
    ENDDO
    buf_prod3 = buf_prod1
    
    
!       buf_prod2(:,:) = matmul( transpose(rotloc(:,:,jatom)),buf_prod1(:,:) )
!       buf_prod1(:,:) = matmul( buf_prod2(:,:),rotloc(:,:,jatom) )
! 
!       buf_prod2(:,:) = matmul(transpose(rotij_st(:,:,index_a)),buf_prod1(:,:))
!       buf_prod1(:,:) = matmul(buf_prod2(:,:),rotij_st(:,:,index_a))
     
!-------Below modification is done because of problem in Rutile if A^inv T A------    
      A1(:,:) = rotloc(:,:,jatom)
      detinv = 1/(A1(1,1)*A1(2,2)*A1(3,3) - A1(1,1)*A1(2,3)*A1(3,2)&
              - A1(1,2)*A1(2,1)*A1(3,3) + A1(1,2)*A1(2,3)*A1(3,1)&
              + A1(1,3)*A1(2,1)*A1(3,2) - A1(1,3)*A1(2,2)*A1(3,1))
      IF(detinv .eq. 0.0) STOP 'inverse of rotloc does not exist -- EC stress'    
      B1(1,1) = +detinv * (A1(2,2)*A1(3,3) - A1(2,3)*A1(3,2))
      B1(2,1) = -detinv * (A1(2,1)*A1(3,3) - A1(2,3)*A1(3,1))
      B1(3,1) = +detinv * (A1(2,1)*A1(3,2) - A1(2,2)*A1(3,1))
      B1(1,2) = -detinv * (A1(1,2)*A1(3,3) - A1(1,3)*A1(3,2))
      B1(2,2) = +detinv * (A1(1,1)*A1(3,3) - A1(1,3)*A1(3,1))
      B1(3,2) = -detinv * (A1(1,1)*A1(3,2) - A1(1,2)*A1(3,1))
      B1(1,3) = +detinv * (A1(1,2)*A1(2,3) - A1(1,3)*A1(2,2))
      B1(2,3) = -detinv * (A1(1,1)*A1(2,3) - A1(1,3)*A1(2,1))
      B1(3,3) = +detinv * (A1(1,1)*A1(2,2) - A1(1,2)*A1(2,1))
    
    buf_prod2(:,:) = matmul( B1(:,:),buf_prod1(:,:) )
    buf_prod1(:,:) = matmul( buf_prod2(:,:),A1(:,:) )

      A1(:,:) = rotij_st(:,:,index_a)
      detinv = 1/(A1(1,1)*A1(2,2)*A1(3,3) - A1(1,1)*A1(2,3)*A1(3,2)&
              - A1(1,2)*A1(2,1)*A1(3,3) + A1(1,2)*A1(2,3)*A1(3,1)&
              + A1(1,3)*A1(2,1)*A1(3,2) - A1(1,3)*A1(2,2)*A1(3,1))
      IF(detinv .eq. 0.0) STOP 'inverse of rotij does not exist -- -- EC stress'    
      B1(1,1) = +detinv * (A1(2,2)*A1(3,3) - A1(2,3)*A1(3,2))
      B1(2,1) = -detinv * (A1(2,1)*A1(3,3) - A1(2,3)*A1(3,1))
      B1(3,1) = +detinv * (A1(2,1)*A1(3,2) - A1(2,2)*A1(3,1))
      B1(1,2) = -detinv * (A1(1,2)*A1(3,3) - A1(1,3)*A1(3,2))
      B1(2,2) = +detinv * (A1(1,1)*A1(3,3) - A1(1,3)*A1(3,1))
      B1(3,2) = -detinv * (A1(1,1)*A1(3,2) - A1(1,2)*A1(3,1))
      B1(1,3) = +detinv * (A1(1,2)*A1(2,3) - A1(1,3)*A1(2,2))
      B1(2,3) = -detinv * (A1(1,1)*A1(2,3) - A1(1,3)*A1(2,1))
      B1(3,3) = +detinv * (A1(1,1)*A1(2,2) - A1(1,2)*A1(2,1))    
    
    buf_prod2(:,:) = matmul(B1(:,:),buf_prod1(:,:))
    buf_prod1(:,:) = matmul(buf_prod2(:,:),A1(:,:))
    
!------------end the modification block--------------------------------    
    
!     WRITE(*,*)'-----ROTLOC--------',jatom
!     DO beta = 1, 3
!        WRITE(*,'(3f10.5)') (rotloc(alpha,beta,jatom),alpha=1,3)
!     ENDDO
!     WRITE(*,*)
! 
!     WRITE(*,*)'-----ROTIJ--------',index_a
!     DO beta = 1, 3
!        WRITE(*,'(3f10.5)') (rotij_st(alpha,beta,index_a),alpha=1,3)
!     ENDDO
         
!   print*
   index2 = 0
   DO alpha = 1, 3
      DO beta = 1,3
         index2 = index2+1
         ele_tensor(index2) = buf_prod1(alpha,beta)
!         print*, ele_tensor(index2), buf_prod3(alpha,beta) 
      ENDDO
   ENDDO
!   print*   
   
RETURN
end subroutine symmetry_stress_rotij

