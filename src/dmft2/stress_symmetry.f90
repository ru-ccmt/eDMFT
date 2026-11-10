! Stress-tensor rotation and crystal-symmetry averaging.
! Adapted from WIEN2k 24 SRC_lapw2; ndif is named natm in dmft2.

! Print a diagnostic decomposition before crystal symmetrization.  All input
! tensors have already been reduced over MPI ranks, but remain in the raw
! extensive form accumulated by the LAPW kernels.  Dividing here by VOL makes
! every row directly comparable to :STR_VAL in Ry/Bohr^3.  None of these
! quantities is fed back into the production stress.
subroutine print_dmft_stress_diagnostics(main,window,norm,response_overlap, &
     sigma_moment, &
     sigma_match,sigma_radial,sigma_angular,sigma_interstitial, &
     sigma_coeff_fd_mt, &
     sigma_coeff_analytic_mt, &
     paired_window_matching, &
     paired_norm_matching, &
     paired_main_correction_matching, &
     paired_window_complete, &
     paired_norm_complete, &
     paired_main_correction_complete, &
     paired_window_kinetic_matching, &
     paired_norm_kinetic_matching, &
     paired_main_correction_kinetic_matching, &
     paired_window_kinetic_complete, &
     paired_norm_kinetic_complete, &
     paired_main_correction_kinetic_complete, &
     paired_sigma_matching, &
     paired_sigma_complete, &
     paired_diagonal_reference, &
     density_diagonal, &
     density_offdiagonal,projector,rela_corr,rela_surf,vol)
  implicit none
  real*8,intent(in) :: main(3,3,4),window(3,3,4),norm(3,3,4)
  real*8,intent(in) :: response_overlap(3,3,3)
  real*8,intent(in) :: sigma_moment(3,3),projector(3,3)
  real*8,intent(in) :: sigma_match(3,3),sigma_radial(3,3)
  real*8,intent(in) :: sigma_angular(3,3)
  real*8,intent(in) :: sigma_interstitial(3,3)
  real*8,intent(in) :: sigma_coeff_fd_mt(3,3)
  real*8,intent(in) :: sigma_coeff_analytic_mt(3,3)
  real*8,intent(in) :: paired_window_matching(3,3)
  real*8,intent(in) :: paired_norm_matching(3,3)
  real*8,intent(in) :: paired_main_correction_matching(3,3)
  real*8,intent(in) :: paired_window_complete(3,3)
  real*8,intent(in) :: paired_norm_complete(3,3)
  real*8,intent(in) :: paired_main_correction_complete(3,3)
  real*8,intent(in) :: paired_window_kinetic_matching(3,3)
  real*8,intent(in) :: paired_norm_kinetic_matching(3,3)
  real*8,intent(in) :: paired_main_correction_kinetic_matching(3,3)
  real*8,intent(in) :: paired_window_kinetic_complete(3,3)
  real*8,intent(in) :: paired_norm_kinetic_complete(3,3)
  real*8,intent(in) :: paired_main_correction_kinetic_complete(3,3)
  real*8,intent(in) :: paired_sigma_matching(3,3)
  real*8,intent(in) :: paired_sigma_complete(3,3)
  real*8,intent(in) :: paired_diagonal_reference(3,3,4)
  real*8,intent(in) :: density_diagonal(3,3,4)
  real*8,intent(in) :: density_offdiagonal(3,3,4)
  real*8,intent(in) :: rela_corr(3,3),rela_surf(3,3),vol
  real*8 :: main_total(3,3),window_total(3,3),norm_total(3,3)
  real*8 :: dft_like(3,3),assembled(3,3)

  main_total=sum(main,dim=3)
  window_total=sum(window,dim=3)
  norm_total=sum(norm,dim=3)
  dft_like=main_total-sigma_moment
  assembled=main_total+window_total+norm_total+projector

  write(6,'(A)') '=== DMFT analytic-stress diagnostic decomposition ==='
  write(6,'(A)') 'Rows are MPI-reduced but unsymmetrized; units are Ry/Bohr^3.'
  call print_one('MAIN_TOTAL',main_total)
  call print_one('MAIN_DFTLIKE',dft_like)
  call print_one('SIGMA_MOMENT_OVERLAP',sigma_moment)
  call print_one('SIGMA_OVERLAP_MATCHING',sigma_match)
  call print_one('SIGMA_OVERLAP_RADIAL',sigma_radial)
  call print_one('SIGMA_OVERLAP_ANGULAR',sigma_angular)
  call print_one('SIGMA_OVERLAP_INTERSTITIAL',sigma_interstitial)
  call print_one('SIGMA_OVERLAP_FIXED_SPHERE', &
       sigma_match+sigma_interstitial)
  call print_one('SIGMA_OVERLAP_COEFF_FD_MT',sigma_coeff_fd_mt)
  call print_one('SIGMA_OVERLAP_COEFF_FD_TOTAL', &
       sigma_coeff_fd_mt+sigma_interstitial)
  call print_one('SIGMA_OVERLAP_COEFF_ANALYTIC_MT', &
       sigma_coeff_analytic_mt)
  call print_one('SIGMA_OVERLAP_COEFF_ANALYTIC_TOTAL', &
       sigma_coeff_analytic_mt+sigma_interstitial)
  call print_one('PAIRED_WINDOW_MATCHING',paired_window_matching)
  call print_one('PAIRED_NORM_MATCHING',paired_norm_matching)
  call print_one('PAIRED_MAIN_CORRECTION_MATCHING', &
       paired_main_correction_matching)
  call print_one('PAIRED_WINDOW_COMPLETE',paired_window_complete)
  call print_one('PAIRED_NORM_COMPLETE',paired_norm_complete)
  call print_one('PAIRED_MAIN_CORRECTION_COMPLETE', &
       paired_main_correction_complete)
  call print_one('PAIRED_WINDOW_KINETIC_MATCHING', &
       paired_window_kinetic_matching)
  call print_one('PAIRED_NORM_KINETIC_MATCHING', &
       paired_norm_kinetic_matching)
  call print_one('PAIRED_MAIN_CORRECTION_KINETIC_MATCHING', &
       paired_main_correction_kinetic_matching)
  call print_one('PAIRED_WINDOW_KINETIC_COMPLETE', &
       paired_window_kinetic_complete)
  call print_one('PAIRED_NORM_KINETIC_COMPLETE', &
       paired_norm_kinetic_complete)
  call print_one('PAIRED_MAIN_CORRECTION_KINETIC_COMPLETE', &
       paired_main_correction_kinetic_complete)
  call print_one('PAIRED_SIGMA_MATCHING',paired_sigma_matching)
  call print_one('PAIRED_SIGMA_COMPLETE',paired_sigma_complete)
  call print_one('PAIRED_DIAGONAL_REFERENCE', &
       sum(paired_diagonal_reference,dim=3))
  call print_one('PAIRED_RECONSTRUCTED_MAIN', &
       sum(paired_diagonal_reference,dim=3)+ &
       paired_main_correction_complete)
  call print_one('PAIRED_RECONSTRUCTED_MAIN_KINETIC', &
       sum(paired_diagonal_reference,dim=3)+ &
       paired_main_correction_kinetic_complete)
  call print_one('DFTLIKE_DIAGONAL',sum(density_diagonal,dim=3))
  call print_one('DFTLIKE_OFFDIAGONAL',sum(density_offdiagonal,dim=3))
  call print_one('WINDOW_TOTAL',window_total)
  call print_one('NORM_WINDOW_TOTAL',norm_total)
  call print_one('MAIN_OVERLAP_PART',response_overlap(:,:,1))
  call print_one('WINDOW_OVERLAP_PART',response_overlap(:,:,2))
  call print_one('WINDOW_HAMILTONIAN_PART', &
       window_total-response_overlap(:,:,2))
  call print_one('NORM_OVERLAP_PART',response_overlap(:,:,3))
  call print_one('NORM_HAMILTONIAN_PART', &
       norm_total-response_overlap(:,:,3))
  call print_one('EXPLICIT_PROJECTOR',projector)
  ! These two tensors are included in MAIN_KIN_MT and are repeated here
  ! only to make the CALC=RELA audit explicit.
  call print_one('RELA_VOLUME_IN_KIN_MT',rela_corr)
  call print_one('RELA_SURFACE_IN_KIN_MT',rela_surf)
  call print_one('ASSEMBLED_VALENCE',assembled)

  call print_one('MAIN_VAL_CORR',main(:,:,1))
  call print_one('MAIN_KIN_MT',main(:,:,2))
  call print_one('MAIN_KIN_INT',main(:,:,3))
  call print_one('MAIN_APW',main(:,:,4))
  call print_one('WINDOW_VAL_CORR',window(:,:,1))
  call print_one('WINDOW_KIN_MT',window(:,:,2))
  call print_one('WINDOW_KIN_INT',window(:,:,3))
  call print_one('WINDOW_APW',window(:,:,4))
  call print_one('NORM_VAL_CORR',norm(:,:,1))
  call print_one('NORM_KIN_MT',norm(:,:,2))
  call print_one('NORM_KIN_INT',norm(:,:,3))
  call print_one('NORM_APW',norm(:,:,4))

contains

  subroutine print_one(label,tensor)
    character(len=*),intent(in) :: label
    real*8,intent(in) :: tensor(3,3)
    character(len=96) :: tag
    real*8,parameter :: rybohr3_to_gpa=14710.5164057d0
    integer :: alpha
    do alpha=1,3
       write(tag,'(A,A,A,I3.3)') ':DMFT_STRESS_DIAG_',trim(label), &
            '_RY',alpha
       write(6,'(A,3(1x,ES20.11))') trim(tag),tensor(alpha,:)/vol
    enddo
    write(tag,'(A,A,A)') ':DMFT_STRESS_DIAG_',trim(label),'_MEAN_GPA'
    write(6,'(A,1x,ES20.11)') trim(tag), &
         (tensor(1,1)+tensor(2,2)+tensor(3,3))* &
         rybohr3_to_gpa/(3.d0*vol)
  end subroutine print_one

end subroutine print_dmft_stress_diagnostics


! Selectively remove DMFT projector-response pieces from the physical stress
! after their MPI reduction.  This routine is an audit mechanism, not a
! silent empirical correction: FULL is the default and reproduces the normal
! assembly bit for bit.  BASE is useful for testing whether the response
! pieces are already represented by the stationary A/AE-weighted LAPW
! contractions.  The four one-at-a-time modes localize any disagreement.
subroutine apply_dmft_stress_response_mode(val_corr,kin_mt,apw,kin_int, &
     projector,window,norm,sigma_moment)
  implicit none
  real*8,intent(inout) :: val_corr(3,3),kin_mt(3,3),apw(3,3)
  real*8,intent(inout) :: kin_int(3,3),projector(3,3)
  real*8,intent(in) :: window(3,3,4),norm(3,3,4)
  real*8,intent(in) :: sigma_moment(3,3)
  character(len=32) :: mode
  integer :: status,i,letter
  logical :: omit_sigma,omit_window,omit_norm,omit_projector

  mode='FULL'
  call get_environment_variable('EDMFT_STRESS_RESPONSE_MODE',mode, &
       status=status)
  if (status.ne.0 .or. len_trim(mode).eq.0) mode='FULL'
  mode=adjustl(mode)
  do i=1,len_trim(mode)
     letter=iachar(mode(i:i))
     if (letter.ge.iachar('a') .and. letter.le.iachar('z')) &
          mode(i:i)=achar(letter-iachar('a')+iachar('A'))
  enddo

  omit_sigma=.false.
  omit_window=.false.
  omit_norm=.false.
  omit_projector=.false.
  select case(trim(mode))
  case('FULL')
     continue
  case('BASE')
     omit_sigma=.true.
     omit_window=.true.
     omit_norm=.true.
     omit_projector=.true.
  case('NO_SIGMA')
     omit_sigma=.true.
  case('NO_WINDOW')
     omit_window=.true.
  case('NO_NORM')
     omit_norm=.true.
  case('NO_PROJECTOR')
     omit_projector=.true.
  case default
     write(6,'(A,A)') 'ERROR: unknown EDMFT_STRESS_RESPONSE_MODE=', &
          trim(mode)
     stop 'invalid EDMFT stress response mode'
  end select

  if (omit_sigma) val_corr=val_corr-sigma_moment
  if (omit_window) then
     val_corr=val_corr-window(:,:,1)
     kin_mt=kin_mt-window(:,:,2)
     kin_int=kin_int-window(:,:,3)
     apw=apw-window(:,:,4)
  endif
  if (omit_norm) then
     val_corr=val_corr-norm(:,:,1)
     kin_mt=kin_mt-norm(:,:,2)
     kin_int=kin_int-norm(:,:,3)
     apw=apw-norm(:,:,4)
  endif
  if (omit_projector) projector=0.d0

  write(6,'(A,A)') ':DMFT_STRESS_RESPONSE_MODE: ',trim(mode)
end subroutine apply_dmft_stress_response_mode


! Apply a reference-only additive correction assembled outside dmft2 from
! complete LAPW H/O matrix derivatives and, when stable, a joint derivative
! of the normalized DMFT projector.  This is intentionally not a fitted or
! unconditional correction: ordinary calculations never enter this routine
! unless EDMFT_STRESS_REFERENCE_CORRECTION names a valid file.
!
! File format (extensive, before symmetry and division by the cell volume):
!
!   EDMFT_STRESS_REFERENCE_CORRECTION_V1
!   RAW_RY
!   dT11 dT12 dT13
!   dT21 dT22 dT23
!   dT31 dT32 dT33
!
! The correction is accumulated in PROJECTOR only as a convenient additive
! carrier.  It may contain both the replacement of the incomplete M_Sigma*dO
! kernel and the replacement of the split projector/window/normalization
! response.  Adding it after apply_dmft_stress_response_mode prevents an
! audit selector from silently removing part of the supplied reference.
subroutine apply_dmft_stress_reference_correction(projector,vol)
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  real*8,intent(inout) :: projector(3,3)
  real*8,intent(in) :: vol
  character(len=1024) :: filename
  character(len=64) :: header,units
  real*8 :: correction(3,3),mean_gpa
  real*8,parameter :: rybohr3_to_gpa=14710.5164057d0
  integer :: status,fh,ios,row

  filename=''
  call get_environment_variable('EDMFT_STRESS_REFERENCE_CORRECTION', &
       filename,status=status)
  if (status.ne.0 .or. len_trim(filename).eq.0) return

  open(newunit=fh,file=trim(filename),status='old',action='read',iostat=ios)
  if (ios.ne.0) then
     write(6,'(A,A)') 'ERROR: cannot open EDMFT stress reference correction ', &
          trim(filename)
     stop 'invalid EDMFT stress reference correction'
  endif
  read(fh,'(A)',iostat=ios) header
  if (ios.ne.0 .or. trim(header).ne. &
       'EDMFT_STRESS_REFERENCE_CORRECTION_V1') then
     close(fh)
     write(6,'(A,A)') 'ERROR: invalid reference-correction header in ', &
          trim(filename)
     stop 'invalid EDMFT stress reference correction'
  endif
  read(fh,'(A)',iostat=ios) units
  if (ios.ne.0 .or. trim(units).ne.'RAW_RY') then
     close(fh)
     write(6,'(A,A)') 'ERROR: reference correction must use RAW_RY in ', &
          trim(filename)
     stop 'invalid EDMFT stress reference units'
  endif
  do row=1,3
     read(fh,*,iostat=ios) correction(row,:)
     if (ios.ne.0) then
        close(fh)
        write(6,'(A,A)') 'ERROR: truncated reference correction in ', &
             trim(filename)
        stop 'invalid EDMFT stress reference tensor'
     endif
  enddo
  close(fh)
  if (.not.all(ieee_is_finite(correction))) then
     write(6,'(A,A)') 'ERROR: non-finite reference correction in ', &
          trim(filename)
     stop 'invalid EDMFT stress reference tensor'
  endif
  if (maxval(abs(correction-transpose(correction))).gt.1.d-12* &
       max(1.d0,maxval(abs(correction)))) then
     write(6,'(A,A)') 'ERROR: nonsymmetric reference correction in ', &
          trim(filename)
     stop 'invalid EDMFT stress reference tensor'
  endif

  projector=projector+correction
  mean_gpa=(correction(1,1)+correction(2,2)+correction(3,3))* &
       rybohr3_to_gpa/(3.d0*vol)
  write(6,'(A,1x,A)') ':DMFT_STRESS_REFERENCE_CORRECTION_FILE', &
       trim(filename)
  do row=1,3
     write(6,'(A,I3.3,A,3(1x,ES20.11))') &
          ':DMFT_STRESS_REFERENCE_CORRECTION_RAW_RY',row,':',correction(row,:)
  enddo
  write(6,'(A,1x,ES20.11)') &
       ':DMFT_STRESS_REFERENCE_CORRECTION_MEAN_GPA',mean_gpa
end subroutine apply_dmft_stress_reference_correction

! subroutine symmetrization_stress(in_stress)
subroutine symmetrization_stress(val_corr,kin_rmt,kin_rmt_apw, &
     kin_int,projector,norm_window,nsp_match,nsp_radial,nsp_angular)
! subroutine symmetrization_stress(val_corr,kin_rmt,kin_rmt_apw,kin_int,rela_corr, rel)
    use sym2
    use com  , only : rel
    use structure, only : lattic,br2,aa,bb,cc, vol    
    use dmft_stress, only: stress_g2_mode
    implicit none
    integer               :: i,j,k,alpha,beta, index2
!     real*8, intent(inout) :: in_stress(1:3,1:3)
    real*8              :: val_corr(1:3,1:3), kin_rmt(1:3,1:3), kin_int(1:3,1:3),  &
                           kin_rmt_apw(1:3,1:3),projector(1:3,1:3), &
                           norm_window(1:3,1:3), &
                           nsp_match(1:3,1:3),nsp_radial(1:3,1:3), &
                           nsp_angular(1:3,1:3), &
                           rela_corr(1:3,1:3)
    
    real*8              :: buf_prod1(3,3,iord),buf_prod2(3,3,iord), buf_out(1:3,1:3)
    
    real*8              :: val_corr_fn(3,3), kin_rmt_fn(3,3), kin_int_fn(3,3),       &
                           kin_rmt_apw_fn(3,3),projector_fn(3,3), &
                           norm_window_fn(3,3), &
                           nsp_match_fn(3,3),nsp_radial_fn(3,3), &
                           nsp_angular_fn(3,3), &
                           rela_corr_fn(3,3), tot_val_stress(1:9),     & 
                           con_fac, BR4(3,3), sym_struct(3,3,iord), A1(3,3), B1(3,3), &
                           detinv, mat_prod1(3,3), mat_prod2(3,3),fin_val_str(3,3), SQRT3, &
                           PIA3(3), PI
                         
        !BR4 is same as BR2 of SRC_symmetry/dirlat.f                 
    
    PI=     3.1415926535897932d0
    SQRT3 = SQRT(3.D0)
    PIA3(1) = 2.D0*PI/AA
    PIA3(2) = 2.D0*PI/BB
    PIA3(3) = 2.D0*PI/CC
        
    buf_out(:,:) = 0.d0
    
    val_corr_fn(:,:)    = 0.d0
    kin_rmt_fn(:,:)     = 0.d0
    kin_int_fn(:,:)     = 0.d0
    kin_rmt_apw_fn(:,:) = 0.d0
    projector_fn(:,:)   = 0.d0
    norm_window_fn(:,:) = 0.d0
    nsp_match_fn(:,:) = 0.d0
    nsp_radial_fn(:,:) = 0.d0
    nsp_angular_fn(:,:) = 0.d0
    fin_val_str(:,:)    = 0.d0
    rela_corr_fn(:,:)   = 0.d0
    con_fac          = 147105.164057
            
    IF(LATTIC(1:1) .EQ. 'H') THEN
      BR4(1,1)=SQRT(3.d0)/2.d0                                              
      BR4(1,2)=-.5d0     !-0.5d0 is original                                                 
      BR4(1,3)=0.0d0                                                      
      BR4(2,1)=0.0d0                                                     
      BR4(2,2)=1.0d0                                                      
      BR4(2,3)=0.0d0                                                      
      BR4(3,1)=0.0d0                                                      
      BR4(3,2)=0.0d0                                                      
      BR4(3,3)=1.d0  
            
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
    
    sym_struct = iz          !iz integer and sym_struct is real
    
    IF( (LATTIC(1:1) .eq. 'H' ) .OR. (LATTIC(1:1) .eq. 'R' ) ) THEN
       A1(:,:) = BR4(:,:)
       DO i = 1, iord                    
           
           detinv = 1/(A1(1,1)*A1(2,2)*A1(3,3) - A1(1,1)*A1(2,3)*A1(3,2)&
              - A1(1,2)*A1(2,1)*A1(3,3) + A1(1,2)*A1(2,3)*A1(3,1)&
              + A1(1,3)*A1(2,1)*A1(3,2) - A1(1,3)*A1(2,2)*A1(3,1))
           
           IF(detinv .eq. 0.0) STOP 'inverse of sym_struct does not exist'
           ! Calculate the inverse of the matrix
           B1(1,1) = +detinv * (A1(2,2)*A1(3,3) - A1(2,3)*A1(3,2))
           B1(2,1) = -detinv * (A1(2,1)*A1(3,3) - A1(2,3)*A1(3,1))
           B1(3,1) = +detinv * (A1(2,1)*A1(3,2) - A1(2,2)*A1(3,1))
           B1(1,2) = -detinv * (A1(1,2)*A1(3,3) - A1(1,3)*A1(3,2))
           B1(2,2) = +detinv * (A1(1,1)*A1(3,3) - A1(1,3)*A1(3,1))
           B1(3,2) = -detinv * (A1(1,1)*A1(3,2) - A1(1,2)*A1(3,1))
           B1(1,3) = +detinv * (A1(1,2)*A1(2,3) - A1(1,3)*A1(2,2))
           B1(2,3) = -detinv * (A1(1,1)*A1(2,3) - A1(1,3)*A1(2,1))
           B1(3,3) = +detinv * (A1(1,1)*A1(2,2) - A1(1,2)*A1(2,1))
           
                                                     
            mat_prod1(:,:) = matmul( B1(:,:),iz(:,:,i) )            
            mat_prod2(:,:) = matmul( mat_prod1(:,:),A1(:,:) )
            
            sym_struct(:,:,i) = mat_prod2(:,:)            
            
               
           
        ENDDO
    ENDIF
    
!write(6,*) 'vc',val_corr(1,1:3)  
!write(6,*) val_corr(2,1:3)  
!write(6,*) val_corr(3,1:3)  
    
    DO i = 1, iord
       !=======calculating the inverse of all szmmetry operations presented in the structure file=======
          A1(:,:) = sym_struct(:,:,i)
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
        
       !=======any 3 by 3 stress tensor
!        buf_prod1(:,:,i) = matmul(B1(:,:),in_stress(:,:) )                           !for A^-1*Tensor*A
!        buf_prod2(:,:,i) = matmul(buf_prod1(:,:,i),sym_struct(:,:,i))           
!        buf_out(:,:) = buf_out(:,:) + buf_prod2(:,:,i)
       
       !=======for valence correction stress       
!      buf_prod1(:,:,i) = matmul(transpose(sym_struct(:,:,i)),val_corr(:,:) )    !for A^T*Tensor*A
       buf_prod1(:,:,i) = matmul(B1(:,:),val_corr(:,:) )                           !for A^-1*Tensor*A
       buf_prod2(:,:,i) = matmul(buf_prod1(:,:,i),sym_struct(:,:,i))
              
       val_corr_fn(:,:) = val_corr_fn(:,:) + buf_prod2(:,:,i)
       
       !=======for kinetic energy stress inside an atomic sphere             
!      buf_prod1(:,:,i) = matmul(transpose(sym_struct(:,:,i)),kin_rmt(:,:) )
         buf_prod1(:,:,i) = matmul(B1(:,:),kin_rmt(:,:) )
       buf_prod2(:,:,i) = matmul(buf_prod1(:,:,i),sym_struct(:,:,i))
                     
       kin_rmt_fn(:,:) = kin_rmt_fn(:,:) + buf_prod2(:,:,i)
       
       !=======kinetic energy stress in the interstitial region                
!      buf_prod1(:,:,i) = matmul(transpose(sym_struct(:,:,i)),kin_int(:,:) )
         buf_prod1(:,:,i) = matmul(B1(:,:),kin_int(:,:) )
       buf_prod2(:,:,i) = matmul(buf_prod1(:,:,i),sym_struct(:,:,i))
       
       kin_int_fn(:,:) = kin_int_fn(:,:) + buf_prod2(:,:,i)
       
       !=======extra surface term from the kinetic energy term============       
!      buf_prod1(:,:,i) = matmul(transpose(sym_struct(:,:,i)),kin_rmt_apw(:,:) )
         buf_prod1(:,:,i) = matmul(B1(:,:),kin_rmt_apw(:,:) )
       buf_prod2(:,:,i) = matmul(buf_prod1(:,:,i),sym_struct(:,:,i))
                     
       kin_rmt_apw_fn(:,:) = kin_rmt_apw_fn(:,:) + buf_prod2(:,:,i)

       ! Explicit derivative of the normalized DMFT projector.
       buf_prod1(:,:,i)=matmul(B1(:,:),projector(:,:))
       buf_prod2(:,:,i)=matmul(buf_prod1(:,:,i),sym_struct(:,:,i))
       projector_fn(:,:)=projector_fn(:,:)+buf_prod2(:,:,i)

       ! Diagnostic copy of the normalization/window cross response.  The
       ! same contribution has already been accumulated into the four LAPW
       ! tensors above, so it must not be added to fin_val_str a second time.
       buf_prod1(:,:,i)=matmul(B1(:,:),norm_window(:,:))
       buf_prod2(:,:,i)=matmul(buf_prod1(:,:,i),sym_struct(:,:,i))
       norm_window_fn(:,:)=norm_window_fn(:,:)+buf_prod2(:,:,i)

       ! Diagnostic decomposition of the nonspherical Eq. (G2) block.
       ! The selected combination is already present exactly once in
       ! VAL_CORR; these three tensors are symmetrized only for reporting.
       buf_prod1(:,:,i)=matmul(B1(:,:),nsp_match(:,:))
       buf_prod2(:,:,i)=matmul(buf_prod1(:,:,i),sym_struct(:,:,i))
       nsp_match_fn(:,:)=nsp_match_fn(:,:)+buf_prod2(:,:,i)
       buf_prod1(:,:,i)=matmul(B1(:,:),nsp_radial(:,:))
       buf_prod2(:,:,i)=matmul(buf_prod1(:,:,i),sym_struct(:,:,i))
       nsp_radial_fn(:,:)=nsp_radial_fn(:,:)+buf_prod2(:,:,i)
       buf_prod1(:,:,i)=matmul(B1(:,:),nsp_angular(:,:))
       buf_prod2(:,:,i)=matmul(buf_prod1(:,:,i),sym_struct(:,:,i))
       nsp_angular_fn(:,:)=nsp_angular_fn(:,:)+buf_prod2(:,:,i)
       
       !===============Relativistic correction from kinetic energy operator=======
!        buf_prod1(:,:,i) = matmul(transpose(sym_struct(:,:,i)),rela_corr(:,:) )
         buf_prod1(:,:,i) = matmul(B1(:,:),rela_corr(:,:) )
         buf_prod2(:,:,i) = matmul(buf_prod1(:,:,i),sym_struct(:,:,i))
       
       rela_corr_fn(:,:) = rela_corr_fn(:,:) + buf_prod2(:,:,i)
    ENDDO
            
        val_corr_fn(:,:) = val_corr_fn(:,:)/iord
                      ! val_corr_fn(:,:) = val_corr_fn(:,:)*147105.164057/iord  !if pressure need to print out
        kin_rmt_fn(:,:) = kin_rmt_fn(:,:)/iord
        kin_int_fn(:,:) = kin_int_fn(:,:)/iord
        kin_rmt_apw_fn(:,:) = kin_rmt_apw_fn(:,:)/iord
        projector_fn(:,:) = projector_fn(:,:)/iord
        norm_window_fn(:,:)= norm_window_fn(:,:)/iord
        nsp_match_fn(:,:)=nsp_match_fn(:,:)/iord
        nsp_radial_fn(:,:)=nsp_radial_fn(:,:)/iord
        nsp_angular_fn(:,:)=nsp_angular_fn(:,:)/iord
        rela_corr_fn(:,:) = rela_corr_fn(:,:)/iord 
        
!         in_stress = buf_out / iord
        
       WRITE(6,*) '==============INDIVIDUAL VALENCE STRESS TNESOR===================='
!           
      WRITE(6,*)    
      WRITE(6,*)'---SPHERICAL PART OF THE "valence correction" ---'
      DO alpha = 1,3
        ! WRITE(21,112) ( val_corr_fn(alpha,beta), beta = 1,3 )
         WRITE(6,*) ( val_corr_fn(alpha,beta), beta = 1,3 )
      ENDDO
!       
        WRITE(6,*)    
        WRITE(6,*)'---Valence kinetic stress inside atomic sphere ---'
        DO alpha = 1,3
          ! WRITE(21,112) (kin_rmt_fn(alpha,beta),beta = 1,3)
           WRITE(6,*) (kin_rmt_fn(alpha,beta),beta = 1,3)
        ENDDO
        
        WRITE(6,*)
        WRITE(6,*)'---Valence kinetic stress in the intersititial region ---'
        DO alpha = 1,3
          ! WRITE(21,112) (kin_int_fn(alpha,beta),beta = 1,3)
           WRITE(6,*) (kin_int_fn(alpha,beta),beta = 1,3)
        ENDDO
        
        WRITE(6,*)
        WRITE(6,*)'---Extra surface term from the kinetic energy ---'
        DO alpha = 1,3        
           WRITE(6,*) (kin_rmt_apw_fn(alpha,beta),beta = 1,3)
        ENDDO
        WRITE(6,*)
        WRITE(6,*)'---Explicit normalized-DMFT-projector stress ---'
        DO alpha=1,3
           WRITE(6,*) (projector_fn(alpha,beta),beta=1,3)
        ENDDO
        WRITE(6,*)
        WRITE(6,*)'---Normalization-window response stress ---'
        DO alpha=1,3
           WRITE(6,*) (norm_window_fn(alpha,beta),beta=1,3)
        ENDDO
        WRITE(6,*)
        WRITE(6,*) '---Nonspherical G2 terms; selected mode: ', &
             trim(stress_g2_mode),' ---'
        DO alpha=1,3
           WRITE(6,*) 'MATCH  ',(nsp_match_fn(alpha,beta),beta=1,3)
           WRITE(6,*) 'RADIAL ',(nsp_radial_fn(alpha,beta),beta=1,3)
           WRITE(6,*) 'ANGULAR',(nsp_angular_fn(alpha,beta),beta=1,3)
           WRITE(21,2231) ':DMFT_G2_MATCH_RY',alpha, &
                (nsp_match_fn(alpha,beta)/vol,beta=1,3)
           WRITE(21,2231) ':DMFT_G2_RADIAL_RY',alpha, &
                (nsp_radial_fn(alpha,beta)/vol,beta=1,3)
           WRITE(21,2231) ':DMFT_G2_ANGULAR_RY',alpha, &
                (nsp_angular_fn(alpha,beta)/vol,beta=1,3)
        ENDDO
        WRITE(6,*)
! ! 
!       IF(REL) THEN
!         WRITE(21,*)    
!         WRITE(21,*)'---Extra correction in the case of Relativity ---'
!         DO alpha = 1,3        
!            WRITE(21,*) (rela_corr_fn(alpha,beta),beta = 1,3)
!         ENDDO
!         WRITE(21,*)
!          WRITE(21,'(A,1x,ES19.12)') ':REL_CORR_001',rela_corr_fn(1,1)
!          WRITE(21,'(A,1x,ES19.12)') ':REL_CORR_002',rela_corr_fn(1,2)
!          WRITE(21,'(A,1x,ES19.12)') ':REL_CORR_003',rela_corr_fn(1,3)
!          WRITE(21,'(A,1x,ES19.12)') ':REL_CORR_004',rela_corr_fn(2,1)
!          WRITE(21,'(A,1x,ES19.12)') ':REL_CORR_005',rela_corr_fn(2,2)
!          WRITE(21,'(A,1x,ES19.12)') ':REL_CORR_006',rela_corr_fn(2,3)
!          WRITE(21,'(A,1x,ES19.12)') ':REL_CORR_007',rela_corr_fn(3,1)
!          WRITE(21,'(A,1x,ES19.12)') ':REL_CORR_008',rela_corr_fn(3,2)
!          WRITE(21,'(A,1x,ES19.12)') ':REL_CORR_009',rela_corr_fn(3,3)
!  
!       ENDIF
                                      
!       index2 = 0
!       DO alpha = 1, 3
!          DO beta = 1, 3
!             index2 = index2 + 1
!             tot_val_stress(index2) = val_corr_fn(alpha,beta)+kin_rmt_fn(alpha,beta)+ &
!                                  kin_int_fn(alpha,beta)+kin_rmt_apw_fn(alpha,beta)
!          ENDDO
!       ENDDO

      WRITE(21,*)
      WRITE(21,*) '==============TOTAL STRESS FOR VALENCES ELECTRONS================='
      fin_val_str(:,:) = val_corr_fn(:,:)+kin_rmt_fn(:,:)+ &
                          kin_int_fn(:,:)+kin_rmt_apw_fn(:,:)+ &
                          projector_fn(:,:)
      fin_val_str = fin_val_str/vol                    
      DO alpha = 1, 3
         WRITE(21,2225) alpha, (fin_val_str(alpha,beta),beta=1,3)
      ENDDO
      WRITE(21,*)      
!       WRITE(21,'(A,1x,ES19.12)') ':TOT_VAL_001',fin_val_str(1,1)
!       WRITE(21,'(A,1x,ES19.12)') ':TOT_VAL_002',fin_val_str(1,2)
!       WRITE(21,'(A,1x,ES19.12)') ':TOT_VAL_003',fin_val_str(1,3)      
!       WRITE(21,'(A,1x,ES19.12)') ':TOT_VAL_004',fin_val_str(2,1)
!       WRITE(21,'(A,1x,ES19.12)') ':TOT_VAL_005',fin_val_str(2,2)
!       WRITE(21,'(A,1x,ES19.12)') ':TOT_VAL_006',fin_val_str(2,3)      
!       WRITE(21,'(A,1x,ES19.12)') ':TOT_VAL_007',fin_val_str(3,1)
!       WRITE(21,'(A,1x,ES19.12)') ':TOT_VAL_008',fin_val_str(3,2)
!       WRITE(21,'(A,1x,ES19.12)') ':TOT_VAL_009',fin_val_str(3,3)
!       
! !       DO alpha = 1,3
! !          DO beta = 1, 3
! !             WRITE(21,*) (fin_val_str(alpha,beta),beta=1,3)
! !          ENDDO   
! !       ENDDO
! !       WRITE(21,*)

         
111   FORMAT(3(3I2,/))
112   FORMAT(3(3ES19.12,/))
2225  FORMAT(':STR_VAL',i3.3,':',2x,3f20.12)
2231  FORMAT(A,i3.3,':',2x,3es20.12)

return
end subroutine symmetrization_stress


subroutine symmetry_stress_rotij(index_a,jatom,buf_sum1p)

use structure, only : rotij, rotloc, lattic, natm
implicit none
    integer, intent(in)    :: index_a,jatom
    real*8 , intent(inout) :: buf_sum1p(1:3,1:3)
    
    real*8              :: buf_prod1(1:3,1:3),               &
                           buf_prod2(1:3,1:3), A1(3,3),      &                           
                           detinv,B1(3,3), final_stress(3,3), &
                           rotij_st(3,3,natm), BR4(3,3)!, B1_rot(3,3), B1_ij(3,3)
    integer             :: alpha, beta, ii, index2
    
    rotij_st = rotij        
    
    IF(LATTIC(1:1) .EQ. 'H') THEN
      BR4(1,1)=SQRT(3.d0)/2.d0                                              
      BR4(1,2)=-.5d0     !-0.5d0 is original                                                 
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
   
!     index2 = 0
!     DO alpha = 1, 3
!        DO beta = 1,3
!           index2 = index2 + 1
!           buf_prod1(alpha,beta) = ele_tensor(index2)
!        ENDDO
!     ENDDO

      buf_prod1 = buf_sum1p
                    
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
    
    buf_sum1p = buf_prod1
!    index2 = 0
!    DO alpha = 1, 3
!       DO beta = 1,3
!          index2 = index2+1
!          ele_tensor(index2) = buf_prod1(alpha,beta)
! !         print*, ele_tensor(index2), buf_prod3(alpha,beta) 
!       ENDDO
!    ENDDO
!   print*   
   
RETURN
end subroutine symmetry_stress_rotij


