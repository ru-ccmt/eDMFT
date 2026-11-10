subroutine Reduce_stress(val_corr,kin_mt,apw,kin_int,projector,norm_window, &
     nsp_match,nsp_radial,nsp_angular,rela_corr,rela_surf, &
     diag_main,diag_window,diag_norm,diag_response_overlap, &
     diag_sigma,diag_sigma_match, &
     diag_sigma_radial,diag_sigma_angular,diag_sigma_interstitial, &
     diag_sigma_coeff_fd_mt, &
     diag_sigma_coeff_analytic_mt, &
     diag_paired_window_matching, &
     diag_paired_norm_matching, &
     diag_paired_main_correction_matching, &
     diag_paired_window_complete, &
     diag_paired_norm_complete, &
     diag_paired_main_correction_complete, &
     diag_paired_window_kinetic_matching, &
     diag_paired_norm_kinetic_matching, &
     diag_paired_main_correction_kinetic_matching, &
     diag_paired_window_kinetic_complete, &
     diag_paired_norm_kinetic_complete, &
     diag_paired_main_correction_kinetic_complete, &
     diag_paired_sigma_matching, &
     diag_paired_sigma_complete, &
     diag_paired_diagonal_reference, &
     diag_density_diagonal, &
     diag_density_offdiagonal,trace_rhoE,trace_erho, &
     trace_sigma,radial_hsph,radial_hnsp,radial_overlap,radial_by_pass, &
     profile_times)
  use com_mpi, only: myrank,master,ierr
  implicit none
  real*8,intent(inout) :: val_corr(3,3),kin_mt(3,3),apw(3,3), &
       kin_int(3,3),projector(3,3),norm_window(3,3), &
       nsp_match(3,3),nsp_radial(3,3),nsp_angular(3,3), &
       rela_corr(3,3),rela_surf(3,3), &
       diag_main(3,3,4),diag_window(3,3,4),diag_norm(3,3,4), &
       diag_response_overlap(3,3,3), &
       diag_sigma(3,3),diag_sigma_match(3,3), &
       diag_sigma_radial(3,3),diag_sigma_angular(3,3), &
       diag_sigma_interstitial(3,3), &
       diag_sigma_coeff_fd_mt(3,3), &
       diag_sigma_coeff_analytic_mt(3,3), &
       diag_paired_window_matching(3,3), &
       diag_paired_norm_matching(3,3), &
       diag_paired_main_correction_matching(3,3), &
       diag_paired_window_complete(3,3), &
       diag_paired_norm_complete(3,3), &
       diag_paired_main_correction_complete(3,3), &
       diag_paired_window_kinetic_matching(3,3), &
       diag_paired_norm_kinetic_matching(3,3), &
       diag_paired_main_correction_kinetic_matching(3,3), &
       diag_paired_window_kinetic_complete(3,3), &
       diag_paired_norm_kinetic_complete(3,3), &
       diag_paired_main_correction_kinetic_complete(3,3), &
       diag_paired_sigma_matching(3,3), &
       diag_paired_sigma_complete(3,3), &
       diag_paired_diagonal_reference(3,3,4), &
       diag_density_diagonal(3,3,4), &
       diag_density_offdiagonal(3,3,4),trace_rhoE,trace_erho, &
       trace_sigma,radial_hsph,radial_hnsp,radial_overlap,profile_times(8)
  real*8,intent(inout) :: radial_by_pass(3,3)
#ifdef _MPI
  include 'mpif.h'
  real*8 :: work(3,3),work3(3,3,3),work4(3,3,4),send_scalar(1),work_scalar(1), &
       work_profile(8)
  call MPI_REDUCE(val_corr,work,9,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_WORLD,ierr)
  if (myrank.eq.master) val_corr=work
  call MPI_REDUCE(kin_mt,work,9,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_WORLD,ierr)
  if (myrank.eq.master) kin_mt=work
  call MPI_REDUCE(apw,work,9,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_WORLD,ierr)
  if (myrank.eq.master) apw=work
  call MPI_REDUCE(kin_int,work,9,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_WORLD,ierr)
  if (myrank.eq.master) kin_int=work
  call MPI_REDUCE(projector,work,9,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_WORLD,ierr)
  if (myrank.eq.master) projector=work
  call MPI_REDUCE(norm_window,work,9,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_WORLD,ierr)
  if (myrank.eq.master) norm_window=work
  call MPI_REDUCE(nsp_match,work,9,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_WORLD,ierr)
  if (myrank.eq.master) nsp_match=work
  call MPI_REDUCE(nsp_radial,work,9,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_WORLD,ierr)
  if (myrank.eq.master) nsp_radial=work
  call MPI_REDUCE(nsp_angular,work,9,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_WORLD,ierr)
  if (myrank.eq.master) nsp_angular=work
  call MPI_REDUCE(rela_corr,work,9,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_WORLD,ierr)
  if (myrank.eq.master) rela_corr=work
  call MPI_REDUCE(rela_surf,work,9,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_WORLD,ierr)
  if (myrank.eq.master) rela_surf=work
  call MPI_REDUCE(diag_main,work4,36,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_WORLD,ierr)
  if (myrank.eq.master) diag_main=work4
  call MPI_REDUCE(diag_window,work4,36,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_WORLD,ierr)
  if (myrank.eq.master) diag_window=work4
  call MPI_REDUCE(diag_norm,work4,36,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_WORLD,ierr)
  if (myrank.eq.master) diag_norm=work4
  call MPI_REDUCE(diag_response_overlap,work3,27,MPI_DOUBLE_PRECISION, &
       MPI_SUM,master,MPI_COMM_WORLD,ierr)
  if (myrank.eq.master) diag_response_overlap=work3
  call MPI_REDUCE(diag_sigma,work,9,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_WORLD,ierr)
  if (myrank.eq.master) diag_sigma=work
  call MPI_REDUCE(diag_sigma_match,work,9,MPI_DOUBLE_PRECISION, &
       MPI_SUM,master,MPI_COMM_WORLD,ierr)
  if (myrank.eq.master) diag_sigma_match=work
  call MPI_REDUCE(diag_sigma_radial,work,9,MPI_DOUBLE_PRECISION, &
       MPI_SUM,master,MPI_COMM_WORLD,ierr)
  if (myrank.eq.master) diag_sigma_radial=work
  call MPI_REDUCE(diag_sigma_angular,work,9,MPI_DOUBLE_PRECISION, &
       MPI_SUM,master,MPI_COMM_WORLD,ierr)
  if (myrank.eq.master) diag_sigma_angular=work
  call MPI_REDUCE(diag_sigma_interstitial,work,9,MPI_DOUBLE_PRECISION, &
       MPI_SUM,master,MPI_COMM_WORLD,ierr)
  if (myrank.eq.master) diag_sigma_interstitial=work
  call MPI_REDUCE(diag_sigma_coeff_fd_mt,work,9,MPI_DOUBLE_PRECISION, &
       MPI_SUM,master,MPI_COMM_WORLD,ierr)
  if (myrank.eq.master) diag_sigma_coeff_fd_mt=work
  call MPI_REDUCE(diag_sigma_coeff_analytic_mt,work,9, &
       MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_WORLD,ierr)
  if (myrank.eq.master) diag_sigma_coeff_analytic_mt=work
  call MPI_REDUCE(diag_paired_window_matching,work,9, &
       MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_WORLD,ierr)
  if (myrank.eq.master) diag_paired_window_matching=work
  call MPI_REDUCE(diag_paired_norm_matching,work,9, &
       MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_WORLD,ierr)
  if (myrank.eq.master) diag_paired_norm_matching=work
  call MPI_REDUCE(diag_paired_main_correction_matching,work,9, &
       MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_WORLD,ierr)
  if (myrank.eq.master) diag_paired_main_correction_matching=work
  call MPI_REDUCE(diag_paired_window_complete,work,9, &
       MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_WORLD,ierr)
  if (myrank.eq.master) diag_paired_window_complete=work
  call MPI_REDUCE(diag_paired_norm_complete,work,9, &
       MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_WORLD,ierr)
  if (myrank.eq.master) diag_paired_norm_complete=work
  call MPI_REDUCE(diag_paired_main_correction_complete,work,9, &
       MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_WORLD,ierr)
  if (myrank.eq.master) diag_paired_main_correction_complete=work
  call MPI_REDUCE(diag_paired_window_kinetic_matching,work,9, &
       MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_WORLD,ierr)
  if (myrank.eq.master) diag_paired_window_kinetic_matching=work
  call MPI_REDUCE(diag_paired_norm_kinetic_matching,work,9, &
       MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_WORLD,ierr)
  if (myrank.eq.master) diag_paired_norm_kinetic_matching=work
  call MPI_REDUCE(diag_paired_main_correction_kinetic_matching,work,9, &
       MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_WORLD,ierr)
  if (myrank.eq.master) diag_paired_main_correction_kinetic_matching=work
  call MPI_REDUCE(diag_paired_window_kinetic_complete,work,9, &
       MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_WORLD,ierr)
  if (myrank.eq.master) diag_paired_window_kinetic_complete=work
  call MPI_REDUCE(diag_paired_norm_kinetic_complete,work,9, &
       MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_WORLD,ierr)
  if (myrank.eq.master) diag_paired_norm_kinetic_complete=work
  call MPI_REDUCE(diag_paired_main_correction_kinetic_complete,work,9, &
       MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_WORLD,ierr)
  if (myrank.eq.master) diag_paired_main_correction_kinetic_complete=work
  call MPI_REDUCE(diag_paired_sigma_matching,work,9, &
       MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_WORLD,ierr)
  if (myrank.eq.master) diag_paired_sigma_matching=work
  call MPI_REDUCE(diag_paired_sigma_complete,work,9, &
       MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_WORLD,ierr)
  if (myrank.eq.master) diag_paired_sigma_complete=work
  call MPI_REDUCE(diag_paired_diagonal_reference,work4,36, &
       MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_WORLD,ierr)
  if (myrank.eq.master) diag_paired_diagonal_reference=work4
  call MPI_REDUCE(diag_density_diagonal,work4,36,MPI_DOUBLE_PRECISION, &
       MPI_SUM,master,MPI_COMM_WORLD,ierr)
  if (myrank.eq.master) diag_density_diagonal=work4
  call MPI_REDUCE(diag_density_offdiagonal,work4,36,MPI_DOUBLE_PRECISION, &
       MPI_SUM,master,MPI_COMM_WORLD,ierr)
  if (myrank.eq.master) diag_density_offdiagonal=work4
  send_scalar(1)=trace_rhoE
  call MPI_REDUCE(send_scalar,work_scalar,1,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_WORLD,ierr)
  if (myrank.eq.master) trace_rhoE=work_scalar(1)
  send_scalar(1)=trace_erho
  call MPI_REDUCE(send_scalar,work_scalar,1,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_WORLD,ierr)
  if (myrank.eq.master) trace_erho=work_scalar(1)
  send_scalar(1)=trace_sigma
  call MPI_REDUCE(send_scalar,work_scalar,1,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_WORLD,ierr)
  if (myrank.eq.master) trace_sigma=work_scalar(1)
  send_scalar(1)=radial_hsph
  call MPI_REDUCE(send_scalar,work_scalar,1,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_WORLD,ierr)
  if (myrank.eq.master) radial_hsph=work_scalar(1)
  send_scalar(1)=radial_hnsp
  call MPI_REDUCE(send_scalar,work_scalar,1,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_WORLD,ierr)
  if (myrank.eq.master) radial_hnsp=work_scalar(1)
  send_scalar(1)=radial_overlap
  call MPI_REDUCE(send_scalar,work_scalar,1,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_WORLD,ierr)
  if (myrank.eq.master) radial_overlap=work_scalar(1)
  call MPI_REDUCE(radial_by_pass,work,9,MPI_DOUBLE_PRECISION,MPI_SUM, &
       master,MPI_COMM_WORLD,ierr)
  if (myrank.eq.master) radial_by_pass=work
  ! Wall-clock phase timings are reduced with MAX because MPI completion time
  ! is controlled by the slowest rank, not by the sum over ranks.
  call MPI_REDUCE(profile_times,work_profile,8,MPI_DOUBLE_PRECISION, &
       MPI_MAX,master,MPI_COMM_WORLD,ierr)
  if (myrank.eq.master) profile_times=work_profile
#endif
end subroutine Reduce_stress
