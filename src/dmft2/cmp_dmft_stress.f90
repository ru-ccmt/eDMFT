! Contract Blaha's LAPW strain kernels with the two stationary one-particle
! objects of the DFT+eDMFT Klein functional.
!
! For an energy-dependent LAPW basis the stationary derivative contains
!
!   -2 Tr[rho H_{,epsilon}] + 2 Tr[rhoE O_{,epsilon}],
!
! not a derivative of the non-stationary total-energy expression.  rho and
! rhoE are Hermitian matrices in the Kohn-Sham band window.  l2main
! diagonalizes each matrix and passes its eigenvectors as Aweight/AEweight;
! the eigenvalues are weight/eweight.  Consequently the existing diagonal
! WIEN2k kernels can be used exactly, once for each pseudo-band set.
!
! VAL_CORR_STRESS is linear in the band energy e:
!
!       K(e) = H_epsilon - e O_epsilon.
!
! Evaluating the unchanged, well-tested kernel at e=0 gives H_epsilon, and
! K(0)-K(1) gives O_epsilon.  This algebraic split avoids manually editing
! hundreds of angular-momentum terms and is independently testable in the
! DFT limit, where rhoE=e*rho must recover the original WIEN2k result.

subroutine cmp_dmft_stress(Aweight,AEweight,weight,zwe,nemin,nemax, &
     DM_nemin,DM_nemax,DM_nemaxx,nbands,nbandsx,n0,nnlo, &
     time_bl,time_bl_w,time_lo,time_lo_w)
  use param, only: lmax2,lomax,nume,nloat,nmat
  use com, only: nat
  use dmfts, only: iso
  use structure, only: mult,rmt
  use xa, only: E,fj,dfj
  use xa3, only: bk3
  use w_atpar, only: retrieve_w_atpar
  use alm_blm_clm, only: cmp_alm_blm_clm
  use stress_kb, only: load_stress_integrals
  use stress_pw, only: prepare_stress_pw
  use dmft_stress, only: stress_val_corr,stress_kin_mt,stress_apw, &
       stress_kin_int,apw_surface_coeff,stress_nsp_match, &
       stress_nsp_radial,stress_nsp_angular,stress_g2_mode, &
       stress_rela_corr,stress_rela_surf
  implicit none
  integer,intent(in) :: nemin,nemax,DM_nemin,DM_nemax,DM_nemaxx
  integer,intent(in) :: nbands,nbandsx,n0,nnlo
  complex*16,intent(in) :: Aweight(nbands,nbands),AEweight(nbands,nbands)
  real*8,intent(in) :: weight(nume),zwe(nbands)
  real*8,intent(inout) :: time_bl,time_bl_w,time_lo,time_lo_w

  integer :: lmfull,lmlo,jatom,lfirst,lmmax,latom,mu,is,num
  integer :: alpha,beta,jlo,isize
  real*8 :: esave(nume),eweight(nume)
  real*8 :: buf_kin(3,3),buf_corr(3,3),buf_apw(3,3)
  real*8 :: buf_match(3,3),buf_radial(3,3),buf_angular(3,3)
  real*8 :: buf_rela_corr(3,3),buf_rela_surf(3,3)
  real*8 :: g2_trace
  complex*16,allocatable :: alm(:,:),blm(:,:)
  complex*16,allocatable :: clm_compact(:,:,:),clm_stress(:,:,:)
  complex*16,allocatable :: aalm(:,:,:),bblm(:,:,:),cclm(:,:,:,:)
  complex*16,allocatable :: kh(:,:,:),ke0(:,:,:),ke1(:,:,:)
  complex*16,allocatable :: knsp_match(:,:,:),knsp_radial(:,:,:)
  complex*16,allocatable :: knsp_angular(:,:,:)
  complex*16,allocatable :: kkin(:,:,:),kapw(:,:,:)
  complex*16,allocatable :: krela_corr(:,:,:),krela_sr(:,:,:)
  complex*16,allocatable :: krela_old(:,:,:),krela_surf(:,:,:)

  lmfull=(lmax2+1)*(lmax2+1)
  lmlo=(lomax+1)*(lomax+1)
  isize=n0-nnlo
  allocate(alm(lmfull,nume),blm(lmfull,nume))
  ! dmft2 stores local-orbital coefficients only through LOMAX, whereas the
  ! historical WIEN2k stress kernels declare CLM with the full LMAX2 angular
  ! leading dimension.  Fortran explicit-shape arrays use that leading
  ! dimension to compute the band stride.  Passing the compact dmft2 array
  ! directly would therefore make every band after the first address the
  ! wrong memory (and, in practice, produced NaNs).  Keep the native compact
  ! array for CMP_ALM_BLM_CLM and explicitly pad it before every stress call.
  allocate(clm_compact(lmlo,nume,nloat),clm_stress(lmfull,nume,nloat))
  ! Force-response outputs are not referenced when Qforce_j=.false., but the
  ! allocatable dummy arguments must still be allocated for a portable call.
  allocate(aalm(1,1,1),bblm(1,1,1),cclm(1,1,1,1))
  allocate(kh(3,3,nume),ke0(3,3,nume),ke1(3,3,nume))
  allocate(knsp_match(3,3,nume),knsp_radial(3,3,nume), &
       knsp_angular(3,3,nume))
  allocate(kkin(3,3,nume),kapw(3,3,nume))
  allocate(krela_corr(3,3,nume),krela_sr(3,3,nume), &
       krela_old(3,3,nume),krela_surf(3,3,nume))

  eweight=0.d0
  do num=nemin,DM_nemin-1
     eweight(num)=E(num)*weight(num)
  enddo
  do num=DM_nemin,DM_nemaxx
     eweight(num)=zwe(num-DM_nemin+1)
  enddo

  do jatom=1,nat
     call retrieve_w_atpar(jatom,lfirst,lmmax)
     call load_stress_integrals(jatom)
     ! DFJ_ST/DDFJ_ST depend on RMT and on this atom's Bessel functions;
     ! F1P and KplusG are refreshed in the same call for clarity.
     call prepare_stress_pw(bk3,fj(:,:,jatom),dfj(:,:,jatom),rmt(jatom),isize)

     latom=lfirst-1
     do mu=1,mult(jatom)
        latom=latom+1
        do is=1,iso
           ! -------- rho pseudo-bands: kinetic, APW surface, and H kernel.
           call cmp_alm_blm_clm(alm,blm,clm_compact,aalm,bblm,cclm,Aweight, &
                latom,lfirst,jatom,nemin,nemax,n0,nnlo,is,DM_nemin, &
                DM_nemax,DM_nemaxx,nbands,nbandsx,.false., &
                time_bl,time_bl_w,time_lo,time_lo_w)

           clm_stress=0.d0
           clm_stress(1:lmlo,:,:)=clm_compact
           call val_kinstrmt(nemin,DM_nemaxx,alm,blm,clm_stress,kkin,rmt(jatom))
           ! VAL_KINSTRMT is the historical nonrelativistic p_alpha p_beta
           ! operator.  In CALC=RELA the radial functions instead solve the
           ! Koelling-Harmon equation.  Replace the old muffin-tin metric by
           ! Q_SR and repair the APW flux at R_MT.  The routine returns
           ! Q_SR-Q_old, so its occupation prefactor has the opposite sign
           ! to BUF_KIN, exactly as in patched WIEN2k LAPW2.
           call val_rela_kinstrmt(nemin,DM_nemaxx,jatom,alm,blm, &
                clm_stress,krela_corr,krela_sr,krela_old,krela_surf)
           call stress_APW_surface(nemin,DM_nemaxx,jatom,rmt(jatom), &
                alm,blm,clm_stress,apw_surface_coeff,kapw)

           esave=E
           E(nemin:DM_nemaxx)=0.d0
           call val_corr_stress(isize,isize,isize,nemin,DM_nemaxx, &
                alm,blm,clm_stress,kh,knsp_match,knsp_radial, &
                knsp_angular,.true.,jatom,rmt(jatom),latom,lfirst)
           E=esave

           buf_kin=0.d0; buf_corr=0.d0; buf_apw=0.d0
           buf_match=0.d0; buf_radial=0.d0; buf_angular=0.d0
           buf_rela_corr=0.d0; buf_rela_surf=0.d0
           do num=nemin,DM_nemaxx
              do beta=1,3
                 do alpha=1,3
                    buf_kin(alpha,beta)=buf_kin(alpha,beta)-2.d0* &
                         weight(num)*dble(kkin(alpha,beta,num))
                    buf_corr(alpha,beta)=buf_corr(alpha,beta)-2.d0* &
                         weight(num)*dble(kh(alpha,beta,num))
                    buf_apw(alpha,beta)=buf_apw(alpha,beta)+2.d0* &
                         weight(num)*dble(kapw(alpha,beta,num))
                    buf_rela_corr(alpha,beta)=buf_rela_corr(alpha,beta)+ &
                         2.d0*weight(num)*dble(krela_corr(alpha,beta,num))
                    buf_rela_surf(alpha,beta)=buf_rela_surf(alpha,beta)+ &
                         2.d0*weight(num)*dble(krela_surf(alpha,beta,num))
                    ! These are nonspherical pieces of H_,epsilon.  They
                    ! carry the same -2 Tr[rho H_,epsilon] prefactor as the
                    ! spherical valence-correction kernel, and no rhoE
                    ! counterpart because the overlap contains no V_ns.
                    buf_match(alpha,beta)=buf_match(alpha,beta)-2.d0* &
                         weight(num)*dble(knsp_match(alpha,beta,num))
                    buf_radial(alpha,beta)=buf_radial(alpha,beta)-2.d0* &
                         weight(num)*dble(knsp_radial(alpha,beta,num))
                    buf_angular(alpha,beta)=buf_angular(alpha,beta)-2.d0* &
                         weight(num)*dble(knsp_angular(alpha,beta,num))
                 enddo
              enddo
           enddo

           ! The interstitial term has no atom label.  At this point A still
           ! contains the rho pseudo-bands, so evaluate it once per spin and
           ! k point before the rhoE transformation overwrites A.
           if (jatom.eq.1 .and. mu.eq.1) then
              ! WEIGHT is the occupation/eigenvalue array for this particular
              ! contraction.  It is the physical occupation in the main pass
              ! and a signed response eigenvalue in either response pass.
              call val_kin_interstitial(isize,stress_kin_int,nemin, &
                   DM_nemaxx,weight)
           endif

           ! -------- rhoE pseudo-bands: only the overlap kernel is needed.
           call cmp_alm_blm_clm(alm,blm,clm_compact,aalm,bblm,cclm,AEweight, &
                latom,lfirst,jatom,nemin,nemax,n0,nnlo,is,DM_nemin, &
                DM_nemax,DM_nemaxx,nbands,nbandsx,.false., &
                time_bl,time_bl_w,time_lo,time_lo_w)
           clm_stress=0.d0
           clm_stress(1:lmlo,:,:)=clm_compact
           esave=E
           E(nemin:DM_nemaxx)=0.d0
           call val_corr_stress(isize,isize,isize,nemin,DM_nemaxx, &
                alm,blm,clm_stress,ke0,knsp_match,knsp_radial, &
                knsp_angular,.false.,jatom,rmt(jatom),latom,lfirst)
           E(nemin:DM_nemaxx)=1.d0
           call val_corr_stress(isize,isize,isize,nemin,DM_nemaxx, &
                alm,blm,clm_stress,ke1,knsp_match,knsp_radial, &
                knsp_angular,.false.,jatom,rmt(jatom),latom,lfirst)
           E=esave
           do num=nemin,DM_nemaxx
              do beta=1,3
                 do alpha=1,3
                    buf_corr(alpha,beta)=buf_corr(alpha,beta)+2.d0* &
                         eweight(num)*dble(ke0(alpha,beta,num)-ke1(alpha,beta,num))
                 enddo
              enddo
           enddo

           ! Keep every G2 component available for diagnostics, but include
           ! exactly one selected combination in the physical stress.  TRACE
           ! is invariant under the equivalent-atom rotation and is the
           ! current conservative default in patched WIEN2k-24.  FULL is the
           ! literal affine-coordinate formula (match+radial+angular).
           select case(trim(stress_g2_mode))
           case('NONE')
              continue
           case('MATCH')
              buf_corr=buf_corr+buf_match
           case('RADIAL')
              buf_corr=buf_corr+buf_radial
           case('ANGULAR')
              buf_corr=buf_corr+buf_angular
           case('FULL')
              buf_corr=buf_corr+buf_match+buf_radial+buf_angular
           case default                 ! TRACE
              g2_trace=(buf_match(1,1)+buf_match(2,2)+buf_match(3,3)+ &
                   buf_radial(1,1)+buf_radial(2,2)+buf_radial(3,3)+ &
                   buf_angular(1,1)+buf_angular(2,2)+ &
                   buf_angular(3,3))/3.d0
              do alpha=1,3
                 buf_corr(alpha,alpha)=buf_corr(alpha,alpha)+g2_trace
              enddo
           end select

           ! Rotate the equivalent atom's local tensor into the global frame
           ! exactly as lapw2 does before summing over all atoms and k points.
           call symmetry_stress_rotij(latom,jatom,buf_kin)
           call symmetry_stress_rotij(latom,jatom,buf_corr)
           call symmetry_stress_rotij(latom,jatom,buf_apw)
           call symmetry_stress_rotij(latom,jatom,buf_match)
           call symmetry_stress_rotij(latom,jatom,buf_radial)
           call symmetry_stress_rotij(latom,jatom,buf_angular)
           call symmetry_stress_rotij(latom,jatom,buf_rela_corr)
           call symmetry_stress_rotij(latom,jatom,buf_rela_surf)
           ! Include the relativistic replacement exactly once in the
           ! physical muffin-tin kinetic tensor.  The two named accumulators
           ! below are diagnostic copies only.
           buf_kin=buf_kin+buf_rela_corr+buf_rela_surf
           stress_kin_mt=stress_kin_mt+buf_kin
           stress_val_corr=stress_val_corr+buf_corr
           stress_apw=stress_apw+buf_apw
           stress_nsp_match=stress_nsp_match+buf_match
           stress_nsp_radial=stress_nsp_radial+buf_radial
           stress_nsp_angular=stress_nsp_angular+buf_angular
           stress_rela_corr=stress_rela_corr+buf_rela_corr
           stress_rela_surf=stress_rela_surf+buf_rela_surf

        enddo
     enddo
  enddo

  deallocate(alm,blm,clm_compact,clm_stress,aalm,bblm,cclm, &
       kh,ke0,ke1,knsp_match,knsp_radial,knsp_angular,kkin,kapw, &
       krela_corr,krela_sr,krela_old,krela_surf)
end subroutine cmp_dmft_stress


! Diagnostic-only contraction of an arbitrary Hermitian energy-moment matrix
! with the LAPW overlap derivative.  The physical stationary stress uses
!
!   rhoE = 1/2 {epsilon,rho} + M_Sigma,
!
! where M_Sigma is the frequency-summed embedded-self-energy moment.  This
! routine evaluates only 2 Tr[M_Sigma O_,epsilon].  It deliberately has no
! access to the global stress accumulators, so printing this decomposition
! cannot change the stress reported in case.scf2.
subroutine cmp_dmft_overlap_stress_diag(moment_vectors,moment_values, &
     sigma_stress,nemin,nemax,DM_nemin,DM_nemax,DM_nemaxx,nbands, &
     nbandsx,n0,nnlo,time_bl,time_bl_w,time_lo,time_lo_w)
  use param, only: lmax2,lomax,nume,nloat
  use com, only: nat
  use dmfts, only: iso
  use structure, only: mult,rmt
  use xa, only: E,fj,dfj
  use xa3, only: bk3
  use w_atpar, only: retrieve_w_atpar
  use alm_blm_clm, only: cmp_alm_blm_clm
  use stress_kb, only: load_stress_integrals
  use stress_pw, only: prepare_stress_pw
  implicit none
  integer,intent(in) :: nemin,nemax,DM_nemin,DM_nemax,DM_nemaxx
  integer,intent(in) :: nbands,nbandsx,n0,nnlo
  complex*16,intent(in) :: moment_vectors(nbands,nbands)
  real*8,intent(in) :: moment_values(nbands)
  real*8,intent(inout) :: sigma_stress(3,3)
  real*8,intent(inout) :: time_bl,time_bl_w,time_lo,time_lo_w

  integer :: lmfull,lmlo,jatom,lfirst,lmmax,latom,mu,is,num
  integer :: alpha,beta,isize
  real*8 :: esave(nume),eweight(nume),buf(3,3)
  complex*16,allocatable :: alm(:,:),blm(:,:)
  complex*16,allocatable :: clm_compact(:,:,:),clm_stress(:,:,:)
  complex*16,allocatable :: aalm(:,:,:),bblm(:,:,:),cclm(:,:,:,:)
  complex*16,allocatable :: ke0(:,:,:),ke1(:,:,:)
  complex*16,allocatable :: nsp0(:,:,:),nsp1(:,:,:),nsp2(:,:,:)

  lmfull=(lmax2+1)*(lmax2+1)
  lmlo=(lomax+1)*(lomax+1)
  isize=n0-nnlo
  allocate(alm(lmfull,nume),blm(lmfull,nume))
  allocate(clm_compact(lmlo,nume,nloat),clm_stress(lmfull,nume,nloat))
  allocate(aalm(1,1,1),bblm(1,1,1),cclm(1,1,1,1))
  allocate(ke0(3,3,nume),ke1(3,3,nume))
  allocate(nsp0(3,3,nume),nsp1(3,3,nume),nsp2(3,3,nume))

  eweight=0.d0
  do num=DM_nemin,DM_nemaxx
     eweight(num)=moment_values(num-DM_nemin+1)
  enddo

  do jatom=1,nat
     call retrieve_w_atpar(jatom,lfirst,lmmax)
     call load_stress_integrals(jatom)
     call prepare_stress_pw(bk3,fj(:,:,jatom),dfj(:,:,jatom), &
          rmt(jatom),isize)

     latom=lfirst-1
     do mu=1,mult(jatom)
        latom=latom+1
        do is=1,iso
           call cmp_alm_blm_clm(alm,blm,clm_compact,aalm,bblm,cclm, &
                moment_vectors,latom,lfirst,jatom,nemin,nemax,n0,nnlo, &
                is,DM_nemin,DM_nemax,DM_nemaxx,nbands,nbandsx,.false., &
                time_bl,time_bl_w,time_lo,time_lo_w)
           clm_stress=0.d0
           clm_stress(1:lmlo,:,:)=clm_compact

           esave=E
           E(nemin:DM_nemaxx)=0.d0
           call val_corr_stress(isize,isize,isize,nemin,DM_nemaxx, &
                alm,blm,clm_stress,ke0,nsp0,nsp1,nsp2,.false., &
                jatom,rmt(jatom),latom,lfirst)
           E(nemin:DM_nemaxx)=1.d0
           call val_corr_stress(isize,isize,isize,nemin,DM_nemaxx, &
                alm,blm,clm_stress,ke1,nsp0,nsp1,nsp2,.false., &
                jatom,rmt(jatom),latom,lfirst)
           E=esave

           buf=0.d0
           do num=DM_nemin,DM_nemaxx
              do beta=1,3
                 do alpha=1,3
                    buf(alpha,beta)=buf(alpha,beta)+2.d0*eweight(num)* &
                         dble(ke0(alpha,beta,num)-ke1(alpha,beta,num))
                 enddo
              enddo
           enddo
           call symmetry_stress_rotij(latom,jatom,buf)
           sigma_stress=sigma_stress+buf
        enddo
     enddo
  enddo

  deallocate(alm,blm,clm_compact,clm_stress,aalm,bblm,cclm, &
       ke0,ke1,nsp0,nsp1,nsp2)
end subroutine cmp_dmft_overlap_stress_diag
