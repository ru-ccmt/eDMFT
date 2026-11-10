! Contract Blaha's LAPW strain kernels with the two stationary one-particle
! objects of the DFT+eDMFT Klein functional.
!
! For an energy-dependent LAPW basis the stationary derivative contains
!
!   -2 Tr[rho H_{,epsilon}] + 2 Tr[rhoE O_{,epsilon}],
!
! in the stored WIEN convention for these one-sided local kernels; the
! physical whole-matrix derivative is Re Tr(rho H'-rhoE O'). This is
! not a derivative of the non-stationary total-energy expression.  rho and
! rhoE are Hermitian matrices in the Kohn-Sham band window.  l2main
! diagonalizes each matrix and passes its eigenvectors as Aweight/AEweight;
! the eigenvalues are weight/eweight. Each scalar kernel must therefore be
! the complete quadratic form on arbitrary complex vectors, not a formula
! simplified with the on-shell KS equation. The local G2 matrices below
! meet that requirement; the explicit metric/surface kernels are separate.
!
! The complete local G2 matrices are cached by stress_g2_matrices.
! Each one-sided kernel is linear in the test band energy E, so K(0)-K(1)
! extracts its overlap part exactly for arbitrary complex pseudo-bands.
! Unlike the legacy diagonal-only expansion, this includes all LO primitives,
! the shared LAPW/APW matching derivative, and strain-allowed V_ns channels.

subroutine cmp_dmft_stress(Aweight,AEweight,weight,zwe,nemin,nemax, &
     DM_nemin,DM_nemax,DM_nemaxx,nbands,nbandsx,n0,nnlo, &
     time_bl,time_bl_w,time_lo,time_lo_w,radial_pass)
  use param, only: lmax2,lomax,nume,nloat,nmat
  use com, only: nat
  use dmfts, only: iso
  use structure, only: mult,rmt
  use xa, only: E,fj,dfj
  use xa3, only: bk3
  use w_atpar, only: retrieve_w_atpar
  use alm_blm_clm, only: cmp_alm_blm_clm,cmp_alm_blm_clm_strain
  use stress_g2_matrices, only: evaluate_g2
  use stress_kb, only: load_stress_integrals
  use stress_pw, only: prepare_stress_pw
  use dmft_stress, only: stress_val_corr,stress_kin_mt,stress_apw, &
       stress_kin_int,apw_surface_coeff,stress_nsp_match, &
       stress_nsp_radial,stress_nsp_angular,stress_g2_mode, &
       stress_rela_corr,stress_rela_surf
  use dmft_stress, only: stress_diag_response_overlap
  use stress_radial_fd, only: radial_fd_enabled
  implicit none
  integer,intent(in) :: nemin,nemax,DM_nemin,DM_nemax,DM_nemaxx
  integer,intent(in) :: nbands,nbandsx,n0,nnlo
  complex*16,intent(in) :: Aweight(nbands,nbands),AEweight(nbands,nbands)
  real*8,intent(in) :: weight(nume),zwe(nbands)
  real*8,intent(inout) :: time_bl,time_bl_w,time_lo,time_lo_w
  integer,intent(in) :: radial_pass

  integer :: lmfull,lmlo,jatom,lfirst,lmmax,latom,mu,is,num
  integer :: alpha,beta,jlo,isize
  real*8 :: esave(nume),eweight(nume)
  real*8 :: buf_kin(3,3),buf_corr(3,3),buf_apw(3,3)
  real*8 :: buf_match(3,3),buf_radial(3,3),buf_angular(3,3)
  real*8 :: buf_rela_corr(3,3),buf_rela_surf(3,3)
  real*8 :: buf_overlap(3,3)
  real*8 :: g2_trace
  complex*16,allocatable :: salm(:,:,:),sblm(:,:,:),sclm(:,:,:,:)
  complex*16,allocatable :: alm(:,:),blm(:,:)
  complex*16,allocatable :: clm_compact(:,:,:),clm_stress(:,:,:)
  complex*16,allocatable :: aalm(:,:,:),bblm(:,:,:),cclm(:,:,:,:)
  complex*16,allocatable :: kh(:,:,:),ke0(:,:,:),ke1(:,:,:)
  complex*16,allocatable :: ksph0(:,:,:),ksph1(:,:,:),ksph2(:,:,:)
  complex*16,allocatable :: knsp_match(:,:,:),knsp_radial(:,:,:)
  complex*16,allocatable :: knsp_angular(:,:,:)
  complex*16,allocatable :: kkin(:,:,:),kapw(:,:,:)
  complex*16,allocatable :: krela_corr(:,:,:),krela_sr(:,:,:)
  complex*16,allocatable :: krela_old(:,:,:),krela_surf(:,:,:)

  lmfull=(lmax2+1)*(lmax2+1)
  lmlo=(lomax+1)*(lomax+1)
  isize=n0-nnlo
  allocate(alm(lmfull,nume),blm(lmfull,nume))
  allocate(salm(lmfull,nume,6),sblm(lmfull,nume,6),sclm(lmlo,nume,nloat,6))
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
  allocate(ksph0(3,3,nume),ksph1(3,3,nume),ksph2(3,3,nume))
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

           call cmp_alm_blm_clm_strain(salm,sblm,sclm,alm,blm, &
                clm_compact,Aweight,latom,lfirst,jatom,n0,nnlo,is, &
                DM_nemin,DM_nemax,DM_nemaxx,nbands,nbandsx,first_band=nemin)
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
           if (radial_fd_enabled) call radial_fd_band_contraction( &
                alm,blm,clm_stress,weight,nemin,DM_nemaxx,jatom,1, &
                radial_pass)

           esave=E
           E(nemin:DM_nemaxx)=0.d0
           call evaluate_g2(jatom,nemin,DM_nemaxx,alm,blm,clm_compact, &
                salm,sblm,sclm,E,kh,knsp_match,knsp_radial,knsp_angular, &
                ksph0,ksph1,ksph2,.true.)
           E=esave

           buf_kin=0.d0; buf_corr=0.d0; buf_apw=0.d0
           buf_match=0.d0; buf_radial=0.d0; buf_angular=0.d0
           buf_rela_corr=0.d0; buf_rela_surf=0.d0
           buf_overlap=0.d0
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
           call cmp_alm_blm_clm_strain(salm,sblm,sclm,alm,blm, &
                clm_compact,AEweight,latom,lfirst,jatom,n0,nnlo,is, &
                DM_nemin,DM_nemax,DM_nemaxx,nbands,nbandsx,first_band=nemin)
           clm_stress=0.d0
           clm_stress(1:lmlo,:,:)=clm_compact
           if (radial_fd_enabled) call radial_fd_band_contraction( &
                alm,blm,clm_stress,eweight,nemin,DM_nemaxx,jatom,2, &
                radial_pass)
           esave=E
           E(nemin:DM_nemaxx)=0.d0
           call evaluate_g2(jatom,nemin,DM_nemaxx,alm,blm,clm_compact, &
                salm,sblm,sclm,E,ke0,knsp_match,knsp_radial,knsp_angular, &
                ksph0,ksph1,ksph2,.false.)
           E(nemin:DM_nemaxx)=1.d0
           call evaluate_g2(jatom,nemin,DM_nemaxx,alm,blm,clm_compact, &
                salm,sblm,sclm,E,ke1,knsp_match,knsp_radial,knsp_angular, &
                ksph0,ksph1,ksph2,.false.)
           E=esave
           do num=nemin,DM_nemaxx
              do beta=1,3
                 do alpha=1,3
                    buf_corr(alpha,beta)=buf_corr(alpha,beta)+2.d0* &
                         eweight(num)*dble(ke0(alpha,beta,num)-ke1(alpha,beta,num))
                    buf_overlap(alpha,beta)=buf_overlap(alpha,beta)+2.d0* &
                         eweight(num)*dble(ke0(alpha,beta,num)-ke1(alpha,beta,num))
                 enddo
              enddo
           enddo

           ! Keep every G2 component available for diagnostics, but include
           ! exactly one selected combination in the physical stress. FULL
           ! uses the complete G2 matrices; TRACE is retained as an explicit
           ! historical approximation, never as an automatic shear filter.
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
           call symmetry_stress_rotij(latom,jatom,buf_overlap)
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
           if (radial_pass.ge.1 .and. radial_pass.le.3) &
                stress_diag_response_overlap(:,:,radial_pass)= &
                stress_diag_response_overlap(:,:,radial_pass)+buf_overlap

        enddo
     enddo
  enddo

  deallocate(salm,sblm,sclm)
  deallocate(alm,blm,clm_compact,clm_stress,aalm,bblm,cclm, &
       kh,ke0,ke1,ksph0,ksph1,ksph2,knsp_match,knsp_radial, &
       knsp_angular,kkin,kapw, &
       krela_corr,krela_sr,krela_old,krela_surf)
end subroutine cmp_dmft_stress


! Build matrix elements of the measured radial-basis response.  This routine
! is called once, before the k loop, only when
! EDMFT_STRESS_RADIAL_RESPONSE points to a finite-difference response file.
! It is intentionally a diagnostic reference: no matrix below is inserted in
! the production tensor.
subroutine prepare_stress_radial_fd
  use param, only: nrad,lmax2,lomax,nloat,ncom,ngau,filename_V_vns
  use com, only: nat,rel
  use structure, only: jri,dx,r0,mult
  use atspdt, only: e_store,el
  use lo, only: elo_store,elo,loor,rlo,lapw,nlo,nlov,nlon,ilo, &
       alo,blo,clo
  use p_project, only: P_rfk,drix_mat,n_al_ucase,l_al_ucase, &
       j_al_ucase
  use w_atpar, only: w_RRAD1,w_RRAD2,w_RADE1,w_RADE2, &
       w_RADP1,w_RADP2,w_RADPE1,w_RADPE2,w_a1lo,w_b1lo, &
       w_lm,w_lmmax,Vts,lv,lpv,llv,mv,mpv,mmv,Nhmx, &
       retrieve_w_atpar
  use readPotential, only: init_V_vns,read_V_vns,close_V_vns
  use stress_radial_fd
  implicit none
  external :: get_ilo
  integer :: jatom,l,jlo,ir,i,j,ih,np,nq,nr,ivlm,ilm,best_ivlm, &
       iind,j_apw,n_apw,lfirst,lmmax
  real*8 :: det,rhs1,rhs2,da,db,scale,integral0,residual, &
       max_residual,denom_int,denom_vts,candidate_residual, &
       best_candidate_residual
  ! Keep the atom index explicitly.  The nonspherical pass is separate from
  ! the spherical pass below; atom-less work arrays would otherwise contain
  ! only the last atom when VNS is assembled.
  real*8 :: central1(nrad,2+nloat,0:lmax2,nat), &
       central2(nrad,2+nloat,0:lmax2,nat)
  real*8 :: response1(nrad,2+nloat,0:lmax2,nat), &
       response2(nrad,2+nloat,0:lmax2,nat)
  real*8 :: hket1(nrad),hket2(nrad),tmp1(nrad),tmp2(nrad)
  real*8 :: Vlm(nrad,0:ncom-1)
  real*8 :: radial0(2+nloat,2+nloat)
  complex*16 :: angular_factor,candidate_factor,numerator
  real*8 :: projector_response(2+nloat),ju,judot,jlo_response, &
       q1_kernel
  real*8,external :: Rint13n
  interface
     function RINT13(REL,A,B,X,Y,NRAD_,DX_,JRI_,R0_) result(S)
       logical,intent(in) :: REL
       integer,intent(in) :: NRAD_,JRI_
       real*8,intent(in) :: A(NRAD_),B(NRAD_),X(NRAD_),Y(NRAD_),DX_,R0_
       real*8 :: S
     end function RINT13
  end interface

  if (.not.radial_fd_enabled) return
  call allocate_stress_radial_fd_matrices(nat)
  central1=0.d0; central2=0.d0
  response1=0.d0; response2=0.d0

  ! First form response functions whose interstitial matching data are held
  ! fixed.  LAPW has two boundary conditions; APW+lo has only the value
  ! condition.  Complete LO responses already contain the differentiated ABC
  ! coefficients and therefore need no further boundary repair.
  do jatom=1,nat
     call get_ilo(e_store,elo_store,jatom,mult,nat,lmax2,lomax,nloat, &
          el,elo,loor,rlo,lapw,nlov,nlo,nlon,ilo)
     nr=jri(jatom)
     do l=0,lmax2
        if (.not.all(radial_fd_seen(1:nr,l,jatom))) then
           write(6,*) 'Incomplete :DR records for atom/l ',jatom,l
           stop 'incomplete EDMFT radial response'
        endif
        central1(:,1,l,jatom)=w_RRAD1(:,l,jatom)
        central2(:,1,l,jatom)=w_RRAD2(:,l,jatom)
        central1(:,2,l,jatom)=w_RADE1(:,l,jatom)
        central2(:,2,l,jatom)=w_RADE2(:,l,jatom)

        if (lapw(l)) then
           det=w_RRAD1(nr,l,jatom)*w_RADPE1(nr,l,jatom)- &
                w_RADE1(nr,l,jatom)*w_RADP1(nr,l,jatom)
           if (abs(det).lt.1.d-14) stop 'singular LAPW radial boundary matrix'
           rhs1=-radial_fd_du(1,nr,l,jatom)
           rhs2=-radial_fd_dup(1,nr,l,jatom)
           da=(rhs1*w_RADPE1(nr,l,jatom)- &
                w_RADE1(nr,l,jatom)*rhs2)/det
           db=(w_RRAD1(nr,l,jatom)*rhs2-rhs1* &
                w_RADP1(nr,l,jatom))/det
           response1(:,1,l,jatom)=radial_fd_du(1,:,l,jatom)+ &
                da*w_RRAD1(:,l,jatom)+db*w_RADE1(:,l,jatom)
           response2(:,1,l,jatom)=radial_fd_du(2,:,l,jatom)+ &
                da*w_RRAD2(:,l,jatom)+db*w_RADE2(:,l,jatom)

           rhs1=-radial_fd_dudot(1,nr,l,jatom)
           rhs2=-radial_fd_dudotp(1,nr,l,jatom)
           da=(rhs1*w_RADPE1(nr,l,jatom)- &
                w_RADE1(nr,l,jatom)*rhs2)/det
           db=(w_RRAD1(nr,l,jatom)*rhs2-rhs1* &
                w_RADP1(nr,l,jatom))/det
           response1(:,2,l,jatom)=radial_fd_dudot(1,:,l,jatom)+ &
                da*w_RRAD1(:,l,jatom)+db*w_RADE1(:,l,jatom)
           response2(:,2,l,jatom)=radial_fd_dudot(2,:,l,jatom)+ &
                da*w_RRAD2(:,l,jatom)+db*w_RADE2(:,l,jatom)
        else
           if (abs(w_RRAD1(nr,l,jatom)).lt.1.d-14) &
                stop 'singular APW radial boundary value'
           scale=radial_fd_du(1,nr,l,jatom)/w_RRAD1(nr,l,jatom)
           response1(:,1,l,jatom)=radial_fd_du(1,:,l,jatom)- &
                scale*w_RRAD1(:,l,jatom)
           response2(:,1,l,jatom)=radial_fd_du(2,:,l,jatom)- &
                scale*w_RRAD2(:,l,jatom)
        endif

        if (l.le.lomax) then
           do jlo=1,ilo(l)
              if (.not.all(radial_fd_lo_seen(1:nr,jlo,l,jatom))) then
                 write(6,*) 'Incomplete :DL records for atom/l/jlo ', &
                      jatom,l,jlo
                 stop 'incomplete EDMFT local-orbital response'
              endif
              central1(:,2+jlo,l,jatom)=w_a1lo(:,jlo,l,jatom)
              central2(:,2+jlo,l,jatom)=w_b1lo(:,jlo,l,jatom)
              response1(:,2+jlo,l,jatom)=radial_fd_dlo(1,:,jlo,l,jatom)
              response2(:,2+jlo,l,jatom)=radial_fd_dlo(2,:,jlo,l,jatom)
           enddo
        endif

        do i=1,2+ilo(l)
           do j=1,2+ilo(l)
              radial_fd_overlap(i,j,l,jatom)=RINT13(rel, &
                   response1(:,i,l,jatom),response2(:,i,l,jatom), &
                   central1(:,j,l,jatom),central2(:,j,l,jatom), &
                   nrad,dx(jatom),nr,r0(jatom))
              select case(j)
              case(1)
                 hket1=el(l)*central1(:,1,l,jatom)
                 hket2=el(l)*central2(:,1,l,jatom)
              case(2)
                 ! ATPARN stores udot=d u/d E_Ry, hence
                 ! (H_Ry-E_l)udot=u with coefficient one.
                 hket1=el(l)*central1(:,2,l,jatom)+central1(:,1,l,jatom)
                 hket2=el(l)*central2(:,2,l,jatom)+central2(:,1,l,jatom)
              case default
                 ! The ket column is the primitive u_LO,j and its
                 ! coefficient is CLM=q_j*C_j.  The response-bra column is
                 ! the derivative of the complete local orbital L_j and its
                 ! coefficient is q_j.  Do not apply the Hamiltonian of the
                 ! complete L_j here: doing so would count A_j/B_j twice.
                 hket1=elo(l,j-2)*central1(:,j,l,jatom)
                 hket2=elo(l,j-2)*central2(:,j,l,jatom)
              end select
              radial_fd_hsph(i,j,l,jatom)=RINT13(rel, &
                   response1(:,i,l,jatom),response2(:,i,l,jatom), &
                   hket1,hket2, &
                   nrad,dx(jatom),nr,r0(jatom))
           enddo
        enddo
     enddo
  enddo

  ! Optional paired projector-5 response.  The radial H/O diagnostic above
  ! is not basis covariant by itself: the same radial change also changes the
  ! raw local projector P=<psi|phi> and the global normalization
  ! X=(sum_k w_k P_k^+P_k)^(-1/2).  Build_DMFT_Projector and cmp_overlap
  ! already propagate d(rix_mat) through both P_,q X and P X_,q.  What they
  ! need here is the derivative expressed in the *central primitive radial
  ! coefficients* (ALM,BLM,CLM), even though the finite-difference response
  ! was constructed in the boundary-homogenized basis
  ! (a,b,q_j)=(APW residuals, complete-LO amplitudes).
  !
  ! For an ordinary LAPW channel,
  !
  !   q_j=C_j/CLO_j,
  !   a=A-sum_j q_j ALO_j,  b=B-sum_j q_j BLO_j,
  !
  ! so the effective primitive LO derivative is
  !
  !   J'_j=(J_Lj-ALO_j J_u-BLO_j J_udot)/CLO_j.
  !
  ! APW+lo has CLO_1=0.  Its first local-orbital amplitude is carried by B,
  ! and the equivalent linear transformation is written explicitly below.
  ! The input file is a hydrostatic derivative
  ! d/deta=d/d eps_xx+d/d eps_yy+d/d eps_zz; distributing it equally over
  ! the three diagonal channels makes the trace recover that derivative and
  ! leaves every shear component zero.  This path is diagnostic-only and is
  ! activated separately by EDMFT_STRESS_RADIAL_PROJECTOR.
  if (radial_fd_projector_enabled) then
     drix_mat=0.d0
     do iind=1,n_al_ucase
        l=l_al_ucase(iind); jatom=j_al_ucase(iind); nr=jri(jatom)
        ! GET_ILO reconstructs the LO index flags, but ALO/BLO/CLO are
        ! atom-dependent cached ATPAR results.  Load the complete atom cache
        ! before transforming the boundary-homogenized response back to the
        ! primitive ALM/BLM/CLM coefficients.
        call retrieve_w_atpar(jatom,lfirst,lmmax)
        projector_response=0.d0
        do i=1,2+ilo(l)
           ! P_rfk is allocated to the projector file's actual radial
           ! length, not PARAM%NRAD.  Use the same finite-length contraction
           ! as p_cmp_rixmat; passing it to the full-size RINT13 interface
           ! can read beyond the allocated projector array.
           projector_response(i)=Rint13n(nr, &
                response1(:nr,i,l,jatom),response2(:nr,i,l,jatom), &
                P_rfk(:nr,1,iind),P_rfk(:nr,2,iind), &
                r0(jatom),dx(jatom),rel)
        enddo
        ju=projector_response(1); judot=projector_response(2)
        if (lapw(l)) then
           drix_mat(1,iind,1:3)=ju/3.d0
           drix_mat(2,iind,1:3)=judot/3.d0
           do jlo=1,ilo(l)
              if (.not.loor(jlo,l)) cycle
              if (abs(clo(l,jlo)).lt.1.d-12) &
                   stop 'zero CLO in LAPW radial-projector response'
              jlo_response=(projector_response(2+jlo)- &
                   alo(l,jlo)*ju-blo(l,jlo)*judot)/clo(l,jlo)
              drix_mat(2+jlo,iind,1:3)=jlo_response/3.d0
           enddo
        else
           j_apw=0; n_apw=0
           do jlo=1,ilo(l)
              if (abs(clo(l,jlo)).lt.1.d-12 .and. &
                  abs(blo(l,jlo)).gt.1.d-12) then
                 j_apw=jlo; n_apw=n_apw+1
              endif
           enddo
           if (n_apw.ne.1) then
              write(6,'(A,4I5)') ':RADIAL_PROJECTOR_APW_LO_FAILURE ', &
                   iind,jatom,l,ilo(l)
              do jlo=1,ilo(l)
                 write(6,'(A,I5,3ES22.12,2L3)') &
                      ':RADIAL_PROJECTOR_LO_COEFF ',jlo,alo(l,jlo), &
                      blo(l,jlo),clo(l,jlo),loor(jlo,l),lapw(l)
              enddo
              stop 'unsupported APW+lo radial-projector transformation'
           endif
           q1_kernel=(projector_response(2+j_apw)- &
                alo(l,j_apw)*ju)/blo(l,j_apw)
           drix_mat(1,iind,1:3)=ju/3.d0
           drix_mat(2,iind,1:3)=q1_kernel/3.d0
           do jlo=1,ilo(l)
              if (.not.loor(jlo,l)) cycle
              if (jlo.eq.j_apw) cycle
              if (abs(clo(l,jlo)).lt.1.d-12) &
                   stop 'zero higher CLO in APW+lo radial-projector response'
              jlo_response=(projector_response(2+jlo)- &
                   alo(l,jlo)*ju-blo(l,jlo)*q1_kernel)/clo(l,jlo)
              drix_mat(2+jlo,iind,1:3)=jlo_response/3.d0
           enddo
        endif
     enddo
     write(6,'(A,1x,ES24.15)') &
          ':DMFT_STRESS_RADIAL_PROJECTOR_DRIX_NORM',sqrt(sum(drix_mat**2))
     ! The response file stores the endpoint convention
     ! (u(+h)-u(-h))/(2h).  Keep the historical sign as the default, but
     ! permit a projector-only closure test of the convention used when the
     ! measured radial overlaps are inserted into Build_DMFT_Projector.  The
     ! switch is diagnostic-only because this entire branch requires
     ! EDMFT_STRESS_RADIAL_PROJECTOR explicitly.
     drix_mat=radial_fd_projector_sign*drix_mat
     write(6,'(A,1x,ES24.15)') &
          ':DMFT_STRESS_RADIAL_PROJECTOR_SIGNED_DRIX_NORM', &
          sqrt(sum(drix_mat**2))
  endif

  ! Nonspherical response.  Recover each record's angular factor from the
  ! already validated central Vts matrix, then replace its bra radial
  ! function by the measured response.  This avoids reimplementing WIEN's
  ! real/cubic-harmonic phase conventions.
  call init_V_vns(filename_V_vns,9918)
  max_residual=0.d0
  do jatom=1,nat
     call get_ilo(e_store,elo_store,jatom,mult,nat,lmax2,lomax,nloat, &
          el,elo,loor,rlo,lapw,nlov,nlo,nlon,ilo)
     nr=jri(jatom); Vlm=0.d0
     call read_V_vns(Vlm,jatom,nrad,nr,w_lmmax(jatom))
     do ih=1,Nhmx(jatom)
        np=2+ilo(lpv(ih,jatom)); nq=2+ilo(lv(ih,jatom))
        ! The (llv,mmv) labels in case.nsh are complex-harmonic coupling
        ! labels, whereas case.vns is stored in WIEN's real/cubic harmonic
        ! order.  Do not identify the two index sets directly.  Instead find
        ! the radial VLM component whose complete central radial matrix is
        ! proportional to the already trusted Vts record.  The proportionality
        ! constant is precisely its angular Gaunt/phase factor.
        best_ivlm=0; best_candidate_residual=huge(1.d0)
        angular_factor=(0.d0,0.d0); denom_vts=0.d0
        do i=1,np
           do j=1,nq
              denom_vts=denom_vts+abs(Vts(i,j,ih,jatom))**2
           enddo
        enddo
        if (denom_vts.lt.1.d-28) cycle
        do ilm=1,w_lmmax(jatom)-1
           denom_int=0.d0; numerator=(0.d0,0.d0)
           do i=1,np
              tmp1=central1(:,i,lpv(ih,jatom),jatom)*Vlm(:,ilm)
              tmp2=central2(:,i,lpv(ih,jatom),jatom)*Vlm(:,ilm)
              do j=1,nq
                 radial0(i,j)=RINT13(rel,tmp1,tmp2, &
                      central1(:,j,lv(ih,jatom),jatom), &
                      central2(:,j,lv(ih,jatom),jatom), &
                      nrad,dx(jatom),nr,r0(jatom))
                 denom_int=denom_int+radial0(i,j)**2
                 numerator=numerator+radial0(i,j)*Vts(i,j,ih,jatom)
              enddo
           enddo
           if (denom_int.lt.1.d-30) cycle
           candidate_factor=numerator/denom_int
           candidate_residual=0.d0
           do i=1,np
              do j=1,nq
                 candidate_residual=candidate_residual+ &
                      abs(Vts(i,j,ih,jatom)- &
                      candidate_factor*radial0(i,j))**2
              enddo
           enddo
           candidate_residual=sqrt(candidate_residual/denom_vts)
           if (candidate_residual.lt.best_candidate_residual) then
              best_candidate_residual=candidate_residual
              best_ivlm=ilm; angular_factor=candidate_factor
           endif
        enddo
        if (best_ivlm.eq.0) then
           write(6,*) 'Cannot reconstruct VNS record ',jatom,ih
           stop 'radial-response VNS reconstruction failed'
        endif
        ivlm=best_ivlm
        max_residual=max(max_residual,best_candidate_residual)
        do i=1,np
           tmp1=response1(:,i,lpv(ih,jatom),jatom)*Vlm(:,ivlm)
           tmp2=response2(:,i,lpv(ih,jatom),jatom)*Vlm(:,ivlm)
           do j=1,nq
              radial_fd_vnsp(i,j,ih,jatom)=angular_factor*RINT13(rel, &
                   tmp1,tmp2,central1(:,j,lv(ih,jatom),jatom), &
                   central2(:,j,lv(ih,jatom),jatom), &
                   nrad,dx(jatom),nr,r0(jatom))
              integral0=RINT13(rel, &
                   central1(:,i,lpv(ih,jatom),jatom)*Vlm(:,ivlm), &
                   central2(:,i,lpv(ih,jatom),jatom)*Vlm(:,ivlm), &
                   central1(:,j,lv(ih,jatom),jatom), &
                   central2(:,j,lv(ih,jatom),jatom), &
                   nrad,dx(jatom),nr,r0(jatom))
              ! Per-entry residuals can be ill-conditioned at symmetry
              ! zeros.  The normalized Frobenius residual accumulated above
              ! is the stable reconstruction check reported to the user.
           enddo
        enddo
     enddo
  enddo
  call close_V_vns()
  radial_fd_prepared=.true.
  write(6,'(A,1x,ES24.15)') &
       ':DMFT_STRESS_RADIAL_FD_VNS_RECON_MAX_REL',max_residual
end subroutine prepare_stress_radial_fd


! Contract one pseudo-band representation of rho or rhoE with the measured
! radial response.  WHICH=1 accumulates 2 Re <dpsi|H|psi>; WHICH=2 accumulates
! 2 Re <dpsi|psi>.  Their stationary Klein combination is H_rho-O_rhoE.
subroutine radial_fd_band_contraction(alm,blm,clm,weights,nemin,nemax, &
     jatom,which,radial_pass)
  use param, only: lmax2,lomax,nloat,nume,ngau
  use lo, only: ilo,lapw,alo,blo,clo
  use w_atpar, only: Nhmx,lv,lpv,mv,mpv
  use stress_radial_fd
  implicit none
  integer,intent(in) :: nemin,nemax,jatom,which,radial_pass
  complex*16,intent(in) :: alm((lmax2+1)*(lmax2+1),nume), &
       blm((lmax2+1)*(lmax2+1),nume), &
       clm((lmax2+1)*(lmax2+1),nume,nloat)
  real*8,intent(in) :: weights(nume)
  integer :: num,l,m,lm,ih,lmb,lmk,i,j,np,nq,jlo
  complex*16 :: central(2+nloat),response(2+nloat)
  complex*16 :: ket(2+nloat),bra_response(2+nloat)
  complex*16 :: zoverlap,zhsph,zhnsp

  if (.not.radial_fd_prepared) stop 'radial FD matrices were not prepared'
  do num=nemin,nemax
     if (weights(num).eq.0.d0) cycle
     zoverlap=(0.d0,0.d0); zhsph=(0.d0,0.d0); zhnsp=(0.d0,0.d0)
     do l=0,lmax2
        np=2+ilo(l)
        do m=-l,l
           lm=l*(l+1)+m+1
           call radial_fd_coefficients(l,lm,num,alm,blm,clm,central,response)
           do j=1,np
              do i=1,np
                 zoverlap=zoverlap+conjg(response(i))* &
                      radial_fd_overlap(i,j,l,jatom)*central(j)
                 zhsph=zhsph+conjg(response(i))* &
                      radial_fd_hsph(i,j,l,jatom)*central(j)
              enddo
           enddo
        enddo
     enddo
     do ih=1,Nhmx(jatom)
        lmb=lpv(ih,jatom)*(lpv(ih,jatom)+1)+mpv(ih,jatom)+1
        lmk=lv(ih,jatom)*(lv(ih,jatom)+1)+mv(ih,jatom)+1
        np=2+ilo(lpv(ih,jatom)); nq=2+ilo(lv(ih,jatom))
        call radial_fd_coefficients(lpv(ih,jatom),lmb,num,alm,blm,clm, &
             central,bra_response)
        call radial_fd_coefficients(lv(ih,jatom),lmk,num,alm,blm,clm, &
             ket,response)
        do j=1,nq
           do i=1,np
              zhnsp=zhnsp+conjg(bra_response(i))* &
                   radial_fd_vnsp(i,j,ih,jatom)*ket(j)
           enddo
        enddo
     enddo
     select case(which)
     case(1)
        radial_fd_rho_hsph=radial_fd_rho_hsph+ &
             2.d0*weights(num)*dble(zhsph)
        radial_fd_rho_hnsp=radial_fd_rho_hnsp+ &
             2.d0*weights(num)*dble(zhnsp)
        if (radial_pass.ge.1 .and. radial_pass.le.3) then
           radial_fd_by_pass(radial_pass,1)= &
                radial_fd_by_pass(radial_pass,1)+ &
                2.d0*weights(num)*dble(zhsph)
           radial_fd_by_pass(radial_pass,2)= &
                radial_fd_by_pass(radial_pass,2)+ &
                2.d0*weights(num)*dble(zhnsp)
        endif
     case(2)
        radial_fd_rhoe_overlap=radial_fd_rhoe_overlap+ &
             2.d0*weights(num)*dble(zoverlap)
        if (radial_pass.ge.1 .and. radial_pass.le.3) &
             radial_fd_by_pass(radial_pass,3)= &
             radial_fd_by_pass(radial_pass,3)+ &
             2.d0*weights(num)*dble(zoverlap)
     end select
  enddo
contains
  subroutine radial_fd_coefficients(l,lm,num,alm,blm,clm,central,response)
    integer,intent(in) :: l,lm,num
    complex*16,intent(in) :: alm((lmax2+1)*(lmax2+1),nume), &
         blm((lmax2+1)*(lmax2+1),nume), &
         clm((lmax2+1)*(lmax2+1),nume,nloat)
    complex*16,intent(out) :: central(2+nloat),response(2+nloat)
    complex*16 :: amplitude(nloat),apw_a,apw_b
    integer :: jlo,j_apw,n_apw
    central=(0.d0,0.d0); response=(0.d0,0.d0)
    amplitude=(0.d0,0.d0)
    central(1)=alm(lm,num); central(2)=blm(lm,num)
    do jlo=1,ilo(l)
       central(2+jlo)=clm(lm,num,jlo)
       if (abs(clo(l,jlo)).gt.1.d-12) &
            amplitude(jlo)=clm(lm,num,jlo)/clo(l,jlo)
    enddo
    if (.not.lapw(l)) then
       j_apw=0; n_apw=0
       do jlo=1,ilo(l)
          if (abs(clo(l,jlo)).le.1.d-12 .and. &
              abs(blo(l,jlo)).gt.1.d-12) then
             j_apw=jlo; n_apw=n_apw+1
          endif
       enddo
       if (n_apw.ne.1) stop 'cannot identify unique APW+lo amplitude'
       apw_b=blm(lm,num)
       do jlo=1,ilo(l)
          if (jlo.eq.j_apw) cycle
          apw_b=apw_b-amplitude(jlo)*blo(l,jlo)
       enddo
       amplitude(j_apw)=apw_b/blo(l,j_apw)
    endif
    apw_a=alm(lm,num); apw_b=blm(lm,num)
    do jlo=1,ilo(l)
       apw_a=apw_a-amplitude(jlo)*alo(l,jlo)
       apw_b=apw_b-amplitude(jlo)*blo(l,jlo)
       response(2+jlo)=amplitude(jlo)
    enddo
    response(1)=apw_a; response(2)=apw_b
  end subroutine radial_fd_coefficients
end subroutine radial_fd_band_contraction


! Default-off audit of the matrix generalization used by CMP_DMFT_STRESS.
! RHO and RHOE are arbitrary Hermitian matrices in the correlated KS-band
! window.  They are diagonalized into signed pseudo-bands and contracted by
! the unchanged analytic kernels.  Crucially, every module accumulator that
! CMP_DMFT_STRESS touches is restored before return.  The only observable
! result is the increment of DIAGNOSTIC, so this routine cannot perturb the
! production stress even when the audit is enabled.
subroutine cmp_dmft_stress_matrix_diagnostic(rho,rhoE,diagnostic,label, &
     ikp,nemin,nemax,DM_nemin,DM_nemax,nbands,n0,nnlo, &
     time_bl,time_bl_w,time_lo,time_lo_w)
  use param, only: nume
  use dmft_stress, only: stress_val_corr,stress_kin_mt,stress_apw, &
       stress_kin_int,stress_nsp_match,stress_nsp_radial, &
       stress_nsp_angular,stress_rela_corr,stress_rela_surf
  use stress_radial_fd, only: radial_fd_rho_hsph,radial_fd_rho_hnsp, &
       radial_fd_rhoe_overlap
  implicit none
  integer,intent(in) :: ikp,nemin,nemax,DM_nemin,DM_nemax,nbands,n0,nnlo
  complex*16,intent(in) :: rho(nbands,nbands),rhoE(nbands,nbands)
  real*8,intent(inout) :: diagnostic(3,3,4)
  character(len=*),intent(in) :: label
  real*8,intent(inout) :: time_bl,time_bl_w,time_lo,time_lo_w
  complex*16,allocatable :: rhoVectors(:,:),rhoEVectors(:,:)
  real*8,allocatable :: rhoValues(:),rhoEValues(:),weights(:)
  real*8 :: before(3,3,9),before_radial(3)
  integer :: audit_nemaxx,nbandsx

  allocate(rhoVectors(nbands,nbands),rhoEVectors(nbands,nbands), &
       rhoValues(nbands),rhoEValues(nbands),weights(nume))
  rhoVectors=rho
  rhoEVectors=rhoE
  audit_nemaxx=DM_nemax
  call Diagonalize_DMFT_WEIGHTS(rhoValues,rhoVectors,nbands,DM_nemin, &
       audit_nemaxx,.False.,trim(label)//'_rho',ikp)
  audit_nemaxx=DM_nemax
  call Diagonalize_DMFT_WEIGHTS(rhoEValues,rhoEVectors,nbands,DM_nemin, &
       audit_nemaxx,.False.,trim(label)//'_rhoE',ikp)
  nbandsx=audit_nemaxx-DM_nemin+1
  weights=0.d0
  weights(DM_nemin:audit_nemaxx)=rhoValues(:nbandsx)

  before(:,:,1)=stress_val_corr
  before(:,:,2)=stress_kin_mt
  before(:,:,3)=stress_kin_int
  before(:,:,4)=stress_apw
  before(:,:,5)=stress_nsp_match
  before(:,:,6)=stress_nsp_radial
  before(:,:,7)=stress_nsp_angular
  before(:,:,8)=stress_rela_corr
  before(:,:,9)=stress_rela_surf
  before_radial=(/radial_fd_rho_hsph,radial_fd_rho_hnsp, &
       radial_fd_rhoe_overlap/)

  call cmp_dmft_stress(rhoVectors,rhoEVectors,weights,rhoEValues, &
       nemin,nemax,DM_nemin,DM_nemax,audit_nemaxx,nbands,nbandsx, &
       n0,nnlo,time_bl,time_bl_w,time_lo,time_lo_w,0)

  diagnostic(:,:,1)=diagnostic(:,:,1)+stress_val_corr-before(:,:,1)
  diagnostic(:,:,2)=diagnostic(:,:,2)+stress_kin_mt-before(:,:,2)
  diagnostic(:,:,3)=diagnostic(:,:,3)+stress_kin_int-before(:,:,3)
  diagnostic(:,:,4)=diagnostic(:,:,4)+stress_apw-before(:,:,4)

  stress_val_corr=before(:,:,1)
  stress_kin_mt=before(:,:,2)
  stress_kin_int=before(:,:,3)
  stress_apw=before(:,:,4)
  stress_nsp_match=before(:,:,5)
  stress_nsp_radial=before(:,:,6)
  stress_nsp_angular=before(:,:,7)
  stress_rela_corr=before(:,:,8)
  stress_rela_surf=before(:,:,9)
  ! The matrix-decomposition pass is observational only.  The radial-response
  ! accumulators are global for the physical main/window/normalization passes,
  ! so restore them together with the ordinary stress tensors.  Omitting this
  ! restoration made an opt-in radial diagnostic depend on whether an
  ! unrelated matrix-decomposition audit was also enabled.
  radial_fd_rho_hsph=before_radial(1)
  radial_fd_rho_hnsp=before_radial(2)
  radial_fd_rhoe_overlap=before_radial(3)

  deallocate(rhoVectors,rhoEVectors,rhoValues,rhoEValues,weights)
end subroutine cmp_dmft_stress_matrix_diagnostic


! On-shell diagonal reference used by the paired reconstruction.  The lower
! bands retain their ordinary irreducible-k/spin weight; inside the DMFT
! window the occupation is diag(rho), and rhoE_nn=epsilon_n rho_nn.  Thus the
! unchanged WIEN/Blaha kernel is used only in the situation for which it was
! derived.  All global accumulators are restored before return.
subroutine cmp_dmft_stress_diagonal_reference(rho,wkp,diagnostic,nemin, &
     nemax,DM_nemin,DM_nemax,nbands,n0,nnlo,time_bl,time_bl_w,time_lo, &
     time_lo_w)
  use param,only: nume
  use xa,only: E
  use dmft_stress,only: stress_val_corr,stress_kin_mt,stress_apw, &
       stress_kin_int,stress_nsp_match,stress_nsp_radial, &
       stress_nsp_angular,stress_rela_corr,stress_rela_surf
  use stress_radial_fd,only: radial_fd_rho_hsph,radial_fd_rho_hnsp, &
       radial_fd_rhoe_overlap
  implicit none
  integer,intent(in) :: nemin,nemax,DM_nemin,DM_nemax,nbands,n0,nnlo
  complex*16,intent(in) :: rho(nbands,nbands),wkp
  real*8,intent(inout) :: diagnostic(3,3,4)
  real*8,intent(inout) :: time_bl,time_bl_w,time_lo,time_lo_w
  complex*16,allocatable :: identity(:,:)
  real*8,allocatable :: weights(:),energy_weights(:)
  real*8 :: before(3,3,9),before_radial(3)
  integer :: i,DM_nemaxx

  allocate(identity(nbands,nbands),weights(nume),energy_weights(nbands))
  identity=0.d0
  do i=1,nbands
     identity(i,i)=1.d0
  enddo
  weights=0.d0
  weights(nemin:DM_nemin-1)=dble(wkp)
  do i=1,nbands
     weights(DM_nemin+i-1)=dble(rho(i,i))
     energy_weights(i)=E(DM_nemin+i-1)*dble(rho(i,i))
  enddo
  DM_nemaxx=DM_nemax

  before(:,:,1)=stress_val_corr
  before(:,:,2)=stress_kin_mt
  before(:,:,3)=stress_kin_int
  before(:,:,4)=stress_apw
  before(:,:,5)=stress_nsp_match
  before(:,:,6)=stress_nsp_radial
  before(:,:,7)=stress_nsp_angular
  before(:,:,8)=stress_rela_corr
  before(:,:,9)=stress_rela_surf
  before_radial=(/radial_fd_rho_hsph,radial_fd_rho_hnsp, &
       radial_fd_rhoe_overlap/)

  call cmp_dmft_stress(identity,identity,weights,energy_weights,nemin, &
       nemax,DM_nemin,DM_nemax,DM_nemaxx,nbands,nbands,n0,nnlo, &
       time_bl,time_bl_w,time_lo,time_lo_w,0)
  diagnostic(:,:,1)=diagnostic(:,:,1)+stress_val_corr-before(:,:,1)
  diagnostic(:,:,2)=diagnostic(:,:,2)+stress_kin_mt-before(:,:,2)
  diagnostic(:,:,3)=diagnostic(:,:,3)+stress_kin_int-before(:,:,3)
  diagnostic(:,:,4)=diagnostic(:,:,4)+stress_apw-before(:,:,4)

  stress_val_corr=before(:,:,1)
  stress_kin_mt=before(:,:,2)
  stress_kin_int=before(:,:,3)
  stress_apw=before(:,:,4)
  stress_nsp_match=before(:,:,5)
  stress_nsp_radial=before(:,:,6)
  stress_nsp_angular=before(:,:,7)
  stress_rela_corr=before(:,:,8)
  stress_rela_surf=before(:,:,9)
  radial_fd_rho_hsph=before_radial(1)
  radial_fd_rho_hnsp=before_radial(2)
  radial_fd_rhoe_overlap=before_radial(3)
  deallocate(identity,weights,energy_weights)
end subroutine cmp_dmft_stress_diagonal_reference


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
     sigma_stress,sigma_match,sigma_radial,sigma_angular, &
     sigma_interstitial,sigma_coeff_fd_mt,sigma_coeff_analytic_mt, &
     nemin,nemax,DM_nemin,DM_nemax,DM_nemaxx, &
     nbands,nbandsx,n0,nnlo,time_bl,time_bl_w,time_lo,time_lo_w)
  use param, only: lmax2,lomax,nume,nloat
  use com, only: nat
  use dmfts, only: iso
  use structure, only: mult,rmt
  use xa, only: E,fj,dfj
  use xa3, only: bk3
  use w_atpar, only: retrieve_w_atpar
  use alm_blm_clm, only: cmp_alm_blm_clm,cmp_alm_blm_clm_strain
  use stress_kb, only: load_stress_integrals
  use stress_pw, only: prepare_stress_pw
  use dmft_stress, only: stress_coeff_fd_enabled, &
       stress_coeff_analytic_enabled
  use stress_g2_matrices,only: evaluate_g2
  implicit none
  integer,intent(in) :: nemin,nemax,DM_nemin,DM_nemax,DM_nemaxx
  integer,intent(in) :: nbands,nbandsx,n0,nnlo
  complex*16,intent(in) :: moment_vectors(nbands,nbands)
  real*8,intent(in) :: moment_values(nbands)
  real*8,intent(inout) :: sigma_stress(3,3),sigma_match(3,3), &
       sigma_radial(3,3),sigma_angular(3,3),sigma_interstitial(3,3)
  real*8,intent(inout) :: sigma_coeff_fd_mt(3,3)
  real*8,intent(inout) :: sigma_coeff_analytic_mt(3,3)
  real*8,intent(inout) :: time_bl,time_bl_w,time_lo,time_lo_w

  integer :: lmfull,lmlo,jatom,lfirst,lmmax,latom,mu,is,num
  integer :: alpha,beta,isize
  real*8 :: esave(nume),eweight(nume),buf(3,3)
  real*8 :: buf_match(3,3),buf_radial(3,3),buf_angular(3,3)
  real*8 :: buf_interstitial(3,3)
  complex*16,allocatable :: alm(:,:),blm(:,:)
  complex*16,allocatable :: clm_compact(:,:,:),clm_stress(:,:,:)
  complex*16,allocatable :: aalm(:,:,:),bblm(:,:,:),cclm(:,:,:,:)
  complex*16,allocatable :: ke0(:,:,:),ke1(:,:,:)
  complex*16,allocatable :: kmatch0(:,:,:),kmatch1(:,:,:), &
       kradial0(:,:,:),kradial1(:,:,:),kangular0(:,:,:), &
       kangular1(:,:,:)
  complex*16,allocatable :: nsp0(:,:,:),nsp1(:,:,:),nsp2(:,:,:)
  complex*16,allocatable :: salm(:,:,:),sblm(:,:,:),sclm(:,:,:,:)
  real*8 :: buf_coeff_analytic(3,3)

  lmfull=(lmax2+1)*(lmax2+1)
  lmlo=(lomax+1)*(lomax+1)
  isize=n0-nnlo
  allocate(alm(lmfull,nume),blm(lmfull,nume))
  allocate(clm_compact(lmlo,nume,nloat),clm_stress(lmfull,nume,nloat))
  allocate(aalm(1,1,1),bblm(1,1,1),cclm(1,1,1,1))
  allocate(ke0(3,3,nume),ke1(3,3,nume))
  allocate(kmatch0(3,3,nume),kmatch1(3,3,nume), &
       kradial0(3,3,nume),kradial1(3,3,nume), &
       kangular0(3,3,nume),kangular1(3,3,nume))
  allocate(nsp0(3,3,nume),nsp1(3,3,nume),nsp2(3,3,nume))
  allocate(salm(lmfull,nume,6), &
       sblm(lmfull,nume,6),sclm(lmlo,nume,nloat,6))

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

           call cmp_alm_blm_clm_strain(salm,sblm,sclm,alm,blm, &
                   clm_compact,moment_vectors,latom,lfirst,jatom,n0, &
                   nnlo,is,DM_nemin,DM_nemax,DM_nemaxx,nbands,nbandsx,first_band=nemin)
           if (stress_coeff_analytic_enabled) then
              call sigma_overlap_coeff_analytic_atom(salm,sblm,sclm, &
                   alm,blm,clm_compact,eweight,buf_coeff_analytic, &
                   jatom,DM_nemin,DM_nemaxx)
              sigma_coeff_analytic_mt=sigma_coeff_analytic_mt+ &
                   buf_coeff_analytic
           endif

           ! The isolated overlap derivative has an interstitial
           ! step-function piece which is absent from the on-shell Blaha IBS
           ! kernel.  It is global rather than atom-resolved, so evaluate it
           ! exactly once per spin component and k point after A has been
           ! transformed to the M_Sigma pseudo-band basis.
           if (jatom.eq.1 .and. mu.eq.1) then
              call val_overlap_interstitial_stress(isize,buf_interstitial, &
                   DM_nemin,DM_nemaxx,eweight)
              sigma_interstitial=sigma_interstitial+buf_interstitial
           endif

           esave=E
           E(nemin:DM_nemaxx)=0.d0
           call evaluate_g2(jatom,nemin,DM_nemaxx,alm,blm,clm_compact, &
                salm,sblm,sclm,E,ke0,nsp0,nsp1,nsp2,kmatch0, &
                kradial0,kangular0,.false.)
           E(nemin:DM_nemaxx)=1.d0
           call evaluate_g2(jatom,nemin,DM_nemaxx,alm,blm,clm_compact, &
                salm,sblm,sclm,E,ke1,nsp0,nsp1,nsp2,kmatch1, &
                kradial1,kangular1,.false.)
           E=esave

           buf=0.d0; buf_match=0.d0; buf_radial=0.d0; buf_angular=0.d0
           do num=DM_nemin,DM_nemaxx
              do beta=1,3
                 do alpha=1,3
                    buf(alpha,beta)=buf(alpha,beta)+2.d0*eweight(num)* &
                         dble(ke0(alpha,beta,num)-ke1(alpha,beta,num))
                    buf_match(alpha,beta)=buf_match(alpha,beta)+ &
                         2.d0*eweight(num)*dble( &
                         kmatch0(alpha,beta,num)-kmatch1(alpha,beta,num))
                    buf_radial(alpha,beta)=buf_radial(alpha,beta)+ &
                         2.d0*eweight(num)*dble( &
                         kradial0(alpha,beta,num)-kradial1(alpha,beta,num))
                    buf_angular(alpha,beta)=buf_angular(alpha,beta)+ &
                         2.d0*eweight(num)*dble( &
                         kangular0(alpha,beta,num)-kangular1(alpha,beta,num))
                 enddo
              enddo
           enddo
           call symmetry_stress_rotij(latom,jatom,buf)
           call symmetry_stress_rotij(latom,jatom,buf_match)
           call symmetry_stress_rotij(latom,jatom,buf_radial)
           call symmetry_stress_rotij(latom,jatom,buf_angular)
           sigma_stress=sigma_stress+buf
           sigma_match=sigma_match+buf_match
           sigma_radial=sigma_radial+buf_radial
           sigma_angular=sigma_angular+buf_angular

           ! Independent fixed-coefficient oracle for the muffin-tin part of
           ! 2 Tr[M_Sigma O_,epsilon].  This deliberately bypasses CAPD and
           ! the K(0)-K(1) rearrangement: it perturbs the geometry entering
           ! the actual ALM/BLM/CLM builder while holding eigenvectors and
           ! radial functions fixed.  The diagnostic is expensive and is
           ! enabled only by EDMFT_STRESS_COEFF_FD.
           if (stress_coeff_fd_enabled) call sigma_overlap_coeff_fd_atom( &
                moment_vectors,eweight,buf,latom,lfirst,jatom,is,nemin, &
                nemax,DM_nemin,DM_nemax,DM_nemaxx,nbands,nbandsx,n0, &
                nnlo,time_bl,time_bl_w,time_lo,time_lo_w)
           if (stress_coeff_fd_enabled) &
                sigma_coeff_fd_mt=sigma_coeff_fd_mt+buf
        enddo
     enddo
  enddo

  deallocate(alm,blm,clm_compact,clm_stress,aalm,bblm,cclm, &
       ke0,ke1,kmatch0,kmatch1,kradial0,kradial1,kangular0, &
       kangular1,nsp0,nsp1,nsp2)
  if (allocated(salm)) deallocate(salm,sblm,sclm)
end subroutine cmp_dmft_overlap_stress_diag


! Contract the literal physical matching-coefficient derivative with the
! fixed radial overlap metric.  SALM/SBLM/SCLM are one-sided derivatives, so
! the factor two supplies the Hermitian bra+ket partner.  This is the
! analytic counterpart of SIGMA_OVERLAP_COEFF_FD_MT below.
subroutine sigma_overlap_coeff_analytic_atom(salm,sblm,sclm,alm,blm, &
     clm,eweight,buf,jatom,DM_nemin,DM_nemaxx)
  use param,only: lmax2,lomax,nume,nloat
  use lo,only: ilo,loor
  use w_atpar,only: ri_mat
  implicit none
  complex*16,intent(in) :: salm((lmax2+1)*(lmax2+1),nume,6), &
       sblm((lmax2+1)*(lmax2+1),nume,6), &
       sclm((lomax+1)*(lomax+1),nume,nloat,6), &
       alm((lmax2+1)*(lmax2+1),nume), &
       blm((lmax2+1)*(lmax2+1),nume), &
       clm((lomax+1)*(lomax+1),nume,nloat)
  real*8,intent(in) :: eweight(nume)
  real*8,intent(out) :: buf(3,3)
  integer,intent(in) :: jatom,DM_nemin,DM_nemaxx
  integer,parameter :: strain_a(6)=(/1,2,3,2,1,1/), &
       strain_b(6)=(/1,2,3,3,3,2/)
  integer :: iq,ia,ib,num,l,m,ly,jlo,ip,ip2
  complex*16 :: coeff(2+nloat),dcoeff(2+nloat),metric_value

  buf=0.d0
  do iq=1,6
     ia=strain_a(iq); ib=strain_b(iq)
     do num=DM_nemin,DM_nemaxx
        do l=0,lmax2
           do m=-l,l
              ly=l*(l+1)+m+1
              coeff=0.d0; dcoeff=0.d0
              coeff(1)=alm(ly,num); coeff(2)=blm(ly,num)
              dcoeff(1)=salm(ly,num,iq); dcoeff(2)=sblm(ly,num,iq)
              if (l.le.lomax) then
                 do jlo=1,ilo(l)
                    if (loor(jlo,l)) then
                       coeff(2+jlo)=clm(ly,num,jlo)
                       dcoeff(2+jlo)=sclm(ly,num,jlo,iq)
                    endif
                 enddo
              endif
              metric_value=0.d0
              do ip2=1,2+nloat
                 do ip=1,2+nloat
                    metric_value=metric_value+conjg(dcoeff(ip))* &
                         ri_mat(ip,ip2,l,jatom)*coeff(ip2)
                 enddo
              enddo
              buf(ia,ib)=buf(ia,ib)+2.d0*eweight(num)* &
                   dble(metric_value)
           enddo
        enddo
     enddo
     buf(ib,ia)=buf(ia,ib)
  enddo
end subroutine sigma_overlap_coeff_analytic_atom


!=============================================================================
! Paired Force1-like coefficient response for arbitrary rho/rhoE matrices
!=============================================================================
!
! This diagnostic keeps rho and rhoE in one common band basis and evaluates
!
!   STR_VAL^match = -2 Re Tr[(dC)^dagger(H_sph C rho-R C rhoE)].
!
! Writing K=2 H_sph gives the implementation below,
!
!   STR_VAL^match = -Re Tr[(dC)^dagger(K C rho-2 R C rhoE)].
!
! It is the homogeneous-strain analogue of Force1_DMFT.  Unlike two
! independent pseudo-band diagonalizations, the large H/O cancellation is
! performed before summation.  The routine contains only the spherical
! coefficient/matching channel; operator-metric, interstitial, surface and
! nonspherical terms remain separate.
subroutine cmp_dmft_paired_matching_matrix(rho,rhoE,band_nemin,band_nemax, &
     n0,nnlo,paired_tensor,complete_tensor,kinetic_tensor, &
     kinetic_complete_tensor,time_bl,time_bl_w,time_lo,time_lo_w)
  use param,only: lmax2,lomax,nume,nloat
  use com,only: nat
  use dmfts,only: iso
  use structure,only: mult,rmt
  use w_atpar,only: retrieve_w_atpar,Vts,lv,lpv,mv,mpv,Nhmx, &
       ul_Rmt,dul_Rmt
  use alm_blm_clm,only: cmp_alm_blm_clm,cmp_alm_blm_clm_strain
  implicit none
  integer,intent(in) :: band_nemin,band_nemax,n0,nnlo
  complex*16,intent(in) :: rho(band_nemax-band_nemin+1, &
       band_nemax-band_nemin+1),rhoE(band_nemax-band_nemin+1, &
       band_nemax-band_nemin+1)
  real*8,intent(inout) :: paired_tensor(3,3)
  real*8,intent(inout) :: complete_tensor(3,3)
  real*8,intent(inout) :: kinetic_tensor(3,3)
  real*8,intent(inout) :: kinetic_complete_tensor(3,3)
  real*8,intent(inout) :: time_bl,time_bl_w,time_lo,time_lo_w
  integer :: nb,lmfull,lmlo,jatom,lfirst,lmmax,latom,mu,is,i
  complex*16,allocatable :: identity(:,:),alm(:,:),blm(:,:),clm(:,:,:), &
       salm(:,:,:),sblm(:,:,:),sclm(:,:,:,:),aalm(:,:,:), &
       bblm(:,:,:),cclm(:,:,:,:)
  real*8 :: buf(3,3),buf_kinetic(3,3),buf_surface(3,3),buf_nsp(3,3), &
       buf_interstitial(3,3)

  nb=band_nemax-band_nemin+1
  lmfull=(lmax2+1)*(lmax2+1)
  lmlo=(lomax+1)*(lomax+1)
  allocate(identity(nb,nb)); identity=0.d0
  do i=1,nb
     identity(i,i)=1.d0
  enddo
  allocate(alm(lmfull,nume),blm(lmfull,nume),clm(lmlo,nume,nloat), &
       salm(lmfull,nume,6),sblm(lmfull,nume,6), &
       sclm(lmlo,nume,nloat,6))
  allocate(aalm(1,1,1),bblm(1,1,1),cclm(1,1,1,1))

  do jatom=1,nat
     call retrieve_w_atpar(jatom,lfirst,lmmax)
     latom=lfirst-1
     do mu=1,mult(jatom)
        latom=latom+1
        do is=1,iso
           call cmp_alm_blm_clm(alm,blm,clm,aalm,bblm,cclm,identity, &
                latom,lfirst,jatom,band_nemin,band_nemax,n0,nnlo,is, &
                band_nemin,band_nemax,band_nemax,nb,nb,.false., &
                time_bl,time_bl_w,time_lo,time_lo_w)
           call cmp_alm_blm_clm_strain(salm,sblm,sclm,alm,blm,clm, &
                identity,latom,lfirst,jatom,n0,nnlo,is,band_nemin, &
                band_nemax,band_nemax,nb,nb)
           call paired_matching_contract(salm,sblm,sclm,alm,blm,clm, &
                rho,rhoE,buf,buf_kinetic,jatom,band_nemin,band_nemax)
           paired_tensor=paired_tensor+buf
           kinetic_tensor=kinetic_tensor+buf_kinetic
           call paired_surface_contract(salm,sblm,sclm,alm,blm,clm,rho, &
                buf_surface,jatom,rmt(jatom),ul_Rmt(:,:,jatom), &
                dul_Rmt(:,:,jatom),band_nemin,band_nemax)
           call paired_nonspherical_contract(salm,sblm,sclm,alm,blm,clm, &
                rho,buf_nsp,jatom,band_nemin,band_nemax)
           buf_interstitial=0.d0
           if (jatom.eq.1 .and. mu.eq.1) call paired_interstitial_contract( &
                rho,rhoE,buf_interstitial,band_nemin,band_nemax,n0-nnlo)
           complete_tensor=complete_tensor+buf+buf_surface+buf_nsp+ &
                buf_interstitial
           kinetic_complete_tensor=kinetic_complete_tensor+buf_kinetic+ &
                buf_surface+buf_interstitial
        enddo
     enddo
  enddo

  deallocate(identity,alm,blm,clm,salm,sblm,sclm,aalm,bblm,cclm)
end subroutine cmp_dmft_paired_matching_matrix


subroutine paired_matching_contract(salm,sblm,sclm,alm,blm,clm,rho, &
     rhoE,buf,buf_kinetic,jatom,band_nemin,band_nemax)
  use param,only: lmax2,lomax,nume,nloat
  use atspdt,only: el
  use lo,only: ilo,loor,elo
  use w_atpar,only: ri_mat
  implicit none
  integer,intent(in) :: jatom,band_nemin,band_nemax
  complex*16,intent(in) :: salm((lmax2+1)*(lmax2+1),nume,6), &
       sblm((lmax2+1)*(lmax2+1),nume,6), &
       sclm((lomax+1)*(lomax+1),nume,nloat,6), &
       alm((lmax2+1)*(lmax2+1),nume), &
       blm((lmax2+1)*(lmax2+1),nume), &
       clm((lomax+1)*(lomax+1),nume,nloat), &
       rho(band_nemax-band_nemin+1,band_nemax-band_nemin+1), &
       rhoE(band_nemax-band_nemin+1,band_nemax-band_nemin+1)
  real*8,intent(out) :: buf(3,3),buf_kinetic(3,3)
  integer,parameter :: strain_a(6)=(/1,2,3,2,1,1/), &
       strain_b(6)=(/1,2,3,3,3,2/)
  integer :: nb,nr,l,m,ly,jlo,jlop,n,iq,ia,ib
  real*8 :: Rmat(2+nloat,2+nloat),Kmat(2+nloat,2+nloat), &
       Vmat(2+nloat,2+nloat),Tmat(2+nloat,2+nloat),zene
  complex*16,allocatable :: C(:,:),S(:,:,:),HC(:,:),TC(:,:),OC(:,:), &
       B(:,:),Bkinetic(:,:)

  nb=band_nemax-band_nemin+1
  allocate(C(2+nloat,nb),S(2+nloat,nb,6),HC(2+nloat,nb), &
       TC(2+nloat,nb),OC(2+nloat,nb),B(2+nloat,nb), &
       Bkinetic(2+nloat,nb))
  buf=0.d0; buf_kinetic=0.d0

  do l=0,lmax2
     nr=2
     if (l.le.lomax) nr=2+ilo(l)
     Rmat=0.d0; Kmat=0.d0
     Rmat(1:nr,1:nr)=ri_mat(1:nr,1:nr,l,jatom)
     Kmat(1,1)=2.d0*el(l)
     Kmat(1,2)=1.d0; Kmat(2,1)=1.d0
     Kmat(2,2)=2.d0*el(l)*Rmat(2,2)
     if (l.le.lomax) then
        do jlo=1,ilo(l)
           zene=el(l)+elo(l,jlo)
           Kmat(1,2+jlo)=zene*Rmat(1,2+jlo)
           Kmat(2+jlo,1)=Kmat(1,2+jlo)
           Kmat(2,2+jlo)=Rmat(1,2+jlo)+zene*Rmat(2,2+jlo)
           Kmat(2+jlo,2)=Kmat(2,2+jlo)
           do jlop=1,ilo(l)
              Kmat(2+jlo,2+jlop)=(elo(l,jlo)+elo(l,jlop))* &
                   Rmat(2+jlo,2+jlop)
           enddo
        enddo
     endif
     call paired_spherical_potential_matrix(Vmat,l,jatom)
     ! KMAT is twice the spherical Hamiltonian volume matrix because the
     ! Hermitian bra partner is supplied by the outer real part.  VMAT is in
     ! Rydbergs, so twice the kinetic volume matrix is KMAT-2*VMAT.
     Tmat=Kmat-2.d0*Vmat

     do m=-l,l
        ly=l*(l+1)+m+1
        C=0.d0; S=0.d0
        do n=1,nb
           C(1,n)=alm(ly,band_nemin+n-1)
           C(2,n)=blm(ly,band_nemin+n-1)
           S(1,n,:)=salm(ly,band_nemin+n-1,:)
           S(2,n,:)=sblm(ly,band_nemin+n-1,:)
           if (l.le.lomax) then
              do jlo=1,ilo(l)
                 if (loor(jlo,l)) then
                    C(2+jlo,n)=clm(ly,band_nemin+n-1,jlo)
                    S(2+jlo,n,:)=sclm(ly,band_nemin+n-1,jlo,:)
                 endif
              enddo
           endif
        enddo
        HC=matmul(Kmat,C)
        TC=matmul(Tmat,C)
        OC=matmul(Rmat,C)
        B=matmul(HC,rho)-2.d0*matmul(OC,rhoE)
        Bkinetic=matmul(TC,rho)-2.d0*matmul(OC,rhoE)
        do iq=1,6
           ia=strain_a(iq); ib=strain_b(iq)
           buf(ia,ib)=buf(ia,ib)-dble(sum(conjg(S(:,:,iq))*B))
           buf(ib,ia)=buf(ia,ib)
           buf_kinetic(ia,ib)=buf_kinetic(ia,ib)- &
                dble(sum(conjg(S(:,:,iq))*Bkinetic))
           buf_kinetic(ib,ia)=buf_kinetic(ia,ib)
        enddo
     enddo
  enddo
  deallocate(C,S,HC,TC,OC,B,Bkinetic)
end subroutine paired_matching_contract


! Spherical potential matrix in the primitive radial basis used by ALM,
! BLM and CLM.  W_VR stores r*V_sph/2 after ATPARN's Ry-to-Hartree
! conversion; consequently 2*W_VR/r is the physical potential in Rydbergs.
! Subtracting this matrix from the radial-equation Hamiltonian is preferable
! to differentiating a standalone kinetic differential operator: it keeps
! the scalar-relativistic and local-orbital conventions exactly identical to
! the generalized eigenproblem.
subroutine paired_spherical_potential_matrix(Vmat,l,jatom)
  use param,only: nrad,lmax2,nloat
  use com,only: rel
  use structure,only: jri,dx,r0
  use lo,only: ilo
  use w_atpar,only: w_RRAD1,w_RRAD2,w_RADE1,w_RADE2,w_a1lo,w_b1lo,w_VR
  implicit none
  integer,intent(in) :: l,jatom
  real*8,intent(out) :: Vmat(2+nloat,2+nloat)
  integer :: nr,ir,i,j
  real*8 :: radius,potential
  real*8 :: radial1(nrad,2+nloat),radial2(nrad,2+nloat)
  real*8 :: ket1(nrad),ket2(nrad)
  interface
     function RINT13(REL,A,B,X,Y,NRAD_,DX_,JRI_,R0_) result(S)
       logical,intent(in) :: REL
       integer,intent(in) :: NRAD_,JRI_
       real*8,intent(in) :: A(NRAD_),B(NRAD_),X(NRAD_),Y(NRAD_),DX_,R0_
       real*8 :: S
     end function RINT13
  end interface

  nr=2+ilo(l)
  radial1=0.d0; radial2=0.d0; Vmat=0.d0
  radial1(:,1)=w_RRAD1(:,l,jatom)
  radial2(:,1)=w_RRAD2(:,l,jatom)
  radial1(:,2)=w_RADE1(:,l,jatom)
  radial2(:,2)=w_RADE2(:,l,jatom)
  do i=1,ilo(l)
     radial1(:,2+i)=w_a1lo(:,i,l,jatom)
     radial2(:,2+i)=w_b1lo(:,i,l,jatom)
  enddo
  do j=1,nr
     ket1=0.d0; ket2=0.d0
     do ir=1,jri(jatom)
        radius=r0(jatom)*exp(dx(jatom)*dble(ir-1))
        potential=2.d0*w_VR(ir,jatom)/radius
        ket1(ir)=potential*radial1(ir,j)
        ket2(ir)=potential*radial2(ir,j)
     enddo
     do i=1,nr
        Vmat(i,j)=RINT13(rel,radial1(:,i),radial2(:,i),ket1,ket2, &
             nrad,dx(jatom),jri(jatom),r0(jatom))
     enddo
  enddo
  ! The central potential is Hermitian.  Explicit symmetrization removes the
  ! insignificant quadrature asymmetry before it is amplified by the large
  ! H-V cancellation in deep core-like radial channels.
  Vmat(1:nr,1:nr)=0.5d0*(Vmat(1:nr,1:nr)+ &
       transpose(Vmat(1:nr,1:nr)))
end subroutine paired_spherical_potential_matrix


! The radial Hamiltonian matrix used by PAIRED_MATCHING is the volume part
! obtained from the radial Schroedinger equation.  LAPW1's symmetric
! Hamiltonian has in addition the kinetic surface matrix
!
!   Hsurf_ij = R_MT^2 [u_i(R) u'_j(R)+u'_i(R) u_j(R)]/2.
!
! This is the homogeneous-strain counterpart of Force3.  The sphere radius
! and radial functions are fixed in the present convention; only the LAPW
! matching coefficients are differentiated.
subroutine paired_surface_contract(salm,sblm,sclm,alm,blm,clm,rho,buf, &
     jatom,radius,ul,dul,band_nemin,band_nemax)
  use param,only: lmax2,lomax,nume,nloat
  use lo,only: ilo,loor
  implicit none
  integer,intent(in) :: jatom,band_nemin,band_nemax
  real*8,intent(in) :: radius,ul(nloat+2,0:lmax2), &
       dul(nloat+2,0:lmax2)
  complex*16,intent(in) :: salm((lmax2+1)*(lmax2+1),nume,6), &
       sblm((lmax2+1)*(lmax2+1),nume,6), &
       sclm((lomax+1)*(lomax+1),nume,nloat,6), &
       alm((lmax2+1)*(lmax2+1),nume), &
       blm((lmax2+1)*(lmax2+1),nume), &
       clm((lomax+1)*(lomax+1),nume,nloat), &
       rho(band_nemax-band_nemin+1,band_nemax-band_nemin+1)
  real*8,intent(out) :: buf(3,3)
  integer,parameter :: strain_a(6)=(/1,2,3,2,1,1/), &
       strain_b(6)=(/1,2,3,3,3,2/)
  integer :: nb,nr,l,m,ly,jlo,n,iq,ia,ib,i,j
  real*8 :: hsurf(2+nloat,2+nloat)
  complex*16,allocatable :: C(:,:),S(:,:,:),HC(:,:),B(:,:)

  nb=band_nemax-band_nemin+1
  allocate(C(2+nloat,nb),S(2+nloat,nb,6),HC(2+nloat,nb), &
       B(2+nloat,nb))
  buf=0.d0
  do l=0,lmax2
     nr=2
     if (l.le.lomax) nr=2+ilo(l)
     hsurf=0.d0
     do j=1,nr
        do i=1,nr
           hsurf(i,j)=0.5d0*radius**2*(ul(i,l)*dul(j,l)+ &
                dul(i,l)*ul(j,l))
        enddo
     enddo
     do m=-l,l
        ly=l*(l+1)+m+1
        C=0.d0; S=0.d0
        do n=1,nb
           C(1,n)=alm(ly,band_nemin+n-1)
           C(2,n)=blm(ly,band_nemin+n-1)
           S(1,n,:)=salm(ly,band_nemin+n-1,:)
           S(2,n,:)=sblm(ly,band_nemin+n-1,:)
           if (l.le.lomax) then
              do jlo=1,ilo(l)
                 if (loor(jlo,l)) then
                    C(2+jlo,n)=clm(ly,band_nemin+n-1,jlo)
                    S(2+jlo,n,:)=sclm(ly,band_nemin+n-1,jlo,:)
                 endif
              enddo
           endif
        enddo
        HC=matmul(hsurf,C)
        B=matmul(HC,rho)
        do iq=1,6
           ia=strain_a(iq); ib=strain_b(iq)
           buf(ia,ib)=buf(ia,ib)-2.d0*dble(sum(conjg(S(:,:,iq))*B))
           buf(ib,ia)=buf(ia,ib)
        enddo
     enddo
  enddo
  deallocate(C,S,HC,B)
end subroutine paired_surface_contract


! Coefficient derivative of the complete nonspherical muffin-tin potential.
! Vts is the same full (u,udot,LO_1,...) matrix used by Force2.  The factor
! two supplies the Hermitian bra partner, so the result is exactly
! -Tr[rho d(C^+ V_ns C)] in WIEN's stored STR_VAL convention.
subroutine paired_nonspherical_contract(salm,sblm,sclm,alm,blm,clm,rho, &
     buf,jatom,band_nemin,band_nemax)
  use param,only: lmax2,lomax,nume,nloat
  use lo,only: ilo,loor
  use w_atpar,only: Vts,lv,lpv,mv,mpv,Nhmx
  implicit none
  integer,intent(in) :: jatom,band_nemin,band_nemax
  complex*16,intent(in) :: salm((lmax2+1)*(lmax2+1),nume,6), &
       sblm((lmax2+1)*(lmax2+1),nume,6), &
       sclm((lomax+1)*(lomax+1),nume,nloat,6), &
       alm((lmax2+1)*(lmax2+1),nume), &
       blm((lmax2+1)*(lmax2+1),nume), &
       clm((lomax+1)*(lomax+1),nume,nloat), &
       rho(band_nemax-band_nemin+1,band_nemax-band_nemin+1)
  real*8,intent(out) :: buf(3,3)
  integer,parameter :: strain_a(6)=(/1,2,3,2,1,1/), &
       strain_b(6)=(/1,2,3,3,3,2/)
  integer :: nb,ih,l,lp,ly,lpy,n,jlo,np,nq,iq,ia,ib
  complex*16,allocatable :: Cbra(:,:),Cket(:,:),Sket(:,:,:), &
       VC(:,:),dH(:,:)

  nb=band_nemax-band_nemin+1
  allocate(Cbra(2+nloat,nb),Cket(2+nloat,nb), &
       Sket(2+nloat,nb,6),VC(2+nloat,nb),dH(nb,nb))
  buf=0.d0
  do ih=1,Nhmx(jatom)
     lp=lpv(ih,jatom); l=lv(ih,jatom)
     if (lp.gt.lmax2 .or. l.gt.lmax2) cycle
     lpy=lp*(lp+1)+mpv(ih,jatom)+1
     ly=l*(l+1)+mv(ih,jatom)+1
     np=2
     if (lp.le.lomax) np=2+ilo(lp)
     nq=2
     if (l.le.lomax) nq=2+ilo(l)
     Cbra=0.d0; Cket=0.d0; Sket=0.d0
     do n=1,nb
        Cbra(1,n)=alm(lpy,band_nemin+n-1)
        Cbra(2,n)=blm(lpy,band_nemin+n-1)
        Cket(1,n)=alm(ly,band_nemin+n-1)
        Cket(2,n)=blm(ly,band_nemin+n-1)
        Sket(1,n,:)=salm(ly,band_nemin+n-1,:)
        Sket(2,n,:)=sblm(ly,band_nemin+n-1,:)
        if (lp.le.lomax) then
           do jlo=1,ilo(lp)
              if (loor(jlo,lp)) Cbra(2+jlo,n)= &
                   clm(lpy,band_nemin+n-1,jlo)
           enddo
        endif
        if (l.le.lomax) then
           do jlo=1,ilo(l)
              if (loor(jlo,l)) then
                 Cket(2+jlo,n)=clm(ly,band_nemin+n-1,jlo)
                 Sket(2+jlo,n,:)=sclm(ly,band_nemin+n-1,jlo,:)
              endif
           enddo
        endif
     enddo
     do iq=1,6
        VC=0.d0
        VC(1:np,:)=matmul(Vts(1:np,1:nq,ih,jatom), &
             Sket(1:nq,:,iq))
        dH=matmul(conjg(transpose(Cbra)),VC)
        ia=strain_a(iq); ib=strain_b(iq)
        buf(ia,ib)=buf(ia,ib)-2.d0*dble(sum(rho*transpose(dH)))
        buf(ib,ia)=buf(ia,ib)
     enddo
  enddo
  deallocate(Cbra,Cket,Sket,VC,dH)
end subroutine paired_nonspherical_contract


! Direct interstitial derivative in the same plane-wave representation used
! to build LAPW1's generalized eigenproblem:
!
!   O^I_GG' = Theta_G-G',
!   T^I_GG' = (k+G).(k+G') Theta_G-G'.
!
! WARPIN_G-G' is invariant in the frozen reduced-coordinate potential test,
! hence it has no explicit basis derivative.  Differentiating both the
! reciprocal vectors and the fixed-sphere step function gives the complete
! interstitial contribution without invoking an on-shell band energy.
subroutine paired_interstitial_contract(rho,rhoE,buf,band_nemin, &
     band_nemax,isize)
  use xa3,only: a,k3,bk3
  use structure,only: br1
  implicit none
  integer,intent(in) :: band_nemin,band_nemax,isize
  complex*16,intent(in) :: rho(band_nemax-band_nemin+1, &
       band_nemax-band_nemin+1),rhoE(band_nemax-band_nemin+1, &
       band_nemax-band_nemin+1)
  real*8,intent(out) :: buf(3,3)
  integer :: nb,npair,g,gp,index,alpha,beta
  integer,allocatable :: gx(:),gy(:),gz(:)
  real*8,allocatable :: p(:,:)
  complex*16,allocatable :: theta(:),dtheta(:,:,:),Avec(:,:), &
       Drho(:,:),DrhoE(:,:)
  complex*16 :: dO,dH
  real*8 :: kinetic,dkinetic

  nb=band_nemax-band_nemin+1
  npair=isize*isize
  allocate(gx(npair),gy(npair),gz(npair),theta(npair), &
       dtheta(3,3,npair),p(3,isize),Avec(isize,nb), &
       Drho(isize,isize),DrhoE(isize,isize))
  index=0
  do g=1,isize
     p(:,g)=matmul(br1,bk3(:,g))
     do gp=1,isize
        index=index+1
        gx(index)=k3(1,g)-k3(1,gp)
        gy(index)=k3(2,g)-k3(2,gp)
        gz(index)=k3(3,g)-k3(3,gp)
     enddo
  enddo
  call step_fun(gx,gy,gz,npair,theta)
  call step_fun_stress_derivative(gx,gy,gz,npair,dtheta)
  Avec=a(1:isize,band_nemin:band_nemax)
  Drho=matmul(matmul(Avec,rho),conjg(transpose(Avec)))
  DrhoE=matmul(matmul(Avec,rhoE),conjg(transpose(Avec)))
  buf=0.d0
  index=0
  do g=1,isize
     do gp=1,isize
        index=index+1
        kinetic=dot_product(p(:,g),p(:,gp))
        do alpha=1,3
           do beta=1,3
              dkinetic=-(p(alpha,g)*p(beta,gp)+ &
                   p(beta,g)*p(alpha,gp))
              dO=dtheta(alpha,beta,index)
              dH=dkinetic*theta(index)+kinetic*dO
              ! D(gp,g) is the coefficient-space density paired with the
              ! operator matrix element (g,gp).
              buf(alpha,beta)=buf(alpha,beta)- &
                   dble(dH*Drho(gp,g))+dble(dO*DrhoE(gp,g))
           enddo
        enddo
     enddo
  enddo
  deallocate(gx,gy,gz,theta,dtheta,p,Avec,Drho,DrhoE)
end subroutine paired_interstitial_contract


!=============================================================================
! Fixed-coefficient finite-difference oracle for the muffin-tin overlap
!=============================================================================
!
! For one atom/equivalent-atom/spin component, form
!
!   D^MT_q = 2 sum_a m_a d/d eps_q <psi_a|Theta_MT|psi_a>,
!
! where m_a are the signed eigenvalues of M_Sigma and psi_a are its
! pseudo-band eigenvectors.  Only BR1, VOL, |k+G| and the Bessel matching
! data are perturbed.  The plane-wave eigenvectors, radial basis functions,
! muffin-tin radii, and radial metric RI_MAT remain fixed.  Consequently this
! is exactly the object that the spherical overlap part of the analytic
! stress must reproduce under WIEN2k's fixed-sphere finite-strain convention.
!
! The routine changes no production accumulator and restores every global
! array byte-for-byte before returning.
subroutine sigma_overlap_coeff_fd_atom(moment_vectors,eweight,buf_fd, &
     latom,lfirst,jatom,is,nemin,nemax,DM_nemin,DM_nemax,DM_nemaxx, &
     nbands,nbandsx,n0,nnlo,time_bl,time_bl_w,time_lo,time_lo_w)
  use param, only: lmax2,lomax,nume,nloat
  use com, only: nat
  use structure, only: BR1,VOL,RMT,rotij
  use xa, only: fj,dfj
  use xa3, only: bk3,aK
  use lo, only: ilo,loor
  use w_atpar, only: ri_mat
  use alm_blm_clm, only: cmp_alm_blm_clm
  implicit none
  integer,intent(in) :: latom,lfirst,jatom,is,nemin,nemax,DM_nemin, &
       DM_nemax,DM_nemaxx,nbands,nbandsx,n0,nnlo
  complex*16,intent(in) :: moment_vectors(nbands,nbands)
  real*8,intent(in) :: eweight(nume)
  real*8,intent(out) :: buf_fd(3,3)
  real*8,intent(inout) :: time_bl,time_bl_w,time_lo,time_lo_w

  integer,parameter :: strain_a(6)=(/1,2,3,2,1,1/), &
       strain_b(6)=(/1,2,3,3,3,2/)
  real*8,parameter :: h=1.d-6
  integer :: isize,lmfull,lmlo,iq,ia,ib,i
  real*8 :: BR10(3,3),VOL0,Kn(3),plus_value,minus_value
  real*8,allocatable :: aK0(:),fj0(:,:,:),dfj0(:,:,:)
  complex*16,allocatable :: alm(:,:),blm(:,:),clm(:,:,:)
  complex*16,allocatable :: aalm(:,:,:),bblm(:,:,:),cclm(:,:,:,:)

  isize=n0-nnlo
  lmfull=(lmax2+1)*(lmax2+1)
  lmlo=(lomax+1)*(lomax+1)
  allocate(alm(lmfull,nume),blm(lmfull,nume),clm(lmlo,nume,nloat))
  allocate(aalm(1,1,1),bblm(1,1,1),cclm(1,1,1,1))
  allocate(aK0(isize),fj0(0:lmax2,isize,1), &
       dfj0(0:lmax2,isize,1))

  BR10=BR1
  VOL0=VOL
  aK0=aK(:isize)
  fj0(:,:,1)=fj(:,:isize,jatom)
  dfj0(:,:,1)=dfj(:,:isize,jatom)
  buf_fd=0.d0

  do iq=1,6
     call set_geometry(+h,iq)
     call cmp_alm_blm_clm(alm,blm,clm,aalm,bblm,cclm, &
          moment_vectors,latom,lfirst,jatom,nemin,nemax,n0,nnlo,is, &
          DM_nemin,DM_nemax,DM_nemaxx,nbands,nbandsx,.false., &
          time_bl,time_bl_w,time_lo,time_lo_w)
     plus_value=weighted_mt_overlap(alm,blm,clm)

     call set_geometry(-h,iq)
     call cmp_alm_blm_clm(alm,blm,clm,aalm,bblm,cclm, &
          moment_vectors,latom,lfirst,jatom,nemin,nemax,n0,nnlo,is, &
          DM_nemin,DM_nemax,DM_nemaxx,nbands,nbandsx,.false., &
          time_bl,time_bl_w,time_lo,time_lo_w)
     minus_value=weighted_mt_overlap(alm,blm,clm)

     ia=strain_a(iq); ib=strain_b(iq)
     buf_fd(ia,ib)=(plus_value-minus_value)/(2.d0*h)
     buf_fd(ib,ia)=buf_fd(ia,ib)
  enddo

  BR1=BR10
  VOL=VOL0
  aK(:isize)=aK0
  fj(:,:isize,jatom)=fj0(:,:,1)
  dfj(:,:isize,jatom)=dfj0(:,:,1)
  deallocate(alm,blm,clm,aalm,bblm,cclm,aK0,fj0,dfj0)

contains

  subroutine set_geometry(delta,iq_local)
    real*8,intent(in) :: delta
    integer,intent(in) :: iq_local
    integer :: ia_local,ib_local

    ia_local=strain_a(iq_local)
    ib_local=strain_b(iq_local)
    BR1=BR10
    if (ia_local.eq.ib_local) then
       ! F=I+delta e_a e_a is used exactly.  The distinction from the
       ! linearized inverse is O(h^2), but the exact form improves the
       ! central-difference conditioning at no cost.
       BR1(ia_local,:)=BR10(ia_local,:)/(1.d0+delta)
       VOL=VOL0*(1.d0+delta)
    else
       ! Off-diagonal q denotes eps_ab=eps_ba=delta/2.  The following is
       ! sufficient through first order, which is all the central derivative
       ! retains as h -> 0.
       BR1(ia_local,:)=BR10(ia_local,:)-0.5d0*delta*BR10(ib_local,:)
       BR1(ib_local,:)=BR10(ib_local,:)-0.5d0*delta*BR10(ia_local,:)
       VOL=VOL0*(1.d0-0.25d0*delta*delta)
    endif

    ! At a symmetry-equivalent atom the matching Bessel argument contains
    ! R_a(k+G).  At the unstrained structure its magnitude equals that of the
    ! representative vector; after the diagnostic strain the equality no
    ! longer holds, so construct the rotated magnitude explicitly.
    do i=1,isize
       Kn=matmul(BR1,matmul(rotij(:,:,latom),bk3(:,i)))
       aK(i)=sqrt(dot_product(Kn,Kn))
    enddo
    call HARMON(isize,aK(:isize),lmax2,fj(:,:isize,jatom), &
         dfj(:,:isize,jatom),RMT(jatom))
  end subroutine set_geometry

  function weighted_mt_overlap(alm_local,blm_local,clm_local) result(value)
    complex*16,intent(in) :: alm_local(lmfull,nume), &
         blm_local(lmfull,nume),clm_local(lmlo,nume,nloat)
    real*8 :: value,band_norm
    complex*16 :: coeff(2+nloat),metric_value
    integer :: num,l,m,ly,ip,ip2,jlo

    value=0.d0
    do num=DM_nemin,DM_nemaxx
       band_norm=0.d0
       do l=0,lmax2
          do m=-l,l
             ly=l*(l+1)+m+1
             coeff=0.d0
             coeff(1)=alm_local(ly,num)
             coeff(2)=blm_local(ly,num)
             if (l.le.lomax) then
                do jlo=1,ilo(l)
                   if (loor(jlo,l)) coeff(2+jlo)=clm_local(ly,num,jlo)
                enddo
             endif
             metric_value=0.d0
             do ip2=1,2+nloat
                do ip=1,2+nloat
                   metric_value=metric_value+conjg(coeff(ip))* &
                        ri_mat(ip,ip2,l,jatom)*coeff(ip2)
                enddo
             enddo
             band_norm=band_norm+dble(metric_value)
          enddo
       enddo
       ! BAND_NORM is the complete bra+ket norm at the endpoint.  Its
       ! central finite difference therefore already contains the 2*Re
       ! associated with the two one-sided coefficient derivatives.  Do not
       ! apply the explicit factor two used by the analytic one-sided CAPD
       ! contraction here.
       value=value+eweight(num)*band_norm
    enddo
  end function weighted_mt_overlap

end subroutine sigma_overlap_coeff_fd_atom
