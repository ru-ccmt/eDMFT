! State and radial-integral storage for the analytic LAPW stress tensor.
!
! WIEN2k lapw2 has an outer atom loop and may keep one atom's radial stress
! integrals in module variables.  dmft2 has an outer k-point loop.  Recomputing
! those radial integrals at every k point would be wasteful, so we compute them
! once in precmp_w_atpar, save one copy per inequivalent atom, and load the
! appropriate copy immediately before applying the unchanged WIEN2k kernels.

module stress_kb
  use param
  implicit none
  logical :: STRESS=.false.

  real*8,allocatable :: Jllp00(:,:), Jllp01(:,:), Jllp10(:,:), Jllp11(:,:)
  real*8,allocatable :: Illp00(:,:), Illp01(:,:), Illp10(:,:), Illp11(:,:)
  real*8,allocatable :: jr3upu_20(:,:,:),jr3upu_21(:,:,:),jr3upu_22(:,:,:,:)
  real*8,allocatable :: jr3upu_02(:,:,:),jr3upu_12(:,:,:)
  real*8,allocatable :: jlaplo1_02(:,:),jlaplo1_12(:,:)

  real*8,allocatable :: TC1S00(:,:),TC1S10(:,:),TC1S01(:,:),TC1S11(:,:)
  real*8,allocatable :: TC2S00(:,:),TC2S10(:,:),TC2S01(:,:),TC2S11(:,:)
  real*8,allocatable :: TC3S00(:,:),TC3S10(:,:),TC3S01(:,:),TC3S11(:,:)
  real*8,allocatable :: juddu_02(:,:,:),juddu_20(:,:,:),juddu_12(:,:,:)
  real*8,allocatable :: juddu_21(:,:,:),juddu_22(:,:,:,:)
  real*8,allocatable :: judu_02(:,:,:),judu_20(:,:,:),judu_12(:,:,:)
  real*8,allocatable :: judu_21(:,:,:),judu_22(:,:,:,:)
  real*8,allocatable :: juu_02(:,:,:),juu_20(:,:,:),juu_12(:,:,:)
  real*8,allocatable :: juu_21(:,:,:),juu_22(:,:,:,:)

  type stress_radial_cache
     real*8,allocatable :: J00(:,:),J01(:,:),J10(:,:),J11(:,:)
     real*8,allocatable :: I00(:,:),I01(:,:),I10(:,:),I11(:,:)
     real*8,allocatable :: r20(:,:,:),r21(:,:,:),r22(:,:,:,:),r02(:,:,:),r12(:,:,:)
     real*8,allocatable :: jl02(:,:),jl12(:,:)
     real*8,allocatable :: t100(:,:),t110(:,:),t101(:,:),t111(:,:)
     real*8,allocatable :: t200(:,:),t210(:,:),t201(:,:),t211(:,:)
     real*8,allocatable :: t300(:,:),t310(:,:),t301(:,:),t311(:,:)
     real*8,allocatable :: jd02(:,:,:),jd20(:,:,:),jd12(:,:,:),jd21(:,:,:),jd22(:,:,:,:)
     real*8,allocatable :: j102(:,:,:),j120(:,:,:),j112(:,:,:),j121(:,:,:),j122(:,:,:,:)
     real*8,allocatable :: j002(:,:,:),j020(:,:,:),j012(:,:,:),j021(:,:,:),j022(:,:,:,:)
  end type stress_radial_cache
  type(stress_radial_cache),allocatable :: radial_cache(:)

contains

  subroutine get_variable_integrals
    if (.not.allocated(Jllp00)) then
       allocate(Jllp00(0:lmax2,0:lmax2),Jllp01(0:lmax2,0:lmax2), &
                Jllp10(0:lmax2,0:lmax2),Jllp11(0:lmax2,0:lmax2), &
                Illp00(0:lmax2,0:lmax2),Illp01(0:lmax2,0:lmax2), &
                Illp10(0:lmax2,0:lmax2),Illp11(0:lmax2,0:lmax2))
       allocate(jr3upu_20(nloat,0:lomax,0:lmax2),jr3upu_21(nloat,0:lomax,0:lmax2), &
                jr3upu_22(nloat,nloat,0:lomax,0:lomax), &
                jr3upu_02(nloat,0:lomax,0:lmax2),jr3upu_12(nloat,0:lomax,0:lmax2), &
                jlaplo1_02(nloat,0:lomax),jlaplo1_12(nloat,0:lomax))
       allocate(TC1S00(0:lmax2,0:lmax2),TC1S10(0:lmax2,0:lmax2), &
                TC1S01(0:lmax2,0:lmax2),TC1S11(0:lmax2,0:lmax2), &
                TC2S00(0:lmax2,0:lmax2),TC2S10(0:lmax2,0:lmax2), &
                TC2S01(0:lmax2,0:lmax2),TC2S11(0:lmax2,0:lmax2), &
                TC3S00(0:lmax2,0:lmax2),TC3S10(0:lmax2,0:lmax2), &
                TC3S01(0:lmax2,0:lmax2),TC3S11(0:lmax2,0:lmax2))
       allocate(juddu_02(nloat,0:lomax,0:lmax2),juddu_20(nloat,0:lomax,0:lmax2), &
                juddu_12(nloat,0:lomax,0:lmax2),juddu_21(nloat,0:lomax,0:lmax2), &
                juddu_22(nloat,nloat,0:lomax,0:lomax), &
                judu_02(nloat,0:lomax,0:lmax2),judu_20(nloat,0:lomax,0:lmax2), &
                judu_12(nloat,0:lomax,0:lmax2),judu_21(nloat,0:lomax,0:lmax2), &
                judu_22(nloat,nloat,0:lomax,0:lomax), &
                juu_02(nloat,0:lomax,0:lmax2),juu_20(nloat,0:lomax,0:lmax2), &
                juu_12(nloat,0:lomax,0:lmax2),juu_21(nloat,0:lomax,0:lmax2), &
                juu_22(nloat,nloat,0:lomax,0:lomax))
    endif
  end subroutine get_variable_integrals

  subroutine init_stress_cache(nat)
    integer,intent(in) :: nat
    call get_variable_integrals
    if (.not.allocated(radial_cache)) allocate(radial_cache(nat))
  end subroutine init_stress_cache

  subroutine save_stress_integrals(jatom)
    integer,intent(in) :: jatom
    radial_cache(jatom)%J00=Jllp00; radial_cache(jatom)%J01=Jllp01
    radial_cache(jatom)%J10=Jllp10; radial_cache(jatom)%J11=Jllp11
    radial_cache(jatom)%I00=Illp00; radial_cache(jatom)%I01=Illp01
    radial_cache(jatom)%I10=Illp10; radial_cache(jatom)%I11=Illp11
    radial_cache(jatom)%r20=jr3upu_20; radial_cache(jatom)%r21=jr3upu_21
    radial_cache(jatom)%r22=jr3upu_22; radial_cache(jatom)%r02=jr3upu_02
    radial_cache(jatom)%r12=jr3upu_12
    radial_cache(jatom)%jl02=jlaplo1_02; radial_cache(jatom)%jl12=jlaplo1_12
    radial_cache(jatom)%t100=TC1S00; radial_cache(jatom)%t110=TC1S10
    radial_cache(jatom)%t101=TC1S01; radial_cache(jatom)%t111=TC1S11
    radial_cache(jatom)%t200=TC2S00; radial_cache(jatom)%t210=TC2S10
    radial_cache(jatom)%t201=TC2S01; radial_cache(jatom)%t211=TC2S11
    radial_cache(jatom)%t300=TC3S00; radial_cache(jatom)%t310=TC3S10
    radial_cache(jatom)%t301=TC3S01; radial_cache(jatom)%t311=TC3S11
    radial_cache(jatom)%jd02=juddu_02; radial_cache(jatom)%jd20=juddu_20
    radial_cache(jatom)%jd12=juddu_12; radial_cache(jatom)%jd21=juddu_21
    radial_cache(jatom)%jd22=juddu_22
    radial_cache(jatom)%j102=judu_02; radial_cache(jatom)%j120=judu_20
    radial_cache(jatom)%j112=judu_12; radial_cache(jatom)%j121=judu_21
    radial_cache(jatom)%j122=judu_22
    radial_cache(jatom)%j002=juu_02; radial_cache(jatom)%j020=juu_20
    radial_cache(jatom)%j012=juu_12; radial_cache(jatom)%j021=juu_21
    radial_cache(jatom)%j022=juu_22
  end subroutine save_stress_integrals

  subroutine load_stress_integrals(jatom)
    integer,intent(in) :: jatom
    Jllp00=radial_cache(jatom)%J00; Jllp01=radial_cache(jatom)%J01
    Jllp10=radial_cache(jatom)%J10; Jllp11=radial_cache(jatom)%J11
    Illp00=radial_cache(jatom)%I00; Illp01=radial_cache(jatom)%I01
    Illp10=radial_cache(jatom)%I10; Illp11=radial_cache(jatom)%I11
    jr3upu_20=radial_cache(jatom)%r20; jr3upu_21=radial_cache(jatom)%r21
    jr3upu_22=radial_cache(jatom)%r22; jr3upu_02=radial_cache(jatom)%r02
    jr3upu_12=radial_cache(jatom)%r12
    jlaplo1_02=radial_cache(jatom)%jl02; jlaplo1_12=radial_cache(jatom)%jl12
    TC1S00=radial_cache(jatom)%t100; TC1S10=radial_cache(jatom)%t110
    TC1S01=radial_cache(jatom)%t101; TC1S11=radial_cache(jatom)%t111
    TC2S00=radial_cache(jatom)%t200; TC2S10=radial_cache(jatom)%t210
    TC2S01=radial_cache(jatom)%t201; TC2S11=radial_cache(jatom)%t211
    TC3S00=radial_cache(jatom)%t300; TC3S10=radial_cache(jatom)%t310
    TC3S01=radial_cache(jatom)%t301; TC3S11=radial_cache(jatom)%t311
    juddu_02=radial_cache(jatom)%jd02; juddu_20=radial_cache(jatom)%jd20
    juddu_12=radial_cache(jatom)%jd12; juddu_21=radial_cache(jatom)%jd21
    juddu_22=radial_cache(jatom)%jd22
    judu_02=radial_cache(jatom)%j102; judu_20=radial_cache(jatom)%j120
    judu_12=radial_cache(jatom)%j112; judu_21=radial_cache(jatom)%j121
    judu_22=radial_cache(jatom)%j122
    juu_02=radial_cache(jatom)%j002; juu_20=radial_cache(jatom)%j020
    juu_12=radial_cache(jatom)%j012; juu_21=radial_cache(jatom)%j021
    juu_22=radial_cache(jatom)%j022
  end subroutine load_stress_integrals

  subroutine deallocate_stress_integrals
    if (allocated(Jllp00)) then
       deallocate(Jllp00,Jllp01,Jllp10,Jllp11,Illp00,Illp01,Illp10,Illp11)
       deallocate(jr3upu_20,jr3upu_21,jr3upu_22,jr3upu_02,jr3upu_12)
       deallocate(jlaplo1_02,jlaplo1_12)
       deallocate(TC1S00,TC1S10,TC1S01,TC1S11,TC2S00,TC2S10,TC2S01,TC2S11)
       deallocate(TC3S00,TC3S10,TC3S01,TC3S11)
       deallocate(juddu_02,juddu_20,juddu_12,juddu_21,juddu_22)
       deallocate(judu_02,judu_20,judu_12,judu_21,judu_22)
       deallocate(juu_02,juu_20,juu_12,juu_21,juu_22)
    endif
    if (allocated(radial_cache)) deallocate(radial_cache)
  end subroutine deallocate_stress_integrals
end module stress_kb

module stress_pw
  use param, only: lmax2
  implicit none
  real*8,allocatable :: DFJ_ST(:,:),DDFJ_ST(:,:),KplusG(:),F1P(:,:)
  real*8,allocatable :: DFJ1(:,:),DDFJ(:,:)
contains
  subroutine get_variable_pw(nw)
    integer,intent(in) :: nw
    if (allocated(KplusG)) deallocate(KplusG,F1P,DFJ_ST,DDFJ_ST,DFJ1,DDFJ)
    allocate(KplusG(nw),F1P(nw,3),DFJ_ST(0:lmax2,nw),DDFJ_ST(0:lmax2,nw))
    allocate(DFJ1(0:lmax2,nw),DDFJ(0:lmax2,nw))
  end subroutine get_variable_pw

  subroutine prepare_stress_pw(bk3,fj,dfj,rmt,isize)
    use structure, only: br1
    implicit none
    integer,intent(in) :: isize
    real*8,intent(in)  :: bk3(3,*),fj(0:lmax2,*),dfj(0:lmax2,*),rmt
    real*8 :: avec(3),ak,xa,dj(0:lmax2),ddj(0:lmax2)
    integer :: i,l
    do i=1,isize
       ! bk3 contains k+G in the reciprocal-lattice coordinates used by
       ! case.vector.  BR1 converts it to Cartesian inverse-Bohr units.
       avec=matmul(br1,bk3(:,i))
       ak=sqrt(dot_product(avec,avec)); xa=rmt*ak
       F1P(i,:)=avec; KplusG(i)=ak*ak
       if (ak.gt.1.d-14) then
          dj(:)=dfj(:,i)/ak
       else
          dj(:)=0.d0
       endif
       call ddvbes1_stress(fj(:,i),dj,ddj,xa,lmax2+1)
       do l=0,lmax2
          DFJ_ST(l,i)=dj(l)
          DDFJ_ST(l,i)=xa*ddj(l)
          DDFJ(l,i)=xa*ak*ddj(l)
          DFJ1(l,i)=rmt*dfj(l,i)
       enddo
    enddo
  end subroutine prepare_stress_pw

  subroutine ddvbes1_stress(fj,dj,ddj,x,nt)
    implicit none
    integer,intent(in) :: nt
    real*8,intent(in) :: fj(nt),dj(nt),x
    real*8,intent(out):: ddj(nt)
    integer :: i,l
    if (x.le.1.d-5) then
       ddj=0.d0; ddj(1)=-1.d0/3.d0
       if (nt.ge.3) ddj(3)=2.d0/15.d0
    else
       do i=1,nt
          l=i-1
          ddj(i)=((l+1)*l/(x*x)-1.d0)*fj(i)-2.d0*dj(i)/x
       enddo
    endif
  end subroutine ddvbes1_stress
end module stress_pw

module dmft_stress
  implicit none
  logical :: Qstress=.false.
  ! The isolated M_Sigma overlap contraction is useful when decomposing the
  ! stress, but the production FULL tensor already contains it through rhoE.
  ! Re-evaluating the complete LAPW overlap kernel only to print that number
  ! is expensive.  Keep it automatic for response-mode audits and opt-in for
  ! detailed diagnostics; leave it off in ordinary FULL calculations.
  logical :: stress_sigma_diag_enabled=.false.
  logical :: stress_detailed_diagnostics=.false.
  ! Accumulated wall-clock decomposition of the scalar-relativistic kernel.
  ! These counters are diagnostic only and are reduced with MPI_MAX together
  ! with the outer stress-pass timings in l2main.
  real*8 :: stress_rela_setup_wall=0.d0,stress_rela_band_wall=0.d0
  real*8 :: stress_val_corr(3,3),stress_kin_mt(3,3),stress_apw(3,3)
  real*8 :: stress_kin_int(3,3),stress_projector(3,3)
  ! Diagnostic selector for the direct projector contraction before mixer.
  ! CURRENT preserves the historical dmft2 assembly (+1).  WIEN inserts the
  ! physical dGamma/depsilon with the opposite raw sign (-1), because mixer
  ! later reports sigma=-Tstress.  The default remains CURRENT until the
  ! finite-difference closure tests settle all compensating terms.
  real*8 :: stress_projector_raw_sign=1.d0
  ! Diagnostic copies of the two scalar-relativistic additions.  They are
  ! already included in STRESS_KIN_MT; keeping them separate permits a
  ! direct comparison with LAPW2's :DBG_RELA_CORR and :DBG_RELA_SURF terms
  ! without risking a second addition to the physical stress.
  real*8 :: stress_rela_corr(3,3),stress_rela_surf(3,3)
  ! The three nonspherical lines of Blaha et al. Eq. (G2) are kept
  ! separately for term-by-term comparison with LAPW2.  A selected
  ! combination is also included once in stress_val_corr.  TRACE is the
  ! conservative production default used by the patched WIEN2k-24 code;
  ! FULL/AFFINE exposes the complete, still more delicate deviatoric tensor.
  real*8 :: stress_nsp_match(3,3),stress_nsp_radial(3,3)
  real*8 :: stress_nsp_angular(3,3)
  character(len=16) :: stress_g2_mode='TRACE'
  ! Part of the complementary-band response which enters through the
  ! k-global normalization X=(sum_k w_k P_k^+P_k)^(-1/2).  This tensor is
  ! also contained in the ordinary LAPW components after the second response
  ! pass; the separate copy is retained only for diagnostics and validation.
  real*8 :: stress_norm_window(3,3)
  ! Diagnostic decomposition only.  These arrays do not enter STR_VAL.
  ! The last index is (valence correction, muffin-tin kinetic,
  ! interstitial kinetic, APW surface).  Keeping the four kernels separate
  ! is important: a response pass must use its own signed pseudo-band
  ! weights in every kernel, and this makes an accidental reuse of the
  ! physical density immediately visible.
  real*8 :: stress_diag_main(3,3,4),stress_diag_window(3,3,4)
  real*8 :: stress_diag_norm(3,3,4),stress_diag_sigma_moment(3,3)
  complex*16,allocatable :: apw_surface_coeff(:,:)
  ! Explicit derivative of the (unnormalized on entry, normalized later in
  ! l2main) LAPW-to-DMFT projector.  The last index uses Voigt-like ordering
  ! (xx,yy,zz,yz,xz,xy), but stores tensor derivatives rather than engineering
  ! shear derivatives.
  complex*16,allocatable :: projector_strain(:,:,:,:)
contains
  subroutine init_dmft_stress(lmmax,nume)
    integer,intent(in) :: lmmax,nume
    integer :: env_status
    character(len=32) :: env_value
    stress_val_corr=0.d0; stress_kin_mt=0.d0; stress_apw=0.d0
    stress_kin_int=0.d0; stress_projector=0.d0
    stress_rela_corr=0.d0; stress_rela_surf=0.d0
    stress_nsp_match=0.d0; stress_nsp_radial=0.d0
    stress_nsp_angular=0.d0
    stress_norm_window=0.d0
    stress_diag_main=0.d0; stress_diag_window=0.d0
    stress_diag_norm=0.d0; stress_diag_sigma_moment=0.d0
    stress_sigma_diag_enabled=.false.
    stress_detailed_diagnostics=.false.
    stress_projector_raw_sign=1.d0
    stress_rela_setup_wall=0.d0; stress_rela_band_wall=0.d0
    if (.not.allocated(apw_surface_coeff)) allocate(apw_surface_coeff(lmmax,nume))
    apw_surface_coeff=0.d0

    ! Share the mode name with the patched LAPW2 implementation so DFT and
    ! eDMFT runs can be compared without translating input conventions.
    stress_g2_mode='TRACE'
    env_value=''
    call get_environment_variable('WIEN_STRESS_G2_MODE',env_value, &
         status=env_status)
    if (env_status.eq.0 .and. len_trim(env_value).gt.0) then
       select case(trim(adjustl(env_value)))
       case('NONE','none','0')
          stress_g2_mode='NONE'
       case('MATCH','match')
          stress_g2_mode='MATCH'
       case('RADIAL','radial')
          stress_g2_mode='RADIAL'
       case('ANGULAR','angular')
          stress_g2_mode='ANGULAR'
       case('TRACE','trace')
          stress_g2_mode='TRACE'
       case('AFFINE','affine','FULL','full')
          stress_g2_mode='FULL'
       case default
          write(6,*) 'Unknown WIEN_STRESS_G2_MODE: ',trim(env_value)
          stop 'G2 mode must be NONE, MATCH, RADIAL, ANGULAR, TRACE, or FULL'
       end select
    endif

    env_value='CURRENT'
    call get_environment_variable('EDMFT_STRESS_PROJECTOR_SIGN',env_value, &
         status=env_status)
    if (env_status.ne.0 .or. len_trim(env_value).eq.0) env_value='CURRENT'
    select case(trim(adjustl(env_value)))
    case('CURRENT','current','PLUS','plus','+1')
       stress_projector_raw_sign=1.d0
    case('WIEN','wien','PHYSICAL','physical','MINUS','minus','-1')
       stress_projector_raw_sign=-1.d0
    case default
       write(6,*) 'Unknown EDMFT_STRESS_PROJECTOR_SIGN: ',trim(env_value)
       stop 'projector sign must be CURRENT or WIEN'
    end select

    ! BASE and NO_SIGMA need the isolated M_Sigma tensor in order to remove
    ! it from the already assembled main stress.  Other response modes do not.
    env_value='FULL'
    call get_environment_variable('EDMFT_STRESS_RESPONSE_MODE',env_value, &
         status=env_status)
    if (env_status.ne.0 .or. len_trim(env_value).eq.0) env_value='FULL'
    select case(trim(adjustl(env_value)))
    case('BASE','base','NO_SIGMA','no_sigma')
       stress_sigma_diag_enabled=.true.
    end select

    ! Explicit opt-in for users who want the expensive printed decomposition
    ! while retaining the normal FULL production tensor.
    env_value=''
    call get_environment_variable('EDMFT_STRESS_DETAILED_DIAGNOSTICS', &
         env_value,status=env_status)
    if (env_status.eq.0 .and. len_trim(env_value).gt.0) then
       select case(trim(adjustl(env_value)))
       case('1','YES','yes','TRUE','true','ON','on')
          stress_sigma_diag_enabled=.true.
          stress_detailed_diagnostics=.true.
       end select
    endif

    ! The kernel-independent M_Sigma*dO finite-difference check needs the
    ! un-diagonalized M_Sigma matrix.  Reuse the same construction as the
    ! printed analytic decomposition, but keep the dump itself separately
    ! gated in WriteSigmaMomentOverlapDump.  An unset/empty/false value leaves
    ! the ordinary production path unchanged.
    env_value=''
    call get_environment_variable('EDMFT_MSIGMA_DUMP',env_value, &
         status=env_status)
    if (env_status.eq.0 .and. len_trim(env_value).gt.0) then
       select case(trim(adjustl(env_value)))
       case('0','OFF','off','FALSE','false','NO','no')
          continue
       case default
          stress_sigma_diag_enabled=.true.
       end select
    endif
  end subroutine init_dmft_stress
  subroutine fini_dmft_stress
    if (allocated(apw_surface_coeff)) deallocate(apw_surface_coeff)
    if (allocated(projector_strain)) deallocate(projector_strain)
  end subroutine fini_dmft_stress
end module dmft_stress
