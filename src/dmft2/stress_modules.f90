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

! Default-off reference for the radial part which is deliberately discarded
! by Blaha's frozen-augmentation approximation.  A center/minus/plus LAPW2
! triplet supplies finite-difference responses of u, udot, and each complete
! local-orbital basis function.  The arrays are consumed only by the
! diagnostic contraction in cmp_dmft_stress.f90; they never enter the normal
! production stress tensor.
module stress_radial_fd
  use param, only: nrad,lmax2,lomax,nloat,ngau
  implicit none
  logical :: radial_fd_enabled=.false.,radial_fd_prepared=.false., &
       radial_fd_projector_enabled=.false.
  real*8 :: radial_fd_projector_sign=1.d0
  character(len=512) :: radial_fd_path=''
  real*8,allocatable :: radial_fd_du(:,:,:,:),radial_fd_dudot(:,:,:,:), &
       radial_fd_dup(:,:,:,:),radial_fd_dudotp(:,:,:,:), &
       radial_fd_dlo(:,:,:,:,:)
  logical,allocatable :: radial_fd_seen(:,:,:),radial_fd_lo_seen(:,:,:,:)
  complex*16,allocatable :: radial_fd_overlap(:,:,:,:), &
       radial_fd_hsph(:,:,:,:),radial_fd_vnsp(:,:,:,:)
  real*8 :: radial_fd_rho_hsph=0.d0,radial_fd_rho_hnsp=0.d0, &
       radial_fd_rhoe_overlap=0.d0
  ! Pass-resolved copies of the three measured radial contractions.
  ! First index: 1=physical rho/rhoE, 2=fixed-X window response,
  ! 3=normalization response.  Second index: spherical H, nonspherical H,
  ! and overlap O.  These values are diagnostics only; they are never added
  ! to the production stress tensor.
  real*8 :: radial_fd_by_pass(3,3)=0.d0
contains
  subroutine init_stress_radial_fd(nat)
    integer,intent(in) :: nat
    integer :: env_status,projector_env_status,projector_sign_status, &
         unit,ios,atom,l,ir,jlo
    character(len=1024) :: line
    character(len=16) :: tag
    character(len=64) :: projector_environment,projector_sign_environment
    real*8 :: radius,values8(8),values4(4)

    radial_fd_enabled=.false.; radial_fd_prepared=.false.
    radial_fd_projector_enabled=.false.
    radial_fd_projector_sign=1.d0
    radial_fd_rho_hsph=0.d0; radial_fd_rho_hnsp=0.d0
    radial_fd_rhoe_overlap=0.d0
    radial_fd_by_pass=0.d0
    radial_fd_path=''
    call get_environment_variable('EDMFT_STRESS_RADIAL_RESPONSE', &
         radial_fd_path,status=env_status)
    if (env_status.ne.0 .or. len_trim(radial_fd_path).eq.0) return
    projector_environment=''
    call get_environment_variable('EDMFT_STRESS_RADIAL_PROJECTOR', &
         projector_environment,status=projector_env_status)
    radial_fd_projector_enabled=projector_env_status.eq.0 .and. &
         len_trim(projector_environment).gt.0 .and. &
         trim(projector_environment).ne.'0'
    projector_sign_environment='CURRENT'
    call get_environment_variable('EDMFT_STRESS_RADIAL_PROJECTOR_SIGN', &
         projector_sign_environment,status=projector_sign_status)
    if (projector_sign_status.ne.0 .or. &
        len_trim(projector_sign_environment).eq.0) &
         projector_sign_environment='CURRENT'
    select case(trim(adjustl(projector_sign_environment)))
    case('CURRENT','current','PLUS','plus','+1')
       radial_fd_projector_sign=1.d0
    case('PHYSICAL','physical','MINUS','minus','-1')
       radial_fd_projector_sign=-1.d0
    case default
       write(6,*) 'Unknown EDMFT_STRESS_RADIAL_PROJECTOR_SIGN: ', &
            trim(projector_sign_environment)
       stop 'radial projector sign must be CURRENT or PHYSICAL'
    end select

    allocate(radial_fd_du(2,nrad,0:lmax2,nat), &
         radial_fd_dudot(2,nrad,0:lmax2,nat), &
         radial_fd_dup(2,nrad,0:lmax2,nat), &
         radial_fd_dudotp(2,nrad,0:lmax2,nat), &
         radial_fd_dlo(2,nrad,nloat,0:lomax,nat), &
         radial_fd_seen(nrad,0:lmax2,nat), &
         radial_fd_lo_seen(nrad,nloat,0:lomax,nat))
    radial_fd_du=0.d0; radial_fd_dudot=0.d0
    radial_fd_dup=0.d0; radial_fd_dudotp=0.d0
    radial_fd_dlo=0.d0; radial_fd_seen=.false.
    radial_fd_lo_seen=.false.

    open(newunit=unit,file=trim(radial_fd_path),status='old', &
         action='read',iostat=ios)
    if (ios.ne.0) then
       write(6,*) 'Cannot open EDMFT_STRESS_RADIAL_RESPONSE=', &
            trim(radial_fd_path)
       stop 'missing finite-difference radial-response diagnostic'
    endif
    do
       read(unit,'(a)',iostat=ios) line
       if (ios.ne.0) exit
       if (len_trim(line).eq.0 .or. line(1:1).eq.'#') cycle
       tag=''
       read(line,*,iostat=ios) tag
       if (ios.ne.0) cycle
       select case(trim(tag))
       case(':DR')
          read(line,*,iostat=ios) tag,atom,l,ir,radius,values8
          if (ios.ne.0 .or. atom.lt.1 .or. atom.gt.nat .or. &
              l.lt.0 .or. ir.lt.1 .or. ir.gt.nrad) &
               stop 'invalid :DR radial-response record'
          ! WIEN2k may tabulate radial functions beyond dmft2's LMAX2.
          ! Those channels are not present in ALM/BLM and must simply be
          ! ignored; rejecting them made a perfectly valid dump unreadable.
          if (l.gt.lmax2) cycle
          radial_fd_du(:,ir,l,atom)=values8(1:2)
          radial_fd_dudot(:,ir,l,atom)=values8(3:4)
          radial_fd_dup(:,ir,l,atom)=values8(5:6)
          radial_fd_dudotp(:,ir,l,atom)=values8(7:8)
          radial_fd_seen(ir,l,atom)=.true.
       case(':DL')
          read(line,*,iostat=ios) tag,atom,l,jlo,ir,radius,values4
          if (ios.ne.0 .or. atom.lt.1 .or. atom.gt.nat .or. &
              l.lt.0 .or. jlo.lt.1 .or. ir.lt.1 .or. ir.gt.nrad) &
               stop 'invalid :DL radial-response record'
          ! As for :DR, a producer is allowed to include channels and LO
          ! indices not represented by this dmft2 build.
          if (l.gt.lomax .or. jlo.gt.nloat) cycle
          radial_fd_dlo(:,ir,jlo,l,atom)=values4(1:2)
          radial_fd_lo_seen(ir,jlo,l,atom)=.true.
       case default
          write(6,*) 'Unknown radial-response record: ',trim(tag)
          stop 'unsupported EDMFT radial-response file version'
       end select
    enddo
    close(unit)
    radial_fd_enabled=.true.
    write(6,'(A,1x,A)') ':DMFT_STRESS_RADIAL_FD_INPUT',trim(radial_fd_path)
    write(6,'(A,1x,L1)') ':DMFT_STRESS_RADIAL_PROJECTOR_DIAGNOSTIC', &
         radial_fd_projector_enabled
    write(6,'(A,1x,F6.1)') ':DMFT_STRESS_RADIAL_PROJECTOR_SIGN', &
         radial_fd_projector_sign
    write(6,'(A)') ':DMFT_STRESS_RADIAL_FD_DIAGNOSTIC_ONLY T'
  end subroutine init_stress_radial_fd

  subroutine allocate_stress_radial_fd_matrices(nat)
    integer,intent(in) :: nat
    if (.not.radial_fd_enabled) return
    allocate(radial_fd_overlap(2+nloat,2+nloat,0:lmax2,nat), &
         radial_fd_hsph(2+nloat,2+nloat,0:lmax2,nat), &
         radial_fd_vnsp(2+nloat,2+nloat,ngau,nat))
    radial_fd_overlap=0.d0; radial_fd_hsph=0.d0
    radial_fd_vnsp=0.d0
  end subroutine allocate_stress_radial_fd_matrices

  subroutine fini_stress_radial_fd
    if (allocated(radial_fd_du)) deallocate(radial_fd_du,radial_fd_dudot, &
         radial_fd_dup,radial_fd_dudotp,radial_fd_dlo,radial_fd_seen, &
         radial_fd_lo_seen)
    if (allocated(radial_fd_overlap)) deallocate(radial_fd_overlap, &
         radial_fd_hsph,radial_fd_vnsp)
    radial_fd_enabled=.false.; radial_fd_prepared=.false.
    radial_fd_projector_enabled=.false.
  end subroutine fini_stress_radial_fd
end module stress_radial_fd

module dmft_stress
  implicit none
  logical :: Qstress=.false.
  ! The isolated M_Sigma overlap contraction is useful when decomposing the
  ! stress, but the production FULL tensor already contains it through rhoE.
  ! Re-evaluating the complete LAPW overlap kernel only to print that number
  ! is expensive.  Keep it automatic for response-mode audits and opt-in for
  ! detailed diagnostics; leave it off in ordinary FULL calculations.
  logical :: stress_sigma_diag_enabled=.false.
  ! Default-off selector for the off-shell overlap derivative carried by
  ! M_Sigma.  CURRENT retains K(0)-K(1) extracted from Blaha's on-shell,
  ! volume-compensated kernel.  FIXED_SPHERE replaces only that M_Sigma
  ! channel by the literal derivative of the physical ALM/BLM/CLM metric plus
  ! the exact interstitial step-function derivative.  The analytic
  ! coefficient derivative has an independent in-code central-difference
  ! oracle and does not alter the ordinary DFT/on-shell CAPD path.
  logical :: stress_fixed_sphere_overlap_enabled=.false.
  logical :: stress_detailed_diagnostics=.false.
  ! Default-off fixed-coefficient oracle for the isolated M_Sigma overlap
  ! channel.  It central-differences the actual ALM/BLM/CLM construction at
  ! fixed LAPW eigenvectors and radial functions, and therefore tests the
  ! complete muffin-tin metric derivative without using Blaha's rearranged
  ! on-shell kernel.  The result is diagnostic only until analytic closure
  ! has been demonstrated on all reference materials.
  logical :: stress_coeff_fd_enabled=.false.
  logical :: stress_coeff_analytic_enabled=.false.
  ! Default-off paired fixed-band diagnostic.  It contracts rho and rhoE in
  ! one common band basis with the full physical coefficient derivative, as
  ! Force1_DMFT does for atomic displacements.  This preserves the H/O
  ! cancellation before summation and is used to audit the complementary-
  ! window response without changing the production tensor.
  logical :: stress_paired_diag_enabled=.false.
  ! Default-off audit of the generalized Blaha contraction.  When enabled,
  ! l2main separates the interacting density matrix into its diagonal and
  ! off-diagonal KS-band parts and sends each part, together with
  ! rhoE=1/2{epsilon,rho}, through exactly the same analytic LAPW kernels.
  ! The resulting tensors are diagnostic copies only: every physical stress
  ! accumulator is restored after each audit call.
  logical :: stress_matrix_decomposition_enabled=.false.
  ! Default-off covariant assembly audit.  In this representation the
  ! ordinary LAPW stress sees only the diagonal interacting band
  ! occupations.  The complete non-diagonal energy moment is carried by the
  ! derivative of the normalized projector, including its within-window
  ! band connection and generalized-metric connection.  It is algebraically
  ! equivalent to the fixed-coefficient H/O representation when every term
  ! is exact, and is substantially better conditioned for a truncated DMFT
  ! band window.
  logical :: stress_moving_band_enabled=.false.
  ! Opt-in end-to-end assembly of the validated paired matrix derivative.
  ! This is deliberately separate from PAIRED_DIAG: the latter only prints
  ! reference tensors, whereas PAIRED_ASSEMBLY replaces the historical
  ! off-shell pseudo-band contractions in STR_VAL.  It remains false by
  ! default until the total frozen-Sigma finite differences close on all
  ! reference materials.
  logical :: stress_paired_assembly_enabled=.false.
  ! Accumulated wall-clock decomposition of the scalar-relativistic kernel.
  ! These counters are diagnostic only and are reduced with MPI_MAX together
  ! with the outer stress-pass timings in l2main.
  real*8 :: stress_rela_setup_wall=0.d0,stress_rela_band_wall=0.d0
  real*8 :: stress_val_corr(3,3),stress_kin_mt(3,3),stress_apw(3,3)
  real*8 :: stress_kin_int(3,3),stress_projector(3,3)
  ! Direct +dGamma/depsilon must enter raw STR_VAL with a minus sign,
  ! since mixer reports -STR_VAL. CURRENT/+1 remains an explicit legacy test.
  real*8 :: stress_projector_raw_sign=-1.d0
  ! H'/O' response matrices already pass through CMP_DMFT_STRESS's WIEN sign
  ! conversion. PHYSICAL therefore means +1; LEGACY_REVERSED is diagnostic.
  real*8 :: stress_norm_response_sign=1.d0
  ! Diagnostic copies of the two scalar-relativistic additions.  They are
  ! already included in STRESS_KIN_MT; keeping them separate permits a
  ! direct comparison with LAPW2's :DBG_RELA_CORR and :DBG_RELA_SURF terms
  ! without risking a second addition to the physical stress.
  real*8 :: stress_rela_corr(3,3),stress_rela_surf(3,3)
  ! The three nonspherical lines of Blaha et al. Eq. (G2) are kept
  ! separately for term-by-term comparison with LAPW2.  A selected
  ! combination is included once in stress_val_corr. FULL uses the complete
  ! local matrices. TRACE remains available only as an explicit approximation
  ! for comparing historical hydrostatic calculations.
  real*8 :: stress_nsp_match(3,3),stress_nsp_radial(3,3)
  real*8 :: stress_nsp_angular(3,3)
  character(len=16) :: stress_g2_mode='FULL'
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
  ! The overlap part +2 Tr[rhoE O_,epsilon] of each CMP_DMFT_STRESS
  ! pass.  Index 1 is the physical/main pass, 2 the fixed-window response,
  ! and 3 the global-normalization response.  This diagnostic lets the
  ! finite-difference coefficient-space oracle distinguish an error in
  ! O_,epsilon from one in H_,epsilon without changing the assembled stress.
  real*8 :: stress_diag_response_overlap(3,3,3)
  ! Spherical pieces of the isolated M_Sigma overlap diagnostic.  These sum
  ! exactly to stress_diag_sigma_moment and never enter the physical tensor.
  real*8 :: stress_diag_sigma_match(3,3),stress_diag_sigma_radial(3,3)
  real*8 :: stress_diag_sigma_angular(3,3)
  real*8 :: stress_diag_sigma_interstitial(3,3)
  real*8 :: stress_diag_sigma_coeff_fd_mt(3,3)
  real*8 :: stress_diag_sigma_coeff_analytic_mt(3,3)
  real*8 :: stress_diag_paired_window_matching(3,3)
  real*8 :: stress_diag_paired_norm_matching(3,3)
  real*8 :: stress_diag_paired_main_correction_matching(3,3)
  ! Complete fixed-sphere matrix derivative associated with the three
  ! paired diagnostics above.  MATCHING contains only the derivative of the
  ! augmented radial coefficients.  COMPLETE additionally contains the
  ! muffin-tin kinetic surface matrix, the nonspherical potential matrix,
  ! and the interstitial kinetic/overlap derivative.  These remain
  ! diagnostic-only until they have closed against the independent LAPW1
  ! coefficient-space finite difference on every reference material.
  real*8 :: stress_diag_paired_window_complete(3,3)
  real*8 :: stress_diag_paired_norm_complete(3,3)
  real*8 :: stress_diag_paired_main_correction_complete(3,3)
  ! Potential-free counterpart of COMPLETE.  It is the derivative of the
  ! kinetic-overlap combination
  !
  !   -Tr[rho T_,epsilon] + Tr[rhoE O_,epsilon],
  !
  ! and is therefore the object which may replace WIEN's complete STR_VAL
  ! channel.  The full-H COMPLETE tensors above are indispensable operator
  ! diagnostics, but inserting them into STR_VAL would double count the
  ! potential derivative which is already cancelled by STR_DENS.
  real*8 :: stress_diag_paired_window_kinetic_matching(3,3)
  real*8 :: stress_diag_paired_norm_kinetic_matching(3,3)
  real*8 :: stress_diag_paired_main_correction_kinetic_matching(3,3)
  real*8 :: stress_diag_paired_window_kinetic_complete(3,3)
  real*8 :: stress_diag_paired_norm_kinetic_complete(3,3)
  real*8 :: stress_diag_paired_main_correction_kinetic_complete(3,3)
  ! Exact full-matrix contraction of the many-body energy moment
  !
  !   M_Sigma = rhoE - 1/2 {epsilon,rho}
  !
  ! with the LAPW metric derivative.  This is evaluated before any
  ! pseudo-band diagonalization.  That distinction is essential when the
  ! signed moment is ill-conditioned or truncated: FeSe loses most of this
  ! contribution in the older eigenvector reconstruction, although the
  ! direct band-space matrix is perfectly regular.  With rho=0 in the
  ! paired H/O kernel, MATCHING is the muffin-tin coefficient response and
  ! COMPLETE additionally contains the interstitial overlap response.
  real*8 :: stress_diag_paired_sigma_matching(3,3)
  real*8 :: stress_diag_paired_sigma_complete(3,3)
  ! Ordinary on-shell LAPW stress of the fully occupied lower bands plus the
  ! diagonal part of the interacting density in the DMFT window.  Together
  ! with PAIRED_MAIN_CORRECTION_COMPLETE this is a non-overlapping
  ! reconstruction of the complete one-particle valence channel.
  real*8 :: stress_diag_paired_diagonal_reference(3,3,4)
  real*8 :: stress_diag_density_diagonal(3,3,4)
  real*8 :: stress_diag_density_offdiagonal(3,3,4)
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
    stress_diag_response_overlap=0.d0
    stress_diag_sigma_match=0.d0; stress_diag_sigma_radial=0.d0
    stress_diag_sigma_angular=0.d0
    stress_diag_sigma_interstitial=0.d0
    stress_diag_sigma_coeff_fd_mt=0.d0
    stress_diag_sigma_coeff_analytic_mt=0.d0
    stress_diag_paired_window_matching=0.d0
    stress_diag_paired_norm_matching=0.d0
    stress_diag_paired_main_correction_matching=0.d0
    stress_diag_paired_window_complete=0.d0
    stress_diag_paired_norm_complete=0.d0
    stress_diag_paired_main_correction_complete=0.d0
    stress_diag_paired_window_kinetic_matching=0.d0
    stress_diag_paired_norm_kinetic_matching=0.d0
    stress_diag_paired_main_correction_kinetic_matching=0.d0
    stress_diag_paired_window_kinetic_complete=0.d0
    stress_diag_paired_norm_kinetic_complete=0.d0
    stress_diag_paired_main_correction_kinetic_complete=0.d0
    stress_diag_paired_sigma_matching=0.d0
    stress_diag_paired_sigma_complete=0.d0
    stress_diag_paired_diagonal_reference=0.d0
    stress_diag_density_diagonal=0.d0
    stress_diag_density_offdiagonal=0.d0
    stress_sigma_diag_enabled=.false.
    stress_fixed_sphere_overlap_enabled=.false.
    stress_detailed_diagnostics=.false.
    stress_coeff_fd_enabled=.false.
    stress_coeff_analytic_enabled=.false.
    stress_paired_diag_enabled=.false.
    stress_matrix_decomposition_enabled=.false.
    stress_moving_band_enabled=.false.
    stress_paired_assembly_enabled=.false.
    stress_projector_raw_sign=-1.d0
    stress_norm_response_sign=1.d0
    stress_rela_setup_wall=0.d0; stress_rela_band_wall=0.d0
    if (.not.allocated(apw_surface_coeff)) allocate(apw_surface_coeff(lmmax,nume))
    apw_surface_coeff=0.d0

    ! Share the mode name with the patched LAPW2 implementation so DFT and
    ! eDMFT runs can be compared without translating input conventions.
    stress_g2_mode='FULL'
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

    env_value=''
    call get_environment_variable('EDMFT_STRESS_PAIRED_DIAG',env_value, &
         status=env_status)
    if (env_status.eq.0 .and. len_trim(env_value).gt.0) then
       select case(trim(adjustl(env_value)))
       case('1','YES','yes','TRUE','true','ON','on')
          stress_paired_diag_enabled=.true.
       end select
    endif

    env_value=''
    call get_environment_variable('EDMFT_STRESS_PAIRED_ASSEMBLY', &
         env_value,status=env_status)
    if (env_status.eq.0 .and. len_trim(env_value).gt.0) then
       select case(trim(adjustl(env_value)))
       case('1','YES','yes','TRUE','true','ON','on')
          stress_paired_assembly_enabled=.true.
          stress_paired_diag_enabled=.true.
          stress_sigma_diag_enabled=.true.
          stress_coeff_analytic_enabled=.true.
          ! The exact paired main correction replaces one inseparable block:
          ! the DFT-like off-diagonal density plus M_Sigma.  Request the
          ! corresponding legacy off-diagonal tensor so that block can be
          ! removed without touching the on-shell diagonal WIEN stress.
          stress_matrix_decomposition_enabled=.true.
       end select
    endif

    ! This audit is intentionally expensive: for each k point and atom it
    ! rebuilds the LAPW matching coefficients at +/-h for all six strain
    ! components.  It is never activated implicitly by ordinary detailed
    ! diagnostics.
    env_value=''
    call get_environment_variable('EDMFT_STRESS_COEFF_FD',env_value, &
         status=env_status)
    if (env_status.eq.0 .and. len_trim(env_value).gt.0) then
       select case(trim(adjustl(env_value)))
       case('1','YES','yes','TRUE','true','ON','on')
          stress_coeff_fd_enabled=.true.
          stress_coeff_analytic_enabled=.true.
          stress_sigma_diag_enabled=.true.
       end select
    endif

    env_value=''
    call get_environment_variable('EDMFT_STRESS_MOVING_BAND',env_value, &
         status=env_status)
    if (env_status.eq.0 .and. len_trim(env_value).gt.0) then
       select case(trim(adjustl(env_value)))
       case('1','YES','yes','TRUE','true','ON','on')
          stress_moving_band_enabled=.true.
          ! The assembly needs the otherwise diagnostic separation of the
          ! DFT-like off-diagonal density and M_Sigma overlap channels.
          stress_matrix_decomposition_enabled=.true.
          stress_sigma_diag_enabled=.true.
       end select
    endif

    env_value='CURRENT'
    call get_environment_variable('EDMFT_STRESS_OVERLAP_MODE',env_value, &
         status=env_status)
    if (env_status.ne.0 .or. len_trim(env_value).eq.0) env_value='CURRENT'
    select case(trim(adjustl(env_value)))
    case('CURRENT','current')
       stress_fixed_sphere_overlap_enabled=.false.
    case('FIXED_SPHERE','fixed_sphere','FIXED-SPHERE','fixed-sphere')
       stress_fixed_sphere_overlap_enabled=.true.
       stress_coeff_analytic_enabled=.true.
       ! The replacement tensors are produced by the otherwise optional
       ! M_Sigma diagnostic pass, so enable that pass automatically.
       stress_sigma_diag_enabled=.true.
    case default
       write(6,*) 'Unknown EDMFT_STRESS_OVERLAP_MODE: ',trim(env_value)
       stop 'overlap mode must be CURRENT or FIXED_SPHERE'
    end select

    ! This decomposition is intentionally expensive (two additional analytic
    ! kernel contractions per k point), so it must never be activated merely
    ! by requesting ordinary detailed output.
    env_value=''
    call get_environment_variable('EDMFT_STRESS_MATRIX_DECOMPOSITION', &
         env_value,status=env_status)
    if (env_status.eq.0 .and. len_trim(env_value).gt.0) then
       select case(trim(adjustl(env_value)))
       case('1','YES','yes','TRUE','true','ON','on')
          stress_matrix_decomposition_enabled=.true.
       end select
    endif

    env_value=''
    call get_environment_variable('EDMFT_STRESS_PROJECTOR_SIGN',env_value, &
         status=env_status)
    ! The literal projector derivative enters STR_VAL with a minus sign;
    ! mixer supplies the final minus sign for every assembly, not just PAIRED.
    if (env_status.ne.0 .or. len_trim(env_value).eq.0) env_value='WIEN'
    select case(trim(adjustl(env_value)))
    case('CURRENT','current','PLUS','plus','+1')
       stress_projector_raw_sign=1.d0
    case('WIEN','wien','PHYSICAL','physical','MINUS','minus','-1')
       stress_projector_raw_sign=-1.d0
    case default
       write(6,*) 'Unknown EDMFT_STRESS_PROJECTOR_SIGN: ',trim(env_value)
       stop 'projector sign must be CURRENT or WIEN'
    end select

    env_value='CURRENT'
    call get_environment_variable('EDMFT_STRESS_NORM_RESPONSE_SIGN', &
         env_value,status=env_status)
    if (env_status.ne.0 .or. len_trim(env_value).eq.0) env_value='CURRENT'
    select case(trim(adjustl(env_value)))
    case('CURRENT','current','PLUS','plus','+1','PHYSICAL','physical','WIEN','wien')
       stress_norm_response_sign=1.d0
    case('LEGACY_REVERSED','legacy_reversed','MINUS','minus','-1')
       stress_norm_response_sign=-1.d0
    case default
       write(6,*) 'Unknown EDMFT_STRESS_NORM_RESPONSE_SIGN: ', &
            trim(env_value)
       stop 'normalization response sign must be PHYSICAL or LEGACY_REVERSED'
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
