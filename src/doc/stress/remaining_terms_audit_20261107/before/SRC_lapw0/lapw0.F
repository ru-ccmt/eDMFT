     PROGRAM LAPW0
!
!              P O T E N T I A L   G E N E R A T O R
!              ACTUAL VERSION
!
!     LAPW0 CALCULATES THE TOTAL POTENTIAL BY SOLVING THE POISSON-
!     EQUATION  AND USING THE LSD, GGA OR MGGA APPROXIMATIONS
!     FULL NON-MUFFIN-TIN TOTAL POTENTIAL COMPLEX VERSION
!
      USE parallel
      USE defs
      USE ams
      USE struct
      USE efeld
      USE symetr
      USE densit
      USE work
      USE vresp1
      USE xcparam
      use norm_kub,   only: c_kub, init_c_kub
      USE PW_mult,    only: DECOUL_POT_PW,CFFT_STR,SGGP,CFFT8,SGGP_GGA, &          !CKAMAL
                            CFFT_STR_GGA,META_FFT
      use gga_str_kb,  only: gga_int !, sigma_str
      use str_switch,  only: SW_STRESS, SW_TOT_STR
      use stress_time
      use, intrinsic :: ISO_C_BINDING, only: C_INT
      
      USE fftw_parallel_base
      USE TestCases,only        : CoulEM
#ifdef LIBXC
      use libxc_mod, only: libxc_init, libxc_end
#endif
#if defined Parallel
      USE fftw_parallel, only: prepare_parallel_ffts, cleanup_ffts
#else
      USE fftw_serial, only: cleanup_ffts
#endif
!      USE fftw_para, only : iff3_g, local_last_start
      implicit none
      INCLUDE 'param.inc'

#ifdef Parallel
      include 'mpif.h'
#endif

!     
! THE PARAMETERS MUST SATISFY THE FOLLOWING EQUATIONS :
! NCOM :     NCOM+3 >= SUM(L,M=0) + 2*SUM(L,M>0)
! DIMENSION OF YKA, QQ,..
!     
      LOGICAL          LR2V, force,force_buf,lrhopw,xcpot1qq,isopen, &
           vsp_only, core_add, core_sub, stress_int_delta_diag
      LOGICAL,ALLOCATABLE :: logefg(:)
      CHARACTER*19     nkktext
      CHARACTER*11     STATUS,FORM
      CHARACTER*4      SWITCH,R2V
      CHARACTER*5      switc1,switch1,switch2,switch3,idstr
      CHARACTER*67     ERRMSG
      CHARACTER*77     RMARGI
      CHARACTER*180     fname,fnavn(99)
      CHARACTER*80     DEFFN, ERRFN, parseline
      CHARACTER*32     stress_int_delta_env
      CHARACTER*60     EX_SWITCH_TMP, EC_SWITCH_TMP, VX_SWITCH_TMP, VC_SWITCH_TMP, VCVDISP_SWITCH_TMP
      integer            :: I,J,k,JJ,J1,J2,iunit,iunitr2v,irecl,il,lmma,nuc, &
                            iunit_vtau, iunit_vlapl
      integer            :: JSPIN,ISCF,JATOM,IND,LMMULT1
      integer            :: INDEX,JROW,JCOL,MU,JX,LM1,L,M,NKK,NKK1,NKK2,NKK3,nkknonloc,nkkr2v,LLMM
      integer            :: LL,LMX,NCON,NCON1,IBES,JRJ,LMOD,IRAD,IK,ISPIN
      integer            :: JRE,ITAPS,ITAPN,INFO
      integer            :: IFF1,IFF2,IFF3,II1,II2,II3,IFF1T,IFF2T,IFF3T
      integer            :: iffx2,iffy2,iffz2
      integer            :: jrjstart,jrjend,jrj1
      integer            :: jff1,jff2,jff3
      integer            :: l_max,l_min
      integer            :: ib_max
      integer            :: ios,iprint,luse2,inod,numbeece,nxcprint
      integer            :: stress_int_delta_status
      integer            :: i1,i2,i3
      integer            :: it1, it2, ibbuff, myid1
      REAL*8             :: MP,av(3),dfak(42),pi,sqfp,sqrt2,pi4
      REAL*8             :: SUPN,DELCC,QIN,ARG,RK,RK1,d1,rnuc
      REAL*8             :: T947,TE947,TS947,TIMTOT,TMULTI,PMULTI,TCOULI
      REAL*8             :: PCOULI,TXCSPH,PXCSPH,TXCPW,PXCPW,TENERG
      REAL*8             :: PENERG,TREANA,PREANA,TIMEY,TIME0,TIME1,TIME2
      REAL*8             :: TIME3,TIME4,TIME5,TIME6,TCOUL,PCOUL,time4a,time4b
      REAL*8             :: RJ1,RJ2,ACORR,VOLINT,VLIN1,VLIN2,SURFIN
      REAL*8             :: ZSHIFT,VINF,VOUTF,EKVOL1,UST0
      REAL*8             :: ETOTAL,ekinint,excint,ekin,exc2
      real*8             :: etot_bg,espxc,vspxc,fexck,fvxck
      real*8             :: pxcpww,penergw
      REAL*8             :: VC,vxcrhoc,vtautauc,vxcrhoc_buf,vtautauc_buf
      real*8             :: tmp,tmp1,potktest,faciff
      real*8             :: time0w,time1w,time2w,time3w,time4aw,time4bw,time4w,time5w,&
                            time6w,timeyw,tcoulw,txcsphw,txcpww,pxcsphw,treanaw,&
                            preanaw,tenergw,pcoulw,pcouliw,tcouliw,pmultiw,tmultiw,&
                            timtotw
      real*8             :: t947w,ts947w,te947w
      REAL*8             ::  vfield(ncom+3,nrad)
      real*8,allocatable :: clm_tmp(:)
      real*8             :: bextene, vxmix1,dis
      integer            :: imfield,nmod,nsp,natorb
      REAL*8,ALLOCATABLE :: sigma(:,:),bes(:),bespw(:,:)
      COMPLEX*16       IMAG,IMAG1,QOUT
      COMPLEX*16       RHOKIN,RHOKT
      COMPLEX*16       CP,CQ(ncom+3),CVOUT(ncom+3)
!      COMPLEX*16 CrhoOUT(ncom+3,5)
      COMPLEX*16,ALLOCATABLE   ::  rhopw00(:,:)
      COMPLEX*16,ALLOCATABLE   ::  IMAGA(:,:)
      COMPLEX*16,ALLOCATABLE   ::  cfft(:,:,:),cfft2(:,:,:),cfft3(:,:,:)
      COMPLEX*16,ALLOCATABLE   ::  fft(:,:,:,:)
!      COMPLEX*16,ALLOCATABLE   ::  meta_fft(:,:,:,:)
      COMPLEX*16,ALLOCATABLE   ::  ust(:,:,:,:)
      
      COMPLEX*16,ALLOCATABLE   ::  potk(:),buf1(:),vkcoul(:),vkexc(:,:), &
           vkvxc(:,:),cvalue(:),taup(:),vkvxcold(:,:),vkvxcnonloc(:,:)
      REAL*8,ALLOCATABLE       ::  absk(:)
      REAL*8,ALLOCATABLE       ::  fhf(:,:)
      REAL*8,ALLOCATABLE       ::  vtot(:,:,:)
      REAL*8,ALLOCATABLE       ::  rhoc(:,:,:), tauc(:,:,:)

      COMPLEX*16,ALLOCATABLE :: fc(:,:)
      INTEGER,ALLOCATABLE    :: lmmtmx(:),lmmult(:,:,:),lmmax1(:)
      INTEGER,ALLOCATABLE    :: lmmax(:),lm(:,:,:)
      REAL*8,ALLOCATABLE     :: qel(:),vcoul0(:),rldm1(:),rldm2(:),wldm(:) ! *** LDM ADDITION
      COMPLEX*16,ALLOCATABLE :: qq(:,:),qqnew(:,:)
      
!=============defining the variables for the stress tensor calculation==============
      real*8,allocatable       :: sigma_str(:,:,:,:)
      REAL*8                   :: fac1_p,E_xc_m_V_xc,E_xc_m_V_xc_rmt(1:9)      
      REAL*8                   :: electrostatic_int_jacobian
      INTEGER                  :: index_p,alpha1,beta1, indexp
      COMPLEX*16               :: coulomb_str_int(1:9),gga_str_int(1:9)      
      
!!$      TYPE efg_type
!!$        REAL*8 :: v20,v21,v22
!!$        REAL*8 :: v21m,v22m
!!$        REAL*8 :: v20srf,v21srf,v22srf
!!$        REAL*8 :: v21msr,v22msr,vzz
!!$        REAL*8 :: winkel,eta
!!$        REAL*8 :: qmat(3,3)
!!$        REAL*8 :: tens(3)
!!$        REAL*8 :: eivec(3,3)
!!$      END TYPE efg_type      
!!$      TYPE (efg_type),ALLOCATABLE,TARGET :: efgb(:)


        REAL*8, allocatable :: efg_v20(:),efg_v21(:),efg_v22(:)
        REAL*8, allocatable :: efg_v21m(:),efg_v22m(:)
        REAL*8, allocatable :: efg_v20srf(:),efg_v21srf(:),efg_v22srf(:)
        REAL*8, allocatable :: efg_v21msr(:),efg_v22msr(:),efg_vzz(:)
        REAL*8, allocatable :: efg_winkel(:),efg_eta(:)
        REAL*8, allocatable :: efg_qmat(:,:,:)
        REAL*8, allocatable :: efg_tens(:,:)
        REAL*8, allocatable :: efg_eivec(:,:,:)

        COMPLEX*16, allocatable  ::     vk_vtau(:,:),vk_vlapl(:,:), &
             vk_vtau_nonrean(:,:)
        
#ifdef Parallel
!!$      INTEGER               :: efgb_blocklength(16),efgb_displace(16)
!!$      INTEGER               :: efgb_address(17),efgb_typelist(16)
!!$      INTEGER               :: typempi_efgb
#endif

#ifdef FFTW_OMP
      integer(C_INT) omp_error
      integer omp_get_max_threads
      interface
         subroutine fftw_plan_with_nthreads(nthreads) bind(C, name='fftw_plan_with_nthreads')
            import :: C_INT
            integer(C_INT), value :: nthreads
         end subroutine fftw_plan_with_nthreads
         integer(C_INT) function fftw_init_threads() bind(C, name='fftw_init_threads')
            import :: C_INT
         end function fftw_init_threads
      end interface
#endif

! additional variables for EECE
      real*8 vcoef
      real*8 ,allocatable:: elnum(:),elnumsp(:,:),fk(:,:)
      integer ifil
      integer, allocatable:: jeece(:)
!
! **** LDM ADDITIONAL WORKSPACE ***
!     For Kahan Sums
      COMPLEX*16 CWK(ncom+3), TC
!
!     For extrapolation
      complex*16 Extrapol(26,ncom+3)
! **** END LDM CHANGE
      
      DATA         IMAG    /(0.D0,1.D0)/
!-----------------------------------------------------------------------
!     Global initializer of stack setting (and signal-handlers for mpi modes)
      CALL W2kinit
!     initialise atomic masses      
      CALL init_ams
      call init_c_kub
      ALLOCATE(yka(ncom+3),r(nrad),value(nrad))
      ALLOCATE(RLDM1(nrad),rldm2(nrad),bes(0:max(lmax2+1,20)),bespw(0:max(lmax2+1,20),1:nrad+10))
      ALLOCATE(WLDM(NRAD) ) ! *** LDM ADDITION
      force_buf = .FALSE.

!     Equation (G5) contains two interstitial electrostatic terms:
!
!       1/2 int_I rho [ -delta_ab V_C + dV_C/d epsilon_ab ].
!
!     The historical stress path constructs the second term in
!     CFFT_STR.  The first, purely diagonal Jacobian term is tested below
!     behind an environment switch until finite-difference validation is
!     complete.  It costs no extra FFT: EKVOL1 is the interstitial integral
!     of rho*(V_xc-e_xc+V_C/2), while E_xc_m_V_xc is the integral of
!     rho*(e_xc-V_xc); their sum is exactly int_I rho*V_C/2 in the same
!     normalization used by REAN3.
      stress_int_delta_env = ''
      call get_environment_variable('WIEN_STRESS_INTERSTITIAL_DELTA', &
                                    stress_int_delta_env, &
                                    status=stress_int_delta_status)
      stress_int_delta_diag = stress_int_delta_status.eq.0 .and. &
           len_trim(stress_int_delta_env).gt.0 .and. &
           trim(stress_int_delta_env).ne.'0'


! Double factorial
      dfak(1)=1.d0
      DO i=2,41
        dfak(i)= dfak(i-1)*(2*i-1)
      ENDDO

      pi=4.0d+0*ATAN(1.0d+0)
      sqfp=SQRT(4.d0*pi)
      sqrt2=SQRT(2.D0)
      pi4=4.D0*pi

      GRR_SWITCH=.false.
      read_nlvdw=.false.
      read_vhalf=.false.
      lmbj_switch=.false.
      vsp_only=.false.
      core_add=.false.
      core_sub=.false.

      CALL GTFNAM(DEFFN,ERRFN)
      CALL ERRFLG(ERRFN,'Error in LAPW0')
#ifdef FFTW_OMP
     omp_error = fftw_init_threads()
     if (omp_error.eq.0_C_INT) then
        write(6,*) 'ERROR: FFTW thread initialization failed in LAPW0'
        stop 1
     endif
     call fftw_plan_with_nthreads(int(omp_get_max_threads(),C_INT))
#endif      
!print*, myid,deffn
!      if(myid.eq.1) then
         OPEN (1,FILE=DEFFN,STATUS='OLD',ERR=910)
!         print*, 'open with myid',myid
!      endif
!call  MPI_BARRIER(MPI_COMM_WORLD,ierr)
!if(myid.ne.1) then
!         OPEN (1,FILE=DEFFN,STATUS='OLD',ERR=910)
!         print*, 'open with myid',myid
!      endif
8000  CONTINUE
      READ (1,*,END=8001,ERR=960) IUNIT,FNAME,STATUS,FORM,IRECL
      if(iunit.eq.5) then
         if (fname((LEN_TRIM(fname)-2):LEN_TRIM(fname)) .eq. 'grr') then
            GRR_SWITCH=.true.
         else if (fname((LEN_TRIM(fname)-2):LEN_TRIM(fname)) .eq. 'vsp') then
            vsp_only=.true.
         end if
      endif 
      if(iunit.eq.66) then
        read_nlvdw=.true.
      endif 
      if(iunit.eq.68) then
        lmbj_switch=.true.
      endif
      if(iunit.eq.69) then
         read_vhalf=.true.
      endif
      if(iunit.eq.22) then
         if (fname((LEN_TRIM(fname)-2):LEN_TRIM(fname)) .eq. 'add') then
            core_add=.true.
         else if(fname((LEN_TRIM(fname)-2):LEN_TRIM(fname)) .eq. 'sub') then
            core_sub=.true.
         else
            call OUTERR('LAPW0','WRONG CEE SWITCH USED: ONLY ADD/SUB POSSIBLE')
            goto 999
         end if
      endif
#ifdef Parallel
      IF(iunit.EQ.6) THEN 
        il=LEN_TRIM(fname)
        i=myid
        k=1000
        DO j=1,3
          k=k/10
          fname(il+j:il+j)=CHAR(ICHAR("0")+i/k)
          i=i-i/k*k
        ENDDO
        if(restrict_output.le.myid) then
           FNAME='/dev/null'
           STATUS='OLD'
        endif
      ENDIF
#endif

      fnavn(IUNIT)=FNAME
      OPEN (IUNIT,FILE=FNAME,STATUS=STATUS,FORM=FORM,ERR=920)
      GOTO 8000
8001  CONTINUE
      CLOSE (1)
      CLOSE (10)
! test for .xcpot1q file
      xcpot1qq=.false.
      inquire(file='.xcpot1q',exist=xcpot1qq)
!print*,'xcpot1q',xcpot1qq

!call date_and_time(datexx,timexx)
!write(6,*) 'time',timexx(1:2),timexx(3:4),timexx(5:6) 
      CALL CPUTIM(TIME0)
      CALL walltim(TIME0w)

      supn=ZERO
      vxmix=0.5d0
      dis=0.d0
      READ(4,'(a)') parseline       
      read(parseline(6:80),*,err=3029) supn
!      read(4,'(a5,f8.2)',end=3029,err=3029) switc1,supn
      read(14,*,end=3029,err=3029)
      read(14,'(a)',end=3029,err=3029) parseline 
      do i=1,80
        if(parseline(i:i+5).eq.'GREED:') then
          READ(PARSELINE(i+6:i+12),*,IOSTAT=IOS) vxmix
          exit
        endif
      enddo
      read(14,'(a)',end=3029,err=3029) parseline 
      do i=1,80
        if(parseline(i:i).eq.'(') then    ! :dis
          READ(PARSELINE(i+1:i+8),*,IOSTAT=IOS) dis
          vxmix1=min(0.8,sqrt(0.04/dis))
          exit
        endif
      enddo
          vxmix=min(sqrt(vxmix*2.d0),vxmix1)
          vxmix=max(vxmix,0.24d0)
3029 CONTINUE
      
!     spin-polarized or not?
      READ(12,2030,END=3031)
      READ(13,2030,END=3031)
      JSPIN=2
      GOTO 3032
3031  JSPIN=1
3032  CONTINUE

      NCON=9
      call inputpars(PARSELINE,SWITCH,FACIFF,IFFX2,IFFY2,IFFZ2,IPRINT,SWITCH1, &
             SWITCH2,SWITCH3,R2V,LUSE2,LR2V,NCON,jspin)


      imfield=0
      read(7,*,end=5556,err=5556)nmod,nsp,natorb,Bextene
      if (nmod.eq.3) imfield=1
5556  continue

!     Insert Informational output
      if(myid .eq. 0)then
        call TestPrint()
      endif

! READ SCFDATA UNTIL BOTTOM AND INCREASE CYCLENUMBER BY ONE
      if(myid.eq.0) then
         READ(8,2044)  ISCF
         read(8,2043) lmmaxx
      endif
#ifdef Parallel
         call MPI_Bcast(iscf,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
         call MPI_Bcast(lmmaxx,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif
      ISCF=ISCF+1
      if(lmmaxx.eq.0) lmmaxx=ncom
      IF(ISCF.GT.999) ISCF=1
!     if((switch1.ne.' EECE') .and. ((iex .lt. 51) .or. (iex .gt. 55))) then     !THORUH
   !   if((switch1.ne.' EECE') .and. (EX_SWITCH(4:4).eq.'S') .and. &
   !      (VX_SWITCH(5:5) .ne. 'O')) then

      if((switch1.ne.' EECE') .and. (GRR_SWITCH .neqv. .true.) .and. .not. vsp_only) then
        IF(myid.EQ.0)  WRITE(21,705) ISCF,ISCF
      endif

      if(imfield.ne.0) then
        IF(myid.EQ.0)  write(21,5557) Bextene
        write(6,5557) Bextene
5557    format('        Interaction with Bext included, muB*Bext=',f12.7,' Ry')
      endif 

! Initialize E-field
      call EPOT0  (myid)

      CALL readstruct
      call haveinversion

#ifdef Parallel
      CALL begend(nat)
#else
      nstart(0)=1
      nstop(0)=nat
      coul=.TRUE.
      xc=.TRUE.
      icoul=0
#endif

      ALLOCATE(kkk(3,nsym))
      ALLOCATE(fc(ncom+3,nat),lmmtmx(nat),lmmult(2,ncom+3,nat))
      ALLOCATE(qel(nat),qq(ncom+3,nat),qqnew(ncom+3,nat))
      ALLOCATE(lmmax1(nat),lmmax(nat),lm(2,ncom+3,nat))
      ALLOCATE(vsph(nat))
      if (nstart(myid).le.nstop(myid)) &
        ALLOCATE(v(1:nrad,1:lmmaxx,nstart(myid):nstop(myid)))
      ALLOCATE(vcoul2(1:nrad,1:nat))
      ALLOCATE(taup(nsym))
      if (nstart(myid).le.nstop(myid)) &
        ALLOCATE(clm(1:nrad,1:lmmaxx,nstart(myid):nstop(myid)),&
        clmsp(1:nrad,1:lmmaxx,nstart(myid):nstop(myid),1:2))
      if (nstart(myid).le.nstop(myid)) &
        ALLOCATE(rhopw0(1:nrad+10,1:lmmaxx,nstart(myid):nstop(myid)),rhopw0r(1:nrad+10,nstart(myid):nstop(myid)))
      ALLOCATE(rhopw00(1:nrad+10,1:lmmaxx))
      allocate(clm_tmp(1:nrad))
      ALLOCATE(logefg(1:nat))

      allocate(elnum(nat),elnumsp(nat,2),fk(0:6,nat))
      fk=0.d0
      if (switch1 .eq. ' EECE') then
        allocate(jeece(nat))
        jeece = 0
        read(5,*) inod, numbeece
        read(5,'(a)') errmsg
        do j=1,numbeece
          read(5,*) i1, i2, i3
          do jatom=1,nat
            if (jatom .eq. i1) jeece(jatom) = 1
          enddo
        enddo
      endif

      logefg=.false.
      if (coul) v=zero

      DO jatom=1,nat
        vsph(jatom)=pi4/3.D0*rmt(jatom)**3
      ENDDO

      CALL LATGEN(NAT)

      IF(myid.EQ.0) THEN  
        WRITE(6,800)
        WRITE(6,805) TITLE
        WRITE(6,810) LATTIC
        nxcprint = 0
        if ((EX_SWITCH .eq. 'EX_PBE-ALPHA') .or. (VX_SWITCH .eq. 'VX_PBE-ALPHA') .or. &
            (VX_SWITCH .eq. 'VX_BR89') .or. (VX_SWITCH .eq. 'VX_MBJ') .or. &
            (VX_SWITCH .eq. 'VX_LMBJ')) then
          nxcprint = 1
        elseif ((EX_SWITCH .eq. 'EX_MPBE') .or. (VX_SWITCH .eq. 'VX_MPBE') .or. &
                (EX_SWITCH .eq. 'EX_MGGA_MS') .or. (VX_SWITCH .eq. 'VX_GBJ')) then
          nxcprint = 3
        endif
        IF (trim(EX_SWITCH) .eq. 'EX_LIBXC') THEN
          write(idstr,'(I5)') id_func_x
          EX_SWITCH_TMP = 'EX_LIBXC('//trim(X_SWITCH_LIBXC)//',id='//trim(adjustl(idstr))//')'
        ELSE
          EX_SWITCH_TMP = trim(EX_SWITCH)
        ENDIF
        IF (trim(EC_SWITCH) .eq. 'EC_LIBXC') THEN
          write(idstr,'(I5)') id_func_c
          EC_SWITCH_TMP = 'EC_LIBXC('//trim(C_SWITCH_LIBXC)//',id='//trim(adjustl(idstr))//')'
        ELSE
          EC_SWITCH_TMP = trim(EC_SWITCH)
        ENDIF
        IF (trim(VX_SWITCH) .eq. 'VX_LIBXC') THEN
          write(idstr,'(I5)') id_func_x
          VX_SWITCH_TMP = 'VX_LIBXC('//trim(X_SWITCH_LIBXC)//',id='//trim(adjustl(idstr))//')'
        ELSE
          VX_SWITCH_TMP = trim(VX_SWITCH)
        ENDIF
        IF (trim(VC_SWITCH) .eq. 'VC_LIBXC') THEN
          write(idstr,'(I5)') id_func_c
          VC_SWITCH_TMP = 'VC_LIBXC('//trim(C_SWITCH_LIBXC)//',id='//trim(adjustl(idstr))//')'
          if (VDISP_SWITCH .NE. 'VDISP_NONE') then
            VCVDISP_SWITCH_TMP = trim(VC_SWITCH_TMP)//' '//trim(VDISP_SWITCH)
          else
            VCVDISP_SWITCH_TMP = trim(VC_SWITCH_TMP)
          endif
        ELSE
          if (VDISP_SWITCH .NE. 'VDISP_NONE') then
             VCVDISP_SWITCH_TMP = trim(VC_SWITCH)//' '//trim(VDISP_SWITCH)
          else
             VCVDISP_SWITCH_TMP = trim(VC_SWITCH)
          endif
        ENDIF
        if (nxcprint .eq. 0) then
          WRITE(6,819) trim(EX_SWITCH_TMP), trim(EC_SWITCH_TMP), trim(VX_SWITCH_TMP), trim(VCVDISP_SWITCH_TMP)
        elseif (nxcprint .eq. 1) then
          WRITE(6,823) trim(EX_SWITCH_TMP), trim(EC_SWITCH_TMP), trim(VX_SWITCH_TMP), trim(VCVDISP_SWITCH_TMP), XC1
        elseif (nxcprint .eq. 2) then
          WRITE(6,824) trim(EX_SWITCH_TMP), trim(EC_SWITCH_TMP), trim(VX_SWITCH_TMP), trim(VCVDISP_SWITCH_TMP), XC1, XC2
        elseif (nxcprint .eq. 3) then
          WRITE(6,825) trim(EX_SWITCH_TMP), trim(EC_SWITCH_TMP), trim(VX_SWITCH_TMP), trim(VCVDISP_SWITCH_TMP), XC1, XC2, XC3
       endif
        IF((EX_SWITCH .eq. 'EX_GRR') .or. (VX_SWITCH .eq. 'VX_GRR')) THEN
          write(6,*) 'xcconst,ycconst for IXC=GRR (Average for Gradient Rho over Rho)',xcconst,ycconst
        ENDIF
        IF((VX_SWITCH .eq. 'VX_MBJ') .or. (VX_SWITCH .eq. 'VX_SMBJ') .or. &
           (VX_SWITCH .eq. 'VX_GBJ') .or. (id_func_x .eq. 208)) THEN
          write(6,*) 'xca,xcb,xcp,v-mixing,<grad rho/rho>,c for IXC=modified Becke-Johnson:'
          write(6,'(6f10.3)') xca,xcb,xcp,vxmix,grr,XCA + XCB*GRR**XCP
        ELSEIF (VX_SWITCH .eq. 'VX_LMBJ') THEN
          write(6,*) 'xca,xcb,xcp,v-mixing,r_th for IXC=local modified Becke-Johnson:'
          write(6,'(5f10.3)') xca,xcb,xcp,vxmix,r_th
        ENDIF
        IF (dostagfield .eqv. .true.) THEN
          write(6,*) 'A staggered field is applied'
        ENDIF
      ENDIF
      if(supn.ne.0.d0) then
         if(GRR_SWITCH .eqv. .true.) then
             supn=0.d0
             write(6,*) 'CHARGED CELL switched off for a -grr calculation'
         else 
             write(6,*) ' CHARGED CELL with ',supn
         endif
      endif
      if(supn.ne.0.d0 .and. SWITCH3.eq.'YES' ) then
         write(6,*) 'Charged cells AND E-tot decomposition not possible'
         CALL OUTERR('LAPW0','Charged cell and SWITCH3 not possible')
         call abort_parallel
         STOP 'Charged cell and SWITCH3 not possible'
      endif

      if(lattic(1:1).eq.'H'.or.lattic(1:1).eq.'R') then
         WRITE(6,820) AA,BB,CC
      else
         WRITE(6,820) AA,BB,CC,ALPHA
      endif
      WRITE(6,840) NAT,ndif
      WRITE(6,850) IREL
      WRITE(6,860) ISCF
      if(myid .eq. 0)then
        WRITE(21,840) NAT,ndif
        WRITE(21,805) TITLE
        WRITE(21,810) LATTIC
        if(lattic(1:1).eq.'H'.or.lattic(1:1).eq.'R') then
           WRITE(21,820) AA,BB,CC
        else
           WRITE(21,820) AA,BB,CC,ALPHA
        endif
        WRITE(21,821) VOL
        nxcprint = 0
        if ((EX_SWITCH .eq. 'EX_PBE-ALPHA') .or. (VX_SWITCH .eq. 'VX_PBE-ALPHA') .or. &
            (VX_SWITCH .eq. 'VX_BR89') .or. (VX_SWITCH .eq. 'VX_MBJ') .or. &
            (VX_SWITCH .eq. 'VX_LMBJ')) then
          nxcprint = 1
        elseif ((EX_SWITCH .eq. 'EX_MPBE') .or. (VX_SWITCH .eq. 'VX_MPBE') .or. &
                (EX_SWITCH .eq. 'EX_MGGA_MS') .or. (VX_SWITCH .eq. 'VX_GBJ')) then
          nxcprint = 3
        endif
        IF (trim(EX_SWITCH) .eq. 'EX_LIBXC') THEN
          write(idstr,'(I5)') id_func_x
          EX_SWITCH_TMP = 'EX_LIBXC('//trim(X_SWITCH_LIBXC)//',id='//trim(adjustl(idstr))//')'
        ELSE
          EX_SWITCH_TMP = trim(EX_SWITCH)
        ENDIF
        IF (trim(EC_SWITCH) .eq. 'EC_LIBXC') THEN
          write(idstr,'(I5)') id_func_c
          EC_SWITCH_TMP = 'EC_LIBXC('//trim(C_SWITCH_LIBXC)//',id='//trim(adjustl(idstr))//')'
        ELSE
          EC_SWITCH_TMP = trim(EC_SWITCH)
        ENDIF
        IF (trim(VX_SWITCH) .eq. 'VX_LIBXC') THEN
          write(idstr,'(I5)') id_func_x
          VX_SWITCH_TMP = 'VX_LIBXC('//trim(X_SWITCH_LIBXC)//',id='//trim(adjustl(idstr))//')'
        ELSE
          VX_SWITCH_TMP = trim(VX_SWITCH)
        ENDIF
        IF (trim(VC_SWITCH) .eq. 'VC_LIBXC') THEN
          write(idstr,'(I5)') id_func_c
          VC_SWITCH_TMP = 'VC_LIBXC('//trim(C_SWITCH_LIBXC)//',id='//trim(adjustl(idstr))//')'
          if (VDISP_SWITCH .NE. 'VDISP_NONE') then
            VCVDISP_SWITCH_TMP = trim(VC_SWITCH_TMP)//' '//trim(VDISP_SWITCH)
          else
            VCVDISP_SWITCH_TMP = trim(VC_SWITCH_TMP)
          endif
        ELSE
          if (VDISP_SWITCH .NE. 'VDISP_NONE') then
             VCVDISP_SWITCH_TMP = trim(VC_SWITCH)//' '//trim(VDISP_SWITCH)
          else
             VCVDISP_SWITCH_TMP = trim(VC_SWITCH)
          endif
        ENDIF
        IF(MGGA_TAU_SWITCH .eqv. .TRUE.) WRITE(21,826)
        if (nxcprint .eq. 0) then
          WRITE(21,819) trim(EX_SWITCH_TMP), trim(EC_SWITCH_TMP), trim(VX_SWITCH_TMP), trim(VCVDISP_SWITCH_TMP)
        elseif (nxcprint .eq. 1) then
          WRITE(21,823) trim(EX_SWITCH_TMP), trim(EC_SWITCH_TMP), trim(VX_SWITCH_TMP), trim(VCVDISP_SWITCH_TMP), XC1
        elseif (nxcprint .eq. 2) then
          WRITE(21,824) trim(EX_SWITCH_TMP), trim(EC_SWITCH_TMP), trim(VX_SWITCH_TMP), trim(VCVDISP_SWITCH_TMP), XC1, XC2
        elseif (nxcprint .eq. 3) then
          WRITE(21,825) trim(EX_SWITCH_TMP), trim(EC_SWITCH_TMP), trim(VX_SWITCH_TMP), trim(VCVDISP_SWITCH_TMP), XC1, XC2, XC3
        endif
        IF((EX_SWITCH .eq. 'EX_GRR') .or. (VX_SWITCH .eq. 'VX_GRR')) THEN
          write(21,*) 'xcconst,ycconst for IXC=GRR (Average for Gradient Rho over Rho)',xcconst,ycconst
        ENDIF
        IF((VX_SWITCH .eq. 'VX_MBJ') .or. (VX_SWITCH .eq. 'VX_SMBJ') .or. &
           (VX_SWITCH .eq. 'VX_GBJ') .or. (id_func_x .eq. 208)) THEN
          write(21,'(":MBJ  : xca,xcb,xcp,v-mixing,<grad rho/rho>,c for IXC=modified Becke-Johnson:")')
          write(21,'(6f10.3)') xca,xcb,xcp,vxmix,grr,XCA + XCB*GRR**XCP
        ELSEIF (VX_SWITCH .eq. 'VX_LMBJ') THEN
          write(21,'(":LMBJ  : xca,xcb,xcp,v-mixing,r_th for IXC=local modified Becke-Johnson:")')
          write(21,'(5f10.3)') xca,xcb,xcp,vxmix,r_th
        ENDIF
        IF (dostagfield .eqv. .true.) THEN
          write(21,*) 'A staggered field is applied'
        ENDIF
 822    format( ':INFO  :        CHARGED CELL with ',f7.3,/,8x, &
        'an energy correction like C Q**2/(L eps) is not included (PRB51,4014; PRB73,35215)')
        if(supn.gt.1.d-6) write(21,822) supn
        if(supn.lt.-1.d-6) write(21,822) supn
        WRITE(21,850) IREL
      endif
      delcc=supn/vol

      INDEX=0
      if(myid .eq. 0)then
         DO JATOM=1,NAT
            WRITE(6,1060) ANAME(JATOM)
            DO JROW=1,3
               WRITE(6,1070) ( ROTLOC(JCOL,JROW,JATOM),JCOL=1,3 )
            ENDDO
            DO MU=1,MULT(JATOM)
               INDEX=INDEX+1
               WRITE(6,1080) MU,( POS(JX,INDEX),JX=1,3 )
               DO JROW=1,3
                  WRITE(6,1070) (ROTIJ(JCOL,JROW,INDEX),JCOL=1,3)
               ENDDO
            ENDDO
         ENDDO
      endif
!
! read total charge density CLM from CLMSUM
! norm of CLM(0,0)=CLM(0,0)*R*R*SQRT(4.*PI), else CLM(0,0)*R*R
!

      DO JATOM=1,NAT
        if(myid.eq.0) READ(8,2000) LMMAX(JATOM)
#ifdef Parallel
         call MPI_Bcast(LMMAX(JATOM),1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif
        LMMAX1(JATOM)=LMMAX(JATOM)
        IF(LMMAX1(JATOM).GT.NCOM+3) GOTO 970
        DO LM1=1,LMMAX(JATOM)
          if(myid.eq.0) READ(8,2010) L,M
#ifdef Parallel
         call MPI_Bcast(L,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
         call MPI_Bcast(M,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif
          if(l.gt.lmax2) goto 971
          LM(1,LM1,JATOM)=L
          LM(2,LM1,JATOM)=M
          IF(M.NE.0) THEN
            LMMAX1(JATOM)=LMMAX1(JATOM)+1
            IF(LMMAX1(JATOM).GT.NCOM+3) GOTO 970
            LM(1,LMMAX1(JATOM),JATOM)=L
            LM(2,LMMAX1(JATOM),JATOM)=-M
          endif
          if(myid.eq.0) then
            if (jatom.ge.nstart(myid).and.jatom.le.nstop(myid)) then
              READ(8,2021) ( CLM(J,LM1,JATOM), J=1,JRI(JATOM))
            else
              read(8,2021) ( CLM_tmp(J), J=1,JRI(JATOM))
#ifdef Parallel
               DO i=1,icoul
                 DO j=nstart(i),nstop(i)
                   if(jatom.eq.j) then
                      myid1=i
                      goto 3478
                   endif
                 enddo
               enddo
               print*, 'error: myid1 not defined', jatom,myid
 3478          continue
!print*, 'sending clm to', myid1, jatom,lm1
              CALL MPI_SEND(clm_tmp, jri(jatom), MPI_DOUBLE_PRECISION, myid1, 0, comm_coul, ierr)
#endif
            endif
          endif
#ifdef Parallel
            IF(myid.nE.0) THEN
              IF(jatom.ge.nstart(myid).and.jatom.le.nstop(myid)) THEN
!print*, 'should be receiving clm ', myid, jatom,lm1
                CALL MPI_RECV(CLM(1,LM1,JATOM),jri(jatom) , MPI_DOUBLE_PRECISION, 0, &
                  MPI_ANY_TAG, comm_coul, statusmpi, ierr)
!print*, 'received clm ', myid, jatom
              endif
            endif
#endif
          if(myid.eq.0) READ(8,2031)
        ENDDO
        if(myid.eq.0) READ(8,1980)
        if(myid.eq.0) READ(8,2030)
      ENDDO

      IF(jspin.EQ.1) THEN
        IF(myid.EQ.0)  WRITE(6,*) '      NON-SPINPOLARIZED CALCULATION'
        IF(myid.EQ.0)  WRITE(21,*) '      NON-SPINPOLARIZED CALCULATION'
      ELSE
        IF(myid.EQ.0)  WRITE(6,*) '      SPINPOLARIZED CALCULATION'
        IF(myid.EQ.0)  WRITE(21,*) '      SPINPOLARIZED CALCULATION'
      endif

      DO jatom=1,nat

        IF(jspin.EQ.2) THEN

          if(myid.eq.0) READ(12,2000) lmma                                             
#ifdef Parallel
         call MPI_Bcast(LMMA,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif
          lmmax(jatom)=lmma                  !eece
          LMMAX1(JATOM)=LMMAX(JATOM)         !eece
          if(myid.eq.0) READ(13,2000) lmma                                             
          DO lm1=1 ,lmma                                             
            if(myid.eq.0) READ(12,2010) l,m                                              
            if(myid.eq.0) READ(13,2010) l,m                                              
#ifdef Parallel
         call MPI_Bcast(L,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
         call MPI_Bcast(M,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif
            if(l.gt.lmax2) goto 971          !eece   use lm-list from unit 12
            LM(1,LM1,JATOM)=L
            LM(2,LM1,JATOM)=M
            IF(M.NE.0) THEN
              LMMAX1(JATOM)=LMMAX1(JATOM)+1
              IF(LMMAX1(JATOM).GT.NCOM+3) GOTO 970
              LM(1,LMMAX1(JATOM),JATOM)=L
              LM(2,LMMAX1(JATOM),JATOM)=-M
            endif
            if(myid.eq.0) then
              if (jatom.ge.nstart(myid).and.jatom.le.nstop(myid)) then
                READ(12,2021) ( CLMsp(J,LM1,JATOM,1), J=1,JRI(JATOM))
                READ(13,2021) ( CLMsp(J,LM1,JATOM,2), J=1,JRI(JATOM))
              else
                read(12,2021) ( CLM_tmp(J), J=1,JRI(JATOM))
#ifdef Parallel
               DO i=1,icoul
                 DO j=nstart(i),nstop(i)
                   if(jatom.eq.j) then
                      myid1=i
                      goto 3479
                   endif
                 enddo
               enddo
               print*, 'error: myid1 not defined', jatom,myid
 3479          continue
!print*, 'sending clm to', myid1, jatom,lm1
              CALL MPI_SEND(clm_tmp, jri(jatom), MPI_DOUBLE_PRECISION, myid1, 0, comm_coul, ierr)
#endif
                read(13,2021) ( CLM_tmp(J), J=1,JRI(JATOM))
#ifdef Parallel
               DO i=1,icoul
                 DO j=nstart(i),nstop(i)
                   if(jatom.eq.j) then
                      myid1=i
                      goto 3477
                   endif
                 enddo
               enddo
               print*, 'error: myid1 not defined', jatom,myid
 3477          continue
!print*, 'sending clm to', myid1, jatom,lm1
              CALL MPI_SEND(clm_tmp, jri(jatom), MPI_DOUBLE_PRECISION, myid1, 0, comm_coul, ierr)
#endif
              endif
            endif
#ifdef Parallel
            IF(myid.nE.0) THEN
              IF(jatom.ge.nstart(myid).and.jatom.le.nstop(myid)) THEN
!print*, 'should be receiving clm ', myid, jatom,lm1
                CALL MPI_RECV(CLMsp(1,LM1,JATOM,1),jri(jatom) , MPI_DOUBLE_PRECISION, 0, &
                  MPI_ANY_TAG, comm_coul, statusmpi, ierr)
                CALL MPI_RECV(CLMsp(1,LM1,JATOM,2),jri(jatom) , MPI_DOUBLE_PRECISION, 0, &
                  MPI_ANY_TAG, comm_coul, statusmpi, ierr)
!print*, 'received clm ', myid, jatom
              endif
            endif
#endif
            if(myid.eq.0) READ(13,2031)                                                  
            if(myid.eq.0) READ(12,2031)                                                  
            if(switch1.eq.' EECE')then
              if (jatom.ge.nstart(myid).and.jatom.le.nstop(myid)) then
                do j=1,jri(jatom)
                  if (l.eq.0)then
!!  if (j.eq.0) print *,'stupid print to fix ifort bug' ! fix for stupid overoptimization
                    clmsp(j,lm1,jatom,1)=clmsp(j,lm1,jatom,1)*sqfp
                    clmsp(j,lm1,jatom,2)=clmsp(j,lm1,jatom,2)*sqfp
                  endif
                    clm(j,lm1,jatom)=clmsp(j,lm1,jatom,1)+clmsp(j,lm1,jatom,2)
                enddo
              endif
           endif
          ENDDO
          if(myid.eq.0) then
            READ(12,1980)                                                  
            READ(12,2030)                                                  
            READ(13,1980)                                                  
            READ(13,2030)
          endif 
                                                 
        ELSE
                                                            
          if (jatom.ge.nstart(myid).and.jatom.le.nstop(myid)) then
            DO  lm1=1,lmmax(jatom)                                      
              DO  j=1,jri(jatom)                                          
                clmsp(j,lm1,jatom,1)=clm(j,lm1,jatom)*0.5D0       
                clmsp(j,lm1,jatom,2)=clm(j,lm1,jatom)*0.5D0      
              ENDDO
            ENDDO
          endif

        endif

      ENDDO
! Only for functionals using the kinetic-energy density
!      IF (IEX.eq.12.or.IEX.eq.27.or.IEX.eq.28.or.IEX.eq.29.or.IEX.eq.30.or.IEX.eq.34.or. &         !THORUH
!          IEX.eq.38.or.IEX.eq.39.or.IEX.eq.40.or.IEX.eq.41.or.IEX.eq.42.or.IEX.eq.43) THEN         !THORUH
      IF (VRESP_SWITCH .eqv. .true.) THEN
         if (nstart(myid).le.nstop(myid)) &
              allocate(vclmsp(1:nrad,1:lmmaxx,nstart(myid):nstop(myid),1:2))
         IF(myid.EQ.0)  write(21,*) '      Calculation includes Kinetic Energy Density'
         
         DO jatom=1,nat
            
            IF(JSPIN.EQ.2) THEN
               
               if(jatom.eq.1) then
                  READ(29,2030) 
                  READ(30,2030)
               endif
               
               READ(29,2000) LMMA                                             
               READ(30,2000) LMMA                                             
               DO LM1=1 ,LMMA                                             
                  READ(29,2010) L,M
                  READ(30,2010) L,M
                  if (jatom.ge.nstart(myid).and.jatom.le.nstop(myid)) then
                     READ(29,2021) ( VCLMSP(J,LM1,JATOM,1), J=1,JRI(JATOM))
                     READ(30,2021) ( VCLMSP(J,LM1,JATOM,2), J=1,JRI(JATOM))
                  else
                     read(30,2021) ( CLM_tmp(J), J=1,JRI(JATOM))
                     read(29,2021) ( CLM_tmp(J), J=1,JRI(JATOM))
                  endif
                  READ(30,2031)
                  READ(29,2031)
               ENDDO
               READ(29,1980)                                                  
               READ(29,2030)                                                  
               READ(30,1980)                                                  
               READ(30,2030)
               
            ELSE
               
               if(jatom.eq.1) then
                  READ(28,2030,IOSTAT=IOS) 
                  if(IOS.ne.0) then
                     STOP 'ERROR reading vrespsum'
                  endif
               endif
               READ(28,2000) LMMA                                             
               DO LM1=1 ,LMMA                                             
                  READ(28,2010) L,M
                  if (jatom.ge.nstart(myid).and.jatom.le.nstop(myid)) then
                     READ(28,2021) ( VCLMSP(J,LM1,JATOM,1), J=1,JRI(JATOM))
                  else
                     read(28,2021) ( CLM_tmp(J), J=1,JRI(JATOM))
                  endif
                  READ(28,2031)
               ENDDO
               READ(28,1980)                                                  
               READ(28,2030)                                                  
               if (jatom.ge.nstart(myid).and.jatom.le.nstop(myid)) then
                  DO LM1=1,LMMAX(JATOM)
                     DO J=1,JRI(JATOM)                                          
                        VCLMSP(J,LM1,JATOM,1)=VCLMSP(J,LM1,JATOM,1)*0.5D0 
                        VCLMSP(J,LM1,JATOM,2)=VCLMSP(J,LM1,JATOM,1)
                     ENDDO
                  ENDDO
               endif
               
            endif
            
! Simpler Normalization  for -eece -vresp ?
         if (jatom.ge.nstart(myid).and.jatom.le.nstop(myid)) then
            if (switch1 .EQ. ' EECE') &
                 VCLMSP(1:JRI(JATOM),1,JATOM,1:2)= VCLMSP(1:JRI(JATOM),1,JATOM,1:2)*sqfp
         endif
         
      ENDDO
      
!     end if(vresp_switch)
   endif


      IF (read_tau_switch) THEN
         if (nstart(myid).le.nstop(myid)) &
              allocate(slater_tau_lm(1:nrad,1:lmmaxx,nstart(myid):nstop(myid),1:2))
         slater_tau_lm = 0.d0
         IF(myid.EQ.0)  write(21,*) '      Calculation includes Kinetic Energy Density'
         IF(myid.EQ.0)  write(21,*) '      Using Slater form of KED from file case.tausum' 

         DO jatom=1,nat
            
            IF(JSPIN.EQ.2) THEN
               
               if(jatom.eq.1) then
                  READ(39,2030) 
                  READ(40,2030)
               endif
               READ(39,2000) LMMA                                             
               READ(40,2000) LMMA                                             
               DO LM1=1 ,LMMA                                             
                  READ(39,2010) L,M
                  READ(40,2010) L,M
                  if (jatom.ge.nstart(myid).and.jatom.le.nstop(myid)) then
                     READ(39,2021) ( slater_tau_lm(J,LM1,JATOM,1), J=1,JRI(JATOM))
                     READ(40,2021) ( slater_tau_lm(J,LM1,JATOM,2), J=1,JRI(JATOM))
                  else
                     read(40,2021) ( CLM_tmp(J), J=1,JRI(JATOM))
                     read(39,2021) ( CLM_tmp(J), J=1,JRI(JATOM))
                  endif
                  READ(40,2031)
                  READ(39,2031)
               ENDDO
               READ(39,1980)                                                  
               READ(39,2030)                                                  
               READ(40,1980)                                                  
               READ(40,2030)
               
            ELSE
               
               if(jatom.eq.1) then
                  READ(38,2030,IOSTAT=IOS) 
                  if(IOS.ne.0) then
                     STOP 'ERROR reading tausum'
                  endif
               endif
               READ(38,2000) LMMA                                             
               DO LM1=1 ,LMMA                                             
                  READ(38,2010) L,M
                  if (jatom.ge.nstart(myid).and.jatom.le.nstop(myid)) then
                     READ(38,2021) ( slater_tau_lm(J,LM1,JATOM,1), J=1,JRI(JATOM))
                  else
                     read(38,2021) ( CLM_tmp(J), J=1,JRI(JATOM))
                  endif
                  READ(38,2031)
               ENDDO
               READ(38,1980)                                                  
               READ(38,2030)                                                  
               if (jatom.ge.nstart(myid).and.jatom.le.nstop(myid)) then
                  DO LM1=1,LMMAX(JATOM)
                     DO J=1,JRI(JATOM)                                          
                        slater_tau_lm(J,LM1,JATOM,1)=slater_tau_lm(J,LM1,JATOM,1)*0.5D0 
                        slater_tau_lm(J,LM1,JATOM,2)=slater_tau_lm(J,LM1,JATOM,1)
                     ENDDO
                  ENDDO
               endif
               
            endif
            
! Simpler Normalization  for -eece -vresp ?
         if (jatom.ge.nstart(myid).and.jatom.le.nstop(myid)) then
            if (switch1 .EQ. ' EECE') &
                 slater_tau_lm(1:JRI(JATOM),1,JATOM,1:2)= slater_tau_lm(1:JRI(JATOM),1,JATOM,1:2)*sqfp
         endif
         
      ENDDO
! Convert Ry -> Ha
      if(allocated(slater_tau_lm)) &
           slater_tau_lm(:,:,nstart(myid):nstop(myid),:) = 0.5d0*slater_tau_lm(:,:,nstart(myid):nstop(myid),:)
!     end if(read_tau_switch)
   endif


!!$      ALLOCATE(efgb(1:nat))
!!$      DO jatom=1,nat
!!$        efgb(jatom)%v20=zero ;efgb(jatom)%v21=zero; efgb(jatom)%v22=zero
!!$        efgb(jatom)%v21m=zero ; efgb(jatom)%v22m=zero
!!$        efgb(jatom)%v20srf=zero ; efgb(jatom)%v21srf=zero ; efgb(jatom)%v22srf=zero
!!$        efgb(jatom)%v21msr=zero ; efgb(jatom)%v22msr=zero ; efgb(jatom)%vzz=zero
!!$        efgb(jatom)%winkel=zero ; efgb(jatom)%eta=zero
!!$        efgb(jatom)%qmat(1:3,1:3)=zero
!!$        efgb(jatom)%tens(3)=zero
!!$        efgb(jatom)%eivec(1:3,1:3)=zero
!!$      ENDDO
      ALLOCATE(efg_v20(1:nat),efg_v21(1:nat),efg_v22(1:nat),efg_v21m(1:nat),&
           efg_v22m(1:nat),efg_v20srf(1:nat),efg_v21srf(1:nat),efg_v22srf(1:nat),&
           efg_v21msr(1:nat),efg_v22msr(1:nat),efg_vzz(1:nat),efg_winkel(1:nat),&
           efg_eta(1:nat),efg_qmat(1:3,1:3,1:nat),efg_tens(1:3,1:nat),efg_eivec(1:3,1:3,1:nat))



        efg_v20=zero 
        efg_v21=zero
        efg_v22=zero
        efg_v21m=zero 
        efg_v22m=zero
        efg_v20srf=zero
        efg_v21srf=zero
        efg_v22srf=zero
        efg_v21msr=zero
        efg_v22msr=zero
        efg_vzz=zero
        efg_winkel=zero
        efg_eta=zero
        efg_qmat=zero
        efg_tens=zero
        efg_eivec=zero

      ALLOCATE(sigma(2,1:nat))
      sigma=zero

      if (nstart(myid).le.nstop(myid)) &
        ALLOCATE(exc(1:nrad,1:lmmaxx,nstart(myid):nstop(myid),1:2),&
        vxc(1:nrad,1:lmmaxx,nstart(myid):nstop(myid),1:2))

      if (nstart(myid).le.nstop(myid) .and. MGGA_TAU_SWITCH) &
           ALLOCATE(vtau(1:nrad,1:lmmaxx,nstart(myid):nstop(myid),1:2))

      if (nstart(myid).le.nstop(myid) .and. MGGA_LAPL_SWITCH) &
           ALLOCATE(vlapl(1:nrad,1:lmmaxx,nstart(myid):nstop(myid),1:2))

      IF((VX_SWITCH.eq.'VX_BR89') .or. (VX_SWITCH.eq.'VX_SLATER') .or. &
         (VX_SWITCH.eq.'VX_MBJ') .or. (VX_SWITCH.eq.'VX_SMBJ') .or. &
         (VX_SWITCH.eq.'VX_GBJ') .or. (VX_SWITCH.eq.'VX_KLI') .or. &
         (VX_SWITCH.eq.'VX_R2V') .or. (VDISP_SWITCH.eq.'VDISP_NLVDW') .or. &
         (VHALF_SWITCH.ne.'VHALF_NONE') .or. (VX_SWITCH.eq.'VX_LMBJ') .or. &
         (id_func_x .eq. 207) .or. (id_func_x .eq. 208) .or. (id_func_x .eq. 209)) THEN
        if (nstart(myid).le.nstop(myid)) then
          ALLOCATE(vxcold(1:nrad,1:lmmaxx,nstart(myid):nstop(myid),1:2))
          vxcold=zero
        endif
      endif
      IF((VX_SWITCH.eq.'VX_SLATER') .or. (VX_SWITCH.eq.'VX_SMBJ') .or. &
         (VX_SWITCH.eq.'VX_KLI') .or. (EX_SWITCH(1:12).eq.'EX_DELTAGLLB') .or. &
         (VX_SWITCH.eq.'VX_LMBJ')) THEN
        if (nstart(myid).le.nstop(myid)) then
          ALLOCATE(vxcnonloc(1:nrad,1:lmmaxx,nstart(myid):nstop(myid),1:2))
          vxcnonloc=zero
        endif
      endif

      IF(myid.EQ.0)  write(6,'("allocation of V,clm,clmsp,exc,vxc:",i5," MB")') &
        nrad*lmmaxx*(nstop(myid)-nstart(myid)+1)*8*8/1000000

!        print*,'after clm alloc'
!        do j1=1,100000
!        do j2=1,100000
!        do i=1,10000
!        if(j1*j2.eq.1.and.j2.eq.j1.and.i.eq.1) print*,j1
!        enddo
!        enddo
!        enddo
      if (coul) then
        vxc=zero
        exc=zero
      endif
!
! read plane wave density RHOK(K) from CLMSUM and average YLM
!
      if(myid.eq.0) then
         read(8,'(/,a19)') nkktext
         read(nkktext,'(9x,i10)',err=6767) nkk
         goto 6768
6767     read(nkktext,'(13x,i6)') nkk
6768  continue
      endif
#ifdef Parallel
         call MPI_Bcast(nkk,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif
      if(switch1 .EQ. ' EECE') nkk=1        !Cleanup multipole calculation
      allocate (absk(nkk))   

!====allocating variables for stress=================
     IF(SW_STRESS .eq. 'YES') THEN
      allocate( DECOUL_POT_PW(1:nkk,1:9) ) !contains equation 6.15
      allocate( SGGP(1:nkk,1:9) )          !contains (k+G)a*(k+G)b      
      DECOUL_POT_PW  = zeroc            
      SGGP           = zero      
      
      allocate( SGGP_GGA(nkk,1:9) )
      allocate( gga_int(1:nkk,1:9) )       !contains gga correction term
      gga_int        = zeroc
      SGGP_GGA       = zero
      
      time_electro   = zero
      
     ENDIF 
!====================================================      
         if (nstart(myid).le.nstop(myid).and.GGA_SWITCH .eqv. .TRUE. &
             .and. SW_STRESS .eq. 'YES') then
             ALLOCATE(sigma_str(1:nrad,1:lmmaxx,nstart(myid):nstop(myid),1:3))         
         else
             ALLOCATE(sigma_str(1,1,1,1))         
         endif
      
      IF((VX_SWITCH.eq.'VX_BR89') .or. (VX_SWITCH.eq.'VX_SLATER') .or. &
         (VX_SWITCH.eq.'VX_MBJ') .or. (VX_SWITCH.eq.'VX_SMBJ') .or. &
         (VX_SWITCH.eq.'VX_GBJ') .or. (VX_SWITCH.eq.'VX_KLI') .or. &
         (VX_SWITCH.eq.'VX_R2V') .or. (VDISP_SWITCH.eq.'VDISP_NLVDW') .or. &
         (VHALF_SWITCH.ne.'VHALF_NONE') .or. (VX_SWITCH.eq.'VX_LMBJ') .or. &
         (id_func_x .eq. 207) .or. (id_func_x .eq. 208) .or. (id_func_x .eq. 209)) THEN
        allocate (vkvxcold(1:nkk,1:2))
        VKVXCOLD(1:nkk,1:2)=ZEROC
        nkk2=nkk
      else
        allocate (vkvxcold(1,1:2))
        nkk2=1
      endif

      IF((VX_SWITCH.eq.'VX_SLATER') .or. (VX_SWITCH.eq.'VX_SMBJ') .or. &
         (VX_SWITCH.eq.'VX_KLI') .or. (EX_SWITCH(1:12).eq.'EX_DELTAGLLB') .or. &
         (VX_SWITCH.eq.'VX_LMBJ')) THEN
        allocate (vkvxcnonloc(1:nkk,1:2))
        VKVXCNONLOC(1:nkk,1:2)=ZEROC
        nkk3=nkk
      else
        allocate (vkvxcnonloc(1,1:2))
        nkk3=1
      endif

      allocate (kzz(3,nkk),inst(nkk))
      allocate (rhok(1:nkk))
      IF(myid.EQ.0)  write(6,'("allocation of X(nkk) arrays:",i10,i5," MB")') &
        nkk,nkk*16*10/1000000
  
      if(myid.eq.0) READ(8,2071)  ( KZZ(JX,1), JX=1,3 ), RHOKIN
#ifdef Parallel
         call MPI_Bcast(kzz(1,1),3,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
         call MPI_Bcast(rhokin,1,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
#endif

      RHOK(1)=RHOKIN
      ABSK(1)=0.0
      INST(1)=1

      IFF1=0
      IFF2=0
      IFF3=0
      DO 40 J=2,NKK
        if(myid.eq.0) READ(8,2071)  ( KZZ(JX,J), JX=1,3 ), RHOKIN
#ifdef Parallel
         call MPI_Bcast(kzz(1,J),3,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
         call MPI_Bcast(rhokin,1,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
#endif
        RHOK(J)=RHOKIN
        AV(1)=BR1(1,1)*KZZ(1,J)+BR1(1,2)*KZZ(2,J)+BR1(1,3)*KZZ(3,J)
        AV(2)=BR1(2,1)*KZZ(1,J)+BR1(2,2)*KZZ(2,J)+BR1(2,3)*KZZ(3,J)
        AV(3)=BR1(3,1)*KZZ(1,J)+BR1(3,2)*KZZ(2,J)+BR1(3,3)*KZZ(3,J)
        ABSK(J)=SQRT(AV(1)*AV(1)+AV(2)*AV(2)+AV(3)*AV(3))
        CALL STERN(J,IND,KKK,TAUP)
        INST(J)=IND
        DO 42 JJ=1,IND
!
! transfer carthesian into internal coord.
! find max k-vector
!
          IF(ORTHO) THEN
            II1=NINT((AVEC(1,1)*KKK(1,JJ)+AVEC(1,2)*KKK(2,JJ)+ &
              AVEC(1,3)*KKK(3,JJ))/AA)
            II2=NINT((AVEC(2,1)*KKK(1,JJ)+AVEC(2,2)*KKK(2,JJ)+ &
              AVEC(2,3)*KKK(3,JJ))/BB)
            II3=NINT((AVEC(3,1)*KKK(1,JJ)+AVEC(3,2)*KKK(2,JJ)+ &
              AVEC(3,3)*KKK(3,JJ))/CC)
          ELSE
            II1=KKK(1,JJ)
            II2=KKK(2,JJ)
            II3=KKK(3,JJ)
          endif
          IF(II1.GT.IFF1) IFF1=II1
          IF(II2.GT.IFF2) IFF2=II2
          IF(II3.GT.IFF3) IFF3=II3
          IF(KKK(1,JJ).GT.IFF1) IFF1=KKK(1,JJ)
          IF(KKK(2,JJ).GT.IFF2) IFF2=KKK(2,JJ)
          IF(KKK(3,JJ).GT.IFF3) IFF3=KKK(3,JJ)
  42    CONTINUE
  40  CONTINUE

! Set E-field values if needed
      CALL EPOTINIT(NKK,myid)

! set IFF values for Fast Fourier calculation
! check dimensions
!
      IFF1=(IFF1+1)*2
      IFF2=(IFF2+1)*2
      IFF3=(IFF3+1)*2
      if(iffx2.gt.0) iff1=iffx2
      if(iffy2.gt.0) iff2=iffy2
      if(iffz2.gt.0) iff3=iffz2
      IFF1=IFF1 * faciff
      IFF2=IFF2 * faciff
      IFF3=IFF3 * faciff
      IFF1T=IFF1 
      IFF2T=IFF2
      IFF3T=IFF3
      CALL IFFLIM(IFF1T,IFF1)
      CALL IFFLIM(IFF2T,IFF2)
      CALL IFFLIM(IFF3T,IFF3)

      iff3_g=0

      if (xc) then

!        call init_fftw_para(iff1,iff2,iff3)
!#ifdef Parallel
!        if (myid_coul.eq.0.and.local_last_start.ne.0) then
!          write(0,*) 'FFTW problem: myid_coul.eq.0.and.local_last_start.ne.0)'
!          call mpi_abort(MPI_COMM_WORLD,1,ierr)
!          stop
!        endif
!#endif
#ifdef Parallel
        call prepare_parallel_ffts (iff1,iff2,iff3)
        iff3_g=iff3 
        ifft3_g=iff3 
        iff3=local_nlast        

        if (myid_coul.eq.0.and.local_last_start.ne.0) then
          write(0,*) 'FFTW problem: myid_coul.eq.0.and.local_last_start.ne.0)'
          call mpi_abort(MPI_COMM_WORLD,1,ierr)
          stop
        endif
#else
        local_last_start=0
        total_local_size=iff3*iff2*iff1
        local_nlast=iff3
        iff3_g=iff3
        ifft3_g=iff3
#endif
!        print*,'before fft'
!        do j1=1,100000
!        do j2=1,100000
!        if(j1*j2.eq.1.and.j2.eq.1) print*,j1
!        enddo
!        enddo
!        print*,'ok'
        allocate (cfft(iff1,iff2,iff3),fft(iff1,iff2,iff3,2))
        allocate (ust(iff1,iff2,iff3,2))
        cfft=zero
        fft=zero
        ust=zero
        write(6,'(":FFT    global FFT-parameters:",3i5," Factor:",f5.2)') &
          iff1,iff2,iff3_g,faciff
        write(6,'(":FFT    local  FFT-parameters:",3i5)') iff1,iff2,iff3
        write(6,'(":FFT    allocation of cfft,cfft1,ust,fft:",f10.1," MB")') &
          iff1*iff2*iff3*6.*16/1000000
        IF(myid.EQ.0)  write(21,'(":IFFT  : FFT-parameters:",3i5," Factor:",f5.2)') &
          iff1,iff2,iff3_g,faciff
!        print*,'after fft'
!        do j1=1,100000
!        do j2=1,100000
!        if(j1*j2.eq.1.and.j2.eq.1) print*,j1
!        enddo
!        enddo

      endif

!Check for task parallelism
      IF(.NOT.coul) GOTO 777

!
! calculate multipolmoments of true charge density inside spheres from CLM
!
      DO jatom=nstart(myid),nstop(myid)
        DO J=1,JRI(JATOM)
          R(J)=R0(JATOM)*EXP(DX(JATOM)*dble(J-1))
        ENDDO
        DO LM1=1,LMMAX(JATOM)
          L=IABS( LM(1,LM1,JATOM) )
          DO J=1,JRI(JATOM)
            VALUE(J)=CLM(J,LM1,JATOM)*R(J)**L
          ENDDO
          CALL CHARGE(R,DX(JATOM),VALUE,1,JRI(JATOM),QIN)
          QQ(LM1,JATOM)=-QIN
          IF(L.EQ.0) THEN
            QEL(JATOM) =QIN
            QQ(1,JATOM)=QQ(1,JATOM)/SQFP + ZZ(JATOM)/SQFP 
! for charged cells
            delcc=supn/vol
            QEL(JATOM)=QEL(JATOM)-delcc*vsph(jatom)
          ENDIF
        ENDDO
      ENDDO
!
! write out multipolmoments
!
      WRITE(6,870)

      INDEX=1
      do jatom=1,NSTART(myid)-1
        INDEX=INDEX+MULT(JATOM)
      enddo

      DO JATOM=NSTART(myid),NSTOP(myid)

        WRITE(6,2040) JATOM,ANAME(JATOM),ZZ(JATOM),LMMAX(JATOM), &
          ( POS(J,INDEX), J=1,3 )
        INDEX=INDEX+MULT(JATOM)
        IF(IATNR(JATOM).GT.0) THEN
          LMMAX1(JATOM)=LMMAX(JATOM)
        ELSE
          IND=LMMAX(JATOM)
          DO 36 LM1=1,LMMAX(JATOM)
           IF(LM(2,LM1,JATOM).NE.0) THEN
             IND=IND+1
             QQ(IND,JATOM)=QQ(LM1,JATOM)
           endif
  36      CONTINUE
        endif
        CALL MULTFC(FC,JATOM,LMMAX1(JATOM),LM(1,1,JATOM),IATNR(JATOM))
        CALL MULTSU(QQ,fc,JATOM,LMMAX1(JATOM),IATNR(JATOM),lm)
        LMMULT1=0

        DO 38 LM1=1,LMMAX1(JATOM)
!
! combines qq of l,m and -l,m and defines new indices lmmult
!
          LMMULT1=LMMULT1+1
          QQ(LMMULT1,JATOM)=QQ(LM1,JATOM)
          LMMULT(1,LMMULT1,JATOM)=LM(1,LM1,JATOM)
          LMMULT(2,LMMULT1,JATOM)=LM(2,LM1,JATOM)
          IF (LMMULT1.EQ.1) GOTO 38
          IF ((IABS(LMMULT(1,LMMULT1,JATOM)) .EQ. &
            IABS(LMMULT(1,LMMULT1-1,JATOM))) .AND. &
            (LMMULT(2,LMMULT1,JATOM) .EQ. &
            LMMULT(2,LMMULT1-1,JATOM))) THEN
            QQ(LMMULT1-1,JATOM)=QQ(LMMULT1-1,JATOM)+QQ(LMMULT1,JATOM)
            LMMULT1=LMMULT1-1
          endif
  38    CONTINUE

        LMMTMX(JATOM)=LMMULT1

        WRITE(6,2050) LMMULT(1,1,JATOM), LMMULT(2,1,JATOM),QQ(1,JATOM)
        if(iprint.gt.0) then
           DO LMMULT1=2,LMMTMX(JATOM)
              WRITE(6,2050) LMMULT(1,LMMULT1,JATOM), &
                            LMMULT(2,LMMULT1,JATOM),QQ(LMMULT1,JATOM)
           ENDDO
        endif
        IF(IATNR(JATOM).GT.0) THEN
          DO LM1=1,LMMAX(JATOM)
            IF(lm(2,lm1,jatom).NE.0) THEN
              LMMAX1(JATOM)=LMMAX1(JATOM)+1
              IF(LMMAX1(JATOM).GT.NCOM+3) GOTO 970
              LM(1,LMMAX1(JATOM),JATOM)=lm(1,lm1,jatom)
              LMMULT(1,LMMAX1(JATOM),JATOM)=lm(1,lm1,jatom)
              LM(2,LMMAX1(JATOM),JATOM)=-lm(2,lm1,jatom)
              LMMULT(2,LMMAX1(JATOM),JATOM)=-lm(2,lm1,jatom)
            endif
          ENDDO
        ENDIF

      ENDDO
!
! calculate multipolmoments of plane waves in spheres
!
      INDEX=1
      do jatom=1,NSTART(myid)-1
        INDEX=INDEX+MULT(JATOM)
      enddo

      DO JATOM=NSTART(myid),NSTOP(myid)

        WRITE(6,2040) JATOM,ANAME(JATOM),ZZ(JATOM),LMMAX(JATOM), &
          ( POS(J,INDEX),J=1,3 )
        LLMM=LMMTMX(JATOM)
        ib_max=-999
        DO LM1=1,LLMM
!          CQ(LM1)=(0.d0,0.d0)
          ib_max=max(ib_max,IABS( lmmult(1,LM1,JATOM) ))
        ENDDO
!
! *** LDM CHANGES ***
        CWK(:) = (0.D0,0.D0)
        CQ(:) = (0.D0,0.D0)
! jdoumont 30/7/20 omp statements cause floating invalid during runtime???
! TODO: try allocating private arrays in parallel section!
!$omp parallel do &
!$omp& private(arg,yka,rhokt,lm1,ll,tc,bes) &
!$omp& reduction(+:cq,cwk)  schedule(guided)
        DO J=NKK,2,-1
         if(abs(rhok(j)) .gt. 1D-14)then
          ARG=RMT(JATOM)*ABSK(J)
          CALL SPHBES(ib_max+1,ARG,BES)
          CALL YKAV(J,JATOM,YKA,LMMTMX,LMMULT,LMMAX(JATOM))
          RHOKT=RHOK(J)/ABSK(J)
          DO LM1=1,LLMM
            LL=IABS( lmmult(1,LM1,JATOM) )
            TC=RHOKT*YKA(LM1)*BES(LL+1)
! ** Kahan Summation
!            Yw = TC - CWK(LM1)
!            TW = CQ(LM1) + YW
!            CWK(LM1) = (TW-CQ(LM1)) - YW
!            CQ(LM1) = TW
            call KSumC(CQ(LM1),TC,CWK(LM1))
! normal summation
!            CQ(lm1)= TC + cq(lm1) 
          ENDDO
         endif
        ENDDO
        CQ(1:LLMM) = CQ(1:LLMM) - 0.5D0*CWK(1:LLMM)
! *** END OF LDM CHANGES ***

        DO LM1=1,LLMM
          LL=IABS( LMMULT(1,LM1,JATOM) )
          QOUT=PI4*RMT(JATOM)**(LL+2)
          IF(LL.EQ.0) CQ(LM1)=CQ(LM1)+RHOK(1)*RMT(JATOM)/3.d0/SQFP
          CQ(LM1)=CQ(LM1)*QOUT*IMAG**LL
          QQNEW(LM1,JATOM)=QQ(LM1,JATOM) + CQ(LM1)
          if(iprint.gt.0) then
             WRITE(6,996)LMMULT(1,LM1,JATOM),LMMULT(2,LM1,JATOM),CQ(LM1)
          else if(lm1.eq.1) then
             WRITE(6,996)LMMULT(1,LM1,JATOM),LMMULT(2,LM1,JATOM),CQ(LM1)
          endif
        ENDDO
        INDEX=INDEX+MULT(JATOM)

      ENDDO

!
! calculate pseudo multipolmoments and define QQs for M<0
!
      INDEX=1
      do jatom=1,NSTART(myid)-1
        INDEX=INDEX+MULT(JATOM)
      enddo

      DO JATOM=NSTART(myid),NSTOP(myid)

        LLMM=LMMTMX(JATOM)
        IF(IATNR(JATOM).GT.0) THEN
          DO LM1=1,LMMAX(JATOM)
            IF(lm(2,lm1,jatom).NE.0) THEN
              LLMM=LLMM+1
              IF(llmm.GT.NCOM+3) GOTO 970
              QQNEW(llmm,JATOM)=dCONJG( QQNEW(lm1,JATOM) )
            endif
          ENDDO
        endif
        WRITE(6,2040) JATOM,ANAME(JATOM),ZZ(JATOM),LMMAX(JATOM), &
          (POS(J,INDEX),J=1,3)
        WRITE(6,994) ( LMMULT(J,1,JATOM),J=1,2 ),QQNEW(1,JATOM)
        if(iprint.gt.0) then
           DO LMX=2,LLMM
              WRITE(6,994) ( LMMULT(J,LMX,JATOM),J=1,2 ),QQNEW(LMX,JATOM)
           ENDDO
        endif
        INDEX=INDEX+MULT(JATOM)

      ENDDO

      CALL CPUTIM(TIME1)
      CALL wallTIM(TIME1w)

      IF(SWITCH.EQ.'MULT') then
          call abort_parallel
          STOP 'MULT OK'
      endif
      IF(SWITCH.EQ.'EXCH') GOTO 445
      RK=0.0
      DO JATOM=1,NAT
        RK1=RMT(JATOM)*ABSK(NKK)
        RK=MAX(RK,RK1)
      ENDDO

1212  WRITE(6,880)
      WRITE(6,1520) NCON,RK
      if(myid.eq.0) WRITE(21,1520) NCON,RK

! Buffer for parallel
      allocate(buf1(nkk),potk(nkk))
      buf1(1:nkk)=zeroc
      potk(1:nkk)=zeroc
!
! calculate plane wave Coulomb potential
!
! Find some parameters for the next loop

      DO JATOM=NSTART(myid),NSTOP(myid)
        IB_MAX=-999
        L_MAX=-999
        L_MIN=999
        LLMM=LMMTMX(JATOM)
        IF(IATNR(JATOM).GT.0) LLMM=LMMAX1(JATOM)
        DO LM1=1,LLMM
          LL=IABS( LMMULT(1,LM1,JATOM) )
          NCON1=NCON-LL
          IBES=LL+NCON1+1
          IF(IBES.GT.20) GOTO 901
          IF(IBES.GT.IB_MAX)IB_MAX=IBES
          IF(LL.GT.L_MAX)L_MAX=LL
          IF(LL.LT.L_MIN)L_MIN=LL
        ENDDO
      ENDDO

      L_MIN=MIN(0,L_MIN)
      L_MAX=max(L_MAX+2,20)
      allocate(IMAGA(L_MIN:L_MAX,NAT))

      DO JATOM=NSTART(myid),NSTOP(myid)
        DO LL=L_MIN,L_MAX
          NCON1=NCON-LL
          IBES=LL+NCON1+1
          MP=DFAK(IBES+1)/DFAK(LL+1)
          IMAGA(LL,JATOM)=((-IMAG/RMT(JATOM))**LL)*MP*MULT(JATOM)
        ENDDO
      ENDDO
      
      !This had to initialise here, otherwise EST in the interstitial become 2 time bigger.
     IF(SW_STRESS .eq. 'YES') THEN 
      SGGP     = 0.0  
      SGGP_GGA = 0.0
     ENDIF 
!*****************************************************             
! ! ! $omp parallel do private(arg,llmm,ll,ncon1,ibes,cp,jj,bes,jatom,lm1,yka) &
! ! ! $omp& schedule(guided)
!$omp parallel do private(arg,llmm,ll,ncon1,ibes,cp,jj,bes,jatom,lm1,yka) &
!$omp& schedule(guided)
      DO J=2,NKK
        DO JATOM=NSTART(myid),NSTOP(myid)

!!!#if defined(Parallel) || defined(_OPENMP)
!!!          ARG=RMT(JATOM)*ABSK(J)
!!!          CALL SPHBES(IB_MAX,ARG,BES)
!!!#else
! Do not recalculate bessel functions unless we have to
          JJ=max(1,JATOM-1)
!!          if( (JATOM.EQ.1) .OR. (ABS(RMT(JJ)-RMT(JATOM)).GT.1D-6))then
          if( (JATOM.EQ.NSTART(myid)) .OR. (ABS(RMT(JJ)-RMT(JATOM)).GT.1D-6))then
            ARG=RMT(JATOM)*ABSK(J)
            CALL SPHBES(IB_MAX,ARG,BES)
          endif
!!!#endif

          LLMM=LMMTMX(JATOM)
          IF(IATNR(JATOM).GT.0) LLMM=LMMAX1(JATOM)
          CALL YKAV(J,JATOM,YKA,LMMTMX,LMMULT,LLMM)
          DO LM1=1,LLMM
            LL=IABS( LMMULT(1,LM1,JATOM) )
            NCON1=NCON-LL
            IBES=LL+NCON1+1
            CP=  IMAGA(LL,JATOM)/(ARG**(NCON1+1))* &
              BES(IBES)*QQNEW(LM1,JATOM)*dCONJG(YKA(LM1))
            buf1(J)=buf1(J)+CP
          ENDDO
          
        ENDDO   !jatom loop
      ENDDO
     
#ifdef Parallel
      call MPI_ALLREDUCE(buf1,potk,nkk,MPI_DOUBLE_COMPLEX,MPI_SUM,comm_coul,ierr)
#else
! writing pseudo-densities into a file for tests
! write(888,2071) 0,0,0,(0.d0,0.d0)
! write(889,2071) 0,0,0,rhok(1)
      DO j=2, nkk
        potk(j)=buf1(j)
! write(888,2071)  ( KZZ(JX,j), JX=1,3 ), -potk(j)*inst(J)/vol*4.d0*pi
! write(889,2071)  ( KZZ(JX,j), JX=1,3 ), -potk(j)*inst(J)/vol*4.d0*pi+rhok(j)
      ENDDO
#endif
     
     IF(SW_STRESS .eq. 'YES') THEN
      call walltim(time1_st)
      DO J=2,NKK
          !========col stress int, this is here becoz SGGP comes from ykav which depends on jatom
          fac1_p = 2.d0*pi4/( ABSK(J)*ABSK(J) )
          index_p = 0
          DO alpha1 = 1,3
             DO beta1 = 1,3
                index_p = index_p + 1
                DECOUL_POT_PW(J,index_p) = &
                -SGGP(J,index_p)*( RHOK(J)/dble(INST(J))-PI4/VOL*potk(J) )*fac1_p*0.5d0!/NSTOP(myid) 
!                -SGGP(J,index_p)*( RHOK(J)/dble(INST(J))-PI4/VOL*buf1(J) )*fac1_p*0.5d0!/NSTOP(myid) 
                gga_int(j,index_p) = SGGP_GGA(j,index_p)*RHOK(J)/dble(INST(j))                
             ENDDO
          ENDDO           
      ENDDO
      call walltim(time2_st)
      time_electro = time_electro + (time1_st-time2_st)
     ENDIF      
          !****************************************************************
          !================stress block====================================
            

      TMP=4.D0*PI/VOL
      TMP1=-8.D0*PI
      DO j =2, nkk
        POTK(J)=(TMP*POTK(J)*dble(INST(J))-RHOK(J))*TMP1/(ABSK(J)*ABSK(J))
      ENDDO
! sum largest POTK and check convergence
      potktest=0.d0
      if(nkk.gt.10) then
        do j=nkk-9,nkk
          potktest=potktest+abs(potk(j))
        enddo
      endif
      if(myid.eq.0) then
        write(21,'(":VKCOUL :  VK-COUL convergence:",e10.3)') potktest/10.d0
        if(potktest/10.d0.gt.1.d-5) then
          write(21,'(":WARN     VK-COUL not well converged: Increase GMAX or decrease NCON")') 
        endif
      endif
      deallocate(buf1)
      CALL CPUTIM(TIMEy)
      CALL wallTIM(TIMEyw)
      t947=0.d0
      t947w=0.d0
!
! calculate Coulomb potential inside muffin-tin spheres
!
      allocate(fhf(0:3,1:nat),vcoul0(1:nat))
      fhf = zero ; vcoul0 = zero

      NKK1=NKK-1
      INDEX=1
      do jatom=1,NSTART(myid)-1
        INDEX=INDEX+MULT(JATOM)
      enddo
      rhopw0=0.d0
      lrhopw=.false.
      IF(SWITCH.EQ.'COUL') lrhopw=.true.
      DO 50 JATOM=NSTART(myid),NSTOP(myid)

        jrjstart=jri(jatom)
        jrjend=jri(jatom)
        IF(SWITCH.EQ.'COUL') jrjstart=jri(jatom)-3
        IF(SWITCH.EQ.'COUL') jrjend=jri(jatom)+3
        LLMM=LMMAX(JATOM)
        JRJ=JRI(JATOM)
        DO 54 J=1,JRJ
          R(J)=R0(JATOM)*EXP(dble(J-1)*DX(JATOM))
          Rhopw0r(J,jatom)=R(J)
  54    CONTINUE
        DO J=jrj+1,JRJ+10
          Rhopw0r(J,jatom)=R0(JATOM)*EXP(dble(J-1)*DX(JATOM))
        enddo
        CALL CPUTIM(ts947)
        CALL wallTIM(ts947w)
!
! Coulomb pot at R-Mt
!
       if(nint(zz(jatom)).gt.118)  stop 'atomic mass: Z gt 118'
        rnuc=2.2677D-05*(atom_mass(nint(zz(jatom)))**(1.d0/3.d0))
        write(6,*)'amass, r0:',atom_mass(nint(zz(jatom))),rnuc
        DO J=1,JRJ
          if (r(j).ge.rnuc) exit
        ENDDO
        nuc=J
        rnuc=r(j)
        write(6,*)'nuc:',nuc
        rhopw00=(0.d0,0.d0)
        ib_max=-999
        DO LM1=1,LLMM
          CVOUT(LM1)=(0.d0,0.d0)
          ib_max=max(ib_max,IABS( lmmult(1,LM1,JATOM) ))
 !         ib_max=max(ib_max,IABS(LM(1,LM1,JATOM)),IABS( lmmult(1,LM1,JATOM) )) ! LDM ?
        ENDDO

!        CrhoOUT=(0.d0,0.d0)
! *** LDM CHANGES ***
        CWK(:)=(0.D0,0.D0)
        CVOUT(:)=(0.D0,0.D0)
        RHOPW00(:,:)=(0.D0,0.D0)
! jdoumont 30/7/20 omp statements cause floating invalid/segmentation error (big case) during runtime???
!$omp parallel do reduction(+:rhopw00,cwk,cvout) &
!$omp& private(jrj1,arg,BESpw,bes,lm1,L,M,LMMULT1,YKA,TC) &
!$omp& schedule(guided)
        DO J=NKK,2,-1
         if(abs(Potk(j)).gt.1D-15)then
          do jrj1=jrjstart,jrjend
            arg=rhopw0r(jrj1,jatom)*absk(j)
            CALL SPHBES(ib_max+1,ARG,BESpw(0,jrj1))
          enddo
          bes(:)=bespw(:,jrj)
          CALL YKAV(J,JATOM,YKA,LMMTMX,LMMULT,LLMM)
          LMMULT1=0

          DO LM1=1,LLMM
            L=IABS( LM(1,LM1,JATOM) )
            M=LM(2,LM1,JATOM)
            LMMULT1=LMMULT1+1
            IF (LMMULT1.NE.1) THEN
              IF ((L.EQ.IABS(LMMULT(1,LMMULT1-1,JATOM))).AND. &
                (M.EQ.LMMULT(2,LMMULT1-1,JATOM))) THEN
                LMMULT1=LMMULT1-1
              endif
            endif
!            CVOUT(LM1)=CVOUT(LM1) + POTK(J)*BES(L)*YKA(LMMULT1)
! Kahan Summation
            TC = POTK(J)*BES(L)*YKA(LMMULT1)
 !           YW = TC - CWK(LM1)
 !           TW = CVOUT(LM1) + YW
 !           CWK(LM1)  = (TW - CVOUT(LM1)) - YW
 !           CVOUT(LM1) = TW
            call KSumC(CVOUT(LM1),TC,CWK(LM1))

!            CrhoOUT(LM1,1)=CrhoOUT(LM1,1) + rhoK(J)*BES(L)*YKA(LMMULT1)
!            if(lrhopw) call rhopw(rhok(j),yka(lmmult1),ABSK(J),rhopw0r(1,jatom), &
!                       jrjstart,jrj,rhopw00(1,lm1),l,jrjend-jrj)
            rhopw00(jrjstart:jrjend,lm1)=rhopw00(jrjstart:jrjend,lm1)+rhok(j)*yka(lmmult1)*bespw(l,jrjstart:jrjend)
          ENDDO
         endif
      ENDDO
        CVOUT(1:LLMM)=CVOUT(1:LLMM)-CWK(1:LLMM)*0.5D0
! *** END OF CHANGES

          CALL SPHBES(2,0.d0,BES)
          do j=jrjstart,jrjend
            rhopw0(j,1,jatom)=(rhopw00(j,1)+rhok(1)/sqfp*bes(0))*sqfp*4.d0*pi*rhopw0r(j,jatom)**2
          enddo
          do lm1=2,llmm
            CALL potfac (vinf,voutf,lm(1,lm1,jatom), &
              lm(2,lm1,jatom),iatnr(jatom))
            M=LM(2,LM1,JATOM)
            IMAG1=(1.d0,0.d0)
            IF(M.NE.0) THEN
              IF(LM(1,LM1,JATOM).LT.0) IMAG1=-IMAG
              IF(MOD(M,2).EQ.1) IMAG1=-IMAG1
              imag1=1.d0/IMAG1*SQRT2
            ENDIF
            rhopw0(jrjstart:jrjend,lm1,jatom)=rhopw00(jrjstart:jrjend,lm1)*voutf*imag1* &
               IMAG**IABS(LM(1,LM1,JATOM))*4.d0*pi*rhopw0r(jrjstart:jrjend,jatom)**2
         enddo

        IF(SWITCH.EQ.'COUL') then
          write(3,'("ATOM",i4,/,(50f15.8))') jatom
          do j=1,jrjend
             write(3,'((50f15.8))')  rhopw0r(j,jatom),&
               (rhopw0(j,lm1,jatom),lm1=1,llmm)
          enddo
          do j=1,jrj
             write(333,'((50f15.8))')  rhopw0r(j,jatom),&
               (clm(j,lm1,jatom),lm1=1,llmm)
          enddo
        endif

        write(6,*) 'TEST of density continuity (clm,rhopw) for atom',jatom
        DO LM1=1,LLMM
          write(6,'(i4,2e15.5)') lm1,clm(jrj,lm1,jatom),rhopw0(jrj,lm1,jatom)
          M=LM(2,LM1,JATOM)
          IF(M.NE.0) THEN
            IMAG1=(1.d0,0.d0)
            IF(LM(1,LM1,JATOM).LT.0) IMAG1=-IMAG
            IF(MOD(M,2).EQ.1) IMAG1=-IMAG1
            CVOUT(lm1)=CVOUT(lm1)/IMAG1*SQRT2
          ENDIF
        ENDDO
!!!        endif
!end of rhopw test

        CALL CPUTIM(te947)
        CALL wallTIM(te947w)
        t947=t947+(te947-ts947)
        t947w=t947w+(te947w-ts947w)

        if (switch1 .eq. ' EECE') then

          if (jeece(jatom) .eq. 1) then

            do ispin=1,2
              DO J=1,JRJ
                VALUE(J)=CLMSP(J,1,JATOM,ISPIN)
              ENDDO
              CALL CHARGE (R,DX(JATOM),VALUE,1,jrj,VLIN1) ! Was charg2
!              ACORR=VALUE(1)*R(1)*0.5D0
              elnumsp(jatom,ispin)=vlin1  ! +acorr
            enddo
            elnum(jatom)=elnumsp(jatom,1)+elnumsp(jatom,2)
            write(6,553) jatom,(elnumsp(jatom,ispin),ispin=1,2),elnum(jatom)
 553        format(i4,3f12.6,' atom, in sphere electron numbers: up, down, total')

            if(CoulEM)then
             fk(:,jatom)=0.D0
             if(elnum(jatom) .gt. 1D-8)then
              do L=0,6,2
                DO J=1,JRJ
                  vaLue(J) = cLm(J,1,jatom)*r(J)**L     !/r(irad)**(L+1)
                  WLDM (J) = cLm(J,1,jatom)/r(J)**(L+1) ! *r(irad)**L
                enddo
!
                call CoulInt(R,DX(JATOM),JRJ,VALUE)
                call CoulInt(R,DX(JATOM),JRJ,WLDM )
!
                do irad=1,jrj
                  VLin1 = Value(IRAD)/r(irad)**(L+1)
                  VLin2 = (WLDM(JRJ)-WLDM(IRAD))*r(irad)**L         

                  RLDM1(irad)=(VLIN1+VLIN2)*cLm(irad,1,jatom)
                enddo
                CALL CHARGE (R,DX(JATOM),RLDM1,1,jrj,VLIN1) 
                fk(L,jatom)=vLin1/eLnum(jatom)**2

                write(6,551)jatom,L,fk(L,jatom)*2.d0
              enddo
             endif
            else
             do L=0,6,2
              do irad=1,jrj
                do j=1,irad-1
                  vaLue(j)=cLm(j,1,jatom)*r(j)**L/r(irad)**(L+1)
                enddo
                do j=irad,jrj
                  vaLue(j)=cLm(j,1,jatom)*r(irad)**L/r(j)**(L+1)
                enddo

                CALL CHARG2 (R,DX(JATOM),VALUE,1,IRAD,VLIN1)
                CALL CHARG3 (R,DX(JATOM),VALUE,IRAD,JRJ,VLIN2)
!                ACORR=VALUE(1)*R(1)*0.5D0
                RLDM1(irad)=VLIN1+VLIN2 !+ACORR

              enddo

              if (eLnum(jatom) .gt. 1d-8) then
                DO J=1,JRJ
                  VALUE(J)=RLDM1(J)*cLm(j,1,jatom)/eLnum(jatom)**2
                ENDDO
              else
                value = 0d0
              endif

              CALL CHARGE (R,DX(JATOM),VALUE,1,jrj,VLIN1)      ! Was charg2
!              ACORR=VALUE(1)*R(1)*0.5D0
              fk(L,jatom)=vLin1 ! +acorr
            
              write(6,551)jatom,L,fk(L,jatom)*2.d0
 551          format(2i4,f12.6,'  atom, L, SLater integraL F^2L, Ry')
             enddo
            endif

            DO LM1=1,LLMM

              CALL potfac (vinf,voutf,Lm(1,Lm1,jatom),Lm(2,Lm1,jatom),iatnr(jatom))
              L=IABS( LM(1,LM1,JATOM) )
              M=LM(2,LM1,JATOM)
              LMOD=0
              IF(L.EQ.0) LMOD=1

! EECE caLcuLation of Vc in the sphere using 'in sphere rho' onLy
              if(CoulEM)then
                 do J=1,JRJ
                   vaLue(J)=cLm(J,1,jatom)*r(j)**L !/r(irad)**(L+1)
                   WLDM(J) =cLm(J,1,jatom)/r(j)**(L+1) ! *r(irad)**L
                 enddo
!
                 call CoulInt(R,DX(JATOM),JRJ,VALUE)
                 call CoulInt(R,DX(JATOM),JRJ,WLDM )
              endif
              DO IRAD=1,JRJ
                if(CoulEM)then
                   VLin1 = Value(IRAD)/r(irad)**(L+1)
                   VLin2 = (WLDM(JRJ)-WLDM(IRAD))*r(irad)**L     
                else
                  do j=1,irad-1
                    vaLue(j)=cLm(j,Lm1,jatom)*r(j)**L/r(irad)**(L+1)
                  enddo
                  do j=irad,jrj
                    value(j)=clm(j,lm1,jatom)*r(irad)**l/r(j)**(l+1)
                  enddo

                  CALL CHARG2 (R,DX(JATOM),VALUE,1,IRAD,VLIN1)
                  CALL CHARG3 (R,DX(JATOM),VALUE,IRAD,JRJ,VLIN2)
!                 ACORR=VALUE(1)*R(1)*0.5D0
                endif
                 VC=VLIN1+VLIN2   !+ACORR
! EECE: the coefficient is the same as for VOLINT
! for 00 comp. vcoef=2*sqrt(4*pi), coef 2 because v is in Ry
! sqrt(4*pi) because integral over angles (4pi) is multiplied
! by norm of clm_00 (1/sqrt(4pi))
                vcoef =8.0d0*PI/(2*L+1)*VINF
                v(irad,lm1,jatom) =VC*vcoef
              enddo

            ENDDO

          endif

        else !EECE end, original LAPW0 part follows

!
!     v-coul inside spheres
!
!$omp  parallel do &
!$omp& private(L,M,LMOD,J,RLDM1,rldm2,WLDM,IRAD,rj1,rj2,volint,surfin,zshift, &
!$omp&         value,acorr,vlin1,vlin2,vinf,voutf) schedule(guided)
          DO LM1=1,LLMM

            CALL potfac (vinf,voutf,lm(1,lm1,jatom), &
              lm(2,lm1,jatom),iatnr(jatom))
            L=IABS( LM(1,LM1,JATOM) )
            M=LM(2,LM1,JATOM)
            LMOD=0
            IF(L.EQ.0) LMOD=1
            DO J=1,JRJ
              RLDM1(J)=R(J)**L
            ENDDO
            DO J=1,JRJ
              RLDM2(J)=( 1.0D0-(R(J)/RMT(JATOM))**(2*L+1))/R(J)**(L+1)
            ENDDO

! *** LDM CHANGE HERE
            if( (L.eq.0) .and. (abs(delcc).gt.1D-12) )then
                VALUE(1:JRJ)=(CLM(1:JRJ,LM1,JATOM)-delcc*r(1:JRJ)*r(1:JRJ)*pi4 ) !*RLDM1(1:JRJ)
                WLDM (1:JRJ)=(CLM(1:JRJ,LM1,JATOM)-delcc*r(1:JRJ)*r(1:JRJ)*pi4 )*RLDM2(1:JRJ)
            else if(L.eq.0)then
                VALUE(1:JRJ)=CLM(1:JRJ,LM1,JATOM)
                WLDM (1:JRJ)=CLM(1:JRJ,LM1,JATOM)*RLDM2(1:JRJ)
            else
                VALUE(1:JRJ)=CLM(1:JRJ,LM1,JATOM)*RLDM1(1:JRJ)
                WLDM (1:JRJ)=CLM(1:JRJ,LM1,JATOM)*RLDM2(1:JRJ)
            endif
!
        if(CoulEM)then
!               Calculate terms for Euler-Maclaurian
!               Coulomb integral
                call CoulInt(R,DX(JATOM),JRJ,VALUE)
                call CoulInt(R,DX(JATOM),JRJ,WLDM )
        endif
!       Revision to this loop, LDM, April 2014
            DO 53 IRAD=1,JRJ
              rj1=rldm2(irad)
              rj2=RLDM1(IRAD)
!
!     integration in two parts, simpson + 3/8 formula at start(end)
!     for even points
!
!              CALL CHARG2 (R,DX(JATOM),VALUE,1,IRAD,VLIN1)
!              CALL CHARG3 (R,DX(JATOM),VALUE,IRAD,JRJ,VLIN2)
!              ACORR=VALUE(1)*R(1)*0.5D0
              if(CoulEM)then
                VLin1 = Value(IRAD)*RJ1
                VLin2 = (WLDM(JRJ)-WLDM(IRAD))*RJ2
              else
                CALL CHARGE (R,DX(JATOM),VALUE,1,IRAD,VLIN1)    ; VLIN1=VLIN1*RJ1
                CALL CHARGE (R,DX(JATOM),WLDM ,IRAD,JRJ,VLIN2)  ; VLIN2=VLIN2*RJ2
              endif
              VOLINT=VLIN1+VLIN2 ! +ACORR
              VOLINT =VOLINT*8.0*PI/dble(2*L+1)*VINF
              SURFIN=CVOUT(LM1)*(-4.D0*PI)*(R(IRAD)/RMT(JATOM))**L* &
                IMAG**L*VOUTF
              ZSHIFT=(-2.D0*ZZ(JATOM)/R(IRAD)+ &
                2.D0*ZZ(JATOM)/RMT(JATOM))*SQFP
              V(IRAD,LM1,JATOM)=(VOLINT-SURFIN+LMOD*ZSHIFT)
!if(l.eq.0.and.irad.gt.jrj-5) print*, volint,surfin,zshift*lmod
!if(l.eq.0.and.irad.lt.20) print*, r(irad),V(IRAD,LM1,JATOM)*r(irad)
              IF(IRAD.EQ.nuc) THEN

                IF(L.EQ.2 .AND. M.EQ.0) THEN
                  efg_V20SRF(jatom)=SURFIN
                  efg_V20(jatom)=-V(nuc,LM1,JATOM)
                ELSE IF(LM(1,LM1,JATOM).EQ. 2.AND. M.EQ.2) THEN
                  efg_V22SRF(jatom)=SURFIN
                  efg_V22(jatom)=-V(nuc,LM1,JATOM)
                ELSE IF(LM(1,LM1,JATOM).EQ.-2.AND. M.EQ.2) THEN
                  efg_V22MSR(jatom)=SURFIN
                  efg_V22M(jatom)=-V(nuc,LM1,JATOM)
                ELSE IF(LM(1,LM1,JATOM).EQ.2.AND. M.EQ.1) THEN
                  efg_V21SRF(jatom)=SURFIN
                  efg_V21(jatom) =-V(nuc,LM1,JATOM)
                ELSE IF(LM(1,LM1,JATOM).EQ.-2.AND. M.EQ.1) THEN
                  efg_V21MSR(jatom)=SURFIN
                  efg_V21M(jatom)=-V(nuc,LM1,JATOM)
                ENDIF

              ENDIF

              IF(LM1.EQ.1) VCOUL2(IRAD,JATOM)=(VOLINT-SURFIN)/2.D0 + &
                LMOD*(ZZ(JATOM)/RMT(JATOM))*SQFP

  53        CONTINUE

! calculation of ik. component of Hellman-Feynman force
! from:  lim(r-->0) V_1m / r * grad( r * y_lmp)
!
! component of force: y
! ll,m --> ik
! +1,1     1
! -1,1     2
! +1,0     3

            IF (l.EQ.1) THEN
              force_buf=.TRUE.
              ll=lm(1,lm1,jatom)
              ik=MOD(ll+m+1,3)+1
              fhf(ik,jatom)=v(1,lm1,jatom)/r(1)*zz(jatom)*SQRT(3.d0/4.d0/pi)
            ENDIF

          ENDDO

          vcoul0(jatom)=v(jri(jatom),1,jatom)/sqfp
          WRITE(6,1335)  jatom,jatom,aname(jatom),vcoul0(jatom)

          IF(lmmax(jatom).GT.1.AND.iatnr(jatom).LT.0) THEN
            CALL EFG(efg_v20(jatom),efg_v21(jatom),efg_v22(jatom) &
              ,efg_v21m(jatom),efg_v22m(jatom) &
              ,efg_v20srf(jatom),efg_v21srf(jatom),efg_v22srf(jatom) &
              ,efg_v21msr(jatom),efg_v22msr(jatom),efg_vzz(jatom) &
              ,efg_qmat(1,1,jatom),efg_tens(1,jatom),efg_eivec(1,1,jatom) &
              ,efg_winkel(jatom),efg_eta(jatom),rnuc,jatom )
            logefg(jatom)=.TRUE.
          ENDIF

        endif !EECE the original part remained

  50  CONTINUE

!
! Serial output of coulomb potential
!
      IF (LR2V) THEN

#ifdef Parallel            
        IF(myid.NE.0) THEN
          CALL MPI_RECV(i,1,MPI_INTEGER,myid-1,MPI_ANY_TAG,comm_coul &
            ,statusmpi,ierr)
        ENDIF
#endif
        
        IF(myid.EQ.0) THEN
          OPEN(10,FILE=fnavn(10))
          WRITE(10,*) '    COULOMB POTENTIAL'
          WRITE(10,*) '    NORM OF VC(L,M)=VC(L,M)*r**2'
          WRITE(10,*)
        else
          OPEN(10,FILE=fnavn(10),POSITION='append')
        ENDIF
        
        DO jatom=nstart(myid),nstop(myid)
          WRITE(10,1990) jatom
          WRITE(10,2001) lmmax(jatom)
          DO j=1,jri(jatom)
            r(j)=r0(jatom)*EXP(dx(jatom)*(j-1))
          ENDDO
          DO lm1=1,lmmax(jatom)
            WRITE(10,2011) lm(1,lm1,jatom),lm(2,lm1,jatom)
            WRITE(10,2021) (v(j,lm1,jatom)*r(j)*r(j),j=1,jri(jatom))
            WRITE(10,2031)
          ENDDO
          WRITE(10,2030)
        ENDDO

        IF(myid.EQ.icoul) THEN
          WRITE(10,*) '    COULOMB POTENTIAL IN INTERSTITIAL'
          WRITE(10,*)
          WRITE(10,2061) nkk
          if(inversion)then
            DO j=1,nkk
              WRITE(10,2071) (kzz(jx,j),jx=1,3), dble(potk(j)),0.D0
            ENDDO
          else
            DO j=1,nkk
              WRITE(10,2071) (kzz(jx,j),jx=1,3), potk(j)
            ENDDO
          endif
        ENDIF

        CLOSE(10)

#ifdef Parallel
        IF(myid.LT.icoul) THEN
          CALL MPI_SEND(i,1,MPI_INTEGER,myid+1, 0,comm_coul,ierr)
        ENDIF
#endif

      ENDIF

#ifdef Parallel
!!$! Reduce efgb for later output
!!$      efgb_typelist(1:16) = MPI_DOUBLE_PRECISION
!!$      efgb_blocklength(1:13) = 1
!!$      efgb_blocklength(14) = 9
!!$      efgb_blocklength(15) = 3
!!$      efgb_blocklength(16) = 9
!!$!MPI1_
!!$      CALL MPI_ADDRESS(efgb(1)%v20,   efgb_address(2),  ierr)
!!$      CALL MPI_ADDRESS(efgb(1)%v21,   efgb_address(3),  ierr)
!!$      CALL MPI_ADDRESS(efgb(1)%v22,   efgb_address(4),  ierr)
!!$      CALL MPI_ADDRESS(efgb(1)%v21m,  efgb_address(5),  ierr)
!!$      CALL MPI_ADDRESS(efgb(1)%v22m,  efgb_address(6),  ierr)
!!$      CALL MPI_ADDRESS(efgb(1)%v20srf,efgb_address(7),  ierr)
!!$      CALL MPI_ADDRESS(efgb(1)%v21srf,efgb_address(8),  ierr)
!!$      CALL MPI_ADDRESS(efgb(1)%v22srf,efgb_address(9),  ierr)
!!$      CALL MPI_ADDRESS(efgb(1)%v21msr,efgb_address(10), ierr)
!!$      CALL MPI_ADDRESS(efgb(1)%v22msr,efgb_address(11), ierr)
!!$      CALL MPI_ADDRESS(efgb(1)%vzz,   efgb_address(12), ierr)
!!$      CALL MPI_ADDRESS(efgb(1)%winkel,efgb_address(13), ierr)
!!$      CALL MPI_ADDRESS(efgb(1)%eta,   efgb_address(14), ierr)
!!$      CALL MPI_ADDRESS(efgb(1)%qmat,  efgb_address(15), ierr)
!!$      CALL MPI_ADDRESS(efgb(1)%tens,  efgb_address(16), ierr)
!!$      CALL MPI_ADDRESS(efgb(1)%eivec, efgb_address(17), ierr)
!!$!MPI-2
!!$!     CALL MPI_GET_ADDRESS(efgb(1)%v20,   efgb_address(2),  ierr)
!!$!     CALL MPI_GET_ADDRESS(efgb(1)%v21,   efgb_address(3),  ierr)
!!$!     CALL MPI_GET_ADDRESS(efgb(1)%v22,   efgb_address(4),  ierr)
!!$!     CALL MPI_GET_ADDRESS(efgb(1)%v21m,  efgb_address(5),  ierr)
!!$!     CALL MPI_GET_ADDRESS(efgb(1)%v22m,  efgb_address(6),  ierr)
!!$!     CALL MPI_GET_ADDRESS(efgb(1)%v20srf,efgb_address(7),  ierr)
!!$!     CALL MPI_GET_ADDRESS(efgb(1)%v21srf,efgb_address(8),  ierr)
!!$!     CALL MPI_GET_ADDRESS(efgb(1)%v22srf,efgb_address(9),  ierr)
!!$!     CALL MPI_GET_ADDRESS(efgb(1)%v21msr,efgb_address(10), ierr)
!!$!     CALL MPI_GET_ADDRESS(efgb(1)%v22msr,efgb_address(11), ierr)
!!$!     CALL MPI_GET_ADDRESS(efgb(1)%vzz,   efgb_address(12), ierr)
!!$!     CALL MPI_GET_ADDRESS(efgb(1)%winkel,efgb_address(13), ierr)
!!$!     CALL MPI_GET_ADDRESS(efgb(1)%eta,   efgb_address(14), ierr)
!!$!     CALL MPI_GET_ADDRESS(efgb(1)%qmat,  efgb_address(15), ierr)
!!$!     CALL MPI_GET_ADDRESS(efgb(1)%tens,  efgb_address(16), ierr)
!!$!     CALL MPI_GET_ADDRESS(efgb(1)%eivec, efgb_address(17), ierr)
!!$     
!!$      DO i=1,16
!!$        efgb_displace(i)=efgb_address(i+1)-efgb_address(2)
!!$      ENDDO
!!$
!!$      CALL MPI_TYPE_STRUCT(16, efgb_blocklength, efgb_displace, efgb_typelist &
!!$        , typempi_efgb, ierr)
!!$
!!$      CALL MPI_TYPE_COMMIT(typempi_efgb, ierr)
!!$
!!$      IF(myid.NE.0) THEN
!!$        DO j=nstart(myid),nstop(myid)
!!$          CALL MPI_SEND(efgb(j)%v20, 1, typempi_efgb, 0, 0, comm_coul, ierr)
!!$          CALL MPI_SEND(vcoul0(j), 1, MPI_DOUBLE_PRECISION, 0, 0, comm_coul, ierr)
!!$          CALL MPI_SEND(fhf(0,j), 4, MPI_DOUBLE_PRECISION, 0, 0, comm_coul, ierr)
!!$          CALL MPI_SEND(logefg(j), 1, MPI_LOGICAL, 0, 0, comm_coul, ierr) 
!!$          CALL MPI_SEND(fk(0,j), 7, MPI_DOUBLE_PRECISION, 0, 0, comm_coul, ierr)
!!$          CALL MPI_SEND(elnum(j), 1, MPI_DOUBLE_PRECISION, 0, 0, comm_coul, ierr)
!!$        ENDDO
!!$      ELSE
!!$        DO i=1,icoul
!!$          DO j=nstart(i),nstop(i)
!!$            CALL MPI_RECV(efgb(j)%v20, 1, typempi_efgb, i, MPI_ANY_TAG, comm_coul, statusmpi, ierr)
!!$            CALL MPI_RECV(vcoul0(j), 1, MPI_DOUBLE_PRECISION, i, MPI_ANY_TAG, comm_coul, statusmpi, ierr)
!!$            CALL MPI_RECV(fhf(0,j), 4, MPI_DOUBLE_PRECISION, i, MPI_ANY_TAG, comm_coul, statusmpi, ierr)
!!$            CALL MPI_RECV(logefg(j), 1, MPI_LOGICAL, i, MPI_ANY_TAG, comm_coul, statusmpi, ierr)
!!$            CALL MPI_RECV(fk(0,j), 7, MPI_DOUBLE_PRECISION, i, MPI_ANY_TAG, comm_coul, statusmpi, ierr)
!!$            CALL MPI_RECV(elnum(j), 1, MPI_DOUBLE_PRECISION, i, MPI_ANY_TAG, comm_coul, statusmpi, ierr)
!!$          ENDDO
!!$        ENDDO
!!$      ENDIF

      IF(myid.NE.0) THEN
        DO j=nstart(myid),nstop(myid)
          CALL MPI_SEND(efg_v20(j)   ,    1, MPI_DOUBLE_PRECISION, 0, 0, comm_coul, ierr)
          CALL MPI_SEND(efg_v21(j)   ,    1, MPI_DOUBLE_PRECISION, 0, 0, comm_coul, ierr)
          CALL MPI_SEND(efg_v22(j)   ,    1, MPI_DOUBLE_PRECISION, 0, 0, comm_coul, ierr)
          CALL MPI_SEND(efg_v21m(j)  ,    1, MPI_DOUBLE_PRECISION, 0, 0, comm_coul, ierr)
          CALL MPI_SEND(efg_v22m(j)  ,    1, MPI_DOUBLE_PRECISION, 0, 0, comm_coul, ierr)
          CALL MPI_SEND(efg_v20srf(j),    1, MPI_DOUBLE_PRECISION, 0, 0, comm_coul, ierr)
          CALL MPI_SEND(efg_v21srf(j),    1, MPI_DOUBLE_PRECISION, 0, 0, comm_coul, ierr)
          CALL MPI_SEND(efg_v22srf(j),    1, MPI_DOUBLE_PRECISION, 0, 0, comm_coul, ierr)
          CALL MPI_SEND(efg_v21msr(j),    1, MPI_DOUBLE_PRECISION, 0, 0, comm_coul, ierr)
          CALL MPI_SEND(efg_v22msr(j),    1, MPI_DOUBLE_PRECISION, 0, 0, comm_coul, ierr)
          CALL MPI_SEND(efg_vzz(j),       1, MPI_DOUBLE_PRECISION, 0, 0, comm_coul, ierr)
          CALL MPI_SEND(efg_winkel(j),    1, MPI_DOUBLE_PRECISION, 0, 0, comm_coul, ierr)
          CALL MPI_SEND(efg_eta(j),       1, MPI_DOUBLE_PRECISION, 0, 0, comm_coul, ierr)
          CALL MPI_SEND(efg_qmat(1,1,j),  9, MPI_DOUBLE_PRECISION, 0, 0, comm_coul, ierr)
          CALL MPI_SEND(efg_tens(1,j),    3, MPI_DOUBLE_PRECISION, 0, 0, comm_coul, ierr)
          CALL MPI_SEND(efg_eivec(1,1,j), 9, MPI_DOUBLE_PRECISION, 0, 0, comm_coul, ierr)

          CALL MPI_SEND(vcoul0(j), 1, MPI_DOUBLE_PRECISION, 0, 0, comm_coul, ierr)
          CALL MPI_SEND(fhf(0,j), 4, MPI_DOUBLE_PRECISION, 0, 0, comm_coul, ierr)
          CALL MPI_SEND(logefg(j), 1, MPI_LOGICAL, 0, 0, comm_coul, ierr) 
          CALL MPI_SEND(fk(0,j), 7, MPI_DOUBLE_PRECISION, 0, 0, comm_coul, ierr)
          CALL MPI_SEND(elnum(j), 1, MPI_DOUBLE_PRECISION, 0, 0, comm_coul, ierr)
        ENDDO
      ELSE
        DO i=1,icoul
          DO j=nstart(i),nstop(i)
            CALL MPI_RECV(efg_v20(j)   ,     1, MPI_DOUBLE_PRECISION, i, MPI_ANY_TAG, comm_coul, statusmpi, ierr)
            CALL MPI_RECV(efg_v21(j)   ,     1, MPI_DOUBLE_PRECISION, i, MPI_ANY_TAG, comm_coul, statusmpi, ierr)
            CALL MPI_RECV(efg_v22(j)   ,     1, MPI_DOUBLE_PRECISION, i, MPI_ANY_TAG, comm_coul, statusmpi, ierr)
            CALL MPI_RECV(efg_v21m(j)  ,     1, MPI_DOUBLE_PRECISION, i, MPI_ANY_TAG, comm_coul, statusmpi, ierr)
            CALL MPI_RECV(efg_v22m(j)  ,     1, MPI_DOUBLE_PRECISION, i, MPI_ANY_TAG, comm_coul, statusmpi, ierr)
            CALL MPI_RECV(efg_v20srf(j),     1, MPI_DOUBLE_PRECISION, i, MPI_ANY_TAG, comm_coul, statusmpi, ierr)
            CALL MPI_RECV(efg_v21srf(j),     1, MPI_DOUBLE_PRECISION, i, MPI_ANY_TAG, comm_coul, statusmpi, ierr)
            CALL MPI_RECV(efg_v22srf(j),     1, MPI_DOUBLE_PRECISION, i, MPI_ANY_TAG, comm_coul, statusmpi, ierr)
            CALL MPI_RECV(efg_v21msr(j),     1, MPI_DOUBLE_PRECISION, i, MPI_ANY_TAG, comm_coul, statusmpi, ierr)
            CALL MPI_RECV(efg_v22msr(j),     1, MPI_DOUBLE_PRECISION, i, MPI_ANY_TAG, comm_coul, statusmpi, ierr)
            CALL MPI_RECV(efg_vzz(j),        1, MPI_DOUBLE_PRECISION, i, MPI_ANY_TAG, comm_coul, statusmpi, ierr)
            CALL MPI_RECV(efg_winkel(j),     1, MPI_DOUBLE_PRECISION, i, MPI_ANY_TAG, comm_coul, statusmpi, ierr)
            CALL MPI_RECV(efg_eta(j),        1, MPI_DOUBLE_PRECISION, i, MPI_ANY_TAG, comm_coul, statusmpi, ierr)
            CALL MPI_RECV(efg_qmat(1,1,j),   9, MPI_DOUBLE_PRECISION, i, MPI_ANY_TAG, comm_coul, statusmpi, ierr)
            CALL MPI_RECV(efg_tens(1,j),     3, MPI_DOUBLE_PRECISION, i, MPI_ANY_TAG, comm_coul, statusmpi, ierr)
            CALL MPI_RECV(efg_eivec(1,1,j),  9, MPI_DOUBLE_PRECISION, i, MPI_ANY_TAG, comm_coul, statusmpi, ierr)

            CALL MPI_RECV(vcoul0(j),1, MPI_DOUBLE_PRECISION, i, MPI_ANY_TAG, comm_coul, statusmpi, ierr)
            CALL MPI_RECV(fhf(0,j), 4, MPI_DOUBLE_PRECISION, i, MPI_ANY_TAG, comm_coul, statusmpi, ierr)
            CALL MPI_RECV(logefg(j), 1, MPI_LOGICAL, i, MPI_ANY_TAG, comm_coul, statusmpi, ierr)
            CALL MPI_RECV(fk(0,j), 7, MPI_DOUBLE_PRECISION, i, MPI_ANY_TAG, comm_coul, statusmpi, ierr)
            CALL MPI_RECV(elnum(j), 1, MPI_DOUBLE_PRECISION, i, MPI_ANY_TAG, comm_coul, statusmpi, ierr)
          ENDDO
        ENDDO
      ENDIF


#endif 


      IF(SWITCH.EQ.'COUL') THEN
        CALL POISSN (LM,LMMAX)
          call abort_parallel
          STOP ' COUL OK'
      ENDIF

      if (switch1 .eq. ' EECE') then
        if (myid .EQ. 0) then
          do jatom=1,nat
            if (jeece(jatom) .eq. 1) then
              WRITE(50,*) JATOM,LMMAX(JATOM),' Atom, number of LM components'
              WRITE(51,*) JATOM,LMMAX(JATOM),' Atom, number of LM components'
              do L=0,6,2
                if (myid .EQ. 0) write(21,551) jatom,L,fk(L,jatom)*2.d0
                write(50,555) jatom,L,fk(L,jatom)*2.d0
                write(51,555) jatom,L,fk(L,jatom)*2.d0
 555            format(2i4,es20.12,'  atom, L, SLater integraL F^2L, Ry')
              enddo
            endif
          enddo
        endif
        close(50)
        close(51)
      endif

 445  CONTINUE

      WRITE(6,2031)
      WRITE(6,1540)

      CALL CPUTIM(TIME2)
      CALL wallTIM(TIME2w)
!
! calculate exchange potential for MT-spheres and for PW
!
!Task paralism
      IF(.NOT.xc) GOTO 888

 777  CONTINUE

      IF(.NOT.xc) goto 889
!      if (.not.coul) goto 888

!
!Only for Becke-Johnson potential or for potential fixed to case.r2v
!
      IF((VX_SWITCH.eq.'VX_BR89') .or. (VX_SWITCH.eq.'VX_SLATER') .or. &
         (VX_SWITCH.eq.'VX_MBJ') .or. (VX_SWITCH.eq.'VX_SMBJ') .or. &
         (VX_SWITCH.eq.'VX_GBJ') .or. (VX_SWITCH.eq.'VX_KLI') .or. &
         (VX_SWITCH.eq.'VX_R2V') .or. (VDISP_SWITCH.eq.'VDISP_NLVDW') .or. &
         (VHALF_SWITCH.ne.'VHALF_NONE') .or. (VX_SWITCH.eq.'VX_LMBJ') .or. &
         (id_func_x .eq. 207) .or. (id_func_x .eq. 208) .or. (id_func_x .eq. 209)) THEN

        ALLOCATE(KZZR2V(3,NKK))
        KZZR2V = 0

        if (VDISP_SWITCH.eq.'VDISP_NLVDW') then
           iunitr2v = 66
        else if (VHALF_SWITCH.ne.'VHALF_NONE') then
           iunitr2v = 69
        else
           iunitr2v = 11
        endif

        DO ISPIN = 1, JSPIN
            IF (ISPIN .EQ. 1) THEN
              OPEN(iunitr2v,FILE=FNAVN(iunitr2v))
            ELSE
              OPEN(iunitr2v,FILE=FNAVN(iunitr2v)(1:LEN_TRIM(FNAVN(iunitr2v)))//'dn')
            ENDIF
            READ(iunitr2v,'(/,/)',err=2826,end=2826)

!          DO JATOM = NSTART(MYID), NSTOP(MYID)
          DO JATOM = 1, nat

            READ(iunitr2v,'(/,/,/)')

            DO LM1 = 1, LMMAX(JATOM)

              READ(iunitr2v,'(/)')
              if(jatom.ge.NSTART(MYID).and.jatom.le.NSTOP(MYID)) then
                 READ(iunitr2v,2021) (VXCOLD(J,LM1,JATOM,ISPIN), J=1, JRI(JATOM))
              else
                 READ(iunitr2v,2021) (clm_tmp(j), J=1, JRI(JATOM))
              endif
              READ(iunitr2v,'(/)')

            ENDDO

            READ(iunitr2v,'(/,/,/)')

          ENDDO

            READ(iunitr2v,'(/)')
   !         READ(iunitr2v,'(1X,"        ",I10)') NKKR2V
            READ(iunitr2v,'(9X,I10)') NKKR2V
            IF (NKKR2V .GT. NKK) NKKR2V = NKK
            DO J = 1, NKKR2V
              READ(iunitr2v,2071) (KZZR2V(JX,J),JX=1,3), VKVXCOLD(J,ISPIN)
            ENDDO
        ENDDO

        IF (JSPIN .EQ. 1) THEN

        if (nstart(myid).le.nstop(myid))then
          VXCOLD(:,:,:,2) = VXCOLD(:,:,:,1)
        endif
          VKVXCOLD(:,2) = VKVXCOLD(:,1)

        ENDIF
        CLOSE(iunitr2v)
        goto 2825

 2826   continue
        CALL OUTERR('LAPW0','R2V files not present, which are required for mBJ, lmBJ, DFT-1/2, or NLvdW')
        GOTO 999

      ENDIF
 2825 continue
!
!Only for the Slater/KLI potential
!
      IF((VX_SWITCH.eq.'VX_SLATER') .or. (VX_SWITCH.eq.'VX_SMBJ') .or. &
         (VX_SWITCH.eq.'VX_KLI') .or. (EX_SWITCH(1:12).eq.'EX_DELTAGLLB') .or. &
         (VX_SWITCH.eq.'VX_LMBJ')) THEN

        ALLOCATE(KZZNONLOC(3,NKK))
        KZZNONLOC = 0

        DO ISPIN = 1, JSPIN
            IF (ISPIN .EQ. 1) THEN
              OPEN(64,FILE=FNAVN(64))
            ELSE
              OPEN(64,FILE=FNAVN(64)(1:LEN_TRIM(FNAVN(64)))//'dn')
            ENDIF
            READ(64,'(/,/)')

!          DO JATOM = NSTART(MYID), NSTOP(MYID)
          DO JATOM = 1, nat

            READ(64,'(/,/,/)')

            DO LM1 = 1, LMMAX(JATOM)

              READ(64,'(/)')
              if(jatom.ge.NSTART(MYID).and.jatom.le.NSTOP(MYID)) then
                 READ(64,2021) (VXCNONLOC(J,LM1,JATOM,ISPIN), J=1, JRI(JATOM))
              else
                 READ(64,2021) (clm_tmp(j), J=1, JRI(JATOM))
              endif
              READ(64,'(/)')

            ENDDO

            READ(64,'(/,/,/)')

          ENDDO

            READ(64,'(/)')
            READ(64,'(9X,I10)') NKKNONLOC
            IF (NKKNONLOC .GT. NKK) NKKNONLOC = NKK
            DO J = 1, NKKNONLOC
              READ(64,2071) (KZZNONLOC(JX,J),JX=1,3), VKVXCNONLOC(J,ISPIN)
            ENDDO
        ENDDO

        IF (JSPIN .EQ. 1) THEN

        if (nstart(myid).le.nstop(myid))then
          VXCNONLOC(:,:,:,2) = VXCNONLOC(:,:,:,1)
        endif
          VKVXCNONLOC(:,2) = VKVXCNONLOC(:,1)

        ENDIF
        CLOSE(64)

        goto 2835

 2836   continue
        CALL OUTERR('LAPW0','R2V_nonloc files not present, which are required for Slater/KLI')
        GOTO 999

      ENDIF
 2835 continue

      if (.not.coul) goto 888

!
! Exchange-correlation inside the spheres
!  
      DO JATOM=NSTART(myid),NSTOP(myid)

        JRJ=JRI(JATOM)
        DO J=1,JRJ
          R(J)=R0(JATOM)*EXP(dble(J-1)*DX(JATOM))
        ENDDO

        if (switch1 .eq. ' EECE') then
          if (jeece(jatom) .eq. 1) then
            if(.not.xcpot1qq) then
              CALL XCPOT1(luse2,LM,JATOM,LMMAX,LR2V,JSPIN, &
                  sigma(1,jatom),imfield,bextene,iprint,switch1,sigma_str)
            else
              if(MGGA_SWITCH) then      
                 call outerr("LAPW0", "xcpot1q not supported with MGGA")
                 goto 999
              endif   
              CALL XCPOT1Q(luse2,LM,JATOM,LMMAX,LR2V,JSPIN, &
                   sigma(1,jatom),imfield,bextene,iprint,switch1)
            endif
            IF(myid.EQ.0) WRITE(21,699) JATOM,luse2,JATOM,SIGMA(1,jatom)
            IF(myid.EQ.0) WRITE(21,699) JATOM,luse2,JATOM,SIGMA(2,jatom)
          endif
        else
            if(.not.xcpot1qq) then
              CALL XCPOT1(luse2,LM,JATOM,LMMAX,LR2V,JSPIN, &
                   sigma(1,jatom),imfield,bextene,iprint,switch1,sigma_str)
            else
              if(MGGA_SWITCH) then      
                 call outerr("LAPW0", "xcpot1q not supported with MGGA")
                 goto 999
              endif   
               CALL XCPOT1Q(luse2,LM,JATOM,LMMAX,LR2V,JSPIN, &
                    sigma(1,jatom),imfield,bextene,iprint,switch1)
            endif
        endif

      ENDDO
#ifdef Parallel
      call mpi_reduce(force_buf,force,1,MPI_LOGICAL,MPI_LOR,0,comm_coul,ierr)
      IF(myid.NE.0) THEN
        DO j=nstart(myid),nstop(myid)
          CALL MPI_SEND(sigma(1,j), 2, MPI_DOUBLE_PRECISION, 0, 0, comm_coul, ierr)
        ENDDO
      ELSE
        DO i=1,icoul
          DO j=nstart(i),nstop(i)
            CALL MPI_RECV(sigma(1,j), 2, MPI_DOUBLE_PRECISION, i, &
              MPI_ANY_TAG, comm_coul, statusmpi, ierr)
          ENDDO
        ENDDO
      ENDIF
#else
      force=force_buf
#endif

!    COMPUTE VXC*RHO CORE
     if (core_add .or. core_sub) then

!       MPI task 0 reads all data and distributes it
        if (myid.eq.0) then
           allocate(rhoc(nrad,nat,jspin))
           rhoc = 0.0d0
           do ispin=1,jspin
              do jatom=1,nat
                 read(92+ispin,'(8/)')
                 read(92+ispin,'(3x,4e19.12)') (rhoc(i,jatom,ispin), i=1, jri(jatom))
                 read(92+ispin,'(2/)')
              enddo
           enddo
           rhoc = rhoc/sqrt(4d0*pi)/dble(jspin)
        end if

!       TASK 0 sends; other ids receive their rhoc of "their" atoms
#ifdef Parallel

        if (myid .ne. 0) then
           allocate(rhoc(nrad,nstart(myid):nstop(myid),jspin))
           rhoc=0.d0
        end if
        
        do ispin=1,jspin
           do jatom=1,nat

              if (myid .eq. 0) then
!             find which MPI task should receive jatom
                 j = findloc( jatom .ge. nstart(:) .and. jatom .le. nstop(:), .true., 1) - 1
                 if (j .eq. 0) cycle
                 
                 call mpi_send(rhoc(:,jatom,ispin), jri(jatom), MPI_DOUBLE_PRECISION, j, 0, comm_coul, ierr)
                 if (ierr .ne. MPI_SUCCESS) STOP 'ERROR in mpi_send of rhoc'
                 
              else if (jatom .ge. nstart(myid) .and. jatom .le. nstop(myid)) then

                 call mpi_recv(rhoc(:,jatom,ispin), jri(jatom), MPI_DOUBLE_PRECISION, 0, 0, comm_coul, statusmpi, ierr)
                 if (ierr .ne. MPI_SUCCESS) STOP 'ERROR in mpi_recv of rhoc'

              end if
              
           end do
        end do
#endif
        
        vxcrhoc = 0.0d0
        vxcrhoc_buf = 0.0d0
        do ispin=1,jspin
           do jatom=nstart(myid),nstop(myid)

              jrj = jri(jatom)
              do i=1,jrj
                 r(i)=r0(jatom)*exp(dble(i-1)*dx(jatom))
              enddo

              value(1:jrj) = vxc(1:jrj,1,jatom,ispin)*rhoc(1:jrj,jatom,ispin)
              call charge(r,dx(jatom),value,1,jrj,qin)
              vxcrhoc_buf = vxcrhoc_buf + mult(jatom)*qin
              
           enddo
        enddo
        
#ifdef Parallel
        call mpi_reduce(vxcrhoc_buf,vxcrhoc,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,comm_coul,ierr)
        if (ierr  .ne. MPI_SUCCESS) STOP 'ERROR im mpi_reduce of vxcrhoc'
#else
        vxcrhoc = vxcrhoc_buf
#endif
        deallocate(rhoc)

        if(myid .eq. 0) then
           if (core_sub) then
              write(21,'(":VXCRHOC2 : ",1x,f22.8,3X,"(Ry)")') vxcrhoc
           elseif (core_add) then
              write(21,'(":VXCRHOC1 : ",1x,f22.8,3X,"(Ry)")') vxcrhoc
           endif
        end if

     endif

!    COMPUTE VTAU*TAU CORE
     if (core_sub .and. MGGA_TAU_SWITCH) then

        if (myid.eq.0) then
           allocate(tauc(nrad,nat,jspin))
           tauc = 0.0d0
           do ispin=1,jspin
              do jatom=1,nat
                 read(94+ispin,'(8/)')
                 read(94+ispin,'(3x,4e19.12)') (tauc(i,jatom,ispin), i=1, jri(jatom))
                 read(94+ispin,'(2/)')
              enddo
           enddo
           tauc = 0.5d0*tauc/sqrt(4d0*pi)/dble(jspin) ! Ry -> Ha
        end if

#ifdef Parallel
        
        if (myid .ne. 0) then
           allocate(tauc(nrad,nstart(myid):nstop(myid),jspin))
           tauc=0.d0
        end if
        
        do ispin=1,jspin
           do jatom=1,nat

              if (myid .eq. 0) then
!             find which MPI task should receive jatom
                 j = findloc( jatom .ge. nstart(:) .and. jatom .le. nstop(:), .true., 1) - 1
                 if (j .eq. 0) cycle

                 call mpi_send(tauc(:,jatom,ispin), jri(jatom), MPI_DOUBLE_PRECISION, j, 0, comm_coul, ierr)
                 if (ierr .ne. MPI_SUCCESS) STOP 'ERROR in mpi_send of tauc'

              else if (jatom .ge. nstart(myid) .and. jatom .le. nstop(myid)) then
                 
                 call mpi_recv(tauc(:,jatom,ispin), jri(jatom), MPI_DOUBLE_PRECISION, 0, 0, comm_coul, statusmpi, ierr)
                 if (ierr .ne. MPI_SUCCESS) STOP 'ERROR in mpi_recv of tauc'

              end if

           end do
        end do
#endif
        
        vtautauc = 0.0d0
        vtautauc_buf = 0.0d0
        do ispin=1,jspin
           do jatom=nstart(myid),nstop(myid)

              jrj = jri(jatom)
              do i=1,jrj
                 r(i)=r0(jatom)*exp(dble(i-1)*dx(jatom))
              enddo

              value(1:jrj) = vtau(1:jrj,1,jatom,ispin)*tauc(1:jrj,jatom,ispin)
              call charge(r,dx(jatom),value,1,jrj,qin)
              vtautauc_buf = vtautauc_buf + mult(jatom)*qin
              
           enddo
        enddo

        
#ifdef Parallel
        call mpi_reduce(vtautauc_buf,vtautauc,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,comm_coul,ierr)
        if (ierr .ne. MPI_SUCCESS) STOP 'ERROR in mpi_reduce of vtautauc'
#else
        vtautauc = vtautauc_buf
#endif
        deallocate(tauc)

        if (myid .eq. 0) write(21,'(":VTAUTAUC : ",1x,f22.8,3X,"(Ry)")') vtautauc

     endif

     
 888  CONTINUE

#ifdef Parallel
      if (.not.coul) allocate(potk(nkk))
       ibbuff=nkk     !  set it to lower value if you have a lousy mpi, 32768 
       i2=NKK/ibbuff
       it1=1
       do i=1,i2
         call MPI_Bcast(potk(it1),ibbuff,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
         it1=it1+ibbuff
       enddo
       if(it1 .lt. nkk)then
         it2=nkk-it1+1
         call MPI_Bcast(potk(it1),it2,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
       endif
       call MPI_BARRIER(MPI_COMM_WORLD,ierr) 
!      call MPI_Bcast(potk, NKK, MPI_DOUBLE_COMPLEX, 0,MPI_COMM_WORLD , ierr)
#endif

      CALL CPUTIM(TIME3)
      CALL wallTIM(TIME3w)
  
      if(switch1.eq.' EECE' .or. vsp_only) goto 3456  !EECE jump over interstitial

!
! Exchange-correlation in the interstitial
!
      deallocate(absk)
      allocate(vkexc(1:nkk,1:2),vkvxc(1:nkk,1:2))
      VKEXC(1:nkk,1:2)=ZEROC
      VKVXC(1:nkk,1:2)=ZEROC


      if (switch3 .eq. 'YES' .or. abs(delcc) .gt. 1.d-14) then
        allocate(cfft2(iff1,iff2,iff3),cfft3(iff1,iff2,iff3))
        jff1 = iff1
        jff2 = iff2
        jff3 = iff3
        cfft2 = ZEROC
        cfft3 = ZEROC
      else
        allocate(cfft2(1,1,1),cfft3(1,1,1))
        jff1 = 1
        jff2 = 1
        jff3 = 1
     endif
     if (MGGA_SWITCH) then
        allocate (META_FFT(IFF1,IFF2,IFF3,2))
        allocate (vk_vtau_nonrean(nkk,2))
     else
        allocate (META_FFT(1,1,1,1))
        allocate (vk_vtau_nonrean(1,1))
     end if
     meta_fft = zeroc
     vk_vtau_nonrean = zeroc
     
     !stress tensor
     IF(SW_STRESS .eq. 'YES') then
       allocate( CFFT_STR(iff1,iff2,iff3,1:9) )    !for coulomb stress int
       allocate( CFFT8(iff1,iff2,iff3) )           !exc_vxc in int
       IF(GGA_SWITCH .eqv. .TRUE.) allocate( CFFT_STR_GGA(iff1,iff2,iff3,1:9) )
        write(6,'(":FFT    allocation of stress for GGA ",i5," MB")') &
          iff1*iff2*iff3*19*16/1000000
     ELSE
       allocate( CFFT_STR(1,1,1,1:9) )
       allocate( CFFT8(1,1,1) )
       IF(GGA_SWITCH .eqv. .TRUE.) allocate( CFFT_STR_GGA(1,1,1,1:9) )
     ENDIF
     CFFT_STR = ZEROC
     CFFT8    = ZEROC
     IF(GGA_SWITCH .eqv. .TRUE.) CFFT_STR_GGA = ZEROC              
     !end stress tensor
     CALL XCPOT3 (NKK,NKK2,NKK3,NKKR2V,NKKNONLOC,&
          VKEXC,VKVXC,VKVXCOLD,VKVXCNONLOC,vk_vtau_nonrean,LR2V,JSPIN, &
          IFF1,IFF2,IFF3,JFF1,JFF2,JFF3,CFFT,CFFT2,CFFT3, &
          UST,FFT,POTK,imfield,bextene,switch3,delcc)
!
!     OUTPUT OF XCPOT3: NOT REANALYSED!!!
!
!     1. CFFT(:,:,:) holds total energy (rho*v)
!     2. CFFT2(:,:,:) only for charged cells
!     3. CFFT3(:,:,:) only for switch3
!     4. VKEXC and VKVXC star coefficients of exc/vxc 
!     5. FFT(:,:,:,1/2) holds potential (vkvxc+potk) up/down
!     6. UST(:,:,:,1/2) holds total energy (vkvxc + potk - exc)*rho up/down
!     7. META_FFT: VTAU or VLAPL
!
      SW_TOT_STR = .true.  !added bz Kamal to limit print options in rean3
      
      IF(SWITCH.EQ.'EXCH') then
          call abort_parallel
          STOP 'EXCH OK'
      endif

      CALL CPUTIM(TIME4)
      CALL wallTIM(TIME4w)

      ! generate the step function
      CALL REAN0(NKK,KZZ,IFF1,IFF2,IFF3,UST)

      ISPIN=1
      CALL CPUTIM(TIME4a)
      CALL wallTIM(TIME4aw)
      allocate(cvalue(nkk))

      if (switch3 .eq. 'YES') then
        CALL REAN3(NKK,KZZ,CVALUE,IFF1,IFF2,IFF3,CFFT2,UST,UST0)
        EKININT = CVALUE(1)*UST0
        CALL REAN3(NKK,KZZ,CVALUE,IFF1,IFF2,IFF3,CFFT3,UST,UST0)
        EXCINT = CVALUE(1)*UST0
      endif

      if (abs(delcc) .gt. 1.d-14) then
        CALL REAN3(NKK,KZZ,CVALUE,IFF1,IFF2,IFF3,CFFT2,UST,UST0)
        EKININT = CVALUE(1)*UST0
        ! print*,'ekinint',ekinint
      endif

      deallocate(cfft2,cfft3)      

!     step function for k=0
      CALL REAN3(NKK,KZZ,CVALUE,IFF1,IFF2,IFF3,CFFT,UST,UST0)

      EKVOL1=CVALUE(1)*UST0
      CALL CPUTIM(TIME4b)
      CALL wallTIM(TIME4bw)
      IF(SWITCH.EQ.'TOT ') then
        if (coul .and. (VX_SWITCH.ne.'VX_LGRR') .and. .not. vsp_only) CALL ENERGY(VCOUL0,QEL,LM,LMMAX, &
          ekvol1,ekinint,excint,etotal,ekin,exc2,etot_bg,switch3,iprint,delcc)
      endif
      CALL CPUTIM(TIME5)
      CALL wallTIM(TIME5w)
      
      
      !=========Electrostatic stress================            
      SW_TOT_STR = .false.
      IF( SW_STRESS .eq. 'YES' ) THEN      
      call walltim(time1_st)
         CALL REAN3(NKK,KZZ,CVALUE,IFF1,IFF2,IFF3,CFFT8,UST,UST0)
         E_xc_m_V_xc  = CVALUE(1)*UST0     !Exc_vxc in interstitial                        
         
         DO index_p = 1, 9
            CALL REAN3(NKK,KZZ,CVALUE,IFF1,IFF2,IFF3,CFFT_STR(1,1,1,index_p),UST,UST0)
            coulomb_str_int(index_p) = CVALUE(1)*UST0                     
                   
           IF(GGA_SWITCH .eqv. .TRUE.) THEN           
!               DO index_p = 1, 9
                 CALL REAN3(NKK,KZZ,CVALUE,IFF1,IFF2,IFF3,CFFT_STR_GGA(1,1,1,index_p),UST,UST0)
                 gga_str_int(index_p) = CVALUE(1)*UST0        
!               ENDDO            
           ENDIF         
         ENDDO   

!        Add the missing -delta_ab/2 * int_I rho V_C part of Eq. (G5)
!        only in diagnostic mode.  Indices 1,5,9 are xx,yy,zz in the
!        flattened tensor convention used throughout LAPW0.
         electrostatic_int_jacobian = EKVOL1 + E_xc_m_V_xc
         if(myid.eq.0) write(21,2091) electrostatic_int_jacobian, &
              stress_int_delta_diag
         if(stress_int_delta_diag) then
            coulomb_str_int(1) = coulomb_str_int(1) - &
                                 electrostatic_int_jacobian
            coulomb_str_int(5) = coulomb_str_int(5) - &
                                 electrostatic_int_jacobian
            coulomb_str_int(9) = coulomb_str_int(9) - &
                                 electrostatic_int_jacobian
         endif
         !columb+Madelung stress+exchange-correlation for lda
         call exchange_str(lm,lmmax,E_xc_m_V_xc_rmt)
         IF(GGA_SWITCH .eqv. .TRUE.) then
           call stress_gga(lm,lmmax, jspin,gga_str_int,sigma_str,lmmaxx)      
           if(allocated(sigma_str)) deallocate( sigma_str)
         ENDIF  
         call finl_elect_str(lm,lmmax,coulomb_str_int,E_xc_m_V_xc,E_xc_m_V_xc_rmt)  
         
         deallocate( CFFT_STR,CFFT8 )   !from stress
         deallocate( SGGP )             !from stress
         deallocate( DECOUL_POT_PW )          
         deallocate( SGGP_GGA )
         deallocate( gga_int ) 
           
      call walltim(time2_st)
      time_electro = time_electro + (time2_st-time1_st)
      
      ENDIF
      SW_TOT_STR = .true.
      !*********************************************
      !==========Electrostatic Kamal================      
      
!
! transform plane wave total potential to zero inside spheres
! (stepfunction)
! CVALUE(K) is now the reanalyzed  PW  potential  coefficient
!
      CALL REAN4 (NKK,IFF1,IFF2,IFF3,UST,FFT,VKeXC)
      deallocate(cfft,fft)

      if (MGGA_LAPL_SWITCH) then
         allocate (vk_vlapl(nkk,2))
         CALL REAN_DIRECT_FFT (NKK,IFF1,IFF2,IFF3,UST,VK_VLAPL)         
         deallocate(meta_fft)
         vk_vlapl = vk_vlapl/vol
      else if (MGGA_TAU_SWITCH) then
         allocate (vk_vtau(nkk,2))
!         REAN_DIRECT_FFT COULD BE USED IF NO VTAU PLOTTING IS NECESSARY
!         CALL REAN_DIRECT_FFT (NKK,IFF1,IFF2,IFF3,UST,META_FFT,VK_VTAU)


! regenerate the step function
         CALL REAN0(NKK,KZZ,IFF1,IFF2,IFF3,UST)
         CALL REAN4(NKK,IFF1,IFF2,IFF3,UST,META_FFT,VK_VTAU)
         
         deallocate(meta_fft)
         vk_vtau = vk_vtau/vol
      end if
      
      deallocate(ust)

      IF(.NOT.coul) goto 889
!  
! sum up the potentials and do the output
!
      SCFWRITE: IF(myid.EQ.0) THEN

        IF(myid.EQ.0) then

          DO jatom=1,nat
            WRITE(21,1335) jatom,jatom,aname(jatom),vcoul0(jatom)
            IF(logefg(jatom)) THEN
              WRITE(21,1340)  jatom,efg_vzz(jatom)
              WRITE(21,1341)  efg_V20(jatom),efg_V20SRF(jatom)
              WRITE(21,1342)  efg_V22(jatom),efg_V22SRF(jatom)
              WRITE(21,1343)  efg_V22M(jatom),efg_V22MSR(jatom)
              WRITE(21,1344)  efg_V21(jatom),efg_V21SRF(jatom)
              WRITE(21,1345)  efg_V21M(jatom),efg_V21MSR(jatom)
              WRITE(21,1350)  efg_qmat(1,1,jatom),efg_qmat(1,2,jatom), &
                efg_qmat(1,3,jatom),efg_TENS(1,jatom),0.,0.
              WRITE(21,1350)   efg_qmat(2,1,jatom),efg_qmat(2,2,jatom), &
                efg_qmat(2,3,jatom),0.,efg_TENS(2,jatom),0.
              WRITE(21,1350)   efg_qmat(3,1,jatom),efg_qmat(3,2,jatom), &
                efg_qmat(3,3,jatom),0.,0.,efg_TENS(3,jatom)
              WRITE(21,1355)   ((efg_eivec(i,j,jatom),j=1,3),i=1,3), &
                jatom,efg_winkel(jatom)
              WRITE(21,1360)  jatom,efg_eta(jatom)
            ENDIF
          ENDDO

          IF (force) THEN
! Efeld term
            if(abs(iefeld) .gt. 0) call EFORCE(FHF)
            DO jatom=1,nat
              WRITE(21,710)
               fhf(0,jatom)=SQRT (fhf(1,jatom)**2 + fhf(2,jatom)**2 &
                 + fhf(3,jatom)**2)
              WRITE(21,711)jatom,jatom,(fhf(ik,jatom)*1.d3,ik=0,3)
            ENDDO
          ENDIF

          DO jatom=1,nat
            WRITE(21,699) JATOM,luse2,JATOM,SIGMA(1,jatom)
            WRITE(21,699) JATOM,luse2,JATOM,SIGMA(2,jatom)
          ENDDO
     
! write(6,*) 'total energy correction due to charge cell (q*Z/2)',etotcorr/2.d0
! write(21,*) '          total energy correction due to charge cell (q*Z/2)',etotcorr/2.d0
!          if((iex .lt. 50) .or. (iex .gt. 55)) then
          if((GRR_SWITCH .eqv. .false.) .and. (vsp_only .eqv. .false.)) THEN
             WRITE(6,1121) ETOTAL      
             WRITE(21,1121) ETOTAL 
          elseif((GRR_SWITCH.eqv..true.) .and. &
                 (EX_SWITCH.eq.'EX_GRR') .and. (VX_SWITCH.eq.'VX_GRR')) then 
             WRITE(6,1120) ETOTAL/vol      
             WRITE(21,1120) ETOTAL/vol 
          endif
          if (switch3 .eq. 'YES') then
            if((GRR_SWITCH.eqv..true.) .and. &
               (EX_SWITCH(4:4).eq.'S') .and. (VX_SWITCH(4:4).eq.'S') .and. &
               (EC_SWITCH(4:7).eq.'NONE') .and. (VC_SWITCH(4:7).eq.'NONE')) then
              WRITE(6,1124) alphax*EXC2
              WRITE(21,1124) alphax*EXC2
            else
              WRITE(6,1122) EKIN
              WRITE(21,1122) EKIN
              if (EX_SWITCH(1:12).eq.'EX_DELTAGLLB') then
                WRITE(6,1125) EXC2
                WRITE(21,1125) EXC2
              else
                WRITE(6,1123) EXC2
                WRITE(21,1123) EXC2
              endif
            endif
         endif
          if(abs(delcc) .gt. 1.d-14) then
            write(21,'(":EBG   :  Background correction",f20.6)') etot_bg
          endif
! Add efeld terms for convenience here to potk (not used in fact!)
          if(abs(iefeld).gt.0) call epot3(potk,nkk)
! calculate potential at z=0 for workfunction (requires origin in interstitial)
           if(.not. (GRR_SWITCH )) call workf1(kzz,potk,vkvxc,nkk)

        endif

      ENDIF SCFWRITE

3456  continue !EECE

      IF(.NOT.coul) goto 889
 
      if (nstart(myid).le.nstop(myid)) &
        allocate(VTOT(NRAD,lmmaxx,nstart(myid):nstop(myid)))

      SPIN: DO ISPIN=1,JSPIN

        if(switch1.eq.' EECE') goto 3458

        if (.not. vsp_only) then
           CVALUE(1:NKK)=(0.d0,0.d0)
           TMP=1.D0/VOL
           DO JRE=1,NKK
              CVALUE(JRE)=VKeXC(JRE,ISPIN)*tmp
           ENDDO
        end if
           
        DO JATOM=NSTART(myid),NSTOP(myid)
          DO J=1,JRI(JATOM)
            R(J)=R0(JATOM)*EXP(DX(JATOM)*dble(J-1))
          ENDDO
          if(abs(iefeld) .gt. 0)call epot1c(LM,luse2,JATOM,LMMAX,vfield)
          DO LM1=1,LMMAX(JATOM)
            DO J=1,JRI(JATOM)
! E-field, added completely after the energy terms
              if(abs(iefeld) .gt. 0)then
                VTOT(J,LM1,JATOM) = V(J,LM1,JATOM) + &
                  VXC(J,LM1,JATOM,ISPIN) +vfield(LM1,J)
              else
                VTOT(J,LM1,JATOM) = V(J,LM1,JATOM) + &
                  VXC(J,LM1,JATOM,ISPIN)
              endif

            ENDDO
          ENDDO
        ENDDO

3458    continue

        close(9)
        close(10)
        close(11)
        close(16)
        close(17)
        close(18)
        close(19)
        
        close(15)

        inquire(unit=88, opened=isopen)
        if (isopen) close(88)
        inquire(unit=89, opened=isopen)
        if (isopen) close(89)
        inquire(unit=91, opened=isopen)
        if (isopen) close(91)
        inquire(unit=92, opened=isopen)
        if (isopen) close(92)
     
#ifdef Parallel            
        IF(myid.NE.0) THEN
          CALL MPI_RECV(i,1,MPI_INTEGER,myid-1, MPI_ANY_TAG,comm_coul &
            ,statusmpi,ierr)
        ENDIF
#endif

        if(switch1.eq.' EECE' .or. vsp_only) goto 3457

        IF (LR2V) THEN

          IF(myid.EQ.0) THEN
            if(ispin.eq.1) then
              OPEN(9,FILE=fnavn(9))
              OPEN(11,file=fnavn(11))
            else
              OPEN(9,FILE=fnavn(9)(1:LEN_TRIM(fnavn(9)))//'dn')
              OPEN(11,file=fnavn(11)(1:LEN_TRIM(fnavn(11)))//'dn')
            endif
            WRITE(9,*) '    TOTAL POTENTIAL'
            WRITE(9,*) '    NORM OF VC(L,M)=VC(L,M)*r**2'
            WRITE(9,*)
            WRITE(11,1970) ISCF
            WRITE(11,78)
            WRITE(11,77)
          else
            if(ispin.eq.1) then
              OPEN(9,FILE=fnavn(9),POSITION='append')
              OPEN(11,file=fnavn(11),POSITION='append')
            else
              OPEN(9,FILE=fnavn(9)(1:LEN_TRIM(fnavn(9)))//'dn',POSITION='append')
              OPEN(11,file=fnavn(11)(1:LEN_TRIM(fnavn(11)))//'dn',POSITION='append')
            endif
          ENDIF
        
          DO jatom=nstart(myid),nstop(myid)
            WRITE(9,1990) JATOM
            WRITE(9,2001) LMMAX(JATOM)
            WRITE(11,1990) JATOM                                        
            WRITE(11,2001) LMMAX(JATOM)                                 
            DO J=1,JRI(JATOM)
              R(J)=R0(JATOM)*EXP(DX(JATOM)*dble(J-1))
            ENDDO
!if(jatom.eq.1) then
! DO j=1,jri(jatom)
! if(ispin.eq.1)write(878,'(3f27.20)') r(j), vxc(j,1,jatom,1), vxc(j,1,jatom,1)*r(j)*r(j)
! if(ispin.eq.2)write(879,'(3f27.20)') r(j), vxc(j,1,jatom,2), vxc(j,1,jatom,2)*r(j)*r(j)
! enddo
!endif
              DO LM1=1,LMMAX(JATOM)
                 WRITE(9,2011) LM(1,LM1,JATOM),LM(2,LM1,JATOM)
                 WRITE(9,2021) (VTOT(J,LM1,JATOM)*R(J)*R(J),J=1,JRI(JATOM))
                 WRITE(9,2031)
                 WRITE(11,2011) LM(1,LM1,JATOM),LM(2,LM1,JATOM)
                 WRITE(11,2021) (VXC(J,LM1,JATOM,ispin)*r(j)*r(j),J=1,JRI(JATOM))
                 WRITE(11,2031)                                           
              ENDDO
              WRITE(9,2030)
              WRITE(11,2030)
           ENDDO
           
           IF(myid.EQ.icoul) THEN
              WRITE(9,*) '    TOTAL POTENTIAL IN INTERSTITIAL'
              WRITE(9,*)
              WRITE(9,2061) nkk
              WRITE(11,*) '    EXCHANGE-CORRELATION POTENTIAL IN INTERSTITIAL'
              WRITE(11,*)
              WRITE(11,2061) nkk
              if(inversion)then
                 DO j=1,nkk
                    WRITE(9,2071) (kzz(jx,j),jx=1,3), dble(potk(j)+vkvxc(j,ispin)),0.D0
                    WRITE(11,2071) (kzz(jx,j),jx=1,3), dble(vkvxc(j,ispin)),0.D0
                 ENDDO
              else
                 DO j=1,nkk
                    WRITE(9,2071) (kzz(jx,j),jx=1,3), potk(j)+vkvxc(j,ispin)
                    WRITE(11,2071) (kzz(jx,j),jx=1,3), vkvxc(j,ispin)
                 ENDDO
              endif
           ENDIF
        
           close(9)
           close(11)
           
           if (MGGA_TAU_SWITCH) then
              IF(myid.EQ.0) THEN
                 if(ispin.eq.1) then
                    OPEN(15,file=fnavn(15))
                 else
                    OPEN(15,file=fnavn(15)(1:LEN_TRIM(fnavn(15)))//'dn')
                 endif
                 WRITE(15,1970) ISCF
                 WRITE(15,78)
                 WRITE(15,77)
              else
                 if(ispin.eq.1) then
                    OPEN(15,file=fnavn(15),POSITION='append')
                 else
                    OPEN(15,file=fnavn(15)(1:LEN_TRIM(fnavn(15)))//'dn',POSITION='append')
                 endif
              ENDIF
              
              DO jatom=nstart(myid),nstop(myid)
                 WRITE(15,1990) JATOM                                        
                 WRITE(15,2001) LMMAX(JATOM)                                 
                 DO J=1,JRI(JATOM)
                    R(J)=R0(JATOM)*EXP(DX(JATOM)*dble(J-1))
                 ENDDO
                 
                 DO LM1=1,LMMAX(JATOM)
                    WRITE(15,2011) LM(1,LM1,JATOM),LM(2,LM1,JATOM)
                    WRITE(15,2021) (VTAU(J,LM1,JATOM,ispin)*r(j)*r(j),J=1,JRI(JATOM))
                    WRITE(15,2031)                                           
                 ENDDO
                 WRITE(15,2030)
              ENDDO
              
              IF(myid.EQ.icoul) THEN
                 WRITE(15,*) '    EXCHANGE-CORRELATION POTENTIAL IN INTERSTITIAL'
                 WRITE(15,*)
                 WRITE(15,2061) nkk
                 if(inversion)then
                    DO j=1,nkk
                       WRITE(15,2071) (kzz(jx,j),jx=1,3), dble(vk_vtau_nonrean(j,ispin)),0.D0
                    ENDDO
                 else
                    DO j=1,nkk
                       WRITE(15,2071) (kzz(jx,j),jx=1,3), vk_vtau_nonrean(j,ispin)
                    ENDDO
                 endif
              ENDIF
              close(15)
           endif
        endif

        close(9)
        close(11)
!
! Write out spherical part of total potential to TAPE 16,17 = VSP
! PW and nonspherical part of total potential to TAPE 18,19 = VNS
! PW and nonspherical part of laplacian-dependent part of potential to TAPE 88,89 = vlapl
! PW and nonspherical part of tau-dependent part of potential to TAPE 98,99 = vtau
!

3457    continue

        if(switch1.ne.' EECE')then

          ITAPS=15+ISPIN
          ITAPN=17+ISPIN

          IF(myid.EQ.0) THEN
            open(itaps,file=fnavn(itaps))
            open(itapn,file=fnavn(itapn))
            WRITE(itaps,787) iscf
            WRITE(itaps,*) '    NORM:  V*R'
            WRITE(itaps,*)
            WRITE(itapn,*) '    TOTAL NONSPHERICAL POTENTIAL OF MT-SPHERES'
            WRITE(itapn,*) '    NORM  OF VLM(R) = VLM(R)'
            WRITE(itapn,*)
          else
            open(itaps,file=fnavn(itaps),POSITION='append')
            open(itapn,file=fnavn(itapn),POSITION='append')
          ENDIF

        endif

        DO jatom=nstart(myid),nstop(myid)

          JRJ=JRI(JATOM)
          if(switch1.ne.' EECE')then
            WRITE(itaps,1990) jatom
            WRITE(itapn,1990) jatom
            WRITE(itaps,2001) 1
            WRITE(itapn,2001) lmmax(jatom)-1
          endif

!EECE potential output
          if (switch1 .eq. ' EECE') then

            if (jeece(jatom) .eq. 1) then

              if (myid.EQ.0) then
                open(49+ispin,file=fnavn(49+ispin),POSITION='append')
              else
                open(49+ispin,file=fnavn(49+ispin),POSITION='append')
              endif

              DO J=1,JRI(JATOM)
                R(J)=R0(JATOM)*EXP(DX(JATOM)*(J-1))
              ENDDO

              IF (IATNR(JATOM) .LT. 0) THEN !non-cubic case

                espxc = 0d0
                vspxc = 0d0

                DO LM1=1,LMMAX(JATOM)

!espxc = integral of exc*rho inside the sphere
                  do j=1,jri(jatom)
                    value(j) = clmsp(j,lm1,jatom,ispin)*exc(j,lm1,jatom,ispin)
                    if (lm1.eq.1) value(j) = value(j)/sqfp
                  enddo
                  CALL CHARGE (R,DX(JATOM),VALUE,1,jrj,VLIN1)   ! Was charg2
!                  ACORR = VALUE(1)*R(1)*0.5d0
                  espxc = espxc + vlin1 !+ acorr

!vspxc = integral of vxc*rho inside the sphere
                  do j=1,jri(jatom)
                    value(j) = clmsp(j,lm1,jatom,ispin)*vxc(j,lm1,jatom,ispin)
                    if (lm1.eq.1) value(j) = value(j)/sqfp
                  enddo
                  CALL CHARGE (R,DX(JATOM),VALUE,1,jrj,VLIN1)   ! Was charg2
!                  ACORR = VALUE(1)*R(1)*0.5d0
                  vspxc = vspxc + vlin1 !+ acorr

!calculation of <g|exc|g> radial part
                  do j=1,jri(jatom)
                    value(j) = exc(j,lm1,jatom,ispin)*clmsp(j,1,jatom,ispin)
                  enddo
                  CALL CHARGE (R,DX(JATOM),VALUE,1,jrj,VLIN1)   ! Was charg2
!                  ACORR = VALUE(1)*R(1)*0.5D0
                  if (elnumsp(jatom,ispin) .gt. 1D-12) then
!                    fexck = (vlin1+acorr)/elnumsp(jatom,ispin)
                    fexck = (vlin1)/elnumsp(jatom,ispin)
                  else
                    fexck = 0d0
                  endif

!calculation of <g|vxc|g> radial part
                  do j=1,jri(jatom)
                    value(j) = vxc(j,lm1,jatom,ispin)*clmsp(j,1,jatom,ispin)
                  enddo
                  CALL CHARGE (R,DX(JATOM),VALUE,1,jrj,VLIN1)   ! Was charg2
!                  ACORR = VALUE(1)*R(1)*0.5D0
                  if (elnumsp(jatom,ispin) .gt. 1D-12) then
!                    fvxck = (vlin1+acorr)/elnumsp(jatom,ispin)
                    fvxck = (vlin1)/elnumsp(jatom,ispin)
                  else
                    fvxck = 0d0
                  endif

                  write(6,652) jatom, lm1, (lm(i,lm1,jatom),i=1,2), fexck, fvxck
 652              format(4i4,2f12.6,' atom, lm1, l, m, fexck, fvxck')

                  write(49+ispin,653) jatom, lm1, (lm(i,lm1,jatom),i=1,2), fexck, fvxck
 653              format(4i4,2es20.12)

                ENDDO

                write(6,654) jatom, espxc, vspxc
 654            format(i4,2f12.6,' atom, integral of Exc, integral of Vxc ')

                write(49+ispin,655) jatom, espxc, vspxc
 655            format(i4,2es20.12)

              ELSEIF (IATNR(JATOM) .GT. 0) THEN !cubic case

                ESPXC = 0D0
                VSPXC = 0D0

                DO LM1=1, LMMAX(JATOM)

                  L = IABS(LM(1,LM1,JATOM))
                  M = LM(2,LM1,JATOM)

                  IF ((L .EQ. 0) .OR. (L .EQ. 3)) THEN

!espxc = integral of exc*rho inside the sphere
                    DO J=1, JRI(JATOM)
                      VALUE(J) = CLMSP(J,LM1,JATOM,ISPIN)*EXC(J,LM1,JATOM,ISPIN)
                      IF (L .EQ. 0) VALUE(J) = VALUE(J)/SQFP
                    ENDDO
                    CALL CHARGE(R,DX(JATOM),VALUE,1,JRJ,VLIN1)   ! Was charg2
!                    ACORR = VALUE(1)*R(1)*0.5D0
                    ESPXC = ESPXC + VLIN1 !+ ACORR

!vspxc = integral of vxc*rho inside the sphere
                    DO J=1, JRI(JATOM)
                      VALUE(J) = CLMSP(J,LM1,JATOM,ISPIN)*VXC(J,LM1,JATOM,ISPIN)
                      IF (L .EQ. 0) VALUE(J) = VALUE(J)/SQFP
                    ENDDO
                    CALL CHARGE(R,DX(JATOM),VALUE,1,JRJ,VLIN1)   ! Was charg2
!                    ACORR = VALUE(1)*R(1)*0.5D0
                    VSPXC = VSPXC + VLIN1 !+ ACORR

!calculation of <g|exc|g> radial part
                    DO J=1, JRI(JATOM)
                      VALUE(J) = EXC(J,LM1,JATOM,ISPIN)*CLMSP(J,1,JATOM,ISPIN)
                    ENDDO
                    CALL CHARGE(R,DX(JATOM),VALUE,1,JRJ,VLIN1)   ! Was charg2
!                    ACORR = VALUE(1)*R(1)*0.5D0
                    IF (ELNUMSP(JATOM,ISPIN) .GT. 1D-12) THEN
!                      FEXCK = (VLIN1+ACORR)/ELNUMSP(JATOM,ISPIN)
                      FEXCK = (VLIN1)/ELNUMSP(JATOM,ISPIN)
                    ELSE
                      FEXCK = 0d0
                    ENDIF

!calculation of <g|vxc|g> radial part
                    DO J=1, JRI(JATOM)
                      VALUE(J) = VXC(J,LM1,JATOM,ISPIN)*CLMSP(J,1,JATOM,ISPIN)
                    ENDDO
                    CALL CHARGE(R,DX(JATOM),VALUE,1,JRJ,VLIN1)   ! Was charg2
!                    ACORR = VALUE(1)*R(1)*0.5D0
                    IF (ELNUMSP(JATOM,ISPIN) .GT. 1D-12) THEN
!                      FVXCK = (VLIN1+ACORR)/ELNUMSP(JATOM,ISPIN)
                      FVXCK = (VLIN1)/ELNUMSP(JATOM,ISPIN)
                    ELSE
                      FVXCK = 0d0
                    ENDIF

                    WRITE(6,652) JATOM, LM1, (LM(I,LM1,JATOM),I=1,2), FEXCK, FVXCK
                    WRITE(49+ISPIN,653) JATOM, LM1, (LM(I,LM1,JATOM),I=1,2), FEXCK, FVXCK

                  ELSEIF ((L .NE. 9) .AND. (M .EQ. 4)) THEN

                    FEXCK = 0d0
                    FVXCK = 0d0

                    WRITE(6,652) JATOM, LM1, (LM(I,LM1,JATOM),I=1,2), FEXCK, FVXCK
                    WRITE(49+ISPIN,653) JATOM, LM1, (LM(I,LM1,JATOM),I=1,2), FEXCK, FVXCK

                  ELSEIF ((M .EQ. 6) .OR. (M .EQ. 8) .OR. (M .EQ. 10)) THEN

                    FEXCK = 0d0
                    FVXCK = 0d0

                    WRITE(6,652) JATOM, LM1, (LM(I,LM1,JATOM),I=1,2), FEXCK, FVXCK
                    WRITE(49+ISPIN,653) JATOM, LM1, (LM(I,LM1,JATOM),I=1,2), FEXCK, FVXCK

                  ELSEIF ((L .EQ. 4) .OR. (L .EQ. 6) .OR. (L .EQ. 7) .OR. (L .EQ. 9)) THEN

!espxc=Exc integrated over the sphere
                    DO J=1, JRI(JATOM)
                      VALUE(J) = (C_KUB(L,M)*CLMSP(J,LM1,JATOM,ISPIN) + &
                        C_KUB(L,M+4)*CLMSP(J,LM1+1,JATOM,ISPIN))* &
                        (C_KUB(L,M)*EXC(J,LM1,JATOM,ISPIN) + &
                        C_KUB(L,M+4)*EXC(J,LM1+1,JATOM,ISPIN))
                    ENDDO
                    CALL CHARGE(R,DX(JATOM),VALUE,1,JRJ,VLIN1)   ! Was charg2
!                    ACORR = VALUE(1)*R(1)*0.5D0
                    ESPXC = ESPXC + VLIN1 !+ ACORR

!vspxc=Vxc integrated over the sphere
                    DO J=1, JRI(JATOM)
                      VALUE(J) = (C_KUB(L,M)*CLMSP(J,LM1,JATOM,ISPIN) + &
                        C_KUB(L,M+4)*CLMSP(J,LM1+1,JATOM,ISPIN))* &
                        (C_KUB(L,M)*VXC(J,LM1,JATOM,ISPIN) + &
                        C_KUB(L,M+4)*VXC(J,LM1+1,JATOM,ISPIN))
                    ENDDO
                    CALL CHARGE(R,DX(JATOM),VALUE,1,JRJ,VLIN1)   ! Was charg2
!                    ACORR = VALUE(1)*R(1)*0.5D0
                    VSPXC = VSPXC + VLIN1 !+ ACORR

!calculation of <g|exc|g> radial part
                    DO J=1, JRI(JATOM)
                      VALUE(J) = (C_KUB(L,M)*EXC(J,LM1,JATOM,ISPIN) + &
                        C_KUB(L,M+4)*EXC(J,LM1+1,JATOM,ISPIN))*CLMSP(J,1,JATOM,ISPIN)
                    ENDDO
                    CALL CHARGE(R,DX(JATOM),VALUE,1,JRJ,VLIN1)   ! Was charg2
!                    ACORR = VALUE(1)*R(1)*0.5D0
                    IF (ELNUMSP(JATOM,ISPIN) .GT. 1D-12) THEN
!                      FEXCK = (VLIN1+ACORR)/ELNUMSP(JATOM,ISPIN)
                      FEXCK = (VLIN1)/ELNUMSP(JATOM,ISPIN)
                    ELSE
                      FEXCK = 0d0
                    ENDIF

!calculation of <g|vxc|g> radial part
                    DO J=1, JRI(JATOM)
                      VALUE(J) = (C_KUB(L,M)*VXC(J,LM1,JATOM,ISPIN) + &
                        C_KUB(L,M+4)*VXC(J,LM1+1,JATOM,ISPIN))*CLMSP(J,1,JATOM,ISPIN)
                    ENDDO
                    CALL CHARGE(R,DX(JATOM),VALUE,1,JRJ,VLIN1)   ! Was charg2
!                    ACORR = VALUE(1)*R(1)*0.5D0
                    IF (ELNUMSP(JATOM,ISPIN) .GT. 1D-12) THEN
!                      FVXCK = (VLIN1+ACORR)/ELNUMSP(JATOM,ISPIN)
                      FVXCK = (VLIN1)/ELNUMSP(JATOM,ISPIN)
                    ELSE
                      FVXCK = 0d0
                    ENDIF

                    WRITE(6,652) JATOM, LM1, (LM(I,LM1,JATOM),I=1,2), FEXCK, FVXCK
                    WRITE(49+ISPIN,653) JATOM, LM1, (LM(I,LM1,JATOM),I=1,2), FEXCK, FVXCK

                  ENDIF

                ENDDO

                WRITE(6,654) JATOM, ESPXC, VSPXC
                WRITE(49+ISPIN,655) JATOM, ESPXC, VSPXC

              ENDIF

              close(49+ispin)

            endif

          endif
!EECE end of potential and energy output

          DO J=1,JRI(JATOM)
            R(J)=R0(JATOM)*EXP(DX(JATOM)*(J-1))
            VTOT(J,1,JATOM)=VTOT(J,1,JATOM)*R(J)/SQFP
          ENDDO

          if(switch1.ne.' EECE')then
            WRITE(ITAPS,2011) 0,0
            WRITE(ITAPS,2021) ( VTOT(J,1,JATOM), J=1,JRI(JATOM) )
            WRITE(ITAPS,2031)
            WRITE(ITAPS,2030)

            DO LM1=2,LMMAX(JATOM)
              WRITE(ITAPN,2011) LM(1,LM1,JATOM),LM(2,LM1,JATOM)
              WRITE(ITAPN,2021) ( VTOT(J,LM1,JATOM), J=1,JRI(JATOM) )
              WRITE(ITAPN,2031)
            ENDDO

            WRITE(ITAPN,2030)
          endif

        ENDDO !jatom

        if(switch1.ne.' EECE' .and. .not. vsp_only)then
          IF(myid.EQ.icoul) THEN
            WRITE(itapn,*)    '    TOTAL POTENTIAL IN INTERSTITIAL'
            WRITE(itapn,*)
            WRITE(itapn,2061) nkk
            if(inversion)then
              DO j=1,nkk
                WRITE(itapn,2071) (kzz(jx,j),jx=1,3), DBLE(cvalue(j)),0.D0
              ENDDO
            else
              DO j=1,nkk
                WRITE(itapn,2071) (kzz(jx,j),jx=1,3), dCONJG(cvalue(j))
              ENDDO
            endif
          ENDIF
       endif
       CLOSE(itaps)
       CLOSE(itapn)

! Write MGGA potential to files 88/89 (vlapl)
!                            or 91/92 (vtau)
! Note that each MPI task writes its own data (atoms) to the same
! shared file. It seems the parts are written in order of ascending
! task number; hence the interstitial part is written by the last tasknumber.
! I do not know why it works, but this approach has been reliable for years for
! the VNS/VSP output.
! If bugs are found here with other mpi's/filesystems... this may need to be
! changed to use MPI IO (MPI_FILE_OPEN, .. ), or copy all data to task 0 and
! do sequential io using "if (myid.eq.0)"
! jdoumont 29/10/2019
!
       iunit_vlapl = 87 + ispin
       iunit_vtau = 90 + ispin
       
       if(switch1.ne.' EECE' .and. .not. vsp_only ) then
          if (MGGA_LAPL_SWITCH) then
             if (myid.eq.0) then
                open(iunit_vlapl,file=fnavn(iunit_vlapl))
                WRITE(iunit_vlapl,*) '   LAPLACIAN DEPENDENT POTENTIAL OF MT-SPHERES'
                WRITE(iunit_vlapl,*) '   NORM  OF VLAPL(R) = VLAPL(R)'
                WRITE(iunit_vlapl,*)
             else
                open(iunit_vlapl,file=fnavn(iunit_vlapl), position='append')
             endif
             
!           loop over atoms/per MPI task
             do jatom = nstart(myid), nstop(myid)
                WRITE(iunit_vlapl,1990) jatom
                WRITE(iunit_vlapl,2001) lmmax(jatom)

!               loop over LM's of the atom
                DO LM1=1,LMMAX(JATOM)
                   WRITE(iunit_vlapl,2011) LM(1,LM1,JATOM),LM(2,LM1,JATOM)
                   WRITE(iunit_vlapl,2021) ( VLAPL(J,LM1,JATOM,ISPIN), J=1,JRI(JATOM) )
                   WRITE(iunit_vlapl,2031)
                ENDDO
                WRITE(iunit_vlapl,2030)
             end do
             
!            last MPI task writes interstitial
             IF(myid.EQ.icoul) THEN
                WRITE(iunit_vlapl,*)    'LAPLACIAN POTENTIAL IN INTERSTITIAL'
                WRITE(iunit_vlapl,*)
                WRITE(iunit_vlapl,2061) nkk
                if(inversion)then
                   DO j=1,nkk
                      WRITE(iunit_vlapl,2071) (kzz(jx,j),jx=1,3), DBLE(vk_vlapl(j,ispin)),0.D0
                   ENDDO
                else
                   DO j=1,nkk
                      WRITE(iunit_vlapl,2071) (kzz(jx,j),jx=1,3), dCONJG(vk_vlapl(j,ispin))
                   ENDDO
                endif
             ENDIF
             CLOSE(iunit_vlapl)
          else
!               delete the file if potential has no laplacian dependence
!               this is important as the precence of the file triggers inclusion
!               of matrix elements in lapw1.
             if(myid.eq.0)  then
                open(iunit_vlapl,file=fnavn(iunit_vlapl))
                CLOSE(iunit_vlapl, status='DELETE')
             end if
          end if
       
          if (MGGA_TAU_SWITCH) then
             if (myid.eq.0) then
                open(iunit_vtau,file=fnavn(iunit_vtau))
                WRITE(iunit_vtau,*) '   TAU DEPENDENT POTENTIAL OF MT-SPHERES'
                WRITE(iunit_vtau,*) '   NORM  OF VTAU(R) = VTAU(R)'
                WRITE(iunit_vtau,*)
             else
                open(iunit_vtau,file=fnavn(iunit_vtau), position='append')
             endif
             
!           loop over atoms/per MPI task
             do jatom = nstart(myid), nstop(myid)
                WRITE(iunit_vtau,1990) jatom
                WRITE(iunit_vtau,2001) lmmax(jatom)

!               loop over LM's of the atom
                DO LM1=1,LMMAX(JATOM)
                   WRITE(iunit_vtau,2011) LM(1,LM1,JATOM),LM(2,LM1,JATOM)
                   WRITE(iunit_vtau,2021) ( VTAU(J,LM1,JATOM,ISPIN), J=1,JRI(JATOM) )
                   WRITE(iunit_vtau,2031)
                ENDDO
                WRITE(iunit_vtau,2030)
             end do
             
!            last MPI task writes interstitial
             IF(myid.EQ.icoul) THEN
                WRITE(iunit_vtau,*)    'TAUACIAN POTENTIAL IN INTERSTITIAL'
                WRITE(iunit_vtau,*)
                WRITE(iunit_vtau,2061) nkk
                if(inversion)then
                   DO j=1,nkk
                      WRITE(iunit_vtau,2071) (kzz(jx,j),jx=1,3), DBLE(vk_vtau(j,ispin)),0.D0
                   ENDDO
                else
                   DO j=1,nkk
                      WRITE(iunit_vtau,2071) (kzz(jx,j),jx=1,3), dCONJG(vk_vtau(j,ispin))
                   ENDDO
                endif
             ENDIF
             CLOSE(iunit_vtau)
          else
!               delete the file if potential has no tau dependence
!               this is important as the precence of the file triggers inclusion
!               of matrix elements in lapw1.
             if(myid.eq.0) then
                open(iunit_vtau,file=fnavn(iunit_vtau))
                CLOSE(iunit_vtau, status='DELETE')
             end if
          end if
       endif
       
#ifdef Parallel
        IF(myid.LT.icoul) THEN
          CALL MPI_SEND(i,1,MPI_INTEGER,myid+1, 0,comm_coul,ierr)
        ENDIF
#endif

      ENDDO SPIN

#ifdef LIBXC
      if(VX_SWITCH.eq.'VX_LIBXC' .or. VC_SWITCH.eq.'VC_LIBXC' .or. &
           EX_SWITCH.eq.'EX_LIBXC' .or. EC_SWITCH.eq.'EC_LIBXC') then
         call libxc_end
      endif
#endif
      if(allocated(vkvxc)) deallocate (vkvxc,vkexc)
      if(allocated(vk_vtau)) deallocate (vk_vtau)
      if(allocated(vk_vlapl)) deallocate (vk_vlapl)
!
! calculate used cpu-time
!
!      call walltim(time6w)
      CALL CPUTIM(time6)
      CALL wallTIM(time6w)
      timtot=time6-time0
      tmulti=time1-time0
      pmulti=tmulti/timtot*100.
      tcouli=timey-time1
      pcouli=tcouli/timtot*100.
      tcoul=time2-time1-tcouli-t947
      pcoul=tcoul/timtot*100.
      txcsph=time3-time2
      pxcsph=txcsph/timtot*100.
      timtotw=time6w-time0w
      tmultiw=time1w-time0w
      pmultiw=tmultiw/timtotw*100.
      tcouliw=timeyw-time1w
      pcouliw=tcouliw/timtotw*100.
      tcoulw=time2w-time1w-tcouliw-t947w
      pcoulw=tcoulw/timtotw*100.
      txcsphw=time3w-time2w
      pxcsphw=txcsphw/timtotw*100.
      if(switch1.ne.' EECE' .and. .not. vsp_only)then
        txcpw=time4-time3
        pxcpw=txcpw/timtot*100.
        tenerg=time5-time4
        penerg=tenerg/timtot*100.
        treana=time6-time5
        preana=treana/timtot*100.
        txcpww=time4w-time3w
        pxcpww=txcpww/timtotw*100.
        tenergw=time5w-time4w
        penergw=tenergw/timtotw*100.
        treanaw=time6w-time5w
        preanaw=treanaw/timtotw*100.
      endif
!      write(6,*) 'total walltime',time6w-time0w
      WRITE(6,3000) timtot,100.,timtotw,100.
!      write(6,*) 'multi walltime (reading clmsum)',time1w-time0w
      WRITE(6,3010) tmulti,pmulti,tmultiw,pmultiw
      WRITE(6,3021) tcouli,pcouli,tcouliw,pcouliw
      WRITE(6,3022) t947,t947/timtot*100.,t947w,t947w/timtotw*100.
      WRITE(6,3020) tcoul,pcoul,tcoulw,pcoulw
      WRITE(6,3030) txcsph,pxcsph,txcsphw,pxcsphw
      if(switch1.ne.' EECE' .and. .not. vsp_only)then
        WRITE(6,3040) txcpw,pxcpw,txcpww,pxcpww
        WRITE(6,3050) tenerg,penerg,tenergw,penergw
        WRITE(6,3070) time4a-time4,time4b-time4a,time4aw-time4w,time4bw-time4aw
        WRITE(6,3060) treana,preana,treanaw,preanaw
        IF(SW_STRESS .eq. 'YES') then
          WRITE(6,3666)  time_electro, time_electro/timtotw*100
        ENDIF
      endif      
      IF(myid.EQ.0) THEN
        if((VX_SWITCH.ne.'VX_RS').and.(VX_SWITCH.ne.'VX_S').and. &
           (VX_SWITCH.ne.'VX_LAPRHO').and.(VX_SWITCH.ne.'VX_TAU').and. &
           (VX_SWITCH.ne.'VX_TAUTF').and.(VX_SWITCH.ne.'VX_TAUW').and. &
           (VX_SWITCH.ne.'VX_TAU-TAUW').and.(VX_SWITCH.ne.'VX_Z').and. &
           (VX_SWITCH.ne.'VX_ALPHA').and.(VX_SWITCH.ne.'VX_ELF')) then
          CALL errclr(errfn)
        else
          CALL OUTERR('LAPW0','finishes okay, but SCF not allowed')
        endif
      ENDIF

 889  CONTINUE

      call cleanup_ffts()

#ifdef Parallel
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call MPI_FINALIZE(ierr)
#endif

      if (myid.eq.0) then
            STOP ' LAPW0 END'
      else
            STOP 
      endif

!
! error handling
!
 910  INFO = 1
!
! lapw0.def couldnt be opened
!
      WRITE (ERRMSG,9000) FNAME
      CALL OUTERR('LAPW0',ERRMSG)
      GOTO 999
 920  INFO = 2
!
! file FNAME couldn''t be opened
!
      WRITE (ERRMSG,9010) IUNIT
      CALL OUTERR('LAPW0',ERRMSG)
      WRITE (ERRMSG,9020) FNAME
      CALL OUTERR('LAPW0',ERRMSG)
      WRITE (ERRMSG,9030) STATUS, FORM
      CALL OUTERR('LAPW0',ERRMSG)
      GOTO 999
 930  INFO = 3
!
! illegal number of equivalent atoms
!
      CALL OUTERR('LAPW0','MULT .EQ. 0')
      GOTO 999
 960  INFO = 7
!
!        Error reading file 'lapw2.def'
!
      WRITE (ERRMSG,9040) FNAME
      CALL OUTERR('LAPW0',ERRMSG)
      GOTO 999
 970  INFO = 8
      WRITE (ERRMSG,9050) JATOM,LMMAX1(JATOM),NCOM+3
      CALL OUTERR('LAPW0',ERRMSG)
      GOTO 999
 971  INFO = 8
      WRITE (ERRMSG,9051) JATOM, l,lmax2
      CALL OUTERR('LAPW0',ERRMSG)
      GOTO 999
 901  INFO = 11
      CALL OUTERR('LAPW0','IBES.GT.20')
      GOTO 999

 999  write(0,*) myid, ' LAPW0-Error. Check file lapw0.error and case.output0*.'
      call abort_parallel
      STOP 'LAPW0 - Error. Check file lapw0.error.'

  77  FORMAT(1X,9(F9.7))
  78  FORMAT(1X,A20,A4,3F10.6,I5)
 699  FORMAT(':FIT',i3.3,':',1x,'LUSE=',i3,' SIGMA OF V-XC FIT FOR ATOM',i3,e15.7)
 700  FORMAT(I3,A77)
 702  FORMAT(' IFFT PARAMETERS:',3I3)
 705  FORMAT(/,/,12X,9('-'),/,':ITE',i3.3,':',I3,'. ITERATION',/,12X,9('-'),/)
 710  FORMAT (5x,'HELLMANN-FEYNMAN-FORCE IN mRy/a.u. = |F|',2x,'Fx',7x,'Fy',13x,'Fz')
 711  FORMAT (':FHF',i3.3,':',1x,i3,'.ATOM',4f20.12,/)
 787  FORMAT('     TOTAL SPHERICAL POTENTIAL IN MT SPHERES',5X,I3,'.ITERATION')
 800  FORMAT(/,/,/,/,30X,50(1H-),/,33X,'S T R U C T U R A L   ','I N F O R M A T I O N',/,30X,50(1H-),/,/)
 805  FORMAT(7X,'SUBSTANCE: ',A60,/)
 810  FORMAT(7X,'LATTICE',22X,'= ',A4)
 819  FORMAT(':POT  : POTENTIAL OPTION ',A,x,A,x,A,x,A)
 820  FORMAT(':LAT  :',1x,'LATTICE CONSTANTS=',3F9.5,3F9.3)
 821  FORMAT(':VOL  :',1x,'UNIT CELL VOLUME =',F14.5)
 823  FORMAT(':POT  : POTENTIAL OPTION ',A,x,A,x,A,x,A,' (with xc1 =',F13.8,')')
 824  FORMAT(':POT  : POTENTIAL OPTION ',A,x,A,x,A,x,A,' (with xc1, xc2 =',2F13.8,')')
 825  FORMAT(':POT  : POTENTIAL OPTION ',A,x,A,x,A,x,A,' (with xc1, xc2, xc3 =',3F13.8,')')
 826  FORMAT(':GKS  : POTENTIAL IS A GKS MGGA POTENTIAL')
 840  FORMAT(':NATO :',I5,' INDEPENDENT AND',I5,' TOTAL ATOMS IN UNITCELL')
 850  FORMAT(7X,'MODE OF CALCULATION IS',7X,'= ',A4)
 860  FORMAT(7X,'SELFCONSISTENT CYCLE-NUMBER  = ',I3,/)
 870  FORMAT(1H1,/,30X,30(1H-),/,30X,'M U L T I P O L M O M E N T S',/,30X,30(1H-))
 880  FORMAT(1H1,/,30X,40(1H-),/,33X,'C O N V E R G E N C E   PARA','METERS',/,30X,40(1H-),/,/)
 991  FORMAT(2E15.7,3I5)
 994  FORMAT(3X,'L=',I3,3X,'M=',I2,3X,'PSEUDO MULTIPOLMOMENT     = ',2F13.9)
 996  FORMAT(3X,'L=',I3,3X,'M=',I2,3X,'PLANE WAVE MULTIPOLMOMENT = ',2F13.9)
 997  FORMAT(A4,2A5)
1011  FORMAT(A10,5X,I5,5X,F10.9,5X,F10.5,5X,F10.5)
1022  FORMAT(5E14.7)
1120  FORMAT(':GRR  :',' <AVERAGE GRAD RHO / RHO>  =',1x,F22.8)           
1121  FORMAT(':DEN  :',' DENSITY INTEGRAL  =',1x,F22.8,3X,'(Ry)')           
1122  FORMAT(':INTRV:',' RHO*V_eff INTEGRAL=',1x,F22.8,3X,'(Ry)')
1123  FORMAT(':EXC  :',' XC ENERGY         =',1x,F22.8,3X,'(Ry)')
1124  FORMAT(':AEXSL:',' AX*EX             =',1x,F22.8,3X,'(Ry)')
1125  FORMAT(':DELTAXC:',' XC DISCONTINUITY=',1x,F22.8,3X,'(Ry)')
1333  FORMAT(2I5,A10)
1334  FORMAT(6E12.5)
1335  FORMAT(':VCOUL',I3.3,1X,'ATOMNUMBER=',I3,1X,A7,1X,'VCOUL-ZERO = ',E12.5)
1340  FORMAT(':EFG',i3.3,':',24X,'EFG         =',F12.5,'   *10**21  V / m**2')         
1341  FORMAT(31X,'V20  TOT/SRF=',2F12.5)                                 
1342  FORMAT(31X,'V22  TOT/SRF=',2F12.5)                               
1343  FORMAT(31X,'V22M TOT/SRF=',2F12.5)                               
1344  FORMAT(31X,'V21  TOT/SRF=',2F12.5)                               
1345  FORMAT(31X,'V21M TOT/SRF=',2F12.5,/)                               
1350  FORMAT(10X,3F11.5,4X,3F11.5)                                            
1355  FORMAT(/,9X,'MAIN DIRECTIONS OF THE EFG ',3F8.4,/,36X,3F8.4,/,36X,3F8.4,/,&
        & ':ANG',i3.3,':',2X,'ANGLE WITH OLD X-AXIS =    ',f8.1)
1360  FORMAT(/,':ETA',i3.3,':',25X,'ASYMM. ETA =',F12.5,/,/) 
1520  FORMAT(/,/,3X,'CONVERGENCE PARAMETER FOR PSEUDOCHARGE: NCON= ', &
        I2,/,3X,'MAXIMAL VALUE OF RMT(JATOM)*ABSK(NKK) : RK  =', &
        F8.5,/,/)
1540  FORMAT(3X,'XC-potentials inside spheres (XCPOT1)')
1550  FORMAT(/,7X,'VALUE OF V0-EXCHANGE',11X,'V0 = ',2E13.5, &
        / ,7X,'VALUE OF V -BAR     ',8X,'V-BAR = ',2E13.5)
1560  FORMAT(7X,'VALUE OF MUFFIN-TIN ZERO',10X,'= ',2E13.5)
1040  FORMAT(8I2)
1060  FORMAT(/,/,3X,'NOT EQUIV ATOM ',A10,'  LOCAL ROTATION MATRIX')
1070  FORMAT(30X,3F10.5)
1080  FORMAT(13X,'EQUIV ATOM ',I3,3X,'POSITION: ',3F8.3)
1970  FORMAT(3X,'             exchange correlation potential ', &
        'BY',I3,'. ITERATION (NORM: CLM=CLM*R*R)')
1980  FORMAT(3X)
1990  FORMAT(3X,'ATOMNUMBER =',I3,5X,10A4)
2000  FORMAT(15X,I3,/,/)
2001  FORMAT(3X,'NUMBER OF LM',I3,/,/)
2010  FORMAT(15X,I3,5X,I2,/)
2011  FORMAT(3X,'VLM(R) FOR L',I3,3X,'M=',I2,/)
2020  FORMAT(3X,5E13.7)
2021  FORMAT(3X,4ES19.12)
2030  FORMAT(/,/,/)
2031  FORMAT(/)
2040  FORMAT(/,/,3X,'ATOM=',I3,3X,A10,3X,'Z=',F5.2,3X, &
        'LM=',I3,5X,'POSITION=',3F7.3,/)
2043  FORMAT(8X,I4,/)
2044  FORMAT(49X,I3,/)
2050  FORMAT(3X,'L=',I3,3X,'M=',I2,3X, &
        'SPHERE MM = ',2F13.9,3x,2f6.2)
2061  FORMAT(1X,'        ',I10,1X,'NUMBER OF PW')
2060  FORMAT(/,13X,I6)
2070  FORMAT(3X,3I5,2E15.7)
2071  FORMAT(3X,3I5,2ES19.12)
2091  FORMAT(':DBG_STR_INT_COUL_JAC:',1X,ES23.15,1X,L1)
3000  FORMAT(/,/,3X,'=====>>> CPU-TIME SUMMARY',/, &
             12X,'TOTAL CPU/WALL-TIME USED  : ',2(F7.1,5X,F4.0,' PERCENT'))
3010  FORMAT(12X,'TIME MULTIPOLMOMENTS      :      ',2(F7.1,5X,F4.0,' PERCENT'))
3020  FORMAT(12X,'TIME COULOMB POT SPH      :      ',2(F7.1,5X,F4.0,' PERCENT'))
3021  FORMAT(12X,'TIME COULOMB POT INT      :      ',2(F7.1,5X,F4.0,' PERCENT'))
3022  FORMAT(12X,'TIME COULOMB POT RMT      :      ',2(F7.1,5X,F4.0,' PERCENT'))
3030  FORMAT(12X,'TIME XCPOT SPHERES        :      ',2(F7.1,5X,F4.0,' PERCENT'))
3040  FORMAT(12X,'TIME XCPOT INTERST        :      ',2(F7.1,5X,F4.0,' PERCENT'))
3050  FORMAT(12X,'TIME TOTAL ENERGY         :      ',2(F7.1,5X,F4.0,' PERCENT'))
3060  FORMAT(12X,'TIME REANALYSE            :      ',2(F7.1,5X,F4.0,' PERCENT'))
3070  FORMAT(12X,'TIME REAN0, REAN3         :      ',2(F7.1,5X,F4.0,' PERCENT'))
3666  FORMAT(12X,'TIME ELECTROSTATIC STRESS :      ',2(F7.1,5X,F4.0,' PERCENT'))
9000  FORMAT('can''t open definition file ',A40)
9010  FORMAT('can''t open unit: ',I2)
9020  FORMAT('       filename: ',A50)
9030  FORMAT('         status: ',A,'  form: ',A)
9040  FORMAT('Error reading file: ',A47)
9050  FORMAT('lmmax1 for jatom',i3,' of ',i4,' gt ',i3)
9051  FORMAT('l-value for jatom',i3,' =',i3,' gt lmax2=',i2)
9060  FORMAT('NWAV TOO SMALL COMPARED WITH NKK',2i7)
9070  FORMAT('NEEDED VALUES:',3I4)
9080  FORMAT('PRESENT PARAMETERS:',3I4)

END
