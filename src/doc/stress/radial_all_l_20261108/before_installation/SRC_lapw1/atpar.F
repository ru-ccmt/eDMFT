      SUBROUTINE ATPAR(NT,REL,NAT,LNSMAX,ZZ,force,efermi)
!
      use reallocate, only: doreallocate
      use loabc, only : ALO, BLO, CLO, ELO, PLO, DPLO, PELO, DPELO, PEILO, &
                        PI12LO, PE12LO, pilolo,init_loabc,nloat_new,LOSECDER
      use radial, only: A,B,AP,BP,AE,BE,VLM,A1,B1,AE1,BE1,A1LO,B1LO, &
                        init_radial_mgga, APE, BPE, AP1, BP1, APE1, BPE1,&
                        AP1LO, BP1LO, VTAU, init_radial, end_radial
      use lolog, only : nlo, loor, lapw, ilo, init_lolog
      use atspdt, only: INS, NL, LMMAX, LQIND, LM, LQNS, LQMM, DP, DPE, E, P, PE, &
           PEI, VNS1, VNS2, VNS3, WNS1, WNS2, WNS3, WNS4, GFAC, init_atspdt, NGAU, &
           INT_U___VTAU_U__,INT_U___VTAU_UD_, INT_UD__VTAU_UD_, INT_U___VTAU_UP_,& 
           INT_U___VTAU_UPD, INT_UP__VTAU_UD_, INT_UD__VTAU_UPD,&
           INT_UP__VTAU_UP_,INT_UP__VTAU_UPD, INT_UPD_VTAU_UPD,&
           init_atspdt_mgga, calculate_mgga_fac, atom_gaunt,&
           mgga_fac
      use char, only  : LATTIC, NAME
      use loint, only : VNS1LO, VNS2LO, VNS3LO, WNS1LO, WNS2LO, WNS3LO, &
           WNS4LO, WNS5LO, &
           init_loint, &
           init_loint_mgga,INT_Ulo__VTAU_U__,INT_Ulo__VTAU_UD_,&
           INT_Ulo__VTAU_UP_,INT_Ulo__VTAU_UPD,INT_UloP_VTAU_U__,&
           INT_UloP_VTAU_UD_,INT_UloP_VTAU_UP_,INT_UloP_VTAU_UPD,&
           INT_Ulo__VTAU_Ulo,INT_UloP_VTAU_Ulo,INT_UloP_VTAU_UloP
      use potnlc, only: JRJ, DH, R0, VR
      use rotmat, only: ROTIJ, ROTLOC
      use struk, only : POS, ALAT, ALPHA, RMT, V, PIA, VI, IATNR, MULT
      use parallel, only: MYID, abort_parallel
      use comi, only: NSLMAX
      use coml, only: mgga, mgga_central_correction
      use orb, only     :  nlorb,LORB,iatom,natorb
      IMPLICIT NONE
      INCLUDE 'param.inc'
#ifdef Parallel
      include 'mpif.h'
#endif
!
!        Arguments
!
      INTEGER            LNSMAX, NAT, NT, IERR, L_next,iapw_next
      DOUBLE PRECISION   ZZ(NAT), TRAPLEVEL, efermi, EI_next,DE_next,EMAIN_next
      parameter (traplevel=1D-4) 
      LOGICAL            REL, DOTRAP
!
!     ..................................................................
!
!        ATPAR calculates the solutions u(l) of the radial Schroedinger
!        Equation, the energy derivatives ue(l), the radial non muffin-
!        tin integrals uvu, uevu, and the Gaunt-coefficients.
!        Spin-orbitless relativistic equations are used.
!
!     ..................................................................
!
!        Local Scalars
!
      INTEGER            I, ICOUNT, IDUMMY, INDEX, J, JATOM, JC, JCOL,jat
      INTEGER            JRI, JROW, L, L0, L0BEFO, LALT, LATOM, LL
      INTEGER            LLBEFO, LLL, LLMM, LM1, LM2, LMBEFO, LMX, LP
      INTEGER            LPBEFO, LQX, M, MINU, MM, MP, NLR, NODEL
      INTEGER            NODES, NODEU, IAPW, JLO, LMMIN
      DOUBLE PRECISION   AAINT, CFEIN1, CFEIN2, CROSS, DE, DE1, DELE, DELEI, DELES
      DOUBLE PRECISION   DUV, DUVB, DUVE, DX, E1, EI, FL, OVLP, RNOT
      DOUBLE PRECISION   TRX, TRY, UV, UVB, UVE,  RNET(NRAD)
      DOUBLE PRECISION   AW(NRAD),BW(NRAD),AEW(NRAD),BEW(NRAD)
!!      COMPLEX*16         IMAG, IMAG1
      CHARACTER*4        EMAIN
      CHARACTER*67       ERRMSG
      integer            jlo1,klo,klo1
      INTEGER            MGGA_GAUNTS, MCOUNT
      logical, dimension(0:4) ::  logic_mgga_gaunt
      character*80  ::   fname
      character*80  ::   u_print_fname_format
      character*32  ::   g2_diag_env
      integer       ::   g2_diag_status,g2_bad_count
      logical       ::   g2_diag_enabled,ham_gaunt_allowed
      character*512 ::   basis_vsp_file
      integer       ::   basis_vsp_status,basis_vsp_ios,basis_vsp_unit
      logical       ::   basis_vsp_enabled,basis_vsp_exists
      real*8, allocatable :: VRB(:,:)
      real*8 :: basis_vr_solve(NRAD)
      character(len=16) :: basis_low_l_env
      integer :: basis_low_l_status
      logical :: basis_low_l_only
      real*8, parameter ::    PI = 4.D0*DATAN(1.D0)
!fb
!     variables for [yu91]-force formalism
      logical     force
!fe
!
!        External Functions
!
      DOUBLE PRECISION   GAUNT
      EXTERNAL           GAUNT
!
!        External Subroutines
!
      EXTERNAL           ABC,CBCOMB,FORFHS,OUTERR,OUTWINb,RINT13,SELECT
!
!        Intrinsic Functions
!
      INTRINSIC          ABS, DBLE, MOD, SQRT
!
!        Data statements
!
      DATA DELE /2.0D-3/
!!      DATA IMAG /(0.0D+0,1.0D+0)/
!
!fb
!      force=.false.
!fe
      NGAU=2350
      LMMIN = 1
      basis_low_l_env=''
      call get_environment_variable('WIEN_LAPW1_BASIS_LOW_L_ONLY', &
           basis_low_l_env,status=basis_low_l_status)
      basis_low_l_only=(basis_low_l_status.eq.0.and.trim(basis_low_l_env).eq.'1')
      nloat=nloat_new
      mgga_gaunts = -1
      g2_diag_env=''
      call get_environment_variable('WIEN_STRESS_G2_DIAG',g2_diag_env, &
                                    status=g2_diag_status)
      g2_diag_enabled=(g2_diag_status.eq.0 .and. &
                       len_trim(g2_diag_env).gt.0 .and. &
                       trim(g2_diag_env).ne.'0')
      call init_loabc(LOMAX, NLOAT, NAT)
      call init_radial(NRAD, LOMAX, NSLMAX, NLOAT, LMMX)
      call init_radial_mgga(NRAD, LOMAX, NSLMAX, NLOAT, LMMX)
      call init_lolog(LOMAX, NAT, LMAX)
      call init_atspdt(NAT, LMMX, LMAX, NSLMAX)
      call init_loint(LOMAX+1, LMMX, NSLMAX, NAT, nloat)
      allocate(VRB(NRAD,NAT))
!        initialize constants
!        CFEIN1 and CFEIN2 are used in RINT13 to calculate the
!        integrals of radialfunctions and are derivates of the
!        Feinstruktur-constant (in Hartrees)
!
      IF (REL) THEN
         CFEIN1 = 1.0D+0
         CFEIN2 = 1.d0/137.0359895d0**2
      ELSE
         CFEIN1 = 1.0D+0
         CFEIN2 = 1.0D-22
      ENDIF
!
!        read total spherical potential V(0,0) of TAPE18=VSP
!        norm of V(0,0)=V00(R)*R/SQRT(4.0D+0*PI)
!
      if(myid.eq.0)  READ (18,5070)
      DO 20 JATOM = 1, NAT
      if(myid.eq.0) then
         READ (18,5010)
         READ (18,5020) IDUMMY
         READ (18,5060)
         READ (18,5040) (VR(J,JATOM),J=1,JRJ(JATOM))
         READ (18,5060)
         READ (18,5050)
         DO J = 1, JRJ(JATOM)
            VR(J,JATOM) = VR(J,JATOM)/2.0D+0
         ENDDO
      endif
#ifdef Parallel
           CALL MPI_BCAST(vr(1,jatom), JRJ(JATOM) , MPI_DOUBLE_PRECISION, 0,  &
                     MPI_COMM_WORLD, ierr)
#endif
   20 CONTINUE

!     Optional radial-basis response diagnostic.
!
!     In ordinary LAPW1 the same spherical potential generates u_l,
!     udot_l and all local orbitals and is also the spherical part of the
!     physical Hamiltonian.  This hides the distinction between a change of
!     the radial trial space and a change of the operator.
!
!     If WIEN_LAPW1_BASIS_VSP names another VSP file, only the radial
!     functions are solved with that potential (VRB).  The regular unit-18
!     VSP (VR) remains the physical potential.  Farther below, VR-VRB is
!     added back as an L=0 matrix element.  Thus the generalized eigenproblem
!     is the center physical Hamiltonian represented in the perturbed radial
!     basis.  Its central-difference eigenvalue derivative directly measures
!
!          2 Re < d_phi/d_epsilon | H-epsilon | phi >,
!
!     including the matching-coefficient and local-orbital response.  This
!     path is strictly opt-in: without the environment variable VRB=VR.
      VRB=VR
      basis_vsp_file=''
      call get_environment_variable('WIEN_LAPW1_BASIS_VSP', &
           basis_vsp_file,status=basis_vsp_status)
      basis_vsp_enabled=(basis_vsp_status.eq.0 .and. &
           len_trim(basis_vsp_file).gt.0)
      basis_vsp_ios=0
      basis_vsp_unit=98
      if (basis_vsp_enabled.and.myid.eq.0) then
         inquire(file=trim(basis_vsp_file),exist=basis_vsp_exists)
         if (.not.basis_vsp_exists) then
            basis_vsp_ios=-1
         else
            open(unit=basis_vsp_unit,file=trim(basis_vsp_file), &
                 status='old',form='formatted',action='read', &
                 iostat=basis_vsp_ios)
         endif
         if (basis_vsp_ios.eq.0) then
            read(basis_vsp_unit,5070,iostat=basis_vsp_ios)
            do jat=1,nat
               if (basis_vsp_ios.ne.0) exit
               read(basis_vsp_unit,5010,iostat=basis_vsp_ios)
               if (basis_vsp_ios.eq.0) &
                    read(basis_vsp_unit,5020,iostat=basis_vsp_ios) idummy
               if (basis_vsp_ios.eq.0) &
                    read(basis_vsp_unit,5060,iostat=basis_vsp_ios)
               if (basis_vsp_ios.eq.0) &
                    read(basis_vsp_unit,5040,iostat=basis_vsp_ios) &
                    (VRB(j,jat),j=1,jrj(jat))
               if (basis_vsp_ios.eq.0) &
                    read(basis_vsp_unit,5060,iostat=basis_vsp_ios)
               if (basis_vsp_ios.eq.0) &
                    read(basis_vsp_unit,5050,iostat=basis_vsp_ios)
               if (basis_vsp_ios.eq.0) &
                    VRB(1:jrj(jat),jat)=VRB(1:jrj(jat),jat)/2.0d0
            enddo
            close(basis_vsp_unit)
         endif
      endif
#ifdef Parallel
      call MPI_BCAST(basis_vsp_ios,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif
      if (basis_vsp_enabled.and.basis_vsp_ios.ne.0) then
         if (myid.eq.0) write(6,*) &
              'ERROR reading WIEN_LAPW1_BASIS_VSP=',trim(basis_vsp_file)
         call OUTERR('ATPAR','cannot read alternate radial-basis VSP')
         call abort_parallel
         stop 'ATPAR - alternate basis VSP error'
      endif
#ifdef Parallel
      if (basis_vsp_enabled) then
         do jat=1,nat
            call MPI_BCAST(VRB(1,jat),JRJ(jat),MPI_DOUBLE_PRECISION,0, &
                 MPI_COMM_WORLD,ierr)
         enddo
      endif
#endif
      if (basis_vsp_enabled.and.myid.eq.0) write(6,*) &
           ':BASIS_VSP: radial basis from ',trim(basis_vsp_file), &
           '; physical Hamiltonian remains unit-18 VSP'
!
     if (mgga) then
        if(myid.eq.0) then
           read (15,5050,IOSTAT=MGGA_GAUNTS)
           rewind 15
           READ (15,5070,IOSTAT=MGGA_GAUNTS)
        endif
#ifdef Parallel
           CALL MPI_BCAST(MGGA_GAUNTS, 1, MPI_INTEGER, 0, &
                      MPI_COMM_WORLD, ierr)
#endif
     endif
                
     if(myid.eq.0) then
        READ (19,5050,END=30,IOSTAT=INS)
        REWIND 19
        READ (19,5070)
!fb
!     force-output only 
!     if non-spherical potential exists
!      IF(myid.eq.0) then
        if(force)  write(71) &
          &" non-spherical Matrixelements x Gaunt-Coefficient (d=du/dE):", &
          &" l0 ll lp m0 mm mp     uVu dVd", &
          &"                       uVd dVu"
!fe
   30 CONTINUE
     endif

     if (MGGA_GAUNTS .EQ. 0) then
        mgga = .true.
        call init_atspdt_mgga(NAT, LMMX, LMAX, NSLMAX)
        call init_loint_mgga(LOMAX+1, LMMX, NSLMAX, NAT, nloat)
     else
         mgga = .false.
         mgga_central_correction = .false.
     end if
     if (basis_vsp_enabled.and.mgga) then
        if (myid.eq.0) write(6,*) &
             'WIEN_LAPW1_BASIS_VSP is not implemented for MGGA'
        call OUTERR('ATPAR','alternate radial-basis VSP with MGGA')
        call abort_parallel
        stop 'ATPAR - alternate basis VSP/MGGA error'
     endif
!
!        start loop over atoms in unitcell
!
      LATOM = 0
      INDEX = 0
      NLO = 0
      DO 340 JATOM = 1, NAT
!re-initialize radial functions
         A=0.0D0
         AP=0.0D0
         A1=0.0D0
         A1LO=0.0D0
         AE=0.0D0
         AE1=0.0D0
         B=0.0D0
         BP=0.0D0
         B1=0.0D0
         B1LO=0.0D0
         BE=0.0D0
         BE1=0.0D0
         VLM=0.0D0
         if (mgga) then
            APE = 0.0D0
            BPE = 0.0D0
            AP1 = 0.0D0
            BP1 = 0.0D0
            APE1 = 0.0D0
            BPE1 = 0.0D0
            AP1LO = 0.0D0
            BP1LO = 0.0D0
            VTAU = 0.0D0
         end if

         do l=0,lmax-1
            lapw(l,jatom)=.true.
         enddo
         RNOT = R0(JATOM)
         LATOM = LATOM + 1
         JRI = JRJ(JATOM)
         DX = DH(JATOM)
!        Setup RMT mesh. Moved out of loop for efficiency
         do i=1,NRAD
           RNET(i)=RNOT*(exp(DX*(i-1.d0)))
         enddo
         iapw=0
         READ (5,*,err=41) EI, NLR,iapw
         if(abs(EI-0.3d0).lt.1.d-6) EI=efermi - 0.2d0
         if(iapw.eq.1) then
            do l=0,lmax-1
               lapw(l,jatom)=.false.
            enddo
         endif
 41      continue
!         WRITE (6,6040) IATNR(JATOM), (POS(I,LATOM),I=1,3),
!     ,        MULT(JATOM)
         LATOM = LATOM + MULT(JATOM) - 1
         IF(myid.eq.0) WRITE (6,6080) jatom,NAME(JATOM)
         DO 40 JCOL = 1, 3
!            WRITE (6,6090) (ROTLOC(JROW,JCOL,JATOM),JROW=1,3)
 40      CONTINUE
         DO 60 M = 1, MULT(JATOM)
            INDEX = INDEX + 1
!            WRITE (6,6100) M, (POS(JC,INDEX),JC=1,3)
            DO 50 JCOL = 1, 3
!               WRITE (6,6110) (ROTIJ(JROW,JCOL,INDEX),JROW=1,3)
 50         CONTINUE
 60      CONTINUE
         if(iapw.eq.1) then
!            WRITE (6,6000) NAME(JATOM), jatom, EI,'  APW'
            IF(myid.eq.0) WRITE (21,6000) NAME(JATOM), jatom, EI,'  APW'
         else
!            WRITE (6,6000) NAME(JATOM), jatom, EI,' LAPW'
            IF(myid.eq.0) WRITE (21,6000) NAME(JATOM), jatom, EI,' LAPW'
         endif
         DO 70 L = 1, LMAX
            IF ((L-1) .LE. (LOMAX)) THEN
               LOOR(L-1,JATOM) = .FALSE.
               ilo(l-1,jatom)=0
if(nloat.gt.3) then
               do i=1,nloat
                  ELO(L-1,i,JATOM) = 99997.0D+0
                  if(lapw(l-1,jatom)) ELO(L-1,i,JATOM) = 99999.0D+0
               enddo
else
               do i=1,nloat
                  ELO(L-1,i,JATOM) = 997.0D+0
                  if(lapw(l-1,jatom)) ELO(L-1,i,JATOM) = 999.0D+0
               enddo
endif
            ENDIF
            E(L,JATOM) = EI
 70      CONTINUE
         L = -1
         DO 80 J = 1, NLR
            LALT = L
            iapw=0
            READ (5,5000,err=42) L, EI, DE, EMAIN, iapw
            if(l.gt.LOMAX) then
!            if(l.gt.LOMAX.and. l.eq.lalt) then ! this would be a good fix, but
!            incompatible with lapw2
              write(*,*) 'l-value of LOs.gt.LOMAX, increase LOMAX in param.inc'
              goto 901
            endif
            if (de.lt.-1.0d-8) de=-de/10000.d0
            L_next=-1
            iapw_next=3
            EI_next=999.
            if(j.ne.NLR) then
               READ (5,5000,err=42) L_next,EI_next,DE_next,EMAIN_next,iapw_next
               backspace (5)
            endif
!            if(abs(EI-0.3d0).lt.1.d-6.and.l.ne.l_next) EI=efermi - 0.2d0
            if(abs(EI-0.3d0).lt.1.d-6.and.l.eq.l_next.and.EI_next.lt.-2.d0) &
               EI=efermi - 0.2d0
            if(abs(EI-0.3d0).lt.1.d-6.and.l.eq.l_next.and.EI_next.gt.-2.d0 &
               .and.EI_next.lt.efermi.and.iapw_next.ne.2) EI=efermi + 0.2d0
            if(abs(EI-0.3d0).lt.1.d-6.and.l.eq.l_next.and.EI_next.gt.-2.d0 &
               .and.EI_next.gt.efermi.and.iapw_next.ne.2) EI=efermi - 0.2d0
            if(abs(EI-0.3d0).lt.1.d-6.and.l.ne.lalt) EI=efermi - 0.2d0
            if(abs(EI-0.3d0).lt.1.d-6.and.l.eq.lalt) EI=efermi + 0.2d0
 42         continue
            if (iapw.eq.2) then
               if (ilo(l,jatom).ge.1) then
                  EI=ELO(L,1,JATOM)
               else
                  EI=E(L,JATOM)
               endif
            endif
            IF (DE .GT. 1.0D-9.and.(iapw.ne.2)) THEN
               CALL SELECT(L,EI,DE,EMAIN,REL, &
                    VRB(1,JATOM),RNOT,DX,JRI,ZZ(JATOM),jatom,efermi)
            ELSE
               IF(myid.eq.0) WRITE (6,6010) l,jatom, L, EI
               IF(myid.eq.0) WRITE (21,6010) l,jatom, L, EI
            ENDIF
            IF (L .EQ. LALT) THEN
               ilo(l,jatom)=ilo(l,jatom)+1
               if (ilo(l,jatom).gt.nloat) then
                  goto 902
               endif
               ELO(L,ilo(l,jatom),JATOM) = EI
               NLO = NLO + ((2*L+1))*MULT(JATOM)
               if (iapw.eq.2) then
                  do jat=1,NATORB
                     if(jatom.eq.iatom(jat)) then
                       if(lorb(1,jat).eq.l) then
                      print*,'hdlo+orbital potential not supported, STOP for atom,l',iatom(jat),l
                      call outerr('ATPAR','case.in1: HDLO + VORB not supported')
                      call abort_parallel
                      STOP 'ATPAR - hdlo+orb Error'
                        endif     
                     endif     
                  enddo   
                  LOSECDER(L,(ilo(l,jatom)),JATOM) = .true.
                  IF(myid.eq.0) write(6,*)  '            LOCAL ORBITAL(SECDER)'
                  IF(myid.eq.0) write(21,*) '            LOCAL ORBITAL(SECDER)'
                  if(lapw(l,jatom)) then
                       print*, "ERROR in case.in1: second derivative LO only with APW, not LAPW"
                       call outerr('ATPAR','case.in1: second derivative LO only with APW')
                       call abort_parallel
                       STOP 'ATPAR - HDLO-LAPW Error'
                  endif     
               else
                  if(ilo(l,jatom).gt.1) then
                     if(LOSECDER(L,(ilo(l,jatom)-1),JATOM)) then
                       print *,"ERROR in case.in1: second derivative LO must be last"
                       call outerr('ATPAR','case.in1: second derivative LO must be last')
                       call abort_parallel
                       STOP 'ATPAR - HDLO-order Error'
                     endif
                  endif
                  IF(myid.eq.0) write(6,*)  '            LOCAL ORBITAL'
                  IF(myid.eq.0) write(21,*) '            LOCAL ORBITAL'
               endif
               loor(l,jatom)=.true.
            ELSE
               if(iapw.eq.1) then
                  lapw(l,jatom)=.false.
                  NLO = NLO + ((2*L+1))*MULT(JATOM)
                  IF(myid.eq.0) write(6,*) '            APW+lo'
                  IF(myid.eq.0) write(21,*) '            APW+lo'
                  ilo(l,jatom)=1
                  elo(l,1,jatom)=EI
if(nloat.gt.3) then
                  elo(l,2,jatom)=99997.0D+0
else
                  elo(l,2,jatom)=997.0D+0
endif
               else
                  IF(myid.eq.0) write(6,*) '            LAPW'
                  IF(myid.eq.0) write(21,*) '            LAPW'
                  lapw(l,jatom)=.true.
               endif
               E(L+1,JATOM) = EI
            ENDIF
 80      CONTINUE
         
         IF(myid.eq.0) WRITE (6,6020)
         DELEI = .250D+0/DELE
!         write(6,*) e(:,jatom)
         DO 150 J = 1, NT
!  Experimental Trap location
            DOTRAP=.true.
1500        CONTINUE
            L = J - 1
            FL = L
            basis_vr_solve(:)=VRB(:,JATOM)
            ! Diagnostic: the alternate-basis L=0 operator restoration is
            ! currently in HNS, whose orbital legs end at LNSMAX.  Freeze
            ! all higher-l basis functions to remove the uncovered leg.
            ! This is a restricted response, not an all-l correction.
            if (basis_low_l_only.and.L.gt.LNSMAX) &
                 basis_vr_solve(:)=VR(:,JATOM)
            EI = E(J,JATOM)/2.D0
!
!        calculate energy-derivative by finite difference
!        DELE is the up and downward energy-shift in Hartrees
!
            E1 = EI - DELE
            CALL OUTWINB(REL,basis_vr_solve,RNeT,DX,JRI,E1,FL,UVB,DUVB, &
                 NODEL,ZZ(JATOM))
            CALL RINT13(CFEIN1,CFEIN2,A,B,A,B,OVLP,JATOM)
            TRX = 1.0D+0/SQRT(OVLP)
            AE(1:JRI) = TRX*A(1:JRI)
            BE(1:JRI) = TRX*B(1:JRI)
            APE(1:JRI) = TRX*AP(1:JRI)
            BPE(1:JRI) = TRX*BP(1:JRI)
            UVB = TRX*UVB
            DUVB = TRX*DUVB

            E1 = EI + DELE
            CALL OUTWINB(REL,basis_vr_solve,RNeT,DX,JRI,E1,FL,UVE,DUVE, &
                 NODEU,ZZ(JATOM))
            CALL RINT13(CFEIN1,CFEIN2,A,B,A,B,OVLP,JATOM)
            TRX = 1.0D+0/SQRT(OVLP)
            UVE = DELEI*(TRX*UVE - UVB)
            DUVE = DELEI*(TRX*DUVE - DUVB)
            AE(1:JRI) = DELEI*(TRX*A(1:JRI) - AE(1:JRI))
            BE(1:JRI) = DELEI*(TRX*B(1:JRI) - BE(1:JRI))
            APE(1:JRI) = DELEI*(TRX*AP(1:JRI) - APE(1:JRI))
            BPE(1:JRI) = DELEI*(TRX*BP(1:JRI) - BPE(1:JRI))
!     
!     calculate function at EI
!
            CALL OUTWINB(REL,basis_vr_solve,RNeT,DX,JRI,EI,FL,UV,DUV,NODES &
                 ,ZZ(JATOM))
            CALL RINT13(CFEIN1,CFEIN2,A,B,A,B,OVLP,JATOM)
            TRX = 1.0D+0/SQRT(OVLP)
            P(J,JATOM) = TRX*UV
            DP(J,JATOM) = TRX*DUV
            A(1:JRI) = TRX*A(1:JRI)
            B(1:JRI) = TRX*B(1:JRI)
            AP(1:JRI) = TRX*AP(1:JRI)
            BP(1:JRI) = TRX*BP(1:JRI)
!
!     insure orthogonalization
!
            CALL RINT13(CFEIN1,CFEIN2,A,B,AE,BE,CROSS,JATOM)
            TRY = -CROSS
            IF (TRY .LT. -0.05D+0) WRITE (6,6070) L, TRY, OVLP
            AE(1:JRI) = AE(1:JRI) + TRY*A(1:JRI)
            BE(1:JRI) = BE(1:JRI) + TRY*B(1:JRI)
!              The radial derivatives are not explicitly normalized
!              The same factor is subtracted to keep them consistent with their
!              primitives
            APE(1:JRI) = APE(1:JRI) + TRY*AP(1:JRI)
            BPE(1:JRI) = BPE(1:JRI) + TRY*BP(1:JRI)
            PE(J,JATOM) = UVE + TRY*P(J,JATOM)
            DPE(J,JATOM) = DUVE + TRY*DP(J,JATOM)
            CALL RINT13(CFEIN1,CFEIN2,AE,BE,AE,BE,PEI(J,JATOM),JATOM)
            IF (J .LE. NSLMAX) THEN
               A1(1:JRI,J) = A(1:JRI)
               B1(1:JRI,J) = B(1:JRI)
               AP1(1:JRI,J) = AP(1:JRI)
               BP1(1:JRI,J) = BP(1:JRI)
               AE1(1:JRI,J) = AE(1:JRI)
               BE1(1:JRI,J) = BE(1:JRI)
               APE1(1:JRI,J) = APE(1:JRI)
               BPE1(1:JRI,J) = BPE(1:JRI)
            ENDIF

! Trap P(J,JATOM) too small
        IF((Abs(P(J,JATOM)) .lt. traplevel).and.DOTRAP)then
                DOTRAP=.false.
! An alternative might be to set the level to LAPW....?
                E(J,JATOM)=E(J,JATOM)-0.01
                IF(myid.eq.0) write(21,1501)P(J,JATOM),E(J,JATOM)
                IF(myid.eq.0) write(6,1501)P(J,JATOM),E(J,JATOM)
1501 format(':INFO  : P(J,JATOM) almost zero (',e15.3,'). Shifted Energy by -0.01 down to ',F10.4)
                goto 1500
        endif
            IF(myid.eq.0) WRITE (6,6030) L, P(J,JATOM), DP(J,JATOM), PE(J,JATOM), &
                 DPE(J,JATOM), PEI(J,JATOM), NODEL, NODES, &
                 NODEU
 150     CONTINUE
!     
!     and now for local orbitals
!     
         IF(myid.eq.0) write(6,6021)
         DO L = 0, LOMAX
            IF (LOOR(L,JATOM)) THEN
               J = L + 1
               FL = L
               jlo=0

               do jlo=1,ilo(l,jatom)

                  if ((.not.lapw(l,jatom)).and.jlo.eq.1) cycle

                  EI = ELO(L,jlo,JATOM)/2.D0
                  if (LOSECDER(L,jlo,JATOM)) then
                     DELES=DELE*2.0d0
!     calculate second derivative at EI (3 point method)
!begin     
!                     CALL OUTWINB(REL,VR(1,JATOM),RNeT,DX,JRI,EI,FL,UV,DUV,NODES &
!                                  ,ZZ(JATOM))
!                     CALL RINT13(CFEIN1,CFEIN2,A,B,A,B,OVLP,JATOM)
!                     TRX = 1.0D+0/SQRT(OVLP)
!                     UV = -2.0d0*TRX*UV
!                     DUV = -2.0d0*TRX*DUV
!                     DO M = 1, JRI
!                        AE(M) = -2.0d0*TRX*A(M)
!                        BE(M) = -2.0d0*TRX*B(M)
!                     ENDDO
!
!                     E1 = EI - DELES
!                     CALL OUTWINB(REL,VR(1,JATOM),RNeT,DX,JRI,E1,FL,UVB,DUVB, &
!                          NODEL,ZZ(JATOM))
!                     CALL RINT13(CFEIN1,CFEIN2,A,B,A,B,OVLP,JATOM)
!                     TRX = 1.0D+0/SQRT(OVLP)
!                     UV = UV+TRX*UVB
!                     DUV = DUV+TRX*DUVB
!                     DO M = 1, JRI
!                        AE(M) = AE(M)+TRX*A(M)
!                        BE(M) = BE(M)+TRX*B(M)
!                     ENDDO
!
!                     E1 = EI + DELES
!                     CALL OUTWINB(REL,VR(1,JATOM),RNeT,DX,JRI,E1,FL,UVE,DUVE, &
!                          NODEU,ZZ(JATOM))
!                     CALL RINT13(CFEIN1,CFEIN2,A,B,A,B,OVLP,JATOM)
!                     TRX = 1.0D+0/SQRT(OVLP)
!                     PLO(L,jlo,JATOM) = (UV+TRX*UVE)/(4.0d0*DELES*DELES)
!                     DPLO(L,jlo,JATOM) = (DUV+TRX*DUVE)/(4.0d0*DELES*DELES)
!                     DO M = 1, JRI
!                        A(M) = (AE(M)+TRX*A(M))/(4.0d0*DELES*DELES)
!                        B(M) = (BE(M)+TRX*B(M))/(4.0d0*DELES*DELES)
!                     ENDDO
!     calculate second derivative at EI (3 point method)
!end     
!     calculate second derivative at EI (5 point method)
!begin
                     CALL OUTWINB(REL,VRB(1,JATOM),RNeT,DX,JRI,EI,FL,UV,DUV,NODES &
                                  ,ZZ(JATOM))
                     CALL RINT13(CFEIN1,CFEIN2,A,B,A,B,OVLP,JATOM)
                     TRX = 1.0D+0/SQRT(OVLP)
                     UV = -30.0d0*TRX*UV
                     DUV = -30.0d0*TRX*DUV
                     AE(1:JRI) = -30.0d0*TRX*A(1:JRI)
                     BE(1:JRI) = -30.0d0*TRX*B(1:JRI)
                     APE(1:JRI) = -30.0d0*TRX*AP(1:JRI)
                     BPE(1:JRI) = -30.0d0*TRX*BP(1:JRI)

                     E1 = EI - 2.0d0*DELES
                     CALL OUTWINB(REL,VRB(1,JATOM),RNeT,DX,JRI,E1,FL,UVB,DUVB, &
                          NODEL,ZZ(JATOM))
                     CALL RINT13(CFEIN1,CFEIN2,A,B,A,B,OVLP,JATOM)
                     TRX = 1.0D+0/SQRT(OVLP)
                     UV = UV-TRX*UVB
                     DUV = DUV-TRX*DUVB
                     AE(1:JRI) = AE(1:JRI)-TRX*A(1:JRI)
                     BE(1:JRI) = BE(1:JRI)-TRX*B(1:JRI)
                     APE(1:JRI) = APE(1:JRI)-TRX*AP(1:JRI)
                     BPE(1:JRI) = BPE(1:JRI)-TRX*BP(1:JRI)

                     E1 = EI - DELES
                     CALL OUTWINB(REL,VRB(1,JATOM),RNeT,DX,JRI,E1,FL,UVB,DUVB, &
                          NODEL,ZZ(JATOM))
                     CALL RINT13(CFEIN1,CFEIN2,A,B,A,B,OVLP,JATOM)
                     TRX = 1.0D+0/SQRT(OVLP)
                     UV = UV+16.0d0*TRX*UVB
                     DUV = DUV+16.0d0*TRX*DUVB
                     AE(1:JRI) = AE(1:JRI)+16.0d0*TRX*A(1:JRI)
                     BE(1:JRI) = BE(1:JRI)+16.0d0*TRX*B(1:JRI)
                     APE(1:JRI) = APE(1:JRI)+16.0d0*TRX*AP(1:JRI)
                     BPE(1:JRI) = BPE(1:JRI)+16.0d0*TRX*BP(1:JRI)

                     E1 = EI + 2.0d0*DELES
                     CALL OUTWINB(REL,VRB(1,JATOM),RNeT,DX,JRI,E1,FL,UVE,DUVE, &
                          NODEU,ZZ(JATOM))
                     CALL RINT13(CFEIN1,CFEIN2,A,B,A,B,OVLP,JATOM)
                     TRX = 1.0D+0/SQRT(OVLP)
                     UV = UV-TRX*UVE
                     DUV = DUV-TRX*DUVE
                     AE(1:JRI) = AE(1:JRI)-TRX*A(1:JRI)
                     BE(1:JRI) = BE(1:JRI)-TRX*B(1:JRI)
                     APE(1:JRI) = APE(1:JRI)-TRX*AP(1:JRI)
                     BPE(1:JRI) = BPE(1:JRI)-TRX*BP(1:JRI)

                     E1 = EI + DELES
                     CALL OUTWINB(REL,VRB(1,JATOM),RNeT,DX,JRI,E1,FL,UVE,DUVE, &
                          NODEU,ZZ(JATOM))
                     CALL RINT13(CFEIN1,CFEIN2,A,B,A,B,OVLP,JATOM)
                     TRX = 1.0D+0/SQRT(OVLP)
                     PLO(L,jlo,JATOM) = (UV+16.0d0*TRX*UVE)/(12.0d0*4.0d0*DELES*DELES)
                     DPLO(L,jlo,JATOM) = (DUV+16.0d0*TRX*DUVE)/(12.0d0*4.0d0*DELES*DELES)
                     A(1:JRI) = (AE(1:JRI)+16.0d0*TRX*A(1:JRI))/(12.0d0*4.0d0*DELES*DELES)
                     B(1:JRI) = (BE(1:JRI)+16.0d0*TRX*B(1:JRI))/(12.0d0*4.0d0*DELES*DELES)
                     AP(1:JRI) = (APE(1:JRI)+16.0d0*TRX*AP(1:JRI))/(12.0d0*4.0d0*DELES*DELES)
                     BP(1:JRI) = (BPE(1:JRI)+16.0d0*TRX*BP(1:JRI))/(12.0d0*4.0d0*DELES*DELES)
!end     
!     calculate second derivative at EI (5 point method)
!
!     insure orthogonalization of second derivative
!
!                     CALL RINT13(CFEIN1,CFEIN2,A1(1,J),B1(1,J),A,B,CROSS,JATOM)
!                     TRY = -CROSS
!                     DO M = 1, JRI
!                        A(M) = A(M) + TRY*A1(M,J)
!                        B(M) = B(M) + TRY*B1(M,J)
!                     ENDDO
!                     PLO(L,jlo,JATOM) = PLO(L,jlo,JATOM) + TRY*P(J,JATOM)
!                     DPLO(L,jlo,JATOM) = DPLO(L,jlo,JATOM) + TRY*DP(J,JATOM)

                  else
!     
!     calculate function at EI
!     
                     CALL OUTWINB(REL,VRB(1,JATOM),RNeT,DX,JRI,EI,FL,UV,DUV,&
                               NODES,ZZ(JATOM))
                     CALL RINT13(CFEIN1,CFEIN2,A,B,A,B,OVLP,JATOM)

                     TRX = 1.0D+0/SQRT(OVLP)
                     PLO(L,jlo,JATOM) = TRX*UV
                     DPLO(L,jlo,JATOM) = TRX*DUV
                     A(1:JRI) = TRX*A(1:JRI)
                     B(1:JRI) = TRX*B(1:JRI)
                     AP(1:JRI) = TRX*AP(1:JRI)
                     BP(1:JRI) = TRX*BP(1:JRI)
                  endif


                  IF (J .LE. NSLMAX) THEN
                     A1LO(1:JRI,jlo,L) = A(1:JRI)
                     B1LO(1:JRI,jlo,L) = B(1:JRI)
                     AP1LO(1:JRI,jlo,L) = AP(1:JRI)
                     BP1LO(1:JRI,jlo,L) = BP(1:JRI)
                  else
                     print*,'nslmax for LOs too small'     
                     goto 903
                  ENDIF

                  CALL RINT13(CFEIN1,CFEIN2,A1(1,J),B1(1,J),A,B,&
                              PI12LO(L,jlo,JATOM),JATOM)
                  CALL RINT13(CFEIN1,CFEIN2,AE1(1,J),BE1(1,J),A,B,&
                              PE12LO(L,jlo,JATOM),JATOM)

                  IF(myid.eq.0) WRITE(6,6031) L,PLO(L,jlo,JATOM),DPLO(L,jlo,JATOM), &
                       PI12LO(L,jlo,JATOM),PE12LO(L,jlo,JATOM), NODEL, NODES, NODEU

               enddo

               do jlo=1,ilo(l,jatom)
                  if ((.not.lapw(l,jatom)).and.jlo.eq.1) cycle
                  do klo=1,ilo(l,jatom)
                     if ((.not.lapw(l,jatom)).and.klo.eq.1) cycle

                     CALL RINT13(CFEIN1,CFEIN2,A1LO(1,jlo,L),B1LO(1,jlo,L),A1LO(1,klo,L),B1LO(1,klo,L),PILOLO(L,jlo,klo,JATOM),JATOM)
                     IF(myid.eq.0) write(6,'(2i5,3f15.6)') jlo,klo,PILOLO(L,jlo,klo,JATOM),ELO(L,jlo,JATOM),ELO(L,klo,JATOM)
                  enddo
               enddo

            endif
         enddo
         do l=0,lomax
            do jlo=1,ilo(l,jatom)
               CALL ABC(JATOM,L,jlo,lapw(l,jatom))
            enddo
         enddo

         IF (INS .EQ. 0) THEN
! fb
            IF(myid.eq.0) then
              if (force) write(71) jatom
!fe
!
!        read total nonspherical potential of TAPE19=VNS
!        norm of VLM=VLM*1
!
              READ (19,5010)
              READ (19,5020) LMMAX(JATOM)
           endif
#ifdef Parallel
           CALL MPI_BCAST(LMMAX(JATOM), 1, MPI_INTEGER, 0,  &
                     MPI_COMM_WORLD, ierr)
#endif
            LM2 = 0
            DO 230 LM1 = 1, LMMAX(JATOM)
               IF(myid.eq.0) READ (19,5030) LLL, MM
#ifdef Parallel
                  CALL MPI_BCAST(LLL, 1, MPI_INTEGER, 0,  &
                     MPI_COMM_WORLD, ierr)
                  CALL MPI_BCAST(MM, 1, MPI_INTEGER, 0,  &
                     MPI_COMM_WORLD, ierr)
#endif
               LL = ABS(LLL)
               LM2 = LM2 + 1
               IF (LM2 .GT. LMMX) THEN
!
!        Error: LMMX too small - stop execution
!
                  GOTO 900
               ENDIF
               IF(myid.eq.0) READ (19,5040) (VLM(J,LM2),J=1,JRJ(JATOM))
#ifdef Parallel
                CALL MPI_BCAST(VLM(1,LM2), JRJ(JATOM), MPI_DOUBLE_PRECISION, 0, &
                     MPI_COMM_WORLD, ierr)
#endif
!     
               IF (LL .GT. LNSMAX*2) THEN
                  LM2 = LM2 - 1
               ELSE
                  LM(1,LM2,JATOM) = LLL
                  LM(2,LM2,JATOM) = MM
               ENDIF
               IF(myid.eq.0) READ (19,5060)
 230        CONTINUE
            IF(myid.eq.0) THEN
              READ (19,5050)
            endif

!           Restore the *physical* spherical Hamiltonian after generating
!           the radial basis with VRB.  VR and VRB store r*V_00/sqrt(4*pi)
!           in Hartree after the factor-of-two conversion above.  VLM is a
!           physical radial potential in Rydberg; hence the conversion
!
!             Delta V_0(r) = 2 [VR(r)-VRB(r)] sqrt(4*pi) / r .
!
!           Index zero is already reserved by the MGGA central-correction
!           implementation for exactly such a spherical matrix element.
!           Adding it here makes H = T+VRB+(VR-VRB) while all matching
!           functions continue to be those generated from VRB.
            if (basis_vsp_enabled) then
               LMMIN=0
               LM(1,0,JATOM)=0
               LM(2,0,JATOM)=0
               VLM(1:JRI,0)=2.0d0*(VR(1:JRI,JATOM)- &
                    VRB(1:JRI,JATOM))*SQRT(4.0d0*PI)/RNET(1:JRI)
            endif
!
!        read total kinetic-energy dependent part of potential (only for mgga)
!        (UNIT 15 = case.vtau)
!        norm of VTAU=VTAU*1
!
           if (mgga) then
              IF(myid.eq.0) then
                 READ (15,5010)
                 READ (15,5020) LM2
              endif
#ifdef Parallel
                 CALL MPI_BCAST(LM2, 1, MPI_INTEGER, 0, &
                      MPI_COMM_WORLD, ierr)
#endif
              if(LM2 .ne.  LMMAX(JATOM) + 1) goto 905 ! error in LM list
              LM2 = -1
              DO LM1 = 0, LMMAX(JATOM)
                 IF(myid.eq.0) READ (15,5030) LLL, MM
#ifdef Parallel
                 CALL MPI_BCAST(LLL, 1, MPI_INTEGER, 0, &
                      MPI_COMM_WORLD, ierr)
                 CALL MPI_BCAST(MM, 1, MPI_INTEGER, 0, &
                      MPI_COMM_WORLD, ierr)
#endif
                 LL = ABS(LLL)
                 LM2 = LM2 + 1
                 IF (LM2 .GT.  LMMX) THEN
!
!        Error: LMMX too small - stop execution   
                    GOTO 900
                 ENDIF
                 IF(myid.eq.0) READ (15,5040) (VTAU(J,LM2),J=1,JRJ(JATOM))
     
#ifdef Parallel
                 CALL MPI_BCAST(VTAU(1,LM2), JRJ(JATOM), MPI_DOUBLE_PRECISION, 0, &
                      MPI_COMM_WORLD, ierr)
#endif
!     
                 IF (LL .GT. LNSMAX*2) THEN
                    LM2 = LM2 - 1
                 ELSE
                    if(LM(1,LM2,JATOM) .ne. LLL) goto 905 ! error in LM list
                    if(LM(2,LM2,JATOM) .ne. MM) goto 905 !  error in LM list
                 ENDIF
                 IF(myid.eq.0) READ (15,5060)
              END DO
              IF(myid.eq.0) READ (15,5050)
           endif
!
!          Read VSP associated with local part of MGGA potential
!          ONLY FOR MGGA CENTRAL CORRECTION
           if(mgga_central_correction) then
              LMMIN = 0 ! there is a spherical matrix element!
!          First, SUBTRACT VSP of fully local potential from VLM(R,0)
!          Note the different normalisation of VR compared to VSP
!          this will be renormalized after VSP_MGGA is added
               VLM(1:JRI,0) = -1.0d0*VR(1:JRI,JATOM)
               VR(:,JATOM) = 0.0d0

!          Store VSP of MGGA potential in VR, then ADD to VLM(R,0)
              if(myid.eq.0) then
                 if(jatom.eq.1)  READ (17,5070)
                 READ (17,5010)
                 READ (17,5020) IDUMMY
                 READ (17,5060)
                 READ (17,5040) (VR(J,JATOM),J=1,JRJ(JATOM))
                 READ (17,5060)
                 READ (17,5050)
! Rydberg to Hartree
                                                                                                                                                            VR(:JRJ(JATOM),JATOM) = VR(1:JRJ(JATOM),JATOM)/2.0d0
              endif
#ifdef Parallel
              CALL MPI_BCAST(vr(1,jatom), JRJ(JATOM) , MPI_DOUBLE_PRECISION, 0,  &
                   MPI_COMM_WORLD, ierr)
#endif
              VLM(1:JRI,0) = VLM(1:JRI,0) + VR(1:JRI,JATOM)
!          Normalise & Hartree to Rydberg
              VLM(1:JRI,0) = 2.d0*VLM(1:JRI,0)*SQRT(4.0d0*PI)/RNET(1:JRI)
           end if
!        calculate possible nonspherical contributions to H
!
           LQIND(JATOM) = 0
            LLMM = LM2
!     
            DO 300 L0 = 0, LNSMAX
               DO 290 LP = 0, LNSMAX
                  DO 280 LMX = LMMIN, LLMM    ! 1, LLMM
                    LL = ABS(LM(1,LMX,JATOM))
                    LLL= LM(1,LMX,JATOM)
                    MM = LM(2,LMX,JATOM)
!     
!     which gaunts are not zero ?
!
                     IF (MOD((L0+LP+LL),2) .EQ. 1) GOTO 280

                     ! The ordinary Hamiltonian needs only triples obeying
                     ! the triangle rule for <Y_l0|Y_LL|Y_lp>.  The last two
                     ! lines of stress Eq. (G2), however, insert a rank-two
                     ! strain tensor.  They can therefore require radial
                     ! V/W matrices for triples whose *ordinary* Gaunt is
                     ! zero but which become allowed after coupling to
                     ! ranks 0 or 2.  Preserve the production list exactly;
                     ! append those auxiliary radial channels only for the
                     ! explicit G2 diagnostic.
                     ham_gaunt_allowed = (L0+LP-LL).GE.0 .AND. &
                                          (L0-LP+LL).GE.0 .AND. &
                                          (LP-L0+LL).GE.0
                     IF (.NOT.ham_gaunt_allowed .AND. &
                         .NOT.g2_diag_enabled) GOTO 280
!     
 240                 CONTINUE
                     DO 270 M = -L0, L0
                        DO 260 MP = -LP, LP
                           IF (-M+MM+MP .ne. 0) CYCLE
                           LQIND(JATOM) = LQIND(JATOM) + 1
                           IF (LQIND(JATOM) .GT. NGAU) THEN
                              ngau=ngau*2
                              call doreallocate(GFAC, NGAU, NAT)
                              call doreallocate(LQNS,6,NGAU,NAT)
                              call doreallocate(LQMM,NGAU,NAT)
!     Error: more then NGAU gaunt-factors - stop execution
!
!                              GOTO 910
                           ENDIF
                           IF (ham_gaunt_allowed) THEN
                              GFAC(LQIND(JATOM),JATOM) = &
                                  atom_gaunt(L0,LLL,LP,M,MM,MP,iatnr(jatom))
                           ELSE
                              ! Never let an auxiliary stress channel enter
                              ! the Hamiltonian or force contraction.
                              GFAC(LQIND(JATOM),JATOM) = (0.d0,0.d0)
                           ENDIF
!!                       GAUNT(L0,LL,LP,M,MM,MP)
!!                           IF (IATNR(JATOM) .GT. 0) THEN
!!			      IF (LL.EQ.3.OR.LL.EQ.7.OR.LL.EQ.9) THEN
!!				 IF(MM.LT.0) THEN
!!                                    GFAC(LQIND(JATOM),JATOM) = &
!!                                      GFAC(LQIND(JATOM),JATOM)*IMAG
!!				 ELSE
!!                                    GFAC(LQIND(JATOM),JATOM) = &
!!                                      -GFAC(LQIND(JATOM),JATOM)*IMAG
!!				 ENDIF
!!                              ENDIF
!!                           ELSE
!!                              IF (MM .NE. 0) THEN
!!                                 MINU = 1
!!                                 IMAG1 = (1.0D+0,0.0D+0)
!!                                 IF (LM(1,LMX,JATOM) .LT. 0) THEN
!!                                    IMAG1 = -IMAG
!!                                    MINU = -1
!!                                 ENDIF
!!                                 IF (MOD(MM,2) .EQ. 1) THEN
!!                                    IMAG1 = -IMAG1
!!                                    MINU = -MINU
!!                                 ENDIF
!!                                 IF (MM .GT. 0) MINU = 1
!!                                 GFAC(LQIND(JATOM),JATOM) = &
!!                                      GFAC(LQIND(JATOM),JATOM)*IMAG1* &
!!                                      DBLE(MINU)/SQRT(2.0D+0)
!!                              ENDIF
!!                           ENDIF
                           LQNS(1,LQIND(JATOM),JATOM) = L0 + 1
                           LQNS(2,LQIND(JATOM),JATOM) = LL
                           LQNS(3,LQIND(JATOM),JATOM) = LP + 1
                           LQNS(4,LQIND(JATOM),JATOM) = LMX
                           LQNS(5,LQIND(JATOM),JATOM) = M
                           LQNS(6,LQIND(JATOM),JATOM) = MP
                           ! LM(2,LMX,JATOM) contains only |M|.  The loop
                           ! below explicitly visits the conjugate -M
                           ! channel, whose sign must survive until FORFHS.
                           LQMM(LQIND(JATOM),JATOM) = MM
!                           IF(myid.eq.0) WRITE(6,6050) L0, LL, LP, M, MM, MP,
!     ,                                   GFAC(LQIND(JATOM),JATOM)
 260                    CONTINUE
 270                 CONTINUE
!
                     IF (MM.GT.0) THEN
                        MM = -MM
                        GOTO 240
                     ENDIF
 280              CONTINUE
 290           CONTINUE
 300        CONTINUE
              if (g2_diag_enabled.and.myid.eq.0) then
                 g2_bad_count=0
                 do lqx=1,lqind(jatom)
                    if (-lqns(5,lqx,jatom)+lqmm(lqx,jatom)+ &
                         lqns(6,lqx,jatom).ne.0) then
                       g2_bad_count=g2_bad_count+1
                       if (g2_bad_count.le.20) write(6,*) &
                            ':DBG_G2_ATPAR_M:',jatom,lqx, &
                            lqns(5,lqx,jatom),lqmm(lqx,jatom), &
                            lqns(6,lqx,jatom)
                    endif
                 enddo
                 write(6,*) ':DBG_G2_ATPAR_M_SUM:',jatom, &
                              lqind(jatom),g2_bad_count
              endif
              WRITE (6,*) LQIND(JATOM),' POSSIBLE NONSPHERICAL CONTRIBUTIONS TO H'
!
            IF ((IATNR(JATOM) .GT. 0) .AND. (LMMAX(JATOM) .GT. 1)) THEN
!
!        combine the radial total potential coefficients according to
!        cubic harmonic
!
!!              CALL CBCOMB(JRI,LLMM,VLM(1,1),LM(1,1,1),JATOM)
              CALL CBCOMB(JRI,LLMM,VLM(1:NRAD,1:LMMX),LM(1:2,1:LMMX,1:NAT),JATOM)
              if(mgga) CALL CBCOMB(JRI,LLMM,Vtau(1:NRAD,1:LMMX),LM(1:2,1:LMMX,1:NAT),JATOM)
!!              if(mgga) CALL CBCOMB(JRI,LLMM,Vtau(1,1),LM(1,1,1),JATOM)
            ENDIF
!
!        integrate product of radialfunctions (and derivatives) and
!        total potential coefficients (combinations) from 0 to R-MT
!
            ICOUNT = 0
            L0BEFO = 999
            LLBEFO = 999
            LPBEFO = 999
            LMBEFO = 999
            DO 330 LQX = 1, LQIND(JATOM)
               L0 = LQNS(1,LQX,JATOM) - 1
               LL = LQNS(2,LQX,JATOM)
               LP = LQNS(3,LQX,JATOM) - 1
               LM1 = LQNS(4,LQX,JATOM)
               IF ((L0 .EQ. L0BEFO) .AND. (LL .EQ. LLBEFO) .AND. &
                   (LP .EQ. LPBEFO) .AND. (LM1 .EQ. LMBEFO)) THEN
                  M = LQNS(5,LQX,JATOM)
                  MP = LQNS(6,LQX,JATOM)
!                  WRITE (6,6050) L0, LL, LP, M, LM(2,LM1,JATOM), MP,
!     ,                 GFAC(LQX,JATOM)
               ELSE
!     
                  A(1:JRI) = A1(1:JRI,L0+1)*VLM(1:JRI,LM1)
                  B(1:JRI) = B1(1:JRI,L0+1)*VLM(1:JRI,LM1)
                  AE(1:JRI) = AE1(1:JRI,L0+1)*VLM(1:JRI,LM1)
                  BE(1:JRI) = BE1(1:JRI,L0+1)*VLM(1:JRI,LM1)
                  ICOUNT = ICOUNT + 1
!
                  CALL RINT13(CFEIN1,CFEIN2,A,B,A1(1,LP+1),B1(1,LP+1), &
                       AAINT,JATOM)
                  VNS1(L0+1,LM1,LP+1,JATOM) = AAINT
                  CALL RINT13(CFEIN1,CFEIN2,AE,BE,AE1(1,LP+1), &
                       BE1(1,LP+1),AAINT,JATOM)
                  VNS2(L0+1,LM1,LP+1,JATOM) = AAINT
                  CALL RINT13(CFEIN1,CFEIN2,A,B,AE1(1,LP+1),BE1(1,LP+1), &
                       AAINT,JATOM)
                  VNS3(L0+1,LM1,LP+1,JATOM) = AAINT

                  ! The fifth line of stress Eq. (G2) contains
                  !
                  !   W_{l,L,l'} = int r^3 dr [d u_l(r)/dr] V_L(r) u_l'(r).
                  !
                  ! LAPW1 stores A1=r*u and AP1=r*du/dr.  RINT13 integrates
                  ! products of these stored radial functions over dr, so the
                  ! bra supplied to it must be r*AP1*V, not r*A1*V.  The
                  ! leading RNET supplies the remaining power of r in the
                  ! definition above.  This distinction is essential: W is a
                  ! radial-derivative matrix element, not a radial moment of V.
                  !
                  ! The derivative is on the bra only.  Consequently W is not
                  ! symmetric under bra/ket interchange, unlike the ordinary
                  ! VNS matrix elements.  Store all four independent APW
                  ! channels explicitly: u'-u, udot'-udot, u'-udot, and
                  ! udot'-u.
                  AW(1:JRI)  = RNET(1:JRI)*AP1(1:JRI,L0+1)*VLM(1:JRI,LM1)
                  BW(1:JRI)  = RNET(1:JRI)*BP1(1:JRI,L0+1)*VLM(1:JRI,LM1)
                  AEW(1:JRI) = RNET(1:JRI)*APE1(1:JRI,L0+1)*VLM(1:JRI,LM1)
                  BEW(1:JRI) = RNET(1:JRI)*BPE1(1:JRI,L0+1)*VLM(1:JRI,LM1)
                  CALL RINT13(CFEIN1,CFEIN2,AW,BW,A1(1,LP+1),B1(1,LP+1), &
                       AAINT,JATOM)
                  WNS1(L0+1,LM1,LP+1,JATOM) = AAINT
                  CALL RINT13(CFEIN1,CFEIN2,AEW,BEW,AE1(1,LP+1), &
                       BE1(1,LP+1),AAINT,JATOM)
                  WNS2(L0+1,LM1,LP+1,JATOM) = AAINT
                  CALL RINT13(CFEIN1,CFEIN2,AW,BW,AE1(1,LP+1),BE1(1,LP+1), &
                       AAINT,JATOM)
                  WNS3(L0+1,LM1,LP+1,JATOM) = AAINT
                  CALL RINT13(CFEIN1,CFEIN2,AEW,BEW,A1(1,LP+1),B1(1,LP+1), &
                       AAINT,JATOM)
                  WNS4(L0+1,LM1,LP+1,JATOM) = AAINT

                  ! These are the independent APW-bra/LO-ket W channels.
                  ! They cannot be reconstructed by transposing the LO-bra
                  ! integrals below because the radial derivative remains on
                  ! the APW bra.
                  IF (LP .LE. LOMAX) THEN
                     do klo=1,ilo(lp,jatom)
                        if (.not.lapw(lp,jatom).and.klo.eq.1) cycle
                        if (.NOT. loor(lp,jatom)) cycle
                        CALL RINT13(CFEIN1,CFEIN2,AW,BW, &
                             A1LO(1,klo,LP),B1LO(1,klo,LP),AAINT,JATOM)
                        WNS4LO(L0+1,LM1,LP+1,klo,JATOM) = AAINT
                        CALL RINT13(CFEIN1,CFEIN2,AEW,BEW, &
                             A1LO(1,klo,LP),B1LO(1,klo,LP),AAINT,JATOM)
                        WNS5LO(L0+1,LM1,LP+1,klo,JATOM) = AAINT
                     enddo
                  ENDIF
!
!     integrals needed for local orbitals
!
                  IF (L0 .LE. LOMAX) THEN

                     do jlo=1,ilo(l0,jatom)

                        VNS1LO(L0+1,LM1,LP+1,jlo,JATOM) = 0.0D0
                        VNS2LO(L0+1,LM1,LP+1,jlo,JATOM) = 0.0D0
                        VNS3LO(L0+1,LM1,LP+1,jlo,:,JATOM) = 0.0d0
                        WNS1LO(L0+1,LM1,LP+1,jlo,JATOM) = 0.0D0
                        WNS2LO(L0+1,LM1,LP+1,jlo,JATOM) = 0.0D0
                        WNS3LO(L0+1,LM1,LP+1,jlo,:,JATOM) = 0.0d0

                        if (.not.lapw(l0,jatom).and.jlo.eq.1) cycle
                        if( .NOT. loor(l0,jatom)) cycle

                        A(1:JRI) = A1LO(1:JRI,jlo,L0)*VLM(1:JRI,LM1)
                        B(1:JRI) = B1LO(1:JRI,jlo,L0)*VLM(1:JRI,LM1)

                        CALL RINT13(CFEIN1,CFEIN2,A,B,A1(1,LP+1),B1(1,LP+1),AAINT,JATOM)
                        VNS1LO(L0+1,LM1,LP+1,jlo,JATOM) = AAINT
                        CALL RINT13(CFEIN1,CFEIN2,A,B,AE1(1,LP+1),BE1(1,LP+1),AAINT,JATOM)
                        VNS2LO(L0+1,LM1,LP+1,jlo,JATOM) = AAINT
                        ! The LO is now the bra, hence AP1LO (r*u'_LO) is the
                        ! differentiated radial function required by W.
                        AW(1:JRI) = RNET(1:JRI)*AP1LO(1:JRI,jlo,L0)*VLM(1:JRI,LM1)
                        BW(1:JRI) = RNET(1:JRI)*BP1LO(1:JRI,jlo,L0)*VLM(1:JRI,LM1)
                        CALL RINT13(CFEIN1,CFEIN2,AW,BW,A1(1,LP+1),B1(1,LP+1),AAINT,JATOM)
                        WNS1LO(L0+1,LM1,LP+1,jlo,JATOM) = AAINT
                        CALL RINT13(CFEIN1,CFEIN2,AW,BW,AE1(1,LP+1),BE1(1,LP+1),AAINT,JATOM)
                        WNS2LO(L0+1,LM1,LP+1,jlo,JATOM) = AAINT
                        IF (LP .LE. LOMAX) THEN
                           do klo=1,ilo(lp,jatom)
                           if (.not.lapw(lp,jatom).and.klo.eq.1) cycle
                              CALL RINT13(CFEIN1,CFEIN2,A,B,A1LO(1,klo,LP),B1LO(1,klo,LP),AAINT,JATOM)
                              VNS3LO(L0+1,LM1,LP+1,jlo,klo,JATOM) = AAINT
                              CALL RINT13(CFEIN1,CFEIN2,AW,BW,A1LO(1,klo,LP),B1LO(1,klo,LP),AAINT,JATOM)
                              WNS3LO(L0+1,LM1,LP+1,jlo,klo,JATOM) = AAINT
                           enddo
                        M = LQNS(5,LQX,JATOM)
                        MP = LQNS(6,LQX,JATOM)
!                     WRITE (6,6060) (L0, LL, LP, M, LM(2,LM1,JATOM), MP,jlo,klo, &
!                         GFAC(LQX,JATOM), &
!                         VNS1LO(L0+1,LM1,LP+1,jlo,JATOM), &
!                         VNS2LO(L0+1,LM1,LP+1,jlo,JATOM), &
!                         VNS3LO(L0+1,LM1,LP+1,jlo,klo,JATOM), LM1,klo=1,ilo(lp,jatom))
                        ENDIF
                     enddo

                  ENDIF
 320              continue
!        end of lo-integrals
!
                  L0BEFO = L0
                  LLBEFO = LL
                  LPBEFO = LP
                  LMBEFO = LM1
                  M = LQNS(5,LQX,JATOM)
                  MP = LQNS(6,LQX,JATOM)
!                  WRITE (6,6050) L0, LL, LP, M, LM(2,LM1,JATOM), MP, &
!                                GFAC(LQX,JATOM),&
!                                VNS1(L0+1,LM1,LP+1,JATOM),&
!                                VNS2(L0+1,LM1,LP+1,JATOM),&
!                                VNS3(L0+1,LM1,LP+1,JATOM), LM1
               ENDIF
  330       CONTINUE
!
!fb
            if (force) then
                call forfhs(jatom) 
!
!               final line of [yu91]-output
!               final line for each atom = 6 x 0 and 8 x 0.0d0
                IF(myid.eq.0) write(71) &
                      ! Four matrix blocks: original V, original W,
                      ! phase-weighted bare V, and phase-weighted bare W.
                      (0, i=1,6),(0.d0,i=1,32)
            endif
!fe

!           MGGA contributions
!           Here, <u_l | VTAU(LM) | u_l'> couple to
!           GAUNT(L0+-1,LL,LP+-1,M+-1,MM,MP+-1) !!
            IF (mgga) THEN
               LQIND(JATOM) = 0
               LLMM = LM2
               MGGA_GAUNTS = 0
               !              reallocate in case NGAU was raised above
               call doreallocate(MGGA_FAC,4,NGAU,NAT)
!              Find the additional L0,LL,LP combinations, add them to LQNS,GFAC
!              All previous Gaunts give an MGGA contribution as well
               DO L0 = 0, LNSMAX
                  DO LP = 0, LNSMAX
!                    note that LMX = 0, .. (not 1 ,) because
!                    VTAU starts with L=0
!                    whereas VNS starts with L=1
                     DO LMX = 0, LLMM
                        LL = ABS(LM(1,LMX,JATOM))
                        LLL= LM(1,LMX,JATOM)
                        MM = (LM(2,LMX,JATOM)) ! ABS?

                        IF( MOD((L0+LP+LL),2) .NE. 0)  CYCLE

!                       check if this combination is already included above
!                       (has a VNS contribution)
                        logic_mgga_gaunt(0) = ( ((L0+LP-LL) .GE. 0) .AND. &
                             & ((L0-LP+LL) .GE. 0) .AND. ((LP-L0+LL) .GE. 0) )

!                       check if this combination has a mgga contribution
!                                    Contribution of gaunt(l0+1,ll,lp+1,...)
                        logic_mgga_gaunt(4) = (IABS(L0-LP) .LE. LL .AND. L0+LP-LL .GE. -2)
!                                    Contribution of gaunt(l0+1,ll,lp-1,...)
                        logic_mgga_gaunt(3) = (IABS(L0-LP+2) .LE. LL .AND. L0+LP-LL .GE. 0)
!                                    Contribution of gaunt(l0-1,ll,lp+1,...)
                        logic_mgga_gaunt(2) = (IABS(L0-LP-2) .LE. LL .AND. L0+LP-LL .GE. 0)
!                                    Contribution of gaunt(l0-1,ll,lp-1,...)
                        logic_mgga_gaunt(1) = (IABS(L0-LP) .LE. LL .AND.  L0+LP-LL .GE. 2)

                        IF( .not. any(logic_mgga_gaunt)) CYCLE
!                       Check M's
                        DO M = -L0, L0
                           DO MP = -LP, LP
                              if ( M .eq.  MM + MP) then
                                 continue
                              elseif ( M .eq.  MP - MM) then
                                 MM = -MM 
                                 continue
                              else
                                 cycle
                              end if

                              if (.not.  logic_mgga_gaunt(0)) then
!                                this combination has a mgga contribution but no 
!                                vns contribution
                                 MGGA_GAUNTS = MGGA_GAUNTS + 1
                              endif

                             LQIND(JATOM) = LQIND(JATOM) + 1

                             IF (LQIND(JATOM) .GT.  NGAU) THEN
                                ngau=ngau*2
                                call doreallocate(GFAC, NGAU, NAT)
                                call doreallocate(MGGA_FAC, 4, NGAU, NAT)
                                call doreallocate(LQNS,6,NGAU,NAT)
                                call doreallocate(LQMM,NGAU,NAT)
                             ENDIF
                             
                             LQNS(1,LQIND(JATOM),JATOM) = L0 + 1
                             LQNS(2,LQIND(JATOM),JATOM) = LL
                             LQNS(3,LQIND(JATOM),JATOM) = LP + 1
                             LQNS(4,LQIND(JATOM),JATOM) = LMX
                             LQNS(5,LQIND(JATOM),JATOM) = M
                             LQNS(6,LQIND(JATOM),JATOM) = MP
                             LQMM(LQIND(JATOM),JATOM) = MM
!                              Calculate Gaunt and mgga_fac for this combination
                              call calculate_mgga_fac(logic_mgga_gaunt,&
                                   IATNR(JATOM),                       &
                                   GFAC(LQIND(JATOM),JATOM), &
                                   MGGA_FAC(:,LQIND(JATOM),JATOM), &
                                   L0,LLL,LP,M,MM,MP)
                           end DO
                        end DO
                     end DO
                  end DO
               end DO

               IF ((IATNR(JATOM) .GT. 0) .AND. (LMMAX(JATOM) .GT. 1)) THEN
!
!        combine the radial total potential (real harmonic) coefficients
!        according to cubic harmonic for vtau
!
!                  CALL CBCOMB(JRI,LLMM,VTAU,LM,JATOM)
               ENDIF

!              Calculate integrals for MGGA matrix elements <u|vtau|u>
               MCOUNT = 0
               L0BEFO = -1
               LLBEFO = -1
               LPBEFO = -1
               LMBEFO = -1
               DO LQX = 1, LQIND(JATOM)
                  L0 = LQNS(1,LQX,JATOM) - 1
                  LL = LQNS(2,LQX,JATOM)
                  LP = LQNS(3,LQX,JATOM) - 1
                  LM1 = LQNS(4,LQX,JATOM)
                  IF ((L0 .EQ. L0BEFO) .AND. (LL .EQ. LLBEFO) .AND.  &
                       (LP .EQ. LPBEFO) .AND. (LM1 .EQ. LMBEFO)) THEN
                     M = LQNS(5,LQX,JATOM)
                     MP = LQNS(6,LQX,JATOM)
                     CYCLE
                  END IF
                  MCOUNT = MCOUNT + 1

!                 In the following section A/B store first function times VTAU
!                 and AE/BE used as needed to store temporary A1/R arrays

                  A(1:JRI) = A1(1:JRI,L0+1)*VTAU(1:JRI,LM1)/RNET(1:JRI)
                  B(1:JRI) = B1(1:JRI,L0+1)*VTAU(1:JRI,LM1)/RNET(1:JRI)

                  AE(1:JRI) = A1(1:JRI,LP+1)/RNET(1:JRI)
                  BE(1:JRI) = B1(1:JRI,LP+1)/RNET(1:JRI)
                  CALL RINT13(CFEIN1,CFEIN2,A,B,AE,BE,AAINT,JATOM)
                  INT_U___VTAU_U__(L0+1,LM1,LP+1,JATOM) = AAINT

                  AE(1:JRI) = AE1(1:JRI,LP+1)/RNET(1:JRI)
                  BE(1:JRI) = BE1(1:JRI,LP+1)/RNET(1:JRI)
                  CALL RINT13(CFEIN1,CFEIN2,A,B,AE,BE,AAINT,JATOM)
                  INT_U___VTAU_UD_(L0+1,LM1,LP+1,JATOM) = AAINT

                  CALL RINT13(CFEIN1,CFEIN2,A,B,AP1(:,LP+1),BP1(:,LP+1),&
                          AAINT,JATOM)
                  INT_U___VTAU_UP_(L0+1,LM1,LP+1,JATOM) = AAINT

                  CALL RINT13(CFEIN1,CFEIN2,A,B,APE1(:,LP+1),BPE1(:,LP+1),&
                          AAINT,JATOM)
                  INT_U___VTAU_UPD(L0+1,LM1,LP+1,JATOM) = AAINT
                  A(1:JRI) = AE1(1:JRI,L0+1)*VTAU(1:JRI,LM1)/RNET(1:JRI)
                  B(1:JRI) = BE1(1:JRI,L0+1)*VTAU(1:JRI,LM1)/RNET(1:JRI)
                  CALL RINT13(CFEIN1,CFEIN2,A,B,AE,BE,AAINT,JATOM)
                  INT_UD__VTAU_UD_(L0+1,LM1,LP+1,JATOM) = AAINT

                  CALL RINT13(CFEIN1,CFEIN2,A,B,APE1(:,LP+1),BPE1(:,LP+1),&
                              AAINT,JATOM)
                  INT_UD__VTAU_UPD(L0+1,LM1,LP+1,JATOM) = AAINT

                  A(1:JRI) = AP1(1:JRI,L0+1)*VTAU(1:JRI,LM1)
                  B(1:JRI) = BP1(1:JRI,L0+1)*VTAU(1:JRI,LM1)
                  CALL RINT13(CFEIN1,CFEIN2,A,B,AE,BE,AAINT,JATOM)
                  INT_UP__VTAU_UD_(L0+1,LM1,LP+1,JATOM) = AAINT

                  CALL RINT13(CFEIN1,CFEIN2,A,B,AP1(:,LP+1),BP1(:,LP+1),&
                              AAINT,JATOM)
                  INT_UP__VTAU_UP_(L0+1,LM1,LP+1,JATOM) = AAINT

                  CALL RINT13(CFEIN1,CFEIN2,A,B,APE1(:,LP+1),BPE1(:,LP+1),&
                              AAINT,JATOM)
                  INT_UP__VTAU_UPD(L0+1,LM1,LP+1,JATOM) = AAINT

                  A(1:JRI) = APE1(1:JRI,L0+1)*VTAU(1:JRI,LM1)
                  B(1:JRI) = BPE1(1:JRI,L0+1)*VTAU(1:JRI,LM1)
                  CALL RINT13(CFEIN1,CFEIN2,A,B,APE1(:,LP+1),BPE1(:,LP+1),&
                              AAINT,JATOM)
                  INT_UPD_VTAU_UPD(L0+1,LM1,LP+1,JATOM) = AAINT

!                 And for the local orbitals
                  IF (L0 .LE. LOMAX) THEN
                     DO JLO = 1, ILO(L0,JATOM)
                        IF( .NOT.  LAPW(L0,JATOM) .AND. JLO .EQ. 1) CYCLE
                        IF( .NOT.  LOOR(L0,JATOM)) CYCLE

!                       In this section, A/B holds Ulo/RNET*VTAU
!                       AP/BP holds UloP*VTAU
!                       AE/BP is used as needed to the second rad. function

                        A(1:JRI) = A1LO(1:JRI,JLO,L0)*VTAU(1:JRI,LM1)/RNET(1:JRI)
                        B(1:JRI) = B1LO(1:JRI,JLO,L0)*VTAU(1:JRI,LM1)/RNET(1:JRI)
                        AP(1:JRI) = AP1LO(1:JRI,JLO,L0)*VTAU(1:JRI,LM1)
                        BP(1:JRI) = BP1LO(1:JRI,JLO,L0)*VTAU(1:JRI,LM1)

                        AE(1:JRI) = A1(1:JRI,LP+1)/RNET(1:JRI)
                        BE(1:JRI) = B1(1:JRI,LP+1)/RNET(1:JRI)

                        CALL RINT13(CFEIN1,CFEIN2,A,B,AE,BE,AAINT,JATOM)
                        INT_Ulo__VTAU_U__(L0+1,LM1,LP+1,jlo,JATOM) = AAINT

                        CALL RINT13(CFEIN1,CFEIN2,AP,BP,AE,BE,AAINT,JATOM)
                        INT_UloP_VTAU_U__(L0+1,LM1,LP+1,jlo,JATOM) = AAINT

                        AE(1:JRI) = AE1(1:JRI,LP+1)/RNET(1:JRI)
                        BE(1:JRI) = BE1(1:JRI,LP+1)/RNET(1:JRI)

                        CALL RINT13(CFEIN1,CFEIN2,A,B,AE,BE,AAINT,JATOM)
                        INT_Ulo__VTAU_UD_(L0+1,LM1,LP+1,jlo,JATOM) = AAINT

                        CALL RINT13(CFEIN1,CFEIN2,AP,BP,AE,BE,AAINT,JATOM)
                        INT_UloP_VTAU_UD_(L0+1,LM1,LP+1,jlo,JATOM) = AAINT

                        CALL RINT13(CFEIN1,CFEIN2,A,B,AP1(:,LP+1),BP1(:,LP+1)&
                                    ,AAINT,JATOM)
                        INT_Ulo__VTAU_UP_(L0+1,LM1,LP+1,jlo,JATOM) = AAINT

                        CALL RINT13(CFEIN1,CFEIN2,AP,BP,AP1(:,LP+1),BP1(:,LP+1)&
                                    ,AAINT,JATOM)
                        INT_UloP_VTAU_UP_(L0+1,LM1,LP+1,jlo,JATOM) = AAINT

                        CALL RINT13(CFEIN1,CFEIN2,A,B,APE1(:,LP+1),BPE1(:,LP+1)&
                                    ,AAINT,JATOM)
                        INT_Ulo__VTAU_UPD(L0+1,LM1,LP+1,jlo,JATOM) = AAINT

                        CALL RINT13(CFEIN1,CFEIN2,AP,BP,APE1(:,LP+1),BPE1(:,LP+1)&
                                    ,AAINT,JATOM)
                        INT_UloP_VTAU_UPD(L0+1,LM1,LP+1,jlo,JATOM) = AAINT

                        IF (LP .LE. LOMAX) THEN
                           DO KLO = 1,ILO(LP,JATOM)
                              IF( .NOT.  LAPW(LP,JATOM) .AND.  KLO .EQ.  1) CYCLE
!                              WHY IS THIS CHECK NOT THERE IN THE VNS3 SECTION???
                              IF( .NOT.  LOOR(LP,JATOM)) CYCLE

                              AE(1:JRI) = A1LO(1:JRI,KLO,LP)/RNET(1:JRI)
                              BE(1:JRI) = B1LO(1:JRI,KLO,LP)/RNET(1:JRI)

                              CALL RINT13(CFEIN1,CFEIN2,A,B,AE,BE,AAINT,JATOM)
                              INT_Ulo__VTAU_Ulo(L0+1,LM1,LP+1,jlo,klo,JATOM) = AAINT

                              CALL RINT13(CFEIN1,CFEIN2,AP,BP,AE,BE,AAINT,JATOM)
                              INT_UloP_VTAU_Ulo(L0+1,LM1,LP+1,jlo,klo,JATOM) = AAINT

                              CALL RINT13(CFEIN1,CFEIN2,AP,BP,&
                                          AP1LO(:,KLO,LP),BP1LO(:,KLO,LP),AAINT,JATOM)
                              INT_UloP_VTAU_UloP(L0+1,LM1,LP+1,jlo,klo,JATOM) = AAINT

                           END DO
                        end IF
                     END DO
                  end IF
                  L0BEFO = L0
                  LLBEFO = LL
                  LPBEFO = LP
                  LMBEFO = LM1
                  M = LQNS(5,LQX,JATOM)
                  MP = LQNS(6,LQX,JATOM)
               end DO
            end IF

!fb
!!            if (force) then
!!                call forfhs(jatom) 
!!!
!!!               final line of [yu91]-output
!!!               final line for each atom = 6 x 0 and 8 x 0.0d0
!!                IF(myid.eq.0) write(71) &
!!                      (0, i=1,6),(0.d0,i=1,8)
!!            endif
!fe
            IF(myid.eq.0) WRITE (6,*) '   NUMBER OF RADIAL INTEGRALS FOR ATOM', &
                        JATOM,' = ',ICOUNT
         ENDIF
!
!        end of loop over atoms in unitcell
!
  340 CONTINUE
      
      CALL END_RADIAL
!
      RETURN
!
  900 CONTINUE
      WRITE (ERRMSG,9000) 'LMMX too small', LMMX, LM2
      CALL OUTERR('ATPAR',ERRMSG)
      STOP 'ATPAR - Error'
!
  901 CONTINUE
      WRITE (ERRMSG,9001) 'LOMAX too small', LOMAX, L
      CALL OUTERR('ATPAR',ERRMSG)
      STOP 'ATPAR - Error'
!
  902 CONTINUE
      WRITE (ERRMSG,9002) 'NLOAT too small', NLOAT,l,jatom
      CALL OUTERR('ATPAR',ERRMSG)
      STOP 'ATPAR - Error'
!
  903 CONTINUE
      WRITE (ERRMSG,9003) 'NSLMAX too small', NSLMAX,l,jatom
      CALL OUTERR('ATPAR',ERRMSG)
      STOP 'ATPAR - Error'
!
  904 CONTINUE
      WRITE (ERRMSG,9010) 'ERROR reading VTAU (unit 15)'
      CALL OUTERR('ATPAR',ERRMSG)
      STOP 'ATPAR - Error'
!
  905 CONTINUE
      WRITE (ERRMSG,9010) 'ERROR in VTAU (unit 15): error in LM list'
      CALL OUTERR('ATPAR',ERRMSG)
      STOP 'ATPAR - Error'
!
            910 CONTINUE
      CALL OUTERR('ATPAR','more than NGAU gaunts')
      WRITE (ERRMSG,9010) '  NGAU,   L0,   LP,   LL,    M,   MP,   MM'
      CALL OUTERR('ATPAR',ERRMSG)
      WRITE (ERRMSG,9020) NGAU, L0, LP, LL, M, MP, MM
      CALL OUTERR('ATPAR',ERRMSG)
      STOP 'ATPAR - Error'
!
 5000 FORMAT (1X,I1,2F10.5,A4,i2)
 5010 FORMAT (3X)
 5020 FORMAT (15X,I3,/,/)
 5030 FORMAT (15X,I3,5X,I2,/)
#ifdef Extended
 5040 FORMAT (3X,4D27.20)
#else
 5040 FORMAT (3X,4es19.12)
#endif
 5050 FORMAT (/,/,/)
 5060 FORMAT (/)
 5070 FORMAT (/,/)
 6000 FORMAT (/,10X,'ATOMIC SPHERE DEPENDENT PARAMETERS FOR ATOM  ',A10, &
              /,':e__',i4.4,':',1X,'OVERALL ENERGY PARAMETER IS',F10.4 &
              /,10X,'OVERALL BASIS SET ON ATOM IS',A5)
 6010 FORMAT (':E',i1,'_',i4.4,':',1X,'E(',I2,2H)=,F10.4)
 6020 FORMAT (/,10X,'POTENTIAL PARAMETERS ',/,11X,'L',5X,'U(R)',10X, &
              'U''(R)',9X,'DU/DE',8X,'DU''/DE',6X,'NORM-DU/DE')
 6021 FORMAT (/,10X,'LOCAL ORBITAL POTENTIAL PARAMETERS',/, &
           11X,'L',5X,'U(R)',10X, &
           'U''(R)',5X,'NORM U1U2',8X,'NORM UE1U2')
 6030 FORMAT (10X,I2,5E14.6,5X,3I2)
 6031 FORMAT (10X,I2,4E14.6,5X,3I2)
 6040 FORMAT (/,/,/,' AT.NR.:',I2,5X,'POSITION: ',3F8.3,5X, &
              'MULTIPLICITY:',I5)
 6050 FORMAT (3X,6I3,5F15.8,I3)
 6060 FORMAT ('lo ',8I3,5F15.8,I3)
 6070 FORMAT (10X,'FOR L=',I2,' CORRECTION=',E14.5,' OVERLAP=',E14.5)
 6080 FORMAT (/,3X,'NOT EQUIV ATOM ',i4,3x,A10)
 6090 FORMAT (30X,3F10.5)
 6100 FORMAT (13X,'EQUIV ATOM ',I3,3X,'POSITION: ',3F8.3)
 6110 FORMAT (30X,3F10.5)
 9000 FORMAT (A,' (LMMX =',I5,', LM2 =',I5,')')
 9001 FORMAT (A,' (LOMAX =',I5,', L  =',I5,')')
 9002 FORMAT (A,' (NLOAT =',I3,', L =',I3,' JATOM =',i3,')')
 9003 FORMAT (A,' (NSLMAX=',I3,', L =',I3,' JATOM =',i3,')')
 9010 FORMAT (A)
 9020 FORMAT (I6,',',6(I5,','))
!
!        End of 'ATPAR'
!
      END






