SUBROUTINE read_and_initialize(NUMKPT,DEFFN,INFO,nmat_only)
  use param, only   : filename_v_sph, lmax, lmmx, lomax, ngau, nloat, nrad, nslmax
  use loabc, only   : ALO, BLO, CLO, ELO, nloat_new, PLO, DPLO, PI12LO, PE12LO, init_loabc
  use lapw_timer, only : READ_CPU_TIME, time_g1, time_atpar,START_TIMER, STOP_TIMER
  use lolog, only   : lapw, ilo, loor, nlo, init_lolog, nlov
  use orb, only     : init_orb
  use atspdt, only  : E, DP, DPE, P, PE, PEI, init_atspdt
  use comc, only    : IPGR, KNAME
  use comi, only    : NT, LNSMAX, NAT,nemax_it
  use coml, only    : SPRSWF, PRNTWF, WFTAPE, REL, ITER, NOHNS
  use comr, only    : RKM, ELow, EUp, WEIGHT, sproj_limit,eshift_iter,readHinv
  use kpts, only    : K3
  use potnlc, only  : VR
  use mgaunt, only  : gaunt_init
  use structure, only: lattic, aname, alpha, alat, jri, r0, dx, POS, RMT, V, VOL, MULT, ndf, ZZ, BR1, BR2, ORTHO, ReadStructure, WriteInfoStructure, iord, imat, tau, rotij, tauij, spaceGroupInit, ReadStructureGroup
  use readPotential, ONLY: init_V_vsp, read_V_vsp, read_V_vns_interstital_nk, read_V_vns_interstital_pot_real, read_V_vns_interstital_pot_cmplx
  use radial, only  : init_radial, end_radial
  use loint, only   : init_loint
  use nmr, only : dokorig, korig, writeham
  use mpi, only : Qprint, mpi_bcast_V_vsp, vector_para, myrank, master
  IMPLICIT NONE
  !INCLUDE 'param.inc'
  !
  INTEGER,      intent(out) :: numkpt, info
  CHARACTER*(*), intent(in) :: DEFFN
  LOGICAL,      intent(out) :: nmat_only
  ! locals
  INTEGER        :: ITAPE, JATOM, J, jlo
  INTEGER        :: I, ios
  INTEGER        :: L, K
  REAL*8         :: E1n, E2n, E1, E2, eorb_ind, efermi
  REAL*8         :: PI
  CHARACTER*67   :: ERRMSG
  LOGICAL        :: RUNORB
  LOGICAL        :: FORCE
  INTEGER        :: lxdos, iscf, tnlov
  !REAL*8         :: BRNN(3,3)
  !  External Subroutines
  EXTERNAL           OUTERR, LATGEN, NN, ATPAR
  EXTERNAL           RDSWAR
  !
  CALL START_TIMER(time_g1)  !        initialize timing for 'INILPW'
  !
  PI=ACOS(-1.0D0)
  !        open all files listed in 'lapw1.def'
  call read_def_file(deffn, nloat_new, nohns, iter, readHinv, force, runorb, nmat_only, dokorig, writeham, ERRMSG, info)
  if (info.ne.0) then
     CALL OUTERR('read_d',ERRMSG); return
  endif
  !
  if (Qprint) then
#ifdef _REAL
      WRITE(6,*) 'real lapw1 version'
#else
      WRITE(6,*) 'complex lapw1c version'
#endif
      WRITE (6,*) 'max number of local orbitals/atom',nloat_new
  endif
  !
  NT = LMAX
  !
  ! read crystal-structure (atompositions, symmetry-parameters,muffin-tin radius, ...)
  ! allocates: rmt(nat),v(nat),iatnr(nat),mult(nat),isplit(nat),rotloc(3,3,nat),r0(nat), dx(nat), jri(nat),aname(nat),ZZ(nat),pos(3,ndf)
  call ReadStructure(20,nat,rel,lxdos,ERRMSG,info)
  if (info.ne.0) then
     CALL OUTERR('structure.f90',ERRMSG); return
  endif
  !
  ! read the orbital potential
  call ReadOrbitalPotential(nat,runorb)
  !
  ! Reading case.in1
  call ReadFirstPart_of_in1(PRNTWF, SPRSWF, WFTAPE, Efermi, NT, RKM, LNSMAX, NSLMAX, LMAX)
  if (Qprint) WRITE (6,6020) RKM, NT-1, LNSMAX, NAT
  !
  if (Qprint) call WriteInfoStructure(6, nat)
  !
  CALL generate_lattice(alpha, BR1, BR2, VOL, ORTHO, ALAT(1), ALAT(2), ALAT(3), LATTIC)
  if (Qprint) then
     WRITE(6,*) 'BR1='
     WRITE(6,'(3F12.6)') BR1
     WRITE(6,*) 'BR2='
     WRITE(6,'(3F12.6)') BR2
     !call direct_lattice(BRNN, nat, alpha, LATTIC)
     !WRITE(6,*) 'DIRECT LATTICE='
     !WRITE(6,'(3F12.6)'), BRNN
     !        define spherevolume/cellvolume
  endif
  V(1:nat) = 4.0d0*PI*Rmt(1:nat)**3/(3.0D0*VOL)
  !
  !KH: Warning -- check if high-precision atoms are correctly read. Not sure....
  CALL ReadStructureGroup(20, nat)
  !
  ! define rotation matrices if required
  allocate( rotij(3,3,ndf) )
  allocate( tauij(3,ndf) )
  call get_rotij_tauij(rotij,tauij,pos,alat,imat,tau,iord,nat,ndf,mult,lattic)
  ! Transforms rotij to cartesian coordinates by transformation : rotij_cartesian = BR2 * rotij * ((BR2^T)^{-1})^T
  CALL correct_rotij(rotij,lattic,BR2,ORTHO,mult,nat,ndf)
  !  check that atomic spheres are not overlapping
  CALL NN(NAT)
  !  initialize GAUNT-routine
  call gaunt_init
  CALL START_TIMER(time_atpar)
  !  solve radial Schroedinger-equation, calculate gaunts and
  !  radial integrals UVU, UVUE of nonspherical potential
  nloat=nloat_new
  call init_loabc(lomax,nloat,nat) ! local orbitals arrays
  call init_lolog(lomax,nat,lmax)  ! 
  call init_atspdt(nat,lmmx,ngau,lmax,nslmax)
  call init_loint(lomax+1,lmmx,nslmax,nat,nloat)
  ! initialization of storage for the radial wave functions
  call init_radial(nrad,lomax,nslmax,nloat,lmmx)
  ! spherical radial potential allocated
  allocate( VR(NRAD,NAT) )
  ! Opens file with spherically symmetric potential
  if (myrank.eq.master) then
     CALL init_V_vsp(filename_V_sph, 18, ISCF)
     ! Reads the spherically symmetric potential
     Vr(:,:)=0.d0
     do jatom=1,nat
        CALL read_V_vsp(Vr(:jri(jatom),jatom),jatom,jri(jatom))  ! This potential is in Rydbergs
     enddo
     Vr(:,:) = Vr(:,:)/2.d0    ! Converst from Rydbergs to Hartrees
  endif
  call mpi_bcast_V_vsp(Vr, NRAD,NAT)
  
  ! case.in1 contains information of how to choose linearization energies.
  ! The linearization energies E & Elo are determined here inside on the fly.
  CALL ReadSecondPart_of_in1(E,Elo,lapw,loor,ilo,nlo, Vr,ZZ,jri,mult,dx,r0,Efermi,aname,REL,5,6,21,nat,nrad,nloat,lmax,lomax,Qprint)

  !print *, 'Came out of ReadSecondPart_of_in1 and have E='
  !do jatom=1,nat
  !   do l=0,3
  !      print *, jatom, l, E(l+1,jatom)
  !   enddo
  !enddo
  ! call exit


  ! Compute how many localized orbitals appear before each atom nlov(jatom)
  allocate(nlov(nat+1))
  tnlov=0   ! all localized orbitals before jatom
  do jatom=1,nat
     nlov(jatom)=tnlov  ! number of local orbitals before this atom
     DO l = 0,lomax
        DO jlo=1,ilo(l,jatom)
           IF (Elo(l,jlo,jatom).LT.(995.d+0)) tnlov=tnlov+((2*l+1))*mult(jatom)
        ENDDO
     ENDDO
  enddo
  nlov(nat+1)=tnlov  ! # of all local orbitals
  
  ! The radial wave functions are computed, and all necessary radial integrals are computed and stored
  CALL ATPAR(NT,REL,NAT,LNSMAX,ZZ,force,lmmx,lomax,ngau,nrad,nslmax)
  ! Frees the memory for the radial wave functions, as they will no longer be needed.
  CALL END_RADIAL  
  ! Coefficients for local orbitals are computed such that the real local orbitals vanish at the MT-boundary
  DO jatom=1,NAT
     if (Qprint) then
        WRITE(6, '(/,10X,"LO PARAMETERS FOR JATOM=",I3," ",A10)') jatom, aname(jatom)
     endif
     do l=0,lomax
        do jlo=1,ilo(l,jatom)
           CALL ABC0(alo(l,jlo,jatom),blo(l,jlo,jatom),clo(l,jlo,jatom),lapw(l,jatom),p(l+1,jatom),dp(l+1,jatom),pe(l+1,jatom),dpe(l+1,jatom),pei(l+1,jatom),plo(l,jlo,jatom),dplo(l,jlo,jatom),pi12lo(l,jlo,jatom),pe12lo(l,jlo,jatom),Rmt(JATOM),l,jlo,Qprint)
        enddo
     enddo
  ENDDO
  
  CALL STOP_TIMER(time_atpar)

  if (Qprint) then
     WRITE (6,6030) READ_CPU_TIME(time_atpar)
     WRITE (6,6040) (RMT(J), J = 1, NAT)
     WRITE (6,6050) (V(J),   J = 1, NAT)
     WRITE (6,6060) (1./VOL)
     WRITE (6,6070) ALAT
  endif
  
  ! Reads the file handle for k-points and min,max energy, and maximum number of bands
  CALL ReadThirdPart_of_in1(ITAPE,E1n,E2n,nemax_it,sproj_limit,eshift_iter,5)
  e2n=e2n+efermi

  CALL findKlistLength(ITAPE,numkpt,ERRMSG,info)
  if (info.ne.0) then
     CALL OUTERR('findKlistLength',ERRMSG); return
  endif
  
  allocate( IPGR(numkpt), KNAME(numkpt), WEIGHT(numkpt) )
  allocate( K3(3,numkpt) )
  
  CALL readKlist(ITAPE,numkpt,KNAME,K3,WEIGHT,IPGR,E1,E2,ERRMSG,info)
  if (info.ne.0) then
     CALL OUTERR('readKlist',ERRMSG); return
  endif

  ! Set up energy cutoffs, L-cutoffs, and Kmax-cutoffs.
  !   set up energy range for eigenvalues
  EUp = -1.0D+0    ! cutoff energy for eigenvalues
  ELow =  0.0D+0   ! cutoff energy for eigenvalues
  IF (E1 .NE. E2) THEN  ! first set it to the values from case.klist if they are OK
     ELow=E1
     EUp=E2
  ENDIF
  IF (E1n .NE. E2n) THEN ! overwrite with the values from case.in1 if these are also OK
     ELow=E1n
     EUp=E2n
  ENDIF
  !
  !        set up Kmax for list of plane waves for hamilton-matrix
  RKM = RKM/minval(RMT)           ! Rmt*Kmax/min(Rmt) => maxK
  NT = min(NT,LMAX)               ! NT if maxL from case.in1, which should be smaller than maximum alloved
  !
  
  korig(1:3)=0.0d0
  if (dokorig) then ! maybe the center of k-points is not at (0,0,0)
     read (3,*,iostat=ios) korig(1:3)
     if (ios.eq.0) then
        K3(1,1:numkpt)=K3(1,1:numkpt)+korig(1)
        K3(2,1:numkpt)=K3(2,1:numkpt)+korig(2)
        K3(3,1:numkpt)=K3(3,1:numkpt)+korig(3)
     endif
  end if

  ! write energy parameter to vector-file
  eorb_ind=0.d0
  if(runorb) eorb_ind=-99.d0
  DO i=1,nat
     do j=1,lmax
        if(.not.lapw(j-1,i)) E(j,i)=E(j,i)+200.D+0
     enddo

     if (vector_para .or. myrank.eq.master) then
        WRITE(10) (E(j,i),j=1,lmax)
        WRITE(10) ((Elo(j,k,i),j=0,lomax),k=1,nloat_new)
        !keep compatibility with old scheme
        if(nloat_new.gt.3) then
           WRITE(11,'(1000(f12.5))') (E(j,i),j=1,lmax),eorb_ind
           WRITE(11,'(1000(f12.5))') ((Elo(j,k,i),j=0,lomax),k=1,nloat_new)
        else
           WRITE(11,'(100(f9.5))') (E(j,i),j=1,lmax),eorb_ind
           WRITE(11,'(100(f9.5))') ((Elo(j,k,i),j=0,lomax),k=1,nloat_new)
        endif
     endif
     do j=1,lmax
        if(.not.lapw(j-1,i)) E(j,i)=E(j,i)-200.D+0
     enddo
  ENDDO
  !

  CALL STOP_TIMER(time_g1)
  deallocate( ZZ )
  info=0
  RETURN
  !
6020 FORMAT(10X,'R-MT TIMES K-MAX IS',F5.2,/10X,'MAX L IS',I3,5X,'MAX L IN NONSPHERICAL MATRIXELEMENTS:',I3,/ ' NUMBER OF ATOMS IS',I5)
6030 FORMAT(40X,'CPTIME ATPAR  =',F8.1)
6040 FORMAT(///10X,'R-MT=',4F10.7/15X,4F10.7/15X,4F10.7)
6050 FORMAT(10X,' FRACTIONAL VOLUME WITHIN MT=',4F10.7,2(/39X,4F10.7))
6060 FORMAT(10X,' ONE/UNIT CELL VOLUME=',E16.9)
6070 FORMAT(/10X,'LATTICE CONSTANTS ARE:',3F9.5)
END SUBROUTINE read_and_initialize
