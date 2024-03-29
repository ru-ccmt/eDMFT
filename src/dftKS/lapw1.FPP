! @Copyright 2007 Kristjan Haule 

#ifdef _REAL
#define tpcast(X) dble(X)  !! typecasting in fortran
#define XSWAP DSWAP
#define Ttype REAL*8
#else
#define tpcast(X) X
#define XSWAP ZSWAP
#define Ttype COMPLEX*16
#endif

PROGRAM LAPW1
  use lapw_timer, only : init_all_timer, start_timer, time_total, stop_timer, time_g1, read_wall_time, read_cpu_time, time_setwar, time_pkp, READ_CPU_TIME, time_hamilt, time_hns, time_diag, time_horb,START_TIMER, STOP_TIMER
  use cut1, only  : CUT
  use kpts, only    : K3
  use matrices, only : HS, HSDIAG, HSROWS, Z, EIGVAL, init_matrices, END_MATRICES, END_MATRICES_PART, KZZ! , Kn, RK
  use lolog, only   : nlo
  use orb, only     : NMOD
  use atspdt, only  : LQIND, LQNS, DP, DPE, E, P, PE, PEI, VNS1, VNS2, VNS3, GFAC
  use comi, only    : NVAA, NAT, NT, LNSMAX, NE
  use coml, only    : NOHNS
  use comr, only    : RKM, ELow, EUp
  use lstapw, only  : NV
  use structure, only: POS, RMT, V, VOL, MULT, ROTLOC, rotij, DeallocateStructure
  use out, only: end_out, WARP, kmax, NKK
  use nmr, only: writeham, dokorig
  use totpot, only  : KPxyz, POTK
  use structure, only: iord
  use readPotential, ONLY: read_V_vns_interstital_nk, read_V_vns_interstital_pot_real, read_V_vns_interstital_pot_cmplx
  use param, ONLY: nume, nume2, filename_energy, filename_vector, nmatmax
  use mpi, ONLY: Qprint, start_MPI, stop_MPI, nprocs, myrank, master, vector_para, mpi_bcast_nkk, mpi_bcast_V_vns_interstitial, mpi_SendReceive, WriteProcesses
  IMPLICIT NONE
  CHARACTER*5        MNAM
  PARAMETER          (MNAM = 'LAPW1')
  !        Local Scalars
  logical            nmat_only
  INTEGER            NUMKPT, INFO, KPOINT, numemin, nmatmax1, first_kp, pr_proc, pr_procr, ikp, i
  CHARACTER*80       DEFFN, ERRFN
  REAL*8             r1, r2, r3
  !
  INTEGER, allocatable  :: KP_equiv(:,:), kpoints(:)
  Ttype,   allocatable  :: taup(:)
  
  CALL start_MPI()

  cut=.false.
  
  CALL GTFNAM(DEFFN,ERRFN)
  CALL init_all_timer
  CALL START_TIMER(time_total)
  CALL ERRFLG(ERRFN,'Error in LAPW1')
  !
  CALL read_and_initialize(NUMKPT,DEFFN,INFO,nmat_only)
  IF (INFO .NE. 0) THEN
     CALL OUTERR(MNAM,'INILPW aborted unsuccessfully.')
     GOTO 999
  ENDIF

  if (nprocs.EQ.1 .and. Qprint) then
     write(6,*)'Running DMFT-lapw1 in single processor mode'
     write(6,*)' '
  else if (myrank.EQ.master .and. nprocs.gt.1) then
     write(6,*)'Running DMFT-lapw1 in mpi mode'
     write(6,*)' '
  endif

  CALL START_TIMER(time_setwar)

  if (myrank.eq.master) then
     ! Finds how many kreciprocal vectors are listed in case.vns for interstital
     NKK = read_V_vns_interstital_nk()
  endif
  call mpi_bcast_nkk(nkk)
  ! Now allocate arrays to hold interstitial potential
  allocate( KPxyz(3,NKK+1), POTK(NKK+1) )
  ! Read the interstitial potential from case.vns
  if (myrank.eq.master) then
#ifdef _REAL  
     CALL read_V_vns_interstital_pot_real(nkk,KPxyz,POTK)
#else
     CALL read_V_vns_interstital_pot_cmplx(nkk,KPxyz,POTK)
#endif
  endif
  call mpi_bcast_V_vns_interstitial(KPxyz, POTK, nkk)
  
  ! Finds maximum value of G in each direction kmax=(Gx,Gy,Gz)
  allocate(KP_equiv(3,iord), taup(iord))
  CALL HashTableSize(kmax, KP_equiv, taup)
  ! Allocated the hash-table to hold V(G1)
  allocate(warp(-kmax(1):kmax(1),-kmax(2):kmax(2),-kmax(3):kmax(3)))
  if (Qprint) write(6,'(a20,f12.1,a,a,2i6)') 'alloc Warp',(kmax(1)*kmax(2)*kmax(3))*128.0/(1024.0*1024.0),' MB '
  ! Fills in the hash table for all G-vectors, even those that are generated with any group operation.
  CALL SetHashTable(WARP, kmax(1), kmax(2), kmax(3), KP_equiv, taup)
  ! No need to keep these beyond this point
  deallocate(KP_equiv, taup)
  deallocate( KPxyz, POTK )  !     This is the only routine which uses potk, so deallocate it 

  CALL STOP_TIMER(time_setwar)

  pr_proc  = floor(numkpt/DBLE(nprocs)+0.999)  ! The maximum number of points calculated per processor     
  ! correcting pr_proc is the number of k-points is not divisible
  pr_procr = pr_proc
  first_kp = pr_proc*myrank+1

  if ((myrank+1)*pr_proc > numkpt) then
     if (numkpt-myrank*pr_proc > 0) then
        pr_procr = numkpt-myrank*pr_proc
     else
        pr_procr = 0
     endif
  endif

  allocate( kpoints(pr_proc) )
  if (vector_para) then
     ! If we output into separate files, we can compute first pr_proc points by the first core, next pr_proc by the second, etc...
     do i=0,pr_proc-1
        kpoints(i+1) = first_kp+i
     enddo
     if (Qprint) WRITE(6,'(A,I3,1x,A,I4,1x,A,I4,1x,A,I4)') 'first_kp=', first_kp, 'pr_procr=', pr_procr, 'pr_proc=', pr_proc, 'tot-k=', numkpt
     !WRITE(6,*) kpoints(:)
  else
     ! If we output into a common file, the second core needs to calculate kpoint number 2, and third processor kpoint 3, etc...
     do i=0,pr_proc-1
        kpoints(i+1) = 1 + myrank + i*nprocs
     enddo
     if (Qprint) WRITE(6,'(A,I4,1x,A,I4)') 'pr_proc=', pr_proc, 'tot-k=', numkpt
     !WRITE(6,*) kpoints(:)
  endif
  

  !DO kpoint=first_kp,first_kp+pr_proc-1
  DO ikp=1,pr_proc
     kpoint = kpoints(ikp)
     ! We could use pr_procr, but we rather use pr_prc due to mpi communication. All processes have to be enter equal number of times so that we can call MPI_Barrier. The work is done only if needed.
     if (kpoint .le. numkpt) then 
        !call INIT_MATRICES(0, 0,nmatmax1, Qprint)
        nmatmax1=nmatmax
#ifndef _REAL
        nmatmax1=int(nmatmax1/sqrt(2.d0))
#endif
        !  set up array warp and generate plane wave list
        call generate_reciprocal(RKM,K3(:,kpoint),nat,nmatmax1,nmat_only)
        if(nmat_only) then
           CALL ERRCLR(ERRFN)
           call stop_MPI
           stop 'NMAT_ONLY'
        endif
        numemin=min(nume,(nv+nlo))
        if(dokorig) numemin=nv+nlo     ! fix for NMR
        call INIT_MATRICES(nv+nlo, numemin,Qprint)
        if (nlo.GT.0) CALL LOPW(nat)
        NVAA = NV
        !  construct secular equation
        CALL HAMILTON(NT,NAT,NV,P,PE,DP,DPE,MULT,POS,V,Vol,PEI,E)
        !  Diagonal elements of Hamilt in HSDIAG, of Overlap in HS
        ! Swap matrix diagonal so that HSDIAG contains diagonal part of overlap, while HS(:,:) contains Hamiltonian on the lower part and diagonal
        !       This is needed for HNS and HORB
        CALL XSWAP(HSROWS,HS,HSROWS+1,HSDIAG,1)
        !  Diagonal elements of Hamilt in HS, of Overlap in HSDIAG
        if(.not.nohns) CALL HNS(NV,NAT,LNSMAX+1,POS,Vol,RMT,MULT,LQIND,LQNS,GFAC,ROTLOC,ROTIJ,VNS1(1,1,1,1),VNS2(1,1,1,1),VNS3(1,1,1,1))
        ! Add orbital potential to the hamiltonian
        if (Nmod.ne.0)  CALL  horb(NV,NAT,POS,Vol,RMT,MULT,PEI,ROTLOC,ROTIJ)
        !
        call END_MATRICES_PART
        if (numkpt .eq. 1) call end_out
        if (writeham) call Write_Ham(hs,hsdiag,hsrows)
        ! Swap back so that HSDIAG contains diagonal part of Hamiltonian, while HS(:,:) contains Overlap on the upper part and diagonal
        CALL XSWAP(HSROWS,HS,HSROWS+1,HSDIAG,1)
        !.....full diagonalization
        CALL START_TIMER(time_diag)
        CALL Diagonalize(Elow, Eup, NE, HS, HSDIAG, Z, EIGVAL, HSROWS, nume2, Qprint)
        CALL STOP_TIMER(time_diag)
        if (Qprint) call print_times()
        CALL PRTKPT(kpoint, KZZ, Eigval, Z, HSrows, NE, Qprint)
     else
        NE=0
        HSrows=0
	!allocate( KZZ(3,1),  eigval(1), Z(1,1) )
	!allocate( HS(1,1), HSDIAG(1) )
     endif

     ! If the output is not parallel (each processor does not produce its own vector/energy file), we always need to
     ! send data to the master node, and the master node needs to print

     if (.not. vector_para) call mpi_SendReceive(kpoint,NE,HSrows,Eigval,KZZ,Z,numkpt)
     if (kpoint .le. numkpt) CALL END_MATRICES
  ENDDO


  if (vector_para) call WriteProcesses(pr_procr,first_kp,filename_vector,filename_energy)
  
  deallocate( kpoints )
  call DeallocateStructure()
  if (Qprint) then
     !        print summary
     WRITE(21,'(":KPT   :      NUMBER OF K-POINTS:",i6)') numkpt
     WRITE(6,*)  '      NUMBER OF K-POINTS:', numkpt
     if(cut) then
        WRITE(21,99) 
99      format(':WARN :      WARNING: RKmax reduced due to NMATMAX')
        WRITE(6,*)  '      WARNING: RKmax reduced due to NMATMAX'
     endif
     if(ne.eq.NUME) then
        WRITE(21,98) 
98      format(':WARN :      WARNING: NE limited by NUME in param.inc')
        WRITE(6,*)  '      WARNING: NE limited by NUME in param.inc'
     endif
     !
     !        timing results
     r1= READ_CPU_TIME(time_setwar) +  READ_CPU_TIME(time_g1) + READ_CPU_TIME(time_hamilt) + READ_CPU_TIME(time_hns) + READ_CPU_TIME(time_diag) + READ_CPU_TIME(time_pkp)
     r2 = READ_CPU_TIME(time_setwar) + READ_CPU_TIME(time_g1)
     r3 = READ_CPU_TIME(time_hamilt) + READ_CPU_TIME(time_hns) + READ_CPU_TIME(time_diag) +  READ_CPU_TIME(time_pkp)  
     Write(6,6000) r1, r2, r3
6000 FORMAT(3X,'===> TOTAL CPU       TIME: ',F8.1, ' (INIT = ',F8.1,' + K-POINTS = ',F8.1, ')' )
  
     r1 = READ_WALL_TIME(time_setwar) +   READ_WALL_TIME(time_g1) + READ_WALL_TIME(time_hamilt) +   READ_WALL_TIME(time_hns) + READ_WALL_TIME(time_diag) +  READ_WALL_TIME(time_pkp) 
     r2 = READ_WALL_TIME(time_setwar) +  READ_WALL_TIME(time_g1)
     r3 =  READ_WALL_TIME(time_hamilt) +  READ_WALL_TIME(time_hns) + READ_WALL_TIME(time_diag) +  READ_WALL_TIME(time_pkp)
     Write(6,6010)r1,r2,r3
6010 FORMAT(3X,   '> SUM OF WALL CLOCK TIMES: ',F8.1,' (INIT = ',F8.1, ' + K-POINTS = ',F8.1,')')
  endif
  !
  CALL STOP_TIMER(time_total)
  if (Qprint) then
     WRITE(6,*) '     Maximum WALL clock time: ',READ_WALL_TIME(time_total)
     WRITE(6,*) '     Maximum CPU time:        ',READ_CPU_TIME(time_total)
  endif
  CALL ERRCLR(ERRFN)

  CALL stop_MPI()
  
  STOP !' LAPW1 END'
  
999 CONTINUE
  call stop_MPI
  STOP 'LAPW1 - Error'
  
END PROGRAM LAPW1


subroutine Write_Ham(hs,hsdiag,hsrows)
  IMPLICIT NONE
  integer, intent(in) :: hsrows
  Ttype,   intent(in) :: hsdiag(hsrows), hs(hsrows,hsrows)
  ! locals
  INTEGER :: i,j
  write(12) hsrows
  do i=1,hsrows
     write(12) (hs(j,i),j=1,i-1),hsdiag(i)
  enddo
  do i=1,hsrows
     write(12) (hs(j,i),j=i,hsrows)
  enddo
end subroutine Write_Ham

subroutine print_times()
  use lapw_timer, only: READ_CPU_TIME, READ_WALL_TIME, time_hamilt, time_hns, time_diag, time_horb
  IMPLICIT NONE
  REAL*8 :: rr1, rr2, rr3,rr4
  rr1 = READ_CPU_TIME(time_hamilt)
  rr2 = READ_CPU_TIME(time_hns)
  rr3 =  READ_CPU_TIME(time_horb)
  rr4 =  READ_CPU_TIME(time_diag)
  WRITE(6,6020) rr1, rr2, rr3,rr4
  rr1 = READ_WALL_TIME(time_hamilt)
  rr2 = READ_WALL_TIME(time_hns)
  rr3 =  READ_WALL_TIME(time_horb)
  rr4 =  READ_WALL_TIME(time_diag)
  WRITE(6,6021) rr1, rr2, rr3,rr4
6020 FORMAT(7X,'TIME HAMILT (CPU)  =',F8.1,', HNS =',F8.1,', HORB =',F8.1,', DIAG =',F8.1)
6021 FORMAT(7X,'TIME HAMILT (WALL) =',F8.1,', HNS =',F8.1,', HORB =',F8.1,', DIAG =',F8.1)
end subroutine print_times

