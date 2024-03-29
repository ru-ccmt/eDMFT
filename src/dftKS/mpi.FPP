!**************************************************!
!******* § Copyright by Kristjan Haule, 2016 ******!
!**************************************************!
#ifdef _REAL
#define Ttype REAL*8
#define tsize 8
#else
#define Ttype COMPLEX*16
#define tsize 16
#endif

MODULE mpi
!!! This module should contain everything connected with paralelezation 
!!! of the code with MPI

#ifdef _MPI
  include 'mpif.h'
#endif

  INTEGER :: myrank    ! processor ID
  INTEGER :: nprocs    ! # of all processors awailable
  INTEGER :: ierr      ! returned error code
  INTEGER :: master    ! # of master processor
  CHARACTER*3 :: cpuID ! number of cpu in string representation
  INTEGER, PARAMETER :: clen = 100      ! Length of characters in command-line
  CHARACTER*100, ALLOCATABLE :: argv(:) ! Command-line arguments
  INTEGER      :: nargs                 ! Number of command line arguments
  LOGICAL      :: vector_para           ! vector files will be created by each processor.
  LOGICAL      :: Qprint                ! Should we open stdout (6) output file
  !INTEGER      :: pr_proc, pr_procr
  INTEGER, ALLOCATABLE :: pr_procs(:)
  LOGICAL      :: fastFilesystem       ! Should each core print info?
CONTAINS

#ifdef _MPI
  ! What needs to be done for parallel job
  SUBROUTINE start_MPI
    IMPLICIT NONE
    INTEGER :: iargc ! external function gives number of command line arguments
    INTEGER :: j
    ! getting in contact with MPI
    CALL MPI_INIT( ierr )
    CALL MPI_COMM_SIZE( MPI_COMM_WORLD, nprocs, ierr)
    CALL MPI_COMM_RANK( MPI_COMM_WORLD, myrank, ierr)
!    PRINT *,'nprocs=',nprocs,'myrank =',myrank 
    master = 0
    write(cpuID,'(I3)') myrank
    ! Get command-line arguments
    IF (myrank .EQ. master) THEN
       nargs = iargc()
       IF (nargs .GT. 4) nargs = nargs-4  ! Seems that MPI adds 4 additional arguments which we
       ALLOCATE (argv(nargs))                                 ! wouldn't like to parse
       !WRITE(*,'(A,I2)') 'nargs=', nargs
       DO j=1,nargs
          CALL getarg(j, argv(j))
          !WRITE(*,'(A,A)') 'argi=', TRIM(argv(j))
       ENDDO
    ENDIF
    ! Send the number of arguments to other nodes
    CALL MPI_BCAST(nargs, 1, MPI_INTEGER, master, MPI_COMM_WORLD,ierr)
    IF (myrank .NE. master) THEN
       ALLOCATE (argv(nargs))  ! Only now we can allocate correct size of array
    ENDIF
    ! Send all arguments to other nodes
    CALL MPI_BCAST(argv, nargs*clen, MPI_CHARACTER, master, MPI_COMM_WORLD,ierr)

    fastFilesystem = .false. !.true.
    Qprint = (myrank .EQ. master) .OR. fastFilesystem
  END SUBROUTINE start_MPI

  SUBROUTINE stop_MPI
    CALL MPI_FINALIZE(ierr)
  ENDSUBROUTINE stop_MPI

  SUBROUTINE FilenameMPI(infout)
    CHARACTER(LEN=*) :: infout
    infout    = TRIM(infout)//"."//trim(ADJUSTL(cpuID))
  ENDSUBROUTINE FilenameMPI

  SUBROUTINE FilenameMPI2(infout)
    CHARACTER(LEN=*) :: infout
    infout    = TRIM(infout)//"_"//trim(ADJUSTL(cpuID))
  ENDSUBROUTINE FilenameMPI2


  subroutine mpi_bcast_V_vsp(Vr, NRAD,NAT)
    IMPLICIT NONE
    INTEGER, intent(in)   :: NRAD, NAT
    REAL*8, intent(inout) :: Vr(NRAD,NAT)
    !
    CALL MPI_BCAST(Vr, NRAD*NAT, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
    if (ierr.ne.0) print *, 'ERROR in MPI_BCAST mpi_bcast_V_vsp', ierr
  end subroutine mpi_bcast_V_vsp

  subroutine mpi_bcast_V_vns(Vlm,lmmax,lm_stored,LM,nrad,lmmx)
    IMPLICIT NONE
    REAL*8, intent(inout):: Vlm(nrad,1:lmmx)
    INTEGER, intent(inout):: lmmax, lm_stored, LM(2,lmmx)
    INTEGER,intent(in)    :: nrad, lmmx
    ! locals
    INTEGER :: itmp(2)
    !
    if (myrank.eq.master) then
       itmp(1) = lmmax
       itmp(2) = lm_stored
    endif
    CALL MPI_BCAST(itmp, 2, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
    if (ierr.ne.0) print *, 'ERROR in MPI_BCAST mpi_bcast 1', ierr
    if (myrank.ne.master) then
       lmmax = itmp(1)
       lm_stored = itmp(2)
    endif
    CALL MPI_BCAST(LM, 2*lm_stored, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
    if (ierr.ne.0) print *, 'ERROR in MPI_BCAST mpi_bcast 2', ierr
    CALL MPI_BCAST(Vlm, nrad*lm_stored, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
    if (ierr.ne.0) print *, 'ERROR in MPI_BCAST mpi_bcast 3', ierr
  end subroutine mpi_bcast_V_vns
  
  subroutine mpi_bcast_nkk(nkk)
    IMPLICIT NONE
    INTEGER, intent(inout) :: nkk
    CALL MPI_BCAST(nkk, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
    if (ierr.ne.0) print *, 'ERROR in MPI_BCAST mpi_bcast 4', ierr
  end subroutine mpi_bcast_nkk

  subroutine mpi_bcast_V_vns_interstitial(KPxyz, POTK, nkk)
    IMPLICIT NONE
    INTEGER, intent(in)    :: nkk
    INTEGER, intent(inout) :: KPxyz(3,nkk+1)
    Ttype,   intent(inout) :: POTK(nkk+1)
    ! locals
    CALL MPI_BCAST(KPxyz, 3*nkk, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
    if (ierr.ne.0) print *, 'ERROR in MPI_BCAST mpi_bcast 5', ierr
    CALL MPI_BCAST(POTK, nkk*tsize/8, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
    if (ierr.ne.0) print *, 'ERROR in MPI_BCAST mpi_bcast 6', ierr
  end subroutine mpi_bcast_V_vns_interstitial

  subroutine mpi_SendReceive(kpoint,NE,HSrows,Eigval,KZZ,Z,numkpt)
    IMPLICIT NONE
    INTEGER, intent(in)    :: kpoint, NE, HSrows, numkpt
    INTEGER, POINTER,   intent(inout) :: KZZ(:,:)   ! (3,HSrows+1)
    Ttype, allocatable, intent(inout) :: Z(:,:)     ! (HSrows,NE)
    REAL*8, allocatable,intent(inout) :: Eigval(:)  ! Eigval(NE)
    ! locals
    INTEGER :: itmp(3), kpt, HSr, NEi, ipr
    INTEGER :: status(MPI_STATUS_SIZE)
    !
    itmp(1) = kpoint
    itmp(2) = NE
    itmp(3) = HSrows

    !if (myrank.eq.master) print *, 'On master we have', itmp(1:3)
    do ipr=1,nprocs-1 ! Notice that all processors need to send integer data, even if they contain no data.
       ! This is because otherwise master will still wait to receive data, even when tere is nothing to send.
       if (myrank.eq.master) then
          call MPI_Recv ( itmp, 3, MPI_INTEGER, ipr,      0, MPI_COMM_WORLD, status, ierr )
          if (ierr.ne.0) print *, 'ERROR in MPI_Recv 0', ierr
          kpt = itmp(1)
          NEi = itmp(2)
          HSr = itmp(3)
	  !print *, 'I am receiving', myrank, itmp(1:3)
          if (kpt.le.numkpt) then  ! Now even master knows that it might be empty k-point, and can skip it.
	     
	     if (HSr.gt.size(Z,1) .or. NEi.gt.size(Z,2)) then
                !print *, 'reallocating Z increasing', size(Z,1), HSr, size(Z,2), NEi
                deallocate(Z)
		allocate( Z(1:HSr,1:NEi) )
             end if
             if (HSr.gt.size(KZZ,2)) then
                !print *, 'reallocating KZZ increasing', size(KZZ,2), HSr
                deallocate( KZZ )
                allocate( KZZ(1:3,1:HSr) )
             endif
             if (NEi.gt.size(Eigval)) then
                !print *, 'reallocaying Eigval increasing ', size(Eigval), NEi
                deallocate(Eigval)
                allocate( Eigval(NEi) )
             endif
             !
             call MPI_Recv ( Eigval, NEi, MPI_DOUBLE_PRECISION, ipr,1, MPI_COMM_WORLD, status, ierr )
             if (ierr.ne.0) print *, 'ERROR in MPI_Recv 1', ierr
             call MPI_Recv ( KZZ, 3*HSr, MPI_INTEGER, ipr,      2, MPI_COMM_WORLD, status, ierr )
             if (ierr.ne.0) print *, 'ERROR in MPI_Recv 2', ierr
             call MPI_Recv ( Z, HSr*NEi*tsize/8, MPI_DOUBLE_PRECISION, ipr,    3, MPI_COMM_WORLD, status, ierr )
             if (ierr.ne.0) print *, 'ERROR in MPI_Recv 3', ierr
             !
             CALL PRTKPT(kpt,KZZ,Eigval,Z,HSr,NEi,Qprint)  ! Finally printing this k-point on master
          endif
       else if (myrank.eq.ipr) then
          !print *, 'I am sending', myrank, itmp(1:3)
          call MPI_Send ( itmp, 3, MPI_INTEGER, master, 0, MPI_COMM_WORLD, ierr )
          if (ierr.ne.0) print *, 'ERROR in MPI_Send 0', ierr
          if (kpoint.le.numkpt) then
             call MPI_Send ( Eigval, NE, MPI_DOUBLE_PRECISION, master, 1, MPI_COMM_WORLD,     ierr )
             if (ierr.ne.0) print *, 'ERROR in MPI_Send 1', ierr
             call MPI_Send ( KZZ, 3*HSrows, MPI_INTEGER, master,       2, MPI_COMM_WORLD,     ierr )
             if (ierr.ne.0) print *, 'ERROR in MPI_Send 2', ierr
             call MPI_Send ( Z, HSrows*NE*tsize/8, MPI_DOUBLE_PRECISION, master, 3, MPI_COMM_WORLD,     ierr )
             if (ierr.ne.0) print *, 'ERROR in MPI_Send 3', ierr
          endif
       endif
    enddo

  end subroutine mpi_SendReceive
  
#else

! What needs to be done for serial job

  SUBROUTINE start_MPI
    IMPLICIT NONE
    INTEGER :: iargc ! external function gives number of command line arguments
    INTEGER :: j
    myrank=0
    master=0
    nprocs=1
    ! Get command-line arguments
    nargs = iargc()
    ALLOCATE (argv(nargs))
    DO j=1,nargs
       CALL getarg(j, argv(j))
    ENDDO
    cpuID='0'
    Qprint = .true.
    !print *, 'nprocs=', nprocs
  END SUBROUTINE start_MPI

  SUBROUTINE stop_MPI
  ENDSUBROUTINE stop_MPI

  SUBROUTINE FilenameMPI(infout)
    CHARACTER(LEN=*) :: infout
  ENDSUBROUTINE FilenameMPI

  SUBROUTINE FilenameMPI2(infout)
    CHARACTER(LEN=*) :: infout
  ENDSUBROUTINE FilenameMPI2

  subroutine mpi_bcast_V_vsp(Vr, NRAD,NAT)
    IMPLICIT NONE
    INTEGER, intent(in)   :: NRAD, NAT
    REAL*8, intent(inout) :: VR(NRAD,NAT)
  end subroutine mpi_bcast_V_vsp
  
  subroutine mpi_bcast_V_vns(Vlm,lmmax,lm_stored,LM,nrad,lmmx)
    IMPLICIT NONE
    REAL*8, intent(inout):: Vlm(nrad,1:lmmx)
    INTEGER, intent(inout):: lmmax, lm_stored, LM(2,lmmx)
    INTEGER,intent(in)    :: nrad, lmmx
  end subroutine mpi_bcast_V_vns
  
  subroutine mpi_bcast_nkk(nkk)
    IMPLICIT NONE
    INTEGER, intent(inout) :: nkk
  end subroutine mpi_bcast_nkk

  subroutine mpi_bcast_V_vns_interstitial(KPxyz, POTK, nkk)
    IMPLICIT NONE
    INTEGER, intent(in)    :: nkk
    INTEGER, intent(inout) :: KPxyz(3,nkk+1)
    Ttype,   intent(inout) :: POTK(nkk+1)
  end subroutine mpi_bcast_V_vns_interstitial
  
  subroutine mpi_SendReceive(kpoint,NE,HSrows,Eigval,KZZ,Z,numkpt)
    IMPLICIT NONE
    INTEGER, intent(in)    :: kpoint, NE, HSrows, numkpt
    INTEGER, intent(inout) :: KZZ(3,HSrows+1)
    Ttype,   intent(inout) :: Z(HSrows,NE)
    REAL*8,  intent(inout) :: Eigval(NE)
  end subroutine mpi_SendReceive
  
#endif

subroutine WriteProcesses(pr_procr,first_kp,filename_vector,filename_energy)
  IMPLICIT NONE
  INTEGER,       intent(in) :: pr_procr, first_kp
  CHARACTER*180, intent(in) :: filename_vector, filename_energy
  ! locals
  INTEGER       :: i
  CHARACTER*180 :: filename_vectordn, filename_energydn
  INTEGER       :: idat(3)
  CHARACTER*180 :: sdat(2)
  INTEGER,      allocatable :: idat_all(:,:)
  CHARACTER*180,allocatable :: sdat_all(:,:)

  idat(1)=myrank
  idat(2)=pr_procr
  idat(3)=first_kp-1
  sdat(1)=filename_vector
  sdat(2)=filename_energy
  if (myrank.eq.master) then
     allocate( idat_all(3,nprocs), sdat_all(2,nprocs) )
  endif

#ifdef _MPI     
  call MPI_GATHER (idat,3,MPI_INTEGER,idat_all,3,MPI_INTEGER,master,MPI_COMM_WORLD, ierr)
  call MPI_GATHER (sdat,2*180,MPI_CHAR,sdat_all,2*180,MPI_CHAR,master,MPI_COMM_WORLD, ierr)
#else
  idat_all(:,1) = idat(:)
  sdat_all(:,1) = sdat(:)
#endif
  if (myrank.eq.master) then
     open(999,FILE='_processes_',status='unknown')
     do i=1,nprocs
        idat(:)=idat_all(:,i)
        sdat(:)=sdat_all(:,i)
        call give_dn_filename(filename_vectordn,sdat(1))
        call give_dn_filename(filename_energydn,sdat(2))
        WRITE(999,'(I6,1x,I6,1x,I6,1x,A,A,A,1x,A,A,A,1x,A,A,A,1x,A,A,A)') idat(1),idat(2),idat(3),'"',TRIM(sdat(1)),'"','"',TRIM(filename_vectordn),'"','"',TRIM(filename_energydn),'"','"',TRIM(sdat(2)),'"'
     enddo
     close(999)
  endif
end subroutine WriteProcesses
  
end MODULE mpi

SUBROUTINE FilenameMPI3(infout, ii)
  IMPLICIT NONE
  CHARACTER(LEN=*) :: infout
  INTEGER, intent(in) :: ii
  CHARACTER*10 :: strii
  write(strii,'(I10)') ii
  infout    = TRIM(infout)//trim(ADJUSTL(strii))
END SUBROUTINE FilenameMPI3

subroutine give_dn_filename(filenamedn,filename)
  IMPLICIT NONE
  CHARACTER*180, intent(in) :: filename
  CHARACTER*180, intent(out) :: filenamedn
  ! locals
  INTEGER :: i, ln
  ln=len_trim(filename)
  DO i=len_trim(filename),1,-1
     IF (filename(i:i) .EQ. '_') THEN
        ln = i-1
        exit
     ENDIF
  ENDDO
  if (filename(ln-1:ln) .EQ. 'dn') then
     filenamedn(:) = filename(:ln-2)//filename(ln+1:)
  else
     filenamedn(:) = filename(1:ln)//'dn'//filename(ln+1:180-2)
  endif
end subroutine give_dn_filename

subroutine add_xend(filename,xend)
  IMPLICIT NONE
  CHARACTER*180, intent(inout) :: filename
  CHARACTER*180, intent(in) :: xend
  ! locals                                                                                                                                                                  
  INTEGER :: ln
  !print *, 'add_xend: input filename=', TRIM(filename)
  ln=len_trim(filename)
  if (filename(ln-1:ln).eq.'_x') then
     filename(:) = filename(1:ln-2)//xend
  else
     filename(:) = filename(1:ln)//xend
  endif
  !print *, 'add_xend: xend=', TRIM(xend)
  !print *, 'add_xend: result =', filename
end subroutine add_xend


