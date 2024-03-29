!**************************************************!
!******* § Copyright by Kristjan Haule, 2002 ******!
!**************************************************!
MODULE com_mpi
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
    !PRINT *,'nprocs=',nprocs,'myrank =',myrank 
    master = 0
    write(cpuID,'(I3)') myrank
    ! Get command-line arguments
    IF (myrank .EQ. master) THEN
       nargs = iargc()
       IF (nargs .GT. 4) nargs = nargs-4  ! Seems that MPI adds 4 additional arguments which we
       ALLOCATE (argv(nargs))                                 ! wouldn't like to parse
       WRITE(*,'(A,I2)') 'nargs=', nargs
       DO j=1,nargs
          CALL getarg(j, argv(j))
          WRITE(*,'(A,A)') 'argi=', TRIM(argv(j))
       ENDDO
    ENDIF
    ! Send the number of arguments to other nodes
    CALL MPI_BCAST(nargs, 1, MPI_INTEGER, master, MPI_COMM_WORLD,ierr)
    IF (myrank .NE. master) THEN
       ALLOCATE (argv(nargs))  ! Only now we can allocate correct size of array
    ENDIF
    ! Send all arguments to other nodes
    CALL MPI_BCAST(argv, nargs*clen, MPI_CHARACTER, master, MPI_COMM_WORLD,ierr)
  END SUBROUTINE start_MPI

  SUBROUTINE stop_MPI
    CALL MPI_FINALIZE(ierr)
  ENDSUBROUTINE stop_MPI

  SUBROUTINE FilenameMPI(infout)
    CHARACTER(LEN=*) :: infout
    infout    = TRIM(infout)//"."//trim(ADJUSTL(cpuID))
  ENDSUBROUTINE FilenameMPI

  SUBROUTINE FindMax_MPI(max_bands, nbandsk, pr_proc)
    INTEGER, intent(out) :: max_bands
    INTEGER, intent(in)  :: pr_proc
    INTEGER, intent(in)  :: nbandsk(pr_proc)
    ! locals
    INTEGER :: maxb
    maxb=1
    DO i=1,pr_proc
       maxb = max(maxb,nbandsk(i))
    ENDDO
    CALL MPI_ALLREDUCE(maxb, max_bands, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
  END SUBROUTINE FindMax_MPI
  
  SUBROUTINE Reduce_MPI(conduc, COmega, Nw, Ndirection, renorm_wgh, Nd) 
    IMPLICIT NONE
    REAL*8, intent(inout) :: conduc(Nw,Ndirection), COmega(2*Nd+1), renorm_wgh
    INTEGER, intent(in)   :: Nw, Ndirection, Nd
    ! locals
    REAL*8 :: wconduc(Nw,Ndirection)
    REAL*8 :: wrenorm_wgh
    REAL*8 :: wCOmega(2*Nd+1)

    !print *, '*Nw=', Nw, 'myrank=', myrank
    !print *, '*renorm_wgh=', renorm_wgh
    !print *, '*myrank=', myrank, '*conduc=', conduc

    CALL MPI_REDUCE(renorm_wgh, wrenorm_wgh,  1, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_WORLD, ierr)
    if (ierr.ne.0) print *, 'ERROR in MPI_REDUCE 1', ierr

    if (myrank.eq.master) print *, '**wrenorm_wgh=', wrenorm_wgh

    CALL MPI_REDUCE(conduc,  wconduc, Nw*Ndirection, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_WORLD, ierr)
    if (ierr.ne.0) print *, 'ERROR in MPI_REDUCE 2', ierr
    
    if (myrank.eq.master) print *, '**wconduc=', wconduc

    CALL MPI_REDUCE(COmega,  wCOmega, 2*Nd+1, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_WORLD, ierr)
    if (ierr.ne.0) print *, 'ERROR in MPI_REDUCE 3', ierr
    
    if (myrank.eq.master) then
       renorm_wgh = wrenorm_wgh
       conduc(:,:) = wconduc(:,:)
       COmega(:) = wCOmega(:)
    endif
  END SUBROUTINE Reduce_MPI
  
  SUBROUTINE Reduce_MPI_dos(gc,Nw)
    IMPLICIT NONE
    COMPLEX*16, intent(inout) :: gc(2*Nw)
    INTEGER, intent(in)       :: Nw
    ! locals
    COMPLEX*16 :: wgc(2*Nw)
    INTEGER    :: ierr
    CALL MPI_REDUCE(gc,  wgc, 2*2*Nw, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_WORLD, ierr)
    if (ierr.ne.0) print *, 'ERROR in MPI_REDUCE 5', ierr
    if (myrank.eq.master) then
       gc(:) = wgc(:)
    endif
  END SUBROUTINE Reduce_MPI_dos
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
    print *, 'Non parallel calculation'
  END SUBROUTINE start_MPI

  SUBROUTINE stop_MPI
  ENDSUBROUTINE stop_MPI

  SUBROUTINE FilenameMPI(infout)
    CHARACTER(LEN=*) :: infout
  ENDSUBROUTINE FilenameMPI

  SUBROUTINE FindMax_MPI(max_bands, nbandsk, pr_proc)
    INTEGER, intent(out) :: max_bands
    INTEGER, intent(in)  :: pr_proc
    INTEGER, intent(in)  :: nbandsk(pr_proc)
    ! locals
    INTEGER :: maxb
    maxb=1
    DO i=1,pr_proc
       maxb = max(maxb,nbandsk(i))
    ENDDO
    max_bands = maxb
  END SUBROUTINE FindMax_MPI

  SUBROUTINE Reduce_MPI(conduc, COmega, Nw, Ndirection, renorm_wgh, Nd) 
    IMPLICIT NONE
    REAL*8, intent(inout) :: conduc(Nw,Ndirection), renorm_wgh, COmega(2*Nd+1)
    INTEGER, intent(in)   :: Nw, Ndirection, Nd
    ! locals
  END SUBROUTINE Reduce_MPI

  SUBROUTINE Reduce_MPI_dos(gc,Nw)
    IMPLICIT NONE
    COMPLEX*16, intent(inout) :: gc(2*Nw)
    INTEGER, intent(in)       :: Nw
    ! locals
  END SUBROUTINE Reduce_MPI_dos
#endif
END MODULE com_mpi
