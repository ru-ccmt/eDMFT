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
  LOGICAL      :: vector_para           ! is this parallel Wien2K run?
  LOGICAL      :: Qprint
  INTEGER      :: vectors(20,3)        ! which vector files should be read
  INTEGER      :: nvector
  CHARACTER*200:: VECFN(4)
  CHARACTER*200:: fvectors(20,4)
  INTEGER      :: pr_proc, pr_procr
  INTEGER, ALLOCATABLE :: pr_procs(:)
  LOGICAL,parameter :: fastFilesystem=.FALSE.
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
       WRITE(*,'(A,I2)') 'nargs=', nargs
       DO j=1,nargs
          CALL getarg(j, argv(j))
          WRITE(*,'(A,A)') 'argi=', TRIM(argv(j))
       ENDDO
    ENDIF
    ! Send the number of arguments to other nodes
    CALL MPI_BCAST(nargs, 1, MPI_INTEGER, master, MPI_COMM_WORLD,ierr)
    if (ierr.ne.0) WRITE(*,*) 'ERROR in MPI_BCAST 1', ierr
    IF (myrank .NE. master) THEN
       ALLOCATE (argv(nargs))  ! Only now we can allocate correct size of array
    ENDIF
    ! Send all arguments to other nodes
    CALL MPI_BCAST(argv, nargs*clen, MPI_CHARACTER, master, MPI_COMM_WORLD,ierr)
    if (ierr.ne.0) WRITE(*,*) 'ERROR in MPI_BCAST 2', ierr
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
    infout    = TRIM(infout)//trim(ADJUSTL(cpuID))
  ENDSUBROUTINE FilenameMPI2

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
    !print *, 'max_bands on', myrank, '=', maxb
    CALL MPI_ALLREDUCE(maxb, max_bands, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
    !print *, 'total number of maxbands=', max_bands
    ! MPI_ALLREDUCE(sendbuf, recbuf, count, MPI_TYPE, MPI_OP, comm, ierr)
  END SUBROUTINE FindMax_MPI

  SUBROUTINE FindMaxK_MPI(kmax)
    INTEGER, intent(inout) :: kmax(3)
    INTEGER :: tkmax(3)
    CALL MPI_ALLREDUCE(kmax, tkmax, 3, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
    kmax(:) = tkmax(:)
  END SUBROUTINE FindMaxK_MPI
  
  SUBROUTINE Reduce_MPI0(Olapm, maxdim, ncix, max_nbands)
    IMPLICIT NONE
    COMPLEX*16, intent(inout) :: Olapm(maxdim,maxdim,ncix)
    INTEGER, intent(inout) :: max_nbands
    INTEGER, intent(in) :: maxdim, ncix
    !
    COMPLEX*16 :: w_Olapm(maxdim,maxdim,ncix)
    INTEGER :: w_max_nbands
    
    CALL MPI_ALLREDUCE(Olapm, w_Olapm, maxdim*maxdim*ncix, MPI_DOUBLE_COMPLEX,  MPI_SUM, MPI_COMM_WORLD, ierr)
    if (ierr.ne.0) print *, 'ERROR in MPI_REDUCE 1', ierr
    Olapm(:,:,:) = w_Olapm(:,:,:)

    CALL MPI_ALLREDUCE(max_nbands, w_max_nbands, 1, MPI_INTEGER,  MPI_MAX, MPI_COMM_WORLD, ierr)
    if (ierr.ne.0) print *, 'ERROR in MPI_REDUCE 2', ierr
    max_nbands = w_max_nbands
  END SUBROUTINE Reduce_MPI0

  SUBROUTINE Reduce_maxbands(max_nbands)
    IMPLICIT NONE
    INTEGER, intent(inout) :: max_nbands
    !
    INTEGER :: w_max_nbands
    CALL MPI_ALLREDUCE(max_nbands, w_max_nbands, 1, MPI_INTEGER,  MPI_MAX, MPI_COMM_WORLD, ierr)
    if (ierr.ne.0) print *, 'ERROR in MPI_REDUCE maxbands', ierr
    max_nbands = w_max_nbands
  END SUBROUTINE Reduce_maxbands

  SUBROUTINE Reduce_maxnume(max_nume, nume)
    IMPLICIT NONE
    INTEGER, intent(out) :: max_nume
    INTEGER, intent(in)  :: nume
    !
    CALL MPI_ALLREDUCE(nume, max_nume, 1, MPI_INTEGER,  MPI_MAX, MPI_COMM_WORLD, ierr)
    if (ierr.ne.0) print *, 'ERROR in MPI_REDUCE maxnume', ierr
  END SUBROUTINE Reduce_maxnume
  
  SUBROUTINE Reduce_MPI2(numkpt, nkpt, sumw0, DM_EMIN, DM_EMAX)
    IMPLICIT NONE
    INTEGER, intent(out)  :: numkpt
    INTEGER, intent(in)   :: nkpt
    REAL*8, intent(inout) :: sumw0
    REAL*8, intent(inout) :: DM_EMIN, DM_EMAX
    !
    REAL*8  :: sumw0_final
    REAL*8  :: w_DM_EMIN, w_DM_EMAX
    
    if (vector_para) then
        CALL MPI_ALLREDUCE(nkpt, numkpt, 1, MPI_INTEGER,  MPI_SUM, MPI_COMM_WORLD, ierr)
	if (ierr.ne.0) print *, 'ERROR in MPI_ALLREDUCE 0', ierr
    
        CALL MPI_ALLREDUCE(sumw0, sumw0_final, 1, MPI_REAL8,  MPI_SUM, MPI_COMM_WORLD, ierr)
	if (ierr.ne.0) print *, 'ERROR in MPI_ALLREDUCE 1', ierr
	sumw0 = sumw0_final

	CALL MPI_ALLREDUCE(DM_EMIN, w_DM_EMIN, 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, ierr)
	if (ierr.ne.0) print *, 'ERROR in MPI_ALLREDUCE 2', ierr
	DM_EMIN = w_DM_EMIN
	
	CALL MPI_ALLREDUCE(DM_EMAX, w_DM_EMAX, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
	if (ierr.ne.0) print *, 'ERROR in MPI_ALLREDUCE 3', ierr
	DM_EMAX = w_DM_EMAX
    else
        numkpt=nkpt
    endif

  END SUBROUTINE Reduce_MPI2

  SUBROUTINE Gather_procs(pr_procr, pr_procs, nprocs)
    IMPLICIT NONE
    INTEGER, intent(in) :: pr_procr, nprocs
    INTEGER, intent(out):: pr_procs(nprocs)
    INTEGER :: ier
    CALL MPI_GATHER(pr_procr, 1, MPI_INTEGER, pr_procs, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
    if (ierr.ne.0) print *, 'ERROR in MPI_GATHER pr_procs', ierr
  END SUBROUTINE Gather_procs


  SUBROUTINE Reduce_MPI1(pr_procr, nomega, max_nume)
    use muzero, ONLY: nkp, wgh, nemm, zEk, Ek, max_nbands
    use param, ONLY: numkpt
    IMPLICIT NONE
    INTEGER, intent(in)   :: pr_procr, nomega, max_nume
    !
    INTEGER, allocatable    :: w_nemm(:,:), offs(:), offsl(:), pr_procsl(:)
    COMPLEX*16, allocatable :: w_zEk(:,:,:)
    REAL*8, allocatable     :: w_Ek(:,:), w_wgh(:)
    INTEGER :: i, ikp, iom
    ! MPI_Gatherv(*sendbuf, int sendcnt, MPI_Datatype sendtype, *recvbuf, int recvcnts[], int displs[], MPI_Datatype recvtype, int root, MPI_Comm comm)
    !WRITE(6,*) 'CAME into MPI1'
    
    if (myrank.eq.master)  nkp = sum(pr_procs)
    
    ALLOCATE( offs(nprocs) )
    offs(1)=0
    do i=2,nprocs
       offs(i)= offs(i-1)+pr_procs(i-1)
    enddo
    ALLOCATE( offsl(nprocs), pr_procsl(nprocs) )

    pr_procsl(:) = pr_procs(:)
    offsl(:) = offs(:)

    if (myrank.eq.master)  then
       ALLOCATE( w_wgh(nkp) )
       !write(6,*) 'Just allocated w_wgh of size', nkp
       !write(6,*) 'pr_procs=', pr_procs
    endif
    !write(6,*) 'pr_procr=', pr_procr
    
    CALL MPI_GATHERV(wgh, pr_procr, MPI_REAL8, w_wgh, pr_procsl, offsl, MPI_REAL8, master, MPI_COMM_WORLD, ierr)
    
    if (ierr.ne.0) print *, 'ERROR in MPI_GATHER 4', ierr
    if (myrank.eq.master) then
       DEALLOCATE( wgh )
       ALLOCATE( wgh(nkp) )
       wgh(:) = w_wgh(:)
       DEALLOCATE( w_wgh )
    endif

    if (myrank.eq.master) ALLOCATE( w_nemm(3,nkp) )
    
    pr_procsl(:) = pr_procs(:)*3
    offsl(:) = offs(:)*3
    CALL MPI_GATHERV(nemm, pr_procr*3, MPI_INTEGER, w_nemm, pr_procsl, offsl, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
    if (ierr.ne.0) print *, 'ERROR in MPI_GATHER 1', ierr
    if (myrank.eq.master) then
       DEALLOCATE( nemm )
       ALLOCATE( nemm(3,nkp) )
       nemm(:,:) = w_nemm(:,:)
       DEALLOCATE( w_nemm )
    endif
    
    if (myrank.eq.master) print *, 'Data-size pr_processor[MBy]=', 16.*max_nbands*nomega*pr_proc/(1024**2)
    if (myrank.eq.master) print *, 'Total datasize is[MBy]=', 16.*max_nbands*nomega*pr_proc*nprocs/(1024**2)
    if (myrank.eq.master) WRITE(*,'(A,I4,1x,A,I4,1x,A,I4,1x,A,I4,1x)') 'max_nbands=', max_nbands, 'nomega=', nomega, 'pr_proc=', pr_proc, 'nprocs=', nprocs
    
    if (myrank.eq.master) ALLOCATE( w_zEk(max_nbands,nomega,nkp) )
    pr_procsl(:) = pr_procs(:)*(max_nbands*nomega)
    offsl(:) = offs(:)*(max_nbands*nomega)
    CALL MPI_GATHERV(zEk, max_nbands*nomega*pr_procr, MPI_DOUBLE_COMPLEX, w_zEk, pr_procsl, offsl, MPI_DOUBLE_COMPLEX, master, MPI_COMM_WORLD, ierr)
    if (ierr.ne.0) print *, 'ERROR in MPI_GATHER 2', ierr
    if (myrank.eq.master) then
       DEALLOCATE( zEk )
       ALLOCATE( zEk(max_nbands,nomega,nkp) )
       zEk(:,:,:) = w_zEk(:,:,:)
       DEALLOCATE( w_zEk )
    endif
    
    if (myrank.eq.master) ALLOCATE( w_Ek(max_nume,nkp) )
    pr_procsl(:) = pr_procs(:)*(max_nume)
    offsl(:) = offs(:)*(max_nume)
    CALL MPI_GATHERV(Ek, max_nume*pr_procr, MPI_REAL8, w_Ek, pr_procsl, offsl, MPI_REAL8, master, MPI_COMM_WORLD, ierr)
    if (ierr.ne.0) print *, 'ERROR in MPI_GATHER 3', ierr
    if (myrank.eq.master) then
       DEALLOCATE( Ek )
       ALLOCATE( Ek(max_nume,nkp) )
       Ek(:,:nkp) = w_Ek(:,:nkp)
       DEALLOCATE( w_Ek )
    endif
    CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

    DEALLOCATE( offs, offsl, pr_procsl)
  END SUBROUTINE Reduce_MPI1


  SUBROUTINE Reduce_MPI1b(nomega, Ebmax, Ebmin, weib, gtot, numex)
    IMPLICIT NONE
    INTEGER, intent(in)   :: nomega, numex
    REAL*8, intent(inout) :: Ebmin(numex), Ebmax(numex), weib(numex)
    COMPLEX*16, intent(inout) :: gtot(nomega)
    !
    COMPLEX*16, allocatable :: w_gtot(:)
    REAL*8, allocatable     :: w_Ebmin(:), w_Ebmax(:), w_eib(:)
    INTEGER :: i, ikp, iom
    
    if (myrank.eq.master) ALLOCATE( w_gtot(nomega) )
    CALL MPI_REDUCE(gtot,w_gtot, nomega, MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, ierr)
    if (ierr.ne.0) print *, 'ERROR in MPI_REDUCE 4', ierr
    if (myrank.eq.master) then
       gtot(:) = w_gtot(:)
       DEALLOCATE( w_gtot )
    endif
    
    if (myrank.eq.master) ALLOCATE( w_Ebmax(numex) )
    CALL MPI_REDUCE(Ebmax,w_Ebmax, numex, MPI_DOUBLE_PRECISION, MPI_MAX, master, MPI_COMM_WORLD, ierr)
    if (ierr.ne.0) print *, 'ERROR in MPI_REDUCE 5', ierr
    if (myrank.eq.master) then
       Ebmax(:) = w_Ebmax(:)
       DEALLOCATE( w_Ebmax )
    endif
    
    if (myrank.eq.master) ALLOCATE( w_Ebmin(numex) )
    CALL MPI_REDUCE(Ebmin,w_Ebmin, numex, MPI_DOUBLE_PRECISION, MPI_MIN, master, MPI_COMM_WORLD, ierr)
    if (ierr.ne.0) print *, 'ERROR in MPI_REDUCE 6', ierr
    if (myrank.eq.master) then
       Ebmin(:) = w_Ebmin(:)
       DEALLOCATE( w_Ebmin )
    endif
    
    if (myrank.eq.master) ALLOCATE( w_eib(numex) )
    CALL MPI_REDUCE(weib,w_eib, numex, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_WORLD, ierr)
    if (ierr.ne.0) print *, 'ERROR in MPI_REDUCE 7', ierr
    if (myrank.eq.master) then
       weib(:) = w_eib(:)
       DEALLOCATE( w_eib )
    endif
    
  END SUBROUTINE Reduce_MPI1b



  
  SUBROUTINE Bcast_MPI(EF)
    IMPLICIT NONE
    REAL*8, intent(inout) :: EF
    CALL MPI_BCAST(EF, 1, MPI_REAL8, master, MPI_COMM_WORLD, ierr)
    if (ierr.ne.0) print *, 'ERROR in MPI_BCAST 1', ierr
  END SUBROUTINE Bcast_MPI

  SUBROUTINE Reduce_MPI(xwt, etot, ftot, ftot0, gloc, w_RHOLM, DM, Gdloc, fsph, fnsp, fsph2, Edimp, Nds, w_xwt1, w_xwteh, w_xwtel, w_xwt1h, w_xwt1l, sumfft, tloc, nomega, NRAD, LM_MAX, nat, natm, iff1, iff2, iff3, ift1, ift2, ift3, maxdim, ncix, ntcix, nip, Qforce)
    IMPLICIT NONE
    REAL*8, intent(inout)     :: xwt, etot, ftot, ftot0
    COMPLEX*16, intent(inout) :: gloc(nomega)
    REAL*8, intent(inout)     :: w_RHOLM(1:NRAD,1:LM_MAX,1:nat)!, w_vRHOLM(1:NRAD,1:LM_MAX,1:nat)
    REAL*8, intent(inout)     :: w_xwt1(0:21,1:nat), w_xwteh(0:3,1:nat), w_xwtel(0:3,1:nat), w_xwt1h(0:3,1:nat), w_xwt1l(0:3,1:nat)
    COMPLEX*16, intent(inout) :: sumfft(iff1,iff2,iff3), tloc(ift1,ift2,ift3)
    INTEGER, intent(in)       :: nomega, NRAD, LM_MAX, nat, natm, iff1, iff2, iff3, ift1, ift2, ift3, maxdim, ncix, ntcix, nip
    COMPLEX*16, intent(inout) :: DM(maxdim,maxdim,ncix)
    REAL*8, intent(inout)     :: Edimp(ntcix), Nds(nip,ntcix)
    COMPLEX*16, intent(inout) :: Gdloc(ntcix,nomega,nip)
    REAL*8, intent(inout)     :: fsph(3,natm), fnsp(3,natm), fsph2(3,natm)
    LOGICAL, intent(in)       :: Qforce
    !
    REAL*8 :: w_xwt, w_etot, w_ftot, w_ftot0
    COMPLEX*16, allocatable :: w_gloc(:)
    REAL*8, allocatable :: ww_RHOLM(:,:,:)!, ww_vRHOLM(:,:,:)
    REAL*8, allocatable :: ww_xwt1(:,:), ww_xwteh(:,:), ww_xwtel(:,:), ww_xwt1h(:,:), ww_xwt1l(:,:)
    COMPLEX*16, allocatable :: w_sumfft(:,:,:), w_tloc(:,:,:), w_Gdloc(:,:,:)
    COMPLEX*16, allocatable :: w_DM(:,:,:)
    REAL*8, allocatable :: w_Edimp(:), w_Nds(:,:)
    REAL*8, allocatable :: w_force(:,:)

    CALL MPI_REDUCE(xwt, w_xwt, 1,               MPI_REAL8,          MPI_SUM, master, MPI_COMM_WORLD, ierr)
    if (ierr.ne.0) print *, 'ERROR in MPI_REDUCE 1', ierr
    CALL MPI_REDUCE(etot,w_etot,1,               MPI_REAL8,          MPI_SUM, master, MPI_COMM_WORLD, ierr)
    if (ierr.ne.0) print *, 'ERROR in MPI_REDUCE 2', ierr
    CALL MPI_REDUCE(ftot,w_ftot,1,               MPI_REAL8,          MPI_SUM, master, MPI_COMM_WORLD, ierr)
    if (ierr.ne.0) print *, 'ERROR in MPI_REDUCE 2', ierr
    CALL MPI_REDUCE(ftot0,w_ftot0,1,             MPI_REAL8,          MPI_SUM, master, MPI_COMM_WORLD, ierr)
    if (ierr.ne.0) print *, 'ERROR in MPI_REDUCE 2c', ierr

    if (myrank.eq.master) then
       xwt = w_xwt
       etot=w_etot
       ftot=w_ftot
       ftot0=w_ftot0
    endif

    if (myrank.eq.master) ALLOCATE( w_gloc(nomega) )
    CALL MPI_REDUCE(gloc,w_gloc, nomega,          MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, ierr)
    if (ierr.ne.0) print *, 'ERROR in MPI_REDUCE 3', ierr
    if (myrank.eq.master) then
       gloc(:) = w_gloc(:)
       DEALLOCATE( w_gloc )
    endif

    if (myrank.eq.master) ALLOCATE( ww_RHOLM(1:NRAD,1:LM_MAX,1:nat) )
    CALL MPI_REDUCE(w_RHOLM, ww_RHOLM, NRAD*LM_MAX*nat, MPI_REAL8,          MPI_SUM, master, MPI_COMM_WORLD, ierr)
    if (ierr.ne.0) print *, 'ERROR in MPI_REDUCE 4', ierr
    if (myrank.eq.master) then
       w_RHOLM(:,:,:) = ww_RHOLM(:,:,:)
       DEALLOCATE( ww_RHOLM )
    endif

    !if (myrank.eq.master) ALLOCATE( ww_vRHOLM(1:NRAD,1:LM_MAX,1:nat) )
    !CALL MPI_REDUCE(w_vRHOLM, ww_vRHOLM, NRAD*LM_MAX*nat, MPI_REAL8,          MPI_SUM, master, MPI_COMM_WORLD, ierr)
    !if (ierr.ne.0) print *, 'ERROR in MPI_REDUCE 5', ierr
    !if (myrank.eq.master) then
    !   w_vRHOLM(:,:,:) = ww_vRHOLM(:,:,:)
    !   DEALLOCATE( ww_vRHOLM )
    !endif

    if (myrank.eq.master) ALLOCATE( ww_xwt1(0:21,1:nat) )
    CALL MPI_REDUCE(w_xwt1, ww_xwt1, 22*nat, MPI_REAL8,          MPI_SUM, master, MPI_COMM_WORLD, ierr)
    if (ierr.ne.0) print *, 'ERROR in MPI_REDUCE 6', ierr
    if (myrank.eq.master) then
       w_xwt1(:,:) = ww_xwt1(:,:)
       DEALLOCATE( ww_xwt1 )
    endif

    if (myrank.eq.master) ALLOCATE( ww_xwteh(0:3,1:nat) )
    CALL MPI_REDUCE(w_xwteh, ww_xwteh, 4*nat, MPI_REAL8,          MPI_SUM, master, MPI_COMM_WORLD, ierr)
    if (ierr.ne.0) print *, 'ERROR in MPI_REDUCE 7', ierr
    if (myrank.eq.master) then
       w_xwteh(:,:) = ww_xwteh(:,:)
       DEALLOCATE( ww_xwteh )
    endif

    if (myrank.eq.master) ALLOCATE( ww_xwtel(0:3,1:nat) )
    CALL MPI_REDUCE(w_xwtel, ww_xwtel, 4*nat, MPI_REAL8,          MPI_SUM, master, MPI_COMM_WORLD, ierr)
    if (ierr.ne.0) print *, 'ERROR in MPI_REDUCE 8', ierr
    if (myrank.eq.master) then
       w_xwtel(:,:) = ww_xwtel(:,:)
       DEALLOCATE( ww_xwtel )
    endif

    if (myrank.eq.master) ALLOCATE( ww_xwt1h(0:3,1:nat) )
    CALL MPI_REDUCE(w_xwt1h, ww_xwt1h, 4*nat, MPI_REAL8,          MPI_SUM, master, MPI_COMM_WORLD, ierr)
    if (ierr.ne.0) print *, 'ERROR in MPI_REDUCE 9', ierr
    if (myrank.eq.master) then
       w_xwt1h(:,:) = ww_xwt1h(:,:)
       DEALLOCATE( ww_xwt1h )
    endif

    if (myrank.eq.master) ALLOCATE( ww_xwt1l(0:3,1:nat) )
    CALL MPI_REDUCE(w_xwt1l, ww_xwt1l, 4*nat, MPI_REAL8,          MPI_SUM, master, MPI_COMM_WORLD, ierr)
    if (ierr.ne.0) print *, 'ERROR in MPI_REDUCE 10', ierr
    if (myrank.eq.master) then
       w_xwt1l(:,:) = ww_xwt1l(:,:)
       DEALLOCATE( ww_xwt1l )
    endif

    if (myrank.eq.master) ALLOCATE( w_sumfft(iff1,iff2,iff3) )
    CALL MPI_REDUCE(sumfft, w_sumfft,   iff1*iff2*iff3, MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, ierr)
    if (ierr.ne.0) print *, 'ERROR in MPI_REDUCE 11', ierr
    if (myrank.eq.master) then
       sumfft(:,:,:) = w_sumfft(:,:,:)
       DEALLOCATE( w_sumfft )
    endif

    if (Qforce) then
        if (myrank.eq.master) ALLOCATE( w_tloc(iff1,iff2,iff3) )
        CALL MPI_REDUCE(tloc, w_tloc,   iff1*iff2*iff3, MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, ierr)
        if (ierr.ne.0) print *, 'ERROR in MPI_REDUCE 12', ierr
        if (myrank.eq.master) then
           tloc(:,:,:) = w_tloc(:,:,:)
           DEALLOCATE( w_tloc )
        endif
    endif

    if (myrank.eq.master) ALLOCATE( w_DM(maxdim,maxdim,ncix) )
    CALL MPI_REDUCE(DM, w_DM,  maxdim*maxdim*ncix, MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, ierr)
    if (ierr.ne.0) print *, 'ERROR in MPI_REDUCE 13', ierr
    if (myrank.eq.master) then
       DM(:,:,:) = w_DM(:,:,:)
       DEALLOCATE( w_DM )
    endif

    if (myrank.eq.master) ALLOCATE( w_Edimp(ntcix) )
    CALL MPI_REDUCE(Edimp, w_Edimp,  ntcix, MPI_REAL8, MPI_SUM, master, MPI_COMM_WORLD, ierr)
    if (ierr.ne.0) print *, 'ERROR in MPI_REDUCE 14', ierr
    if (myrank.eq.master) then
       Edimp(:) = w_Edimp(:)
       DEALLOCATE( w_Edimp )
    endif


    if (myrank.eq.master) ALLOCATE( w_Gdloc(ntcix,nomega,nip) )
    CALL MPI_REDUCE(Gdloc, w_Gdloc,  ntcix*nomega*nip, MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, ierr)
    if (ierr.ne.0) print *, 'ERROR in MPI_REDUCE 15', ierr
    if (myrank.eq.master) then
       Gdloc(:,:,:) = w_Gdloc(:,:,:)
       DEALLOCATE( w_Gdloc )
    endif


    if (myrank.eq.master) ALLOCATE( w_Nds(nip,ntcix) )
    CALL MPI_REDUCE(Nds, w_Nds,  nip*ntcix, MPI_REAL8, MPI_SUM, master, MPI_COMM_WORLD, ierr)
    if (ierr.ne.0) print *, 'ERROR in MPI_REDUCE 16', ierr
    if (myrank.eq.master) then
       Nds(:,:) = w_Nds(:,:)
       DEALLOCATE( w_Nds )
    endif

    if (Qforce) then
       if (myrank.eq.master) ALLOCATE( w_force(3,natm) )
       CALL MPI_REDUCE(fsph, w_force,  3*natm, MPI_REAL8, MPI_SUM, master, MPI_COMM_WORLD, ierr)
       if (ierr.ne.0) print *, 'ERROR in MPI_REDUCE 17', ierr
       if (myrank.eq.master) fsph(:,:) = w_force(:,:)
       
       CALL MPI_REDUCE(fnsp, w_force,  3*natm, MPI_REAL8, MPI_SUM, master, MPI_COMM_WORLD, ierr)
       if (ierr.ne.0) print *, 'ERROR in MPI_REDUCE 17', ierr
       if (myrank.eq.master) fnsp(:,:) = w_force(:,:)
       
       CALL MPI_REDUCE(fsph2, w_force,  3*natm, MPI_REAL8, MPI_SUM, master, MPI_COMM_WORLD, ierr)
       if (ierr.ne.0) print *, 'ERROR in MPI_REDUCE 17', ierr
       if (myrank.eq.master) fsph2(:,:) = w_force(:,:)
       
       if (myrank.eq.master) DEALLOCATE( w_force )
    endif
  END SUBROUTINE Reduce_MPI

  


  SUBROUTINE Bcast_Size(size)
    IMPLICIT NONE
    INTEGER, intent(inout) :: size
    CALL MPI_BCAST(size, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
    if (ierr.ne.0) print *, 'ERROR in MPI_BCAST 3', ierr
  END SUBROUTINE Bcast_Size
  
  SUBROUTINE Bcast_Projector(P_rfk, P_rfi, dri, max_nrad, n_al_ucase, Nri)
    IMPLICIT NONE
    REAL*8, intent(inout)  :: P_rfk(max_nrad,2,n_al_ucase)
    REAL*8, intent(inout)  :: P_rfi(Nri,n_al_ucase)
    REAL*8, intent(inout)  :: dri(n_al_ucase)
    INTEGER, intent(inout) ::  max_nrad, n_al_ucase, Nri
    
    CALL MPI_BCAST(P_rfk, max_nrad*2*n_al_ucase, MPI_REAL8, master, MPI_COMM_WORLD, ierr)
    if (ierr.ne.0) print *, 'ERROR in MPI_BCAST 4', ierr
    CALL MPI_BCAST(P_rfi, Nri*n_al_ucase, MPI_REAL8, master, MPI_COMM_WORLD, ierr)
    if (ierr.ne.0) print *, 'ERROR in MPI_BCAST 5', ierr
    CALL MPI_BCAST(dri, n_al_ucase, MPI_REAL8, master, MPI_COMM_WORLD, ierr)
    if (ierr.ne.0) print *, 'ERROR in MPI_BCAST 6', ierr
  END SUBROUTINE Bcast_Projector
  
  subroutine Scatter_Vector_data()
    IMPLICIT NONE
    INTEGER       :: idat(2)
    CHARACTER*180 :: sdat(4), xend
    INTEGER,       allocatable :: idat_all(:,:)
    CHARACTER*180, allocatable :: sdat_all(:,:)
    INTEGER       :: ierr, irank, inkp, ik_start, i
    CHARACTER*180 :: vec_up, vec_dn, ene_dn, ene_up
    CHARACTER*67  :: ERRMSG
    !
    allocate( idat_all(2,nprocs), sdat_all(4,nprocs) )
    if (myrank.eq.master) then
       DO i=1,nprocs
          READ (999,*,END=970,ERR=970) irank,inkp,ik_start,vec_up,vec_dn,ene_dn,ene_up  ! can jump to 20
          idat_all(1,irank+1) = inkp       ! number of k-points
          idat_all(2,irank+1) = ik_start   ! what is the first k-point in this vector file
          sdat_all(1,irank+1) = vec_up     ! filenames....
          sdat_all(2,irank+1) = vec_dn
          sdat_all(3,irank+1) = ene_dn
          sdat_all(4,irank+1) = ene_up
       ENDDO
       !print *, 'idat_all, sdat_all='
       !do i=1,nprocs
       !   print *, i, idat_all(1:2,i)
       !   print *, i, TRIM(sdat_all(1,i)), TRIM(sdat_all(2,i)), TRIM(sdat_all(3,i)), TRIM(sdat_all(4,i))
       !end do
    end if

    call MPI_Scatter(idat_all, 2,  MPI_INTEGER, idat, 2, MPI_INT,  master, MPI_COMM_WORLD, ierr)
    if (ierr.NE.0) then
       WRITE(6,*) 'ERROR in MPI_Scatter 1'
       call stop_MPI
       STOP 'MPI ERROR 1'
    endif
    call MPI_Scatter(sdat_all, 4*180, MPI_CHAR, sdat, 4*180, MPI_CHAR, master, MPI_COMM_WORLD, ierr)
    if (ierr.NE.0) then
       WRITE(6,*) 'ERROR in MPI_Scatter 2'
       call stop_MPI
       STOP 'MPI ERROR 2'
    endif
    deallocate( idat_all, sdat_all )
    
    nvector=1
    vectors(nvector,1) = 1         ! the successive number of vector file
    vectors(nvector,2) = idat(1)   ! number of k-points
    vectors(nvector,3) = idat(2)   ! what is the first k-point in this vector file
    ! THIS IS NEW PART FROM SUMMER 2016
    ! Here we decided to take filenames from the dmft2.def-file, so that user
    ! has flexibility to change the order of files from the input.
    call get_xend(xend,sdat(1))  ! This will extract the suffix from _processes_ files, i.e., xend=='_myrank'
    ! Now we add xend to the filenames from 
    call add_xend(sdat(1), vecfn(1), xend) ! case.vectorso
    call add_xend(sdat(2), vecfn(2), xend) ! case.vectorsodn
    call add_xend(sdat(3), vecfn(3), xend) ! case.energysodn
    call add_xend(sdat(4), vecfn(4), xend) ! case.energyso
    ! THIS IS NEW PART FROM SUMMER 2016
    !
    fvectors(nvector,1) = sdat(1)  ! filenames....
    fvectors(nvector,2) = sdat(2)
    fvectors(nvector,3) = sdat(3)
    fvectors(nvector,4) = sdat(4)
    !print*, myrank, vectors(1,1), vectors(1,2), vectors(1,3)
    !print*, myrank, TRIM(fvectors(1,1)), TRIM(fvectors(1,2)), TRIM(fvectors(1,3)), TRIM(fvectors(1,4))
    return
970 CONTINUE
    WRITE (ERRMSG,'(A,A)')  'read error _processes_'
    CALL OUTERR('DMFT',ERRMSG)
    STOP 'DMFT - Error'
  end subroutine Scatter_Vector_data
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
    Qprint = .True.
  END SUBROUTINE start_MPI

  SUBROUTINE stop_MPI
  ENDSUBROUTINE stop_MPI

  SUBROUTINE FilenameMPI(infout)
    CHARACTER(LEN=*) :: infout
  ENDSUBROUTINE FilenameMPI

  SUBROUTINE FilenameMPI2(infout)
    CHARACTER(LEN=*) :: infout
  ENDSUBROUTINE FilenameMPI2
  
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

  SUBROUTINE FindMaxK_MPI(kmax)
    INTEGER, intent(inout) :: kmax(3)
  END SUBROUTINE FindMaxK_MPI


  SUBROUTINE Reduce_MPI2(numkpt, nkpt, sumw0, DM_EMIN, DM_EMAX)
    IMPLICIT NONE
    INTEGER, intent(out)  :: numkpt
    INTEGER, intent(in)   :: nkpt
    REAL*8, intent(inout) :: sumw0
    REAL*8, intent(inout) :: DM_EMIN, DM_EMAX
    numkpt = nkpt
  END SUBROUTINE Reduce_MPI2

  SUBROUTINE Reduce_maxbands(max_nbands)
    IMPLICIT NONE
    INTEGER, intent(inout) :: max_nbands
    !
  END SUBROUTINE Reduce_maxbands

  SUBROUTINE Reduce_maxnume(max_nume, nume)
    IMPLICIT NONE
    INTEGER, intent(out) :: max_nume
    INTEGER, intent(in)  :: nume
    !
    max_nume = nume
  END SUBROUTINE Reduce_maxnume
  
  SUBROUTINE Reduce_MPI0(Olapm, maxdim, ncix, max_nbands)
    IMPLICIT NONE
    COMPLEX*16, intent(inout) :: Olapm(maxdim,maxdim,ncix)
    INTEGER, intent(in) :: maxdim, ncix, max_nbands
  END SUBROUTINE Reduce_MPI0

  SUBROUTINE Reduce_MPI1(pr_procr, nomega, max_nume)
    use muzero, ONLY: nkp, wgh, nemm, zEk, Ek, max_nbands
    IMPLICIT NONE
    INTEGER, intent(in)   :: pr_procr, nomega, max_nume
    nkp = sum(pr_procs)
  END SUBROUTINE Reduce_MPI1

  SUBROUTINE Reduce_MPI1b(nomega, Ebmax, Ebmin, weib, gtot, numex)
    IMPLICIT NONE
    INTEGER, intent(in)   :: nomega, numex
    REAL*8, intent(inout) :: Ebmin(numex), Ebmax(numex), weib(numex)
    COMPLEX*16, intent(inout) :: gtot(nomega)
    !
  END SUBROUTINE Reduce_MPI1b


  SUBROUTINE Reduce_MPI(xwt, etot, ftot, ftot0, gloc, w_RHOLM, DM, Gdloc, fsph, fnsp, fsph2, Edimp, Nds, w_xwt1, w_xwteh, w_xwtel, w_xwt1h, w_xwt1l, sumfft, tloc, nomega, NRAD, LM_MAX, nat, natm, iff1, iff2, iff3, ift1, ift2, ift3, maxdim, ncix, ntcix, nip, Qforce)
    IMPLICIT NONE
    REAL*8, intent(inout)     :: xwt, etot, ftot, ftot0
    COMPLEX*16, intent(inout) :: gloc(nomega)
    REAL*8, intent(inout)     :: w_RHOLM(1:NRAD,1:LM_MAX,1:nat)!, w_vRHOLM(1:NRAD,1:LM_MAX,1:nat)
    REAL*8, intent(inout)     :: w_xwt1(0:21,1:nat), w_xwteh(0:3,1:nat), w_xwtel(0:3,1:nat), w_xwt1h(0:3,1:nat), w_xwt1l(0:3,1:nat)
    COMPLEX*16, intent(inout) :: sumfft(iff1,iff2,iff3), tloc(ift1,ift2,ift3)
    INTEGER, intent(in)       :: nomega, NRAD, LM_MAX, nat, natm, iff1, iff2, iff3, ift1, ift2, ift3, maxdim, ncix, ntcix, nip
    COMPLEX*16, intent(inout) :: DM(maxdim,maxdim,ncix)
    REAL*8, intent(inout)     :: Edimp(ntcix), Nds(nip,ntcix)
    COMPLEX*16, intent(inout) :: Gdloc(ntcix,nomega,nip)
    REAL*8, intent(inout)     :: fsph(3,natm), fnsp(3,natm), fsph2(3,natm)
    LOGICAL, intent(in)       :: Qforce
  END SUBROUTINE Reduce_MPI

  SUBROUTINE Bcast_MPI(EF)
    IMPLICIT NONE
    REAL*8, intent(inout) :: EF
  END SUBROUTINE Bcast_MPI

  SUBROUTINE Gather_procs(pr_procr, pr_procs, nprocs)
    IMPLICIT NONE
    INTEGER, intent(in) :: pr_procr, nprocs
    INTEGER, intent(out):: pr_procs(nprocs)
    pr_procs(1) = pr_procr
  END SUBROUTINE Gather_procs

  SUBROUTINE Bcast_Size(size)
    IMPLICIT NONE
    INTEGER, intent(inout) :: size
  END SUBROUTINE Bcast_Size
  
  SUBROUTINE Bcast_Projector(P_rfk, P_rfi, dri, max_nrad, n_al_ucase, Nri)
    IMPLICIT NONE
    REAL*8, intent(inout)  :: P_rfk(max_nrad,2,n_al_ucase)
    REAL*8, intent(inout)  :: P_rfi(Nri,n_al_ucase)
    REAL*8, intent(inout)  :: dri(n_al_ucase)
    INTEGER, intent(inout) ::  max_nrad, n_al_ucase, Nri
  END SUBROUTINE Bcast_Projector

  subroutine Scatter_Vector_data()
  end subroutine Scatter_Vector_data

#endif
END MODULE com_mpi

SUBROUTINE FilenameMPI3(infout, ii)
  IMPLICIT NONE
  CHARACTER(LEN=*) :: infout
  INTEGER, intent(in) :: ii
  CHARACTER*10 :: strii
  write(strii,'(I10)') ii
  infout    = TRIM(infout)//trim(ADJUSTL(strii))
END SUBROUTINE FilenameMPI3

subroutine get_xend(xend,filename)
  IMPLICIT NONE
  CHARACTER*180, intent(in) :: filename
  CHARACTER*180, intent(out) :: xend
  ! locals                                                                                                                                                                  
  INTEGER :: i, ln, ln0
  ln0=len_trim(filename)
  ln=ln0
  DO i=ln0,1,-1
     IF (filename(i:i) .EQ. '_') THEN
        ln = i-1
        exit
     ENDIF
  ENDDO
  xend = filename(ln+1:ln0)
end subroutine get_xend

subroutine add_xend(fileout,filein,xend)
  IMPLICIT NONE
  CHARACTER*180, intent(out) :: fileout
  CHARACTER*200, intent(in) :: filein
  CHARACTER*180, intent(in) :: xend
  ! locals                                                                                                                                                                  
  INTEGER :: ln
  ln=len_trim(filein)
  if (filein(ln-1:ln).eq.'_x') then
     fileout(:) = filein(1:ln-2)//xend
  else
     fileout(:) = filein(1:ln)//xend
  endif
end subroutine add_xend
