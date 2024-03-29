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
    if (ierr.ne.0) WRITE(*,*) 'ERROR in MPI_INIT ', ierr
    CALL MPI_COMM_SIZE( MPI_COMM_WORLD, nprocs, ierr)
    if (ierr.ne.0) WRITE(*,*) 'ERROR in MPI_COMM_SIZE ', ierr
    CALL MPI_COMM_RANK( MPI_COMM_WORLD, myrank, ierr)
    if (ierr.ne.0) WRITE(*,*) 'ERROR in MPI_COMM_RANK ', ierr
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

  SUBROUTINE FindMax_MPI(max_bands, maxbands)
    INTEGER, intent(out) :: max_bands
    INTEGER, intent(in)  :: maxbands
    ! locals
    CALL MPI_ALLREDUCE(maxbands, max_bands, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
  END SUBROUTINE FindMax_MPI
  
  SUBROUTINE Reduce_MPI(gloc, gtot, gmloc, Olapm, Eimpm, norbitals, nomega, maxdim, ncix)                                                                     
    IMPLICIT NONE                                                                                                                                             
    COMPLEX*16, allocatable, intent(inout) :: gloc(:,:)  ! gloc(norbitals,nomega)                                                                                                       
    COMPLEX*16, intent(inout) :: gtot(nomega)                                                                                                                 
    COMPLEX*16, intent(inout) :: gmloc(maxdim, maxdim, ncix, nomega)                                                                                          
    COMPLEX*16, intent(inout) :: Olapm(maxdim,maxdim,ncix), Eimpm(maxdim,maxdim,ncix)                                                                         
    INTEGER, intent(in)       :: norbitals, nomega, maxdim, ncix                                                                                              
    ! locals                                                                                                                                                  
    COMPLEX*16, allocatable :: wgtot(:)                                                                                                                       
    COMPLEX*16, allocatable :: wgloc(:,:)                                                                                                                     
    COMPLEX*16, allocatable :: wOlapm(:,:,:), wEimpm(:,:,:)                                                                                                   
    COMPLEX*16, allocatable :: wgmloc(:,:,:,:)                                                                                                                
                                                                                                                                                              
    if (myrank.eq.master) then                                                                                                                                
       ALLOCATE( wgtot(nomega), wOlapm(maxdim,maxdim,ncix), wEimpm(maxdim,maxdim,ncix) )
       if (allocated(gloc)) ALLOCATE( wgloc(norbitals,nomega) )
    endif
    CALL MPI_REDUCE(gtot,  wgtot,  nomega,           MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, ierr)                                               
    if (ierr.ne.0) print *, 'ERROR in MPI_REDUCE 1', ierr

    if (allocated(gloc)) then
       CALL MPI_REDUCE(gloc,  wgloc,  norbitals*nomega, MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, ierr)                                               
       if (ierr.ne.0) print *, 'ERROR in MPI_REDUCE 2', ierr
    endif
    
    CALL MPI_REDUCE(Olapm, wOlapm, maxdim*maxdim*ncix, MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, ierr)                                             
    if (ierr.ne.0) print *, 'ERROR in MPI_REDUCE 3', ierr
    CALL MPI_REDUCE(Eimpm, wEimpm, maxdim*maxdim*ncix, MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, ierr)                                             
    if (myrank.eq.master) then                                                                                                                                
       gtot = wgtot                                                                                                                                           
       Olapm = wOlapm                                                                                                                                         
       Eimpm = wEimpm                                                                                                                                         
       DEALLOCATE( wgtot, wOlapm, wEimpm )
       if (allocated(gloc)) then
          gloc = wgloc
          DEALLOCATE( wgloc )
       endif
    endif
                                                                                                                                                              
    if (myrank.eq.master) then                                                                                                                                
       ALLOCATE( wgmloc(maxdim, maxdim, ncix, nomega) )                                                                                                       
    endif
    CALL MPI_REDUCE(gmloc, wgmloc, maxdim*maxdim*ncix*nomega, MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, ierr)                                      
    if (ierr.ne.0) print *, 'ERROR in MPI_REDUCE 4', ierr
    if (myrank.eq.master) then                                                                                                                                
       gmloc = wgmloc                                                                                                                                         
       DEALLOCATE( wgmloc )                                                                                                                                   
    endif
  END SUBROUTINE Reduce_MPI

  SUBROUTINE Reduce_dm_MPI(g_inf, g_ferm, maxdim, ncix, nomega)
    IMPLICIT NONE
    complex*16, intent(inout) :: g_inf(maxdim,maxdim,ncix,nomega), g_ferm(maxdim,maxdim,ncix)
    INTEGER, intent(in) :: maxdim, ncix, nomega    
    ! locals
    complex*16, allocatable :: w_g_inf(:,:,:,:), w_g_ferm(:,:,:)
    
    if (myrank.eq.master) then                                                                                                                                
       ALLOCATE( w_g_inf(maxdim,maxdim,ncix,nomega), w_g_ferm(maxdim,maxdim,ncix) )
    endif
    CALL MPI_REDUCE(g_inf ,  w_g_inf,  maxdim*maxdim*ncix*nomega, MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, ierr)
    CALL MPI_REDUCE(g_ferm,  w_g_ferm, maxdim*maxdim*ncix,        MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, ierr)
    if (ierr.ne.0) print *, 'ERROR in MPI_REDUCE 10', ierr
    if (myrank.eq.master) then
       g_inf(:,:,:,:) = w_g_inf(:,:,:,:)
       g_ferm(:,:,:) = w_g_ferm(:,:,:)
       DEALLOCATE( w_g_inf, w_g_ferm )
    endif
  END SUBROUTINE Reduce_dm_MPI
        

  SUBROUTINE AllReduce_MPI(Olapm, maxdim, ncix)                                                                     
    IMPLICIT NONE                                                                                                                                             
    COMPLEX*16, intent(inout) :: Olapm(maxdim,maxdim,ncix)
    INTEGER, intent(in)       :: maxdim, ncix                                                                                              
    ! locals                                                                                                                                                  
    COMPLEX*16, allocatable :: wOlapm(:,:,:)
    
    ALLOCATE( wOlapm(maxdim,maxdim,ncix) )
    
    CALL MPI_ALLREDUCE(Olapm, wOlapm, maxdim*maxdim*ncix, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ierr)
    
    if (ierr.ne.0) print *, 'ERROR in MPI_REDUCE 3', ierr
    Olapm = wOlapm
    DEALLOCATE( wOlapm )
    
  END SUBROUTINE AllReduce_MPI


  SUBROUTINE Gather_MPI(tEk, tnbands, tnemin, tn_ik, Ekp, nbandsk, nemink, n_ik, pr_procr, pr_procs, nprocs, numkpt, max_bands, nomega, nsymop)
    IMPLICIT NONE
    INTEGER, intent(in)     :: pr_procr, nprocs, numkpt, max_bands, nomega, nsymop, pr_procs(nprocs)
    !COMPLEX*16, intent(out) :: tEk(nomega,nsymop,max_bands,numkpt)
    COMPLEX*16, intent(out) :: tEk(*)
    INTEGER, intent(out)    :: tnbands(*), tnemin(*), tn_ik(*)  ! tnbands(numkpt), tnemin(numkpt), tn_ik(numkpt)
    COMPLEX*16, intent(in)  :: Ekp(nomega,nsymop,max_bands,pr_procr)
    INTEGER, intent(in)     :: nbandsk(pr_procr), nemink(pr_procr), n_ik(pr_procr)
    ! locals
    INTEGER :: i
    INTEGER, allocatable :: offs(:), pr_procs_local(:)
    
    ALLOCATE( offs(nprocs), pr_procs_local(nprocs) )
    offs(1)=0
    do i=2,nprocs
       offs(i)= offs(i-1)+pr_procs(i-1)
    enddo

    CALL MPI_GATHERV(nbandsk, pr_procr, MPI_INTEGER, tnbands, pr_procs, offs, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
    if (ierr.ne.0) print *, 'ERROR in MPI_GATHER 1', ierr
    
    CALL MPI_GATHERV(nemink,  pr_procr, MPI_INTEGER, tnemin,  pr_procs, offs, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
    if (ierr.ne.0) print *, 'ERROR in MPI_GATHER 2', ierr

    CALL MPI_GATHERV(n_ik,  pr_procr, MPI_INTEGER, tn_ik,  pr_procs, offs, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
    if (ierr.ne.0) print *, 'ERROR in MPI_GATHER 2', ierr

    do i=1,nprocs
       pr_procs_local(i) = pr_procs(i)*nomega*nsymop*max_bands
       offs(i) = offs(i)*nomega*nsymop*max_bands
    enddo
    
    !print *, 'sending from: ',  myrank, nomega*nsymop*max_bands*pr_procr, 'totalsize=', nomega*nsymop*max_bands*numkpt
    !if (myrank.EQ.master) then
    !    print *, 'pr_procs_local:', pr_procs_local(:)
    !	print *, 'offs:', offs
    !endif
	    
    CALL MPI_GATHERV(Ekp,  nomega*nsymop*max_bands*pr_procr, MPI_DOUBLE_COMPLEX, tEk,  pr_procs_local, offs, MPI_DOUBLE_COMPLEX, master, MPI_COMM_WORLD, ierr)
    if (ierr.ne.0) print *, 'ERROR in MPI_GATHER 3', ierr

    DEALLOCATE( offs, pr_procs_local )
    
  END SUBROUTINE Gather_MPI

  SUBROUTINE Gather_procs(pr_procr, pr_procs, nprocs)
    IMPLICIT NONE
    INTEGER, intent(in) :: pr_procr, nprocs
    INTEGER, intent(out):: pr_procs(nprocs)
    CALL MPI_GATHER(pr_procr, 1, MPI_INTEGER, pr_procs, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
    if (ierr.ne.0) print *, 'ERROR in MPI_GATHER pr_procs', ierr
  END SUBROUTINE Gather_procs
  

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
    ! Here we decided to take filenames from the dmft1.def-file, so that user
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
    print*, myrank, TRIM(fvectors(1,1)), ' ', TRIM(fvectors(1,2)), ' ', TRIM(fvectors(1,3)), ' ', TRIM(fvectors(1,4))
    return
970 CONTINUE
    WRITE(99,*) 'ERROR dmft1 - read error _processes_'
    call stop_MPI
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

  SUBROUTINE FindMax_MPI(max_bands, maxbands)
    INTEGER, intent(out) :: max_bands
    INTEGER, intent(in)  :: maxbands
    max_bands = maxbands
  END SUBROUTINE FindMax_MPI

  SUBROUTINE Reduce_MPI(gloc, gtot, gmloc, Olapm, Eimpm, norbitals, nomega, maxdim, ncix)
    IMPLICIT NONE
    COMPLEX*16, intent(inout) :: gloc(norbitals,nomega)
    COMPLEX*16, intent(inout) :: gtot(nomega)
    COMPLEX*16, intent(inout) :: gmloc(maxdim, maxdim, ncix, nomega)
    COMPLEX*16, intent(inout) :: Olapm(maxdim,maxdim,ncix), Eimpm(maxdim,maxdim,ncix)
    INTEGER, intent(in)       :: norbitals, nomega, maxdim, ncix
    ! locals
  END SUBROUTINE Reduce_MPI

  SUBROUTINE Reduce_dm_MPI(g_inf, g_ferm, maxdim2, ncix, nomega)
    IMPLICIT NONE
    complex*16, intent(inout) :: g_inf(maxdim2,maxdim2,ncix,nomega), g_ferm(maxdim2,maxdim2,ncix)
    INTEGER, intent(in) :: maxdim2, ncix, nomega
  END SUBROUTINE Reduce_dm_MPI
  
  SUBROUTINE AllReduce_MPI(Olapm, maxdim, ncix)                                                                     
    IMPLICIT NONE                                                                                                                                             
    COMPLEX*16, intent(inout) :: Olapm(maxdim,maxdim,ncix)
    INTEGER, intent(in)       :: maxdim, ncix                                                                                              
  END SUBROUTINE AllReduce_MPI


  SUBROUTINE Gather_MPI(tEk, tnbands, tnemin, tn_ik, Ekp, nbandsk, nemink, n_ik, pr_procr, pr_procs, nprocs, numkpt, max_bands, nomega, nsymop)
    IMPLICIT NONE
    INTEGER, intent(in)     :: pr_procr, nprocs, numkpt, max_bands, nomega, nsymop, pr_procs(nprocs)
    COMPLEX*16, intent(out) :: tEk(nomega,nsymop,max_bands,numkpt)
    INTEGER, intent(out)    :: tnbands(numkpt), tnemin(numkpt), tn_ik(numkpt)
    COMPLEX*16, intent(in)  :: Ekp(nomega,nsymop,max_bands,pr_procr)
    INTEGER, intent(in)     :: nbandsk(pr_procr), nemink(pr_procr), n_ik(pr_procr)
    tnbands(:numkpt) = nbandsk(:)
    tnemin(:numkpt) = nemink(:)
    tEk(:,:,:,:) = Ekp(:,:,:,:)
  END SUBROUTINE Gather_MPI

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


!SUBROUTINE FNAME_PROCESS(FNAME, myrank)
!  IMPLICIT NONE
!  CHARACTER*200,intent(out) :: FNAME
!  INTEGER,intent(in) :: myrank
!
!  if (myrank.LT.10) then
!     WRITE(FNAME,'(A,I1)') '.processes_',myrank
!  else if (myrank.LT.100) then
!     WRITE(FNAME,'(A,I2)') '.processes_',myrank
!  else if (myrank.LT.1000) then
!     WRITE(FNAME,'(A,I3)') '.processes_',myrank
!  else if (myrank.LT.10000) then
!     WRITE(FNAME,'(A,I4)') '.processes_',myrank
!  else if (myrank.LT.10000) then
!     WRITE(FNAME,'(A,I5)') '.processes_',myrank
!  else if (myrank.LT.10000) then
!     WRITE(FNAME,'(A,I6)') '.processes_',myrank
!  endif
!END SUBROUTINE FNAME_PROCESS
!
!SUBROUTINE vector_filename(FNAME, name, ii)
!  IMPLICIT NONE
!  CHARACTER*200,intent(out) :: FNAME
!  CHARACTER*200,intent(in) :: name
!  INTEGER,intent(in) :: ii
!  INTEGER:: lng
!  lng = Len_Trim(name)
!  if (ii.LT.10) then
!     WRITE(FNAME,'(A,A1,I1)') name, '_', ii
!  else if (ii.LT.100) then
!     WRITE(FNAME,'(A,A1,I2)') name, '_', ii
!  else if (ii.LT.1000) then
!     WRITE(FNAME,'(A,A1,I3)') name, '_', ii
!  else if (ii.LT.10000) then
!     WRITE(FNAME,'(A,A1,I4)') name, '_', ii
!  else if (ii.LT.10000) then
!     WRITE(FNAME,'(A,A1,I5)') name, '_', ii
!  else if (ii.LT.10000) then
!     WRITE(FNAME,'(A,A1,I6)') name, '_', ii
!  endif
!END SUBROUTINE vector_filename

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
