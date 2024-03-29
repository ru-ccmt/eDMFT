!**************************************************!
!******* § Copyright by Kristjan Haule, 2016 ******!
!**************************************************!
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
  INTEGER      :: pr_proc       ! maximum number of k-points, which must be done on any processor
  LOGICAL      :: fastFilesystem        ! Should each core print info?
  INTEGER, parameter :: tsize=16
  !
  CHARACTER*180     :: filename_vector(2), filename_energy(2), filename_vectorso(2), filename_energyso(3), filename_norm(2)
  INTEGER           :: ikps(2)
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

  subroutine mpi_bcast_V_vsp(Vr, NRAD,NAT,jspin)
    IMPLICIT NONE
    INTEGER, intent(in)   :: NRAD, NAT,jspin
    REAL*8, intent(inout) :: Vr(NRAD,NAT,jspin)
    !
    CALL MPI_BCAST(Vr, NRAD*NAT*jspin, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
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

  subroutine mpi_SendReceive12(ss,weight,bname, en,vt, nv_plus_nnrlo,neig,kt,isi,jspin,must_compute)
    IMPLICIT NONE
    real*8, intent(in)  :: ss(3), weight
    character*10, intent(in) :: bname
    integer, intent(in) :: nv_plus_nnrlo, neig, isi, jspin
    integer, intent(in) :: kt(3,nv_plus_nnrlo)
    logical, intent(in) :: must_compute
    !
    real*8, intent(in)  :: en(neig)
    complex*16, intent(in) :: vt(nv_plus_nnrlo,neig)
    ! locals
    LOGICAL :: was_computed
    INTEGER :: status(MPI_STATUS_SIZE)
    real*8  :: realdata(4)
    integer :: intdata(3)
    !
    real*8  :: tss(3), tweight
    integer :: ipr, tisi, tnv_plus_nnrlo, tneig, j
    character*10 :: tname
    integer, allocatable :: tkt(:,:)
    real*8, allocatable :: ten(:)
    complex*16, allocatable :: tvt(:,:)
    !
    do ipr=1,nprocs-1 ! Notice that all processors need to tell the master if they have something to send or not.
       ! Master can not know whether it should expect something or no.
       if (myrank.eq.master) then
          call MPI_Recv ( was_computed, 1, MPI_LOGICAL, ipr,      0, MPI_COMM_WORLD, status, ierr )
          if (ierr.ne.0) print *, 'ERROR in MPI_Recv 0', ierr
          if (was_computed) then  ! Now even master knows that it might be empty k-point, and can skip it.
             call MPI_Recv ( intdata, 3, MPI_INTEGER, ipr,      1, MPI_COMM_WORLD, status, ierr )
             if (ierr.ne.0) print *, 'ERROR in MPI_Recv 1', ierr
             tisi=intdata(1)
             tneig=intdata(2)
             tnv_plus_nnrlo=intdata(3)
             call MPI_Recv ( realdata, 4, MPI_DOUBLE_PRECISION, ipr, 2, MPI_COMM_WORLD, status, ierr )
             if (ierr.ne.0) print *, 'ERROR in MPI_Recv 2', ierr
             tss(1:3) = realdata(1:3)
             tweight = realdata(4)
             call MPI_Recv ( tname, 10, MPI_CHAR, ipr,    3, MPI_COMM_WORLD, status, ierr )
             if (ierr.ne.0) print *, 'ERROR in MPI_Recv 3', ierr
             allocate(tkt(3,tnv_plus_nnrlo))
             call MPI_Recv ( tkt, 3*tnv_plus_nnrlo, MPI_INTEGER, ipr,    4, MPI_COMM_WORLD, status, ierr )
             if (ierr.ne.0) print *, 'ERROR in MPI_Recv 4', ierr
             allocate( ten(tneig) )
             call MPI_Recv ( ten, tneig, MPI_DOUBLE_PRECISION, ipr, 5, MPI_COMM_WORLD, status, ierr )
             if (ierr.ne.0) print *, 'ERROR in MPI_Recv 5', ierr
             allocate(tvt(tnv_plus_nnrlo,tneig))
             call MPI_Recv ( tvt, 2*tnv_plus_nnrlo*tneig, MPI_DOUBLE_PRECISION, ipr, 6, MPI_COMM_WORLD, status, ierr )
             if (ierr.ne.0) print *, 'ERROR in MPI_Recv 6', ierr
             
             ! Finally printing this k-point on master
             call Print_Kpoint1(41+isi-1,51+isi-1, isi.eq.2, tkt, tss, tname, tweight, tnv_plus_nnrlo, tneig)
             do j=1,tneig
                call Print_Kpoint2(41+tisi-1, 51+tisi-1, (tisi.eq.2.and.j.eq.tneig), j, ten(j), tvt(:,j), tnv_plus_nnrlo, tneig)
             enddo
             
             deallocate( tkt )
             deallocate( ten, tvt )
          endif
       else if (myrank.eq.ipr) then
          !print *, 'I am sending', myrank, itmp(1:3)
          call MPI_Send ( must_compute, 1, MPI_LOGICAL, master, 0, MPI_COMM_WORLD, ierr )
          if (ierr.ne.0) print *, 'ERROR in MPI_Send 0', ierr
          if (must_compute) then
             realdata(1:3) = ss(1:3)
             realdata(4) = weight
             intdata(1)=isi
             intdata(2)=neig
             intdata(3)=nv_plus_nnrlo
             call MPI_Send ( intdata, 3, MPI_INTEGER, master, 1, MPI_COMM_WORLD, ierr )
             if (ierr.ne.0) print *, 'ERROR in MPI_Send 1', ierr
             call MPI_Send ( realdata, 4, MPI_DOUBLE_PRECISION, master, 2, MPI_COMM_WORLD, ierr )
             if (ierr.ne.0) print *, 'ERROR in MPI_Send 2', ierr
             call MPI_Send ( bname, 10, MPI_CHAR, master, 3, MPI_COMM_WORLD, ierr )
             if (ierr.ne.0) print *, 'ERROR in MPI_Send 3', ierr
             call MPI_Send (kt, 3*nv_plus_nnrlo, MPI_INTEGER, master, 4, MPI_COMM_WORLD, ierr )
             if (ierr.ne.0) print *, 'ERROR in MPI_Send 4', ierr
             call MPI_Send ( en, neig, MPI_DOUBLE_PRECISION, master, 5, MPI_COMM_WORLD, ierr )
             if (ierr.ne.0) print *, 'ERROR in MPI_Send 5', ierr
             call MPI_Send ( vt, 2*nv_plus_nnrlo*neig, MPI_DOUBLE_PRECISION, master, 6, MPI_COMM_WORLD, ierr )
             if (ierr.ne.0) print *, 'ERROR in MPI_Send 6', ierr
          endif
       endif
    enddo
  end subroutine mpi_SendReceive12

  
  subroutine mpi_SendReceive3(vnorm, neig, must_compute)
    IMPLICIT NONE
    real*8, intent(in)  :: vnorm(neig,2)
    integer, intent(in) :: neig
    logical, intent(in) :: must_compute
    ! locals
    INTEGER :: status(MPI_STATUS_SIZE)
    integer :: tneig, ipr
    real*8, allocatable :: tnorm(:,:)
    !
    tneig=neig
    if (.not.must_compute) tneig=0
    
    do ipr=1,nprocs-1 ! Notice that all processors need to tell the master if they have something to send or not.
       ! Master can not know whether it should expect something or no.
       if (myrank.eq.master) then
          call MPI_Recv ( tneig, 1, MPI_INTEGER, ipr,      0, MPI_COMM_WORLD, status, ierr )
          if (ierr.ne.0) print *, 'ERROR in MPI_Recv 0', ierr
          if (tneig.gt.0) then  ! Now even master knows that it might be empty k-point, and can skip it.
             allocate( tnorm(tneig,2) )
             call MPI_Recv ( tnorm, tneig*2, MPI_DOUBLE_PRECISION, ipr, 1, MPI_COMM_WORLD, status, ierr )
             if (ierr.ne.0) print *, 'ERROR in MPI_Recv 1', ierr
             call Print_Kpoint3( tnorm, tneig)
             deallocate( tnorm )
          endif
       else if (myrank.eq.ipr) then
          call MPI_Send ( tneig, 1, MPI_INTEGER, master, 0, MPI_COMM_WORLD, ierr )
          if (ierr.ne.0) print *, 'ERROR in MPI_Send 0', ierr
          if (tneig.gt.0) then
             call MPI_Send ( vnorm, neig*2, MPI_DOUBLE_PRECISION, master, 1, MPI_COMM_WORLD,     ierr )
             if (ierr.ne.0) print *, 'ERROR in MPI_Send 1', ierr
          endif
       endif
    enddo
  end subroutine mpi_SendReceive3
  
  subroutine Scatter_Vector_data(ierr0)
    IMPLICIT NONE
    INTEGER, intent(out) :: ierr0
    !
    INTEGER       :: idat(3), i, ios
    CHARACTER*180 :: sdat(4)
    INTEGER,       allocatable :: idat_all(:,:)
    CHARACTER*180, allocatable :: sdat_all(:,:)
    INTEGER       :: ierr, irank, inkp, ik_start
    CHARACTER*180 :: vec_up, vec_dn, ene_dn, ene_up
    CHARACTER*67  :: ERRMSG
    CHARACTER*180 :: xend, FNAME
    !
    FNAME = '_processes_'
    if (myrank.eq.master) open(999,FILE=FNAME,STATUS='OLD',ERR=77,IOSTAT=ios) ! opens _processes_, hence this parallel dmft_lapwso run, and not parallel w2k_lapwso run
    
    allocate( idat_all(3,nprocs), sdat_all(4,nprocs) )
    if (myrank.eq.master) then
       idat_all(:,:)=0
       DO i=1,nprocs
          READ (999,*,END=970,ERR=970) irank,inkp,ik_start,vec_up,vec_dn,ene_dn,ene_up  ! can jump to 20
          idat_all(1,irank+1) = inkp       ! number of k-points
          idat_all(2,irank+1) = ik_start   ! what is the first k-point in this vector file
          sdat_all(1,irank+1) = vec_up     ! filenames....
          sdat_all(2,irank+1) = vec_dn
          sdat_all(3,irank+1) = ene_dn
          sdat_all(4,irank+1) = ene_up
       ENDDO
       pr_proc = maxval(idat_all(1,:))
       do i=1,nprocs
          idat_all(3,i) = pr_proc
       enddo
       !print *, 'idat_all, sdat_all='
       !do i=1,nprocs
       !   print *, i, idat_all(1:2,i)
       !   print *, i, TRIM(sdat_all(1,i)), TRIM(sdat_all(2,i)), TRIM(sdat_all(3,i)), TRIM(sdat_all(4,i))
       !end do
    end if

    call MPI_Scatter(idat_all, 3,  MPI_INTEGER, idat, 3, MPI_INT,  master, MPI_COMM_WORLD, ierr)
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

    call get_xend(xend,sdat(1))  ! This will extract the suffix from _processes_ files, i.e., xend=='_myrank'
    ! Now we take the filenames from lapwso.def def and attach the correct suffix on each core
    call add_xend(xend, filename_vector(1))    !'case.vector_x'
    call add_xend(xend, filename_vector(2))    !'case.vectorup_x'
    call add_xend(xend, filename_energy(1))    !'case.energy_x'
    call add_xend(xend, filename_energy(2))    !'case.energyup_x'
    call add_xend(xend, filename_vectorso(1))  !'case.vectorsodn'
    call add_xend(xend, filename_vectorso(2))  !'case.vectorso'
    call add_xend(xend, filename_energyso(1))  !'case.energysodn'
    call add_xend(xend, filename_energyso(2))  !'case.energyso'
    call add_xend(xend, filename_energyso(3))  !'case.energydum'
    call add_xend(xend, filename_norm(1))      !'case.normsodn'
    call add_xend(xend, filename_norm(2))      !'case.normsoup'
    ikps(2) = idat(1)   ! number of k-points
    ikps(1) = idat(2)   ! what is the first k-point in this vector file
    pr_proc = idat(3)   ! maximum number of k-points which should be done on any core.
    !print*, myrank, TRIM(fvectors(1,1)), TRIM(fvectors(1,2)), TRIM(fvectors(1,3)), TRIM(fvectors(1,4))

    if (Qprint) print *, 'Vector files are read in parallel mode'
    if (myrank.eq.master) close(999)

    ierr0=0
    return 
77  CONTINUE
    ierr0=1
    return
970 CONTINUE
    WRITE (ERRMSG,'(A,A)')  'read error _processes_'
    CALL OUTERR('DMFT',ERRMSG)
    STOP 'DMFT - Error'
  end subroutine Scatter_Vector_data


  subroutine MPI_Bcast_dims(nmat,nume,nkp)
    IMPLICIT NONE
    INTEGER, intent(inout) :: nmat, nume, nkp
    INTEGER :: data(3)
    data(:)=0
    if (myrank.eq.master) then ! pack the data into array to use single broadcat
       data(1)=nmat
       data(2)=nume
       data(3)=nkp
    endif
    CALL MPI_BCAST(data, 3, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
    if (ierr.ne.0) print *, 'ERROR in MPI_BCAST mpi_bcast 10', ierr
    if (myrank.ne.master) then ! unpack the data
       nmat=data(1)
       nume=data(2)
       nkp=data(3)
    endif
  end subroutine MPI_Bcast_dims

  subroutine MPI_Bcast_dims_para(nmat,nume,nkp)
    IMPLICIT NONE
    INTEGER, intent(inout) :: nmat, nume, nkp
    INTEGER :: data(2), data_(2), nkp_
    
    nkp_=nkp;
    call MPI_Allreduce(nkp_, nkp, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD,ierr)

    ! This might not be necessary, as it is enough to have maximum nmat,nume for the k-points being done on this processor
    !data_(:)=0
    !if (myrank.eq.master) then ! pack the data into array to use single broadcat
    !   data_(1)=nmat
    !   data_(2)=nume
    !endif
    !call MPI_Allreduce(data_, data, 2, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD,ierr)
    !if (ierr.ne.0) print *, 'ERROR in MPI_bcast_dims_para', ierr
    !nmat=data(1)
    !nume=data(2)
  end subroutine MPI_Bcast_dims_para

  subroutine WriteProcesses()
    IMPLICIT NONE
    ! locals
    INTEGER       :: i
    CHARACTER*180 :: filename_vectordn, filename_energydn
    INTEGER       :: idat(3)
    CHARACTER*180 :: sdat(4)
    INTEGER,      allocatable :: idat_all(:,:)
    CHARACTER*180,allocatable :: sdat_all(:,:)

    idat(1)=myrank
    idat(2)=ikps(2) ! number of k-points
    idat(3)=ikps(1) ! what is the first k-point in this vector file
    sdat(1)=filename_vectorso(2)
    sdat(2)=filename_vectorso(1)
    sdat(3)=filename_energyso(1)
    sdat(4)=filename_energyso(2)
    if (myrank.eq.master) then
       allocate( idat_all(3,nprocs), sdat_all(4,nprocs) )
    endif
    call MPI_GATHER (idat,3,MPI_INTEGER,idat_all,3,MPI_INTEGER,master,MPI_COMM_WORLD, ierr)
    call MPI_GATHER (sdat,4*180,MPI_CHAR,sdat_all,4*180,MPI_CHAR,master,MPI_COMM_WORLD, ierr)
    if (myrank.eq.master) then
       open(999,FILE='_processes_',status='unknown')
       do i=1,nprocs
          idat(:)=idat_all(:,i)
          sdat(:)=sdat_all(:,i)
          WRITE(999,'(I6,1x,I6,1x,I6,1x,A,A,A,1x,A,A,A,1x,A,A,A,1x,A,A,A)') idat(1),idat(2),idat(3),'"',TRIM(sdat(1)),'"','"',TRIM(sdat(2)),'"','"',TRIM(sdat(3)),'"','"',TRIM(sdat(4)),'"'
       enddo
       close(999)
    endif
  end subroutine WriteProcesses
  
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

  subroutine mpi_bcast_V_vsp(Vr, NRAD,NAT,jspin)
    IMPLICIT NONE
    INTEGER, intent(in)   :: NRAD, NAT,jspin
    REAL*8, intent(inout) :: VR(NRAD,NAT,jspin)
  end subroutine mpi_bcast_V_vsp
  
  subroutine mpi_bcast_V_vns(Vlm,lmmax,lm_stored,LM,nrad,lmmx)
    IMPLICIT NONE
    REAL*8, intent(inout):: Vlm(nrad,1:lmmx)
    INTEGER, intent(inout):: lmmax, lm_stored, LM(2,lmmx)
    INTEGER,intent(in)    :: nrad, lmmx
  end subroutine mpi_bcast_V_vns
  
  subroutine mpi_SendReceive12(ss,weight,bname, en,vt, nv_plus_nnrlo,neig,kt,isi,jspin,must_compute)
    IMPLICIT NONE
    real*8, intent(in)  :: ss(3), weight
    character*10, intent(in) :: bname
    integer, intent(in) :: nv_plus_nnrlo, neig, isi, jspin
    integer, intent(in) :: kt(3,nv_plus_nnrlo)
    logical, intent(in) :: must_compute
    !
    real*8, intent(in)  :: en(neig)
    complex*16, intent(in) :: vt(nv_plus_nnrlo,neig)
  end subroutine mpi_SendReceive12

  subroutine mpi_SendReceive3(vnorm, neig, must_compute)
    IMPLICIT NONE
    real*8, intent(in)  :: vnorm(neig,2)
    integer, intent(in) :: neig
    logical, intent(in) :: must_compute
  end subroutine mpi_SendReceive3
  
  subroutine Scatter_Vector_data(ierr)
    IMPLICIT NONE
    INTEGER, intent(out) :: ierr
    ierr=0
  end subroutine Scatter_Vector_data
  
  subroutine MPI_Bcast_dims(nmat,nume,nkp)
    IMPLICIT NONE
    INTEGER, intent(inout) :: nmat, nume, nkp
  end subroutine MPI_Bcast_dims

  subroutine WriteProcesses()
  end subroutine WriteProcesses
  
  subroutine MPI_Bcast_dims_para(nmat,nume,nkp)
    IMPLICIT NONE
    INTEGER, intent(inout) :: nmat, nume, nkp
  end subroutine MPI_Bcast_dims_para
#endif

  
end MODULE mpi

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

subroutine add_xend(xend,filename)
  IMPLICIT NONE
  CHARACTER*180, intent(inout) :: filename
  CHARACTER*180, intent(in) :: xend
  ! locals                                                                                                                                                                  
  INTEGER :: ln
  ln=len_trim(filename)
  if (filename(ln-1:ln).eq.'_x') then
     filename(:) = filename(1:ln-2)//xend
  else
     filename(:) = filename(1:ln)//xend
  endif
end subroutine add_xend
