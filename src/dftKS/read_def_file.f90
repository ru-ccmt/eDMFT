subroutine read_def_file(deffn, nloat, nohns, iter, readHinv, force, runorb, nmat_only, dokorig, writeham, ERRMSG, info)
  !        open all files listed in 'lapw1.def'
  use param, only   : filename_V_sph, filename_V_vns, filename_vector, filename_energy
  use mpi, only : Qprint, FilenameMPI2, FilenameMPI, vector_para, myrank, master, cpuID, nprocs
  IMPLICIT NONE
  CHARACTER*(*), intent(in) :: deffn
  INTEGER, intent(out) :: nloat, info
  logical, intent(out) :: nohns, iter, readHinv, force, runorb, nmat_only, dokorig, writeham
  CHARACTER*67, intent(out) :: ERRMSG
  !
  INTERFACE 
     FUNCTION find_nloat(fh)
       INTEGER :: find_nloat
       INTEGER, INTENT(IN) :: fh
     END FUNCTION find_nloat
  END INTERFACE
  !locals
  INTEGER       :: IUNIT, IRECL, myid, lngth
  CHARACTER*180 :: FNAME, xend
  CHARACTER*11  :: STATUS, FORM
  LOGICAL       :: Hinv_open
  CHARACTER*80  :: myids
  !
  myid=0
  info=0
  !
  vector_para=.false.
  !
  nohns=.false.
  iter=.false.
  readHinv=.false.
  force=.false.
  runorb=.false.
  nmat_only=.false.
  dokorig=.false.
  filename_V_vns=''
  filename_V_sph=''
  OPEN (1,FILE=DEFFN,STATUS='OLD',ERR=910)
  DO
     READ (1,*,END=20,ERR=960) IUNIT,FNAME,STATUS,FORM,IRECL
     lngth = Len_Trim(FNAME)
     
     if (iunit.eq.10) filename_vector=FNAME  ! vector file
     if (iunit.eq.11) filename_energy=FNAME  ! energy file

     if ((vector_para .or. FNAME(lngth-1:lngth).EQ.'_x') .and. (iunit.eq.10 .or. iunit.eq.11)) then
        ! It is enugh that one of the two files has '_x' and we will try vector_para.
        xend=''
        if (nprocs > 1) then
           xend = '_'//trim(ADJUSTL(cpuID))
           vector_para = .True.
        endif
        call add_xend(FNAME, xend ) ! remove '_x' if present, and then add _{myrank}
        if (iunit.eq.10) filename_vector=FNAME  ! correct name of vector file
        if (iunit.eq.11) filename_energy=FNAME  ! correct name of energy file
        OPEN (IUNIT,FILE=FNAME,STATUS=STATUS,FORM=FORM,ERR=920)
        cycle
     endif

     if (iunit.EQ.6 .or. iunit.eq.21) then
        if (Qprint) then
           if (myrank.ne.master) CALL FilenameMPI(FNAME)    ! Each processor will have different information file
           ! master will use usual filename, the rest will append .myrank
           OPEN(IUNIT,FILE=FNAME,STATUS=STATUS,FORM=FORM,ERR=920)
        endif
        cycle
     endif
     if (iunit.eq.10 .or. iunit.eq.11 .or. iunit.eq.71 .or. iunit.eq.55) then
        ! These are output files and should be opened only on master node
        if (iunit.eq.71) force=.true.
        if (myrank.EQ.master) OPEN (IUNIT,FILE=FNAME,STATUS=STATUS,FORM=FORM,ERR=920)
        cycle
     endif
     if (iunit.eq.18 .or. iunit.eq.19) then
        ! These are input files that will be opened later
        if (iunit.eq.18) filename_V_sph =FNAME
        if (iunit.eq.19) filename_V_vns =FNAME
        cycle
     endif
     if (iunit.eq.201) then
        readHinv=.true.
        CYCLE
     endif
     if (iunit.eq.200) then   ! Seems to be something for iterative deiagonalization in jacdavblock.FP
        write(myids,*) myid
        FNAME=trim(FNAME)//"_proc_"//trim(adjustl(myids))
        if (status.eq.'unknown') readHinv=.true.
     endif
     
     OPEN (IUNIT,FILE=FNAME,STATUS=STATUS,FORM=FORM,ERR=920)
     if(iunit.eq.72) nmat_only=.true.
     if(iunit.eq.7) runorb=.true.
     if(iunit.eq.97)   then
        close(97)
        nohns=.true.
     end if
     if(iunit.eq.98) then
        iter=.true.
     endif
     if(iunit.eq.3) then
        dokorig=.true.
     endif
     if(iunit.eq.12) then
        writeham=.true.
     endif
     if(iunit.eq.5) then
        nloat = find_nloat(5)
     endif
  ENDDO
20 CONTINUE
  CLOSE (1)
  
  !if(.not.readHinv) then
  !   inquire(unit=200, opened=Hinv_open)
  !   if(Hinv_open) then
  !      close(unit=200, status='DELETE')
  !   endif
  !endif
  
  return
  
910 CONTINUE  ! 'lapw1.def' couldn't be opened
  info = 1
  WRITE (ERRMSG,9000) trim(DEFFN)
  return
920 CONTINUE  !  file FNAME couldn't be opened
  info = 2
  WRITE(ERRMSG,'(A,A)') 'cant open file:',trim(FNAME)
  WRITE(6,*) ERRMSG
  return
960 CONTINUE
  info = 7 !  Error reading file 'lapw1.def'
  WRITE (ERRMSG,9040) trim(FNAME)
  return
9000 FORMAT('can''t open definition file ',A40)
9040 FORMAT('Error reading file: ',A47)
end subroutine read_def_file
