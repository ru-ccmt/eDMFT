!***************************
! Reading lapwso.def  file *
!***************************
SUBROUTINE read_def_file(deffn, Qcomplex, iorbpot,ERRMSG,info)
  !        open all files listed in 'lapw1.def'
  USE param, ONLY: filename_V_sp, filename_V_vns
  use mpi, only : Qprint, FilenameMPI2, FilenameMPI, vector_para, myrank, master, stop_MPI, filename_vector, filename_energy, filename_vectorso, filename_energyso, filename_norm
  IMPLICIT NONE
  CHARACTER*(*), intent(in) :: deffn
  LOGICAL, intent(out) :: Qcomplex
  INTEGER, intent(out) :: iorbpot
  CHARACTER*67, intent(out):: ERRMSG
  INTEGER, intent(out)     :: info
  !locals
  INTEGER       :: IUNIT, IRECL, lngth
  CHARACTER*180 :: FNAME
  CHARACTER*11  :: STATUS, FORM
  !
  !print *, 'myrank=', myrank
  Qcomplex = .false.    ! for inversion symmetry breaking, complex case
  vector_para = .false.
  iorbpot=0          ! for oprbital potential
  OPEN(1,FILE=DEFFN,STATUS='OLD',ERR=910)
  DO
     READ(1,*,END=20,ERR=910) IUNIT,FNAME,STATUS,FORM,IRECL
     lngth = Len_Trim(FNAME)
     !        check for in1c file for complex case
     if ( FNAME(lngth-3:lngth).EQ.'in1c' ) Qcomplex=.true. ! no inversion symmetry. Eigenvectors are complex.
     
     ! files to read
     if (iunit.eq.9) filename_vector(1)=FNAME    ! vector
     if (iunit.eq.10) filename_vector(2)=FNAME   ! vectorup
     if (iunit.eq.54) filename_energy(1)=FNAME   ! energy
     if (iunit.eq.55) filename_energy(2)=FNAME   ! energyup
     ! files to write
     if (iunit.eq.41) filename_vectorso(1)=FNAME ! vectorsodn
     if (iunit.eq.42) filename_vectorso(2)=FNAME ! vectorso
     if (iunit.eq.51) filename_energyso(1)=FNAME ! energysodn
     if (iunit.eq.52) filename_energyso(2)=FNAME ! energyso
     if (iunit.eq.53) filename_energyso(3)=FNAME ! energydum
     if (iunit.eq.45) filename_norm(1)=FNAME     ! normsodn
     if (iunit.eq.46) filename_norm(2)=FNAME     ! normsoup

     if (vector_para .and. (iunit.eq.9 .or. iunit.eq.10 .or. iunit.eq.54 .or. iunit.eq.55 .or. iunit.eq.41 .or. iunit.eq.42 .or. iunit.eq.51 .or. iunit.eq.52 .or. iunit.eq.53 .or. iunit.eq.45 .or. iunit.eq.46))then
        ! If we are storing vector files parallely, we will open these files later, when we read _processes_
        cycle
     endif
     if ( FNAME(lngth-1:lngth).EQ.'_x') then  ! This is how we recognize that we have vector_para
        ! We will rea/write vector and energy files for each processor separately. Each mpi process has his own vector file
        if (iunit.eq.9 .or. iunit.eq.10 .or. iunit.eq.54 .or. iunit.eq.55) then
           vector_para = .True.
           if (Qprint) WRITE(6,*) 'Found para=', FNAME, vector_para
        endif
        cycle
     endif
     if (iunit.EQ.6 .or. iunit.eq.8) then ! case.outputso or case.scfso
        if (Qprint) then
           if (myrank.ne.master) CALL FilenameMPI(FNAME)    ! Each processor will have different information file
           ! master will use usual filename, the rest will append .myrank
           OPEN(IUNIT,FILE=FNAME,STATUS=STATUS,FORM=FORM,ERR=920)
        endif
        cycle
     endif
     if (iunit.eq.41.or.iunit.eq.42.or.iunit.eq.44.or.iunit.eq.45.or.iunit.eq.46.or.iunit.eq.51.or.iunit.eq.52.or.iunit.eq.53) then
        ! These are output files and should be opened only on master node
        if (iunit.eq.44) cycle                      ! vec is empty
        if (myrank.EQ.master) OPEN (IUNIT,FILE=FNAME,STATUS=STATUS,FORM=FORM,ERR=920)
        cycle
     endif
     if (iunit.eq.18 .or. iunit.eq.19 .or. iunit.eq.22 .or. iunit.eq.23) then
        ! These are input files that will be opened later
        if (iunit.eq.18) filename_V_sp(1) = FNAME
        if (iunit.eq.19) filename_V_sp(2) = FNAME
        if (iunit.eq.22) filename_V_vns(1) =FNAME
        if (iunit.eq.23) filename_V_vns(2) =FNAME
        if (iunit.eq.18 .or. iunit.eq.19) cycle ! we will open it later
        OPEN(IUNIT,FILE=FNAME,STATUS=STATUS,FORM=FORM,ERR=920)  ! eventually we would like to read also V_vns from atpar library, so this shoudl be commented out
        cycle
     endif
     if(iunit.eq.12)then    
        ! find out whether orb potential should be added (iorbpot=1)
        iorbpot=1
     endif
     
     OPEN(IUNIT,FILE=FNAME,STATUS=STATUS,FORM=FORM,ERR=920)
  ENDDO
20 CONTINUE
  CLOSE (1)
  info=0
  ERRMSG=''
  !call stop_MPI
  !STOP 'debug'
  return
  
910 continue
  WRITE(6,*) ' ERROR IN OPENING LAPWSO.DEF !!!!'
  WRITE (6,*) ' ERROR IN OPENING LAPWSO.DEF !!!!' 
  info = 3
  ERRMSG = 'In read_def_file.f90 the file named lapwso.def could not be opened.'
  return
920 continue
  WRITE(6,*) ' ERROR IN OPENING UNIT:',IUNIT
  WRITE(6,*) '       FILENAME: ',FNAME,'  STATUS: ',STATUS,'  FORM:',FORM
  info = 4
  ERRMSG = 'In read_def_file.f90 was an error reading file lapwso.def'
  return
END SUBROUTINE read_def_file
