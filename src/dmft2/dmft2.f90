! @Copyright 2007 Kristjan Haule
! 

PROGRAM DMFT2
  USE param, only: nsym, filename_V_nsh, filename_V_sph, filename_V_vns, lxdos, nemax0, nemin0, nkpt, nmat, nume, numkpt, nwave
  USE defs, only: Ry2eV
  USE reallocate, only: 
  USE ams, only: init_ams
  USE structure, only : natm, mult, rotij, tauij, pos, rotloc, BR1, ReadStructure,  WriteInfoStructure, DeallocateStructure, rot_spin_quantization, lattic
  USE xa2, only: init_xa2
  USE com, only: ef, elecn, emax, emin, nat, rel
  USE char, only: efmod, modus, modus1
  USE dmfts, only: Qcomplex, iso, natom, crotloc, DM_EF, DM_Emin, DM_Emax, shft, iatom, Read_indmf2, Read_indmfl, DeallocateDmf
  USE sym2,  only: init_sym2, iord, iz, tau
  USE com_mpi
  IMPLICIT NONE
  REAL*8        :: GMAX
  CHARACTER*4   :: RCFILE
  CHARACTER*5   :: COORD, CHARS, DUMMY
  CHARACTER*11  :: STATUS,FORM                                      
  CHARACTER*67  :: ERRMSG
  CHARACTER*80  :: DEFFN, ERRFN
  CHARACTER*180 :: FNAME,FNAME1,fnamehelp
  CHARACTER*200 :: vec_up, vec_dn, ene_up, ene_dn
  LOGICAL       :: helpfiles!, Qident
  INTEGER       :: i, loro, latom
  REAL*8        :: sw(3), sumw, BR1inv(3,3)
  INTEGER       :: INDEX, INDEX1, MI, JR, JC, fh_dos, lngth, ios, ivector, inkp, ik_start
  REAL*8        :: POST(3), BKRLOC(3,3)
  REAL*8        :: TTOTAL_w, TFOUR_w, PFOUR_w, TCLM_w, PCLM_w, TFERMI_w, PFERMI_w
  REAL*8        :: cost, cosf, esepermin, eseper, eseper0, fi, sint, sinf, theta
  REAL*8        :: Tstart_w, ttime, tstart, tfour, tclm, tfermi, ttotal, PFERMI, PCLM, PFOUR
  INTEGER       :: IDUMMY, nspin1, icase, INFO, IRECL, IUNIT, jatom, kxmax, kymax, kzmax
  INTEGER       :: identity3(3,3), diff, ig, j1, j2
  save ttotal,tfermi,tclm,tfour,tstart,errfn
  
  CALL start_MPI()

  fh_dos = 500
  Qcomplex=.FALSE.
  nspin1=2
  RCFILE='NOFI'
  fnamehelp=''
  helpfiles=.false.
  iso=1
  !-----------------------------------------------------------------------  
  !                                                                       
  CALL init_ams

  CALL GTFNAM(DEFFN,ERRFN)

  if (myrank.EQ.master) CALL ERRFLG(ERRFN,'Error in DMFT2')

  OPEN (1,FILE=DEFFN,STATUS='OLD',ERR=910)
  
  DO
     READ (1,*,END=20,ERR=960) IUNIT,FNAME,STATUS,FORM,IRECL
     lngth = Len_Trim(FNAME)
     if ( FNAME(lngth-1:lngth).EQ.'_x' ) then
        FNAME = FNAME(:lngth-2)
        if(iunit.eq.9 .or. iunit.eq.10) vector_para = .True.
     else if (iunit.EQ.6) then
        if (myrank.eq.master) then
           print *, 'Opening', FNAME
           OPEN(IUNIT,FILE=FNAME,STATUS=STATUS,FORM=FORM,ERR=920)
        else if (fastFilesystem) then
           CALL FilenameMPI(FNAME)    ! Each processor will have different information file
           OPEN(IUNIT,FILE=FNAME,STATUS=STATUS,FORM=FORM,ERR=920)
        endif
        cycle
     else if (iunit.eq.8 .or. iunit.eq.21 .or. iunit.eq.fh_dos .or. iunit.eq.15) then 
        if (myrank.eq.master) OPEN(IUNIT,FILE=FNAME,STATUS=STATUS,FORM=FORM,ERR=920)
     else if (iunit.eq.29 .or. iunit.eq.30 .or. iunit .eq. 31 .or. iunit.eq.18 .or. iunit.eq.9902 .or. iunit.eq.9919) then
        ! Energy file and potential file! Do not do anything yet
        if (iunit.eq.18) filename_V_sph = FNAME
        if (iunit.eq.9902) filename_V_nsh = FNAME
        if (iunit.eq.9919) filename_V_vns = FNAME
     else
        OPEN(IUNIT,FILE=FNAME,STATUS=STATUS,FORM=FORM,ERR=920)
     endif
     
     if(iunit.eq.9)  VECFN(1)=FNAME
     if(iunit.eq.10) VECFN(2)=FNAME
     if(iunit.eq.30) VECFN(3)=FNAME
     if(iunit.eq.31) VECFN(4)=FNAME
     
     if(iunit.eq.18) then
        do i=180,5,-1
           if(fname(i:i).ne.' ') then
              if(fname(i-2:i).eq.'vsp') nspin1=1
              CYCLE
           endif
        enddo
     endif
     if(iunit.eq.12) then
        do i=180,6,-1
           if(fname(i:i).ne.' ') then
              if(fname(i-5:i).eq.'normso') then
                 Qcomplex=.TRUE.
                 iso=2
              endif
           endif
        enddo
     endif
     if(iunit.eq.13) fname1=fname
     close(13)
     if(iunit.eq.3) then
        read(iunit,'(A5)',end=12) CHARS
        Qcomplex=.true. ! switches Qcomplex  on if case.in1c exists and is nonzero
12      continue
     end if
  ENDDO
20 CONTINUE
  CLOSE (1)
  ! finished dmtft2.def
  
  if (nprocs.EQ.1) then
     write(6,*)'Running DMFT in single processor mode'
     !write(6,*)' '
  else if (Qprint) then
     write(6,*)'Running DMFT in mpi mode'
     !write(6,*)' '
  endif
  
  call ReadStructure(20,nat,rel,lxdos)
  ! natm = sum(mult)  ! ndif==natm
  !                                                                       
  !.....READ IN SYMMETRY OPERATIONS AND NONPRIMITIVE TRANSLATIONS         
  READ(20,'(i4)') IORD                                                  
  CALL init_sym2(iord)
  DO  ig=1,iord
     READ(20,111) ( (iz(J1,J2,ig),J1=1,3),TAU(J2,ig), J2=1,3 )
  enddo
111 FORMAT(3(3I2,F11.8/))                                              
  !
  ! The symmetry operations iz have just been read. Now correct it if necessary:
  identity3(:,:)=0
  identity3(1,1)=1
  identity3(2,2)=1
  identity3(3,3)=1
  do ig=1,iord
     diff = sum(abs(iz(:,:,ig)-identity3(:,:)))
     if (diff.eq.0) exit
  enddo
  !print *, 'Found identity at ig=', ig
  if (ig.ne.1) then
     if (ig.eq.iord+1) then
        WRITE(6,*) 'ERROR : Could not find indenity among symmetry operations'
     endif
     CALL swapGroup(iz(:,:,ig),tau(:,ig),iz(:,:,1),tau(:,1))
  endif
  
  if (Qprint) CALL WriteInfoStructure(6, nat)
  
  ALLOCATE (rotij(3,3,natm),tauij(3,natm))
  nsym=iord
  !.....DEFINE ROTATION MATRICES IN NONSYMMORPHIC CASE                    
  CALL ROTDEF(iz,tau,iord,nat,pos,natm,rotij,tauij,mult,lattic)
  
  !.....READING case.in2 file
  READ(5,1003)  MODUS, MODUS1, coord                                        
  if (Qprint) write(6,*)  ' Modus: ', MODUS
  esepermin=0.5d0
  eseper0=0.05d0

  READ(5,*,err=4716)  EMIN, ELECN, esepermin, eseper0
  if (Qprint) write(6,*) 'no read error'
4716 continue
  if(esepermin.lt.1.d-3)   esepermin=0.5d0

  READ(5,1004) efmod,ef
  DO I=1,NAT
     READ(5,*) IDUMMY
  ENDDO
  READ(5,*,END=4715) GMAX

  READ(5,'(A4)',END=4715) RCFILE
4715 CONTINUE
  
  EMAX=ef+100.  ! Just take very large emax for now. It will reduce it later on with DMFT cutoff

  REWIND 5
  READ(5,1234) DUMMY
  if (myrank.EQ.master) write(21,720)gmax
  
  if (Qprint) WRITE(6,*) 'RECPR LIST: ',RCFILE
  if (Qprint) WRITE(6,870)  COORD
  
  CALL LATGEN(nat)
  
  CALL CPUTIM(ttime)
  tstart = ttime
  call walltim(ttime)
  tstart_w = ttime

!!!! This generates reciprocal vectors for dmft2 calculations
!!!! These are in general different than what is needed in lapw1 for the hamiltonian.
  CALL RECFIL(FNAME1,GMAX,RCFILE,KXMAX,KYMAX,KZMAX,NWAVE)
  
  if (vector_para) then
     FNAME = '_processes_'
     !WRITE(6,*) 'Trying to open ', FNAME
     open(999,FILE=FNAME,STATUS='OLD',ERR=77,IOSTAT=ios) ! opens _processes_, hence this parallel dmft_lapw1 run, and not parallel w2k_lapw1 run
     call Scatter_Vector_data()
     goto 88
77   CONTINUE     
     
     CALL FilenameMPI2(FNAME)
     !WRITE(6,*) 'Trying to open ', FNAME
     open(999,FILE=FNAME,STATUS='OLD',ERR=88,IOSTAT=ios) ! opens _processes_x  
     nvector=0
     vectors=0
     DO 
        READ (999,*,END=88,ERR=970) ivector,inkp,ik_start,vec_up,vec_dn,ene_up,ene_dn  ! can jump to 20
        nvector = nvector+1
        vectors(nvector,1) = ivector       ! the successive number of vector file          
        vectors(nvector,2) = inkp          ! number of k-points                            
        vectors(nvector,3) = ik_start      ! what is the first k-point in this vector file 
        fvectors(nvector,1) = vec_up       ! filenames....                                 
        fvectors(nvector,2) = vec_dn
        fvectors(nvector,3) = ene_up
        fvectors(nvector,4) = ene_dn
        if (Qprint) WRITE(6,*) 'FNAMES=', fvectors(nvector, :)
     ENDDO
88   CONTINUE
     if (Qprint) print *, 'Vector files are read in parallel mode'
     close(999)
  else
     nvector=1
  endif

!!!! Begining of DMFT specific reading !!!!!
  ! Reading case.indmf2
  call Read_indmf2(2)
  ! Reading case.indmfl
  call Read_indmfl(7, 6, nemin0, nemax0, loro, rotloc, mult, nat, natm, (modus.eq.'FOR '), Qprint )

  ! For spin-orbit coupling crates transformation "rot_spin_quantization", which rotates from spin quantization axis to global cartesian coordinate system
  !....reading *.inso
  if (iso.eq.2) then
     do i=1,3
        read(4,*)
     end do
     read(4,*)sw(1),sw(2),sw(3)
     if (Qprint) WRITE(6,*) 'so-dir=', sw(:)
     call angle(sw,theta,fi)
     if (Qprint) write(6,126) theta,fi
     !.... test of collinearity of the z-axis with spin quantization axis
     cost=dcos(theta)
     cosf=dcos(fi)
     sint=dsin(theta)
     sinf=dsin(fi)
     rot_spin_quantization(1,1)=cosf*cost
     rot_spin_quantization(1,2)=-sinf
     rot_spin_quantization(1,3)=cosf*sint
     rot_spin_quantization(2,1)=sinf*cost
     rot_spin_quantization(2,2)=cosf
     rot_spin_quantization(2,3)=sinf*sint
     rot_spin_quantization(3,1)=-sint
     rot_spin_quantization(3,2)=0
     rot_spin_quantization(3,3)=cost
  else
     rot_spin_quantization(:,:)=0.d0
     rot_spin_quantization(1,1)=1.d0
     rot_spin_quantization(2,2)=1.d0
     rot_spin_quantization(3,3)=1.d0
  end if  ! spin-orbit
  
  ! If 'EF.dat' exists, we should use it instead of case.indm2
  OPEN(1000,FILE='EF.dat',STATUS='OLD',ERR=991)
  READ(1000,*) DM_EF
  CLOSE(1000)
  if (myrank.EQ.master) WRITE(6,*) 'Found EF.dat and EF is ', DM_EF
  DM_EF = DM_EF/Ry2eV
991 CONTINUE


  CALL INIT_ENERGY_FILE(sumw, nkpt, numkpt, nmat, nume)  ! This computes weights and also reads all k-points
  
  if (Qprint) write(6, 877) DM_Emin,DM_Emax  ! Energy window
  if (Qprint) write(6, 887) natom      ! Projected density of states calculation for xxx atoms

  if (.not. vector_para) then
     vectors(1,2)=numkpt
  endif

  
  CALL init_xa2(nume,nkpt)  !  allocates WEIGHT, E, NE

  eseper = DM_EMIN
  
  if (Qprint) WRITE(6,*) 'sumw=', sumw, 'eseper=', eseper
  
  if (Qprint)  WRITE(6,1060)  ELECN,EF
  
  if (myrank.EQ.master) WRITE(21,1060) ELECN,0.0
  
  if (myrank.EQ.master) then
     !call INVERSSYMDEF(BR1,BR1inv)
     call inv_3x3(BR1,BR1inv)
     INDEX = 0
     DO JATOM=1, NAT
        INDEX1 = INDEX+1
        DO MI=1, MULT(JATOM)
           INDEX = INDEX + 1

           POST = matmul(POS(:,index1), ROTIJ(:,:,index))
           POST = POST + TAUIJ(:,index) + shft(INDEX,:)
           
           write(6,'(A,I2,A)',advance='no') 'Actual position of atom ', INDEX, ' is:'
           write(6,'(f10.6,2x,f10.6,2x,f10.6)') POST
           write(6,'(A)') 'The local coordinate system (withoth user rotation "locrot") is:'
           do i=1,3
              write(6,'(3F8.2,3x,F7.3)') ROTIJ(i,:,index), TAUIJ(i,index)
           enddo
        ENDDO
     ENDDO
     
     write(6,*)
     write(6,'(A)',advance='no') 'The local coordinate systems on correlated atoms '
     write(6,'(A)') 'including users rotation ("locrot") are'
     DO icase=1,natom
        latom = iatom(icase)
        write(6, '(A,I3,1x,A,1x,I3,1x,A)') 'catom', icase, '( atom', latom, ')'
        BKRLOC = matmul(crotloc(:,:,icase),matmul(BR1, matmul(rotij(:,:,latom),BR1inv) ))
        DO JR=1,3
           WRITE(6, '(3F10.5)') (BKRLOC(JR,JC),JC=1,3) 
        ENDDO
        WRITE(6,*)
     ENDDO
     write(6,*)
  endif
  
  CALL CPUTIM(TTIME)
  TFERMI=TTIME                                        
  call walltim(TTIME)
  TFermi_w=TTIME

  !.....CALCULATE CHARGE DENSITY CLM(R) IN SPHERES,  PARTIAL CHARGES      
  call l2main(coord,NSPIN1,sumw,TCLM,TCLM_w,TFOUR,TFOUR_w) !!!! SHOULD BE CHANGED

  call DeallocateStructure()
  call DeallocateDmf()
  
  !.....CALCULATE CPUTIME REQUIRED                                        
  TTOTAL=TFOUR-TSTART                                               
  TFOUR=TFOUR-TCLM                                                  
  PFOUR=TFOUR/TTOTAL*100.                                           
  TCLM=TCLM-TFERMI                                                  
  PCLM=TCLM/TTOTAL*100.                                             
  TFERMI=TFERMI-TSTART                                              
  PFERMI=TFERMI/TTOTAL*100.                                         
  if (Qprint) then
     WRITE(6,2000)                                                     
     WRITE(6,2010) TTOTAL,100.0
     WRITE(6,2020) TFERMI,PFERMI
     WRITE(6,2030) TCLM,PCLM
     WRITE(6,2040) TFOUR,PFOUR
  endif
  TTOTAL_w=TFOUR_w-TSTART_w                                               
  TFOUR_w=TFOUR_w-TCLM_w                                                  
  PFOUR_w=TFOUR_w/TTOTAL_w*100.                                           
  TCLM_w=TCLM_w-TFERMI_w                                                  
  PCLM_w=TCLM_w/TTOTAL_w*100.                                             
  TFERMI_w=TFERMI_w-TSTART_w                                              
  PFERMI_w=TFERMI_w/TTOTAL_w*100.                                         
  if (Qprint) then
     WRITE(6,2001)                                                     
     WRITE(6,2010) TTOTAL_w,100.0
     WRITE(6,2020) TFERMI_w,PFERMI_w
     WRITE(6,2030) TCLM_w,PCLM_w
     WRITE(6,2040) TFOUR_w,PFOUR_w
  endif
  if (myrank.EQ.master) CALL ERRCLR(ERRFN)

  CALL stop_MPI()
  
  STOP !' DMFT2 END'                                                 
  !
  !        error handling
  !
910 INFO = 1
  !
  !        'lapw2.def' couldn't be opened
  !
  WRITE (ERRMSG,9000) FNAME
  CALL OUTERR('DMFT2',ERRMSG)
  GOTO 999
920 INFO = 2
  !
  !        file FNAME couldn't be opened
  !
  WRITE (ERRMSG,9010) IUNIT
  CALL OUTERR('DMFT2',ERRMSG)
  WRITE (ERRMSG,9020) FNAME
  CALL OUTERR('DMFT2',ERRMSG)
  WRITE (ERRMSG,9030) STATUS, FORM
  CALL OUTERR('DMFT2',ERRMSG)
  GOTO 999
960 INFO = 7
  !
  !        Error reading file 'lapw2.def'
  !
  WRITE (ERRMSG,9040) FNAME
  CALL OUTERR('DMFT2',ERRMSG)
  GOTO 999

970 INFO = 8
  !
  ! Error in parallel: could not find .processes.x
  !
  WRITE (ERRMSG,'(A,A)')  'file open error:', trim(ADJUSTL(FNAME))
  CALL OUTERR('DMFT',ERRMSG)
  GOTO 999
999 STOP 'DMFT2 - Error. Check file dmft2.error'
  !                                                                       
720 FORMAT(':GMA  :',' POTENTIAL AND CHARGE CUT-OFF',f7.2,' Ry**.5')
870 FORMAT(3X,'TYPE OF COORDINATES IN DSPLIT= ',A5)                   
1003 FORMAT(A5,a5,a5)                                                        
1004 FORMAT(A5,f10.5)                                                        
1060 FORMAT(//,':NLDA  :',1X,'NUMBER OF ELECTRONS',10X,'=',F8.3,//   ,':FLDA  :',1X,'F E R M I - ENERGY',11X,'= ',F9.5)            
1234 FORMAT(//,1A)
2000 FORMAT(//,3X,'=====>>> CPU TIME SUMMARY',/)                       
2001 FORMAT(//,3X,'=====>>> WALL TIME SUMMARY',/)                       
2010 FORMAT(12X,'TOTAL       : ',F8.1,5X,'... ',F4.0,' PERCENT')       
2020 FORMAT(12X,'PART FERMI  : ',F8.1,5X,'... ',F4.0,' PERCENT')       
2030 FORMAT(12X,'PART CLM    : ',F8.1,5X,'... ',F4.0,' PERCENT')       
2040 FORMAT(12X,'PART FOURIR : ',F8.1,5X,'... ',F4.0,' PERCENT')       
9000 FORMAT('can''t open definition file ',A40)
9010 FORMAT('can''t open unit: ',I2)
9020 FORMAT('       filename: ',A50)
9030 FORMAT('         status: ',A,'  form: ',A)
9040 FORMAT('Error reading file: ',A47)
877 FORMAT(' Energy/band window:',f10.4,' < E < ',f10.3)
887 FORMAT(' Greens function calculation for',i3,' atoms')
126 FORMAT(' Magnetic system with s-o coupling; M theta, phi:',2f8.4)
END PROGRAM DMFT2


subroutine swapGroup(iz1,tau1,iz2,tau2)
  IMPLICIT NONE
  INTEGER, intent(inout) :: iz1(3,3), iz2(3,3)
  REAL*8,  intent(inout) :: tau1(3), tau2(3)
  ! locals
  INTEGER :: izx(3,3)
  REAL*8  :: taux(3)
  izx(:,:) = iz1(:,:)
  iz1(:,:) = iz2(:,:)
  iz2(:,:) = izx(:,:)
  taux(:) = tau1(:)
  tau1(:) = tau2(:)
  tau2(:) = taux(:)
end subroutine swapGroup
