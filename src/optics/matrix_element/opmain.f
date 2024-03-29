PROGRAM optic
  use struk
  use potnlc
  use xrpar
  use fermi
  use comi
  !
  !               O P T I C A L   M A T R I X  E L E M E N T S
  !                         
  !     written by Robert Ab
  !
  !     LAST UPDATES: March    1996   Robert Abt: local orbital
  !                   December 1998   Claudia Ambrosch-Draxl: SO
  !                   May      1999   Jan Kunes, spin-polarized SO
  !                   August   1999   cad: modifications
  !ole                November 2001   !ole Bernd Olejnik : p^a_nm NLO
  !
  !          FILE  3         MATRIX ELEMENTS
  !          FILE  5         INPUT
  !          FILE  6         OUTPUT
  !          FILE 10         VECTOR
  !          FILE 18         VSP
  !          FILE 20         STRUCT
  !ole       FILE 24         case.mme       simple momentum matrix elements
  !ole       FILE 25         case.symop     symmetry matrices
  !
  !
  !
  INCLUDE 'param.inc' 
  PARAMETER (NCASE=9)                     
  IMPLICIT REAL*8 (A-H,O-Z)
  !      COMPLEX*16       OPMATX,OPMATY,OPMATZ,MX_,MY_,MZ_,SX_,SY_,SZ_
  CHARACTER*4      LATTIC,IREL,RELA
  CHARACTER*3      MODUS,OUTME,edge    ! edge added LP
  CHARACTER*9      HEADER(NCASE),amodus
  CHARACTER*11     STATUS,FORM                                
  CHARACTER*67     ERRMSG
  CHARACTER*80     TITLE                                    
  CHARACTER*180    FNAME                                    
  CHARACTER*90     DEFFN,ERRFN
  character*300    title1
  !ole ##### Begin #####
  !ole  LOGICAL          REL,LSO,SPIN
  LOGICAL          REL,LSO,SPIN,MME_FLAG 
  !ole #####  End  #####
  INTEGER(4) :: atom_num,xlo	!LP
  real*8           xms(3)
  
  COMMON /outop/ Ncol,icol(NCASE)                                  
  COMMON /LEAD/ K1,K2,KEND,KSTEP,KOUT,KSTOP   
  COMMON /COM/  EMIN,EMAX,ELECN,EULIMIT,EDLIMIT, NK,IOUT,NSPIN,NAT,NBAND,ix,NB(NKPT),MINWAV,MAXWAV,ifpmat
  COMMON /MIM/  MIMA(2)         
  COMMON /CHAR/   TITLE,LATTIC,MODUS,OUTME                   
  COMMON /GENER/  BR1(3,3),BR2(3,3)                                 
  COMMON /RADFU/  RRAD1(NRAD,LMAX1),RADE1(NRAD,LMAX1),RRAD2(NRAD,LMAX1),RADE2(NRAD,LMAX1)
  COMMON /UHELP/  UDUM(NRAD,6)                      
  COMMON /SYM2/   TAU(3,NSYM),IORD,IMAT(3,3,NSYM)
  COMMON /SYMo/   opimat(3,3,NSYM)
  COMMON /so/     theta,fi
  COMMON /SYMd/   det(NSYM)
  COMMON /CLOGIC/ LSO,SPIN,REL,MME_FLAG
  EXTERNAL GTFNAM                                                     
  DATA   RELA/'RELA'/                                     
  !-----------------------------------------------------------------------       
  !
  LSOUP=.FALSE.
  LSO=.FALSE.
  SPIN=.FALSE.
  !ad
  !ad
  !ad_____________________________ OPEN FILES ____________________________
  !ad
  
  call GTFNAM(DEFFN,ERRFN)
  call ERRFLG(ERRFN,'Error in OPTIC')
  OPEN (1,FILE=DEFFN,STATUS='OLD',ERR=910)
10 CONTINUE
  READ (1,*,END=20,ERR=960) IUNIT,FNAME,STATUS,FORM,IRECL
  OPEN (IUNIT,FILE=FNAME,STATUS=STATUS,FORM=FORM,ERR=920)
  
  !ad
  !jan
  !LO
  if(iunit.eq.10.) then
     do i=180,8,-1
        if(fname(i-5:i).eq.'ctorup') then
           LSOUP=.true.
           SPIN=.true.
        else 
           if(fname(i-5:i).eq.'ctordn') then
              LSOUP=.false.
              SPIN=.true.
           endif
        endif
     enddo
  endif
  !LO
  if((iunit.eq.10.).or.(iunit.eq.11)) then
     do i=180,8,-1
        if(fname(i-7:i).eq.'ctorsoup') then
           LSO=.true.
           SPIN=.TRUE.
           goto 10
        endif
     enddo
  endif
  !jan
  !ad
  if(iunit.eq.10) then
     do i=80,8,-1
        if(fname(i-7:i).eq.'vectorso') then
           LSO=.true.
           goto 10
        endif
     enddo
  endif
  !ad
  
  GOTO 10  
20 CONTINUE
  CLOSE (1)
  !ad
  !ad
  !ad _________________ INITIALIZE VARIABLES AND ARRAYS __________________
  !ad
  EEF=0.0
  DEEF=0.0
  TSTART=0.0
  TOPMAT=0.0
  !ad
  !ad _____________________ READ STRUCTURAL INFORMATION __________________
  !ad
  READ(20,1000) TITLE
  READ(20,1010) LATTIC,NAT,IREL    
  !
  nato=nat
  natti=nat   ! LP: natti in "comi" module 
  allocate (ANAME(nato),POS(3,48*nato),RMT(NATO),V(NATO),ZZ(NATO))
  allocate (IATNR(NATO),MULT(NATO),ISPLIT(NATO))
  allocate (Rnot(NATO),DX(NATO),JRI(NATO),imodus(NATO))
  DO I=1,NATO ! 100                                                    
     MULT(I)=0                                                         
     V(I)=0.D0                                                         
     RMT(I)=0.D0
  ENDDO
  !                    
  !      if(nat.gt.nato) stop 'NATO too small!'
  !.....READ IN LATTICE CONSTANTS                                         
  READ(20,1020) AA,BB,CC,ALPHA(1),ALPHA(2),ALPHA(3)
  IF(ABS(ALPHA(1)).LT.1.D-5) ALPHA(1)=90.0D0
  IF(ABS(ALPHA(2)).LT.1.D-5) ALPHA(2)=90.0D0
  IF(ABS(ALPHA(3)).LT.1.D-5) ALPHA(3)=90.0D0
  !C
  IF(IREL.EQ.'RELA') REL=.TRUE.                                     
  IF(IREL.EQ.'NREL') REL=.FALSE.                                    
  WRITE(6,800)                                                      
  WRITE(6,805)  TITLE                                               
  WRITE(6,810)  LATTIC                                              
  WRITE(6,820)  AA,BB,CC                                            
  WRITE(6,840)  NAT                                                 
  WRITE(6,850)  IREL                                                
  REL=.FALSE.                                                       
  IF(IREL.EQ.RELA) REL=.TRUE.   
  !     INDEX COUNTS ALL ATOMS IN THE UNIT CELL, JATOM COUNTS ONLY THE    
  !     NOTEQUIVALENT ATOMS                                               
  INDEX=0                                                           
  DO JATOM = 1,NAT  ! 50
     INDEX=INDEX+1
     READ(20,1030) IATNR(JATOM),( POS(J,INDEX),J=1,3 ),MULT(JATOM),ISPLIT(JATOM)                                    
     IF (MULT(JATOM).EQ.0) THEN                                     
        WRITE(6,1040) JATOM,INDEX,MULT(JATOM)                       
        STOP ' OPTIC: MULT EQ 0'                                    
     ENDIF
     DO M=1,MULT(JATOM)-1   ! 55
        INDEX=INDEX+1                                            
        READ(20,1031) IATNR(JATOM),( POS(J,INDEX),J=1,3)         
     ENDDO
     READ(20,1050) ANAME(JATOM),JRI(JATOM),Rnot(JATOM),RMT(JATOM),ZZ(JATOM)
     DX(JATOM)=LOG(RMT(JATOM)/Rnot(JATOM)) / (JRI(JATOM)-1)    
     RMT(JATOM)=Rnot(JATOM)*EXP( DX(JATOM)*(JRI(JATOM)-1) ) 
     do i=1,3      
        READ(20,*)  
     end do
  ENDDO
  !.....reading symmetry operations
  READ(20,12) IORD
  do j=1,IORD
     READ(20,11) ( (IMAT(j1,j2,j),j1=1,3),TAU(j2,j), j2=1,3 )
  end do
12 FORMAT(I4)
11 FORMAT(3(3I2,F15.8/))                                                                   
  MODUS='ALL'
  !ad
  !ad   KUPLIMIT gives maximal k-points to calculate 
  !ad   KFIRST gives the k-point to start calculation 
  !ad
  READ(5,*) KUPLIMIT,KFIRST
  kuplimit=kuplimit+kfirst-1
  !ad
  EDLIMIT=EEF-DEEF                                           
  EULIMIT=EEF+DEEF                                           
  nbvalmax=9999999        
  read(5,'(a)') errmsg                               
  READ(errmsg,*,end=1,err=1) EMIN,EMAX,nbvalmax
1 print *, 'emin,emax,nbvalmax', emin,emax,nbvalmax
  errmsg=''
  !      READ(5,*) EMIN,EMAX
  !      READ(5,*)
  xmcd=0
  read(5,'(a)') errmsg   
  do i=1,10
     if(errmsg(i:i+3).eq.'XMCD') then
        READ(errmsg(i+4:67),*,end=1,err=1) atom_num,edge
        xmcd=1
        print*, 'XMCD selected for atom',atom_num,edge
     endif
  enddo
  if(xmcd.eq.1)then
     select case (edge)
     case ('K  ')
        core_name1='1s'
        xlo=1
     case ('L1 ')
        core_name1='2s'
        xlo=1
     case ('L23')
        core_name1='2p'
        core_name2='2ps'
        xlo=2
     case ('M1 ')
        core_name1='3s'
        xlo=1
     case ('M23')
        core_name1='3p'
        core_name2='3ps'
        xlo=2
     case ('M45')
        core_name1='3d'
        core_name2='3ds'
        xlo=2
     case default
        write(6,*)'Error: check edge!'
        stop
     end select
     !ad
     !jan
     !jan  isp=0 non spinpol isp=1 spinpol
     !jan
     READ(5,*) NCOL
  else
     read(errmsg,*) ncol
  endif
  !ole ##### Start #####
  !ole  output matrix elements p^a_nm (a=x,y,z)
  !ole
  header(1)='Re <x>'
  header(2)='Im <x>'
  header(3)='Re <y>'
  header(4)='Im <y>'
  header(5)='Re <z>'
  header(6)='Im <z>'
  
  DO i=1,NCASE
     icol(i)=i
  END DO

  SELECT CASE (NCOL)
     CASE (0) ! ### get simple momentum matrix elements
        MME_FLAG=.TRUE.
        NCOL=6
        WRITE(24,330) ncol,(header(icol(i)),i=1,ncol)
     CASE (1,2,3,4,5,6,7,8,9) 
        MME_FLAG=.FALSE.   ! This is for optics
        do i=1,NCOL
           read(5,*) icol(i)  ! icol(:)=1,2,3
           if (icol(i).gt.NCASE) stop 'col number too large'
        end do
        
        if(xmcd.eq.1) then
           write(*,*)'LSO= ',LSO
           if(.not.LSO) then 
              write(6,*) 'XMCD calculation requires spin-polarized AND spin-orbit setup'
              stop
           endif
           NCOL=2
           do i=1,NCOL
              icol(i)=i
           enddo
        endif
        outme='OFF'
        read(5,'(A3,a9)',end=777) OUTME,amodus
        if(OUTME.eq.'on ' ) OUTME = 'ON '    ! In our case, OUTME -> 'ON'
        moddo = outme  
        nmodus=0
        read(amodus,*,err=777,end=777) nmodus
        if(nmodus.gt.0.and.nmodus.le.nat+1) then
           read(5,*) (imodus(i),i=1,nmodus)
           print*, 'optical matrix elements for atoms:',(imodus(i),i=1,nmodus)
           modus='ato'
        endif
        !ad
        !ad _______________________________ HEADERS ____________________________
        !ad
777     write(6,*) '  MODUS',modus
        if(xmcd.eq.1) then
           modus='ALL'
        endif
        if(modus.eq.'ATO') then
           write(6,*) 'optical matrix elements for atoms:',(imodus(i),i=1,nmodus)
        endif
        header(1)='Re <x><x>'
        header(2)='Re <y><y>'
        header(3)='Re <z><z>'
        header(4)='Re <x><y>'
        header(5)='Re <x><z>'
        header(6)='Re <y><z>'
        header(7)='Im <x><y>'
        header(8)='Im <x><z>'
        header(9)='Im <y><z>'
        if(modus.eq.'ALL') then
           WRITE(3,330) ncol,(header(icol(i)),i=1,ncol),modus
           WRITE(9,330) ncol,(header(icol(i)),i=1,ncol),modus
           if(xmcd.eq.1) then
              WRITE(13,332) LSO,edge,modus
              WRITE(14,332) LSO,edge,modus
           endif
        else
           j=11+13*(ncol+1)+1
           nmodus1=min(80/3,nmodus)
           write(title1,330) ncol,(header(icol(i)),i=1,ncol),modus
           write(title,'(30i3)') (imodus(i),i=1,nmodus1)
           WRITE(3,'(A,a)') title1(1:j),title(1:3*nmodus1)
           WRITE(9,'(A,a)') title1(1:j),title(1:3*nmodus1)
        endif
        if(OUTME.EQ.'ON ') WRITE(4,9044)
        !
        !Clas-start: I have also included ifpmat in COMMON block /COM/
        read(5,*,end=333,err=333) ifpmat
        goto 334
333     ifpmat=0
334     continue
        if(ifpmat.eq.1) write(6,*) " Outputs <j,k|-i*grad|j',k> to case.pmat"
        !Clas-end
        !
        !ole ##### Start #####
        ! default error
     CASE DEFAULT
        GOTO 970
  END SELECT
  !ole ##### End #####
  !
  !ad ____________________________________________________________________
  !ad
  iout=1
  ix=1
  !ad
  if(SPIN) write(6,*) '  spin-polarized calculation'
  if(LSO)  write(6,*) '  spin-orbit coupling included'
  !ad
  !jan
  !jan   read *.inso
  !jan
  call LATGEN (NAT)                                                 
  if (lso.AND.spin) then
     do i=1,3
        read(28,*)
     end do
     read(28,*) xms(1),xms(2),xms(3)
     call angle(xms,theta,fi)
     
     !      read(28,*) theta,fi
111  format(2F10.3)
     write(6,*)'read inso:'
     write(6,*)'theta=',theta
     write(6,*)'phi=',fi
  end if
  !jan
  !
  !      call LATGEN (NAT)                                                 
  call SYMOP (NAT,OUTME)
  !jan
  if (lso.AND.spin) then
     call SYM
  else
     do i=1,IORD
        det(i)=1
     end do
  end if
  !ad
  call CPUTIM(TTIME)
  TSTART0=TTIME
  !ad 
  !ad...CALCULATE OPTICAL MATRIXELEMENTS ................     
  !ad
  if(xmcd.eq.0) then
     call MOM_MAT(KUPLIMIT,KFIRST,IEF,topmat,nbvalmax) 
  else
     do i=1,xlo
        if(i.eq.1) core_name=core_name1
        if(i.eq.2) core_name=core_name2
        !        write(91,*) 'i=',i,'atom_num=',atom_num,'KUPLIMIT=',KUPLIMIT,'KFIRST=',KFIRST,'IEF=',IEF
        call COR_MAT(i,atom_num,KUPLIMIT,KFIRST,IEF)
     enddo
     TOPMAT=TOPMAT-TSTART
  endif
  !ad
  call CPUTIM(TTIME)
  TOUTP=TTIME
  !                                                                       
  !.....CALCULATE CPUTIME REQUIRED                                        
  TTOTAL=TOUTP -TSTART0 
  !      TOPMAT=TOPMAT-TSTART
  WRITE(6,2000) 
  WRITE(6,2020) TOPMAT,TOPMAT/TTOTAL*100.
  WRITE(6,2010) TTOTAL,100.0                                        
  call ERRCLR(ERRFN)
  STOP ' OPTIC END'                
  
  !
  !        error handling
  !
910 INFO = 1
  !
  !     optic.def couldn't be opened
  !
  WRITE (ERRMSG,9000) DEFFN
  call OUTERR('OPTIC',ERRMSG)
  GOTO 999
  !
920 INFO = 2
  !
  !     file FNAME couldn't be opened
  !
  WRITE (ERRMSG,9010) IUNIT
  call OUTERR('OPTIC',ERRMSG)
  WRITE (ERRMSG,9020) FNAME
  call OUTERR('OPTIC',ERRMSG)
  WRITE (ERRMSG,9030) STATUS, FORM
  call OUTERR('OPTIC',ERRMSG)
  GOTO 999
  !
960 INFO = 7
  !
  !     error reading file *.def
  !
  WRITE (ERRMSG,9040) FNAME
  call OUTERR('OPTIC',ERRMSG)
  !ole ##### Start #####
  GOTO 999
  !
  !     bad number of choices in <case>.inop 
  !
970 WRITE (ERRMSG,9050)  
  CALL OUTERR('OPTIC',ERRMSG)
  !ole ##### End #####
999 STOP 'OPTIC - ERROR'
  !ad
  
43 FORMAT(3X,A77)
44 FORMAT(I3,A77)
330 format(10X,I1,10(2X,A9,2X))
  ! 331  format(L5,3X,A3,10(2X,A9,2X))
332 format(L5,3X,A3,'     LEFT        RIGHT       ',A5)
444 format(I4)	 
700 FORMAT(I3,A77)
701 FORMAT('CPU TIME (',2I3,') -> ',F9.2)
800 FORMAT(////,30X,50(1H-),/,33X,'S T R U C T U R A L   ', 'I N F O R M A T I O N',/,30X,50(1H-),//)
805 FORMAT(3X,'SUBSTANCE',20X,'= ',A80,/)
810 FORMAT(3X,'LATTICE',22X,'= ',A4)
820 FORMAT(3X,'LATTICE CONSTANTS ARE',8X,'= ',3F12.7)
830 FORMAT(3X,'SYMMETRY ATTITUDE IS',9X,'= ',A4)
840 FORMAT(3X,'NUMBER OF ATOMS IN UNITCELL  = ',I3)
850 FORMAT(3X,'MODE OF CALCULATION IS',7X,'= ',A4)
860 FORMAT(3X,'SELFCONSISTENT CYCLE-NUMBER  = ',I3,/)
1000 FORMAT(A80)
1010 FORMAT(A4,23X,I3,/,13X,A4)
1020 FORMAT(6F10.5)
1030 FORMAT(4X,I4,4X,F10.7,3X,F10.7,3X,F10.7,/,15X,I2,17X,I2)          
1031 FORMAT(4X,I4,4X,F10.7,3X,F10.7,3X,F10.7)                          
1040 FORMAT(///,3X,'ERROR IN LAPW2 : MULT(JATOM)=0 ...',/, 20X,'JATOM=',I3,3X,'INDEX=',I3,3X,'MULT=',I3)           
1050 FORMAT(A10,5X,I5,5X,F10.9,5X,F10.5,5X,F10.5)                      
2000 FORMAT(//,3X,'=====>>> CPU TIME SUMMARY',/)                       
2010 FORMAT(12X,'TOTAL          : ',F8.1,5X,'... ',F4.0,' PERCENT')       
2020 FORMAT(12X,'Interst. MatEl : ',F8.1,5X,'... ',F4.0,' PERCENT')       
9000 FORMAT(' can''t open definition file ',A40)
9010 FORMAT(' can''t open unit: ',I2)
9020 FORMAT(' filename: ',A50)
9030 FORMAT(' status: ',A,'  form: ',A)
9040 FORMAT('Error reading file: ',A47)
9044 FORMAT(3X,'Complex momentum matrix elements: (Re_x,Im_x),(Re_y,Im_y),(Re_z,Im_z)')
  !ole ##### Start #####
9050 FORMAT('error : bad number of choices in <case>.inop')
  !ole ##### End #####
END PROGRAM optic


