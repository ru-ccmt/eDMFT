!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! The code reads two W2K-charge densities      !!
!! and writes the average of the two densities  !!
!! The argumens is                              !!
!!        case                                  !!
!! The input files are                          !!
!!        case.clmval, case.clmvaldn            !!
!! The output files is                          !!
!!        case.clmval_aver                      !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE struct
  LOGICAL                  :: rel
  REAL*8                   :: AA,BB,CC,VOL,pia(3),alpha(3)
  REAL*8,ALLOCATABLE       :: R0(:),DX(:),RMT(:),zz(:),rotloc(:,:,:),v(:)
  REAL*8,ALLOCATABLE       :: tau(:,:)
  REAL*8,POINTER           :: pos(:,:)
  CHARACTER*4              :: lattic,irel,cform
  CHARACTER*80             :: title
  CHARACTER*10,ALLOCATABLE :: aname(:)
  INTEGER                  :: nat,iord
  INTEGER,ALLOCATABLE      :: mult(:),jrj(:),iatnr(:),isplit(:)
  INTEGER,ALLOCATABLE      :: iz(:,:,:),inum(:)
  REAL*8,ALLOCATABLE       :: rotij(:,:,:),tauij(:,:)

CONTAINS
  
  SUBROUTINE dealc_struct
    DEALLOCATE(aname, mult, jrj, r0, dx, rmt)
    deallocate(zz, rotloc, iatnr, isplit, v)
    DEALLOCATE (pos)
    DEALLOCATE(iz, tau, inum)
  END SUBROUTINE dealc_struct
  
  SUBROUTINE init_struct(struct_file)
    !USE reallocate
    IMPLICIT NONE
    CHARACTER*100, intent(in) :: struct_file
    INTEGER                   :: ios
    REAL*8                    :: test,ninety
    INTEGER                   :: index,i,j,j1,j2,m,jatom

    test=1.D-5
    ninety=90.0D0

    open (20, file=struct_file, status='old')
    read (20,1000) title
    read (20,1010)  lattic,nat,cform,irel
    REL=.TRUE.                                     
    IF(IREL.EQ.'NREL') REL=.FALSE.                                    
    ALLOCATE(aname(nat),mult(nat),jrj(nat),r0(nat),dx(nat),rmt(nat))
    allocate(zz(nat),rotloc(3,3,nat),iatnr(nat),isplit(nat),v(nat))
    v=0.0d0
    ALLOCATE (pos(3,48*nat))
    read (20,1020) aa,bb,cc,alpha(1),alpha(2),alpha(3)
    IF(ABS(ALPHA(1)).LT.test) ALPHA(1)=ninety
    IF(ABS(ALPHA(2)).LT.test) ALPHA(2)=ninety
    IF(ABS(ALPHA(3)).LT.test) ALPHA(3)=ninety
    INDEX=0
    
    DO jatom=1,NAT
       INDEX=INDEX+1
       READ(20,1030,iostat=ios) iatnr(jatom),( pos(j,index),j=1,3 ), mult(jatom),isplit(jatom) 
       IF(ios /= 0) THEN
          WRITE(*,*) iatnr(jatom),( pos(j,index),j=1,3 ), mult(jatom),isplit(jatom) 
          WRITE(*,*) 'ERROR IN STRUCT FILE READ'
          STOP
       ENDIF
       IF (mult(jatom) .EQ. 0) THEN
          WRITE (*,6000) jatom, index, mult(jatom)
          STOP
       ENDIF
       DO m=1,mult(jatom)-1                                     
          index=index+1                                            
          READ(20,1031) iatnr(jatom),( pos(j,index),j=1,3)   ! pos -- position inside the unit cell read from case.struct
       ENDDO
       READ(20,1050) aname(jatom),jrj(jatom),r0(jatom),rmt(jatom),zz(jatom) ! zz-nuclear charge, jrj--number of radial data points
       dx(jatom)=LOG(rmt(jatom)/r0(jatom)) / (jrj(jatom)-1)           
       rmt(jatom)=r0(jatom)*EXP( dx(jatom)*(jrj(jatom)-1) )           
       READ(20,1051) ((rotloc(i,j,jatom),i=1,3),j=1,3)                
    ENDDO
    !CALL doreallocate(pos, 3, index)
    READ(20,1151) iord
    ALLOCATE(iz(3,3,iord),tau(3,iord),inum(iord))
    DO j=1,iord  ! iz(:,:,iord) - all symmetry transformations
                 ! tau(:,iord)  - translations
       READ(20,1101) ( (iz(j1,j2,j),j1=1,3),tau(j2,j),j2=1,3 ),inum(j)
    ENDDO

1000 FORMAT(A80)                                                       
1010 FORMAT(A4,23X,I3,1x,a4,/,13X,A4,18X,A4)                                 
1020 FORMAT(6F10.7,10X,F10.7)                                          
1030 FORMAT(4X,I4,4X,F10.7,3X,F10.7,3X,F10.7,/,15X,I2,17X,I2)          
1031 FORMAT(4X,I4,4X,F10.7,3X,F10.7,3X,F10.7)                          
1050 FORMAT(A10,5X,I5,5X,F10.9,5X,F10.5,5X,F10.5)                      
1051 FORMAT(20X,3F11.8)
1101 FORMAT(3(3I2,F11.8/),I8)
1151 FORMAT(I4)
6000 FORMAT(///,3X,'ERROR IN LAPW0 : MULT(JATOM)=0 ...', &
          /, 20X,'JATOM=',I3,3X,'INDEX=',I3,3X,'MULT=',I3)
  END SUBROUTINE init_struct

END MODULE struct

subroutine ReadCLMVAL(fclmval1, fclmval2, fclmvala, pfc)
  use struct, ONLY: nat, jrj
  IMPLICIT NONE
  CHARACTER*100, intent(in) :: fclmval1, fclmval2, fclmvala
  ! locals
  INTEGER :: fhv1, fhv2, fhv3
  INTEGER :: ISCF, jatom, jatom1, jatom2, iat, IMAX, LMMAX, ILM, Li, Mi, I, J, JX, NWAVE
  CHARACTER*100 :: str1, str2
  COMPLEX*16 :: RHOK1,RHOK2
  INTEGER :: KZZ(3)
  REAL*8  :: pfc
  INTEGER, ALLOCATABLE :: JRI(:)
  REAL*8, ALLOCATABLE :: RHOLM1(:), RHOLM2(:)
  
  fhv1 = 10
  fhv2 = 11
  fhv3 = 12
  
  IMAX = jrj(1)
  do iat=2,nat
     if (jrj(iat).gt.IMAX) IMAX = jrj(iat)
  enddo
  
  ALLOCATE( RHOLM1(IMAX), RHOLM2(IMAX) )
  
  OPEN(fhv1,file=fclmval1, status='old')
  OPEN(fhv2,file=fclmval2, status='old')
  OPEN(fhv3,file=fclmvala, status='unknown')
  
  
  READ(fhv1,'(A43,5X,I3,A)') str1, ISCF, str2
  READ(fhv2,'(A43,5X,I3,A)') str1, ISCF, str2
  READ(fhv1,*)
  READ(fhv2,*)
  READ(fhv1,*)
  READ(fhv2,*)
  
  WRITE(fhv3,'(A43,5X,I3,A)') '    VALENCE CHARGE DENSITY   IN  MT SPHERES', ISCF, ' ITERATION'
  WRITE(fhv3,'(1x,A)') '   NORM OF CLM(R) =          '
  WRITE(fhv3,*)
  
  DO jatom=1,nat
     IMAX=jrj(JATOM)
     
     READ(fhv1,'(3X,A13,I3,5X)') str1, jatom1
     READ(fhv2,'(3X,A13,I3,5X)') str1, jatom2
     if (jatom1.ne.jatom .or. jatom2.ne.jatom) print *, 'ERROR: atom number not correct! Something wrong reading clmval'
     READ(fhv1,'(3X,A12,I3//)') str1, LMMAX
     READ(fhv2,'(3X,A12,I3//)') str1, LMMAX
     
     WRITE(fhv3,'(3X,A13,I3,5X)')'ATOMNUMBER   ', jatom
     WRITE(fhv3,'(3X,A12,I3//)') 'NUMBER OF LM', LMMAX
     
     DO ILM=1,LMMAX
        READ(fhv1,'(3X,A12,I3,3X,A2,I2/)') str1, Li, str2, Mi
        READ(fhv2,'(3X,A12,I3,3X,A2,I2/)') str1, Li, str2, Mi
        READ(fhv1,'(3X,4E19.12)') ( RHOLM1(I), I=1,IMAX )
        READ(fhv2,'(3X,4E19.12)') ( RHOLM2(I), I=1,IMAX )
        READ(fhv1,'(/)')
        READ(fhv2,'(/)')
        
        WRITE(fhv3,'(3X,A12,I3,3X,A2,I2/)') 'CLM(R) FOR L', Li, 'M=', Mi
        WRITE(fhv3,'(3X,4E19.12)') ( 0.5*pfc*(RHOLM1(I)+RHOLM2(I)), I=1,IMAX )
        WRITE(fhv3,'(/)')
     end DO
     READ(fhv1,'(///)')
     READ(fhv2,'(///)')
     WRITE(fhv3,'(///)')
  ENDDO
  
  READ(fhv1,*)
  READ(fhv2,*)
  READ(fhv1,*)
  READ(fhv2,*)
  READ(fhv1,'(9X,I10,1X,A12)')  NWAVE, str1
  READ(fhv2,'(9X,I10,1X,A12)')  NWAVE, str1

  WRITE(fhv3,'(1x,A)') '   VALENCE CHARGE DENSITY IN INTERSTITIAL '
  WRITE(fhv3,*)
  WRITE(fhv3,'(9X,I10,1X,A12)')  NWAVE, 'NUMBER OF PW'

  DO J=1,NWAVE
     READ(fhv1,'(3X,3I5,2E19.12)')  (KZZ(JX),JX=1,3),RHOK1
     READ(fhv2,'(3X,3I5,2E19.12)')  (KZZ(JX),JX=1,3),RHOK2
     WRITE(fhv3, '(3X,3I5,2E19.12)')   (KZZ(JX),JX=1,3),0.5*pfc*(RHOK1+RHOK2)
  ENDDO
    
  close(fhv1)
  close(fhv2)
  close(fhv3)
  
  DEALLOCATE( RHOLM1, RHOLM2 )
  
end subroutine ReadCLMVAL


program CombineUpDn
  use struct
  IMPLICIT NONE
  !DEFINE BUFFER HOLDS THE COMMAND LINE ARGUMENT
  CHARACTER*100 :: case
  !
  CHARACTER*100 :: struct_file, fclmval1, fclmval2, fclmvala, cpfc
  INTEGER :: iargc
  REAL*8  :: pfc
  !GET THE PARAMETERS FROM THE COMMAND LINE ARGUMENT
  
  if (iargc().LT.1) then
     print *, 'ERROR: Give case as the argument!'
     CALL EXIT(1)
  endif
  
  CALL GETARG(1,case)
  
  pfc=1.0d0
  if (iargc().GT.1) then
     CALL GETARG(2,cpfc)
     READ(cpfc,*) pfc
  endif
  
  struct_file = TRIM(case)//'.struct'
  fclmval1 = TRIM(case)//'.clmval'
  fclmval2 = TRIM(case)//'.clmvaldn'
  fclmvala = TRIM(case)//'.clmval_aver'
  
  !print *, struct_file
  !print *, fclmval1
  !print *, fclmval2
  !print *, fclmvala
  print *, 'pfc=', pfc
  
  CALL init_struct(struct_file)
  CALL ReadCLMVAL(fclmval1, fclmval2, fclmvala, pfc)
  CALL dealc_struct()
end program CombineUpDn
