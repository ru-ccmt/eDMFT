! @Copyright 2007 Kristjan Haule
! 

MODULE structure
  INTEGER,ALLOCATABLE      :: iatnr(:),mult(:),isplit(:),jri(:)
  REAL*8,ALLOCATABLE       :: pos(:,:)
  REAL*8,ALLOCATABLE       :: rotij(:,:,:),tauij(:,:)
  REAL*8,ALLOCATABLE       :: rmt(:),v(:),rotloc(:,:,:)
  REAL*8,ALLOCATABLE       :: r0(:),dx(:)
  CHARACTER*10,ALLOCATABLE :: aname(:)
  REAL*8, allocatable      :: ZZ(:)
  CHARACTER*4              :: irel, cform
  LOGICAL                  :: ortho
  INTEGER                  :: natm
  REAL*8                   :: aa,bb,cc,alpha(3),pia(3),vol
  REAL*8                   :: BR1(3,3), BR2(3,3)
  CHARACTER                :: title*80,lattic*4
  REAL*8                   :: rot_spin_quantization(3,3)
CONTAINS
  
  SUBROUTINE ReadStructure(fh_str,nat,rel,lxdos)
    IMPLICIT NONE
    INTEGER, intent(in)      :: fh_str ! should be 20
    INTEGER, intent(out)     :: nat
    LOGICAL, intent(out)     :: rel
    INTEGER, intent(inout)   :: lxdos
    !
    REAL*8,PARAMETER          :: test = 1.D-12
    REAL*8, allocatable       :: tpos(:,:)
    INTEGER :: index, jatom, i, j, m
    !
    READ(fh_str,1000) title
    READ(fh_str,1010) lattic,nat,cform,irel
    allocate( rmt(nat),v(nat),iatnr(nat),mult(nat),isplit(nat) )
    allocate( rotloc(3,3,nat) )
    allocate( r0(nat), dx(nat), jri(nat) )
    allocate( aname(nat) )
    allocate( ZZ(nat) )
    allocate ( tpos(3,48*nat) ) ! temporary allocate large array
    
    !.....READ IN LATTICE CONSTANTS                                         
    READ(fh_str,1020) aa,bb,cc,alpha(1),alpha(2),alpha(3)
    if (abs(alpha(1)).LT.test) alpha(1)=90.0d0
    if (abs(alpha(2)).LT.test) alpha(2)=90.0d0
    if (abs(alpha(3)).LT.test) alpha(3)=90.0d0
    rel=.FALSE.
    IF(IREL.EQ.'RELA') rel=.TRUE.
    !IF(IREL.EQ.'NREL') REL=.FALSE.
    !
    !  read crystal-structure (atompositions, symmetry-parameters,muffin-tin radius, ...)
    !  'INDEX' counts all atoms in the unit cell,
    !  'JATOM' counts only the non-equivalent atoms
    index = 0                                                          
    DO jatom = 1,NAT                                               
       index = index+1
       READ(fh_str,1030) iatnr(jatom),( tpos(j,index),j=1,3), mult(jatom),isplit(jatom)
       if(isplit(jatom).eq.99)    lxdos=   3
       if(isplit(jatom).eq.88)    lxdos=   3
       if ( mult(jatom) .EQ. 0 .OR. mult(jatom).gt.(48*nat)) THEN
          !...illegal number of equivalent atoms
          WRITE (6,6000) jatom, index, mult(jatom)
          !INFO = 3
          !        illegal number of equivalent atoms
          CALL OUTERR('structure.f90','MULT .EQ. 0')
          STOP 'DMFT2 - Error. Check file dmft2.error'
       ENDIF
       DO m=1,mult(jatom)-1
          index = index+1                                            
          READ(fh_str,1031) iatnr(jatom),( tpos(j,index),j=1,3)
       ENDDO
       READ(fh_str,1050) aname(jatom), jri(jatom), r0(jatom), rmt(jatom), ZZ(jatom)
       dx(jatom) = log(rmt(jatom)/r0(jatom)) / (jri(jatom)-1)
       rmt(jatom) = r0(jatom)*exp(dx(jatom)*(jri(jatom)-1) )
       READ(fh_str,1051) ((rotloc(i,j,jatom),i=1,3),j=1,3)
    ENDDO
    natm = sum(mult)
    allocate( pos(3,natm) )
    pos(:,:) = tpos(:,:natm)
    deallocate( tpos )
    
1000 FORMAT(A80)
1010 FORMAT(A4,23X,I3,1x,a4,/,13X,A4,18X,A4)                                 
1020 FORMAT(6F10.7,10X,F10.7)                                          
1030 FORMAT(4X,I4,4X,F10.7,3X,F10.7,3X,F10.7,/,15X,I2,17X,I2)          
1031 FORMAT(4X,I4,4X,F10.7,3X,F10.7,3X,F10.7)                          
1040 FORMAT(///,3X,'ERROR IN DMFT2 : MULT(JATOM)=0 ...',/, 20X,'JATOM=',I3,3X,'INDEX=',I3,3X,'MULT=',I3)           
1050 FORMAT(A10,5X,I5,5X,F10.9,5X,F10.5,5X,F10.5)
1051 FORMAT(20X,3F10.8)
6000 FORMAT(///,3X,'ERROR IN DMFT2 : MULT(JATOM)=0 ...',/, 20X,'JATOM=',I3,3X,'INDEX=',I3,3X,'MULT=',I3)
  END SUBROUTINE ReadStructure

  SUBROUTINE WriteInfoStructure(fh_stdout, nat)
    INTEGER, intent(in) :: fh_stdout  ! should be 6
    INTEGER, intent(in) :: nat
    WRITE(fh_stdout,800)                                                      
    WRITE(fh_stdout,805)  title
    WRITE(fh_stdout,810)  lattic
    WRITE(fh_stdout,820)  aa,bb,cc                                            
    WRITE(fh_stdout,840)  nat
    WRITE(fh_stdout,850)  irel
800 FORMAT(////,30X,50(1H-),/,33X,'S T R U C T U R A L   ','I N F O R M A T I O N',/,30X,50(1H-),//)                  
805 FORMAT(3X,'SUBSTANCE',20X,'= ',A80,/)                             
810 FORMAT(3X,'LATTICE',22X,'= ',A4)                                  
820 FORMAT(3X,'LATTICE CONSTANTS ARE',8X,'= ',3F12.7)                 
840 FORMAT(3X,'NUMBER OF ATOMS IN UNITCELL  = ',I3)                   
850 FORMAT(3X,'MODE OF CALCULATION IS',7X,'= ',A4)                    
  END SUBROUTINE WriteInfoStructure
  
  SUBROUTINE DeallocateStructure()
    IMPLICIT NONE
    deallocate( rmt, v, iatnr, mult, isplit )
    deallocate( rotloc )
    deallocate( r0, dx, jri )
    deallocate( aname )
    deallocate( ZZ )
    deallocate( pos )
  END SUBROUTINE DeallocateStructure
  
END MODULE structure

