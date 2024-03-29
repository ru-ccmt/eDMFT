! @Copyright 2007 Kristjan Haule
! 

MODULE structure
  LOGICAL                  :: rel
  REAL*8                   :: aa,bb,cc,alpha(3),pia(3),vol
  REAL*8,ALLOCATABLE       :: r0(:),dx(:)
  REAL*8, allocatable      :: rmt(:), ZZ(:), rotloc(:,:,:), v(:)
  REAL*8,ALLOCATABLE       :: tau(:,:)
  REAL*8,ALLOCATABLE       :: pos(:,:)
  CHARACTER*4              :: lattic, irel, cform
  CHARACTER*80             :: title
  CHARACTER*10,ALLOCATABLE :: aname(:)
  INTEGER,ALLOCATABLE      :: mult(:),jri(:),iatnr(:),isplit(:)
  INTEGER,ALLOCATABLE      :: iz(:,:,:),inum(:)
  REAL*8,ALLOCATABLE       :: rotij(:,:,:),tauij(:,:)
  INTEGER                  :: nat,iord, ndif 
  LOGICAL                  :: ortho
  REAL*8                   :: BR1(3,3), BR2(3,3)
  REAL*8                   :: rot_spin_quantization(3,3)
CONTAINS
  
  SUBROUTINE ReadStructure(fh_str) ! fh_str should be 20
    use error, ONLY: outerr
    IMPLICIT NONE
    INTEGER, intent(in)      :: fh_str
    !
    REAL*8,PARAMETER          :: test = 1.D-12
    REAL*8, allocatable       :: tpos(:,:)
    INTEGER :: index, jatom, i, j, m, ios, j1, j2
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
       READ(fh_str,1030,iostat=ios) iatnr(jatom),( tpos(j,index),j=1,3), mult(jatom),isplit(jatom)
       if ( ios /= 0 .OR. mult(jatom) .EQ. 0 .OR. mult(jatom).gt.(48*nat)) THEN
          !...illegal number of equivalent atoms
          WRITE (6,6000) jatom, index, mult(jatom)
          WRITE(6,*) 'ERROR IN STRUCT FILE READ'
          CALL OUTERR('structure.f90','MULT .EQ. 0')
          STOP 'DMFT2 - Error. Check file dmft1.error'
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
    ndif = sum(mult)
    allocate( pos(3,ndif) )
    pos(:,:) = tpos(:,:ndif)
    deallocate( tpos )

    READ(fh_str,1151) iord
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
!1040 FORMAT(///,3X,'ERROR IN DMFT2 : MULT(JATOM)=0 ...',/, 20X,'JATOM=',I3,3X,'INDEX=',I3,3X,'MULT=',I3)           
1050 FORMAT(A10,5X,I5,5X,F10.9,5X,F10.5,5X,F10.5)
1051 FORMAT(20X,3F10.8)
1101 FORMAT(3(3I2,F10.8/),I8)
1151 FORMAT(I4)
6000 FORMAT(///,3X,'ERROR IN DMFT2 : MULT(JATOM)=0 ...',/, 20X,'JATOM=',I3,3X,'INDEX=',I3,3X,'MULT=',I3)
  END SUBROUTINE ReadStructure


  ! mprint = myrank.eq.master
  SUBROUTINE LATGEN(mprint)                                            
    !     LATGEN GENERATES TWO BRAVAIS MATRICES, DEFINES THE VOLUME OF THE UNIT CELL.
    !     BR1(3,3)  : TRANSFORMS INTEGER RECIPROCAL LATTICE VECTORS AS      
    !                 GIVEN IN THE VECTORLIST OF LAPW1  ( GENERATED IN      
    !                 COORS, TRANSFORMED IN BASISO, AND WRITTEN OUT IN      
    !                 WFTAPE) INTO CARTESIAN SYSTEM                         
    !     BR2(3,3) :  TRANSFORMS A RECIPROCAL LATTICE VECTOR OF A SPE-      
    !                 CIAL COORDINATE SYSTEM ( IN UNITS OF 2 PI / A )       
    !                 TO CARTESIAN SYSTEM                                   
    !
    use error, ONLY: outerr
    IMPLICIT NONE
    LOGICAL, intent(in) :: mprint
    ! local variables
    REAL*8  :: pi, cosab, cosac, cosbc, sinab, sinbc, sqrt3, wurzel, rvfac
    INTEGER :: i, j
    !---------------------------------------------------------------------  
    !                      
    PI=ACOS(-1.0D0)                                                   
    SQRT3=SQRT(3.D0)
    ALPHA(1)=ALPHA(1)*PI/180.0D0                                             
    ALPHA(2)=ALPHA(2)*PI/180.0D0                                             
    ALPHA(3)=ALPHA(3)*PI/180.0D0                                             
    PIA(1)=2.D0*PI/AA                                                 
    PIA(2)=2.D0*PI/BB                                                 
    PIA(3)=2.D0*PI/CC                                                 

    SELECT CASE (LATTIC(1:1))
    CASE ('H')  !.....HEXAGONAL LATTICE   
       BR1(1,1)=2.D0/SQRT3*PIA(1)                                        
       BR1(1,2)=1.D0/SQRT3*PIA(1)                                        
       BR1(1,3)=0.0D0                                                    
       BR1(2,1)=0.0D0                                                    
       BR1(2,2)=PIA(2)                                                   
       BR1(2,3)=0.0D0                                                    
       BR1(3,1)=0.0D0                                                    
       BR1(3,2)=0.0D0                                                    
       BR1(3,3)=PIA(3)                                                   
       !                                                                       
       BR2(1,1)=2.D0/SQRT3*PIA(1)                                        
       BR2(1,2)=1.D0/SQRT3*PIA(1)                                        
       BR2(1,3)=0.0D0                                                    
       BR2(2,1)=0.0D0                                                    
       BR2(2,2)=PIA(2)                                                   
       BR2(2,3)=0.0D0                                                    
       BR2(3,1)=0.0D0                                                    
       BR2(3,2)=0.0D0                                                    
       BR2(3,3)=PIA(3)                                                   
       !                                                                       
       RVFAC=2.D0/SQRT(3.D0)                                             
       ORTHO=.FALSE.                                             
    CASE ('S', 'P')    !.....PRIMITIVE LATTICE    
       SINBC=SIN(ALPHA(1))
       COSAB=COS(ALPHA(3))
       COSAC=COS(ALPHA(2))
       COSBC=COS(ALPHA(1))
       WURZEL=SQRT(SINBC**2-COSAC**2-COSAB**2+2*COSBC*COSAC*COSAB)
       BR2(1,1)= SINBC/WURZEL*PIA(1)
       BR2(1,2)= (-COSAB+COSBC*COSAC)/(SINBC*WURZEL)*PIA(2)
       BR2(1,3)= (COSBC*COSAB-COSAC)/(SINBC*WURZEL)*PIA(3)
       BR2(2,1)= 0.0
       BR2(2,2)= PIA(2)/SINBC
       BR2(2,3)= -PIA(3)*COSBC/SINBC
       BR2(3,1)= 0.0
       BR2(3,2)= 0.0
       BR2(3,3)= PIA(3)
       !
       BR1(1,1)= SINBC/WURZEL*PIA(1)
       BR1(1,2)= (-COSAB+COSBC*COSAC)/(SINBC*WURZEL)*PIA(2)
       BR1(1,3)= (COSBC*COSAB-COSAC)/(SINBC*WURZEL)*PIA(3)
       BR1(2,1)= 0.0
       BR1(2,2)= PIA(2)/SINBC
       BR1(2,3)= -PIA(3)*COSBC/SINBC
       BR1(3,1)= 0.0
       BR1(3,2)= 0.0
       BR1(3,3)= PIA(3)
       !
       RVFAC= 1.d0/WURZEL
       ORTHO=.TRUE.
       if(abs(alpha(1)-pi/2.d0).gt.0.0001) ortho=.false.
       if(abs(alpha(2)-pi/2.d0).gt.0.0001) ortho=.false.
       if(abs(alpha(3)-pi/2.d0).gt.0.0001) ortho=.false.
    CASE ('F')  !.....FC LATTICE
       BR1(1,1)=PIA(1)                                                   
       BR1(1,2)=0.0D0                                                    
       BR1(1,3)=0.0D0                                                    
       BR1(2,1)=0.0D0                                                    
       BR1(2,2)=PIA(2)                                                   
       BR1(2,3)=0.0D0                                                    
       BR1(3,2)=0.0D0                                                    
       BR1(3,1)=0.0D0                                                    
       BR1(3,3)=PIA(3)                                                   
       !     definitions according to column, rows convention for BR2
       BR2(1,1)=-PIA(1)                                                  
       BR2(1,2)= PIA(1)                                                  
       BR2(1,3)= PIA(1)                                                  
       BR2(2,1)= PIA(2)                                                  
       BR2(2,2)=-PIA(2)                                                  
       BR2(2,3)= PIA(2)                                                  
       BR2(3,1)= PIA(3)                                                  
       BR2(3,2)= PIA(3)                                                  
       BR2(3,3)=-PIA(3)                                                  
       !                                                                       
       RVFAC=4.D0                                                        
       ORTHO=.TRUE.                                             
    CASE ('B')  !.....BC LATTICE                                                        
       BR1(1,1)=PIA(1)                                                   
       BR1(1,2)=0.0D0                                                    
       BR1(1,3)=0.0D0                                                    
       BR1(2,1)=0.0D0                                                    
       BR1(2,2)=PIA(2)                                                   
       BR1(2,3)=0.0D0                                                    
       BR1(3,1)=0.0D0                                                    
       BR1(3,2)=0.0D0                                                    
       BR1(3,3)=PIA(3)                                                   
       !                                                                       
       BR2(1,1)= 0.0D0                                                    
       BR2(1,2)= PIA(1)                                                  
       BR2(1,3)= PIA(1)                                                  
       BR2(2,1)= PIA(2)                                                  
       BR2(2,2)= 0.0D0                                                    
       BR2(2,3)= PIA(2)                                                  
       BR2(3,1)= PIA(3)                                                  
       BR2(3,2)= PIA(3)                                                  
       BR2(3,3)= 0.0D0
       !
       RVFAC=2.D0
       ORTHO=.TRUE.                                             
    CASE ('R')   !.....RHOMBOHEDRAL LATTICE
       BR1(1,1)=1.D0/SQRT(3.D0)*PIA(1)                                          
       BR1(1,2)=1.D0/SQRT(3.D0)*PIA(1)                                          
       BR1(1,3)=-2.d0/sqrt(3.d0)*PIA(1)                                         
       BR1(2,1)=-1.0d0*PIA(2)                                                  
       BR1(2,2)=1.0d0*PIA(2)                                                    
       BR1(2,3)=0.0d0*PIA(2)                                                    
       BR1(3,1)=1.0d0*PIA(3)                                                    
       BR1(3,2)=1.0d0*PIA(3)                                                    
       BR1(3,3)=1.0d0*PIA(3)                                                    
       !
       BR2(1,1)=1.D0/SQRT(3.D0)*PIA(1)                                          
       BR2(1,2)=1.D0/SQRT(3.D0)*PIA(1)                                          
       BR2(1,3)=-2.d0/sqrt(3.d0)*PIA(1)                                         
       BR2(2,1)=-1.0d0*PIA(2)                                                   
       BR2(2,2)=1.0d0*PIA(2)                                                    
       BR2(2,3)=0.0d0*PIA(2)                                                    
       BR2(3,1)=1.0d0*PIA(3)                                                    
       BR2(3,2)=1.0d0*PIA(3)                                                    
       BR2(3,3)=1.0d0*PIA(3)                                                    
       RVFAC=6.D0/SQRT(3.D0)
       ORTHO=.FALSE.                                             
    CASE ('C')
       IF(LATTIC(2:3).EQ.'XY') then !.....CXY ORTHOROMBIC CASE
          BR1(1,1)=PIA(1)                                                   
          BR1(1,2)=0.0D0                                                    
          BR1(1,3)=0.0D0                                                    
          BR1(2,1)=0.0D0                                                    
          BR1(2,2)=PIA(2)                                                   
          BR1(2,3)=0.0D0                                                    
          BR1(3,1)=0.0D0                                                    
          BR1(3,2)=0.0D0                                                    
          BR1(3,3)=PIA(3)                                                   
          !                                                                       
          BR2(1,1)= PIA(1)                                                    
          BR2(1,2)= PIA(1)                                                  
          BR2(1,3)= 0.0D0                                                  
          BR2(2,1)=-PIA(2)                                                  
          BR2(2,2)= PIA(2)                                                    
          BR2(2,3)= 0.0D0                                                 
          BR2(3,1)= 0.0D0                                                  
          BR2(3,2)= 0.0D0                                                 
          BR2(3,3)= PIA(3)                                                    
          !                                                                       
          RVFAC=2.D0                                                        
          ORTHO=.TRUE.
       ELSE IF(LATTIC(2:3).EQ.'YZ') then !.....CYZ ORTHOROMBIC
          BR1(1,1)=PIA(1)                                                   
          BR1(1,2)=0.0D0                                                    
          BR1(1,3)=0.0D0                                                    
          BR1(2,1)=0.0D0                                                    
          BR1(2,2)=PIA(2)                                                   
          BR1(2,3)=0.0D0                                                    
          BR1(3,1)=0.0D0                                                    
          BR1(3,2)=0.0D0                                                    
          BR1(3,3)=PIA(3)                                                   
          !                                                                       
          BR2(1,1)= PIA(1)                                                      
          BR2(1,2)= 0.0                                                   
          BR2(1,3)= 0.0                                                      
          BR2(2,1)= 0.0                                                      
          BR2(2,2)= PIA(2)                                                     
          BR2(2,3)= PIA(2)                                                     
          BR2(3,1)= 0.0                                                     
          BR2(3,2)=-PIA(3)                                                     
          BR2(3,3)= PIA(3)                                                     
          !                                                                       
          RVFAC=2.0                                                         
          ORTHO=.TRUE.                                             
       ELSE IF(LATTIC(2:3).EQ.'XZ') then !.....CXZ ORTHOROMBIC OR MONOCLINIC CASE 
          IF(ABS(ALPHA(3)-PI/2.0D0).LT.0.0001) then ! ORTHOROMBIC
             BR1(1,1)=PIA(1)                                                   
             BR1(1,2)=0.0D0                                                    
             BR1(1,3)=0.0D0                                                    
             BR1(2,1)=0.0D0                                                    
             BR1(2,2)=PIA(2)                                                   
             BR1(2,3)=0.0D0                                                    
             BR1(3,1)=0.0D0                                                    
             BR1(3,2)=0.0D0                                                    
             BR1(3,3)=PIA(3)                                                   
             !                                                                       
             BR2(1,1)= PIA(1)                                                   
             BR2(1,2)= 0.0                                                   
             BR2(1,3)= PIA(1)                                                      
             BR2(2,1)= 0.0                                                      
             BR2(2,2)= PIA(2)                                                     
             BR2(2,3)= 0.0                                                     
             BR2(3,1)=-PIA(3)                                                     
             BR2(3,2)= 0.0                                                     
             BR2(3,3)= PIA(3)                                                     
             !                                                                       
             RVFAC=2.0                                                         
             ORTHO=.TRUE.                                             
          ELSE                     !.....CXZ MONOCLINIC CASE 
             !         write(*,*) '  gamma not equal 90'
             SINAB=SIN(ALPHA(3))
             COSAB=COS(ALPHA(3))
             !                                                                       
             BR1(1,1)= PIA(1)/SINAB 
             BR1(1,2)= -PIA(2)*COSAB/SINAB
             BR1(1,3)= 0.0                                                   
             BR1(2,1)= 0.0                                                      
             BR1(2,2)= PIA(2)                                                     
             BR1(2,3)= 0.0                                                     
             BR1(3,1)= 0.0                                                     
             BR1(3,2)= 0.0                                                     
             BR1(3,3)= PIA(3)                                                     
             !                                                                       
             BR2(1,1)= PIA(1)/SINAB 
             BR2(1,2)= -PIA(2)*COSAB/SINAB
             BR2(1,3)= PIA(1)/SINAB 
             BR2(2,1)= 0.0                                                      
             BR2(2,2)= PIA(2)                                                     
             BR2(2,3)= 0.0                                                     
             BR2(3,1)=-PIA(3)                                                     
             BR2(3,2)= 0.0                                                     
             BR2(3,3)= PIA(3)                                                     
             !                                                                       
             RVFAC=2.0/SINAB                                                   
             ORTHO=.FALSE.                                             
          ENDIF
       ELSE
          GOTO 900 !        Error: wrong lattice, stop execution
       ENDIF
    CASE DEFAULT
       !        Error: wrong lattice, stop execution
       GOTO 900
    END SELECT

    if (mprint) then
       write(6,*)' BR1,  BR2'
       do i=1,3
          write(6,654)(br1(i,j),j=1,3),(br2(i,j),j=1,3)
       enddo
    endif
    VOL=AA*BB*CC/RVFAC                                                
    RETURN
    !        Error messages
900 CALL OUTERR('LATGEN','wrong lattice.')
    STOP 'LATGEN - Error'
    !        End of 'LATGEN'
654 format(3f10.5,3x,3f10.5)
  END SUBROUTINE LATGEN
  
  SUBROUTINE ROTDEF(mprint)
    use error, ONLY: outerr
    !                                                                       
    !     ROTDEF GENERATES THE ROTATION-MATRICES ROTIJ(3,3,ntot) TAUIJ(3,ntot) FOR      
    !     NONSYMMORPHIC STRUCTURES. THESE MATRICES PERFORM ROTATIONS        
    !     FROM THE GENERAL COORDINATION-SYSTEM INTO LOCAL SYSTEMS WITH      
    !     SPECIFIED POINTSYMMETRY.                                          
    !     THE MATRICES  ROTIJ(3,3,INDEX)  TRANSFORM THE POSITION OF AN      
    !     ATOM TO IT'S CORRESPONDING POSITION OF AN EQUIVALENT              
    !     ATOM.                                                             
    !                                                                       
    IMPLICIT NONE
    LOGICAL, intent(in) :: mprint
    !
    CHARACTER*67  :: ERRMSG
    REAL*8        :: toler, toler2, ONE, x(3), x1(3)
    INTEGER       :: INDEX, NCOUNT, JATOM, M, INDEX1, I, J, I1
    LOGICAL       :: FOUND
    !
    ! Additions for averaging
    integer         :: jtimes,nthis
    real*8          :: X00, Y00, Z00, X0, Y0, Z0,XA,YA,ZA
    real*8          :: XB, YB, ZB, TT
    !
    toler = 1.D-7
    toler2=1.5d0*toler            
    ONE = 1.D0
    !
    ALLOCATE (rotij(3,3,ndif),tauij(3,ndif))
    !ndif = sum(mult) ! all atoms
    !
    !---------------------------------------------------------------------  
    !  Patch up symmetry issues with positions in limited precision
    DO jtimes=1,2
       DO JATOM=1,ndif
          NTHIS=0
          X00=POS(1,JATOM)
          Y00=POS(2,JATOM)
          Z00=POS(3,JATOM)
          X0=0
          Y0=0
          Z0=0                                           
          DO I=1,IORD                                              
             XA=TAU(1,I)+1.D0
             YA=TAU(2,I)+1.D0
             ZA=TAU(3,I)+1.D0                                                      
             DO J=1,3 
                XA=XA+dble(IZ(J,1,I))*POS(J,JATOM)                               
                YA=YA+dble(IZ(J,2,I))*POS(J,JATOM) 
                ZA=ZA+dble(IZ(J,3,I))*POS(J,JATOM)
             ENDDO
             XA = mod(XA,1.D0)
             YA = mod(YA,1.D0)
             ZA = mod(ZA,1.D0)
             if(xa.gt.0.999999)XA=XA-1.D0
             if(ya.gt.0.999999)YA=YA-1.D0
             if(Za.gt.0.999999)ZA=ZA-1.D0
             XB=ABS(XA-X00)*aa                                      
             YB=ABS(YA-Y00)*bb                 
             ZB=ABS(ZA-Z00)*cc                                      

             IF(XB.LT.1d-2.AND.YB.LT.1d-2.AND.ZB.LT.1d-2) THEN        
                NTHIS=NTHIS+1                                             
                X0=X0+XA
                Y0=Y0+YA
                Z0=Z0+ZA
             ENDIF
          ENDDO
          if(nthis.gt.1)then
             TT=1.D0/dble(nthis)                         
             POS(1,JATOM)=X0*TT
             POS(2,JATOM)=Y0*TT
             POS(3,JATOM)=Z0*TT
             !       write(88,88)'OLD     ',X00,Y00,Z00
             !       write(88,88)'Average ',X0*TT,Y0*TT,Z0*TT,NTHIS,JATOM
          endif
!88        format(a,3f16.13,2i3)
       ENDDO
    ENDDO
    !---------------------------------------------------------------------  
    INDEX=0                                                           
    NCOUNT=0                                                          
    DO JATOM=1,NAT          ! 20
       INDEX1=INDEX+1                                                 
       DO M=1,MULT(JATOM)   ! 30
          INDEX=INDEX+1                                               
          FOUND = .FALSE.
          DO I=1,IORD       ! 25
             x(1:3)=0.d0
             x=MATMUL(TRANSPOSE(iz(1:3,1:3,i)),pos(1:3,index1))
             x(1:3)=x(1:3)+tau(1:3,i)
             x(1) = MOD(x(1) + 1.d0, ONE)
             x(2) = MOD(x(2) + 1.d0, ONE)
             x(3) = MOD(x(3) + 1.d0, ONE)
             X1(1:3)=MOD(ABS(X(1:3)-POS(1:3,INDEX))+toler,one)-toler
             IF(MAXVAL(ABS(X1)).LT.TOLER2) THEN 
                NCOUNT=NCOUNT+1
                TAUIJ(1:3,INDEX)=TAU(1:3,I)
                ROTIJ(1:3,1:3,INDEX)=IZ(1:3,1:3,I) 
                FOUND = .TRUE.
                EXIT !GOTO 30                                                     
             END IF
             !....check positions for centered lattices
             if(lattic(1:1).eq.'B') then
                x1(1:3)=mod(x1(1:3)+0.5d0+toler,one)-toler
                IF (MAXVAL(ABS(X1)).LT.TOLER2) THEN
                   NCOUNT=NCOUNT+1
                   TAUIJ(1:3,INDEX)=TAU(1:3,I)
                   ROTIJ(1:3,1:3,INDEX)=IZ(1:3,1:3,I) 
                   FOUND = .TRUE.
                   EXIT !GOTO 30                                                     
                END IF
             endif
             if(lattic(1:1).eq.'F'.or.lattic(1:3).eq.'CXY') then
                x1(1:2)=mod(x1(1:2)+0.5d0+toler,one)-toler
                IF(MAXVAL(ABS(X1)).LT.TOLER2) THEN   
                   NCOUNT=NCOUNT+1
                   TAUIJ(1:3,INDEX)=TAU(1:3,I)
                   ROTIJ(1:3,1:3,INDEX)=IZ(1:3,1:3,I) 
                   FOUND = .TRUE.
                   EXIT !GOTO 30                                                     
                END IF
                x1(1:2)=mod(x1(1:2)+0.5d0,one)
             endif
             if(lattic(1:1).eq.'F'.or.lattic(1:3).eq.'CXZ') then
                x1(1)=mod(x1(1)+0.5d0+toler,one)-toler
                x1(3)=mod(x1(3)+0.5d0+toler,one)-toler
                IF(MAXVAL(ABS(X1)).LT.TOLER2) THEN
                   NCOUNT=NCOUNT+1
                   TAUIJ(1:3,INDEX)=TAU(1:3,I)
                   ROTIJ(1:3,1:3,INDEX)=IZ(1:3,1:3,I) 
                   FOUND = .TRUE.
                   EXIT  !GOTO 30                                                     
                END IF
                x1(1)=mod(x1(1)+0.5d0,one)
                x1(3)=mod(x1(3)+0.5d0,one)
             endif
             if(lattic(1:1).eq.'F'.or.lattic(1:3).eq.'CYZ') then
                x1(2:3)=mod(x1(2:3)+0.5d0+toler,one)-toler
                IF(MAXVAL(ABS(X1)).LT.TOLER2) THEN        
                   NCOUNT=NCOUNT+1                                             
                   TAUIJ(1:3,INDEX)=TAU(1:3,I)
                   ROTIJ(1:3,1:3,INDEX)=IZ(1:3,1:3,I) 
                   FOUND = .TRUE.
                   EXIT  !GOTO 30                                                     
                END IF
             end if
          ENDDO

          IF (.NOT.FOUND) THEN
             !
             ! Error: no symmetry operation found
             !
             CALL OUTERR('ROTDEF','no symmetry operation found.')
             WRITE(ERRMSG,9000) jatom, index          
             CALL OUTERR('ROTDEF',ERRMSG)
             WRITE(ERRMSG,9010) (POS(I1,JATOM),I1=1,3) 
             CALL OUTERR('ROTDEF',ERRMSG)
             WRITE(ERRMSG,9020) (POS(I1,INDEX),I1=1,3) 
             CALL OUTERR('ROTDEF',ERRMSG)
             STOP 'ROTDEF - Error'
          ELSE
             if (mprint) WRITE(6,'(A,I3,1x,A,I3,1x,A,I3,1x)') 'Operation', I, 'transforms atom ', index, 'from first atom ', index1
          ENDIF
       ENDDO       !  30  
    ENDDO          !  20
    IF(NCOUNT.NE.INDEX) THEN
       CALL OUTERR('ROTDEF','ROTIJ not defined for all atoms of basis.')
       WRITE (ERRMSG,9030) NCOUNT
       CALL OUTERR('ROTDEF',ERRMSG)
       STOP 'ROTDEF - Error'
    ENDIF

    RETURN                                                            
9000 FORMAT ('for jatom, index',I2,i2)
9010 FORMAT ('atomposition of jatom',3F12.7)
9020 FORMAT ('atomposition of index',3F12.7)
9030 FORMAT ('NCOUNT=',I2)
  END SUBROUTINE ROTDEF

  
  SUBROUTINE WriteInfoStructure(fh_stdout)
    INTEGER, intent(in) :: fh_stdout  ! should be 6
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
    deallocate( iz, tau, inum )
    if (allocated(rotij)) deallocate(rotij)
    if (allocated(tauij)) deallocate(tauij)
  END SUBROUTINE DeallocateStructure
  
END MODULE structure

