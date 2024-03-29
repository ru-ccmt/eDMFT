! @Copyright 2007 Kristjan Haule
SUBROUTINE struct_sizes(nat,nsym,ndif,lattic,AA,BB,CC,alpha,structf)
  IMPLICIT NONE
  CHARACTER*80, intent(in):: structf
  INTEGER, intent(out)    :: nat, nsym, ndif
  CHARACTER*4, intent(out):: lattic
  REAL*8, intent(out)      :: AA,BB,CC,alpha(3)
  !----------- local variables ---------------
  CHARACTER*4  :: irel,cform
  CHARACTER*80 :: title
  CHARACTER*10 :: aname
  REAL*8       :: r0,rmt,zz,rotloc(3,3),pos(3)
  INTEGER      :: iord,ios,mult,jrj,iatnr,isplit,index,i,j,m,jatom

  open(20,FILE=structf,STATUS='OLD')
  read (20,1000) title
  read (20,1010) lattic,nat,cform,irel
  read (20,1020) aa,bb,cc,alpha(1),alpha(2),alpha(3)

  IF(ABS(ALPHA(1)).LT.1.D-5) ALPHA(1)=90.
  IF(ABS(ALPHA(2)).LT.1.D-5) ALPHA(2)=90.
  IF(ABS(ALPHA(3)).LT.1.D-5) ALPHA(3)=90.
  
  INDEX=0
  DO jatom=1,NAT
     INDEX=INDEX+1
     READ(20,1030,iostat=ios) iatnr,( pos(j),j=1,3 ), mult,isplit 
     IF(ios /= 0 ) THEN
        WRITE(6,*) iatnr,( pos(j),j=1,3 ), mult,isplit
        WRITE(6,*) 'ERROR IN STRUCT FILE READ'
        STOP
     ENDIF
     IF (mult .EQ. 0) THEN
        WRITE (6,6000) jatom, index, mult
        STOP
     ENDIF
     DO m=1,mult-1                                     
        index=index+1                                            
        READ(20,1031) iatnr,( pos(j),j=1,3)   ! pos -- position inside the unit cell read from case.struct
     ENDDO
     READ(20,1050) aname,jrj,r0,rmt,zz ! zz-nuclear charge, jrj--number of radial data points
     READ(20,1051) ((rotloc(i,j),i=1,3),j=1,3)
  ENDDO
  ndif=index
  READ(20,1151) iord
  nsym=iord

  CLOSE(20)
1000 FORMAT(A80)                                                       
1010 FORMAT(A4,23X,I3,1x,a4,/,13X,A4,18X,A4)                                 
1020 FORMAT(6F10.7,10X,F10.7)                                          
1030 FORMAT(4X,I4,4X,F10.7,3X,F10.7,3X,F10.7,/,15X,I2,17X,I2)          
1031 FORMAT(4X,I4,4X,F10.7,3X,F10.7,3X,F10.7)                          
1050 FORMAT(A10,5X,I5,5X,F10.9,5X,F10.5,5X,F10.5)                      
1051 FORMAT(20X,3F10.8)
1151 FORMAT(I4)
6000 FORMAT(///,3X,'ERROR IN struct-read : MULT(JATOM)=0 ...',/, 20X,'JATOM=',I3,3X,'INDEX=',I3,3X,'MULT=',I3)
END SUBROUTINE struct_sizes



SUBROUTINE init_struct(lattic,aname,AA,BB,CC,alpha,tau,pos,rel,r0,dx,RMT,zz,rotloc,mult,jrj,&
&iatnr,isplit,iz,inum, structf,nat,nsym,ndif)
  IMPLICIT NONE
  ! input
  CHARACTER*80, intent(in) :: structf
  INTEGER, intent(in)      :: nat, nsym, ndif
  ! output
  CHARACTER*4, intent(out) :: lattic
  CHARACTER*10,intent(out) :: aname(nat)
  REAL*8, intent(out)      :: AA,BB,CC,alpha(3)
  REAL*8, intent(out)      :: tau(3,nsym), pos(3,ndif)
  LOGICAL, intent(out)     :: rel
  REAL*8, intent(out)      :: r0(nat),dx(nat),RMT(nat),zz(nat),rotloc(3,3,nat)
  INTEGER, intent(out)     :: mult(nat),jrj(nat),iatnr(nat),isplit(nat)
  INTEGER, intent(out)     :: iz(3,3,nsym),inum(nsym)
  !----------- local variables ---------------
  CHARACTER*80     :: title
  CHARACTER*4      :: irel,cform
  INTEGER          :: ios,iord
!loop indexs
  INTEGER          :: index,i,j,j1,j2,m,jatom,nato

  open(20,FILE=structf,STATUS='OLD')
  read (20,1000) title
  read (20,1010) lattic,nato,cform,irel
  if (nato.NE.nat) WRITE(6,*) 'ERROR init_struct: nat(1)!=nat(2)'
  REL=.TRUE.                                     
  IF(IREL.EQ.'NREL') REL=.FALSE.                                    
  read (20,1020) aa,bb,cc,alpha(1),alpha(2),alpha(3)
  IF(ABS(ALPHA(1)).LT.1.D-5) ALPHA(1)=90.0D0
  IF(ABS(ALPHA(2)).LT.1.D-5) ALPHA(2)=90.0D0
  IF(ABS(ALPHA(3)).LT.1.D-5) ALPHA(3)=90.0D0
  INDEX=0
  DO jatom=1,nat
     INDEX=INDEX+1
     READ(20,1030,iostat=ios) iatnr(jatom),( pos(j,index),j=1,3 ), mult(jatom),isplit(jatom) 
     IF(ios /= 0 ) THEN
        WRITE(6,*) iatnr(jatom),( pos(j,index),j=1,3 ), mult(jatom),isplit(jatom) 
        WRITE(6,*) 'ERROR IN STRUCT FILE READ'
        STOP
     ENDIF
     IF (mult(jatom) .EQ. 0) THEN
        WRITE (6,6000) jatom, index, mult(jatom)
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
  if (ndif.NE.index) WRITE(6,*) 'ERROR init_struct: ndif(1)!=ndif(2)' 
  READ(20,1151) iord
  if (nsym.NE.iord) WRITE(6,*) 'ERROR init_struct: nsym(1)!=nsym(2)'
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
1051 FORMAT(20X,3F10.8)
1101 FORMAT(3(3I2,F10.8/),I8)
1151 FORMAT(I4)
6000 FORMAT(///,3X,'ERROR IN LAPW0 : MULT(JATOM)=0 ...',/, 20X,'JATOM=',I3,3X,'INDEX=',I3,3X,'MULT=',I3)
END SUBROUTINE init_struct



SUBROUTINE LATGEN(BR1, BR2, vol, ORTHO, lattic, AA, BB, CC, alpha)
  !*******************************************************************
  !*  LATGEN GENERATES TWO BRAVAIS MATRICES, DEFINES THE VOLUME OF   *
  !*  THE UNIT CELL AND CALLS ROTDEF                                 *
  !*  BR1(3,3)  : TRANSFORMS INTEGER RECIPROCAL LATTICE VECTORS AS   *
  !*              GIVEN IN THE VECTORLIST OF LAPW1  ( GENERATED IN   *
  !*              COORS, TRANSFORMED IN BASISO, AND WRITTEN OUT IN   *
  !*              WFTAPE) INTO CARTESIAN SYSTEM                      *
  !*  BR2(3,3) :  TRANSFORMS A RECIPROCAL LATTICE VECTOR OF A SPE-   *
  !*              CIAL COORDINATE SYSTEM ( IN UNITS OF 2 PI / A )    *
  !*              TO CARTESIAN SYSTEM                                *
  !*******************************************************************
  IMPLICIT NONE
  REAL*8, intent(out)    :: BR1(3,3),BR2(3,3), vol
  LOGICAL, intent(out)   :: ORTHO
  CHARACTER*4, intent(in):: lattic
  REAL*8, intent(in)     :: AA,BB,CC,alpha(3)
  !-----------------------------------------
  ! local variables
  REAL*8 :: PI, SQRT3, RVFAC, SINAB, SINBC, COSAB, COSBC, WURZEL, COSAC
  REAL*8 :: PIA(3)
  REAL*8 :: talpha(3)
  !---------------------------------------------------------------------  
  PI=ACOS(-1.0D0)                                                   
  SQRT3=SQRT(3.D0)
  tALPHA(:)=ALPHA(:)*PI/180.0D0                                             
  PIA(1)=2.D0*PI/AA                                                 
  PIA(2)=2.D0*PI/BB                                                 
  PIA(3)=2.D0*PI/CC
  
  IF ( LATTIC(1:1).EQ.'H') THEN !.....HEXAGONAL LATTICE
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
  ELSE IF(LATTIC(1:1).EQ.'R') THEN  !.....RHOMBOHEDRAL CASE                                                    
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
  ELSE IF( LATTIC(1:1).EQ.'S' .OR. LATTIC(1:1).EQ.'P' ) THEN !.....PRIMITIVE LATTICE
     SINBC=SIN(tALPHA(1))
     COSAB=COS(tALPHA(3))
     COSAC=COS(tALPHA(2))
     COSBC=COS(tALPHA(1))
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
     if(abs(talpha(1)-pi/2.d0).gt.0.0001) ortho=.false.
     if(abs(talpha(2)-pi/2.d0).gt.0.0001) ortho=.false.
     if(abs(talpha(3)-pi/2.d0).gt.0.0001) ortho=.false.
     !
  ELSE IF(LATTIC(1:1).EQ.'F') THEN !.....FC LATTICE
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
  ELSE IF(LATTIC(1:1).EQ.'B') THEN  !.....BC LATTICE                                                        
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
     !GOTO 100                                                    
  ELSE IF(LATTIC(1:1).EQ.'C') THEN
     !             
     IF(LATTIC(2:3).EQ.'XZ') THEN
        !.....CXZ CASE (CXZ LATTICE BUILD UP)                                     
        !.....CXZ ORTHOROMBIC CASE 
        IF(ABS(tALPHA(3)-PI/2.0D0).LT.0.0001) then
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
        ELSE
           !.....CXZ MONOCLINIC CASE 
           SINAB=SIN(tALPHA(3))
           COSAB=COS(tALPHA(3))
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
     ELSE IF(LATTIC(2:3).EQ.'YZ') THEN
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
     ELSE
        !.....CXY LATTICE                                                          
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
     ENDIF
  ELSE
     ! Error: wrong lattice, stop execution
     !
     WRITE(6,*) 'LATGEN: wrong lattice.'
     STOP 'LATGEN - Error'
  ENDIF
  
  !.....DEFINE VOLUME OF UNIT CELL                                        
  VOL=AA*BB*CC/RVFAC                                                
  RETURN
  !
654 format(3f10.5,3x,3f10.5)
END SUBROUTINE LATGEN
