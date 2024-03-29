! @Copyright 2007 Kristjan Haule
SUBROUTINE LATGEN(BR1, BR2, VOL, ORTHO, ALPHA_IN, AA, BB, CC, LATTIC)
  !                                                                       
  !     LATGEN GENERATES TWO BRAVAIS MATRICES, DEFINES THE VOLUME OF      
  !     THE UNIT CELL AND CALLS ROTDEF                                    
  !     BR1(3,3)  : TRANSFORMS INTEGER RECIPROCAL LATTICE VECTORS AS      
  !                 GIVEN IN THE VECTORLIST OF LAPW1  ( GENERATED IN      
  !                 COORS, TRANSFORMED IN BASISO, AND WRITTEN OUT IN      
  !                 WFTAPE) INTO CARTESIAN SYSTEM                         
  !     BR2(3,3) :  TRANSFORMS A RECIPROCAL LATTICE VECTOR OF A SPE-      
  !                 CIAL COORDINATE SYSTEM ( IN UNITS OF 2 PI / A )       
  !                 TO CARTESIAN SYSTEM                                   
  !      
  IMPLICIT NONE
  REAL*8, intent(out)     :: BR1(3,3), BR2(3,3)
  REAL*8, intent(out)     :: VOL
  LOGICAL, intent(out)    :: ORTHO
  REAL*8, intent(in)      :: ALPHA_IN(3), AA, BB, CC
  CHARACTER*4, intent(in) :: LATTIC
  ! locals
  REAL*8 :: Pi, SQRT3, ALPHA(3), PIA(3), RVFAC
  REAL*8 :: SINBC, COSAB, COSAC, COSBC, WURZEL, SINAB
  !
  PI=ACOS(-1.0D0)                                                   
  SQRT3=SQRT(3.D0)
  ! From degrees to radians
  ALPHA(1)=ALPHA_IN(1)*PI/180.0D0
  ALPHA(2)=ALPHA_IN(2)*PI/180.0D0
  ALPHA(3)=ALPHA_IN(3)*PI/180.0D0
  ! 2*pi/a, 2*pi/b, 2*pi/c
  PIA(1)=2.D0*PI/AA
  PIA(2)=2.D0*PI/BB
  PIA(3)=2.D0*PI/CC
  
  RVFAC=1.0
  BR1(:,:) = 0.0D0
  BR2(:,:) = 0.0D0
  IF(LATTIC(1:1).EQ.'H') THEN         !.....HEXAGONAL LATTICE                                                 
     BR1(1,1)=2.D0/SQRT3*PIA(1)
     BR1(1,2)=1.D0/SQRT3*PIA(1)
     BR1(2,2)=PIA(2)                                                   
     BR1(3,3)=PIA(3)
     !
     BR2(1,1)=2.D0/SQRT3*PIA(1)                                        
     BR2(1,2)=1.D0/SQRT3*PIA(1)                                        
     BR2(2,2)=PIA(2)                                                   
     BR2(3,3)=PIA(3)                                                   
     !                                                                       
     RVFAC=2.D0/SQRT(3.D0)
     ORTHO=.FALSE.                                             
  ELSEIF (LATTIC(1:1).EQ.'R') THEN  !.....RHOMBOHEDRAL CASE                                                    
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
  ELSEIF (LATTIC(1:1).EQ.'S' .OR. LATTIC(1:1).EQ.'P') THEN  !.....PRIMITIVE LATTICE
     SINBC=SIN(ALPHA(1))
     COSAB=COS(ALPHA(3))
     COSAC=COS(ALPHA(2))
     COSBC=COS(ALPHA(1))
     WURZEL=SQRT(SINBC**2-COSAC**2-COSAB**2+2*COSBC*COSAC*COSAB)
     BR2(1,1)= SINBC/WURZEL*PIA(1)
     BR2(1,2)= (-COSAB+COSBC*COSAC)/(SINBC*WURZEL)*PIA(2)
     BR2(1,3)= (COSBC*COSAB-COSAC)/(SINBC*WURZEL)*PIA(3)
     BR2(2,2)= PIA(2)/SINBC
     BR2(2,3)= -PIA(3)*COSBC/SINBC
     BR2(3,3)= PIA(3)
     !
     BR1(1,1)= SINBC/WURZEL*PIA(1)
     BR1(1,2)= (-COSAB+COSBC*COSAC)/(SINBC*WURZEL)*PIA(2)
     BR1(1,3)= (COSBC*COSAB-COSAC)/(SINBC*WURZEL)*PIA(3)
     BR1(2,2)= PIA(2)/SINBC
     BR1(2,3)= -PIA(3)*COSBC/SINBC
     BR1(3,3)= PIA(3)
     !
     RVFAC= 1.d0/WURZEL
     ORTHO=.TRUE.
     if(abs(alpha(1)-pi/2.d0).gt.0.0001) ortho=.false.
     if(abs(alpha(2)-pi/2.d0).gt.0.0001) ortho=.false.
     if(abs(alpha(3)-pi/2.d0).gt.0.0001) ortho=.false.
  ELSEIF (LATTIC(1:1).EQ.'F') THEN       !.....FC LATTICE                                                        
     BR1(1,1)=PIA(1)                                                   
     BR1(2,2)=PIA(2)                                                   
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
  ELSEIF (LATTIC(1:1).EQ.'B') THEN    !.....BC LATTICE
     BR1(1,1)=PIA(1)                                                   
     BR1(2,2)=PIA(2)                                                   
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
  ELSEIF (LATTIC(1:3).EQ.'CXZ') THEN           !     CXZ case
     IF(ABS(ALPHA(3)-PI/2.0D0).LT.0.0001) then !.....CXZ ORTHOROMBIC CASE 
        BR1(1,1)=PIA(1)                                                   
        BR1(2,2)=PIA(2)                                                   
        BR1(3,3)=PIA(3)                                                   
        !                                                                       
        BR2(1,1)= PIA(1)                                                   
        BR2(1,3)= PIA(1)                                                      
        BR2(2,2)= PIA(2)                                                     
        BR2(3,1)=-PIA(3)                                                     
        BR2(3,3)= PIA(3)                                                     
        !                                                                       
        RVFAC=2.D0
        ORTHO=.TRUE.                                             
     ELSE             !.....CXZ MONOCLINIC CASE 
        SINAB=SIN(ALPHA(3))
        COSAB=COS(ALPHA(3))
        !                                                                       
        BR1(1,1)= PIA(1)/SINAB 
        BR1(1,2)= -PIA(2)*COSAB/SINAB
        BR1(2,2)= PIA(2)                                                     
        BR1(3,3)= PIA(3)                                                     
        !                                                                       
        BR2(1,1)= PIA(1)/SINAB 
        BR2(1,2)= -PIA(2)*COSAB/SINAB
        BR2(1,3)= PIA(1)/SINAB 
        BR2(2,2)= PIA(2)                                                     
        BR2(3,1)=-PIA(3)                                                     
        BR2(3,3)= PIA(3)                                                     
        !                                                                       
        RVFAC=2.0/SINAB                                                   
        ORTHO=.FALSE.                                             
     ENDIF
  ELSEIF(LATTIC(1:3).EQ.'CYZ') THEN  !.....CYZ CASE (CYZ LATTICE BUILD UP)
     BR1(1,1)=PIA(1)                                                   
     BR1(2,2)=PIA(2)                                                   
     BR1(3,3)=PIA(3)                                                   
     !                                                                       
     BR2(1,1)= PIA(1)                                                      
     BR2(2,2)= PIA(2)                                                     
     BR2(2,3)= PIA(2)                                                     
     BR2(3,2)=-PIA(3)                                                     
     BR2(3,3)= PIA(3)                                                     
     !                                                                       
     RVFAC=2.0                                                         
     ORTHO=.TRUE.                                             
  ELSEIF (LATTIC(1:1).EQ.'C') THEN    !.....CXY LATTICE
     BR1(1,1)=PIA(1)                                                   
     BR1(2,2)=PIA(2)                                                   
     BR1(3,3)=PIA(3)                                                   
     !                                                                       
     BR2(1,1)= PIA(1)                                                    
     BR2(1,2)= PIA(1)                                                  
     BR2(2,1)=-PIA(2)                                                  
     BR2(2,2)= PIA(2)                                                    
     BR2(3,3)= PIA(3)                                                    
     !                                                                       
     RVFAC=2.D0                                                        
     ORTHO=.TRUE.
  ELSE
     print *, 'ERROR: LATGEN can not find lattice. Problem with structure file'
  ENDIF
  !.....DEFINES VOLUME OF UNIT CELL                                        
  VOL=AA*BB*CC/RVFAC                                                
  RETURN
END SUBROUTINE LATGEN
