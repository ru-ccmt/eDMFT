SUBROUTINE generate_lattice(alpha, BR1, BR2, VOL, ORTHO, AA,BB,CC, LATTIC)
  IMPLICIT NONE
  REAL*8, intent(out)   :: VOL
  REAL*8, intent(out)   :: BR1(3,3)   ! conventional reciprocal
  REAL*8, intent(out)   :: BR2(3,3)   ! BR2(:,i) is primitive reciprocal vector b_i
  LOGICAL, intent(out)  :: ORTHO
  REAL*8, intent(in) :: alpha(3)
  REAL*8, intent(in)    :: AA, BB, CC
  CHARACTER*4, intent(in) :: LATTIC
  ! locals
  REAL*8 :: PIA(3)
  REAL*8 :: PI, SQRT3, RVFAC
  REAL*8 :: SINAB, SINBC, COSAB, COSAC, COSBC, WURZEL
  
  PI=ACOS(-1.0D0)
  SQRT3=SQRT(3.0D+0)
  ! It was already divided!
  !ALPHA(1)=ALPHA(1)*PI/180.0D0
  !ALPHA(2)=ALPHA(2)*PI/180.0D0
  !ALPHA(3)=ALPHA(3)*PI/180.0D0
  PIA(1)=2.d0*PI/AA
  PIA(2)=2.d0*PI/BB
  PIA(3)=2.d0*PI/CC

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
     RVFAC=2.D0/SQRT3
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
     BR1(1,1)=1.D0/SQRT3*PIA(1)                                          
     BR1(1,2)=1.D0/SQRT3*PIA(1)                                          
     BR1(1,3)=-2.d0/SQRT3*PIA(1)                                         
     BR1(2,1)=-1.0d0*PIA(2)                                                  
     BR1(2,2)=1.0d0*PIA(2)                                                    
     BR1(2,3)=0.0d0*PIA(2)                                                    
     BR1(3,1)=1.0d0*PIA(3)                                                    
     BR1(3,2)=1.0d0*PIA(3)                                                    
     BR1(3,3)=1.0d0*PIA(3)                                                    
     !
     BR2(1,1)=1.D0/SQRT3*PIA(1)                                          
     BR2(1,2)=1.D0/SQRT3*PIA(1)                                          
     BR2(1,3)=-2.d0/SQRT3*PIA(1)                                         
     BR2(2,1)=-1.0d0*PIA(2)                                                   
     BR2(2,2)=1.0d0*PIA(2)                                                    
     BR2(2,3)=0.0d0*PIA(2)                                                    
     BR2(3,1)=1.0d0*PIA(3)                                                    
     BR2(3,2)=1.0d0*PIA(3)                                                    
     BR2(3,3)=1.0d0*PIA(3)                                                    
     RVFAC=6.D0/SQRT3
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
        GOTO 900 !  ERROR: stop the program : wrong lattice
     ENDIF
  CASE DEFAULT
     !        
     GOTO 900  !  ERROR: stop the program : wrong lattice
  END SELECT
  
  VOL=AA*BB*CC/RVFAC
  
  !print *, 'PIA=', PIA
  !print *, 'ALPHA=', ALPHA
  !print *, 'lattic=', lattic
  !WRITE(6,*) 'BR1='
  !WRITE(6,'(3F12.6)'), BR1
  !WRITE(6,*) 'BR2='
  !WRITE(6,'(3F12.6)'), BR2
  !print *, 'VOL=', VOL
  RETURN
  
900 CONTINUE
  WRITE(6,*) 'Error in generate_lattice : wrong lattice.'
  FLUSH(6)
  STOP 'LATGEN - Error'
end SUBROUTINE generate_lattice

SUBROUTINE direct_lattice(brnn, NAT,alpha, lattic)
  IMPLICIT NONE
  DOUBLE PRECISION, intent(out) :: BRnn(3,3)
  INTEGER,          intent(in)  :: NAT
  DOUBLE PRECISION, intent(in)  :: ALPHA(3)
  CHARACTER*4,      intent(in)  :: lattic
  ! locals
  DOUBLE PRECISION :: PI, GAMMA, BETA, ALPH, COSG1, GAMMA1, SINAB, COSAB, SQRT3
  !-----------------------------------------------------------------------
  !
  PI=ACOS(-1.0D0)
  SQRT3=SQRT(3.0D+0)

  gamma=alpha(3)
  beta=alpha(2)
  alph=alpha(1)
  
  cosg1=(cos(gamma)-cos(alph)*cos(beta))/sin(alph)/sin(beta)
  gamma1=acos(cosg1)
  
  SELECT CASE (LATTIC(1:1))
  CASE ('H') !.....HEXAGONAL CASE
     BRnn(1,1)=SQRT3/2.0d0
     BRnn(1,2)=-0.5d0
     BRnn(1,3)= 0.0d0
     BRnn(2,1)= 0.0d0
     BRnn(2,2)= 1.0d0
     BRnn(2,3)= 0.0d0
     BRnn(3,1)= 0.0d0
     BRnn(3,2)= 0.0d0
     BRnn(3,3)= 1.d0
  CASE ('S', 'P')  !.....PRIMITIVE LATTICE CASE
     BRnn(1,1)=1.0d0*sin(gamma1)*sin(beta)
     BRnn(1,2)=1.0d0*cos(gamma1)*sin(beta)
     BRnn(1,3)=1.0d0*cos(beta)
     BRnn(2,1)=0.0d0
     BRnn(2,2)=1.0d0*sin(alph)
     BRnn(2,3)=1.0d0*cos(alph)
     BRnn(3,1)=0.0d0
     BRnn(3,2)=0.0d0
     BRnn(3,3)=1.0d0
  CASE ('F') !.....FC CASE (DIRECT LATTICE)
     BRnn(1,1)=0.0d0
     BRnn(1,2)=0.5d0
     BRnn(1,3)=0.5d0
     BRnn(2,1)=0.5d0
     BRnn(2,2)=0.0d0
     BRnn(2,3)=0.5d0
     BRnn(3,1)=0.5d0
     BRnn(3,2)=0.5d0
     BRnn(3,3)=0.0d0
  CASE ('B') !.....BC CASE (DIRECT LATTICE)
     BRnn(1,1)=-0.5d0
     BRnn(1,2)=0.5d0
     BRnn(1,3)=0.5d0
     BRnn(2,1)=0.5d0
     BRnn(2,2)=-0.5d0
     BRnn(2,3)=0.5d0
     BRnn(3,1)=0.5d0
     BRnn(3,2)=0.5d0
     BRnn(3,3)=-0.5d0
  CASE ('R')   !.....RHOMBOHEDRAL LATTICE
     BRnn(1,1)=1/(2.d0*SQRT3)
     BRnn(1,2)=-1/2.d0
     BRnn(1,3)=1/3.d0
     BRnn(2,1)=1/(2.d0*SQRT3)
     BRnn(2,2)=1/2.D0
     BRnn(2,3)=1/3.d0
     BRnn(3,1)=-1/SQRT3
     BRnn(3,2)=0.d0
     BRnn(3,3)=1/3.d0
  CASE ('C')   
     IF(LATTIC(2:3).EQ.'XY') then      !.....CXZ ORTHOROMBIC CASE
        BRnn(1,1)=0.5d0
        BRnn(1,2)=-0.5d0
        BRnn(1,3)=0.0d0
        BRnn(2,1)=0.5d0
        BRnn(2,2)=0.5d0
        BRnn(2,3)=0.0d0
        BRnn(3,1)=0.0d0
        BRnn(3,2)=0.0d0
        BRnn(3,3)=1.0d0
     ELSE IF(LATTIC(2:3).EQ.'YZ') then  !.....CYZ ORTHOROMBIC
        BRnn(1,1)=1.0d0
        BRnn(1,2)=0.0d0
        BRnn(1,3)=0.0d0
        BRnn(2,1)=0.0d0
        BRnn(2,2)=0.5d0
        BRnn(2,3)=0.5d0
        BRnn(3,1)=0.0d0
        BRnn(3,2)=-0.5d0
        BRnn(3,3)=0.5d0
     ELSE IF(LATTIC(2:3).EQ.'XZ') then      !.....CXZ ORTHOROMBIC OR MONOCLINIC
        IF(ABS(ALPHA(3)-PI/2.0D0).LT.0.0001) then !.....CXZ ORTHOROMBIC CASE
           BRnn(1,1)=0.5d0
           BRnn(1,2)=0.0d0
           BRnn(1,3)=-0.5d0
           BRnn(2,1)=0.0d0
           BRnn(2,2)=1.0d0
           BRnn(2,3)=0.0d0
           BRnn(3,1)=0.5d0
           BRnn(3,2)=0.0d0
           BRnn(3,3)=0.5d0
        ELSE    !.....CXZ MONOCLINIC CASE
           ! write(6,*) 'gamma not equal 90'
           SINAB=SIN(ALPHA(3))
           COSAB=COS(ALPHA(3))
           BRNN(1,1)=0.5d0*sinab
           BRNN(1,2)=0.5d0*cosab
           BRNN(1,3)=-0.5d0
           BRNN(2,1)=0.0d0
           BRNN(2,2)=1.0d0
           BRNN(2,3)=0.0d0
           BRNN(3,1)=0.5d0*sinab
           BRNN(3,2)=0.5d0*cosab
           BRNN(3,3)=0.5d0
        ENDIF
     ELSE
        GOTO 900 !  ERROR: stop the program : wrong lattice
     ENDIF
  CASE DEFAULT
     GOTO 900  !  ERROR: stop the program : wrong lattice
  END SELECT
  !      write(6,*) 'Bravais Matrix:'
  !      write(6,999) brnn
  !999  format(3f15.5)
  !
  RETURN
  !        Error messages
900 CONTINUE
  WRITE(6,*) 'Error in DIRECT_LATTICE : wrong lattice.'
  FLUSH(6)
  STOP 'LATGEN - Error'
END SUBROUTINE direct_lattice
