SUBROUTINE LATGEN(NAT)                                           
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
  use structure, ONLY: alpha, pia, lattic, BR1, BR2, aa, bb, cc, mult, natm, ortho, pos, vol!, rotij, tauij
  USE sym2,      ONLY: !init_sym2, iord, iz, tau
  USE com_mpi,   ONLY: myrank, master
  !IMPLICIT REAL*8 (A-H,O-Z)
  IMPLICIT NONE
  INTEGER, intent(in) :: nat
  ! locals
  REAL*8  :: pi, cosab, cosac, cosbc, sinab, sinbc, sqrt3, wurzel, rvfac
  INTEGER :: i, j, j1, j2
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
  
  if (myrank.EQ.master) then
     write(6,*)' BR1,  BR2'
     do i=1,3
        write(6,654)(br1(i,j),j=1,3),(br2(i,j),j=1,3)
     enddo
  endif
  
  VOL=AA*BB*CC/RVFAC


  
!  !                                                                       
!  !.....DEFINE ROTATION MATRICES IN NONSYMMORPHIC CASE                    
!  !      CALL ROTDEF (NAT)                                            
!  READ(20,'(i4)') IORD                                                  
!  nsym=iord
!  CALL init_sym2(nsym)
!  DO  J=1,IORD                                                     
!     READ(20,11) ( (iz(J1,J2,J),J1=1,3),TAU(J2,J), J2=1,3 )          
!  enddo
!
!  CALL ROTDEF(iz,tau,iord,nat,pos,natm,rotij,tauij,mult,lattic)
!  !.....READ IN SYMMETRY OPERATIONS AND NONPRIMITIVE TRANSLATIONS         
!11 FORMAT(3(3I2,F11.8/))                                              
!  !


  
  RETURN
  !
  !        Error messages
  !
900 CALL OUTERR('LATGEN','wrong lattice.')
  STOP 'LATGEN - Error'
  !
  !        End of 'LATGEN'
  !
654 format(3f10.5,3x,3f10.5)
END SUBROUTINE LATGEN
