SUBROUTINE ATPAR (REL,NAT,JATOM,LATOM,is) 
  use struk
  use potnlc
  use lolog1, only : loor1
  use xrpar
  !adpi12
  !ad   last updates:
  !ad   09/01  cad: implementation of APW+LO
  !ad                           
  !                                                                       
  !     LINEAR EXPANDED ORBITAL APW ROUTINE                               
  !     SET UP THE RADIAL INTEGRAL PARAMETERS                             
  !     SPIN-ORBITLESS RELATIVISTIC EQUATIONS USED                        
  !     L IS 7  (LMAX)                                                    
  !     REQUIRED SUBROUTINES ARE OUTWIN AND RINT13                        
  !     SET UP FOR 4 ATOMS, AND CALL 4 TIMES                              
  !                                                                       
  !
  IMPLICIT REAL*8 (A-H,O-Z)
  !                                                                       
  INCLUDE 'param.inc'
  CHARACTER*4      LATTIC                                           
  CHARACTER*3      MODUS                                            
  CHARACTER*80     TITLE                                            
  LOGICAL          REL,lapw(0:lmax2),loor(0:lomax),lloor(0:lmax2)
  dimension        emist(0:lomax,nloat)
  !                                                                       
  COMMON /ATSPDT/  E(0:LMAX2),P(0:LMAX2),DP(0:LMAX2),PE(0:LMAX2),DPE(0:LMAX2),PEI(0:LMAX2)              
  COMMON /RADFU/   RRAD1(NRAD,0:LMAX2),RADE1(NRAD,0:LMAX2),RRAD2(NRAD,0:LMAX2),RADE2(NRAD,0:LMAX2)                                    
  COMMON /CHAR/   TITLE,LATTIC,MODUS                   
  real*8 VR(NRAD)
  COMMON /UHELP/   A(NRAD),B(NRAD),AP(NRAD),BP(NRAD),AE(NRAD),BE(NRAD)                                         
  common /loabc/   alo(0:lomax,nloat),blo(0:lomax,nloat),clo(0:lomax,nloat),elo(0:lomax,nloat),plo(0:lomax),dplo(0:lomax),pelo(0:lomax),dpelo(0:lomax),peilo(0:lomax),pi12lo(0:lomax),pe12lo(0:lomax),a1lo(nrad,0:lomax),b1lo(nrad,0:lomax)
  common /lolog/   nlo,nlov,nlon,lapw,ilo(0:lomax),loor,lloor
  !---------------------------------------------------------------------  
  !      if(cform.eq.'NEW ') then
  !        assign 2022 to iform1
  !      else
  !        assign 2021 to iform1
  !      end if       
2021 FORMAT(3X,5E13.7)                                                 
2022 FORMAT(3X,4E19.12)  
  IF (LATOM.EQ.1) THEN
     READ(17+is,2032) ISCF
  END IF
  !                                                                       
  !.....READ TOTAL SPHERICAL POTENTIAL V(0,0) OF TAPE18=VSP               
  !     NORM OF V(0,0)=V00(R)*R/SQRT(4.D0*PI)    
  !     write(6,*)'REL',REL                         
  READ(17+is,1980)                                                     
  READ(17+is,2000) IDUMMY                                               
  READ(17+is,2031)                                                   
  READ(17+is,2022) ( VR(J), J=1,JRI(JATOM) )                          
  READ(17+is,2031)                                                     
  READ(17+is,2030)

  DO J=1,JRI(JATOM)                                             
     VR(J)=VR(J)/2.0D0   ! Converting potential from Rydberg to Hartree                                              
  ENDDO
  !   
  nlo=0
  nlov=0
  nlon=0
  !ad
  ! GM nlo #LO on this atom, nlov #LO up til now, nlop #LO left
  !      lapw=.false.
  do l=0,lomax
     ilo(l)=0
  enddo
  !ad
  DO I=1,JATOM
     READ(9+is) E
     READ(9+is) elo
     !      write(*,*) 'atpar_op: atom,E,Elo ',i,E,Elo
     !ad    check   lapw /apw by E(l)
     if(i.eq.jatom) then
        do l=0,lmax2
           lapw(l)=.true.
           if(e(l).gt.150.) then
              e(l)=e(l)-200.d+0
              lapw(l)=.false.
           endif
        enddo
     endif
     !ad    check LO's by ELO value
     do l = 0,lomax
        loor(l)=.false.
        do k=1,nloat
           if (i.eq.jatom) then 
              if (elo(l,k).lt.(995.D+0)) then
                 ilo(l)=k
                 nlo=nlo+((2*l+1))*mult(i)
              endif
           else
              if (elo(l,k).lt.(995.D+0)) nlov=nlov+((2*l+1))*mult(i)
           endif
        enddo
        if((lapw(l).and.ilo(l).eq.1).or.ilo(l).eq.2) then
           loor(l)=.true.
           !adr       write(*,*) 'l,lapw(l),ilo(l)',l,lapw(l),ilo(l)
        endif
     enddo
  enddo
  !ad
  ! LO
  if (xmcd.eq.1) loor1(0:lomax)=loor(0:lomax)
  ! LO
  IF(JATOM.EQ.NAT) GOTO 30                                          
  DO I=JATOM+1,NAT  
     READ(9+is) EMIST                                                    
     READ(9+is) EMIST                                                    
     do l = 0,lomax
        do k=1,nloat
           if (emist(l,k).lt.(995.0D+0)) nlon=nlon+((2*l+1))*mult(i)
        enddo
     enddo
  enddo
  !ad
30 CONTINUE     
  WRITE(6,13)  JATOM,IATNR(JATOM),(POS(I,LATOM),I=1,3),MULT(JATOM)        
  WRITE(6,1060) ANAME(JATOM)                                    
  INDEX=LATOM-1                                                 
  DO M=1,MULT(JATOM) ! 210
     INDEX=INDEX+1
     WRITE(6,1080) M,( POS(JC,INDEX),JC=1,3 )
  ENDDO
  WRITE(6,7) ANAME(JATOM)                                           
  WRITE(6,5) E                                                      
  WRITE(6,14)                                                       
  !                                                                       
  DO l=0,LMAX2    ! 70                                              
     DELE=2.0D-3                                                       
     DELEI=0.25D0/DELE                                                 
     FL=L                                                              
     EI=E(l)/2                                                         
     !                                                                       
     !     CALCULATE ENERGY-DERIVATIVE BY FINITE DIFFERENCE                  
     !     DELE IS THE UPWARD AND DOWNWARD ENERGY SHIFT IN HARTREES          
     !                                                                       
     E1=EI-DELE                                                        
     CALL OUTWIN(REL,VR,RNOT(JATOM),DX(JATOM),JRI(JATOM),E1,FL,UVB,DUVB,NODEL,ZZ(jatom))                                 
     CALL RINT13(REL,A,B,A,B,OVLP,JATOM)                               
     TRX=1.0D0/SQRT(OVLP)                                              
     IMAX=JRI(JATOM)                                                   
     DO M=1,IMAX                                                    
        AE(M)=TRX*A(M)                                                    
        BE(M)=TRX*B(M)                                                    
     ENDDO
     UVB=TRX*UVB                                                       
     DUVB=TRX*DUVB                                                     
     E1=EI+DELE                                                        
     CALL OUTWIN(REL,VR,RNOT(JATOM),DX(JATOM),JRI(JATOM),E1,FL,UVE,DUVE,NODE,ZZ(jatom))                                  
     CALL RINT13(REL,A,B,A,B,OVLP,JATOM)                               
     TRX=1.0/SQRT(OVLP)                                                
     UVE=DELEI*(TRX*UVE-UVB)                                           
     DUVE=DELEI*(TRX*DUVE-DUVB)                                        
     IMAX=JRI(JATOM)                                                   
     DO M=1,IMAX  ! 50
        AE(M)=DELEI*(TRX*A(M)-AE(M))                                      
        BE(M)=DELEI*(TRX*B(M)-BE(M))
     ENDDO        ! 50
     !     CALCULATE FUNCTION AT EI                                          
     CALL OUTWIN(REL,VR(1),RNOT(JATOM),DX(JATOM),JRI(JATOM),EI,FL,UV,DUV,NODES,ZZ(jatom))                                
     CALL RINT13(REL,A,B,A,B,OVLP,JATOM)                               
     TRX=1.0/SQRT(OVLP)                                                
     P(l)=TRX*UV                                                       
     DP(l)=TRX*DUV                                                     
     IMAX=JRI(JATOM) 
     DO M=1,IMAX                                                    
        A(M)=TRX*A(M)                                                     
        B(M)=TRX*B(M)
     ENDDO
     !                                                                       
     !     INSURE ORTHOGONALIZATION                                          
     !                                                                       
     CALL RINT13(REL,A,B,AE,BE,CROSS,JATOM)                            
     TRY=-CROSS                                                        
     IF(TRY.LT.(-0.05)) WRITE(6,9) L,TRY,OVLP                          
     IMAX=JRI(JATOM)                                                   
     DO M=1,IMAX        ! 55
        AE(M)=(AE(M)+TRY*A(M))                                            
        BE(M)=(BE(M)+TRY*B(M)) 
     ENDDO              ! 55
     IMAX=JRI(JATOM)                                                   
     DO I=1,IMAX        ! 80                                            
        RRAD1(I,l)=A(I)                                                   
        RRAD2(I,l)=B(I)                                                   
        RADE1(I,l)=AE(I)                                                  
        RADE2(I,l)=BE(I)                                                 
     ENDDO              ! 80    CONTINUE                                                          
     PE(l)=UVE+TRY*P(l)                                                
     DPE(l)=DUVE+TRY*DP(l)                                             
     CALL RINT13(REL,AE,BE,AE,BE,PEI(l),JATOM)                         
     WRITE(6,8) L,P(l),DP(l),PE(l),DPE(l),PEI(l),NODEL,NODES,NODE                         
  ENDDO
  !
  !   now for LO's
  !                         
  DO l=0,lomax                ! 170                                       
     !ad
     pi12lo(l)=0.d0 
     pe12lo(l)=0.d0      
     !ad
     if (.not.loor(l)) CYCLE  !goto 170
     DELE=2.0D-3                                                       
     DELEI=0.25D0/DELE                                                 
     FL=L                        
     !ad                      
     EI=elo(l,ilo(l))/2     
     !adr    write(*,*) 'atpar_op: EI for LO ',EI                                 
     !ad
     !                                                                       
     !     CALCULATE ENERGY-DERIVATIVE BY FINITE DIFFERENCE                  
     !     DELE IS THE UPWARD AND DOWNWARD ENERGY SHIFT IN HARTREES          
     !                                                                       
     E1=EI-DELE                                                        
     CALL OUTWIN(REL,VR,RNOT(JATOM),DX(JATOM),JRI(JATOM),E1,FL,UVB,DUVB,NODEL,ZZ(jatom))                             
     CALL RINT13(REL,A,B,A,B,OVLP,JATOM)                               
     TRX=1.0D0/SQRT(OVLP)                                              
     IMAX=JRI(JATOM)                                                   
     DO M=1,IMAX      ! 145                                              
        AE(M)=TRX*A(M)                                                    
        BE(M)=TRX*B(M)                                                    
     ENDDO            ! 145
     UVB=TRX*UVB                                                       
     DUVB=TRX*DUVB                                                     
     E1=EI+DELE                                                        
     CALL OUTWIN(REL,VR,RNOT(JATOM),DX(JATOM),JRI(JATOM),E1,FL,UVE,DUVE,NODE,ZZ(jatom))                          
     CALL RINT13(REL,A,B,A,B,OVLP,JATOM)                               
     TRX=1.0/SQRT(OVLP)                                                
     UVE=DELEI*(TRX*UVE-UVB)                                           
     DUVE=DELEI*(TRX*DUVE-DUVB)                                        
     IMAX=JRI(JATOM)                                                   
     DO M=1,IMAX     ! 150
        AE(M)=DELEI*(TRX*A(M)-AE(M))                                      
        BE(M)=DELEI*(TRX*B(M)-BE(M))                                      
     ENDDO           ! 150
     !                                                                       
     !     CALCULATE FUNCTION AT EI                                          
     !                                                                       
     CALL OUTWIN(REL,VR(1),RNOT(JATOM),DX(JATOM),JRI(JATOM),EI,FL,UV,DUV,NODES,ZZ(jatom))                 
     CALL RINT13(REL,A,B,A,B,OVLP,JATOM)                               
     TRX=1.0/SQRT(OVLP)                                                
     Plo(l)=TRX*UV                                                       
     DPlo(l)=TRX*DUV                                                     
     IMAX=JRI(JATOM)                                                   
     DO M=1,IMAX       ! 160                                             
        A(M)=TRX*A(M)                                                     
        B(M)=TRX*B(M)
     ENDDO             ! 160
     !ad
     !ad   maybe not necessary
     !ad
     !                                                                       
     !     INSURE ORTHOGONALIZATION                                          
     !                                                                       
     CALL RINT13(REL,A,B,AE,BE,CROSS,JATOM)                            
     TRY=-CROSS                                                        
     IF(TRY.LT.(-0.05)) WRITE(6,9) L,TRY,OVLP                          
     IMAX=JRI(JATOM)                                                   
     DO M=1,IMAX                ! 155                                    
        AE(M)=(AE(M)+TRY*A(M))                                            
        BE(M)=(BE(M)+TRY*B(M))  
     ENDDO                      ! 155
     Pelo(l)=UVE+TRY*Plo(l)                                                
     DPElo(l)=DUVE+TRY*DPlo(l)                                             
     CALL RINT13(REL,AE,BE,AE,BE,peilo(l),JATOM)         
     !ad
     !ad   not necessary end
     !ad
     CALL RINT13(REL,rrad1(1,l),rrad2(1,l),A,B,pi12lo(l),JATOM)         
     CALL RINT13(REL,rade1(1,l),rade2(1,l),A,B,pe12lo(l),JATOM)         
     DO I=1,IMAX  ! 180                                                   
        a1lo(I,l)=A(I)                                                   
        b1lo(I,l)=B(I)
     ENDDO        ! 180
     WRITE(6,8) L,Plo(l),DPlo(l),PElo(l),DPElo(l),PEIlo(l),NODEL,NODES,NODE                         
  ENDDO  ! 170  continue

     
  do l=0,lomax
     do jlo=1,ilo(l)
        call abc(l,jlo,jatom)
     enddo
  enddo
  RETURN                                                            
  !                                                                       
2 FORMAT(5E14.7)                                                    
3 FORMAT(16I5)                                                      
4 FORMAT(I4,E16.7)                                                  
5 FORMAT(10X,' ENERGY PARAMETERS ARE',7F7.2)                        
6 FORMAT(10X,'E(',I2,2H)=,F10.4)                                    
7 FORMAT(/10X,'ATOMIC PARAMETERS FOR ',A10/)                        
14 FORMAT(/11X,1HL,5X,4HU(R),10X,5HU'(R),9X,5HDU/DE,8X,6HDU'/DE,6X,7HNORM-U')               
8 FORMAT(10X,I2,5E14.6,5X,3I2)                                      
9 FORMAT(10X,'FOR L=',I2,' CORRECTION=',E14.5,' OVERLAP=',E14.5)    
11 FORMAT(7F10.5)                                                    
13 FORMAT(////,':POS',i2.2,':',1x,'AT.NR.',I3,2X,'POSITION =',3F8.5,2X,'MULTIPLICITY =',I3)
1040 FORMAT(8I2)                                                       
1050 FORMAT(A10,I5,5X,2F10.9,I5,F5.2)                                  
1060 FORMAT(//,3X,'NOT EQUIV ATOM ',A10)
1080 FORMAT(13X,'EQUIV ATOM ',I3,3X,'POSITION: ',3F8.3)                
1980 FORMAT(3X)                                                        
2000 FORMAT(15X,I3//)                                                  
2030 FORMAT(///)                                                       
2031 FORMAT(/)                                                         
2032 FORMAT(49X,I3,//)
END SUBROUTINE ATPAR
