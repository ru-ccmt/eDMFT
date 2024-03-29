SUBROUTINE csplit(nemin,nemax,l,jatom,mu,alm,blm,clm,coord)
  !
  ! performs charge analysis
  !                                                                       
  use param, ONLY: lomax, lxdos
  use defs, ONLY: PI
  use structure, ONLY: isplit, mult
  use lo, ONLY: ilo, pi12lo, pe12lo, pr12lo
  USE lohelp, ONLY: TC22, TC12, TC21, TCe12, TCe21, sum12, sum21, sume12, sume21, sum22
  USE xa, ONLY: TC100, TCA100, TCB100, LMAX2, NUME, nloat, suma, sumb, E 
  USE atspdt, only: p,dp,pe,dpe,pei
  IMPLICIT REAL*8 (A-H,O-Z)
  !
  CHARACTER*5:: COORD 
  REAL*8     :: tc21h,tc12h,tce21h
  REAL*8     :: tce12h,tc22h
  REAL*8     :: wtca100, wtcb100
  COMPLEX*16 :: BLM((LMAX2+1)*(LMAX2+1),NUME),ALM((LMAX2+1)*(LMAX2+1),NUME),cLM((lomax+1)*(lomax+1),NUME,nloat)
  !                                                                       
  ! DEFINE CONSTANTS, FACTORIALS                                      
  !
  ! need: SUMA, SUMB, SUM12, SUM21, SUME12, SUME21, SUM22
  !       
  ! creates: TCxxx
  ! nested subroutines produce:
  !          app variables in charp, chard, charf
  ! all are functions of L,num

  SQFP=SQRT(4.D0*PI)
  
  !print *, 'L=', L, 'PEI=', PEI(l)
  !do num=nemin,nemax
  !   wtca100 = suma(num)*SQFP * 100
  !   wtcb100 = sumb(num)*SQFP*PEI(l) * 100
  !   WRITE(*, '("num=",I3,1x,"E=",F8.4,2x,"wtca100=",F12.6,2x,"wtcb100=",F12.6)') num, E(num), wtca100, wtcb100
  !end do
  
  num_loop:  DO NUM=NEMIN,NEMAX
     TCA=SUMA(NUM)*SQFP
     TCB=SUMB(NUM)*SQFP*PEI(L)
     TC=TCA+TCB
     tc22h=0.d0                         
     tc12h=0.d0                          
     tc21h=0.d0                         
     tce12h=0.d0                         
     tce21h=0.d0                         
     IF(l.LE.lomax) THEN
        DO jlo=1,ilo(l)
           tc12h=sum12(num,jlo)*sqfp*pi12lo(jlo,l)
           tc21h=sum21(num,jlo)*sqfp*pi12lo(jlo,l)
           tce12h=sume12(num,jlo)*sqfp*pe12lo(jlo,l)
           tce21h=sume21(num,jlo)*sqfp*pe12lo(jlo,l)
           TC=TC+tc12h+tc21h+tce12h+tce21h   
           DO jlop=1,ilo(l)
              tc22h=sum22(num,jlop,jlo)*sqfp*pr12lo(jlop,jlo,l)
              tc=tc+tc22h
           ENDDO
        ENDDO
     ENDIF
     TC100(L,NUM)=TC*100.0D0+TC100(L,NUM)
     !                                                                       
     ! Create "CROSS"-partial charges 
     !                           
     ipip=max(ilo(l),1)
     !
     ! SPLIT P IN 3 components
     !                           
     IF( ISPLIT(JATOM).NE.2 .AND. L.EQ.1 ) CALL P3SPLT (ALM(1,NUM),BLM(1,NUM),cLM(1,NUM,ipip),MULT(JATOM),PEI(L),num,coord)
     !
     ! SPLIT D IN 5 components              
     !
     IF(L.EQ.2) CALL D5SPLT (ALM(1,NUM),BLM(1,NUM),cLM(1,NUM,ipip),MULT(JATOM),PEI(L),num,coord,jatom)        
     !                                                                       
     ! SPLIT F IN 7 components                
     !
     IF(L.EQ.3 .and. isplit(jatom).eq.15) CALL F7SPLT (ALM(1,NUM),BLM(1,NUM),MULT(JATOM),PEI(L),num,coord)
     !
     !
     !                                                                       
     TCA100(L,NUM)=TCA*100.0D0+TCA100(L,NUM)                         
     TCB100(L,NUM)=TCB*100.0D0+TCB100(L,NUM)                         
     TC22(L,NUM)=tc22h*100.0D0+TC22(L,NUM)                         
     TC12(L,NUM)=tc12h*100.0D0+TC12(L,NUM)                         
     TC21(L,NUM)=tc21h*100.0D0+TC21(L,NUM)                         
     TCe12(L,NUM)=tce12h*100.0D0+TCe12(L,NUM)                         
     TCe21(L,NUM)=tce21h*100.0D0+TCe21(L,NUM)                         
  ENDDO num_loop
  
  IF(MU.NE.MULT(JATOM)) RETURN !GOTO 55                                     
  !
  !.....WRITE OUT L,M DECOMPOSITIONS                                      
  !                                                                       
  DO NUM=NEMIN,NEMAX
     TC100(L,NUM)=TC100(L,NUM)/MULT(JATOM)                           
     TCA100(L,NUM)=TCA100(L,NUM)/MULT(JATOM)                         
     TCB100(L,NUM)=TCB100(L,NUM)/MULT(JATOM)                         
     TC22(L,NUM)=TC22(L,NUM)/MULT(JATOM)                         
     TC12(L,NUM)=TC12(L,NUM)/MULT(JATOM)                         
     TC21(L,NUM)=TC21(L,NUM)/MULT(JATOM)                         
     TCe12(L,NUM)=TCe12(L,NUM)/MULT(JATOM)                         
     TCe21(L,NUM)=TCe21(L,NUM)/MULT(JATOM)                         
  ENDDO
  return ! 55
END SUBROUTINE csplit
