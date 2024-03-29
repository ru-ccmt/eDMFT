SUBROUTINE psplit(XWT1,xwteh,xwtel,xwt1h,xwt1l,jatom,nemin,nemax,test,eqbad,jatombad,lbad) 
  USE char,  ONLY:
  USE param, ONLY: lxdos
  USE bandm, ONLY: eseper
  !USE com,   ONLY: weigh
  USE charp, ONLY: apx, bpx, cpx, capx, cbpx, acpx, bcpx, apy, bpy, cpy, capy, cbpy, acpy, bcpy, apz, bpz, cpz, capz, cbpz, acpz, bcpz
  USE chard, ONLY: adxy, bdxy, cdxy, acdxy, bcdxy, adz2, bdz2, cdz2, acdz2, cadz2, cbdz2, bcdz2, adx2y2, bdx2y2, cdx2y2, cadx2y2, acdx2y2, bcdx2y2, cbdx2y2, cadxy, cbdxy, adxy, bdxy, adxz, bdxz, cdxz, cadxz, acdxz, bcdxz, cbdxz, adyz, bdyz, cdyz, cadyz, acdyz, bcdyz, cbdyz, bcdyz
  USE charf, ONLY: af00, bf00, af11, bf11, af22, bf22, af33, bf33, af1m, bf1m, af2m, bf2m, af3m, bf3m
  USE lohelp,ONLY: tc22, tc12, tce12
  USE structure, ONLY: isplit
  USE xa,    ONLY: E, tc100, tca100, tcb100, weight
  !USE xdos
  IMPLICIT REAL*8 (A-H,O-Z)
  !
  REAL*8, intent(inout) :: XWT1(0:21), xwteh(0:3),xwtel(0:3),xwt1h(0:3),xwt1l(0:3) 
  !common /xmean/  xwteh(0:3),xwtel(0:3),xwt1h(0:3),xwt1l(0:3) 
  lxdos2=(lxdos+1)*(lxdos+1)
  !
  num_loop: DO NUM=NEMIN,NEMAX                                           
     DO L=0,6
        XWT1(L)=XWT1(L)+TC100(L,NUM)/100.D0*WEIGHT(NUM)
        if(l.lt.4) then
           if(e(num).lt.eseper) then
              xwt1l(l)=xwt1l(l)+TC100(L,NUM)/100.D0*WEIGHT(NUM) 
              xwtel(l)=xwtel(l)+e(num)*TC100(L,NUM)/100.D0*WEIGHT(NUM) 
           else
              xwt1h(l)=xwt1h(l)+TC100(L,NUM)/100.D0*WEIGHT(NUM) 
              xwteh(l)=xwteh(l)+e(num)*TC100(L,NUM)/100.D0*WEIGHT(NUM) 
           endif
        endif
        !
        IF((ISPLIT(JATOM).EQ.15).AND.L.EQ.3) THEN    
           f00=Af00(NUM) + Bf00(NUM)                                        
           f11=Af11(NUM) + Bf11(NUM)                                
           f22=Af22(NUM) + Bf22(NUM)                                      
           f33=Af33(NUM) + Bf33(NUM)                                      
           f1m=Af1m(NUM) + Bf1m(NUM)                                      
           f2m=Af2m(NUM) + Bf2m(NUM)                                      
           f3m=Af3m(NUM) + Bf3m(NUM)                                      
           XWT1(15)=XWT1(15)+f00/100.D0*WEIGHT(NUM) 
           XWT1(16)=XWT1(16)+f11/100.D0*WEIGHT(NUM) 
           XWT1(17)=XWT1(17)+f22/100.D0*WEIGHT(NUM) 
           XWT1(18)=XWT1(18)+f33/100.D0*WEIGHT(NUM) 
           XWT1(19)=XWT1(19)+f1m/100.D0*WEIGHT(NUM) 
           XWT1(20)=XWT1(20)+f2m/100.D0*WEIGHT(NUM) 
           XWT1(21)=XWT1(21)+f3m/100.D0*WEIGHT(NUM) 
        ENDIF
        !
        !
        IF((ISPLIT(JATOM).EQ.5.OR.ISPLIT(JATOM).EQ.8.or.isplit(jatom).eq.15).AND.L.EQ.2) THEN    
           DZ2=ADZ2(NUM)+BDZ2(NUM)+cdz2(num)+cadz2(num)+acdz2(num)+bcdz2(num)+cbdz2(num)
           DX2Y2=ADX2Y2(NUM)+BDX2Y2(NUM)+cdx2y2(num)+cadx2y2(num)+acdx2y2(num)+bcdx2y2(num)+cbdx2y2(num)                        
           DXY=ADXY(NUM)+BDXY(NUM)+cdxy(num)+cadxy(num)+acdxy(num)+bcdxy(num)+cbdxy(num)                                      
           DXZ=ADXZ(NUM)+BDXZ(NUM)+cdxz(num)+cadxz(num)+acdxz(num)+bcdxz(num)+cbdxz(num)                                      
           DYZ=ADYZ(NUM)+BDYZ(NUM)+cdyz(num)+cadyz(num)+acdyz(num)+bcdyz(num)+cbdyz(num)                                      
           XWT1(10)=XWT1(10)+DZ2/100.D0*WEIGHT(NUM)                       
           XWT1(11)=XWT1(11)+DX2Y2/100.D0*WEIGHT(NUM)                     
           XWT1(12)=XWT1(12)+DXY/100.D0*WEIGHT(NUM)                       
           XWT1(13)=XWT1(13)+DXZ/100.D0*WEIGHT(NUM)                       
           XWT1(14)=XWT1(14)+DYZ/100.D0*WEIGHT(NUM)                       
        ENDIF
        !
        !
        !                                                                  
        IF((ISPLIT(JATOM).EQ.6.OR.ISPLIT(JATOM).EQ.8.or.isplit(jatom).eq.15).AND.L.EQ.1 ) THEN   
           PX=APX(NUM)+BPX(NUM)+cpx(num)+capx(num)+acpx(num)+cbpx(num)+bcpx(num)                                        
           PY=APY(NUM)+BPY(NUM)+cpy(num)+capy(num)+acpy(num)+cbpy(num)+bcpy(num)                                        
           PZ=APZ(NUM)+BPZ(NUM)+cpz(num)+capz(num)+acpz(num)+cbpz(num)+bcpz(num)                                        
           XWT1(7)=XWT1(7)+PX/100.D0*WEIGHT(NUM)                             
           XWT1(8)=XWT1(8)+PY/100.D0*WEIGHT(NUM)                             
           XWT1(9)=XWT1(9)+PZ/100.D0*WEIGHT(NUM)                             
        ENDIF
        !   
        !     SPLIT D INTO 2-4 COMP (SEE ABOVE)  
        !
        IF( (ISPLIT(JATOM).EQ.-2.OR.ISPLIT(JATOM).EQ.3.OR.ISPLIT(JATOM).EQ. 2.OR.ISPLIT(JATOM).EQ.4).AND.L.EQ.2 ) THEN
           DZ2=ADZ2(NUM)+BDZ2(NUM)+cdz2(num)+cadz2(num)+acdz2(num)+bcdz2(num)+cbdz2(num)
           DX2Y2=ADX2Y2(NUM)+BDX2Y2(NUM)+cdx2y2(num)+cadx2y2(num)+acdx2y2(num)+bcdx2y2(num)+cbdx2y2(num)                        
           DXY=ADXY(NUM)+BDXY(NUM)+cdxy(num)+cadxy(num)+acdxy(num)+bcdxy(num)+cbdxy(num)                                      
           ADXZYZ=ADXZ(NUM)+ADYZ(NUM)                                        
           BDXZYZ=BDXZ(NUM)+BDYZ(NUM)                                        
           cDXZYZ=cDXZ(NUM)+cDYZ(NUM)                                        
           acDXZYZ=acDXZ(NUM)+acDYZ(NUM)                                    
           BcDXZYZ=BcDXZ(NUM)+BcDYZ(NUM)                                      
           DXZYZ=ADXZYZ+BDXZYZ+cDXZYZ+2*(acDXZYZ+bcDXZYZ)                   
           IF(ISPLIT(JATOM).EQ.-2) THEN                                      
              XWT1(10)=XWT1(10)+DZ2/100.D0*WEIGHT(NUM)                       
              XWT1(11)=XWT1(11)+DXY/100.D0*WEIGHT(NUM)                       
              XWT1(12)=XWT1(12)+DX2Y2/100.D0*WEIGHT(NUM)                     
              XWT1(13)=XWT1(13)+DXZYZ/100.D0*WEIGHT(NUM)                     
           ELSE IF(ISPLIT(JATOM).EQ.3.OR.ISPLIT(JATOM).EQ.4) THEN            
              DXY=DXY+DX2Y2                                                  
              ADXY(NUM)=ADXY(NUM)+ADX2Y2(NUM)                                
              BDXY(NUM)=BDXY(NUM)+BDX2Y2(NUM)                                
              cDXY(NUM)=cDXY(NUM)+cDX2Y2(NUM)                                
              acDXY(NUM)=acDXY(NUM)+acDX2Y2(NUM)            
              BcDXY(NUM)=BcDXY(NUM)+BcDX2Y2(NUM)          
              XWT1(10)=XWT1(10)+DZ2/100.D0*WEIGHT(NUM)                       
              XWT1(11)=XWT1(11)+DXY/100.D0*WEIGHT(NUM)                       
              XWT1(12)=XWT1(12)+DXZYZ/100.D0*WEIGHT(NUM)                     
           ELSE IF(ISPLIT(JATOM).EQ.2) THEN                                  
              DEG=DZ2+DX2Y2                                                  
              ADEG=ADZ2(NUM)+ADX2Y2(NUM)                                     
              BDEG=BDZ2(NUM)+BDX2Y2(NUM)                                     
              cDEG=cDZ2(NUM)+cDX2Y2(NUM)                                     
              AcDEG=AcDZ2(NUM)+AcDX2Y2(NUM)                                  
              BcDEG=BcDZ2(NUM)+BcDX2Y2(NUM)                                   
              DT2G=DXY+DXZYZ                                                 
              ADT2G=ADXY(NUM)+ADXZYZ                                         
              BDT2G=BDXY(NUM)+BDXZYZ                                         
              cDT2G=cDXY(NUM)+cDXZYZ                                         
              AcDT2G=AcDXY(NUM)+AcDXZYZ                                       
              BcDT2G=BcDXY(NUM)+BcDXZYZ                                       
              XWT1(10)=XWT1(10)+DEG/100.D0*WEIGHT(NUM)                       
              XWT1(11)=XWT1(11)+DT2G/100.D0*WEIGHT(NUM)                      
           END IF
        END IF
        !
        !     SPLIT P IN PZ, (PX,PY) 
        !
        IF((ISPLIT(JATOM).EQ.1.OR.ISPLIT(JATOM).EQ.4.or.ISPLIT(JATOM).EQ.-2).AND.L.EQ.1) THEN    
           PZ=APZ(NUM)+BPZ(NUM)+cpz(num)+capz(num)+acpz(num)+cbpz(num)+bcpz(num)                                        
           APXY=APX(NUM)+APY(NUM)                                            
           BPXY=BPX(NUM)+BPY(NUM)
           cPXY=cPX(NUM)+cPY(NUM)                                            
           AcPXY=AcPX(NUM)+AcPY(NUM)
           BcPXY=BcPX(NUM)+BcPY(NUM)
           PXY=APXY+BPXY+cPXY+2.d0*AcPXY+2.d0*BcPXY
           XWT1(7)=XWT1(7)+PZ/100.D0*WEIGHT(NUM)                             
           XWT1(8)=XWT1(8)+PXY/100.D0*WEIGHT(NUM)                            
        END IF
        BTEST=5.D0                                                       
        IF(TCB100(L,NUM).GT.TEST) then
           TEST=TCB100(L,NUM)
           EQBAD=E(NUM)
           jatombad=jatom
           lbad=l
        endif
     ENDDO
  ENDDO num_loop
  return

3 FORMAT(1X,' X2Y2:',F10.5,2X,5F10.3)                               
6 FORMAT(1X,'   Z2:',F10.5,2X,5F10.3)                               
7 FORMAT(1X,'   XY:',F10.5,2X,5F10.3)                               
8 FORMAT(1X,' XZYZ:',F10.5,2X,5F10.3)                               
11 FORMAT(1X,'XYX2Y2',F10.5,2X,5F10.3)                               
12 FORMAT(1X,' D-EG:',F10.5,2X,5F10.3)                               
13 FORMAT(1X,'D-T2G:',F10.5,2X,5F10.3)                               
61 FORMAT(1X,'  F00:',F10.5,2X,5F10.3)                               
62 FORMAT(1X,'  F11:',F10.5,2X,5F10.3)                               
63 FORMAT(1X,'  F22:',F10.5,2X,5F10.3)                               
64 FORMAT(1X,'  F33:',F10.5,2X,5F10.3)                               
65 FORMAT(1X,'  F1M:',F10.5,2X,5F10.3)                               
66 FORMAT(1X,'  F2M:',F10.5,2X,5F10.3)                               
67 FORMAT(1X,'  F3M:',F10.5,2X,5F10.3)                               
70 FORMAT(1X,'  DZ2:',F10.5,2X,5F10.3)                               
71 FORMAT(1X,'DX2Y2:',F10.5,2X,5F10.3)                               
72 FORMAT(1X,'  DXY:',F10.5,2X,5F10.3)                               
73 FORMAT(1X,'  DXZ:',F10.5,2X,5F10.3)                               
74 FORMAT(1X,'  DYZ:',F10.5,2X,5F10.3)                               
80 FORMAT(1X,'   PX:',F10.5,2X,5F10.3)                               
81 FORMAT(1X,'   PY:',F10.5,2X,5F10.3)                               
82 FORMAT(1X,'   PZ:',F10.5,2X,5F10.3)                               
204 FORMAT(1X,' BAND#',I4,'  E=',F9.5,'  WEIGHT=',F10.7)             
205 FORMAT(/,1X,' K-POINT:',3F8.4,2X,2I4,2X,A10)                      
206 FORMAT(2X,'L=',I2,1X,F10.5,2X,5F10.3)                             
207 FORMAT(2X,'XDOS',1X,(10F10.5))                             
END SUBROUTINE psplit
