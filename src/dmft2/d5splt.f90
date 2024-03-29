SUBROUTINE D5SPLT (ALM,BLM,clm,MULT,UENORM,num,coord,jatom)                 
!                                                                       
!     DECOMPOSITION OF D CHARGE IN DZ2,DX2-Y2,DXY,DXZ,DYZ               
!     new (and reversed) definitions for dxz, dyz (6.2.89., Blaha)
!                                                                      
  use param
  USE chard
  use lo
  IMPLICIT REAL*8 (A-H,O-Z)
  COMPLEX*16       ALM,BLM,clm,CSUMa,csumb,csumc
  DIMENSION      ALM((lmax2+1)*(lmax2+1)),BLM((lmax2+1)*(lmax2+1))
  DIMENSION        cLM((lomax+1)*(lomax+1))
  CHARACTER*5 COORD                                                 
!---------------------------------------------------------------------  
!                                                                       
  ipip=max(ilo(2),1)
  SQRT2=SQRT(2.D0)                                                  
  SQRT3=SQRT(3.D0)                                                  
  IF (COORD.EQ.'OCTAH')  THEN                                       
     !
     ! non standard dsplit  
     ! classical octahedral coordinates                          
     !                                                                       
     IF (JATOM.EQ.1) THEN                                              
        !     ===========================                                       
        !     1. METALLATOM                                                     
        !                                                                       
        !
        CSUMa=( SQRT3*(ALM(9)+ALM(5))/SQRT2 - ALM(7) )*0.5                 
        CSUMb=( SQRT3*(bLM(9)+bLM(5))/SQRT2 - bLM(7) )*0.5                 
        CSUMc=( SQRT3*(cLM(9)+cLM(5))/SQRT2 - cLM(7) )*0.5                 
        ADz2(num)=ADz2(num)+CSUMa*CONJG(CSUMa)*100.D0/MULT                      
        BDz2(num)=BDz2(num)+CSUMb*CONJG(CSUMb)*UENORM*100.D0/MULT               
        cdz2(num)=cdz2(num)+CSUMc*CONJG(CSUMc)*100.D0/MULT
        cadz2(num)=cadz2(num)+CSUMc*CONJG(CSUMa)*pi12lo(ipip,2)*100.D0/MULT
        acdz2(num)=acdz2(num)+CSUMa*CONJG(CSUMc)*pi12lo(ipip,2)*100.D0/MULT
        cbdz2(num)=cbdz2(num)+CSUMc*CONJG(CSUMb)*pe12lo(ipip,2)*100.D0/MULT
        bcdz2(num)=bcdz2(num)+CSUMb*CONJG(CSUMc)*pe12lo(ipip,2)*100.D0/MULT
!                                                                       
        CSUMa=(-1)*(ALM(8)+ALM(6))/SQRT2                                   
        CSUMb=(-1)*(BLM(8)+BLM(6))/SQRT2                                   
        CSUMc=(-1)*(cLM(8)+cLM(6))/SQRT2                                   
        ADX2Y2(num)=ADX2Y2(num)+CSUMa*CONJG(CSUMa)*100.D0/MULT
        BDX2Y2(num)=BDX2Y2(num)+CSUMb*CONJG(CSUMb)*UENORM*100.D0/MULT
        cdx2y2(num)=cdx2y2(num)+CSUMc*CONJG(CSUMc)*100.D0/MULT
        cadx2y2(num)=cadx2y2(num)+ &
             CSUMc*CONJG(CSUMa)*pi12lo(ipip,2)*100.D0/MULT
        acdx2y2(num)=acdx2y2(num)+ &
             CSUMa*CONJG(CSUMc)*pi12lo(ipip,2)*100.D0/MULT
        cbdx2y2(num)=cbdx2y2(num)+ &
             CSUMc*CONJG(CSUMb)*pe12lo(ipip,2)*100.D0/MULT
        bcdx2y2(num)=bcdx2y2(num)+ &
             CSUMb*CONJG(CSUMc)*pe12lo(ipip,2)*100.D0/MULT
!                                                                       
        CSUMa=( (ALM(9)+ALM(5))/SQRT2 + SQRT3*ALM(7) )*0.5                 
        CSUMb=( (BLM(9)+BLM(5))/SQRT2 + SQRT3*BLM(7) )*0.5                 
        CSUMc=( (cLM(9)+cLM(5))/SQRT2 + SQRT3*cLM(7) )*0.5                 
        ADXY(num)=ADXY(num)+CSUMa*CONJG(CSUMa)*100.D0/MULT                      
        BDXY(num)=BDXY(num)+CSUMb*CONJG(CSUMb)*UENORM*100.D0/MULT               
        cdxy(num)=cdxy(num)+CSUMc*CONJG(CSUMc)*100.D0/MULT
        cadxy(num)=cadxy(num)+CSUMc*CONJG(CSUMa)*pi12lo(ipip,2)*100.D0/MULT
        acdxy(num)=acdxy(num)+CSUMa*CONJG(CSUMc)*pi12lo(ipip,2)*100.D0/MULT
        cbdxy(num)=cbdxy(num)+CSUMc*CONJG(CSUMb)*pe12lo(ipip,2)*100.D0/MULT
        bcdxy(num)=bcdxy(num)+CSUMb*CONJG(CSUMc)*pe12lo(ipip,2)*100.D0/MULT
!                                                                       
        CSUMa=((ALM(9)-ALM(5))/SQRT2-(ALM(8)-ALM(6))/SQRT2)/SQRT2      
        CSUMb=((BLM(9)-BLM(5))/SQRT2-(BLM(8)-BLM(6))/SQRT2)/SQRT2      
        CSUMc=((cLM(9)-cLM(5))/SQRT2-(cLM(8)-cLM(6))/SQRT2)/SQRT2      
        ADXz(num)=ADXz(num)+CSUMa*CONJG(CSUMa)*100.D0/MULT                      
        BDXz(num)=BDXz(num)+CSUMb*CONJG(CSUMb)*UENORM*100.D0/MULT               
        cdxz(num)=cdxz(num)+CSUMc*CONJG(CSUMc)*100.D0/MULT
        cadxz(num)=cadxz(num)+CSUMc*CONJG(CSUMa)*pi12lo(ipip,2)*100.D0/MULT
        acdxz(num)=acdxz(num)+CSUMa*CONJG(CSUMc)*pi12lo(ipip,2)*100.D0/MULT
        cbdxz(num)=cbdxz(num)+CSUMc*CONJG(CSUMb)*pe12lo(ipip,2)*100.D0/MULT
        bcdxz(num)=bcdxz(num)+CSUMb*CONJG(CSUMc)*pe12lo(ipip,2)*100.D0/MULT
        !                                                                       
        CSUMa=((-1.)*(ALM(9)-ALM(5))/SQRT2-(ALM(8)-ALM(6))/SQRT2)/SQRT2  
        CSUMb=((-1.)*(BLM(9)-BLM(5))/SQRT2-(BLM(8)-BLM(6))/SQRT2)/SQRT2  
        CSUMc=((-1.)*(cLM(9)-cLM(5))/SQRT2-(cLM(8)-cLM(6))/SQRT2)/SQRT2  
        ADYz(num)=ADYz(num)+CSUMa*CONJG(CSUMa)*100.D0/MULT                      
        BDYz(num)=BDYz(num)+CSUMb*CONJG(CSUMb)*UENORM*100.D0/MULT               
        cdyz(num)=cdyz(num)+CSUMc*CONJG(CSUMc)*100.D0/MULT
        cadyz(num)=cadyz(num)+CSUMc*CONJG(CSUMa)*pi12lo(ipip,2)*100.D0/MULT
        acdyz(num)=acdyz(num)+CSUMa*CONJG(CSUMc)*pi12lo(ipip,2)*100.D0/MULT
        cbdyz(num)=cbdyz(num)+CSUMc*CONJG(CSUMb)*pe12lo(ipip,2)*100.D0/MULT
        bcdyz(num)=bcdyz(num)+CSUMb*CONJG(CSUMc)*pe12lo(ipip,2)*100.D0/MULT
        !                                                                       
        
                                                                        
                                                                        
     ELSE                                                              
        !     ===========================                                       
        !     2.METALLATOM                                                      
        !                                                                       
        CSUMa=(-SQRT3*(ALM(9)+ALM(5))/SQRT2-ALM(7) )*0.5                
        CSUMb=(-SQRT3*(BLM(9)+BLM(5))/SQRT2-BLM(7) )*0.5                
        CSUMc=(-SQRT3*(cLM(9)+cLM(5))/SQRT2-cLM(7) )*0.5                
        ADz2(num)=ADz2(num)+CSUMa*CONJG(CSUMa)*100.D0/MULT                      
        BDz2(num)=BDz2(num)+CSUMb*CONJG(CSUMb)*UENORM*100.D0/MULT               
        cdz2(num)=cdz2(num)+CSUMc*CONJG(CSUMc)*100.D0/MULT
        cadz2(num)=cadz2(num)+CSUMc*CONJG(CSUMa)*pi12lo(ipip,2)*100.D0/MULT
        acdz2(num)=acdz2(num)+CSUMa*CONJG(CSUMc)*pi12lo(ipip,2)*100.D0/MULT
        cbdz2(num)=cbdz2(num)+CSUMc*CONJG(CSUMb)*pe12lo(ipip,2)*100.D0/MULT
        bcdz2(num)=bcdz2(num)+CSUMb*CONJG(CSUMc)*pe12lo(ipip,2)*100.D0/MULT
        !                                                                       
        CSUMa=(-1.)*(ALM(8)-ALM(6))/SQRT2                                  
        CSUMb=(-1.)*(BLM(8)-BLM(6))/SQRT2                                  
        CSUMc=(-1.)*(cLM(8)-cLM(6))/SQRT2                                  
        ADX2Y2(num)=ADX2Y2(num)+CSUMa*CONJG(CSUMa)*100.D0/MULT         
        BDX2Y2(num)=BDX2Y2(num)+CSUMb*CONJG(CSUMb)*UENORM*100.D0/MULT
        cdx2y2(num)=cdx2y2(num)+CSUMc*CONJG(CSUMc)*100.D0/MULT
        cadx2y2(num)=cadx2y2(num)+ &
             CSUMc*CONJG(CSUMa)*pi12lo(ipip,2)*100.D0/MULT
        acdx2y2(num)=acdx2y2(num)+ &
             CSUMa*CONJG(CSUMc)*pi12lo(ipip,2)*100.D0/MULT
        cbdx2y2(num)=cbdx2y2(num)+ &
             CSUMc*CONJG(CSUMb)*pe12lo(ipip,2)*100.D0/MULT
        bcdx2y2(num)=bcdx2y2(num)+ &
             CSUMb*CONJG(CSUMc)*pe12lo(ipip,2)*100.D0/MULT
        !                                                                       
        CSUMa=( (-1.)*(ALM(9)+ALM(5))/SQRT2 + SQRT3*ALM(7) )*0.5           
        CSUMb=( (-1.)*(BLM(9)+BLM(5))/SQRT2 + SQRT3*BLM(7) )*0.5           
        CSUMc=( (-1.)*(cLM(9)+cLM(5))/SQRT2 + SQRT3*cLM(7) )*0.5           
        ADXY(num)=ADXY(num)+CSUMa*CONJG(CSUMa)*100.D0/MULT                      
        BDXY(num)=BDXY(num)+CSUMb*CONJG(CSUMb)*UENORM*100.D0/MULT               
        cdxy(num)=cdxy(num)+CSUMc*CONJG(CSUMc)*100.D0/MULT
        cadxy(num)=cadxy(num)+CSUMc*CONJG(CSUMa)*pi12lo(ipip,2)*100.D0/MULT
        acdxy(num)=acdxy(num)+CSUMa*CONJG(CSUMc)*pi12lo(ipip,2)*100.D0/MULT
        cbdxy(num)=cbdxy(num)+CSUMc*CONJG(CSUMb)*pe12lo(ipip,2)*100.D0/MULT
        bcdxy(num)=bcdxy(num)+CSUMb*CONJG(CSUMc)*pe12lo(ipip,2)*100.D0/MULT
        !                                                                       
        CSUMa=((ALM(9)-ALM(5))/SQRT2-(ALM(8)+ALM(6))/SQRT2 )/SQRT2      
        CSUMb=((BLM(9)-BLM(5))/SQRT2-(BLM(8)+BLM(6))/SQRT2 )/SQRT2      
        CSUMc=((cLM(9)-cLM(5))/SQRT2-(cLM(8)+cLM(6))/SQRT2 )/SQRT2      
        ADXz(num)=ADXz(num)+CSUMa*CONJG(CSUMa)*100.D0/MULT                      
        BDXz(num)=BDXz(num)+CSUMb*CONJG(CSUMb)*UENORM*100.D0/MULT               
        cdxz(num)=cdxz(num)+CSUMc*CONJG(CSUMc)*100.D0/MULT
        cadxz(num)=cadxz(num)+CSUMc*CONJG(CSUMa)*pi12lo(ipip,2)*100.D0/MULT
        acdxz(num)=acdxz(num)+CSUMa*CONJG(CSUMc)*pi12lo(ipip,2)*100.D0/MULT
        cbdxz(num)=cbdxz(num)+CSUMc*CONJG(CSUMb)*pe12lo(ipip,2)*100.D0/MULT
        bcdxz(num)=bcdxz(num)+CSUMb*CONJG(CSUMc)*pe12lo(ipip,2)*100.D0/MULT
        !                                                                       
        CSUMa=((-1.)*(ALM(9)-ALM(5))/SQRT2-(ALM(8)+ALM(6))/SQRT2)/SQRT2  
        CSUMb=((-1.)*(BLM(9)-BLM(5))/SQRT2-(BLM(8)+BLM(6))/SQRT2)/SQRT2  
        CSUMc=((-1.)*(cLM(9)-cLM(5))/SQRT2-(cLM(8)+cLM(6))/SQRT2)/SQRT2  
        ADYz(num)=ADYz(num)+CSUMa*CONJG(CSUMa)*100.D0/MULT                      
        BDYz(num)=BDYz(num)+CSUMb*CONJG(CSUMb)*UENORM*100.D0/MULT               
        cdyz(num)=cdyz(num)+CSUMc*CONJG(CSUMc)*100.D0/MULT
        cadyz(num)=cadyz(num)+CSUMc*CONJG(CSUMa)*pi12lo(ipip,2)*100.D0/MULT
        acdyz(num)=acdyz(num)+CSUMa*CONJG(CSUMc)*pi12lo(ipip,2)*100.D0/MULT
        cbdyz(num)=cbdyz(num)+CSUMc*CONJG(CSUMb)*pe12lo(ipip,2)*100.D0/MULT
        bcdyz(num)=bcdyz(num)+CSUMb*CONJG(CSUMc)*pe12lo(ipip,2)*100.D0/MULT
        !                                                                       
     END IF
     !                                                                        
  ELSE if (COORD.EQ.'TRIGO')  THEN                                       
     !
     CSUMa=(1/SQRT3)*((ALM(8)-ALM(6))+(ALM(9)+ALM(5))/SQRT2)                 
     CSUMb=(1/SQRT3)*((bLM(8)-bLM(6))+(bLM(9)+bLM(5))/SQRT2)             
     CSUMc=(1/SQRT3)*((cLM(8)-cLM(6))+(cLM(9)+cLM(5))/SQRT2)               
     ADz2(num)=ADz2(num)+CSUMa*CONJG(CSUMa)*100.D0/MULT                      
     BDz2(num)=BDz2(num)+CSUMb*CONJG(CSUMb)*UENORM*100.D0/MULT               
     cdz2(num)=cdz2(num)+CSUMc*CONJG(CSUMc)*100.D0/MULT
     cadz2(num)=cadz2(num)+CSUMc*CONJG(CSUMa)*pi12lo(ipip,2)*100.D0/MULT
     acdz2(num)=acdz2(num)+CSUMa*CONJG(CSUMc)*pi12lo(ipip,2)*100.D0/MULT
     cbdz2(num)=cbdz2(num)+CSUMc*CONJG(CSUMb)*pe12lo(ipip,2)*100.D0/MULT
     bcdz2(num)=bcdz2(num)+CSUMb*CONJG(CSUMc)*pe12lo(ipip,2)*100.D0/MULT
     !                                                                       
     CSUMa=(1/SQRT3)*((aLM(8)+aLM(6))-(aLM(9)-aLM(5))/SQRT2)                 
     CSUMb=(1/SQRT3)*((bLM(8)+bLM(6))-(bLM(9)-bLM(5))/SQRT2)                 
     CSUMc=(1/SQRT3)*((cLM(8)+cLM(6))-(cLM(9)-cLM(5))/SQRT2)                 
     ADX2Y2(num)=ADX2Y2(num)+CSUMa*CONJG(CSUMa)*100.D0/MULT                      
     BDX2Y2(num)=BDX2Y2(num)+CSUMb*CONJG(CSUMb)*UENORM*100.D0/MULT               
     cdx2y2(num)=cdx2y2(num)+CSUMc*CONJG(CSUMc)*100.D0/MULT
     cadx2y2(num)=cadx2y2(num)+ &
          CSUMc*CONJG(CSUMa)*pi12lo(ipip,2)*100.D0/MULT
     acdx2y2(num)=acdx2y2(num)+ &
          CSUMa*CONJG(CSUMc)*pi12lo(ipip,2)*100.D0/MULT
     cbdx2y2(num)=cbdx2y2(num)+ &
          CSUMc*CONJG(CSUMb)*pe12lo(ipip,2)*100.D0/MULT
     bcdx2y2(num)=bcdx2y2(num)+ &
          CSUMb*CONJG(CSUMc)*pe12lo(ipip,2)*100.D0/MULT
     !                                                                       
     CSUMa=ALM(7)                  
     CSUMb=BLM(7)                 
     CSUMc=cLM(7)                  
     ADXY(num)=ADXY(num)+CSUMa*CONJG(CSUMa)*100.D0/MULT                      
     BDXY(num)=BDXY(num)+CSUMb*CONJG(CSUMb)*UENORM*100.D0/MULT               
     cdxy(num)=cdxy(num)+CSUMc*CONJG(CSUMc)*100.D0/MULT
     cadxy(num)=cadxy(num)+CSUMc*CONJG(CSUMa)*pi12lo(ipip,2)*100.D0/MULT
     acdxy(num)=acdxy(num)+CSUMa*CONJG(CSUMc)*pi12lo(ipip,2)*100.D0/MULT
     cbdxy(num)=cbdxy(num)+CSUMc*CONJG(CSUMb)*pe12lo(ipip,2)*100.D0/MULT
     bcdxy(num)=bcdxy(num)+CSUMb*CONJG(CSUMc)*pe12lo(ipip,2)*100.D0/MULT
     !                                                                       
     CSUMa=(1/SQRT3)*((ALM(8)+ALM(6))/SQRT2+(ALM(9)-ALM(5)))
     CSUMb=(1/SQRT3)*((bLM(8)+bLM(6))/SQRT2+(bLM(9)-bLM(5)))
     CSUMc=(1/SQRT3)*((cLM(8)+cLM(6))/SQRT2+(cLM(9)-cLM(5)))         
     ADXz(num)=ADXz(num)+CSUMa*CONJG(CSUMa)*100.D0/MULT                      
     BDXz(num)=BDXz(num)+CSUMb*CONJG(CSUMb)*UENORM*100.D0/MULT               
     cdxz(num)=cdxz(num)+CSUMc*CONJG(CSUMc)*100.D0/MULT
     cadxz(num)=cadxz(num)+CSUMc*CONJG(CSUMa)*pi12lo(ipip,2)*100.D0/MULT
     acdxz(num)=acdxz(num)+CSUMa*CONJG(CSUMc)*pi12lo(ipip,2)*100.D0/MULT
     cbdxz(num)=cbdxz(num)+CSUMc*CONJG(CSUMb)*pe12lo(ipip,2)*100.D0/MULT
     bcdxz(num)=bcdxz(num)+CSUMb*CONJG(CSUMc)*pe12lo(ipip,2)*100.D0/MULT
     !                                                                       
     CSUMa=(1/SQRT3)*((ALM(8)-ALM(6))/SQRT2-(ALM(9)+ALM(5)))
     CSUMb=(1/SQRT3)*((bLM(8)-bLM(6))/SQRT2-(bLM(9)+bLM(5)))            
     CSUMc=(1/SQRT3)*((cLM(8)-cLM(6))/SQRT2-(cLM(9)+cLM(5)))       
     ADYz(num)=ADYz(num)+CSUMa*CONJG(CSUMa)*100.D0/MULT                      
     BDYz(num)=BDYz(num)+CSUMb*CONJG(CSUMb)*UENORM*100.D0/MULT               
     cdyz(num)=cdyz(num)+CSUMc*CONJG(CSUMc)*100.D0/MULT
     cadyz(num)=cadyz(num)+CSUMc*CONJG(CSUMa)*pi12lo(ipip,2)*100.D0/MULT
     acdyz(num)=acdyz(num)+CSUMa*CONJG(CSUMc)*pi12lo(ipip,2)*100.D0/MULT
     cbdyz(num)=cbdyz(num)+CSUMc*CONJG(CSUMb)*pe12lo(ipip,2)*100.D0/MULT
     bcdyz(num)=bcdyz(num)+CSUMb*CONJG(CSUMc)*pe12lo(ipip,2)*100.D0/MULT
     
  ELSE                                        
     
     !
     ! standard-split
     !                                                             
     ADZ2(num)=ADZ2(num)+ALM(7)*CONJG(ALM(7))*100.D0/MULT                     
     BDZ2(num)=BDZ2(num)+BLM(7)*CONJG(BLM(7))*UENORM*100.D0/MULT              
     cdz2(num)=cdz2(num)+cLM(7)*CONJG(cLM(7))*100.D0/MULT               
     cadz2(num)=cadz2(num)+Clm(7)*CONJG(alm(7))*pi12lo(ipip,2)*100.D0/MULT
     acdz2(num)=acdz2(num)+alm(7)*CONJG(Clm(7))*pi12lo(ipip,2)*100.D0/MULT
     cbdz2(num)=cbdz2(num)+Clm(7)*CONJG(blm(7))*pe12lo(ipip,2)*100.D0/MULT
     bcdz2(num)=bcdz2(num)+blm(7)*CONJG(Clm(7))*pe12lo(ipip,2)*100.D0/MULT
     !                                                                       
     CSUMa=(ALM(9)+ALM(5))/SQRT2                                        
     CSUMb=(BLM(9)+BLM(5))/SQRT2                                        
     CSUMc=(cLM(9)+cLM(5))/SQRT2                                        
     ADX2Y2(num)=ADX2Y2(num)+CSUMa*CONJG(CSUMa)*100.D0/MULT               
     BDX2Y2(num)=BDX2Y2(num)+CSUMb*CONJG(CSUMb)*UENORM*100.D0/MULT      
     cdx2y2(num)=cdx2y2(num)+CSUMc*CONJG(CSUMc)*100.D0/MULT
     cadx2y2(num)=cadx2y2(num)+ &
          CSUMc*CONJG(CSUMa)*pi12lo(ipip,2)*100.D0/MULT
     acdx2y2(num)=acdx2y2(num)+ &
          CSUMa*CONJG(CSUMc)*pi12lo(ipip,2)*100.D0/MULT
     cbdx2y2(num)=cbdx2y2(num)+ &
          CSUMc*CONJG(CSUMb)*pe12lo(ipip,2)*100.D0/MULT
     bcdx2y2(num)=bcdx2y2(num)+ &
          CSUMb*CONJG(CSUMc)*pe12lo(ipip,2)*100.D0/MULT
     !                                                                       
     CSUMa=(ALM(9)-ALM(5))/SQRT2                                        
     CSUMb=(BLM(9)-BLM(5))/SQRT2                                        
     CSUMc=(cLM(9)-cLM(5))/SQRT2                                        
     ADXY(num)=ADXY(num)+CSUMa*CONJG(CSUMa)*100.D0/MULT                      
     BDXY(num)=BDXY(num)+CSUMb*CONJG(CSUMb)*UENORM*100.D0/MULT               
     cdxy(num)=cdxy(num)+CSUMc*CONJG(CSUMc)*100.D0/MULT
     cadxy(num)=cadxy(num)+CSUMc*CONJG(CSUMa)*pi12lo(ipip,2)*100.D0/MULT
     acdxy(num)=acdxy(num)+CSUMa*CONJG(CSUMc)*pi12lo(ipip,2)*100.D0/MULT
     cbdxy(num)=cbdxy(num)+CSUMc*CONJG(CSUMb)*pe12lo(ipip,2)*100.D0/MULT
     bcdxy(num)=bcdxy(num)+CSUMb*CONJG(CSUMc)*pe12lo(ipip,2)*100.D0/MULT
     !                                                                       
     CSUMa=(ALM(8)-ALM(6))/SQRT2                                        
     CSUMb=(BLM(8)-BLM(6))/SQRT2                                        
     CSUMc=(cLM(8)-cLM(6))/SQRT2                                        
     ADXz(num)=ADXz(num)+CSUMa*CONJG(CSUMa)*100.D0/MULT                      
     BDXz(num)=BDXz(num)+CSUMb*CONJG(CSUMb)*UENORM*100.D0/MULT               
     cdxz(num)=cdxz(num)+CSUMc*CONJG(CSUMc)*100.D0/MULT
     cadxz(num)=cadxz(num)+CSUMc*CONJG(CSUMa)*pi12lo(ipip,2)*100.D0/MULT
     acdxz(num)=acdxz(num)+CSUMa*CONJG(CSUMc)*pi12lo(ipip,2)*100.D0/MULT
     cbdxz(num)=cbdxz(num)+CSUMc*CONJG(CSUMb)*pe12lo(ipip,2)*100.D0/MULT
     bcdxz(num)=bcdxz(num)+CSUMb*CONJG(CSUMc)*pe12lo(ipip,2)*100.D0/MULT
     !                                                                       
     CSUMa=(ALM(8)+ALM(6))/SQRT2                                        
     CSUMb=(BLM(8)+BLM(6))/SQRT2                                        
     CSUMc=(cLM(8)+cLM(6))/SQRT2                                        
     ADYz(num)=ADYz(num)+CSUMa*CONJG(CSUMa)*100.D0/MULT                      
     BDYz(num)=BDYz(num)+CSUMb*CONJG(CSUMb)*UENORM*100.D0/MULT               
     cdyz(num)=cdyz(num)+CSUMc*CONJG(CSUMc)*100.D0/MULT
     cadyz(num)=cadyz(num)+CSUMc*CONJG(CSUMa)*pi12lo(ipip,2)*100.D0/MULT
     acdyz(num)=acdyz(num)+CSUMa*CONJG(CSUMc)*pi12lo(ipip,2)*100.D0/MULT
     cbdyz(num)=cbdyz(num)+CSUMc*CONJG(CSUMb)*pe12lo(ipip,2)*100.D0/MULT
     bcdyz(num)=bcdyz(num)+CSUMb*CONJG(CSUMc)*pe12lo(ipip,2)*100.D0/MULT
  endif
  RETURN                                                            
END SUBROUTINE D5SPLT
