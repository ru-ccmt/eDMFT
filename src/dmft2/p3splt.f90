SUBROUTINE P3SPLT (ALM,BLM,clm,MULT,UENORM,num,coord)    
  !                                                                       
  !     DECOMPOSITION OF P CHARGE IN PX,PY,PZ                             
  !                                                                       
  use param
  USE charp
  use lo
  IMPLICIT REAL*8 (A-H,O-Z)
  !
  COMPLEX*16     ALM,BLM,clm,CSUMa,csumb,csumc
  DIMENSION      ALM((lmax2+1)*(lmax2+1)),BLM((lmax2+1)*(lmax2+1))
  DIMENSION      cLM((lomax+1)*(lomax+1))
  CHARACTER*5 COORD                        
  ! uses:
  !   alm, blm, clm
  ! output:
  !   APX, BPX, CPX, CAPX,....
  !   APY, BPY, CPY, CAPY,....
  !   APZ, BPZ, CPZ, CAPZ,....
  !---------------------------------------------------------------------  
  !     new (and reversed) definitions for p-x and p-y (6.2.89, Blaha)
  !
  ipip=max(ilo(1),1)
  SQRT2=SQRT(2.D0)                                                  
  CSUMa=(ALM(4)-ALM(2))/SQRT2                                        
  CSUMb=(bLM(4)-bLM(2))/SQRT2                                        
  CSUMc=(cLM(4)-cLM(2))/SQRT2                                        
  APX(num)=APX(num)+CSUMa*CONJG(CSUMa)*100.D0/MULT
  BPX(num)=BPX(num)+CSUMb*CONJG(CSUMb)*UENORM*100.D0/MULT
  cPX(num)=cPX(num)+CSUMc*CONJG(CSUMc)*100.D0/MULT
  caPX(num)=caPX(num)+CSUMc*CONJG(CSUMa)*pi12lo(ipip,1)*100.D0/MULT
  acPX(num)=acPX(num)+CSUMa*CONJG(CSUMc)*pi12lo(ipip,1)*100.D0/MULT
  cbPX(num)=cbPX(num)+CSUMc*CONJG(CSUMb)*pe12lo(ipip,1)*100.D0/MULT
  bcPX(num)=bcPX(num)+CSUMb*CONJG(CSUMc)*pe12lo(ipip,1)*100.D0/MULT
!                                                                       
  CSUMa=(ALM(4)+ALM(2))/SQRT2                                        
  CSUMb=(bLM(4)+bLM(2))/SQRT2                                        
  CSUMc=(cLM(4)+cLM(2))/SQRT2                                        
  APy(num)=APy(num)+CSUMa*CONJG(CSUMa)*100.D0/MULT
  BPy(num)=BPy(num)+CSUMb*CONJG(CSUMb)*UENORM*100.D0/MULT                
  cPy(num)=cPy(num)+CSUMc*CONJG(CSUMc)*100.D0/MULT
  caPy(num)=caPy(num)+CSUMc*CONJG(CSUMa)*pi12lo(ipip,1)*100.D0/MULT
  acPy(num)=acPy(num)+CSUMa*CONJG(CSUMc)*pi12lo(ipip,1)*100.D0/MULT
  cbPy(num)=cbPy(num)+CSUMc*CONJG(CSUMb)*pe12lo(ipip,1)*100.D0/MULT
  bcPy(num)=bcPy(num)+CSUMb*CONJG(CSUMc)*pe12lo(ipip,1)*100.D0/MULT
!                                                                       
  APZ(num)=APZ(num)+ALM(3)*CONJG(ALM(3))*100.D0/MULT 
  BPZ(num)=BPZ(num)+BLM(3)*CONJG(BLM(3))*UENORM*100.D0/MULT
  cPZ(num)=cPZ(num)+cLM(3)*CONJG(cLM(3))*100.D0/MULT
  caPz(num)=caPz(num)+Clm(3)*CONJG(alm(3))*pi12lo(ipip,1)*100.D0/MULT
  acPz(num)=acPz(num)+alm(3)*CONJG(Clm(3))*pi12lo(ipip,1)*100.D0/MULT
  cbPz(num)=cbPz(num)+Clm(3)*CONJG(blm(3))*pe12lo(ipip,1)*100.D0/MULT
  bcPz(num)=bcPz(num)+blm(3)*CONJG(Clm(3))*pe12lo(ipip,1)*100.D0/MULT
  RETURN                                                            
END SUBROUTINE P3SPLT
