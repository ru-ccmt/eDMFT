      SUBROUTINE F7SPLT (ALM,BLM,MULT,UENORM,num,coord)
!                                                                       
!     DECOMPOSITION OF F CHARGE                
!                                                                      
      USE param
      USE charf
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16    ALM,BLM,CSUM                                        
      DIMENSION  ALM((lmax2+1)*(lmax2+1)),bLM((lmax2+1)*(lmax2+1))
      CHARACTER*5 COORD                                                 
!---------------------------------------------------------------------  
!                                                                       
      SQRT2=SQRT(2.D0)                                                  
      Af00(num)=Af00(num)+ALM(13)*CONJG(ALM(13))*100.D0/MULT                     
      Bf00(num)=Bf00(num)+BLM(13)*CONJG(BLM(13))*UENORM*100.D0/MULT              
!                                                                       
      CSUM=(ALM(14)+ALM(12))/SQRT2                                        
      Af11(num)=Af11(num)+CSUM*CONJG(CSUM)*100.D0/MULT                      
      CSUM=(BLM(14)+BLM(12))/SQRT2                                        
      Bf11(num)=Bf11(num)+CSUM*CONJG(CSUM)*UENORM*100.D0/MULT               
!                                                                       
      CSUM=(ALM(14)-ALM(12))/SQRT2                                        
      Af1m(num)=Af1m(num)+CSUM*CONJG(CSUM)*100.D0/MULT                          
      CSUM=(BLM(14)-BLM(12))/SQRT2                                        
      Bf1m(num)=Bf1m(num)+CSUM*CONJG(CSUM)*UENORM*100.D0/MULT                   
!                                                                       
      CSUM=(ALM(15)-ALM(11))/SQRT2                                        
      Af2m(num)=Af2m(num)+CSUM*CONJG(CSUM)*100.D0/MULT                          
      CSUM=(BLM(15)-BLM(11))/SQRT2                                        
      Bf2m(num)=Bf2m(num)+CSUM*CONJG(CSUM)*UENORM*100.D0/MULT                   
!                                                                       
      CSUM=(ALM(15)+ALM(11))/SQRT2                                        
      Af22(num)=Af22(num)+CSUM*CONJG(CSUM)*100.D0/MULT                          
      CSUM=(BLM(15)+BLM(11))/SQRT2                                        
      Bf22(num)=Bf22(num)+CSUM*CONJG(CSUM)*UENORM*100.D0/MULT                   
!                                                                       
      CSUM=(ALM(16)-ALM(10))/SQRT2                                        
      Af3m(num)=Af3m(num)+CSUM*CONJG(CSUM)*100.D0/MULT                          
      CSUM=(BLM(16)-BLM(10))/SQRT2                                        
      Bf3m(num)=Bf3m(num)+CSUM*CONJG(CSUM)*UENORM*100.D0/MULT                   
!                                                                       
      CSUM=(ALM(16)+ALM(10))/SQRT2                                        
      Af33(num)=Af33(num)+CSUM*CONJG(CSUM)*100.D0/MULT                          
      CSUM=(BLM(16)+BLM(10))/SQRT2                                        
      Bf33(num)=Bf33(num)+CSUM*CONJG(CSUM)*UENORM*100.D0/MULT                   
      RETURN                                                            
      END                                                               
