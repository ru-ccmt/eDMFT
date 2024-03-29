SUBROUTINE  IFFLIM(IFFT,IFF)                                      
  !     SETS FFT PARAMETERS TO MULTIPLES OF 2, 3, AND 5                   
   10 IF((IFFT/2)*2.NE.IFFT) GO TO 20                                   
      IFFT=IFFT/2                                                       
      IF(IFFT.GT.1) GO TO 10                                            
      RETURN                                                            
   20 IF((IFFT/3)*3.NE.IFFT) GO TO 30                                   
      IFFT=IFFT/3                                                       
      IF(IFFT.GT.1) GO TO 20                                            
      RETURN                                                            
   30 IF((IFFT/5)*5.NE.IFFT) GO TO 40                                   
      IFFT=IFFT/5                                                       
      IF(IFFT.GT.1) GO TO 30                                            
      RETURN                                                            
   40 CONTINUE                                                          
      IFF=IFF+1                                                         
      IFFT=IFF                                                          
      GO TO 10                                                          
END SUBROUTINE IFFLIM
