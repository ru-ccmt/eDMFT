REAL*8 FUNCTION GAUNT(L1,L2,L3,MM1,MM2,MM3)
  !IMPLICIT REAL*8 (A-H,O-Z)
  IMPLICIT NONE
  INTEGER, intent(in) :: L1,L2,L3,MM1,MM2,MM3
  !
  REAL*8  :: T3J0, T3J
  REAL*8  :: FAC, FPI
  INTEGER :: J1,J2,J3,M1,M2,M3
  DATA FPI/12.56637061436D0/                                        
  J1=2*L1                                                           
  J2=2*L2                                                           
  J3=2*L3                                                           
  M1=-2*MM1                                                         
  M2=2*MM2                                                          
  M3=2*MM3                                                          
  FAC=(-1.)**MM1*SQRT((J1+1)*(J2+1)*(J3+1)/FPI)                     
  GAUNT=FAC*T3J0(J1,J2,J3)*T3J(J1,J2,J3,M1,M2,M3)                   
  RETURN                                                            
END FUNCTION GAUNT
