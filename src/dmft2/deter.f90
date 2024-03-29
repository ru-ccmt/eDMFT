SUBROUTINE DETER(GMAX,pia,RBAS,kxmax,kymax,kzmax,lattic)              
  !     **                                                              **
  !     **  CALCULATE RECIPROCAL LATTICE VECTORS FROM REAL SPACE        **
  !     **  LATTICE VECTORS OR VICE VERSA                               **
  !     **                                                              **
  IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
  character*4 lattic
  DIMENSION RBAS(3,3),GBAS(3,3),pia(3)                                    
  GBAS(1,1)=RBAS(2,2)*RBAS(3,3)-RBAS(3,2)*RBAS(2,3)                 
  GBAS(2,1)=RBAS(3,2)*RBAS(1,3)-RBAS(1,2)*RBAS(3,3)                 
  GBAS(3,1)=RBAS(1,2)*RBAS(2,3)-RBAS(2,2)*RBAS(1,3)                 
  GBAS(1,2)=RBAS(2,3)*RBAS(3,1)-RBAS(3,3)*RBAS(2,1)                 
  GBAS(2,2)=RBAS(3,3)*RBAS(1,1)-RBAS(1,3)*RBAS(3,1)                 
  GBAS(3,2)=RBAS(1,3)*RBAS(2,1)-RBAS(2,3)*RBAS(1,1)                 
  GBAS(1,3)=RBAS(2,1)*RBAS(3,2)-RBAS(3,1)*RBAS(2,2)                 
  GBAS(2,3)=RBAS(3,1)*RBAS(1,2)-RBAS(1,1)*RBAS(3,2)                 
  GBAS(3,3)=RBAS(1,1)*RBAS(2,2)-RBAS(2,1)*RBAS(1,2)                 
  DET=0.D0                                                          
  DO I=1,3                                                      
     DET=DET+GBAS(I,1)*RBAS(I,1)                                       
  ENDDO
  kxmax0=DABS(GBAS(1,1)/DET*GMAX)
  kymax0=DABS(GBAS(1,2)/DET*GMAX)
  kzmax0=DABS(GBAS(1,3)/DET*GMAX)
  xtest=DABS(GBAS(2,1)/DET*GMAX)
  if (xtest.gt.kxmax0) kxmax0= xtest
  xtest=DABS(GBAS(3,1)/DET*GMAX)
  if (xtest.gt.kxmax0) kxmax0= xtest
  ytest=DABS(GBAS(2,2)/DET*GMAX)
  if (ytest.gt.kymax0) kymax0= ytest
  ytest=DABS(GBAS(3,2)/DET*GMAX)
  if (ytest.gt.kymax0) kymax0= ytest
  ztest=DABS(GBAS(2,3)/DET*GMAX)
  if (ztest.gt.kzmax0) kzmax0= ztest
  ztest=DABS(GBAS(3,3)/DET*GMAX)
  if (ztest.gt.kzmax0) kzmax0= ztest
  if(lattic(1:1).eq.'B') then
     kxmax1=gmax/sqrt(dble(pia(1)**2+pia(2)**2))
     kymax1=kxmax1
     kxmax2=gmax/sqrt(dble(pia(1)**2+pia(3)**2))
     kzmax1=kxmax2
     kymax2=gmax/sqrt(dble(pia(2)**2+pia(3)**2))
     kzmax2=kymax2
     kxmax=max(kxmax0,kxmax1,kxmax2)
     kymax=max(kymax0,kymax1,kymax2)
     kzmax=max(kzmax0,kzmax1,kzmax2)
     !          write(*,*) kxmax,kymax,kzmax,kxmax0,kymax0,kzmax0 &
     !         ,kxmax1,kymax1,kzmax1,kxmax2,kymax2,kzmax2
  else if(lattic(1:1).eq.'F') then
     kxmax1=gmax/sqrt(dble(pia(1)**2+pia(2)**2+pia(3)**2))
     kymax1=kxmax1
     kzmax1=kxmax1
     kxmax=max(kxmax0,kxmax1)
     kymax=max(kymax0,kymax1)
     kzmax=max(kzmax0,kzmax1)
     !          write(*,*) kxmax,kymax,kzmax,kxmax0,kymax0,kzmax0 &
     !         ,kxmax1,kymax1,kzmax1
  else if(lattic(1:3).eq.'CXY') then
     kxmax1=gmax/sqrt(dble(pia(1)**2+pia(2)**2))
     kymax1=kxmax1
     kzmax1=1
     kxmax=max(kxmax0,kxmax1)
     kymax=max(kymax0,kymax1)
     kzmax=max(kzmax0,kzmax1)
     !          write(*,*) kxmax,kymax,kzmax,kxmax0,kymax0,kzmax0 &
     !         ,kxmax1,kymax1,kzmax1
  else if(lattic(1:3).eq.'CXZ') then
     kxmax1=gmax/sqrt(dble(pia(1)**2+pia(3)**2))
     kzmax1=kxmax1
     kymax1=1
     kxmax=max(kxmax0,kxmax1)
     kymax=max(kymax0,kymax1)
     kzmax=max(kzmax0,kzmax1)
     !          write(*,*) kxmax,kymax,kzmax,kxmax0,kymax0,kzmax0 &
     !         ,kxmax1,kymax1,kzmax1
  else if(lattic(1:3).eq.'CYZ') then
     kzmax1=gmax/sqrt(dble(pia(2)**2+pia(3)**2))
     kymax1=kzmax1
     kxmax1=1
     kxmax=max(kxmax0,kxmax1)
     kymax=max(kymax0,kymax1)
     kzmax=max(kzmax0,kzmax1)
     !          write(6,*) kxmax,kymax,kzmax,kxmax0,kymax0,kzmax0 &
     !         ,kxmax1,kymax1,kzmax1
  else 
     kxmax=max(kxmax0,1)
     kymax=max(kymax0,1)
     kzmax=max(kzmax0,1)
  endif
  RETURN                                                            
END SUBROUTINE DETER
      
