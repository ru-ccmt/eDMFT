       subroutine rintuout(jatom,is)
       use intu
       include 'param.inc'
       IMPLICIT REAL*8 (A-H,O-Z)
!      COMMON /INTU/ Duu1(LMAX1,NATO,2),Duu2(LMAX1,NATO,2), &
!                    Duup1(LMAX1,NATO,2),Duup2(LMAX1,NATO,2), &
!                    Dupu1(LMAX1,NATO,2),Dupu2(LMAX1,NATO,2), &
!                    Dupup1(LMAX1,NATO,2),Dupup2(LMAX1,NATO,2), &
!                    Ruu(LMAX1,NATO,2),Ruup(LMAX1,NATO,2), &
!                    Rupu(LMAX1,NATO,2),Rupup(LMAX1,NATO,2)
!      COMMON /INTUL/Duul1(LOMAX1,NATO,2), Duul2(LOMAX1,NATO,2), &
!                    Dulup1(LOMAX1,NATO,2),Dulup2(LOMAX1,NATO,2), &
!                    Dupul1(LOMAX1,NATO,2),Dupul2(LOMAX1,NATO,2), &
!                    Dulu1(LOMAX1,NATO,2), Dulu2(LOMAX1,NATO,2), &
!                    Dulul1(LOMAX1,NATO,2),Dulul2(LOMAX1,NATO,2), &
!                    Ruul(LOMAX1,NATO,2),  Rulu(LOMAX1,NATO,2), &
!                    Rupul(LOMAX1,NATO,2), Rulup(LOMAX1,NATO,2), &
!                    Rulul(LOMAX1,NATO,2)
!       do 10 jatom=1,nat
       write(6,*) ' ATOMNUMBER: ',jatom
       write(6,*) &
       '   L   Duu1          Dupu1         Duup1         Dupup1'
       do 11 i=1,lmax
       write(6,999) i,Duu1(i,jatom,is),Dupu1(i,jatom,is), &
                     Duup1(i,jatom,is),Dupup1(i,jatom,is)   
 11    continue
       write(6,*) &
       '   L   Duu2          Dupu2         Duup2         Dupup2'
       do 12 i=1,lmax
       write(6,999) i,Duu2(i,jatom,is),Dupu2(i,jatom,is),  &
                     Duup2(i,jatom,is),Dupup2(i,jatom,is)
 12    continue
       write(6,*)  &
       '   L   Ruu           Rupu          Ruup          Rupup'
       do 13 i=1,lmax
       write(6,999) i,Ruu(i,jatom,is),Ruup(i,jatom,is),Rupu(i,jatom,is), &
                      Rupup(i,jatom,is)
 13    continue
!ad
!ad   now for the LO's
!ad
       write(6,*) &
       '   L   Duul1         Duul2        Dulup1        Dulup2'
      do 14 i=1,lomax
       write(6,999) i,Duul1(i,jatom,is),Dupul2(i,jatom,is), &
                     Dulup1(i,jatom,is),Dulup2(i,jatom,is)
 14    continue
       write(6,*) &
       '   L   Dupul1        Dupul2       Dulu1         Dulu2'
      do 15 i=1,lomax
       write(6,999) i,Dupul1(i,jatom,is),Dupul2(i,jatom,is), &
                     Dulu1(i,jatom,is),Dulu2(i,jatom,is)
 15    continue
       write(6,*) &
       '   L   Dulul1        Dulul2       Ruul          Rulu '
      do 16 i=1,lomax
       write(6,999) i,Dulul1(i,jatom,is),Dulul2(i,jatom,is), &
                     Ruul(i,jatom,is),Rulu(i,jatom,is)
 16    continue
       write(6,*) &
       '   L   Rupul         Rulup        Rulul '              
      do 17 i=1,lomax
       write(6,999) i,Rupul(i,jatom,is),Rulup(i,jatom,is), &
                     Rulul(i,jatom,is)
 17    continue
 10    continue
       return
 999   format(i4,4(2x,e12.5))
       end

