SUBROUTINE LOHNS(JATOM,MULT,I1,I2,L0,jlo)
  !
  use lolog, only : ILO
  use comi, only  : NAT
  use param
  IMPLICIT NONE
  !INCLUDE 'param.inc'
  !        Arguments
  INTEGER, intent(in) :: jatom, l0, jlo
  INTEGER, intent(in) :: MULT(NAT)
  INTEGER, intent(out) :: i1, i2
  !..................................................................
  !   LOHNS calculates indices of loop in hns
  !..................................................................
  ! Local Scalars
  INTEGER            I, L
  INTEGER            JLO1
  I2 = 1
  DO I=1,JATOM-1
     DO L=0,LOMAX
        do jlo1=1,ilo(l,i)
           I2 = I2 + (2*L+1)*MULT(I)
        enddo
     ENDDO
  ENDDO
  
  DO L=0,L0-1
     do jlo1=1,ilo(l,jatom)
        I2 = I2 + (2*L+1)*MULT(jatom)
     enddo
  ENDDO

  do jlo1=1,jlo
     I1 = I2
     I2 = I2 + (2*L+1)*MULT(jatom)
  enddo
  !
  I2 = I2 - 1
  !
  RETURN
END SUBROUTINE LOHNS
