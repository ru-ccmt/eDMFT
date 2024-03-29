SUBROUTINE CBCOMB(JRI,LLMM,VLM,lm,jatom,nat,nrad,lmmx)
  !
  !use comi, only : NAT
  !use param, only: 
  IMPLICIT NONE
  !INCLUDE 'param.inc'
  !        Arguments
  INTEGER, intent(in)  :: jri, LLMM, LM(2,LMMX,NAT), jatom, nat,nrad, lmmx
  REAL*8, intent(inout):: VLM(NRAD,LMMX)
  !..................................................................
  !   CBCOMB combines the radial total potential coefficients  of
  !   the cubic harmonics in the following manner:
  !      VLM40 = (VLM40*SQRT(7/12) + VLM44*SQRT(5/12))*SQRT(7/12)
  !      VLM44 = (VLM40*SQRT(7/12) + VLM44*SQRT(5/12))*SQRT(5/12/2)
  !      VLM60 = (VLM60*SQRT(2)/4  - VLM64*SQRT(14)/4)*SQRT(2)/4
  !      VLM64 = (VLM60*SQRT(2)/4  - VLM64*SQRT(14)/4)*-SQRT(14/2)/4
  !      VLM32 =  VLM32/SQRT(2)
  !..................................................................
  !   Local Scalars
  INTEGER            I, J
  REAL*8 c_kub(0:10,0:10)
  DOUBLE PRECISION ::       SQ1, SQRT2, C1, C2, C3
  !        Intrinsic Functions
  INTRINSIC          SQRT
  !
  c_kub(0,0)=1.d0
  c_kub(3,2)=1.d0
  c_kub(4,0)=.5d0*SQRT(7.d0/3.d0)
  c_kub(4,4)=.5*SQRT(5.d0/3.d0)
  c_kub(6,0)=.5d0*SQRT(.5d0)
  c_kub(6,2)=.25d0*SQRT(11.d0)
  c_kub(6,4)=-.5d0*SQRT(7.d0/2.d0)
  c_kub(6,6)=-.25d0*SQRT(5.d0)
  c_kub(7,2)=.5d0*SQRT(13.d0/6.d0)
  c_kub(7,6)=.5d0*SQRT(11.d0/6.d0)
  c_kub(8,0)=.125d0*SQRT(33.d0)
  c_kub(8,4)=.25d0*SQRT(7.d0/3.d0)
  c_kub(8,8)=.125d0*SQRT(65.d0/3.d0)
  c_kub(9,2)=.25d0*SQRT(3.d0)
  c_kub(9,4)=.5d0*SQRT(17.d0/6.d0)
  c_kub(9,6)=-.25d0*SQRT(13.d0)
  c_kub(9,8)=-.5d0*SQRT(7.d0/6.d0)
  c_kub(10,0)=.125d0*SQRT(65.D0/6.D0)
  c_kub(10,2)=.125d0*SQRT(247.D0/6.D0)
  c_kub(10,4)=-.25d0*SQRT(11.D0/2.D0)
  c_kub(10,6)=0.0625d0*SQRT(19.D0/3.D0)
  c_kub(10,8)=-.125d0*SQRT(187.D0/6.D0)
  c_kub(10,10)=-.0625d0*SQRT(85.d0)
  
  sqrt2=SQRT(2.d0)

  i=1
  DO
     IF(i.gt.llmm) exit
     !         print*,lm(1,i,jatom),lm(2,i,jatom)
     IF(lm(1,i,jatom).EQ.0.AND.lm(2,i,jatom).EQ.0) THEN
        i=i+1
     ELSEIF (lm(1,i,jatom).EQ.-3.AND.lm(2,i,jatom).EQ.2) THEN  
        DO J = 1, JRI
           VLM(J,i) = VLM(J,i)/SQRT2
        ENDDO
        i=i+1
     ELSEIF (lm(1,i,jatom).EQ.4.OR.lm(1,i,jatom).EQ.6.OR.lm(1,i,jatom).EQ.-7.OR.lm(1,i,jatom).EQ.-9) THEN
        IF (lm(2,i,jatom).EQ.0) THEN
           sq1=1.d0
        ELSE
           sq1=sqrt2
        ENDIF
        c1=c_kub(ABS(lm(1,i,jatom)),lm(2,i,jatom))
        c2=c_kub(ABS(lm(1,i,jatom)),lm(2,i,jatom)+4)
        DO J = 1, JRI
           VLM(J,i) = VLM(J,i)*C1 + VLM(J,i+1)*C2
           VLM(J,i+1) = VLM(J,i)*C2/SQRT2
           VLM(J,i) = VLM(J,i)*C1/SQ1
        ENDDO
        i=i+2
     ELSEIF (lm(1,i,jatom).EQ.8.OR.lm(1,i,jatom).EQ.10) THEN 
        IF (lm(2,i,jatom).EQ.0) THEN
           sq1=1.d0
        ELSE
           sq1=sqrt2
        ENDIF
        c1=c_kub(ABS(lm(1,i,jatom)),lm(2,i,jatom))
        c2=c_kub(ABS(lm(1,i,jatom)),lm(2,i,jatom)+4)
        c3=c_kub(ABS(lm(1,i,jatom)),lm(2,i,jatom)+8)
        DO J = 1, JRI
           vlm(j,i) = vlm(j,i)*c1 + vlm(j,i+1)*c2 + vlm(j,i+2) *c3
           VLM(J,i+1) = VLM(J,i)*C2/sqrt2
           VLM(J,i+2) = VLM(J,i)*C3/sqrt2
           VLM(J,i) = VLM(J,i)*C1/sq1
        ENDDO
        i=i+3
     ELSE
        WRITE(6,*) 'UNCORRECT LM LIST FOR CUBIC STRUCTURE'
        STOP 'Incorrect LM in CBCOMB'
     ENDIF
  END DO
  RETURN
END SUBROUTINE CBCOMB
