subroutine abc_r (l,alor_l,blor_l,clor_l,p_l,dp_l,pe_l,dpe_l,pei_l,plor_l,dplor_l,pi2lor_l,pe2lor_l,rmt,lapw)
  use mpi, ONLY: Qprint
  IMPLICIT NONE
  REAL*8, intent(out) :: alor_l, blor_l, clor_l
  INTEGER, intent(in) :: l
  REAL*8,  intent(in) :: Rmt
  logical, intent(in) :: lapw
  REAL*8, intent(in)  :: DP_l,DPE_l, P_l,PE_l,PEI_l
  REAL*8, intent(in)  :: plor_l, dplor_l, pi2lor_l, pe2lor_l
  ! locals
  REAL*8 :: CUTOFF
  PARAMETER          (CUTOFF = 200.0D+0)
  REAL*8 :: XAC, XBC
  INTRINSIC          SQRT

  if (lapw) then
     XAC =  PLOR_l*DPE_l - DPLOR_l* PE_l
     XAC = XAC*RMT*RMT
     XBC =  PLOR_l* DP_l - DPLOR_l*  P_l
     XBC = -XBC*RMT*RMT
     CLOR_l = XAC*(XAC + 2.0D0*PI2LOR_l) + XBC*(XBC*PEI_L + 2.0D0*PE2LOR_l) + 1.0D0
     CLOR_l = 1.0D0/SQRT(CLOR_l)
     CLOR_l = MIN(CLOR_l,CUTOFF)
     ALOR_l = CLOR_l*XAC
     BLOR_l = CLOR_l*XBC
  else
     xbc=-P_l/PLOR_l
     xac=sqrt(1+xbc**2+2*xbc*PI2LOR_l)
     ALOR_l = 1.d0/xac
     BLOR_l = 0.d0
     CLOR_l = xbc/xac
  end if
  if (Qprint) write (6,10)l,alor_l,blor_l,clor_l
  RETURN
  !
10 FORMAT ('RLO COEFFICIENT: l,A,B,C  ',i2,5X,3F12.5)
END subroutine abc_r
