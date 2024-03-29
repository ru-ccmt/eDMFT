SUBROUTINE abc0(alo,blo,clo,lapw,p,dp,pe,dpe,pei,plo,dplo,pi12lo,pe12lo,rmt_,l,jlo,qprint_)
  IMPLICIT NONE
  REAL*8, intent(out) :: alo, blo, clo
  logical, intent(in) :: lapw
  REAL*8,  intent(in) :: p, dp, pe, dpe, pei, plo, dplo, pi12lo, pe12lo, rmt_
  INTEGER, intent(in) :: l, jlo
  LOGICAL, intent(in) :: qprint_
  !..................................................................
  !   ABC calculates the cofficients A,B,C of the local orbitals (LO)
  !..................................................................
  !        Local Parameters
  DOUBLE PRECISION   CUTOFF
  PARAMETER          (CUTOFF = 200.0D+0)
  !        Local Scalars
  DOUBLE PRECISION   XAC, XBC, ALONORM
  !        Intrinsic Functions
  
  if(lapw) then
     ! DPE are indexed from 0 to LMAX2 in lapw2 usw
     ! DPE are indexed from 1 to LMAX2+1 in lapw1 usw
     !   P(l)   = ul_Rmt(1,l,jatom)
     !   PE(l)  = ul_Rmt(2,l,jatom)
     !   DP(l)  = dul_Rmt(1,l,jatom)
     !   DPE(l) = dul_Rmt(2,l,jatom)
     !   PLO(l) = ul_Rmt(2+jlo,l,jatom)
     !   DPLO(l)= dul_Rmt(2+jlo,l,jatom)
     !   ri_mat(2,2,l,jatom)=pei(l) ! <udot|udot>
     !   ri_mat(1,2+jlo,l,jatom) = pi12lo(jlo,l)  ! <u | u_lo>
     !   ri_mat(2+jlo,1,l,jatom) = pi12lo(jlo,l)  ! <u_lo | u>
     !   ri_mat(2,2+jlo,l,jatom) = pe12lo(jlo,l)  ! <udot | u_lo>
     !   ri_mat(2+jlo,2,l,jatom) = pe12lo(jlo,l)  ! <u_lo | udot>
     !   ri_mat(2+jlo,2+jlop,l,jatom) = pr12lo(jlo,jlop,l)  ! <u_lo | u_lo >
     !
     !  We construct the LO orbtial as   u_new = ALO*u + BLO*dotu + CLO*u_LO
     !    and require       u_new(R) = 0
     !                  du_new/dr(R) = 0
     !         and     <u_new|u_new> = 1
     !  which leads to:
     ! xac =  (u_LO*ddotu-du_LO*dotu)*R^2
     ! xbc = -(u_LO*du - du_LO*u)*R^2
     ! clo =   1/sqrt( 1 + xac*(xac + 2*<u_LO|u>)+xbc*(xbc*<dotu|dotu>+2*<dotu|u_LO> )
     ! alo = xac/sqrt( 1 + xac*(xac + 2*<u_LO|u>)+xbc*(xbc*<dotu|dotu>+2*<dotu|u_LO> )
     ! blo = xbc/sqrt( 1 + xac*(xac + 2*<u_LO|u>)+xbc*(xbc*<dotu|dotu>+2*<dotu|u_LO> )
     xac = (plo*dpe - dplo*pe)*rmt_**2
     xbc = -(plo*dp - dplo*p)*rmt_**2
     clo = 1.0d0/sqrt(1.d0 + xac*(xac + 2.0D0*pi12lo) + xbc*(xbc*pei + 2.0D0*pe12lo))
     clo = min(clo,CUTOFF)
     alo = clo*xac
     blo = clo*xbc
  else
     !.....APW definitions
     if(jlo.eq.1) then  ! APW+lo, hence ony u and dotu are used
        ! We construct the lo orbtial as   u_new = ALO*u + BLO*dotu
        !   and require       u_new(R) = 0
        !        and     <u_new|u_new> = 1
        ! which leads to:
        ! ALO = 1/sqrt(1 + (u/dotu)**2 * <dotu|dotu>)
        ! BLO = -u/dotu/sqrt(1 + (u/dotu)**2 * <dotu|dotu>)
        ! CLO = 0.d0
        alonorm=sqrt(1.d0+(P/PE)**2*pei)
        alo = 1.d0 /alonorm 
        blo = -P/PE/alonorm
        clo = 0.d0
     else ! apw+LO
        ! We construct the LO orbital as u_new = ALO*u + CLO*u_LO
        !   and require       u_new(R) = 0
        !       and      <u_new|u_new> = 1
        ! which leads to:
        ! xac = sqrt(1 + (u/u_LO)**2 - 2*(u/u_LO)*<u|u_LO>)
        ! ALO = 1/sqrt(1 + (u/u_LO)**2 - 2*(u/u_LO)*<u|u_LO>)
        ! CLO = -u/u_LO /sqrt(1 + (u/u_LO)**2 - 2*(u/u_LO)*<u|u_LO>)
        xbc=-P/PLO
        xac=sqrt(1+xbc**2+2*xbc*pi12lo)
        alo = 1.d0/xac
        blo = 0.d0
        clo = xbc/xac
     endif
  endif
  if (qprint_) WRITE(6,6000) l,alo,blo,clo
  !
  RETURN
  !
6000 FORMAT ('LO COEFFICIENT: l,A,B,C  ',i2,5X,3F12.5)
END SUBROUTINE ABC0
