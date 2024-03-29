!qprint == myrank.EQ.master .OR. fastFilesystem
SUBROUTINE abc(alo,blo,clo,l,jlo,lapw,p,dp,pe,dpe,pei,plo,dplo,pi12lo,pe12lo, nrad,lomax,nloat,lmax2,rmt_,qprint_)
  !
  ! Ouput: alo(l,jlo), blo(l,jlo), clo(l,jlo)
  !                                                                       
  IMPLICIT NONE
  REAL*8, intent(out) :: alo(0:lomax,nloat), blo(0:lomax,nloat), clo(0:lomax,nloat)
  INTEGER, intent(in) :: l,jlo
  LOGICAL, intent(in) :: lapw(0:lmax2)
  REAL*8, intent(in)  :: p(0:lmax2),dp(0:lmax2),pe(0:lmax2),DPE(0:lmax2),pei(0:lmax2)
  REAL*8, intent(in)  :: plo(nloat,0:lomax),dplo(nloat,0:lomax),pi12lo(nloat,0:lomax),pe12lo(nloat,0:lomax)
  INTEGER, intent(in) :: nrad, lomax, nloat, lmax2 
  REAL*8, intent(in)  :: rmt_
  LOGICAL, intent(in) :: qprint_
  ! locals
  REAL*8              :: xac,xbc,alonorm
  !---------------------------------------------------------------------  
  !                                                                       
  if(lapw(l)) then
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
     xac=(plo(jlo,l)*dpe(l)-dplo(jlo,l)*pe(l))*rmt_**2
     xbc=-(plo(jlo,l)*dp(l)-dplo(jlo,l)*p(l))*rmt_**2
     clo(l,jlo)=1.0d0/sqrt(1.d0 + xac*(xac+2.0D0*pi12lo(jlo,l))+xbc*(xbc*pei(l)+2.0D0*pe12lo(jlo,l)))
     if (clo(l,jlo).gt.2.0D2) clo(l,jlo)=2.0d2
     alo(l,jlo)=clo(l,jlo)*xac
     blo(l,jlo)=clo(l,jlo)*xbc
  else
     if(jlo.eq.1) then  ! APW+lo, hence ony u and dotu are used
        ! We construct the lo orbtial as   u_new = ALO*u + BLO*dotu
        !   and require       u_new(R) = 0
        !        and     <u_new|u_new> = 1
        ! which leads to:
        ! ALO = 1/sqrt(1 + (u/dotu)**2 * <dotu|dotu>)
        ! BLO = -u/dotu/sqrt(1 + (u/dotu)**2 * <dotu|dotu>)
        ! CLO = 0.d0
        alonorm=sqrt(1.d0+(P(l)/PE(l))**2*PEI(l))
        ALO(l,jlo) = 1.d0 /alonorm 
        BLO(l,jlo) = -P(l)/PE(l)/alonorm
        CLO(l,jlo) = 0.d0
     else ! apw+LO
        ! We construct the LO orbital as u_new = ALO*u + CLO*u_LO
        !   and require       u_new(R) = 0
        !       and      <u_new|u_new> = 1
        ! which leads to:
        ! xac = sqrt(1 + (u/u_LO)**2 - 2*(u/u_LO)*<u|u_LO>)
        ! ALO = 1/sqrt(1 + (u/u_LO)**2 - 2*(u/u_LO)*<u|u_LO>)
        ! CLO = -u/u_LO /sqrt(1 + (u/u_LO)**2 - 2*(u/u_LO)*<u|u_LO>)
        xbc=-P(l)/PLO(jlo,l)
        xac=sqrt(1+xbc**2+2*xbc*PI12LO(jlo,l))
        !         xac=1.d0
        ALO(l,jlo) = 1.d0/xac 
        BLO(l,jlo) = 0.d0
        CLO(l,jlo) = xbc/xac
     endif
  endif
  if (qprint_) WRITE(6,10) l,alo(l,jlo),blo(l,jlo),clo(l,jlo)
10 FORMAT ('LO COEFFICIENT: l,A,B,C  ',i2,5X,3F12.5)
  return
END SUBROUTINE abc

SUBROUTINE abc0(alo,blo,clo,l,jlo,lapw,p,dp,pe,dpe,pei,plo,dplo,pi12lo,pe12lo, nrad,lomax,nloat,lmax2,rmt_,qprint_)
  !
  ! Ouput: alo(l,jlo), blo(l,jlo), clo(l,jlo)
  !                                                                       
  IMPLICIT NONE
  REAL*8, intent(out) :: alo, blo, clo
  INTEGER, intent(in) :: l,jlo
  LOGICAL, intent(in) :: lapw
  REAL*8, intent(in)  :: p,dp,pe,DPE,pei
  REAL*8, intent(in)  :: plo,dplo,pi12lo,pe12lo
  INTEGER, intent(in) :: nrad, lomax, nloat, lmax2 
  REAL*8, intent(in)  :: rmt_
  LOGICAL, intent(in) :: qprint_
  ! locals
  REAL*8              :: xac,xbc,alonorm
  !---------------------------------------------------------------------  
  !                                                                       
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
     xac=(plo*dpe-dplo*pe)*rmt_**2
     xbc=-(plo*dp-dplo*p)*rmt_**2
     clo = 1.0d0/sqrt(1.d0 + xac*(xac+2.0D0*pi12lo)+xbc*(xbc*pei+2.0D0*pe12lo))
     if (clo.gt.2.0D2) clo=2.0d2
     alo = clo*xac
     blo = clo*xbc
  else
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
        xac=sqrt(1+xbc**2+2*xbc*PI12LO)
        !         xac=1.d0
        ALO = 1.d0/xac 
        BLO = 0.d0
        CLO = xbc/xac
     endif
  endif
  if (qprint_) WRITE(6,10) l,alo,blo,clo
10 FORMAT ('LO COEFFICIENT: l,A,B,C  ',i2,5X,3F12.5)
  return
END SUBROUTINE abc0
