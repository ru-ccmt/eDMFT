      subroutine abc (l,jlo,jatom)                           
      use struk
!                                                                       
!     abc calculates the cofficients a,b,c of the lo    
!                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)
! ---
      INCLUDE 'param.inc'
! ---
      integer l,jatom,jlo
      logical loor(0:lomax),lloor(0:lmax2)
      logical lapw(0:lmax2)
      common /loabc/   alo(0:lomax,nloat),blo(0:lomax,nloat), &
                       clo(0:lomax,nloat), &
                       elo(0:lomax,nloat),plo(0:lomax), &
                       dplo(0:lomax),pelo(0:lomax), &
                       dpelo(0:lomax),peilo(0:lomax),      &
                       pi12lo(0:lomax),pe12lo(0:lomax), &
                       a1lo(nrad,0:lomax),b1lo(nrad,0:lomax)
      common /lolog/   nlo,nlov,nlon,lapw,ilo(0:lomax),loor,lloor
!
      COMMON /ATSPDT/  E(0:LMAX2),P(0:LMAX2),DP(0:LMAX2),PE(0:LMAX2), &
                       DPE(0:LMAX2),PEI(0:LMAX2)              
!      COMMON /STRUK/   POS(3,NDIF),ALATT(3),ALPHA(3),RMT(NATO),V(NATO), &
!                       PIA(3),VOL,IATNR(NATO),MULT(NATO),ISPLIT(NATO)         
!---------------------------------------------------------------------  
!      
      if(lapw(l)) then
! DPE are indexed from 0 to LMAX2 in lapw2 usw
! DPE are indexed from 1 to LMAX2+1 in lapw1 usw
! therefore difference in indexing
         xac=plo(l)*dpe(l)-dplo(l)*pe(l)
         xac=xac*rmt(jatom)*rmt(jatom)
         xbc=plo(l)*dp(l)-dplo(l)*p(l)
         xbc=-xbc*rmt(jatom)*rmt(jatom)
         clo(l,jlo)=xac*(xac+2.0D0*pi12lo(l))+xbc* &
              (xbc*pei(l)+2.0D0*pe12lo(l))+1.0D0  
         clo(l,jlo)=1.0D0/sqrt(clo(l,jlo))
         if (clo(l,jlo).gt.2.0D2) clo(l,jlo)=2.0d2
         alo(l,jlo)=clo(l,jlo)*xac
         blo(l,jlo)=clo(l,jlo)*xbc
      else
      if(jlo.eq.1) then
         alonorm=sqrt(1.d0+(P(l)/PE(l))**2*PEI(l))
         ALO(l,jlo) = 1.d0 /alonorm 
         BLO(l,jlo) = -P(l)/PE(l)/alonorm
         CLO(l,jlo) = 0.d0
      else 
        xbc=-P(l)/PLO(l)
        xac=sqrt(1+xbc**2+2*xbc*PI12LO(l))
        !         xac=1.d0
        ALO(l,jlo) = 1.d0/xac 
        BLO(l,jlo) = 0.d0
        CLO(l,jlo) = xbc/xac
      endif
      endif
      write (6,10) l,alo(l,jlo),blo(l,jlo),clo(l,jlo)
 10   FORMAT ('LO COEFFICIENT: l,A,B,C  ',i2,5X,3F12.5)
      return
      end








