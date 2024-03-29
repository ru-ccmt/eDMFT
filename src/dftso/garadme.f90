subroutine garadme(e,vru,p,dp,pe,dpe,ri_mat,jspin,kpot,ipr)
  USE param,  ONLY: lmax, nrad, nato, labc, nloat, clight, lomax
  USE loabc,  ONLY: Elo
  USE loabcr, ONLY: Elor2, pi2lor, pe2lor, elor, elor2
  USE lolog,  ONLY: mrf, lapw, ilo
  USE rlolog, ONLY: nrlo, loorext
  USE structure, ONLY: dx, jri, Rmt, ZZ
  USE hexpt,  ONLY: hexl, allocate_hexpt
  USE peic,   ONLY: pei, allocate_peic
  USE radovlp,ONLY: ruu, allocate_radovlp
  IMPLICIT NONE
  REAL*8, intent(in)   :: E(0:lmax,NATO,2)
  REAL*8, intent(in)   :: VRU(NRAD,NATO,jspin)
  REAL*8, intent(inout):: P(Labc+1,NATO,2), DP(Labc+1,NATO,2), PE(Labc+1,NATO,2), DPE(Labc+1,NATO,2)
  REAL*8, intent(out)  :: ri_mat(nloat,nloat,0:labc,nato,2,2)
  INTEGER, intent(in)  :: jspin, kpot, ipr
  ! external functions
  REAL*8 :: hscalc
  !
  !   The overlapp matrix and the radial matrix elements of the spin-orbit
  !   term as given in Koelling, Harmon (J. Phys. C (1976) ????) are calculated.
  !   A basis extension for an improved representation of p_{1/2} has been included.  
  !
  !   Input:
  !        e
  !        vru          potential * r in Rydberg.  
  !        p
  !        dp
  !        dpe
  !        zz(nato)     nuclear charge of the different atoms
  !        jspin        = 1 for unspinpolarized calculations
  !                     = 2 for spin-polarized calculations
  !         (Note that the present version is not intended for spin-polarized calculations,
  !          as some of the programming for this case might be wrong. )
  !        kpot           
  !        ipr          = debugging option
  !
  !   Output
  !        ri_mat      matrix containing the radial overlap integrals.
  !
  real*8, allocatable:: radf(:,:,:,:,:)
  REAL*8 :: vrr(nrad,2),vder(nrad,2),rx(nrad)
  REAL*8 :: ame(nrad), uu(nrad), tota(3)
  REAL*8 :: enei, enej, TOT
  INTEGER:: iradf, irfi, irfj, jlo, jradf
  integer i,ic
  integer isj,isi
  integer jric,jatom
  integer ir,l,lpot
  integer lspin
  
  real*8 coc,coe,cin
  real*8 rmt2,dx2
  real*8 dnomi,dnomj,dnom
  real*8 vde
  
  allocate(radf(nrad,0:labc,2,2,nloat))
  
  CALL allocate_peic(labc,nato)
  CALL allocate_hexpt(lomax,nato)
  CALL allocate_radovlp(labc,nato,nloat)
  !    coc  : velocity of light in Rydberg units
  !    cin  : inverse velocity of light squared. (Hartree units)
  COC=2.d0*CLIGHT
  COE=2.D0
  CIN=1.d0/CLIGHT**2
  
  ri_mat=0.0d0
  !  jatom is a loop over all non-equivalent atoms.
  DO jatom=1,nato
     DX2  = DX(jatom)
     JRIc = JRI(jatom)
     RMT2= RMT(jatom)
     !  The logarithmic mesh is setup in RX for the respective atom.
     do i = 1, jric
        rx(i) = rmt2 * dexp(dfloat(i-jric)*dx2)
     enddo
     ! if lpot=3, averaged potential is used in spinpolarized calculations
     ! when calculating dV/dr (kpot=1, jspin=2)
     lpot=jspin+kpot
     if (lpot.eq.3) then
        lspin=1
     else
        lspin=jspin
     endif
     !  Setup potential ( V*r -> V ), do averaging and calculate
     !  derivative of potential
     call vderiv(lpot,jatom,lspin,rx,vru,vrr,vder,jric,jspin)
     
     do isi=1,jspin
        call atpar(e(:,jatom,isi),elo(:,:,jatom,isi),elor(:,jatom,isi),elor2(:,jatom,isi),vru(:,:,jspin),ZZ(jatom),ilo(:,jatom),labc+1,jatom,radf,p(:,jatom,isi),dp(:,jatom,isi),pe(:,jatom,isi),dpe(:,jatom,isi),pei(:,jatom,isi),isi)
        ! contribution of V(2,m) to <p1/2|H|p1/2>
        if (nrlo(jatom).gt.0) call vnsrint(jatom,isi,radf,mrf(1,jatom))
     enddo
     !  Calculate some radial integrals needed for the construction
     !  of the overlap matrix.
     do isi=1,jspin
        do isj = 1,jspin
           do l=0,labc
              do iradf=1,mrf(l,jatom)
                 do jradf=1,iradf
                    do ic=1,2
                       do ir=1,jric
                          uu(ir)=radf(ir,l,ic,isi,iradf)*radf(ir,l,ic,isj,jradf)
                       enddo
                       call cali(uu,tota(ic),rx(1),dx2,jric)
                    enddo
                    !   Note that ic = 2 is not the whole expression for
                    !   the small component. If eq. 10 of Koelling, Harmon
                    !   is used for the overlapp matrix and S_l is discarded
                    !   an additional angular momentum dependent contribution
                    !   has to be added.
                    !   This contribution amouts to approx. 0.05 mRyd for p states
                    !   in Au. (It has not been taken into account in other parts
                    !   of the code.) For excited states there are contributions
                    !   up to 0.1 mRyd.  
                    if (l .le. lomax) then
                       if (loorext(l,jatom)) then 
                          if (iradf.lt.mrf(l,jatom).and.jradf.lt.mrf(l,jatom)) then
                             ruu(l,jatom,isi,isj,iradf,jradf)=tota(1)+(cin)*tota(2)
                          else
                             ruu(l,jatom,isi,isj,iradf,jradf)=tota(1)+sqrt(cin)*tota(2)
                          endif
                          if (iradf.eq.mrf(l,jatom).and.jradf.eq.mrf(l,jatom)) then 
                             ruu(l,jatom,isi,isj,iradf,jradf)=tota(1)+tota(2)  
                          endif
                          if (iradf.eq.mrf(l,jatom).and.jradf.lt.3) then 
                             ruu(l,jatom,isi,isj,iradf,jradf)=tota(1)+sqrt(cin)*tota(2)  
                          endif
                       else
                          ruu(l,jatom,isi,isj,iradf,jradf)=tota(1)+(cin)*tota(2)
                       endif
                    else
                       ruu(l,jatom,isi,isj,iradf,jradf)=tota(1)+(cin)*tota(2)
                    endif
                    ruu(l,jatom,isi,isj,jradf,iradf)=ruu(l,jatom,isi,isj,iradf,jradf)
                    !write(300,'(6i4,f15.8)') l,jatom,isi,isj,iradf,jradf,ruu(l,jatom,isi,isj,iradf,jradf) 
                 enddo
              enddo
           enddo
        enddo
     enddo
     
     do isi=1,jspin
        do isj=1,jspin
           !.........................................................................
           !.....CALCULATE THE MATRIX ELEMENT
           !
           DO L=0,LABC
              ic=1
              do irfi=1,mrf(l,jatom)
                 do irfj=1,mrf(l,jatom)
                    if (irfi.le.2) then
                       enei=E(L,jatom,isi)
                    else
                       if (loorext(l,jatom).and.irfi.eq.mrf(l,jatom)) then
                          !p1/2 local
                          enei=elor2(l,jatom,isi)                              
                       else
                          if (lapw(l,jatom)) then
                             jlo=irfi-2
                          else
                             jlo=irfi-1
                          endif
                          enei=ELO(L,jlo,jatom,isi)
                       endif
                    endif
                    
                    if (irfj.le.2) then
                       enej=E(L,jatom,isj)
                    else
                       if (loorext(l,jatom).and.irfj.eq.mrf(l,jatom)) then
                          !p1/2 local
                          enej=elor2(l,jatom,isj)                              
                       else
                          if (lapw(l,jatom)) then
                             jlo=irfj-2
                          else
                             jlo=irfj-1
                          endif
                          enej=ELO(L,jlo,jatom,isj)
                       endif
                    endif
                    
                    DO IR=1,JRIc
                       DNOMi=(COC+(enei-VRR(IR,isi))/COC)
                       DNOMj=(COC+(enej-VRR(IR,isj))/COC)
                       dnom=dnomi*dnomj
                       
                       if (isi.eq.isj) then
                          vde=vder(ir,isi)
                       else
                          vde=(vder(ir,isi)+vder(ir,isj))/2.0d0
                       endif
                       
                       AME(IR)=RADF(IR,L,IC,isi,irfi)*RADF(IR,L,IC,isj,irfj)*vde/RX(IR)/DNOM
                    enddo
                    
                    CALL CALI(AME,TOT,RX,DX2,JRIc)                      
                    
                    ri_mat(irfi,irfj,L,jatom,isi,isj)=COE*TOT
                 enddo
              enddo
           enddo
           !
           !   Calculate some integrals required for computation of 
           !   < \phi_sc | H_dirac | \phi_sc >.
           do l=0,lomax
              if (isi.eq.isj.and.loorext(l,jatom)) then
                 !honza
                 hexl(l,jatom,isi,1)=e(l,jatom,isi)
                 hexl(l,jatom,isi,2)=1.d0
                 hexl(l,jatom,isi,3)=pi2lor(l,jatom,isi)*elor2(l,jatom,isi)+hscalc(radf(1,l,1,isi,1),radf(1,l,1,isi,mrf(l,jatom)),vrr(1,isi),elor2(l,jatom,isi),l,jric,dx2,rx(1))
                 hexl(l,jatom,isi,4)=0.d0
                 hexl(l,jatom,isi,5)=0.d0
                 hexl(l,jatom,isi,6)=pe2lor(l,jatom,isi)*elor2(l,jatom,isi)+hscalc(radf(1,l,1,isi,2),radf(1,l,1,isi,mrf(l,jatom)),vrr(1,isi),elor2(l,jatom,isi),l,jric,dx2,rx(1))
                 hexl(l,jatom,isi,7)=e(l,jatom,isi)*pi2lor(l,jatom,isi)
                 hexl(l,jatom,isi,8)=pi2lor(l,jatom,isi)
                 hexl(l,jatom,isi,9)=elor2(l,jatom,isi)+hscalc(radf(1,l,1,isi,mrf(l,jatom)),radf(1,l,1,isi,mrf(l,jatom)),vrr(1,isi),elor2(l,jatom,isi),l,jric,dx2,rx(1))
              end if
           enddo
!!! Only the upper part of the spin-matrix has yet been calculated.
        enddo
     enddo
  enddo
  deallocate(radf)
END subroutine garadme
