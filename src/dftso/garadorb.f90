subroutine garadorb(e,vru,p,dp,pe,dpe,ri_orb,jspin,kpot,ipr)
  USE param, ONLY: lmax, nrad, nato, labc, clight, nloat
  USE lolog, ONLY: mrf, ilo
  USE orb, ONLY: natorb, iat, ll
  USE structure, ONLY: dx, jri, Rmt, ZZ
  USE peic, ONLY: pei
  USE loabc, ONLY: elo
  USE loabcr, ONLY: elor, elor2
  IMPLICIT NONE
  REAL*8, intent(in)   :: E(0:lmax,NATO,2)
  REAL*8, intent(in)   :: VRU(NRAD,NATO,jspin)
  REAL*8, intent(inout):: P(Labc+1,NATO,2),DP(Labc+1,NATO,2),PE(Labc+1,NATO,2),DPE(Labc+1,NATO,2)
  REAL*8, intent(out)  :: ri_orb(4,4,0:labc,nato,2,2)
  INTEGER, intent(in)  :: jspin, kpot, ipr
  !
  ! subroutine based on garadme, calculate radial elements of orb potential
  ! in this version vorb(r)=1 is assumed
  ! local orbitals included, but p1/2 probably is not correct
  REAL*8, allocatable:: radf(:,:,:,:,:)
  REAL*8 :: rx(nrad)
  REAL*8 :: uu(nrad),tota(2)
  INTEGER :: i,ic, index, iradf, jradf
  INTEGER :: isj,isi
  INTEGER :: jric,jatom
  INTEGER :: ir,l
  REAL*8 :: coe, cin
  REAL*8 :: rmt2,dx2
  !
  !    coc  : velocity of light in Rydberg units
  !    cin  : inverse velocity of light squared. (Hartree units)
  COE=1.D0
  CIN=1.d0/CLIGHT**2
  
  allocate(radf(nrad,0:labc,2,2,nloat))
  
  ri_orb(:,:,:,:,:,:)=0.0d0

  !  jatom runs over all non-equivalent atoms.
  DO index=1,natorb 
     jatom=iat(index)
     DX2 = DX(JATOM)
     JRIc = JRI(JATOM)
     RMT2= RMT(JATOM)
     !  The logarithmic mesh is setup in RX for the respective atom.
     do i = 1, jric
        rx(i) = rmt2 * dexp(dfloat(i-jric)*dx2)
     enddo
     
     do isi=1,jspin
        call atpar(e(:,jatom,isi),elo(:,:,jatom,isi),elor(:,jatom,isi),elor2(:,jatom,isi),vru(:,:,jspin),zz(jatom),ilo(:,jatom),labc+1,jatom,radf,p(:,jatom,isi),dp(:,jatom,isi),pe(:,jatom,isi),dpe(:,jatom,isi),pei(:,jatom,isi),isi)
     enddo
     
     do isi=1,jspin
        do isj=1,jspin
           !.........................................................................
           !.....CALCULATE THE MATRIX ELEMENT
           !
           L=ll(index)
           do iradf=1,mrf(l,jatom)
              do jradf=1,mrf(l,jatom)   
                 DO IC=1,2
                    DO IR=1,JRIc
                       uu(ir)=radf(ir,l,ic,isi,iradf)*radf(ir,l,ic,isj,jradf)                           
                    enddo
                    call cali(uu,tota(ic),rx(1),dx2,jric)
                 enddo
                 ri_orb(iradf,jradf,l,jatom,isi,isj)=coe*(tota(1)+cin*tota(2))
              enddo
           enddo
        enddo
     enddo
  enddo
  deallocate(radf)
end subroutine garadorb
