subroutine atpar(e,elo,elor,elor2,vru,zz,ilo,nt,jatom,radf,p,dp,pe,dpe,pei,isi)
  USE param, ONLY: lmax, nato, nrad, labc, nloat, clight, lomax
  USE loabc, ONLY: plo, dplo, pi12lo, pe12lo, alo, blo, clo !, elo
  USE loabcr, ONLY: plor, dplor, pi2lor, pe2lor, alor, blor, clor !, elor, elor2
  USE lolog, ONLY: mrf, lapw!, ilo
  USE rlolog, ONLY: loorext
  USE rpars, ONLY: extde, extei
  USE structure, ONLY: dx, jri, Rmt
  USE mpi, ONLY: Qprint
  IMPLICIT NONE
  REAL*8, intent(out) :: p(0:labc),dp(0:labc),pe(0:labc),dpe(0:labc),pei(0:labc)
  REAL*8, intent(in)  :: e(0:lmax), zz
  REAL*8, intent(inout):: elo(0:lomax,nloat), elor(0:lomax), elor2(0:lomax)
  INTEGER, intent(in) :: jatom, isi, nt, ilo(0:lomax)
  interface
     Function RINT13(REL,A,B,X,Y,NRAD,DX_,JRI_,R0_) result(S) ! Calculates overlap between psi_1=(A,B) and psi_2=(X,Y) functions
       REAL*8 :: S
       LOGICAL, intent(in) :: REL   ! relativistic or not
       REAL*8, intent(in)  :: A(NRAD),B(NRAD),X(NRAD),Y(NRAD)
       INTEGER, intent(in) :: NRAD, JRI_
       REAL*8, intent(in)  :: DX_, R0_
     End Function RINT13
     Function RINT13g(C1,C2,A,B,X,Y,NRAD,DX_,JRI_,R0_) result(S)
       ! instead of (C1,C2) = (1,1/137.0359895d0**2) we allow arbitrary C1 and C2 
       REAL*8 :: S
       REAL*8, intent(in) :: C1, C2
       REAL*8, intent(in) :: A(NRAD), B(NRAD), X(NRAD), Y(NRAD)
       INTEGER, intent(in):: NRAD, JRI_
       REAL*8, intent(in) :: DX_, R0_
     End Function RINT13g
  end interface
  !
  !  Note that vru is  v*r in Rydberg unit on input.
  !  It is transformed to v*r in Hartree units.
  !rschmid
  REAL*8  :: vru(nrad,nato), vr(nrad)
  LOGICAL :: REL
  REAL*8  :: cross, de, delei, r22, r_m, rmt2
  REAL*8  :: uv, uvb, uve, duv, duvb, duve, dx_, EI, eparam, L_, ovlp, R0_, TRX
  INTEGER :: jlo, jri_, kappa, l, ll, m, nl, nodel, nodes, nodeu
  !..................................................................
  !   ATPAR calculates the solutions u(l) of the radial Schroedinger
  !   Equation, the energy derivatives ue(l), the radial non muffin-
  !   tin integrals uvu, uevu, and the Gaunt-coefficients.
  !   Spin-orbitless relativistic equations are used.
  !..................................................................
  !        Common blocks
  real*8  :: radf(nrad,0:labc,2,2,nloat)
  real*8  :: a(nrad),b(nrad),ap(nrad),bp(nrad),ae(nrad),be(nrad)
  common  /work/     a,b,ap,bp,ae,be
  save    /work/
  character*4 emain
  !        Data statements
  REAL*8, PARAMETER     :: DELE = 2.0D-3
  COMPLEX*16, PARAMETER :: IMAG = (0.0D+0,1.0D+0)
  !        initialize constants
  !        CFEIN1 and CFEIN2 are used in RINT13 to calculate the
  !        integrals of radialfunctions and are derivates of the
  !        Feinstruktur-constant (in Hartrees)
  !
  REL=.True. ! relativistic nby default
  !
  DX_ = DX(JATOM)
  Jri_ = JRI(JATOM)
  RMT2=RMT(JATOM)
  R0_ = RMT2*DEXP(DX_*(1.D0-Jri_))
  !  The potential*r in Hartree units is calculated.
  if (Qprint) then
     WRITE(6,'(/10X,A,I2,A,I1)') 'ATOMIC PARAMETERS FOR ',jatom,'  SPIN=',isi
     WRITE(6,'(10X,A,7F7.2)') ' LINEARIZATION ENERGYS ARE', (E(ll),ll=0,6)
     WRITE(6,'(10X,A)') ' LINEARIZATION ENERGIES FOR LOs ARE'
     do nl=1,nloat-1
        write(6,'(10X,4F9.2)')(elo(ll,nl),ll=0,lomax)
     end do
     WRITE(6,14)
  end if
  
  ! Converting from Rydberg to Hartree
  vr(:)=0.d0
  vr(1:jri_)=vru(1:jri_,jatom)/2.d0
  delei = 0.25d0/dele
  DO L=0,NT-1
     !J=L+1
     mrf(l,jatom)=2
     L_ = L
     EI = e(l)/2.D0 ! el is in Ry, want Hartree's
     !     CALCULATE ENERGY-DERIVATIVE BY FINITE DIFFERENCE                  
     !     DELE IS THE UPWARD AND DOWNWARD ENERGY SHIFT IN HARTREES          
     !     CALCULATE FUNCTION AT EI-DELE
     CALL OUTWIN(A,B,NODEL,UVB,DUVB,REL,VR,R0_,DX_,JRI_,EI-DELE,L_,ZZ,NRAD)
     OVLP = RINT13(REL,A,B,A,B,NRAD,DX_,JRI_,R0_) ! OVLP = A**2+B**2
     TRX=1.0D0/SQRT(OVLP)                ! normalization=1/sqrt(A**2+B**2)
     AE(1:JRI_) = A(1:JRI_)*TRX
     BE(1:JRI_) = B(1:JRI_)*TRX
     UVB  = UVB*TRX
     DUVB = DUVB*TRX
     !     CALCULATE FUNCTION AT EI+DELE
     CALL OUTWIN(A,B,NODEU,UVE,DUVE, REL,VR,R0_,DX_,JRI_,EI+DELE,L_,ZZ,NRAD)
     OVLP = RINT13(REL,A,B,A,B,NRAD,DX_,JRI_,R0_)
     TRX=1.0d0/SQRT(OVLP)
     UVE=DELEI*(TRX*UVE-UVB)                                           
     DUVE=DELEI*(TRX*DUVE-DUVB)
     AE(1:JRI_)=DELEI*( A(1:JRI_)*TRX-AE(1:JRI_))
     BE(1:JRI_)=DELEI*( B(1:JRI_)*TRX-BE(1:JRI_))
     !     CALCULATE FUNCTION AT EI
     CALL OUTWIN(A,B,NODES,UV,DUV,REL,VR,R0_,DX_,JRI_,EI,L_,ZZ,NRAD)
     OVLP = RINT13(REL,A,B,A,B,NRAD,DX_,JRI_,R0_)
     TRX=1.0d0/SQRT(OVLP)
     
     P(l)=TRX*UV
     DP(l)=TRX*DUV
     A(1:JRI_) = A(1:JRI_)*TRX
     B(1:JRI_) = B(1:JRI_)*TRX
     !     INSURE ORTHOGONALIZATION                                          
     CROSS = RINT13(REL,A,B,AE,BE,NRAD,DX_,JRI_,R0_)
     AE(1:JRI_) = AE(1:JRI_)-CROSS*A(1:JRI_)
     BE(1:JRI_) = BE(1:JRI_)-CROSS*B(1:JRI_)                                            
     ! stores into : 
     !    RRAD1(r,l), RRAD2(r,l), RRADE1(r,l), RRADE2(r,l)
     !    P(l), DP(l), PE(l), DPE(l), PEI(l)
     RADF(1:JRI_,L,1,isi,1) = A(1:JRI_)
     RADF(1:JRI_,L,2,isi,1) = B(1:JRI_)
     RADF(1:JRI_,L,1,isi,2) = AE(1:JRI_)
     RADF(1:JRI_,L,2,isi,2) = BE(1:JRI_)
     PE(l) = UVE - CROSS*P(l)
     DPE(l) = DUVE - CROSS*DP(l)
     PEI(l) = RINT13(REL,AE,BE,AE,BE,NRAD,DX_,JRI_,R0_)
     if(l.le.6 .and. Qprint) then
        write(6,8) l,p(l),dp(l),pe(l),dpe(l),pei(l)
     end if
  enddo
  !
  !        and now for local orbitals
  !
  DO l=0,lomax
     do jlo=1,ilo(l)
        if ((.not.lapw(l,jatom)).and.jlo.eq.1) cycle
        mrf(l,jatom)=mrf(l,jatom)+1
        if (mrf(l,jatom).gt.nloat) then
           write(0,*) 'increase nloat in module  (atpar --- lo part)',mrf(l,jatom),nloat
           stop
        endif
        L_ = L
        EI = elo(l,jlo)/2.D0  ! From Rydbergs to Hartree's
        !     CALCULATE FUNCTION AT EI       
        CALL OUTWIN(A,B,NODES,UV,DUV,REL,VR,R0_,DX_,jri_,ei,L_,zz,NRAD)
        OVLP = RINT13(rel,a,b,a,b,NRAD,DX_,JRI_,R0_)
        ! TRX is used to normalize the wave functions.
        TRX=1.0d0/SQRT(OVLP)
        A(1:jri_) = A(1:jri_)*TRX
        B(1:jri_) = B(1:jri_)*TRX
        RADF(1:jri_,l,1,isi,mrf(l,jatom)) = A(1:JRI_)  ! this is called rf1 in other atpar
        RADF(1:jri_,l,2,isi,mrf(l,jatom)) = B(1:JRI_)   ! this is called rf2 in other atpar
        plo(l,jatom,isi,jlo)  = UV*trx
        dplo(l,jatom,isi,jlo) = DUV*trx
        pi12lo(l,jatom,isi,jlo) = RINT13(REL,RADF(:,l,1,isi,1),RADF(:,l,2,isi,1),A,B,NRAD,DX_,JRI_,R0_)
        pe12lo(l,jatom,isi,jlo) = RINT13(REL,RADF(:,l,1,isi,2),RADF(:,l,2,isi,2),A,B,NRAD,DX_,JRI_,R0_)
        if (Qprint) WRITE(6,'(10X,I2,2E14.6)') l,Plo(l,jatom,isi,jlo),DPlo(l,jatom,isi,jlo)
     enddo
  enddo
  
  DO l=0,lomax
     DO jlo=1,ilo(l)
        CALL abc0(alo(l,jlo,jatom,isi),blo(l,jlo,jatom,isi),clo(l,jlo,jatom,isi),lapw(l,jatom),P(l),dP(l),PE(l),dPE(l),PEI(l),PLO(l,jatom,isi,jlo),dPLO(l,jatom,isi,jlo),pi12lo(l,jatom,isi,jlo),pe12lo(l,jatom,isi,jlo),Rmt(jatom),l,jlo,Qprint)
     ENDDO
  ENDDO
  !    Now for the additional local orbital used in second
  !    diagonalisation.
  DO l=0,lomax
     IF (LOOREXT(L,JATOM)) THEN
        mrf(l,jatom)=mrf(l,jatom)+1
        if (mrf(l,jatom).gt.nloat) then
           write(0,*) 'increase nloat in module  (atpar --- rlo part)',mrf(l,jatom),nloat
           stop
        endif
        L_ = L   
        kappa = l
        !        calculate energy-derivative by finite difference
        !        DELE is the up and downward energy-shift in hartrees
        de = extde(jatom,l)
        emain = 'CONT'
        ! Test for RLO energy in  valence band
        ei = extei(jatom, l)
        CALL SELECT(L,ei,DE,EMAIN,REL,VR(1),R0_,DX_,Jri_,ZZ)
        elor(l) = ei
        elo(l,nloat-1)=ei

        if (Qprint) then
           write(6,*)'relativistic lo:',l,elor(l),ilo(l)+1,nloat
           write(8,*)'RELATIVISTIC LOs:'
           write(8,*)'on atom ',jatom,'e(',l,')=',elor(l)
           write(8,*)
        endif
        eparam = elor(l)
        ei = ei/2.d0
        CALL diracout_old(REL,VR,R0_,DX_,Jri_,EI,L,kappa,UV,DUV,NODES,ZZ)
        
        call dergl_old(a,b,r0_,dx_,jri_)
        
        do m = 1, jri_
           r_m = r0_ * exp(dx_*(m-1))
           b(m) =  b(m)*r_m/(2.d0*clight+(eparam-2.d0*vr(m)/r_m)/(2.d0*clight)) 
        enddo
        
        !CALL RINT13_old(1.d0,1.d0,A,B,A,B,OVLP,JATOM,R0_,DX_,Jri_)
        OVLP = RINT13g(1.d0,1.d0,A,B,A,B,NRAD,DX_,JRI_,R0_)
        
        TRX = 1.0D+0/SQRT(OVLP)
        PLOR(L,JATOM,isi) = TRX*UV
        DPLOR(L,JATOM,isi) = TRX*DUV
        
        A(1:jri_) = TRX*A(1:jri_)
        B(1:jri_) = TRX*B(1:jri_)

        !CALL RINT13_old(1.d0,1.d0/clight,RADF(:,L,1,isi,1),RADF(:,L,2,isi,1),a,b,PI2LOR(L,JATOM,isi),JATOM,R0_,DX_,Jri_)
        !CALL RINT13_old(1.d0,1.d0/clight,RADF(:,L,1,isi,2),RADF(:,L,2,isi,2),a,b,PE2LOR(L,JATOM,isi),JATOM,R0_,DX_,Jri_)
        PI2LOR(L,JATOM,isi) = RINT13g(1.d0,1.d0/clight,RADF(:,L,1,isi,1),RADF(:,L,2,isi,1),a,b,NRAD,DX_,JRI_,R0_)
        PE2LOR(L,JATOM,isi) = RINT13g(1.d0,1.d0/clight,RADF(:,L,1,isi,2),RADF(:,L,2,isi,2),a,b,NRAD,DX_,JRI_,R0_)
        
        
        do jlo=1,ilo(l)
           if ((.not.lapw(l,jatom)).and.jlo.eq.1) cycle       
           !call RINT13_old(1.d0,1.d0/clight,RADF(1,l,1,isi,jlo+2),RADF(1,l,2,isi,jlo+2),a,b,r22,jatom,r0_,dx_,jri_)
           r22 = RINT13g(1.d0,1.d0/clight,RADF(1,l,1,isi,jlo+2),RADF(1,l,2,isi,jlo+2),a,b,NRAD,DX_,JRI_,R0_)
        enddo

        do m = 1, jri_
           radf(m,l,1,isi,mrf(l,jatom)) =a(m)
           radf(m,l,2,isi,mrf(l,jatom)) =b(m)
        enddo
        
        elor2(l) = elor(l)
        if (Qprint) write(6,*)p(l)
        call abc_r(l,alor(l,jatom,isi),blor(l,jatom,isi),clor(l,jatom,isi),p(l),dp(l),pe(l),dpe(l),pei(l),plor(l,jatom,isi),dplor(l,jatom,isi),pi2lor(l,jatom,isi),pe2lor(l,jatom,isi),rmt(jatom),lapw(l,jatom))
     else
        elor(l) = 1.d+6
     ENDIF
     
  enddo
  !        end of loop over atoms in unitcell
  RETURN
14 FORMAT(/11X,1HL,5X,4HU(R),10X,5HU'(R),9X,5HDU/DE,8X,6HDU'/DE,6X,7HNORM-U')
8 FORMAT(10X,I2,5E21.10,5X,3I2)
END subroutine atpar
