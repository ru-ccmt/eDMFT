SUBROUTINE ATPARN(RRAD1,RADE1,RRAD2,RADE2,a1lo,b1lo,ul_Rmt,dul_Rmt,ri_mat,alo,blo,clo, Vr_input,nr,el,elo,loor,rlo,lapw,nlov,nlo,nlon,ilo,rel,ZZ,qprint_,jri_,r0_,dx_,rmt_,lmax2,lomax,nloat,nrad)
  !USE param,  ONLY: lmax2, lomax, nloat, nrad
  !use struk,  ONLY: jri, r0, dx, rmt ! V
  !use lo,     ONLY: plo, dplo, pi12lo, pe12lo, pr12lo, a1lo, b1lo, alo, blo, clo
  !USE atspdt, ONLY: p,dp,pe,dpe,pei
  !USE com_mpi
  IMPLICIT NONE
  ! Solves the radial scalar relativistic Dirac equation and returns the wave functions and their values at MT-sphere and their overlaps.
  ! Inputs are:
  !        Vr_input(nr)  -- spherical part of the KS radial potential in Rydbergs
  !        el, elo       -- linearization energies
  !        lapw          -- which atom/l is treated by LAPW versus APW+lo basis
  !        loor,rlo,nlov,nlo,nlon,ilo  -- information on local orbitals
  !        rel           -- relativistic or not
  !        qprint_       -- should we print extra info
  !        ZZ            -- nucleous charge
  !        jri_          -- numer of radial point for this atom
  !        r0_,dx_       -- radial mesh r0="first point", dx="log distance in points"
  !        rmt_          -- R_mt for this atom
  !        lmax2         -- maximum l for LAPW
  !        lomax         -- maximum l for local orbitals
  !        nloat         -- maximum number of local orbitals
  !        nrad          -- maximum number of radial points
  ! Output are:
  !  1) Radial wave functions:
  !        RRAD1(r,l), RRAD2(r,l), RADE1(r,l), RADE2(r,l)
  !  2)  Values of the wave functions at the MT-sphere:
  !        ul_Rmt(1,l)      = P(l) 
  !        ul_Rmt(2,l)      = PE(l)
  !        dul_Rmt(1,l)     = DP(l)
  !        dul_Rmt(2,l)     = DPE(l)
  !        ul_Rmt(2+jlo,l)  = PLO(jlo,l)
  !        dul_Rmt(2+jlo,l) = DPLO(jlo,l)
  !  3) Overlaps of the above wave functions:
  !        ri_mat(2,2,l)          = pei(l)             ! <udot|udot>
  !        ri_mat(1,2+jlo,l)      = pi12lo(jlo,l)      ! <u | u_lo>
  !        ri_mat(2,2+jlo,l)      = pe12lo(jlo,l)      ! <udot | u_lo>
  !        ri_mat(2+jlo,2+jlop,l) = pr12lo(jlo,jlop,l) ! <u_lo | u_lo >
  !  4) Coefficients for local orbitals:
  !        a1lo(nrad,nloat,0:lomax), b1lo(nrad,nloat,0:lomax)
  !        alo, blo, clo                                       
  !
  REAL*8, intent(out) :: RRAD1(nrad,0:lmax2),RADE1(nrad,0:lmax2),RRAD2(nrad,0:lmax2),RADE2(nrad,0:lmax2)
  REAL*8, intent(out) :: a1lo(nrad,nloat,0:lomax),b1lo(nrad,nloat,0:lomax)
  REAL*8, intent(out) :: ul_Rmt(nloat+2,0:lmax2), dul_Rmt(nloat+2,0:lmax2), ri_mat(2+nloat,2+nloat,0:lmax2)
  REAL*8, intent(out) :: alo(0:lomax,nloat), blo(0:lomax,nloat), clo(0:lomax,nloat)
  !
  REAL*8,  intent(in) :: Vr_input(nr)
  INTEGER, intent(in) :: nr
  REAL*8,  intent(in) :: el(0:lmax2), elo(0:lomax,1:nloat)
  LOGICAL, intent(in) :: loor(nloat,0:lomax), rlo(nloat,0:lomax), lapw(0:lmax2)
  INTEGER, intent(in) :: nlov, nlo, nlon
  INTEGER, intent(in) :: ilo(0:lmax2)
  LOGICAL, intent(in) :: rel
  REAL*8,  intent(in) :: ZZ
  LOGICAL, intent(in) :: qprint_
  INTEGER, intent(in) :: jri_
  REAL*8,  intent(in) :: r0_, dx_, rmt_
  INTEGER, intent(in) :: lmax2, lomax, nloat, nrad
  !CHARACTER*4, intent(in):: cform
  interface
     Function RINT13(REL,A,B,X,Y,NRAD,DX_,JRI_,R0_) result(S) ! Calculates overlap between psi_1=(A,B) and psi_2=(X,Y) functions
       REAL*8 :: S
       LOGICAL, intent(in) :: REL   ! relativistic or not
       REAL*8, intent(in)  :: A(NRAD),B(NRAD),X(NRAD),Y(NRAD)
       INTEGER, intent(in) :: NRAD, JRI_
       REAL*8, intent(in)  :: DX_, R0_
     End Function RINT13
  end interface
  ! locals
  REAL*8,PARAMETER  :: CLIGHT= 137.0359895d0
  REAL*8,PARAMETER  :: PI=     3.1415926535897932d0
  REAL*8 :: P(0:lmax2), PE(0:lmax2), DP(0:lmax2), DPE(0:lmax2), PLO(nloat,0:lomax), DPLO(nloat,0:lomax)
  REAL*8 :: pei(0:lmax2), pi12lo(nloat,0:lomax), pe12lo(nloat,0:lomax), pr12lo(nloat,nloat,0:lomax)
  REAL*8     :: A(NRAD), B(NRAD), AE(NRAD), BE(NRAD)
  real*8     :: VR(NRAD), dele, delei, L_, EI, E1, UVB, DUVB, OVLP, TRX, UVE, DUVE, UV, DUV, CROSS, TRY, r_m
  INTEGER    :: IDUMMY, I, J, K, L, M, JC, JR, IDEST, JROW, JCOL, INDEX, NODEL, IMAX, NODE, NODES, jlo, kappa, jlop, lp
  !---------------------------------------------------------------------  
  !
  VR(:)=0.0
  VR(:nr) = VR_input(:nr)/2.d0  ! Converting from Rydberg to Hartree

  DO l=0,lmax2
     DELE=2.0D-3
     DELEI=0.25D0/DELE
     L_=L
     EI=el(l)/2.0d0  ! el is in Ry, want Hartree's
     !     CALCULATE ENERGY-DERIVATIVE BY FINITE DIFFERENCE                  
     !     DELE IS THE UPWARD AND DOWNWARD ENERGY SHIFT IN HARTREES          
     ! outputs: UVB, DUVB, NODEL
     !     CALCULATE FUNCTION AT EI-DELE
     !print *, 'EI=', EI, 'L=', L_
     CALL OUTWIN(A,B,NODEL,UVB,DUVB,REL,VR,R0_,DX_,JRI_,EI-DELE,L_,ZZ,NRAD)
     OVLP = RINT13(REL,A,B,A,B,NRAD,DX_,JRI_,R0_) ! OVLP = A**2+B**2
     TRX=1.0D0/SQRT(OVLP)                ! normalization=1/sqrt(A**2+B**2)
     AE(1:JRI_) = A(1:JRI_)*TRX
     BE(1:JRI_) = B(1:JRI_)*TRX
     UVB  = UVB*TRX
     DUVB = DUVB*TRX
     ! results are: AE(:), BE(:), UVB, DUVB, E1
     ! outputs: UVE, DUVE, NODE, OVLP
     !     CALCULATE FUNCTION AT EI+DELE
     CALL OUTWIN(A,B,NODE,UVE,DUVE, REL,VR,R0_,DX_,JRI_,EI+DELE,L_,ZZ,NRAD)
     OVLP = RINT13(REL,A,B,A,B,NRAD,DX_,JRI_,R0_)
     TRX=1.0d0/SQRT(OVLP)
     UVE=DELEI*(TRX*UVE-UVB)                                           
     DUVE=DELEI*(TRX*DUVE-DUVB)
     AE(1:JRI_)=DELEI*( A(1:JRI_)*TRX-AE(1:JRI_))
     BE(1:JRI_)=DELEI*( B(1:JRI_)*TRX-BE(1:JRI_))
     !     CALCULATE FUNCTION AT EI
     ! results are: A(:), B(:), P(:), DP(:), UV, DUV
     CALL OUTWIN(A,B,NODES,UV,DUV,REL,VR,R0_,DX_,JRI_,EI,L_,ZZ,NRAD)
     OVLP = RINT13(REL,A,B,A,B,NRAD,DX_,JRI_,R0_)
     TRX=1.0d0/SQRT(OVLP)
     P(l)=TRX*UV
     DP(l)=TRX*DUV
     A(1:JRI_) = A(1:JRI_)*TRX
     B(1:JRI_) = B(1:JRI_)*TRX
     !                                                                       
     !     INSURE ORTHOGONALIZATION                                          
     !                                                                       
     CROSS = RINT13(REL,A,B,AE,BE,NRAD,DX_,JRI_,R0_)
     AE(1:JRI_) = AE(1:JRI_)-CROSS*A(1:JRI_)
     BE(1:JRI_) = BE(1:JRI_)-CROSS*B(1:JRI_)                                            
     ! stores into : 
     !    RRAD1(r,l), RRAD2(r,l), RRADE1(r,l), RRADE2(r,l)
     !    P(l), DP(l), PE(l), DPE(l), PEI(l)
     RRAD1(1:JRI_,l) = A(1:JRI_)
     RRAD2(1:JRI_,l) = B(1:JRI_)
     RADE1(1:JRI_,l) = AE(1:JRI_)
     RADE2(1:JRI_,l) = BE(1:JRI_)
     PE(l)=UVE-CROSS*P(l)
     DPE(l)=DUVE-CROSS*DP(l)
     PEI(l) = RINT13(REL,AE,BE,AE,BE,NRAD,DX_,JRI_,R0_)                         
     if (qprint_) then
        !WRITE(6,8) L,P(l),DP(l),PE(l),DPE(l),PEI(l),NODEL,NODES,NODE                         
        WRITE(6,8) L,P(l),DP(l),PE(l),DPE(l),PEI(l),NODEL,NODES,NODE                         
     endif
  ENDDO
  !                         
  ! only fur lo
  !
  pr12lo(1:nloat,1:nloat,0:lomax)=0.0d0
  pi12lo(1:nloat,0:lomax)=0.0d0
  pe12lo(1:nloat,0:lomax)=0.0d0
  a1lo(1:nrad,1:nloat,0:lomax)=0.0d0
  b1lo(1:nrad,1:nloat,0:lomax)=0.0d0
  plo(:,:)  = 0.d0
  dplo(:,:) = 0.d0
  DO l=0,lomax 
     DO jlo=1,ilo(l)
        if (.not.loor(jlo,l)) CYCLE
        DELE=2.0D-3
        DELEI=0.25D0/DELE
        L_=L
        EI=elo(l,jlo)/2.d0   ! From Rydbergs to Hartree's
        !     CALCULATE FUNCTION AT EI                                          
        IF(rlo(jlo,l)) THEN
           ei=elo(l,nloat)/2.d0
           kappa=l
           ! output: A, B, UV, DUV, nodes
           CALL diracout(A,B,rel,VR,R0_,DX_,jri_,ei,l,kappa,UV,DUV,nodes,ZZ,NRAD)
           ! output: b
           CALL dergl(A,B,R0_,DX_,jri_,nrad)
           DO m = 1, jri_
              r_m = R0_*exp(DX_*(m-1))
              B(m) = B(m)*r_m/(2.d0*clight+(elo(l,jlo)-2.d0*VR(m)/r_m)/(2.d0*clight))
           ENDDO
        ELSE
           CALL OUTWIN(A,B,NODES,UV,DUV,REL,VR,R0_,DX_,jri_,ei,L_,zz,NRAD) 
        ENDIF
        ovlp = rint13(rel,a,b,a,b,NRAD,DX_,JRI_,R0_)
        ! TRX is used to normalize the wave functions.
        TRX=1.0d0/SQRT(OVLP)
        a1lo(1:JRI_,jlo,l) = A(1:JRI_)*TRX   ! this is called rf1 in other atpar
        b1lo(1:JRI_,jlo,l) = B(1:JRI_)*TRX   ! this is called rf2 in other atpar
        plo(jlo,l)  = UV*trx
        dplo(jlo,l) = DUV*trx
        pi12lo(jlo,l) = RINT13(REL,rrad1(1,l),rrad2(1,l),a1lo(1,jlo,l),b1lo(1,jlo,l),NRAD,DX_,JRI_,R0_)
        pe12lo(jlo,l) = RINT13(REL,rade1(1,l),rade2(1,l),a1lo(1,jlo,l),b1lo(1,jlo,l),NRAD,DX_,JRI_,R0_)
        if (qprint_) WRITE(6,800) L,Plo(jlo,l),DPlo(jlo,l),NODEL,NODES,NODE
     ENDDO
     DO jlo=1,ilo(l)
        CALL abc(alo,blo,clo,l,jlo,lapw,p,dp,pe,dpe,pei,plo,dplo,pi12lo,pe12lo,nrad,lomax,nloat,lmax2,rmt_,qprint_)
     ENDDO
     DO jlo=1,ilo(l)
        IF (.NOT.loor(jlo,l)) CYCLE
        DO jlop=1,ilo(l)
           IF (.NOT.loor(jlop,l)) CYCLE
           pr12lo(jlop,jlo,l) = RINT13(REL,a1lo(:,jlop,l),b1lo(:,jlop,l),a1lo(:,jlo,l),b1lo(:,jlo,l),NRAD,DX_,JRI_,R0_)
        ENDDO
     ENDDO
  ENDDO
  
  ! Finally store data into output variables
  ul_Rmt(:,:)=0.d0
  dul_Rmt(:,:)=0.d0
  do l=0,lmax2
     ul_Rmt(1,l) = P(l) 
     ul_Rmt(2,l) = PE(l)
     dul_Rmt(1,l) = DP(l)
     dul_Rmt(2,l) = DPE(l)
  enddo
  do l=0,lomax
     DO jlo=1,ilo(l)
        ul_Rmt(2+jlo,l) = PLO(jlo,l)
        dul_Rmt(2+jlo,l) = DPLO(jlo,l)
     ENDDO
  enddo
  
  ri_mat(:,:,:)=0.0
  do l=0,lmax2
     ri_mat(1,1,l)=1.0    ! <u|u>
     ri_mat(1,2,l)=0.0    ! <udot|u>
     ri_mat(2,1,l)=0.0    ! <u|udot>
     ri_mat(2,2,l)=pei(l) ! <udot|udot>
  enddo
  do l=0,lomax
     DO jlo=1,ilo(l)
        ri_mat(1,2+jlo,l) = pi12lo(jlo,l)  ! <u | u_lo>
        ri_mat(2+jlo,1,l) = pi12lo(jlo,l)  ! <u_lo | u>
        ri_mat(2,2+jlo,l) = pe12lo(jlo,l)  ! <udot | u_lo>
        ri_mat(2+jlo,2,l) = pe12lo(jlo,l)  ! <u_lo | udot>
     ENDDO
     DO jlo=1,ilo(l)
        DO jlop=1,ilo(l)
           ri_mat(2+jlo,2+jlop,l) = pr12lo(jlo,jlop,l)  ! <u_lo | u_lo >
        ENDDO
     ENDDO
  enddo
  RETURN
8   FORMAT(10X,I2,5E14.6,5X,3I2)
800 FORMAT(10X,I2,2E14.6,42X,5X,3I2)                                      
END SUBROUTINE ATPARN

subroutine get_ilo(e_store,elo_store,jatom,mult,nat,lmax2,lomax,nloat, el,elo,loor,rlo,lapw,nlov,nlo,nlon,ilo)
  ! Uses the linearization energies e_store and elo_store (from vector file) to determine the precise type 
  !  of basis function for each atom and l (LAPW, APW+lo, LAPW+LO, APW+lo+LO).
  !  In the output, it computes the following quantities:
  !
  !     el   -- actual linearization energies
  !     elo  -- actial linearization energies for LO and lo
  !     loor -- signals presence of LO
  !     lapw -- does this l have LAPW basis (otherwise is APW+lo)
  !     nlov -- number of LO's before this atom
  !     nlo  -- number of LO's on this atom
  !     nlon -- number of LO's after this atom
  !     ilo  -- number of LO's on this atom and this l.
  !     rlo  -- is this real LO or just lo
  !
  IMPLICIT NONE
  INTEGER, intent(in)  :: lmax2, jatom, mult(nat), nat, lomax, nloat
  REAL*8, intent(in)   :: e_store(0:lmax2,nat), elo_store(0:lomax,nloat,nat)
  REAL*8, intent(out)  :: el(0:lmax2), elo(0:lomax,1:nloat)
  INTEGER, intent(out) :: nlo, nlov, nlon
  INTEGER, intent(out) :: ilo(0:lmax2)
  LOGICAL, intent(out) :: loor(nloat,0:lomax), rlo(nloat,0:lomax), lapw(0:lmax2)
  ! locals
  INTEGER :: i, k, l
  ! nlo  :  #LO on this atom
  ! nlov :  #LO for atoms before this atom
  ! nlon :  #LO for atoms after this atom
  ! for this atom we also need:  ilo, rlo, loor
  nlov=0   ! all localized orbitals before jatom
  DO i=1,jatom-1
     elo(0:lomax,1:nloat)=elo_store(0:lomax,1:nloat,i)
     DO l = 0,lomax
        DO k=1,nloat
           IF (elo(l,k).LT.(995.d+0)) nlov=nlov+((2*l+1))*mult(i)
        ENDDO
     ENDDO
  ENDDO
  nlon=0  ! all localized orbitals after jatom
  DO i=jatom+1,nat  ! all localized orbitals after jatom
     elo(0:lomax,1:nloat)=elo_store(0:lomax,1:nloat,i)
     DO l = 0,lomax
        DO k=1,nloat
           IF (elo(l,k).LT.(995.0d+0)) nlon=nlon+((2*l+1))*mult(i)
        ENDDO
     ENDDO
  ENDDO
  nlo=0  ! nlo -- How many localized orbitals on this atom?
  elo(0:lomax,1:nloat) = elo_store(0:lomax,1:nloat,jatom)
  DO l=0,lomax
     DO k=1,nloat
        IF (elo(l,k).LT.(995.d+0)) THEN  ! localized orbital for this l,k
           nlo=nlo+((2*l+1))*mult(jatom) ! nlo is the number of localized orbitals for this atom
        ENDIF
     ENDDO
  ENDDO
  !
  el(0:lmax2) = e_store(0:lmax2,jatom)
  DO l=0,lmax2
     lapw(l)=.TRUE.  ! Do we have LAPW or APW+lo?
     IF(el(l).GT.150.) THEN
        el(l)=el(l)-200.d+0
        lapw(l)=.FALSE.
     ENDIF
  ENDDO

  ilo(0:lmax2)=0 ! Which l has localized orbital?
  DO l=0,lomax
     DO k=1,nloat
        loor(k,l) = .FALSE.
        rlo(k,l)  = .FALSE.
        IF (elo(l,k).LT.(995.d+0)) THEN  ! localized orbital for this l,k
           ilo(l)=ilo(l)+1
           IF(.NOT.lapw(l).AND.k.EQ.1) CYCLE    ! in this case, apw+lo, rather than lapw+LO
           IF(k.EQ.nloat) rlo(ilo(l),l)=.TRUE.  ! We have localized orbital LO
           loor(ilo(l),l) = .TRUE.              ! This is LO, rather than apw+lo
        ENDIF
     ENDDO
  ENDDO
end subroutine get_ilo

