SUBROUTINE ATPAR(NT,REL,NAT,LNSMAX,ZZ,force,lmmx,lomax,ngau,nrad,nslmax)
  use param, only : filename_V_vns
  USE readPotential, ONLY: read_V_vns2, init_V_vns, CheckFileExists
  use loabc, only : ELO, PLO, DPLO, PI12LO, PE12LO, pilolo
  use radial, only: A,B, AE, BE, VLM, A1, B1, AE1, BE1, A1LO, B1LO
  use lolog, only : loor, lapw, ilo
  use atspdt, only: LMMAX, LQIND, LM, LQNS, DP, DPE, E, P, PE,PEI, VNS1, VNS2, VNS3, GFAC!, NL
  use structure, only  : jri, dx, r0, aname
  use loint, only : VNS1LO, VNS2LO, VNS3LO
  use potnlc, only: VR 
  use structure, only : IATNR
  use mgaunt, only: gaunt1
  use mpi, only : myrank, master, mpi_bcast_V_vns, Qprint, stop_MPI
  IMPLICIT NONE
  !        Arguments
  INTEGER, intent(in) :: NT     ! number of L's which we need for hamiltonian
  INTEGER, intent(in) :: LNSMAX !
  INTEGER, intent(in) :: NAT    ! number of sorts
  REAL*8,  intent(in) :: ZZ(NAT)! nulcear charges
  LOGICAL, intent(in) :: REL    ! relativistic Dirac or no
  LOGICAL, intent(in) :: force  ! are we calculating force?
  INTEGER, intent(in) :: lmmx,lomax,ngau,nrad,nslmax ! various constants from params
  ! locals
  REAL*8 ::   TRAPLEVEL
  parameter (traplevel=1D-3)
  LOGICAL            DOTRAP
  !..................................................................
  !   ATPAR calculates the solutions u(l) of the radial Schroedinger
  !   Equation, the energy derivatives ue(l), the radial non muffin-
  !   tin integrals uvu, uevu, and the Gaunt-coefficients.
  !   Spin-orbitless relativistic equations are used.
  !..................................................................
  !   Local Scalars
  INTEGER     :: I, ICOUNT, jatom
  INTEGER     :: JRI_, L, L0, L0BEFO, LL
  INTEGER     :: LLBEFO, lm_stored, LMBEFO, ilm, LP
  INTEGER     :: LPBEFO, LQX, M, MINU, MM, MP, NODEL
  INTEGER     :: NODES, NODEU, JLO
  REAL*8      :: CROSS, DELE, DELEI
  REAL*8      :: DUV, DUVB, DUVE, DX_, E1, EI, FL, OVLP, R0_
  REAL*8      :: TRX, UV, UVB, UVE,  RNET(NRAD)
  COMPLEX*16  :: IMAG, IMAG1
  CHARACTER*67:: ERRMSG
  integer     :: klo
  !
  interface
     Function RINT13(REL,A,B,X,Y,NRAD,DX_,JRI_,R0_) result(S) ! Calculates overlap between psi_1=(A,B) and psi_2=(X,Y) functions
       REAL*8 :: S
       LOGICAL, intent(in) :: REL   ! relativistic or not
       REAL*8, intent(in)  :: A(NRAD),B(NRAD),X(NRAD),Y(NRAD)
       INTEGER, intent(in) :: NRAD, JRI_
       REAL*8, intent(in)  :: DX_, R0_
     End Function RINT13
  end interface
  !        External Functions
  EXTERNAL           CBCOMB,FORFHS,OUTERR,OUTWINC,SELECT!,ABC0
  !
  DATA IMAG /(0.0D+0,1.0D+0)/
  !
  ! If non-spherical potential exists, than we might want to calculate forces and store the
  !   matrix elements of the non-spherical potential <Y_{l'm'}u|V|Y_{lm}u>
  if (myrank.eq.master) then
     if ( (filename_V_vns/="") .and. force) then !CheckFileExists(filename_V_vns)) then
        !     force-output only 
        write(71,'( " non-spherical Matrixelements x Gaunt-Coefficient (d=du/dE):",/, " l0 ll lp m0 mm mp     uVu dVd",/, "                       uVd dVu")')
     endif
     CALL init_V_vns(filename_V_vns, 19)  ! Start reading non-spherical potential
  endif
  
  !        start loop over atoms in unitcell
  DO jatom = 1, NAT
     JRI_ = JRI(jatom)
     DX_ = dx(jatom)
     R0_ = r0(jatom)
     do i=1,jri(jatom)  ! Radial mesh created just once at the beginning
        RNET(i)=r0(jatom)*exp(dx(jatom)*(i-1.d0))
     enddo

     if (Qprint) WRITE (6,6020) jatom, aname(jatom) ! Starts calculating potential parameters
     DELE = 2.0D-3
     DELEI = 0.250D0/DELE
     DO l=0,NT-1
        !  Experimental Trap location
        DOTRAP=.true.
1500    CONTINUE
        FL = L
        EI = E(l+1,jatom)/2.D0
        !        calculate energy-derivative by finite difference
        !        DELE is the up and downward energy-shift in Hartrees
        E1 = EI - DELE
        CALL OUTWINC(A,B,NODEL,UVB,DUVB,REL,VR(:,jatom),Rnet,DX_,JRI_,EI-DELE,FL,ZZ(jatom),nrad)
        OVLP = RINT13(REL,A,B,A,B,NRAD,DX_,JRI_,R0_) ! OVLP = A**2+B**2
        TRX = 1.0D+0/SQRT(OVLP)
        AE(1:JRI_) = A(1:JRI_)*TRX
        BE(1:JRI_) = B(1:JRI_)*TRX
        UVB  = UVB*TRX
        DUVB = DUVB*TRX
        ! results are: AE(:), BE(:), UVB, DUVB, E1
        ! outputs: UVE, DUVE, NODE, OVLP
        !     CALCULATE FUNCTION AT EI+DELE
        E1 = EI + DELE
        CALL OUTWINC(A,B,NODEU,UVE,DUVE,REL,VR(:,jatom),Rnet,DX_,JRI_,EI+DELE,FL,ZZ(jatom),nrad)
        OVLP = RINT13(REL,A,B,A,B,NRAD,DX_,JRI_,R0_)
        TRX = 1.0D+0/SQRT(OVLP)
        UVE = DELEI*(TRX*UVE-UVB)
        DUVE = DELEI*(TRX*DUVE-DUVB)
        AE(1:JRI_) = DELEI*( A(1:JRI_)*TRX - AE(1:JRI_))
        BE(1:JRI_) = DELEI*( B(1:JRI_)*TRX - BE(1:JRI_))
        ! results are: A(:), B(:), P(:), DP(:), UV, DUV
        !     calculate function at EI
        CALL OUTWINC(A,B,NODES,UV,DUV,REL,VR(:,jatom),Rnet,DX_,JRI_,EI,FL,ZZ(jatom),nrad)
        OVLP = RINT13(REL,A,B,A,B,NRAD,DX_,JRI_,R0_)
        TRX = 1.0D+0/SQRT(OVLP)
        P(l+1,jatom) = TRX*UV
        DP(l+1,jatom) = TRX*DUV
        A(1:JRI_) = A(1:JRI_)*TRX
        B(1:JRI_) = B(1:JRI_)*TRX
        !     insure orthogonalization
        CROSS = RINT13(REL,A,B,AE,BE,NRAD,DX_,JRI_,R0_)
        IF (Qprint .and. CROSS .GT. 0.05D+0) WRITE (6,6070) L, CROSS, OVLP
        AE(1:JRI_) = AE(1:JRI_) - CROSS*A(1:JRI_)
        BE(1:JRI_) = BE(1:JRI_) - CROSS*B(1:JRI_)
        ! stores radial functions into A1, B1, AE1, BE1
        IF ((l+1) .LE. NSLMAX) THEN
           A1(1:JRI_,l+1) = A(1:JRI_)
           B1(1:JRI_,l+1) = B(1:JRI_)
           AE1(1:JRI_,l+1) = AE(1:JRI_)
           BE1(1:JRI_,l+1) = BE(1:JRI_)
        ENDIF
        PE(l+1,jatom) = UVE - CROSS*P(l+1,jatom)
        DPE(l+1,jatom) = DUVE - CROSS*DP(l+1,jatom)
        PEI(l+1,jatom) = RINT13(REL,AE,BE,AE,BE,NRAD,DX_,JRI_,R0_)

        ! Trap P(J,jatom) too small
        IF((Abs(P(l+1,jatom)) .lt. traplevel).and.DOTRAP)then
           DOTRAP=.false.
           ! An alternative might be to set the level to LAPW....?
           E(l+1,jatom)=E(l+1,jatom)-0.05
           write(21,1501) E(l+1,jatom)
           write(6,1501) E(l+1,jatom)
1501       format(':WARN  : P(J,JATOM) almost zero. Shifted Energy by -0.05 down to ',F10.4)
           goto 1500
        endif
        if (Qprint) WRITE (6,6030) l, P(l+1,jatom), DP(l+1,jatom), PE(l+1,jatom),DPE(l+1,jatom), PEI(l+1,jatom), NODEL, NODES,NODEU
     ENDDO
     !     
     !     and now for local orbitals
     !
     if (Qprint) write(6,6021)
     DO L = 0, LOMAX
        IF (LOOR(L,jatom)) THEN
           FL = L
           do jlo=1,ilo(l,jatom)
              if ((.not.lapw(l,jatom)).and.jlo.eq.1) cycle  ! apw+lo
              EI = Elo(L,jlo,jatom)/2.D0  ! from Ry2Hartree
              !     calculate function at EI
              CALL OUTWINC(A,B,NODES,UV,DUV,REL,VR(:,jatom),Rnet,DX_,jri_,Ei,FL,ZZ(jatom),nrad) 
              OVLP = rint13(REL,A,B,A,B,NRAD,DX_,JRI_,R0_)
              TRX = 1.0D+0/SQRT(OVLP)
              A(1:JRI_) = A(1:JRI_)*TRX
              B(1:JRI_) = B(1:JRI_)*TRX
              PLO(l,jlo,jatom) = UV*TRX
              DPLO(l,jlo,jatom) = DUV*TRX
              pi12lo(l,jlo,jatom) = RINT13(REL,A1(:,L+1),B1(:,L+1),A,B,NRAD,DX_,JRI_,R0_)
              pe12lo(l,jlo,jatom) = RINT13(REL,AE1(:,L+1),BE1(:,L+1),A,B,NRAD,DX_,JRI_,R0_)
              IF ((l+1) .LE. NSLMAX) THEN
                 A1LO(1:JRI_,jlo,l) = A(1:JRI_)
                 B1LO(1:JRI_,jlo,l) = B(1:JRI_)
              ENDIF
              if (Qprint) WRITE(6,6031) l,PLO(l,jlo,jatom),DPLO(l,jlo,jatom),PI12LO(l,jlo,jatom),PE12LO(l,jlo,jatom), NODEL, NODES, NODEU
           enddo
           do jlo=1,ilo(l,jatom)
              if ((.not.lapw(l,jatom)).and.jlo.eq.1) cycle ! apw+lo
              do klo=1,ilo(l,jatom)
                 if ((.not.lapw(l,jatom)).and.klo.eq.1) cycle ! apw+lo
                 pilolo(l,jlo,klo,jatom) = RINT13(REL,a1lo(:,jlo,l),b1lo(:,jlo,l),a1lo(:,klo,l),b1lo(:,klo,l),NRAD,DX_,JRI_,R0_)
                 if (Qprint) write(6,'(2i5,3f15.6)') jlo,klo,PILOLO(L,jlo,klo,jatom),ELO(L,jlo,jatom),ELO(L,klo,jatom)
              enddo
           enddo
        endif
     enddo

!!!!  Below we calculate the matrix elements of potential in the following way
!!!!   <Y_{l0,m}| V(ll,mm) | Y_{lp,mp}>
!!!! where V(l,m) is non-spherical potential in the muffin-thin sphere     
     JRI_ = JRI(JATOM)
     DX_ = dx(JATOM)
     R0_ = r0(jatom)
!!!  We need to do this in the same loop because the wave functions, such as a1lo(:,:,:), b1lo(:,:,:) and c1lo(:,:,:), are not saved for all atoms, but just for the current atom.

     if (myrank.eq.master) then
        if (force) write(71,'(i3," .Atom ")') jatom
        !        read total nonspherical potential of TAPE19=VNS
        !        norm of VLM=VLM*1
        CALL read_V_vns2(Vlm,lmmax(jatom),lm_stored,LM(:,:,jatom),jatom,nrad,jri(jatom),lmmx,lnsmax)
     endif
     !call mpi_bcast_V_vns(Vlm,lmmax(jatom),lm_stored,LM(:,:,jatom),jatom,nrad,jri(jatom),lmmx,lnsmax)
     call mpi_bcast_V_vns(Vlm,lmmax(jatom),lm_stored,LM(:,:,jatom),nrad,lmmx)
     !        calculate possible nonspherical contributions to H
     if (Qprint) then
        WRITE (6,*) ' POSSIBLE NONSPHERICAL CONTRIBUTIONS TO H: '
        WRITE (6,*) '(L0,LL,LP,M, MM,MP, GNT,  U V U,    UE V UE,','   U V UE     )'
     endif
     LQIND(jatom) = 0
     !     
     DO l0=0,LNSMAX
        DO lp=0,LNSMAX
           DO ilm=1,lm_stored
              ll = ABS(lm(1,ilm,jatom))
              mm = lm(2,ilm,jatom)
              !     which gaunts are not zero ?
              ! Wigner-Eckart for <(l0,m)| O(ll,mm) |(lp,mp)>
              !    ll-l0 < lp < ll+l0
              IF (MOD((l0+lp+ll),2) .EQ. 1) CYCLE
              IF ((l0+lp-ll) .LT. 0) CYCLE
              IF ((l0-lp+ll) .LT. 0) CYCLE
              IF ((lp-l0+ll) .LT. 0) CYCLE
              !     
240           CONTINUE
              DO m=-l0,l0
                 DO mp=-lp,lp
                    IF (-m+mm+mp .NE. 0) CYCLE       ! for matrix element <(l0,m)| O(ll,mm) |(lp,mp)> should work.
                    LQIND(jatom) = LQIND(jatom) + 1  ! we have one more nonzero matrix element
                    IF (LQIND(jatom) .GT. NGAU) THEN ! there is an upper limit to the number of these matrix elements
                       !     Error: more then NGAU gaunt-factors - stop execution
                       GOTO 910
                    ENDIF
                    ! We need Gaunt coefficients for real harmonics, which is a bit more cumbersome...
                    GFAC(LQIND(jatom),jatom) = GAUNT1(l0,ll,lp,m,mm,mp)   ! store the gaunt coefficient
                    IF (IATNR(jatom) .GT. 0) THEN                         ! the number in structure file in front of atoms position
                       IF (ll.EQ.3.OR.ll.EQ.7.OR.ll.EQ.9) THEN            ! positive means cubic symmetry.?
                          IF(mm.LT.0) THEN
                             GFAC(LQIND(jatom),jatom) = GFAC(LQIND(jatom),jatom)*IMAG
                          ELSE
                             GFAC(LQIND(jatom),jatom) = -GFAC(LQIND(jatom),jatom)*IMAG
                          ENDIF
                       ENDIF
                    ELSE
                       IF (mm .NE. 0) THEN
                          MINU = 1
                          IMAG1 = (1.0D+0,0.0D+0)
                          IF (lm(1,ilm,jatom).LT.0) THEN  ! real harmonics with sin
                             IMAG1 = -IMAG
                             MINU = -1
                          ENDIF
                          IF (MOD(mm,2) .EQ. 1) THEN   ! mm is odd
                             IMAG1 = -IMAG1
                             MINU = -MINU
                          ENDIF
                          IF (mm .GT. 0) MINU = 1
                          GFAC(LQIND(jatom),jatom) = GFAC(LQIND(jatom),jatom)*IMAG1*dble(MINU)/SQRT(2.0D+0)
                       ENDIF
                    ENDIF
                    LQNS(1,LQIND(jatom),jatom) = l0 + 1
                    LQNS(2,LQIND(jatom),jatom) = ll
                    LQNS(3,LQIND(jatom),jatom) = lp + 1
                    LQNS(4,LQIND(jatom),jatom) = ilm
                    LQNS(5,LQIND(jatom),jatom) = m
                    LQNS(6,LQIND(jatom),jatom) = mp
                 ENDDO
              ENDDO
              !
              IF (mm.GT.0) THEN
                 mm = -mm
                 GOTO 240
              ENDIF
           ENDDO
        ENDDO
     ENDDO
     !
     IF ((IATNR(jatom).GT.0) .AND. (LMMAX(jatom).GT.1)) THEN
        !        combine the radial total potential coefficients according to cubic harmonic
        CALL CBCOMB(JRI_,lm_stored,VLM,LM,jatom,nat,nrad,lmmx)
     ENDIF
     !   integrate product of radialfunctions (and derivatives) and
     !   total potential coefficients (combinations) from 0 to R-MT
     ICOUNT = 0
     L0BEFO = 999
     LLBEFO = 999
     LPBEFO = 999
     LMBEFO = 999
     DO LQX = 1, LQIND(jatom)
        L0  = LQNS(1,LQX,jatom) - 1
        LL  = LQNS(2,LQX,jatom)
        LP  = LQNS(3,LQX,jatom) - 1
        ilm = LQNS(4,LQX,jatom)
        IF ((L0 .EQ. L0BEFO) .AND. (LL .EQ. LLBEFO) .AND. (LP .EQ. LPBEFO) .AND. (ilm .EQ. LMBEFO)) THEN
           M = LQNS(5,LQX,jatom)
           MP = LQNS(6,LQX,jatom)
        ELSE
           A(1:JRI_)  = A1(1:JRI_,L0+1) *VLM(1:JRI_,ilm)
           B(1:JRI_)  = B1(1:JRI_,L0+1) *VLM(1:JRI_,ilm)
           AE(1:JRI_) = AE1(1:JRI_,L0+1)*VLM(1:JRI_,ilm)
           BE(1:JRI_) = BE1(1:JRI_,L0+1)*VLM(1:JRI_,ilm)
           !
           ICOUNT = ICOUNT + 1
           !
           VNS1(L0+1,ilm,LP+1,jatom) = RINT13(REL,A,B,A1(:,lp+1),B1(:,lp+1),NRAD,DX_,JRI_,R0_)      ! <A1|V|A1>
           VNS2(L0+1,ilm,LP+1,jatom) = RINT13(REL,AE,BE,AE1(:,lp+1),BE1(:,lp+1),NRAD,DX_,JRI_,R0_)  ! <AE1|V|AE1>
           VNS3(L0+1,ilm,LP+1,jatom) = RINT13(REL,A,B,AE1(:,lp+1),BE1(:,lp+1),NRAD,DX_,JRI_,R0_)    ! <A1|V|AE1>
           !     integrals needed for local orbitals
           IF (L0 .LE. LOMAX) THEN
              do jlo=1,ilo(l0,jatom)
                 VNS1LO(L0+1,ilm,LP+1,jlo,jatom) = 0.0D0
                 VNS2LO(L0+1,ilm,LP+1,jlo,jatom) = 0.0D0
                 VNS3LO(L0+1,ilm,LP+1,jlo,:,jatom) = 0.0d0
                 if (.not.lapw(l0,jatom).and.jlo.eq.1) cycle
                 if(loor(l0,jatom)) THEN
                    A(1:JRI_) = A1LO(1:JRI_,jlo,L0)*VLM(1:JRI_,ilm)
                    B(1:JRI_) = B1LO(1:JRI_,jlo,L0)*VLM(1:JRI_,ilm)
                 else
                    exit
                 endif
                 !
                 VNS1LO(L0+1,ilm,LP+1,jlo,jatom) = RINT13(REL,A,B,A1(:,lp+1),B1(:,lp+1),NRAD,DX_,JRI_,R0_)
                 VNS2LO(L0+1,ilm,LP+1,jlo,jatom) = RINT13(REL,A,B,AE1(:,lp+1),BE1(:,lp+1),NRAD,DX_,JRI_,R0_)
                 !
                 IF (LP .LE. LOMAX) THEN
                    do klo=1,ilo(lp,jatom)
                       if (.not.lapw(lp,jatom).and.klo.eq.1) cycle
                       VNS3LO(L0+1,ilm,LP+1,jlo,klo,jatom) = RINT13(REL,A,B,A1lo(:,klo,lp),B1lo(:,klo,lp),NRAD,DX_,JRI_,R0_)
                    enddo
                 ENDIF
                 M = LQNS(5,LQX,jatom)
                 MP = LQNS(6,LQX,jatom)
              enddo
           ENDIF !        end of lo-integrals
           L0BEFO = L0
           LLBEFO = LL
           LPBEFO = LP
           LMBEFO = ilm
           M = LQNS(5,LQX,jatom)
           MP = LQNS(6,LQX,jatom)
        ENDIF
     ENDDO
     if (force .and. myrank.eq.master) then
        call forfhs(jatom) 
        !               final line for each atom = 6 x 0 and 8 x 0.0d0
        write(71,'(6i3,4es19.12,/,6(3x),4es19.12)') (0, i=1,6),(0.d0,i=1,8)
     endif
     if (Qprint) WRITE (6,*) '   NUMBER OF RADIAL INTEGRALS FOR ATOM', JATOM,' = ',ICOUNT
  ENDDO
  !
  RETURN
  !
910 CONTINUE
  CALL OUTERR('ATPAR','more than NGAU gaunts')
  WRITE (ERRMSG,9010) '  NGAU,   L0,   LP,   LL,    M,   MP,   MM'
  CALL OUTERR('ATPAR',ERRMSG)
  WRITE (ERRMSG,9020) NGAU, L0, LP, LL, M, MP, MM
  CALL OUTERR('ATPAR',ERRMSG)
  call stop_MPI
  STOP 'ATPAR - Error'
  !
!6000 FORMAT (/,10X,'ATOMIC SPHERE DEPENDENT PARAMETERS FOR ATOM  ',A10)
6020 FORMAT (/,10X,'POTENTIAL PARAMETERS FOR JATOM=',I3,' ',A10,/,11X,'L',5X,'U(R)',10X,'U''(R)',9X,'DU/DE',8X,'DU''/DE',6X,'NORM-U''')
6021 FORMAT (/,10X,'LOCAL ORBITAL POTENTIAL PARAMETERS',/,11X,'L',5X,'U(R)',10X,'U''(R)',5X,'NORM U1U2',8X,'NORM UE1U2')
6030 FORMAT (10X,I2,5E14.6,5X,3I2)
6031 FORMAT (10X,I2,4E14.6,5X,3I2)
6070 FORMAT (10X,'FOR L=',I2,' CORRECTION=',E14.5,' OVERLAP=',E14.5)
9010 FORMAT (A)
9020 FORMAT (I6,',',6(I5,','))
END SUBROUTINE ATPAR
