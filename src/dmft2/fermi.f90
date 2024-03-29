
SUBROUTINE fermi_tetra(EF,weigh,nbmax,ELECN,Eb,nkpt,nume,nspin,qprint_)
  !     EF      -- The Fermi energy to be calculated
  !     weigh   -- tetrahedron weights
  !     nbmax   -- last band which needs to be considered for density
  !
  !     ELECN   -- The number of valence electrons
  !     Eb      -- all bands
  !     qprint_ -- controls printing
  IMPLICIT NONE
  REAL*8, intent(out)    :: EF
  REAL*8, intent(out)    :: weigh(nkpt*nspin,nume)
  INTEGER, intent(out)   :: nbmax
  REAL*8, intent(in)     :: ELECN
  REAL*8, intent(in)     :: Eb(nume,nkpt,nspin)
  INTEGER, intent(in)    :: nkpt, nume, nspin
  LOGICAL, intent(in)    :: qprint_
  ! locals
  REAL*8,  PARAMETER  :: TEST1=2.D-8                                                      
  REAL*8, allocatable :: Ebmax(:), Ebmin(:)
  REAL*8, allocatable :: en2(:,:), weight2(:,:)
  INTEGER          :: ik, inum!, iw(nw)
  REAL*8           :: elecn_, cordeg, emax
  INTEGER          :: icor, i, iloop, ispin, itap, n, nemax, num, inum_spin, ik_spin
  !     icor switches non-linear correction  on/off
  !     cordeg is a switch to correct occupancy of degenerate states
  icor=1
  if(abs(ef).ge.100.d0) icor=0
  cordeg=-1.d-6
  if(ef.gt.0.d0)  cordeg=-ef   
  if(ef-100.d0.gt.0.d0) cordeg=-ef+100.d0
  if(ef.lt.0.d0)        cordeg=ef
  if(ef.lt.-100.d0)     cordeg=ef+100.d0
  if(abs(cordeg).gt.0.01d0)  cordeg=-1.d-6
  !
  if (qprint_) write(6,'(" BZ-integration with TETRA-program.   icor=:",I2  )') icor
  if(cordeg.lt.0.d0 .and. qprint_) write(6,'(" Equal occupancy of degenerate states, tol=:",E8.1)') -cordeg
  
  nemax=nume
  ELECN_ = ELECN-1.D-10

  ! finds maximum and minimum energy for each band
  allocate( Ebmax(nume), Ebmin(nume) )
  Ebmin(:)=999.d0; Ebmax(:)=-999.d0
  do ispin=1,nspin
     DO ik=1,nkpt
        DO inum=1,nume
           if (Eb(inum,ik,ispin).LT.3.d0) Ebmax(inum) = max(Ebmax(inum),Eb(inum,ik,ispin))
           Ebmin(inum) = min(Ebmin(inum),Eb(inum,ik,ispin))
        ENDDO
     ENDDO
  enddo

  ! below we need different memory arrangement of energies and weights
  ! we create en2(nkpt,nemax*nspin) and weight2(nkpt,nemax*nspin)
  allocate(en2(nkpt,nemax*nspin))
  if (qprint_) write(6,*)'call eord...' 
  do ik=1,nkpt
     en2(ik,:nemax) = eb(:nemax,ik,1)
     if (nspin.eq.2) en2(ik,nemax:2*nemax) = eb(:nemax,ik,2)
  enddo
  
  allocate(weight2(nkpt,nemax*nspin))
  weight2=0.0d0
  
  if (qprint_) write(6,*)'call dos...'
  call dos(nemax*nspin,nkpt,en2,weight2,elecn_/2.d0*nspin,ef,icor,qprint_) 
  if (qprint_) write(6,*)'call correctw...'
  call correctw(weight2,nemax,nspin,Eb,nkpt,nume,cordeg,qprint_)
  
  nbmax=0  
  emax=-10.d0
  DO ik=1,nkpt
     DO inum=1,nemax
        do ispin=1,nspin
           inum_spin = inum + (ispin-1)*nemax
           ik_spin = ik + (ispin-1)*nkpt
           if(abs(weight2(ik,inum_spin)).gt.test1) then
              nbmax=max(nbmax,inum)
              !emax=max(emax,Eb(inum,ik,ispin))
              emax = max(emax,en2(ik,inum_spin))
           end if
           WEIGH(ik_spin,inum) = WEIGHT2(ik,inum_spin)*2.0d0/nspin
        enddo
     ENDDO
  ENDDO
  if (qprint_) write(6,*) '  number of occupied bands:',nbmax
  if (qprint_) write(6,*) '  highest energy:',emax
  if(ef.gt.emax.and.abs(emax-ebmin(nbmax+1)).gt.1.d-3) then
     if (qprint_) write(6,*) 'insulator !'
     if((EF-emax).lt.1.d-4.or.(EF-emax).ge.1.d-4) then
        if (qprint_) write(6,*) 'EF-inconsistency corrected'
        if (qprint_) write(21,'(A)') '       Insulator, EF-inconsistency corrected'
        EF=emax
     endif
  endif
  deallocate( en2 )
  deallocate( weight2 )
  deallocate( Ebmax, Ebmin )
END SUBROUTINE fermi_tetra

subroutine correctw(weight2,nemax,nspin,Eb,nkpt,nume,cordeg,qprint_)
  IMPLICIT NONE
  REAL*8,  intent(inout):: weight2(nkpt,nemax*nspin)
  INTEGER, intent(in)   :: nkpt, nume
  REAL*8,  intent(in)   :: cordeg
  INTEGER, intent(in)   :: nemax, nspin
  REAL*8,  intent(in)   :: Eb(nume,nkpt,2)
  LOGICAL, intent(in)   :: qprint_
  ! local
  REAL*8  :: emax, wecp
  INTEGER :: icp, ifcp, ifcpt, ispin, ik, ndeg, nn
  DO ik=1,nkpt
     !
     ! ### 2007-08-31, Clas Persson, Juergen Spitaler, and Claudia Ambrosch-Draxl
     ! We suggest to change on page 90 in the usersguide.pdf [note, now abs(eval)]:
     !
     !                                     ...    abs(eval).gt.100 specifies
     ! the use of the standard tetrahedron method instead of the modified one
     ! (see above). Using eval.lt.0 in combination with TETRA forces the
     ! occupancy of degenerate states to be equal in order to avoid incorrect
     ! split of these states, which otherwise may occur for partially filled
     ! states. Degeneracy is here determined by a tolerance of abs(eval) for
     ! -1.gt.eval.lt.0, and of 1e-6 Ry for eval.le.-1.
     !
     if(cordeg.lt.0.d0) then
        if(ik.eq.1) ifcpt=1
        do ispin=1,nspin
           nn=1
           do while(nn.le.nemax)
              ifcp=0
              ndeg=1
              wecp=weight2(ik,nn+(ispin-1)*nemax)
              !
              ! check degeneracy
              do while((nn+ndeg).le.nemax)
                 if(abs(eb(nn+ndeg,ik,ispin)-eb(nn,ik,ispin)).gt.abs(cordeg)) exit ! not degenerate
                 if(abs(weight2(ik,nn+ndeg+(ispin-1)*nemax)-weight2(ik,nn+(ispin-1)*nemax)).ge.1e-6) ifcp=1
                 wecp=wecp+weight2(ik,(nn+ndeg)+(ispin-1)*nemax)
                 ndeg=ndeg+1
              enddo
              !
              ! equalizes occupancy and outputs to case.output2
              if(ifcp.eq.1) then
                 if(ifcpt.eq.1) then
                    if (qprint_) write(6,'("k-pnt spin band  energy          old/new occ.    ")')
                    ifcpt=0
                 endif
                 do icp=0,ndeg-1
                    if (qprint_) write(6,'(3i5,f10.6,2x,f9.6,"/",$)') ik,ispin,nn,eb(nn+icp,ik,ispin),weight2(ik,(nn+icp)+(ispin-1)*nemax)
                    weight2(ik,(nn+icp)+(ispin-1)*nemax) = wecp/dble(ndeg)
                    if (qprint_) write(6,'(f9.6)') weight2(ik,(nn+icp)+(ispin-1)*nemax)
                 enddo
              endif
              nn=nn+ndeg
           enddo
        enddo
     endif
  ENDDO
end subroutine correctw

!-----------------------------------------------------------------
!-----------------------------------------------------------------
!----                                                         ----
!----  BLOCK KINTR                                            ----
!----  CALCULATION OF SAMPLING WEIGHTS                        ----
!----                                                         ----
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!.....................................................DOS.........
SUBROUTINE DOS(NB,NKP,EB,WGHT,RNTOT,EF,icor,qprint_)
  ! **                                                              *
  ! **  CALCULATES THE SAMPLING WEIGHTS FROM TETRAHEDRON INTEGRATION*
  ! **                                                              *
  ! **  INPUT :                                                     *
  ! **    NB          NUMBER OF BANDS                               *
  ! **    NKP         NUMBER OF IRREDUCIBLE K-POINTS                *
  ! **    EB          ENERGIES ( DETERMINE FERMI SURFACE )          *
  ! **    RNTOT       NUMBER OF OCCUPIED STATES                     *
  ! **    W           (INTEGER) WORK ARRAY                          *
  ! **    NWX         LENGTH OF WROK ARRAY W                        *
  ! **  OUTPUT :                                                    *
  ! **    WGHT        SAMPLING WEIGHTS                              *
  ! **    EF          FERMI LEVEL                                   *
  ! **                                                              *
  ! **  AUTHOR : PETER E. BLOECHL                                   *
  ! **                                                              *
  ! **  
  IMPLICIT NONE
  REAL*8, intent(in)     :: EB(NKP,NB), RNTOT
  REAL*8, intent(out)    :: WGHT(NKP,NB), EF
  INTEGER, intent(in)    :: NKP, NB, icor
  LOGICAL, intent(in)    :: qprint_
  ! locals
  REAL*8 :: wSUM
  CHARACTER *67 ::   ERRMSG
  REAL*8,PARAMETER :: TOLMAX = 1.D-5
  ! a larger TOLMAX may not lead to the correct number of electrons!
  ! eventually one could issue just a warning, but continue if normalization
  ! is only slightly violated.
  !
  !-----------------------------------------------------------------
  !--  CALCULATE FERMI LEVEL                                       -
  !-----------------------------------------------------------------
  CALL EFERMI(RNTOT,EF,TOLMAX,NB,NKP,EB,qprint_)
  if (qprint_) write(6,*) '  FERMI ENERGY AT ',EF                                   
  !-----------------------------------------------------------------
  !--  CALCULATE WEIGHTS                                           -
  !-----------------------------------------------------------------
  CALL SAMFAC(NB,NKP,EB,EF,WGHT,icor)
  !-----------------------------------------------------------------
  !--  CHECK WHETHER SUMRULE IS FULLFILLED                         -
  !-----------------------------------------------------------------
  wSUM=sum(WGHT)                                                       
  IF(DABS(wSUM-RNTOT).GT.TOLMAX) THEN
     if (qprint_) WRITE (6,9000) wSUM,RNTOT
     if (qprint_) write(21,9001) wSUM,RNTOT
  ELSE
     if (qprint_) write(6,*) '  SUM RULE OF TETRAHEDRON INTEGRATION CHECKED '
  END IF
  RETURN
9000 FORMAT(' RESULT OF INTEGRATION: ',f10.5, '; SHOULD BE: ',f10.5)
9001 FORMAT(':WARN : RESULT OF INTEGRATION: ',f10.5, '; SHOULD BE: ',f10.5)
END SUBROUTINE DOS


module ReadTetra
  INTEGER :: ipos, ntet, Tinit, mwrit, nrec
  INTEGER, allocatable :: iwork(:)
contains
  SUBROUTINE T_init_(VOL,ntet_,nkp)
    IMPLICIT NONE
    REAL*8, intent(out) :: VOL
    INTEGER, intent(out):: ntet_
    INTEGER, intent(in) :: nkp
    ! locals
    INTEGER :: nkp_
    REWIND 14
    READ(14,'(2i10,e20.12,2i10)') nkp_,ntet,VOL,mwrit,nrec
    if(nkp_.ne.nkp) then
       WRITE(6,*) 'FERMI','number of k-points inconsistent when reading kgen'
       WRITE(6,*) 'FERMI','check case.in1 case.energy and case.kgen files!'
       STOP  'FERMI - Error'
    endif
    ntet_ = ntet
    ipos=0
    Tinit=0
    if (allocated(iwork)) deallocate(iwork)
    allocate(iwork(5*mwrit))
    !print *, 'Just allocated mwrit=', mwrit, 'ntet=', ntet, 'nkp=', nkp, 'Vol=', Vol
  END SUBROUTINE T_init_
  
  SUBROUTINE Tetra1(ITET,IWGHT,IKP)
    IMPLICIT NONE
    INTEGER, intent(out)   :: IKP(4), IWGHT
    INTEGER, intent(in)    :: ITET
    ! locals
    !INTEGER :: ipos, ntet
    REAL*8  :: V
    INTEGER :: i, IP, irec, nkp, nleft
    ! each tetragedron is written with 5 ints : (weight,k1,k2,k3,k4)
    !IF(Tinit.EQ.1) THEN  ! Start reading the file from the beginning
    !   REWIND 14
    !   READ(14,'(2i10,e20.12,2i10)') NKP,NTET,V,mwrit,NREC                             
    !   Tinit=0
    !   ipos=0
    !   if (allocated(iwork)) deallocate(iwork)
    !   allocate(iwork(5*mwrit))
    !END IF
    IF(ITET.GT.NTET) THEN
       WRITE(6,*) 'FERMI','ASK FOR NONEXISTING TETRAHEDRON'
       WRITE(6,*) 'FERMI','STOP IN TETR1'
       STOP  'FERMI - Error'
    ENDIF
    irec=(ITET-1)/MWRIT+1
    IF(irec.NE.ipos) THEN
       READ(14,'(6i10)') iwork
       !---------------------------------------------------------------
       NLEFT=NTET-(IREC-1)*MWRIT
       NLEFT=5*min(MWRIT,NLEFT)
       !---------------------------------------------------------------
       ipos=irec
    END IF
    IP=5*(ITET-1-(IPOS-1)*MWRIT)                                     
    IWGHT=IWORK(IP+1)
    DO I=1,4
       IKP(I)=IWORK(IP+1+I)                                             
    ENDDO
    !print *, 'current irec=', irec, 'ipos=', ipos, 'nrec=', nrec, 'ikp=', ikp, 'w=', iwght
  END SUBROUTINE Tetra1
  
  SUBROUTINE T_finish_()
    IMPLICIT NONE
    if (allocated(iwork)) deallocate( iwork )
  END SUBROUTINE T_finish_
end module ReadTetra

!.....................................................EFERMI......
SUBROUTINE EFERMI(RNTOT,EF,TOLMAX,NB,NKP,EB,qprint_)               
  ! **                                                              *
  ! **  CALCUALTES THE FERMILEVEL BY INTEGRATION OF THE TOTAL       *
  ! **  DENSITY OF STATES                                           *
  ! **                                                              *
  ! **  INPUT :                                                     *
  ! **    RNTOT       NUMBER OF OCCUPIED STATES                     *
  ! **    NB          NUMBER OF BANDS                               *
  ! **    NKP         NUMBER OF IRREDUCIBLE K-POINTS                *
  ! **    EB          ENERGIES ( DETERMINE FERMI SURFACE )          *
  ! **    W           (INTEGER) WORK ARRAY                          *
  ! **    NWX         LENGTH OF WROK ARRAY W                        *
  ! **    TOLMAX      TOLERANCE IN THE NUMBER OF STATES AT EF       *
  ! **  OUTPUT :                                                    *
  ! **    EF          FERMI LEVEL                                   *
  ! **                                                              *
  use ReadTetra, ONLY: T_init_, Tetra1, T_finish_
  IMPLICIT NONE
  REAL*8, intent(in)     :: EB(NKP,NB)
  INTEGER, intent(in)    :: NB, NKP
  REAL*8, intent(out)    :: EF
  REAL*8, intent(in)     :: RNTOT, TOLMAX
  LOGICAL, intent(in)    :: qprint_
  ! locals
  INTEGER, PARAMETER     :: NP=1000  ! Number of energy divisions in the interval
  REAL*8                 :: dW_ONOS(NP)
  REAL*8                 :: dW_OSOS(NP)
  CHARACTER*67           :: ERRMSG
  REAL*8  :: E(4), DE, EMAX, EMIN, ESTEP, wSUM, tol, vol, vol1
  INTEGER :: IKP(4), I, IB, IK, ILOOP, INIT, IP, ITET, IWGHT, MWRIT, NTET, ONOS, OSOS, OWORK
  !-----------------------------------------------------------------
  !--  FIND EMIN EMAX     (ENERGYBANDS ARE ASSUMED                 -
  !--                      TO BE ORDERED WITH RESPECT TO SYMMETRY  -
  !-----------------------------------------------------------------
  Emin = minval(Eb(:,:))
  Emax = maxval(Eb(:,:))
  !print*, 'Emin=', Emin, 'Emax=', Emax
  DE=(EMAX-EMIN)/DBLE(NP-1)
  Emax = Emax+DE
  Emin = Emin-DE
  DO ILOOP=1,6
     dW_ONOS(:)=0.d0
     dW_OSOS(:)=0.d0
     CALL T_init_(VOL1,ntet,nkp)
     wSUM=0.D0
     DO ITET=1,NTET
        CALL Tetra1(itet, iwght, ikp)
        VOL=VOL1*DBLE(IWGHT)
        wSUM=wSUM+VOL
        DO IB=1,NB           
           DO I=1,4           
              E(I)=EB(IKP(I),IB)                                               
           ENDDO
           CALL TOTNOS(VOL,E,EMIN,EMAX,NP,dW_ONOS,dW_OSOS)                  
        ENDDO
     ENDDO
     
     IF(DABS(wSUM-1.D0).GT.1.D-5) then
        WRITE(6,*) 'FERMI',' TETRAHEDRA DO NOT FILL VOLUME wsum=', wSUM
        WRITE(6,*) 'FERMI',' STOP IN EFERMI'
        STOP  'FERMI - Error'
     ENDIF
     !-----------------------------------------------------------------
     !--  GET FERMI LEVEL                                             -
     !-----------------------------------------------------------------
     tol=tolmax
     CALL EFI(RNTOT,TOL,EF,EMIN,EMAX,NP,dW_ONOS,dW_OSOS)
     !-----------------------------------------------------------------
     !--  CHECK ACCURACY AND RESTART IF NECCESARY                     -
     !-----------------------------------------------------------------
     IF(TOL.GT.TOLMAX) THEN                                           
        ESTEP=(EMAX-EMIN)/DBLE(NP-1)
        IP=1+(EF-EMIN)/ESTEP                                           
        EMIN=EMIN+ESTEP*DBLE(IP-1)                                     
        EMAX=EMIN+ESTEP                                                
        if(estep.lt.1.d-8) then
           if (qprint_) write(*,*) 'WARNING: EF not accurate, new emin,emax,NE-min,', 'NE-max',emin,emax,dW_ONOS(IP),dW_ONOS(IP+1)
           ef=(emin+emax)/2.d0
           EXIT
        endif
        IF(RNTOT-dW_ONOS(IP).LE.TOLMAX) THEN
           EF=EMIN                                                      
           EXIT
        ELSE IF(dW_ONOS(IP+1)-RNTOT.LE.TOLMAX) THEN              
           EF=EMAX
           EXIT
        END IF
        IF(ILOOP.GT.5) THEN
           WRITE(6,*) 'FERMI',' CANNOT FIND FERMI LEVEL tol=', tol, 'emin,emax=', emin,emax
           WRITE(6,*) 'w(ip-1)=', dW_ONOS(IP), 'w(ip)=', dW_ONOS(IP+1)
           WRITE(6,*) 'FERMI',' STOP IN EFERMI'
           STOP  'FERMI - Error'
        ENDIF
        if (qprint_) write(6,*) 'TETRA ILOOP ',ILOOP
     END IF
  ENDDO
  CALL T_finish_()
  RETURN                                                           
END SUBROUTINE EFERMI

!.....................................................EFI.........
SUBROUTINE EFI(RNTOT,TOL,EFERMI,EMIN,EMAX,NP,NOS,SOS)
  IMPLICIT NONE
  REAL*8,  intent(out)  :: EFERMI
  REAL*8,  intent(inout):: TOL
  REAL*8,  intent(inout):: NOS(NP)
  REAL*8,  intent(in)   :: SOS(NP)
  INTEGER, intent(in)   :: NP
  REAL*8,  intent(in)   :: EMIN, EMAX
  REAL*8,  intent(in)   :: RNTOT
  ! locals
  REAL*8  :: NOSUP, NOSLOW, NOSIP, ADD, DNOS, ELOW, ESTEP
  INTEGER :: i, IFIND, IP, IPLOW, IPUP
  CHARACTER*67 :: ERRMSG
  ADD=0.D0
  DO I=1,NP
     ADD=ADD+SOS(I)
     NOS(I)=NOS(I)+ADD
  ENDDO
  IF(NOS(1).GT.RNTOT+.5D0*TOL.OR.NOS(NP).LT.RNTOT-.5D0*TOL) THEN
     WRITE(6,*) 'FERMI','EFERMI OUT OF ENERGY RANGE'
     WRITE(6,*) 'FERMI','STOP IN EFI'
     WRITE(6,9000) EMIN
     WRITE(6,9010) NOS(1)
     WRITE(6,9020) EMAX
     WRITE(6,9030) NOS(NP)
     WRITE(6,9040) ADD
     WRITE(6,9050) (SOS(I),I=100,1000,100)
     WRITE(6,9060) (NOS(I),I=100,1000,100)
     STOP  'FERMI - Error'
  ENDIF
  IPUP=NP
  IPLOW=1
  nosup=nos(ipup)
  noslow=nos(iplow)
  DO IFIND=1,NP
     IP=IPLOW+0.5*(IPUP-IPLOW)
     NOSIP=NOS(IP)
     IF(RNTOT-NOSIP.GT.0.d0) THEN
        IPLOW=IP
        NOSLOW=NOSIP
     ELSE
        IPUP=IP
        NOSUP=NOSIP
     END IF
     IF(IPUP.EQ.IPLOW+1) THEN
        TOL=NOSUP-NOSLOW                                                 
        ESTEP=(EMAX-EMIN)/DBLE(NP-1)                                     
        ELOW=EMIN+DBLE(IPLOW-1)*ESTEP                                    
        DNOS=NOSUP-NOSLOW                                                
        IF(DNOS.NE.0.D0) THEN                                            
           EFERMI=ELOW+(RNTOT-NOSLOW)/(NOSUP-NOSLOW)*ESTEP                
        ELSE                                                             
           EFERMI=ELOW                                                    
        END IF
        IF(EFERMI-ELOW.LT.0) WRITE(6,*) 'ERROR IN EFI '
        RETURN                                                           
     ENDIF
  ENDDO
  WRITE(6,*) 'FERMI','EFERMI NOT FOUND'
  WRITE(6,*) 'FERMI','STOP IN EFI'
  STOP  'FERMI - Error'
  
9000 FORMAT('ENERGY OF LOWER BOUND                 :',f10.5)
9010 FORMAT('NUMBER OF STATES AT THE LOWER BOUND   :',f10.5)
9020 FORMAT('ENERGY OF UPPER BOUND                 :',f10.5)
9030 FORMAT('NUMBER OF STATES AT THE UPPER BOUND   :',f10.5)
9040 FORMAT('ADD ',f10.5)
9050 FORMAT('SOS ',10f5.3)
9060 FORMAT('NOS ',10f5.5)
END SUBROUTINE EFI

!.....................................................TOTNOS......
SUBROUTINE TOTNOS(VOL,E,EMIN,EMAX,NP,NOS,SOS)                    
  ! **                                                              *
  ! **  CALCULATES THE INTEGRATED DOS                               *
  ! **  FOR ONE TETRAHEDRON ON THE ENERGYMESH                       *
  ! **  INPUT :                                                     *
  ! **    VOL         WEIGHT OF THIS TETRAHEDRON                    *
  ! **    E           ENERGIES AT THE EDGEPOINTS                    *
  ! **    EMIN        MINIMUM VALUE OF ENERGY MESH FOR NOS AND SOS  *
  ! **    EMAX        MAXIMUM VALUE OF ENERGY MESH FOR NOS AND SOS  *
  ! **    NP          NUMBER OF POINTS ON THE ENERGY MESH           *
  ! **  OUTPUT:                                                     *
  ! **    SOS         > NOS(E)+ SUM OVER E: SOS(E) =                *
  ! **    NOS         > NUMBER OF STATES BELOW E                    *
  ! **                                                              *
  IMPLICIT NONE
  REAL*8, intent(inout) :: NOS(NP), SOS(NP)
  REAL*8, intent(inout) :: E(4)
  REAL*8, intent(in)    :: VOL, EMIN, EMAX
  INTEGER, intent(in)   :: NP
  ! locals
  REAL*8  :: A, B, C, D, DE, E21, E31, E41, E32, E42, E43, EN, ESTEP, SVAR1, SVAR2, X
  INTEGER :: i, j, imax, imin
  !-----------------------------------------------------------------
  !--  INTEGRATION WITHOUT FERMISURFACE                            -
  !-----------------------------------------------------------------
  X=DMIN1(E(1),E(2),E(3),E(4))
  IF(X.GE.EMAX) THEN                                               
     RETURN                                                         
  END IF
  X=DMAX1(E(1),E(2),E(3),E(4))
  IF(X.LE.EMIN) THEN
     SOS(1)=SOS(1)+VOL
     RETURN
  END IF
  !-----------------------------------------------------------------
  !--  ORDER ENERGIES                                              -
  !-----------------------------------------------------------------
  DO I=1,3
     DO J=I+1,4
        SVAR1=DMIN1(E(I),E(J))                                           
        SVAR2=DMAX1(E(I),E(J))                                           
        E(I)=SVAR1
        E(J)=SVAR2
     ENDDO
  ENDDO
  !-----------------------------------------------------------------
  !--  CALCULATE UNCORRECTED INTEGRAL AS MEANVALUE                 -
  !-----------------------------------------------------------------
  E21=E(2)-E(1)                                                    
  if(e21.lt.1.d-10) e21=1.d-10
  E31=E(3)-E(1)
  E41=E(4)-E(1)
  E32=E(3)-E(2)
  E42=E(4)-E(2)
  E43=E(4)-E(3)
  ESTEP=(EMAX-EMIN)/DBLE(NP-1)
  !print*, 'Emin=', Emin, 'Emax=', Emax, 'Estep=', Estep, 'E(1)=', E(1), 'E(2)=', E(2), 'E(4)=', E(4)
  IMIN=IDINT(2.D0+(E(1)-EMIN)/ESTEP)                               
  IMIN=MAX0(1,IMIN)                                                
  IMAX=IDINT(1.D0+(E(2)-EMIN)/ESTEP)                               
  IMAX=MIN0(NP,IMAX)                                               
  EN=EMIN+ESTEP*(IMIN-1)
  !print*, 'IMIN=', IMIN, 'IMAX=', IMAX, 'EN=', EN
  IF(IMAX.GE.IMIN) THEN                                            
     A=VOL/(E21*E31*E41)                                            
     DO I=IMIN,IMAX
        NOS(I) = NOS(I)+A*(EN-E(1))**3                                   
        EN=EN+ESTEP
     ENDDO
  END IF
  IMIN=MAX0(1,IMAX+1)                                              
  IMAX=INT(1.D0+(E(3)-EMIN)/ESTEP)                                 
  IMAX=MIN0(NP,IMAX)                                               
  IF(IMAX.GE.IMIN) THEN                                            
     A=VOL*E21**2/(E31*E41)                                         
     B=3.D0*VOL*E21/(E31*E41)                                       
     C=3.D0*VOL/(E31*E41)                                           
     D=-VOL/(E32*E41*E31*E42)*(E31+E42)                             
     DO I=IMIN,IMAX
        DE=EN-E(2)                                                     
        NOS(I)=NOS(I)+A+DE*(B+DE*(C+D*DE))                             
        EN=EN+ESTEP                                                    
     ENDDO
  END IF
  IMIN=MAX0(1,IMAX+1)                                              
  IMAX=INT(1.D0+(E(4)-EMIN)/ESTEP)                                 
  IMAX=MIN0(NP,IMAX)                                               
  IF(E43.GT.0.D0) THEN                                             
     A=VOL                                                          
     D=VOL/(E41*E42*E43)                                            
     DO I=IMIN,IMAX
        NOS(I)=NOS(I)+A+D*(EN-E(4))**3                                 
        EN=EN+ESTEP
     ENDDO
  END IF
  IMIN=MAX0(1,IMAX+1)                                              
  IF(IMIN.GT.NP) RETURN                                            
  SOS(IMIN)=SOS(IMIN)+VOL                                          
  RETURN                                                           
END SUBROUTINE TOTNOS

! .....................................................SAMFAC......
SUBROUTINE SAMFAC(NB,NKP,EB,EF,WGHT,icor)                       
  !     **                                                              *
  !     **  CALCULATES SAMPLING WEIGHTS                                 *
  !     **  INPUT :                                                     *
  !     **    NB          NUMBER OF BANDS                               *
  !     **    NKP         NUMBER OF K-POINTS                            *
  !     **    EF          FERMI LEVEL                                   *
  !     **    W           INTEGER WORK ARRAY                            *
  !     **    NWX         LENTH OF WORK ARRAY                           *
  !     **  OUTPUT :                                                    *
  !     **    WGHT        SAMPLING WEIGHTS                              *
  !     **                                                              *
  use ReadTetra, ONLY: T_init_, Tetra1, T_finish_
  IMPLICIT NONE
  REAL*8, intent(in)    :: EB(NKP,NB), EF
  REAL*8, intent(out)   :: WGHT(NKP,NB)
  INTEGER, intent(in)   :: icor, NB, NKP
  !
  REAL*8  :: E(4), WGHT0(4), VOL, VOL0
  INTEGER :: IKP(4), I, IB, INIT, ITET, IWGHT, MWRIT, NTET, OWORK
  WGHT(:,:)=0
  CALL T_init_(VOL0,ntet,nkp)
  DO ITET=1,NTET                                               
     CALL Tetra1(itet, iwght, ikp)
     VOL=VOL0*DBLE(IWGHT)
     DO IB=1,NB
        DO I=1,4                                                     
           E(I)=EB(IKP(I),IB)                                               
        ENDDO
        WGHT0(:)=0
        CALL WEIGHT(VOL,E,EF,WGHT0,icor)
        DO I=1,4
           WGHT(IKP(I),IB)=WGHT(IKP(I),IB) + WGHT0(I)
        ENDDO
     ENDDO
  ENDDO
  CALL T_finish_()
  RETURN                                                           
END SUBROUTINE SAMFAC

!.....................................................WHEIGT......
SUBROUTINE WEIGHT(VOL,E,EF,WGHT,icor)                                 
  !**                                                              *
  !**  CALCULATES THE WEIGHTS FOR TETRAHEDRON-SAMPLING             *
  !**  CORRESPONDING TO INTEGRATION OVER ONE TETRAHEDRON           *
  !**                                                              *
  !**  CORRECTION FOR THE NONLINEAR SHAPE INCLUDED IF ICOR=1       *
  !**                                                              *
  !**  AUTHOR : P.BLOECHL                                          *
  !**                                                              *
  !**    VOL.........VOLUME OF THIS TETRAHEDRON                    *
  !**    EF..........FERMI ENERGY                                  *
  !**    D...........KT (NOT USED)                                 *
  !**    E...........ENERGIES AT THE EDGEPOINTS                    *
  !**                                                              *
  IMPLICIT NONE
  REAL*8, intent(in)    :: VOL, EF
  REAL*8, intent(inout) :: E(4)
  REAL*8, intent(out)   :: WGHT(4)
  INTEGER, intent(in)   :: icor
  ! locals
  DOUBLE PRECISION  :: FA(4),FB(4)!, cordeg
  INTEGER :: INDEX(4)
  DOUBLE PRECISION :: DA, DB, DC, DE, DE1, DE2, DE3, DE4, DOS, E21, E31, E32, E41, E42, E43, VOL14, vprime, X
  INTEGER :: i, IP, J, K, M, N
  !-----------------------------------------------------------------
  !--  INTEGRATION WITHOUT FERMISURFACE                            -
  !-----------------------------------------------------------------
  !print *, 'icor=', icor, 'cordeg=', cordeg
  X=DMIN1(E(1),E(2),E(3),E(4))                                     
  IF(X.GE.EF) THEN                                                 
     WGHT(:)=0
     RETURN                                                         
  END IF
  X=DMAX1(E(1),E(2),E(3),E(4))                                     
  IF(X.LE.EF) THEN                                                 
     VPRIME=.25D0*VOL                                               
     DO I=1,4  !10                                                  
        WGHT(I)=VPRIME                                                 
     ENDDO
     RETURN                                                         
  END IF
  !-----------------------------------------------------------------
  !--  ORDER ENERGIES                                              -
  !-----------------------------------------------------------------
  !-- INDEX HOLDS THE ORIGINAL POSITION OF THE ENERGIES AND WEIGHTS 
  DO I=1,4
     INDEX(I)=I
  ENDDO
  DO I=1,3
     IP=I                                                             
     DO J=I+1,4
        IF(E(IP).GT.E(J)) IP=J
     ENDDO
     IF(IP.GT.I) THEN
        X=E(IP)
        E(IP)=E(I)
        E(I)=X
        K=INDEX(IP)
        INDEX(IP)=INDEX(I)
        INDEX(I)=K
     END IF
  ENDDO
  !-----------------------------------------------------------------
  !--  CALCULATE UNCORRECTED INTEGRAL AS MEANVALUE                 -
  !-----------------------------------------------------------------
  E21=E(2)-E(1)                                                    
  E31=E(3)-E(1)                                                    
  E41=E(4)-E(1)                                                    
  E32=E(3)-E(2)                                                    
  E42=E(4)-E(2)                                                    
  E43=E(4)-E(3)                                                    
  WGHT(:)=0.d0
  IF(EF.GT.E(1).AND.EF.LE.E(2)) THEN                               
     DE=EF-E(1)                                                     
     VPRIME=.25D0*VOL*DE**3/(E21*E31*E41)                           
     WGHT(1)=VPRIME*(4.D0-DE/E21-DE/E31-DE/E41)                     
     WGHT(2)=VPRIME*DE/E21                                          
     WGHT(3)=VPRIME*DE/E31                                          
     WGHT(4)=VPRIME*DE/E41                                          
     !       ------  PARAMETERS FOR CORRECION                               
     DOS=3.D0*VPRIME*4.D0/(EF-E(1))                                 
  ELSE IF(EF.GT.E(2).AND.EF.LT.E(3)) THEN                          
     DE1=EF-E(1)                                                    
     DE2=EF-E(2)                                                    
     DE3=E(3)-EF                                                    
     DE4=E(4)-EF                                                    
     !       ------  TETRAHEDRON X1,X2,X13',X14'                            
     VPRIME=VOL*DE1**2/(E41*E31)*.25D0                              
     WGHT(2)=VPRIME                                                 
     WGHT(3)=VPRIME*(DE1/E31)                                       
     WGHT(4)=VPRIME*(DE1/E41)                                       
     WGHT(1)=VPRIME*(3.D0-DE1/E41-DE1/E31)                          
     !       ------  TETRAHEDRON X2,X13',X23',X14'                          
     VPRIME=.25D0*VOL*DE2*DE3*DE1/(E32*E31*E41)                     
     WGHT(1)=WGHT(1)+VPRIME*(2.D0-DE1/E31-DE1/E41)                  
     WGHT(2)=WGHT(2)+VPRIME*(2.D0-DE2/E32)                          
     WGHT(3)=WGHT(3)+VPRIME*(DE2/E32+DE1/E31)                       
     WGHT(4)=WGHT(4)+VPRIME*(DE1/E41)                               
     !       ------  TETRAHEDRON X2,X23',X24',X14'                          
     VPRIME=.25D0*VOL*DE2**2*DE4/(E42*E32*E41)                      
     WGHT(1)=WGHT(1)+VPRIME*(1.D0-DE1/E41)                          
     WGHT(2)=WGHT(2)+VPRIME*(3.D0-DE2/E32-DE2/E42)                  
     WGHT(3)=WGHT(3)+VPRIME*(DE2/E32)                               
     WGHT(4)=WGHT(4)+VPRIME*(DE2/E42+DE1/E41)                       
     !       ------  DOS=A+B*(EF-E2)+C*(EF-E2)**2                           
     DA=3.D0*VOL*E21/(E31*E41)                                      
     DB=6.D0*VOL/(E31*E41)                                          
     DC=-3.D0*VOL/(E32*E41*E31*E42)*(E31+E42)                       
     DOS=DA+DB*DE2+DC*DE2**2                                        
  ELSE IF(EF.GE.E(3).AND.EF.LT.E(4)) THEN                          
     DE=E(4)-EF                                                     
     VPRIME=.25D0*VOL*DE**3/(E41*E42*E43)                           
     VOL14=.25D0*VOL                                                
     WGHT(1)=VOL14-VPRIME*DE/E41                                    
     WGHT(2)=VOL14-VPRIME*DE/E42                                    
     WGHT(3)=VOL14-VPRIME*DE/E43                                    
     WGHT(4)=VOL14-VPRIME*(4.D0-DE/E41-DE/E42-DE/E43)               
     !       ------  PARAMETERS FOR CORRECION                               
     DOS=3.D0*VPRIME*4.D0/(E(4)-EF)                                 
  ELSE                                                             
     WRITE(6,*) 'FERMI','ERROR IN TETINT'
     STOP  'FERMI - Error'
  END IF
  !-----------------------------------------------------------------
  !--  ADD CORRECTION FOR QUADRATIC DEVIATION                      -
  !-----------------------------------------------------------------
  IF(ICOR.EQ.1) THEN                                               
     DO M=1,4                                                   
        DO N=1,4                                                   
           WGHT(M)=WGHT(M)+.25D0*(E(N)-E(M))*DOS*.1D0                     
        ENDDO
     ENDDO
  END IF
  !-----------------------------------------------------------------
  !--  REORDER WEIGHTS                                             -
  !-----------------------------------------------------------------
  DO I=1,4                                                     
     FA(INDEX(I))=WGHT(I)                                             
     FB(INDEX(I))=E(I)                                                
  ENDDO
  WGHT(:)=FA(:)
  E(:)=FB(:)                                                       
  RETURN                                                           
END SUBROUTINE WEIGHT

SUBROUTINE Read_energy_dont_store(energy_filenames, sumw0, nkpt, nmat, nume, nspin, nat)
  IMPLICIT NONE
  !
  REAL*8, intent(out) :: sumw0
  INTEGER, intent(out):: nkpt, nmat, nume, nspin
  INTEGER, intent(in) :: nat
  CHARACTER*200, intent(in) :: energy_filenames(2)
  !
  CHARACTER*10  :: KNAME
  INTEGER       :: ik, itape, ios
  INTEGER       :: N, NEn, NUM, I, ii
  REAL*8        :: EMIST, SS, TT, ZZ, wgh, E1
  !--------------------------------------------------------------------- 
  !find nkpt, nmat and nume in energy file
  open(30, FILE=energy_filenames(1),STATUS='old')
  DO I=1,NAT 
     READ(30,'(f9.5)') EMIST
     READ(30,'(f9.5)') EMIST
  ENDDO
  
  sumw0=0
  ik=0  
  nmat=0
  nume=0
  ios=0
  DO WHILE (ios == 0)
     READ(30,'(3e19.12,a10,2i6,F5.1)',IOSTAT=ios) SS,TT,ZZ,KNAME,N,NEn,wgh
     IF (ios /= 0) CYCLE
     ik=ik+1
     nmat=MAX(n,nmat)
     nume=MAX(nen,nume)
     DO ii=1,nen
        READ(30,*) NUM,E1
     ENDDO
  ENDDO
  nkpt=ik
  close(30)
  
  !***************************
  ! Check if both vector and vectordn exists. Set nspin accordingly
  NSPIN=1                                                         
  open(29, FILE=energy_filenames(2),STATUS='unknown')
  DO I=1,NAT                                                  
     READ(29,'(f9.5)',END=1004) EMIST
     READ(29,'(f9.5)',END=1004) EMIST
  ENDDO
  NSPIN=2
1004 CONTINUE

  IF(nspin.EQ.2) THEN
     ios=0
     DO WHILE (ios == 0)
        READ(29,'(3e19.12,a10,2i6,F5.1)',IOSTAT=ios) SS,TT,ZZ,KNAME,N,NEn,wgh
        nmat=MAX(n,nmat)
        nume=MAX(nen,nume)
        DO ii=1,nen
           READ(29,*) NUM,E1
        ENDDO
     END DO
     sumw0=sumw0+wgh
  END IF
  close(29)
  !print *, 'nkpt=', nkpt, 'nume=', nume
END SUBROUTINE Read_energy_dont_store

SUBROUTINE Read_energy_file(energy_filenames, Eb, nkpt, nume, nspin, nat)
  IMPLICIT NONE
  !
  REAL*8, intent(out) :: Eb(nume,nkpt,nspin)
  !INTEGER,intent(out) :: Ne(nkpt,nspin)
  CHARACTER*200, intent(in) :: energy_filenames(2)
  INTEGER, intent(in) :: nkpt, nume, nspin, nat
  !
  CHARACTER*10  :: KNAME
  INTEGER       :: ik, itape, ios, nen, N, NUM, inum, ispin, i
  REAL*8        :: EMIST, SS, TT, ZZ, wei, E1
  !--------------------------------------------------------------------- 
  
  Eb(:,:,:)=3.0d0
  do ispin=1,nspin
     itape=30-ispin+1
     open(itape, FILE=energy_filenames(ispin),STATUS='old')
     DO I=1,NAT 
        READ(itape,'(f9.5)') EMIST
        READ(itape,'(f9.5)') EMIST
     ENDDO
     
     DO ik=1,nkpt
        READ(itape,'(3e19.12,a10,2i6,F5.1)',IOSTAT=ios) SS,TT,ZZ,KNAME,N,nen,wei
        !Ne(ik,ispin) = nen
        DO inum=1,nen
           READ(itape,*) num,E1
           Eb(num,ik,ispin)=E1
        ENDDO
     ENDDO
  enddo
END SUBROUTINE Read_energy_file
