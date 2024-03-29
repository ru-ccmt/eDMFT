SUBROUTINE LOPW(NAT)
  use matrices, only: KZZ, Kn
  use lolog, only : ilo
  use lstapw, only  : NV
  use structure, only: pos, mult, ndf, ROTLOC, rotij
  use mpi, only: Qprint, stop_MPI
  use param
  IMPLICIT NONE
  !INCLUDE 'param.inc'
  !        Scalar Arguments
  INTEGER            NAT
  !..................................................................
  !   generates the LAPW (K+G)-vector for local orbitals
  !..................................................................
  !        Locals
  INTEGER            IA1, IEQ, IIX, INDEX, J, K, KOFF, L, LM, LMDN
  INTEGER            LMUP, LMX, N, NATX, NATXX, NB, NBM
  INTEGER            JLO,ipass
  DOUBLE PRECISION   HL, RKGM, SX, TPI, check
  DOUBLE PRECISION   ROTV1(3), ROTV2(3), VEC(3)
  COMPLEX*16         CC
  COMPLEX*16         HH((2*LOMAX+1)*48,(2*LOMAX+1)*48)
  COMPLEX*16         SF(NDF), YL(0:(LOMAX+1)**2) !,nv:HSROWS)
  !        External Subroutines
  EXTERNAL           ROTATE, YLM
  !        Intrinsic Functions
  INTRINSIC          ATAN, DCMPLX, DCONJG, EXP, SQRT
  !     ** Maybe Experiment **
  DOUBLE PRECISION  VEC2(3), TMP1, TMP2
  !
  TPI = 8.0D+0*ATAN(1.0D+0)
  !
  check=2.0D-2
  ipass=0
  
1 continue
  
  check=check/2.d0
  KOFF = NV
  IA1 = 0
  DO N = 1, NAT   ! 140
     DO L = 0, LOMAX
        do jlo=1,ilo(l,n)
           LMDN = L*L + 1
           LMUP = (L+1)*(L+1)
           INDEX = 0
           NB = 0
           NBM = MULT(N)*(1+LMUP-LMDN)
           DO IEQ = 1, MULT(N)     ! 120
              DO LM = LMDN, LMUP   ! 110
                 NB = NB + 1
                 K = KOFF + NB
10               CONTINUE
                 INDEX = INDEX + 1
                 IF (INDEX .GT. NV) GOTO 900
                 !                  WRITE (6,*) 'INDEX,K,N,L,IEQ,LM',INDEX,K,N,L,IEQ,LM
                 KZZ(1,K) = KZZ(1,INDEX)
                 KZZ(2,K) = KZZ(2,INDEX)
                 KZZ(3,K) = KZZ(3,INDEX)
                 Kn(:,K) = Kn(:,index)
                 RKGM = sqrt(Kn(1,K)**2+Kn(2,K)**2+Kn(3,K)**2)
                 IF (NBM .NE. 1) THEN
                    DO NATX = 1, MULT(N)
                       NATXX = IA1 + NATX
                       SX = KZZ(1,K)*POS(1,NATXX) + KZZ(2,K)*POS(2,NATXX) + KZZ(3,K)*POS(3,NATXX)
                       SF(NATX) = DCMPLX(DCOS(TPI*SX),DSIN(TPI*SX))
                    ENDDO
                    IIX = 0
                    DO NATX = 1, MULT(N)
                       IF (RKGM .LE. 1.0D-5) THEN
                          DO LMX = LMDN, LMUP
                             YL(LMX-1) = 0.0D0
                          ENDDO
                          YL(0) = 1.D0 
                       ELSE
                          ROTV1 = matmul(ROTIJ(:,:,IA1+NATX), Kn(:,K))
                          ROTV2 = matmul(ROTLOC(:,:,N), ROTV1)
                          CALL YLM(ROTV2,LOMAX,YL(0))
                       ENDIF
                       DO LMX = LMDN, LMUP
                          IIX = IIX + 1
                          HH(IIX,NB) = SF(NATX)*YL(LMX-1) !,K)
                       ENDDO
                    ENDDO
                    IF (NB .NE. 1) THEN
                       DO J = 1, NB - 1
                          CC = (0.0D+0,0.0D+0)
                          DO IIX = 1, NBM
                             CC = CC + HH(IIX,NB)*DCONJG(HH(IIX,J))
                          ENDDO
                          DO IIX = 1, NBM
                             HH(IIX,NB) = HH(IIX,NB) - CC*HH(IIX,J)
                          ENDDO
                       ENDDO
                    ENDIF
                    HL = 0.0D+0
                    DO IIX = 1, NBM
                       HL = HL + dble( DCONJG(HH(IIX,NB))*HH(IIX,NB) )
                    ENDDO
                    !
                    !       Change here to increase test so now it is RMS > 0.1D0
                    !
                    IF (HL .LE. dble(NBM)*check) GOTO 10  !! PB Change here, was 1.0D-3
                    !        WRITE (6,6001) n,l,ieq,index,K, RKGM, (KZZ(J,K),J=1,3),hl
!6001                format(' atom',i3,' L',i2,' equiv',i2,' PW',i5,i5,f10.4,3i4,e15.5)
                    !
                    !       *** Maybe Experiment **
                    !       Ensure that the K vectors themselves are not linearly dependent
                    !       This can be done quicker by storing the VEC2 values....
                    IF(NB .GT. 2) THEN
                       TMP1=KZZ(1,K)*KZZ(1,K)+KZZ(2,K)*KZZ(2,K)+KZZ(3,K)*KZZ(3,K)
                       IF(TMP1 .lt. 1D-15)TMP1=1.D0
                       VEC(1:3)=KZZ(1:3,K)/sqrt(TMP1)
                       DO J = 1, NB-1
                          IIX=KOFF+J
                          TMP2=KZZ(1,IIX)*KZZ(1,IIX)+KZZ(2,IIX)*KZZ(2,IIX)+KZZ(3,IIX)*KZZ(3,IIX)
                          IF(TMP2 .lt. 1D-15)TMP2=1.D0
                          VEC2(1:3)=KZZ(1:3,IIX)/sqrt(TMP2)
                          if((abs(VEC(1)-VEC2(1)).lt. 1D-10).and.( abs(VEC(2)-VEC2(2)).lt. 1D-10).and.( abs(VEC(3)-VEC2(3)).lt. 1D-10) ) GOTO 10
                          if((abs(VEC(1)+VEC2(1)).lt. 1D-10).and.( abs(VEC(2)+VEC2(2)).lt. 1D-10) .and.( abs(VEC(3)+VEC2(3)).lt. 1D-10) ) GOTO 10
                       ENDDO
                    ENDIF
                    !       *** End of Experiment **
                 ELSE
                    CYCLE !GOTO 110
                 ENDIF
                 HL = 1.0D+0/SQRT(HL)
                 DO IIX = 1, NBM  ! 100
                    HH(IIX,NB) = HH(IIX,NB)*HL
                 ENDDO            ! 100
              ENDDO  !110                 CONTINUE
           ENDDO  ! 120          CONTINUE
           KOFF = KOFF + MULT(N)*(2*L + 1)
        ENDdo
     ENDDO
     IA1 = IA1 + MULT(N)
  ENDDO !140      CONTINUE
  
!6000 FORMAT (5X,I5,F15.8,3I5,10X,3F10.5,5X)
  !
  RETURN
900 if(ipass.lt.10) then
     ipass=ipass+1
     if (Qprint) write(6,901) n,ipass,check
     if (Qprint) write(21,901) n,ipass,check
901  format(':WAR   : LOPW-exhausted for atom',i5,' PASS',i2,'  had to reduce check',f9.6) 
     goto 1
  endif
  CALL OUTERR('LOPW','Plane waves exhausted ')
  call stop_MPI
  STOP 'LOPW - Error'
END SUBROUTINE LOPW
