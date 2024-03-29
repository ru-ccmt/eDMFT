SUBROUTINE NN(NAT)
  use reallocate, only: doreallocate  
  use structure, only: lattic, ALAT, ALPHA, ortho, pos, Rmt, mult
  use param
  use mpi, only: stop_MPI
  IMPLICIT NONE
  INTEGER, intent(in) ::  NAT
  !INCLUDE 'param.inc'
  INTEGER            MM, NNN, N1, N2, N3, NC, nnn1
  INTEGER            K, I3, JAT, L, JATOM, INDEX, M, I2, I1, J
  DOUBLE PRECISION   DIST,rlarge
  DOUBLE PRECISION   PI, SINGAM, COSGAM, SUMRAD
  CHARACTER*67       ERRMSG
  PARAMETER          (NNN=10000)
  INTEGER, pointer :: NR(:), NNAT(:)      ! nnn
  DOUBLE PRECISION, pointer :: DISTS(:),PNN(:,:)    !(3,NNN)
  DOUBLE PRECISION   DIF(3),XX(3),PP(3),P(3), HELP(3), BRNN(3,3)
  !-----------------------------------------------------------------------
  !     INDEX COUNTS ALL ATOMS IN THE UNIT CELL, JATOM COUNTS ONLY THE
  !     NONEQUIVALENT ATOMS
  !.....SET UP LATTICE, AND IF REQUIRED ROTATION MATRICES

  CALL DIRECT_LATTICE(brnn, nat, alpha, lattic)

  !print *, 'alat=', alat
  !print *, 'alpha=', alpha
  !print *, 'ortho=', ortho
  !print *, 'ndf=', ndf
  !print *, 'pos=', pos
  !print *, 'Rmt=', Rmt
  !print *, 'mult=', mult
  pi=4.d0*atan(1.d0)
  cosgam=cos(alpha(3))
  singam=sin(alpha(3))
  INDEX=0
  rlarge=min(alat(1),alat(2))
  rlarge=min(rlarge,alat(3))+0.01
  !
  allocate( NR(NNN),nnat(nnn),dists(nnn),pnn(3,nnn) )
  nnn1=nnn
  INDEX=0
  DO JATOM=1,NAT         ! 200
     DO M=1,MULT(JATOM)  ! 190
        INDEX=INDEX+1
        DO J=1,3
           XX(J)=POS(J,INDEX)
        ENDDO
        NC=0
        DO I1=-2,2        ! 180
           DO I2=-2,2     ! 180
              DO I3=-2,2  ! 180
                 IF(ortho) THEN
                    P(1)=I1*BRnn(1,1)+I2*BRnn(2,1)+I3*BRnn(3,1)
                    P(2)=I1*BRnn(1,2)+I2*BRnn(2,2)+I3*BRnn(3,2)
                    P(3)=I1*BRnn(1,3)+I2*BRnn(2,3)+I3*BRnn(3,3)
                 ELSE
                    P(1)=I1
                    P(2)=I2
                    P(3)=I3
                    IF(LATTIC(1:3).eq.'CXZ') THEN
                       P(1)=I1*0.5d0+i3*0.5d0
                       P(2)=I2
                       P(3)=-I1*0.5d0+i3*0.5d0
                    END IF
                 ENDIF
                 K=0
                 DO JAT=1,NAT   ! 120
                    DO MM=1,MULT(JAT)
                       K=K+1
                       DIST=0.d0
                       DO L=1,3  ! 100
                          PP(L)=POS(L,K)+P(L)
                          DIF(L)=XX(L)-PP(L)
                       ENDDO  ! 100
                       IF (.not.ortho) THEN
                          help(1)=dif(1)
                          help(2)=dif(2)
                          help(3)=dif(3)
                          if(lattic(1:1).eq.'R') then
                             dif(1)=help(1)*BRnn(1,1)+help(2)*BRnn(2,1)+help(3)*BRnn(3,1)
                             dif(2)=help(1)*BRnn(1,2)+help(2)*BRnn(2,2)+help(3)*BRnn(3,2)
                             dif(3)=help(1)*BRnn(1,3)+help(2)*BRnn(2,3)+help(3)*BRnn(3,3)
                          elseif(lattic(1:3).eq.'CXZ') then
                             dif(1)=help(1)*singam
                             dif(2)=(help(1)*cosgam*alat(1)+help(2)*alat(2))/alat(2)
                             dif(3)=help(3)
                          else
                             dif(1)=(help(1)*BRnn(1,1)*ALAT(1)+help(2)*BRnn(2,1)*ALAT(2)+help(3)*BRnn(3,1)*ALAT(3))/ALAT(1)
                             dif(2)=(help(1)*BRnn(1,2)*ALAT(1)+help(2)*BRnn(2,2)*ALAT(2)+help(3)*BRnn(3,2)*ALAT(3))/ALAT(2)
                             dif(3)=(help(1)*BRnn(1,3)*ALAT(1)+help(2)*BRnn(2,3)*ALAT(2)+help(3)*BRnn(3,3)*ALAT(3))/ALAT(3)
                          endif
                       ENDIF
                       DO L=1,3  ! 103
                          DIST=DIST+DIF(L)*DIF(L)*ALAT(L)*ALAT(L)
                       ENDDO     ! 103
                       DIST=SQRT(DIST)
                       IF(DIST.GT.rlarge) CYCLE !GO TO 110
                       IF(DIST.LT..001) CYCLE   !GO TO 110
                       NC=NC+1
                       if(nc.gt.nnn1) then
                          !       goto 900
                          nnn1=nnn1*5
                          call doreallocate( nr, nnn1)
                          call doreallocate( nnat, nnn1)
                          call doreallocate( dists, nnn1)
                          call doreallocate( pnn, 3, nnn1)
                       endif
                       DISTS(NC)=DIST
                       NNAT(NC)=JAT
                       DO L=1,3    ! 105
                          PNN(L,NC)=PP(L)
                       ENDDO    ! 105
                    ENDDO  !110                    CONTINUE
                 ENDDO ! 120
              ENDDO ! 180
           ENDDO ! 180
        ENDDO ! 180
        !180     CONTINUE
        CALL ORD2(DISTS,NR,NC)
        N1=1
        N2=NR(N1)
        N3=NNAT(N2)
        SUMRAD=RMT(JATOM)+RMT(N3)
        IF(M.EQ.1) THEN
           IF(SUMRAD.GT.DISTS(N1)) GOTO 910
        ENDIF
        !
        DO N1=1,NC   ! 185
           N2=NR(N1)
           N3=NNAT(N2)
           SUMRAD=RMT(JATOM)+RMT(N3)
           IF(SUMRAD.GT.DISTS(N1)) GOTO 910
        ENDDO        ! 185
     ENDDO  ! 190     CONTINUE
  ENDDO ! 200     CONTINUE
  !
  deallocate(nr,nnat,dists,pnn)
  RETURN
  !        Error messages
  !
!900 CALL OUTERR('NN','nnn too small')
!  STOP 'NN - Error'
910 CALL OUTERR('NN','overlapping spheres')
  WRITE (ERRMSG,9000) JATOM,RMT(JATOM),N3,RMT(N3)
  CALL OUTERR('NN',ERRMSG)
  WRITE (ERRMSG,9010) SUMRAD, DISTS(N1)
  CALL OUTERR('NN',ERRMSG)
  print*, 'NN','overlapping spheres'
  call stop_MPI
  STOP 'NN - Error'
9000 FORMAT('RMT(',I2,')=',F7.5,' AND RMT(',I2,')=',F7.5)
9010 FORMAT('SUMS TO',F8.5,' GT NNN-DIST=',F8.5)
END SUBROUTINE NN
!


SUBROUTINE ORD2(A,NR,IMAX)
  !     ORDERS ELEMENTS IN ARRAY A INCREASING IN SIZE
  !       REORDERS CORRESPONDING INDICES (NR)
  IMPLICIT NONE
  DOUBLE PRECISION A(*)
  INTEGER NR(*)
  INTEGER IMAX
  LOGICAL CONT
  INTEGER I, NHLP
  DOUBLE PRECISION HLP
  
  DO I=1,IMAX 
     NR(I)=I
  ENDDO
  
100 CONTINUE
  I=1
  CONT=.FALSE.
110 I=I+1
  IF(A(I).LT.A(I-1))THEN
     !       INTERCHANGE I AND (I-1) ELEMENT IN ARRAYS A AND NR
     HLP=A(I)
     A(I)=A(I-1)
     A(I-1)=HLP
     NHLP=NR(I)
     NR(I)=NR(I-1)
     NR(I-1)=NHLP
     CONT=.TRUE.
  ENDIF
  IF(I.LT.IMAX)GO TO 110
  IF(CONT) GO TO 100
  RETURN
END SUBROUTINE ORD2
