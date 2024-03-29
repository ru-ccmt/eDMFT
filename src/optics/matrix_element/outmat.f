SUBROUTINE OUTMAT(OUTME,Nbvalmax)
  use opme
  use bindex
  use xa3
  !ad ___________________ CALCULATION OF MATRIX ELEMENTS _________________
  INCLUDE 'param.inc'
  IMPLICIT REAL*8 (A-H,O-Z) 
  LOGICAL     INFO_FLAG
  LOGICAL     REL,IM,LSO,SPIN,MME_FLAG
  CHARACTER*3  OUTME
  CHARACTER*10 KNAME
  COMPLEX*16  O(3),OX1,OX2,OY1,OY2,OZ2,OZ1     
  REAL*8  :: ss1, ss2
  INTEGER :: iss
  COMMON /outop/ Ncol,icol(9)          
  COMMON /KPOI/ S,T,Z,NEMIN,NEMAX,KKZ,N,NNLO,KNAME
  COMMON /LEAD/ KFIRST,KLAST,KEND,KSTEP,KOUT,KSTOP
  COMMON /COM/  EMIN,EMAX,ELECN,EULIMIT,EDLIMIT,NK,IOUT,NSPIN,NAT,NBAND,ix,NB(NKPT),MINWAV,MAXWAV
  COMMON /MIM / MIMA(2) 
  COMMON /CLOGIC/ LSO,SPIN,REL,MME_FLAG
  NEMIN=MIMA(1)
  NEMAX=MIMA(2)
  INFO_FLAG=.FALSE.
  IF (MME_FLAG) THEN
     WRITE(24,9010) NK,NEMIN,NEMAX,EMIN,EMAX,KNAME
  ELSE
     WRITE(3,9010) NK,NEMIN,NEMAX,EMIN,EMAX,KNAME
     if(OUTME.EQ.'ON ') then
        WRITE(4,9010) NK,NEMIN,NEMAX,EMIN,EMAX,KNAME
     endif
     WRITE(9,9010) NK,NEMIN,NEMAX,EMIN,EMAX,KNAME
  END IF
  NBINDEX=0

  ss1=0.d0
  ss2=0.d0
  iss=0
  IF (LSO.AND.(.NOT.SPIN)) THEN
     IF (MME_FLAG) THEN ! This is only for simple optics
        write(6,*) 'SPIN-orbit coupling for systems without inversion requires spinpolarized setup. Contact authors'
        stop 'SO requires spinpol.setup'
        INFO_FLAG=.TRUE.
        DO  NB1=NEMIN,Nbvalmax
           DO  NB2=NB1,NEMAX
              N1=NIN(NB1,NB2)
              IF (NB2.NE.NB1+1) THEN
                 WRITE(24,9030) NB1,NB2,DBLE(OPMATX(N1)),AIMAG(OPMATX(N1)),DBLE(OPMATY(N1)),AIMAG(OPMATY(N1)),DBLE(OPMATZ(N1)),AIMAG(OPMATZ(N1))  
              ELSE
                 WRITE(24,9030) NB1,NB2,0.,0.,0.,0.,0.,0.
              END IF
           END DO
        END DO
     ELSE
        !kh This symmetrization is because for time invariant case and SO we set iso to 1 in the main part of the program,
        !kh  hence only half of the bands were properly computed. We need to symmetrized the matrix elements here.
        !kh  We would not need this if iso was set to 2 when SO is present.
        DO NB1=NEMIN,Nbvalmax,2
           DO NB2=NB1,NEMAX,2
              if(nb2.gt.nb1) then
                 N1=NIN(NB1,NB2)
                 N2=NIN(NB1+1,NB2+1)
                 N3=NIN(NB1,NB2+1)
                 N4=NIN(NB1+1,NB2)
                 OX1=(OPMATX(n1)+conjg(OPMATX(n2)))
                 OY1=(OPMATY(n1)+conjg(OPMATY(n2)))
                 OZ1=(OPMATZ(n1)+conjg(OPMATZ(n2)))
                 OX2=(OPMATX(n3)-conjg(OPMATX(n4)))
                 OY2=(OPMATY(n3)-conjg(OPMATY(n4)))
                 OZ2=(OPMATZ(n3)-conjg(OPMATZ(n4)))
                 !ad
                 OPMATX(n1)=OX1
                 OPMATX(n2)=OX1
                 OPMATX(n3)=OX2
                 OPMATX(n4)=OX2
                 OPMATY(n1)=OY1
                 OPMATY(n2)=OY1
                 OPMATY(n3)=OY2
                 OPMATY(n4)=OY2
                 OPMATZ(n1)=OZ1
                 OPMATZ(n2)=OZ1
                 OPMATZ(n3)=OZ2
                 OPMATZ(n4)=OZ2
                 !ad
              else
                 !ad
                 N1=NIN(NB1,NB2)
                 N2=NIN(NB1+1,NB2+1)
                 N3=NIN(NB1,NB2+1)
                 OX1=OPMATX(n1)+OPMATX(n2)
                 OY1=OPMATY(n1)+OPMATY(n2)
                 OZ1=OPMATZ(n1)+OPMATZ(n2)
                 !ad

                 ss1 = ss1 + abs(OPMATX(n1)/OX1)
                 ss2 = ss2 + abs(OPMATX(n2)/OX1)
                 iss = iss + 1
                 
                 OPMATX(n1)=OX1
                 OPMATX(n2)=OX1
                 OPMATX(n3)=0.0d0
                 OPMATY(n1)=OY1
                 OPMATY(n2)=OY1
                 OPMATY(n3)=0.0d0
                 OPMATZ(n1)=OZ1
                 OPMATZ(n2)=OZ1
                 OPMATZ(n3)=0.0d0
              endif
           END DO
        END DO
     END IF
  END IF
  
  DO NB1=NEMIN,nemax   ! 119
     DO NB2=NB1,NEMAX  ! 119
        NBINDEX=NBINDEX+1  
        O(1) = OPMATX(nbindex)
        O(2) = OPMATY(nbindex)
        O(3) = OPMATZ(nbindex)
        !ad   write out complex matrix elements before symmetrization
        if(OUTME.EQ.'ON ') write(4,9040) NB1,NB2,O(1),O(2),O(3),E(NB2)-E(NB1)
        IF ((MME_FLAG).AND.(.NOT.INFO_FLAG)) THEN
           WRITE(24,9030) NB1,NB2,(DBLE(O(i)),AIMAG(O(i)),i=1,Ncol/2)
        END IF
        IF (.NOT.MME_FLAG) THEN
           call outsym(o,nb1,nb2)
        END IF
     ENDDO
  ENDDO
  
  RETURN
9010 FORMAT(/,2X,' KP:',I6,' NEMIN NEMAX : ',2I5,' dE:',2f5.1,' K:',a10 /)
9030 FORMAT(3X,2I4,6E13.6)
9040 FORMAT(3X,2I4,6E13.6,F13.8)
END SUBROUTINE OUTMAT

SUBROUTINE OUTSYM(O,nb1,nb2)
  !ad  ________________ SYMMETRIZATION OF MATRIX ELEMENTS ________________                  
  !ad   May 1999: updated by Jan Kunes
  IMPLICIT REAL*8 (A-H,O-Z)
  INCLUDE 'param.inc'
  COMMON /outop/  Ncol,icol(9)
  COMMON /SYM2/   TAU(3,NSYM),IORD,IMAT(3,3,NSYM)
  COMMON /SYMo/   opimat(3,3,NSYM)
  COMMON /SYMd/   det(NSYM)
  COMPLEX*16  O(3),o2(6),o1(3)
  DIMENSION outm(9)
  syfac=1.D0/REAL(iord)
  do ii=1,6
     o2(ii)=0.0
  end do
  do l=1,iord   ! 772
     do ii=1,3
        o1(ii)=0.0
     end do
     do   i=1,3
        do  ii=1,3
           o1(i)=o1(i)+opIMAT(ii,i,l)*O(ii)
        end do
     end do
     o2(1)=o2(1)+o1(1)*conjg(o1(1))
     o2(2)=o2(2)+o1(2)*conjg(o1(2))
     o2(3)=o2(3)+o1(3)*conjg(o1(3))
     if (det(l).LT.0) then
        o2(4)=o2(4)+conjg(o1(2))*o1(1)
        o2(5)=o2(5)+conjg(o1(3))*o1(1)
        o2(6)=o2(6)+conjg(o1(3))*o1(2)
     else
        o2(4)=o2(4)+o1(2)*conjg(o1(1))
        o2(5)=o2(5)+o1(3)*conjg(o1(1))
        o2(6)=o2(6)+o1(3)*conjg(o1(2))
     end if
  enddo
  !____________________________________________________
  outm(1)=Dble(o2(1))*syfac
  outm(2)=Dble(o2(2))*syfac
  outm(3)=Dble(o2(3))*syfac        
  outm(4)=Dble(o2(4))*syfac
  outm(5)=Dble(o2(5))*syfac
  outm(6)=Dble(o2(6))*syfac
  outm(7)=aIMAG(o2(4))*syfac
  outm(8)=aIMAG(o2(5))*syfac
  outm(9)=aIMAG(o2(6))*syfac
  !____________________________________________________
  if(nb1.eq.nb2) then
     write(9,9030)  NB1,NB2,(outm(icol(i)),i=1,Ncol)
  endif
  WRITE(3,9030)  NB1,NB2,(outm(icol(i)),i=1,Ncol)
9030 FORMAT(3X,2I4,9E13.6)
23 continue
END SUBROUTINE OUTSYM
