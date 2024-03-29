SUBROUTINE fopen(fhp, filename, nkp, nsymop, norbitals)
  INTEGER, intent(in) :: fhp
  CHARACTER*100, intent(in) :: filename
  INTEGER, intent(out) :: nkp, nsymop, norbitals
  open(fhp, file=filename, status='old', form='unformatted')
  READ(fhp) nkp, nsymop, norbitals
END SUBROUTINE fopen

SUBROUTINE fclose(fhp)
  INTEGER, intent(in) :: fhp
  close(fhp)
END SUBROUTINE fclose

SUBROUTINE Read1(fhp, norbitals, nindo)
  INTEGER, intent(in)  :: fhp, norbitals
  INTEGER, intent(out) :: nindo(norbitals)

  INTEGER :: iorb
  READ(fhp) (nindo(iorb), iorb=1,norbitals)
  !print *, nindo
END SUBROUTINE Read1


SUBROUTINE Read2(fhp, iikp, nbands, tmaxdim2, tnorbitals, nemin)
  INTEGER, intent(in) :: fhp
  INTEGER, intent(out) :: iikp, nbands, tmaxdim2, tnorbitals, nemin

  READ(fhp) iikp, nbands, tmaxdim2, tnorbitals, nemin
  
END SUBROUTINE Read2

SUBROUTINE Read3(fhp, iisym)
  INTEGER, intent(in) :: fhp
  INTEGER, intent(out) :: iisym
  READ(fhp) iisym
END SUBROUTINE Read3

SUBROUTINE Read4(fhp, nbands, DMFTU)
  INTEGER, intent(in) :: fhp, nbands
  COMPLEX*16, intent(out) :: DMFTU(nbands)
  !locals
  INTEGER :: i
  READ(fhp) (DMFTU(i),i=1,nbands)
END SUBROUTINE Read4

SUBROUTINE fEopen(fh, filename, nat, parseline1, parseline2)
  IMPLICIT NONE
  INTEGER, intent(in) :: fh, nat
  CHARACTER*100, intent(in) :: filename
  CHARACTER*400, intent(out) :: parseline1(nat), parseline2(nat)
  ! locals
  INTEGER :: i
  open(fh,FILE=filename,STATUS='old',form='formatted')
  DO I=1,NAT
     READ(fh,'(a)') parseline1(i)
     READ(fh,'(a)') parseline2(i)
  ENDDO
END SUBROUTINE fEopen

SUBROUTINE fEclose(fh)
  IMPLICIT NONE
  INTEGER, intent(in) :: fh
  close(fh)
END SUBROUTINE fEclose

SUBROUTINE fERead1(fh, K, KNAME, wgh, ios, n0, nen)
  IMPLICIT NONE
  INTEGER, intent(in) :: fh
  REAL*8, intent(out) :: K(3), wgh
  CHARACTER*10, intent(out) :: KNAME
  INTEGER, intent(out):: ios, n0, nen
  !
  READ(fh,'(3e19.12,a10,2i6,f5.2)',IOSTAT=ios) K(1),K(2),K(3),KNAME,n0,nen,wgh
END SUBROUTINE fERead1

SUBROUTINE fERead2(fh, nen, Ek)
  INTEGER, intent(in):: fh, nen
  REAL*8, intent(out):: Ek(nen)
  INTEGER :: ii
  DO ii=1,nen
     READ(fh,*) NUM, Ek(ii)
  ENDDO
END SUBROUTINE fERead2

SUBROUTINE fWEopen(fh, filename, nat, parseline1, parseline2)
  IMPLICIT NONE
  INTEGER, intent(in) :: fh, nat
  CHARACTER*100, intent(in) :: filename
  CHARACTER*400, intent(in) :: parseline1(nat), parseline2(nat)
  ! locals
  INTEGER :: i
  open(fh,FILE=filename,STATUS='replace',form='formatted')
  DO I=1,NAT
     WRITE(fh,'(a)') parseline1(i)
     WRITE(fh,'(a)') parseline2(i)
  ENDDO
END SUBROUTINE fWEopen

SUBROUTINE fEWrite1(fh, K, KNAME, wgh, ios, n0, nen)
  IMPLICIT NONE
  INTEGER, intent(in) :: fh, n0, nen
  REAL*8, intent(in) :: K(3), wgh
  CHARACTER*10, intent(in) :: KNAME
  INTEGER, intent(out):: ios
  !
  WRITE(fh,'(3e19.12,a10,2i6,f5.2)',IOSTAT=ios) K(1),K(2),K(3),KNAME,n0,nen,wgh
END SUBROUTINE fEWrite1

SUBROUTINE fEWrite2(fh, nen, Ek)
  INTEGER, intent(in):: fh
  INTEGER, intent(in):: nen
  REAL*8, intent(in):: Ek(nen)
  INTEGER :: ii
  DO ii=1,nen
     WRITE(fh,*) ii, Ek(ii)
  ENDDO
END SUBROUTINE fEWrite2


SUBROUTINE fVopen(fh, filename, nat,Elinear)
  IMPLICIT NONE
  INTEGER, intent(in) :: fh, nat
  CHARACTER*100, intent(in) :: filename
  REAL*8, intent(out) :: Elinear(2)
  ! locals
  INTEGER :: i
  print *, filename
  open(fh,FILE=filename,STATUS='old',FORM='unformatted')
  DO I=1,NAT
     READ(fh) Elinear(1)  ! linearization energy
     READ(fh) Elinear(2)  ! linearization energy for lo
  ENDDO
END SUBROUTINE fVopen

SUBROUTINE fVclose(fh)
  IMPLICIT NONE
  INTEGER, intent(in) :: fh
  ! locals
  close(fh)
END SUBROUTINE fVclose


SUBROUTINE fVRead1(fh, K, KNAME, wgh, ios, n0, ne)
  INTEGER, intent(in) :: fh
  REAL*8, intent(out) :: K(3), wgh
  CHARACTER*10, intent(out) :: KNAME
  INTEGER, intent(out):: ios, n0,ne
  READ(fh,IOSTAT=ios) K(1),K(2),K(3),KNAME,n0,ne,wgh
END SUBROUTINE fVRead1


SUBROUTINE fVRead2(fh, n0, GS)
  INTEGER, intent(in) :: fh, n0
  INTEGER, intent(out) :: GS(3,n0)
  INTEGER :: i, ios
  READ(fh,IOSTAT=ios) (GS(1,i),GS(2,i),GS(3,i),i=1,n0)
END SUBROUTINE fVRead2

SUBROUTINE fVRead3(fh, n0, NUM, ek, A)
  INTEGER, intent(in) :: fh, n0
  INTEGER, intent(out) :: NUM
  REAL*8, intent(out) :: ek
  REAL*8, intent(out) :: A(n0)
  INTEGER :: i, ios
  READ(fh,IOSTAT=ios) NUM,ek
  READ(fh,IOSTAT=ios) (A(i),i=1,n0)
END SUBROUTINE fVRead3

SUBROUTINE fVRead3c(fh, n0, NUM, ek, A)
  INTEGER, intent(in) :: fh, n0
  INTEGER, intent(out) :: NUM
  REAL*8, intent(out) :: ek
  COMPLEX*16, intent(out) :: A(n0)
  INTEGER :: i, ios
  READ(fh,IOSTAT=ios) NUM,ek
  READ(fh,IOSTAT=ios) (A(i),i=1,n0)
END SUBROUTINE fVRead3c


!------------------

SUBROUTINE fWopen(fh, filename, nat,Elinear)
  IMPLICIT NONE
  INTEGER, intent(in) :: fh, nat
  CHARACTER*100, intent(in) :: filename
  REAL*8, intent(in) :: Elinear(2)
  ! locals
  INTEGER :: i
  open(fh,FILE=filename,STATUS='replace',FORM='unformatted')
  DO I=1,NAT
     WRITE(fh) Elinear(1)  ! linearization energy
     WRITE(fh) Elinear(2)  ! linearization energy for lo
  ENDDO
END SUBROUTINE fWopen

SUBROUTINE fWclose(fh)
  IMPLICIT NONE
  INTEGER, intent(in) :: fh
  ! locals
  close(fh)
END SUBROUTINE fWclose

SUBROUTINE fVWrite1(fh, K, KNAME, wgh, n0, ne, ios)
  INTEGER, intent(in) :: fh
  REAL*8, intent(in) :: K(3), wgh
  CHARACTER*10, intent(in) :: KNAME
  INTEGER, intent(in):: n0,ne
  INTEGER, intent(out):: ios
  WRITE(fh,IOSTAT=ios) K(1),K(2),K(3),KNAME,n0,ne,wgh
END SUBROUTINE fVWrite1

SUBROUTINE fVWrite2(fh, n0, GS)
  INTEGER, intent(in) :: fh, n0
  INTEGER, intent(in) :: GS(3,n0)
  INTEGER :: i, ios
  WRITE(fh,IOSTAT=ios) (GS(1,i),GS(2,i),GS(3,i),i=1,n0)
END SUBROUTINE fVWrite2

SUBROUTINE fVWrite3(fh, n0, NUM, ek, A)
  INTEGER, intent(in) :: fh, n0
  INTEGER, intent(in) :: NUM
  REAL*8, intent(in) :: ek
  REAL*8, intent(in) :: A(n0)
  INTEGER :: i, ios
  WRITE(fh,IOSTAT=ios) NUM,ek
  WRITE(fh,IOSTAT=ios) (A(i),i=1,n0)
END SUBROUTINE fVWrite3

SUBROUTINE fVWrite3c(fh, n0, NUM, ek, A)
  INTEGER, intent(in) :: fh, n0
  INTEGER, intent(in) :: NUM
  REAL*8, intent(in) :: ek
  COMPLEX*16, intent(in) :: A(n0)
  INTEGER :: i, ios
  WRITE(fh,IOSTAT=ios) NUM,ek
  WRITE(fh,IOSTAT=ios) (A(i),i=1,n0)
END SUBROUTINE fVWrite3c




