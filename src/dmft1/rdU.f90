SUBROUTINE fopen(fhp, filename, nkp, nsymop, norbitals)
  INTEGER, intent(in) :: fhp
  CHARACTER*100, intent(in) :: filename
  INTEGER, intent(out) :: nkp, nsymop, norbitals
  
  
  !print *, fhp, filename
  open(fhp, file=filename, status='old', form='unformatted')

  READ(fhp) nkp, nsymop, norbitals
  !print *, nkp, nsymop, norbitals

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

SUBROUTINE fEopen(fh, filename, nat)
  IMPLICIT NONE
  INTEGER, intent(in) :: fh, nat
  CHARACTER*100, intent(in) :: filename
  ! locals
  INTEGER :: i
  REAL*8 :: EMIST
  open(fh,FILE=filename,STATUS='old')
  DO I=1,NAT
     READ(fh,'(f9.5)') EMIST
     READ(fh,'(f9.5)') EMIST
  ENDDO
  !print *, EMIST
END SUBROUTINE fEopen

SUBROUTINE fEclose(fh)
  IMPLICIT NONE
  INTEGER, intent(in) :: fh
  close(fh)
END SUBROUTINE fEclose

SUBROUTINE fERead1(fh, K, KNAME, wgh, ios, nen)
  IMPLICIT NONE
  INTEGER, intent(in) :: fh
  REAL*8, intent(out) :: K(3), wgh
  CHARACTER*10, intent(out) :: KNAME
  INTEGER, intent(out):: ios, nen
  !
  INTEGER :: n
  READ(fh,'(3e19.12,a10,2i6,f5.2)',IOSTAT=ios) K(1),K(2),K(3),KNAME,n,nen,wgh
  !print *, 'wgh=', wgh, 'ios=', ios
END SUBROUTINE fERead1

SUBROUTINE fERead2(fh, nen, Ek)
  INTEGER, intent(in):: fh, nen
  REAL*8, intent(out):: Ek(nen)
  INTEGER :: ii
  DO ii=1,nen
     READ(fh,*) NUM, Ek(ii)
  ENDDO
END SUBROUTINE fERead2

SUBROUTINE fVopen(fh, filename, nat)
  IMPLICIT NONE
  INTEGER, intent(in) :: fh, nat
  CHARACTER*100, intent(in) :: filename
  ! locals
  INTEGER :: i
  REAL*8 :: EMIST
  print *, filename
  open(fh,FILE=filename,STATUS='old',FORM='unformatted')
  DO I=1,NAT
     READ(fh) EMIST
     READ(fh) EMIST
  ENDDO
  !print *, EMIST
END SUBROUTINE fVopen


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
