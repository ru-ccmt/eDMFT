SUBROUTINE TestXqtlDiff(xqtl2, E, l1, nbands, iso, nume, maxdim2)
  ! Checks that xqtl computed through DMFT and the original calculation gives the same result
  ! BE CAREFULL : This works only if split=88 (no local transformation)
  IMPLICIT NONE
  COMPLEX*16, intent(in) :: xqtl2(nbands,maxdim2,maxdim2)
  REAL*8, intent(in)     :: E(nume)
  INTEGER, intent(in)    :: l1, nbands, iso, nume, maxdim2
  ! local variables
  INTEGER :: lcase, nind, is1, m1, lms1, ind1, is2, m2, lms2, ind2, N1, num
  REAL*8  :: small
  INTEGER :: ks(2)
  ks(1)=1
  ks(2)=-1
  small=1e-8
  N1=2*l1+1
  nind=N1*iso
  do num=1,nbands
     WRITE(210,904)  num, E(num), 0.0
     do is1=1,iso           ! over spin-1
        do m1=-l1,l1        ! over m-1
           lms1=l1+1+m1
           ind1=l1+1+m1+N1*(is1-1)
           do is2=1,iso        ! over spin-2
              do m2=-l1,l1     ! over m-2
                 lms2=l1+1+m2
                 ind2=l1+1+m2+N1*(is2-1)
                 if (abs(xqtl2(num,ind1,ind2))>1e-3) then
                    WRITE(210,931) xqtl2(num,ind1,ind2), l1, m1, ks(is1), l1, m2, ks(is2), ind1, ind2
                 endif
              enddo
           enddo
        enddo
     enddo
  enddo
931 FORMAT(1X,'*****:',2f17.10,8i4)                               
904 FORMAT(1X,' BAND #',I3,'  E=',F9.5,'  WEIGHT=',F10.7)             
END SUBROUTINE TestXqtlDiff
