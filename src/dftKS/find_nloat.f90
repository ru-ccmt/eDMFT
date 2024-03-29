INTEGER Function find_nloat(fh)  
  ! Finds maximum number of local orbitals from file case.in1
  IMPLICIT NONE
  INTEGER, intent(in) :: fh
  !locals
  INTEGER :: nloat_max, i, l, l_old, il, nl
  REAL*8  :: a
  !
  nloat_max=0
  read(fh,*)
  read(fh,*)
  DO   ! atom loop
     i=0
     l_old=9999
     read(fh,*,err=20,end=20) a,nl  ! a==GLOBAL E-PARAMETER, nl == n OTHER CHOICES
     do il=1,nl
        read(fh,*) l
        if(l.eq.l_old) then
           i=i+1
        else
           nloat_max=max(nloat_max,i)
           i=1
           l_old=l
        endif
     enddo
     nloat_max=max(nloat_max,i)
  ENDDO
20 continue
  nloat_max=nloat_max+1
  if(nloat_max.lt.3) nloat_max=3
  rewind fh
  find_nloat = nloat_max
  return
end Function find_nloat
