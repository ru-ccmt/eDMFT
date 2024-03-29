SUBROUTINE Cmp_ri_mat(pei, pi12lo, pe12lo, pr12lo, ilo, lmax2, nloat, lomax, ri_mat)
  IMPLICIT NONE
  !  rf1    rf2       ri_mat
  !   1      1         1.0                 RRADx * RRADx
  !   1      2         0.0                 RRADx * RADEx
  !   1      2+jlo     pi12lo(jlo,l)       RRADx * a1lo
  !   2      2         PEI(l)              RADEx * RADEx
  !   2      2+jlo     pe12lo(jlo,l)       RADEx * a1lo
  !   2+jlo  2+jlop    pr12lo(jlo,jlop,l)  a1lo  * a1lo
  INTEGER, intent(in)  :: lmax2, nloat, lomax, ilo(0:lmax2)
  REAL*8 , intent(in)  :: pei(0:lmax2)
  REAL*8 , intent(in)  :: pi12lo(1:nloat,0:lomax)
  REAL*8 , intent(in)  :: pe12lo(1:nloat,0:lomax)
  REAL*8 , intent(in)  :: pr12lo(1:nloat,1:nloat,0:lomax)
  REAL*8, intent(out)  :: ri_mat(1:2+nloat,1:2+nloat,0:lomax)
  ! locals
  INTEGER :: jlo, jlop, l
  ri_mat(:,:,:)=0.0
  do l=0,lmax2
     ri_mat(1,1,l)=1.0 ! <u|u>
     ri_mat(1,2,l)=0.0 ! <udot|u>
     ri_mat(2,1,l)=0.0 ! <u|udot>
     ri_mat(2,2,l)=pei(l) ! <udot|udot>
  enddo
  do l=0,lomax
     DO jlo=1,ilo(l)
        ri_mat(1,2+jlo,l) = pi12lo(jlo,l)  ! <u | u_lo>
        ri_mat(2+jlo,1,l) = pi12lo(jlo,l)  ! <u_lo | u>
        ri_mat(2,2+jlo,l) = pe12lo(jlo,l) ! <udot | u_lo>
        ri_mat(2+jlo,2,l) = pe12lo(jlo,l) ! <u_lo | udot>
     ENDDO
     DO jlo=1,ilo(l)
        DO jlop=1,ilo(l)
           ri_mat(2+jlo,2+jlop,l) = pr12lo(jlo,jlop,l)
        ENDDO
     ENDDO
  enddo

END SUBROUTINE Cmp_ri_mat
