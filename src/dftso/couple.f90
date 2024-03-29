subroutine couple(indj,a,b,c)
  !
  !.....THIS IS THE ROUTINE CALCULATES THE MATRIX ELEMENT OF
  !     <Ylm' Xs'|s*l|Ylm Xs>. IT GENERATES THE ARRAY COUP WHICH WILL
  !     BE USED IN THE SPIN-ORIBIT PROGRAM.   L=0 DOESN'T COUPLE.
  !
  USE param, ONLY: labc
  USE couples, ONLY: couplo
  IMPLICIT NONE
  INTEGER, intent(in) :: indj
  REAL*8, intent(in)  :: a, b, c
  ! locals
  COMPLEX*16 :: COUP(LABC,-LABC:LABC,-LABC:LABC,2,2), T(2,2), sum
  INTEGER :: i, is1, is2, isf, isi, j, l, mf, mi
  COMPLEX*16, PARAMETER :: imag = (0.d0,1.d0)
  !!  spin rotation has the form     
  !!
  !! T =  | e^{ i(a+c)/2}*cos(b/2) , e^{ i(a-c)/2}*sin(b/2) |
  !!      |-e^{-i(a-c)/2}*sin(b/2) , e^{-i(a+c)/2}*cos(b/2) |
  !! which is the following sequence of rotations:
  !!
  !!                T =  Rz(a) * Ry(b) * Rz(c)
  !!
  T(1,1)=exp(imag*(a+c)/2.d0)*cos(b/2.d0)
  T(1,2)=exp(imag*(a-c)/2.d0)*sin(b/2.d0)
  T(2,1)=-dconjg(T(1,2))
  T(2,2)=dconjg(T(1,1))

  !write(6,*)'COUPLE:',a,b,c
  !write(6,*)
  !write(6,*)t(1,1),t(1,2)
  !write(6,*)t(2,1),t(2,2)
  !write(6,*)

  do isi=1,2
     do isf=1,2
        do l=1,labc
           do mi=-labc,labc
              do mf=-labc,labc
                 coup(l,mi,mf,isi,isf)=(0.d0,0.d0)
                 couplo(indj,l,mi,mf,isi,isf)=(0.d0,0.d0)
              enddo
           enddo
        enddo
     enddo
  enddo
  
  !
  DO L=1,LABC
     DO MF=-L,L
        !
        ! <mf||mf>
        !!PB 15.9.2004          a = dfloat(mf) * cos(theta)/2.d0
        coup(l,mf,mf,1,1) = -dfloat(mf)/2.d0
        coup(l,mf,mf,2,2)  = dfloat(mf)/2.d0
        coup(l,mf,mf,1,2) = 0.d0
        coup(l,mf,mf,2,1) = 0.d0
        do is1=1,2
           do is2=1,2
              sum=0.d0
              do i=1,2
                 do j=1,2
                    sum=sum+dconjg(T(i,is1))*T(j,is2)*coup(l,mf,mf,i,j)
                 end do
              end do
              couplo(indj,l,mf,mf,is1,is2)=sum
           end do
        end do
        ! <mf||mf-1>
        if (mf.gt.-l) then
           coup(l,mf,mf-1,1,1) = 0.d0
           coup(l,mf,mf-1,2,2) = 0.d0
           coup(l,mf,mf-1,2,1) = 0.d0
           coup(l,mf,mf-1,1,2) = sqrt(dfloat((l+mf)*(l-mf+1)))/2.d0
           do is1=1,2
              do is2=1,2
                 sum=0.d0
                 do i=1,2
                    do j=1,2
                       sum=sum+dconjg(T(i,is1))*T(j,is2)*coup(l,mf,mf-1,i,j)
                    end do
                 end do
                 couplo(indj,l,mf,mf-1,is1,is2)=sum
              end do
           end do
        endif
        ! <mf||mf+1>
        if (mf.lt.l) then 
           coup(l,mf,mf+1,1,1) = 0.d0
           coup(l,mf,mf+1,2,2) = 0.d0 
           coup(l,mf,mf+1,2,1) = sqrt(dfloat((l-mf)*(l+mf+1)))/2.d0
           coup(l,mf,mf+1,1,2) = 0.d0
           do is1=1,2
              do is2=1,2
                 sum=0.d0
                 do i=1,2
                    do j=1,2
                       sum=sum+dconjg(T(i,is1))*T(j,is2)*coup(l,mf,mf+1,i,j)
                    end do
                 end do
                 couplo(indj,l,mf,mf+1,is1,is2)=sum
              end do
           end do
        endif
        ! <mf-1||mf>
        if (mf.gt.-l) then
           coup(l,mf-1,mf,1,1) = 0.d0
           coup(l,mf-1,mf,2,2) = 0.d0 
           coup(l,mf-1,mf,2,1) = sqrt(dfloat((l+mf)*(l-mf+1)))/2.d0
           coup(l,mf-1,mf,1,2) = 0.d0
           do is1=1,2
              do is2=1,2
                 sum=0.d0
                 do i=1,2
                    do j=1,2
                       sum=sum+dconjg(T(i,is1))*T(j,is2)*coup(l,mf-1,mf,i,j)
                    end do
                 end do
                 couplo(indj,l,mf-1,mf,is1,is2)=sum
              end do
           end do
        endif
        ! <mf+1||mf>
        if (mf.lt.l) then
           coup(l,mf,mf+1,1,1) = 0.d0
           coup(l,mf,mf+1,2,2) = 0.d0
           coup(l,mf,mf+1,2,1) = 0.d0
           coup(l,mf,mf+1,1,2) = sqrt(dfloat((l-mf)*(l+mf+1)))/2.d0
           do is1=1,2
              do is2=1,2
                 sum=0.d0
                 do i=1,2
                    do j=1,2
                       sum=sum+dconjg(T(i,is1))*T(j,is2)*coup(l,mf+1,mf,i,j)
                    end do
                 end do
                 couplo(indj,l,mf+1,mf,is1,is2)=sum
              end do
           end do
        endif
     enddo
  enddo
  
  RETURN
END subroutine couple



