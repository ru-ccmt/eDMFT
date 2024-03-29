!-*- f90 -*- - file contains F90 code in free format

subroutine trafo(l,pseudo,T)
  IMPLICIT NONE
  integer, intent(in) :: L
  logical, intent(in) :: pseudo
  complex*16, intent(out) :: T(2*L+1,2*L+1)
  ! trafo provides transformation matrices for eigenvectors
  ! of an octahedral potential from |l,m> (l=0,1,2,3) basis to
  ! basis of triplets, doublets and singlets (both 
  ! pseudotriplets and real functions) see appendices 2.1 and 4.2 of
  ! H. Watanabe 'Operator Methods in Ligand Field Theory'
  ! Prentice Hall, 1966.
  complex*16 :: img
  complex*16 :: TF(1:2*L+1,1:2*L+1), sum
  logical    :: unitary
  integer    :: ipr, m, m1, m2, m3, sh
  real*8     :: s2, s3, s5, s8, a, b
  ipr=-1
  img=(0.,1.)
  s2=sqrt(2.)
  T(:,:)=0
  sh=L+1
! real functions are ordered: x,y,z (T1 triplet), ksi,eta,dzeta (T2)
  if(L.eq.0)then
     T(sh+0,sh+0)=1.
  else if(L.eq.1)then
     if(pseudo)then
        do m=-L,L
           T(sh+m,sh+m)=1.
        enddo
     else
        T(sh-1,sh-1)=  1./s2    !x
        T(sh-1,sh+1)= -1./s2    !x
        T(sh+0,sh-1)=  img/s2   !y
        T(sh+0,sh+1)=  img/s2   !y
        T(sh+1,sh+0) = 1.       !z
     endif
  else if(L.eq.2)then 
! ordering: doublet z2, doublet x2-y2, triplet yz, triplet xz, triplet xy
     T(sh-2,sh+0) = 1.          !z2
     T(sh-1,sh-2) = 1./s2       !x2-y2
     T(sh-1,sh+2) = 1./s2       !x2-y2
     if(pseudo)then
        T(sh+0,sh+1)=  1.
        T(sh+1,sh-2)= -1./s2
        T(sh+1,sh+2)=  1./s2
        T(sh+2,sh-1)= -1.
     else
        T(sh+0,sh-1) =  1./s2    !ksi=yz
        T(sh+0,sh+1) = -1./s2    !ksi
        T(sh+1,sh-1) =  img/s2   !eta=xz
        T(sh+1,sh+1) =  img/s2   !eta
        T(sh+2,sh-2) = -img/s2   !zeta=xy
        T(sh+2,sh+2) =  img/s2   !zeta
     endif
  else if(L.eq.3)then
     s3=sqrt(3.)
     s5=sqrt(5.)
     s8=sqrt(8.)
! f-states ordering: singlet, triplet T1, triplet T2
! singlet
     T(sh-3,sh-2)=-1./s2
     T(sh-3,sh+2)= 1./s2
! triplet T1
     TF(sh-2,sh+3) = -s5/s8
     TF(sh-2,sh-1) = -s3/s8
     TF(sh-1,sh+0) =  1.
     TF(sh+0,sh+1) = -s3/s8
     TF(sh+0,sh-3) = -s5/s8
! triplet T2
     TF(sh+1,sh+1) = s5/s8
     TF(sh+1,sh-3) =-s3/s8
     TF(sh+2,sh-2) = 1./s2
     TF(sh+2,sh+2) = 1./s2
     TF(sh+3,sh+3) =-s3/s8
     TF(sh+3,sh-1) = s5/s8
     if(.not.pseudo)then
        a=s5/4.
        b=s3/4.
        ! real singlet A2
        T(sh-3,sh-2) =img/s2
        T(sh-3,sh+2) =-img/s2
! real triplet T1
        T(sh-2,sh-3)= a        !x
        T(sh-2,sh-1)=-b
        T(sh-2,sh+1)= b
        T(sh-2,sh+3)=-a
        T(sh-1,sh-3)=-img*a    !y
        T(sh-1,sh-1)=-img*b
        T(sh-1,sh+1)=-img*b
        T(sh-1,sh+3)=-img*a
        T(sh+0,sh+0)= 1.       !z
! real triplet T2
        T(sh+1,sh-3)=-b         !ksi
        T(sh+1,sh-1)=-a
        T(sh+1,sh+1)= a
        T(sh+1,sh+3)= b
        T(sh+2,sh-3)=-img*b     !eta
        T(sh+2,sh-1)= img*a
        T(sh+2,sh+1)= img*a
        T(sh+2,sh+3)=-img*b
        T(sh+3,sh-2)= 1./s2     !zeta
        T(sh+3,sh+2)= 1./s2
     else
! pseudo triplets
        do m1=-2,3
           do m2=-3,3
              T(sh+m1,sh+m2)=TF(sh+m1,sh+m2)
           enddo
        enddo
     endif
  endif
  if(ipr.gt.0)then
     if(pseudo)then
        write(6,100)L
100     format(' L=',i3,'. Unitary transformation to pseudotriplets')
     else
        write(6,101)L
101     format(' L=',i3,'. Unitary transformation real basis')
     endif
     write(6,*)' Real part of unitary matrix'
     do m1=-l,l
        write(6,555)(dble(T(m1+sh,m2+sh)),m2=-l,l)
     enddo
     write(6,*)' Imaginary part of unitary matrix'
     do m1=-l,l
        write(6,555)(imag(T(m1+sh,m2+sh)),m2=-l,l)
     enddo
     if(ipr.gt.2)then
        write(6,*)' Real part of herm. conjugate matrix'
        do m1=-l,l
           write(6,555)(dble(conjg(T(m2+sh,m1+sh))),m2=-l,l)
        enddo
        write(6,*)' Imaginary part of herm. conjugate matrix'
        do m1=-l,l
           write(6,555)(imag(conjg(T(m2+sh,m1+sh))),m2=-l,l)
        enddo
     endif
  endif
! test of unitarity
  unitary=.true.
  TF = matmul(T, conjg(transpose(T)))
  do m1=1,2*l+1
     do m2=1,2*l+1
        if(m1.eq.m2)then
           if(abs(TF(m1,m2)-1.).gt.0.000001)then
              write(6,666) m1-sh,m2-sh,sum
              unitary=.false.
666           format(' incorrect normalization',2i4,2f12.5)
           endif
        else
           if(abs(TF(m1,m2)).gt.0.000001)then
              write(6,667)m1-sh,m2-sh,sum
              unitary=.false.
667           format(' incorrect orthogonality',2i4,2f12.5)
           endif
        endif
     enddo
  enddo
  if(ipr.gt.1)then
     if(unitary)write(6,*)' Unitarity of transformation checked'
  endif
555 format(7f9.4)
  !print *, 'Trafor finished'
  !print *, T
  return
end subroutine trafo
