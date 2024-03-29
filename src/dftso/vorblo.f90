subroutine vorblo(index,iatom,l,isi,vorb)
  USE orb,    ONLY: iat, iiat,ll, vv
  USE structure, ONLY: mult
  USE rotmat, ONLY: det, phase
  ! implicit real*8 (a-h,o-z)
  IMPLICIT NONE
  INTEGER, intent(in)    :: iatom, l, isi
  INTEGER, intent(inout) :: index
  COMPLEX*16, intent(in) :: vorb(-3:3,-3:3)
  ! local
  COMPLEX*16 :: imag
  INTEGER    :: indj, iph, jatom, mu, mi, mf
  imag=(0.d0,1.d0)
  indj=0
  iph=0
  do jatom=1,iatom-1
     do mu=1,mult(jatom)
        indj=indj+1
     end do
  end do

  do mu=1,mult(iatom)
     index=index+1
     iat(index)=iatom
     ll(index)=l
     iiat(index)=indj+mu
     if (isi.eq.3) iph=1
     do mi=-l,l
        do mf=-l,l
           if (det(indj+mu).gt.0) then
              vv(index,mf,mi,isi)=vorb(mf,mi)*exp(imag*iph*phase(indj+mu))
           else
              vv(index,mf,mi,isi)=dconjg(vorb(-mf,-mi))*(-1)**(mf+mi)*exp(imag*iph*phase(indj+mu))
           end if
        end do
     end do
  enddo
end subroutine vorblo
	     

			

