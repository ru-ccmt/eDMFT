SUBROUTINE CmpXqtl(xqtl2, DMFTrans, nbands, nind, maxdim2)
  ! Creates xqtl from transformation DMFTrans
  IMPLICIT NONE
  COMPLEX*16, intent(inout):: xqtl2(nbands,maxdim2,maxdim2)
  COMPLEX*16, intent(in)   :: DMFTrans(nbands,nbands,maxdim2,maxdim2)
  INTEGER, intent(in)      :: nbands, nind, maxdim2
  ! local variables
  INTEGER    :: ind1, ind2, i
  do ind1=1,nind
     do ind2=1,nind
        do i=1,nbands ! over bands
           xqtl2(i,ind1,ind2) = xqtl2(i,ind1,ind2) + DMFTrans(i,i,ind1,ind2)
        enddo
     enddo
  enddo
END SUBROUTINE CmpXqtl
