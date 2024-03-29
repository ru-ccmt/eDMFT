subroutine Find_nlos(lapw, loor, ilo, nlo, nlov, nlon, nrlo, nrlov, nnrlo, loorext, jspin, nato, lmax, lomax, nloat, e, elo)
  ! nlo -- how many for this atom
  ! nlov , nrlov -- first for this atom
  ! nlon -- how many after this atom
  use structure, ONLY: mult
  IMPLICIT NONE
  LOGICAL, intent(out) :: loor(0:lomax,nato), lapw(0:lmax,nato)
  INTEGER, intent(out) :: ilo(0:lomax,nato), nlov(nato), nrlov(nato), nlo(nato), nrlo(nato), nlon(nato), nnrlo
  LOGICAL, intent(in)  :: loorext(0:lomax,nato)
  INTEGER, intent(in)  :: jspin, nato, lmax, lomax, nloat
  REAL*8, intent(inout):: e(0:lmax,nato,2), elo(0:lomax,nloat,nato,2)
  ! locals
  REAL*8  :: ecut
  INTEGER :: isi, jatom, k, l
  INTEGER :: nnlo
  
  do isi=1,jspin
     nnlo=0
     nnrlo=0
     do jatom=1,nato
        !     Counting of rlo's included.
        !        nnlo    total number of local orbitals
        !        nlov    index-1 of the first local orbital belonging
        !                to atom type jatom
        !        nnrlo   is nnlo for rlo's
        !        nrlov   is nlov for rlo's
        DO l=0,lmax
           lapw(l,jatom)=.true.
           IF(e(l,jatom,isi).gt.150) THEN
              !print *, 'e-removing 200 : ', jatom, l, e(l,jatom,isi)
              e(l,jatom,isi)=e(l,jatom,isi)-200.0D0
              lapw(l,jatom)=.False.
           ENDIF
        ENDDO
        DO l=0,lomax
           ilo(l,jatom)=0
        ENDDO
        nlov(jatom)=nnlo
        nrlov(jatom) = nnrlo
  
        nlo(jatom)  = 0
        nrlo(jatom) = 0
        do l = 0,lomax
           loor(l,jatom)=.false.
           DO k=1,nloat-1
              if ((nloat-1).gt.3) then
                 ecut=99995.D+0
              else
                 ecut=995.D+0
              endif
              if (elo(l,k,jatom,isi).lt.ecut) then
                 ilo(l,jatom) = k
                 IF((lapw(l,jatom).AND.ilo(l,jatom).EQ.1).OR.ilo(l,jatom).eq.2) loor(l,jatom)=.true.
                 nnlo        =  nnlo       +((2*l+1))*mult(jatom)
                 nlo(jatom)    =  nlo(jatom)   +((2*l+1))*mult(jatom)
              endif
           ENDDO
           if (loorext(l, jatom)) then
              nnrlo          = nnrlo    +((2*l+1))*mult(jatom)
              nrlo(jatom)      = nrlo(jatom)+((2*l+1))*mult(jatom)
           endif
        enddo
     enddo
     !  Determines the position of the first k-vector in klist associated
     !  with the local orbitals belonging to atom type jatom.
     do jatom=1,nato
        nlon(jatom)=nnlo-nlo(jatom)-nlov(jatom)
     enddo
  enddo
end subroutine Find_nlos
