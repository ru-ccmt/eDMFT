subroutine hsocalc(h,ri_mat,ne,nnrlo,jspin,indj,ity)
  USE abcd, ONLY: abcdlm
  USE couples, ONLY: couplo
  USE param, ONLY: nume2, nloat, nato, labc
  USE lolog, ONLY: mrf
  IMPLICIT NONE
  complex*16, intent(out) :: h(nume2,nume2)
  real*8, intent(in)  :: ri_mat(nloat,nloat,0:labc,nato,2,2)   
  INTEGER, intent(in) :: ne(2), nnrlo, jspin, indj, ity
  ! locals
  complex*16  :: hso
  complex*16  :: dd
  logical     :: t_check
  integer     :: ibf, ibi, ief, iei, ilmf, ilmi, index, isd, isf, isi, isu, l, l_mat1, l_mat2, mf, mi, n_scr, nban2, nmatr
  complex*16,allocatable ::  t(:,:)   
  !
  nban2=ne(1)+ne(2)+nnrlo*2
  n_scr=ne(1)+ne(2)
  nmatr=nban2*(nban2+1)/2
  allocate(t((labc+1)**2*(nloat),2))   
  t=(0.0d0,0.0d0)

  do ibi = 1, nban2
     t_check=.FALSE.
     do ibf = 1, ibi
        if (ibi.le.n_scr) then
           if (ibi.gt.ne(1)) then
              iei=ibi-ne(1)
              isi=2
           else
              iei=ibi
              isi=1
           end if
        else
           if ((ibi-n_scr).gt.nnrlo) then
              iei=ibi-ne(1)-nnrlo
              isi=2
           else
              iei=ibi-ne(2)
              isi=1
           end if
        end if
        
        if (jspin.eq.1) then
           isd=1
        else
           isd=isi
        endif
        if (t_check) then
        else
           do isf=1,2
              if (jspin.eq.1) then
                 isu=1
              else
                 isu=isf
              endif
              index=0
              DO L = 1,LABC
                 DO MF = -L,L
                    ILMF = L*L +L +MF +1
                    DO l_mat2=1,mrf(l,ity)
                       index=index+1
                       t(index,isf)=0.d0
                       do mi=MF-1,MF+1
                          ILMI = L*L +L +MI +1
                          IF(ABS(mi).le.L) THEN
                             dd=0.d0
                             do l_mat1=1,mrf(l,ity)
                                dd=dd+abcdlm(l_mat1,ilmi,iei,isd)*ri_mat(l_mat1,l_mat2,l,ity,isd,isu)
                             end do
                             t(index,isf)=t(index,isf)+dd*couplo(indj,l,mf,mi,isf,isi)
                          END IF
                       end do
                    ENDDO
                 ENDDO
              ENDDO
           enddo
           t_check=.TRUE.
        endif
        
        if (ibf.le.n_scr) then
           if (ibf.gt.ne(1)) then
              ief=ibf-ne(1)
              isf=2
           else
              ief=ibf
              isf=1
           end if
        else
           if ((ibf-n_scr).gt.nnrlo) then
              ief=ibf-ne(1)-nnrlo
              isf=2
           else
              ief=ibf-ne(2)
              isf=1
           end if
        end if
        
        if (jspin.eq.1) then
           isd=1
           isu=1
        else
           isd=isi
           isu=isf
        endif
        
        index=0
        hso=0.d0
        DO L = 1,LABC
           DO MF = -L,L
              ILMF = L*L +L +MF +1
              DO l_mat2=1,mrf(l,ity)
                 index=index+1
                 hso=hso+dconjg(abcdlm(l_mat2,ilmf,ief,isu))*t(index,isf)
              ENDDO
           ENDDO
        ENDDO

        if (ibi.eq.ibf) then
           h(ibf,ibi)=h(ibf,ibi)+dconjg(hso)
        else
           h(ibf,ibi)=h(ibf,ibi)+hso
           h(ibi,ibf)=h(ibi,ibf)+dconjg(hso)
        endif
     enddo
  enddo
  deallocate (t) 
END subroutine hsocalc
