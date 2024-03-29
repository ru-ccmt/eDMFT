subroutine horbcalc(h,ri_orb,ne,jspin,indj,ity)
  USE abcd,  ONLY: abcdlm
  USE orb,   ONLY: natorb, iat, iiat, ll, vv
  USE lolog, ONLY: mrf
  USE rlolog,ONLY: nnrlo
  USE param,ONLY: nume2, labc, nato, ndif, nloat
  IMPLICIT NONE
  COMPLEX*16, intent(inout) :: h(nume2,nume2)
  REAL*8,     intent(in)    :: ri_orb(4,4,0:labc,nato,2,2)
  INTEGER,    intent(in)    :: Ne(2), jspin, indj, ity
  ! locals
  complex*16, allocatable ::  t(:,:)   
  complex*16 :: czero,dd
  INTEGER    :: nban2,n_scr,nmatr, ibf, ibi, ief, iei, ii, ilmf, ilmi, index, isd, isf, isi, iso, isu, l, l_mat1, l_mat2, mf, mi
  logical    :: t_check,calcyes
  complex*16 :: hso
  
  czero=(0.0d0,0.0d0)
  nban2=ne(1)+ne(2)+nnrlo*2
  n_scr=ne(1)+ne(2)
  nmatr=nban2*(nban2+1)/2
  allocate(t(ndif*(labc+1)**2*(nloat),2))
  t=czero
  
  calcyes=.FALSE.
  DO index=1,natorb
     if (iiat(index).eq.indj.and.iat(index).eq.ity) then
        calcyes=.TRUE.
        l=ll(index)
        goto 225
     endif
  enddo 

225 continue
  
  if (calcyes) then 
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
              DO isf=1,2
                 if (isf.eq.isi) then
                    iso=isi
                 else
                    iso=3
                 end if
                 if (jspin.eq.1) then
                    isu=1
                 else
                    isu=isf
                 endif
                 ii=0
                 DO MF = -L,L
                    ILMF = L*L +L +MF +1
                    DO l_mat2=1,mrf(l,ity)
                       ii=ii+1
                       t(ii,isf)=czero
                       do mi= -l,l
                          ILMI = L*L +L +MI +1
                          dd=czero
                          do l_mat1=1,mrf(l,ity)
                             dd=dd+abcdlm(l_mat1,ilmi,iei,isd)*ri_orb(l_mat1,l_mat2,l,ity,isd,isu)
                          enddo
                          t(ii,isf)=t(ii,isf)+dd*vv(index,mf,mi,iso)
                       enddo
                    ENDDO
                 ENDDO
              ENDDO
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
        
           ii=0
           hso=czero
           DO MF = -L,L
              ILMF = L*L +L +MF +1
              DO l_mat2=1,mrf(l,ity)
                 ii=ii+1
                 hso=hso+dconjg(abcdlm(l_mat2,ilmf,ief,isu))*t(ii,isf)
              enddo
           enddo
        
           if (ibf.eq.ibi) then
              h(ibi,ibf)=h(ibf,ibi)+dconjg(hso)
           else
              h(ibf,ibi)=h(ibf,ibi)+hso
              h(ibi,ibf)=h(ibi,ibf)+dconjg(hso)
           endif
        enddo
     enddo
  endif
  deallocate (t) 
END subroutine horbcalc
