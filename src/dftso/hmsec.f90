SUBROUTINE HMSEC(FL,emm,ne,nv,ee,ri_mat,ri_orb,jspin,ipr,ss,kv,p,dp,pe,dpe)
  USE param,  ONLY: nume, labc, nato, nmat, hblock, num2, nume2
  USE lolog,  ONLY: mrf, lso
  USE rlolog, ONLY: nnrlo, loorext
  USE orb,    ONLY: natorb
  USE structure, ONLY: mult
  USE hexpt,  ONLY: hexl
  USE abcd,   ONLY: abcdlm
  USE radovlp,ONLY: ruu
  USE hmsout, ONLY: vect, vnorm, en, neig, allocate_hmsout
  USE couples,ONLY: couplo
  USE vns,    ONLY: lvns
  USE store_vec, ONLY: vec
  USE meigve_mod, ONLY: meigve
  USE mpi, ONLY: Qprint
  IMPLICIT NONE
  logical, intent(in) :: FL
  REAL*8,  intent(in) :: emm(2), ee(nume,2), SS(3)
  INTEGER, intent(in) :: Ne(2), Nv(2), jspin, ipr
  INTEGER, intent(in) :: KV(3,NMAT,2)
  REAL*8,  intent(in) :: P(Labc+1,NATO,2),DP(Labc+1,NATO,2),PE(Labc+1,NATO,2),DPE(Labc+1,NATO,2)
  !
  !     THIS IS THE MAIN PROGRAM FOR SPIN-ORBIT INTERACTION. A modified
  !     version of the original B.N. Harmon program.
  !
  !     The input is the wavefunction outside Mt which is read in from file
  !     #10. Other inputs include file5(main input) and file11(potential
  !     which is *ONLY* spherical part of self-consistent flapw potential) and
  !     file20(the master structure file for flapw package)
  !
  !     RADIAL PART OF MATRIX ELEMENT OF Hs.o. IS CALCULATED IN GARADME.
  !     ANGULAR PART IS DONE BY COUPPLE. WAVEFUNCTION INSIDE MT (Alm and Blm)
  !     IS CALCULATED BY ABLM.   ............Jun Ye 1991
  !**********************************************************************
  !TESTING parallel
  !Arrays and descriptors For testing of h_ and s_ in parallel mode
  !      complex*16,allocatable :: s_test(:,:)
  !      integer :: DESC_s_test(9)
  !      complex*16 :: testsum
  !TESTING parallel
  !----------------------------------------------
  COMPLEX*16,allocatable :: hso(:)
  COMPLEX*16,allocatable :: h_(:,:),s_(:,:)
  REAL*8,    allocatable :: rwork(:)
  INTEGER,   allocatable :: iwork(:), ifail(:)
  COMPLEX*16,allocatable :: dwork(:)
  COMPLEX*16,allocatable :: ovrful(:,:)
  COMPLEX*16,allocatable :: s_vec(:,:)
  REAL*8 :: ri_mat(4,4,0:labc,nato,2,2)
  REAL*8 :: ri_orb(4,4,0:labc,nato,2,2)
  CHARACTER*1:: job2
  COMPLEX*16 :: ovrlc,hsoc,vnstot
  INTEGER    :: info, i, ibf, ibi, ief, iei, il, ilm, ilmi, ina, index, indj, irf1, irf2, is
  INTEGER    :: isd, isf, isi, ism, isu, ity, iu, j, k, l, lda, ldc, lfirst, m, m1, mi, ml
  INTEGER    :: mlp, ms, n_scr, nban2, nbelw, nf, nlow, nmatr, nnum, nrf
  REAL*8     :: cp(7)
  REAL*8     :: pi, abstol, tol, dtime1, dtime2, dtime3, dtime4, dtime5, dtime6, dtime7
  REAL*8     :: htim1,htim2
  LOGICAL    :: hs_check
  COMPLEX*16, parameter :: imag = (0.0d0,1.0d0)
  DTIME4=0.0d0
  DTIME5=0.0d0
  DTIME6=0.0d0
  DTIME7=0.0d0
  !**********************************************************************
  !  If required output for diagnostic purposes is provided.
  !  Diagnostic for abcdlm only for sequential version available
  
  if (ipr.gt.2 .and. Qprint) then
     write(6,*)'  ************ HMSEC **************'
     write(6,*)' l mlp ml    -,-    -,+    +,-    +,+'
     do l=1,2 
        do mlp=-l,l
           do ml =-l,l
              write(6,529) l,mlp,ml,(couplo(1,l,mlp,ml,1,ms),ms=1,2),(couplo(1,l,mlp,ml,2,ms),ms=1,2)
           enddo
        enddo
     enddo
     !Diagnostic for ALM,BLM,CLM in this part ONLY works for non parallel case
     write(6,*)'*********** ALM,BLM,CLM **********'
     write(6,*)'TEST:',abcdlm(4,2,1,1)
     do i=1,ne(1)+nnrlo
        write(6,602) i,ee(i,1)
        l=2
        do m=-l,l
           ilm=l*l+l+m+1
           write(6,601) m,abcdlm(1,ilm,i,1),abcdlm(2,ilm,i,1),abcdlm(4,ilm,i,1)
        enddo
     enddo
  endif

529 format(3i3,4(2f7.3,1X))
602 format(' level',i3,' energy',f20.14)
601 format(i3,3(2f12.6,2x))
  
  PI=3.1415926535898D0
  TOL=1.D-6
  CALL CPUTIM(DTIME1)
  
  nban2=ne(1)+ne(2)+nnrlo*2
  n_scr=ne(1)+ne(2)
  nmatr=nban2*(nban2+1)/2
  
  allocate(h_(nume2,nume2))

  if (nnrlo.eq.0) then
     allocate      (s_(1,1))
     allocate      (ovrful(1,1))
  else
     allocate      (s_(nume2,nume2))
     allocate      (ovrful(2*nnrlo,2*nnrlo))
     s_=(0.0d0,0.0d0)
  endif
  h_=(0.0d0,0.0d0)
  
  hs_check=.TRUE.
  do ity=1,nato
     !lfirst=lfirst+mult(ity-1)
     lfirst = 1+sum(mult(1:ity-1))
     do ina=0,mult(ity)-1
        indj=lfirst+ina
        call cputim(htim1)            
        do isi=1,jspin
           call abclm(meigve(:,:,isi),ss,ne(isi),nv(isi),kv,p,dp,pe,dpe,isi,ity,indj,lfirst)
        enddo
        call cputim(htim2) 
        DTIME4=DTIME4+htim2-htim1           
        do ibi = 1, nban2
           do ibf = 1, ibi
              hsoc  = 0.d0
              ovrlc = 0.d0
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
              
              if (isf.ne.isi) goto 700
              
              if (ibi.le.n_scr.and.hs_check) then
                 
                 if (ibi.eq.ibf)  then
                    hsoc  = ee(iei,isd)
                    ovrlc = (1.d0, 0.d0)
                 end if
                 
              else if (ibi.gt.n_scr) then
                 do  l=1,1
                    if (.not.loorext(l,ity)) cycle
                    do mi=-l,l
                       ilmi = l*l +l +mi +1
                       nrf=mrf(l,ity)
                       do irf1=1,nrf
                          do irf2=1,nrf
                             ovrlc=ovrlc+abcdlm(irf1,ilmi,iei,isd)*dconjg(abcdlm(irf2,ilmi,ief,isu))*ruu(l,ity,isu,isd,irf1,irf2)
                          enddo
                       enddo
                    enddo
                 enddo
                 
                 if (ibf.le.n_scr.and.hs_check) then
                    hsoc = ee(ief,isu)*ovrlc
                 else if (ibi.gt.n_scr.and.ibf.le.n_scr.and.(.not.hs_check)) then
                    hsoc = ee(ief,isu)*ovrlc
                 else if (ibf.gt.n_scr) then
                    do l=1,1
                       if (.not.loorext(l,ity)) cycle
                       do m1=-l,l
                          ilmi = l*l +l + m1 +1
                          nrf=mrf(l,ity)
                          hsoc = hsoc                            &
                               +   abcdlm(1,ilmi,iei,isd)*dconjg(abcdlm(1,ilmi,ief,isu))*hexl(l,ity,isu,1)         &
                               +   abcdlm(1,ilmi,iei,isd)*dconjg(abcdlm(2,ilmi,ief,isu))*hexl(l,ity,isu,2)         &
                               +   abcdlm(1,ilmi,iei,isd)*dconjg(abcdlm(nrf,ilmi,ief,isu))*hexl(l,ity,isu,3)         &
                               +   abcdlm(2,ilmi,iei,isd)*dconjg(abcdlm(nrf,ilmi,ief,isu))*hexl(l,ity,isu,6)         &
                               +   abcdlm(nrf,ilmi,iei,isd)*dconjg(abcdlm(1,ilmi,ief,isu))*hexl(l,ity,isu,7)         &
                               +   abcdlm(nrf,ilmi,iei,isd)*dconjg(abcdlm(2,ilmi,ief,isu))*hexl(l,ity,isu,8)         &
                               +   abcdlm(nrf,ilmi,iei,isd)*dconjg(abcdlm(nrf,ilmi,ief,isu))*hexl(l,ity,isu,9)
                       enddo
                       if (lvns(ity)) then
                          call hns(ity,isu,iei,ief,vnstot)
                          hsoc = hsoc +vnstot
                       end if
                    enddo
                 else
                 end if
              else
              end if
              if (ibf.eq.ibi) then
                 h_(ibf,ibi)=h_(ibf,ibi)+dconjg(hsoc)
              else
                 h_(ibf,ibi)=h_(ibf,ibi)+hsoc
                 h_(ibi,ibf)=h_(ibi,ibf)+dconjg(hsoc)
              endif
              if (nnrlo.eq.0) cycle
              if (ibf.eq.ibi) then
                 s_(ibi,ibf)=s_(ibi,ibf)+dconjg(ovrlc)
              else
                 s_(ibf,ibi)=s_(ibf,ibi)+ovrlc
                 s_(ibi,ibf)=s_(ibi,ibf)+dconjg(ovrlc)
              endif
700           continue 
           enddo   !end ibf
        enddo      !end ibi
        
            
        call cputim(htim1)
        if (lso(ity)) then
           call hsocalc(h_,ri_mat,ne,nnrlo,jspin,indj,ity)
        endif
        call cputim(htim2)
        DTIME5=DTIME5+htim2-htim1
        call cputim(htim1)
        if (natorb.gt.0) call horbcalc(h_,ri_orb,ne,jspin,indj,ity)
        call cputim(htim2)
        DTIME6=DTIME6+htim2-htim1
        hs_check=.FALSE.
     enddo
  enddo
  DEALLOCATE(abcdlm)
  
  call cputim(htim1)

  CALL allocate_hmsout(nmat,nume2)
  call cputim(htim2)
  DTIME7=DTIME7+htim2-htim1

  if (nnrlo.ne.0) then
     do i=1,2*nnrlo
        do j=1,2*nnrlo
           ovrful(i,j)=s_(n_scr+i,n_scr+j)
        end do
     end do
  endif
  !rschmid
  !   The hamiltonmatrix has been set up now.
  !rschmid
  
  CALL CPUTIM(DTIME2)
  
  !rschmid
  !    Provides diagnostic ouput of the Hamiltonmatrix (only non-parallel version)
  !rschmid
  if (ipr.gt.3 .and. Qprint) then
     write(6,*)' ******* real. part of h_ *****'
     do i=1,nban2
        write(6,568)(real(h_(i,j)),j=1,i)
     enddo
     write(6,*)' ******* imag. part of h_ *****'
     do i=1,nban2
        write(6,568)(dimag(h_(i,j)),j=1,i)
     enddo
     if (nnrlo.ne.0) then
        write(6,*)' ******* real  part of s_ *****'
        do i=n_scr+1,nban2
           write(6,568)(real(s_(i,j)),j=1,i)
        enddo
        write(6,*)' ******* imag. part of s_ *****'
        do i=n_scr+1,nban2
           write(6,568)(dimag(s_(i,j)),j=1,i)
        enddo
     endif
  endif
568 format(12f10.6)
  call cputim(cp(2))
  call cputim(cp(1))
  
  
  if (nnrlo.eq.0) goto 555
  
  !**************************************************************
  !.....SOLVE THE SECULAR EQUATION HSO*VEC=EN*VEC
  !
  !....computes the subblock of Cholesky transformation of the
  !....overlap matrix
  M=2*nnrlo
  K=n_scr

  !>>>>>>>>> start CHOLESKY
  
  !      cp=0
  call cputim(cp(1))
  call cputim(cp(2))
  call cputim(cp(3))
  call cputim(cp(4))
  call cputim(cp(5))
  
  lda=nume2
  call zgemm('N','C',M,M,K,(-1.d0,0.d0),s_(n_scr+1,1),lda,s_(n_scr+1,1),lda,(1.d0,0.d0),s_(n_scr+1,n_scr+1),lda)
  if(.not.fl) then
     job2='N'
  else
     job2='V'
  endif
  il=0
  iu=0
  abstol=0.0d0
  call zpotrf('L',M,s_(n_scr+1,n_scr+1),lda,info)
  !>>>>>>>>> end CHOLESKY
  
  call cputim(cp(2))
  
  !>>>>>  H~_ba=H_ba-S_ba*H_aa
  call zgemm('N','N',M,K,K,(-1d0,0.d0),s_(n_scr+1,1),lda,h_,lda,(1.d0,0.d0),h_(n_scr+1,1),lda)
  !>>>>>
  call zgemm('C','C',M,M,K,(-1.d0,0.d0),h_(1,n_scr+1),lda,s_(n_scr+1,1),lda,(1.d0,0.d0),h_(n_scr+1,n_scr+1),lda)
  call zgemm('N','C',M,M,K,(-1.d0,0.d0),s_(n_scr+1,1),lda,h_(n_scr+1,1),lda,(1.d0,0.d0),h_(n_scr+1,n_scr+1),lda)

  call ztrtri('L','N',M,s_(n_scr+1,n_scr+1),lda,info)

  !>>>>>  H`_ba=(U^+_bb)^-1*H~_ba
  call ztrmm('L','L','N','N',M,K,(1.d0,0.d0),s_(n_scr+1,n_scr+1),lda,h_(n_scr+1,1),lda)
  !>>>>>

  !>>>>>  H`_bb=(U^+_bb)^-1*(H_bb-S_ba*H~^+_ba-H_ba*S_ba^+)*U^-1
  call ztrmm('L','L','N','N',M,M,(1.d0,0.d0),s_(n_scr+1,n_scr+1),lda,h_(n_scr+1,n_scr+1),lda)
  call ztrmm('R','L','C','N',M,M,(1.d0,0.d0),s_(n_scr+1,n_scr+1),lda,h_(n_scr+1,n_scr+1),lda)
  !>>>>>

555 continue
  if(.not.fl) then
     job2='N'
  else
     job2='V'
  endif
  il=0
  iu=0
  abstol=0.0d0
  
  allocate (hso(num2))
  index=0
  do ibi=1,nban2
     do ibf=1,ibi
        index=index+1
        hso(index)=dconjg(h_(ibi,ibf))
     end do
  end do
  deallocate (h_)
  allocate (vec(nume2,nume2))
  allocate (rwork(7*nume2),iwork(5*nume2),ifail(nume2),dwork((2+hblock)*nume2))
  call cputim(cp(3))
  
  !TESTING sequential
  ! Official lapack routine zheevx as alternative to Kvasnicka routine, zhhevx
  !call ZHEEVX( job2,'V','L',nban2, h_,nume2,emm(1), emm(2),il,iu,abstol,neig,en, vec, nume2,dwork,(2+hblock)*nume2,rwork,iwork,ifail,info)
  !TESTING sequential
  call zhhevx(job2,'V','U',nban2,hso,emm(1),emm(2),il,iu,abstol,neig,en,vec,nume2,dwork,rwork,iwork,ifail,hblock,nume2,nbelw,info)                                               
         
  deallocate   (hso)
  deallocate   (rwork,iwork,ifail,dwork)
  
  call cputim(cp(4))
  
  if (nnrlo.eq.0) goto 1555
  !  backsubstitution  
  call ztrmm('L','L','C','N',M,neig,(1.d0,0.d0),s_(n_scr+1,n_scr+1),lda,vec(n_scr+1,1),lda)
  ldc=nume2
  call zgemm('C','N',k,neig,m,(-1.d0,0.d0),s_(n_scr+1,1),lda,vec(n_scr+1,1),lda,(1.d0,0.d0),vec(1,1),lda)

1555 continue

      
!-> Getting lapwso back
  do is = 1, 2
     do i = 1, nv(is)+nnrlo
        do ibi = 1,neig
           vect(i,ibi,is) = (0.d0,0.d0)
        enddo
     enddo
  enddo
  
  do is=1,2
     if (jspin.eq.2) then
        ism=is
     else
        ism=1
     end if
     
     call zgemm('N','N',nv(is),neig,ne(is),(1.d0,0.d0),meigve(1,1,ism),nmat,vec(1+ne(1)*(is -1),1),nume2,(1.d0,0.d0),vect(1,1,is),nmat)
     
     ! This part is calculated in kptout.F for the parallel version:
     do j=1,nnrlo
        do ibi=1, neig
           nf = ne(1)+ne(2)+(is-1)*nnrlo
           vect(nv(is)+j,ibi,is) = vect(nv(is)+j,ibi,is)+(vec(nf+j,ibi))
        enddo
     enddo
  enddo
  deallocate(meigve)

!-> Getting lapwso back
  
  call cputim(cp(5))
      
  !  calculate the norm of the spin parts of vector
  !      if (ipr.ge.1) then
  if (nnrlo.ne.0) then
     do i=1,2*nnrlo
        do j=1,2*nnrlo
           s_(n_scr+i,n_scr+j)=ovrful(i,j)
        end do
     end do
  endif
  

  vnorm(j,1) = 0.d0
  vnorm(j,2) = 0.d0
  if (nnrlo.ne.0) then
     allocate(s_vec(nume2,nume2))
     !-------------s_ times vec --------------------------------
     s_vec=(0.0d0,0.0d0)
     call zhemm('L','U',nban2,neig,(1.0d0,0.0d0),s_,nume2,vec,nume2,(1.0d0,0.0d0),s_vec,nume2)
     !-------------s_ times vec --------------------------------
     !----------up spin------------------
     s_=(0.0d0,0.0d0)
     call zgemm('C','N',neig,neig,ne(1),(1.0d0,0.0d0),vec,nume2,s_vec,nume2,(0.0d0,0.0d0),s_,nume2)
     nlow=n_scr+1
     nnum=n_scr+nnrlo-nlow+1
     call zgemm('C','N',neig,neig,nnum,(1.0d0,0.0d0),vec(nlow,1),nume2,s_vec(nlow,1),nume2,(1.0d0,0.0d0),s_,nume2)
     do j=1,neig
        vnorm(j,1)=vnorm(j,1)+s_(j,j)
     enddo
     !----------up spin------------------
     !----------dn spin------------------
     s_=(0.0d0,0.0d0)
     nlow=ne(1)+1
     nnum=n_scr-nlow+1
     call zgemm('C','N',neig,neig,nnum,(1.0d0,0.0d0),vec(nlow,1),nume2,s_vec(nlow,1),nume2,(0.0d0,0.0d0),s_,nume2)
     nlow=n_scr+nnrlo+1
     nnum=nban2-nlow+1
     call zgemm('C','N',neig,neig,nnum,(1.0d0,0.0d0),vec(nlow,1),nume2,s_vec(nlow,1),nume2,(1.0d0,0.0d0),s_,nume2)
     do j=1,neig
        vnorm(j,2)=vnorm(j,2)+s_(j,j)
     enddo
     !----------dn spin------------------
     deallocate(s_vec)
  else   !nnrlo.eq.0
     allocate(h_(nume2,nume2))
     !-----------up spin-----------------------
     h_=(0.0d0,0.0d0)
     call zgemm('C','N',neig,neig,ne(1),(1.0d0,0.0d0),vec,nume2,vec,nume2,(0.0d0,0.0d0),h_,nume2)
     do j=1,neig
        vnorm(j,1)=vnorm(j,1)+h_(j,j)
     enddo
     !-----------up spin-----------------------
     !-----------dn spin-----------------------
     h_=(0.0d0,0.0d0)
     nlow=ne(1)+1
     nnum=n_scr-nlow+1
     call zgemm('C','N',neig,neig,nnum,(1.0d0,0.0d0),vec(nlow,1),nume2,vec(nlow,1),nume2,(0.0d0,0.0d0),h_,nume2)
     do j=1,neig
        vnorm(j,2)=vnorm(j,2)+h_(j,j)
     enddo
     deallocate(h_)
  endif
  deallocate  (s_,ovrful)
  
  call cputim(cp(6))

  do j = 1,5
     cp(j) = cp(j+1)-cp(j)
  enddo

  if (Qprint) then
     WRITE(6,"(A)") 'Computation time on process 0'
     WRITE(6,1001) 'Cholesky complete' , CP(1) 
     WRITE(6,1001) 'Transform to eigenwertproblem' , CP(2) 
     WRITE(6,1001) 'Compute eigenvalues' , CP(3) 
     WRITE(6,1001) 'Backtransform' , CP(4) 
     WRITE(6,1001) 'norm' , CP(5) 
  endif
  
1001 FORMAT (1X,'Seclr4(',A,') :', t50, f9.3)
  
  CALL CPUTIM(DTIME3)
  if (Qprint) then
     write(6,6020)dtime2-dtime1,dtime3-dtime2
     write(6,6021)DTIME5,DTIME6
     write(6,6023)DTIME4
     write(6,6022)DTIME7
  endif
  RETURN
6020 FORMAT(7X,'CPTIME HAMILT=',F9.3,', DIAG=',F9.3/)
6021 FORMAT(7X,'CPTIME Hsoset=',F9.3,', Horb=',F9.3/)
6023 FORMAT(7X,'CPTIME HAMILTabclm=',F9.3,/)
6022 FORMAT(7X,'CPTIME HAMILTextra=',F9.3,/)
END SUBROUTINE HMSEC
