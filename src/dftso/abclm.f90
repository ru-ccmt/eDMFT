SUBROUTINE ABCLM(meigve,SS,NE,NV,KV,P,DP,PE,DPE,isi,JA,indj,lfirst) 
  USE param, ONLY: nmat, nume, labc, nato, hblock, lmax, labc2, nloat
  USE lolog, ONLY: nlo, nlov, nlon, lapw
  USE rlolog, ONLY: nnrlo, nrlo
  USE abcd, ONLY: abcdlm
  USE structure, ONLY: VOL, Rmt, Rotloc, pos, BR1, rotij, tauij
  IMPLICIT NONE
  INTEGER, intent(in) :: NE, NV, isi, ja, indj, lfirst, KV(3,NMAT,2)
  REAL*8, intent(in)  :: SS(3)
  REAL*8, intent(in)  :: P(Labc+1,NATO,2),DP(Labc+1,NATO,2),PE(Labc+1,NATO,2), DPE(Labc+1,NATO,2)
  !
  !-----X CALCULATE RADIAL ENERGY DERIVATIVE ANALYTICALLY.
  !     THIS PROGRAM CALCULATE THE WAVEFUNCTION INSIDE MT SPHERE
  !     (TRANSFERED FROM OUTSIDE MT). THE WAVEFUNCTION INSIDE MT IS
  !     EXPRESSED IN TERMS OF Alm Blm Clm FOR EACH ATOM.
  !
  !**********************************************************************
  complex*16,allocatable ::  yl(:),h_yl(:,:), h_alyl(:,:),h_blyl(:,:)
  complex*16,allocatable ::  alm(:,:),blm(:,:)
  !
  complex*16:: meigve(nmat,nume)
  REAL*8    :: al(hblock),bl(hblock)
  REAL*8    :: BK(3),BK0(3),BKROT(3),bkrloc(3)
  REAL*8    :: ARG1, ARG2, ARG3, ARGT
  complex*16:: phshel,cfac
  integer   :: i, i3, ibb, ii, index, ix, l, lda, ldb, m, num
  REAL*8    :: BKX(NMAT),BKY(NMAT),BKZ(NMAT)
  REAL*8    :: dt0, fac, pi, twopi, y
  REAL*8    :: FJ(0:LMAX,NMAT),DFJ(0:LMAX,NMAT), FCT(100)
  REAL*8,     PARAMETER :: tol1 = 1.0D-6
  COMPLEX*16, PARAMETER :: IMAG = (0.0D0,1.0D0)
  allocate     ( yl(labc2),h_yl(labc2,hblock), h_alyl(labc2,hblock),h_blyl(labc2,hblock)) 
  allocate(alm(labc2,nume),blm(labc2,nume))
  PI=ACOS(-1.0D0)                                                   
  TWOPI=2.D0*PI                                                     
  Y=1.0D0                                                           
  I=1
  
  Y=1.0
  do i=1,50
     FCT(2*i-1)=Y
     Y=Y*i
  enddo
  
  DO I=1,NV
     BKX(I)=(SS(1)+KV(1,I,isi))                                               
     BKY(I)=(SS(2)+KV(2,I,isi))                                               
     BKZ(I)=(SS(3)+KV(3,I,isi))         
  enddo

  call cputim(dt0)

  FAC=4.0D0*PI*RMT(JA)**2/SQRT(VOL)
  CALL HARMON2(NV,BKX,BKY,BKZ,LMAX,FJ,DFJ,RMT(JA))

  do ix=1,labc2                                           
     do num=1,ne+nnrlo 
        alm(ix,num)=0.d0
        blm(ix,num)=0.d0
        abcdlm(1:nloat,ix,num,isi)=0.d0
     enddo
  enddo

  !print *, 'NV=', NV, 'nlo=', nlo, 'nlon=', nlon, 'nlov=', nlov, 'hblock=', hblock
  
  DO II=1,NV-(nlo(1)+nlon(1)+nlov(1)),hblock
     i3=0

     DO I=ii,min(ii+hblock-1,NV-(nlo(1)+nlon(1)+nlov(1)))
        i3=i3+1
        BK0(1)=BKX(I)                                               
        BK0(2)=BKY(I)                                             
        BK0(3)=BKZ(I)
        !
        BKROT = matmul(ROTIJ(:,:,indj),BK0)
        BK = matmul(BR1, BKROT)
        BKRLOC = matmul(ROTLOC(:,:,ja),BK)
        !
        CALL YLM(BKRLOC,LABC,YL)
        !
        ARG1 = dot_product(BKROT,POS(:,lfirst))*twopi
        ARGT = dot_product(BK0,TAUIJ(:,indj))*twopi
        !
        PHSHEL=EXP(IMAG*(ARG1+ARGT))
        do index=1,LABC2
           h_yl(index,i3)=conjg(yl(index))*phshel
        end do
     ENDDO


     index=0
     do L=0,LABC
        i3=0
        do i=ii,min(ii+hblock-1,NV-(nlo(1)+nlon(1)+nlov(1)))
           i3=i3+1
           IF(lapw(l,ja)) THEN
              al(i3)=dfj(l,i)*pe(l+1,ja,isi)- fj(l,i)*dpe(l+1,ja,isi)
              bl(i3)= fj(l,i)*dp(l+1,ja,isi)-dfj(l,i)*  p(l+1,ja,isi)
           ELSE
              al(i3)= fj(l,i)/p(l+1,ja,isi)/rmt(ja)**2
              bl(i3)= 0.d0
           ENDIF
        end do

        do m=1,2*l+1
           index=index+1
           i3=0
           do i=ii,min(ii+hblock-1,NV-(nlo(1)+nlon(1)+nlov(1)))
              i3=i3+1
              h_alyl(index,i3)=AL(i3)*h_YL(INDEX,i3)
              h_blyl(index,i3)=BL(i3)*h_YL(INDEX,i3)
           enddo
        enddo
     enddo

     ibb=min(hblock,NV-(nlo(1)+nlon(1)+nlov(1))-ii+1)
     lda=labc2
     ldb=nmat
     call zgemm('N','N',index,ne,ibb,(1.d0,0.d0),h_alyl,lda,meigve(ii,1),ldb,(1.d0,0.d0),alm,lda)
     call zgemm('N','N',index,ne,ibb,(1.d0,0.d0),h_blyl,lda,meigve(ii,1),ldb,(1.d0,0.d0),blm,lda)
  ENDDO

  index=0
  do l=0,LABC
     cfac=fac*(imag**l)
     do m=1,2*l+1
        index=index+1
        do num=1,ne
           abcdlm(1,index,num,isi)=alm(index,num)*cfac
           abcdlm(2,index,num,isi)=blm(index,num)*cfac
        enddo
     enddo
  enddo

  if (nlo(ja).ne.0) call lomain(meigve,bkx,bky,bkz,fac,ne,lfirst,indj,JA,nv,isi)
  if (nrlo(ja).ne.0) call rlomain(bkx,bky,bkz,fac,ne,lfirst,indj,JA,nv,isi,kv)

  deallocate   (yl,h_yl,h_alyl,h_blyl)
  deallocate   (alm,blm)

  RETURN
END SUBROUTINE ABCLM

