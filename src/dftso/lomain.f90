subroutine lomain(a,bkx,bky,bkz,fac,ne,lfirst,indj,jatom,n,isi)
  USE param, ONLY: nmat,nume, labc, labc2, lomax
  USE abcd, ONLY: abcdlm
  USE loabc, ONLY: alo, blo, clo
  USE lolog, ONLY: nlon, lapw, nlo, ilo
  USE structure, ONLY: rotij, rotloc, POS, tauij, BR1, mult
  IMPLICIT NONE
  REAL*8, intent(in) :: BKX(NMAT),BKY(NMAT),BKZ(NMAT)
  REAL*8, intent(in) :: FAC
  INTEGER, intent(in):: NE, lfirst, indj, jatom, n, isi
  COMPLEX*16, intent(in):: A(NMAT,NUME)
  !
  ! intrinsic
  REAL*8 :: ddot
  COMPLEX*16:: PHSHEL,CFAC,IMAG,PT
  COMPLEX*16:: YL(LABC2)!,ALM,BLM,CLM
  COMPLEX*16:: PHS(NUME)
  REAL*8 :: arg1, arg2, arg3, argt, tol1, pi, twopi
  REAL*8 :: BK(3),BK0(3),BKROT(3),BKRLOC(3)
  INTEGER:: i, index, jlo, irf, jneq, l, m, m1, num
  DATA           IMAG/(0.0D0,1.0D0)/
  DATA           TOL1/1.0D-6/
  !------------------------------------------------------------------     
  !                                                                       
  !.initiales a,b,c of lo                                      
  !	
  PI=ACOS(-1.0D0)                                                   
  TWOPI=2.D0*PI
  i=n-(nlo(jatom)+nlon(jatom)) 
  DO L=0,LOMAX
     DO jlo=1,ilo(l,jatom)
        CFAC=FAC*(IMAG**L)
        do jneq=1,mult(jatom)
           DO M1=-l,+l                                                    
              i=i+1                                   
              BK0(1)=BKX(i)                                                      
              BK0(2)=BKY(i)                                                      
              BK0(3)=BKZ(i)
              !
              BKROT = matmul(ROTIJ(:,:,indj), BK0)
              BK = matmul(BR1, BKROT)
              BKRLOC = matmul(ROTLOC(:,:,jatom),BK)
              !
              CALL YLM(BKRLOC,LABC,YL)
              !
              ARG1 = dot_product(BKROT,POS(:,lfirst))*twopi
              ARGT = dot_product(BK0,TAUIJ(:,indj))*twopi
              !
              PHSHEL=EXP(IMAG*(ARG1+ARGT))
              DO NUM=1,NE
                 PHS(NUM)=0.d0
                 PT=PHSHEL*A(I,NUM)
                 if(dabs(dreal(pt)).gt.tol1)phs(num)=phs(num)+dreal(pt)
                 if(dabs(dimag(pt)).gt.tol1)phs(num)=phs(num)+dimag(pt)*imag
              ENDDO
              do m=-l,+l                                                    
                 index=l*(l+1)+m+1 
                 do num=1,ne
                    abcdlm(1,index,num,isi)=abcdlm(1,index,num,isi)+alo(l,jlo,jatom,isi)*conjg(yl(index))*phs(num)*cfac
                    abcdlm(2,index,num,isi)=abcdlm(2,index,num,isi)+blo(l,jlo,jatom,isi)*conjg(yl(index))*phs(num)*cfac
                    if (lapw(l,jatom)) then
                       irf=jlo+2
                    else
                       irf=jlo+1
                    endif
                    abcdlm(irf,index,num,isi)=abcdlm(irf,INDEX,num,isi)+clo(l,jlo,jatom,isi)*conjg(yl(index))*phs(num)*cfac
                 enddo
              enddo
           enddo
        enddo
     enddo
  enddo
END subroutine lomain
