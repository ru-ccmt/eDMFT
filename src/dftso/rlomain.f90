subroutine rlomain(bkx,bky,bkz,fac,ne,lfirst,indj,jatom,n,isi,kv)
  USE param, ONLY: nmat, labc, labc2, lomax, nume
  USE abcd, ONLY: abcdlm
  USE loabcr, ONLY: alor, blor, clor
  USE lolog, ONLY: nlo, nlon, ilo, mrf
  USE rlolog, ONLY: loorext, nrlov
  USE structure, ONLY: mult, rotij, tauij, BR1, rotloc, POS
  implicit none
  real*8, intent(in) :: bkx(nmat),bky(nmat),bkz(nmat)
  real*8, intent(in) :: fac
  integer, intent(in):: ne, lfirst, indj, jatom, n, isi
  integer, intent(inout) :: kv(3,nmat,2)
  !
  complex*16 :: yl(labc2)
  complex*16 :: phs(nume),phshel,cfac,imag,pt
  integer    :: ii,iat
  real*8     :: argt
  real*8     :: bk(3), bk0(3), bkrot(3), bkrloc(3)
  real*8     :: tol1
  real*8 :: arg1, pi, twopi!, arg2, arg3
  integer:: i,l,m,m1,jneq,num, index, irf
  data   imag/(0.0d0,1.0d0)/
  data   tol1/1.0d-6/
  !
  !.initiales a,b,c of rlo                                      
  !
  pi=acos(-1.0d0)                                                   
  twopi=2.d0*pi
  ii=0
  do iat=1,jatom-1
     do l=0,lomax
        if (loorext(l,iat)) then
           ii=ii+(2*l+1)*mult(iat)
        end if
     end do
  end do
  
  i=n-(nlo(jatom)+nlon(jatom))
  
  do l=0,lomax
     if (.not.loorext(l,jatom)) then
        i=i+(2*l+1)*mult(jatom)*ilo(l,jatom)
     else
        cfac=fac*(imag**l)
        do jneq=1,mult(jatom) 
           do m1=-l,+l                                                    
              i=i+1 
              ii=ii+1 
              BK0(1)=BKX(i)
              BK0(2)=BKY(i)
              BK0(3)=BKZ(i)
              kv(1,n+ii,isi)=kv(1,i,isi)
              kv(2,n+ii,isi)=kv(2,i,isi)
              kv(3,n+ii,isi)=kv(3,i,isi)
              !
              BKROT = matmul(ROTIJ(:,:,indj), BK0)
              BK = matmul(BR1(:,:), BKROT)
              BKRLOC = matmul(ROTLOC(:,:,jatom), BK)
              !
              CALL YLM(BKRLOC,LABC,YL)
              !
              arg1 = dot_product(BKROT, POS(:,lfirst) )*twopi
              argt = dot_product(BK0, TAUIJ(:,indj))*twopi
              !
              PHSHEL=EXP(IMAG*(ARG1+ARGT))
              num = ne + nrlov(jatom) + (m1+l+1) + (jneq-1)*(2*l+1)
              phs(num) = 0.d0
              pt = phshel*1.d0
              
              if (dabs(dreal(pt)).gt.tol1) phs(num) = phs(num) + dreal(pt)
              if (dabs(dimag(pt)).gt.tol1) phs(num) = phs(num) + dimag(pt)*imag
              
              do m = -l,l
                 index = l*(l+1)+m+1
                 abcdlm(1,index,num,isi) = abcdlm(1,index,num,isi)+alor(l,jatom,isi)*conjg(yl(index))*phs(num)*cfac
                 abcdlm(2,index,num,isi) = abcdlm(2,index,num,isi)+blor(l,jatom,isi)*conjg(yl(index))*phs(num)*cfac
                 irf=mrf(l,jatom)
                 abcdlm(irf,index,num,isi) = abcdlm(irf,index,num,isi)+clor(l,jatom,isi)*conjg(yl(index))*phs(num)*cfac
              enddo
           enddo
        enddo
     endif
  enddo
END subroutine rlomain
