subroutine make_albl(jneq,r,ltop,NV,sub,imin,imax)
  use atspdt, only : P,DP,PE, DPE
  use albl
  use matrices, only : RK
  use lolog, only : lapw
  !INCLUDE 'param.inc'
  use param
  character*3 :: sub
  integer     :: jneq, ltop, nv, imin, imax, ltop2, N,j
  real*8      :: r, rkn
  DOUBLE PRECISION, allocatable ::   DFJ(:), FJ(:)
      
  ! fix for ltop=0,1 for sphbes
  ltop2=max(ltop,2)
  allocate( DFJ(0:ltop2), FJ(0:ltop2) )
  if (sub.eq.'HAM') then
     DO N=1,NV
        RKN = RK(N)
        CALL SPHBES(ltop2,R*RKN,FJ)
        CALL DVBES1(FJ,DFJ,RKN,R,ltop+1)
        DO L= 0,ltop
           if(lapw(l,jneq)) then
              AL_r(N,L) = RKN*DFJ(L)*PE(L+1,JNEQ) - FJ(L)*DPE(L+1,JNEQ)
              BL_r(N,L) = FJ(L)*DP(L+1,JNEQ) - RKN*DFJ(L)*P(L+1,JNEQ)
           else
              AL_r(N,L) = FJ(L)/P(L+1,JNEQ)/R**2
              BL_r(N,L) = 0.0d0
           endif
        enddo
     enddo
     !j=0
     do N=imin,imax
        !j=j+1
        j=N-imin+1
        AL_c(j,:)=AL_r(N,:)
        BL_c(j,:)=BL_r(N,:)
     enddo
  else
     j=0
     DO N=1,NV
        j=j+1
        RKN = RK(N)
        CALL SPHBES(ltop2,R*RKN,FJ)
        CALL DVBES1(FJ,DFJ,RKN,R,ltop+1)
        DO L = 0, ltop
           if(lapw(l,jneq)) then
              AL(j,L,1) = RKN*DFJ(L)*PE(L+1,JNEQ) - FJ(L)*DPE(L+1,JNEQ)
              BL(j,L,1) = FJ(L)*DP(L+1,JNEQ) - RKN*DFJ(L)*P(L+1,JNEQ)
           else
              AL(j,L,1) = FJ(L)/P(L+1,JNEQ)/R**2
              BL(j,L,1) = 0.0d0
           endif
        enddo
     enddo
  endif
  deallocate(dfj,fj)
  
end subroutine make_albl
