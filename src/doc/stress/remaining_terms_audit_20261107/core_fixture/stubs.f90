module TestCases
 logical :: SimpInt=.false., ConOpt=.true.
end module
module struct
 integer :: nat=1,ndif=2,iatnr(1)=-1,mult(1)=1,jri(1)=881
 real(8) :: rotij(3,3,2),rotloc(3,3,1)
 character(4) :: lattic='P'
end module
module core_nsp
 integer,allocatable :: lm11_st(:,:,:),lmmax_st(:)
 real(8),allocatable :: vns_st(:,:,:)
 real(8) :: sum_tot(9,1)
end module
module core_sp
 real(8) :: dv1p(1022)
end module
module norm_kub
 real(8) :: c_kub(0:10,0:10)=0.d0
end module
