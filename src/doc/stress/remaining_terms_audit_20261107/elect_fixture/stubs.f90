module TestCases
 logical :: SimpInt=.false.,ConOpt=.true.
end module
module struct
 integer :: nat=1,ndif=1,iatnr(1)=-1,mult(1)=1,jri(1)=881
 real(8) :: zz(1)=0,vol=1,rmt(1)=1,r0(1)=1d-5,dx(1),rotij(3,3,1),rotloc(3,3,1)
 character(4) :: irel='NREL',lattic='P'
end module
module densit
 integer :: lmmaxx=264
 real(8) :: clm(881,264,1),clmsp(881,264,1,2)
end module
module work
 real(8) :: v(881,264,1),exc(881,264,1,2),vxc(881,264,1,2)
end module
module gga_str_kb
 real(8) :: gga_fin_str(9)=0
end module
module parallel
 integer :: myid=0,nstart(0:0)=1,nstop(0:0)=1
end module
module norm_kub
 real(8) :: c_kub(0:10,0:10)=0
end module
