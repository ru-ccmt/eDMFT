! Small all-channel fixture: real radial functions with exact gamma-function
! integrals, two distinct LO energies for every l, and a Y20 potential.
module param
 integer,parameter :: nrad=801,nloat=2,lmax2=2,lomax=2,ncom=2
 character(180) :: filename_V_vns='unused'
end module
module com
 integer,parameter :: nat=1
 logical :: rel=.false.
end module
module com_mpi
 integer :: myrank=0,master=0
end module
module structure
 use param
 integer :: jri(1)=[nrad],iatnr(1)=[-1]
 real(8) :: r0(1),dx(1)=[0.012d0],rotloc(3,3,1)
end module
module lo
 use param
 integer :: ilo(0:lmax2)=nloat
 logical :: loor(nloat,0:lomax)=.true.
 real(8) :: elo(0:lomax,nloat)
end module
module atspdt
 use param
 real(8) :: el(0:lmax2)
end module
module w_atpar
 use param
 real(8) :: w_rrad1(nrad,0:lmax2,1),w_rrad2(nrad,0:lmax2,1), &
 w_rade1(nrad,0:lmax2,1),w_rade2(nrad,0:lmax2,1), &
 w_radp1(nrad,0:lmax2,1),w_radp2(nrad,0:lmax2,1), &
 w_radpe1(nrad,0:lmax2,1),w_radpe2(nrad,0:lmax2,1), &
 w_a1lo(nrad,nloat,0:lomax,1),w_b1lo(nrad,nloat,0:lomax,1), &
 w_ap1lo(nrad,nloat,0:lomax,1),w_bp1lo(nrad,nloat,0:lomax,1), &
 ri_mat(2+nloat,2+nloat,0:lmax2,1)
 integer :: w_lm(2,ncom,1),w_lmmax(1)=[2],Nhmx(1)=[0]
 integer :: lpv(1,1),lv(1,1),mpv(1,1),mv(1,1)
 complex(8) :: Vts(2+nloat,2+nloat,1,1)
contains
 subroutine retrieve_w_atpar(jatom,lfirst,lmmax)
 integer,intent(in) :: jatom
 integer,intent(out) :: lfirst,lmmax
 lfirst=1;lmmax=2
 end subroutine
end module
module readPotential
 use structure,only: r0,dx
contains
 subroutine init_V_vns(fname,unit)
 character(*),intent(in) :: fname
 integer,intent(in) :: unit
 end subroutine
 subroutine close_V_vns()
 end subroutine
 subroutine read_V_vns(v,jatom,nrad,nr,nlm)
 integer,intent(in) :: jatom,nrad,nr,nlm
 real(8),intent(out) :: v(nrad,0:nlm-1)
 integer :: i
 v=0
 do i=1,nr
 v(i,1)=0.31d0*(r0(1)*exp((i-1)*dx(1)))**2
 enddo
 end subroutine
end module
