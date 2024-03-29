subroutine kptout(ss,bname,weight,ikp,kv,jspin,nv,ne,must_compute)
  USE param, only: nmat, lomax!, nato
  USE lolog, only: nlov, nlon, nlo
  USE rlolog, only: nnrlo
  USE hmsout, only: neig, en, vect, vnorm
  USE mpi, only: Qprint, vector_para, myrank, master, mpi_SendReceive12, mpi_SendReceive3
  IMPLICIT NONE
  real*8, intent(in) :: ss(3)
  character*10, intent(in) :: bname
  real*8, intent(in) :: weight
  integer, intent(in):: ikp, jspin
  integer, intent(in):: kv(3,nmat,2)
  integer, intent(in):: nv(2),ne(2)
  logical, intent(in):: must_compute
  ! locals
  integer   :: j, ie, isi, is, nv0, nv_plus_nnrlo
  integer,   allocatable :: kt(:,:)
  complex*16, allocatable:: vt(:,:)
  !
  if (Qprint .and. must_compute) then
     !  Write output to case.outputso and case.scfso
     WRITE(6,"(3x,'K=',3f10.5,1x,a10,2x,'MATRIX SIZE=',i5,2x,'WEIGHT=',f5.2/,5x,'EIGENVALUES ARE:')") SS(1),SS(2),SS(3),BNAME,ne(1)+ne(2)+2*nnrlo,WEIGHT
     if(ikp.eq.1) then
        write(8,*)
        write(8,*) '       SPIN-ORBIT EIGENVALUES:'
        WRITE(8,"(3x,'K=',3f10.5,1x,a10,2x,'MATRIX SIZE=',i5,2x,'WEIGHT=',f5.2/,5x,'EIGENVALUES ARE:')") SS(1),SS(2),SS(3),BNAME,ne(1)+ne(2)+2*nnrlo,WEIGHT
     end if
     write(6,531) (en(ie),ie=1,neig)
     write(6,6010) 0
     write(6,6030)  
     if(ikp.eq.1) then
        write(8,530) (en(ie),ie=1,neig)
        write(8,6030)
     endif
  endif
  
  ! Writing to vector and energy file (case.vectorso and case.energyso) in parallel
  ! When vector_para=True, each processor just writes into its file
  ! When vector_para=False, all prcocessor need to send their data to master, which eventually prints
  ! the data to the file, opened only on master proces.
  do isi=1,2
     is=isi
     if (jspin.eq.1) is=1
     if (must_compute) then                    ! if must_compute=False, there is no data on such node.
        nv0 = nv(isi)-(nlo(1)+nlon(1)+nlov(1)) ! G's withouth local orbitals
        nv_plus_nnrlo = nv(isi)+nnrlo          ! all G's incluing local orbitals and p1/2 states
        allocate(kt(3,nv_plus_nnrlo))          ! rearanging, so that p1/2 orbitals follow lo's for each atom
        allocate(vt(nv_plus_nnrlo,neig) )
        call Rearange_kpts1(kt, kv(:,:,is), nv0, nv(is), nmat, nv_plus_nnrlo)
        call Rearange_kpts2(vt, vect(:,1:neig,isi), nv0, nv(isi), nmat, nv_plus_nnrlo, neig)
        ! now we write what can be written independently
        if (vector_para .or. myrank.eq.master) then
           call Print_Kpoint1(41+isi-1,51+isi-1, isi.eq.2, kt, ss, bname, weight, nv_plus_nnrlo, neig)
           do j=1,neig
              call Print_Kpoint2(41+isi-1, 51+isi-1, (isi.eq.2.and.j.eq.neig), j, en(j), vt(:,j), nv_plus_nnrlo, neig)
           enddo
        endif
     endif
     if (.not.vector_para) then
        if (.not.must_compute) then 
           nv_plus_nnrlo = nv(1)+nnrlo          ! all G's incluing local orbitals and p1/2 states
           allocate(kt(3,nv_plus_nnrlo))          ! so that send-receive will work
           allocate(vt(nv_plus_nnrlo,neig) )           
        endif
        ! if vector_para=False, all but master node did not write anythning yet.
        ! It needs to communicate its results with master node
        !call mpi_SendReceive1(ss,weight,bname,nv_plus_nnrlo,neig,kt,isi,jspin,must_compute)
        !call mpi_SendReceive2(en(:neig),vt,nv_plus_nnrlo,isi,neig,must_compute)
        call mpi_SendReceive12(ss,weight,bname, en(:neig),vt, nv_plus_nnrlo,neig,kt,isi,jspin,must_compute)
        if (.not.must_compute) then 
           deallocate(kt)
           deallocate(vt)           
        endif
     endif
     if (must_compute) then
        deallocate( kt )
        deallocate( vt )
     endif
  enddo
  ! The norms respective norms of the eigenstates are written
  !  to unit 45 and unit 46.
  ! write norms of spin down part of vectors on file 45
  ! write norms of spin up part of vectors on file 46
  if (must_compute) then
     if (vector_para .or. myrank.eq.master) then
        call Print_Kpoint3(vnorm(1:neig,1:2), neig)
     endif
  endif
  if (.not.vector_para) then
     call mpi_SendReceive3(vnorm(1:neig,:),neig,must_compute)
  endif
  RETURN
  
530 FORMAT(8(7X,5F13.7/))
531 FORMAT(8(2X,5F13.7/))
6010 FORMAT(I13,' EIGENVALUES BELOW THE ENERGY ',F10.5)
6030 FORMAT(7X,14('****'))
end subroutine kptout


subroutine Print_Kpoint1(fh_vector, fh_energy, qprint_53, kt, ss, bname, weight, nv_plus_nnrlo, neig)
  IMPLICIT NONE
  integer, intent(in):: fh_vector, fh_energy  ! fh_vector=41+isi-1, fh_energy=51+isi-1
  logical, intent(in):: qprint_53
  REAL*8, intent(in) :: ss(3)
  character*10, intent(in) :: bname
  real*8, intent(in) :: weight
  integer, intent(in):: nv_plus_nnrlo,neig
  integer, intent(in):: kt(3,nv_plus_nnrlo)
  ! locals
  integer :: i
  write(fh_energy,'(3e19.12,a10,2i6,f5.1)') ss(1),ss(2),ss(3),bname,nv_plus_nnrlo,neig,weight
  write(fh_vector) ss(1),ss(2),ss(3),bname,nv_plus_nnrlo,neig,weight
  write(fh_vector)(kt(1,i),kt(2,i),kt(3,i),i=1,nv_plus_nnrlo)
  ! write also on dummy vector file 43 to be used in LAPW2 (see comment).
  if (qprint_53) write(53,'(3e19.12,a10,2i6,f5.1)') SS(1),SS(2),SS(3),BNAME,1,neig,weight
end subroutine Print_Kpoint1

subroutine Print_Kpoint2(fh_vector, fh_energy, qprint_53, j, ene, vt, nv_plus_nnrlo, neig)
  IMPLICIT NONE
  integer, intent(in)  :: fh_vector, fh_energy
  logical, intent(in)  :: qprint_53
  integer, intent(in)  :: j, nv_plus_nnrlo, neig
  real*8, intent(in)   :: ene
  complex*16,intent(in):: vt(nv_plus_nnrlo)
  ! locals
  real*8 :: tene
  integer:: i, k
  !
  write(fh_vector) j,ene
  write(fh_energy,*) j,ene
  write(fh_vector) (vt(i),i=1,nv_plus_nnrlo)
  
  if (qprint_53) then
     tene = 0.5 + ene
     do k=1,neig
        write(53,*) k,tene
        tene=tene+0.001
     enddo
  endif
end subroutine Print_Kpoint2

subroutine Print_Kpoint3(vnorm, neig)
  IMPLICIT NONE
  integer, intent(in):: neig
  real*8, intent(in) :: vnorm(neig,2)
  ! locals
  integer :: j
  write(45,'(4e20.12)') (vnorm(j,1),j=1,neig)
  write(46,'(4e20.12)') (vnorm(j,2),j=1,neig)
end subroutine Print_Kpoint3

subroutine Rearange_kpts1(kt, kv, nv0, nv_is, nmat, nv_plus_nnrlo)
  use param, only : nato, lomax
  use structure, only: mult
  USE lolog, only: ilo
  USE rlolog, only: loorext
  IMPLICIT NONE
  integer, intent(out) :: kt(3,nv_plus_nnrlo)
  integer, intent(in)  :: kv(3,nmat)
  INTEGER, intent(in)  :: nv0, nv_is, nmat, nv_plus_nnrlo
  ! locals
  INTEGER :: n, nn, nl, jatom, num, l
  !...rearrange so that p1/2 orbitals are between local orbitals of this atom and k-points of the next atom
  kt(1:3,1:nv0) = kv(1:3,1:nv0)
  ! Below is the local orbital part and p1/2 orbitals
  n=nv0
  nn=nv0
  nl=0
  do jatom=1,nato
     do l=0,lomax
        num = (2*l+1)*mult(jatom)*ilo(l,jatom)  ! how many local orbitals used in lapw1
        ! first local orbitals for this atom
        kt(1:3,n+1:n+num) = kv(1:3,nn+1:nn+num)
        n = n+num
        nn = nn+num
        if (loorext(l,jatom)) then
           ! next we add p1/2 for this atom, if exist
           num = (2*l+1)*mult(jatom) ! how many p1/2 orbitals?
           kt(1:3,n+1:n+num)=kv(1:3,nv_is+nl+1:nv_is+nl+num) ! copy into kt, followd directly by lo's of this atom
           n = n+num
           nl = nl+num
        end if
     end do
  end do
end subroutine Rearange_kpts1

subroutine Rearange_kpts2(vt, vect, nv0, nv_is, nmat, nv_plus_nnrlo, neig)
  use param, only : nato, lomax
  use structure, only: mult
  USE lolog, only: ilo
  USE rlolog, only: loorext
  IMPLICIT NONE
  complex*16, intent(out) :: vt(nv_plus_nnrlo,neig)
  complex*16, intent(in)  :: vect(nmat,neig)
  INTEGER, intent(in)     :: nv0, nv_is, nmat, nv_plus_nnrlo, neig
  ! locals
  INTEGER :: n, nn, nl, jatom, num, l
  vt(1:nv0,1:neig) = vect(1:nv0,1:neig)
  n=nv0
  nn=nv0
  nl=0
  do jatom=1,nato
     do l=0,lomax
        num = (2*l+1)*mult(jatom)*ilo(l,jatom)
        ! first local orbitals for this atom
        vt(n+1:n+num,1:neig)=vect(nn+1:nn+num,1:neig)
        n = n+num
        nn = nn+num
        if (loorext(l,jatom)) then
           ! next we add p1/2 for this atom, if exist
           num = (2*l+1)*mult(jatom)
           vt(n+1:n+num,1:neig)=vect(nv_is+nl+1:nv_is+nl+num,1:neig)
           n = n+num
           nl = nl+num
        endif
     end do
  end do
end subroutine Rearange_kpts2
  
