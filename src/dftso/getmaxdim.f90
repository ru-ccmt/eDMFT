subroutine getMaxDim(nmat,nume,nkp,jspin,elin,elo, Qcomplex)
  ! Gives maximum number of bands, plane waves, nspins, and linearization energies
  use MPI, only: filename_energy, MPI_Bcast_dims, myrank, master, vector_para, MPI_Bcast_dims_para, stop_MPI
  use param, only: nato, nloat, lomax, lmax
  IMPLICIT NONE
  INTEGER, intent(out) :: nmat, nume, nkp, jspin
  REAL*8, intent(out)  :: elo(0:lomax,nloat,nato,2), elin(0:lmax,nato,2)
  LOGICAL, intent(in)  :: Qcomplex
  ! locals
  LOGICAL      :: energy_file_exists
  real*8       :: s, t, z, E1, wgh
  character*10 :: bname
  integer      :: i, jatom, nn, l, ii, ios, n, nen, num, wkx, wky, wkz, isi
  real*8       :: wAsr
  complex*16   :: wAsc
  
  ! Read precise linearization energies from the vector file
  ! All nodes have to read the header of vector file, so that they can later read up to their needed k-point
  do jatom=1,nato
     read(9)(elin(l,jatom,1),l=0,lmax)            ! linearization energy from vector file
     read(9)((elo(l,nn,jatom,1),l=0,lomax),nn=1,nloat-1)
  enddo

  jspin=1
  do jatom=1,nato
     read(10,IOSTAT=ios)(elin(l,jatom,2),l=0,lmax)            ! linearization energy from vector file
     read(10,IOSTAT=ios)((elo(l,nn,jatom,2),l=0,lomax),nn=1,nloat-1)
  enddo

  if (ios.eq.0) then
     jspin=2
  else
     if (myrank.eq.master) then
        close(10, status='delete')
        close(55, status='delete')
        close(23, status='delete')
     else
        close(10)
        close(55)
        close(23)
     endif
  endif

  if (vector_para .or. myrank.eq.master) then
     ! If every node has his own file, all nodes read. If only one file exists, only master reads.
     
     INQUIRE(FILE=filename_energy(1), EXIST=energy_file_exists) ! Due to slownes of the filesystem, the energy file might not exist
     !energy_file_exists = .False.
     !print *, 'Searcing for ', TRIM(filename_energy(1)), ' returns ', energy_file_exists
     
     nkp=0
     nmat=0
     nume=0
     ! Reading energy file to find nkp,nmat,nume
     if (energy_file_exists) then  ! Due to slownes of the filesystem, the energy file might not exist
        do isi=1,jspin
           DO i=1,nato
              READ(54+isi-1,*)  ! linearization energies, which are less precise
              READ(54+isi-1,*)  ! lo linearization energies
           ENDDO
           DO
              READ(54+isi-1,'(3e19.12,a10,2i6)',IOSTAT=ios) s,t,z,bname,n,nen
              IF (ios /= 0) EXIT
              if (isi.eq.1) nkp = nkp+1
              nmat=MAX(n,nmat)
              nume=MAX(nen,nume)
              DO ii=1,nen
                 READ(54+isi-1,*) num, e1
              ENDDO
           ENDDO
           rewind(54+isi-1)
        enddo
     else
        ! Reading vector file to find nkp,nmat,nume
        do isi=1,jspin
           DO 
              READ(9+isi-1,IOSTAT=ios) s,t,z,bname,n,nen,wgh
              IF (ios /= 0) EXIT
              if (isi.eq.1) nkp=nkp+1
              nmat=MAX(n,nmat)
              nume=MAX(nen,nume)
              READ(9+isi-1) (wkx,wky,wkz,i=1,N)
              DO ii=1,nen 
                 READ(9+isi-1) num,e1
                 if (Qcomplex) then
                    READ(9+isi-1) (wAsc,i=1,N)
                 else
                    READ(9+isi-1) (wAsr,i=1,N)
                 endif
              ENDDO
           ENDDO
           rewind(9+isi-1)
           ! Now put the reading to the right place
           DO i=1,nato
              READ(9+isi-1) ! linearization energy
              READ(9+isi-1) ! linearization energy
           END DO
        end do
     end if
  endif


  !print *, myrank, 'nkp=', nkp, 'nmat=', nmat, 'nume=', nume
  ! If only master reads, we need to broadcast information at the end
  if (vector_para) then
     call MPI_Bcast_dims_para(nmat,nume,nkp)
  else
     call MPI_Bcast_dims(nmat,nume,nkp)
  endif
  !print *, myrank, 'nkp=', nkp, 'nmat=', nmat, 'nume=', nume
END subroutine getMaxDim
