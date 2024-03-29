SUBROUTINE LDA_fermi(EF,weigh,nbmax,nkpt,nume,nspin,nat)
  !     NELEC IS THE NUMBER OF ELECTRONS IN THIS COMPOUND                
  !     EF IS THE FERMI-ENERGY TO BE CALCULATED                          
  !USE param,   only: nkpt, nume, nmat
  USE com_mpi, only: myrank, master, VECFN
  USE com,     only: ELECN!, ef,elecn,xwt,nspin,nat,nband,rel,nk
  IMPLICIT NONE
  !
  REAL*8, intent(out) :: EF
  REAL*8,  intent(out):: weigh(nkpt*nspin,nume)
  INTEGER, intent(out):: nbmax
  INTEGER, intent(in) :: nkpt, nume, nat, nspin
  !
  REAL*8, allocatable  :: Eb(:,:,:)
  LOGICAL :: qprint_
  !--------------------------------------------------------------------- 
  !find nkpt, nmat and nume in energy file
  qprint_ = myrank.eq.master
  
  !CALL Read_energy_dont_store(VECFN(3:4), sumw0, nkpt, nmat, nume, nspin, nat)
  allocate( Eb(nume,nkpt,nspin) )
  CALL Read_energy_file(VECFN(3:4), Eb, nkpt, nume, nspin, nat)
  
  !ALLOCATE(weigh(nkpt*nspin,nume))
  weigh(:,:)=0.0d0
  call fermi_tetra(EF,weigh,nbmax,ELECN,Eb,nkpt,nume,nspin,qprint_)
  
  WRITE(6,*) 'Finally EF=', EF
  
  !open(999,file='checkwght.dat',status='unknown')
  !wsum=0.d0
  !print*, 'nkpt=', nkpt, 'nume=', nume
  !do ispin=1,nspin
  !   do ik=1,nkpt
  !      k = ik + nkpt*(ispin-1)
  !      do ib=1,nbmax
  !         wsum = wsum + weigh(k,ib)
  !         WRITE(999,'(I3,1x,I4,1x,g25.12,1x,g25.12,1x,g25.12)') ik,ib,weigh(k,ib)/2.d0, wsum/2.d0, Eb(ib,ik,ispin)
  !      enddo
  !   enddo
  !enddo
  !close(999)
  
  deallocate( Eb )
  !DEALLOCATE( weigh )
  !STOP
  !return	
  !
5001 FORMAT(3e19.12,a10,2i6,F5.1)
END SUBROUTINE LDA_fermi

