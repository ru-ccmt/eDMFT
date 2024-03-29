! @Copyright 2007 Kristjan Haule
module RealBubble
  !  ###############################
  !  # Computing real axis Bubble  #
  !  ###############################
  IMPLICIT NONE
  REAL*8, allocatable :: Ome(:)
  REAL*8, allocatable :: chi0r0(:,:,:)
  INTEGER :: norb, nOme
CONTAINS

  SUBROUTINE RealBubble__Init__()
    use LinLogMesh
    use Qlist
    use greenk    
    IMPLICIT NONE
    
    nOme=nom-zero_ind+1
    allocate( Ome(nOme) )
    Ome(:) = om(zero_ind:)
    norb=cixdm
    ! This is the zero frequency value
    allocate( chi0r0(nQ,norb,norb) )
  END SUBROUTINE RealBubble__Init__
  
  INTEGER FUNCTION find_index(list, element)
    INTEGER, dimension(:) :: list
    INTEGER, intent(in) :: element
    !
    INTEGER :: i
    do i=1,size(list)
       if (list(i).EQ.element) EXIT
    enddo
    find_index =i
    return
  END FUNCTION find_index
  
  SUBROUTINE CmpChiQ0_real(fileout)
    ! This routine actaually computes Bubble on real axis using formula
    ! P''(Om)_{q,ab} = \sum_k \int_0^{Om}  G''(x)_{k,ab} G''(x-Om)_{k-q,ba} /pi
    ! here G'' = (G - G^+)/(2*i)  and is a complex function for off-diagonal components
    ! We compute P'(Om) by generalized Kramar's-Kronig relation, i.e.,
    !
    !    1/2[chi(om+i*delta)+chi(om-i*delta)] =-1/pi Principal \int_{-\infty}^{\infty} [chi(x+i*delta)-chi(x-i*delta)]/(2*i)/(om-x)
    !
    use LinLogMesh
    use Qlist
    use Klist
    use greenk
    IMPLICIT NONE
    CHARACTER*200, intent(in) :: fileout
    !
    CHARACTER*10 :: skii
    CHARACTER*210:: FNAME
    INTEGER,PARAMETER :: CHUNK= 4
    INTEGER,PARAMETER :: fh_out= 99
    INTEGER :: Q(3)
    INTEGER :: iQ, Qi, ik, ikq, iOm, im, iom1, iom2, izero, dOm, level, iomfin, iorb, jorb
    REAL*8     :: PI, t1c, t1w, t2c, t2w, time_Q_c, time_Q_w, rf, dsum, dsum2
    COMPLEX*16 :: csum, irhok, irhoq, cf
    COMPLEX*16 :: ImBubO(norb,norb)
    INTEGER, allocatable:: om_idx(:)
    COMPLEX*16, allocatable :: ImBub(:,:,:), Bub(:,:,:)
    REAL*8, allocatable :: dx(:), Imf(:), Ref(:), OmKK(:)
    COMPLEX*16, allocatable :: Imc(:), Rec(:)
    COMPLEX*16, PARAMETER  :: IMAG = (0.0D0,1.0D0)
    
    PI=ACOS(-1.0D0)
    
    allocate( ImBub(norb,norb,nOme), Bub(norb,norb,nOme) )
    allocate( om_idx(nom), dx(nom) )

    ! Arrays for Kramars-Kronig
    allocate( Imf(2*nOme-1), Ref(nOme) )
    allocate( Imc(2*nOme-1), Rec(nOme) )
    allocate( OmKK(2*nOme-1) )
    ! Mesh for Kramars-Kronig
    do iOm=1,nOme
       OmKK(nOme+iOm-1) = Ome(iOm)
       OmKK(nOme-iOm+1) = -Ome(iOm)
    enddo
    
    time_Q_c=0
    time_Q_w=0
    do iQ=1,nQ
       Q = Qp(:3,iQ)
       Qi = k_index(Q)
       !print *, 'Qi=', Qi, Q
       
       call cputim(t1c)
       call walltim(t1w)

       ImBub(:,:,1) = 0
       
       do iOm=2,nOme
          level = (iOm-2)/Nd+1                      !# which linear mesh should we use?
          om_idx(:) = idx(level,:)                  !# index for the linear mesh on this level
          izero= nidx(level)/2+1                    !# zero_ind on this level

          iomfin = find_index(om_idx(:nidx(level)), zero_ind+iOm-1) ! finds index such that om(om_idx(iomfin))==Ome(iOm)
          dOm = iomfin - izero   !# om-Om in integer notation is im-dOm

          ! dx for trapezoid integration
          dx(izero)=0.5*(om(om_idx(izero+1))-om(om_idx(izero)))
          do im=izero+1,iomfin-1
             dx(im) = 0.5*(om(om_idx(im+1))-om(om_idx(im-1)))
          enddo
          dx(iomfin) = 0.5*(om(om_idx(iomfin))-om(om_idx(iomfin-1)))
          
          ImBubO=0.0
          !$OMP  PARALLEL DO SHARED(gk,nkp,norb,k_m_q,Qi,om_idx,dx)&
          !$OMP& PRIVATE(ik,iorb,jorb,im,ikq,iom1,iom2,csum,irhok,irhoq)&
          !$OMP& SCHEDULE(STATIC,CHUNK)&
          !$OMP& REDUCTION(+:ImBubO)
          do ik=1,nkp
             ikq  = k_m_q(ik,Qi)
             do iorb=1,norb
                do jorb=1,norb
                   csum=0.0
                   do im=izero,iomfin
                      iom1 = om_idx(im)
                      iom2 = om_idx(im-dOm)
                      irhok=(gk(iom1,iorb,jorb,ik) -conjg(gk(iom1,jorb,iorb,ik)))/2.0
                      irhoq=(gk(iom2,jorb,iorb,ikq)-conjg(gk(iom2,iorb,jorb,ikq)))/2.0
                      csum = csum - irhok*irhoq*dx(im) !  ImG_k * ImG_{k+q}
                   enddo
                   ImBubO(iorb,jorb) = ImBubO(iorb,jorb) + csum/(nkp*PI)
                enddo
             enddo
          enddo
          !$OMP END PARALLEL DO
          ImBub(:,:,iOm)=ImBubO(:,:)
       enddo

       do iorb=1,norb
          do jorb=1,norb
             dsum=0.
             do iOm=1,nOme
                dsum = dsum + abs(dimag(ImBub(iorb,jorb,iOm)))
             enddo
             ! Check if we need complex Kramars-Kronig
             if (dsum.lt.1e-4) then
                ! creating a function which is defined for both positive and negative frequencies
                Imf(nOme)=0.0
                do iOm=2,nOme
                   Imf(nOme+iOm-1) = dreal(ImBub(iorb,jorb,iOm))
                   Imf(nOme-iOm+1) = -dreal(ImBub(iorb,jorb,iOm))
                enddo
                ! Perform Kramars-Kronig for the real part
                !$OMP  PARALLEL DO SHARED(Imf,nOme,OmKK,Ref) PRIVATE(rf,iOm)
                do iOm=2,nOme
                   CALL kramarskronig(rf, Imf, OmKK, nOme+iOm-1, 2*nOme-1)
                   Ref(iOm) = rf
                enddo
                !$OMP END PARALLEL DO
                !
                ! zero frequency value done here, because it can not be obtained by subroutine
                do iOm=1,nOme-1
                   Imf(iOm) = Imf(iOm)/OmKK(iOm)
                enddo
                do iOm=nOme+1,2*nOme-1
                   Imf(iOm) = Imf(iOm)/OmKK(iOm)
                enddo
                Imf(nOme)=0.5*(Imf(nOme-1)+Imf(nOme+1))
                CALL integrate_trapz(rf, Imf, OmKK, 2*nOme-1)
                Ref(1) = rf/PI
                ! Saving the complex function result into Bub
                do iOm=1,nOme
                   Bub(iorb,jorb,iOm) = Ref(iOm) + IMAG*ImBub(iorb,jorb,iOm)
                enddo
             else
                print *, 'Imaginary Kramars-Kronig', iorb, jorb
                Imc(nOme)=0.0
                do iOm=2,nOme
                   Imc(nOme+iOm-1) = ImBub(iorb,jorb,iOm)
                   Imc(nOme-iOm+1) = -ImBub(iorb,jorb,iOm)
                enddo
                ! Perform Kramars-Kronig for the real part
                !$OMP  PARALLEL DO SHARED(Imc,nOme,OmKK,Rec) PRIVATE(cf,iOm)
                do iOm=2,nOme
                   CALL complex_kramarskronig(cf, Imc, OmKK, nOme+iOm-1, 2*nOme-1)
                   Rec(iOm) = cf
                enddo
                !$OMP END PARALLEL DO
                !
                ! zero frequency value done here, because it can not be obtained by subroutine
                do iOm=1,nOme-1
                   Imc(iOm) = Imc(iOm)/OmKK(iOm)
                enddo
                do iOm=nOme+1,2*nOme-1
                   Imc(iOm) = Imc(iOm)/OmKK(iOm)
                enddo
                Imc(nOme)=0.5*(Imc(nOme-1)+Imc(nOme+1))
                CALL complex_integrate_trapz(cf, Imc, OmKK, 2*nOme-1)
                Rec(1) = cf/PI
                ! Saving the result
                do iOm=1,nOme
                   Bub(iorb,jorb,iOm) = Rec(iOm) + IMAG*ImBub(iorb,jorb,iOm)
                enddo
             endif
          enddo
       enddo
       
       ! Writing out the result
       WRITE(skii,fmt='(I5)') (iQ-1)
       FNAME=TRIM(fileout)//ADJUSTL(TRIM(skii))
       open(fh_out, FILE=FNAME, STATUS='unknown')
       do iOm=1,nOme
          WRITE(fh_out,'(F12.6)',advance='no') Ome(iOm)
          do iorb=1,norb
             do jorb=1,norb
                WRITE(fh_out,'(F14.6,1x,F14.6,3x)',advance='no') dreal(Bub(iorb,jorb,iOm)),dimag(Bub(iorb,jorb,iOm))
             enddo
          enddo
          WRITE(fh_out,*)
       enddo
       close(fh_out)
       ! Timings
       call cputim(t2c)
       call walltim(t2w)
       time_Q_c=time_Q_c+t2c-t1c
       time_Q_w=time_Q_w+t2w-t1w
       WRITE(*,'(A,I4,A,f10.4,A,f10.4)') 'Q=', iQ, ' time[s]=', t2c-t1c, ' walltime[s]', t2w-t1w
    enddo

    deallocate( OmKK )
    deallocate( Imf, Ref )
    deallocate( Imc, Rec )

    
    deallocate( om_idx, dx )
    deallocate( ImBub, Bub )
    WRITE(*,'(A,f15.4,A,f15.4)') 'Total-time [in min] to calculate Bubble=', time_Q_c/60., ' Total-walltime', time_Q_w/60.
    
  end SUBROUTINE CmpChiQ0_real

  SUBROUTINE RealBubble__Destruct__
    DEALLOCATE( Ome , chi0r0 )
  END SUBROUTINE RealBubble__Destruct__
end module RealBubble


program RealAxisBubble
  use LinLogMesh
  use Qlist
  use Klist
  use greenk
  use RealBubble
  IMPLICIT NONE
  CHARACTER*200 :: filemesh, filegkr, fbubbleReal, fileQlist, fileKlist, fileout
  INTEGER :: nargs, j, iom
  INTEGER,PARAMETER :: fhi = 995
  CHARACTER*100, allocatable :: argv(:) ! Command-line arguments  
  
  nargs = iargc()
  if (nargs.LT.1) then
     print *, 'Missing the input file to proceed. Create input file with the following lines'
     print *, 'Ba122.klist    # filename with k-list'
     print *, 'Qlist.dat      # filename with Qlist'
     print *, 'rmesh.dat      # real axis mesh'
     print *, 'G_k1r_         # file with real axis k-dependent Grens function'
     print *, 'G_local1r_     # file with real axis local Grens function'
     print *, 'chi0_real.     # name of the output Bubble on real axis'
     call exit(1)
  endif
  ALLOCATE (argv(nargs))                                 ! wouldn't like to parse                                                                             
  WRITE(*,'(A,I2)') 'nargs=', nargs
  DO j=1,nargs
     CALL getarg(j, argv(j))
     WRITE(*,'(A,A)') 'argi=', TRIM(argv(j))
  ENDDO
  
  open(fhi, FILE=argv(1),status='old', ERR=900, form='formatted')
  READ(fhi,*)   fileKlist   ! k-list filename, Ba122.klist
  READ(fhi,*)   fileQlist   ! Q-list
  READ(fhi,*)   filemesh    ! real axis mesh
  READ(fhi,*)   filegkr     ! Green's function G_k on real axis
  READ(fhi,*)
  READ(fhi,*)   fileout     ! output filename, for example  'fchiQ.'
  close(fhi)
  
  ! Reading real axis mesh 
  CALL LinLogMesh__Init__(filemesh)
  ! Generates real axis mesh
  CALL LinLogMesh_Generate()
  
  ! Reading Q list
  CALL Qlist__Init__(fileQlist)
  CALL Klist__Init__(fileKlist)
  CALL greenk__Init__(filegkr)
  if ( sum(abs(om-oml)).gt.1e-5 ) then
     print *, 'Mesh in greens function and real-mesh are not compatible. Diff=', sum(abs(om-oml))
     do iom=1,noml
        print *, iom, om(iom), oml(iom)
     enddo
  endif
  CALL RealBubble__Init__()
  
  CALL CmpChiQ0_real(fileout)

  CALL RealBubble__Destruct__()
  CALL greenk__Destruct__()
  CALL LinLogMesh__Destruct__()
  CALL Qlist__Destruct__()
  CALL Klist__Destruct__()

  deallocate( argv )
  STOP
900 print *, 'ERROR opening ',TRIM(argv(1)),' file'
end program RealAxisBubble
