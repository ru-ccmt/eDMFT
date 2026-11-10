      module parallel
!       ICTXT     - Context of processors holding matrices H and S
!       ICTXTALL  - Context including all processors
!       DESCHS(9) - BLACS array descriptor of H and S
!       DESCZ(9)  - BLACS array descriptor of Z
!       NPE       - Number of processing elements
!       MYID      - Number of current processor
!       NPROW     - Number of process rows
!       NPCOL     - Number of process columns 
!                   (NPROW x NPCOL chosen 'as square as possible')
!       MYROWHS   - Row number of current processor
!       MYCOLHS   - Column number of current processor
!       BLOCKSIZE - Blocksize for two dimensional distribution of H, S, and Z,
!                   used for both dimensions,
!                   also used as algorithmic block size in ScaLapack
!       LDHS      - Leading dimension of local portion of H, S, and Z
!       LCOLHS    - Number of local columns of H and S
!       LCOLZ     - Number of local columns of Z
        INTEGER :: ICTXT, ICTXTALL
        INTEGER :: DESCHS(9), DESCZ(9), DESCHD(9)
        INTEGER :: DESCLOCAL(9)
        INTEGER :: NPE, MYID, NPROW, NPCOL
        INTEGER :: MYROWHS, MYCOLHS
        INTEGER :: MYROWPANEL, MYCOLPANEL
        INTEGER :: BLOCKSIZE, BLOCKIN, DBLOCKIN, DBLOCKSIZE
        INTEGER :: LDHS, LCOLHS, LCOLZ, DLDHS, DLCOLHS
        INTEGER :: LLD, myrow, mycol, info
        INTEGER :: NBELW_PARA               !! added RUH, July 16
        LOGICAL :: elpa_switch,scala_switch
#ifdef ELPA
        LOGICAL :: elpa_mem_default
#endif
        CHARACTER*3 :: gridshape

      CONTAINS
        SUBROUTINE INIT_PARALLEL
          IMPLICIT NONE
#ifdef Parallel
          include 'mpif.h'
#endif
          INTEGER :: IERR,i,j,p,q
          INTEGER :: iarg, iargc, iunit
          character*180 :: FNAME, DEFFN
          character*80 :: line
#ifdef Parallel
          call MPI_INIT(IERR)
          call MPI_COMM_SIZE( MPI_COMM_WORLD, NPE, IERR)
          call MPI_COMM_RANK( MPI_COMM_WORLD, MYID, IERR)
          CALL BARRIER
          CALL SL_INIT(ICTXTALL, 1, NPE)
#endif
! two-dimensional distribution for global matrix computation
!         NPCOL=SQRT(real(NPE))
!         NPROW=NPE/NPCOL
          iarg=iargc()
          if(iarg.lt.1) then
             STOP 'Exactly one commandline argument has to be given. &
                   Error in GTFNAM'
          endif
          iarg=1
          CALL GETARG(iarg,DEFFN)
          open(1,FILE=DEFFN)
          do 
            read(1,*) iunit, FNAME
            if(iunit.eq.5) exit
          enddo
          close(1)
          open(5,FILE=FNAME)
          gridshape='old'
          elpa_switch=.false.
          scala_switch=.false.
#ifdef ELPA
          elpa_mem_default=.true.
#endif          
          read(5,*) 
          read(5,'(a)') line
          BLOCKIN=0 ; DBLOCKIN=0
          do i=10,80
             if(line(i:i).eq.'(') exit
             if(line(i:i+1).eq.'bl'.or.line(i:i+1).eq.'BL') then
                 read(line(i+2:i+4),'(i3)',err=1212) BLOCKIN
                 if(BLOCKIN.eq.999) goto 1234
                 if(BLOCKIN.lt.16) BLOCKIN=16
                 if(BLOCKIN.gt.384) BLOCKIN=384
             endif
  1234           if(line(i:i+1).eq.'bd'.or.line(i:i+1).eq.'BD') then
                 read(line(i+2:i+4),'(i3)',err=1212) DBLOCKIN
                 if(DBLOCKIN.eq.999) goto 1212
                 if(DBLOCKIN.lt.16) DBLOCKIN=16
                 if(DBLOCKIN.gt.384) DBLOCKIN=384
             endif
   1212      continue    
!             if(line(i:i+3).eq.'ELPA') elpa_switch=.true.
             if(line(i:i+4).eq.'SCALA') scala_switch=.true.
             if(line(i:i).eq.'x') gridshape=line(i-1:i+1)
#ifdef ELPA
             if(line(i:i+3).eq.'ELPA') elpa_switch=.true.
             if(line(i:i).eq.'m') then
                if(line(i-1:i).eq.'lm') then
                   elpa_mem_default=.false.
                elseif(line(i-1:i).eq.'hm') then
                   elpa_mem_default=.true.
                endif
             endif
#endif
          enddo
          if((.not. elpa_switch).and.(.not. scala_switch)) then
#ifndef ELPA
            scala_switch=.true.
#else
            elpa_switch=.true.
#endif
          endif
          rewind 5
          close (5)
          do p=1,1000
             do q=1,p
                if (p*q.eq.npe) then
!                   npcol=i   !j   modification March 2016 Ruh
!                   nprow=j   !i
                   if((gridshape.eq.'pxq').or.(gridshape.eq.'old')) then
                      nprow=p
                      npcol=q
                   endif
                   if(gridshape.eq.'qxp') then
                      nprow=q
                      npcol=p
                   endif
                   goto 300
                endif
             enddo
          enddo
300       continue
#ifdef Parallel
          CALL SL_INIT(ICTXT, NPROW,npcol)
#else
          NPE=1
          MYID=0
          NPCOL = 1
          NPROW = 1
#endif
        END SUBROUTINE INIT_PARALLEL

        SUBROUTINE INIT_PARALLELMATRICES(N,NUME2)
          IMPLICIT NONE
          INTEGER :: N, NUME2
#ifdef Parallel
          include 'mpif.h'
          include 'param.inc'

          INTEGER :: IERR, NPROWHS, NPCOLHS, KONE1
          INTEGER :: NPROWPANEL, NPCOLPANEL
          INTEGER, EXTERNAL :: NUMROC
          integer i,j

          INTEGER I1, I2, II1, II2, IBB, IBL, Excess, EXX
!         INTEGER IB1, IB2
          real*8  ENICE, ENOTNICE
          INTEGER IBB_NICE
!
!         Minimum sizes without blocking
          I1=N/NPROW
          I2=N/NPCOL
          IF(I1*NPROW .lt. N)I1=I1+1
          IF(I2*NPCOL .lt. N)I2=I2+1
!
!         Search for the best Blocksize between 64 and 256
          EXCESS=I1*I2
!_REAL      DO IBL=128,256,2
!_COMPLEX      DO IBL=64,128,2
!          DO IBL=64,384,4
!               How many blocks will fit
                II1=I1/IBL
                if(II1*IBL .lt. I1)II1=II1+1
                II2=I2/IBL
                if(II2*IBL .lt. I2)II2=II2+1
!               How much excess do we have
                EXX=II2*IBL*II1*IBL-I1*I2
                IF(EXX .lt. EXCESS)then
                        IBB=IBL
                        EXCESS=EXX
                endif
        ENDDO
        ENICE=100.d0*EXCESS/dble(I1*I2)
        IBB_NICE=IBB
        write(6,67)IBB_NICE,ENICE
        67      format('Optimum Blocksize for setup',i4,' Excess % ',D10.3 )
!       
!         Try some others between 16 and 48
          EXCESS=I1*I2
          DO IBL=16,48,2
!               How many blocks will fit
                II1=I1/IBL
                if(II1*IBL .lt. I1)II1=II1+1
                II2=I2/IBL
                if(II2*IBL .lt. I2)II2=II2+1
!               How much excess do we have
                EXX=II2*IBL*II1*IBL-I1*I2
                IF(EXX + 100 .lt. EXCESS)then
                        IBB=IBL
                        EXCESS=EXX
                endif
        ENDDO
        ENOTNICE = 100.d0*EXCESS/dble(I1*I2)
        write(6,68)IBB,ENOTNICE
68      format('Optimum Blocksize for diag',i4,' Excess % ',D10.3)
!
!       We could use the nice size if the difference is not much (?)
!       if((IBB .gt. IBB_NICE) .and. (ENICE-ENOTNICE .gt. 0.5))IBB=IBB_NICE
!          BLOCKSIZE = IBB
 
        if(BLOCKIN.eq.0) THEN
!_REAL    BLOCKSIZE = BREAL
!_COMPLEX BLOCKSIZE = BCOMPLEX
        else if(blockin.eq.999) THEN
          BLOCKSIZE=IBB_NICE
        ELSE
          BLOCKSIZE=BLOCKIN
        ENDIF
        IF(DBLOCKIN.eq.0) THEN
          if(elpa_switch)then
!_REAL          DBLOCKSIZE = BELPA
!_COMPLEX       DBLOCKSIZE = BELPAC
          else
!_REAL          DBLOCKSIZE = BREAL
!_COMPLEX       DBLOCKSIZE = BCOMPLEX
          endif
        else if(dblockin.eq.999) THEN
          DBLOCKSIZE=IBB
        else 
          DBLOCKSIZE=DBLOCKIN
        endif

!          BLOCKSIZE = 64
!          DBLOCKSIZE = BLOCKSIZE
!print*,'blocksize',n,dblocksize,mod(n,dblocksize),mod(n,dblocksize*2)
          if(mod(n, DBLOCKSIZE).eq.1)  then
                  DBLOCKSIZE=IBB
                  write(6,*) 'Blocksize for Diagonalization changed'
          endif        
          write(6,69)BLOCKSIZE, DBLOCKSIZE
69        FORMAT('Base Blocksize ',i4,' Diagonalization ',i4)

!do not use IBB, because this leads to wrong eigenvectors for Au-atom, PB 2014
!!         BLOCKSIZE = IBB_nice
#if defined(INTEL_VML)
          call mkl_free_buffers()
#endif
! one-dimensional distribution for local matrix setup
          CALL BLACS_GRIDINFO(ICTXTALL, NPROWPANEL, NPCOLPANEL, &
                              MYROWPANEL, MYCOLPANEL)
          IF(MYROWPANEL.NE.-1) THEN
            CALL DESCINIT(DESCLOCAL, N, N,  &
                          N, N,0,0,ICTXTALL, N, IERR)
          END IF

          call BLACS_GRIDINFO(ICTXT, NPROW, NPCOL, MYROWHS, MYCOLHS)
          LDHS = max(1,NUMROC(N, BLOCKSIZE, MYROWHS, 0, NPROW))
          LCOLHS = max(1,NUMROC(N, BLOCKSIZE, MYCOLHS, 0, NPCOL))
          LCOLZ = max(1,NUMROC(N, BLOCKSIZE, MYCOLHS, 0, NPCOL))
          CALL DESCINIT(DESCHS, N, N, BLOCKSIZE, BLOCKSIZE, 0, 0, ICTXT, LDHS, IERR)
          CALL DESCINIT(DESCZ, N, N, BLOCKSIZE, BLOCKSIZE, 0, 0, ICTXT, LDHS, IERR)
!         For Diagonalization
          DLDHS   = max(1,NUMROC(N, DBLOCKSIZE, MYROWHS, 0, NPROW))
          DLCOLHS = max(1,NUMROC(N, DBLOCKSIZE, MYCOLHS, 0, NPCOL))
          CALL DESCINIT(DESCHD, N, N, DBLOCKSIZE, DBLOCKSIZE, 0, 0, ICTXT, DLDHS, IERR)

#else
        BLOCKSIZE = 64
!use larger blocksize to fix acml problem ?
        BLOCKSIZE = 128
        ICTXT = -1
        LDHS = N
        LCOLHS = N
        LCOLZ = NUME2
        DESCHS(2) = -1
        DESCZ(2) = -1
#endif
        END SUBROUTINE INIT_PARALLELMATRICES

        SUBROUTINE BARRIER
#ifdef Parallel
          IMPLICIT NONE
          include 'mpif.h'
          INTEGER :: IERR
          CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
#endif
        END SUBROUTINE BARRIER


        SUBROUTINE CLOSE_PARALLEL
#ifdef Parallel
          IMPLICIT NONE
          include 'mpif.h'
          INTEGER :: IERR
          CALL MPI_FINALIZE(IERR)
#endif
        END SUBROUTINE CLOSE_PARALLEL

        SUBROUTINE ABORT_PARALLEL
#ifdef Parallel
        IMPLICIT NONE
        include 'mpif.h'
        INTEGER :: IERR
        call flush(6)
!        close(unit=6)
        if(myid .eq. 0)then
                call flush(99)
!                close(unit=99)
        endif
        CALL MPI_ABORT(MPI_COMM_WORLD, IERR)
#endif
        END SUBROUTINE ABORT_PARALLEL

      end module parallel

!--------------------------------------------------------------

      module lapw_timer
!  64 Timers are available in ScaLapack
      double precision :: cpu_time(64)
      double precision :: old_cpu_time(64)
      double precision :: wall_time(64)
      double precision :: old_wall_time(64)
!$omp threadprivate(cpu_time,old_cpu_time,wall_time,old_wall_time)
      integer, parameter :: time_total=1,  &
                   time_hamilt=10, time_hns=11, time_diag=12, time_horb=13,  &
                 time_coor=15,  &
                 time_setwar=17, &
                 time_pkp=18, &
                 time_cholesky=20, time_transform=21, time_sep=22,  &
                 time_backtransform=23, time_transpose_h_s=24, time_pwpw=25, &
                 time_pwlo=26,  time_lolo=27,  time_hns30=28, time_matrixcomplete=29, &
                 time_iter=30,time_H1=31, time_HxF=32, &
                 time_F=33, time_Ffac=34, time_Fleq=35, &
                 time_SxF=36, time_projH=37, time_projH_seq=38, time_backproject=39,&
                 time_albl=40, time_loop260=41, time_phase=42,  &
                 time_legendre=43, time_step_function=44, time_h=45,  &
                 time_us=46, time_loop210=47, time_loop230=48,  &
                 time_loop240=49, time_overlap=50, time_lo=51, &
                 time_g1=55, time_atpar=56, time_distiouter=57,&
                 time_horb_setup=58, time_sync=59
      logical    usecputime, usewalltime
      contains
        subroutine init_all_timer
!!$#ifdef Parallel
!!$          use parallel, only : ICTXTALL
!!$          IMPLICIT NONE
!!$          external slboot
!!$          double precision  :: WTIME
!!$          CALL SLBOOT
!!$          CALL SLCOMBINE( ICTXTALL, 'All', '>', 'C', 1, 1, WTIME )
!!$          if(WTIME.eq.-1.0D0) then
!!$            usecputime=.false.
!!$            CALL SLENABLE
!!$          else
!!$            usecputime=.true.
!!$          end if
!!$          CALL SLCOMBINE( ICTXTALL, 'All', '>', 'W', 1, 1, WTIME )
!!$          if(WTIME.eq.-1.0D0) then
!!$            usewalltime=.false.
!!$            CALL SLENABLE
!!$          else
!!$            usewalltime=.true.
!!$          end if
!!$#endif
          cpu_time(:)=0.0D0
          wall_time(:)=0.0D0
        end subroutine init_all_timer

        subroutine init_timer(n)
          IMPLICIT NONE
          integer :: n
          cpu_time(n) = 0.0D0
          wall_time(n) = 0.0D0
        end subroutine init_timer

        subroutine start_timer(n)
          integer :: n
!!$#ifdef Parallel
!!$          external sltimer
!!$          call SLTIMER(n)
!!$          if(.not.usecputime) then
!!$            call CPUTIM(old_cpu_time(n))
!!$            cpu_time(n)=0.0D0
!!$          end if
!!$          if(.not.usewalltime) then
!!$            call WALLTIM(old_wall_time(n))
!!$            wall_time(n)=0.0D0
!!$          end if
!!$#else
          call CPUTIM(old_cpu_time(n))
          cpu_time(n)=0.0D0
          call WALLTIM(old_wall_time(n))
          wall_time(n)=0.0D0
!!$#endif
        end subroutine start_timer

        subroutine stop_timer(n)
          integer :: n
          double precision :: current_time
!!$#ifdef Parallel
!!$          external sltimer
!!$          call SLTIMER(n)
!!$          if(.not.usecputime) then
!!$            call CPUTIM(current_time)
!!$            cpu_time(n) = cpu_time(n) + current_time - old_cpu_time(n)
!!$          end if
!!$          if(.not.usewalltime) then
!!$            call WALLTIM(current_time)
!!$            wall_time(n) = wall_time(n) + current_time - old_wall_time(n)
!!$          end if
!!$#else
          call CPUTIM(current_time)
          cpu_time(n) = cpu_time(n) + current_time - old_cpu_time(n)
          call WALLTIM(current_time)
          wall_time(n) = wall_time(n) + current_time - old_wall_time(n)
!!$#endif
        end subroutine stop_timer
         
        double precision function read_cpu_time(n)
          use parallel, only : ICTXTALL
          integer :: n
          double precision :: WTIME
!!$#ifdef Parallel
!!$          external slcombine
!!$          if(usecputime) then
!!$            CALL SLCOMBINE( ICTXTALL, 'All', '>', 'C', 1, N, WTIME )
!!$            read_cpu_time=WTIME
!!$          else
!!$            read_cpu_time=cpu_time(n)
!!$          end if
!!$#else
          read_cpu_time=cpu_time(n)
!!$#endif

          return
        end function read_cpu_time

        double precision function read_wall_time(n)
          use parallel, only : ICTXTALL
          integer :: n
          double precision :: WTIME
!!$#ifdef Parallel
!!$          external slcombine
!!$          if(usewalltime) then
!!$            CALL SLCOMBINE( ICTXTALL, 'All', '>', 'W', 1, N, WTIME )
!!$            read_wall_time=WTIME
!!$          else
!!$            read_wall_time=wall_time(n)
!!$          end if
!!$#else
          read_wall_time=wall_time(n)
!!$#endif

        end function read_wall_time
      end module lapw_timer
          
!--------------------------------------------------------------

        module reallocate
          !     nur 1 (generischer) Name wird von aussen angesprochen
          interface doreallocate
            module procedure doreallocate_r8_d1
            module procedure doreallocate_r8_d2
            module procedure doreallocate_r8_d3
            module procedure doreallocate_i4_d1
            module procedure doreallocate_i4_d2
            module procedure doreallocate_i4_d3
            module procedure doreallocate_warp_c16_d2x3
            module procedure doreallocate_c3_d1
            module procedure doreallocate_z16_d2
            module procedure doreallocate_z16_d3
            module procedure doreallocate_z32_d3
            module procedure doreallocate_c10_d1
            module procedure hugo     !   ;)
          end interface
        contains

          !     leider sind mehrere subroutines notwendig fuer verschiedene Typen
          subroutine doreallocate_r8_d1(tf, newdimension)
            implicit none
            real*8, pointer :: hilfsfeld(:), tf(:)
            integer min1, newdimension
            allocate(hilfsfeld(newdimension))
            hilfsfeld=0.0D0
            !     nur 1 mal kopieren reicht
            !     auch fuer mehrdimensionale Felder schaut die Zuweisung gleich aus
            min1=min(newdimension,size(tf,1))
            hilfsfeld(1:min1)=tf(1:min1)
            deallocate(tf)
            !     der Zeiger wird nur auf das neue Feld umgebogen, nicht neu alloziert
            tf=>hilfsfeld
          end subroutine

          subroutine doreallocate_r8_d2(tf, newdimension1, newdimension2)
            implicit none
            real*8, pointer :: hilfsfeld(:,:), tf(:,:)
            integer min1, min2, newdimension1, newdimension2
            allocate(hilfsfeld(newdimension1,newdimension2))
            hilfsfeld=0.0D0
            min1=min(newdimension1,size(tf,1))
            min2=min(newdimension2,size(tf,2))
            hilfsfeld(1:min1,1:min2)=tf(1:min1,1:min2)
            deallocate(tf)
            tf=>hilfsfeld
          end subroutine

          subroutine doreallocate_r8_d3(tf, newdimension1, newdimension2, newdimension3)
            implicit none
            real*8, pointer :: hilfsfeld(:,:,:), tf(:,:,:)
            integer min1, min2, min3, newdimension1, newdimension2, newdimension3
            allocate(hilfsfeld(newdimension1,newdimension2,newdimension3))
            hilfsfeld=0.0D0
            min1=min(newdimension1,size(tf,1))
            min2=min(newdimension2,size(tf,2))
            min3=min(newdimension3,size(tf,3))
            hilfsfeld(1:min1,1:min2,1:min3)=tf(1:min1,1:min2,1:min3)
            deallocate(tf)
            tf=>hilfsfeld
          end subroutine

          subroutine doreallocate_i4_d1(tf, newdimension)
            implicit none
            integer*4, pointer :: hilfsfeld(:), tf(:)
            integer min1, newdimension
            allocate(hilfsfeld(newdimension))
            hilfsfeld=0
            min1=min(newdimension,size(tf,1))
            hilfsfeld(1:min1)=tf(1:min1)
            deallocate(tf)
            tf=>hilfsfeld
          end subroutine

          subroutine doreallocate_i4_d2(tf, newdimension1, newdimension2)
            implicit none
            integer*4, pointer :: hilfsfeld(:,:), tf(:,:)
            integer min1, min2, newdimension1, newdimension2
            allocate(hilfsfeld(newdimension1,newdimension2))
            hilfsfeld=0
            min1=min(newdimension1,size(tf,1))
            min2=min(newdimension2,size(tf,2))
            hilfsfeld(1:min1,1:min2)=tf(1:min1,1:min2)
            deallocate(tf)
            tf=>hilfsfeld
          end subroutine

          subroutine doreallocate_i4_d3(tf, newdimension1, newdimension2, newdimension3)
            implicit none
            integer*4, pointer :: hilfsfeld(:,:,:), tf(:,:,:)
            integer min1, min2, min3, newdimension1, newdimension2, newdimension3
            allocate(hilfsfeld(newdimension1,newdimension2,newdimension3))
            hilfsfeld=0
            min1=min(newdimension1,size(tf,1))
            min2=min(newdimension2,size(tf,2))
            min3=min(newdimension3,size(tf,3))
            hilfsfeld(1:min1,1:min2,1:min3)=tf(1:min1,1:min2,1:min3)
            deallocate(tf)
            tf=>hilfsfeld
          end subroutine

          subroutine doreallocate_warp_c16_d2x3(tf, n1a, n1b, n2a, n2b, n3a, n3b)
            implicit none
            complex*16, pointer :: hilfsfeld(:,:,:), tf(:,:,:)
            integer min1, n1a, n1b, n2a, n2b, n3a, n3b
            allocate( hilfsfeld(-n1b:n1b, -n2b:n2b, -n3b:n3b) )
            hilfsfeld(-n1b:n1b, -n2b:n2b, -n3b:n3b)=(0.0D0,0.0D0)
            hilfsfeld(-n1a:n1a,-n2a:n2a,-n3a:n3a)=tf(-n1a:n1a,-n2a:n2a,-n3a:n3a)
            deallocate(tf)
            tf=>hilfsfeld
          end subroutine

          subroutine doreallocate_c3_d1(tf, newdimension1, dummy)
            implicit none
            character*3, pointer :: hilfsfeld(:), tf(:)
            integer :: min1, newdimension1
            double precision :: dummy
            allocate(hilfsfeld(newdimension1))
            min1=min(newdimension1,size(tf,1))
            hilfsfeld(1:min1)=tf(1:min1)
            deallocate(tf)
            tf=>hilfsfeld
          end subroutine

          subroutine doreallocate_z16_d2(tf, newdimension1, newdimension2)
            implicit none
            complex*16, pointer :: hilfsfeld(:,:), tf(:,:)
            integer :: min1, min2, newdimension1,newdimension2
 !          double precision :: dummy
            allocate(hilfsfeld(newdimension1,newdimension2))
            hilfsfeld=(0.0D0,0.0D0)
            min1=min(newdimension1,size(tf,1))
            min2=min(newdimension2,size(tf,2))
            hilfsfeld(1:min1,1:min2)=tf(1:min1,1:min2)
            deallocate(tf)
            tf=>hilfsfeld
          end subroutine
          
          subroutine doreallocate_z16_d3(tf, newdimension1, newdimension2, newdimension3)
            implicit none
            complex*16, pointer :: hilfsfeld(:,:,:), tf(:,:,:)
            integer :: min1, min2, min3
            integer :: newdimension1, newdimension2, newdimension3
!          double precision :: dummy
            allocate(hilfsfeld(newdimension1,newdimension2,newdimension3))
            hilfsfeld=(0.0D0,0.0D0)
            min1=min(newdimension1,size(tf,1))
            min2=min(newdimension2,size(tf,2))
            min3=min(newdimension2,size(tf,3))
            hilfsfeld(1:min1,1:min2,1:min3)=tf(1:min1,1:min2,1:min3)
            deallocate(tf)
            tf=>hilfsfeld
          end subroutine doreallocate_z16_d3

          subroutine doreallocate_z32_d3(tf, newdimension1, newdimension2, newdimension3)
            implicit none
            complex*32, pointer :: hilfsfeld(:,:,:), tf(:,:,:)
            integer :: min1, min2, min3
            integer :: newdimension1, newdimension2, newdimension3
!          double precision :: dummy
            allocate(hilfsfeld(newdimension1,newdimension2,newdimension3))
            hilfsfeld=(0.0Q0,0.0Q0)
            min1=min(newdimension1,size(tf,1))
            min2=min(newdimension2,size(tf,2))
            min3=min(newdimension2,size(tf,3))
            hilfsfeld(1:min1,1:min2,1:min3)=tf(1:min1,1:min2,1:min3)
            deallocate(tf)
            tf=>hilfsfeld
          end subroutine doreallocate_z32_d3

          subroutine doreallocate_c10_d1(tf, newdimension1)
            implicit none
            character*10, pointer :: hilfsfeld(:), tf(:)
            integer min1, newdimension1
            allocate(hilfsfeld(newdimension1))
            min1=min(newdimension1,size(tf,1))
            hilfsfeld(1:min1)=tf(1:min1)
            deallocate(tf)
            tf=>hilfsfeld
          end subroutine

        !     Es gibt auch Methoden, um das Programm unleserlich zu machen :-)
        !     das sollten wir besser vermeiden ;-)
          subroutine hugo
            write(*,*) " Hier ist Hugo"
          end subroutine
        end module

!--------------------------------------------------------------

      module matrices
        INCLUDE 'param.inc'
        integer hsrows
!
!      sequential code:
!        HS - Hamilton matrix (hermitian matrix, stored below diagonal 
!                                                in conjugated form)
!           - Overlap matrix  (hermitian matrix, stored above diagonal)
!        HSDIAG - Diagonal elements of Overlap matrix
!               - several times swapped with 
!                 diagonal elements of Hamilton matrix
!      parallel code:
!        H - Hamilton matrix (hermitian matrix, stored above diagonal)
!        S - Overlap matrix (hermitian matrix, stored above diagonal)
!      HS, HSDIAG, H, S are set up in HAMILT, HNS 
!                       and solved in SECLR4, SECLR5
#ifdef Parallel
!_REAL      DOUBLE PRECISION, allocatable ::  S(:, :)
!_REAL      DOUBLE PRECISION, allocatable ::  H(:, :)
!_COMPLEX      DOUBLE COMPLEX, allocatable :: S(:, :)
!_COMPLEX      DOUBLE COMPLEX, allocatable :: H(:, :)
#else
!_REAL      DOUBLE PRECISION, allocatable ::  HS(:, :)
!_REAL      DOUBLE PRECISION, allocatable ::  HSDIAG(:)
!_COMPLEX      DOUBLE COMPLEX, allocatable :: HS(:, :)
!_COMPLEX      DOUBLE COMPLEX, allocatable :: HSDIAG(:)
#endif
!_REAL      DOUBLE PRECISION, allocatable ::  Z(:,:)
!_COMPLEX      DOUBLE COMPLEX, allocatable :: Z(:, :)
        DOUBLE PRECISION, allocatable ::  EIGVAL(:)
        DOUBLE PRECISION, pointer ::  XK(:), YK(:), ZK(:)
        DOUBLE PRECISION, pointer ::  RK(:)
        INTEGER, pointer :: KZZ(:,:)
        CONTAINS
        SUBROUTINE INIT_MATRICES( N, NUME2,nmatmax1)
        use parallel, only : DESCHS, DESCZ, NPE, MYID, &
                             NPROW, NPCOL, MYROWHS, MYCOLHS, &
                             BLOCKSIZE, LDHS, LCOLHS, LCOLZ, &
                             INIT_PARALLELMATRICES
        use reallocate, only: doreallocate
        IMPLICIT NONE        
        INTEGER N, NUME2,nmatmax1
        include 'param.inc'
        INTEGER INFO
        nmatmax1=nmatmax
!_COMPLEX       nmatmax1=nmatmax1/sqrt(2.d0)
#ifdef Parallel
        nmatmax1=nmatmax1*sqrt(NPE*1.d0)
#endif          

        if( n.eq.0 .and. nume2.eq.0) then
          allocate( RK(NMATMAX1+1) )
          allocate( XK(NMATMAX1+1) )
          allocate( YK(NMATMAX1+1) )
          allocate( ZK(NMATMAX1+1) )
          allocate( KZZ(3,NMATMAX1+1) )
          RK=0.0D0
          XK=0.0D0
          YK=0.0D0
          ZK=0.0D0
          KZZ=0.0D0
        else
          if(myid.eq.0) write(6,*) 'Matrix size ', N
          HSROWS = N
!          write(6,*) 'Initialize parallelization variables'
          CALL INIT_PARALLELMATRICES(N, NUME2)
          call doreallocate( RK,HSROWS+1 )
          call doreallocate( XK,HSROWS+1 )
          call doreallocate( YK,HSROWS+1 )
          call doreallocate( ZK,HSROWS+1 )
          call doreallocate( KZZ,3,HSROWS+1 )
#ifdef Parallel
          allocate( H(LDHS, LCOLHS))        
          allocate( S(LDHS, LCOLHS))
!_REAL     write(6,'(a20,f12.1,a,a,2i6)') 'allocate H',ldhs*lcolhs*8.0/(1024.0*1024.0),' MB ', '         dimensions',ldhs,lcolhs
!_REAL     write(6,'(a20,f12.1,a,a,2i6)') 'allocate S',ldhs*lcolhs*8.0/(1024.0*1024.0),' MB ', '         dimensions',ldhs,lcolhs
!_COMPLEX  write(6,'(a20,f12.1,a,a,2i6)') 'allocate H',ldhs*lcolhs*16.0/(1024.0*1024.0),' MB ', '         dimensions',ldhs,lcolhs
!_COMPLEX  write(6,'(a20,f12.1,a,a,2i6)') 'allocate S',ldhs*lcolhs*16.0/(1024.0*1024.0),' MB ', '         dimensions',ldhs,lcolhs
          H=0.0D0
          S=0.0D0
#else 
          allocate( HS(HSROWS, HSROWS) )
          allocate( HSDIAG( HSROWS ) )
          allocate( Z(HSROWS, NUME2) )
!_REAL     write(6,'(a20,f12.1,a)') 'allocate HS',HSROWS**2*8.0/(1024.0*1024.0),' MB '
!_COMPLEX  write(6,'(a20,f12.1,a)') 'allocate HS',HSROWS**2*16.0/(1024.0*1024.0),' MB '
!_REAL     write(6,'(a20,f12.1,a)') 'allocate Z ',HSROWS*NUME2*8.0/(1024.0*1024.0),' MB '
!_COMPLEX  write(6,'(a20,f12.1,a)') 'allocate Z ',HSROWS*NUME2*16.0/(1024.0*1024.0),' MB '
          HS=0.0D0
          HSDIAG=0.0D0
          Z=0.0D0
#endif          
          allocate( EIGVAL(HSROWS) )
          EIGVAL=0.0D0
        end if

        END SUBROUTINE INIT_MATRICES
        
        SUBROUTINE END_MATRICE1
          !if(allocated(rk)) deallocate( RK, XK, YK, ZK )
          !deallocate( RK, XK, YK, ZK )
#ifdef Parallel
          if (allocated(H)) deallocate( H)
          if (allocated(S)) deallocate( S)
#else
          deallocate( HS, HSDIAG )
#endif  
        END SUBROUTINE END_MATRICE1
        
        SUBROUTINE END_MATRICES
          deallocate( KZZ )
#ifdef Parallel
          deallocate( Z )
#else
          deallocate( Z )
#endif  
          deallocate( EIGVAL )
        END SUBROUTINE END_MATRICES

      end module matrices

!--------------------------------------------------------------

      module albl
!
!        AL(n,l,i) - coefficient al(kn) for atom i
!                    (including equivalent atoms)
!        BL(n,l,i) - coefficient bl(kn) for atom i
!                    (including equivalent atoms)
!
        DOUBLE PRECISION, allocatable ::   AL_r(:,:)
        DOUBLE PRECISION, allocatable ::   BL_r(:,:)
        DOUBLE PRECISION, allocatable ::   AL_c(:,:)
        DOUBLE PRECISION, allocatable ::   BL_c(:,:)

        DOUBLE PRECISION, allocatable ::   AL(:,:,:)
        DOUBLE PRECISION, allocatable ::   BL(:,:,:)
!$omp threadprivate(AL_r,BL_r,AL_c,BL_c,AL,BL)

      contains

        subroutine init_albl(lmax,SUB)

          use matrices, only : HSrows
          use parallel, only : LDHS, LCOLHS, blocksize

          integer  lmax
          character*3 sub

          if (sub.eq.'HAM') then
#ifdef Parallel
             allocate( AL_r(LDHS, 0:LMAX) )
             allocate( BL_r(LDHS, 0:LMAX) )
!             allocate( AL_c(blocksize, 0:LMAX) )
!             allocate( BL_c(blocksize, 0:LMAX) )
             allocate( AL_c(0:LMAX,blocksize) )
             allocate( BL_c(0:LMAX,blocksize) )
#else
             allocate( AL_r(HSROWS, 0:LMAX) )
             allocate( BL_r(HSROWS, 0:LMAX) )
             allocate( AL_c(0:LMAX,HSROWS) )
             allocate( BL_c(0:LMAX,HSROWS) )
!             allocate( AL_c(HSROWS, 0:LMAX) )
!             allocate( BL_c(HSROWS, 0:LMAX) )
#endif
          else
#ifdef Parallel
             ldm=max(LDHS,lcolhs)
             allocate( AL(ldm, 0:LMAX,2) )
             allocate( BL(ldm, 0:LMAX,2) )
#else
             allocate( AL(HSROWS, 0:LMAX,1) )
             allocate( BL(HSROWS, 0:LMAX,1) )
#endif
             AL=0.0D0
             BL=0.0D0
          endif
        end subroutine init_albl

        subroutine end_albl(sub)
          character*3 sub
          if (sub.eq.'HAM') then
             deallocate( AL_r, BL_r)
             deallocate( AL_c, BL_c)
          else
             deallocate( AL, BL)
          endif
        end subroutine end_albl

      end module albl

!--------------------------------------------------------------

      module loabc
!
!        ALO(l,i)    - local orbital coefficient alm for atom i
!        BLO(l,i)    - local orbital coefficient blm for atom i
!        CLO(l,i)    - local orbital coefficient clm for atom i
!        DPELO(l,i)  - derivative of the energy derivative ul_dot(r,E)
!                      of the local orbital radial-function ul(r,E) with
!                      respect to r
!                      (evaluated at r = muffin tin radius of atom i)
!        DPLO(l,i)   - derivative of the local orbital radial-function
!                      ul(r,E) with respect to radius r
!                      (evaluated at r = muffin tin radius of atom i)
!        ELO(l,i)    - local orbital expansion energy El for atom i
!        PEILO(l,i)  - the norm of the energy derivative ul_dot(r,E)
!                      of the local orbital radial function ul(r,E1)
!                      integrated over the muffin tin sphere
!        PELO(l,i)   - ul_dot(r,E)
!                      (evaluated at r = muffin tin radius of atom i)
!        PLO(l,i)    - local orbital radial-function ul(r,E)
!                      (evaluated at r = muffin tin radius of atom i)
!        PI12LO(l,i) - the product of the local orbital radial functions
!                      ul(r,E1) and ul(r,E2) integrated over the muffin
!                      tin sphere (for atom i)
!        PE12LO(l,i) - the product of the derivative of ul(r,E1) with
!                      respect to r and ul(r,E2) integrated over the
!                      muffin tin sphere (for atom i)
!
        integer nloat_new
        DOUBLE PRECISION, allocatable ::   ALO(:,:,:)
        DOUBLE PRECISION, allocatable ::   BLO(:,:,:)
        DOUBLE PRECISION, allocatable ::   CLO(:,:,:)
        DOUBLE PRECISION, allocatable ::   ELO(:,:,:)
        DOUBLE PRECISION, allocatable ::   PLO(:,:,:)
        DOUBLE PRECISION, allocatable ::   DPLO(:,:,:)
        DOUBLE PRECISION, allocatable ::   PELO(:,:,:)
        DOUBLE PRECISION, allocatable ::   DPELO(:,:,:)
        DOUBLE PRECISION, allocatable ::   PEILO(:,:,:)
        DOUBLE PRECISION, allocatable ::   PI12LO(:,:,:)
        DOUBLE PRECISION, allocatable ::   PE12LO(:,:,:)
        DOUBLE PRECISION, allocatable ::   PILOLO(:,:,:,:)
        LOGICAL, allocatable ::   LOSECDER(:,:,:)
      contains
        SUBROUTINE init_loabc(LOMAX, NLOAT, NATO)
          IMPLICIT NONE
          integer LOMAX, NLOAT, NATO
          allocate( ALO(0:LOMAX,nloat,NATO) )
          allocate( BLO(0:LOMAX,nloat,NATO) )
          allocate( CLO(0:LOMAX,nloat,NATO) )
          allocate( ELO(0:LOMAX,nloat,NATO) )
          allocate( PLO(0:LOMAX,nloat,NATO) )
          allocate( DPLO(0:LOMAX,nloat,NATO) )
          allocate( PELO(0:LOMAX,nloat,NATO) )
          allocate( DPELO(0:LOMAX,nloat,NATO) )
          allocate( PEILO(0:LOMAX,nloat,NATO) )
          allocate( PI12LO(0:LOMAX,nloat,NATO) )
          allocate( PE12LO(0:LOMAX,nloat,NATO) )
          allocate( PILOLO(0:LOMAX,nloat,nloat,NATO) )
          allocate( LOSECDER(0:LOMAX,nloat,NATO))
          ALO=0.0D0
          BLO=0.0D0
          CLO=0.0D0
          ELO=0.0D0
          PLO=0.0D0
          DPLO=0.0D0
          PELO=0.0D0
          DPELO=0.0D0
          PEILO=0.0D0
          PI12LO=0.0D0
          PE12LO=0.0D0
          PILOLO=0.0D0
          LOSECDER=.FALSE.
        END SUBROUTINE init_loabc
      end module loabc

!--------------------------------------------------------------

      module radial


!       A(NRAD)                   - r*u large component
!       A1(NRAD,NSLMAX)           - r*u large component
!       A1LO(NRAD,nloat,0:LOMAX)  - r*u large component for LO
!       AE(NRAD)                  - energy derivative of A
!       AE1(NRAD,NSLMAX)          - energy derivative of A1
!       AE1LO(NRAD,nloat,0:LOMAX) - energy derivative of A1LO
!                                   CURRENTLY NOT USED OR NECESSARY!
!       AP(NRAD)                  - r*u prime large component
!       B(NRAD)                   - r*u small component
!       B1(NRAD,NSLMAX)           - r*u small component
!       B1LO(NRAD,nloat,0:LOMAX)  - r*u small component for LO
!       BE(NRAD)                  -
!       BE1(NRAD,NSLMAX)          -
!       BE1LO(NRAD,nloat,0:LOMAX) -
!                                    CURRENTLY NOT USED OR NECESSARY!
!       BP(NRAD)                  -  r*u prime small component
!       VLM(NRAD,LMMX)            -  Potential from TAPE 19 (case.VNS) / VR(J,JATOM) holds VSP 

!       NEW FOR MGGA IMPLEMENTATION
!
!       AP1(NRAD,NSLMAX)           - r*u' large component
!       AP1LO(NRAD,nloat,0:LOMAX)  - r*u' large component for LO
!       APE(NRAD)                  - energy derivative of AP
!       APE1(NRAD,NSLMAX)          - energy derivative of AP1
!       BP1(NRAD,NSLMAX)           - r*u' small component
!       BP1LO(NRAD,nloat,0:LOMAX)  - r*u' small component for LO
!       BPE(NRAD)                  -
!       BPE1(NRAD,NSLMAX)          -
!       VTAU(NRAD,LMMX)            - Kinetic energy-dependent part of potential (functional derivative
!                                    of energy to non-interacting kinetic energy)
        
        DOUBLE PRECISION, allocatable ::   A(:), A1(:,:)
        DOUBLE PRECISION, allocatable ::   A1LO(:,:,:)
        DOUBLE PRECISION, allocatable ::   AE(:), AE1(:,:)
!       DOUBLE PRECISION, allocatable ::   AE1LO(:,:,:) NOT USED
        DOUBLE PRECISION, allocatable ::   AP(:), B(:), B1(:,:)
        DOUBLE PRECISION, allocatable ::   B1LO(:,:,:)
        DOUBLE PRECISION, allocatable ::   BE(:), BE1(:,:)
!       DOUBLE PRECISION, allocatable ::   BE1LO(:,:,:) NOT USED
        DOUBLE PRECISION, allocatable ::   BP(:), VLM(:,:)

!       New for MGGA implementation
        DOUBLE PRECISION, dimension(:),    allocatable ::   APE, BPE
        DOUBLE PRECISION, dimension(:,:),  allocatable ::   AP1, BP1, APE1, BPE1
        DOUBLE PRECISION, dimension(:,:,:),allocatable ::   AP1LO, BP1LO
        DOUBLE PRECISION, dimension(:,:),  allocatable ::   VTAU
      contains
        SUBROUTINE init_radial(NRAD, LOMAX, NSLMAX, NLOAT, LMMX)
          IMPLICIT NONE
          integer NRAD, NSLMAX, LOMAX, NLOAT, LMMX
          allocate( A(NRAD) )
          allocate( AP(NRAD) )
          allocate( A1(NRAD,NSLMAX) )
          allocate( A1LO(NRAD,nloat,0:LOMAX) )
          allocate( AE(NRAD) )
          allocate( AE1(NRAD,NSLMAX) )
          allocate( B(NRAD) )
          allocate( BP(NRAD) )
          allocate( B1(NRAD,NSLMAX) )
          allocate( B1LO(NRAD,nloat,0:LOMAX) )
          allocate( BE(NRAD) )
          allocate( BE1(NRAD,NSLMAX) )
          allocate( VLM(NRAD,0:LMMX) ) !expanded index to 0 for mgga_central_correction

          A=0.0D0
          AP=0.0D0
          A1=0.0D0
          A1LO=0.0D0
          AE=0.0D0
          AE1=0.0D0
          B=0.0D0
          BP=0.0D0
          B1=0.0D0
          B1LO=0.0D0
          BE=0.0D0
          BE1=0.0D0
          VLM=0.0D0
        END SUBROUTINE init_radial

        SUBROUTINE init_radial_mgga(NRAD, LOMAX, NSLMAX, NLOAT, LMMX)
          implicit none
          integer NRAD, LOMAX, NSLMAX, NLOAT, LMMX
          allocate(APE(NRAD), BPE(NRAD))
          allocate(AP1(NRAD, NSLMAX), BP1(NRAD, NSLMAX), &
               APE1(NRAD, NSLMAX), BPE1(NRAD, NSLMAX))
          allocate(AP1LO(NRAD, NLOAT, 0:LOMAX), BP1LO(NRAD, NLOAT, 0:LOMAX))
          allocate(VTAU(NRAD, 0:LMMX))

          APE = 0.0D0
          BPE = 0.0D0
          AP1 = 0.0D0
          BP1 = 0.0D0
          APE1 = 0.0D0
          BPE1 = 0.0D0
          AP1LO = 0.0D0
          BP1LO = 0.0D0
          VTAU = 0.0D0
        END SUBROUTINE init_radial_mgga
        
        SUBROUTINE end_radial
          deallocate(A, A1, A1LO, AE, AE1, AP, B, B1, B1LO, BE, BE1, BP, VLM)
          ! If VTAU is allocated, so are all other arrays of init_radial_mgga
          if(allocated(VTAU)) deallocate(APE, BPE, AP1, BP1, APE1, BPE1, &
               AP1LO, BP1LO, VTAU)
        END SUBROUTINE end_radial
      end module radial

!--------------------------------------------------------------

      module orb
!   NMOD=1 LDA+U, NMOD=2 OP, NMOD=3 Bext calculation
!   NSP=1 spin up, NSP=-1 spin down
!   NATORB number of atoms for which orbital potential is included
!   Iatom are the indexes of atom types selected for orb.pot.
!   NLORB number of L's of given atom for which OP is considered
!   LORB L's of given atom
!   Bext is energy of ext. magn field muB*Bext in Hartree
!   VORB is orbital potential matrix created by ORB package
!       COMPLEX*16 VORB(NATO,0:3,-3:3,-3:3)
!       DOUBLE PRECISION BEXT
!       INTEGER NMOD, NSP, NATORB, IATOM(NATO), NLORB(NATO), LORB(3,NATO)
!       COMMON /ORB/  Bext,VORB, NMOD,NSP,Natorb,IATOM,NLORB,LORB
!       SAVE   /ORB/
        DOUBLE PRECISION            :: BEXT
        DOUBLE COMPLEX, allocatable :: VORB(:,:,:,:)
        INTEGER                     :: NMOD, NSP, NATORB
        INTEGER, allocatable        :: IATOM(:), NLORB(:), LORB(:,:)
      contains
        SUBROUTINE init_orb(NAT)
          IMPLICIT NONE
          integer NAT
          allocate( VORB(NAT,0:3,-3:3,-3:3) )
          allocate( IATOM(NAT) )
          allocate( NLORB(NAT) )
          allocate( LORB(3,NAT) )
          VORB=(0.0D0,0.0D0)
          IATOM=0
          NLORB=0
          LORB=-1
        END SUBROUTINE init_orb
      end module orb

!--------------------------------------------------------------

      module lolog
!        LOOR(l,i) - selects which local orbitals should be used for
!                    atom i
!        NLO       - total number of used local orbitals
!
!      INTEGER            NLO
!      LOGICAL            LOOR(0:LOMAX,NATO),lapw(0:lmax-1,nato)
!      COMMON  /LOLOG/    NLO, LOOR, lapw, ilo(0:lomax,nato)
!      SAVE    /LOLOG/
        INTEGER              :: NLO
        LOGICAL, allocatable :: LOOR(:,:), LAPW(:,:)
        INTEGER, allocatable :: ILO(:,:)
      contains
        SUBROUTINE init_lolog(LOMAX, NATO, LMAX)
          IMPLICIT NONE
          integer LOMAX, NATO, LMAX
          allocate( LOOR(0:LOMAX, NATO) )
          allocate( LAPW(0:LMAX-1, NATO) )
          allocate( ILO(0:LOMAX, NATO) )
          LOOR=.true.
          LAPW=.true.
          ILO=0
        END SUBROUTINE init_lolog
      end module lolog

!--------------------------------------------------------------

      module atspdt
!        DP(l,i)      - derivative of the radial-function ul(r,E) with
!                       respect to radius r
!                       (evaluated at r = muffin tin radius of atom i)
!        DPE(l,i)     - derivative of the energy derivative ul_dot(r,E)
!                       of the radial-function ul(r,E) with respect to
!                       radius r
!                       (evaluated at r = muffin tin radius of atom i)
!        GFAC(:,i)    - gaunt-factors for atom i
!        INS          - flag if nonspherical correction (warpin) should
!                       take place
!                       if (INS .EQ. 0) then do warpin otherwise don't
!        LM(:,j,i)    - L,M combination for lattice harmonics j for
!                       atom i
!                       LM(1,j,i) ... L for j-th gaunt factor
!                       LM(2,j,i) ... M for j-th gaunt factor
!        LMMAX(i)     - number of L,M combinations for the
!                       nonspherical contribution of atom i
!        LQIND(i)     - number of gaunt factors for atom i
!        LQNS(:,j,i)  - l,m,LM,l',m' combination corresponding to
!                       the gaunt-factor GFAC(j,i) for atom i
!                       LQNS(1,1:LQIND(i),i) ... l'+1 for atom i
!                       LQNS(2,1:LQIND(i),i) ... L
!                       LQNS(3,1:LQIND(i),i) ... l+1 for atom i
!                       LQNS(4,1:LQIND(i),i) ... index of LM for atom i
!                       LQNS(5,1:LQIND(i),i) ... m' for atom i
!                       LQNS(6,1:LQIND(i),i) ... m for atom i
!        E(l,i)       - expansion energy El for atom i
!        P(l,i)       - radial-function ul(r,E)
!                       (evaluated at r = muffin tin radius of atom i)
!        PE(l,i)      - derivative of the radial-function ul(r,E) with
!                       respect to energy E
!                       (evaluated at r = muffin tin radius of atom i)
!        PEI(l,i)     - norm of ul_dot(r,E) integrated over the muffin
!                       tin sphere
!        VNS1(l',L,l,i) - non-spherical matrix elements for atom i
!                         integral(u_{l'}*V_{L}*u_{l} dr)
!        VNS2(l',L,l,i) - non-spherical matrix elements for atom i
!                         integral(u_{l'}_dot*V_{L}*u_{l}_dot dr)
!        VNS3(l',L,l,i) - non-spherical matrix elements for atom i
!                         integral(u_{l'}*V_{L}*u_{l}_dot dr)
!                         VNS3(l,L,l') ...
!                            ... integral(u_{l'}_dot*V_{L}*u_{l} dr)

        INTEGER                       :: INS, NL, NGAU
        INTEGER, allocatable          :: LMMAX(:), LQIND(:), LM(:,:,:)
        DOUBLE PRECISION, allocatable :: DP(:,:), DPE(:,:), E(:,:), P(:,:), PE(:,:), PEI(:,:)
        DOUBLE PRECISION, allocatable :: VNS1(:,:,:,:), VNS2(:,:,:,:), VNS3(:,:,:,:), &
                                         WNS1(:,:,:,:), WNS2(:,:,:,:), WNS3(:,:,:,:), &
                                         WNS4(:,:,:,:)
        DOUBLE COMPLEX,   pointer     :: GFAC(:,:)
        COMPLEX*16,       pointer     :: MGGA_FAC(:,:,:)
        INTEGER,          pointer     :: LQNS(:,:,:)
        ! Signed magnetic quantum number M of the potential harmonic for
        ! every LQNS record.  LQNS(4) stores only the index of WIEN2k's
        ! lattice-harmonic table, whose LM(2,...) entry is nonnegative.
        ! ATPAR also generates the conjugate -M channel, so reconstructing M
        ! from that table later loses its sign.  Keep it explicitly for the
        ! force/stress auxiliary file without changing the long-established
        ! six-field LQNS interface used by HNS/LMSORT.
        INTEGER,          pointer     :: LQMM(:,:)

!       Matrix elements required for mgga calculation        
        DOUBLE PRECISION, allocatable :: INT_U___VTAU_U__(:,:,:,:)  ! like VNS1
        DOUBLE PRECISION, allocatable :: INT_U___VTAU_UD_(:,:,:,:)  ! like VNS3
        DOUBLE PRECISION, allocatable :: INT_UD__VTAU_UD_(:,:,:,:)  ! like VNS2
                                                                   
        DOUBLE PRECISION, allocatable :: INT_U___VTAU_UP_(:,:,:,:)  ! like VNS1
        DOUBLE PRECISION, allocatable :: INT_U___VTAU_UPD(:,:,:,:)  ! like VNS3
        DOUBLE PRECISION, allocatable :: INT_UP__VTAU_UD_(:,:,:,:)  ! like VNS3
        DOUBLE PRECISION, allocatable :: INT_UD__VTAU_UPD(:,:,:,:)  ! like VNS2
                                                                   
        DOUBLE PRECISION, allocatable :: INT_UP__VTAU_UP_(:,:,:,:)  ! like VNS1
        DOUBLE PRECISION, allocatable :: INT_UP__VTAU_UPD(:,:,:,:)  ! like VNS3
        DOUBLE PRECISION, allocatable :: INT_UPD_VTAU_UPD(:,:,:,:)  ! like VNS2

      contains
        SUBROUTINE init_atspdt(NATO, LMMX, LMAX, NSLMAX)
          IMPLICIT NONE
          integer NATO, LMMX, LMAX, NSLMAX
          allocate( LMMAX(NATO) )
          allocate( LQIND(NATO) )
          allocate( LM(2,0:LMMX,NATO) )
          allocate( LQNS(6,NGAU,NATO) )
          allocate( LQMM(NGAU,NATO) )
          allocate( DP(LMAX,NATO) )
          allocate( DPE(LMAX,NATO) )
          allocate( E(LMAX,NATO) )
          allocate( P(LMAX,NATO) )
          allocate( PE(LMAX,NATO) )
          allocate( PEI(LMAX,NATO) )
          allocate( VNS1(NSLMAX,0:LMMX,NSLMAX,NATO) )
          allocate( VNS2(NSLMAX,0:LMMX,NSLMAX,NATO) )
          allocate( VNS3(NSLMAX,0:LMMX,NSLMAX,NATO) )
          ! Stress Eq. (G2): W contains r*u'_bra rather than u_bra.  It is
          ! therefore not symmetric under interchange of bra and ket.  Four
          ! APW channels are required: uu, udot-udot, u-udot, and udot-u.
          ! Keeping them next to VNS guarantees identical potential and
          ! radial-mesh conventions.
          allocate( WNS1(NSLMAX,0:LMMX,NSLMAX,NATO) )
          allocate( WNS2(NSLMAX,0:LMMX,NSLMAX,NATO) )
          allocate( WNS3(NSLMAX,0:LMMX,NSLMAX,NATO) )
          allocate( WNS4(NSLMAX,0:LMMX,NSLMAX,NATO) )
          allocate( GFAC(NGAU,NATO) )
          allocate(MGGA_FAC(1,1,1))
          LMMAX=0
          LQIND=0
          LM=0
          LQNS=0
          LQMM=0
          DP=0.0D0
          DPE=0.0D0
          E=0.0D0
          P=0.0D0
          PE=0.0D0
          PEI=0.0D0
          VNS1=0.0D0
          VNS2=0.0D0
          VNS3=0.0D0
          WNS1=0.0D0
          WNS2=0.0D0
          WNS3=0.0D0
          WNS4=0.0D0
          GFAC=(0.0D0,0.0D0)
        END SUBROUTINE init_atspdt

        SUBROUTINE init_atspdt_mgga(NATO, LMMX, LMAX, NSLMAX)
          implicit none
          integer NATO, LMMX, LMAX, NSLMAX
          allocate(INT_U___VTAU_U__(NSLMAX,0:LMMX,NSLMAX,NATO),&
               INT_U___VTAU_UD_(NSLMAX,0:LMMX,NSLMAX,NATO),&
               INT_UD__VTAU_UD_(NSLMAX,0:LMMX,NSLMAX,NATO),&
               INT_U___VTAU_UP_(NSLMAX,0:LMMX,NSLMAX,NATO),&
               INT_U___VTAU_UPD(NSLMAX,0:LMMX,NSLMAX,NATO),&
               INT_UP__VTAU_UD_(NSLMAX,0:LMMX,NSLMAX,NATO),&
               INT_UD__VTAU_UPD(NSLMAX,0:LMMX,NSLMAX,NATO),&
               INT_UP__VTAU_UP_(NSLMAX,0:LMMX,NSLMAX,NATO),&
               INT_UP__VTAU_UPD(NSLMAX,0:LMMX,NSLMAX,NATO),&
               INT_UPD_VTAU_UPD(NSLMAX,0:LMMX,NSLMAX,NATO))
          deallocate(MGGA_FAC)
          allocate(MGGA_FAC(4,NGAU,NATO))
          !LM = 0
          INT_U___VTAU_U__ = 0.0D0
          INT_U___VTAU_UD_ = 0.0D0
          INT_UD__VTAU_UD_ = 0.0D0
          INT_U___VTAU_UP_ = 0.0D0
          INT_U___VTAU_UPD = 0.0D0
          INT_UP__VTAU_UD_ = 0.0D0
          INT_UD__VTAU_UPD = 0.0D0
          INT_UP__VTAU_UP_ = 0.0D0
          INT_UP__VTAU_UPD = 0.0D0
          INT_UPD_VTAU_UPD = 0.0D0
          MGGA_FAC=(0.0Q0,0.0Q0)
        end SUBROUTINE init_atspdt_mgga

        function atom_gaunt(L0,LLL,LP,M,MM,MP,IATNR)
          implicit none
          complex*16          :: atom_gaunt
          integer, intent(in) :: L0,LLL,LP,M,MM,MP,IATNR

          ! check again the M selection rule
          if(M .ne. MM + MP) then
             atom_gaunt = 0.d0
             return
          end if
          
          IF (IATNR .GT. 0) then
             atom_gaunt = cubic_atom_gaunt(L0,LLL,LP,M,MM,MP)
          else
             atom_gaunt = noncubic_atom_gaunt(L0,LLL,LP,M,MM,MP)
          end IF
        end function atom_gaunt

        function clebsch_gordan_l_one(l,m)
!         Returns clebsch-gordan coefficients of l,m from l+-1 (first index),
!           m+mu (mu={-1,0,1}) (second index) coupled to mu
          implicit none
          real*8        :: l,m
          real*8, dimension(-1:1,-1:1) :: clebsch_gordan_l_one
          clebsch_gordan_l_one = 0.0Q0
              
          clebsch_gordan_l_one(-1,+1) = sqrt( (l+m-1.q0)*(l+m)/((2.q0*l-1.q0)*2.q0*l) )    
          clebsch_gordan_l_one(-1, 0) = sqrt( (l+m)*(l-m)/((2.q0*l-1.q0)*l) )    
          clebsch_gordan_l_one(-1,-1) = sqrt( (l-m-1.q0)*(l-m)/((2.q0*l-1.q0)*2.q0*l) )

          clebsch_gordan_l_one(+1,+1) = sqrt( (l-m+1.q0)*(l-m+2.q0)/((2.q0*l+2.q0)*(2.q0*l+3.q0)) )
          clebsch_gordan_l_one(+1, 0) =-sqrt( (l+m+1.q0)*(l-m+1.q0)/((l+1.q0)*(2.q0*l+3.q0)) )
          clebsch_gordan_l_one(+1,-1) = sqrt( (l+m+1.q0)*(l+m+2.q0)/((2.q0*l+2.q0)*(2.q0*l+3.q0)) )

        end function clebsch_gordan_l_one
        
        function cubic_atom_gaunt(L0,LLL,LP,M,MM,MP)
          implicit none
          
          real*8   gaunt
          external gaunt

          complex*16            :: cubic_atom_gaunt
          integer, intent(in)   :: L0,LLL,LP,M,MM,MP
          integer               :: LL
          complex*16            :: gfac
          complex*16, parameter :: one_imag = (0.0d0,1.0d0)
          LL = ABS(LLL)
          
          gfac = gaunt(L0,LL,LP,M,MM,MP)
          if ((LL .eq. 3) .or. (LL .eq. 5) .or. (LL .eq. 7) .or. (LL .eq. 9)) then
             if (MM .lt. 0) then
                gfac = gfac*one_imag
             else
                gfac = -gfac*one_imag
             end if
          end if
          
          cubic_atom_gaunt = gfac
          
        END function cubic_atom_gaunt
        
        function noncubic_atom_gaunt(L0,LLL,LP,M,MM,MP)
          implicit none
          
          real*8   gaunt
          external gaunt

          complex*16            :: noncubic_atom_gaunt
          integer, intent(in)   :: L0,LLL,LP,M,MM,MP
          complex*16            :: gfac,gaunt_sign
          integer               :: LL
          
          gaunt_sign = (1.0d0,0.0d0)
          LL = ABS(LLL)
          
          gfac = gaunt(L0,LL,LP,M,MM,MP)
          if (MM .ne. 0) then
             if (LLL .lt. 0) then
                gaunt_sign = (0.0d0,1.0d0)
                if (mod(mm,2).eq.0 .and. mm .gt. 0) gaunt_sign = -(0.0d0,1.0d0)
             else
                if (mod(mm,2).eq.1 .and. mm .gt. 0) gaunt_sign = -(1.0d0,0.0d0)
             end if
             gfac = gaunt_sign*gfac/SQRT(2.0D+0)
          end if
        
          noncubic_atom_gaunt = gfac
        END function noncubic_atom_gaunt
        
        subroutine calculate_mgga_fac(logic_mgga_gaunt,iatnr,&
             gfac_slice,mgga_fac_slice,l0,lll,lp,m,mm,mp)
          implicit none
          logical, dimension(0:4), intent(in)      :: logic_mgga_gaunt
          integer, intent(in)                      :: iatnr,l0,lll,lp,m,mm,mp
          complex*16, intent(out)                  :: gfac_slice
!          complex*16, dimension(1:4), intent(out)  :: mgga_fac_slice_out
          complex*16, dimension(1:4)  :: mgga_fac_slice

!         local variables
!         Clebsch-Gordan coefficients C(l+-1,m-mu,1,mu; l,m) 
          real*8, dimension(-1:1,-1:1)   :: C0, Cp
          real*8                         :: prefactor, factor_tmp, prefactor_all
          complex*16                     :: gaunt_tmp
          integer                        :: mu
          real*8                         :: l0_d,lp_d,m0_d,mp_d

          l0_d = real(l0,kind=8)
          lp_d = real(lp,kind=8)
          m0_d = real(m,kind=8)
          mp_d = real(mp,kind=8)

          C0 = clebsch_gordan_l_one(l0_d,m0_d)
          Cp = clebsch_gordan_l_one(lp_d,mp_d)

! We assume RYDBERG UNITS
          prefactor_all = 0.5D0/(SQRT((2.0D0*l0_d+1.0D0)*(2.0D0*lp_d+1.0D0)))
! Should be zerod in doreallocate or init routine, but let's make sure
          gfac_slice = (0.0D0,0.0D0)
          mgga_fac_slice = (0.0D0,0.0D0)
          
           if(logic_mgga_gaunt(0)) then
              gfac_slice = atom_gaunt(l0,lll,lp,m,mm,mp,iatnr)
           end if
          
          if(logic_mgga_gaunt(1)) then
             gaunt_tmp = (0.0D0,0.0D0)
             factor_tmp = 0.0D0
             prefactor = SQRT(l0_d*lp_d)
             do mu = -1, 1
                if ( abs(m-mu) .gt. l0-1 .or. abs(mp-mu) .gt. lp-1) cycle
                gaunt_tmp = atom_gaunt(l0-1,lll,lp-1,m-mu,mm,mp-mu,iatnr)
                factor_tmp = prefactor_all*prefactor*C0(-1,+mu)*Cp(-1,+mu)
                mgga_fac_slice(1) = gaunt_tmp*factor_tmp &
                     + mgga_fac_slice(1)
             end do
          end if
          
          if(logic_mgga_gaunt(2)) then
             gaunt_tmp = (0.0D0,0.0D0)
             factor_tmp = 0.0D0
             prefactor = -1.D0*SQRT(l0_d*(lp_d+1.0D0))
             do mu = -1, 1
                if ( abs(m-mu) .gt. l0-1 ) cycle
                gaunt_tmp = atom_gaunt(l0-1,lll,lp+1,m-mu,mm,mp-mu,iatnr)
                factor_tmp = prefactor_all*prefactor*C0(-1,+mu)*Cp(+1,+mu)
                mgga_fac_slice(2) = gaunt_tmp*factor_tmp &
                     + mgga_fac_slice(2)
             end do
          end if
          
          if(logic_mgga_gaunt(3)) then
             gaunt_tmp = (0.0D0,0.0D0)
             factor_tmp = 0.0D0
             prefactor = -1.d0*SQRT((l0_d+1.0D0)*lp_d)
             do mu = -1, 1
                if ( abs(mp-mu) .gt. lp-1) cycle
                gaunt_tmp = atom_gaunt(l0+1,lll,lp-1,m-mu,mm,mp-mu,iatnr)
                factor_tmp = prefactor_all*prefactor*C0(+1,+mu)*Cp(-1,+mu)
                mgga_fac_slice(3) = gaunt_tmp*factor_tmp &
                     + mgga_fac_slice(3)
             end do
          end if
          
          if(logic_mgga_gaunt(4)) then
             gaunt_tmp = (0.0D0,0.0D0)
             factor_tmp = 0.0D0
             prefactor = SQRT((l0_d+1.0D0)*(lp_d+1.0D0))
             do mu = -1, 1
!               should be satisfied always
!               if ( abs(m-mu) .gt. l0+1 .or. abs(mp-mu) .gt. lp+1) cycle
                gaunt_tmp = atom_gaunt(l0+1,lll,lp+1,m-mu,mm,mp-mu,iatnr)
                factor_tmp = prefactor_all*prefactor*C0(+1,+mu)*Cp(+1,+mu)
                mgga_fac_slice(4) = gaunt_tmp*factor_tmp &
                     + mgga_fac_slice(4)
             end do
          end if

!          mgga_fac_slice_out = mgga_fac_slice
          
!         DEBUG SECTION
!         ENCODE ANGULAR NUMBERS IN THE ENTRIES TO CHECK THEM LATER
!         TODO: REMOVE THIS SECTION LATER
          ! mgga_fac_slice(1) = 1.0d7 + l0_d*1.0d5 + abs(lll)*1.0d4 + &
          !      lp_d*1.0d3 + abs(m)*1.0d2 + abs(mm)*1.0d1 + abs(mp)
          ! mgga_fac_slice(2) = 2.0d7 + l0_d*1.0d5 + abs(lll)*1.0d4 + &
          !      lp_d*1.0d3 + abs(m)*1.0d2 + abs(mm)*1.0d1 + abs(mp)
          ! mgga_fac_slice(3) = 3.0d7 + l0_d*1.0d5 + abs(lll)*1.0d4 + &
          !      lp_d*1.0d3 + abs(m)*1.0d2 + abs(mm)*1.0d1 + abs(mp)
          ! mgga_fac_slice(4) = 4.0d7 + l0_d*1.0d5 + abs(lll)*1.0d4 + &
          !      lp_d*1.0d3 + abs(m)*1.0d2 + abs(mm)*1.0d1 + abs(mp)



          

        end subroutine calculate_mgga_fac
      end module atspdt

!--------------------------------------------------------------

      module char
!        LATTIC  - lattice type
!                  'P' ...... primitive latice (cubic, tetragonal,
!                             orthorhombic, monoclinic, triclin)
!                  'FC' ..... face centered
!                  'BC' ..... body centered
!                  'HEX' .... hexagonal
!                  'CXY' .... c-base centered (only orthorombic)
!                  'CYZ' .... a-base centered (only orthorombic)
!                  'CXZ' .... b-base centered (only orthorombic)
!                  'R' ...... rhombohedral
!        NAME(i) - name of atom i in the unit-cell (e.g. 'Titanium')
!
!      CHARACTER*4        LATTIC
!      CHARACTER*10       NAME(NATO)
!      COMMON  /CHAR/     LATTIC, NAME
!      SAVE    /CHAR/
        CHARACTER*4 :: LATTIC
        CHARACTER*10, allocatable :: NAME(:)
      contains
        SUBROUTINE init_char(NAT)
          IMPLICIT NONE
          integer NAT
          allocate( NAME(NAT) )
        END SUBROUTINE init_char
      end module char

!--------------------------------------------------------------

      module comc
!        IPGR(j)  - point group of j-th k-point (only used for output
!                   (documentation) purposes)
!        KNAME(j) - name of j-th k-point (optional)
!
!      CHARACTER*3        IPGR(NKPT)
!      CHARACTER*10       KNAME(NKPT)
!      COMMON  /COMC/     KNAME, IPGR
!      SAVE    /COMC/
        CHARACTER*3, pointer :: IPGR(:)
        CHARACTER*10, pointer :: KNAME(:)
      contains
        SUBROUTINE init_comc(NKPT, I)
          use reallocate, only: doreallocate
          IMPLICIT NONE
          integer NKPT, I
          double precision :: dummy
          if(I.eq.0) THEN
            allocate( IPGR(NKPT) )
            allocate( KNAME(NKPT) )
          else
            call doreallocate(IPGR, NKPT, dummy)
            call doreallocate(KNAME, NKPT)
          end if
        END SUBROUTINE init_comc
      end module comc

!--------------------------------------------------------------

      module comi
!        LNSMAX - maximum considered azimuthal quantum number l of ul
!                 for non-spherical hamilton matrix contributions
!        NAT    - number of inequivalent atoms
!        NBELW  - number of Eigenvalues below given energy window
!        NE     - number of calculated Eigenvalues
!        NT     - maximum considered azimuthal quantum number + 1
!                 (l+1 of ul) for the spherical contribution
!        NVAA   - number of plane waves
!
!      INTEGER            NVAA, NE, NBELW, NAT, NT, LNSMAX
!      COMMON  /COMI/     NVAA, NE, NBELW, NAT, NT, LNSMAX
!      SAVE    /COMI/
        INTEGER :: NVAA, NE, NBELW, NAT, NT, LNSMAX,nemax_it
!       Explicit opt-in for selecting a fixed number of lowest bands.
!       This must remain separate from nemax_it so input files without the
!       FIXED_NBANDS keyword retain WIEN2k's historical behavior.
        LOGICAL :: fixed_nbands
        INTEGER :: NSLMAX
      end module comi

!--------------------------------------------------------------

      module coml
!     PRNTWF - if .TRUE. print wave functions to output-file
!     SPRSWF - if .TRUE. suppress wave function calculation
!     WFTAPE - if .TRUE. write wave functions to vector-file
!     REL    - if .TRUE. perform relativistic calculations
!     ITER   - if .TRUE. perform iterative diagonalization
!     mgga   - if .TRUE. perform meta-gga calculation
!     mgga_central_correction
!            - if .TRUE. perform mgga central correction
        
        
        
!     NOHNS  - if .TRUE. supress nonspherical matrix elements
        
        
!
!     LOGICAL     SPRSWF,  PRNTWF,  WFTAPE,  REL, ITER,NOHNS,force
!     COMMON  /COML/     SPRSWF,  WFTAPE,  PRNTWF,  REL, ITER,NOHNS
!     SAVE    /COML/
      
        LOGICAL :: SPRSWF,  PRNTWF,  WFTAPE,  REL, ITER,NOHNS
        logical :: mgga, mgga_central_correction
      end module coml

!--------------------------------------------------------------

      module comr
!        EL        - lower boundary of energy window for Eigenvalues
!        EU        - upper boundary of energy window for Eigenvalues
!        EIGVAL(j) - calculated Eigenvalues
!        RKM       - Rmt*Kmax determines matrix size (convergence)
!                    Rmt .... muffin tin radius
!                    Kmax ... plane wave cut-off
!        WEIGHT(j) - weight of the j-th k-vector (order of group of k)
!
!      DOUBLE PRECISION   RKM,  EL,  EU
!      DOUBLE PRECISION   WEIGHT(NKPT)
!      COMMON  /COMR/     RKM,  WEIGHT,  EL,  EU
!      SAVE    /COMR/
        DOUBLE PRECISION :: RKM, EL, EU,sproj_limit, eshift_iter
        DOUBLE PRECISION, pointer :: WEIGHT(:)
        logical readHinv, divide
      contains
        SUBROUTINE init_comr(NKPT, I)
          use reallocate, only: doreallocate
          IMPLICIT NONE
          integer NKPT, I
          if(I.eq.0) THEN
            allocate( WEIGHT(NKPT) )
            WEIGHT(1:NKPT)=0
          else
            call doreallocate(weight, NKPT)
          end if
        END SUBROUTINE init_comr
      end module comr

!--------------------------------------------------------------

      module gener
!
!        BR2(1:3,1:3) - reciprocal Bravais matrix
!
!      DOUBLE PRECISION   BR2(3,3)
!      COMMON  /GENER/    BR2
!      SAVE    /GENER/
        DOUBLE PRECISION ::  BR2(3,3)
      end module gener

!--------------------------------------------------------------

      module kpts
!        SX(j) - x-component of K-point j
!        SY(j) - y-component of K-point j
!        SZ(j) - z-component of K-point j
!
!      DOUBLE PRECISION   SX(NKPT),  SY(NKPT),  SZ(NKPT)
!      COMMON  /KPTS/     SX,  SY,  SZ
!      SAVE    /KPTS/
        DOUBLE PRECISION, pointer :: SX(:),  SY(:),  SZ(:)

      contains
        SUBROUTINE init_kpts(NKPT, I)
          use reallocate, only: doreallocate
          IMPLICIT NONE
          integer NKPT, I
          if(I.eq.0) THEN
            allocate( SX(NKPT) )
            allocate( SY(NKPT) )
            allocate( SZ(NKPT) )
            SX(1:NKPT)=0.0D0
            SY(1:NKPT)=0.0D0
            SZ(1:NKPT)=0.0D0
          else
            call doreallocate( SX, NKPT)
            call doreallocate( SY, NKPT)
            call doreallocate( SZ, NKPT)
          end if
        END SUBROUTINE init_kpts  
      end module kpts

!--------------------------------------------------------------

      module loint
!        VNS1LO(l',L,l) - local orbital non-spherical matrix elements
!                         integral(ulo_{l'}*V_{L}*u_{l} dr)
!        VNS2LO(l',L,l) - local orbital non-spherical matrix elements
!                         integral(ulo_{l'}*V_{L}*u_{l}_dot dr)
!        VNS3LO(l',L,l) - local orbital non-spherical matrix elements
!                         integral(ulo_{l'}*V_{L}*ulo_{l} dr)
!                         VNS3(l,L,l') ...
!                           ... integral(u_{l'}_dot*V_{L}*u_{l} dr)

!        INT_Ulo_VTAU_U(l',L,l)      - local orbital kinetic-energy dependent
!                                        matrix elements
!                        integral(ulo_{l'}*Vtau_{L}*u_{l}/r^2 dr)
!        INT_Ulo_VTAU_UD(l',L,l)     - local orbital kinetic-energy dependent
!                                        matrix elements
!                        integral(ulo_{l'}*Vtau_{L}*u_{l}_dot/r^2 dr)
!        INT_Ulo_VTAU_Ulo(l',L,l)    - local orbital kinetic-energy dependent
!                                        matrix elements
!                        integral(ulo_{l'}*Vtau_{L}*ulo_{l}/r^2 dr)
!        
!        INT_Ulo_VTAU_UP(l',L,l)      - local orbital kinetic-energy dependent
!                                        matrix elements
!                        integral(ulo_{l'}*Vtau_{L}*u'_{l}/r dr)
!        INT_Ulo_VTAU_UPD(l',L,l)     - local orbital kinetic-energy dependent
!                                        matrix elements
!                        integral(ulo_{l'}*Vtau_{L}*u'_{l}_dot/r dr)
!        INT_UloP_VTAU_U(l',L,l)     - local orbital kinetic-energy dependent
!                                        matrix elements
!                        integral(ulo'_{l'}*Vtau_{L}*u_{l}/r dr)
!        INT_UloP_VTAU_UD(l',L,l)    - local orbital kinetic-energy dependent
!                                        matrix elements
!                        integral(ulo'_{l'}*Vtau_{L}*u_{l}_dot/r dr)
!        
!        
!        INT_UloP_VTAU_UP(l',L,l)    - local orbital kinetic-energy dependent
!                                        matrix elements
!                        integral(ulo'_{l'}*Vtau_{L}*u'_{l} dr)
!        INT_UloP_VTAU_UPD(l',L,l)   - local orbital kinetic-energy dependent
!                                        matrix elements
!                        integral(ulo'_{l'}*Vtau_{L}*u'_{l}_dot dr)
!        INT_UloP_VTAU_UloP(l',L,l)  - local orbital kinetic-energy dependent
!                                        matrix elements
!                        integral(ulo'_{l'}_dot*Vtau_{L}*ulo'_{l} dr)
        
        
!      DOUBLE PRECISION   VNS1LO(LOMAX+1,LMMX,NSLMAX,nloat,NATO)
!      DOUBLE PRECISION   VNS2LO(LOMAX+1,LMMX,NSLMAX,nloat,NATO)
!      DOUBLE PRECISION   VNS3LO(LOMAX+1,LMMX,NSLMAX,nloat,nloat,NATO)
        
        ! LOMAX  = max l number for lo's
        ! LMMX   = max lm terms for non-spherical potential VNS
        ! NSLMAX = LNSMAX + 1
        ! LNSMAX - maximum considered azimuthal quantum number l of ul (L)APW
        ! nloat  = Number of (non-equivalent) atoms with lo
        ! NATO   = Number of non-equivalent atoms
        
        !      COMMON  /LOINT/    VNS1LO, VNS2LO, VNS3LO
!      SAVE    /LOINT/
        DOUBLE PRECISION, allocatable :: VNS1LO(:,:,:,:,:), VNS2LO(:,:,:,:,:), VNS3LO(:,:,:,:,:,:), &
                                         WNS1LO(:,:,:,:,:), WNS2LO(:,:,:,:,:), WNS3LO(:,:,:,:,:,:), &
                                         WNS4LO(:,:,:,:,:), WNS5LO(:,:,:,:,:)
        
        DOUBLE PRECISION, allocatable :: INT_Ulo__VTAU_U__(:,:,:,:,:)
        DOUBLE PRECISION, allocatable :: INT_Ulo__VTAU_UD_(:,:,:,:,:)
        DOUBLE PRECISION, allocatable :: INT_Ulo__VTAU_Ulo(:,:,:,:,:,:)

        DOUBLE PRECISION, allocatable :: INT_Ulo__VTAU_UP_(:,:,:,:,:)
        DOUBLE PRECISION, allocatable :: INT_Ulo__VTAU_UPD(:,:,:,:,:)
        DOUBLE PRECISION, allocatable :: INT_UloP_VTAU_U__(:,:,:,:,:)
        DOUBLE PRECISION, allocatable :: INT_UloP_VTAU_UD_(:,:,:,:,:)
        DOUBLE PRECISION, allocatable :: INT_UloP_VTAU_Ulo(:,:,:,:,:,:)

        DOUBLE PRECISION, allocatable :: INT_UloP_VTAU_UP_(:,:,:,:,:)
        DOUBLE PRECISION, allocatable :: INT_UloP_VTAU_UPD(:,:,:,:,:)
        DOUBLE PRECISION, allocatable :: INT_UloP_VTAU_UloP(:,:,:,:,:,:)

      contains
        SUBROUTINE init_loint(LOMAX, LMMX, NSLMAX, NATO,nloat)
          IMPLICIT NONE
          integer LOMAX, LMMX, NSLMAX, NATO,nloat
          allocate( VNS1LO(LOMAX+1, 0:LMMX, NSLMAX, nloat, NATO) )
          allocate( VNS2LO(LOMAX+1, 0:LMMX, NSLMAX, nloat, NATO) )
          allocate( VNS3LO(LOMAX+1, 0:LMMX, NSLMAX, nloat, nloat,NATO) )
          allocate( WNS1LO(LOMAX+1, 0:LMMX, NSLMAX, nloat, NATO) )
          allocate( WNS2LO(LOMAX+1, 0:LMMX, NSLMAX, nloat, NATO) )
          allocate( WNS3LO(LOMAX+1, 0:LMMX, NSLMAX, nloat, nloat,NATO) )
          ! WNS1LO/WNS2LO/WNS3LO have a differentiated LO on the bra.
          ! WNS4LO/WNS5LO are the independent differentiated APW-bra,
          ! LO-ket channels; they cannot be obtained by transposition.
          allocate( WNS4LO(NSLMAX, 0:LMMX, LOMAX+1, nloat, NATO) )
          allocate( WNS5LO(NSLMAX, 0:LMMX, LOMAX+1, nloat, NATO) )
          VNS1LO=0.0D0
          VNS2LO=0.0D0
          VNS3LO=0.0D0
          WNS1LO=0.0D0
          WNS2LO=0.0D0
          WNS3LO=0.0D0
          WNS4LO=0.0D0
          WNS5LO=0.0D0
!todo          ! write(6,'(a20,f12.1,a,a,6i6)') "Allocating VNS3LO ",&
          !      & (LOMAX+1)*LMMX*NSLMAX*nloat*nloat*NATO,' MB ', &
          !      & 'dimensions', LOMAX+1, LMMX, NSLMAX, nloat, nloat,NATO
        END SUBROUTINE init_loint

        subroutine init_loint_mgga(lomax, lmmx, nslmax, nato, nloat)
          implicit none
          integer LOMAX, LMMX, NSLMAX, NATO,nloat
          allocate(INT_Ulo__VTAU_U__(LOMAX+1, 0:LMMX, NSLMAX, nloat, NATO))
          allocate(INT_Ulo__VTAU_UD_(LOMAX+1, 0:LMMX, NSLMAX, nloat, NATO))
          allocate(INT_Ulo__VTAU_UP_(LOMAX+1, 0:LMMX, NSLMAX, nloat, NATO))
          allocate(INT_Ulo__VTAU_UPD(LOMAX+1, 0:LMMX, NSLMAX, nloat, NATO))
          allocate(INT_UloP_VTAU_U__(LOMAX+1, 0:LMMX, NSLMAX, nloat, NATO))
          allocate(INT_UloP_VTAU_UD_(LOMAX+1, 0:LMMX, NSLMAX, nloat, NATO))
          allocate(INT_UloP_VTAU_UP_(LOMAX+1, 0:LMMX, NSLMAX, nloat, NATO))
          allocate(INT_UloP_VTAU_UPD(LOMAX+1, 0:LMMX, NSLMAX, nloat, NATO))
          
          allocate(INT_Ulo__VTAU_Ulo(LOMAX+1, 0:LMMX, NSLMAX, nloat, nloat, NATO))
          allocate(INT_UloP_VTAU_Ulo(LOMAX+1, 0:LMMX, NSLMAX, nloat, nloat, NATO))
          allocate(INT_UloP_VTAU_UloP(LOMAX+1, 0:LMMX, NSLMAX, nloat, nloat, NATO))

          INT_Ulo__VTAU_U__  = 0.0D0
          INT_Ulo__VTAU_UD_  = 0.0D0
          INT_Ulo__VTAU_UP_  = 0.0D0
          INT_Ulo__VTAU_UPD  = 0.0D0
          INT_UloP_VTAU_U__  = 0.0D0
          INT_UloP_VTAU_UD_  = 0.0D0
          INT_UloP_VTAU_UP_  = 0.0D0
          INT_UloP_VTAU_UPD  = 0.0D0

          INT_Ulo__VTAU_Ulo  = 0.0D0
          INT_UloP_VTAU_Ulo  = 0.0D0
          INT_UloP_VTAU_UloP = 0.0D0
          
        end subroutine init_loint_mgga
      end module loint

!--------------------------------------------------------------

      module lstapw
!        NV    - total number of considered reciprocal lattice vectors
!                (number of plane waves)
!        XK(j) - x-coordinate (cartesian) of j-th k-vector
!                (ordered by distance from origin)
!        YK(j) - y-coordinate (cartesian) of j-th k-vector
!                (ordered by distance from origin)
!        ZK(j) - z-coordinate (cartesian) of j-th k-vector
!                (ordered by distance from origin)
!        RKM       - Rmt*Kmax determines matrix size (input)
!
!      INTEGER            NV
!      DOUBLE PRECISION   RKMT
!      DOUBLE PRECISION   XK(HSROWS+1), YK(HSROWS+1), ZK(HSROWS+1)
!      COMMON  /LSTAPW/   RKMT, NV
!      SAVE    /LSTAPW/
        DOUBLE PRECISION :: RKMT
        INTEGER          :: NV
      end module lstapw

!--------------------------------------------------------------

      module onekpt
!        SX,SY,SZ - k-point in the irreducible Brillouin Zone
!
!      DOUBLE PRECISION   SX1, SY1, SZ1
!      COMMON  /ONEKPT/   SX1, SY1, SZ1
!      SAVE    /ONEKPT/
        DOUBLE PRECISION :: SX1, SY1, SZ1
      end module onekpt

!--------------------------------------------------------------

      module orth
!        ORTH - .TRUE. for orthogonal lattice
!               .FALSE. otherwise
!
!      LOGICAL            ORTHO
!      COMMON  /ORTH/     ORTHO
!      SAVE    /ORTH/
        LOGICAL          :: ORTHO
      end module orth

!--------------------------------------------------------------

      module out
!        IMAT(1:3,1:3,j) - matrix representation of (space group)
!                          symmetry operation j
!        IORD            - number of symmetry operations of space group
!        KKK(1:3,j)      - star vectors of reciprocal lattice
!                          vector KZZ(:)
!        KZZ1(1:3)       - reciprocal lattice vector
!        NKK             - number of vectors in KKK (size of the star)
!        TAU(1:3,j)      - non-primitive translation vector for symmetry
!                          operation j
!        TAUP(j)         - tauphase factor of vector KKK(:,j)
!        WARP(i,j,k)     - interstitial warpin for K vector (i,j,k)
!
!      INTEGER            IORD, NKK
!      INTEGER            KKK(3,NSYM), IMAT(3,3,NSYM)
!      DOUBLE PRECISION   TAU(3,NSYM)
!      COMPLEX*16         TAUP(NSYM)
!      COMPLEX*16         WARP(-KMAX1:KMAX1,-KMAX2:KMAX2,-KMAX3:KMAX3)
!      COMMON  /OUT/      WARP, TAU, TAUP, IORD, IMAT, KZZ1, KKK, NKK
!      SAVE    /OUT/
        INTEGER                       :: IORD, NKK
        INTEGER                       :: KZZ1(3)
        INTEGER, allocatable          :: KKK(:,:), IMAT(:,:,:)
        DOUBLE PRECISION, allocatable :: TAU(:,:)
!_COMPLEX        DOUBLE COMPLEX, allocatable   :: TAUP(:)
!_REAL           DOUBLE PRECISION, allocatable :: TAUP(:)
!_COMPLEX        DOUBLE COMPLEX, pointer       :: WARP(:,:,:)
!_REAL           DOUBLE PRECISION, pointer     :: WARP(:,:,:)
!_COMPLEX        DOUBLE COMPLEX, allocatable   :: vtau_warp(:,:,:)
!_REAL           DOUBLE PRECISION, allocatable :: vtau_warp(:,:,:)
        INTEGER                       :: WARP1,WARP2,WARP3
        LOGICAL                       :: kkk_allocated
        LOGICAL                       :: warp_allocated
      contains
        SUBROUTINE init_out(NSYM, KMAX1, KMAX2, KMAX3)
          IMPLICIT NONE
          integer NSYM, KMAX1, KMAX2, KMAX3
          DOUBLE COMPLEX :: CZERO
          PARAMETER ( CZERO=(0.0D0,0.0D0) )
          if(KMAX1.eq.0) THEN
            allocate( KKK(3,NSYM) )
            allocate( IMAT(3,3,NSYM) )
            allocate( TAU(3,NSYM) )
            allocate( TAUP(NSYM) )
            kkk_allocated=.true.
            warp_allocated=.false.
            KKK(1:3,1:nsym)=0
            IMAT(1:3,1:3,1:nsym)=0
            TAU(1:3,1:nsym)=0.0D0
!_COMPLEX            TAUP(1:nsym)=CZERO
!_REAL               TAUP(1:nsym)=0.0D0
          else
            allocate( WARP(-KMAX1:KMAX1,-KMAX2:KMAX2,-KMAX3:KMAX3) )
!_COMPLEX            WARP(-KMAX1:KMAX1,-KMAX2:KMAX2,-KMAX3:KMAX3)=CZERO
!_REAL               WARP(-KMAX1:KMAX1,-KMAX2:KMAX2,-KMAX3:KMAX3)=0.0D0
            WARP1=KMAX1
            WARP2=KMAX2
            WARP3=KMAX3
            warp_allocated=.true.
          end if
        END SUBROUTINE init_out
        SUBROUTINE end_out
          if(kkk_allocated) then
            deallocate( KKK, IMAT, TAU, TAUP )
            kkk_allocated=.false.
          end if
          if(warp_allocated) then
            deallocate( WARP )
            warp_allocated=.false.
          end if
        END SUBROUTINE end_out
      end module out

!--------------------------------------------------------------

      module potnlc
!        DH(i)   - step size of the logarithmic radial mesh for atom i
!        JRJ(i)  - number of radial (logarithmic) mesh points for atom i
!        R0(i)   - first radial mesh point for atom i
!        VR(j,i) - spherical part (l=0, m=0) of the total potential r*V
!                  at mesh point j for atom i
!
!      INTEGER            JRJ(NATO)
!      DOUBLE PRECISION   DH(NATO), R0(NATO), VR(NRAD,NATO)
!      COMMON  /POTNLC/   VR,RO,DH,JRJ
!      SAVE    /POTNLC/
        INTEGER, allocatable :: JRJ(:)
        DOUBLE PRECISION, allocatable :: DH(:), R0(:), VR(:,:)
      contains
        SUBROUTINE init_potnlc(NATO, NRAD)
          IMPLICIT NONE
          integer NATO, NRAD
          allocate( JRJ(NATO) )
          allocate( DH(NATO) )
          allocate( R0(NATO) )
          allocate( VR(NRAD,NATO) )
          JRJ(1:NATO)=0
          DH(1:NATO)=0.0D0
          R0(1:NATO)=0.0D0
          VR(1:NRAD,1:NATO)=0.0D0
        END SUBROUTINE init_potnlc
      end module potnlc

!--------------------------------------------------------------

      module rotmat
!
!        ROTIJ(1:3,1:3,i)  - rotate atom i (including equivalent ones)
!                            according to the symmetry. For atoms with
!                            only one equiv. atom this is always the
!                            identity matrix.
!        ROTLOC(1:3,1:3,i) - local rotation matrix for atom i (from
!                            global to local coordinate system)
!
!      DOUBLE PRECISION   ROTIJ(3,3,NDIF),ROTLOC(3,3,NATO)
!      COMMON  /ROTMAT/   ROTLOC, ROTIJ
!      SAVE    /ROTMAT/
        DOUBLE PRECISION, allocatable :: ROTLOC(:,:,:)
        DOUBLE PRECISION, pointer :: ROTIJ(:,:,:),tauij(:,:)
      contains
        SUBROUTINE init_rotmat(NDIF, NATO, I)
          use reallocate, only: doreallocate
          IMPLICIT NONE
          integer NDIF, NATO, I
          if(I.eq.0) THEN
            allocate( ROTLOC(3,3,NATO) )
            allocate( ROTIJ(3,3,NDIF) )
            allocate( tauIJ(3,NDIF) )
            ROTLOC(1:3,1:3,1:NATO)=0.0D0
            ROTIJ(1:3,1:3,1:NDIF)=0.0D0
            tauIJ(1:3,1:NDIF)=0.0D0
          else
            call doreallocate(ROTIJ, 3,3,NDIF)
          end if
        END SUBROUTINE init_rotmat
      end module rotmat

!--------------------------------------------------------------

      module cut1
!      logical            cut
!      common /cut1/      cut
        LOGICAL :: CUT
      end module cut1

!--------------------------------------------------------------

      module totpot
!        JA(j)   - 1. component of the j-th representative star vector
!        JB(j)   - 2. component of the j-th representative star vector
!        JC(j)   - 3. component of the j-th representative star vector
!        POTK(j) - interstitial total potential for the j-th star
!
!      INTEGER            JA(NWAV),  JB(NWAV),  JC(NWAV)
!      COMPLEX*16         POTK(NWAV)
!      COMMON  /TOTPOT/   POTK,  JA,  JB,  JC
!      SAVE    /TOTPOT/
        INTEGER, allocatable        :: JA(:), JB(:), JC(:)
!_COMPLEX        DOUBLE COMPLEX, allocatable :: POTK(:), vtau_star(:)
!_REAL           DOUBLE PRECISION, allocatable :: POTK(:), vtau_star(:)
      contains
        SUBROUTINE init_totpot(NWAV)
          IMPLICIT NONE
          integer NWAV
          allocate( JA(NWAV) )
          allocate( JB(NWAV) )
          allocate( JC(NWAV) )
          allocate( POTK(NWAV) )
          JA(1:NWAV)=0
          JB(1:NWAV)=0
          JC(1:NWAV)=0
!_COMPLEX          POTK(1:NWAV)=(0.0D0,0.0D0)
!_REAL             POTK(1:NWAV)=0.0D0
        END SUBROUTINE init_totpot

        subroutine init_vtau_star(nwav)
          implicit none
          integer :: nwav

          allocate(vtau_star(nwav))
!_COMPLEX vtau_star = (0.0D0,0.0D0)
!_REAL    vtau_star = 0.0D0
          
        end subroutine init_vtau_star
      end module totpot

!--------------------------------------------------------------

      module struk
!        NDF         - number of atoms
!        ALAT(1:3)   - lattice constants (a,b,c)
!        ALPHA(1:3)  - angles of the unit-cell's unit-vectors
!        IATNR(i)    - atom index of (inequivalent) atom i
!                      also indicates cubic and non-cubic symmetries
!                      IATNR(i) .GT. 0 ... cubic symmetry
!                      IATNR(i) .LT. 0 ... non-cubic symmetry
!        MULT(i)     - number of equivalent atoms for inequiv. atom i
!        PIA(1:3)    - reciprocal lattice constants (2*pi/a, 2*pi/b,
!                      2*pi/c)
!        POS(1:3,i)  - position vector of atom i (including equivalent
!                      atoms) in the unit-cell
!        RMT(i)      - muffin tin radius of atom i
!        V(i)        - relative muffin tin spherevolume for atom i
!        VI          - inverse volume of the direct lattice unit-cell
!
!      INTEGER            IATNR(NATO), MULT(NATO)
!      DOUBLE PRECISION   VI
!      DOUBLE PRECISION   ALAT(3), ALPHA(3), PIA(3), POS(3,NDIF)
!      DOUBLE PRECISION   RMT(NATO), V(NATO)
!      COMMON  /STRUK/    POS, ALAT, ALPHA, RMT, V, PIA, VI, IATNR, MULT
!      SAVE    /STRUK/
        INTEGER              :: NDF,multmax
        INTEGER, allocatable :: IATNR(:), MULT(:)
        DOUBLE PRECISION     :: VI
        DOUBLE PRECISION     :: ALAT(3), ALPHA(3), PIA(3)
        DOUBLE PRECISION, allocatable :: RMT(:), V(:)
        DOUBLE PRECISION, pointer :: POS(:,:)
      contains
        SUBROUTINE init_struk(NATO, NDIF, I)
          use reallocate, only: doreallocate
          IMPLICIT NONE
          integer NATO, NDIF, I
          if(I.eq.0) THEN
            allocate( IATNR(NATO) )
            allocate( MULT(NATO) )
            allocate( RMT(NATO) )
            allocate( V(NATO) )
            allocate( POS(3, NDIF) )
            IATNR(1:NATO)=0
            MULT(1:NATO)=0
            RMT(1:NATO)=0.0D0
            V(1:NATO)=0.0D0
            POS(1:3,1:NDIF)=0.0D0
          else
            call doreallocate(POS, 3, NDIF)
          end if
          NDF = NDIF
        END SUBROUTINE init_struk
      end module struk

!--------------------------------------------------------------
      module kpthold
      integer number_kpt
      end module kpthold

!--------------------------------------------------------------

      module nmr

        logical  writeham,  dokorig
        real*8    korig(3)

      end module nmr

      module stdmath
#if defined (HAVE_LIBMVEC) && !defined(INTEL_VML)
      interface
        pure function dcos(x) bind(C, NAME="cos")
          use,intrinsic :: iso_c_binding
          real(c_double), value, intent(in) :: x
          real(c_double) :: dcos
!$omp declare simd (dcos)
        end function

        pure function dsin(x) bind(C, NAME="sin")
          use,intrinsic :: iso_c_binding
          real(c_double), value, intent(in) :: x
          real(c_double) :: dsin
!$omp declare simd (dsin)
        end function
      end interface
#endif
      end module stdmath
