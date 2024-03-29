! @Copyright 2007 Kristjan Haule
module LinLogMesh
  !##########################
  !# Rreal axis mesh        #
  !##########################
  IMPLICIT NONE
  REAL*8  :: delta, ommax
  INTEGER :: Nd, Nitt, nom, zero_ind
  REAL*8, allocatable :: om(:)             ! om(nom)
  INTEGER, allocatable:: idx(:,:), nidx(:) ! idx(Nitt,nom), nidx(Nitt)
CONTAINS
  SUBROUTINE LinLogMesh__Init__(filename)
    IMPLICIT NONE
    CHARACTER*200, intent(in) :: filename
    !
    INTEGER, PARAMETER :: fh_mesh=91
    INTEGER :: ios
    open(fh_mesh, FILE=filename, STATUS='old', IOSTAT=ios)
    READ(fh_mesh,*) delta, ommax, Nd
    close(fh_mesh)
    if (ios.NE.0) then
       print *, 'ERROR reading mesh ', filename
    endif
    Nitt = int(log(1+ommax/(delta*Nd))/log(2.)+0.5)
    nom = Nitt*Nd*2+1
    ALLOCATE( om(nom), idx(Nitt, nom), nidx(Nitt) )
  END SUBROUTINE LinLogMesh__Init__
  
  SUBROUTINE LinLogMesh__Destruct__()
    DEALLOCATE( om, idx, nidx )
  END SUBROUTINE LinLogMesh__Destruct__

  SUBROUTINE LinLogMesh_Generate()
    ! Creates logarithmic mesh of linear meshes.
    ! In the first level we have Nd*2 points in linear mesh with spacing delta.
    ! In the second level we add Nd*2 points with spacing 2*delta.
    ! In each subsequent level we add Nd*2 points with doubled spacing.
    IMPLICIT NONE
    ! locals
    REAL*8  :: x, y, deltar
    INTEGER :: j, itt, ii, nsubidx
    INTEGER :: subidx(nom)
    deltar = delta
    x=0.0
    
    zero_ind = Nitt*Nd+1
    om(zero_ind)=0.d0
    ii=1
    do itt=1,Nitt
       y = x+deltar
       do j=1,Nd
          om(zero_ind+ii) = y
          om(zero_ind-ii) = -y
          ii = ii + 1
          y = y+deltar
       enddo
       x = x + deltar*Nd
       deltar = deltar*2
    enddo
    
    ! Create an index array such that mesh[idx[i][j]] constitutes a linear mesh with regular spacings.
    ! In the first level  mesh[idx[0][:]] will give linear mesh of points in the interval [-delta*Nd,delta*Nd]
    ! In the second level mesh[idx[1][:]] will give linear mesh of points in the interval [-(1+2)*delta*Nd, (1+2)*delta*Nd] with spacing 2*delta
    ! In the k-th level   mesh[idx[k][:]] will give linear mesh of points in the interval [-(2**k-1)*delta*Nd, (2**k-1)*delta*Nd] with spacing 2**(k-1)*delta
    nsubidx=1
    do itt=1,Nitt
       ii=1
       do j=zero_ind-itt*Nd, zero_ind-(itt-1)*Nd-1
          !print *, 'j0=', j, om(j)
          subidx(ii)=j
          ii=ii+1
       enddo
       if (itt.gt.1) then
          do j=1,(nsubidx+1)/2
             subidx(ii)=idx(itt-1,2*j-1)
             !print *, 'j1=', idx(itt-1,2*j-1), om(idx(itt-1,2*j-1))
             ii=ii+1
          enddo
       else
          subidx(ii)=zero_ind
          ii=ii+1
       end if
       do j=zero_ind+(itt-1)*Nd+1, zero_ind+itt*Nd
          !print *, 'j2=', j, om(j)
          subidx(ii)=j
          ii=ii+1
       enddo
       nsubidx=ii-1
       idx(itt,:) = subidx(:)
       nidx(itt) = nsubidx
       !print *, itt, subidx(:nsubidx)
    enddo
  END SUBROUTINE LinLogMesh_Generate
END module LinLogMesh
