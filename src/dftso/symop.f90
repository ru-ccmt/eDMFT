SUBROUTINE symop(theta,fi)
  !     for .orhto. systems   opimat == imat
  !     for not.ortho. systems opimat == BR1 . imat . BR1^-1 
  USE rotmat, ONLY: det, phase
  USE structure, ONLY: BR1, rotij, mult, rotloc
  USE param, ONLY: nato
  USE mpi, ONLY: Qprint
  IMPLICIT NONE
  REAL*8, intent(in) :: theta, fi
  ! external
  REAL*8 :: determinant
  ! locals
  REAL*8 :: BR1in(3,3), opimat(3,3), rot(3,3), rotinv(3,3), trans(3,3), tmp(3,3)
  REAL*8 :: rloc(3,3), sloc(3,3)!, th(10), ph(10)
  REAL*8 :: pi, ct, cf, st, sf, a, b, c, DD, tt
  INTEGER :: indj, jatom, mu, ii, jj
  pi=acos(-1.d0)
  !.........inverssymdef calculates the inverse of an 3*3 matrix.......
  CALL INVERSSYMDEF(BR1,BR1in)
  ct=dcos(theta)
  cf=dcos(fi)
  st=dsin(theta)
  sf=dsin(fi)
!!! rot is the spin coordinate system, in which case.inso (h,k,l) points in z direction
!!! more precisely, it takes the following form
!!!  ( e_{theta}, e_{phi}, e_r ),
!!! where e_{theta}, e_{phi} and e_r are unit vectors in spherical coordinate system.  
  rot(1,1)=cf*ct
  rot(1,2)=-sf
  rot(1,3)=cf*st
  rot(2,1)=sf*ct
  rot(2,2)=cf
  rot(2,3)=sf*st
  rot(3,1)=-st
  rot(3,2)=0
  rot(3,3)=ct
  call INVERSSYMDEF(rot,rotinv)
  ! sloc  = Det* rotloc * (BR1 * rotij * BR1^{-1}) * rot
  !   my understanding:
  !   sloc :  local_to_first <- from_first_to_each <- global <- spin-coordinate-system
  !
  ! trans =    rot^{-1} * (BR1 * rotij * BR1^{-1}) * rot
  !   my understanding:
  !   trans : spin-coordinate-system <- global <- spin-coordinate-system
  indj=0
  do jatom=1,nato
     do mu=1,mult(jatom)
        indj=indj+1
!!!     rotij specified the transformation between the first atom of certain type
!!!         and all other atoms of the same type 
!!!     opimat is like rotij in cartesian coordinates, i.e.,
!!!         opimat = BR1 * rotij * BR1^{-1}
        tmp(:,:) = matmul(BR1(:,:),rotij(:,:,indj))
        opimat(:,:) = matmul(tmp,BR1in(:,:))
!!!     trans is now rotij transofrmed to spin coordinate system,
!!!         namely, in the axis in which spin is quantizied
        call transform(trans,rotinv,opimat,rot)
!!!     rloc transforms into local coordinate system attached to indj atom.
!!!         To achieve this, we first transformed (using rotij) from the coordinate
!!!         system attached to the first atom, to the other atom, and then into
!!!         local coordinate system (using rotloc).        
        rloc(:,:) = matmul(rotloc(:,:,jatom),opimat(:,:))
!!!     Here we will perform proper rotation only (inversion does not change spin).
!!!     sloc is spin coordinate system written in local coordinate system.
        tt=determinant(rloc)
        sloc = matmul(rloc,rot) * tt

        if (Qprint) then
           write(6,*)'ATOM:',indj
           write(6,*)'local coord. system:'
           do ii=1,3
              write(6,5)(rloc(ii,jj),jj=1,3)
           end do
           write(6,*)'rotij:'
           do ii=1,3
              write(6,5)(rotij(ii,jj,indj),jj=1,3)
           end do
           write(6,*)'opimat:'
           do ii=1,3
              write(6,5)(opimat(ii,jj),jj=1,3)
           end do
           write(6,*)
           write(6,*)'spin with respect to local coord.:'
           do ii=1,3
              write(6,5)(sloc(ii,jj),jj=1,3)
           end do
           write(6,*)
        endif
     
        call euler(sloc,a,b,c)
        !! Now we calculate <Y_{lms}|l*s|Y_{l'm's'}> matrix elements
        !! but we rotate the spin into the local coordinate system of each atom,
        !! in which L is specified.
        call couple(indj,a,b,c)
        
        det(indj)=trans(1,1)*trans(2,2)-trans(1,2)*trans(2,1)
        if (abs(1.-abs(det(indj))).gt.1.d-2) then
           write(6,*)'WRONG SYMMETRY in spin pol. case'
           write(6,555)indj,det(indj)
        endif

        DD=determinant(trans)
        if (DD.lt.-0.5) then
           trans(:,:) = -trans(:,:)
        end if
        
        call euler(trans,a,b,c)
        !       write(6,3)indj,a*180/pi,b*180/pi,c*180/pi
        if (det(indj).gt.0.5) then
           phase(indj)=a+c
        else
           phase(indj)=a-c
        end if
        !       write(6,4)det(indj),phase(indj)
     enddo
  enddo
!3 format(i2,1x,'Euler angles: a,b,c:  ',3f7.1)
!4 format('det=',f3.0,' phase=',f8.4)
5 format(3f12.7)
555 format(' atom',i3,' det',f10.5)
end SUBROUTINE symop

subroutine transform(T,Pinv,A,P)
  IMPLICIT REAL*8 (A-H,O-Z)
  DIMENSION T(3,3),P(3,3),A(3,3),Pinv(3,3), tmp(3,3)
  tmp = matmul(Pinv,A)
  T = matmul(tmp,P)
  return
end subroutine transform

REAL*8 FUNCTION determinant(A)
  REAL* 8 A(3,3)
  determinant=a(1,1)*a(2,2)*a(3,3)+a(1,2)*a(2,3)*a(3,1)+a(1,3)*a(2,1)*a(3,2)-a(3,1)*a(2,2)*a(1,3)-a(1,1)*a(3,2)*a(2,3)-a(2,1)*a(1,2)*a(3,3)
END FUNCTION determinant
