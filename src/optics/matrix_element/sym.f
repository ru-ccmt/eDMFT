      SUBROUTINE Sym
      
      INCLUDE 'param.inc'
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /SYM2/ TAU(3,NSYM),IORD,IMAT(3,3,NSYM)
      COMMON /SYMo/ opimat(3,3,NSYM)
      COMMON /SYMd/ det(NSYM)
      COMMON /so/ theta,fi

      DIMENSION rot(3,3),trans(3,3),rotinv(3,3)
!      parameter (pi=3.141592654d0) 
	PI=ACOS(-1.D0)
      
!      theta=theta/180.d0*pi
!      fi=fi/180.d0*pi
      ct=cos(theta)
      cf=cos(fi)
      st=sin(theta)
      sf=sin(fi)
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
      do 100 i=1,IORD
         write(6,*)'call trans..'
         call transform(trans,rotinv,opimat(1,1,i),rot)
         write(6,*) 'operation',i
         det(i)=trans(1,1)*trans(2,2)-trans(1,2)*trans(2,1)
!         if (abs(1-abs(det(i))).GT.1d-5) stop 'symm. fault'
         write(6,*)'det:', det(i) 
 100  continue
      end
      
      subroutine transform(T,Pinv,A,P)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION T(3,3),P(3,3),A(3,3),Pinv(3,3)
      
      do i=1,3
         do j=1,3
            sum=0
            do k=1,3
               do l=1,3
                  sum=sum+Pinv(i,k)*A(k,l)*P(l,j)
               end do
            end do
            T(i,j)=sum
         end do
       write(6,*)'Transf. matrix:',(T(i,k),k=1,3)
      end do
      return
      end
















