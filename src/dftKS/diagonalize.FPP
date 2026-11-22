! Solver backend shares this compilation unit with the legacy wrapper so the
! existing real/complex build rules and module dependencies remain sufficient.
! This module has no dependencies on LAPW globals; tests exercise the same code.
#ifdef _REAL
#define SCALAR real(8)
#define ADJ(X) (X)
#define GEMM dgemm
#define HEMM dsymm
#define FACTOR dsytrf
#define SOLVE dsytrs
#define CONDITION dsycon
#define DOT ddot
#else
#define SCALAR complex(8)
#define ADJ(X) conjg(X)
#define GEMM zgemm
#define HEMM zhemm
#define FACTOR zhetrf
#define SOLVE zhetrs
#define CONDITION zhecon
#define DOT zdotc
#endif
module dftks_eigensolver_backend
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use, intrinsic :: iso_c_binding, only: c_int, c_char, c_null_char
  implicit none
  private
  public :: dftks_solver_method, dftks_advanced_solve
  SCALAR, parameter :: one=1, zero=0, minusone=-1
#ifdef _REAL
  character, parameter :: trans='T', kindtag='r'
#else
  character, parameter :: trans='C', kindtag='c'
#endif
  interface
     function c_rename(old,new) bind(C,name='rename') result(status)
       import c_char,c_int
       character(c_char), intent(in) :: old(*),new(*)
       integer(c_int) :: status
     end function
     function c_getpid() bind(C,name='getpid') result(pid)
       import c_int
       integer(c_int) :: pid
     end function
  end interface
contains
  function dftks_solver_method() result(method)
    character(16) :: method
    character(80) :: text
    integer :: status,length,i,k
    method='GVX'
    call get_environment_variable('DFTKS_EIGENSOLVER',text,length,status)
    if (status==1.or.length==0) return
    if (status/=0.or.length>len(method)) then
       method='INVALID'
       return
    endif
    text=adjustl(text)
    do i=1,len_trim(text)
       k=iachar(text(i:i))
       if (k>=97.and.k<=122) text(i:i)=achar(k-32)
    enddo
    method=trim(text)
  end function

  ! HS input is the raw LAPW two-triangle storage BEFORE the old wrapper's
  ! complex conjugation: physical H_lower=conjg(HS_lower), O_upper=HS_upper.
  subroutine dftks_advanced_solve(method,n,ldz,raw,hdiag,elow,eup,fixed,nband, &
       ne,z,e,kpoint,labels,kfrac,nlo,namespace,qprint,info)
    character(*), intent(in) :: method,namespace
    integer, intent(in) :: n,ldz,nband,kpoint,nlo,labels(3,n)
    SCALAR, intent(in) :: raw(n,n),hdiag(n)
    real(8), intent(in) :: elow,eup,kfrac(3)
    logical, intent(in) :: fixed,qprint
    integer, intent(out) :: ne,info
    SCALAR, intent(out) :: z(n,ldz)
    real(8), intent(out) :: e(n)
    SCALAR, allocatable :: h(:,:),o(:,:),c(:,:),fac(:,:)
    real(8), allocatable :: values(:)
    integer, allocatable :: piv(:)
    integer :: i,j,p,wanted,below,steps,block,maxiter,maxspace,nout,first,last
    integer :: cache_unit,guard
    real(8) :: sigma,tol,maxres,ortherr
    character(2048) :: cachefile,prefix
    character(160) :: reason
    character(32) :: route
    logical :: cached,iterative,converged,settings_ok
    info=0; ne=0; steps=0; maxres=0; ortherr=0
    if (.not.(method=='GVD'.or.method=='DAVIDSON'.or.method=='AUTO')) then
       info=-100
       if (qprint) write(6,'(A)') ':DFTKS_SOLVER_ERROR: choose GVX, GVD, DAVIDSON, or AUTO'
       return
    endif
    if (n<1.or.ldz<1.or.(fixed.and.(nband<1.or.nband>min(n,ldz)))) then
       info=-101
       return
    endif
    if (.not.fixed.and.(.not.ieee_is_finite(elow).or. &
         .not.ieee_is_finite(eup).or.elow>=eup)) then
       info=-102
       return
    endif
    allocate(h(n,n),o(n,n),values(n),stat=info)
    if (info/=0) return
    do j=1,n
       h(j,j)=real(hdiag(j),8); o(j,j)=real(raw(j,j),8)
       do i=j+1,n
          h(i,j)=ADJ(raw(i,j)); h(j,i)=ADJ(h(i,j))
          o(j,i)=raw(j,i); o(i,j)=ADJ(o(j,i))
       enddo
    enddo
    iterative=method/='GVD'
    reason='dense solver requested'; route='GVD'
    if (iterative) then
       call settings(tol,block,maxiter,maxspace,settings_ok)
       if (.not.settings_ok) then
          info=-103
          if(qprint) write(6,'(A)') ':DFTKS_SOLVER_ERROR: invalid Davidson setting'
          return
       endif
       prefix=namespace
       call get_environment_variable('DFTKS_EIG_CACHE',length=i,status=j)
       if(j==0.and.i>0) then
          if(i>len(prefix)-40) then
             info=-104
             return
          endif
          call get_environment_variable('DFTKS_EIG_CACHE',prefix)
       endif
       if(len_trim(prefix)>len(cachefile)-40) then
          info=-104
          return
       endif
       write(cachefile,'(A,".dftks-",A,".k",I8.8)') trim(prefix),kindtag,kpoint
       call read_cache(trim(cachefile),n,labels,kfrac,nlo,c,values,p,cached)
       if (.not.cached) then
          reason='no compatible cache; dense seed'
       else if (fixed.and.p<=nband) then
          reason='cache lacks guard bands; dense seed'
       else if (method=='AUTO'.and.real(p,8)/n>0.25d0) then
          reason='AUTO: cached band fraction exceeds 0.25'
       else
          ! Inertia of H-sigma*O gives the exact number of eigenvalues below
          ! sigma (O positive definite). Reuse this factor as preconditioner.
          if (fixed) then
             i=min(p-1,nband+max(0,(p-nband)/2))
             sigma=values(i)+0.5d0*(values(i+1)-values(i))
          else
             sigma=eup
          endif
          allocate(fac(n,n),piv(n),stat=info)
          if(info/=0) return
          fac=h-sigma*o
          call shifted_factor(fac,n,piv,below,info)
          if(info/=0) then
             reason='singular or ambiguous inertia factor; dense fallback'
          else if ((fixed.and.below<nband).or.below>=p.or.below<1) then
             reason='inertia outside cached subspace; dense fallback'
          else
             call davidson(h,o,n,p,c,values,fac,piv,tol,block,maxiter, &
                  maxspace,steps,maxres,ortherr,converged)
             if (converged) then
                ! A small residual alone does not certify that no lower band
                ! was missed. Require all states below the inertia threshold.
                i=count(values(1:p)<sigma)
                if(i==below.and.(.not.fixed.or.values(nband)<sigma)) then
                   route='DAVIDSON'; reason='residual, metric and inertia checks passed'
                else
                   reason='inertia certificate mismatch; dense fallback'
                endif
             else
                reason='Davidson did not converge; dense fallback'
             endif
          endif
          deallocate(fac,piv)
          info=0
       endif
    endif
    if(route/='DAVIDSON') then
       if(allocated(c)) deallocate(c)
       ! No copies of the full matrices are needed after a dense fallback.
       ! GVD returns its eigenvectors in H and overwrites O with its factor.
       call dense_all(h,o,n,values,info)
       if(info/=0) return
       call move_alloc(h,c)
       p=n
    endif
    if(fixed) then
       first=1; last=nband
    else
       first=count(values(1:p)<=elow)+1
       last=count(values(1:p)<=eup)
    endif
    nout=max(0,last-first+1)
    if(nout>ldz) then
       info=-105
       if(qprint) write(6,'(A,2I8)') ':DFTKS_SOLVER_ERROR: bands exceed NUME: ',nout,ldz
       return
    endif
    ne=nout
    if(ne>0) then
       e(1:ne)=values(first:last)
       z(:,1:ne)=c(:,first:last)
    endif
    if(iterative) then
       if(fixed) then
          wanted=nband
       else
          wanted=last
       endif
       guard=max(8,(wanted+9)/10)
       cache_unit=min(p,wanted+guard)
       if(cache_unit>0) call write_cache(trim(cachefile),n,cache_unit,labels, &
            kfrac,nlo,c,values,qprint)
    endif
    if(qprint) then
       write(6,*) ne,' Eigenvalues computed '
       if(fixed) write(6,'(A,I8)') ':FIXED_NBANDS: number of lowest eigenpairs =',ne
       write(6,'(A,A,A,I8,A,I8)') ':DFTKS_SOLVER: ',trim(route), &
            ' basis=',n,' returned=',ne
       write(6,'(A,A)') ':DFTKS_SOLVER_REASON: ',trim(reason)
       if(steps>0) write(6,'(A,I5,A,ES12.4,A,ES12.4)') &
            ':DFTKS_DAVIDSON: steps=',steps,' residual=',maxres,' metric_error=',ortherr
    endif
  end subroutine

  subroutine settings(tol,block,maxiter,maxspace,ok)
    real(8), intent(out) :: tol
    integer, intent(out) :: block,maxiter,maxspace
    logical, intent(out) :: ok
    character(80) :: text
    integer :: status,length,ios,i
    character(32), parameter :: names(4)=[character(32):: &
         'DFTKS_DAVIDSON_TOL','DFTKS_DAVIDSON_BLOCK', &
         'DFTKS_DAVIDSON_MAXITER','DFTKS_DAVIDSON_MAXSPACE']
    tol=1d-10; block=16; maxiter=40; maxspace=0; ok=.true.
    do i=1,4
       call get_environment_variable(trim(names(i)),text,length,status)
       if(status==1.or.length==0) cycle
       if(status/=0) then
          ok=.false.; return
       endif
       select case(i)
       case(1)
          read(text,*,iostat=ios) tol
       case(2)
          read(text,*,iostat=ios) block
       case(3)
          read(text,*,iostat=ios) maxiter
       case(4)
          read(text,*,iostat=ios) maxspace
       end select
       if(ios/=0) ok=.false.
    enddo
    ok=ok.and.ieee_is_finite(tol).and.tol>0d0.and.tol<1d-2 &
         .and.block>0.and.maxiter>0.and.maxspace>=0
  end subroutine

  subroutine dense_all(h,o,n,e,info)
    integer, intent(in) :: n
    SCALAR, intent(inout) :: h(n,n),o(n,n)
    real(8), intent(out) :: e(n)
    integer, intent(out) :: info
    SCALAR, allocatable :: work(:)
    SCALAR :: query(1)
    integer :: iwquery(1),lw,liw
    integer, allocatable :: iw(:)
#ifndef _REAL
    real(8) :: rq(1)
    real(8), allocatable :: rw(:)
    integer :: lrw
#endif
#ifdef _REAL
    call dsygvd(1,'V','L',n,h,n,o,n,e,query,-1,iwquery,-1,info)
#else
    call zhegvd(1,'V','L',n,h,n,o,n,e,query,-1,rq,-1,iwquery,-1,info)
#endif
    if(info/=0) return
    lw=max(1,int(real(query(1),8))); liw=max(1,iwquery(1))
    allocate(work(lw),iw(liw),stat=info)
    if(info/=0) return
#ifdef _REAL
    call dsygvd(1,'V','L',n,h,n,o,n,e,work,lw,iw,liw,info)
#else
    lrw=max(1,int(rq(1))); allocate(rw(lrw),stat=info)
    if(info/=0) return
    call zhegvd(1,'V','L',n,h,n,o,n,e,work,lw,rw,lrw,iw,liw,info)
#endif
  end subroutine

  subroutine small_eigh(a,n,e,info)
    integer, intent(in) :: n
    SCALAR, intent(inout) :: a(n,n)
    real(8), intent(out) :: e(n)
    integer, intent(out) :: info
    SCALAR :: query(1)
    SCALAR, allocatable :: work(:)
    integer :: lw,liw,iq(1)
    integer, allocatable :: iw(:)
#ifndef _REAL
    real(8) :: rq(1)
    real(8), allocatable :: rw(:)
    integer :: lrw
#endif
#ifdef _REAL
    call dsyevd('V','L',n,a,n,e,query,-1,iq,-1,info)
#else
    call zheevd('V','L',n,a,n,e,query,-1,rq,-1,iq,-1,info)
#endif
    if(info/=0) return
    lw=max(1,int(real(query(1),8))); liw=max(1,iq(1))
    allocate(work(lw),iw(liw),stat=info)
    if(info/=0) return
#ifdef _REAL
    call dsyevd('V','L',n,a,n,e,work,lw,iw,liw,info)
#else
    lrw=max(1,int(rq(1))); allocate(rw(lrw),stat=info)
    if(info/=0) return
    call zheevd('V','L',n,a,n,e,work,lw,rw,lrw,iw,liw,info)
#endif
  end subroutine

  subroutine shifted_factor(a,n,piv,below,info)
    integer, intent(in) :: n
    SCALAR, intent(inout) :: a(n,n)
    integer, intent(out) :: piv(n),below,info
    SCALAR :: query(1)
    SCALAR, allocatable :: work(:)
    integer :: lw,i
#ifdef _REAL
    integer, allocatable :: iwork(:)
#endif
    real(8) :: x,y,b,scale,det,anorm,rcond
    anorm=maxval(sum(abs(a),dim=1))
    call FACTOR('L',n,a,n,piv,query,-1,info)
    below=0
    if(info/=0) return
    lw=max(1,int(real(query(1),8))); allocate(work(lw),stat=info)
    if(info/=0) return
    call FACTOR('L',n,a,n,piv,work,lw,info)
    if(info/=0) return
    deallocate(work)
    allocate(work(2*n),stat=info)
    if(info/=0) return
#ifdef _REAL
    allocate(iwork(n),stat=info)
    if(info/=0) return
    call CONDITION('L',n,a,n,piv,anorm,rcond,work,iwork,info)
#else
    call CONDITION('L',n,a,n,piv,anorm,rcond,work,info)
#endif
    if(info/=0) return
    if(.not.ieee_is_finite(rcond).or.rcond<1d-12) then
       info=-1; return
    endif
    i=1
    do while(i<=n)
       x=real(a(i,i),8)
       if (.not.ieee_is_finite(x)) then
          info=-1; return
       endif
       if(piv(i)>0) then
          if(x==0d0) then
             info=-1; return
          endif
          if(x<0d0) below=below+1
          i=i+1
       else
          if(i==n) then
             info=-1; return
          endif
          y=real(a(i+1,i+1),8); b=abs(a(i+1,i))
          scale=max(abs(x),abs(y),b)
          if(.not.ieee_is_finite(scale).or.scale==0d0) then
             info=-1; return
          endif
          x=x/scale; y=y/scale; b=b/scale
          det=x*y-b*b
          ! Avoid cancellation when forming a very small eigenvalue of D.
          if(abs(det)<100d0*epsilon(1d0)) then
             info=-1; return
          endif
          if(det<0d0) then
             below=below+1
          else if(x+y<0d0) then
             below=below+2
          endif
          i=i+2
       endif
    enddo
  end subroutine

  ! Twice project in the current metric, then diagonalize the SMALL Gram
  ! matrix. Discard dependent combinations, never modify the physical O.
  subroutine metric_qr(o,n,t,b,v,ov,q,ot,rank,info)
    integer, intent(in) :: n,b,q
    SCALAR, intent(in) :: o(n,n),v(n,q),ov(n,q)
    SCALAR, intent(inout) :: t(n,b)
    SCALAR, intent(out) :: ot(n,b)
    integer, intent(out) :: rank,info
    SCALAR, allocatable :: gram(:,:),projection(:,:),rotated(:,:)
    real(8), allocatable :: eval(:)
    real(8) :: scale
    integer :: pass,i,j
    rank=0; info=0
    call HEMM('L','L',n,b,one,o,n,t,n,zero,ot,n)
    if(q>0) then
       allocate(projection(q,b))
       do pass=1,2
          call GEMM(trans,'N',q,b,n,one,v,n,ot,n,zero,projection,q)
          call GEMM('N','N',n,b,q,minusone,v,n,projection,q,one,t,n)
          call GEMM('N','N',n,b,q,minusone,ov,n,projection,q,one,ot,n)
       enddo
    endif
    ! Normalize column scales before deciding relative numerical rank.
    do j=1,b
       scale=sqrt(max(0d0,real(sum(ADJ(t(:,j))*ot(:,j)),8)))
       if(scale>tiny(1d0)**0.25d0) then
          t(:,j)=t(:,j)/scale; ot(:,j)=ot(:,j)/scale
       else
          t(:,j)=zero; ot(:,j)=zero
       endif
    enddo
    allocate(gram(b,b),eval(b),rotated(n,b))
    call GEMM(trans,'N',b,b,n,one,t,n,ot,n,zero,gram,b)
    gram=(gram+transpose(ADJ(gram)))*0.5d0
    call small_eigh(gram,b,eval,info)
    if(info/=0.or..not.all(ieee_is_finite(eval))) then
       info=-1; return
    endif
    scale=maxval(eval)
    if(scale<=0d0) return
    if(minval(eval)<-1d-9*scale) then
       info=-1; return
    endif
    do i=1,b
       if(eval(i)<=1d-10*scale) cycle
       rank=rank+1
       gram(:,rank)=gram(:,i)/sqrt(eval(i))
    enddo
    if(rank==0) return
    call GEMM('N','N',n,rank,b,one,t,n,gram,b,zero,rotated,n)
    t(:,1:rank)=rotated(:,1:rank)
    ! Reapply O after projection to avoid carrying cancellation errors.
    call HEMM('L','L',n,rank,one,o,n,t,n,zero,ot,n)
  end subroutine

  subroutine davidson(h,o,n,p,c,e,fac,piv,tol,block,maxiter,requested_space, &
       steps,maxres,ortherr,converged)
    integer, intent(in) :: n,p,block,maxiter,requested_space,piv(n)
    SCALAR, intent(in) :: h(n,n),o(n,n),fac(n,n)
    SCALAR, intent(inout) :: c(n,p)
    real(8), intent(out) :: e(n),maxres,ortherr
    real(8), intent(in) :: tol
    integer, intent(out) :: steps
    logical, intent(out) :: converged
    SCALAR, allocatable :: v(:,:),hv(:,:),ov(:,:),hc(:,:),oc(:,:),r(:,:)
    SCALAR, allocatable :: t(:,:),ot(:,:),rr(:,:),metric(:,:)
    real(8), allocatable :: eval(:),res(:),selection(:)
    integer :: space,q,b,rank,info,i,j,idx,iteration
    converged=.false.; steps=0; maxres=huge(1d0); ortherr=huge(1d0)
    space=min(n,max(2*p,p+2*min(block,n-p)))
    if(requested_space>0) space=min(n,requested_space)
    if(space<=p.or.p<1) return
    b=min(block,space-p)
    allocate(v(n,space),hv(n,space),ov(n,space),hc(n,p),oc(n,p),r(n,p), &
         t(n,max(p,b)),ot(n,max(p,b)),res(p),selection(p),metric(p,p))
    t(:,1:p)=c
    call metric_qr(o,n,t(:,1:p),p,v(:,1:0),ov(:,1:0),0,ot(:,1:p),rank,info)
    if(info/=0.or.rank/=p) return
    v(:,1:p)=t(:,1:p); ov(:,1:p)=ot(:,1:p); q=p
    call HEMM('L','L',n,p,one,h,n,v,n,zero,hv,n)
    do iteration=1,maxiter
       steps=iteration
       allocate(rr(q,q),eval(q))
       call GEMM(trans,'N',q,q,n,one,v,n,hv,n,zero,rr,q)
       rr=(rr+transpose(ADJ(rr)))*0.5d0
       call small_eigh(rr,q,eval,info)
       if(info/=0) return
       call GEMM('N','N',n,p,q,one,v,n,rr,q,zero,c,n)
       call GEMM('N','N',n,p,q,one,hv,n,rr,q,zero,hc,n)
       call GEMM('N','N',n,p,q,one,ov,n,rr,q,zero,oc,n)
       e(1:p)=eval(1:p)
       deallocate(rr,eval)
       do i=1,p
          r(:,i)=hc(:,i)-e(i)*oc(:,i)
          res(i)=sqrt(sum(abs(r(:,i))**2))/max(1d0, &
               sqrt(sum(abs(hc(:,i))**2))+abs(e(i))*sqrt(sum(abs(oc(:,i))**2)))
       enddo
       maxres=maxval(res)
       if(.not.all(ieee_is_finite(res))) return
       if(maxres<=tol) then
          call GEMM(trans,'N',p,p,n,one,c,n,oc,n,zero,metric,p)
          do i=1,p
             metric(i,i)=metric(i,i)-one
          enddo
          ortherr=maxval(abs(metric))
          converged=ortherr<=max(1d-10,10d0*tol)
          return
       endif
       if(q+b>space) then
          v(:,1:p)=c; hv(:,1:p)=hc; ov(:,1:p)=oc; q=p
       endif
       selection=res
       rank=0
       do i=1,min(b,count(res>tol))
          idx=maxloc(selection,dim=1); selection(idx)=-1d0
          t(:,i)=r(:,idx); rank=i
       enddo
       call SOLVE('L',n,rank,fac,n,piv,t,n,info)
       if(info/=0) return
       j=rank
       call metric_qr(o,n,t(:,1:j),j,v(:,1:q),ov(:,1:q),q,ot(:,1:j),rank,info)
       if(info/=0.or.rank==0) return
       v(:,q+1:q+rank)=t(:,1:rank)
       ov(:,q+1:q+rank)=ot(:,1:rank)
       call HEMM('L','L',n,rank,one,h,n,t,n,zero,hv(:,q+1:q+rank),n)
       q=q+rank
    enddo
  end subroutine

  subroutine read_cache(file,n,labels,kfrac,nlo,c,e,p,ok)
    character(*), intent(in) :: file
    integer, intent(in) :: n,nlo,labels(3,n)
    real(8), intent(in) :: kfrac(3)
    SCALAR, allocatable, intent(out) :: c(:,:)
    real(8), intent(out) :: e(n)
    integer, intent(out) :: p
    logical, intent(out) :: ok
    integer :: unit,ios,header(4),oldlabels(3,n)
    integer(8) :: bytes,expected
    real(8) :: oldk(3)
    character(8) :: magic
    character :: tag
    ok=.false.; p=0
    open(newunit=unit,file=file,status='old',access='stream',form='unformatted',iostat=ios)
    if(ios/=0) return
    read(unit,iostat=ios) magic,tag,header,oldk
    if(ios==0) then
       if(magic/='DFTKSE01'.or.tag/=kindtag.or.header(1)/=1.or. &
            header(2)/=n.or.header(3)<1.or.header(3)>n.or.header(4)/=nlo) ios=-1
       if(.not.all(ieee_is_finite(oldk))) ios=-1
       if(maxval(abs(oldk-kfrac))>1d-12) ios=-1
    endif
    if(ios==0) then
       p=header(3)
       inquire(unit=unit,size=bytes)
       expected=8+1+4*4+3*8+int(3*n,8)*4+int(p,8)*8 &
            +int(n,8)*p*(storage_size(one)/8)
       if(bytes/=expected) ios=-1
    endif
    if(ios==0) then
       read(unit,iostat=ios) oldlabels
       if(ios==0) then
          if(any(oldlabels/=labels)) ios=-1
       endif
    endif
    if(ios==0) then
       allocate(c(n,p),stat=ios)
       if(ios==0) read(unit,iostat=ios) e(1:p),c
       if(ios==0) then
          if(.not.all(ieee_is_finite(e(1:p)))) ios=-1
#ifdef _REAL
          if(.not.all(ieee_is_finite(c))) ios=-1
#else
          if(.not.all(ieee_is_finite(real(c,8))).or. &
               .not.all(ieee_is_finite(aimag(c)))) ios=-1
#endif
       endif
    endif
    close(unit)
    ok=ios==0
    if(.not.ok) then
       p=0
       if(allocated(c)) deallocate(c)
    endif
  end subroutine

  subroutine write_cache(file,n,p,labels,kfrac,nlo,c,e,qprint)
    character(*), intent(in) :: file
    integer, intent(in) :: n,p,nlo,labels(3,n)
    real(8), intent(in) :: kfrac(3),e(*)
    SCALAR, intent(in) :: c(n,*)
    logical, intent(in) :: qprint
    character(:), allocatable :: temporary
    character(32) :: pidtext
    integer :: unit,ios,close_status
    write(pidtext,'(I0)') c_getpid()
    temporary=file//'.tmp.'//trim(pidtext)
    open(newunit=unit,file=temporary,status='replace',access='stream', &
         form='unformatted',iostat=ios)
    if(ios==0) then
       write(unit,iostat=ios) 'DFTKSE01',kindtag,[1,n,p,nlo],kfrac,labels,e(1:p),c(:,1:p)
       close(unit,iostat=close_status)
       if(ios==0) ios=close_status
       if(ios==0) ios=c_rename(temporary//c_null_char,file//c_null_char)
    endif
    if(ios/=0.and.qprint) write(6,'(A,A)') ':DFTKS_CACHE_WARNING: could not publish ',file
  end subroutine
end module
#undef SCALAR
#undef ADJ
#undef GEMM
#undef HEMM
#undef FACTOR
#undef SOLVE
#undef CONDITION
#undef DOT
#ifdef _REAL
#define tpcast(X) dble(X)  !! typecasting in fortran
#define XSWAP DSWAP
#define Ttype REAL*8
#define tsize 8
#define routine_name 'DSYGVX'
#else
#define tpcast(X) X
#define XSWAP ZSWAP
#define Ttype COMPLEX*16
#define tsize 16
#define routine_name 'ZHEGVX'
#endif
SUBROUTINE Diagonalize(Elow, Eup, NE, HS, HSDIAG, Z, EIGVAL, HSROWS, LCOLZ, Qprint, &
     kpoint,labels,kfrac,nlo,cache_namespace)
  use omp_lib
  use comi, only: nemax_it, fixed_nbands
  use dftks_eigensolver_backend, only: dftks_solver_method,dftks_advanced_solve
  IMPLICIT NONE
  !..................................................................
  ! The interface between matrix set up in LAPW1 and LAPACK (BLAS) matrix subroutines
  !
  !   Solve the secular equations H*Z = EIGVAL*S*Z
  !   H and S are stored in the array HS,
  !   H is stored below the diagonal, whereas S is stored above
  !   the diagonal. The diagonal of H is stored in the 
  !   one dimensional array HSDIAG.
  !   HSROWS is the dimension of the secular matrices
  !   With the explicit FIXED_NBANDS keyword, the lowest nband eigenpairs
  !   are returned exactly.  Otherwise NE is the number of eigenvalues
  !   between Elow and Eup, preserving the historical energy-window behavior.
  !..................................................................
  REAL*8,  intent(in) :: Elow, Eup
  INTEGER, intent(in) :: HSROWS, LCOLZ
  INTEGER, intent(out):: NE
  Ttype, intent(inout):: HS(HSROWS,HSROWS)
  Ttype, intent(inout):: HSDIAG(HSROWS)
  Ttype,   intent(out):: Z(HSROWS,LCOLZ)
  REAL*8,  intent(out):: eigval(HSROWS)
  LOGICAL, intent(in) :: Qprint
  INTEGER, intent(in) :: kpoint,nlo,labels(3,HSROWS)
  REAL*8, intent(in) :: kfrac(3)
  CHARACTER(*), intent(in) :: cache_namespace
  !   Locals
  Ttype,   allocatable :: Olap(:,:)
  Ttype,   allocatable :: work(:)
  INTEGER, allocatable :: ifail(:), iwork(:)
#ifndef _REAL  
  REAL*8,  allocatable :: rwork(:)
#endif  
  REAL*8             :: t1, t2, t1w, t2w, abstol, Flop
  INTEGER            :: i, j, info, lwork, ierr, il, iu
  CHARACTER*67       :: ERRMSG
  CHARACTER(len=1)   :: eig_range
  LOGICAL            :: out, out2, fixed_band_count
  CHARACTER(16)      :: solver_method

  solver_method=dftks_solver_method()
  if(solver_method/='GVX') then
     t1w=omp_get_wtime()
     call dftks_advanced_solve(trim(solver_method),HSROWS,LCOLZ,HS,HSDIAG, &
          Elow,Eup,fixed_nbands,nemax_it,NE,Z,EIGVAL,kpoint,labels,kfrac, &
          nlo,cache_namespace,Qprint,info)
     if(info/=0) then
        write(ERRMSG,'(A,I8)') 'New eigensolver failed, INFO=',info
        call OUTERR('diagonalize',ERRMSG)
        stop 1
     endif
     if(Qprint) write(6,'(A,F12.4)') ':DFTKS_SOLVER_WALL: ',omp_get_wtime()-t1w
     return
  endif
  
  out =.false.
  out2=.false.
  ! A finite nemax_it alone retains the historical energy-window behavior.
  ! Index selection is deliberately opt-in through FIXED_NBANDS in case.in1.
  fixed_band_count=fixed_nbands
  if (fixed_band_count) then
     ! RANGE='I' makes the fourth field of the K-VECTORS line in case.in1 a
     ! genuine fixed band count.  Selecting by eigenvalue index avoids a
     ! discontinuous number of stored bands when a dispersive high-energy
     ! state crosses Eup under strain.
     eig_range='I'
     il=1
     iu=nemax_it
     if (nemax_it.lt.1 .or. nemax_it.gt.HSROWS) then
        WRITE(ERRMSG,'(A,I8,A,I8)') 'Requested nband=',nemax_it, &
             ' is incompatible with matrix size=',HSROWS
        CALL OUTERR('diagonalize',ERRMSG)
        STOP 'invalid fixed band count'
     endif
     if (nemax_it.gt.LCOLZ) then
        WRITE(ERRMSG,'(A,I8,A,I8)') 'Requested nband=',nemax_it, &
             ' exceeds allocated eigenvector columns=',LCOLZ
        CALL OUTERR('diagonalize',ERRMSG)
        STOP 'increase NUME for fixed band count'
     endif
  else
     eig_range='V'
     il=0
     iu=0
  endif
  !out=.true.
  if (out) then
     open(86,FILE='_eigvals_',status='unknown')
     open(89,FILE='_eigvect_',status='unknown')
  endif
  if (out2) then
     open(87,FILE='_h_diag_',status='unknown')
     open(88,FILE='_h_and_o_',status='unknown')
     !
     write(87,*) HSrows!,Nvaa,Nlo
     write(88,'(4e23.15)') HS
     write(87,'(4e23.15)') HSdiag
  endif

#ifndef _REAL  
  !  conjugate entire matrix. When we computed H,O we stored it transposed. Now we conjugate, to get correct O & H.
  HS(:,:) = conjg(HS(:,:))
  HSDIAG(:) = conjg(HSDIAG)
#endif

  
  allocate( Olap(HSROWS, HSROWS), STAT=IERR )
  if(IERR.ne.0) THEN 
     write(6,*) 'not enough memory: ', IERR
     stop 'not enough memory' 
  end if
  Olap(:,:)=0.d0
  do i=1,HSROWS
     Olap(i,1:i) = HS(1:i,i)  ! By transposing here, we effectively conjugated overlap
     HS(i,i) = HSDIAG(i)      ! Diagonal of Ham is set to HS
  enddo  ! Now Olap and HS are lower-diagonal matrices of Overlap and Hamiltonian
  
  allocate( iwork(5*HSROWS), ifail(HSROWS), work(2), STAT=IERR )
  if(IERR.ne.0) THEN 
     write(6,*) 'not enough memory: ', IERR
     stop 'not enough memory' 
  end if
  
  abstol=-1.d0 ! set best tolerance
  ! First we do a query for the best size of lwork
  lwork=-1  ! Set lwork to -1 for query
#ifdef _REAL  
  call DSYGVX(1, 'V', eig_range, 'L', HSROWS, HS, HSROWS, Olap, &
       HSROWS, Elow, Eup, il, iu, abstol, NE, EIGVAL, Z, HSROWS, &
       WORK, lwork, IWORK, ifail, info)
#else
  allocate( rwork(7*HSROWS), STAT=IERR )
  if(IERR.ne.0) THEN
     write(6,*) 'not enough memory for ZHEGVX real workspace: ', IERR
     stop 'not enough memory'
  end if
  call ZHEGVX(1, 'V', eig_range, 'L', HSROWS, HS, HSROWS, Olap, &
       HSROWS, Elow, Eup, il, iu, abstol, NE, EIGVAL, Z, HSROWS, &
       WORK, lwork, RWORK, IWORK, ifail, info)
#endif
  if (info.ne.0) then
     write(6,*) routine_name, ' workspace query failed, INFO = ', info
     stop 'LAPACK workspace query failed'
  endif
  lwork = int(work(1)+0.1) ! This is the best lwork
  !print *, 'lwork=', work(1), lwork
  deallocate( work )
  allocate(work(lwork), STAT=IERR )
  if(IERR.ne.0) THEN 
     write(6,*) 'not enough memory: ', IERR
     stop 'not enough memory' 
  end if
  ! Real diagonalization
  t1w = omp_get_wtime ( )
  call cpu_time ( t1 )
#ifdef _REAL  
  call DSYGVX(1, 'V', eig_range, 'L', HSROWS, HS, HSROWS, Olap, &
       HSROWS, Elow, Eup, il, iu, abstol, NE, EIGVAL, Z, HSROWS, &
       WORK, lwork, IWORK, ifail, info)
#else  
  call ZHEGVX(1, 'V', eig_range, 'L', HSROWS, HS, HSROWS, Olap, &
       HSROWS, Elow, Eup, il, iu, abstol, NE, EIGVAL, Z, HSROWS, &
       WORK, lwork, RWORK, IWORK, ifail, info)
  deallocate( rwork )
#endif
  t2w = omp_get_wtime ( )
  call cpu_time ( t2 )
  
  deallocate(iwork, work)
  deallocate( Olap )
  
  if (INFO.NE.0) THEN
     write(6,*) routine_name, 'INFO = ', INFO
     write(99,*) routine_name, 'INFO = ', INFO
     CALL OUTERR(routine_name,'diagonalize(LAPACK) failed.')
     STOP 'diagonalize - Error'
  ENDIF
  if (fixed_band_count .and. NE.ne.nemax_it) then
     WRITE(ERRMSG,'(A,I8,A,I8)') 'FIXED_NBANDS requested ',nemax_it, &
          ' eigenpairs but solver returned ',NE
     CALL OUTERR('diagonalize',ERRMSG)
     STOP 'FIXED_NBANDS count mismatch'
  endif
  IF (NE>LCOLZ) THEN
     CALL OUTERR('diagonalize','More than NUME eigenvalues.')
     WRITE (ERRMSG,'(A,i4,A,i4,A)') 'Increase NUME from',LCOLZ,'to',NE,'   OR'
     CALL OUTERR('diagonalize',ERRMSG)
     WRITE(ERRMSG, '(A)') 'reduce E2 in *.klist.(see *.output1)'
     CALL OUTERR('diagonalize',ERRMSG)
     STOP 'diagonalize - Error'
  ENDIF
  
  if (Qprint) then
     if (fixed_band_count) WRITE(6,'(A,I8)') &
          ':FIXED_NBANDS: exact number of lowest eigenpairs =',nemax_it
     write(6,*) NE, ' Eigenvalues computed '
     write(6,*)  'info=', info
  endif
  deallocate( ifail )
  
  Flop = 4.d0*dble(HSROWS)**3* (tsize/8.d0)**2/3.d0
  if (Qprint) then
     if (t2.gt.t1) then
        WRITE (6,'(1X,A,A,1x,f10.3,2x,f10.3,2x,f10.2,A)') &
             routine_name,'(Compute eigenvalues (CPU)) :', (t2-t1), &
             (t2w-t1w), (1E-9*Flop/(t2-t1)),' Gflops'
     else
        WRITE (6,'(1X,A,A,1x,f10.3,2x,f10.3)') &
             routine_name,'(Compute eigenvalues (CPU)) :', (t2-t1), (t2w-t1w)
     endif
  endif

  ! Below is printing for debuggig
  if(out) then
     do i=1,ne
        write(86,*) eigval(i)
     enddo
     
     write(89,'(A)',advance='no') '# basis   '
     do i=1,ne
#ifdef _REAL        
        write(89,'(A5,I4,14X)',advance='no') 'ibnd=',i
#else
        write(89,'(A7,I4,12X)',advance='no') 'RE:ibnd=',i
        write(89,'(A7,I4,12X)',advance='no') 'Im:ibnd=',i
#endif        
     enddo
     write(89,*)
     do j=1,hsrows
        write(89,'(I4)',advance='no') j
        do i=1,ne
#ifdef _REAL           
           write(89,'(f23.15)',advance='no') Z(j,i)
#else
           write(89,'(2f23.15)',advance='no') Z(j,i)
#endif           
        enddo
        write(89,*)
     enddo
     !
     close(86)
     close(89)
  endif
  if (out2) then
     close(87)
     close(88)
  endif
END SUBROUTINE Diagonalize
