subroutine bisection ( f, fatol, itmax, xl, xr, xatol, xrtol )
!
!*******************************************************************************
!
!! BISECTION carries out the bisection algorithm.
!
!
!  Modified:
!
!    19 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real FATOL, the absolute error tolerance for F(X).
!
!    Input, integer ITMAX, the maximum number of steps allowed.
!
!    Input, real X1, X2, two points where the function differs in sign.  
!
!    Input, real XATOL, absolute error tolerance for the root.
!
  implicit none
!
  real*8 f
  real*8 fatol
  real*8 fxl
  real*8 fxm
  real*8 fxr
  integer iterate
  integer itmax
  real*8 x_ave
  !real*8 x1
  !real*8 x2
  real*8 xatol
  real*8 xl
  real*8 xm
  real*8 xr
  real*8 xrtol
!
  write ( *, * ) ' '
  write ( *, * ) ' '
  write ( *, * ) 'The bisection method'
  write ( *, * ) ' '
  write ( *, * ) ' '

  !xl = x1
  !xr = x2
  !call f ( xl, fxl, 0 )
  !call f ( xr, fxr, 0 )
  fxl = f(xl)
  fxr = f(xr)
!
!  Check for change of sign.
!
  if ( fxl * fxr > 0.0E+00 ) then
    write ( *, * ) ' '
    write ( *, * ) 'WARNING:'
    write ( *, * ) ' '
    write ( *, * ) '  This method requires a change of sign '
    write ( *, * ) '  interval, but your function does NOT change '
    write ( *, * ) '  sign over the given interval, so this method '
    write ( *, * ) '  cannot be used.'
    write ( *, * ) ' '
    write ( *, * ) 'SUGGESTION:'
    write ( *, * ) ' '
    write ( *, * ) '  Try a different method;'
    write ( *, * ) '  Try changing X1 or X2;'
    write ( *, * ) '  Perhaps your function F(X) is wrong.'
    return
  end if

  write ( *, * ) ' Step           Left          Middle         Right'

  iterate = 0

  write ( *, * ) ' '
  write ( *, '(i4,'' X'',g16.8,16x,f16.8)' ) iterate,  xl,  xr
  write ( *, '(4x,''FX'',g16.8,16x,g16.8)' )          fxl, fxr
!
!  Take another step of the iteration.
!
  do

    iterate = iterate + 1
    x_ave = ( abs ( xr ) + abs ( xl ) ) / 2.0E+00

    if ( iterate > itmax ) then
      write ( *, * ) ' '
      write ( *, * ) '  Maximum number of steps taken.'
      exit
    end if

    if ( abs ( xr - xl ) <= xatol ) then
      write ( *, * ) ' '
      write ( *, * ) '  Absolute convergence of the X interval.'
      exit
    end if

    if ( abs ( xr - xl ) <= xrtol * x_ave ) then
      write ( *, * ) ' '
      write ( *, * ) '  Absolute convergence of the X interval.'
      exit
    end if

    if ( abs ( fxl ) <= fatol ) then
      write ( *, * ) ' '
      write ( *, * ) '  Absolute convergence of |F(X)|.'
      exit   
    end if 

    if ( abs ( fxr ) <= fatol ) then
      write ( *, * ) ' '
      write ( *, * ) '  Absolute convergence of |F(X)|.'
      exit   
    end if 
!
!  Compute the next iterate.
!
    xm = ( xl + xr ) / 2.0E+00
    !call f ( xm, fxm, 0 )
    fxm = f(xm)

    write ( *, * ) ' '
    write ( *, '(i4,'' X'',3f16.8)' ) iterate,  xl,  xm,  xr
    write ( *, '(4x,''FX'',3g16.8)' )          fxl, fxm, fxr
!
!  Prepare for the next step.
!
    if ( fxl * fxm > 0.0E+00 ) then
      xl = xm
      fxl = fxm
    else
      xr = xm
      fxr = fxm
    end if

  end do

  return
end
subroutine brent ( f, fatol, itmax, xb, xc, xatol, xrtol )
!
!*******************************************************************************
!
!! BRENT implements the Brent bisection-based zero finder.
!
!
!  Reference:
!
!    Richard Brent,
!    Algorithms for Minimization without Derivatives,
!    Prentice Hall, 1973.
!
!  Modified:
!
!    20 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real FATOL, an absolute error tolerance for the function
!    value of the root.  If an approximate root X satisfies
!      ABS ( F ( X ) ) <= FATOL, then X will be accepted as the
!    root and the iteration will be terminated.
!
!    Input, integer ITMAX, the maximum number of steps allowed.
!
!    Input, real X1, X2, two points at which the function differs in sign.  
! 
!    Input, real XATOL, XRTOL, absolute and relative error tolerances 
!    for the root.
!
  implicit none
!
  real*8 f
  real*8 d
  real*8 e
  real*8 fatol
  real*8 fxa
  real*8 fxb
  real*8 fxc
  integer iterate
  integer itmax
  real*8 p
  real*8 q
  real*8 r
  real*8 s
  real*8 x_ave
  real*8 x_int
  !real*8 x1
  !real*8 x2
  real*8 xa
  real*8 xb
  real*8 xc
  real*8 xm
  real*8 xatol
  real*8 xrtol
!
!  Initialization.
!
  write ( *, * ) ' '
  write ( *, * ) 'Brent''s Method'
  write ( *, * ) ' '
  write ( *, * ) 'Step      XA            XB             F(XA)         F(XB)'
  write ( *, * ) ' '

  iterate = 0

  xa = xb !x1
  xb = xc !x2
  !call f ( xa, fxa, 0 )
  !call f ( xb, fxb, 0 )
  fxa = f(xa)
  fxb = f(xb)
!
!  Check for change of sign.
!
  if ( fxa * fxb > 0.0E+00 ) then
    write ( *, * ) ' '
    write ( *, * ) 'WARNING:'
    write ( *, * ) ' '
    write ( *, * ) '  This method requires a change of sign '
    write ( *, * ) '  interval, but your function does NOT change '
    write ( *, * ) '  sign over the given interval, so this method '
    write ( *, * ) '  cannot be used.'
    write ( *, * ) ' '
    write ( *, * ) 'SUGGESTION:'
    write ( *, * ) ' '
    write ( *, * ) '  Try a different method;'
    write ( *, * ) '  Try changing X1 or X2;'
    write ( *, * ) '  Perhaps your function F(X) is wrong.'
    return
  end if

  xc = xa
  fxc = fxa
  d = xb - xa
  e = d

  do

    write ( *, '(i4,2x,2g16.8,2g14.6)' ) iterate, xb, xc, fxb, fxc

    iterate = iterate + 1
 
    if ( iterate > itmax ) then
      write ( *, * ) ' '
      write ( *, * ) '  Maximum number of steps taken.'
      exit
    end if
!
!  If necessary, swap data so that |F(XB)| <= |F(XC)|.
!
    if ( abs ( fxc ) < abs ( fxb ) ) then
      call r_swap ( xb, xc )
      call r_swap ( fxb, fxc )
    end if
!
!  XM is the halfwidth of the current change-of-sign interval.
!
    x_int = xc - xb
    x_ave = ( abs ( xc ) + abs ( xb ) ) / 2.0E+00

    xm = 0.5E+00 * ( xc - xb )

    if ( abs ( x_int ) <= xatol ) then
      write ( *, * ) ' '
      write ( *, * ) '  Absolute convergence of the X interval.'
      exit
    end if

    if ( abs ( x_int ) <= xrtol * x_ave ) then
      write ( *, * ) ' '
      write ( *, * ) '  Relative convergence of the X interval.'
      exit
    end if

    if ( abs ( fxb ) <= fatol ) then
      write ( *, * ) ' '
      write ( *, * ) '  Absolute convergence of |F(X)|.'
      exit   
    end if 
!
!  See if a bisection is forced.
!
    if ( abs ( e ) < xatol .or. abs ( fxa ) <= abs ( fxb ) ) then

      d = xm
      e = d
  
    else

      s = fxb / fxa
!
!  Linear interpolation.
!
      if ( xa == xc ) then

        p = 2.0E+00 * xm * s
        q = 1.0E+00 - s
!
!  Inverse quadratic interpolation.
!
      else

        q = fxa / fxc
        r = fxb / fxc
        p = s * ( 2.0E+00 * xm * q * ( q - r ) - ( xb - xa ) * ( r - 1.0E+00 ) )
        q = ( q - 1.0E+00 ) * ( r - 1.0E+00 ) * ( s - 1.0E+00 )

      end if

      if ( p > 0.0E+00 ) then
        q = - q
      else
        p = - p
      end if

      s = e
      e = d

      if ( 2.0E+00 * p >= 3.0E+00 * xm * q - abs ( xatol * q ) .or. &
           p >= abs ( 0.5E+00 * s * q ) ) then
        d = xm
        e = d
      else
        d = p / q
      end if

    end if
!
!  Save in XA, FXA the previous values of XB, FXB.
!
    xa = xb
    fxa = fxb
!
!  Compute the next iterate.
!
    if ( abs ( d ) > xatol ) then
      xb = xb + d
    else if ( xm > 0.0E+00 ) then
      xb = xb + xatol
    else if ( xm <= 0.0E+00 ) then
      xb = xb - xatol
    end if

    !call f ( xb, fxb, 0 )
    fxb = f(xb)
!
!  If the new FXB has the same sign as FXC, then replace XC by XA.
!
    if ( ( fxb > 0.0E+00 .and. fxc > 0.0E+00 ) .or. &
         ( fxb < 0.0E+00 .and. fxc < 0.0E+00 ) ) then
      xc = xa
      fxc = fxa
      d = xb - xa
      e = d
    end if

  end do

  return
end

subroutine muller_r ( f, fatol, itmax, x1, x2, x3, xatol, xrtol )
!
!*******************************************************************************
!
!! MULLER_R carries out Muller's method restricted to real roots.
!
!
!  Modified:
!
!    18 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real FATOL, the absolute error tolerance for F(X).
!
!    Input, integer ITMAX, the maximum number of steps allowed.
!
!    Input, real X1, X2, X3, three points to start the iteration.  
!
!    Input, real XATOL, absolute error tolerance for the root.
!
  implicit none
!
  real*8 f
  double precision a
  double precision b
  double precision c
  double precision discrm
  real*8 fatol
  real*8 fminus
  real*8 fplus
  real*8 fxmid
  real*8 fxnew
  real*8 fxold
  integer iterate
  integer itmax
  real*8 x_ave
  real*8 x_inc
  real*8 x1
  real*8 x2
  real*8 x3
  real*8 xatol
  real*8 xlast
  real*8 xmid
  real*8 xminus
  real*8 xnew
  real*8 xold
  real*8 xplus
  real*8 xrtol
!
  xnew = x1
  xmid = x2
  xold = x3
  !call f ( xnew, fxnew, 0 )
  !call f ( xmid, fxmid, 0 )
  !call f ( xold, fxold, 0 )
  fxnew = f(xnew)
  fxmid = f(xmid)
  fxold = f(xold)
  
  write ( *, * ) ' '
  write ( *, * ) ' '
  write ( *, * ) 'Muller''s method (real root version)'
  write ( *, * ) ' '
  write ( *, * ) 'Note that the existence of a nearby complex root'
  write ( *, * ) 'is suggested when the value of the discriminant'
  write ( *, * ) 'is negative.'
  write ( *, * ) ' '
  write ( *, * ) ' '
  write ( *, * ) 'Iteration          x                 fx            disc'
  write ( *, * ) ' '

  iterate = 0
  discrm = 0.0E+00
  write ( *, '(1x,i6,f20.10,f20.10,f20.10)' ) iterate, xnew, fxnew

  if ( abs ( fxnew ) < fatol ) then
    write ( *, * ) ' '
    write ( *, * ) '|F(X)| is below the tolerance.'
    return
  end if

  do

    if ( abs ( fxnew ) >= abs ( fxmid ) ) then
      call r_swap ( xnew, xmid )
      call r_swap ( fxnew, fxmid )
    end if

    xlast = xnew
    iterate = iterate + 1

    if ( iterate > itmax ) then
      write ( *, * ) ' '
      write ( *, * ) '  Maximum number of steps taken.'
      exit
    end if

    a = dble ( ( xmid - xnew ) * ( fxold - fxnew ) &
             - ( xold - xnew ) * ( fxmid - fxnew ) )

    b = dble ( ( xold - xnew )**2 * ( fxmid - fxnew ) &
             - ( xmid - xnew )**2 * ( fxold - fxnew ) )

    c = dble ( ( xold - xnew ) * ( xmid - xnew ) * ( xold - xmid ) * fxnew )

    xold = xmid
    xmid = xnew
!
!  Apply the quadratic formula to get roots XPLUS and XMINUS.
!
    discrm = b**2 - 4.0E+00 * a * c

    if ( discrm <= 0.0E+00 ) then
      discrm = 0.0E+00
    end if

    if ( a == 0.0E+00 ) then
      write ( *, * ) ' '
      write ( *, * ) 'The algorithm has broken down.'
      write ( *, * ) 'The quadratic coefficient A is zero.'
      exit
    end if

    xplus = xnew + sngl ( ( - b + sqrt ( discrm ) ) / ( 2.0E+00 * a ) )

    !call f ( xplus, fplus, 0 )
    fplus = f(xplus)

    xminus = xnew + sngl ( ( - b - sqrt ( discrm ) ) / ( 2.0E+00 * a ) )

    !call f ( xminus, fminus, 0 )
    fminus = f(xminus)

!
!  Take whichever of the two quadratic roots is closest to a root of 
!  the function.
!
    if ( abs ( fminus ) < abs ( fplus ) ) then
      xnew = xminus
    else
      xnew = xplus
    end if

    fxold = fxmid
    fxmid = fxnew
    !call f ( xnew, fxnew, 0 )
    fxnew = f(xnew)
    
    write ( *, '(1x,i6,f20.10,f20.10,f20.10)' ) iterate, xnew, fxnew, discrm
!
!  Check for convergence.
!
    x_ave = ( abs ( xnew ) + abs ( xmid ) + abs ( xold ) ) / 3.0E+00
    x_inc = xnew - xmid
    
    if ( abs ( x_inc ) <= xatol ) then
      write ( *, * ) ' '
      write ( *, * ) '  Absolute convergence of the X increment.'
      exit
    end if

    if ( abs ( x_inc ) <= xrtol * x_ave ) then
      write ( *, * ) ' '
      write ( *, * ) '  Relative convergence of the X increment.'
      exit
    end if
    
    if ( abs ( fxnew ) <= fatol ) then
      write ( *, * ) ' '
      write ( *, * ) '  Absolute convergence of |F(X)|.'
      exit
    end if   
     
  end do

  return
end


subroutine newton ( f, df, fatol, itmax, x1, xatol, xrtol )
!
!*******************************************************************************
!
!! NEWTON carries out Newton's method.
!
!
!  Modified:
!
!    20 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real FATOL, the absolute error tolerance for F(X).
!
!    Input, integer ITMAX, the maximum number of steps allowed.
!
!    Input, real X1, a point where the iteration begins.  
!
!    Input, real XATOL, absolute error tolerance for the root.
!
  implicit none
!
  integer, parameter :: MAX_ROOT = 20
!
  real*8 f
  real*8 df
  real*8 fatol
  real*8 fprime
  real*8 fx
  integer i
  integer iroot
  character isay
  integer iterat
  integer itmax
  real*8 roots(MAX_ROOT)
  real*8 x_ave
  real*8 x_inc
  real*8 x1
  real*8 xatol
  real*8 xnew
  real*8 xold
  real*8 xrtol
!
  roots(1:MAX_ROOT) = 0.0E+00

  iroot = 0

10 continue

  iroot = iroot + 1
  if ( iroot > MAX_ROOT ) then
    return
  end if

  write ( *, * ) ' '
  write ( *, * ) ' '
  write ( *, * ) ' Newton''s method'
  write ( *, * ) ' '
  write ( *, * ) ' '
  write ( *, * ) 'Iteration            X                    FX'
  write ( *, * ) ' '

  iterat = 0
  xnew = x1
  !call f ( xnew, fx, 0 )
  fx = f(xnew)
  write ( *, '(i6,f20.10,f20.10)' ) iterat, xnew, fx

  if ( abs ( fx ) <= fatol ) then
    write ( *, * ) ' '
    write ( *, * ) '|F(X)| is below the tolerance.'
    return
  end if

  do

    iterat = iterat + 1

    if ( iterat > itmax ) then
      write ( *, * ) ' '
      write ( *, * ) '  Maximum number of steps taken.'
      return
    end if

    xold = xnew
    !call f ( xold, fx, 0 )
    fx = f(xold)
    !call f ( xold, fprime, 1 )
    fprime = df(xold)

    do i = 1, iroot - 1

      if ( xold - roots(i) == 0.0E+00 ) then
        write ( *, * ) ' '
        write ( *, * ) 'The algorithm has broken down.'
        write ( *, * ) 'A zero divisor was computed.'
        exit
      end if

      fprime = fprime - fx / ( xold - roots(i) )

    end do

    if ( fprime == 0.0E+00 ) then
      write ( *, * ) ' '
      write ( *, * ) 'The algorithm has broken down.'
      write ( *, * ) 'A zero divisor was computed.'
      exit
    end if

    x_inc = - fx / fprime
    xnew = xold + x_inc

    !call f ( xnew, fx, 0 )
    fx = f(xnew)

    write ( *, '(i6,f20.10,f20.10)' ) iterat, xnew, fx
!
!  Check for acceptance of function value or interval
!
    x_ave = abs ( xnew ) + abs ( xold )
    
    if ( abs ( x_inc ) <= xatol ) then
      write ( *, * ) ' '
      write ( *, * ) '  Absolute convergence of the X increment.'
      exit
    end if

    if ( abs ( x_inc ) <= xrtol * x_ave ) then
      write ( *, * ) ' '
      write ( *, * ) '  Relative convergence of the X increment.'
      exit
    end if
    
    if ( abs ( fx ) <= fatol ) then
      write ( *, * ) ' '
      write ( *, * ) '  Absolute convergence of |F(X)|.'
      exit
    end if    

  end do

  return
end

subroutine secant ( f, fatol, itmax, x1, x2, xatol, xrtol )
!
!*******************************************************************************
!
!! SECANT carries out the secant algorithm.
!
!
!  Modified:
!
!    18 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real FATOL, the absolute error tolerance for F(X).
!
!    Input, integer ITMAX, the maximum number of steps allowed.
!
!    Input, real X1, X2, two points at which the iteration begins.  
!
!    Input, real XATOL, absolute error tolerance for the root.
!
  implicit none
!
  real*8 f
  real*8 fatol
  real*8 fx
  real*8 fxnew
  real*8 fxold
  integer iterate
  integer itmax
  real*8 x_ave
  real*8 x_inc
  real*8 x1
  real*8 x2
  real*8 xatol
  real*8 x
  real*8 xnew
  real*8 xold
  real*8 xrtol
!
  write ( *, * ) ' '
  write ( *, * ) ' '
  write ( *, * ) 'The secant method:'
  write ( *, * ) ' '
  write ( *, * ) ' '
  write ( *, * ) 'Step       X          FX'
  write ( *, * ) ' '

  iterate = -1
  xold = x1
  !call f ( xold, fxold, 0 )
  fxold = f(xold)
  write ( *, '(i4,2g16.8)' ) iterate, xold, fxold

  if ( abs ( fxold ) <= fatol ) then
    write ( *, * ) ' '
    write ( *, * ) '|F(X)| is below the tolerance.'
    return
  end if

  iterate = 0
  x = x2
  !call f ( x, fx, 0 )
  fx = f(x)
  write ( *, '(i4,2g16.8)' ) iterate, x, fx

  if ( abs ( fx ) <= fatol ) then
    write ( *, * ) ' '
    write ( *, * ) '  Absolute convergence of |F(X)|.'
    return
  end if

  do

    iterate = iterate + 1

    if ( iterate > itmax ) then
      write ( *, * ) ' '
      write ( *, * ) '  Maximum number of steps taken.'
      exit
    end if

    x_inc = x - xold
    x_ave = ( abs ( x ) + abs ( xold ) ) / 2.0E+00
    
    if ( abs ( x_inc ) <= xatol ) then
      write ( *, * ) ' '
      write ( *, * ) '  Absolute convergence of the X increment.'
      exit
    end if

    if ( abs ( x_inc ) <= xrtol * x_ave ) then
      write ( *, * ) ' '
      write ( *, * ) '  Relative convergence of the X increment.'
      exit
    end if
    
    if ( abs ( fx ) <= fatol ) then
      write ( *, * ) ' '
      write ( *, * ) '  Absolute convergence of |F(X)|.'
      exit
    end if
!
!  Compute the next iterate.
!
    xnew = ( fxold * x - fx * xold ) / ( fxold - fx )
    !call f ( xnew, fxnew, 0 )
    fxnew = f(xnew)

    write ( *, '(i4,2g16.8)' ) iterate, xnew, fxnew
!
!  Prepare for the next step.
!
    xold = x
    fxold = fx

    x = xnew
    fx = fxnew

    end do
    
  return
end

subroutine r_swap ( x, y )
!
!*******************************************************************************
!
!! R_SWAP switches two real values.
!
!
!  Modified:
!
!    30 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real X, Y.  On output, the values of X and
!    Y have been interchanged.
!
  implicit none
!
  real*8 x
  real*8 y
  real*8 z
!
  z = x
  x = y
  y = z

  return
end
