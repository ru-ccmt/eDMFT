/*
 ************************************************************************
 *
 *			  Numerical Math Package
 *
 *			    Brent's root finder
 *	       obtains a zero of a function of one variable
 *
 * Synopsis
 *	double zeroin(ax,bx,f,tol=EPSILON)
 *	const double ax 		The root is to be sought within
 *	const double bx  		the interval [ax,bx]
 *	UnivariateFunctor& f		The function under consideration
 *	const double tol		Acceptable tolerance for the root
 *					position. It is an optional parameter
 *					with default value DBL_EPSILON
 *
 *	Zeroin returns an approximate location for the root with accuracy
 *	4*DBL_EPSILON*abs(x) + tol
 *
 * Algorithm
 *	G.Forsythe, M.Malcolm, C.Moler, Computer methods for mathematical
 *	computations. M., Mir, 1980, p.180 of the Russian edition
 *
 * The function makes use of a bissection procedure combined with
 * a linear or quadratic inverse interpolation.
 * At each step the code operates three abscissae - a, b, and c:
 *	b - the last and the best approximation to the root
 *	a - the last but one approximation
 *	c - the last but one or even an earlier approximation such that
 *		1) |f(b)| <= |f(c)|
 *		2) f(b) and f(c) have opposite signs, i.e. b and c encompass
 *		   the root
 * Given these abscissae, the code computes two new approximations, one by the 
 * bissection procedure and the other one from interpolation (if a,b, and c
 * are all different the quadratic interpolation is used, linear otherwise).
 * If the approximation obtained by the interpolation looks
 * reasonable (i.e. falls within the current interval [b,c], not too close
 * to the end points of the interval), the point is accepted as a new
 * approximation to the root. Otherwise, the result of the bissection is used.
 * Therefore, the range of uncertainty is guaranteed to tighten at 
 * least by a factor of 1.6
 *
 *
 ************************************************************************
 */
#ifndef _ZEROIN
#define _ZEROIN
#define DBL_EPSILON_		2.22045e-16

template <class functor>
inline double zeroin(				// An estimate to the root
	const double ax,		// Specify the interval the root
	const double bx,		// to be sought in
	functor& f,		        // Function under investigation
	const double tol)		// Acceptable tolerance
{
  Assert( tol > 0, "Tolerance must be positive");
  Assert( bx > ax, 
	 "Left end point of the interval should be strictly less than the "
	 "right one" );

  double b = bx;		// the last and the best approx to the root
  double fb = f(b);
  double a = ax;		// the last but one approximation
  double fa = f(a);
  double c = a;			// the last but one or even an earlier approx
  double fc = fa;		// see the condition above

  for(;;)		// Main iteration loop
  {
    const double prev_step = b-a;	// Step from the previous iteration
   
    if( std::abs(fc) < std::abs(fb) )
    {                         		// Swap data so that b would be the
      a = b;  b = c;  c = a;          	// best approximation found so far
      fa=fb;  fb=fc;  fc=fa;
    }
					// Estimate the effective tolerance
    const double tol_act = 2*DBL_EPSILON_*std::abs(b) + tol/2;
    double new_step = (c-b)/2;		// Bissection step for this iteration

    if( std::abs(new_step) <= tol_act || fb == 0 )
      return b;				// Acceptable approximation is found

    			// Figuring out if the interpolation can be tried
    if( std::abs(prev_step) >= tol_act	// If prev_step was large enough
	&& std::abs(fa) > std::abs(fb) )		// and was in true direction,
    {					// Interpolatiom may be tried

      double p;      			// Interpolation step is calcu-
      double q;      			// lated in the form p/q; divi-
  					// sion operations is delayed
 					// until the last moment
      const double cb = c-b;

      if( a==c )			// If we've got only two distinct
      {					// points linear interpolation
	/*register*/ const double t1 = fb/fa;	// can only be applied
	p = cb*t1;
	q = 1.0 - t1;
      }
      else				// Quadratic inverse interpolation
      {
	/*register*/ const double t1=fb/fc, t2=fb/fa;
	q = fa/fc;
	p = t2 * ( cb*q*(q-t1) - (b-a)*(t1-1.0) );
	q = (q-1.0) * (t1-1.0) * (t2-1.0);
      }

      if( p > 0 )			// Formulas above computed new_step
	q = -q;				// = p/q with the wrong sign (on purpose).
      else				// Correct this, but in such a way so
	p = -p;				// that p would be positive
      
      if( 2*p < (1.5*cb*q-std::abs(tol_act*q))	// If b+p/q falls in [b,c]
	 && 2*p < std::abs(prev_step*q) )	// and isn't too large
	new_step = p/q;			// it is accepted
					// If p/q is too large then the
					// bissection procedure can
					// reduce [b,c] to a larger
					// extent
    }

    if( std::abs(new_step) < tol_act )	// Adjust the step to be not less
      new_step =  new_step > 0 ?	// than the tolerance
	 tol_act : -tol_act;

    a = b;  fa = fb;			// Save the previous approximation
    b += new_step;  fb = f(b);		// Do step to a new approximation
    if( (fb > 0 && fc > 0) || (fb < 0 && fc < 0) )
    {                 			// Adjust c for it to have the sign
      c = a;  fc = fa;                  // opposite to that of b
    }
  }
}
#endif //_ZEROIN
