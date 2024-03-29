// @Copyright 2007 Kristjan Haule
// 
#include <cstdio>
#include <cstdlib>
#include <gsl/gsl_rng.h>
#include <map>

class RanGSL{
  const gsl_rng_type *T;
  gsl_rng *r;
public:
  RanGSL(int idum)
  {
    gsl_rng_env_setup();
    //    T = gsl_rng_default;
    //    T = gsl_rng_taus;
    T = gsl_rng_ranlux389;
    r = gsl_rng_alloc (T);
    gsl_rng_set(r, idum);
  }
  double operator()()
  { return gsl_rng_uniform (r);}
  ~RanGSL()
  {gsl_rng_free (r);}
};

class dRand48{
public:
  dRand48(int idum){srand48(idum);}
  double operator()(){
    return drand48();
  }
};

class Ran0{
  long idum;
public:
  Ran0(long idum_): idum(idum_){};
  double operator()(){
    // Minimal random number generator of Park and Miller. Returns a uniform random deviate
    // between 0.0 and 1.0. Set or reset idum to any integer value (except the unlikely value MASK)
    // to initialize the sequence; idum must not be altered between calls for successive deviates in
    // a sequence.
    static const int IA=16807;
    static const int IM=2147483647;
    static const double AM=(1.0/IM);
    static const int IQ=127773;
    static const int IR=2836;
    static const int MASK=123459876;
    static long idum;
    long k;
    float ans;
    idum ^= MASK; //XORing with MASK allows use of zero and other simple bit patterns for idum.
    k=idum/IQ; 
    idum = IA*(idum-k*IQ)-IR*k; // Compute idum=(IA*idum) % IM without over
    if (idum < 0) idum += IM; // flows by Schrages method.
    ans=AM*idum; // Convert idum to a floating result.
    idum ^= MASK; // Unmask before return.
    return ans;
  }
};


