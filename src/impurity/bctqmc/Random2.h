/*
  Random.h
  quasi- and pseudo-random-number generation
  K. Haule, 2015
*/
#include <cmath>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <map>

#ifdef _PCG

#include <iostream>
#include <iomanip>
#include <string>
#include <random>
#include <cmath>
#include "pcg/pcg_random.hpp"
using namespace std;
using pcg_extras::operator<<;


class RanPCG{
public:
  pcg64 rng;
  std::uniform_real_distribution<double> uniform_dist;

  // explicit stream = 1 (odd)
  explicit RanPCG(uint64_t initstate=0u, uint64_t initseq=1u)
    : rng(initstate, (initseq|1u)), uniform_dist(0.0, 1.0) {}

  void set_iseed(uint64_t initstate, uint64_t initseq=1u) {
    rng.seed(initstate, (initseq|1u));
  }

  void set_nondeterministic_seed() {
    pcg_extras::seed_seq_from<std::random_device> sseq;
    rng.seed(sseq);  // pcg64 supports seed_seq
  }

  // Reproducible per-rank seeding
  uint64_t seed_repro(uint64_t base_seed, uint64_t rank){
    uint64_t initstate, initseq;
    if (base_seed == 0) {
      pcg_extras::seed_seq_from<std::random_device> sseq;
      sseq.generate(&initstate, &initstate+1);
      sseq.generate(&initseq,   &initseq+1);
    } else {
      initstate = splitmix64(base_seed ^ 0xD2B74407B1CE6E93ULL ^ rank);
      initseq   = splitmix64(base_seed ^ 0xCA5A826395121157ULL ^ (rank*0x9E3779B97F4A7C15ULL));
    }
    rng.seed(initstate, (initseq|1u));
    return initstate;
  }

  double operator()() { return uniform_dist(rng); }

  int Get(double* x, int N){
    for (int i=0; i<N; i++) x[i] = uniform_dist(rng);
    return 0;
  }

  // Distance helpers (must match type + stream)
  std::uint64_t DistanceFromZero(){
    pcg64 rng0(0u, 1u);
    return rng - rng0;                 // distance between two pcg64 with same stream
  }
  void set_iseed_DistanceFromZero(std::uint64_t dist){
    rng.seed(0u, 1u);                  // SAME (state, stream) as above
    rng.advance(dist);
  }

private:
  static inline uint64_t splitmix64(uint64_t x){
    x += 0x9E3779B97F4A7C15ULL;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
  }
};
/*
class RanPCG{
public:
  //pcg64_fast rng;
  pcg64 rng;
  std::uniform_real_distribution<double> uniform_dist; 
  RanPCG(unsigned int idum) : rng(idum), uniform_dist(0, 1) {};
  void set_iseed(unsigned int idum) {rng.seed(idum);}
  void set_nondeterministic_seed() {rng.seed(pcg_extras::seed_seq_from<std::random_device>());}
  // ---------- New: reproducible per-rank seeding ----------
  // Reproducible per-rank seeding (keeps old call sites simple)
  uint64_t set_seed_repro(uint64_t base_seed, uint64_t rank){
    uint64_t initstate, initseq;
    if (base_seed == 0) {
      // nondeterministic base
      pcg_extras::seed_seq_from<std::random_device> sseq;
      sseq.generate(&initstate, &initstate+1);
      sseq.generate(&initseq,   &initseq+1);
    } else {
      // derive per-rank state/stream from base_seed
      initstate = splitmix64(base_seed ^ 0xD2B74407B1CE6E93ULL ^ rank);
      initseq   = splitmix64(base_seed ^ 0xCA5A826395121157ULL ^ (rank*0x9E3779B97F4A7C15ULL));
    }
    rng.seed(initstate, (initseq|1u));          // ensure odd stream
    return initstate;
  }
  double operator()() { return uniform_dist(rng);}
  int Get(double* x, int N){
    for (int i=0; i<N; i++) x[i]=uniform_dist(rng);
    return 0;
  }
  size_t DistanceFromZero(){
    pcg64_fast rng0(0u);
    return rng - rng0;
  }
  void set_iseed_DistanceFromZero(size_t dist){
    rng.seed(0u);
    rng.advance(dist);
  }
private:
  static inline uint64_t splitmix64(uint64_t x){
    x += 0x9E3779B97F4A7C15ULL;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
  }
};
*/
//
//// Good 64-bit mixer (SplitMix64) to decorrelate inputs
//static inline uint64_t splitmix64(uint64_t x) {
//    uint64_t z = x + 0x9E3779B97F4A7C15ull;
//    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
//    z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
//    return z ^ (z >> 31);
//}
//
//// Return a reproducible "base" when nonzero, else time-based
//inline uint64_t resolve_base_seed(uint64_t base_seed_param) {
//    if (base_seed_param != 0) return base_seed_param;
//    auto now = std::chrono::high_resolution_clock::now().time_since_epoch().count();
//    // Mix the clock a bit to avoid weak low bits on some platforms
//    return splitmix64(static_cast<uint64_t>(now));
//}
//
//// Build a per-rank PCG32 with independent stream
//inline pcg32 make_pcg32(uint64_t base_seed_param, uint64_t rank) {
//    const uint64_t base = resolve_base_seed(base_seed_param);
//
//    // Derive initstate (state seed) and initseq (stream selector) from base & rank
//    // initseq must be ODD for PCG
//    uint64_t initstate = splitmix64(base ^ 0xD2B74407B1CE6E93ull ^ rank);
//    uint64_t initseq   = (splitmix64(base ^ 0xCA5A826395121157ull ^ (rank * 0x9E3779B97F4A7C15ull)) << 1) | 1ull;
//
//    // Two-argument ctor seeds and initializes stream selector
//    return pcg32(initstate, initseq);
//}
//
#endif


/*
possible use:
  // SOBOL
  int ndim=10, maxeval=100000;
  Sobol rand(maxeval, ndim);
  vector<double> x(ndim);
  rand.Get(x);

  // Mersenne
  Mersenne rand(iseed);
  vector<double> x(ndim);
  rand.Get(x);

  // Ranlux
  Ranlux rand(iseed);
  vector<double> x(ndim);
  rand.Get(x);

  // MKL random numbers
  RanMKL rand(iseed, type=0,1,2,3,4);
  double x[ndim];
  double x = rand();
  //or
  rand.Get(x,ndim);

  // GSL random numbers
  RanGSL rand(iseed);
  double x = rand();
  // or
  rand.Get(x,ndim);
*/



#define SOBOL_MINDIM 1
#define SOBOL_MAXDIM 40
/*
 Sobol quasi-random-number generator
 adapted from ACM TOMS algorithm 659
*/
/*
class Sobol{
  double norm;
  int v[SOBOL_MAXDIM][30], prev[SOBOL_MAXDIM];
  int seq;
public:
  Sobol(int maxeval, int ndim);
  template<class container>
  void Get(container& x, int size);
  template<class container>
  void Get(container& x) {Get(x,x.size());}
  void Skip(int n, int ndim);
};
*/
/*
 Mersenne Twister pseudo-random-number generator
 adapted from T. Nishimura's and M. Matsumoto's C code at
 http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
*/
#define MERSENNE_M 397
#define MERSENNE_N 624
class Mersenne{
  unsigned int state[MERSENNE_N];
  int next;
  unsigned int Twist(unsigned int a, unsigned int b);
  void Reload();
  unsigned int Int(unsigned int s);
public:
  Mersenne(unsigned int seed);
  template<class container>
  void Get(container& x, int size);
  template<class container>
  void Get(container& x) {Get(x,x.size());}
  void Skip(int n, int ndim);
};

/*
 Ranlux subtract-and-borrow random-number generator 
 proposed by Marsaglia and Zaman, implemented by F. James with 
 the name RCARRY in 1991, and later improved by Martin Luescher 
 in 1993 to produce "Luxury Pseudorandom Numbers".
 Adapted from the CERNlib Fortran 77 code by F. James, 1993.
 
 The available luxury levels are:
 
 level 0  (p = 24): equivalent to the original RCARRY of Marsaglia
          and Zaman, very long period, but fails many tests.
 level 1  (p = 48): considerable improvement in quality over level 0,
          now passes the gap test, but still fails spectral test.
 level 2  (p = 97): passes all known tests, but theoretically still
          defective.
 level 3  (p = 223): DEFAULT VALUE.  Any theoretically possible
          correlations have very small chance of being observed.
 level 4  (p = 389): highest possible luxury, all 24 bits chaotic.
*/
class Ranlux{
  int n24, i24, j24, nskip;
  int carry, state[24];
  int Int(int n);
public:
  Ranlux(int seed, int level=4);
  void Skip(int n, int ndim);
  template<class container>
  void Get(container& x, int size);
  template<class container>
  void Get(container& x) {Get(x,x.size());}
};


////////////////// Implementation of Sobol /////////////////////
/*
Sobol::Sobol(int maxeval,int ndim) // t->
{
  static int ini[9*40] = {
      3,   1,   0,   0,   0,   0,   0,   0,   0,
      7,   1,   1,   0,   0,   0,   0,   0,   0,
     11,   1,   3,   7,   0,   0,   0,   0,   0,
     13,   1,   1,   5,   0,   0,   0,   0,   0,
     19,   1,   3,   1,   1,   0,   0,   0,   0,
     25,   1,   1,   3,   7,   0,   0,   0,   0,
     37,   1,   3,   3,   9,   9,   0,   0,   0,
     59,   1,   3,   7,  13,   3,   0,   0,   0,
     47,   1,   1,   5,  11,  27,   0,   0,   0,
     61,   1,   3,   5,   1,  15,   0,   0,   0,
     55,   1,   1,   7,   3,  29,   0,   0,   0,
     41,   1,   3,   7,   7,  21,   0,   0,   0,
     67,   1,   1,   1,   9,  23,  37,   0,   0,
     97,   1,   3,   3,   5,  19,  33,   0,   0,
     91,   1,   1,   3,  13,  11,   7,   0,   0,
    109,   1,   1,   7,  13,  25,   5,   0,   0,
    103,   1,   3,   5,  11,   7,  11,   0,   0,
    115,   1,   1,   1,   3,  13,  39,   0,   0,
    131,   1,   3,   1,  15,  17,  63,  13,   0,
    193,   1,   1,   5,   5,   1,  27,  33,   0,
    137,   1,   3,   3,   3,  25,  17, 115,   0,
    145,   1,   1,   3,  15,  29,  15,  41,   0,
    143,   1,   3,   1,   7,   3,  23,  79,   0,
    241,   1,   3,   7,   9,  31,  29,  17,   0,
    157,   1,   1,   5,  13,  11,   3,  29,   0,
    185,   1,   3,   1,   9,   5,  21, 119,   0,
    167,   1,   1,   3,   1,  23,  13,  75,   0,
    229,   1,   3,   3,  11,  27,  31,  73,   0,
    171,   1,   1,   7,   7,  19,  25, 105,   0,
    213,   1,   3,   5,   5,  21,   9,   7,   0,
    191,   1,   1,   1,  15,   5,  49,  59,   0,
    253,   1,   1,   1,   1,   1,  33,  65,   0,
    203,   1,   3,   5,  15,  17,  19,  21,   0,
    211,   1,   1,   7,  11,  13,  29,   3,   0,
    239,   1,   3,   7,   5,   7,  11, 113,   0,
    247,   1,   1,   5,   3,  15,  19,  61,   0,
    285,   1,   3,   1,   1,   9,  27,  89,   7,
    369,   1,   1,   3,   7,  31,  15,  45,  23,
    299,   1,   3,   3,   9,   9,  25, 107,  39 };

  int dim, bit, nbits;
  int *pini = ini, max;

  for( nbits = 0, max = maxeval; max; max >>= 1 ) ++nbits;
  this->norm = ldexp(.5, -nbits);

  for( bit = 0; bit <= nbits; ++bit )
    this->v[0][bit] = (int)1 << (nbits - bit);

  for( dim = 1; dim < ndim; ++dim ) {
    int *pv = this->v[dim], *pvv = pv;
    int powers = *pini++, j;
    int inibits = -1, bit;
    for( j = powers; j; j >>= 1 ) ++inibits;

    memcpy(pv, pini, inibits*sizeof *pini);
    pini += 8;

    for( bit = inibits; bit <= nbits; ++bit ) {
      int newv = *pvv, j = powers;
      int b;
      for( b = 0; b < inibits; ++b ) {
        if( j & 1 ) newv ^= pvv[b] << (inibits - b);
        j >>= 1;
      }
      pvv[inibits] = newv;
      ++pvv;
    }

    for( bit = 0; bit < nbits; ++bit )
      pv[bit] <<= nbits - bit;
  }

  this->seq = 0;
  memset(this->prev, 0, ndim*sizeof(*this->prev));
}
template<class container>
void Sobol::Get(container& x, int size)
{
  int seq = this->seq++;
  int zerobit = 0;

  while( seq & 1 ) {
    ++zerobit;
    seq >>= 1;
  }
  for(int dim = 0; dim < size; ++dim ) {
    this->prev[dim] ^= this->v[dim][zerobit];
    x[dim] = this->prev[dim]*this->norm;
  }
}
void Sobol::Skip(int n, int ndim)
{
  while( n-- ) {
    int seq = this->seq++;
    int zerobit = 0, dim;

    while( seq & 1 ) {
      ++zerobit;
      seq >>= 1;
    }

    for( dim = 0; dim < ndim; ++dim )
      this->prev[dim] ^= this->v[dim][zerobit];
  }
}
*/

////////////////// Implementation of Mersenne /////////////////////

inline unsigned int Mersenne::Twist(unsigned int a, unsigned int b)
{
  unsigned int mixbits = (a & 0x80000000) | (b & 0x7fffffff);
  unsigned int matrixA = (-(b & 1)) & 0x9908b0df;
  return (mixbits >> 1) ^ matrixA;
}
inline void Mersenne::Reload()
{
  unsigned int *s = state;
  for(int j = MERSENNE_N - MERSENNE_M + 1; --j; ++s )
    *s = s[MERSENNE_M] ^ Twist(s[0], s[1]);
  for(int j = MERSENNE_M; --j; ++s )
    *s = s[MERSENNE_M - MERSENNE_N] ^ Twist(s[0], s[1]);
  *s = s[MERSENNE_M - MERSENNE_N] ^ Twist(s[0], state[0]);
}

Mersenne::Mersenne(unsigned int seed)
{
  unsigned int *next = this->state;
  for(int j = 1; j <= MERSENNE_N; ++j ) {
    *next++ = seed;
    seed = 0x6c078965*(seed ^ (seed >> 30)) + j;
    /* see Knuth TAOCP Vol 2, 3rd Ed, p. 106 for multiplier */
  }
  Reload();
  this->next = 0;
}


/* 32 or 53 random bits */
#define RANDOM_BITS 32
inline unsigned int Mersenne::Int(unsigned int s)
{
  s ^= s >> 11;
  s ^= (s << 7) & 0x9d2c5680;
  s ^= (s << 15) & 0xefc60000;
  return s ^ (s >> 18);
}


template<class container>
void Mersenne::Get(container& x, int size)
{
  int next = this->next;
  for(int dim = 0; dim < size; ++dim ) {
#if RANDOM_BITS == 53
    unsigned int a, b;
#endif
    
    if( next >= MERSENNE_N ) {
      Reload();
      next = 0;
    }

#if RANDOM_BITS == 53
    a = Int(this->state[next++]) >> 5;
    b = Int(this->state[next++]) >> 6;
    x[dim] = (67108864.*a + b)/9007199254740992.;
#else
    x[dim] = Int(this->state[next++])/4294967296.;
#endif
  }
  this->next = next;
}

void Mersenne::Skip(int n, int ndim)
{
#if RANDOM_BITS == 53
  n = 2*n*ndim + this->next;
#else
  n = n*ndim + this->next;
#endif
  this->next = n % MERSENNE_N;
  n /= MERSENNE_N;
  while( n-- ) Reload();
}

////////////////// Implementation of Ranlux /////////////////////
/* (a < 0) ? -1 : 0 */
#define NegQ(a) ((a) >> (sizeof(a)*8 - 1))
/* (a == 0) ? 0 : -1 */
#define TrueQ(a) NegQ((a) | (-a))

inline int Ranlux::Int(int n)
{
  int s = 0;
  while( n-- ) {
    s = this->state[this->j24] - this->state[this->i24] + this->carry;
    s += (this->carry = NegQ(s)) & (1 << 24);
    this->state[this->i24] = s;
    --this->i24;
    this->i24 += NegQ(this->i24) & 24;
    --this->j24;
    this->j24 += NegQ(this->j24) & 24;
  }
  return s;
}

template<class container>
void Ranlux::Get(container& x, int size)
{
  /* The Generator proper: "Subtract-with-borrow",
     as proposed by Marsaglia and Zaman, FSU, March 1989 */
  for(int dim = 0; dim < size; ++dim ) {
    const int nskip = (--this->n24 >= 0) ? 0 :
      (this->n24 = 24, this->nskip);
    const int s = Int(1 + nskip);
    x[dim] = ldexp(s, -24);
    /* small numbers (with less than 12 significant bits) are "padded" */
    if( s < (1 << 12) )
      x[dim] += ldexp(this->state[this->j24], -48);
  }
}

void Ranlux::Skip(int n, int ndim)
{
  Int(n + this->nskip*(n/24));
  this->n24 = 24 - n % 24;
}

#define Elements(x) (sizeof(x)/sizeof(*x))

Ranlux::Ranlux(int seed, int level)
{
  const int skip[] = {24, 48, 97, 223, 389,
		      223, 223, 223, 223, 223, 223, 223, 223, 223, 223,
		      223, 223, 223, 223, 223, 223, 223, 223, 223, 223};
  
  if( level < Elements(skip) ) level = skip[level];
  this->nskip = level - 24;

  this->i24 = 23;
  this->j24 = 9;
  this->n24 = 24;

  for(int i = 0; i < 24; ++i ) {
    const int k = seed/53668;
    seed = 40014*(seed - k*53668) - k*12211;
    seed += NegQ(seed) & 2147483563;
    this->state[i] = seed & ((1 << 24) - 1);
  }
  this->carry = ~TrueQ(this->state[23]) & (1 << 24);
}


#ifdef _MKL
#include <mkl_vsl.h>

class RanMKL{
  VSLStreamStatePtr stream;
public:
  RanMKL(int seed=777, int type=0){
    switch (type){
    case 0:
      vslNewStream( &stream, VSL_BRNG_MT19937, seed ); // A Mersenne Twister pseudorandom number generator.
      break;
    case 1:
      vslNewStream( &stream, VSL_BRNG_MT2203, seed );  // A set of 6024 Mersenne Twister pseudorandom number generators.
      break;
    case 2:
      vslNewStream( &stream, VSL_BRNG_SFMT19937, seed ); // A SIMD-oriented Fast Mersenne Twister pseudorandom number generator.
      break;
    case 3:
      vslNewStream( &stream, VSL_BRNG_SOBOL, seed ); // A 32-bit Gray code-based generator producing low-discrepancy sequences for dimensions 1 ≤ s ≤ 40; user-defined dimensions are also available.
      break;
    case 4:
      vslNewStream( &stream, VSL_BRNG_NIEDERR, seed); // A 32-bit Gray code-based generator producing low-discrepancy sequences for dimensions 1 ≤ s ≤ 318; user-defined dimensions are also available.
      break;
    default:
      std::cout<<"Not supported random number generator"<<std::endl;
      exit(1);
    }
  };
  double operator()(){
    //vdRngUniform( method, stream, n, r, a, b )
    double r;
    vdRngUniform( VSL_RNG_METHOD_UNIFORM_STD, stream, 1, &r, 0.0, 1.0 );
    return r;
  }
  int Get(double* r, int N){
    return vdRngUniform( VSL_RNG_METHOD_UNIFORM_STD, stream, N, r, 0.0, 1.0 );
  }
  ~RanMKL(){
    vslDeleteStream( &stream );
  }
};
#endif

#ifdef _GSL
#include <gsl/gsl_rng.h>

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
  void set_iseed(int idum)
  {gsl_rng_set(r, idum);}
  double operator()()
  { return gsl_rng_uniform (r);}
  int Get(double* x, int N){
    for (int i=0; i<N; i++) x[i]=gsl_rng_uniform (r);
    return 0;
  }
  ~RanGSL()
  {gsl_rng_free (r);}
};
#endif

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



