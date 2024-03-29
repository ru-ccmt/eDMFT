#ifndef MY_UTIL
#define MY_UTIL

#include <vector>
#include <stdexcept>
#include <blitz/array.h>

namespace bl = blitz;
//using namespace blitz;

inline long double fermi(double x)
{
  if (fabs(x)<300)
    return 1/(exp(x)+1.0);
  else if ( x >=200 ){
    double t=exp(-x);
    return t/(t+1);
  }else{
    return 1;
  }
}
inline double logfermi(double x){
  if (fabs(x)>33){
    if (x>33){
      double t = exp(-x);
      return -x-t*(1-0.5*t);
    }else{
      double t = exp(x);
      return -t*(1-0.5*t*(1-2./3.*t));
    }
  }else{
    return log(1/(exp(x)+1.0));
  }
}

inline long double nbose(double x)
{
  if (fabs(x)<300)
    return 1/(exp(x)-1.0);
  else if ( x >=200 )
    return 0.0;
  else
    return -1.0;
}

// Returns value of Binomial Coefficient C(n, k)
int binomialCoeff(int n, int k)
{
  int res = 1;
  // Since C(n, k) = C(n, n-k)
  if (k > n - k)
    k = n - k;
 
  // Calculate value of
  // [n * (n-1) *---* (n-k+1)] / [k * (k-1) *----* 1]
  for (int i = 0; i < k; ++i) {
    res *= (n - i);
    res /= (i + 1);
  }
  return res;
}
 
//********** Creates Python sign function **************/
template <typename T> int sign(T val) {
    return (T(0) < val) - (val < T(0));
};
template <typename T> double dsign(T val) {
  return (T(0) < val) - (val < T(0));
}
//********** Creates Python range function **************/
template <typename Scalar>
inline bl::Array<Scalar,1> arange(Scalar start, Scalar stop, Scalar step)
{
  if (step == Scalar(0))
  {
    throw std::invalid_argument("step for range must be non-zero");
  }
  int size = static_cast<int>((stop-start+step-0.5*sign(step))/step);
  if (size<=0) return bl::Array<Scalar,1>();
  bl::Array<Scalar,1> result( size );
  Scalar i = start;
  int ii=0;
  while ((step > 0) ? (i < stop) : (i > stop))
  {
    result(ii) = i;
    i += step;
    ++ii;
  }
  return result;
}
template <typename Scalar>
inline bl::Array<Scalar,1> arange(Scalar start, Scalar stop)
{
  return arange(start, stop, Scalar(1));
}
template <typename Scalar>
inline bl::Array<Scalar,1> arange(Scalar stop)
{
  return arange(Scalar(0), stop, Scalar(1));
}

inline int minus_1_2_n(int n)
{ return  n%2 ? -1 : 1;}

//********** TinyVector summation and subtraction **************/
bl::TinyVector<double,3> operator-(const bl::TinyVector<double,3>& a, const bl::TinyVector<double,3>& b)
{ return bl::TinyVector<double,3>(a[0]-b[0],a[1]-b[1],a[2]-b[2]);}
bl::TinyVector<double,3> operator+(const bl::TinyVector<double,3>& a, const bl::TinyVector<double,3>& b)
{ return bl::TinyVector<double,3>(a[0]+b[0],a[1]+b[1],a[2]+b[2]);}
bl::TinyVector<double,3> operator*(const bl::TinyVector<double,3>& a, double c)
{return bl::TinyVector<double,3>(a[0]*c,a[1]*c,a[2]*c);}

template<typename T>
inline T sqr(const T& x){ return x*x;}

//********** Trivial function that has functionality of python bisection, but does not use real bisection **************/
template <int N_rank>
inline int tBisect(double x, const bl::TinyVector<double,N_rank>& Prs){
  for (int i=0; i<N_rank; i++){
    if (x<Prs[i]) return i;
  }
  return N_rank-1;
}
/*
double ipower(double base, int exp){
  if( exp == 0)
    return 1.;
  double temp = ipower(base, exp/2);
  if (exp%2 == 0)
    return temp*temp;
  else {
    if(exp > 0)
      return base*temp*temp;
    else
      return (temp*temp)/base; //negative exponent computation 
  }
} 
*/
template<class T>
inline T ipower(const T& base, int exp)
{
  switch (exp){
  case 0 : return 1.; break;
  case 1 : return base; break;
  case 2 : return base*base; break;
  case 3 : return base*base*base; break;
  default :
    T _base_ = (exp>0) ? base : 1.0/base;
    if (exp<0) exp = -exp;
    T result = 1;
    while (exp){
      if (exp & 1) result *= _base_;
      exp >>= 1;
      _base_ *= _base_;
    }
    return result;
  }
}


double pi = M_PI;

inline bool isPowerOfTwo (int x)
{
  /* First x in the below expression is for the case when x is 0 */
  return x && (!(x&(x-1)));
}
inline double romberg(const bl::Array<double,1>& ff, double dh)
{
  int m = ff.extent(0);
  if ( !isPowerOfTwo(m-1) ) std::cout<<"ERROR : for Romberg method the size of the array should be 2^k+1" << std::endl;
  int n=m-1, Nr=0;
  while (n!=1){
    n = n/2;
    Nr++;
  }
  int N=Nr+1;

  bl::Array<double,1> h(N+1);
  bl::Array<double,2> r(N+1,N+1);
  
  for (int i = 1; i < N + 1; ++i) {
    h(i) = dh / ( 1<<(i-1) ) ;
  }
  r(1,1) = h(1) / 2 * (ff(0) + ff(m-1));
  for (int i = 2; i < N + 1; ++i) {
    double coeff = 0;
    int dr =  1 << (N-i);
    for (int k = 1; k <= (1<<(i-2)); ++k) coeff += ff((2*k-1)*dr);
    r(i,1) = 0.5 * (r(i-1,1) + h(i-1)*coeff);
  }
  for (int i = 2; i < N + 1; ++i) {
    for (int j = 2; j <= i; ++j) {
      r(i,j) = r(i,j - 1) + (r(i,j-1) - r(i-1,j-1))/( (1<<(2*(j-1))) -1 );
    }
  }
  return r(N,N);
}

template<int m, int Nr>
inline double romberg2(const bl::TinyVector<double,m>& ff, double dh){
  const int N=Nr+1;
  int m2 = (1<<Nr) + 1;
  if (m2 != m){
    std::cout<<"ERROR : romberg should have number of points 2**k+1"<<std::endl;
  }
  double h[N+1], r[N+1][N+1];
  for (int i = 1; i < N + 1; ++i) {
    h[i] = dh / ( 1<<(i-1) ) ;
  }
  r[1][1] = h[1] / 2 * (ff[0] + ff[m-1]);
  for (int i = 2; i < N + 1; ++i) {
    double coeff = 0;
    int dr = 1 << (N-i);
    for (int k = 1; k <= (1<<(i-2)); ++k) coeff += ff[(2*k-1)*dr];
    r[i][1] = 0.5 * (r[i-1][1] + h[i-1]*coeff);
  }
  for (int i = 2; i < N + 1; ++i) {
    for (int j = 2; j <= i; ++j) {
      r[i][j] = r[i][j - 1] + (r[i][j-1] - r[i-1][j-1])/( (1<<(2*(j-1))) -1);
    }
  }
  return r[N][N];
}
/*
template<int m, int Nr>
double romberg2(const bl::TinyVector<double,m>& ff, double dh){
  const int N=Nr+1;
  int m2 = (1<<Nr) + 1;
  if (m2 != m){
    std::cout<<"ERROR : romberg should have number of points 2**k+1"<<std::endl;
  }
  double h[N+1], r[N+1][N+1];
  for (int i = 1; i < N + 1; ++i) {
    h[i] = dh / static_cast<int>(pow(2.0,i-1));
  }
  r[1][1] = h[1] / 2 * (ff[0] + ff[m-1]);
  for (int i = 2; i < N + 1; ++i) {
    double coeff = 0;
    int dr = static_cast<int>(pow(2.0,N-i));
    for (int k = 1; k <= static_cast<int>(pow(2.0, i-2)); ++k) coeff += ff[(2*k-1)*dr];
    r[i][1] = 0.5 * (r[i-1][1] + h[i-1]*coeff);
  }
  for (int i = 2; i < N + 1; ++i) {
    for (int j = 2; j <= i; ++j) {
      r[i][j] = r[i][j - 1] + (r[i][j-1] - r[i-1][j-1])/(static_cast<int>(pow(4., j-1))-1);
    }
  }
  return r[N][N];
}
*/

inline bl::Array<int,2> BinomialCoeff(int n){
  // Returns value of Binomial Coefficient C(n, k)
  bl::Array<int,2> C(n+1,n+1);
  // Calculate value of Binomial Coefficient in bottom up manner
  for (int i = 0; i <= n; i++){
    for (int j = 0; j <= i; j++){
      // Base Cases
      if (j == 0 || j == i) C(i,j) = 1;
      // Calculate value using previosly stored values
      else
        C(i,j) = C(i-1,j-1) + C(i-1,j);
    }
  }
  return C;
}

inline static unsigned int log2 (unsigned int val) {
  if (val == 0) return UINT_MAX;
  if (val == 1) return 0;
  unsigned int ret = 0;
  while (val > 1) {
    val >>= 1;
    ret++;
  }
  return ret;
}

template<typename T, int N>
inline double norm2(const bl::TinyVector<T,N>& vector)
{
  double sum = 0.0;
  for (int i=0; i < N; ++i){
    sum += vector[i]*vector[i];
  }
  return sum;
}

#ifndef BZ_VECNORM_CC
#define BZ_VECNORM_CC
template<typename T, int N>
inline double norm(const bl::TinyVector<T,N>& vector)
{
  double sum = 0.0;
  for (int i=0; i < N; ++i){
    sum += vector[i]*vector[i];
  }
  return sqrt(sum);
}
#endif

inline double simpson_nonuniform(const bl::Array<double,1>& x, const bl::Array<double,1>& f)
{ //
  // Simpson rule for irregularly spaced data.
  //
  //    Parameters
  //    ----------
  //    x : array of floats
  //            Sampling points for the function values
  //    f : array of floats
  //            Function values at the sampling points
  //
  //    Returns
  //    -------
  //    float : approximation for the integral
  //
  int N = x.extent(0)-1;
  bl::Array<double,1> h(N+1);
  h=0;
  for (int i=0; i<x.extent(0)-1; i++) h(i) = x(i+1)-x(i);
  
  double result = 0.0;
  for (int i=1; i<N; i+=2){// in range(1, N, 2):
    double hph = h(i) + h(i-1);
    double hi_3  = h(i)*h(i)*h(i);
    double him_3 = h(i-1)*h(i-1)*h(i-1);
    double hip_3 = h(i+1)*h(i+1)*h(i+1);
    //
    result += f(i)   * ( hi_3 + him_3 + 3.* h(i)*h(i-1)*hph       ) / ( 6 * h(i)*h(i-1) );
    result += f(i-1) * ( 2.*him_3 - hi_3  + 3.*h(i)*h(i-1)*h(i-1) ) / ( 6 * h(i-1)*hph);
    result += f(i+1) * ( 2.*hi_3  - him_3 + 3.*h(i-1)*h(i)*h(i)   ) / ( 6 * h(i) * hph );
  }
  if ( (N+1) % 2 == 0){ // N is odd
    result += f(N)   * ( 2*h(N-1)*h(N-1) + 3.*h(N-2)*h(N-1)) / ( 6*( h(N-2) + h(N-1) ) );
    result += f(N-1) * ( h(N-1)*h(N-1) + 3*h(N-1)*h(N-2) )   / ( 6 * h(N-2) );
    result -= f(N-2) * h(N-1)*h(N-1)*h(N-1) / ( 6 * h(N-2) * ( h(N-2) + h(N-1) ) );
  }
  return result;
}

template<class T>
inline T simpson_nonuniform2(const bl::Array<double,1>& x, const bl::Array<T,1>& f)
{ //
  // Simpson rule for irregularly spaced data.
  //
  //    Parameters
  //    ----------
  //    x : array of floats
  //            Sampling points for the function values
  //    f : array of floats
  //            Function values at the sampling points
  //
  //    Returns
  //    -------
  //    float : approximation for the integral
  //
  const double dh_min = 1e6;
  int N = x.extent(0)-1;
  if (N<=1){
    if (N<1) return 0;
    return 0.5*(f(0)+f(1))*(x(1)-x(0));
  }
  bl::Array<double,1> h(N);
  for (int i=0; i<N; i++) h(i) = x(i+1)-x(i);
  double avg_h = sum(h)/N;
  double dh_max = avg_h*dh_min;
  T result = 0.0;
  for (int i=1; i<N; i+=2){
    double h0 = h(i-1);
    double h1 = h(i);
    double hph = h1 + h0;
    double hdh = h1 / h0;
    double hmh = h1 * h0;
    double hpg6 = hph/6;
    double dh0 = hpg6 * (2-hdh);
    double dh1 = hpg6 * (hph*hph/hmh);
    double dh2 = hpg6 * (2-1.0/hdh);
    if (fabs(dh0) < dh_max) result += f(i-1)*dh0;
    if (fabs(dh1) < dh_max) result += f(i)*dh1;
    if (fabs(dh2) < dh_max) result += f(i+1)*dh2;
  }
  if ( N % 2 == 1){ // N is odd
    double h0 = h(N-2);
    double h1 = h(N-1);
    double h12 = h1*h1;
    double h10 = h1*h0;
    double dh0 = (h12*2+3*h10)/(6*(h0+h1));
    double dh1 = (h12 + 3*h10)/(6*h0);
    double dh2 = h1*h1*h1/(6*h0*(h0 + h1));
    if (fabs(dh0) < dh_max) result += f(N)   * dh0;
    if (fabs(dh1) < dh_max) result += f(N-1) * dh1;
    if (fabs(dh2) < dh_max) result -= f(N-2) * dh2;
  }
  return result;
}

template<class functor>
double simpson_recursive(functor F, double a, double b, int& m0/*=2*/, double rprecision=1e-6, double aprecision=1e-10)
{
  // Simpsons integrating routine for given function F (can be lambda function in C++11) and finite interval [a,b]
  // We start with m0 regularly spaced points, and we double the points until desired accuracy is reached.
  // It uses recursive Simpson method, which gives the error estimate.
  // Routine will require at least 128 points in the integral, and one of four conditions met:
  //   a) relative precision is below rprecision
  //   b) absolute precision is below aprecision
  //   c) recursion level is 17, which means that the number of points used is m0*2^17 .
  //   d) number of points exceeds 5*10^5.
  // 
  const int max_recursion=17;
  if (m0<=1){ m0=0;  return 0;}
  if (b<a){
    std::cout << "m0 was " << m0 << " while a=" << a  << " and b=" << b << std::endl;
    m0=-1; return 0;
  }
  double dh = (b-a)/(m0-1.); // starting distance between points
  double s = 0;
  {// trapezoid rule for S_i, from which Simpson will be constructed below
    double f0 = F(a), fN = F(b);
    s = 0.5*(f0+fN);
    for (int i=1; i < m0-1; ++i) s += F(a+i*dh);
  }
  s *= dh;
  int m=m0-1;          // number of subdivisions in completed step
  double oldsm = 0, sm = 0; // simpson approximation is not known yet.
  double olds = s;
  //std::cout << "recursive starting ni="  << m+1 << " and s=" << s << std::endl;
  for (int i=0; i<max_recursion; i++, m*=2) { // m==2^(i-2)
    double x = a + 0.5*dh;       // first point added at
    double sum = 0;
    for (int j=0; j<m; j++, x+=dh) sum += F(x);
    sum *= dh;
    s = 0.5*olds + 0.5*sum;  // The new sum using trapezoid rule
    sm = (4.*s-olds)/3.;     // This is the first fifference with trapezoid
    //std::cout << "recursive ni=" << 2*m+1 << " sm=" << sm << " news=" << s << " olds=" << olds << " sum=" << sum << std::endl;
    double rprec, aprec;
    if (i>=1){
      rprec = std::abs(sm-oldsm)/(std::abs(oldsm));
      aprec = std::abs(sm-oldsm);
    }else{
      rprec = std::abs(sm-olds)/(std::abs(sm));
      aprec = std::abs(sm-olds);
    }
    if ((2*m+1>128) and  (rprec < rprecision or aprec < aprecision )){ //Avoid spurious early convergence.
      //std::clog<<"m=" << m << " prec=" << std::abs(sm-oldsm)/(std::abs(oldsm)) << std::endl;
      m0 = 2*m+1;
      return sm;
    }
    dh /= 2;
    olds=s;
    oldsm=sm;
    if (2*m+1 > 5e5) break;
  }
  //std::cerr<<"Too many steps in routine integrate_simps at a="<<a<<" b="  << b << " m0=" << m0 << " m=" << m << std::endl;
  m0 = 2*m+1;
  return sm; // The best approximation
}

double trapez(const bl::Array<double,1>& ff, const bl::Array<double,1>& x)
{
  int N = ff.extent(0);
  if (N != x.extent(0)) {std::cerr<<"Error in trapez ff.size != x.size "<<std::endl; exit(1);}
  double dsum=0;
  for (int i=0; i<N-1; i++){
    dsum += 0.5*(ff(i+1)+ff(i))*(x(i+1)-x(i));
  }
  return dsum;
}

template <class T>
double power(T p, int k)
{
  T s = 1;
  for (int i=0; i<k; i++) s*=p;
  return s;
}


#ifdef _DEBUG
#define Assert(condition, message)\
{\
  if(!(condition)) std::cerr << (message) << std::endl;\
}
#else /* NO_ARG_CHECK */
#define Assert(condition, message)
#endif /* NO_ARG_CHECK */


#endif //MY_UTIL
