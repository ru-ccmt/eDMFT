// @Copyright 2007 Kristjan Haule
// 
#ifndef MY_UTILS
#define MY_UTILS
#include <cmath>
#include <string>
#include <sstream>
#include "assert.h"

#define HPoffset 8

// template<class T>
// inline T sqr(T x)
// {
//   return x*x;
// }

using namespace std;

int sign(double x)
{
  if (x>0) return 1;
  else return -1;
  //  return 0;
}
template <class T>
T power(T p, int k)
{
  T s = 1;
  if (k>0)
    for (int i=0; i<k; i++) s*=p;
  if (k<0){
    T pm=1/p;
    for (int i=0; i>k; i--) s*=pm;
  }
  return s;
}

inline double _power(double p, double k)
{
  return exp(k*log(p));
}

inline double Power(double p, int k){
  if (abs(k)<50) return power(p,k);
  else return _power(p,k);
}
std::string NameOfFile(const std::string& name, int n, int m=3){
  std::stringstream str;
  str<<name<<".";
  for (int i=1; i<m; i++){
    if (n<power(10,i)) str<<"0";
  }
  str<<n<<std::ends;
  return std::string(str.str());
}
std::string NameOfFile_(const std::string& name, int n1, int n2, int m1=2, int m2=3){
  std::stringstream str;
  str<<name<<".";
  for (int i=1; i<m1; i++){
    if (n1<power(10,i)) str<<"0";
  }
  str<<n1;
  str<<".";
  for (int i=1; i<m2; i++){
    if (n2<power(10,i)) str<<"0";
  }
  str<<n2<<std::ends;
  return std::string(str.str());
}

//************************************************//
// Simple class used to exchange the data between //
// functions and meshes when doing linear	  //
// interpolation  				  //
//************************************************//
class intpar{
public:
  int i;
  double p;
  intpar(int i_, double p_) : i(i_), p(p_){}
  intpar(){};
  intpar(double d) : i(0), p(d) {};
};

enum TypeOfMatrix {_Normal, _Transpose};

template <int Om_gt_zero>
inline double fpferm_f(double f_eps_Om, double f_eps, double f_Om)
{
  return 0;
}
template <>
inline double fpferm_f<1>(double f_eps_Om, double f_eps, double f_Om)
{
  return f_eps_Om * (1-f_eps)/(1-f_Om);
}
template <>
inline double fpferm_f<0>(double f_eps_Om, double f_eps, double f_Om)
{
  return (1-f_eps_Om) * f_eps/f_Om;
}

template <int om_gt_zero>
inline double fpceta_f(double n_ceta_om, double f_ceta, double f_om)
{
  return 0;
}
template <>
inline double fpceta_f<1>(double n_ceta_om, double f_ceta, double f_om)
{
  return n_ceta_om * (1-f_ceta)/(1-f_om);
}
template <>
inline double fpceta_f<0>(double n_ceta_om, double f_ceta, double f_om)
{
  return (1+n_ceta_om) * f_ceta/f_om;
}

inline double ferm_f(double x)
{
  return 1.0/(1+exp(x));
}
inline double nbose_f(double x)
{
  return 1.0/(exp(x)-1);
}

template <class T>
inline T Dilog(const T& x)
{
  double nx = norm(x);
  T z = x;
  int cs = 1;
  if (nx>1.){
    if (x.real()>0.5) {
      z = 1-1/x;
      cs = 2;
    }
    else {
      z = 1-1/(1-x);
      cs = 3;
    }
  }
  
  int i=1;
  T sum = z, px = z, psum=1e10;
  while (norm(sum-psum)>1e-15){
    psum = sum;
    px *= z;
    sum += px/sqr(++i);
  }
  
  if (cs==1) return sum;
  if (cs==2) {
    T ln = log(x);
    return sum + 0.5*ln*(ln-2.*log(1-x)) + sqr(M_PI)/6;
  }
  if (cs==3) return -sum-0.5*sqr(log(1-x));
  else return sum;
}

///////// Simple Newton Method /////////////
template <class functor, class T>
inline T findRoot(const T& x0, functor& f, double tol, int MaxSteps)
{
  Assert( tol > 0, "Tolerance must be positive");
  T a = x0, fa, dfa;

  for (int i=0; i<MaxSteps; i++){
    f(a, fa, dfa);
    a -= fa/dfa;
    if (abs(fa)<tol) return a;
  }
  cerr << "Can't find root after " << MaxSteps;
  cerr << " iterations. The value of the function is " << fa << endl;
  
  if (isnan(fa)) return 0.0;
  else return a;
}

class bounds{
public:
  double lr, ur, li, ui;
};

template <class functor, class T>
inline T findRoot(const T& x0, functor& f, double tol, int MaxSteps, const bounds& m, bool& success)
{
  Assert( tol > 0, "Tolerance must be positive");
  success = true;
  T a = x0, an, da, fa, dfa;
  double alpha;
  for (int i=0; i<MaxSteps; i++){
    f(a, fa, dfa);
    da = -fa/dfa;
    an = a + da;
    alpha = 1.0;
    if (an.real()<m.lr) alpha = min(alpha, 0.75*(m.lr-a.real())/da.real());
    if (an.real()>m.ur) alpha = min(alpha, 0.75*(m.ur-a.real())/da.real());
    if (an.imag()<m.li) alpha = min(alpha, 0.75*(m.li-a.imag())/da.imag());
    if (an.imag()>m.ui) alpha = min(alpha, 0.75*(m.ui-a.imag())/da.imag());
    a += alpha*da;
    if (abs(fa)<tol) return a;
  }
  cerr << "Can't find root after " << MaxSteps;
  cerr << " iterations. The value of the function is " << fa << endl;
  if (abs(fa)>1e-8) success = false;
  if (isnan(fa)) return 0.0;
  else return a;
}

#define RED           "\033[31;1m"
#define GREEN         "\033[32;1m"
#define YELLOW        "\033[33;1m"
#define BLUE          "\033[34;1m"
#define PURPLE        "\033[35;1m"
#define CYAN          "\033[36;1m"
#define WHITE         "\033[37;1m"
#define URED          "\033[31;1;4m"
#define UGREEN        "\033[32;1;4m"
#define UYELLOW       "\033[33;1;4m"
#define UBLUE         "\033[34;1;4m"
#define UPURPLE       "\033[35;1;4m"
#define UCYAN         "\033[36;1;4m"
#define UWHITE        "\033[37;1;4m"
#define DEFAULTCOLOR "\033[0m"
 
#ifndef _PARALEL_
#define COLOR(color, arg) color << arg << DEFAULTCOLOR
#else
#define COLOR(color, arg) arg
#endif

#ifdef _FULL_LOG_
#define LOGOUT(stream, arg) stream<<arg
#else
#define LOGOUT(stream, arg)
#endif

#endif // MY_UTILS
