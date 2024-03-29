// @Copyright 2007 Kristjan Haule
// General purpose interpolating routines
// including 1D and 2D spline interpolation.
// The 2D-spline has a specific implementation with cubic splines in second dimension, and Hermite cubic spline in first dimension.
// It can handle non-uniform meshes, however, there must be  only two meshes, one for each dimension.
//
#include <iostream>
#include <fstream>
#include <assert.h>
#include <fstream>
#include <sstream>
#include <blitz/array.h>
#ifndef _INTERP
#define _INTERP

#include "util.h"
#include "tanmesh.h"
#include "bblas.h"

#ifndef MY_UTILS
class intpar{
public:
  int i;
  double p;
  intpar(int i_, double p_) : i(i_), p(p_){}
  intpar(){};
  bool empty() const{return (i==0 and p==0);}
  intpar(int c) : i(0), p(0) {}
};

class Linintpar{
public:
  int i;
  double p;
  Linintpar(int i_, double p_) : i(i_), p(p_){}
  Linintpar(){};
  bool empty() const{return (i==0 and p==0);}
  Linintpar(int c) : i(0), p(0) {}
};

#endif
/*
#ifndef _SBLAS

#ifdef NO_APPEND_FORTRAN
# define FNAME(x) x
#else
# define FNAME(x) x##_
#endif


extern "C" {
  void FNAME(dptsv)(const int* N, const int* NRHS, double* D, double* E, double* B, const int* LDB, int* INFO);
  void FNAME(zptsv)(const int* N, const int* NRHS, double* D, std::complex<double>* E, std::complex<double>* B, const int* LDB, int* INFO);
}
inline void xptsv_(const int* N, const int* NRHS, double* D, double* E, double* B, const int* LDB, int* INFO)
{  FNAME(dptsv)(N, NRHS, D, E, B, LDB, INFO);}

inline void xptsv_(const int* N, const int* NRHS, double* D, std::complex<double>* E, std::complex<double>* B, const int* LDB, int* INFO)
{  FNAME(zptsv)(N, NRHS, D, E, B, LDB, INFO);}

#endif
*/

//using namespace std;
namespace bl = blitz;

// Class for linear interpolation
class LinInterp{
  int Ntm1;
  double dt, beta;
  intpar p;
public:
  LinInterp(int Nt_, double dt_, double beta_):Ntm1(Nt_-1), dt(dt_), beta(beta_){}
  intpar& operator()(double t){
    p.i = (t/beta)*Ntm1;
    if (p.i>=Ntm1) p.i=Ntm1-1;
    p.p = (t/dt - p.i);
    //    cout<<"Inside LinInterp: "<<t<<" p.i="<<p.i<<" p.p="<<p.p<<" dt="<<dt<<" t/dt="<<t/dt<<" t/beta="<<t/beta<<endl;
    return p;
  }
};

inline int bisection(double x, int& klo, int& khi, int bi, const bl::Array<double,1>& om)
{// Basic routine for searching ordered table is bisection
  int k;
  khi=bi-1;
  while (khi-klo > 1) {
    k=(khi+klo) >> 1;
    if (om(k) > x) khi=k;
    else klo=k;
  }
  return klo;
}
inline int find(double x, const bl::Array<double,1>& om)
{ // If nothing is known about the point x
  if (x<=om(0)) return 0;
  int N = om.extent(0);
  if (x>=om(N-2)) return N-2;
  int ai0=0, ai1=N-1;
  return bisection(x,ai0,ai1,N,om);
}
inline intpar Interp(double x, const bl::Array<double,1>& om)
{
  int i = find(x,om);
  double delta = (i < om.extent(0)-1) ? 1.0/(om(i+1)-om(i)) : 0.0;
  return intpar(i,(x-om(i))*delta);
}
inline Linintpar LInterp(double x, const bl::Array<double,1>& om)
{
  int i = find(x,om);
  double delta = (i < om.extent(0)-1) ? 1.0/(om(i+1)-om(i)) : 0.0;
  return Linintpar(i,(x-om(i))*delta);
}

inline int find_(double x, int& ai0, const bl::Array<double,1>& om)
{ // This is most often called searching routine
  // It is used for searching table in increasing order
  static const int dN = 10; // when searching ordered table, searching is done first between a0 and a0+dN...
  if (x<om(ai0+1)) return ai0; // Checks weather x is stil in [ai0:ai0+1]
  int N = om.extent(0);
  int ai1 = ai0+1;             // Makes a step
  if (ai0>=N-2) return ai0;    // Needs to check for the end of the table
  if (x<om(ai1+1)){            // Checks weather x is in [ai0+1:ai0+2]
    ai0 = ai1;
    ai1 = ai1+1;
    return ai0;
  }
  if (ai1>=N-2) return ai1; // Again checks for the end of the table
  ai0 = ai1+1;              // makes another step
  if (ai1+dN<N){            // First uses bisection is small interval between [ai1:ai1+dN]
    if (x<om(ai1+dN)) return bisection (x, ai0, ai1, ai1+dN+1, om);
  } else return bisection (x, ai0, ai1, N, om);
  if (ai1+dN<N-1) ai0 = ai1+dN;
  else ai0 = ai1+dN-1;
  return bisection (x, ai0, ai1, N, om); // If still not found, use bisection on the rest of the grid
}
inline int _find(double x, int& bi, const bl::Array<double,1>& om)
{ // This routine is used to search ordered table in decreasing order
  static const int dN = 10; // when searching ordered table, searching is done first between a0 and a0+dN...
  int ai0=0;
  if (x>om(bi)) return bi;     // Checks weather x is still in [bi:bi+1]
  if (bi<=0) return bi;        // Checks start of the table
  if (x>om(bi-1)) return bi-1; // Checks weather x is in [bi-1:bi]
  if (bi-2<=ai0) return ai0;   // If [bi-2:bi-1] is first interval (equal to [0:1]) we are done
  if (x>om(bi-2)) return bi-2; // Checks interbal [bi-2:bi-1]
  if (bi-dN>ai0){              // Bisection only between [bi-dN:ai0]
    if (x>om(bi-dN)){
      int ai1, ai00;
      ai00 = bi-dN;
      return bisection(x, ai00, ai1, bi-1, om); 
    }
    bi=bi-dN+1;
  } else bi-=1; 
  {                          // If everything else failed, search everywhere between [ai0:bi]
    int ai1;
    return bisection (x, ai0, ai1, bi, om);
  }
}
inline intpar InterpLeft(const double x, int& a, const bl::Array<double,1>& om)
{
  int i = find_(x,a,om);
  double delta = (i < om.extent(0)-1) ? 1.0/(om(i+1)-om(i)) : 0.0;
  return intpar(i,(x-om(i))*delta); // return value optimization
}
inline intpar InterpRight(const double x, int& a, const bl::Array<double,1>& om)
{
  int i = _find(x,a,om);
  double delta = (i < om.extent(0)-1) ? 1.0/(om(i+1)-om(i)) : 0.0;
  return intpar(i,(x-om(i))*delta); // return value optimization
}
inline intpar InterpBoth(const double x, int& a, const bl::Array<double,1>& om)
{
  int i;
  if (x==om(a)) i = a;
  if (x>om(a)) i = find_(x,a,om);
  else i = _find(x,a,om);
  double delta = (i < om.extent(0)-1) ? 1.0/(om(i+1)-om(i)) : 0.0;
  return intpar(i,(x-om(i))*delta); // return value optimization
}

template <class T>
inline T linInterp(const bl::Array<T,1>& f, const intpar& ip)
{
  int i= ip.i; double p = ip.p, q=1-ip.p;
  return q*f(i) + p*f(i+1);
}
template <class T>
inline T linInterp(const bl::Array<T,2>& f, int l, const intpar& ip)
{
  int i= ip.i; double p = ip.p, q=1-ip.p;
  return q*f(l,i) + p*f(l,i+1);
}

template<class T,int N>
bl::TinyVector<T,N> linInterp(const bl::Array<bl::TinyVector<T,N>,1>& f, const intpar& ip)
{
  int i= ip.i; double p = ip.p, q=1-ip.p;
  return q*f(i) + p*f(i+1);
}

template <class T>
class Spline1D{
public:
  bl::Array<T,1> f;        // function values
  bl::Array<T,1> f2;       // second derivatives
  bl::Array<double,1> dxi; // x_{j+1}-x_j
public:
  Spline1D(){};
  explicit Spline1D(int N) : f(N), f2(N), dxi(N) {};
  void splineIt(const bl::Array<double,1>& om, 
                T df0=std::numeric_limits<double>::max(), // the derivatives at both ends
                T dfN=std::numeric_limits<double>::max()); 
  Spline1D(const Spline1D& m); //copy constructor
  void resize(int N);
  int size() const { return f.extent(0);}
  T& operator[](int i) {if (i<0) i+=f.extent(0); return f(i);}
  const T& operator[](int i) const { if(i<0) i+=f.extent(0); return f(i);}
  T operator()(const intpar& ip) const; // spline interpolation
  T InterpWithLinearExtrpolation(const intpar& ip) const;  
  T operator()(const Linintpar& ip) const; // spline interpolation
  T df(const intpar& ip) const;
  T df(int i) const;
  T d2f(const intpar& ip) const;
  T d2f(int i) const;
  Spline1D& operator=(const Spline1D& m); // copy operator
  T integrate();
  std::complex<double> Fourier(double om, const bl::Array<double,1>& tau) const;
  T Sum_iom(const bl::Array<double,1>& iom) const;
};


template <class T>
inline void Spline1D<T>::resize(int N)
{
  if (N!=f.extent(0)){
    f.resize(N);
    f2.resize(N);
    dxi.resize(N);
  }
}

template <class T>
inline Spline1D<T>::Spline1D(const Spline1D& m)
{
  resize(m.size());
  f = m.f;
  f2 = m.f2;
  dxi = m.dxi;
}

template <class T>
inline Spline1D<T>& Spline1D<T>::operator=(const Spline1D<T>& m)
{
  resize(m.size());
  f = m.f;
  f2 = m.f2;
  dxi = f.dxi;
  return *this;
}

template <class T>
inline T Spline1D<T>::operator()(const intpar& ip) const
{
  int i= ip.i; double p = ip.p, q=1-ip.p;
  return q*f(i) + p*f(i+1) + dxi(i)*dxi(i)*(q*(q*q-1)*f2(i) + p*(p*p-1)*f2(i+1))/6.;
}
template <class T>
inline T Spline1D<T>::InterpWithLinearExtrpolation(const intpar& ip) const
{
  int i= ip.i; double p = ip.p, q=1-ip.p;
  if (p>=0 and p<=1){
    return q*f(i) + p*f(i+1) + dxi(i)*dxi(i)*(q*(q*q-1)*f2(i) + p*(p*p-1)*f2(i+1))/6.;
  }else{
    return q*f(i) + p*f(i+1);
  }
}
template <class T>
inline T Spline1D<T>::operator()(const Linintpar& ip) const
{
  int i= ip.i; double p = ip.p, q=1-ip.p;
  return q*f(i) + p*f(i+1);// + dxi(i)*dxi(i)*(q*(q*q-1)*f2(i) + p*(p*p-1)*f2(i+1))/6.;
}

template <class T>
inline T Spline1D<T>::df(const intpar& ip) const
{
  int i= ip.i; double p = ip.p, q=1-ip.p;
  return (f(i+1)-f(i))/dxi(i) + (1-3*q*q)*dxi(i)*f2(i)/6. - (1-3*p*p)*dxi(i)*f2(i+1)/6.;
}

template <class T>
inline T Spline1D<T>::df(int i) const
{
  if (i==size()-1){
    return (f(i)-f(i-1))/dxi(i-1) + dxi(i-1)*(f2(i-1)+2*f2(i))/6.;
  }else  return (f(i+1)-f(i))/dxi(i) - dxi(i)*(2*f2(i)+f2(i+1))/6.;
}
template <class T>
inline T Spline1D<T>::d2f(const intpar& ip) const
{
  int i= ip.i; double p = ip.p, q=1-ip.p;
  return p*f2(i+1)+q*f2(i);
}

template <class T>
inline T Spline1D<T>::d2f(int i) const
{
  return f2(i);
}



template<class T>
inline void Spline1D<T>::splineIt(const bl::Array<double,1>& om, T df0, T dfN)
{
  int N = f.extent(0);
  if (om.extent(0)!=N) std::cerr<<"Sizes of om and f are different in spline setup"<<std::endl;
  if (N<2){ for (int i=0; i<f.extent(0); i++) f2(i)=0; return;}
  
  if (df0==std::numeric_limits<double>::max() || dfN==std::numeric_limits<double>::max()){
    // tring to get the derivative at the edges from extrapolation.
    double x1 = 0.5*(om(1)+om(0));
    T df1 = (f(1)-f(0))/(om(1)-om(0));
    double x2 = 0.5*(om(2)+om(1));
    T df2 = (f(2)-f(1))/(om(2)-om(1));
    df0 = df1 + (df2-df1)*(om(0)-x1)/(x2-x1);
    x1 = 0.5*(om(N-1)+om(N-2));
    df1 = (f(N-1)-f(N-2))/(om(N-1)-om(N-2));
    x2 = 0.5*(om(N-2)+om(N-3));
    df2 = (f(N-2)-f(N-3))/(om(N-2)-om(N-3));
    dfN = df1 + (df2-df1)*(om(N-1)-x1)/(x2-x1);
  }
  bl::Array<double,1> diag(N);
  bl::Array<T,1> offdiag(N-1); // matrix is stored as diagonal values + offdiagonal values
  // Below, matrix and rhs is setup
  diag(0) = (om(1)-om(0))/3.;
  T dfu0 = (f(1)-f(0))/(om(1)-om(0));
  f2(0) = dfu0-df0;
  for (int i=1; i<N-1; i++){
    diag(i) = (om(i+1)-om(i-1))/3.;
    T dfu1 = (f(i+1)-f(i))/(om(i+1)-om(i));
    f2(i) = dfu1-dfu0;
    dfu0 = dfu1;
  }
  diag(N-1) = (om(N-1)-om(N-2))/3.;
  f2(N-1) = dfN - (f(N-1)-f(N-2))/(om(N-1)-om(N-2));
  for (int i=0; i<N-1; i++) offdiag(i) = (om(i+1)-om(i))/6.;

  // The system of symmetric tridiagonal equations is solved by lapack
  int one=1, info=0;
  xptsv_(&N, &one, diag.data(), offdiag.data(), f2.data(), &N, &info);
  
  if (info!=0) std::cerr<<"dptsv return an error "<<info<<std::endl;
  // Setup of other necessary information for doing splines.
  for (int i=0; i<N-1; i++) dxi(i) = (om(i+1)-om(i));
  dxi(N-1)=0;
}

template <class T>
inline T Spline1D<T>::integrate()
{
  T sum=0;
  for (int i=0; i<size()-1; i++) sum += 0.5*dxi(i)*(f(i+1)+f(i)-(f2(i+1)+f2(i))*dxi(i)*dxi(i)/12.);
  return sum;
}

template <class T>
inline std::complex<double> Spline1D<T>::Fourier(double om, const bl::Array<double,1>& tau) const
{
  std::complex<double> ii(0,1);
  std::complex<double> sum=0;
  std::complex<double> val;
  for (int i=0; i<f.extent(0)-1; i++) {
    double u = om*dxi(i), u2=u*u, u4=u2*u2;
    if (fabs(u)<1e-4){// Taylor expansion for small u
      val = 0.5*(1.+ii*(u/3.))*f(i)+0.5*(1.+ii*(2*u/3.))*f(i+1);
      val -= dxi(i)*dxi(i)/24.*(f2(i)*(1.+ii*7.*u/15.)+f2(i+1)*(1.+ii*8.*u/15.));
    }else{
      std::complex<double> exp(cos(u),sin(u));
      val  = (f(i)*(1.+ii*u-exp) + f(i+1)*((1.-ii*u)*exp-1.))/u2;
      val += dxi(i)*dxi(i)/(6.*u4)*(f2(i)*(exp*(6+u2)+2.*(u2-3.*ii*u-3.))+f2(i+1)*(6+u2+2.*exp*(u2+3.*ii*u-3.)));
    }
    sum += dxi(i)*std::complex<double>(cos(om*tau(i)),sin(om*tau(i)))*val;
  }
  return sum;
}

template <class T>
inline T Spline1D<T>::Sum_iom(const bl::Array<double,1>& om) const
{
  if (om.extent(0)!=f.extent(0)) {std::cerr<<"Sizes of iom and myself are different in spline1D::Sum_iom! Boiling out"<<std::endl; exit(1);}
  if (f.extent(0)<1) return 0;
  if (f.extent(0)<2) return f(0);
  // trapezoid part
  T sum = f(0)*0.5*(om(1)-om(0)+1);
  for (int i=1; i<f.extent(0)-1; i++){
    sum += f(i)*0.5*(om(i+1)-om(i-1));
  }
  sum += f(f.extent(0)-1)*0.5*(om(f.extent(0)-1)-om(f.extent(0)-2)+1);
  // correction due to splines
  T corr=0;
  for (int i=0; i<f.extent(0)-1; i++){
    corr -= (om(i+1)-om(i)+1)*(om(i+1)-om(i))*(om(i+1)-om(i)-1)*(f2(i)+f2(i+1))/24.;
  }
  return sum + corr;
}

template <class T>
std::ostream& operator<<(std::ostream& stream, const Spline1D<T>& f)
{
  for (int i=0; i<f.size(); i++) stream<<i<<" "<<f.f(i)<<std::endl;
  return stream;
}

void inverse_fourier_boson(bl::Array<double,1>& W_tau, double beta, const bl::Array<double,1>& tau, int nom, const bl::Array<int,1>& iom, const bl::Array<double,1>& W_iom, double b=0.0, int order=1, int iq_debug=-1)
{
  // order==1 : at high frequency we expect  : W = c/( om^2 + b^2)
  // order==2 : at high frequency we expect  : W = c/( om^2 + b^2)^2
  
  bool debug=false;
  if (iq_debug>=0) debug = true;
  std::ofstream* log;
  if (debug){
    log = new std::ofstream( std::string("debug.")+std::to_string(iq_debug) );
    (*log) << "# b = "<< b << "  ";
    (*log).precision(15);
  }
  
  bl::Array<double,1> omtail(iom.extent(0)-nom);
  for (int i=0; i<omtail.extent(0); i++) omtail(i) = iom(i+nom);
  
  int n_iom = iom.extent(0);    // integer mesh of imaginary frequency points
  int iw_last = iom(n_iom-1);   // the last integer on matsubara axis, it is usually a very large integer
  bl::Array<double,1> Wq0(iw_last); // this will contain Wq on all Matsubara points up to iw_last
  
  Spline1D<double> Wqtail(n_iom-nom); // first we need to spline the tail. Preparing the container
  double c=0.0;
  if (b>0 && (order==1 || order==2)){
    // subtract something analytic, namely
    //   order==1  : c/(2*b) * ( 1/(iw+b)-1/(iw-b) ) = c/( om^2 + b^2)
    //   order==2  : c/(om^2 + b^2)^2
    
    double om_last = 2*iw_last*M_PI/beta;
    double W_last = W_iom(n_iom-1);
    double om_plast = 2*iom(n_iom-2)*M_PI/beta;
    double W_plast = W_iom(n_iom-2);

    if (order==1){
      c = (om_last*om_last - om_plast*om_plast)/(1/W_last - 1/W_plast);
      b = sqrt(fabs((W_last*om_last*om_last - W_plast*om_plast*om_plast)/(W_plast-W_last)));

      for (int i=0; i<nom; i++){
	double omw = 2*iom(i)*M_PI/beta;
	Wq0(i) = W_iom(i) - c/(omw*omw+b*b);
	if (debug) (*log) << omw << "  " << W_iom(i) << "  " << Wq0(i) << "  " << c/(omw*omw+b*b) << std::endl;
      }
      for (int i=nom; i<n_iom; i++){
	double omw = 2*iom(i)*M_PI/beta;
	Wqtail[i-nom] = W_iom(i) - c/(omw*omw+b*b);
	if (debug) (*log) << omw << "  " << W_iom(i) << "  " << Wqtail[i-nom] << "  " << c/(omw*omw+b*b) << std::endl;
      }
    }else{// order==2
      double signc = sign(W_last);
      W_last = sqrt(fabs(W_last));
      W_plast = sqrt(fabs(W_plast));
      c = (om_last*om_last - om_plast*om_plast)/(1/W_last - 1/W_plast);
      c = c*c*signc;
      b = sqrt(fabs((W_last*om_last*om_last - W_plast*om_plast*om_plast)/(W_plast-W_last)));

      for (int i=0; i<nom; i++){
	double omw = 2*iom(i)*M_PI/beta;
	Wq0(i) = W_iom(i) - c/ipower((omw*omw+b*b),2);
	if (debug) (*log) << omw << "  " << W_iom(i) << "  " << Wq0(i) << "  " << c/ipower((omw*omw+b*b),2) << std::endl;
      }
      for (int i=nom; i<n_iom; i++){
	double omw = 2*iom(i)*M_PI/beta;
	Wqtail[i-nom] = W_iom(i) - c/ipower((omw*omw+b*b),2);
	if (debug) (*log) << omw << "  " << W_iom(i) << "  " << Wqtail[i-nom] << "  " << c/ipower((omw*omw+b*b),2) << std::endl;
      }
    }
    if (debug) (*log) << " om_last="<< om_last << " W_last=" << W_last << " c="<< c << " om_plast=" << om_plast << " W_plast="  << W_plast << " b=" << b << std::endl;
  }else{ // no need to subtract anything
    if (debug) (*log) << std::endl;
    for (int i=0; i<nom; i++){
      Wq0(i) = W_iom(i);
      if (debug) (*log) << 2*iom(i)*M_PI/beta << "  " << W_iom(i) << "  " << Wq0(i) << std::endl;
    }
    for (int i=nom; i<n_iom; i++){
      Wqtail[i-nom] = W_iom(i);
      if (debug) (*log) << 2*iom(i)*M_PI/beta << "  " << W_iom(i) << "  " << Wqtail[i-nom] << std::endl;
    }
  }

  /*
  ofstream debug("debug1.dat");
  for (int i=0; i<omtail.extent(0); i++){
    debug<<omtail(i)<<" "<<Wqtail[i]<<endl;
  }
  debug.close();
  exit(0);
  */
  
  Wqtail.splineIt(omtail); // finally spline the tail

  int ia=0;
  for (int i=nom; i<iw_last; i++){  // Now determine the function at all points on Matsubara axis.
    Wq0(i) = Wqtail( InterpLeft(i, ia, omtail) );
  }
  for (int it=0; it<tau.extent(0); it++){
    double t = tau(it);
    double dsum=0;
    for (int i=1; i<iw_last; i++) dsum += cos( 2*i*M_PI/beta * t) * Wq0(i);
    W_tau(it) = ( 2*dsum + Wq0(0) )/beta;
  }
  if (b>0 && (order==1 || order==2)){
    double cpf = (order==1) ? c/(2*b) : c/(4*b*b*b);
    for (int it=0; it<tau.extent(0); it++){
      double t = tau(it);
      if (order==1)
	// c/(2*b)*(1/(iw+b)-1/(iw-b)) == c/(w**2+b**2)	
	W_tau(it) += cpf * (exp(-b*t)+exp(-b*(beta-t)))/(1.-exp(-beta*b));
      else{
	// c/(w**2+b**2)^2
	double den = 1./ipower( 1.-exp(-beta*b) , 2);
	W_tau(it) += cpf*den*( (1+b*t)*exp(-b*t) + (1+b*(beta-t))*exp(-b*(beta-t)));
	W_tau(it) += cpf*den*( (b*t-1)*exp(-b*(2*beta-t)) + (b*(beta-t)-1)*exp(-b*(beta+t)));
      }
    }
  }
}
void inverse_fourier_boson2(bl::Array<double,1>& W_tau, double beta, const bl::Array<double,1>& tau, int nom, const bl::Array<int,1>& iom, const bl::Array<double,1>& W_iom, bl::Array<double,1> b, int iq_debug=-1)
{
  // order==1 : at high frequency we expect  : W = b[0]/( om^2 + b[1]**2)
  bool debug=false;
  if (iq_debug>=0) debug = true;
  std::ofstream* log;
  if (debug){
    log = new std::ofstream( std::string("debug.")+std::to_string(iq_debug) );
    (*log).precision(15);
    (*log) << "# b = "<< b(0) << "  " << b(1) << " ";
  }
  
  bl::Array<double,1> omtail(iom.extent(0)-nom);
  for (int i=0; i<omtail.extent(0); i++) omtail(i) = iom(i+nom);
  
  int n_iom = iom.extent(0);    // integer mesh of imaginary frequency points
  int iw_last = iom(n_iom-1);   // the last integer on matsubara axis, it is usually a very large integer
  bl::Array<double,1> Wq0(iw_last); // this will contain Wq on all Matsubara points up to iw_last
  
  Spline1D<double> Wqtail(n_iom-nom); // first we need to spline the tail. Preparing the container
  double unit = 2*M_PI/beta;
  double c = b(0)*unit*unit;
  double d = b(1)*unit;
  if (b(0)!=0){
    // subtract something analytic, namely
    //   c/(2*d) * ( 1/(iw+d)-1/(iw-d) ) = c/( om^2 + d^2)
    for (int i=0; i<nom; i++){
      double omw = iom(i)*unit;
      Wq0(i) = W_iom(i) - c/(omw*omw+d*d);
      if (debug) (*log) << omw << "  " << W_iom(i) << "  " << Wq0(i) << "  " << c/(omw*omw+d*d) << std::endl;
    }
    for (int i=nom; i<n_iom; i++){
      double omw = iom(i)*unit;
      Wqtail[i-nom] = W_iom(i) - c/(omw*omw+d*d);
      if (debug) (*log) << omw << "  " << W_iom(i) << "  " << Wqtail[i-nom] << "  " << c/(omw*omw+d*d) << std::endl;
    }
    if (debug) (*log) << " c="<< c << " d=" << d << std::endl;
  }else{ // no need to subtract anything
    if (debug) (*log) << std::endl;
    for (int i=0; i<nom; i++){
      Wq0(i) = W_iom(i);
      if (debug) (*log) << 2*iom(i)*M_PI/beta << "  " << W_iom(i) << "  " << Wq0(i) << std::endl;
    }
    for (int i=nom; i<n_iom; i++){
      Wqtail[i-nom] = W_iom(i);
      if (debug) (*log) << 2*iom(i)*M_PI/beta << "  " << W_iom(i) << "  " << Wqtail[i-nom] << std::endl;
    }
  }

  /*
  ofstream debug("debug1.dat");
  for (int i=0; i<omtail.extent(0); i++){
    debug<<omtail(i)<<" "<<Wqtail[i]<<endl;
  }
  debug.close();
  exit(0);
  */
  
  Wqtail.splineIt(omtail); // finally spline the tail

  int ia=0;
  for (int i=nom; i<iw_last; i++){  // Now determine the function at all points on Matsubara axis.
    Wq0(i) = Wqtail( InterpLeft(i, ia, omtail) );
  }
  for (int it=0; it<tau.extent(0); it++){
    double t = tau(it);
    double dsum=0;
    for (int i=1; i<iw_last; i++) dsum += cos( 2*i*M_PI/beta * t) * Wq0(i);
    W_tau(it) = ( 2*dsum + Wq0(0) )/beta;
  }
  if (b(0)!=0){
    double cpf = c/(2*d);
    for (int it=0; it<tau.extent(0); it++){
      double t = tau(it);
      // c/(2*b)*(1/(iw+b)-1/(iw-b)) == c/(w**2+b**2)	
      W_tau(it) += cpf * (exp(-d*t)+exp(-d*(beta-t)))/(1.-exp(-beta*d));
    }
  }
}

void inverse_fourier_boson(bl::Array<double,1>& W_tau, double beta, const bl::Array<double,1>& tau, int nom, const bl::Array<int,1>& iom, const bl::Array<std::complex<double>,1>& W_iom, double b=0.0)
{
  bl::Array<double,1> omtail(iom.extent(0)-nom);
  for (int i=0; i<omtail.extent(0); i++) omtail(i) = iom(i+nom);
  
  int n_iom = iom.extent(0);    // integer mesh of imaginary frequency points
  int iw_last = iom(n_iom-1);   // the last integer on matsubara axis, it is usually a very large integer
  bl::Array<std::complex<double>,1> Wq0(iw_last); // this will contain Wq on all Matsubara points up to iw_last
  
  Spline1D<double> rWqtail(n_iom-nom); // first we need to spline the tail. Preparing the container
  Spline1D<double> iWqtail(n_iom-nom); // first we need to spline the tail. Preparing the container
  double a=0, eps=1;
  double c=0.0;
  const double very_small=1e-10;
  if (fabs(W_iom(n_iom-1).imag())>very_small){ // We subtract f = a/(iom+eps) = a * (eps-iom)/(eps^2+om^2)
    std::complex<double> W_last = W_iom(n_iom-1);
    double om_last = 2*iw_last*M_PI/beta;
    a  = -om_last*W_last.imag();
    eps = om_last*om_last*W_last.real()/a;
    for (int i=0; i<nom; i++){
      double omw = 2*iom(i)*M_PI/beta;
      Wq0(i) = W_iom(i) - a/std::complex<double>(eps,omw);
    }
    for (int i=nom; i<n_iom; i++){
      double omw = 2*iom(i)*M_PI/beta;
      rWqtail[i-nom] = W_iom(i).real() - a*eps/(omw*omw+eps*eps);
      iWqtail[i-nom] = W_iom(i).imag() + a*omw/(omw*omw+eps*eps);
    }
    b=0;
  }else if (b>0){ // subtract something analytic, namely c/(2*b) * ( 1/(iw+b)-1/(iw-b) )
    double om_last = 2*iw_last*M_PI/beta;
    double W_last = W_iom(n_iom-1).real();
    c = W_last*om_last*om_last;
    double om_plast = 2*iom(n_iom-2)*M_PI/beta;
    double W_plast = W_iom(n_iom-2).real();
    b = sqrt(fabs(-(W_plast*om_plast*om_plast-c)*om_plast*om_plast/c));
    //cout<<iq<<"  c="<<c<<"  b="<<b<<endl;
    for (int i=0; i<nom; i++){
      double omw = 2*iom(i)*M_PI/beta;
      Wq0(i) = W_iom(i) - c/(omw*omw+b*b);
    }
    for (int i=nom; i<n_iom; i++){
      double omw = 2*iom(i)*M_PI/beta;
      rWqtail[i-nom] = W_iom(i).real() - c/(omw*omw+b*b);
      iWqtail[i-nom] = W_iom(i).imag();
    }
  }else{ // no need to subtract anything
    for (int i=0; i<nom; i++) Wq0(i) = W_iom(i);
    for (int i=nom; i<n_iom; i++){
      rWqtail[i-nom] = W_iom(i).real();
      iWqtail[i-nom] = W_iom(i).imag();
    }
  }

  /*
  ofstream debug("debug1.dat");
  for (int i=0; i<omtail.extent(0); i++){
    debug<<omtail(i)<<" "<<Wqtail[i]<<endl;
  }
  debug.close();
  exit(0);
  */
  
  rWqtail.splineIt(omtail); // finally spline the tail
  iWqtail.splineIt(omtail); // finally spline the tail

  int ia=0;
  for (int i=nom; i<iw_last; i++){  // Now determine the function at all points on Matsubara axis.
    intpar p = InterpLeft(i, ia, omtail);
    Wq0(i) = std::complex<double>(rWqtail(p),iWqtail(p));
  }
  for (int it=0; it<tau.extent(0); it++){
    double t = tau(it);
    double dsum=0;
    for (int i=1; i<iw_last; i++) dsum += (cos( 2*i*M_PI/beta * t) * Wq0(i).real() + sin( 2*i*M_PI/beta * t) * Wq0(i).imag());
    W_tau(it) = ( 2*dsum + Wq0(0).real() )/beta;
    // c/(2*b)*(1/(iw+b)-1/(iw-b)) == c/(w**2+b**2)
    if (a!=0){
      if (fabs(beta*eps) > 50){
	W_tau(it) += (eps>0) ? a*exp(-(beta-t)*eps)  : -a*exp(t*eps);
      }else{
	W_tau(it) += a * exp( (t-0.5*beta)*eps )/(2*sinh(0.5*beta*eps));
      }
    }
    if (b>0) W_tau(it) += c/(2*b)*(exp(-b*t)+exp(-b*(beta-t)))/(1.-exp(-beta*b));
  }
}

void inverse_fourier_fermion(bl::Array<double,1>& Gtau, double beta, const bl::Array<double,1>& tau, int nom, const bl::Array<int,1>& iom, const bl::Array<std::complex<double>,1>& Giom, const bl::TinyVector<double,2>& ah)
{
  //std::ofstream debug0("debug0.dat");
  // The difference between G and its analytic approximation  ah[0]/(iom - ah[1])
  bl::Array<std::complex<double>,1> dG(iom.size());
  for (int i=0; i<iom.size(); i++){
    std::complex<double> omw( 0, (2*iom(i)+1)*M_PI/beta );
    dG(i) = Giom(i) - ah[0]/(omw - ah[1]);
    //debug0 << iom(i) << " " << Giom(i).real() << " " << Giom(i).imag() << " " << dG(i).real() << " " << dG(i).imag() << std::endl;
  }
  
  int n_iom = iom.extent(0);    // integer mesh of imaginary frequency points
  int iw_last = iom(n_iom-1);   // the last integer on matsubara axis, it is usually a very large integer

  //std::ofstream debug1("debug1.dat");
  //debug1 << "# " << ah[0] << " " << ah[1] << std::endl;
  //for (int i=0; i<iom.size(); i++){
  //  debug1 << iom(i) << " " << dG(i).real() << " " << dG(i).imag() << std::endl;
  //}
  
  // only the tail stored here
  bl::Array<double,1> omtail(n_iom-nom);
  for (int i=0; i < n_iom-nom; i++) omtail(i) = iom(i+nom);
  Spline1D<double> r_dG(n_iom-nom); // first we need to spline the tail. Preparing the container
  Spline1D<double> i_dG(n_iom-nom); // first we need to spline the tail. Preparing the container
  for (int i=nom; i<n_iom; i++){
    r_dG[i-nom] = dG(i).real();
    i_dG[i-nom] = dG(i).imag();
  }
  r_dG.splineIt(omtail); // finally spline the tail
  i_dG.splineIt(omtail); // finally spline the tail
  bl::Array<std::complex<double>,1> dGa(iw_last); // this will contain Wq on all Matsubara points up to iw_last
  for (int i=0; i<nom; i++) dGa(i) = dG(i);
  int ia=0;
  for (int i=nom; i<iw_last; ++i){  // Now determine the function at all points on Matsubara axis.
    intpar p = InterpLeft(i, ia, omtail);
    dGa(i) = std::complex<double>(r_dG(p),i_dG(p));
  }

  //std::ofstream debug2("debug2.dat");
  //for (int i=0; i<iw_last; i++){
  //  debug2 << i << " " << dGa(i).real() << " " << dGa(i).imag() << std::endl;
  //}
  
  for (int it=0; it<tau.size(); it++){
    double t = tau(it);
    double dsum=0;
    for (int i=0; i<iw_last; i++){
      double omw = (2*i+1)*M_PI/beta;
      dsum += cos(omw*t)*dGa(i).real() + sin(omw*t)*dGa(i).imag();
    }
    Gtau(it) = 2*dsum/beta;
    if (ah[1]>0)
      Gtau(it) -= ah[0]*exp(-ah[1]*tau(it))/(1.+exp(-ah[1]*beta));
    else
      Gtau(it) -= ah[0]*exp( ah[1]*(beta-tau(it)))/(1.+exp(ah[1]*beta));
  }
  //std::ofstream debug3("debug3.dat");
  //for (int it=0; it<tau.size(); it++){
  //  debug3 << tau(it) << " " << Gtau(it) << std::endl;
  //}  
}



template <class T>
class Spline2D{
  // Note that the x-dimension should be equidistant. The y-dimension does not need to be equidistant.
public:
  bl::Array<T,2> f;        // function values
  bl::Array<T,2> f2;       // second derivatives
  bl::Array<double,1> dyi; // y_{j+1}-y_j
  bl::Array<double,1> dxi; // (x_{j+1}-x_j)
public:
  explicit Spline2D(){};
  explicit Spline2D(int N0, int N1) : f(N0,N1), f2(N0,N1), dyi(N1), dxi(N0) {};
  void resize(int N0, int N1){
    f.resize(N0,N1);
    f2.resize(N0,N1);
    dyi.resize(N1);
    dxi.resize(N0);
  }
  void splineIt(const bl::Array<double,1>& x0, const bl::Array<double,1>& om); 
  T& operator()(int i, int j) {if (i<0) i+=f.extent(0); if (j<0) j+=f.extent(1); return f(i,j);}
  const T& operator()(int i, int j) const { if(i<0) i+=f.extent(0); if (j<0) j+=f.extent(1); return f(i,j);}
  T linear(const intpar& ipx, const intpar& ipy) const;
  T operator()(const intpar& ip0, const intpar& ip1) const; // spline interpolation
};

template<class T>
inline void Spline2D<T>::splineIt(const bl::Array<double,1>& x0, const bl::Array<double,1>& om)
{
  if (f.extent(0) != x0.extent(0)) std::cerr<<"Sizes of x0 and f.extent(0) are different in Spline2D.splineIt"<<std::endl;
  int N = f.extent(1);
  if (N!=om.extent(0)) std::cerr<<"Sizes of om and f.extent(1) are different in Spline2D.splineIt"<<std::endl;
  
  double x1 = 0.5*(om(1)+om(0));
  double x2 = 0.5*(om(2)+om(1));
  double xn1 = 0.5*(om(N-1)+om(N-2));
  double xn2 = 0.5*(om(N-2)+om(N-3));

  bl::Array<double,1> diag(N);
  diag(0) = (om(1)-om(0))/3.;
  for (int i=1; i<N-1; i++){
    diag(i) = (om(i+1)-om(i-1))/3.;
  }
  diag(N-1) = (om(N-1)-om(N-2))/3.;
  
  bl::Array<T,1> offdiag(N-1); // matrix is stored as diagonal values + offdiagonal values
  for (int i=0; i<N-1; i++) offdiag(i) = (om(i+1)-om(i))/6.;

  for (int j=0; j<f.extent(0); j++){
    if (N<2){
      for (int i=0; i<N; i++) f2(j,i)=0;
      continue;
    }
    T df0, dfN;
    {
      // tring to get the derivative at the edges from extrapolation.
      double df1 = (f(j,1)-f(j,0))/(om(1)-om(0));
      double df2 = (f(j,2)-f(j,1))/(om(2)-om(1));
      df0 = df1 + (df2-df1)*(om(0)-x1)/(x2-x1);
      double dfn1 = (f(j,N-1)-f(j,N-2))/(om(N-1)-om(N-2));
      double dfn2 = (f(j,N-2)-f(j,N-3))/(om(N-2)-om(N-3));
      dfN = dfn1 + (dfn2-dfn1)*(om(N-1)-xn1)/(xn2-xn1);
    }
    
    // Below, matrix and rhs is setup
    T dfu0 = (f(j,1)-f(j,0))/(om(1)-om(0));
    f2(j,0) = dfu0-df0;
    for (int i=1; i<N-1; i++){
      T dfu1 = (f(j,i+1)-f(j,i))/(om(i+1)-om(i));
      f2(j,i) = dfu1-dfu0;
      dfu0 = dfu1;
    }
    f2(j,N-1) = dfN - (f(j,N-1)-f(j,N-2))/(om(N-1)-om(N-2));
    
    bl::Array<double,1> _diag_(diag.shape());
    _diag_ = diag;
    bl::Array<double,1> _offdiag_(offdiag.shape());
    _offdiag_ = offdiag;
    
    // The system of symmetric tridiagonal equations is solved by lapack
    int one=1, info=0;
    xptsv_(&N, &one, _diag_.data(), _offdiag_.data(), f2.data()+j*N, &N, &info);
    
    if (info!=0){
      std::cerr<<"dptsv return an error "<<info<<std::endl;
      std::ofstream debug("debug.dat");
      debug<<" # j="<<j<<std::endl;
      for (int i=0; i<N; i++){
	debug<<om(i)<<" "<<f(i)<<" "<<f2(i)<<" "<<diag(i)<<" "<<offdiag(i)<<std::endl;
      }
      debug.close();
      exit(0);
    }
  }
  // Setup of other necessary information for doing splines.
  dyi(N-1)=0;
  for (int i=0; i<N-1; i++) dyi(i) = (om(i+1)-om(i));
  dxi(f.extent(0)-1)=0;
  for (int i=0; i<x0.extent(0)-1; i++) dxi(i) = (x0(i+1)-x0(i));
}


template <class T>
inline T Spline2D<T>::linear(const intpar& ipx, const intpar& ipy) const
{
  double py[2];
  int i = ipy.i; double p = ipy.p, q=1-ipy.p;
  double qq = q*(q*q-1)*dyi(i)*dyi(i)/6.;
  double pp = p*(p*p-1)*dyi(i)*dyi(i)/6.;
  for (int k=0; k<2; k++){ // spline interpolation at the ipy point, but for many x-values.
    int j = ipx.i + k;
    py[k] = q*f(j,i) + p*f(j,i+1) + qq*f2(j,i) + pp*f2(j,i+1); // First we interpolate four-splines in ipy point.
  }
  return py[0] + (py[1]-py[0])*ipx.p;
}

bl::TinyVector<double,4> hermiteP(double t)
{
  //double t = ipx.p;
  double t2 = t*t;
  double tm12 = sqr(t-1);
  //double h00 = (1+2*t)*tm12;
  //double h10 = t*tm12;
  //double h01 = t2*(3-2*t);
  //double h11 = t2*(t-1);
  return bl::TinyVector<double,4>( (1+2*t)*tm12, t*tm12, t2*(3-2*t), t2*(t-1) ); // (h00,h10,h01,h11)
}

template <class T>
inline T Spline2D<T>::operator()(const intpar& ipx, const intpar& ipy) const
{
  int N0 = f.extent(0);
  if (ipx.i >= N0-2 && ipx.p > 1.1) return 0.0;
  int i = ipy.i; double p = ipy.p, q=1-ipy.p;
  double qq = q*(q*q-1)*dyi(i)*dyi(i)/6.;
  double pp = p*(p*p-1)*dyi(i)*dyi(i)/6.;
  double py[4];
  int k_start = ipx.i>0    ? 0 : 1;
  int k_end   = ipx.i<N0-2 ? 4 : 3;
  for (int k=k_start; k<k_end; k++){ // spline interpolation at the ipy point, but for many x-values.
    int j = ipx.i + k - 1;
    py[k] = q*f(j,i) + p*f(j,i+1) + qq*f2(j,i) + pp*f2(j,i+1); // First we interpolate four-splines in ipy point.
  }
  bool MONOTONE = false;
  // Now that we have four values at ipy, we can interpolate through these four points a cubic interpolation
  // We will use  Cubic Hermite spline.
  int j = ipx.i;
  if (j>0 && j<N0-2){
    /*
    double dfa = (py[1]-py[0]) * dxi(j)/ dxi(j-1);
    double dfb = (py[2]-py[1]);
    double dfc = (py[3]-py[2]) * dxi(j)/ dxi(j+1);
    double df0 = 0.5*(dfa+dfb);
    double a = 0.5*(dfa + dfc) - dfb;
    double b = 0.5*(dfb - dfa) - a;
    return py[1] + ipx.p*( df0 + ipx.p * ( b + ipx.p * a ));
    */
    // Cubic Hermite spline
    double dfa = (py[1]-py[0])*dxi(j)/dxi(j-1); // df/dx(j-1) * dx(j)
    double dfb = (py[2]-py[1]);                 // df/dx(j)   * dx(j) 
    double dfc = (py[3]-py[2])*dxi(j)/dxi(j+1); // df/dx(j+1) * dx(j)
    double ml_dx = 0.5*(dfa+dfb); // derivative at the beginning * dx
    double mr_dx = 0.5*(dfb+dfc); // derivative at the end point * dx
    bl::TinyVector<double,4> H = hermiteP(ipx.p);
    if (MONOTONE){
      double beta1  = ml_dx/dfa;
      double alpha2 = ml_dx/dfb;
      double beta2  = mr_dx/dfb;
      double alpha3 = mr_dx/dfc;
      if (beta1  > 3) ml_dx = dfa*3;
      if (alpha2 > 3) ml_dx = dfb*3;
      if (beta2  > 3) mr_dx = dfb*3;
      if (alpha3 > 3) mr_dx = dfc*3;
    }
    return py[1]*H[0] + py[2]*H[2] + H[1]*ml_dx + H[3]*mr_dx;
  }
  if (j==0){
    double dfb = (py[2]-py[1]);
    double dfc = (py[3]-py[2]) * dxi(j) / dxi(j+1);
    /*
    double df1_over_df0 = fabs(0.5*dfc/dfb+1.5);
    if (df1_over_df0 > 3. || df1_over_df0 < 1/3.){ // large difference in derivative at both ends
      return py[1] + ipx.p*( dfb + 0.5*(dfb-dfc) * ipx.p * (1-ipx.p));
    }else{
      return py[1] + ipx.p*( dfb + 0.5*(dfb-dfc) * (1-ipx.p));
    }
    */
    double ml_dx = dfb;
    double mr_dx = 0.5*(dfb+dfc);
    bl::TinyVector<double,4> H = hermiteP(ipx.p);
    if (MONOTONE){
      double beta2  = mr_dx/dfb;
      double alpha3 = mr_dx/dfc;
      if (beta2  > 3) mr_dx = dfb*3;
      if (alpha3 > 3) mr_dx = dfc*3;
    }
    return py[1]*H[0] + py[2]*H[2] + H[1]*ml_dx + H[3]*mr_dx;
  }
  if (j==N0-2){
    double dfa = (py[1]-py[0]) * dxi(j)/ dxi(j-1);
    double dfb = (py[2]-py[1]);
    /*
    double df0_over_df1 = 0.5*(1.+dfa/dfb);
    if (df0_over_df1 > 3. || df0_over_df1 < 1/3.){ // large difference in derivative at both ends
      return py[1] + ipx.p*( 0.5*(dfb+dfa) + ipx.p * (dfb-dfa) * (1-0.5*ipx.p));
    }else{
      return py[1] + ipx.p*( 0.5*(dfb+dfa) + ipx.p * (dfb-dfa) * 0.5);
    }
    */
    double ml_dx = 0.5*(dfa+dfb); // derivative at the beginning * dx
    double mr_dx = dfb;           // derivative at the end point * dx
    bl::TinyVector<double,4> H = hermiteP(ipx.p);
    if (MONOTONE){
      double beta1  = ml_dx/dfa;
      double alpha2 = ml_dx/dfb;
      if (beta1  > 3) ml_dx = dfa*3;
      if (alpha2 > 3) ml_dx = dfb*3;
    }
    return py[1]*H[0] + py[2]*H[2] + H[1]*ml_dx + H[3]*mr_dx;
  }
  return 0.0; // should not happen.
}

double my_cinterpolate(const intpar& ipx, const bl::Array<double,1> fx, const bl::Array<double,1> x)
{
  int j = ipx.i;
  int N0 = x.extent(0);
  //dxi(i) = (x0(i+1)-x0(i));
  bool MONOTONE = false;
  if (j>0 && j<N0-2){
    double py[4] = {fx(j-1), fx(j), fx(j+1), fx(j+2)};
    double dfa = (py[1]-py[0]) * (x(j+1)-x(j))/(x(j)-x(j-1));
    double dfb = (py[2]-py[1]);
    double dfc = (py[3]-py[2]) * (x(j+1)-x(j))/(x(j+2)-x(j+1));
    /*
    double df0 = 0.5*(dfa+dfb);
    double a = 0.5*(dfa + dfc) - dfb;
    double b = 0.5*(dfb - dfa) - a;
    return py[1] + ipx.p*( df0 + ipx.p * ( b + ipx.p * a ));
    */
    double ml_dx = 0.5*(dfa+dfb); // derivative at the beginning * dx
    double mr_dx = 0.5*(dfb+dfc); // derivative at the end point * dx
    bl::TinyVector<double,4> H = hermiteP(ipx.p);
    if (MONOTONE){
      double beta1  = ml_dx/dfa;
      double alpha2 = ml_dx/dfb;
      double beta2  = mr_dx/dfb;
      double alpha3 = mr_dx/dfc;
      if (beta1  > 3) ml_dx = dfa*3;
      if (alpha2 > 3) ml_dx = dfb*3;
      if (beta2  > 3) mr_dx = dfb*3;
      if (alpha3 > 3) mr_dx = dfc*3;
    }
    return py[1]*H[0] + py[2]*H[2] + H[1]*ml_dx + H[3]*mr_dx;
  }
  if (j==0){
    double py[4] = {0.0, fx(j), fx(j+1), fx(j+2)};
    double dfb = (py[2]-py[1]);
    double dfc = (py[3]-py[2]) * (x(j+1)-x(j))/(x(j+2)-x(j+1));
    /*
    double df1_over_df0 = fabs(0.5*dfc/dfb+1.5);
    if (df1_over_df0 > 3. || df1_over_df0 < 1/3.){ // large difference in derivative at both ends
      return py[1] + ipx.p*( dfb + 0.5*(dfb-dfc) * ipx.p * (1-ipx.p));
    }else{
      return py[1] + ipx.p*( dfb + 0.5*(dfb-dfc) * (1-ipx.p));
    }
    */
    double ml_dx = dfb;
    double mr_dx = 0.5*(dfb+dfc);
    bl::TinyVector<double,4> H = hermiteP(ipx.p);
    if (MONOTONE){
      double beta2  = mr_dx/dfb;
      double alpha3 = mr_dx/dfc;
      if (beta2  > 3) mr_dx = dfb*3;
      if (alpha3 > 3) mr_dx = dfc*3;
    }
    return py[1]*H[0] + py[2]*H[2] + H[1]*ml_dx + H[3]*mr_dx;
  }
  if (j==N0-2){
    double py[4] = {fx(j-1), fx(j), fx(j+1), 0.0};
    double dfa = (py[1]-py[0]) * (x(j+1)-x(j))/(x(j)-x(j-1));
    double dfb = (py[2]-py[1]);
    /*
    double df0_over_df1 = 0.5*(1.+dfa/dfb);
    if (df0_over_df1 > 3. || df0_over_df1 < 1/3.){ // large difference in derivative at both ends
      return py[1] + ipx.p*( 0.5*(dfb+dfa) + ipx.p * (dfb-dfa) * (1-0.5*ipx.p));
    }else{
      return py[1] + ipx.p*( 0.5*(dfb+dfa) + ipx.p * (dfb-dfa) * 0.5);
    }
    */
    double ml_dx = 0.5*(dfa+dfb); // derivative at the beginning * dx
    double mr_dx = dfb;           // derivative at the end point * dx
    bl::TinyVector<double,4> H = hermiteP(ipx.p);
    if (MONOTONE){
      double beta1  = ml_dx/dfa;
      double alpha2 = ml_dx/dfb;
      if (beta1  > 3) ml_dx = dfa*3;
      if (alpha2 > 3) ml_dx = dfb*3;
    }
    return py[1]*H[0] + py[2]*H[2] + H[1]*ml_dx + H[3]*mr_dx;
  }
  return 0;// should not happen anyway
}

inline void MakeEquidistantMesh(bl::Array<double,1>& tau, int N_, double start, double end)
{ // For building equidistant mesh
  tau.resize(N_);
  double x=start, dh=(end-start)/(N_-1.);
  for (int i=0; i<N_; i++,x+=dh) tau(i) = x;
  //SetUp(N/2);
}

bl::Array<double,1> CalcLogoOnMesh(const bl::Array<double,1>& om)
{
  bl::Array<double,1> f(om.size());
  int N=om.size();
  for (int i=0; i<om.size(); ++i)
    f(i) = (i!=0 && i!=N-1) ? log((om(N-1)-om(i))/(om(i)-om(0))) : 0.0;
  return f;
}

bl::Array<double,1> KramarsKronig(bl::Array<double,1>& fi, const bl::Array<double,1>& om, const bl::Array<double,1>& dh, const bl::Array<double,1>& logo)
{
  if (logo.size()!=om.size() or om.size()!=fi.size()){
    std::cerr<<"In KramarsKronig functions not of equal size"<<std::endl;
    exit(1);
  }
  int N=om.size();
  bl::Array<double,1> fr(N);
  for (int i=0; i<om.size(); i++){
    double fderiv;
    if (i>0 and i<N-1){
      fderiv = 0.5*( (fi(i+1)-fi(i))/(om(i+1)-om(i)) + (fi(i)-fi(i-1))/(om(i)-om(i-1)) );
    }else{
      if (i==0){
	fderiv = (fi(1)-fi(0))/(om(1)-om(0));
      }else{
	fderiv = (fi(N-1)-fi(N-2))/(om(N-1)-om(N-2));
      }
    }
    double fii = fi(i);
    double sum = 0;
    for (int j=0; j<om.size(); j++){
      if (i!=j)
	sum += (fi(j)-fii)*dh(j)/(om(j)-om(i));
      else
	sum += fderiv*dh(j);
    }
    double r = (sum+fii*logo(i))/M_PI;
    //std::cout<<om(i)<<" "<<r<<" "<<fii<<" "<<logo(i)<<std::endl;
    fr(i) = r;
  }
  return fr;
}



template <class T>
class Spline1DLin : public Spline1D<T>{
public:
  int i_start, i_end; // at which i should we start to use Linear interpolation instead of cubic interpolation.
public:
  explicit Spline1DLin(int N) : Spline1D<T>(N) {};
  T operator()(const intpar& ip) const; // spline interpolation
  T operator()(const Linintpar& ip) const{ return (1-ip.p)*this->f(ip.i) + ip.p*this->f(ip.i+1);}; 
  void splineIt(const bl::Array<double,1>& om,
		const bl::Array<T,1>& f_mid,
                T df0=std::numeric_limits<double>::max(), // the derivatives at both ends
		T dfN=std::numeric_limits<double>::max()); 
};

template <class T>
inline T Spline1DLin<T>::operator()(const intpar& ip) const
{
  int i= ip.i; double p = ip.p, q=1-ip.p;
  if (p<0 or p>1 or ip.i<i_start or ip.i>i_end){ // either extrapolation of ip.i is in the regime where linear interpolation should be used
    return q*this->f(i) + p*this->f(i+1);
  }else{
    return q*this->f(i) + p*this->f(i+1) + this->dxi(i)*this->dxi(i)*(q*(q*q-1)*this->f2(i) + p*(p*p-1)*this->f2(i+1))/6.;
  }
}
template <class T>
void Spline1DLin<T>::splineIt(const bl::Array<double,1>& om,
			      const bl::Array<T,1>& f_mid, // we require the exact function knowledge so that we compute it at midpoints.
			      T df0, T dfN)
{
  Spline1D<T>::splineIt(om,df0,dfN);
  // Now we go through all points in the mesh and check in the mid-point between two points 
  i_start=0; i_end=om.size()-2;
  int i_st=0; int i_en = om.size()-2;
  for (; i_st<om.size()-1; ++i_st){
    int i = i_st;
    double x_midpoint = 0.5*(om(i)+om(i+1));
    intpar px(i,0.5);
    Linintpar pxl(i,0.5);
    double f_exact = f_mid(i);
    double f_spline = operator()(px);
    double f_linear = operator()(pxl);
    if (fabs(f_exact-f_spline) <= fabs(f_exact-f_linear)) break;
  }
  for (; i_en>=0; --i_en){
    int i = i_en;
    double x_midpoint = 0.5*(om(i)+om(i+1));
    intpar px(i,0.5);
    Linintpar pxl(i,0.5);
    double f_exact = f_mid(i);
    double f_spline = operator()(px);
    double f_linear = operator()(pxl);
    if (fabs(f_exact-f_spline) <= fabs(f_exact-f_linear)) break;
  }
  i_start = i_st; i_end=i_en;
};



#endif //_INTERP
/*
int main () {
  using namespace std;
  {
    int Nw=40;
    double a = -10.0001*M_PI;
    double b = 10*M_PI;
    bl::Array<double,1> fx(Nw);
    bl::Array<double,1> x(Nw);
    ofstream debug("debug.dat");
    for (int i=0; i<Nw; i++){
      double t = a + (b-a)*i/(Nw-1.);
      x(i) = t;
      fx(i) = (sin(t)/t - cos(t))/(t*t) + 0.1*drand48();
      debug<< t <<" "<< fx(i) <<std::endl;
    }
    debug.close();
    for (int i=0; i<10*Nw; i++){
      double t = a + (b-a)*i/(10*Nw-1.);
      intpar p = Interp(t, x);
      double ft = my_cinterpolate(p, fx, x);
      cout<<t<<" "<<ft<<endl;
    }
    return 0;
  }
  {
    int Nw = 64;
    bl::Array<double,1> x(Nw), fx(Nw);
    double x0 = 0.2;
    double a = 1e-3;
    double b = 32*M_PI;
    GiveTanMesh0(x, x0, b, Nw);
    ofstream debug("debug.dat");
    for (int i=0; i<Nw; i++){
      double t = x(i);
      fx(i) = (sin(t)/t - cos(t))/(t*t);
      debug<< t <<" "<< fx(i) <<std::endl;
    }
    debug.close();
    for (int i=0; i<100*Nw; i++){
      double t = a + (b-a)*i/(100*Nw-1.);
      intpar p = Interp(t, x);
      double ft = my_cinterpolate(p, fx, x);
      cout<<t<<" "<<ft<<endl;
    }
    return 0;
  }
  bl::Array<double,1> tau, qx;
  bl::Array<double,2> Wq;
  ReadWq(tau, qx, Wq, "qx.dat", "Wqt.dat");

  int Nq = qx.extent(0);
  Spline2D<double> Wnq(qx.extent(0),tau.extent(0));
  for (int iq=0; iq<qx.extent(0); iq++)
    for (int it=0; it<tau.extent(0); it++)
      Wnq(iq,it) = Wq(iq,it);
  
  Wnq.splineIt(qx, tau);

  double dq = (qx(qx.extent(0)-1)-qx(0))/(40*Nq-1);
  for (int iq=0; iq<Nq*40; iq++){
    double q = 0.0 + iq*dq;
    intpar pq = Interp(q,qx);
    cout<<q<<" ";
    for (int it=0; it<tau.extent(0)-1; it++){
      double t = (tau(it)+tau(it+1))/2.;
      intpar pt = Interp(t,tau);
      double y = Wnq( pq, pt );
      cout<<y<<" ";
    }
    cout<<endl;
  }
  return 0;
}
*/
