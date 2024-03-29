// @Copyright 2007 Kristjan Haule
// 
#ifndef FUNCTION_
#define FUNCTION_
#include <iostream>
#include <algorithm>
#include <limits>
#include <complex>
#include <cstring>
#include "sutil.h"
#include "sblas.h"

//////////////////////////// Hierarchy of classes: ///////////////////////////////////////
//                             
//                                     base clas  = function<>
//                                               |                                   
//                                          function1D<>                             
//

//*****************************************//
// Classes also used in this header file   //
//*****************************************//
class intpar;

//typedef std::complex<double> dcomplex;

//*******************************************************//
// Classes and functions implemented in this header file //
//*******************************************************//
template<class T> class function;
template<class T> class function1D;
template<class T> class spline1D;
template<class T> class funProxy;
template<class T> class function2D;

//********************************************************************************//
// Base class for two derived classes: function1D<T> and function2D<T>.	  	  //
// It is also used as a proxy class for function2D. Function2D<T> consists of	  //
// arrays of function<T> rather than functions1D<T>.				  //
// Memory is allocated in a fortran-like fashion for better performance.	  //
// Linear interpolation is implemented with the operator() that takes one	  //
// argument (class intpar).							  //
//********************************************************************************//
template<class T>
class functionb{
protected:
  T *f;
  int N0, N;
  functionb() : f(NULL), N0(0), N(0) {};        // constructor is made protected such that mesh can not be instantiated
  explicit functionb(int N_) : N0(N_), N(N_) {};// This class is used only as base class
  ~functionb(){};
  functionb(const functionb&){};
public:
  // OPERATORS
  T& operator[](int i) {if (i<0) i+=N; Assert3(i<N0 && i>=0,"Out of range in functionb[] ", i, N0); return f[i];}
  const T& operator[](int i) const { if(i<0) i+=N; Assert3(i<N0 && i>=0, "Out of range in functionb[] ", i, N0); return f[i];}
  T operator()(const intpar& ip) const; // linear interpolation
  functionb& operator+=(const functionb& m);
  functionb& operator*=(const T& m);
  functionb& operator=(const T& c); // used for initialization
  // SHORT FUNCTIONS
  const T& last() const {return f[N-1];}
  int size() const { return N;}
  int fullsize() const { return N0;}
  T* MemPt() {return f;}
  const T* MemPt() const{return f;}
  // OTHER FUNCTIONS
  T accumulate();
  T* begin() const {return f;}
  T* end() const {return f+N;}
  void AddPart(const functionb& m, double x);

  void xgemv(const string& TN, const function2D<T>& A, const functionb<T>& x, double alpha=1.0, double beta=1.0);
protected:
  template <class U> friend class spline1D;
  template<class U> friend class function2D;
  template<class U> friend class function2Dm;
  template <class U> friend U scalar_product(const functionb<U>& f1, const functionb<U>& f2);
};

//******************************************************************//
// One dimensional functions derived from function<T>. It has it's  //
// own constructors and destructors.				    //
//******************************************************************//
template <class T>
class function1D : public functionb<T>{
public:
  // CONSTRUCTORS AND DESTRUCTORS
  function1D(){};                 // default constructor exists for making arrays
  explicit function1D(int N_);    // basic constructor
  ~function1D();                  // destructor
  void ReleaseMemory();
  function1D(const function1D& m);//copy constructor
  // INITIALIZATION ROUTINES
  void resize(int N_);
  void resize_preserve(int n);
  void resize_virtual(int N_);
  // OPERATORS
  function1D& operator=(const function1D& m); // copy constructor
  function1D& operator=(const T& c) {functionb<T>::operator=(c); return *this;}// used for initialization
  function1D<T>& copy_full(const function1D<T>& m);
};
//************************************************************************//
// One dimensional spline function derived from function<T>. It has it's  //
// own constructors and destructors.				          //
//************************************************************************//
template <class T>
class spline1D : public functionb<T>{
public:
  T* f2;       // second derivatives
  double *dxi; // x_{j+1}-x_j
public:
  // CONSTRUCTORS AND DESTRUCTORS
  spline1D() : f2(NULL), dxi(NULL) {};// default constructor exists for allocating arrays
  // constructor
  explicit spline1D(int N_);
  template <class mesh1D>
  void splineIt(const mesh1D& om, 
		const T& df0=std::numeric_limits<double>::max(), // the derivatives at both ends
		const T& dfN=std::numeric_limits<double>::max()); 
  ~spline1D();                  // destructor
  spline1D(const spline1D& m); //copy constructor
  // INITIALIZATION ROUTINES
  void resize(int N_);
  // OPERATORS
  T operator()(const intpar& ip) const; // spline interpolation
  T df(const intpar& ip) const;
  T df(int i) const;
  spline1D& operator=(const spline1D& m); // copy operator
  // ADVANCED FUNCTIONS
  T integrate();
  template <class mesh1D>
  complex<double> Fourier(double om, const mesh1D& xi) const;
  template <class mesh1D>
  dcomplex Fourier_(double om, const mesh1D& xi) const;
  template <class mesh1D>
  T Sum_iom(const mesh1D& iom) const;
};

template <class T>
class funProxy : public functionb<T>{
public:
  void Initialize(int N_, T* f_);
  void ReInitialize(int N0_, int N_, T* f_);
  void resize(int N_);
  funProxy& operator=(const functionb<T>& m);
  ~funProxy(){};
// private:
//   funProxy(const funProxy& ){};
};

//**********************************************************************//
// Two dimentional function<T> derived from function<T>. It consists	//
// of an array of function<T> rather tham function1D<T>.		//
// Constructor calls operator new and aferwords placement new operator	//
// to allocate the whole memory in one single large peace. 		//
//**********************************************************************//
template<class T>
class function2D{
protected:  
  void *memory;
  T* data;
  funProxy<T> *f;
  int N0, Nd0, N, Nd;
public:
  function2D() : memory(NULL), N0(0), Nd0(0), N(0), Nd(0) {};
  function2D(int N_, int Nd_);
  ~function2D();
  void ReleaseMemory();
  funProxy<T>& operator[](int i) {Assert3(i<N0 && i>=0,"Out of range in function2D[] ", i, N0); return f[i];}
  const funProxy<T>& operator[](int i) const {Assert3(i<N0 && i>=0,"Out of range in function2D[] ", i, N0); return f[i];}

  const T& operator()(int i, int j) const {Assert5(i<N0 && j<Nd0 && i>=0 && j>=0,"Out of range in function2D(i,j) ", i, j, N0, Nd0); return data[i*Nd0+j];}
  T& operator()(int i, int j) {Assert5(i<N0 && j<Nd0 && i>=0 && j>=0,"Out of range in function2D(i,j) ", i, j, N0, Nd0); return data[i*Nd0+j];}
  
  T* MemPt() { return data;}
  const T* MemPt() const { return data;}
  int size_N() const {return N;}
  int size_Nd() const {return Nd;}
  int fullsize_N() const {return N0;}
  int fullsize_Nd() const {return Nd0;}
  int fullsize2() const {return N0*Nd0;}
  int lda() const {return Nd0;}
  void resize(int N_, int Nd_);
  void resize_clear(int N_, int Nd_);
  void resize_virtual(int N_, int Nd_);
  void resize_preserve(int N_, int Nd_);
  
  function2D(const function2D<T>& fun);
  function2D& operator=(const function2D& m);
  function2D& operator+=(double x);
  function2D& operator+=(const function2D& m);
  function2D& operator-=(double x);
  function2D& operator-=(const function2D& m);
  function2D& operator=(const T& u);
  function2D& operator*=(const T& x);
  void SetEqual(const function2D& m);
  
  template <class functor>
  void Set(functor& functn);
  template <class functor, class W>
  void Set(functor& functn, const function2D<W>& Um);

  function2D& minus(const function2D& B, const function2D& C);
  function2D& plus(const function2D& B, const function2D& C);
  // Be careful: Product performs matrix multiplication of  C = A * B^T
  void Product(const function2D& A, const function2D& B, int start=0, int end=-1, const T& alpha=1, const T& beta=0);
  // Be carefu: TProduct performs matrix multiplication of C = B * A^T
  void TProduct(const function2D& A, const function2D& B, const T& alpha=1, const T& beta=0);
  // Be careful: MProduct performs matrix multiplication of C = A * B
  void MProduct(const function2D& A, const function2D& B, const T& alpha=1, const T& beta=0);
  void DotProduct(const function2D& A, const function2D& B, const T& alpha=1, const T& beta=0);
  void SMProduct(const function2D& A, const function2D& B, const T& alpha=1, const T& beta=0);
  void TMProduct(const function2D& A, const function2D& B, const T& alpha, const T& beta);
  void SymmProduct(const function2D& A, const function2D& B, const T& alpha=1, const T& beta=0);
  void TSymmProduct(const function2D& A, const function2D& B, const T& alpha=1, const T& beta=0);
  template <class meshx, class functiond>
  void KramarsKronig(const meshx& om, const functiond& logo, functiond& sum, functiond& odvSig);
  void Product(const std::string& transa, const std::string& transb, const function2D& A, const function2D& B, const T& alpha=1.0, const T& beta=0.0);
  int Inverse(const function2D<T>& A);
  void AddPart(const function2D& m, double x);
  friend void Multiply(function2D<double>& C, const function2D<double>& A, const function2D<double>& B);
  friend void Multiply_N_T(function2D<double>& C, const function2D<double>& A, const function2D<double>& B);
};

template<class T>
class function2DProxy : public function2D<T>{
// A class to take slices of other higher dimensional arrays, and making them appear 2D arrays
  void* _memory_;
public:
  function2DProxy(T* _data_, int _N_, int _Nd_, int _Nd0_)
  {
    this->data = _data_;
    this->N0 = _N_;
    this->N = _N_;
    this->Nd = _Nd_;
    this->Nd0 = _Nd0_;
    _memory_ = operator new (sizeof(funProxy<T>)*this->N0);
    Assert(_memory_!=NULL,"Out of memory");
    this->f = new (_memory_) funProxy<T>[this->N0];
    for (int i=0; i<this->N0; i++) this->f[i].Initialize(this->Nd0,this->data+i*this->Nd0);

    cout<<"N="<<this->N<<" Nd="<<this->Nd<<" N0="<<this->N0<<" Nd0="<<this->Nd0<<endl;
  }
  ~function2DProxy()
  {
    if (_memory_!=NULL){
      for (int i=0; i<this->N0; i++) this->f[i].~funProxy<T>();
      operator delete(_memory_);
      _memory_ = NULL;
    }
  }
  void resize(int N_, int Nd_) {cerr<<"ERROR : Can not resize Proxy class !"<<endl; exit(1);};
  void resize_clear(int N_, int Nd_) {cerr<<"ERROR : Can not resize Proxy class !"<<endl; exit(1);};
  void SetEqual(const function2D<T>& m) {cerr<<"ERROR : Can not resize Proxy class !"<<endl; exit(1);};
};

// functionb ////////////////////////////////////////////////////////////////
// Routine for linear interpolation
template<class T>
inline T functionb<T>::operator()(const intpar& ip) const
{ return  f[ip.i]+ip.p*(f[ip.i+1]-f[ip.i]);}

template<class T> 
inline functionb<T>& functionb<T>::operator+=(const functionb& m)
{
  Assert(N==m.size(),"Functions not of equal length! Can't sum!");
  for (int i=0; i<N; i++) f[i] += m[i];
  return *this;
}

template <class T>
inline void functionb<T>::AddPart(const functionb& m, double x)
{
  for (int i=0; i<N; i++) f[i] += m[i]*x;
}

template<class T>
inline functionb<T>& functionb<T>::operator*=(const T& m)
{
  for (int i=0; i<N; i++) f[i] *= m;
  return *this;
}

template <class T>
inline T functionb<T>::accumulate()
{
  T sum(0);
  for (int i=0; i<N; ++i) sum += f[i];
  return sum;
}

template <class T>
inline functionb<T>& functionb<T>::operator=(const T& c)
{
  _LOG(if (N<=0) std::cerr << "Size of functionb is non positive! "<<N<<std::endl;)
  for (int i=0; i<N; i++) f[i] = c;
  return *this;
}

// function1D ////////////////////////////////////////////////////////////
template<class T>
inline function1D<T>::function1D(int N_) : functionb<T>(N_)
{ this->f = new T[N_];}

template<class T>
inline function1D<T>::~function1D()
{
  if (this->f!=NULL) {
    delete[] this->f;
    this->f = NULL;
  }
  this->N = this->N0 = 0;
}

template<class T>
void inline function1D<T>::ReleaseMemory()
{
  //cout << "Release memory of function1D on pointer "<< this->f << endl;
  if (this->f!=NULL) {
    delete[] this->f;
    this->f = NULL;
  }
  this->N = this->N0 = 0;
}

template<class T>
inline void function1D<T>::resize(int n)
{
  if (n>this->N0){
    if (this->f) delete[] this->f;
    this->f = new T[n];
    this->N0=n;
  }
  this->N = n;
}

template<class T>
inline void function1D<T>::resize_preserve(int n)
{
  if (n>this->N0){
    T* new_f = new T[n]; // The default constructor is called on all the data
    for (int i=0; i<this->N0; i++) new_f[i] = this->f[i];
    if (this->f) delete[] this->f;
    this->f = new_f;
    this->N0=n;
  }
  //this->N = n;
}

template<class T>
inline void function1D<T>::resize_virtual(int n)
{
  if (n>this->N0){
    std::cerr<<"Can not resize_vrtual over the limit of initialized memory!"<<std::endl;
  }
  this->N = n;
}

template<class T>
inline function1D<T>::function1D(const function1D& m)
{
  resize(m.N);
  std::copy(m.f,m.f+this->N,this->f);
}

template <class T>
inline function1D<T>& function1D<T>::operator=(const function1D<T>& m)
{
  resize(m.N);
  std::copy(m.f,m.f+this->N,this->f);
  return *this;
}

template <class T>
inline function1D<T>& function1D<T>::copy_full(const function1D<T>& m)
{
  resize(m.N0);
  this->N = m.N;
  std::copy(m.f,m.f+this->N0,this->f);
  return *this;
}

// spline1D ////////////////////////////////////////////////////////////

template<class T>
inline spline1D<T>::~spline1D()
{
  if (this->f != NULL){
    delete[] this->f;
    this->f = NULL;
    delete[] f2;
    f2 = NULL;
    delete[] dxi;
    dxi = NULL;
  }
  this->N=0;
  this->N0=0;
}

template<class T>
inline void spline1D<T>::resize(int n)
{
  if (n>this->N0){
    if (this->f) delete[] this->f;
    if (f2) delete[] f2;
    if (dxi) delete[] dxi;
    this->f = new T[n];
    f2 = new T[n];
    dxi = new double[n];
    this->N0=n;
  }
  this->N = n;
}

template<class T>
inline spline1D<T>::spline1D(int N_) : functionb<T>(N_)
{
  this->f = new T[N_];
  f2 = new T[N_];
  dxi = new double[N_];
}

template<class T>
inline spline1D<T>::spline1D(const spline1D& m) : functionb<T>(), f2(NULL), dxi(NULL)
{
  resize(m.N);
  std::copy(m.f,m.f+this->N,this->f);
  std::copy(m.f2,m.f2+this->N,f2);
  std::copy(m.dxi,m.dxi+this->N,dxi);
}

template <class T>
inline spline1D<T>& spline1D<T>::operator=(const spline1D<T>& m)
{
  resize(m.N);
  std::copy(m.f,m.f+this->N,this->f);
  std::copy(m.f2,m.f2+this->N,f2);
  std::copy(m.dxi,m.dxi+this->N,dxi);
  return *this;
}

template <class T>
T spline1D<T>::operator()(const intpar& ip) const
{
  int i= ip.i; double p = ip.p, q=1-ip.p;
  return q*this->f[i] + p*this->f[i+1] + dxi[i]*dxi[i]*(q*(q*q-1)*f2[i] + p*(p*p-1)*f2[i+1])/6.;
}

template <class T>
T spline1D<T>::df(const intpar& ip) const
{
  int i= ip.i; double p = ip.p, q=1-ip.p;
  return (this->f[i+1]-this->f[i])/dxi[i] + (1-3*q*q)*dxi[i]*f2[i]/6. - (1-3*p*p)*dxi[i]*f2[i+1]/6.;
}

template <class T>
T spline1D<T>::df(int i) const
{
  if (i==this->N-1){
    return (this->f[i]-this->f[i-1])/dxi[i-1] + dxi[i-1]*(f2[i-1]+2*f2[i])/6.;
  }else  return (this->f[i+1]-this->f[i])/dxi[i] - dxi[i]*(2*f2[i]+f2[i+1])/6.;
}

template <class T>
template <class mesh1D>
inline void spline1D<T>::splineIt(const mesh1D& om, const T& df0, const T& dfN)
{
  if (om.size()!=this->size()) std::cerr<<"Sizes of om and f are different in spline setup"<<std::endl;
  if (this->size()<2){ for (int i=0; i<this->size(); i++) f2[i]=0; return;}
  //  resize(om.size()); // Calling constructor to initialize memory
  //  std::copy(fu.f,fu.f+N,f);
  function1D<double> diag(om.size());
  function1D<T> offdiag(om.size()-1); // matrix is stored as diagonal values + offdiagonal values
  // Below, matrix and rhs is setup
  diag[0] = (om[1]-om[0])/3.;
  T dfu0 = (this->f[1]-this->f[0])/(om[1]-om[0]);
  f2[0] = dfu0-df0;
  for (int i=1; i<om.size()-1; i++){
    diag[i] = (om[i+1]-om[i-1])/3.;
    T dfu1 = (this->f[i+1]-this->f[i])/(om[i+1]-om[i]);
    f2[i] = dfu1-dfu0;
    dfu0 = dfu1;
  }
  diag[this->N-1] = (om[this->N-1]-om[this->N-2])/3.;
  f2[this->N-1] = dfN - (this->f[this->N-1]-this->f[this->N-2])/(om[this->N-1]-om[this->N-2]);
  for (int i=0; i<om.size()-1; i++) offdiag[i] = (om[i+1]-om[i])/6.;
  // The system of symmetric tridiagonal equations is solved by lapack
  int one=1, info=0;
  if (df0==std::numeric_limits<double>::max() || dfN==std::numeric_limits<double>::max()){
    int size = this->N-2;// natural splines	
    xptsv_(&size, &one, diag.MemPt()+1, offdiag.MemPt()+1, f2+1, &(this->N), &info);
    f2[0]=0; f2[this->N-1]=0;
  } else  xptsv_(&(this->N), &one, diag.MemPt(), offdiag.MemPt(), f2, &(this->N), &info);

  if (info!=0) std::cerr<<"dptsv return an error "<<info<<std::endl;
  // Setup of other necessary information for doing splines.
  for (int i=0; i<om.size()-1; i++) dxi[i] = (om[i+1]-om[i]);
}

template <class T>
inline T spline1D<T>::integrate()
{
  T sum=0;
  for (int i=0; i<this->N-1; i++) sum += 0.5*dxi[i]*(this->f[i+1]+this->f[i]-(f2[i+1]+f2[i])*dxi[i]*dxi[i]/12.);
  return sum;
}

template <class T>
template <class mesh1D>
inline complex<double> spline1D<T>::Fourier(double om, const mesh1D& xi) const
{
  complex<double> ii(0,1);
  complex<double> sum=0;
  complex<double> val;
  for (int i=0; i<this->N-1; i++) {
    double u = om*dxi[i], u2=u*u, u4=u2*u2;
    if (fabs(u)<1e-4){// Taylor expansion for small u
      val = 0.5*(1.+ii*(u/3.))*this->f[i]+0.5*(1.+ii*(2*u/3.))*this->f[i+1];
      val -= dxi[i]*dxi[i]/24.*(f2[i]*(1.+ii*7.*u/15.)+f2[i+1]*(1.+ii*8.*u/15.));
    }else{
      complex<double> exp(cos(u),sin(u));
      val  = (this->f[i]*(1.+ii*u-exp) + this->f[i+1]*((1.-ii*u)*exp-1.))/u2;
      val += dxi[i]*dxi[i]/(6.*u4)*(f2[i]*(exp*(6+u2)+2.*(u2-3.*ii*u-3.))+f2[i+1]*(6+u2+2.*exp*(u2+3.*ii*u-3.)));
    }
    sum += dxi[i]*complex<double>(cos(om*xi[i]),sin(om*xi[i]))*val;
    //sum += xi.Dh(i)*complex<double>(cos(om*xi[i]),sin(om*xi[i]))*val;
  }
  return sum;
}

template <class T>
template <class mesh1D>
inline dcomplex spline1D<T>::Fourier_(double om, const mesh1D& xi) const
{
  dcomplex ii(0,1);
  dcomplex sum=0;
  dcomplex val;
  for (int i=0; i<this->N-1; i++) {
    double u = om*dxi[i], u2=u*u, u4=u2*u2;
    if (fabs(u)<1e-4){// Taylor expansion for small u
      val = 0.5*(1.+ii*(u/3.))*this->f[i]+0.5*(1.+ii*(2*u/3.))*this->f[i+1];
      val -= dxi[i]*dxi[i]/24.*(f2[i]*(1.+ii*7.*u/15.)+f2[i+1]*(1.+ii*8.*u/15.));
    }else{
      dcomplex exp(cos(u),sin(u));
      val  = (this->f[i]*(1.+ii*u-exp) + this->f[i+1]*((1.-ii*u)*exp-1.))/u2;
      val += dxi[i]*dxi[i]/(6.*u4)*(f2[i]*(exp*(6+u2)+2.*(u2-3.*ii*u-3.))+f2[i+1]*(6+u2+2.*exp*(u2+3.*ii*u-3.)));
    }
    sum += dxi[i]*dcomplex(cos(om*xi[i]),sin(om*xi[i]))*val;
    //sum += xi.Dh(i)*dcomplex(cos(om*xi[i]),sin(om*xi[i]))*val;
  }
  return sum;
}

template <class T>
template <class mesh1D>
inline T spline1D<T>::Sum_iom(const mesh1D& om) const
{
  if (om.size()!=this->size()) {std::cerr<<"Sizes of iom and myself are different in spline1D::Sum_iom! Boiling out"<<std::endl; exit(1);}
  if (this->size()<1) return 0;
  if (this->size()<2) return this->f[0];
  //  if (size()<3) return (f[0]+f[1]);
  
  // trapezoid part
  T sum = this->f[0]*0.5*(om[1]-om[0]+1);
  for (int i=1; i<this->size()-1; i++){
    sum += this->f[i]*0.5*(om[i+1]-om[i-1]);
  }
  sum += this->f[this->N-1]*0.5*(om[this->N-1]-om[this->N-2]+1);
  // correction due to splines
  T corr=0;
  for (int i=0; i<this->size()-1; i++){
    corr -= (om[i+1]-om[i]+1)*(om[i+1]-om[i])*(om[i+1]-om[i]-1)*(f2[i]+f2[i+1])/24.;
  }
  return sum + corr;
}

template <class T>
std::ostream& operator<<(std::ostream& stream, const functionb<T>& f)
{
  int width = stream.width(); 
  for (int i=0; i<f.size(); i++) stream<<i<<" "<<std::setw(width)<<f[i]<<std::endl;
  return stream;
}

// funProxy ///////////////////////////////////////////////////////////////
template <class T>
inline void funProxy<T>::Initialize(int N_, T* f_)
{
  this->N = this->N0 = N_;// this->f = f_;
  this->f = new (f_) T[N_];
}

template <class T>
inline void funProxy<T>::ReInitialize(int N0_, int N_, T* f_)
{
  this->N0 = N0_;
  this->N = N_;
  this->f = f_;
}

template <class T>
inline void funProxy<T>::resize(int N_)
{
  if (N_>this->N0) std::cerr<<"Can't resize funProxy, to small funProxy!"<<std::endl;
  else this->N=N_;
}

template <class T>
inline funProxy<T>& funProxy<T>::operator=(const functionb<T>& m)
{
  resize(m.size());
  std::copy(m.MemPt(),m.MemPt()+this->N,this->f);
  return *this;
}


// function2D ////////////////////////////////////////////////////////////
template<class T>
function2D<T>::function2D(int N_, int Nd_) : N0(N_), Nd0(Nd_), N(N_), Nd(Nd_) 
{
  memory = operator new (sizeof(funProxy<T>)*N0+sizeof(T)*Nd0*N0+HPoffset);
  Assert(memory!=NULL,"Out of memory");
  f = new (memory) funProxy<T>[N0];
  
  int offset = sizeof(funProxy<T>)*N0+HPoffset;
  data = reinterpret_cast<T*>(static_cast<char*>(memory)+offset);
  
  // First call the default constructor on all newly allocated data
  for (int i=0; i<N0; i++)
    for (int j=0; j<Nd0; j++)
      new(&data[i*Nd0+j])T;
  
  for (int i=0; i<N0; i++) f[i].Initialize(Nd0,data+i*Nd0);
}

template<class T>
function2D<T>::~function2D()
{
  if (memory!=NULL){
    //cout << "And is now executing" <<endl;
    for (int i=0; i<N0; i++) f[i].~funProxy<T>();
    operator delete(memory);
    memory = NULL;
  }
  N0=Nd0=N=Nd=0;
}

template<class T>
void function2D<T>::ReleaseMemory()
{
  //cout << "Release memory of function2D called with memory="<< memory<<endl;
  if (memory!=NULL){
    //cout << "And is now executing" <<endl;
    for (int i=0; i<N0; i++) f[i].~funProxy<T>();
    operator delete(memory);
    memory = NULL;
  }
  N0=Nd0=N=Nd=0;
}


template <class T>
inline function2D<T>& function2D<T>::operator=(const function2D& m)
{
  //if (N0-1<m.N && m.N<=N0 && Nd0-1<m.Nd && m.Nd<=Nd0){
  if (m.N<=N0 && m.Nd<=Nd0){
    N = m.N; Nd = m.Nd;
    for (int i=0; i<N; i++) memcpy(f[i].f, m.f[i].f, sizeof(T)*Nd);
  } else{
    int msize = sizeof(funProxy<T>)*m.N0+sizeof(T)*m.Nd0*m.N0+HPoffset;
    operator delete(memory);
    memory = operator new (msize);
    Assert(memory!=NULL,"Out of memory");
    N0 = m.N0;
    Nd0 = m.Nd0;
    N = m.N;
    Nd = m.Nd;
    f = new (memory) funProxy<T>[N0];
    int offset = sizeof(funProxy<T>)*N0+HPoffset;
    data = reinterpret_cast<T*>(static_cast<char*>(memory)+offset);
    for (int i=0; i<N0; i++) f[i].Initialize(Nd0, data+i*Nd0);
    memcpy(data, m.data, sizeof(T)*Nd0*N0);
  }
  return *this;
}

template <class T>
void function2D<T>::SetEqual(const function2D& m)
{
  if (N0-1<m.N && m.N<=N0 && Nd0-1<m.Nd && m.Nd<=Nd0){
    //if (m.N<=N0 && m.Nd<=Nd0){
    N = m.N; Nd = m.Nd;
    for (int i=0; i<N; i++) memcpy(f[i].f, m.f[i].f, sizeof(T)*Nd);
  } else{
    int msize = sizeof(funProxy<T>)*m.N+sizeof(T)*m.Nd*m.N+HPoffset;
    if (memory!=NULL) operator delete(memory);
    memory = operator new (msize);
    Assert(memory!=NULL,"Out of memory");
    N0 = m.N;
    Nd0 = m.Nd;
    N = m.N;
    Nd = m.Nd;
    f = new (memory) funProxy<T>[N0];
    int offset = sizeof(funProxy<T>)*N0+HPoffset;
    data = reinterpret_cast<T*>(static_cast<char*>(memory)+offset);
    // First call the default constructor on all newly allocated data
    for (int i=0; i<N0; i++)
      for (int j=0; j<Nd0; j++)
	new(&data[i*Nd0+j])T;
    for (int i=0; i<N0; i++)
      f[i].Initialize(Nd0, data+i*Nd0);
    
    for (int i=0; i<N0; i++)
      for (int j=0; j<Nd0; j++)
	data[i*Nd0+j] = m.data[i*m.Nd0+j];
  }
}

template <class T>
inline function2D<T>::function2D(const function2D& m) : memory(NULL)
{
  int msize = sizeof(funProxy<T>)*m.N0+sizeof(T)*m.Nd0*m.N0+HPoffset;
  operator delete(memory);
  memory = operator new (msize);
  Assert(memory!=NULL,"Out of memory");
  N0 = m.N0;
  Nd0 = m.Nd0;
  N = m.N;
  Nd = m.Nd;
  f = new (memory) funProxy<T>[N0];
  int offset = sizeof(funProxy<T>)*N0+HPoffset;
  data = reinterpret_cast<T*>(static_cast<char*>(memory)+offset);
  for (int i=0; i<N0; i++) f[i].Initialize(Nd0, data+i*Nd0);
  //memcpy(data, m.data, sizeof(T)*Nd0*N0);
  for (int i=0; i<N0; i++)
    for (int j=0; j<Nd0; j++)
      data[i*Nd0+j] = m.data[i*m.Nd0+j];
}

template <class T>
inline void function2D<T>::resize(int N_, int Nd_)
{
  if (N_>N0 || Nd_>Nd0){
    //    clog<<"Deleting function2D and resizing from "<<N0<<" "<<Nd0<<" to "<<N_<<" "<<Nd_<<std::endl;
    int msize = sizeof(funProxy<T>)*N_ +sizeof(T)*Nd_*N_+HPoffset;
    operator delete(memory);
    memory = operator new (msize);
    Assert(memory!=NULL,"Out of memory");
    N = N0 = N_; Nd = Nd0 = Nd_;
    f = new (memory) funProxy<T>[N];
    int offset = sizeof(funProxy<T>)*N+HPoffset;
    data = reinterpret_cast<T*>(static_cast<char*>(memory)+offset);
    for (int i=0; i<N; i++) f[i].Initialize(Nd, data+i*Nd);
    // Next call the default constructor on all newly allocated data
    for (int i=0; i<N_; i++)
      for (int j=0; j<Nd_; j++)
	new(&data[i*Nd_+j])T; 
  } else{
    N = N_; Nd = Nd_;
  }
}

template <class T>
inline void function2D<T>::resize_virtual(int N_, int Nd_)
{
  if (N_>N0 || Nd_>Nd0){
    cerr<<"Can not resize virtually beyond boundary "<<N_<<">"<<N0<<" || "<<Nd_<<">"<<Nd0<<endl;
  } else{
    N = N_; Nd = Nd_;
  }
}

template <class T>
inline void function2D<T>::resize_preserve(int N_, int Nd_)
{
  if (N_>N0 || Nd_>Nd0){
    int msize = sizeof(funProxy<T>)*N_ +sizeof(T)*Nd_*N_+HPoffset;
    void* new_memory = operator new (msize);
    Assert(new_memory!=NULL,"Out of memory");

    int offset = sizeof(funProxy<T>)*N_+HPoffset;
    T* new_data = reinterpret_cast<T*>(static_cast<char*>(new_memory)+offset);

    const funProxy<T>* f_old = f;
    f = new (new_memory) funProxy<T>[N_];

    // First call the default constructor on all newly allocated data
    for (int i=0; i<N_; i++)
      for (int j=0; j<Nd_; j++)
	new(&new_data[i*Nd_+j])T;
    
    for (int i=0; i<N_; i++) f[i].ReInitialize(Nd_, f_old[i].N, new_data+i*Nd_);

    // The old part of the matrix needs to be copied by calling operator=
    for (int i=0; i<N0; i++)
      for (int j=0; j<Nd0; j++)
	new_data[i*Nd_+j] = data[i*Nd0+j];
    
    data = new_data;
    //N = N_;
    N0 = N_;
    //Nd = Nd_;
    Nd0 = Nd_;
    operator delete(memory);
    memory = new_memory;
  }
  //else{
  //  N = N_; Nd = Nd_;
  //}
}

template <class T>
inline void function2D<T>::resize_clear(int N_, int Nd_)
{
  if ((N_>N0 || Nd_>Nd0) || (N_<N0-1 || Nd_<Nd0-1)){
    //    clog<<"Deleting function2D and resizing from "<<N0<<" "<<Nd0<<" to "<<N_<<" "<<Nd_<<std::endl;
    int msize = sizeof(funProxy<T>)*N_ +sizeof(T)*Nd_*N_+HPoffset;
    operator delete(memory);
    memory = operator new (msize);
    Assert(memory!=NULL,"Out of memory");
    N = N0 = N_; Nd = Nd0 = Nd_;
    f = new (memory) funProxy<T>[N];
    int offset = sizeof(funProxy<T>)*N+HPoffset;
    data = reinterpret_cast<T*>(static_cast<char*>(memory)+offset);
    for (int i=0; i<N; i++) f[i].Initialize(Nd, data+i*Nd);
    // Next call the default constructor on all newly allocated data
    for (int i=0; i<N_; i++)
      for (int j=0; j<Nd_; j++)
	new(&data[i*Nd_+j])T; 
  } else{
    N = N_; Nd = Nd_;
  }
}


template <class T>
template <class functor>
inline void function2D<T>::Set(functor& functn)
{
  for (int i=0; i<N; i++)
    for (int j=0; j<Nd; j++)
      f[i][j] = functn(f[i][j]);
}

template <class T>
template <class functor, class W>
inline void function2D<T>::Set(functor& functn, const function2D<W>& Um)
{
  if (N0<Um.size_N() || Nd0<Um.size_Nd()){
    std::cerr << "To small function2D for Set functor!" << std::endl;
    return;
  }
  N = Um.size_N();
  if (Nd != Um.size_Nd()){
    Nd = Um.size_Nd();
    for (int i=0; i<N; i++) f[i].resize(Nd);
  }
  for (int i=0; i<N; i++)
    for (int j=0; j<Nd; j++)
      f[i][j] = functn(Um[i][j]);
}

template <class T>
inline function2D<T>& function2D<T>::minus(const function2D& B, const function2D& C)
{
  const T* __restrict__ _B_ = B.MemPt();
  const T* __restrict__ _C_ = C.MemPt();
  T* __restrict__ _data_ = data;
  for (int i=0; i<N; i++){
    int offset = i*Nd0;
    for (int j=0; j<Nd; j++) _data_[j+offset] = _B_[j+offset] - _C_[j+offset];
  }
  return *this;
}
template <class T>
inline function2D<T>& function2D<T>::plus(const function2D& B, const function2D& C)
{
  const T* __restrict__ _B_ = B.MemPt();
  const T* __restrict__ _C_ = C.MemPt();
  T* __restrict__ _data_ = data;
  for (int i=0; i<N; i++){
    int offset = i*Nd0;
    for (int j=0; j<Nd; j++) _data_[j+offset] = _B_[j+offset] + _C_[j+offset];
  }
  return *this;
}

// template <class T>
// inline void function2D<T>::Product(const function2D& A, const function2D& B, const T& alpha=1, const T& beta=0)
// {
//   if (A.Nd != B.Nd || !B.N || !A.N || !A.Nd || !B.Nd || N0<A.N || Nd0<B.N)
//     std::cerr << " Matrix sizes not correct" << std::endl;
//   xgemm("T", "N", B.N, A.N, B.Nd, alpha, B.MemPt(), B.Nd0, A.MemPt(), A.Nd0, beta, MemPt(), Nd0);
//   N = A.N; Nd = B.N;  
// }

template <class T>
inline void function2D<T>::Product(const function2D& A, const function2D& B, int start, int end, const T& alpha, const T& beta)
{
  if (end==-1){end=B.N;}
  if (A.Nd != B.Nd || !B.N || !A.N || !A.Nd || !B.Nd || N0<A.N || Nd0<end-start || end>B.N || start>=B.N) std::cerr << " Matrix sizes not correct" << std::endl;
  //cout<<"A.N="<<A.N<<" B.Nd="<<B.Nd<<" B.Nd0="<<B.Nd0<<" A.Nd0="<<A.Nd0<<" Nd0="<<Nd0<<endl;
  xgemm("T", "N", end-start, A.N, B.Nd, alpha, B[start].MemPt(), B.Nd0, A.MemPt(), A.Nd0, beta, MemPt(), Nd0);
  N = A.N; Nd = end-start;  
}

template <class T>
inline void function2D<T>::TProduct(const function2D& A, const function2D& B, const T& alpha, const T& beta)
{
  if (A.Nd != B.Nd || !B.N || !A.N || !A.Nd || !B.Nd || N0<B.N || Nd0<A.N)
    std::cerr << " Matrix sizes not correct" << std::endl;
//    clog<<"A.MemPt()="<<A.MemPt()<<" B.MemPt()="<<B.MemPt()<<" C.MemPt()="<<MemPt()<<std::endl;
//    clog<<"A.End()="<<A.MemPt()+A.N0*A.Nd0<<" B.End()="<<B.MemPt()+B.N0*B.Nd0<<" C.Endt()="<<MemPt()+N0*Nd0<<std::endl;
  xgemm("T", "N", A.N, B.N, A.Nd, alpha, A.MemPt(), A.Nd0, B.MemPt(), B.Nd0, beta, MemPt(), Nd0);
  N = B.N; Nd = A.N;
}

template <class T>
inline void function2D<T>::MProduct(const function2D& A, const function2D& B, const T& alpha, const T& beta)
{
  if (A.Nd != B.N || !B.Nd || !A.N || !A.Nd || !B.N || N0<A.N || Nd0<B.Nd)
    std::cerr << " Matrix sizes not correct" << std::endl;
  xgemm("N", "N", B.Nd, A.N, B.N, alpha, B.MemPt(), B.Nd0, A.MemPt(), A.Nd0, beta, MemPt(), Nd0);
  N = A.N; Nd = B.Nd;
}
template <class T>
inline void function2D<T>::SMProduct(const function2D& A, const function2D& B, const T& alpha, const T& beta)
{
//   if (A.Nd != B.N || !B.Nd || !A.N || !A.Nd || !B.N || N0<A.N || Nd0<B.Nd)
//     std::cerr << " Matrix sizes not correct" << std::endl;
  if (A.Nd==1 && A.N==1 && B.Nd==1)
    f[0].f[0] = A(0,0)*B(0,0);
  else
    xgemm("N", "N", B.Nd, A.N, B.N, alpha, B.MemPt(), B.Nd0, A.MemPt(), A.Nd0, beta, MemPt(), Nd0);
  N = A.N; Nd = B.Nd;
}
// template <class T>
// inline void function2D<T>::TMProduct(const function2D& A, const function2D& B, const T& alpha=1, const T& beta=0)
// {
//   if (A.N != B.N || !B.Nd || !A.Nd || !A.N || !B.N || N0<A.Nd || Nd0<B.Nd)
//     std::cerr << " Matrix sizes not correct" << std::endl;
//   xgemm("N", "T", B.Nd, A.Nd, B.N, alpha, B.MemPt(), B.Nd0, A.MemPt(), A.N0, beta, MemPt(), Nd0);
//   N = A.Nd; Nd = B.Nd;
// }

template <class T>
inline void function2D<T>::SymmProduct(const function2D& A, const function2D& B, const T& alpha, const T& beta)
{
  if (A.Nd != B.Nd || !B.N || !A.N || !A.Nd || !B.Nd || N0<A.N || Nd0<B.N)
    std::cerr << " Matrix sizes not correct" << std::endl;
  xsymm("R", "L", A.N, B.N, alpha, A.MemPt(), A.Nd0, B.MemPt(), B.Nd0, beta, MemPt(), Nd0);
  N = B.N; Nd = A.N;  
}

template <class T>
inline void function2D<T>::TSymmProduct(const function2D& A, const function2D& B, const T& alpha, const T& beta)
{
  if (A.Nd != B.Nd || !B.N || !A.N || !A.Nd || !B.Nd || N0<A.N || Nd0<B.N)
    std::cerr << " Matrix sizes not correct" << std::endl;
  xsymm("L", "L", A.N, B.N, alpha, A.MemPt(), A.Nd0, B.MemPt(), B.Nd0, beta, MemPt(), Nd0);
  N = B.N; Nd = A.N;  
}

template <class T>
inline void function2D<T>::DotProduct(const function2D& A, const function2D& B, const T& alpha, const T& beta)
{
  if (B.N != A.Nd || !B.N || !A.N || !A.Nd || !B.Nd || Nd0<B.Nd || N0<A.N)
    std::cerr << " Matrix sizes not correct" << std::endl;
  //  xgemm("T", "T", A.N, B.Nd, A.Nd, alpha, A.MemPt(), A.Nd0, B.MemPt(), B.Nd0, beta, MemPt(), Nd0);
  xgemm("N", "N", B.Nd, A.N, B.N, alpha, B.MemPt(), B.Nd0, A.MemPt(), A.Nd0, beta, MemPt(), Nd0);
  Nd = B.Nd; N = A.N;
}

template <class T>
inline function2D<T>& function2D<T>::operator+=(double x)
{
  if (N!=Nd || !Nd || !N) {
    std::cerr << "Can't add number to non-square matrix!" << std::endl;
    return *this;
  }
  for (int i=0; i<Nd; i++) f[i][i] += x;
  return *this;
}

template <class T>
inline function2D<T>& function2D<T>::operator+=(const function2D& m)
{
  if (N!=m.N || Nd!=m.Nd) {
    std::cerr << "Can't sum different matrices!" << std::endl;
    return *this;
  }
  for (int i=0; i<N; i++)
    for (int j=0; j<Nd; j++)
      f[i][j] += m[i][j];
  
  return *this;
}

template <class T>
inline void function2D<T>::AddPart(const function2D& m, double x)
{
  for (int i=0; i<N; i++)
    for (int j=0; j<Nd; j++)
      f[i][j] += m[i][j]*x;
}

template <class T>
inline function2D<T>& function2D<T>::operator-=(double x)
{
  if (N!=Nd || !N || !Nd) {
    std::cerr << "Can't add number to non-square matrix!" << std::endl;
    return *this;
  }
  for (int i=0; i<Nd; i++) f[i][i] -= x;
  return *this;
}

template <class T>
inline function2D<T>& function2D<T>::operator-=(const function2D& m)
{
  if (N!=m.N || Nd!=m.Nd) {
    std::cerr << "Can't sum different matrices!" << std::endl;
    return *this;
  }
  for (int i=0; i<N; i++)
    for (int j=0; j<Nd; j++)
      f[i][j] -= m[i][j];
  
  return *this;
}
template <class T>
inline function2D<T>& function2D<T>::operator=(const T& u)
{
  for (int i=0; i<N; i++) for (int j=0; j<Nd; j++) f[i].f[j]=u;
  return *this;
}
template <class T>
inline function2D<T>& function2D<T>::operator*=(const T& x)
{
  for (int i=0; i<N; i++) for (int j=0; j<Nd; j++) f[i][j] *= x;
  return *this;
}

template <class T>
template <class meshx, class functiond>
void function2D<T>::KramarsKronig(const meshx& om, const functiond& logo, functiond& sum, functiond& odvSig)
{
  sum.resize(size_Nd()); odvSig.resize(size_Nd());
  
  if (logo.size()!=om.size()||om.size()!=N)
    std::cerr<<"Functions not of equal size in KramarsKronig"<<std::endl;

  for (int i=0; i<om.size(); i++){
    const double omom2 = (i>0) ? om.Delta(i-1) : 0.0;
    const int ip1 = (i<om.size()-1) ? i+1 : i, im1 = (i>0) ? i-1 : i;

    for (int l=0; l<size_Nd(); l++)
      odvSig[l] = 0.5*(om.Delta(i)*(f[ip1][l].imag()-f[i][l].imag()) + omom2*(f[i][l].imag()-f[im1][l].imag()));
    
    for (int l=0; l<size_Nd(); l++) sum[l] = 0;
    
    T *g0 = f[i].MemPt();
    
    for (int j=0; j<om.size(); j++){
      double dh = om.Dh(j), dhdom = dh/(om[j]-om[i]);
      T *g = f[j].MemPt();
      
      if (i!=j)
	for (int l=0; l<size_Nd(); l++)  sum[l] += (g[l].imag()-g0[l].imag())*dhdom;
      else
	for (int l=0; l<size_Nd(); l++)  sum[l] += odvSig[l]*dh;
    }
    
    for (int l=0; l<size_Nd(); l++) f[i][l].real()=(sum[l]+g0[l].imag()*logo[i])/M_PI;
  }
}

template <class T>
std::ostream& operator<<(std::ostream& stream, const function2D<T>& f)
{
  int width = stream.width(); 
  for (int i=0; i<f.size_N(); i++){
    for (int j=0; j<f.size_Nd(); j++)
      stream<<std::setw(width)<<f[i][j]<<" ";
    stream<<std::endl;
  }
  return stream;
}

template <class T>
void functionb<T>::xgemv(const string& TN, const function2D<T>& A, const functionb<T>& x, double alpha, double beta)
{
  ::xgemv(TN, A.size_Nd(), A.size_N(), alpha, A.MemPt(), A.lda(), x.MemPt(), beta, MemPt());
}


template<class meshx, class functionx>
void print(std::ostream& stream, const meshx& om, const functionx& f, int width)
{
  if (om.size()!=f.size()) std::cerr<<"Can't print objectc of different size!"<<std::endl;
  for (int i=0; i<om.size(); i++)
    stream <<std::setw(width)<<om[i]<<std::setw(width)<<f[i]<<std::endl;
}
template<class meshx, class functionx, class functiony>
void print(std::ostream& stream, const meshx& om, const functionx& f1, const functiony& f2, int width)
{
  if (om.size()!=f1.size() || om.size()!=f2.size()) std::cerr<<"Can't print objectc of different size!"<<std::endl;
  for (int i=0; i<om.size(); i++)
    stream <<std::setw(width)<<om[i]<<std::setw(width)<<f1[i]<<std::setw(width)<<f2[i]<<std::endl;
}
template<class meshx, class functionx, class functiony, class functionz>
void print(std::ostream& stream, const meshx& om, const functionx& f1, const functiony& f2, const functionz& f3, int width)
{
  if (om.size()!=f1.size() || om.size()!=f2.size() || om.size()!=f3.size()) std::cerr<<"Can't print objectc of different size!"<<std::endl;
  for (int i=0; i<om.size(); i++)
    stream <<std::setw(width)<<om[i]<<std::setw(width)<<f1[i]<<std::setw(width)<<f2[i]<<std::setw(width)<<f3[i]<<std::endl;
}
template<class meshx, class functionx, class functiony, class functionz, class functionw>
void print(std::ostream& stream, const meshx& om, const functionx& f1, const functiony& f2, const functionz& f3, const functionw& f4, int width)
{
  if (om.size()!=f1.size() || om.size()!=f2.size() || om.size()!=f3.size() || om.size()!=f4.size())
    std::cerr<<"Can't print objectc of different size!"<<std::endl;
  for (int i=0; i<om.size(); i++)
    stream <<std::setw(width)<<om[i]<<std::setw(width)<<f1[i]<<std::setw(width)<<f2[i]<<std::setw(width)<<f3[i]<<std::setw(width)<<f4[i]<<std::endl;
}

template <class T>
inline T scalar_product(const functionb<T>& f1, const functionb<T>& f2)
{
  static const int incr=1;
  Assert(f1.size()==f2.size(),"Sizes not the same in scalar_product!");
  return ddot_(&f1.N, f1.f, &incr, f2.f, &incr);
}

template <class funProxy>
inline int SolveSOLA(function2D<funProxy>& A, function2D<funProxy>& B, function1D<int>& ipiv)
{
  return xgesv(A.size_Nd(), B.size_N(), A.MemPt(), A.fullsize_Nd(), ipiv.MemPt(), B.MemPt(), B.fullsize_Nd());
}
template <class funProxy>
inline int Inverse(function2D<funProxy>& A, function2D<funProxy>& B, function1D<int>& ipiv)
{
  if (A.size_N()!=A.size_Nd()) {std::cerr<<"Can not invert nonquadratic matrix! "<<std::endl;return 1;}
  B.resize(A.size_N(),A.size_Nd());
  for (int i=0; i<B.size_N(); i++){
    for (int j=0; j<B.size_Nd(); j++)
      B[i][j] = 0.0;
    B[i][i] = 1.0;
  }
  int rr = SolveSOLA(A,B,ipiv);
  for (int i=0; i<A.size_N(); i++)
    for (int j=0; j<A.size_Nd(); j++)
      A(i,j) = B(j,i);
  for (int i=0; i<A.size_N(); i++)
    for (int j=0; j<A.size_Nd(); j++)
      B(i,j) = A(i,j);
  return rr;
}
 
template <class funProxy, class T>
inline int SolveSOLA(function2D<funProxy>& A, function1D<T>& B, function1D<int>& ipiv)
{
  return xgesv(A.size_Nd(), 1, A.MemPt(), A.fullsize_Nd(), ipiv.MemPt(), B.MemPt(), B.fullsize());
}

template <class T>
int function2D<T>::Inverse(const function2D<T>& A)
{
  static function2D<T> temp;
  static function1D<int> ipiv;

  if (A.size_N()==1){
    f[0][0] = 1/A(0,0);
    return 0;
  }
  if (A.size_N()==2){
    T Det = A[0][0]*A[1][1]-A[1][0]*A[0][1];
    T g00 = A[1][1]/Det;
    T g01 = -A[0][1]/Det;
    T g10 = -A[1][0]/Det;
    T g11 = A[0][0]/Det;
    f[0][0] = g00;
    f[0][1] = g01;
    f[1][0] = g10;
    f[1][1] = g11;
    return 0;
  }
  ipiv.resize(A.size_N());

  temp.resize(A.size_Nd(),A.size_N());
  for (int i=0; i<N; i++)
    for (int j=0; j<Nd; j++)
      temp(i,j) = A(j,i); // temp = transposed(A)
  
  resize(A.size_N(), A.size_Nd());
  for (int i=0; i<N; i++){
    for (int j=0; j<Nd; j++)
      f[i][j] = 0.0;
    f[i][i] = 1.0;
  }
  return SolveSOLA(temp, *this, ipiv);
}

inline double logpart(double a, double b, double x)
{
  return log(fabs((b-x)/(x-a)));
}
template <class complexn>
inline complexn logpart(double a, double b, const complexn& x)
{
  return log((b-x)/(a-x));
}
template <class T, class meshx>
T KramarsKronig(const functionb<double>& fi, const meshx& om, const T& x, int i0, double S0)
{
  T sum=0;
  for (int j=0; j<i0-1; j++) sum += (fi[j]-S0)*om.Dh(j)/(om[j]-x);
  if (i0>0)                  sum += (fi[i0-1]-S0)*(om.Dh(i0-1)+0.5*om.Dh(i0))/(om[i0-1]-x);
  if (i0<om.size()-1)        sum += (fi[i0+1]-S0)*(om.Dh(i0+1)+0.5*om.Dh(i0))/(om[i0+1]-x);
  for (int j=i0+2; j<om.size(); j++) sum += (fi[j]-S0)*om.Dh(j)/(om[j]-x);
  if (x!=om.last() && x!=om[0]) sum += S0*logpart(om[0],om.last(),x);
  return sum/M_PI;
}

template <class fun2D, class fun1D>
inline int Eigensystem(fun2D& H, fun1D& E, fun2D& AL, fun2D& AR)
{
  static fun1D work(4*H.size_N());
  static function1D<double> rwork(2*H.size_N());

  xgeev(H.size_N(), H.MemPt(), H.fullsize_Nd(), E.MemPt(), AL.MemPt(), AL.fullsize_Nd(),
        AR.MemPt(), AR.fullsize_Nd(), work.MemPt(), work.size(), rwork.MemPt());
  return 0;
}

template <class fun2D, class fun1D>
inline int SymEigensystem(fun2D& H, fun1D& E)
{
  static fun1D work(4*H.size_N());
  xsyev(H.size_N(), H.MemPt(), H.fullsize_Nd(), E.MemPt(), work.MemPt(), work.size(), "V");
  return 0;
}

// inline int Eigenvalues(function2D<double>& H, function1D<dcomplex>& E)
// {
//   static function1D<double> work(4*H.size_N());
//   //  static function1D<double> rwork(2*H.size_N());
//   static function1D<double> wr(H.size_Nd()), wi(H.size_Nd());
//   work.resize(4*H.size_N()); wr.resize(H.size_Nd()); wi.resize(H.size_Nd());
//   xgeev_(H.size_N(), H.MemPt(), H.fullsize_Nd(), E.MemPt(), wr.MemPt(), wi.MemPt(), work.MemPt(), work.size());
// }

inline void dger(function2D<double>& A, double alpha, const functionb<double>& x, const functionb<double>& y)
{
  int one=1;
  int ny = y.size(), nx = x.size();
  int lda = A.fullsize_Nd();
  dger_( &ny, &nx, &alpha, y.MemPt(), &one, x.MemPt(), &one, A.MemPt(), &lda);
}
inline void dger(double* A, double alpha, const functionb<double>& x, const functionb<double>& y)
{
  int one=1;
  int ny = y.size(), nx = x.size();
  int lda = ny;
  dger_( &ny, &nx, &alpha, y.MemPt(), &one, x.MemPt(), &one, A, &lda);
}

inline double Det(const function2D<double>& A)
{
  function2D<double> H(A);
  static function1D<double> work(4*H.size_N());
  static function1D<double> wr(H.size_Nd()), wi(H.size_Nd());
  work.resize(4*H.size_N()); wr.resize(H.size_Nd()); wi.resize(H.size_Nd());
  xgeev_(H.size_N(), H.MemPt(), H.fullsize_Nd(), wr.MemPt(), wi.MemPt(), work.MemPt(), work.size());
  dcomplex res=1;
  for (int i=0; i<H.size_N(); i++) res *= dcomplex(wr[i],wi[i]);
  return res.real();
}

/*
inline double Det2(const function2D<double>& A)
{
  // solve Ax=B using A=UDU' factorization, D is placed in A
  // where A represents a qxq matrix, b a 1xq vector
  function1D<int> ipiv(A.size_Nd());
  function2D<double> D(A);
  int info = xgetrf(D.size_Nd(), D.MemPt(), D.fullsize_Nd(), ipiv.MemPt());
  if( info < 0 ) { cerr<<"In Det and dgetrf the "<<-info<<"-th argument had an illegal value"<<endl;}
  if( info > 0 ) { cerr<<"In Det and dgetrf the U("<<info<<") is exactly zero and might cause a problem."; }
  // compute the determinant det = det(A)
  // if ipiv[i] > 0, then D(i,i) is a 1x1 block diagonal
  // if ipiv[i] = ipiv[i-1] < 0, then D(i-1,i-1),
  // D(i-1,i), and D(i,i) form a 2x2 block diagonal
  double det = 1.0;
  for(int i = 0; i < A.size_N(); i++ ){
    if( ipiv[i] > 0){
      det *= D(i,i);
    } else if( i > 0 && ipiv[i] < 0 && ipiv[i-1] == ipiv[i]){
      det *= D(i,i)*D(i-1,i-1) - D(i-1,i)*D(i-1,i);
    }
  }
  return det;
}
*/
template <class container>
inline double Determinant(const container& A0)
{
  if (A0.size_Nd()!=A0.size_N()) {std::cerr<<"Can't compute determinant of nonquadratic matrix!"<<std::endl; return 0;}
  container A(A0);// copy constructor
  int info;
  int n = A.size_N();
  int lda = A.fullsize_Nd();
  function1D<int> ipiv(n);
  dgetrf_(&n, &n, A.MemPt(), &lda, ipiv.MemPt(), &info);
  if (info) {std::cerr<<"LU factorization complains : "<<info<<std::endl; return 0;}
  double det = 1;
  for (int i=0; i<n; i++) det *= ((ipiv[i]==i+1) ? 1 : -1) * A(i,i);
  return det;
}



inline int SymmEigensystem(function2D<dcomplex>& H, function1D<double>& E)
{
  static int N = H.size_N();
  static int lwork = 2*(2*N+N*N), lrwork = 2*(1 + 5*N + 2*N*N), liwork = 2*(3 + 5*N);
  static function1D<dcomplex> work(lwork);
  static function1D<double> rwork(lrwork);
  static function1D<int> iwork(liwork);
  xheevd(H.size_N(), H.MemPt(), H.fullsize_Nd(), E.MemPt(), work.MemPt(), work.size(), rwork.MemPt(), rwork.size(), iwork.MemPt(), iwork.size());
  return 0;
}

// Test for non-qadratic matrices
template <class T>
inline void function2D<T>::Product(const std::string& transa, const std::string& transb, const function2D& A, const function2D& B, const T& alpha, const T& beta)
{
  if (transa!="N" && transa!="T" && transa!="C"){std::cerr<<"Did not recognize your task. Specify how to multiply matrices in dgemm!"<<std::endl; return;}
  if(transa=="N" && transb=="N"){
    if (A.Nd != B.N || !B.Nd || !A.N || !A.Nd || !B.N || N0<A.N || Nd0<B.Nd)
      std::cerr << " Matrix sizes not correct" << std::endl;
    xgemm(transb, transa, B.Nd, A.N, B.N, alpha, B.MemPt(), B.Nd0, A.MemPt(), A.Nd0, beta, MemPt(), Nd0);
    N = A.N; Nd = B.Nd;
  }else if((transa=="T"||transa=="C") && transb=="N"){
    if (A.N != B.N || !B.Nd || !A.Nd || !A.N || !B.N || N0<A.Nd || Nd0<B.Nd){
      std::cerr << " Matrix sizes not correct" << std::endl;
      cout<<"A.N="<<A.N<<" B.N="<<B.N<<" A.Nd="<<A.Nd<<" B.Nd="<<B.Nd<<" N0="<<N0<<" Nd0="<<Nd0<<endl;
    }
    xgemm(transb, transa, B.Nd, A.Nd, B.N, alpha, B.MemPt(), B.Nd0, A.MemPt(), A.Nd0, beta, MemPt(), Nd0);
    N = A.Nd; Nd = B.Nd;
  }else if (transa=="N" && (transb=="T"||transb=="C")){
    if (A.Nd != B.Nd || !B.N || !A.N || !A.Nd || !B.Nd || N0<A.N || Nd0<B.N)
      std::cerr << " Matrix sizes not correct" << std::endl;
    xgemm(transb, transa, B.N, A.N, B.Nd, alpha, B.MemPt(), B.Nd0, A.MemPt(), A.Nd0, beta, MemPt(), Nd0);
    N = A.N; Nd = B.N;
  }else if ((transa=="T"||transa=="C") && (transb=="T"||transb=="C")){
    if (A.N != B.Nd || !B.N || !A.N || !A.Nd || !B.Nd || N0<A.Nd || Nd0<B.N)
      std::cerr << " Matrix sizes not correct" << std::endl;
    xgemm(transb, transa, B.N, A.Nd, B.Nd, alpha, B.MemPt(), B.Nd0, A.MemPt(), A.Nd0, beta, MemPt(), Nd0);
    N = A.Nd; Nd = B.N;
  }
}

template <int N, int K, int M>
inline void Multiply_(function2D<double>& C, const function2D<double>& A, const function2D<double>& B)
{
  C.MProduct(A,B);
}


template <class T>
inline void Multiply(function2D<T>& C, const function2D<T>& A, const function2D<T>& B);


template <class T>
class function2Dm : public function2D<T>{  
public:
  int a0, b0, a1, b1;
  //  function2Dm() : function2Dm<T>, a0(0), b0(0), a1(0), b1(0) {};
  explicit function2Dm(int a0_=0, int b0_=0, int a1_=0, int b1_=0) : 
    function2D<T>(b0_-a0_,b1_-a1_), a0(a0_), b0(b0_), a1(a1_), b1(b1_){};

  ~function2Dm(){};

  funProxy<T>& operator[](int i) {Assert3(i<b0 && i>=a0,"Out of range in function2Dm[] ", i, this->N0); return this->f[i-a0];}
  const funProxy<T>& operator[](int i) const {Assert3(i<b0 && i>=a0,"Out of range in function2Dm[] ", i, this->N0); return this->f[i-a0];}
  const T& operator()(int i, int j) const {Assert5(i<b0 && j<b1 && i>=a0 && j>=a1,"Out of range in function2Dm(i,j) ", i, j, this->N0, this->Nd0); return this->f[i-a0].f[j-a1];}
  T& operator()(int i, int j) {Assert5(i<b0 && j<b1 && i>=a0 && j>=a1,"Out of range in function2Dm(i,j) ", i, j, a0, a1); return this->f[i-a0].f[j-a1];}

  function2Dm& operator=(const function2Dm& m);
  function2Dm& operator=(const T& u);
  
  void resize(int a0_, int b0_, int a1_, int b1_)
  {
    b0 = b0_; a0 = a0_;
    b1 = b1_; a1 = a1_;
    function2D<T>::resize(b0-a0, b1-a1);
  }
};

template <class T>
inline function2Dm<T>& function2Dm<T>::operator=(const T& u)
{
  for (int i=0; i<this->N; i++) for (int j=0; j<this->Nd; j++) this->f[i].f[j]=u;
  return *this;
}
template <class T>
inline function2Dm<T>& function2Dm<T>::operator=(const function2Dm& m)
{
  if (m.N<=this->N0 && m.Nd<=this->Nd0){
    this->N = m.N; this->Nd = m.Nd;
    a0 = m.a0; a1 = m.a1; b0 = m.b0; b1 = m.b1;
    for (int i=0; i<this->N; i++) memcpy(this->f[i].f, m.f[i].f, sizeof(T)*this->Nd);
  } else{
    int msize = sizeof(funProxy<T>)*m.N0+sizeof(T)*m.Nd0*m.N0+HPoffset;
    operator delete(this->memory);
    this->memory = operator new (msize);
    Assert(this->memory!=NULL,"Out of memory");
    this->N0 = m.N0;
    this->Nd0 = m.Nd0;
    this->N = m.N;
    this->Nd = m.Nd;
    a0 = m.a0; a1 = m.a1; b0 = m.b0; b1 = m.b1;
    this->f = new (this->memory) funProxy<T>[this->N0];
    int offset = sizeof(funProxy<T>)*this->N0+HPoffset;
    this->data = reinterpret_cast<T*>(static_cast<char*>(this->memory)+offset);
    for (int i=0; i<this->N0; i++) this->f[i].Initialize(this->Nd0, this->data+i*this->Nd0);
    memcpy(this->data, m.data, sizeof(T)*this->Nd0*this->N0);
  }
  return *this;
}


template <class T>
class function3D{
  T* f;
public:
  int N0, N1, N2;
  function3D() : N0(0), N1(0), N2(0), f(NULL) {};
  function3D(int N0_, int N1_, int N2_) : f(NULL), N0(N0_), N1(N1_), N2(N2_)
  {
    f = new T[N0*N1*N2];
  }
  void resize(int N0_, int N1_, int N2_){
    int size = N0_*N1_*N2_;
    if (size > N0*N1*N2){// old storage too small
      delete[] f;
      f = new T[N0_*N1_*N2_];
      if (f==NULL) std::cerr<<"Could not allocate memory for function3D of size "<<size*sizeof(T)/(1024*1024)<<"MB"<<std::endl;
    }
    N0 = N0_;
    N1 = N1_;
    N2 = N2_;
  }
  ~function3D()
  {
    if (f!=NULL) delete[] f;
    f=NULL;
  }
  T operator()(int i0, int i1, int i2) const
  {
    return f[i2 + N2*(i1 + N1*i0)];
  }
  T& operator()(int i0, int i1, int i2)
  {
    Assert(i0<N0 && i0>=0 && i1<N1 && i1>=0 && i2<N2 && i2>=0,"Out of range in function3D[] ");
    return f[i2 + N2*(i1 + N1*i0)];
  }
  int ind(int i0, int i1, int i2) const
  {return i2 + N2*(i1 + N1*i0);}
  double* operator[](int i0)
  { return &f[N2*N1*i0]; }
  const double* operator[](int i0) const
  { return &f[N2*N1*i0]; }
  function3D& operator*=(const T& c);
  function3D& operator=(const T& c);
};
template <class T>
function3D<T>& function3D<T>::operator*=(const T& c)
{
  for (int i0=0; i0<N0; i0++)
    for (int i1=0; i1<N1; i1++){
      T* __restrict__ f = &(operator()(i0,i1,0));
      for (int i2=0; i2<N2; i2++) f[i2]*=c;
    }
  return *this;
};
template <class T>
function3D<T>& function3D<T>::operator=(const T& c)
{
  for (int i0=0; i0<N0; i0++)
    for (int i1=0; i1<N1; i1++){
      T* __restrict__ f = &(operator()(i0,i1,0));
      for (int i2=0; i2<N2; i2++) f[i2]=c;
    }
  return *this;
};

template <class T>
class function4D{
  T* f;
public:
  int N0, N1, N2, N3;
  function4D() : N0(0), N1(0), N2(0), N3(0), f(NULL) {};
  function4D(int N0_, int N1_, int N2_, int N3_) : f(NULL), N0(N0_), N1(N1_), N2(N2_), N3(N3_)
  {
    f = new T[N0*N1*N2*N3];
  }
  void resize(int N0_, int N1_, int N2_, int N3_){
    int size = N0_*N1_*N2_*N3_;
    if (size > N0*N1*N2*N3){// old storage too small
      delete[] f;
      f = new T[N0_*N1_*N2_*N3_];
      if (f==NULL) std::cerr<<"Could not allocate memory for function4D of size "<<size*sizeof(T)/(1024*1024)<<"MB"<<std::endl;
    }
    N0 = N0_;
    N1 = N1_;
    N2 = N2_;
    N3 = N3_;
  }
  ~function4D()
  {
    if (f!=NULL) delete[] f;
    f=NULL;
  }
  T operator()(int i0, int i1, int i2, int i3) const
  {
    return f[i3 + N3*(i2 + N2*(i1 + N1*i0))];
  }
  T& operator()(int i0, int i1, int i2, int i3)
  {
    Assert(i0<N0 && i0>=0 && i1<N1 && i1>=0 && i2<N2 && i2>=0 && i3<N3 && i3>=0, "Out of range in function4D[] ");
    return f[i3 + N3*(i2 + N2*(i1 + N1*i0))];
  }
  int ind(int i0, int i1, int i2, int i3) const
  {return i3 + N3*(i2 + N2*(i1 + N1*i0));}
  
  function4D& operator*=(const T& c);
  function4D& operator=(const T& c);
};
template <class T>
function4D<T>& function4D<T>::operator*=(const T& c)
{
  for (int i0=0; i0<N0; i0++)
    for (int i1=0; i1<N1; i1++)
      for (int i2=0; i2<N2; i2++){
	T* __restrict__ fc = &(operator()(i0,i1,i2,0));
	for (int i3=0; i3<N3; i3++) fc[i3]*=c;
      }
  return *this;
};
template <class T>
function4D<T>& function4D<T>::operator=(const T& c)
{
  for (int i0=0; i0<N0; i0++)
    for (int i1=0; i1<N1; i1++)
      for (int i2=0; i2<N2; i2++){
	T* __restrict__ fc = &(operator()(i0,i1,i2,0));
	for (int i3=0; i3<N3; i3++) fc[i3]=c;
      }
  return *this;
};


template <class T>
class function5D{
  T* f;
public:
  int N0, N1, N2, N3, N4;
  function5D(int N0_=0, int N1_=0, int N2_=0, int N3_=0, int N4_=0) : f(NULL), N0(N0_), N1(N1_), N2(N2_), N3(N3_), N4(N4_)
  {
    int size = N0*N1*N2*N3*N4;
    if (size!=0) f = new T[size];
  }
  ~function5D()
  {
    if (f!=NULL) delete[] f;
    f=NULL;
  }
  T operator()(int i0, int i1, int i2, int i3, int i4) const
  {
    return f[i4 + N4*(i3 + N3*(i2 + N2*(i1 + N1*i0)))];
  }
  T& operator()(int i0, int i1, int i2, int i3, int i4)
  {
    Assert(i0<N0 && i0>=0 && i1<N1 && i1>=0 && i2<N2 && i2>=0 && i3<N3 && i3>=0 && i4<N4 && i4>=0, "Out of range in function5D[] ");
    return f[i4 + N4*(i3 + N3*(i2 + N2*(i1 + N1*i0)))];
  }
  void resize(int N0_, int N1_, int N2_, int N3_, int N4_)
  {
    int size = N0_*N1_*N2_*N3_*N4_;
    if (size > N0*N1*N2*N3*N4){
      delete[] f;
      f = new T[size];
      if (f==NULL) std::cerr<<"Could not allocate memory for function5D of size "<<size*sizeof(T)/(1024*1024)<<"MB"<<std::endl;
    }
    N0 = N0_; N1 = N1_; N2 = N2_; N3 = N3_; N4 = N4_;
  }
  function5D& operator*=(const T& c);
  int ind(int i0, int i1, int i2, int i3, int i4) const
  {return i4 + N4*(i3 + N3*(i2 + N2*(i1 + N1*i0)));}
  function5D& operator=(const T& c);
};

template <class T>
function5D<T>& function5D<T>::operator*=(const T& c)
{
  for (int i0=0; i0<N0; i0++)
    for (int i1=0; i1<N1; i1++)
      for (int i2=0; i2<N2; i2++)
	for (int i3=0; i3<N3; i3++){
	  T* __restrict__ f = &(operator()(i0,i1,i2,i3,0));
	  for (int i4=0; i4<N4; i4++) f[i4]*=c;
	}
  return *this;
};

template <class T>
function5D<T>& function5D<T>::operator=(const T& c)
{
  int size = N0*N1*N2*N3*N4;
  for (int i=0; i<size; i++) f[i]=c;
  /*
  for (int i0=0; i0<N0; i0++)
    for (int i1=0; i1<N1; i1++)
      for (int i2=0; i2<N2; i2++)
	for (int i3=0; i3<N3; i3++){
	  T* __restrict__ f = &(operator()(i0,i1,i2,i3,0));
	  for (int i4=0; i4<N4; i4++) f[i4]=c;
	}
  */
  return *this;
};

// template <class T>
// class function5Dm{
//   T* f;
// public:
//   int a0, b0, a1, b1, a2, b2, a3, b3, a4, b4;//N0, N1, N2, N3, N4;
//   function5D(int b0_, int b1_, int b2_, int b3_, int b4_,
// 	     int a0_=0, int a1_=0, int a2_=0, int a3_=0, int a4_=0) :
//     b0(b0_), b1(b1_), b2(b2_), b3(b3_), b4(b4_),
//     a0(a0_), a1(a1_), a2(a2_), a3(a3_), a4(a4_), f(NULL)
//   {
//     f = new T[(b0-a0)*(b1-a1)*(b2-a2)*(b3-a3)*(b4-a4)];
//   }
//   ~function5D()
//   {
//     delete[] f;
//   }
//   T operator()(int i0, int i1, int i2, int i3, int i4) const
//   {
//     return f[i4-a4 + b4*(i3-a3 + b3*(i2-a2 + b2*(i1-a1 + b1*(i0-a0))))];
//   }
//   T& operator()(int i0, int i1, int i2, int i3, int i4)
//   {
//     return f[i4-a4 + b4*(i3-a3 + b3*(i2-a2 + b2*(i1-a1 + b1*(i0-a0))))];
//   }
//   int ind(int i0, int i1, int i2, int i3, int i4) const
//   {return i4-a4 + b4*(i3-a3 + b3*(i2-a2 + b2*(i1-a1 + b1*(i0-a0))));}
// };

#endif
