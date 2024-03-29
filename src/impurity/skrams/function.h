// @Copyright 2007 Kristjan Haule
// 
#ifndef FUNCTION_
#define FUNCTION_
#include <iostream>
#include <algorithm>
#include <cstring>
#include "util.h"

//////////////////////////// Hierarchy of classes: ///////////////////////////////////////
//                             
//                                     base clas  = functionb<>
//                   |                                                |
//               function1D<>                                     funProxy<>
//
//                                       function2D<funProxy>
//

//*****************************************//
// Classes also used in this header file   //
//*****************************************//
class intpar;
//class intpari;

//*******************************************************//
// Classes and functions implemented in this header file //
//*******************************************************//
template<class T> class functionb;
template<class T> class function1D;
template<class T> class funProxy;
template<class T> class function2D;


//********************************************************************************//
// Base class for two derived classes: function1D<T> and function2D<T>.	  	  //
// It is also used as a proxy class for function2D. Function2D<T> consists of	  //
// arrays of functionb<T> rather than functions1D<T>.				  //
// Memory is allocated in a fortran-like fashion for better performance.	  //
// Linear interpolation is implemented with the operator() that takes one	  //
// argument (class intpar).							  //
//********************************************************************************//
template<class T>
class functionb{
protected:
  T *f;
  int N0, N;
public:
  T& operator[](int i) {Assert(i<N,"Out of range in functionb[]"); return f[i];}
  const T& operator[](int i) const {Assert(i<N,"Out of range in functionb[]"); return f[i];}
  const T& last() const {return f[N-1];}
  int size() const { return N;}
  int fullsize() const { return N0;}
  T operator()(const intpar& ip) const;
  //  T operator()(const intpari& ip, const functionb& fp) const;
  functionb& operator+=(const functionb& m);
  functionb& operator*=(const T& m);
  T accumulate();
  T* MemPt(){return f;}
  const T* MemPt() const{return f;}
  functionb& operator=(const T& c);
  
  template <class meshx>
  void CalcFermOnMesh(double beta, const meshx& om);
  template <class meshx>
  void CalcLogOnMesh(const meshx& om);
  template <class meshx>
  void CalcTanhOnMesh(double beta, const meshx& om);
  template <class meshx>
  void CalcBoseOnMesh(double beta, const meshx& om);
  template <class meshx, class functiond>
  void KramarsKronig(const meshx& om, const functiond& logo);
  template <class meshx, class functiond, class functiond1D>
  void KramarsKronig(const functiond& Sigt, const meshx& om, const functiond1D& fe, const functiond1D& logo);
  template <class functor>
  void Set(functor& functn);
  void SetProduct(const functionb& m, const T& p);
  
  void SumUp(const functionb<T>& x, const functionb<T>& y);
protected:
  functionb() : f(NULL), N0(0), N(0) {};
  explicit functionb(int N_) : N0(N_), N(N_) {};
  ~functionb(){};
  functionb(const functionb&){};
  template<class U> friend class function2D;
  template <class U> friend U scalar_product(const functionb<U>& f1, const functionb<U>& f2);
};

//******************************************************************//
// One dimensional functions derived from functionb<T>. It has it's  //
// own constructors and destructors.				    //
//******************************************************************//
template <class T>
class function1D : public functionb<T>{
public:
  function1D(){};
  explicit function1D(int N_);
  ~function1D();
  function1D(const function1D& m);
  void resize(int N_);
  // equivalent to f=-f0
  function1D& assignm(const function1D& f0);
  function1D& operator=(const function1D& m);
  function1D& operator=(const T& c) {functionb<T>::operator=(c); return *this;}
  template <class meshx>
  void CalcFermOnMesh(double beta, const meshx& om);
  template <class meshx>
  void CalcLogOnMesh(const meshx& om);
  template <class meshx>
  void CalcBoseOnMesh(double beta, const meshx& om);
  template <class meshx>
  void CalcTanhOnMesh(double beta, const meshx& om);
  template <class meshx, class functiond>
  void KramarsKronig(const functiond& Sigt, const meshx& om, const functiond& fe, const functiond& logo);
  template <class meshx, class functiond>
  void KramarsKronig(const meshx& om, const functiond& logo);
  function1D<double> treshold(const function1D<double>& fe);
};

template <class T>
class funProxy : public functionb<T>{
public:
  void Initialize(int N_, T* f_);
  void ReInitialize(int N_, T* f_);
  void resize(int N_);
  funProxy& operator=(const functionb<T>& m);
  ~funProxy(){};
};

//**********************************************************************//
// Two dimentional function<T> derived from functionb<T>. It consists	//
// of an array of functionb<T> rather tham function1D<T>.		//
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
  funProxy<T>& operator[](int i) {Assert(i<N,"Out of range in function2D[]"); return f[i];}
  const funProxy<T>& operator[](int i) const {Assert(i<N,"Out of range in function2D[]"); return f[i];}
  const T& operator()(int i, int j) const {Assert(i<N && j<Nd,"Out of range in function2D(i,j)"); return f[i].f[j];}
  T& operator()(int i, int j) {Assert(i<N && j<Nd,"Out of range in function2D(i,j)"); return f[i].f[j];}
  T* MemPt() { return data;}
  const T* MemPt() const { return data;}
  const int size_N() const {return N;}
  const int size_Nd() const {return Nd;}
  const int fullsize_N() const {return N0;}
  const int fullsize_Nd() const {return Nd0;}
  const int lda() const {return Nd0;}
  void resize(int N_, int Nd_);
  
  function2D& operator=(const function2D& m);
  function2D& operator+=(double x);
  function2D& operator+=(const function2D& m);
  function2D& operator-=(double x);
  function2D& operator-=(const function2D& m);
  function2D& operator=(const T& u);
  function2D& operator*=(const T& x);
  
  template <class functor>
  inline void Set(functor& functn);
  template <class functor, class W>
  inline void Set(functor& functn, const function2D<W>& Um);

  void Product(const function2D& A, const function2D& B, const T& alpha=1, const T& beta=0);
  void DotProduct(const function2D& A, const function2D& B, const T& alpha=1, const T& beta=0);
  void TProduct(const function2D& A, const function2D& B, const T& alpha=1, const T& beta=0);
  void SymmProduct(const function2D& A, const function2D& B, const T& alpha=1, const T& beta=0);
  void TSymmProduct(const function2D& A, const function2D& B, const T& alpha=1, const T& beta=0);
  //  void TProduct(const function2D& A, const function2D& B, const T& alpha, const T& beta, int offset);
  
  template <class meshx, class functiond>
  void KramarsKronig(const meshx& om, const functiond& logo, functiond& sum, functiond& odvSig);
};

template <class fun, class fun1, class fun2>
void multiply(fun& f, const fun1& a, const fun2& b);

// functionb ////////////////////////////////////////////////////////////////
template<class T>
inline T functionb<T>::operator()(const intpar& ip) const
{
  return  f[ip.i]+ip.p*(f[ip.i+1]-f[ip.i]);
}

template<class T> 
inline functionb<T>& functionb<T>::operator+=(const functionb& m)
{
  CTMA_LOG(if (N!=m.size()) cerr << "Functions not of equal length! Can't sum!" << std::endl;)
  for (int i=0; i<N; i++) f[i] += m[i];
  return *this;
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
  CTMA_LOG(if (N<=0) cerr << "Size of functionb is non positive! "<<N<<std::endl;)
  for (int i=0; i<N; i++) f[i] = c;
  return *this;
}

template <class T>
template <class meshx>
inline void functionb<T>::CalcFermOnMesh(double beta, const meshx& om)
{
  for (int i=0; i<om.size(); i++){
    f[i] = 1.0/(1.0+exp(om[i]*beta));
  }
}

template <class T>
template <class meshx>
inline void functionb<T>::CalcBoseOnMesh(double beta, const meshx& om)
{
  for (int i=0; i<om.size(); i++){
    f[i] = 1.0/(exp(om[i]*beta)-1.0);
  }
}

template <class T>
template <class meshx>
inline void functionb<T>::CalcTanhOnMesh(double beta, const meshx& om)
{
  for (int i=0; i<om.size(); i++){
    double e = exp(beta*om[i]);
    f[i] = (e-1)/(e+1);
  }
}

template <class T>
template <class meshx>
inline void functionb<T>::CalcLogOnMesh(const meshx& om)
{
  for (int i=0; i<om.size(); i++){
    f[i] = (i!=0 && i!=om.size()-1) ? log((om.last()-om[i])/(om[i]-om[0])) : 0.0;
  }
}

template <class T>
template <class meshx, class functiond>
void functionb<T>::KramarsKronig(const meshx& om, const functiond& logo)
{
  if (logo.size()!=om.size()||om.size()!=N)
    cerr<<"Functions not of equal size in KramarsKronig"<<std::endl;
  
  for (int i=0; i<om.size(); i++){
    const double omom2 = (i>0) ? om.Delta(i-1) : 0.0;
    const int ip1 = (i<om.size()-1) ? i+1 : i, im1 = (i>0) ? i-1 : i;
    const double odvSig = 0.5*(om.Delta(i)*(f[ip1].imag()-f[i].imag()) + omom2*(f[i].imag()-f[im1].imag()));
    const double Sigii = f[i].imag();
    double sum = 0;
    for (int j=0; j<om.size(); j++){
      if (i!=j)
	sum += (f[j].imag()-Sigii)*om.Dh(j)/(om[j]-om[i]);
      else
	sum += odvSig*om.Dh(j);
    }
    f[i].real()=(sum+Sigii*logo[i])/M_PI;
  }
}

template <class T>
template <class meshx, class functiond, class functiond1D>
inline void functionb<T>::KramarsKronig(const functiond& Sigi, const meshx& om, const functiond1D& fe, const functiond1D& logo)
{
  if (fe.size()!=om.size()||Sigi.size()!=om.size())
    cerr<<"Functions not of equal size in KramarsKronig"<<std::endl;
  
  for (int i=0; i<om.size(); i++) f[i].imag()=(1-fe[i])*Sigi[i];
  
  KramarsKronig(om, logo);
}

template <class T>
template <class functor>
inline void functionb<T>::Set(functor& functn)
{
  for (int i=0; i<N; i++) f[i]=functn(f[i]);
}

template <class T>
inline void functionb<T>::SetProduct(const functionb& m, const T& p)
{
  CTMA_LOG(if (N<m.N) cerr<<"Size of function too small in SetProduct!"<<std::endl;);
  for (int i=0; i<m.N; i++) f[i] = m.f[i]*p;
}

template <class T>
void functionb<T>::SumUp(const functionb<T>& x, const functionb<T>& y)
{
  if (x.size()!=y.size()) std::cerr<<"Size of functions to sum up is different!"<<std::endl;
  if (N<x.size()) {
    std::cerr<<"The target functionb to small in the SumUp functionb! Can't continue!"<<std::endl;
    return;
  }
  for (int i=0; i<x.size(); i++){
    f[i] = x[i]+y[i];
  }
}

// function1D ////////////////////////////////////////////////////////////
template<class T>
inline function1D<T>::function1D(int N_) : functionb<T>(N_)
{
  this->f = new T[N_];
}

template<class T>
inline function1D<T>::~function1D()
{
  delete[] this->f;
  this->f = NULL;
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

template<class T>
inline function1D<T>& function1D<T>::assignm(const function1D& f0)
{
  resize(f0.N);
  for (int i=0; i<this->N; i++) this->f[i] = -f0.f[i];
  return *this;
}

template <class T>
template <class meshx>
inline void function1D<T>::CalcFermOnMesh(double beta, const meshx& om)
{
  resize(om.size());
  functionb<double>::CalcFermOnMesh(beta,om);
}

template <class T>
template <class meshx>
inline void function1D<T>::CalcBoseOnMesh(double beta, const meshx& om)
{
  resize(om.size());
  functionb<double>::CalcBoseOnMesh(beta,om);
}

template <class T>
template <class meshx>
inline void function1D<T>::CalcTanhOnMesh(double beta, const meshx& om)
{
  resize(om.size());
  functionb<double>::CalcTanhOnMesh(beta,om);
}

template <class T>
template <class meshx>
inline void function1D<T>::CalcLogOnMesh(const meshx& om)
{
  resize(om.size());
  functionb<double>::CalcLogOnMesh(om);
}

template <class T>
template <class meshx, class functiond>
inline void function1D<T>::KramarsKronig(const functiond& Sigi, const meshx& om, const functiond& fe, const functiond& logo)
{
  resize(om.size());
  functionb<T>::KramarsKronig(Sigi, om, fe, logo);
}
template <class T>
template <class meshx, class functiond>
inline void function1D<T>::KramarsKronig(const meshx& om, const functiond& logo)
{
  resize(om.size());
  functionb<T>::KramarsKronig(om, logo);
}

template <>
inline function1D<double> function1D<double>::treshold(const function1D<double>& fe)
{
  CTMA_LOG(if (fe.size()!=N) std::cerr<<"Fermi function not of correct size!"<<std::endl;);
  function1D<double> S(N);
  for (int i=0; i<N; i++) S[i]=f[i]*(1-fe[i]);
  return S;
}


// funProxy ///////////////////////////////////////////////////////////////
template <class T>
inline void funProxy<T>::Initialize(int N_, T* f_)
{
  this->N = this->N0 = N_; this->f = f_;
}

template <class T>
inline void funProxy<T>::ReInitialize(int N_, T* f_)
{
  this->N = N_; this->f = f_;
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
  
  for (int i=0; i<N0; i++) f[i].Initialize(Nd0,data+i*Nd0);
}

template<class T>
function2D<T>::~function2D()
{
  for (int i=0; i<N0; i++){
    f[i].~funProxy<T>();
  }
  operator delete(memory);
  memory = NULL;
}

template <class T>
inline function2D<T>& function2D<T>::operator=(const function2D& m)
{
  if (m.N<=N0 && m.Nd<=Nd0){
    N = m.N; Nd = m.Nd;
    for (int i=0; i<N; i++) memcpy(f[i].f, m.f[i].f, sizeof(T)*Nd);
  } else{
    int msize = sizeof(funProxy<T>)*m.N+sizeof(T)*m.Nd*m.N+HPoffset;
    operator delete(memory);
    memory = operator new (msize);
    Assert(memory!=NULL,"Out of memory");
    memcpy(memory, m.memory, msize);
    N = N0 = m.N; Nd = Nd0 = m.Nd;
    f = new (memory) funProxy<T>[N];
    int offset = sizeof(funProxy<T>)*N+HPoffset;
    data = reinterpret_cast<T*>(static_cast<char*>(memory)+offset);
    for (int i=0; i<N; i++) f[i].Initialize(Nd, data+i*Nd);
  }
  return *this;
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
inline void function2D<T>::Product(const function2D& A, const function2D& B, const T& alpha, const T& beta)
{
  if (A.Nd != B.Nd || !B.N || !A.N || !A.Nd || !B.Nd || N0<A.N || Nd0<B.N)
    std::cerr << " Matrix sizes not correct" << std::endl;
  xgemm("T", "N", B.N, A.N, B.Nd, alpha, B.MemPt(), B.Nd0, A.MemPt(), A.Nd0, beta, MemPt(), Nd0);
  N = A.N; Nd = B.N;  
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

template <class fun, class fun1, class fun2>
inline void multiply(fun& f, const fun1& a, const fun2& b)
{
  CTMA_LOG(if (a.size()!=b.size()) std::cerr << "Functions not of equal length! Can't multiply!" << std::endl;);
  f.resize(a.size());
  for (int i=0; i<f.size(); i++) f[i] = a[i]*b[i];
}

template <class T>
std::ostream& operator<<(std::ostream& stream, const functionb<T>& f)
{
  int width = stream.width(); 
  for (int i=0; i<f.size(); i++) stream<<i<<" "<<std::setw(width)<<f[i]<<std::endl;
  return stream;
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
template<class meshx, class functionx, class functiony, class functionz, class functiona, class functionb>
void print(std::ostream& stream, const meshx& om, const functionx& f1, const functiony& f2, const functionz& f3,
	   const functiona& f4, const functionb& f5, int width)
{
  if (om.size()!=f1.size() || om.size()!=f2.size() || om.size()!=f3.size() || om.size()!=f4.size() || om.size()!=f5.size())
    std::cerr<<"Can't print objectc of different size!"<<std::endl;
  for (int i=0; i<om.size(); i++)
    stream <<std::setw(width)<<om[i]<<std::setw(width)<<f1[i]<<std::setw(width)<<f2[i]<<std::setw(width)<<f3[i]<<std::setw(width)<<f4[i]<<std::setw(width)<<f5[i]<<std::endl;
}
template<class meshx, class functionx, class functiony, class functionz, class functiona, class functionb, class functionc>
void print(std::ostream& stream, const meshx& om, const functionx& f1, const functiony& f2, const functionz& f3,
	   const functiona& f4, const functionb& f5, const functionc& f6, int width)
{
  if (om.size()!=f1.size() || om.size()!=f2.size() || om.size()!=f3.size() || om.size()!=f4.size() || om.size()!=f5.size() || om.size()!=f6.size())
    std::cerr<<"Can't print objectc of different size!"<<std::endl;
  for (int i=0; i<om.size(); i++)
    stream <<std::setw(width)<<om[i]<<std::setw(width)<<f1[i]<<std::setw(width)<<f2[i]<<std::setw(width)<<f3[i]<<std::setw(width)<<f4[i]<<std::setw(width)<<f5[i]<<std::setw(width)<<f6[i]<<std::endl;
}
template<class meshx, class functionx, class functiony, class functionz, class functiona, class functionb, class functionc>
void print(std::ostream& stream, const meshx& om, const functionx& f1, const functiony& f2, const functionz& f3,
	   const functiona& f4, const functionb& f5, const functionc& f6, const functiona& f7, const functionb& f8, const functionc& f9,
	   int width)
{
  if (om.size()!=f1.size() || om.size()!=f2.size() || om.size()!=f3.size() || om.size()!=f4.size() || om.size()!=f5.size() || om.size()!=f6.size())
    std::cerr<<"Can't print objectc of different size!"<<std::endl;
  for (int i=0; i<om.size(); i++)
    stream <<std::setw(width)<<om[i]<<std::setw(width)<<f1[i]<<std::setw(width)<<f2[i]<<std::setw(width)<<f3[i]<<std::setw(width)<<f4[i]<<std::setw(width)<<f5[i]<<std::setw(width)<<f6[i]<<std::setw(width)<<f7[i]<<std::setw(width)<<f8[i]<<std::setw(width)<<f9[i]<<std::endl;
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
 
template <class funProxy, class T>
inline int SolveSOLA(function2D<funProxy>& A, function1D<T>& B, function1D<int>& ipiv)
{
  return xgesv(A.size_Nd(), 1, A.MemPt(), A.fullsize_Nd(), ipiv.MemPt(), B.MemPt(), B.fullsize());
}

inline double logpart(double a, double b, double x)
{
  return log(fabs((b-x)/(x-a)));
}
inline dcomplex logpart(double a, double b, const dcomplex& x)
{
  return log((b-x)/(a-x));
}
template <class T>
T KramarsKronig(const functionb<double>& fi, const mesh& om, const T& x, int i0, double S0)
{
  T sum=0;
  for (int j=0; j<i0-1; j++) sum += (fi[j]-S0)*om.Dh(j)/(om[j]-x);
  if (i0>0)                  sum += (fi[i0-1]-S0)*(om.Dh(i0-1)+0.5*om.Dh(i0))/(om[i0-1]-x);
  if (i0<om.size()-1)        sum += (fi[i0+1]-S0)*(om.Dh(i0+1)+0.5*om.Dh(i0))/(om[i0+1]-x);
  for (int j=i0+2; j<om.size(); j++) sum += (fi[j]-S0)*om.Dh(j)/(om[j]-x);
  if (x!=om.last() && x!=om[0]) sum += S0*logpart(om[0],om.last(),x);
  return sum/M_PI;
}
#endif
