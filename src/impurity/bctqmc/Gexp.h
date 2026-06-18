/*
#include <blitz/array.h>
#include <set>
#include <vector>
#include <unordered_map>
#include "interpolate.h"
using namespace blitz;
using namespace std;
#include "wfile.h"
#include "dynamic_mesh.h"
#include "read_cix.h"
#include "frequency_loops_real_axis.h"
#include "tinyMatmul.h"
*/

class nFexp{
  //  contains the fermi function and its algorithm and the argument.
  //  The use of x is such that log(f(-x)) = log(f(x)) + x
  //  We can generalize this to any product of fermi functions, for example:
  //    log(f(-x1)*f(-x2)/f(-x3)) = log(f(x1)*f(x)/f(x3)) + x1+x2-x3
  //    hence product of logs is the sum, and the arguments also sum.
public:
  double f, lf, x;
  nFexp(){};
  nFexp(double _f_, double _lf_, double _x_) : f(_f_), lf(_lf_), x(_x_){};
  nFexp(double _x_): x(_x_){
    f = fermi(x);
    lf = logfermi(x);
  }
  double mlf() const{ // returns log(fermi(-x)). We store log(fermi(x)) and x.
    //log(fermi(-x)) = log(fermi(x))+x
    return lf+x;
  }
  double plf() const{ // returns log(fermi(x))
    return lf;
  }
  bool empty(){ return (f==0 and lf==0 and x==0); }
  void set(double _f_, double _lf_, double _x_){ f=_f_; lf=_lf_; x=_x_;}
  nFexp& operator*=(const nFexp& F){ f*=F.f; lf+=F.lf; x+=F.x; return *this;}
  nFexp& operator/=(const nFexp& F){ f/=F.f; lf-=F.lf; x-=F.x; return *this;}
  friend nFexp operator*(const nFexp& F1, const nFexp& F2)
  { return nFexp(F1.f*F2.f, F1.lf+F2.lf, F1.x+F2.x); }
  friend nFexp operator/(const nFexp& F1, const nFexp& F2)
  { return nFexp(F1.f/F2.f, F1.lf-F2.lf, F1.x-F2.x); }
  void matmul(const nFexp& F){
    f *= F.f;
    lf += F.lf;
    x += F.x;
  }
  friend double  divide(const nFexp& F1, const nFexp& F2){
    if (F2.f>1e-2) return F1.f/F2.f;
    else return exp(F1.lf-F2.lf);
  }
  friend double divide_minus(const nFexp& F1, const nFexp& F2){
    // f1/f2(-x) = f1/(1-f2)
    if (fabs(1-F2.f)>1e-2){
      return F1.f/(1-F2.f);
    }else{
      return exp(F1.lf-(F2.lf+F2.x));
    }
  }
  void matmul_minus(const nFexp& F){
    // * f(-x)
    f *= (1.-F.f);  // fermi(-x)=1-fermi(x)
    lf += F.lf + F.x; // log(fermi(-x)) = log(fermi(x))+x
    x -= F.x;
  }
  double give_minus() const{ // computes 1-f(x) where f=f(x) and lf=log(f(x))
    // When f~1 f(-x)=1-f would be zero, which is not accurate.
    // We should rather use 1-e^l ~ -l-l^2/2-l^3/3!-.... = -l(1+l/2*(1+l/3))
    if ( fabs(1-f)>1e-2 )
      return  1-f;
    else
      return -lf*(1.+lf/2*(1.+lf/3*(1+lf/4)));
  }
  void combine(const nFexp& f1, const nFexp& f2, const nFexp& f3){
    // This corresponds to operation *this = f1 * f2 / f3;
    lf = f1.lf+f2.lf-f3.lf;
    x = f1.x+f2.x-f3.x;
    if (f3.f>1e-16)
      f = f1.f*f2.f/f3.f;
    else
      f = exp(lf);
  }
  void combinem(const nFexp& f1, const nFexp& f2, const nFexp& f3){
    // This corresponds to operation *this = f1(x) * f2(-x) / f3(-x);
    lf = f1.lf + (f2.lf+f2.x) - (f3.lf+f3.x);
    x = f1.x-f2.x+f3.x;
    if (1-f3.f>1e-10)
      f = f1.f * (1-f2.f) / (1-f3.f);
    else
      f = exp(lf);
  }
  bool check_Accuracy()
  {
    if (fabs(exp(lf)-f)>1e-5){
      cout<<"ERROR in nFexp: df="<<f-exp(lf)<<" f="<<f<<" lf="<<lf<<" x="<<x<<endl;
      return false;
    }
    return true;
  }
  friend ostream& operator<<(ostream& out, const nFexp& F);
};

const bool Compatible_interp=false;
// Forward declarations
class mGexp;

class Gexp{
  // This class contains obects that represent pseudo Green's function but with higher precision.
  // Namely, on real axis the imaginary part of the PP-Green's function falls of exponentially exp(beta*om).
  // But we need quantity like G*exp(-beta*om), which means we need to resolve exponentially small numbers.
  // We will represent G as complex number and in addition a list of prefactors and exponents, i.e.,
  // G''(x) = \sum_i -p[i]*exp(a[i])
public:
  complex<double> c; // complex representation
  double  p, a;      // p--real prefactor for the exponential term, a--exponent
  //
  Gexp() : c(0), p(0), a(0) {};
  Gexp(int _size_) : c(0), p(0), a(0) {};
  void Init2Unity(){ c=1; p=0; a=0; }
  void set2zero(){ c=0; p=0; a=0; }
  void resize(int _size_){
    c=0;
    p=0;
    a=0;
  }
  bool allclose(const Gexp& g){
    if (abs(c)<1e-16 and abs(g.c)<1e-16){
      if ( fabs(p*exp(a)-g.p*exp(g.a))>1e-10 )
	return false;
      else
	return true;
    }
    if (abs(g.c/c-1.0)>1e-6){
      //cout<<" |g.c/c-1|="<<abs(g.c/c-1.0)<<endl;
      return false;
    }
    if ( fabs(p*exp(a)-g.p*exp(g.a))>1e-6 ){
      //cout<<" |p*exp(a)-g.p*exp(g.a)|="<< fabs(p*exp(a)-g.p*exp(g.a)) << endl;
      return false;
    }else{
      return true;
    }
  }
  void set(const complex<double>& Gtc, double mfe, double lmfe){ // mfe=1-f, lmfe==log(1-f)
    constexpr auto small = 1e-16;
    constexpr auto log_small = -36.8413614879047;
    c = complex<double>(Gtc.real(), -Gtc.imag()*mfe); // This is definition of complex Green's function from our stored pair (G', widetilde{G}).
    p = 1.0;                     // log(-G'') = log(Gtc.imag())+log(mfe)=log(Gtc.imag())+lmfe
    if (fabs(Gtc.imag())!=0){
      a = log(fabs(Gtc.imag())) + lmfe; // log(Gtc.imag())+log(mfe)
    }else{
      // pretending that |Gtc.imag|=1e-16 because it was set to 0, which results in nan's
      a = log_small + lmfe;
    }
  }
  //void set(const complex<double>& Gtc, const Fexp& mfe){set(Gtc,mfe.f,mfe.lf);}
  void set(const complex<double>& Gtc, const nFexp& mfe){set(Gtc,mfe.f,mfe.lf);}
  void set(const complex<double>& Gtc, double _a_, const nFexp& mfe){
    //set(Gtc,mfe.f,mfe.lf);
    constexpr auto small = 1e-16;
    constexpr auto log_small = -36.8413614879047;
    c = complex<double>(Gtc.real(), -Gtc.imag()*mfe.f); // This is definition of complex Green's function from our stored pair (G', widetilde{G}).
    p = Gtc.imag()*exp(-_a_);                     // log(-G'') = log(Gtc.imag())+log(mfe.f)=log(Gtc.imag())+lmfe
    a = _a_ + mfe.lf;
  }
  Gexp& operator=(const Gexp& other){
    if (this != &other){
      c = other.c;
      p = other.p;
      a = other.a;
    }
    return *this;
  }
  Gexp& operator*=(const Gexp& Gb){
    if (p==0 and a==0){ // current G is either just constructed or is set to unity
      *this = Gb;      // calls your operator=(…)
    }else{
      if (a>Gb.a and a!=0){ // Decide which exponent is larger
	p = p * Gb.c.real() + c.real() * Gb.p * exp(Gb.a-a);
      }else{
	p = p * Gb.c.real() * exp(a-Gb.a) + c.real() * Gb.p;
	a = Gb.a;
      }
      //c *= Gb.c;
      c = complex<double>(c.real()*Gb.c.real()-c.imag()*Gb.c.imag(), -p*exp(a));
    }
    return *this;
  }
  Gexp& operator*=(double x){
    c *= x;
    p *= x;
    return *this;
  }
  Gexp& operator/=(const Gexp& G){
    if (a > G.a){
      p = (p*G.c.real() - c.real()*G.p*exp(G.a-a))/norm(G.c);
    }else{
      p = (-G.p*c.real() + G.c.real()*p*exp(a-G.a))/norm(G.c);
      a = G.a;
    }
    c /= G.c;
    return *this;
  }
  Gexp& operator*=(const mGexp& Gb);  
  Gexp& operator/=(const mGexp& Gb);
    
  void divide_by_positive(const Gexp& G){
    // We sometimes store G*F, where F is +1 or -1. We than want to divide by G, while we have G*F.
    // Since G is diagonal, we assume G.p>0 when F=1 and G.p<0 when F=-1. If we want to divide by G
    // and not G*F, we simply add -1 if G.p < 0.
    double norm_Gc = norm(G.c);
    double sgn = (G.p<0) ? -1 : 1;
    if (norm_Gc!=0){
      if (a > G.a){
	p = sgn*(p*G.c.real() - c.real()*G.p*exp(G.a-a))/norm_Gc;
      }else{
	p = sgn*(-G.p*c.real() + G.c.real()*p*exp(a-G.a))/norm_Gc;
	a = G.a;
      }
      c /= G.c*sgn;
    }else{ // it can happen that G==0 because ii=0 sometimes.
      if (norm(c)!=0){
	cout<<"WARNING divide_by_positive has x/y, where x="<<(*this)<<" y="<<G<<endl;
      }
    }
  }
  Gexp& matmul_lferm(double lf){
    // G *= exp(lf)
    c *= exp(lf);
    a += lf;
    return *this;
  }
  double getImt(double lmfe){ // gives Gt(x)=-G''(x)/f(-x)=p*exp(a-log(f(-x)))
    return p*exp(a-lmfe);
  }
  double get0(double lf, const nFexp& mfe){ //double lmfe, double mfe){
    // This gives 1/f(-x) * exp(lf) * Im(*this)
    // where lf contains logarithm of product of many fermi functions necessary for self-energy calculation, i.e., lf=log(f(-y1)*f(y2)*...)
    // and where mfe.lf = log(f(-x)) and mfe.f = f(-x)
    // therefore 1/f(-x) = exp(log(1/f(-x))) = exp(-log(f(-x)))= exp(-mfe.lf) 
    return p*exp(a+lf-mfe.lf);
  }
  double get1(double lf, const nFexp& mfe){ //double lmfe, double mfe){
    // This gives  exp(-beta*x) * exp(lf) * Im(*this)
    // where lf contains logarithm of product of many fermi functions necessary for self-energy calculation, i.e., lf=log(f(-y1)*f(y2)*...)
    // and where mfe.lf = log(f(-x)) and mfe.f = f(-x)
    // therefore exp(-beta*x) = f(x)/f(-x) = exp(-log(f(-x)))*(1-f(-x))=exp(-mfe.lf)*(1-mfe.f)
    //return p*exp(a-mfe.lf+lf)*(1-mfe.f);
    return p*exp(a+lf+mfe.x);
  }
  complex<double> getc(double lmfe){
    // This gives [real(G), Gt(x)=-G''(x)/f(-x)=p*exp(a-log(f(-x)))]
    return complex<double>(c.real(), p*exp(a-lmfe));
  }
  complex<double> getc1(double lf, double lmfe){
    // if *this==G, it gives G(x)*exp(lf)
    // It returns exp(lf) * ( Re(G), Im(G)/f(-x) ) where lmfe=log(f(-x))
    // hence we return ( Re(G)*exp(lf),  -Im(G) * exp(-lmfe+lf) )
    return complex<double>(c.real()*exp(lf), p*(exp(a-lmfe+lf)));
  }  
  void linInterp(const bl::Array<Gexp,2>& f, int il, const intpar& ip)
  {
    int i= ip.i; double px = ip.p, qx=1-ip.p;
    if (px > 1.01 or px<-0.01){ // it's extrapolation, which is very dangerous on extrapolating exponents. Newer use it, just set everythng to zero outside!
      p = 0.0;
      a = px>0 ? f(il,i+1).a : f(il,i).a; // keep the exponents at the boundary.
      c = 0.0;
      return;
    }
    if (not Compatible_interp)
      c = f(il,i).c*qx+f(il,i+1).c*px;
    p = f(il,i).p*qx+f(il,i+1).p*px;
    // Note that naively exponent should not be linearly interpolated. However, the numerically unstable part is of the calculation is where
    // G(x) = Gt(x)*f(-x) with Gt(x) of the order of unity and f(-x) exponentially small at very negative frequency. In this case
    // log(G(x)) = log(Gt(x))+log(f(-x))=log(Gt(x))+log(f(x))+beta*x. For negative x log(Gt(x))+log(f(x)) is small, but beta*x is very large and negative
    // hence only linear interpolation can correctly interpolate the large beta*x term.
    a = f(il,i).a*qx+f(il,i+1).a*px;
    if (Compatible_interp)
      c = complex<double>(f(il,i).c.real()*qx+f(il,i+1).c.real()*px, -p*exp(a));
  }
  void linInterp(const bl::Array<Gexp,2>& f, int il, const intpar& ip, double fact)
  {
    int i= ip.i; double px = ip.p, qx=1-ip.p;
    if (px > 1.01 or px<-0.01){ // it's extrapolation, which is very dangerous on extrapolating exponents. Newer use it, just set everythng to zero outside!
      p = 0.0;
      a = px>0 ? f(il,i+1).a : f(il,i).a; // keep the exponents at the boundary.
      c = 0.0;
      return;
    }
    if (not Compatible_interp)
      c = (f(il,i).c*qx+f(il,i+1).c*px)*fact;
    p = (f(il,i).p*qx+f(il,i+1).p*px)*fact;
    // Note that naively exponent should not be linearly interpolated. However, the numerically unstable part is of the calculation is where
    // G(x) = Gt(x)*f(-x) with Gt(x) of the order of unity and f(-x) exponentially small at very negative frequency. In this case
    // log(G(x)) = log(Gt(x))+log(f(-x))=log(Gt(x))+log(f(x))+beta*x. For negative x log(Gt(x))+log(f(x)) is small, but beta*x is very large and negative
    // hence only linear interpolation can correctly interpolate the large beta*x term.
    a = f(il,i).a*qx+f(il,i+1).a*px;
    if (Compatible_interp)
      c = complex<double>((f(il,i).c.real()*qx+f(il,i+1).c.real()*px)*fact, -p*exp(a) );
  }
  void linInterp(const bl::Array<Gexp,1>& f, const intpar& ip)
  {
    int i= ip.i; double px = ip.p, qx=1-ip.p;
    if (px > 1.01 or px<-0.01){ // it's extrapolation, which is very dangerous on extrapolating exponents. Newer use it, just set everythng to zero outside!
      p = 0.0;
      a = px>0 ? f(i+1).a : f(i).a; // keep the exponents at the boundary.
      c = 0.0;
      return;
    }
    if (not Compatible_interp)
      c = f(i).c*qx+f(i+1).c*px;
    p = f(i).p*qx+f(i+1).p*px;
    // Note that naively exponent should not be linearly interpolated. However, the numerically unstable part is of the calculation is where
    // G(x) = Gt(x)*f(-x) with Gt(x) of the order of unity and f(-x) exponentially small at very negative frequency. In this case
    // log(G(x)) = log(Gt(x))+log(f(-x))=log(Gt(x))+log(f(x))+beta*x. For negative x log(Gt(x))+log(f(x)) is small, but beta*x is very large and negative
    // hence only linear interpolation can correctly interpolate the large beta*x term.
    a = f(i).a*qx+f(i+1).a*px;
    if (Compatible_interp)
      c = complex<double>(f(i).c.real()*qx+f(i+1).c.real()*px, -p*exp(a));
  }
  pair<double,double> getI(){ // gives -G''(x)
    return make_pair(c.imag(), -p*exp(a));
  }
  friend ostream& operator<<(ostream& out, const Gexp& G);
};

ostream& operator<<(ostream& out, const Gexp& G)
{
  out<<setprecision(8)<< "["<<PComplex(G.c)<<"; a="<<setprecision(10)<<G.a<<";p="<<G.p<<"]";
  return out;
}
ostream& operator<<(ostream& out, const nFexp& F)
{
  out <<setprecision(6)<<"["<<F.f<<","<<F.lf<<","<<F.x<<"] ";
  return out;
}
ostream& operator<<(ostream& out, const pair<double,double>& x)
{
  out<<"["<<x.first<<","<<x.second<<"] ";
  return out;
}

class mGexp{
  // This class contains obects that represent pseudo Green's function but with higher precision.
  // Namely, on real axis the imaginary part of the PP-Green's function falls of exponentially exp(beta*om).
  // But we need quantity like G*exp(-beta*om), which means we need to resolve exponentially small numbers.
  // We will represent G as complex matrix (c) and in addition, we will also store imaginary part as a prefactor (p) and exponents (a), i.e.,
  // G''_{i,j}(x) = -p_{i,j}*exp(a)
public:
  int N1, N2;
  Array< std::complex<double>,2>  storage_c;
  Array< double,2>  storage_p;
  // --- active “view” into that storage
  Array<std::complex<double>,2>  c;
  Array<double,2>  p;
  double a;
  //
  mGexp() : storage_c(1,1), storage_p(1,1) {
    //cout<<"Default constructor"<<endl;
    setActiveShape(1,1);
    c=0; p=0; a=0;
  };
  mGexp(int _N1_, int _N2_) : storage_c(_N1_,_N2_), storage_p(_N1_,_N2_) {
    // initially point the views to zero-sized subarray
    setActiveShape(_N1_,_N2_);
    c=0; p=0; a=0;
  }
  void reference(mGexp& G){
    c.reference(G.storage_c(Range(0,G.N1-1),Range(0,G.N2-1)));
    p.reference(G.storage_p(Range(0,G.N1-1),Range(0,G.N2-1)));
    a = G.a;
    N1 = G.N1;
    N2 = G.N2;
  }
  // call this whenever you change (N1,N2,a):
  void setActiveShape(int Nnew1, int Nnew2) {
    N1 = Nnew1;  N2 = Nnew2;
    // re-bind the view to the top-left N1×N2 block
    c.reference( storage_c(Range(0,N1-1),Range(0,N2-1)) );
    p.reference( storage_p(Range(0,N1-1),Range(0,N2-1)) );
    //cout<<"setActiveShape to "<<N1<<","<<N2<<endl;
    //cout<<"c="<<c<<endl;
  }
  void resize(int _N1_, int _N2_, bool zero=false){
    if (storage_c.extent(0)<_N1_ or storage_c.extent(1)<_N2_){
      storage_c.resize(_N1_,_N2_);
      storage_p.resize(_N1_,_N2_);
      //cout<<"actively resizing storage to "<<storage_c.shape()<<endl;
    }
    setActiveShape(_N1_,_N2_);
    p(0,0)=0;
    a=0;
    if (zero){
      c=0;
      p=0;
    }
  }
  void Init2Unity(){
    p=0;
    a=0;
    //c=0;
    N1=N2=0;
    //a=0;
    //c=0;
    //for (int i=0; i<min(N1,N2); ++i) c(i,i)=1.0;
  }
  mGexp& operator=(const mGexp& other){
    if (this != &other){
      resize(other.c.extent(0),other.c.extent(1));
      c = other.c;
      p = other.p;
      a = other.a;
    }
    return *this;
  }
  void Set(int j1, int j2, const mGexp& other, double mm1)
  {
    p(j1,j2) = other.p(j1,j2) * mm1;
    c(j1,j2) = other.c(j1,j2) * mm1;
    a = other.a;
  }
  void Set(int j1, int j2, const complex<double>& z){
    c(j1,j2) = z;
    p(j1,j2) = z.imag();
    a = 0;
  }
  mGexp& operator*=(const mGexp& Gb){
    if ((N1==0 and N2==0) or (p(0,0)==0.0 and a==0)){ // current G is either just constructed or is set to unity
      //cout<<"copying Gb with sizes "<<Gb.p.extent(0)<<","<<Gb.p.extent(1)<<" p(0,0)="<<p(0,0)<<" N1="<<N1<<" N2="<<N2<<endl;
      resize(Gb.c.extent(0),Gb.c.extent(1));
      a = Gb.a;
      p = Gb.p;
      c = Gb.c;
      return *this;
    }
    if (storage_c.extent(1)<Gb.c.extent(1)){
      cout<<"ERROR c.extent(1)="<<c.extent(1)<<" is too small, should be "<<Gb.c.extent(1)<<endl;
      return *this;
    }
    TinyVector<int,2> n1(max(Gb.c.extent(0),c.extent(0)), max(Gb.c.extent(1),c.extent(1)));
    TinyVector<int,2> n2(max(c.extent(0), Gb.c.extent(0)), Gb.c.extent(1));
    TinyVector<int,2> n3(c.extent(0), Gb.c.extent(1));
    static thread_local Array<double,2> tmp1, tmp2, tmp3;
    
    if (tmp1.extent(0)<n1[0] or tmp1.extent(1)<n1[1]) tmp1.resize(n1[0],n1[1]);
    if (tmp2.extent(0)<n2[0] or tmp2.extent(1)<n2[1]) tmp2.resize(n2[0],n2[1]);
    if (tmp3.extent(0)<n3[0] or tmp3.extent(1)<n3[1]) tmp3.resize(n3[0],n3[1]);

    //cout<<"tmp1.shape="<<tmp1.shape()<<" tmp2.shape="<<tmp2.shape()<<" tmp3.shape="<<tmp3.shape()<<endl;
    
    static thread_local Array<double,2> tmpd1,tmpd2,tmpd3; // will be used as reference to slices of tmp1,tmp2,tmp3
    tmpd1.reference(tmp1(Range(0,Gb.c.extent(0)-1),Range(0,Gb.c.extent(1)-1))); tmpd1 = real(Gb.c); // tmpd1 = real(Gb.c)
    tmpd2.reference(tmp2(Range(0,c.extent(0)-1),Range(0,Gb.c.extent(1)-1)));

    //cout<<"tmpd1.shape="<<tmpd1.shape()<<" tmpd2.shape="<<tmpd2.shape()<<endl;
    tiny::dot(tmpd2, p, tmpd1);     // tmp2 = p*real(Gb.c)
    //cout<<"p*real(Gb.c)="<<tmpd2<<endl;
    
    tmpd1.reference(tmp1(Range(0,c.extent(0)-1),Range(0,c.extent(1)-1))); tmpd1 = real(c); // tmp1 = real(c)
    //cout<<"tmpd1.shape="<<tmpd1.shape()<<endl;
    tmpd3.reference(tmp3(Range(0,c.extent(0)-1),Range(0,Gb.c.extent(1)-1)));
    //cout<<"tmpd3.shape="<<tmpd3.shape()<<endl;
    tiny::dot(tmpd3, tmpd1, Gb.p); // tmp3 = real(c)*Gp.p
    //cout<<"real(c)*Gp.p="<<tmpd3<<endl;
    
    p.reference( storage_p(Range(0,p.extent(0)-1),Range(0,Gb.c.extent(1)-1)) );
    if (a>Gb.a and a!=0){ // Decide which exponent is larger
      p = tmpd2 + tmpd3*exp(Gb.a-a); // p*real(Gb.c) + real(c)*Gb.p*exp(Gb.a-a)
    }else{
      p = tmpd2*exp(a-Gb.a) + tmpd3; // p*real(Gb.c)*exp(a-Gb.a) + real(c)*Gb.p;
      a = Gb.a;
    }
    tmpd2.reference(tmp2(Range(0,Gb.c.extent(0)-1),Range(0,Gb.c.extent(1)-1))); tmpd2 = real(Gb.c);
    //cout<<"tmpd2.shape="<<tmpd2.shape()<<endl;
    tiny::dot(tmpd3, tmpd1, tmpd2); // tmpd3 = real(c) * real(Gb.c);
    tmpd1 = imag(c);
    tmpd2 = imag(Gb.c);
    c.reference( storage_c(Range(0,c.extent(0)-1),Range(0,Gb.c.extent(1)-1)) );
    c = tmpd3;                // c = real(c) * real(Gb.c)
    tiny::dot(tmpd3, tmpd1, tmpd2); // tmpd3 = imag(c) * imag(Gb.c);
    c -= tmpd3;               // c -= imag(c) * imag(Gb.c)
    c -= p*(exp(a)*complex<double>(0.0,1.0)); // imaginary part

    N1 = c.extent(0); N2 = c.extent(1);
    return *this;
  }
  void matmul(const mGexp& A, const mGexp& B){
    if (storage_c.extent(0)<A.c.extent(0) or storage_c.extent(1)<B.c.extent(1))
      resize(A.c.extent(0),B.c.extent(1));
    else
      setActiveShape(A.c.extent(0), B.c.extent(1));
    
    const double large_number = 200.;
    if (A.c.extent(0)==1 and A.c.extent(1)==1 and B.c.extent(1)==1){
      double t2 = A.p(0,0)        * B.c(0,0).real();
      double t3 = A.c(0,0).real() * B.p(0,0);
      if (A.a==0){ // special case where we can not take A.a>B.a
	if (-B.a < large_number){ // can evaluate nevertheless
	  p(0,0) = t2*exp(A.a-B.a) + t3; // p = A.p*real(B.c)*exp(A.a-B.a) + real(A.c)*B.p;
	}else{
	  // B.a is very negative hence exp_A_B is exploding
	  if (fabs(t2)>1e-100)  // can not evaluate exponent, and tmpd2 is nonzero and can't ignore the term
	    cout<<"ERROR in matmul exp("<<A.a-B.a<<") is too large and tmpd2="<<t2<<" is not small enough"<<endl;
	  p(0,0) = t3;
	}
	a = B.a;
      } else if (B.a==0){ // special case where we can not take A.a>B.a
	if (-A.a < large_number){ // can evaluate nevertheless
	  p(0,0) = t2 + t3 * exp(B.a-A.a);
	}else{
	  if (fabs(t3)>1e-100)  // can not evaluate exponent, and tmpd3 is nonzero and can't ignore the term
	    cout<<"ERROR in matmul exp("<<B.a-A.a<<") is too large and tmpd3="<<t3<<" is not small enough"<<endl;
	  p(0,0) = t2;
	}
	a = A.a;
      } else if (A.a>B.a){ // Decide which exponent is larger
	p(0,0) = t2 + t3 * exp(B.a-A.a);
	a = A.a;
      }else{
	p(0,0) = t2*exp(A.a-B.a) + t3;
	a = B.a;
      }
      c(0,0) = complex<double>( A.c(0,0).real()*B.c(0,0).real()-A.c(0,0).imag()*B.c(0,0).imag(), -p(0,0)*exp(a) );
      return;
    }
    TinyVector<int,2> n1(max(A.c.extent(0),B.c.extent(0)), max(A.c.extent(1),B.c.extent(1)));
    TinyVector<int,2> n2(max(A.c.extent(0),B.c.extent(0)), B.c.extent(1));
    TinyVector<int,2> n3(A.c.extent(0), B.c.extent(1));
    static thread_local Array<double,2> tmp1, tmp2, tmp3;
    
    if (tmp1.extent(0)<n1[0] or tmp1.extent(1)<n1[1]) tmp1.resize(n1[0],n1[1]);
    if (tmp2.extent(0)<n2[0] or tmp2.extent(1)<n2[1]) tmp2.resize(n2[0],n2[1]);
    if (tmp3.extent(0)<n3[0] or tmp3.extent(1)<n3[1]) tmp3.resize(n3[0],n3[1]);

    //cout<<"tmp1.shape="<<tmp1.shape()<<" tmp2.shape="<<tmp2.shape()<<" tmp3.shape="<<tmp3.shape()<<endl;
    static thread_local Array<double,2> tmpd1,tmpd2,tmpd3; // will be used as reference to slices of tmp1,tmp2,tmp3
    tmpd1.reference(tmp1(Range(0,B.c.extent(0)-1),Range(0,B.c.extent(1)-1))); tmpd1 = real(B.c); // tmpd1 = real(B.c)
    tmpd2.reference(tmp2(Range(0,A.c.extent(0)-1),Range(0,B.c.extent(1)-1)));

    //cout<<" matmul A="<<A<<" B="<<B<<endl;
    
    tiny::dot(tmpd2, A.p, tmpd1); // tmp2 = A.p*real(B.c)
    //cout<<" matmul: A.p="<<A.p<<" tmpd1="<<tmpd1<<" result="<<tmpd2<<endl;
    
    tmpd1.reference(tmp1(Range(0,A.c.extent(0)-1),Range(0,A.c.extent(1)-1))); tmpd1 = real(A.c); // tmpd1 = real(A.c)
    tmpd3.reference(tmp3(Range(0,A.c.extent(0)-1),Range(0,B.c.extent(1)-1)));
    tiny::dot(tmpd3, tmpd1, B.p); // tmp3 = real(A.c)*B.p
    //cout<<" matmul: tmpd1="<<tmpd1<<" B.p="<<B.p<<" result="<<tmpd3<<endl;

    if (A.a==0){ // special case where we can not take A.a>B.a
      if (-B.a < large_number) // can evaluate nevertheless
	p = tmpd2*exp(A.a-B.a) + tmpd3; // p = A.p*real(B.c)*exp(A.a-B.a) + real(A.c)*B.p;
      else{
	if (sum(abs(tmpd2))>1e-100)  // can not evaluate exponent, and tmpd2 is nonzero and can't ignore the term
	  cout<<"ERROR in matmul exp("<<A.a-B.a<<") is too large and tmpd2="<<tmpd2<<" is not small enough"<<endl;
	p = tmpd3;
      }
      a = B.a;
    } else if (B.a==0){ // special case where we can not take A.a>B.a
      if (-A.a < large_number) // can evaluate nevertheless
	p = tmpd2 + tmpd3*exp(B.a-A.a); // p = A.p*real(B.c) + real(A.c)*B.p*exp(B.a-A.a)
      else{
	if (sum(abs(tmpd3))>1e-100)  // can not evaluate exponent, and tmpd3 is nonzero and can't ignore the term
	  cout<<"ERROR in matmul exp("<<B.a-A.a<<") is too large and tmpd3="<<tmpd3<<" is not small enough"<<endl;
	p = tmpd2;
      }
      a = A.a;
    } else if (A.a>B.a){ // Decide which exponent is larger
      //cout<<"  matmul: tmpd2="<<tmpd2<<" tmpd3="<<tmpd3<<" B.a-A.a="<<B.a-A.a<<" exp="<<exp(B.a-A.a)<<endl;
      p = tmpd2 + tmpd3*exp(B.a-A.a); // p = A.p*real(B.c) + real(A.c)*B.p*exp(B.a-A.a)
      a = A.a;
    }else{
      //cout<<"  matmul: tmpd2="<<tmpd2<<" tmpd3="<<tmpd3<<" A.a-B.a="<<A.a-B.a<<" exp="<<exp(A.a-B.a)<<endl;
      p = tmpd2*exp(A.a-B.a) + tmpd3; // p = A.p*real(B.c)*exp(A.a-B.a) + real(A.c)*B.p;
      a = B.a;
    }
    tmpd2.reference(tmp2(Range(0,B.c.extent(0)-1),Range(0,B.c.extent(1)-1))); tmpd2 = real(B.c);
    //cout<<"tmpd2.shape="<<tmpd2.shape()<<endl;
    tiny::dot(tmpd3, tmpd1, tmpd2); // tmpd3 = real(A.c) * real(B.c);
    c = tmpd3;                      // c = real(A.c) * real(B.c)
    tmpd1 = imag(A.c);
    tmpd2 = imag(B.c);
    //c.reference( storage_c(Range(0,c.extent(0)-1),Range(0,Gb.c.extent(1)-1)) );
    tiny::dot(tmpd3, tmpd1, tmpd2); // tmpd3 = imag(A.c) * imag(B.c);
    c -= tmpd3;                     // c = real(A.c) * real(B.c) - imag(A.c) * imag(B.c)
    c -= p*(exp(a)*complex<double>(0.0,1.0)); // imaginary part is -p*exp(a)
  }
  void dot_diag(const mGexp& A, const Array<double,1>& Fi){
    assert(Fi.extent(0)==A.c.extent(1) && "ERROR in dot_diag: dimensions of c and Fi dont match.");
    if (storage_c.extent(0)<A.c.extent(0) or storage_c.extent(1)<A.c.extent(1))
      resize(A.c.extent(0),A.c.extent(1));
    else
      setActiveShape(A.c.extent(0), A.c.extent(1));
    a = A.a;
    //cout<<"c.shape="<<c.shape()<<" Fi.shape="<<Fi.shape()<<endl;
    for (int i = 0; i < c.extent(0); ++i) 
      for (int j = 0; j < c.extent(1); ++j){
	p(i,j) = A.p(i,j) * Fi(j);
	c(i,j) = A.c(i,j) * Fi(j);
      }
  }
  void dot_diag(const Array<double,1>& Fi, const mGexp& B){
    assert(Fi.extent(0)==B.c.extent(0) && "ERROR in dot_diag: dimensions of Fi and c dont match.");
    if (storage_c.extent(0)<B.c.extent(0) or storage_c.extent(1)<B.c.extent(1))
      resize(B.c.extent(0),B.c.extent(1));
    else
      setActiveShape(B.c.extent(0), B.c.extent(1));
    a = B.a;
    //cout<<"c.shape="<<c.shape()<<" Fi.shape="<<Fi.shape()<<endl;
    for (int i = 0; i < c.extent(0); ++i){
      p(i,Range::all()) = Fi(i) * B.p(i,Range::all());
      c(i,Range::all()) = Fi(i) * B.c(i,Range::all());
    }
  }
  mGexp& operator*=(const nFexp& F){
    c *= F.f;
    a += F.lf;
    return *this;
  }
  mGexp& matmul_lferm(double lf){
    c *= exp(lf);
    a += lf;
    return *this;
  }
  mGexp& operator*=(double d)
  {
    p *= d;
    c *= d;
    return *this;
  }
  void set(const Array<complex<double>,2>& Gtc, double _a_, double mfe, double lmfe){ // mfe=1-f, lmfe==log(1-f)
    cout<<"Gtc.shape="<<Gtc.shape()<<" while storage_c.shape="<<storage_c.shape()<<endl;
    if (storage_c.extent(0)<Gtc.extent(0) or storage_c.extent(1)<Gtc.extent(1)){
      cout<<"calling resize"<<endl;
      resize(Gtc.extent(0),Gtc.extent(1));
    }else{
      cout<<"skip calling resize, just increase active shape"<<endl;
      setActiveShape(Gtc.extent(0),Gtc.extent(1));
    }
    a = _a_ + lmfe;
    double exp_ma = exp(-_a_);
    for (int i=0; i<Gtc.extent(0); ++i){
      for (int j=0; j<Gtc.extent(1); ++j){
	c(i,j) = complex<double>(Gtc(i,j).real(), -Gtc(i,j).imag()*mfe);
	p(i,j) = Gtc(i,j).imag()*exp_ma;
      }
    }
  }
  //void set(const Array<complex<double>,2>& Gtc, double _a_, const Fexp& fe){set(Gtc,_a_,fe.f,fe.lf);}
  void set(const Array<complex<double>,2>& Gtc, double _a_, const nFexp& fe){set(Gtc,_a_,fe.f,fe.lf);}
  pair<Array<double,2>,Array<double,2>> getI(){ // gives -G''(x)
    Array<double,2> Gc(N1,N2), Gi(N1,N2);
    for (int i=0; i<N1; ++i)
      for (int j=0; j<N2; ++j){
	Gi(i,j) = -p(i,j)*exp(a);
	Gc(i,j) = c(i,j).imag();
      }
    return make_pair(Gc,Gi);
  }
  void get(Array<double,2>& res, double lmfe){ // gives Gt(x)
    res.resize(p.extent(0),p.extent(1));
    res = p*exp(a-lmfe);
  }
  void getc1(Array<complex<double>,2>& res, double lf, double lmfe){
    // if *this==G, it gives G(x)*exp(lf)
    // It returns exp(lf) * ( Re(G), Im(G)/f(-x) ) where lmfe=log(f(-x))
    // hence we return ( Re(G)*exp(lf),  Im(G) * exp(-lmfe+lf) )
    res.resize(p.extent(0),p.extent(1));
    res = real(c)*exp(lf) + p*(exp(a-lmfe+lf)*complex<double>(0.0,1.0));
  }
  void getc(Array<complex<double>,2>& res, double lmfe){
    // if *this==G, it gives Gt(x)
    // It returns Re(G) and Im(G)/f(-x), where lmfe=log(f(-x))
    // hence (Re(G),  Im(G) * exp(-lmfe) )
    res.resize(p.extent(0),p.extent(1));
    res = real(c) + p*(exp(a-lmfe)*complex<double>(0.0,1.0));
  }
  /*
  void linInterp(const bl::Array<mGexp,2>& f, int il, const intpar& ip)
  {
    int i= ip.i; double px = ip.p, qx=1-ip.p;
    const mGexp& f1 = f(il,i);
    const mGexp& f2 = f(il,i+1);
    setActiveShape(f1.c.extent(0),f1.c.extent(1));
    c = f1.c*qx+f2.c*px;
    p = f1.p*qx+f2.p*px;
    // Note that naively exponent should not be linearly interpolated. However, the numerically unstable part is of the calculation is where
    // G(x) = Gt(x)*f(-x) with Gt(x) of the order of unity and f(-x) exponentially small at very negative frequency. In this case
    // log(G(x)) = log(Gt(x))+log(f(-x))=log(Gt(x))+log(f(x))+beta*x. For negative x log(Gt(x))+log(f(x)) is small, but beta*x is very large and negative
    // hence only linear interpolation can correctly interpolate the large beta*x term.
    a = f1.a*qx+f2.a*px;
  }
  void linInterp(const bl::Array<mGexp,1>& f, const intpar& ip)
  {
    int i= ip.i; double px = ip.p, qx=1-ip.p;
    const mGexp& f1 = f(i);
    const mGexp& f2 = f(i+1);
    setActiveShape(f1.c.extent(0),f1.c.extent(1));
    c = f1.c*qx+f2.c*px;
    p = f1.p*qx+f2.p*px;
    // Note that naively exponent should not be linearly interpolated. However, the numerically unstable part is of the calculation is where
    // G(x) = Gt(x)*f(-x) with Gt(x) of the order of unity and f(-x) exponentially small at very negative frequency. In this case
    // log(G(x)) = log(Gt(x))+log(f(-x))=log(Gt(x))+log(f(x))+beta*x. For negative x log(Gt(x))+log(f(x)) is small, but beta*x is very large and negative
    // hence only linear interpolation can correctly interpolate the large beta*x term.
    a = f1.a*qx+f2.a*px;
  }
  */
  void linInterpolate(int j1, int j2, const bl::Array<Gexp,2>& Gms, int ii, const intpar& ip, double pref)
  {
    // Gx(j1,j2) = ii_sign0*linInterpolate(Gms(ii,ip))
    int i= ip.i; double px = ip.p, qx=1-ip.p;
    if (px > 1.01 or px<-0.01){ // it's extrapolation, which is very dangerous on extrapolating exponents. Newer use it, just set everythng to zero outside!
      p(j1,j2) = 0.0;
      a = px>0 ? Gms(ii,i+1).a : Gms(ii,i).a; // keep the exponents at the boundary.
      c(j1,j2) = 0.0;
      return;
    }
    if (not Compatible_interp)
      c(j1,j2) = pref * (Gms(ii,i).c*qx + Gms(ii,i+1).c*px);
    p(j1,j2) = pref * (Gms(ii,i).p*qx + Gms(ii,i+1).p*px);
    double _a_ = Gms(ii,i).a*qx + Gms(ii,i+1).a*px;
    if (a==0){
      a = _a_;
      //cout<<" Setting a="<<a<<" for ("<<j1<<","<<j2<<")="<<ii<<endl;
    }else{
      if (fabs(_a_/a-1.)>1e-7) cout<<"ERROR: previously set a="<<a<<"!="<<_a_<<" currently computed a for ("<<j1<<","<<j2<<")="<<ii<<endl;
      a = _a_;
    }
    if (Compatible_interp)
      c(j1,j2) = complex<double>( pref*(Gms(ii,i).c.real()*qx + Gms(ii,i+1).c.real()*px), -p(j1,j2)*exp(a) );
  }
  double Set_Gt_ferm_1D(const nFexp& mfe) const{
    // GGi = Im(G(x)) * exp(-beta*x) = -Gt(x) * f(-x) * exp(-beta*x) = -Gt(x) * f(x)
    // Here mfe = f(-x) in log form, hence mfe.x == -beta*x
    //GGi.resize(p.extent(0),p.extent(1));
    return p(0,0) * exp(a+mfe.x);
  }
  void Set_Gt_ferm(Array<double,2>& GGi, const nFexp& mfe) const{
    // GGi = Im(G(x)) * exp(-beta*x) = -Gt(x) * f(-x) * exp(-beta*x) = -Gt(x) * f(x)
    // Here mfe = f(-x) in log form, hence mfe.x == -beta*x
    //GGi.resize(p.extent(0),p.extent(1));
    GGi = p;
    double _exp_ = exp(a+mfe.x);
    GGi *= _exp_;
  }
  void Set_Gt_ferm(Array<complex<double>,2>& GGi, const nFexp& mfe) const{
    // GGi = Im(G(x)) * exp(-beta*x) = -Gt(x) * f(-x) * exp(-beta*x) = -Gt(x) * f(x)
    // Here mfe = f(-x) in log form, hence mfe.x == -beta*x
    //GGi.resize(p.extent(0),p.extent(1));
    GGi = p;
    double _exp_ = exp(a+mfe.x);
    GGi *= _exp_;
  }
  friend std::ostream& operator<<(std::ostream& os, const mGexp& g);
};
ostream& operator<<(ostream& out, const mGexp& G)
{
  out<<setprecision(8)<<"[c="<<G.c<<"; a="<<setprecision(10)<<G.a<<"; p="<<G.p<<"]";
  return out;
}

Gexp& Gexp::operator*=(const mGexp& Gb){
  if (Gb.c.extent(0)!=1 or Gb.c.extent(1)!=1){
    cout<<"ERROR Gexp::operator*=(mGexp _) can only use one dimensional mGexp. mGexp shape="<<Gb.c.shape()<<endl;
  }
  if (p==0 and a==0){ // current G is either just constructed or is set to unity
    c = Gb.c(0,0);
    p = Gb.p(0,0);
    a = Gb.a;
  }else{
    if (a>Gb.a and a!=0){ // Decide which exponent is larger
      p = p * Gb.c(0,0).real() + c.real() * Gb.p(0,0) * exp(Gb.a-a);
    }else{
      p = p * Gb.c(0,0).real() * exp(a-Gb.a) + c.real() * Gb.p(0,0);
      a = Gb.a;
    }
    //c *= Gb.c;
    c = complex<double>(c.real()*Gb.c(0,0).real()-c.imag()*Gb.c(0,0).imag(), -p*exp(a));
  }
  return *this;
};

Gexp& Gexp::operator/=(const mGexp& G){
  if (G.c.extent(0)!=1 or G.c.extent(1)!=1){
    cout<<"ERROR Gexp::operator/=(mGexp _) can only use one dimensional mGexp. mGexp shape="<<G.c.shape()<<endl;
  }
  if (a > G.a){
    p = (p*G.c(0,0).real() - c.real()*G.p(0,0)*exp(G.a-a))/norm(G.c(0,0));
  }else{
    p = (-G.p(0,0)*c.real() + G.c(0,0).real()*p*exp(a-G.a))/norm(G.c(0,0));
    a = G.a;
  }
  c /= G.c(0,0);
  return *this;
}

complex<double> trace_dot(const mGexp& A, const mGexp& B){
  // This routine computes Tr(A * B) = Tr( (A'-i*p_A*exp(a_A))*(B'-i*p_B*exp(a_B)) )
  // Re(Tr(A*B)) = Tr( A'*B')-Tr(p_A*p_B)*exp(a_A+a_B)
  // Im(Tr(A*B)) = -Tr(p_A * B') * exp(a_A) - Tr(A' * p_B) * exp(a_B)
  // Expect same shapes for A.c, A.p and B.c, B.p
  //const int M = A.c.extent(0);
  //const int N = A.c.extent(1);
  double reAB   = 0.0; // Tr(Re(A.c) * Re(B.c))-Tr(Im(A.c) * Im(B.c))
  double AreBp  = 0.0; // Tr(Re(A.c) * B.p)
  double ApBre  = 0.0; // Tr(A.p     * Re(B.c))
  // By your earlier trace_dot definition:
  // Tr(X*Y) = sum_{i,j} X(i,j) * Y(j,i)
  for (int i = 0; i < A.c.extent(0); ++i) {
    for (int j = 0; j < A.c.extent(1); ++j) {
      const complex<double> a = A.c(i,j);
      const complex<double> b = B.c(j,i); // note swapped indices for trace
      reAB += a.real() * b.real() - a.imag() * b.imag();
      AreBp  += a.real() * B.p(j,i); // Re(A.c) with B.p transposed
      ApBre  += A.p(i,j) * b.real(); // A.p with Re(B.c) transposed
    }
  }
  // Imaginary contribution factor:
  // Im(Tr(A*B)) = - [ Tr(A.p * ReB)*e^{a_A} + Tr(ReA * B.p)*e^{a_B} ]
  const double Im_p = ApBre * std::exp(A.a) + AreBp * std::exp(B.a);
  return complex<double>(reAB, -Im_p);
}
Gexp trace_dot_Gexp(const mGexp& A, const mGexp& B){
  // This routine computes
  //            Tr(A * B) which is in the form of Gexp, in which Im(Gexp)= -p*exp(a)
  // Tr( (A'-i*p_A*exp(a_A))*(B'-i*p_B*exp(a_B)) )
  // Re(Tr(A*B)) = Tr( A'*B')-Tr(p_A*p_B)*exp(a_A+a_B)
  // Im(Tr(A*B)) = -Tr(p_A * B') * exp(a_A) - Tr(A' * p_B) * exp(a_B)
  //
  // Expect same shapes for A.c, A.p and B.c, B.p
  //const int M = A.c.extent(0);
  //const int N = A.c.extent(1);
  double reAB   = 0.0; // Tr(Re(A.c) * Re(B.c))-Tr(Im(A.c) * Im(B.c))
  double AreBp  = 0.0; // Tr(Re(A.c) * B.p)
  double ApBre  = 0.0; // Tr(A.p     * Re(B.c))
  // By your earlier trace_dot definition:
  // Tr(X*Y) = sum_{i,j} X(i,j) * Y(j,i)
  for (int i = 0; i < A.c.extent(0); ++i) {
    for (int j = 0; j < A.c.extent(1); ++j) {
      const complex<double> a = A.c(i,j);
      const complex<double> b = B.c(j,i); // note swapped indices for trace
      reAB += a.real() * b.real() - a.imag() * b.imag();
      AreBp  += a.real() * B.p(j,i); // Re(A.c) with B.p transposed
      ApBre  += A.p(i,j) * b.real(); // A.p with Re(B.c) transposed
    }
  }
  // Imaginary contribution factor:
  // Im_p = -Im(Tr(A*B)) = [ Tr(A.p * ReB)*e^{a_A} + Tr(ReA * B.p)*e^{a_B} ]
  Gexp r;
  if (A.a > B.a){
    r.a = A.a;
    r.p = ApBre + AreBp * std::exp(B.a-A.a);
    r.c = complex<double>(reAB, -exp(r.a)*r.p);
  }else{
    r.a = B.a;
    r.p = ApBre * std::exp(A.a-B.a) + AreBp;
    r.c = complex<double>(reAB, -exp(r.a)*r.p);
  }
  return r;
}
/*
double trace_dot_get1(const mGexp& A, const mGexp& B, double lf, double lmfe, double mfe){
  // This routine computes
  //                exp(lf-beta*x)*Im(Tr(A * B)) = exp(lf-lmfe)*(1-mfe)*Im(Tr(A * B))
  // Note that exp(-lmfe)*(1-mfe)=(1-mfe)/mfe=1/mfe-1=1+exp(-beta*x)-1=exp(-beta*x)
  // Tr( (A'-i*p_A*exp(a_A))*(B'-i*p_B*exp(a_B)) )
  // Re(Tr(A*B)) = Tr( A'*B')-Tr(p_A*p_B)*exp(a_A+a_B)
  // Im(Tr(A*B)) = -Tr(p_A * B') * exp(a_A) - Tr(A' * p_B) * exp(a_B)
  //
  // Expect same shapes for A.c, A.p and B.c, B.p
  //const int M = A.c.extent(0);
  //const int N = A.c.extent(1);
  //double reAB   = 0.0; // Tr(Re(A.c) * Re(B.c))-Tr(Im(A.c) * Im(B.c))
  double AreBp  = 0.0; // Tr(Re(A.c) * B.p)
  double ApBre  = 0.0; // Tr(A.p     * Re(B.c))
  // By your earlier trace_dot definition:
  // Tr(X*Y) = sum_{i,j} X(i,j) * Y(j,i)
  for (int i = 0; i < A.c.extent(0); ++i) {
    for (int j = 0; j < A.c.extent(1); ++j) {
      const complex<double> a = A.c(i,j);
      const complex<double> b = B.c(j,i); // note swapped indices for trace
      //reAB += a.real() * b.real() - a.imag() * b.imag();
      AreBp  += a.real() * B.p(j,i); // Re(A.c) with B.p transposed
      ApBre  += A.p(i,j) * b.real(); // A.p with Re(B.c) transposed
    }
  }
  // Imaginary contribution factor:
  // Im_p = -Im(Tr(A*B))*exp(-lmfe) = [ Tr(A.p * ReB)*e^{a_A-lmfe} + Tr(ReA * B.p)*e^{a_B-lmfe} ]
  const double Im_p = ApBre * std::exp(A.a+lf-lmfe) + AreBp * std::exp(B.a+lf-lmfe);
  return -Im_p*(1-mfe);
}
*/


/*
ostream& operator<<(ostream& out, const set<int>& s){
  out<<"{";
  for (int x : s) {
    out << x << ",";
  }
  out << "}"<<endl;
  return out;
}
Array<double,2> DeterminPhiR(const Array<complex<double>,2>& Gm, const Array<double,1>& xn, const Array<double,1>& few, const Array<double,1>& dh, const Array<double,1>& logo, double beta)
{
  // This quantity is needed to evaluate Phi on the real axis. We want Phi = 1/Norder Tr(\Sigma_{pseudo} G_{pseudo}) and in time we just need to evaluate time integral.
  // In frequency is harder. We first show that
  // Tr(Sigma_p * G_p) = -e^{-beta*lambda} \int dx/pi e^{-beta*x} Im(\Sigma_p(x) G_p(x)) = -e^{-beta*lambda} \int dx/pi e^{-beta*x} (\Sigma'_p(x) G''_p(x))+\Sigma''_p(x) G'_p(x))
  // Taking into account the Kramars-Kronig property \Sigma'_p(x) = -1/pi \int \Sigma''_p(y)/(x-y)dy (notice no static part to Sigma')
  // we can show that
  // Tr(Sigma_p * G_p) = -e^{-beta*lambda} \int dx/pi [ e^{-beta*x} \Sigma''_p(x) G'_p(x) - \Sigma''_p(x) \widetilde{G}'_p(x) ]
  //
  // where \widetilde{G}'_p(x) = -1/pi \int e^{-beta y} G''_p(y)/(x-y)dy=-1/pi \int f(y) \widetilde{G}_p(y)/(x-y)dy
  // hence we can compute 
  // Tr(Sigma_p * G_p) = -e^{-beta*lambda} \int dx/pi [f(x)G'_p(x)-f(-x)\widetilde{G}'_p(x)]*widetilde{\Sigma}_p(x)
  // hence we need PhiR_p(x) == -[f(x)G'_p(x)-f(-x)\widetilde{G}'_p(x)]/pi and Phi= 1/Q*\int dx PhiR_p(x)*widetilde{\Sigma}_p(x)
  //
  // Alternative derivation:
  //   Phi = -1/<Q> * \int( PhiR_p(w) * \widetilde{\Sigma}_p(w), {w,-\inf,\inf})
  //   This comes from the fact that
  //   Phi = -1/<Q> * Int[ e^{-\beta x} * Im(G_p(x)\Sigma_p(x)) ,{x,-\inf,\inf}]/pi
  //   which we transform to
  //   PhiR(w) = -( Re(G_p(w))*fe(w) + fe(-w)*Int[ widetilde{G}_p(x)*f(x)/(w-x)*(-1/pi),{x,-inf,inf}] )/pi
  //   PhiR(w) = -( Re(G_p(w))*fe(w) + fe(-w)*KK_w(f(y)\widetilde{G}_p(y)) )/pi
  //
  // Note that this element-wise propery works in the absence spin-orbit and other broken symmetries that make G non-symmetric.
  // In general we might need to work with matrices also here.
  // Note that when Greens function is not symmetric (G^T != G) this does not work anymore. We need to work with complex and hermitian quantites, namely,
  // G'(x) =  1/2   *( G(x) + G^\dagger (x) )
  // G''(x) = 1/(2i)*( G(x) - G^\dagger (x) )
  // and both G'(x) and G''(x) are Hermitian, and G = G' + i*G''. However G'!=Re(G) and G''!=Im(G).
  // Kramars Kronig is working for G'(x) and G''(x) but not for Im(G) and Re(G), and hence it does not work component-wise.
  // The issue is that we also need to define \widetilde{G} = G''(x)/f(-x) because we need G''(x)*exp(-beta*x), but not  \widetilde{G} = Im(G)(x)/f(-x)
  // The issue is that all these matrices are now complex, and is more expensive to calculate and store. But it should be generalized when we find
  // that G^T != G.
  int  Npseudo_unique=Gm.extent(0);
  Array<double,2> PhiR(Npseudo_unique,xn.size());
  Array<double,1> mGm(xn.size()), mGr(xn.size());
  Array<double,1> iGm(xn.size());
  Array<double,1> pGr(xn.size());
  for (int ip=0; ip<Npseudo_unique; ++ip){
    iGm(Range::all()) = imag(Gm(ip,Range::all()));
    pGr(Range::all()) = real(Gm(ip,Range::all()));
    mGm = iGm*few;
    mGr = KramarsKronig(mGm, xn, dh, logo);
    PhiR(ip,Range::all()) = -(pGr(Range::all())*few(Range::all()) + mGr(Range::all())*(1.0-few(Range::all())))/M_PI;
  }
  return PhiR;
}


int main()
{
  complex<double> Im(0,1.);
  
  auto[xn, dGm] = loadtxt("PPG.NCA.1");
  Array<complex<double>,2> Gm(dGm.extent(0)/2,xn.size());
  for (int i=0; i<dGm.extent(0)/2; ++i){
    Gm(i,Range::all()) = dGm(2*i,Range::all()) + dGm(2*i+1,Range::all())*Im;
  }

  double beta=250;
  
  Cix cix;
  cix.Read_cix("bctqmc.cix", 8.0, 4.0, beta);
  Diag2Matrix cdiag2matrix(5, cix);
  Determine_Symmetry(cix, cdiag2matrix, false);
  cix.PrintSymmetry(cout);
  
  
  Array<double,1> few(xn.size());
  Spline1DLin<double> femx(xn.size()), lfemx(xn.size());
  //femx.resize(xn.size()); lfemx.resize(xn.size());  // f(-xn), lof(f(-xn))
  for (int i=0; i<xn.size(); ++i){
    femx[i] = fermi(-xn(i)*beta);
    lfemx[i] = logfermi(-xn(i)*beta);
    few(i) = 1-femx[i];
  }
  {
    Array<double,1> fem_mid(xn.size()-1), lfem_mid(xn.size()-1); // evaluation the same functions at the midpoint
    for (int i=0; i<xn.size()-1; ++i){
      double x = 0.5*(xn(i)+xn(i+1));
      fem_mid(i) = fermi(-x*beta);
      lfem_mid(i) = logfermi(-x*beta);
    }
    femx.splineIt(xn,fem_mid); // now using the mid-point knowledge, we can see at which point is linear interpolation becomes more accurate
    lfemx.splineIt(xn,lfem_mid,0,-1);
  }

  Array<double,1> logo(xn.size());
  logo = CalcLogoOnMesh(xn);
  Array<double,1> dhxn(xn.size());
  dhxn = GiveDh(xn);
  Array<double,2> PhiR = DeterminPhiR(Gm,xn,few,dhxn,logo,beta);

  Array<nFexp,1> lmfx(xn.size());    // combination of f(-x) and log(f(-x))
  //lmfx.resize(xn.size()); // Array<nFexp,1> lmfx(xn.size());
  for (int i=0; i<xn.size(); ++i) lmfx(i) = nFexp( femx[i], lfemx[i], -xn(i)*beta );
  
  Array<double,2> common_a(cix.Npseudo_unique,xn.size());
  for (int ix=0; ix<xn.size(); ++ix){
    // Computes average Im(G) in every block, and across equivalent blocks, so that
    // we take average over all blocks that share a common component
    for (int i=0; i<cix.unique_blocks.size(); ++i){
      double average_Gim=0;
      int how_many = 0;
      for (int j=0; j<cix.unique_blocks[i].size(); ++j){
	int ib = cix.unique_blocks[i][j];
	int Ni = cix.Block[ib].size();
	for (int a=0; a<Ni; ++a)
	  if (cix.Sgm_symmetry[ib](a,a)!=0){
	    int ii = abs(cix.Sgm_symmetry[ib](a,a))-1;
	    average_Gim += fabs(Gm(ii,ix).imag());
	    how_many += 1;
	  }
      }
      average_Gim /= how_many;

      // pretending that |Gtc.imag|>=1e-16 because otherwise we get nan's when taking log
      double av = (average_Gim>1e-16) ? log(average_Gim) : log(1e-16);
      // now we can compute a of components in the block and across the blocks.
      // They must have the same exponent a if they share one common matrix element.
      for (int j=0; j<cix.unique_blocks[i].size(); ++j){
	int ib = cix.unique_blocks[i][j];
	int Ni = cix.Block[ib].size();
	for (int a=0; a<Ni; ++a){
	  for (int b=0; b<Ni; ++b){
	    int ii = abs(cix.Sgm_symmetry[ib](a,b))-1;
	    if (ii>=0) common_a(ii,ix) = av;
	  }
	}
      }
      //cout<<"average_Gim["<<ix<<","<<unique_blocks[i]<<"]="<<average_Gim<<" common_a="<< av <<endl;
    }
  }
  Array<Gexp,2> Gms;
  Gms.resize(cix.Npseudo_unique,xn.size()); //Array<nGexp,2> Gms(Npseudo_unique,xn.size());
  for (int ip=0; ip<cix.Npseudo_unique; ++ip){
    for (int i=0; i<xn.size(); ++i){
      Gms(ip,i).set(Gm(ip,i), common_a(ip,i), lmfx(i));
    }
  }

  Array<double,1> dlom(xn);
  Spline1DLin<double> fep, lfep;
  fep.resize(dlom.size()); lfep.resize(dlom.size());   // f(dlom), log(f(dlom))
  for (int i=0; i<dlom.size(); ++i){
    fep[i] = fermi(dlom(i)*beta);
    lfep[i] = logfermi(dlom(i)*beta);
  }
  {
    Array<double,1> fep_mid(dlom.size()-1), lfep_mid(dlom.size()-1); // here we also evaluate mid-points
    for (int i=0; i<dlom.size()-1; ++i){
      double x = 0.5*(dlom(i)+dlom(i+1));
      fep_mid(i) = fermi(x*beta);
      lfep_mid(i) = logfermi(x*beta);
    }
    fep.splineIt(dlom,fep_mid); // so that we can use cubic spline only in the regions where it is more accurate
    lfep.splineIt(dlom,lfep_mid,0,-1);
  }
  
  
  Array<int,2> diag(3,2); diag=0,2,1,4,3,5;
  Array<int,1> dirspin(5); dirspin=11,9,11;
  auto& allowem = cdiag2matrix.Cmp_matrix(diag, dirspin, cix);
  int n = diag.extent(0);

  Array<int,2> wtb;
  DetermineLoops(wtb,
		 diag, dirspin);
  Array<double,1> freq(n); freq=-1.0,-2.0,-3.0;
  Array<int,1> xfsgn(2*n);
  DetermineFermiSigns_fsgn_only(xfsgn,
                                 wtb);
  Array<nFexp,1> fry(n);
  for (int iy=0; iy<n; ++iy){ // 0:Norder for all possible configurations
    int dire = dirspin(iy) % 2;    // direction of hybridization propagator
    int iorb = dirspin(iy) / 2;    // its spin 
    int s = cix.bath_sigind(iorb); // actual entry for hybridization
    double y = freq(1+iy);         // frequency on hybridization
    intpar jpy = Interp(y, dlom);  // 
    fry(iy).set(fep(jpy), lfep(jpy), beta*y ); // fermi function for f(y) stored as logarithm log(f(y))
  }
  Array<double,1> lfpt(2*n);     // -- contains log(f(y1)*f(-y2)..*f(y_n)), the product needed for computing Sigma_pp of pseudoparticle pp.
  for (int j=0; j<2*n; ++j){ // going over the ring
    // Note: lfpt(ipseudo) = log(  \product_{iy} f(\pm_{ipseudo,iy} y)  )
    double lfp=0; // This is the product of fermi functions needed for Self-energy Sigma_j of pseudoparticle j.
    for (int iy=0; iy<n; ++iy){
      if (xfsgn(j) & (1<<iy)){ // xfsgn is represented as binary array in an integer. 
	lfp += fry(iy).plf(); // log(f(y))
      }else{
	lfp += fry(iy).mlf(); // log(f(-y))
      }
    }
    lfpt(j)=lfp; // This is the product of fermi functions needed for Self-energy Sigma_j of pseudoparticle j.
    //fmt::print("fpt[{:d},{:d}]={:g},{:f}\n", ic,j, fpt(ic,j).f, fpt(ic,j).lf);
    //if (DEBUG) cout<<" pp_fr["<<j<<"]="<<pp_freq_t(j)<<" px="<<px<<" Gmt="<<Gmt_t(Range(0,PPs-1),j)<<" mGppt="<<mGppt_t(Range(0,PPs-1))<<endl;
  }
  
  {
    vector<mGexp> Gx;
    vector<Array<double,1>> Fi;
    Array<mGexp,2> Lft, Rgh;

    Array<double,2> dS(cix.Npseudo_unique,2*n); dS=0;
  
    Lft.resize(allowem.size(),2*n); // Array<Array<double,2>,2> Lft(:,:)(N1,N2)
    Rgh.resize(allowem.size(),2*n); // Array<Array<double,2>,2> Rgh(:,:)(N1,N1)
    Gx.resize(2*n); // vector<Array<double,2>>  Gx[2*Norder](N1,N2)
    Fi.resize(2*n); // vector<Array<double,1>>  Fi[2*Norder](N1)

    Array<nFexp,1> frms_(2*n); 
    Array<intpar,1> pxs(2*n);
    for (int j=0; j<2*n; ++j){ // going over all pseudo-particles
      Range all_freq(0,n);
      double x = sum(freq(all_freq)*wtb(n+j,all_freq));
      pxs(j) = Interp(x, xn);     // interpolating point on xn mesh
      frms_(j).set( femx(pxs(j)), lfemx(pxs(j)), -beta*x); // fermi function f(-x) with argument x of pseudo propagator j.
    }
  
  
    int idx=0;
    for (const auto& [ib, alw] : cdiag2matrix.allowed){
      double phiL=0, phiR=0, phi_diag=0;
      int N0 = alw[0].extent(0);
      Lft(idx,0).resize(N0,N0,true);
    
      for (int j=0; j<N0; ++j)  // Lft(:,0) is just set to F
	Lft(idx,0).Set(j,j, sign(alw[0](j,1)));
    
      for (int i=0; i<2*n; ++i){
	int ip1 = (i+1)%(2*n);
	int N1 = alw[i].extent(0), N2 = alw[ip1].extent(0);
	mGexp GxF(N1,N2);
	Gx[i].resize(N1,N2,true);
	Fi[i].resize(N1);
	Array<TinyVector<int,2>,2> iis(N1,N2);
	int ib_a = cix.iBlck(abs(alw[i](0,1)))-1;
	for (int j1=0; j1<N1; ++j1){
	  int _a_ = alw[i](j1,1);
	  int a = cix.inBlock(abs(_a_)); // alw[0] index is given in terms of atomic state index, but a0 is index within the block. cix.inBlock translates this.
	  double mm0 = sign(_a_);         // matrix element of operator op(0) <alw[0](j0,1)|op(0)|alw[0](j0,0)> is stored in the sign of op(0) acting on index, hence sign of alw[0](j0,1) contains op(0) matrix element.
	  Fi[i](j1)=mm0;
	  for (int j2=0; j2<N2; ++j2){
	    int _b_ = alw[ip1](j2,0);
	    int b = cix.inBlock(abs(_b_));
	    //int ibb = cix.iBlck(_b_)-1;
	    int _ii0_ = cix.Sgm_symmetry[ib_a](a,b); // Green's function index for propagator G(a0,b0)==G[|_ii0_|-1]
	    if (_ii0_==0) continue;
	    int ii_sign0 = sign(_ii0_);                 // Green's function value for off-diagonal propagators could be negative/opposite, hence sign
	    int ii = abs(_ii0_)-1;
	    iis(j1,j2) = TinyVector<int,2>(_a_,_b_);
	    //cout<<"i="<<i<<" ib="<<ib_a<<" a="<<_a_<<" b="<<_b_<<" ii="<<ii<<endl;
	    Gx[i].linInterpolate(j1,j2, Gms, ii, pxs(i), ii_sign0); // Gx(j1,j2) = ii_sign0*linInterpolate(Gms(ii,ip)) : this is compatible with uG_(ipp,i)
	    if (i<2*n-1){
	      int mm1 = sign(alw[ip1](j2,1));    // matrix element: <alw[1](j1,1)|op(1)|alw[1](j1,0)>
	      GxF.Set(j1,j2, Gx[i], mm1);        // G*F, i.e., GxF(j1,j2) = Gx(j1,j2)*mm1
	    }
	  }
	}
	//cout<<"Gx ="<<Gx[i]<<endl;
	if (i<2*n-1){
	  cout<<"GxF="<<GxF<<endl;
	  Lft(idx,i+1).matmul( Lft(idx,i), GxF);   //Lft(idx,i+1) = matmul(Lft(idx,i),GxF) : this is similar to uGpp(ipp)
	  cout<<"Lft["<<idx<<","<<i+1<<"]="<<Lft(idx,i+1)<<endl;
	}else{
	  //phiL = trace_dot(Lft(idx,i),Gx[i]);
	}
      }
    
      int N1 = alw[2*n-1].extent(0);
      Rgh(idx,1).resize(N1,N0);
      Rgh(idx,1).dot_diag(Fi[2*n-1], Gx[2*n-1]); // Rgh(idx,1) = diag(Fi[2*n-1]) @ Gx[2*n-1]
      cout<<"Rgh["<<idx<<","<<1<<"]="<<Rgh(idx,1)<<endl;
      for (int i=1; i<2*n-1; ++i){
	int N1 = alw[2*n-1-i].extent(0), N2 = alw[2*n-i].extent(0);
	mGexp FGx(N1,N2);
	FGx.dot_diag(Fi[2*n-1-i], Gx[2*n-1-i]); // FGx = diag(Fi[2*n-1-i]) @ Gx[2*n-1-i]
	//cout<<"Gx ="<<Gx[i]<<endl;
	cout<<"FGx="<<FGx<<endl;
	Rgh(idx,i+1).matmul(FGx, Rgh(idx,i)); // Rgh(idx,i+1) = matmul(FGx,Rgh(idx,i));
	cout<<"Rgh["<<idx<<","<<i+1<<"]="<<Rgh(idx,i+1)<<endl;
      }
      double dphi=0;
      for (int ist=0; ist<2*n; ++ist){
	int ip1 = (ist+1)%(2*n);
	int N1 = alw[ist].extent(0), N2 = alw[ip1].extent(0);
	int ib_a = cix.iBlck(abs(alw[ist](0,1)))-1;
	{
	  mGexp Phis_over_Gm(N2,N1);
	  if (ist<2*n-1)
	    Phis_over_Gm.matmul(Rgh(idx,2*n-1-ist), Lft(idx,ist) ); // Right*Left gives the self-energy contribution, i.e., roughly Phi/G
	  else{
	    Phis_over_Gm = Lft(idx,ist);
	  }
	  cout << "Rgh["<<2*n-1-ist<<"]*Lft["<<ist<<"]="<<Phis_over_Gm<<endl;
	  cout << " *log(f)="<<lfpt(ist)<<" /log(f)="<<frms_(ist).lf<<endl;
	  Array<double,2> _Phis_over_Gm_(N2,N1);
	  Phis_over_Gm.matmul_lferm(lfpt(ist));
	  Phis_over_Gm.get(_Phis_over_Gm_, frms_(ist).lf);
	  cout<<" Phis_over_Gm="<< _Phis_over_Gm_ << endl;
	  //Phis_over_Gm *= (Phi_dd/W)*sgns; // also need sign and hybridization part, and MC-fake-weight.
	  for (int j1=0; j1<N1; ++j1){
	    int _a_ = alw[ist](j1,1);
	    int a = cix.inBlock(abs(_a_)); // alw[0] index is given in terms of atomic state index, but a0 is index within the block. cix.inBlock translates this.
	    for (int j2=0; j2<N2; ++j2){
	      int _b_ = alw[ip1](j2,0);
	      int b = cix.inBlock(abs(_b_));
	      int _ii0_ = cix.Sgm_symmetry[ib_a](b,a); // Green's function index for propagator G(a0,b0)==G[|_ii0_|-1]
	      if (_ii0_==0) continue;
	      int ii_sign0 = sign(_ii0_);              // Green's function value for off-diagonal propagators could be negative/opposite, hence sign
	      int ii = abs(_ii0_)-1;                   // correction to which PP self-energy
	      //if (debug) cout<<"   "<<ccstr<<"P("<<j2<<","<<j1<<")=S("<<_b_<<","<<_a_<<")=S("<<_ii0_<<")="<<Phis_over_Gm(j2,j1)<<endl;
	      double dSt = ii_sign0*_Phis_over_Gm_(j2,j1);
	      dS(ii,ist) += dSt; // ii is PP index, ist is the time-difference index.
	      //Stval(n_mino,ii,Range::all()) += utli(ist,Range::all())*ii_sign0*Phis_over_Gm(j2,j1);
	      double phiR = linInterp(PhiR,ii,pxs(ist));
	      dphi += phiR*dSt;
	      cout<<"dS["<<ii<<","<<ist<<"]="<<ii_sign0*_Phis_over_Gm_(j2,j1)<<endl;
	    }
	  }
	}
      }
      dphi /= (2*n);
      cout<<"Phi["<<idx<<"]= "<<dphi<<endl;
    
      idx += 1;
    }
  }
  return 0;
}
*/
