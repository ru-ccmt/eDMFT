#ifndef _POLYNOM_
#define _POLYNOM_

#include <iostream>
#include <blitz/array.h>
#include "util.h"
namespace bl = blitz;
using namespace std;

bl::Array<int,1> binomialCoeff(int n){
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
  bl::Array<int,1> Coeff(n+1);
  for (int i=0; i<=n; ++i) Coeff(i) = C(n,i);
  return Coeff;
}

void PolyInt_m1(bl::Array<double,1>& In, int i_n, int N)
{
  //   Computes :
  //        Integrate[ t^i /(i_n + t ) , {t,0,1}]
  //    can be shown to satisfy the following recursion relation:
  //        In[n+1] = 1/(n+1) - i_n * In[n]
  In.resize(N); // need In for all nonzero orders
  In = 0;
  if (N==0) return;
  //  Upward recursion is very straighforard:
  //  In[n+1] = 1/(n+1) - i_n * In[n]   and  I[0] = log((i_n+1)/i_n)
  if (i_n==0){ 
    In(0) = 0;
    for (int n=1; n<N; ++n) In(n) = 1.0/n;
  }else{
    In(0) = log((i_n+1.0)/i_n);
    if (N==1) return;
    In(1) = 1 - i_n * In(0);
    for (int n=1; n<N-1; ++n)
      In(n+1) = 1/(n+1.) - i_n * In(n);
  }
}

/*
bl::Array<double,1> multiply(const bl::Array<double,1>& A, const bl::Array<double,1>& B)
{
  bl::Array<double,1> prod(A.extent(0)+B.extent(0)-1);
  prod = 0;
  // Multiply two polynomials term by term
  // Take ever term of first polynomial
  // Multiply the current term of first polynomial
  // with every term of second polynomial.
  for (int i=0; i<A.extent(0); ++i)
    for (int j=0; j<B.extent(0); ++j)
      prod(i+j) += A(i)*B(j);
  return prod;
}
 
// A utility function to print a polynomial
void printPoly(const bl::Array<double,1>& poly)
{
  for (int i=0; i<poly.extent(0); ++i){
    std::cout << poly(i)<<" ";
    if (i != 0) std::cout << "x^" << i;
    if (i != poly.extent(0)-1) std::cout << " + ";
  }
  std::cout<<std::endl;
}


double GetPolyInt(const bl::Array<double,1>& Poly, int i_n, int power=0){
  //    Computes integrals: 
  //        \sum_i a_i * Integrate[ t^i * (i_n + t )^power , {t,0,1}]
  //    For positive powers, the result is just:
  //        \sum_i a_i t^i*(i_n+t)^power =
  //        \sum_{k=0,power} Cobinatorial(power,k) i_n^k \sum_i a_i/(power-k+i+1)
  //    For power == -1 , the corresponding integral
  //        Integrate[ t^i /(i_n + t ) , {t,0,1}]
  //    can be shown to satisfy the following recursion relation:
  //        In[n+1] = 1/(n+1) - i_n * In[n]
  //
  int N = Poly.extent(0);
  bl::Array<double,1> In(N); // need In for all nonzero orders
  if (power>=0){
    bl::Array<int,1> Cn = binomialCoeff(power);
    long i_n_2_k = 1;
    double dsum = 0.0;
    for (int k=0; k<=power; ++k){
      double ds=0;
      for (int i=0; i<N; ++i) ds += Poly(i)/(power-k+i+1.);
      dsum += Cn(k)*i_n_2_k*ds;
      i_n_2_k *= i_n;
    }
    return dsum;
  }else{
    switch(power){
    case -1 :
      PolyInt_m1(In, i_n, N);
      break;
    case -2 :
      if (i_n==0){
	In=0;
	for (int n=2; n<N; ++n) In(n) = 1./(n-1.);
      }else{
	bl::Array<double,1> Kn(N-1); // need In for all nonzero orders
	PolyInt_m1(Kn, i_n, N-1);
	In(0) = 1.0/i_n - 1.0/(i_n+1.0);
	for (int n=1; n<N; ++n) In(n) = -1./(i_n+1.) + n * Kn(n-1);
      }
      break;
    case -3 :
      if (i_n==0){
	In=0;
	for (int n=3; n<N; ++n) In(n) = 1./(n-2.);
      }else{
	bl::Array<double,1> Kn(N-2); // need In for all nonzero orders
	PolyInt_m1(Kn, i_n, N-2);
	In(0) = 0.5*( 1./(i_n*i_n) - 1./((i_n+1.)*(i_n+1.)) );
	In(1) = -0.5/(i_n+1.)*( 1.0/(i_n+1.) + 1 ) + 0.5/i_n;
	for (int n=2; n<N; ++n) In(n) = -0.5/(i_n+1.)*( 1.0/(i_n+1.) + n ) + 0.5*(n-1)*n*Kn(n-2);
      }
      break;
    default : // any negative number
      int k=-power-1;
      if (i_n==0){
	In=0;
	for (int n=k+1; n<N; ++n) In(n) = 1.0/(n-k);
      }else{
	bl::Array<double,1> Kn(N-k); // need In for all nonzero orders
	PolyInt_m1(Kn, i_n, N-k);
	double c0 = 0;
	bl::Array<double,1> cp(k+1), cq(k+1);
	cp(0) = ipower( (i_n+1.), -k)/(k+0.0);
	cq(0) = ipower( (i_n+0.), -k)/(k+0.0);
	// cp(0) = 1/(i+1)^k/k                   cq(0) = 1/i^k/k
	// cp(1) = 1/(i+1)^{k-1}/(k*(k-1))       cq(1) = 1/i^{k-1}/(k*(k-1))
	// cp(2) = 1/(i+1)^{k-2}/(k*(k-1)*(k-2)) cq(2) = 1/i^{k-2}/(k*(k-1)*(k-2))
	// cp(k) = 1/(k!)                        cq(k) = 1/(k!)
	for (int i=1; i<k; ++i){
	  cp(i) = cp(i-1) * (i_n+1.)/(k-i);
	  cq(i) = cq(i-1) * (i_n+0.)/(k-i);
	}
	cp(k) = cp(k-1) * (i_n+1.);
	cq(k) = cq(k-1) * (i_n+0.);
	//In(0) = -[ cp(0) - cq(0)];
	//In(1) = -[ cp(0) + ( cp(1) - cq(1) ) ];
	//In(2) = -[ cp(0) + 2*cp(1) + 2*1*( cp(2) - cq(2) ) ];
	//In(3) = -[ cp(0) + 3*cp(1) + 3*2* cp(2) + 3*2*1*(cp(3)-cq(3)) ];
	//In(i) = -[ cp(0) + i*cp(1) + i*(i-1)*cp(i-1) + i*(i-1)*(i-2)*cp(i-2) ... ]
	for (int i=0; i<k; ++i){
	  double dsm=0;
	  double dd=1;
	  for (int j=0; j<i; ++j){ 
	    dsm += dd*cp(j); 
	    dd *= (i-j);     
	  }
	  dsm += dd*(cp(i)-cq(i));
	  In(i) = -dsm;
	}
	for (int n=k; n<N; ++n){
	  double ds = 0.0;
	  double dd0 = 1.0;
	  for (int i=0; i<k; ++i){
	    ds -= dd0*cp(i);
	    dd0 *= (n-i);
	  }
	  In(n) = ds + dd0*cp(k)*Kn(n-k);
	}
      }
    }
    double dsum = 0.0;
    for (int n=0; n<N; ++n) dsum += In(n)*Poly(n);
    return dsum;
  }
}
*/

class Polynomial{
public:
  bl::Array<double,1> Poly;
  //  
  explicit Polynomial(int n=0) : Poly(n)
  { Poly=0; }
  Polynomial(const bl::Array<double,1>& A) : Poly(A){}
  //    Computes integrals: 
  //        \sum_i a_i * Integrate[ t^i * (i_n + t )^power , {t,0,1}]
  double CmpIntegral(int i_n, int power=0);
  //
  double& operator[](int i){
    return Poly(i);
  }
  Polynomial& operator=(double x){ Poly=x; return *this;}
  Polynomial& operator*=(const Polynomial& A)
  {
    bl::Array<double,1> B = Poly.copy(); // Need to make a copy of this* first
    Poly.resize(A.Poly.extent(0)+B.extent(0)-1);
    Poly=0;
    for (int i=0; i<A.Poly.extent(0); ++i)
      for (int j=0; j<B.extent(0); ++j)
	Poly(i+j) += A.Poly(i)*B(j);
    return *this;
  }
  Polynomial& operator*=(double x){
    for (int i=0; i<Poly.extent(0); ++i) Poly(i) *= x;
    return *this;
  }
  Polynomial& operator+=(const Polynomial& A){
    if (Poly.extent(0)<A.Poly.extent(0)){
      bl::Array<double,1> B = Poly.copy();
      Poly.resize(A.Poly.extent(0));
      Poly=0;
      for (int i=0; i<B.extent(0); ++i) Poly(i) = A.Poly(i) + B(i);
      for (int i=B.extent(0); i<Poly.extent(0); ++i) Poly(i) = A.Poly(i);
    }else{
      for (int i=0; i<A.Poly.extent(0); ++i) Poly(i) += A.Poly(i);
    }
    return *this;
  }
  
  Polynomial operator* (const Polynomial& A)
  {
    Polynomial prod(A.Poly.extent(0)+Poly.extent(0)-1);
    prod.Poly = 0;
    for (int i=0; i<A.Poly.extent(0); ++i)
      for (int j=0; j<Poly.extent(0); ++j)
	prod.Poly(i+j) += A.Poly(i)*Poly(j);
    return prod;
  }
  friend std::ostream& operator<< (std::ostream& stream, const Polynomial& A);
  bool CheckIsLegitimate()
  {
    for (int i=0; i<Poly.extent(0); i++)
      if (std::isnan(Poly(i))){
	cout<<" ERROR Polynomial is not legitimate: "<< *this << endl;
	return false;
      }
    return true;
  }
};

double Polynomial::CmpIntegral(int i_n, int power)
{
  //    Computes integrals: 
  //        \sum_i a_i * Integrate[ t^i * (i_n + t )^power , {t,0,1}]
  //    For positive powers, the result is just:
  //        \sum_i a_i t^i*(i_n+t)^power =
  //        \sum_{k=0,power} Cobinatorial(power,k) i_n^k \sum_i a_i/(power-k+i+1)
  //    For power == -1 , the corresponding integral
  //        Integrate[ t^i /(i_n + t ) , {t,0,1}]
  //    can be shown to satisfy the following recursion relation:
  //        In[n+1] = 1/(n+1) - i_n * In[n]
  //
  int N = Poly.extent(0);
  //cout<<" Size of the polynomial at integration = "<< N <<" and power = "<<power<< endl;
  bl::Array<double,1> In(N); // need In for all nonzero orders
  if (power>=0){
    bl::Array<int,1> Cn = binomialCoeff(power);
    long i_n_2_k = 1;
    double dsum = 0.0;
    for (int k=0; k<=power; ++k){
      double ds=0;
      for (int i=0; i<N; ++i) ds += Poly(i)/(power-k+i+1.);
      dsum += Cn(k)*i_n_2_k*ds;
      i_n_2_k *= i_n;
    }
    return dsum;
  }else{
    switch(power){
    case -1 :
      PolyInt_m1(In, i_n, N);
      break;
    case -2 :
      if (i_n==0){
	In=0;
	for (int n=2; n<N; ++n) In(n) = 1./(n-1.);
      }else{
	bl::Array<double,1> Kn(N-1); // need In for all nonzero orders
	PolyInt_m1(Kn, i_n, N-1);
	In(0) = 1.0/i_n - 1.0/(i_n+1.0);
	for (int n=1; n<N; ++n) In(n) = -1./(i_n+1.) + n * Kn(n-1);
      }
      break;
    case -3 :
      if (i_n==0){
	In=0;
	for (int n=3; n<N; ++n) In(n) = 1./(n-2.);
      }else{
	bl::Array<double,1> Kn(N-2); // need In for all nonzero orders
	PolyInt_m1(Kn, i_n, N-2);
	In(0) = 0.5*( 1./(i_n*i_n) - 1./((i_n+1.)*(i_n+1.)) );
	In(1) = -0.5/(i_n+1.)*( 1.0/(i_n+1.) + 1 ) + 0.5/i_n;
	for (int n=2; n<N; ++n) In(n) = -0.5/(i_n+1.)*( 1.0/(i_n+1.) + n ) + 0.5*(n-1)*n*Kn(n-2);
      }
      break;
    default : // any negative number
      int k=-power-1;
      if (i_n==0){
	In=0;
	for (int n=k+1; n<N; ++n) In(n) = 1.0/(n-k);
      }else{
	bl::Array<double,1> Kn(N-k); // need In for all nonzero orders
	PolyInt_m1(Kn, i_n, N-k);
	bl::Array<double,1> cp(k+1), cq(k+1);
	cp(0) = ipower( (i_n+1.), -k)/(k+0.0);
	cq(0) = ipower( (i_n+0.), -k)/(k+0.0);
	// cp(0) = 1/(i+1)^k/k                   cq(0) = 1/i^k/k
	// cp(1) = 1/(i+1)^{k-1}/(k*(k-1))       cq(1) = 1/i^{k-1}/(k*(k-1))
	// cp(2) = 1/(i+1)^{k-2}/(k*(k-1)*(k-2)) cq(2) = 1/i^{k-2}/(k*(k-1)*(k-2))
	// cp(k) = 1/(k!)                        cq(k) = 1/(k!)
	for (int i=1; i<k; ++i){
	  cp(i) = cp(i-1) * (i_n+1.)/(k-i);
	  cq(i) = cq(i-1) * (i_n+0.)/(k-i);
	}
	cp(k) = cp(k-1) * (i_n+1.);
	cq(k) = cq(k-1) * (i_n+0.);
	//In(0) = -[ cp(0) - cq(0)];
	//In(1) = -[ cp(0) + ( cp(1) - cq(1) ) ];
	//In(2) = -[ cp(0) + 2*cp(1) + 2*1*( cp(2) - cq(2) ) ];
	//In(3) = -[ cp(0) + 3*cp(1) + 3*2* cp(2) + 3*2*1*(cp(3)-cq(3)) ];
	//In(i) = -[ cp(0) + i*cp(1) + i*(i-1)*cp(i-1) + i*(i-1)*(i-2)*cp(i-2) ... ]
	for (int i=0; i<k; ++i){
	  double dsm=0;
	  double dd=1;
	  for (int j=0; j<i; ++j){ 
	    dsm += dd*cp(j); 
	    dd *= (i-j);     
	  }
	  dsm += dd*(cp(i)-cq(i));
	  In(i) = -dsm;
	}
	for (int n=k; n<N; ++n){
	  double ds = 0.0;
	  double dd0 = 1.0;
	  for (int i=0; i<k; ++i){
	    ds -= dd0*cp(i);
	    dd0 *= (n-i);
	  }
	  In(n) = ds + dd0*cp(k)*Kn(n-k);
	}
      }
    }
    double dsum = 0.0;
    for (int n=0; n<N; ++n) dsum += In(n)*Poly(n);
    return dsum;
  }
}

inline std::ostream& operator<< (std::ostream& stream, const Polynomial& A){
  for (int i=0; i<A.Poly.extent(0); ++i){
    stream << A.Poly(i)<<" ";
    if (i != 0) stream << "x^" << i;
    if (i != A.Poly.extent(0)-1) stream << " + ";
  }
  stream<<std::endl;
  return stream;
}


/*
// Driver program to test above functions
int main()
{
  bl::Array<double,1> A(4);
  // The following array represents polynomial 5 + 10x^2 + 6x^3
  A = 5, 0, 10, 6;
  Polynomial pA(A);
  bl::Array<double,1> B(3);
  // The following array represents polynomial 1 + 2x + 4x^2
  B = 1, 2, 4;
  Polynomial pB(B);
  
  
  std::cout << "First polynomial is ";
  //printPoly(A);
  std::cout << pA << std::endl;
  std::cout << "Second polynomial is ";
  //printPoly(B);
  std::cout << pB << std::endl;
  
  //bl::Array<double,1> prod = multiply(A, B);
  std::cout << "Product polynomial is ";
  //printPoly(prod);

  Polynomial pC = pA * pB;
  std:: cout << pC << std::endl;

  int power = -4;
  for (int i_n=0; i_n<10; i_n++)
    std::cout<< i_n << " "<< pC.CmpIntegral(i_n,power) << std::endl;
  std::cout<<std::endl;
  
  //std::cout<<"Int = "<<std::endl;
  //for (int i_n=0; i_n<10; i_n++)
  //    std::cout<<i_n<<" "<< GetPolyInt(prod, i_n, power)<<std::endl;
  
  return 0;
}
*/
#endif
