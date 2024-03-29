// @Copyright 2007 Kristjan Haule
// 
#include <iostream>
#include <cmath>

class Number{
public:
  double mantisa;
  double exponent;
  Number(double mantisa_=0, double exponent_=0) : mantisa(mantisa_), exponent(exponent_){};
  Number(const Number& z){ mantisa=z.mantisa; exponent=z.exponent;}
  Number& operator=(const Number& z){ mantisa=z.mantisa; exponent=z.exponent; return *this;}
  double dbl() const {return mantisa*exp(exponent);}
  double exp_dbl() const { return log(abs(mantisa))+exponent;}
  Number& operator+= (const Number& z);
  Number& operator-= (const Number& z);
  Number& operator*= (const Number& z);
  Number& operator/= (const Number& z);
  
  Number& balance();
  
  friend Number operator* (const Number& x, const Number& y) {return Number(x.mantisa*y.mantisa, x.exponent+y.exponent);};
  friend Number operator* (const Number& x, double y) {return Number(x.mantisa*y, x.exponent);};
  friend Number operator+ (const Number& x, const Number& y);
  friend Number operator- (const Number& x, const Number& y);
  friend Number operator/ (const Number& x, const Number& y) {return Number(x.mantisa/y.mantisa, x.exponent-y.exponent);};
  friend double divide(const Number& x, const Number& y) {return x.mantisa/y.mantisa*exp(x.exponent-y.exponent);};
  
  friend bool operator == (const Number& x, const Number& y) {return x.mantisa==y.mantisa && x.exponent==y.exponent;};
  friend bool operator == (const Number& x, double y) {return x.mantisa*exp(x.exponent) == y;}
  friend bool operator == (double x, const Number& y) {return y==x;};
  
  friend bool operator != (double x, const Number& y) {return !(x==y);}
  friend bool operator != (const Number& x, double y) {return !(x==y);}
  friend bool operator != (const Number& x, const Number& y) { return !(x==y);}
  friend std::ostream& operator<< (std::ostream& stream, const Number& z) { stream<<z.mantisa<<" "<<z.exponent<<" "; return stream;};
  friend Number sqrt(const Number& z) { return Number(sqrt(z.mantisa), z.exponent/2);};
  friend bool isnan(const Number& z) { return ::isnan(z.mantisa) || ::isnan(z.exponent);}
  friend Number abs(const Number& z) { return Number(fabs(z.mantisa),z.exponent);}
  friend double log(const Number& z) { return log(z.mantisa)+z.exponent;};
  friend Number  operator- (const Number& z) { return Number(-z.mantisa,z.exponent);};

  friend bool operator > (const Number& x, const Number& y) {
    return (x.exponent > y.exponent) ? (x.mantisa > y.mantisa * exp(y.exponent-x.exponent)) : (x.mantisa * exp(x.exponent-y.exponent) > y.mantisa);
  }
  friend bool operator < (const Number& x, const Number& y) {
    return (x.exponent > y.exponent) ? (x.mantisa < y.mantisa * exp(y.exponent-x.exponent)) : (x.mantisa * exp(x.exponent-y.exponent) < y.mantisa);
  }
  friend bool operator >= (const Number& x, const Number& y) { return !operator<(x, y); }
  friend bool operator <= (const Number& x, const Number& y) { return !operator>(x, y); }
  
};

Number& Number::operator+= (const Number& z)
{
  if (z.exponent<exponent)
    mantisa += z.mantisa*exp(z.exponent-exponent);
  else {
    mantisa = z.mantisa+mantisa*exp(exponent-z.exponent);
    exponent=z.exponent;
  }
  return *this;
}
Number& Number::operator-= (const Number& z)
{
  if (z.exponent<exponent)
    mantisa -= z.mantisa*exp(z.exponent-exponent);
  else{
    mantisa = -z.mantisa+mantisa*exp(exponent-z.exponent);
    exponent=z.exponent;
  }
  return *this;
} 
Number& Number::operator*= (const Number& z)
{
  mantisa *= z.mantisa;
  exponent += z.exponent;
  return *this;
}
Number& Number::operator/= (const Number& z)
{
  mantisa /= z.mantisa;
  exponent -= z.exponent;
  return *this;
}
Number& Number::balance()
{
  if (mantisa==0) {exponent=0; return *this;}
  exponent += log(fabs(mantisa));
  mantisa = (mantisa>0) ? 1 : -1;
  return *this;
}

Number operator+ (const Number& x, const Number& y)
{ return (x.exponent>y.exponent) ? Number(x.mantisa+y.mantisa*exp(y.exponent-x.exponent), x.exponent) : Number(y.mantisa+x.mantisa*exp(x.exponent-y.exponent), y.exponent);}
Number operator- (const Number& x, const Number& y)
{ return (x.exponent>y.exponent) ? Number(x.mantisa-y.mantisa*exp(y.exponent-x.exponent), x.exponent) : Number(-y.mantisa+x.mantisa*exp(x.exponent-y.exponent), y.exponent);}

Number balance(const Number& x)
{
  if (x.mantisa == 0) return Number(0.0, 0.0);
  return Number( (x.mantisa > 0) ? 1. : -1., x.exponent + log(fabs(x.mantisa)) );
}

Number logDet(const function2D<double>& A)
{
  /*
  ** solve Ax=B using A=UDU' factorization, D is placed in A
  ** where A represents a qxq matrix, b a 1xq vector
  */
  function1D<int> ipiv(A.size_Nd());
  function2D<double> D(A);
  int info = xgetrf(D.size_Nd(), D.MemPt(), D.fullsize_Nd(), ipiv.MemPt());
  if( info < 0 ) { cerr<<"In Det and dgetrf the "<<-info<<"-th argument had an illegal value"<<endl;}
  if( info > 0 ) { cerr<<"In Det and dgetrf the U("<<info<<") is exactly zero and might cause a problem."; }
  /*
  ** compute the determinant det = det(A)
  ** if ipiv[i] > 0, then D(i,i) is a 1x1 block diagonal
  ** if ipiv[i] = ipiv[i-1] < 0, then D(i-1,i-1),
  ** D(i-1,i), and D(i,i) form a 2x2 block diagonal
  */
  double logdet = 0.0;
  int sign=0;
  for (int i=0; i<A.size_N(); i++){
    if (D(i,i)<0){
      ++sign;
      logdet += log(-D(i,i));
    }else{
      logdet += log(D(i,i));
    }
  }
  int change_sign = 0;
  for (int i=0; i<A.size_N(); i++) change_sign += (ipiv[i] != (i+1));
  int overal_sign = 1;
  if ( (sign + change_sign) %2 ) overal_sign = -1;
  return Number(overal_sign,logdet);
}
