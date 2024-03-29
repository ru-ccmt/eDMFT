// @Copyright 2007 Kristjan Haule
// 
#ifndef __MY_DCOMPLEX__
#define __MY_DCOMPLEX__
#include <cmath>
#include <iostream>
#include <iomanip>
#include <complex>
using namespace std;

class dcomplex{
private:
  double re, im;
public:
  dcomplex (double r = 0, double i = 0): re (r), im (i) { }

  dcomplex& operator+= (const double r);
  dcomplex& operator-= (const double r);
  dcomplex& operator*= (const double r);
  dcomplex& operator/= (const double r);
  
  dcomplex& operator+= (const dcomplex& r);
  dcomplex& operator-= (const dcomplex& r);
  dcomplex& operator*= (const dcomplex& r);
  dcomplex& operator/= (const dcomplex& r);

  double& real() { return re; }
  double& imag() { return im; }
  double real() const { return re; }
  double imag() const { return im; }
  //  dcomplex& conj() { im = -im; return *this;}
  dcomplex conj() const {return dcomplex(re,-im);}
  void Set(double r, double i){re=r; im=i;}
  void Add(double r, double i){re+=r, im+=i;};
  void Multiply_zr(const dcomplex& z, const double& r){ re=z.re*r; im=z.im*r;};
  void Multiply_cr(const dcomplex& z, const double& r){ re=z.re*r; im=-z.im*r;};
  // typecasting
  operator complex<double>() {return complex<double>(re,im);}
private:
  friend dcomplex operator + (const dcomplex& x, double y);
  friend dcomplex operator + (double x, const dcomplex& y);
  friend dcomplex operator + (const dcomplex& x, const dcomplex& y);
  friend dcomplex operator - (const dcomplex& x, double y);
  friend dcomplex operator - (double x, const dcomplex& y);
  friend dcomplex operator - (const dcomplex& x, const dcomplex& y);
  friend dcomplex operator * (const dcomplex& x, double y);
  friend dcomplex operator * (double x, const dcomplex& y);
  friend dcomplex operator * (const dcomplex& x, const dcomplex& y);
  friend dcomplex operator / (const dcomplex& x, double y);
  friend dcomplex operator / (double x, const dcomplex& y);
  friend dcomplex operator / (const dcomplex& x, const dcomplex& y);
  friend bool operator == (const dcomplex& x, double y);
  friend bool operator == (double x, const dcomplex& y);
  friend bool operator == (const dcomplex& x, const dcomplex& y);
  friend bool operator != (const dcomplex& x, double y);
  friend bool operator != (double x, const dcomplex& y);
  friend bool operator != (const dcomplex& x, const dcomplex& y);
  friend dcomplex conj (const dcomplex& r);
  friend std::ostream& operator<< (std::ostream& stream, const dcomplex& r);
  friend std::istream& operator>> (std::istream& stream, dcomplex& r);
  friend double real(const dcomplex& r) { return r.re;}
  friend double imag(const dcomplex& r) { return r.im;}
  friend double& real(dcomplex& r) { return r.re;}
  friend double& imag(dcomplex& r) { return r.im;}
  friend dcomplex sqrt(const dcomplex& z);
  friend dcomplex sqrt_(const dcomplex& z);
  friend dcomplex sqrt(const dcomplex& z1, const dcomplex x2);
  friend bool isnan(const dcomplex& z) { return ::isnan(z.re) || ::isnan(z.im);}
  friend double arg(const dcomplex& x) { return atan2(x.imag(), x.real());}
  friend double norm (const dcomplex& r) { return r.re*r.re+r.im*r.im;}
  friend double abs(const dcomplex& r) { return sqrt(norm(r));}
  friend dcomplex log(const dcomplex& x);
  friend dcomplex  operator- (const dcomplex& r);
  friend double real_product(const dcomplex& z1, const dcomplex& z2) {return z1.re*z2.re-z1.im*z2.im;}
  friend double imag_product(const dcomplex& z1, const dcomplex& z2) {return z1.re*z2.im+z1.im*z2.re;}
  friend dcomplex coth(const dcomplex& z);
  friend dcomplex atan(const dcomplex& z);
};

template <class T>
inline T sqr(const T& x)
{  return x*x; }

inline dcomplex& dcomplex::operator+= (const double r)
{
  re += r;
  return *this;
}

inline dcomplex& dcomplex::operator-= (const double r)
{
  re -= r;
  return *this;
}

inline dcomplex& dcomplex::operator*= (const double r)
{
  re *= r;
  im *= r;
  return *this;
}

inline dcomplex& dcomplex::operator/= (const double r)
{
  re /= r;
  im /= r;
  return *this;
}

inline dcomplex& dcomplex::operator+= (const dcomplex& r)
{
  re += r.re;
  im += r.im;
  return *this;
}

inline dcomplex& dcomplex::operator-= (const dcomplex& r)
{
  re -= r.re;
  im -= r.im;
  return *this;
}

inline dcomplex& dcomplex::operator*= (const dcomplex& r)
{
  double a = re * r.re - im * r.im;
  im = re * r.im + im * r.re;
  re = a;
  return *this;
}

inline dcomplex& dcomplex::operator/= (const dcomplex& r)
{
  double norm_ = norm(r);
  double a = re * r.re + im * r.im;
  im = im * r.re - re * r.im;
  re = a/norm_;
  im /= norm_;
  return *this;
}

inline dcomplex operator + (const dcomplex& x, double y)
{
  return dcomplex(x.re + y, x.im);
}

inline dcomplex operator + (double x, const dcomplex& y)
{
  return dcomplex(x + y.re, y.im);
}

inline dcomplex operator + (const dcomplex& x, const dcomplex& y)
{
  return dcomplex(x.re + y.re, x.im + y.im);
}

inline dcomplex operator - (const dcomplex& x, double y)
{
  return dcomplex(x.re - y, x.im);
}

inline dcomplex operator - (double x, const dcomplex& y)
{
  return dcomplex(x - y.re, -y.im);
}

inline dcomplex operator - (const dcomplex& x, const dcomplex& y)
{
  return dcomplex(x.re - y.re, x.im - y.im);
}

inline dcomplex operator * (const dcomplex& x, double y)
{
  return dcomplex(x.re * y, x.im * y);
}

inline dcomplex operator * (double x, const dcomplex& y)
{
  return dcomplex(x * y.re, x * y.im);
}

inline dcomplex operator * (const dcomplex& x, const dcomplex& y)
{
  return dcomplex(x.re * y.re - x.im * y.im, x.re * y.im + x.im * y.re);
}

inline dcomplex operator / (const dcomplex& x, double y)
{
  return dcomplex(x.re/y, x.im/y);
}

inline dcomplex operator / (double x, const dcomplex& y)  
{
  double xnorm = x/norm(y);
  return dcomplex(xnorm * y.re, -xnorm * y.im);
}

inline dcomplex operator / (const dcomplex& x, const dcomplex& y)
{
  double norm_ = norm(y);
  return dcomplex((x.re * y.re + x.im * y.im)/norm_, (x.im * y.re - x.re * y.im)/norm_);
}

inline bool operator == (const dcomplex& x, double y)
{
  return (x.re==y) && (x.im==0);
}

inline bool operator == (double x, const dcomplex& y)
{
  return (x==y.re) && (y.im==0);
}

inline bool operator == (const dcomplex& x, const dcomplex& y)
{
  return (x.re==y.re) && (x.im==y.im);
}

inline bool operator != (const dcomplex& x, double y)
{
  return !(operator==(x,y));
}

inline bool operator != (double x, const dcomplex& y)
{
  return !(operator==(x,y));
}

inline bool operator != (const dcomplex& x, const dcomplex& y)
{
  return !(operator==(x,y));
}

inline dcomplex operator- (const dcomplex& r)
{
  return dcomplex(-r.re,-r.im);
}

inline dcomplex conj (const dcomplex& r)
{
  return dcomplex(r.re, -r.im);
}

inline dcomplex log(const dcomplex& x)
{
  return dcomplex(log(norm(x))*0.5, arg(x));
}

inline dcomplex sqrt(const dcomplex& z1, const dcomplex z2)
{
  
  double fi1 = atan2(z1.im, z1.re);
  double fi2 = atan2(z2.im, z2.re);
  double fi = 0.5*(fi1+fi2);
  double r = sqrt(sqrt((sqr(z1.re)+sqr(z1.im))*(sqr(z2.re)+sqr(z2.im))));
  return dcomplex(r*cos(fi),r*sin(fi));
}

inline dcomplex sqrt(const dcomplex& x)
{
  double r = abs (x);
  double nr, ni;
  if (r == 0.0)
    nr = ni = r;
  else if (real (x) > 0)
    {
      nr = sqrt (0.5 * (r + real (x)));
      ni = imag (x) / nr / 2;
    }
  else
    {
      ni = sqrt (0.5 * (r - real (x)));
      if (imag (x) < 0)
        ni = - ni;
      nr = imag (x) / ni / 2;
    }
  return dcomplex (nr, ni);
}

inline dcomplex sqrt_(const dcomplex& x)
{
  double r = sqrt(abs(x));
  double fi = atan2(x.im, x.re);
  return dcomplex(r*cos(0.5*fi), r*sin(0.5*fi));
}

inline dcomplex coth(const dcomplex& z)
{
  double ex = exp(z.real()), emx = exp(-z.real());
  double cy = cos(z.imag()), sy = sin(z.imag());
  double u1 = ex+emx, u2 = ex-emx;
  double im = 1./(sqr(u2*cy)+sqr(u1*sy));
  return dcomplex (u1*u2*im, cy*sy*(sqr(u2)-sqr(u1))*im);
}

inline dcomplex atan(const dcomplex& z)
{
  double x = z.real(), y = z.imag(), x2 = x*x, y2=y*y;
  double den = 1/(x2+sqr(y+1.));
  double re = 0.5*atan2(2*x*den,(1-x2-y2)*den);
  double im = -0.25*log(den*(x2+sqr(y-1.)));
  return dcomplex(re,im);
}
  
inline std::ostream& operator<< (std::ostream& stream, const dcomplex& r){
  int width = stream.width();
  stream << std::setw(width)<< r.re << " " << std::setw(width) << r.im << " ";
  return stream;
}

inline std::istream& operator>> (std::istream& stream, dcomplex& r){
  stream >> r.re;
  stream >> r.im;
  return stream;
}

#endif
