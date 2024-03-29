// @Copyright 2007 Kristjan Haule
// 
#ifndef _AVERAGE_
#define _AVERAGE_

#include "complex.h"
#include "mesh.h"
#include "function.h"

class cintpar{
public:
  double p[2], q[2], p0;
  int i;
};

/////////////////////////////////////////////////////////
///////////////////////  T3  ////////////////////////////
/////////////////////////////////////////////////////////
template <class T>
class T3{
  T v[3];
public:
  T3(){};
  T3(const T& d)
  { v[0] = v[1] = v[2] = d;}
  T3(const T& v0, const T& v1, const T& v2)
  { v[0] = v0; v[1] = v1; v[2] = v2;}
  T3& operator=(const T3& m)
  { v[0] = m.v[0]; v[1] = m.v[1]; v[2] = m.v[2]; return *this;}
  T3& operator+=(const T3& m)
  { v[0]+= m.v[0]; v[1]+= m.v[1]; v[2]+= m.v[2]; return *this;}
  T& operator[](int i) { return v[i];}
  const T& operator[](int i) const { return v[i];}
  template <class U>
  friend std::ostream& operator<<(std::ostream& stream, const T3<U>& m);
};

template <class T>
std::ostream& operator<<(std::ostream& stream, const T3<T>& m)
{
  stream<<m[0]<<" "<<m[1]<<" "<<m[2];
  return stream;
}

inline void InterpLeft(double x, const mesh& om, tint& a, cintpar& p)
{
  int i = om.find_(x,a);
  p.i = i;
  p.q[0] = (x-om[i]);
  p.p0 = p.q[0]*om.Delta(i);
  p.p[0] = 0.5*p.q[0]*p.q[0]*om.Delta(i);
  p.q[0] = p.q[0]-p.p[0];
  const double x2=x*x, x02=om[i]*om[i], dx2=x2-x02;
  p.q[1] = 0.5*dx2;
  p.p[1] = ((x2*x-x02*om[i])/3-p.q[1]*om[i])*om.Delta(i);
  p.q[1] = p.q[1]-p.p[1];
}

//! Class that stores some constants used to
//! average a function.
class apar{
public:
  //! C1 = \f$ {2\over (\epsilon_{i+1}-\epsilon_{i-1})(\epsilon_i-\epsilon_{i-1})}\f$
  double C1;
  //! C2 = \f$ {2\over (\epsilon_{i+1}-\epsilon_{i-1})(\epsilon_{i+1}-\epsilon_i)}\f$
  double C2;
  //! dy0 = \f$ \epsilon_{i-1}-u_1\f$
  double dy0;
  //! dy1 = \f$ \epsilon_{i+1}-u_1\f$
  double dy1;
  bool SimpleInterpolation;
  void SetUpCsFirst(double u, const mesh& om);
  void SetUpCs(double u, int i, const mesh& om, double dhom);
  void SetUpCsLast(double u, const mesh& om);
};
inline void apar::SetUpCsFirst(double u, const mesh& om)
{
  C1 = dy0 = 0;
  C2 = 2.0*sqr(om.Delta(0));
  dy1 = om[1]-u;
  SimpleInterpolation = false;
}
inline void apar::SetUpCsLast(double u, const mesh& om)
{
  C1 = 2.0*sqr(om.Delta(om.size()-2));
  dy0 = om[om.size()-2]-u;
  C2 = dy1 = 0;
  SimpleInterpolation = false;
}
inline void apar::SetUpCs(double u, int i, const mesh& om, double dhom)
{
  double dh = om.Dh(i);
  if (dh<0.01*dhom){
    SimpleInterpolation = true;
    return;
  }
  
  SimpleInterpolation = false;
  C1 = om.Delta(i-1)/dh;
  C2 = om.Delta(i)/dh;
  dy0 = om[i-1]-u;
  dy1 = om[i+1]-u;
}

//////////////////////////////////////////////////////////
///////////////////// Class   AvFun   ////////////////////
//////////////////////////////////////////////////////////
/**
 * Calculates the average value of the function on the specified mesh in the following way
 * \f[ <G(\epsilon,u_1)> = 
 *   C_i^1 {\int_{\epsilon_{i-1}}^{\epsilon_i} (u-\epsilon_{i-1}) G(u-u_1) du} +
 *   C_i^2 {\int_{\epsilon_i}^{\epsilon_{i+1}} (\epsilon_{i+1}-u) G(u-u_1) du}
 * \f]
 * \f[ <G(\epsilon,u_1)> = 
 *   C_i^1 \left({\int_{\epsilon_{i-1}-u_1}^{\epsilon_i-u_1}G(v) v dv} - dy_i^0 {\int_{\epsilon_{i-1}-u_1}^{\epsilon_i-u_1} G(v) dv}\right)+
 *   C_i^2 \left(dy_i^1 {\int_{\epsilon_{i}-u_1}^{\epsilon_{i+1}-u_1} G(v) dv} - {\int_{\epsilon_i-u_1}^{\epsilon_{i+1}-u_1}G(v)v dv}\right)
 * \f]
 * \arg \f$ dy_i^0 = \epsilon_{i-1}-u_1 \f$
 * \arg \f$ dy_i^1 = \epsilon_{i+1}-u_1 \f$
 * \arg \f$ C_i^1 = {1/(dh_i d\epsilon_i)}\f$
 * \arg \f$ C_i^2 = {1/(dh_i d\epsilon_{i+1})}\f$
 * \arg \f$ dh_i = {1\over 2} (\epsilon_{i+1}-\epsilon_{i-1})\f$
 * \arg \f$ d\epsilon_i = \epsilon_i - \epsilon_{i-1} \f$
 * \arg djg0 = \f$ {\int_{\epsilon_{i-1}-u_1}^{\epsilon_i-u_1} G(v) v dv}\f$
 * \arg dig0 = \f$ {\int_{\epsilon_{i-1}-u_1}^{\epsilon_i-u_1}  G(v) dv}\f$
 * \arg djg1 = \f$ {\int_{\epsilon_i-u_1}^{\epsilon_{i+1}-u_1} G(v) v dv}\f$
 * \arg dig1 = \f$ {\int_{\epsilon_i-u_1}^{\epsilon_{i+1}-u_1}  G(v) dv}\f$
 * \arg dy0 = \f$ dy_i^0 \f$
 * \arg dy1 = \f$ dy_i^1 \f$
 * \arg C1 = \f$ C_i^1 \f$
 * \arg C2 = \f$ C_i^1 \f$
 */
template <class T>
class AvFun{
  function1D<T3<T> > fi;
  //!  ig1 = \f$ \int_{-\infty}^{\epsilon_i    -u_1} G(v) dv \f$
  T ig1;
  //!  jg1 = \f$ \int_{-\infty}^{\epsilon_i    -u_1} G(v) v dv \f$
  T jg1;
  //! dig0 = \f$ {\int_{\epsilon_{i-1}-u_1}^{\epsilon_i-u_1} G(v) dv}\f$
  T dig0;
  //! dig1 = \f$ {\int_{\epsilon_i-u_1}^{\epsilon_{i+1}-u_1} G(v) dv}\f$
  T dig1;
  //! djg0 = \f$ {\int_{\epsilon_{i-1}-u_1}^{\epsilon_i-u_1} G(v) v dv}\f$
  T djg0;
  //! djg1 = \f$ {\int_{\epsilon_i-u_1}^{\epsilon_{i+1}-u_1} G(v) v dv}\f$
  T djg1;
  T simplev;
public:
  AvFun(const functionb<T>& f0, const mesh& om);
  AvFun();
  void SetUp(const functionb<T>& f0, const mesh& om);
  void InterpolateFirst(const cintpar& ip);
  T InterpolateNext(const cintpar& ip, const apar& ap);
  T InterpolateLast(const apar& ap);
};

template <class T>
void AvFun<T>::SetUp(const functionb<T>& f0, const mesh& om)
{
  fi.resize(om.size());
  SUNCA_LOG(if (om.size()!=f0.size()) std::cerr<<"Function and Mesh not compatible in AvFun"<<std::endl;);
  fi[0][0] = f0[0];
  fi[0][1] = fi[0][2]=0;
  for (int i=1; i<om.size(); i++){
    fi[i][0] = f0[i];
    fi[i][1] = fi[i-1][1] + (fi[i][0]+fi[i-1][0])*(om[i]-om[i-1])/2.;
    fi[i][2] = fi[i-1][2] + (fi[i][0]*(2*om[i]+om[i-1])+fi[i-1][0]*(2*om[i-1]+om[i]))*(om[i]-om[i-1])/6.;
  }
}

template <class T>
AvFun<T>::AvFun(const functionb<T>& f0, const mesh& om) : fi(om.size())
{ SetUp(f0,om);}

template <class T>
AvFun<T>::AvFun(){};

template <class T>
inline void AvFun<T>::InterpolateFirst(const cintpar& ip)
{
  dig1 = ig1 = fi[ip.i][1]+ip.q[0]*fi[ip.i][0]+ip.p[0]*fi[ip.i+1][0];
  djg1 = jg1 = fi[ip.i][2]+ip.q[1]*fi[ip.i][0]+ip.p[1]*fi[ip.i+1][0];
}
template <class T>
inline T AvFun<T>::InterpolateNext(const cintpar& ip, const apar& ap)
{
  T f1 = fi[ip.i][1]+ip.q[0]*fi[ip.i][0]+ip.p[0]*fi[ip.i+1][0];
  T f2 = fi[ip.i][2]+ip.q[1]*fi[ip.i][0]+ip.p[1]*fi[ip.i+1][0];
  dig0 = dig1;
  dig1 = f1-ig1;
  ig1 = f1;
  djg0 = djg1;
  djg1 = f2-jg1;
  jg1 = f2;
  if (ap.SimpleInterpolation) return  fi[ip.i][0]+ip.p0*(fi[ip.i+1][0]-fi[ip.i][0]);
  else return ap.C1 * (djg0 - ap.dy0 * dig0) - ap.C2 * (djg1 - ap.dy1 * dig1);
}
template <class T>
inline T AvFun<T>::InterpolateLast(const apar& ap)
{
  return ap.C1 * (djg1 - ap.dy0 * dig1);
}

#endif
