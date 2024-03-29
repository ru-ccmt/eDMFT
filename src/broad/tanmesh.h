// @Copyright 2007-2017 Kristjan Haule
#ifndef _TAN_MESH
#define _TAN_MESH

#include <vector>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

struct qparams{
  double x0, L;
  int Nw;
};

int tanmesh_g(const gsl_vector* x, void* params, gsl_vector* f){
  double x0 = ((struct qparams *) params)->x0;
  double  L = ((struct qparams *) params)->L;
  int Nw = ((struct qparams *) params)->Nw;
  const double d = gsl_vector_get (x, 0);
  const double w = gsl_vector_get (x, 1);
  gsl_vector_set (f, 0, L-w/tan(d) );
  gsl_vector_set (f, 1, x0-w*tan(M_PI/(2*Nw)-d/Nw) );
  return GSL_SUCCESS;
}

typedef int (*gsl_FNC)(const gsl_vector * x, void * params, gsl_vector * f);

class FindRootn{
  const gsl_multiroot_fsolver_type *T;
  gsl_multiroot_fsolver *s;
  size_t n;
  void* params;
  gsl_vector *x;
  gsl_multiroot_function f;
public:
  FindRootn(size_t n_, gsl_FNC fnc, double x_init[], void* params_) : n(n_), params(params_)
  {
    x = gsl_vector_alloc (n);
    f.f = fnc;
    f.n = n;
    f.params = params;
    for (int i=0; i<n; i++) gsl_vector_set(x, i, x_init[i]);
    T = gsl_multiroot_fsolver_hybrids;
    s = gsl_multiroot_fsolver_alloc (T, n);
    gsl_multiroot_fsolver_set (s, &f, x);
  }
  void print_state (size_t iter, gsl_multiroot_fsolver * s){
    std::clog<<"iter = "<<iter<<" x=";
    for (int i=0; i<n; i++) std::clog<<gsl_vector_get(s->x,i)<<" ";
    std::clog<<" f(x)=";
    for (int i=0; i<n; i++) std::clog<<gsl_vector_get(s->f,i)<<" ";
    std::clog<<std::endl;
  };
  std::vector<double> call(){
    size_t iter = 0;
    //print_state (iter, s);
    int status;
    do{
      iter++;
      //print_state (iter, s);
      status = gsl_multiroot_fsolver_iterate(s);
      //print_state (iter, s);
      if (status)   /* check if solver is stuck */
	break;
      status = gsl_multiroot_test_residual(s->f, 1e-13);
    } while (status == GSL_CONTINUE && iter < 100000);
    //clog<<"status = "<<gsl_strerror (status)<<endl;
    std::vector<double> res(n);
    for (int i=0; i<n; i++) res[i] = gsl_vector_get (s->x, i);
    return res;
  }
  ~FindRootn(){
    gsl_multiroot_fsolver_free (s);
    gsl_vector_free (x);
  }
};

template<typename container>
void GiveTanMesh(container& om, container& dom, double& x0, double L, int Nw)
{
  double tnw = tan(M_PI/(2*Nw));
  if (x0 > L*0.25*Nw*tnw*tnw ){
    x0 = L*0.25*Nw*tnw*tnw-1e-15;
  }
  double d0 = Nw/2.*( tnw - sqrt( tnw*tnw - 4*x0/(L*Nw)));
  double w0 = L*d0;
  double x_init[2] = {d0, w0};
  struct qparams p = {x0,L,Nw};
  
  FindRootn fr(2, &tanmesh_g, x_init, &p);
  std::vector<double> dw = fr.call();
  double d = dw[0];
  double w = dw[1];
  
  om.resize(2*Nw+1);
  dom.resize(2*Nw+1);
  double dh = 1.0/static_cast<double>(2*Nw);
  for (int i=0; i<2*Nw+1; i++){
    double t0 = i/static_cast<double>(2*Nw);
    double t = t0*(M_PI-2*d) - M_PI/2 + d;
    om(i)  = w*tan(t);
    dom(i) = dh*w*(M_PI-2*d)/(cos(t)*cos(t));
  }
}
template<typename container>
void GiveTanMesh0(container& om, double& x0, double L, int Nw)
{
  //Nw -> (Nw-1)/2
  double tnw = tan(M_PI/(Nw-1.));
  if (x0 > L*0.25*(Nw-1.)/2.*tnw*tnw ){
    x0 = L*0.25*(Nw-1.)/2.*tnw*tnw-1e-15;
  }
  double d0 = (Nw-1.)/4.*( tnw - sqrt( tnw*tnw - 8*x0/(L*(Nw-1.))));
  double w0 = L*d0;
  double x_init[2] = {d0, w0};
  struct qparams p = {x0,L,(Nw-1)/2};
  
  FindRootn fr(2, &tanmesh_g, x_init, &p);
  std::vector<double> dw = fr.call();
  double d = dw[0];
  double w = dw[1];

  //clog << "x0=" << x0 << " L=" << L << " Nw=" << Nw << " d=" << d << " w=" << w << endl;
  
  om.resize(Nw);
  double dh = 1.0/static_cast<double>(2*Nw);
  for (int i=0; i<Nw; i++){
    double t = (M_PI/2.-d) * i/(Nw-1.);
    if (i==0) t = (M_PI/2.-d) * 1e-5/(Nw-1.);
    om(i)  = w*tan(t);
  }
}

template<class Array_double_1>
void GiveDoubleExpMesh0(Array_double_1& x, double a, double b, int Nx, double al=2.1){
  for (int i=0; i<Nx; i++){
    double t = (Nx-i-1)/(Nx-1.);
    double x0 = M_PI/2.*sinh(al*t);
    x(i) = b + (a-b) * tanh(x0);
    //dx(i) = (b-a)*al*M_PI/2.*cosh(al*t)/ipower(cosh(x0),2);
  }
}

inline double sqr(double x){return x*x;}

template<class Array_double_1>
void GiveDoubleExpMesh(Array_double_1& tau, Array_double_1& tau_dh, double beta, int Ntau_)
{
  double al=2.1;
  tau.resize(Ntau_);
  tau_dh.resize(Ntau_);
  double Ntau = Ntau_-2;
  double dh = beta/(Ntau-1.);
  for (int i=0; i<Ntau; i++){
    double t = 2*i/(Ntau-1.)-1.;
    double x0 = M_PI/2.*sinh(al*t);
    tau(i+1) = (tanh(x0)+1.)*(beta/2.);
    tau_dh(i+1) = dh*al*M_PI/2.*cosh(al*t)/sqr(cosh(x0));
  }
  tau(0)          = 0.0;
  tau(Ntau_-1)    = beta;
  tau_dh(0)       = tau(1)-tau(0);
  tau_dh(Ntau_-1) = beta-tau(Ntau_-2);
  
  //for (int i=0; i<tau.size()-1; i++) tau.Delta(i) = 1/(tau[i+1]-tau[i]);
  //tau.Delta(tau.size()-1) = 0.0;
}

#endif // _TAN_MESH
