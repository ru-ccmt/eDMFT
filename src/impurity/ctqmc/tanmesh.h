// @Copyright 2007-2017 Kristjan Haule
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include "smesh.h"

struct rparams{
  double x0, L;
  int Nw;
};

int tanmesh_f(const gsl_vector* x, void* params, gsl_vector* f){
  double x0 = ((struct rparams *) params)->x0;
  double  L = ((struct rparams *) params)->L;
  int Nw = ((struct rparams *) params)->Nw;
  const double d = gsl_vector_get (x, 0);
  const double w = gsl_vector_get (x, 1);
  gsl_vector_set (f, 0, L-w/tan(d) );
  gsl_vector_set (f, 1, x0-w*tan(M_PI/(2*Nw)-d/Nw) );
  return GSL_SUCCESS;
}

typedef int (*gsl_FNC)(const gsl_vector * x, void * params, gsl_vector * f);

class FindRoot{
  const gsl_multiroot_fsolver_type *T;
  gsl_multiroot_fsolver *s;
  size_t n;
  void* params;
  gsl_vector *x;
  gsl_multiroot_function f;
public:
  FindRoot(size_t n_, gsl_FNC fnc, double x_init[], void* params_) : n(n_), params(params_)
  {
    x = gsl_vector_alloc (n);
    f.f = fnc;
    f.n = n;
    f.params = params;
    for (int i=0; i<n; i++) gsl_vector_set (x, i, x_init[i]);
    T = gsl_multiroot_fsolver_hybrids;
    s = gsl_multiroot_fsolver_alloc (T, 2);
    gsl_multiroot_fsolver_set (s, &f, x);
  }
  void print_state (size_t iter, gsl_multiroot_fsolver * s){
    clog<<"iter = "<<iter<<" x=";
    for (int i=0; i<n; i++) clog<<gsl_vector_get(s->x,i)<<" ";
    clog<<" f(x)=";
    for (int i=0; i<n; i++) clog<<gsl_vector_get(s->f,i)<<" ";
    clog<<endl;
  };
  vector<double> call(){
    size_t i, iter = 0;
    //print_state (iter, s);
    int status;
    do{
      iter++;
      status = gsl_multiroot_fsolver_iterate (s);
      //print_state (iter, s);
      if (status)   /* check if solver is stuck */
	break;
      status = gsl_multiroot_test_residual (s->f, 1e-7);
    } while (status == GSL_CONTINUE && iter < 1000);
    //clog<<"status = "<<gsl_strerror (status)<<endl;
    
    vector<double> res(n);
    for (int i=0; i<n; i++) res[i] = gsl_vector_get (s->x, i);
    return res;
  }
  ~FindRoot(){
    gsl_multiroot_fsolver_free (s);
    gsl_vector_free (x);
  }
};


double sqr(double x){return x*x;}

void GiveTanMesh(mesh1D& om, double& x0, double L, int Nw)
{
  double tnw = tan(M_PI/(2*Nw));
  if (x0 > L*0.25*Nw*sqr(tnw) ){
    x0 = L*0.25*Nw*sqr(tnw)-1e-15;
  }
  double d0 = Nw/2.*( tnw - sqrt( sqr(tnw) - 4*x0/(L*Nw)));
  double w0 = L*d0;
  double x_init[2] = {d0, w0};
  struct rparams p = {x0,L,Nw};

  FindRoot fr(2, &tanmesh_f, x_init, &p);
  vector<double> dw = fr.call();
  double d = dw[0];
  double w = dw[1];
  
  om.resize(2*Nw+1);
  double dh = 1.0/static_cast<double>(2*Nw);
  for (int i=0; i<2*Nw+1; i++){
    double t0 = i/static_cast<double>(2*Nw);
    double t = t0*(M_PI-2*d) - M_PI/2 + d;
    om[i]  = w*tan(t);
    om.Dh(i) = dh*w*(M_PI-2*d)/sqr(cos(t));
  }
  
  om.Delta(0) =  1.0/(om[1]-om[0]);
  for (int i=1; i<om.size()-1; i++) om.Delta(i) = 1/(om[i+1]-om[i]);
  om.Delta(om.size()-1) = 0.0;
}
void GiveDoubleExpMesh(mesh1D& tau, double beta, int Ntau_)
{
  double al=2.1;
  tau.resize(Ntau_);
  
  double Ntau = Ntau_-2;
  double dh = beta/(Ntau-1.);
  for (int i=0; i<Ntau; i++){
    double t = 2*i/(Ntau-1.)-1.;
    double x0 = M_PI/2.*sinh(al*t);
    tau[i+1] = (tanh(x0)+1.)*(beta/2.);
    tau.Dh(i+1) = dh*al*M_PI/2.*cosh(al*t)/sqr(cosh(x0));
  }
  tau[0]          = 0.0;
  tau[Ntau_-1]    = beta;
  tau.Dh(0)       = tau[1]-tau[0];
  tau.Dh(Ntau_-1) = beta-tau[Ntau_-2];
  
  for (int i=0; i<tau.size()-1; i++) tau.Delta(i) = 1/(tau[i+1]-tau[i]);
  tau.Delta(tau.size()-1) = 0.0;
}
