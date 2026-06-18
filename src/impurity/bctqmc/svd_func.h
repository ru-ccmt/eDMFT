// @Copyright 2007-2020 Kristjan Haule
#include <iostream>
#include <fstream>
#include <vector>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <blitz/array.h>
#include <time.h>
#include "util.h"
#include "bblas.h"
#include "interpolate.h"

#ifdef _MPI
#include <mpi.h>
#endif

using namespace std;
namespace bl = blitz;

inline double fermi_kernel(double t, double w, double beta)
{
  double x = beta*w/2;
  double y = 2.*t/beta-1.;
  if (x>100.) return exp(-x*(y+1.));
  if (x<-100) return exp(x*(1.-y));
  return exp(-x*y)/(2*cosh(x));
}
inline double bose_kernel(double t, double w, double beta)
{
  double x = beta*w/2;
  double y = 2.*t/beta-1.;
  if (x>200.) return w*exp(-x*(y+1.));
  if (x<-200.) return -w*exp(x*(1.-y));
  return w*exp(-x*y)/(2*sinh(x));
}
inline void create_log_mesh(bl::Array<int,1>& ind_om, int nom_all, int nom, int ntail_)
{ //Creates logarithmic mesh on Matsubara axis
  //     Takes first istart points from mesh om and the rest of om mesh is replaced by ntail poinst redistribued logarithmically.
  //     Input:
  //         om      -- original long mesh
  //         nom     -- number of points not in the tail
  //         ntail   -- tail replaced by ntail points only
  //     Output:
  //         ind_om  -- index array which conatins index to kept Matsubara points
  int istart = min(nom, nom_all);
  int ntail = min(ntail_, nom_all-istart);
  ind_om.resize(istart+ntail);
  double alpha = log((nom_all-1.)/istart)/(ntail-1.);
  for (int i=0; i< istart; i++) ind_om(i)=i;
  int ii=istart;
  for (int i=0; i < ntail; i++){
    int t = int(istart*exp(alpha*i)+0.5);
    if (t != ind_om(ii-1))
      ind_om(ii++) = t;
  }
  ind_om.resize(ii);
}

inline int SymEigensystem(bl::Array<double,2>& H, bl::Array<double,1>& E)
{
  static bl::Array<double,1> work(4*H.extent(0));
  //cout << "H.size=" << H.extent(0) << " H.lda="<< H.extent(1) << " work.size="<< work.size() << endl;
  xsyev(H.extent(0), H.data(), H.extent(1), E.data(), work.data(), work.size(), "V");
  return 0;
}

class SVDFunc{
public:
  std::vector<Spline1D<double> > fU, fU_bose; // SVD functions in imaginary time 
  bl::Array<double,1> tau;                // imaginary time optimized mesh
  int lmax, lmax_bose;                // cutoff for the singular values
  int k_ovlp;                         // 2^{k_ovlp}+1 is the number of integration points in romberg method
  bl::Array<double,1> om;                        // real axis mesh for SVD
  bl::Array<double,1> S, S_bose;          // eigenvalues of SVD decomposition
  bl::Array<double,2> Vt, Vt_bose;        // real axis basis functions
  bl::Array<double,1> dxi;                // parts of fU splines stored in more efficient way for fast interpolation
  bl::Array<double,3> ff2;                // fU splines stored in different way for fast interpolation
  double beta, L, x0;                 // parameters of the mesh should be remembered for reproducibility
  //
  SVDFunc() : lmax(0)  {  };
  //
  void Init(string,double){
    cout<<"ERROR: svdf.Init should not be called"<<endl;
    exit(1);
  }
  bool operator()() //  was it ever initialized?
  { return lmax>0; }
  double operator()(int l, int i){
    if (l<0) l=lmax-l;
    if (i<0) i=tau.size()+i;
    return fU[l][i];
  }
  double operator()(int l, int i) const{
    if (l<0) l=lmax-l;
    if (i<0) i=tau.size()+i;
    return fU[l][i];
  }
  Spline1D<double>& operator[](int l){
    return fU[l];
  }
  inline vector<int> get_lll(int ii, int lmax_=0)
  { // conversion from a single combined index to the three l indices.
    //int ii = l0 + lmax*l1 + lmax*lmax*lb;
    if (lmax_==0) lmax_ = lmax;
    int lb_ = ii / (lmax_*lmax_);
    int ir = ii % (lmax_*lmax_);
    int l1_ = ir / lmax_;
    int l0_ = ir % lmax_;
    int tmp[] = { l0_, l1_, lb_ };
    return vector<int>(tmp, tmp+3 );
  }

  void _cmp_(double beta, int& lmax, std::vector<Spline1D<double> >& fU, bl::Array<double,1>& S, bl::Array<double,2>& Vt, bl::Array<double,2>& K,
	     const bl::Array<double,1>& tau_dh, ostream& clog, bool Print=true, double smallest=1e-12){
    clock_t t0 = clock();
    int Nl = min(om.size(),tau.size());
    S.resize(Nl);
    Vt.resize(om.size(),Nl);
    bl::Array<double,2> U(Nl,tau.size());
    
    int lwork;
    if (tau.size()>om.size()*11/6 || om.size()>tau.size()*11/6){
      int n = min(tau.size(),om.size());
      lwork = n*(4*n + 8);
    }else{
      int n = min(tau.size(),om.size());
      int mx = max(tau.size(),om.size());
      lwork =  3*n + max(mx, n*(3*n+5));
    }					    
    bl::Array<double,1> work(lwork);
    bl::Array<int,1> iwork(8*min(om.size(),tau.size()));
    int info=0;
    info = xgesdd(true, tau.size(), om.size(), K.data(), K.extent(1), S.data(), U.data(), U.extent(1), Vt.data(), Vt.extent(1), work.data(), lwork, work.data(), iwork.data());
  
    double dt = static_cast<double>( clock() - t0 )/CLOCKS_PER_SEC;
    if (Print) clog<<"svd time="<<dt<<endl;

    for (int l=0; l<lmax; l++)
      if (fabs(S(l))<smallest){
	lmax=l;
	if (Print) clog<<"lmax reduced to "<<lmax<<endl;
	break;
      }
    if (Print) clog<<"last and first singular value="<<S(lmax-1)<<" " << S(0) << endl;
    
    for (int i=0; i<tau.size(); i++){
      for (int l=0; l<Nl; l++)
	U(l,i) *= 1./sqrt(tau_dh(i));
    }
    // Splines functions, which are singular values of the Kernel.
    vector<Spline1D<double> > fu(lmax);
    for (int l=0; l<lmax; l++){
      fu[l].resize(tau.size());
      for (int i=0; i<tau.size(); i++) fu[l][i] = U(l,i);
      fu[l].splineIt(tau);
    }
    // Calculates overlap between spline interpolations
    bl::Array<double,2> overlap(lmax,lmax);
    cmpOverlap(overlap, fu, tau, beta, k_ovlp);
    // Calculates eigensystem of the overlap
    bl::Array<double,1> oE(lmax);
    SymEigensystem(overlap, oE);
    bl::Array<double,2> sou(lmax,lmax);// sou = 1/sqrt(overlap)
    // Computes  1/sqrt(overlap)
    for (int i=0; i<lmax; i++){
      for (int j=0; j<lmax; j++){
	double dsum =0.0;
	for (int k=0; k<lmax; k++)
	  dsum += overlap(k,i)*1/sqrt(oE(k))*overlap(k,j);
	sou(i,j) = dsum;
      }
    }
    // Prepares new spline functions, which are more precisely orthonormal.
    fU.resize(lmax);
    for (int l=0; l<lmax; l++){
      fU[l].resize(tau.size());
      for (int i=0; i<tau.size(); i++){
	double ul=0;
	for (int lp=0; lp<lmax; lp++)  ul += sou(l,lp)*fu[lp][i];
	fU[l][i] = ul;
      }
      fU[l].splineIt(tau);
    }
    
  }
  void cmp(string statistics, double beta_, int lmax_, int Ntau, double x0_, double L_, int Nw, std::ostream& clog, bool Print=true, int k_ovlp_=5){
    k_ovlp=k_ovlp_;
    beta = beta_;
    L = L_;
    x0 = x0_;
    
    bl::Array<double,1> om_dh;
    GiveTanMesh(om, om_dh, x0, L, Nw);
    bl::Array<double,1> tau_dh;
    GiveDoubleExpMesh(tau, tau_dh, beta, Ntau);
    
    if (statistics=="fermi" || statistics=="both"){
      lmax = lmax_;
      clock_t t0 = clock();
      bl::Array<double,2> K(om.size(),tau.size());
      for (int j=0; j<om.size(); j++){
	for (int i=0; i<tau.size(); i++){
	  K(j,i) = fermi_kernel(tau(i),om(j),beta);
	  K(j,i) *= om_dh(j)*sqrt(tau_dh(i));
	}
      }
      double dt = static_cast<double>( clock() - t0 )/CLOCKS_PER_SEC;
      if (Print) clog<<"setup time="<<dt<<endl;
      _cmp_(beta, lmax, fU, S, Vt, K, tau_dh, clog, Print);
      //cout<<"S_fermi="<<S<<endl;
    }
    
    if (statistics=="bose" || statistics=="both"){
      lmax_bose = lmax_;
      clock_t t0 = clock();
      bl::Array<double,2> K(om.size(),tau.size());
      for (int j=0; j<om.size(); j++){
	for (int i=0; i<tau.size(); i++){
	  K(j,i) = bose_kernel(tau(i),om(j),beta);
	  K(j,i) *= om_dh(j)*sqrt(tau_dh(i));
	}
      }
      double dt = static_cast<double>( clock() - t0 )/CLOCKS_PER_SEC;
      if (Print) clog<<"setup time="<<dt<<endl;
      
      _cmp_(beta, lmax_bose, fU_bose, S_bose, Vt_bose, K, tau_dh, clog);
      //cout<<"S_bose="<<S_bose<<endl;
    }
  }
  void Print(const std::string& filename, const std::string& statistics="fermi"){
    std::ofstream pul(filename.c_str());
    pul<<"# "<<lmax<<" "<<tau.size()<<" "<<L<<" "<<x0<<"  "<<(om.size()/2)<<endl;
    pul.precision(16);
    for (int i=0; i<tau.size(); i++){
      pul<< tau(i) <<"  ";
      if (statistics=="fermi")
	for (int l=0; l<lmax; l++) pul<<fU[l][i]<<"  ";
      else
	for (int l=0; l<lmax_bose; l++) pul<<fU_bose[l][i]<<"  ";
      pul<<endl;
    }
  }

  double cdot(int l1, int l2, const vector<Spline1D<double> >& fu, int k=5){
    int M = static_cast<int>(pow(2.0,k))+1;    // Number of points in Romberg integration routine
    bl::Array<double,1> utu(M); // storage for Romberg
    
    int nt = tau.size();
    double odd_1 = fu[l1][0]*fu[l1][nt-1];
    double odd_2 = fu[l2][0]*fu[l2][nt-1];
    if (odd_1*odd_2 < 0) return 0.0; // one is even and one is odd    
    
    double oij=0; // overlap between u_{l1} and u_{l2} functions
    for (int i=0; i<tau.size()-1; i++){
      // here we integrate with a fine mesh of M points using romberg routine
      double a = tau(i);  // integral only between t_i and t_{i+1} points
      double b = tau(i+1);
      double fa = fu[l1][i]   * fu[l2][i];
      double fb = fu[l1][i+1] * fu[l2][i+1];
      utu(0)   = fa;
      utu(M-1) = fb;
      int ia=0;
      for (int j=1; j<M-1; j++){
	intpar p = InterpLeft( a + (b-a)*j/(M-1.0), ia, tau);
	utu(j) = fu[l1](p) * fu[l2](p);
      }
      oij += romberg(utu, b-a);
    }
    return oij;
  }
  void Read(const string& filename, const string& statistics, int lmax_=0, int k_ovlp_=5){
    k_ovlp=k_ovlp_;
    
    ifstream pul(filename.c_str());
    string ch;
    int tau_size, Nw, lmax_read;
    pul>>ch>>lmax_read>>tau_size>>L>>x0>>Nw;
    pul.ignore(2000,'\n');
    int _lmax_ = max(lmax_,lmax_read);
    
    tau.resize(tau_size);
    
    vector<Spline1D<double> > *pfU;
    if (statistics=="fermi"){
      lmax = _lmax_;
      pfU = &fU;
    }else{
      lmax_bose = _lmax_;
      pfU = &fU_bose;
    }
    
    pfU->resize(_lmax_);
    for (int l=0; l<_lmax_; l++){
      (*pfU)[l].resize(tau.size());
      for (int it=0; it<tau.size(); it++) (*pfU)[l][it] = 0;
    }

    for (int it=0; it<tau.size(); it++){
      pul>>tau(it);
      for (int l=0; l<lmax_read; l++) pul>>(*pfU)[l][it];
      pul.ignore(2000,'\n');
    }
    beta = tau(tau.size()-1);
    //clog<<"beta="<<beta<<endl;
    
    bl::Array<double,1> om_dh;
    GiveTanMesh(om, om_dh, x0, L, Nw);
    bl::Array<double,1> tau_dh;
    GiveDoubleExpMesh(tau, tau_dh, beta, tau.size());
    
    // Splining functions from the file
    for (int l=0; l<lmax_read; l++) (*pfU)[l].splineIt(tau);
     // double dsum=0.0;
     // for (int l1=0; l1<lmax_read; l1++){
     //   for (int l2=0; l2<lmax_read; l2++){
     // 	double o12 = cdot(l1, l2, *pfU);
     // 	if (l1!=l2) dsum += fabs(o12);
     // 	else dsum += fabs(o12-1.);
     // 	clog<<"<"<<l1<<"|"<<l2<<">="<<o12<<endl;
     //   }
     // }
     //cout<<"First overlap="<<dsum<<endl;
    if (lmax_ > lmax_read){ // Need to add a few functions, not just those which were read
      bl::Array<double,2> K(om.size(),tau.size());
      for (int j=0; j<om.size(); j++){
	for (int i=0; i<tau.size(); i++){
	  if (statistics=="fermi")
	    K(j,i) = fermi_kernel(tau(i),om(j),beta);
	  else
	    K(j,i) = bose_kernel(tau(i),om(j),beta);
	  K(j,i) *= om_dh(j)*sqrt(tau_dh(i));
	}
      }
      int Nl = min(om.size(),tau.size());
      S.resize(Nl);
      Vt.resize(om.size(),Nl);
      bl::Array<double,2> U(Nl,tau.size());
      int lwork;
      if (tau.size()>om.size()*11/6 || om.size()>tau.size()*11/6){
	int n = min(tau.size(),om.size());
	lwork = n*(4*n + 8);
      }else{
	int n = min(tau.size(),om.size());
	int mx = max(tau.size(),om.size());
	lwork =  3*n + max(mx, n*(3*n+5));
      }					    
      bl::Array<double,1> work(lwork);
      bl::Array<int,1> iwork(8*min(om.size(),tau.size()));
      int info=0;
      info = xgesdd(true, tau.size(), om.size(), K.data(), K.extent(1), S.data(), U.data(), U.extent(1), Vt.data(), Vt.extent(1), work.data(), lwork, work.data(), iwork.data());
      
      clog<<"info="<<info<<endl;
      clog<<"first singular value="<<S(0)<<endl;
      clog<<"last and first singular value="<<S(_lmax_-1)<<" "<< S(0) << endl;
      //clog<<"singular values=" << S << endl;
      
      for (int i=0; i<tau.size(); i++)
	for (int l=0; l<Nl; l++)
	  U(l,i) *= 1./sqrt(tau_dh(i));

      // Additional functions added. Splines functions, which are singular values of the Kernel.
      int n = tau.size();
      for (int l=lmax_read; l<_lmax_; l++){
	for (int i=0; i<tau.size(); i++) (*pfU)[l][i] = U(l,i);
	(*pfU)[l].splineIt(tau);
      }
      // Grahm-Shmidt for the new components
      for (int l=lmax_read; l<_lmax_; l++){
	for (int i=0; i<l; i++){
	  double o_i_l = cdot(i,l, *pfU);
	  //cout<<"l="<<l<<" i="<<i<<" <v_i|v_l>="<<o_i_l<<endl;
	  for (int it=0; it<tau.size(); it++) (*pfU)[l][it] -= (*pfU)[i][it]*o_i_l;
	}
	double o_l_l = cdot(l,l, *pfU);
	//cout<<"l="<<l<<" l="<<l<<" <v_l|v_l>="<<o_l_l<<endl;
	for (int it=0; it<tau.size(); it++) (*pfU)[l][it] *= 1./sqrt(o_l_l);
	(*pfU)[l].splineIt(tau);
      }
      /*
      double dsum=0.0;
      for (int l1=0; l1<_lmax_; l1++){
       	for (int l2=0; l2<_lmax_; l2++){
       	  double o12 = cdot(l1, l2, *pfU);
       	  if (l1!=l2) dsum += fabs(o12);
       	  else dsum += fabs(o12-1.);
       	  clog<<"<"<<l1<<"|"<<l2<<">="<<o12<<endl;
       	}
      }
      cout<<"Second overlap="<<dsum<<endl;
      */
    }
  }
  double CheckOverlap(){
    // Now we check how well is the orthogonality obeyed after reorthogonalization
    bl::Array<double,2> overlap(lmax,lmax);
    cmpOverlap(overlap, fU, tau, beta, k_ovlp);
    double dsum=0;
    for (int l1=0; l1<lmax; l1++){
      for (int l2=0; l2<lmax; l2++){
	if (l1!=l2 && overlap(l1,l2)!=0) dsum += fabs(overlap(l1,l2));
	if (l1==l2) dsum += fabs(overlap(l1,l2)-1.0);
      }
    }
    return dsum;
  }
  void cmpOverlap(bl::Array<double,2>& overlap, const vector<Spline1D<double>>& fu, const bl::Array<double,1>& tau, double beta, int k=5){
    int lmax = fu.size();
    int M = static_cast<int>(pow(2.0,k))+1;    // Number of points in Romberg integration routine
    bl::Array<double,1> utu(M); // storage for Romberg
    overlap.resize(lmax,lmax);
    overlap=0.0;
    int nt = tau.size();
    for (int l1=0; l1<lmax; l1++){
      double odd_1 = fu[l1][0]*fu[l1][nt-1];
      for (int l2=0; l2<=l1; l2++){
	double odd_2 = fu[l2][0]*fu[l2][nt-1];
	if (odd_1*odd_2 > 0){// both are even or odd
	  double oij=0; // overlap between u_{l1} and u_{l2} functions
	  for (int i=0; i<tau.size()-1; i++){
	    // here we integrate with a fine mesh of M points using romberg routine
	    double a = tau(i);  // integral only between t_i and t_{i+1} points
	    double b = tau(i+1);
	    double fa = fu[l1][i]   * fu[l2][i];
	    double fb = fu[l1][i+1] * fu[l2][i+1];
	    utu(0)   = fa;
	    utu(M-1) = fb;
	    int ia=0;
	    for (int j=1; j<M-1; j++){
	      //intpar p = tau.InterpLeft( a + (b-a)*j/(M-1.0), ia);
	      intpar p = InterpLeft(a + (b-a)*j/(M-1.0), ia, tau);
	      utu(j) = fu[l1](p) * fu[l2](p);
	    }
	    oij += romberg(utu, b-a);
	    //inline double romberg(const bl::Array<double,1>& ff, double dh)
	  }
	  overlap(l1,l2) = oij;
	  overlap(l2,l1) = oij;
	}
      }
    }
  }
  void CreateGfromCoeff(bl::Array<double,1>& gf, const bl::Array<double,1>& gl, const string& statistics="fermi")const{
    const vector<Spline1D<double> >* pfU = (statistics=="fermi") ? &fU : &fU_bose;
    gf.resize(tau.size());
    gf = 0;
    for (int l=0; l<lmax; l++)
      for (int i=0; i<tau.size(); i++)
	gf(i) += gl(l)*(*pfU)[l][i];
  }
  void CreateSplineFromCoeff(Spline1D<double>& gf, const bl::Array<double,1>& gl, const string& statistics="fermi") const{
    gf.resize(tau.size());
    CreateGfromCoeff(gf.f, gl, statistics);
    gf.splineIt(tau);
  }
  void MatsubaraFrequency(bl::Array<complex<double>,1>& giom, const Spline1D<double>& gf, const bl::Array<double,1> iom, const string& statistics="fermi")const{
    giom.resize(iom.extent(0));
    double one = (statistics=="fermi") ? 1.0 : 0.0;
    for (int in=0; in<iom.extent(0); in++){
      //double iom = (2*in+one)*M_PI/beta;
      giom(in) = gf.Fourier(iom(in), tau);
    }
  }
  void MatsubaraFrequency(bl::Array<complex<double>,1>& giom, const bl::Array<double,1>& g, const bl::Array<int,1> iw, const string& statistics="fermi")const{
    Spline1D<double> gf(tau.size());
    for (int i=0; i<tau.size(); ++i) gf[i] = g(i);
    gf.splineIt(tau);
    giom.resize(iw.size());
    double one = (statistics=="fermi") ? 1.0 : 0.0;
    for (int in=0; in<iw.size(); in++){
      double iom = (2*iw(in)+one)*M_PI/beta;
      giom(in) = gf.Fourier(iom, tau);
    }
  }
  void MatsubaraFrequencyDebug(bl::Array<complex<double>,1>& giom, const bl::Array<double,1>& g, const bl::Array<int,1> iw, const string& statistics="fermi"){
    Spline1D<double> gf(tau.size());
    for (int i=0; i<tau.size(); ++i) gf[i] = g(i);
    gf.splineIt(tau);
    
    giom.resize(iw.size());
    double one = (statistics=="fermi") ? 1.0 : 0.0;

    ofstream dout("debug1.dat"); dout.precision(12);
    int Nw=500;
    dout << "# ";
    for (int in=0; in<iw.size(); in++){
      double iom = (2*iw(in)+one)*M_PI/beta;
      dout << iom << " ";
    }
    dout << " intg=" << gf.integrate() << " f(0)=" << gf.Fourier(0.0, tau) << " f(w0)=" << gf.Fourier(M_PI/beta, tau) << endl;
    for (int i=0; i<Nw; ++i){
      double t = i/(Nw-1.)*beta;
      intpar p = Interp( t, tau);
      dout << t << " " << gf(p) << endl;
    }
    for (int in=0; in<iw.size(); in++){
      double iom = (2*iw(in)+one)*M_PI/beta;
      giom(in) = gf.Fourier(iom, tau);
    }
  }
  template <class functor>
  void ExpandAnalyticFunction(functor fG, bl::Array<double,1>& gl, const string& statistics="fermi") const
  {
    const vector<Spline1D<double>>* pfU = (statistics=="fermi") ? &fU : &fU_bose;
    gl.resize(lmax);
    int M = static_cast<int>(pow(2.0,k_ovlp))+1;    // Number of points in Romberg integration routine
    bl::Array<double,1> utu(M); // storage for Romberg
    for (int l=0; l<lmax; l++){
      gl(l)=0;
      for (int i=0; i<tau.size()-1; i++){
	double a = tau(i);  // integral only between t_i and t_{i+1} points
	double b = tau(i+1);
	double fa = (*pfU)[l][i]   * fG(a, beta);
	double fb = (*pfU)[l][i+1] * fG(b, beta);
	utu(0)   = fa;
	utu(M-1) = fb;
	int ia=0;
	for (int j=1; j<M-1; j++){
	  double t = a + (b-a)*j/(M-1.0);
	  intpar p = InterpLeft( t, ia, tau);
	  utu(j) = (*pfU)[l](p) * fG(t, beta);
	}
	gl(l) += romberg(utu, b-a);
      }
    }
  }
  void SetUpFastInterp()
  {
    ff2.resize(tau.size(),lmax,2);
    dxi.resize(tau.size());
    for (int i=0; i<tau.size(); i++) dxi(i) = fU[0].dxi(i);
    for (int l=0; l<lmax; l++){
      for (int i=0; i<tau.size(); i++){
	ff2(i,l,0) = fU[l].f(i);
	ff2(i,l,1) = fU[l].f2(i);
      }
    }
  }
  void SetUpFastInterpBose()
  {
    lmax = max(lmax,lmax_bose);
    ff2.resize(tau.size(),lmax_bose,2);
    dxi.resize(tau.size());
    for (int i=0; i<tau.size(); i++) dxi(i) = fU_bose[0].dxi(i);
    for (int l=0; l<lmax; l++){
      for (int i=0; i<tau.size(); i++){
	ff2(i,l,0) = fU_bose[l].f(i);
	ff2(i,l,1) = fU_bose[l].f2(i);
      }
    }
  }
  void FastInterp(bl::Array<double,1>& res, const intpar& ip, double cst) const
  {
    int i= ip.i;
    double p = ip.p, q=1-ip.p;
    double dx26 = dxi(i)*dxi(i)/6.;
    double dq = dx26*q*(q*q-1);
    double dp = dx26*p*(p*p-1);
    q  *= cst;
    dq *= cst;
    p  *= cst;
    dp *= cst;
    double* _res = res.data();
    const double* __restrict__ _ff2 = &ff2(i,0,0);
    // this is fast equivalent of
    // res[l] += q * fU[l].f[i] + dq * fU[l].f2[i];
    for (int l=0; l<lmax; l++) _res[l] += q * _ff2[2*l] + dq * _ff2[2*l+1]; 
    _ff2 += 2*lmax;
    // this is fast equivalent of
    // res[l] += p * fu[l].f[i+1] + dp * fU[l].f2[i+1];
    for (int l=0; l<lmax; l++) _res[l] += p * _ff2[2*l] + dp * _ff2[2*l+1];
  }
  void FastInterp(double* _res, const intpar& ip, double cst) const
  {
    int i= ip.i;
    double p = ip.p, q=1-ip.p;
    double dx26 = dxi(i)*dxi(i)/6.;
    double dq = dx26*q*(q*q-1);
    double dp = dx26*p*(p*p-1);
    q  *= cst;
    dq *= cst;
    p  *= cst;
    dp *= cst;
    const double*  _ff2 = &ff2(i,0,0);
    // this is fast equivalent of
    // res[l] += q * fU[l].f[i] + dq * fU[l].f2[i];
    //cout<<"*ff2="<<_ff2<<" _res="<<_res<<" end_res="<<(_res+lmax)<<" lmax="<<lmax<<endl;
    for (int l=0; l<lmax; l++) _res[l] += q * _ff2[2*l] + dq * _ff2[2*l+1]; 
    _ff2 += 2*lmax;
    // this is fast equivalent of
    // res[l] += p * fu[l].f[i+1] + dp * fU[l].f2[i+1];
    for (int l=0; l<lmax; l++) _res[l] += p * _ff2[2*l] + dp * _ff2[2*l+1];
  }
  void FastInterp(bl::Array<double,1>& res, double t) const
  {
    intpar ip = Interp( t, tau );
    int i= ip.i;
    double p = ip.p, q=1-ip.p;
    double dx26 = dxi(i)*dxi(i)/6.;
    double dq = dx26*q*(q*q-1);
    double dp = dx26*p*(p*p-1);
    double* __restrict__ _res = res.data();
    const double* __restrict__ _ff2 = &ff2(i,0,0);
    // this is fast equivalent of
    // res[l] += q * fU[l].f[i] + dq * fU[l].f2[i];
    for (int l=0; l<lmax; l++) _res[l] = q * _ff2[2*l] + dq * _ff2[2*l+1]; 
    _ff2 += 2*lmax;
    // this is fast equivalent of
    // res[l] += p * fu[l].f[i+1] + dp * fU[l].f2[i+1];
    for (int l=0; l<lmax; l++) _res[l] += p * _ff2[2*l] + dp * _ff2[2*l+1];
  }
  void SlowInterp(bl::Array<double,1>& res, double t) const
  {
    intpar ip = Interp( t, tau );
    for (int l=0; l<lmax; l++) res(l) = fU[l](ip);
  }
};

void Interpolate_SVD_OnNewMesh(Array<double,2>& Sml_new, const SVDFunc& svdf, const Array<double,1>& tau, const Array<double,2>& Sml, double beta)
{
  Sml_new.resize(Sml.extent(0),tau.size());
  Sml_new=0;
  Array<double,1> u_t(svdf.lmax);
  for (int it=0; it<tau.size(); it++){
    intpar p = Interp( tau(it),  svdf.tau);
    for (int l=0; l<svdf.lmax; l++) u_t(l) = svdf.fU[l](p);
    for (int ip=0; ip<Sml.extent(0); ip++){
      double vt=0;
      for (int l=0; l<svdf.lmax; l++) vt += Sml(ip,l) * u_t(l);
      Sml_new(ip,it) = vt;
    }
  }
}

/**/
double fGs(double tau, double beta){
  // On real axis, this corresponds to
  // A(w) = 1/2*(delta(w-x0)+delta(w+x0))
  // with x0=1.
  return -exp(-beta/2.)*cosh(beta/2.-tau);
}
/**/
/*
using namespace std;
int main(){
  double x0=0.005;
  double L=10.;
  int Nw = 500;
  //double beta=9.35887;
  double beta=10;
  int lmax = 40;
  int k_ovlp=5;
  int Ntau = static_cast<int>(5*beta+100);
  clog<<"beta="<<beta<<" Ntau="<<Ntau<<endl;
  
  SVDFunc svdf;
  
  svdf.cmp("both", beta, lmax, Ntau, x0, L, Nw, clog);
  //void cmp(const string& statistics, double beta_, int lmax_, int Ntau, double x0_, double L_, double Nw, std::ofstream& clog, int k_ovlp_=5){
  
  // Print SVD functions
  svdf.Print("uls.dat");
  svdf.Print("ulb.dat", "bose");
  
  // Now we check how well is the orthogonality obeyed after reorthogonalization
  clog<<"Final overlap is smaller than "<<svdf.CheckOverlap()<<endl;

  SVDFunc svdf2;
  svdf2.Read("uls.dat", "fermi");

  ofstream fo("cc.dat");
  bl::Array<double,1>& tau = svdf.tau;
  bl::Array<double,1> Gts(tau.size());
  for (int i=0; i<tau.size(); i++){
    Gts(i) = fGs(tau(i), beta);
    fo<<tau[i]<<" "<<Gts(i)<<" "<<svdf[0][i]<<" "<<svdf[2][i]<<endl;
  }
  fo.close();

  bl::Array<double,1> gl(lmax);
  svdf.ExpandAnalyticFunction(fGs, gl);

  
  bl::Array<double,1> gl2(lmax);
  svdf.SetUpFastInterp();
  {
    k_ovlp=5;
    gl2.resize(lmax);
    gl2 = 0.0;
    int M = pow(2,k_ovlp)+1;    // Number of points in Romberg integration routine
    bl::Array<double,2> utu_T(M,lmax); // storage for Romberg
    double beta = tau(tau.size()-1);
    for (int i=0; i<tau.size()-1; i++){
      utu_T=0;
      double a = tau(i);  // integral only between t_i and t_{i+1} points
      double b = tau(i+1);
      int ia=0;
      for (int j=0; j<M; j++){
	double t = a + (b-a)*j/(M-1.0);
	intpar p = InterpLeft( t, ia, tau );
	double cst;
	cst = fGs(t, beta);
	///// Slow equivalent ////
	//for (int l=0; l<svdf.lmax; l++) utu_T(j,l) += cst * svdf.fU[l](p);
	bl::Array<double,1> uu = utu_T(j,bl::Range::all());
	svdf.FastInterp(uu, p, cst );
      }
      bl::Array<double,1> utu(M);
      for (int l=0; l<lmax; l++){
	for (int j=0; j<M; j++) utu(j) = utu_T(j,l);
	gl2(l) += romberg(utu, b-a);
      }
    }
  }
  cout<<"Should be equal:"<<endl;
  for (int l=0; l<svdf.lmax; l++)
    cout<<l<<" "<<gl(l)<<" "<<gl2(l)<<" "<<gl(l)-gl2(l)<<endl;

  Spline1D<double> gf;
  svdf.CreateSplineFromCoeff(gf, gl);
  
  int nmax=5000;
  bl::Array<complex<double>,1> giom;
  svdf.MatsubaraFrequency(giom, gf, nmax);
  
  ofstream fig("giom.dat");
  fig.precision(16);
  for (int in=0; in<nmax; in++){
    double iom = (2*in+1)*M_PI/beta;
    double gi_exact = -iom/(iom*iom+1.);
    fig<<iom<<" "<< giom(in).real() <<" "<< giom(in).imag() << "  " << gi_exact<<endl;
  }

  
  ofstream fog("gl.dat");
  fog.precision(16);
  for (int l=0; l<svdf.lmax; l++)
    if (fabs(gl(l))>1e-10)
      fog<<l<<" "<<gl(l)<<endl;
  fog.close();

  ofstream foG("Gapp.dat");
  foG.precision(16);
  for (int i=0; i<tau.size(); i++){
    double dval=0;
    for (int l=0; l<svdf.lmax; l++) dval += svdf(l,i)*gl(l);
    foG<<tau(i)<<" "<<dval<<" "<<Gts(i)<<" "<<fabs(dval-Gts(i))<<endl;
  }
  foG.close();

  return 0;
}
*/
