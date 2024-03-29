// @Copyright 2007-2017 Kristjan Haule
#include <iostream>
#include <fstream>
#include <vector>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include "logging.h"
#include "smesh.h"
#include "sfunction.h"
#include "tanmesh.h"
#include "romberg.h"
#include <time.h>

#ifdef _MPI
#include <mpi.h>
#endif

using namespace std;

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
inline void create_log_mesh(function1D<int>& ind_om, int nom_all, int nom, int ntail_)
{
  //Creates logarithmic mesh on Matsubara axis
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
  for (int i=0; i< istart; i++) ind_om[i]=i;
  int ii=istart;
  for (int i=0; i < ntail; i++){
    int t = int(istart*exp(alpha*i)+0.5);
    if (t != ind_om[ii-1])
      ind_om[ii++] = t;
  }
  ind_om.resize(ii);
}



class SVDFunc{
public:
  vector<spline1D<double> > fU, fU_bose; // SVD functions in imaginary time 
  mesh1D tau;                             // imaginary time optimized mesh
  int lmax, lmax_bose;                   // cutoff for the singular values
  int k_ovlp;                             // 2^{k_ovlp}+1 is the number of integration points in romberg method
  mesh1D om;                              // real axis mesh for SVD
  function1D<double> S, S_bose;          // eigenvalues of SVD decomposition
  function2D<double> Vt, Vt_bose;        // real axis basis functions
  function1D<double> dxi;                 // parts of fU splines stored in more efficient way for fast interpolation
  function3D<double> ff2;                 // fU splines stored in different way for fast interpolation
  double beta, L, x0;                     // parameters of the mesh should be remembered for reproducibility
  //
  SVDFunc() : lmax(0)  {  };
  //
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
  spline1D<double>& operator[](int l){
    return fU[l];
  }

  inline vector<int> get_lll(int ii)
  { // conversion from a single combined index to the three l indices.
    //int ii = l0 + lmax*l1 + lmax*lmax*lb;
    int lb_ = ii / (lmax*lmax);
    int ir = ii % (lmax*lmax);
    int l1_ = ir / lmax;
    int l0_ = ir % lmax;
    int tmp[] = { l0_, l1_, lb_ };
    return vector<int>(tmp, tmp+3 );
  }
  void _cmp_(double beta, int& lmax, vector<spline1D<double> >& fU, function1D<double>& S, function2D<double>& Vt, function2D<double>& K, ostream& clog){
    clock_t t0 = clock();
    int Nl = min(om.size(),tau.size());
    S.resize(Nl);
    Vt.resize(om.size(),Nl);
    function2D<double> U(Nl,tau.size());
    
    int lwork;
    {
      int n = min(tau.size(),om.size());
      int mx = max(tau.size(),om.size());
      lwork =  4*n*n + 8*max(mx, n*(n+1));
    }
    function1D<double> work(lwork);
    function1D<int> iwork(8*min(om.size(),tau.size()));
    int info=0;
    info = xgesdd(true, tau.size(), om.size(), K.MemPt(), K.fullsize_Nd(), S.MemPt(), U.MemPt(), U.fullsize_Nd(), Vt.MemPt(), Vt.fullsize_Nd(),work.MemPt(), lwork, work.MemPt(), iwork.MemPt());
    if (info!=0){
      clog << "svd is returning "<<info<<" in svdfunc"<<endl;
    }
    double dt = static_cast<double>( clock() - t0 )/CLOCKS_PER_SEC;
    clog<<"svd time="<<dt<<endl;

    for (int l=0; l<lmax; l++)
      if (fabs(S[l])<1e-13){
	lmax=l;
	clog<<"lmax reduced to "<<lmax<<endl;
	break;
      }
    clog<<"last singular value="<<S[lmax-1]<<endl;
    
    for (int i=0; i<tau.size(); i++)
      for (int l=0; l<Nl; l++)
	U(l,i) *= 1./sqrt(tau.Dh(i));

    // Splines functions, which are singular values of the Kernel.
    vector<spline1D<double> > fu(lmax);
    for (int l=0; l<lmax; l++){
      fu[l].resize(tau.size());
      for (int i=0; i<tau.size(); i++) fu[l][i] = U(l,i);
      int n = tau.size();
      double x1 = 0.5*(tau[1]+tau[0]);
      double df1 = (U(l,1)-U(l,0))/(tau[1]-tau[0]);
      double x2 = 0.5*(tau[2]+tau[1]);
      double df2 = (U(l,2)-U(l,1))/(tau[2]-tau[1]);
      double df0 = df1 + (df2-df1)*(0.0-x1)/(x2-x1);
      x1 = 0.5*(tau[n-1]+tau[n-2]);
      df1 = (U(l,n-1)-U(l,n-2))/(tau[n-1]-tau[n-2]);
      x2 = 0.5*(tau[n-2]+tau[n-3]);
      df2 = (U(l,n-2)-U(l,n-3))/(tau[n-2]-tau[n-3]);
      double dfn = df1 + (df2-df1)*(beta-x1)/(x2-x1);
      fu[l].splineIt(tau, df0, dfn);
    }
    
    // Calculates overlap between spline interpolations
    function2D<double> overlap(lmax,lmax);
    cmpOverlap(overlap, fu, tau, beta, k_ovlp);
    // Calculates eigensystem of the overlap
    function1D<double> oE(lmax);
    SymEigensystem(overlap, oE);
    function2D<double> sou(overlap);// sou = 1/sqrt(overlap)
    // Computes  1/sqrt(overlap)
    for (int i=0; i<lmax; i++){
      for (int j=0; j<lmax; j++){
	double dsum =0.0;
	for (int k=0; k<lmax; k++) dsum += overlap(k,i)*1/sqrt(oE[k])*overlap(k,j);
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
      int n = tau.size();
      double x1 = 0.5*(tau[1]+tau[0]);
      double df1 = (fU[l][1]-fU[l][0])/(tau[1]-tau[0]);
      double x2 = 0.5*(tau[2]+tau[1]);
      double df2 = (fU[l][2]-fU[l][1])/(tau[2]-tau[1]);
      double df0 = df1 + (df2-df1)*(0.0-x1)/(x2-x1);
      x1 = 0.5*(tau[n-1]+tau[n-2]);
      df1 = (fU[l][n-1]-fU[l][n-2])/(tau[n-1]-tau[n-2]);
      x2 = 0.5*(tau[n-2]+tau[n-3]);
      df2 = (fU[l][n-2]-fU[l][n-3])/(tau[n-2]-tau[n-3]);
      double dfn = df1 + (df2-df1)*(beta-x1)/(x2-x1);
      fU[l].splineIt(tau, df0, dfn);
    }
  }
  //
  void cmp(const string& statistics, double beta_, int lmax_, int Ntau, double x0_, double L_, double Nw, ostream& clog, int k_ovlp_=5){
    k_ovlp=k_ovlp_;
    beta = beta_;
    L = L_;
    x0 = x0_;
    
    GiveTanMesh(om, x0, L, Nw);
    GiveDoubleExpMesh(tau, beta, Ntau);

    if (statistics=="fermi" || statistics=="both"){
      lmax = lmax_;
      clock_t t0 = clock();
      function2D<double> K(om.size(),tau.size());
      for (int j=0; j<om.size(); j++){
	for (int i=0; i<tau.size(); i++){
	  K(j,i) = fermi_kernel(tau[i],om[j],beta);
	  K(j,i) *= om.Dh(j)*sqrt(tau.Dh(i));
	}
      }
      double dt = static_cast<double>( clock() - t0 )/CLOCKS_PER_SEC;
      clog<<"setup time="<<dt<<endl;
      
      _cmp_(beta, lmax, fU, S, Vt, K, clog);
    }
    
    if (statistics=="bose" || statistics=="both"){
      lmax_bose = lmax_;
      clock_t t0 = clock();
      function2D<double> K(om.size(),tau.size());
      for (int j=0; j<om.size(); j++){
	for (int i=0; i<tau.size(); i++){
	  K(j,i) = bose_kernel(tau[i],om[j],beta);
	  K(j,i) *= om.Dh(j)*sqrt(tau.Dh(i));
	}
      }
      double dt = static_cast<double>( clock() - t0 )/CLOCKS_PER_SEC;
      clog<<"setup time="<<dt<<endl;
      
      _cmp_(beta, lmax_bose, fU_bose, S_bose, Vt_bose, K, clog);
      //cout<<"S_bose="<<S_bose<<endl;
    }
  }
  
  void Print(const string& filename, const string& statistics="fermi"){
    ofstream pul(filename.c_str());
    pul<<"# "<<lmax<<" "<<tau.size()<<" "<<L<<" "<<x0<<"  "<<(om.size()/2)<<endl;
    pul.precision(16);
    for (int i=0; i<tau.size(); i++){
      pul<<tau[i]<<"  ";
      if (statistics=="fermi")
	for (int l=0; l<lmax; l++) pul<<fU[l][i]<<"  ";
      else
	for (int l=0; l<lmax_bose; l++) pul<<fU_bose[l][i]<<"  ";
      pul<<endl;
    }
  }

  double CheckOverlap(/*double beta*/){
    // Now we check how well is the orthogonality obeyed after reorthogonalization
    function2D<double> overlap(lmax,lmax);
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

  void cmpOverlap(function2D<double>& overlap, const vector<spline1D<double> >& fu, const mesh1D& tau, double beta, int k=5){
    int lmax = fu.size();

    int M = static_cast<int>(pow(2.0,k))+1;    // Number of points in Romberg integration routine
    vector<double> utu(M); // storage for Romberg
    
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
	    double a = tau[i];  // integral only between t_i and t_{i+1} points
	    double b = tau[i+1];
	    double fa = fu[l1][i]   * fu[l2][i];
	    double fb = fu[l1][i+1] * fu[l2][i+1];
	    utu[0]   = fa;
	    utu[M-1] = fb;
	    tint ia=0;
	    for (int j=1; j<M-1; j++){
	      intpar p = tau.InterpLeft( a + (b-a)*j/(M-1.0), ia);
	      utu[j] = fu[l1](p) * fu[l2](p);
	    }
	    oij += romberg(utu, b-a);
	  }
	  overlap(l1,l2) = oij;
	  overlap(l2,l1) = oij;
	}
      }
    }
  }

  template <class container>
  void CreateGfromCoeff(container& gf, const functionb<double>& gl, const string& statistics="fermi"){

    vector<spline1D<double> >* pfU = (statistics=="fermi") ? &fU : &fU_bose;
    
    gf.resize(tau.size());
    for (int i=0; i<tau.size(); i++) gf[i] = 0.0;
    
    for (int l=0; l<lmax; l++)
      for (int i=0; i<tau.size(); i++)
	gf[i] += gl[l]*(*pfU)[l][i];
  }

  void CreateSplineFromCoeff(spline1D<double>& gf, const functionb<double>& gl, const string& statistics="fermi"){
    CreateGfromCoeff(gf,gl, statistics);
    {
      double x1 = 0.5*(tau[1]+tau[0]);
      double df1 = (gf[1]-gf[0])/(tau[1]-tau[0]);
      double x2 = 0.5*(tau[2]+tau[1]);
      double df2 = (gf[2]-gf[1])/(tau[2]-tau[1]);
      double df0 = df1 + (df2-df1)*(0.0-x1)/(x2-x1);
      int n = tau.size();
      x1 = 0.5*(tau[n-1]+tau[n-2]);
      df1 = (gf[n-1]-gf[n-2])/(tau[n-1]-tau[n-2]);
      x2 = 0.5*(tau[n-2]+tau[n-3]);
      df2 = (gf[n-2]-gf[n-3])/(tau[n-2]-tau[n-3]);
      double dfn = df1 + (df2-df1)*(tau[n-1]-x1)/(x2-x1);
      gf.splineIt(tau, df0, dfn);
    }
  }
  
  template <class container>
  void MatsubaraFrequency(container& giom, const spline1D<double>& gf, int nmax, const string& statistics="fermi"){
    giom.resize(nmax);
    //double beta=tau[tau.size()-1];
    double one = (statistics=="fermi") ? 1.0 : 0.0;
    for (int in=0; in<nmax; in++){
      double iom = (2*in+one)*M_PI/beta;
      giom[in] = gf.Fourier_(iom, tau);
    }
  }

  template <class functor>
  void ExpandAnalyticFunction(const functor& fG, function1D<double>& gl, const string& statistics="fermi") const
  {
    const vector<spline1D<double> >* pfU = (statistics=="fermi") ? &fU : &fU_bose;
    
    gl.resize(lmax);
    int M = static_cast<int>(pow(2.0,k_ovlp))+1;    // Number of points in Romberg integration routine
    vector<double> utu(M); // storage for Romberg
    //double beta = tau[tau.size()-1];
    for (int l=0; l<lmax; l++){
      gl[l]=0;
      for (int i=0; i<tau.size()-1; i++){
	double a = tau[i];  // integral only between t_i and t_{i+1} points
	double b = tau[i+1];
	double fa = (*pfU)[l][i]   * fG(a, beta);
	double fb = (*pfU)[l][i+1] * fG(b, beta);
	utu[0]   = fa;
	utu[M-1] = fb;
	tint ia=0;
	for (int j=1; j<M-1; j++){
	  double t = a + (b-a)*j/(M-1.0);
	  intpar p = tau.InterpLeft( t, ia);
	  utu[j] = (*pfU)[l](p) * fG(t, beta);
	}
	gl[l] += romberg(utu, b-a);
      }
    }
  }

  int l_critical(const functionb<double>& gl, double max_ratio=3., const string& statistics="fermi"){
    int nsa=5;

    int _lmax_ = statistics=="fermi" ? lmax : lmax_bose;
    
    if (_lmax_<=nsa) return _lmax_;

    function1D<double> gos(_lmax_);

    function1D<double>* pS = (statistics=="fermi") ? &S : &S_bose;
    
    for (int i=0; i<_lmax_; i++) gos[i] = fabs(gl[i]/(*pS)[i]);

    //cout<<"gos="<<gos<<endl;
    
    double* p = gos.begin() + nsa+1; // pointer to the fifth entry
    
    for (int l=nsa+1; l<_lmax_; ++l && ++p){
      double sa = *max_element(p-nsa,p);  // maximum in the past five coefficients
      //cout<<l<<" g[l]="<<gl[l]<<" gos[l]="<<gos[l]<<" gmax="<<sa<<" g/gmax="<<fabs(gos[l]/sa)<<endl;
      if (fabs(gos[l]/sa)>max_ratio){
	//cout<<"lcritical should be "<<l<<endl;
	return l;
      }
    }
    return _lmax_;
  }

  void GiveRealAxis(function2D<double>& Aw, const functionb<double>& gl, double max_ratio=3., const string& statistics="fermi")
  {
    int lcritical = l_critical(gl,max_ratio,statistics);
    Aw.resize(lcritical,om.size());

    function1D<double> gos(gl.size());

    function1D<double>* pS = (statistics=="fermi") ? &S : &S_bose;
    
    for (int i=0; i<gl.size(); i++) gos[i] = gl[i]/(*pS)[i];

    function2D<double>* pVt = (statistics=="fermi") ? &Vt : &Vt_bose;
      
    for(int i=0; i<om.size(); i++)
      Aw(0,i) = -gos[0]*(*pVt)(i,0);
    
    for (int l=1; l<lcritical; l++){
      for (int i=0; i<om.size(); i++)
	Aw(l,i) = Aw(l-1,i) - gos[l]*(*pVt)(i,l);
    }

    if (statistics=="bose"){
      for (int l=0; l<lcritical; l++)
	for (int iw=0; iw<om.size(); iw++)
	  Aw(l,iw) *= M_PI*om[iw];
    }
  }

  void SetUpFastInterp()
  {
    ff2.resize(tau.size(),lmax,2);
    dxi.resize(tau.size());
    for (int i=0; i<tau.size(); i++) dxi[i] = fU[0].dxi[i];
    for (int l=0; l<lmax; l++){
      for (int i=0; i<tau.size(); i++){
	ff2(i,l,0) = fU[l][i];
	ff2(i,l,1) = fU[l].f2[i];
      }
    }
  }
  void FastInterp(functionb<double>& res, const intpar& ip, double cst) const
  {
    int i= ip.i;
    double p = ip.p, q=1-ip.p;
    double dx26 = dxi[i]*dxi[i]/6.;
    double dq = dx26*q*(q*q-1);
    double dp = dx26*p*(p*p-1);
    q  *= cst;
    dq *= cst;
    p  *= cst;
    dp *= cst;
    double* __restrict__ _res = res.MemPt();
    const double* __restrict__ _ff2 = ff2[i];
    // this is fast equivalent of
    // res[l] += q * fU[l].f[i] + dq * fU[l].f2[i];
    for (int l=0; l<lmax; l++) _res[l] += q * _ff2[2*l] + dq * _ff2[2*l+1]; 
    _ff2 += 2*lmax;
    // this is fast equivalent of
    // res[l] += p * fu[l].f[i+1] + dp * fU[l].f2[i+1];
    for (int l=0; l<lmax; l++) _res[l] += p * _ff2[2*l] + dp * _ff2[2*l+1];
  }
  void FastInterp_(functionb<double>& res, const intpar& ip, double cst) const
  {// This is the same as FastInterp, except res is zeroth first.
    int i= ip.i;
    double p = ip.p, q=1-ip.p;
    double dx26 = dxi[i]*dxi[i]/6.;
    double dq = dx26*q*(q*q-1);
    double dp = dx26*p*(p*p-1);
    q  *= cst;
    dq *= cst;
    p  *= cst;
    dp *= cst;
    double* __restrict__ _res = res.MemPt();
    const double* __restrict__ _ff2 = ff2[i];
    for (int l=0; l<lmax; l++) _res[l] = q * _ff2[2*l] + dq * _ff2[2*l+1]; 
    _ff2 += 2*lmax;
    for (int l=0; l<lmax; l++) _res[l] += p * _ff2[2*l] + dp * _ff2[2*l+1];
  }

#ifdef _MPI
  void ComputeClCoefficients(function2D<double>& Cl, int my_rank, int mpi_size, int Master, int nw=0, int ntail=0, int cutoff_iom=0){
#else  
  void ComputeClCoefficients(function2D<double>& Cl, int nw=0, int ntail=0, int cutoff_iom=0){
#endif    
    bool BRISI = false;
    
    int Nt = tau.size();
    //double L = om[om.size()-1];
    //double beta = tau[tau.size()-1];
    if (nw==0)    nw = L*beta/3;      // number of Matsubara points treated exactly
    if (ntail==0) ntail=150;          // number of points used for the tail
    if (cutoff_iom==0) cutoff_iom=100;// how far the tail will extend. We will go to Matsubara index nw*cutoff_iom
    
    function1D<int> ind_om;
    create_log_mesh(ind_om, nw*cutoff_iom, nw, ntail); // log mesh with nw exact points in [0,...nw-1], and ntail points distributed between [nw,nw*cutoff_iom]
    mesh1D iom(ind_om.size());
    for (int in=0; in<ind_om.size(); in++)  iom[in] = ind_om[in];
    iom.SetUp(0);

    if (BRISI){
      cout.precision(11);
      cout<<"# Starting single loop for u(iw) "<<endl;
    }
    // This creates u_l(iom) from u_l(tau)
    function2D<complex<double> > u_om(lmax,iom.size());
    for (int lf=0; lf<lmax; lf++){
      for (int in=0; in<iom.size(); in++){
	double om = (2*iom[in]+1.)*M_PI/beta;
	u_om(lf,in) = fU[lf].Fourier(om, tau);
      }
    }
    if (BRISI) cout<<"# Starting double loop for w(iw) "<<endl;
    // This creates w_l(iom) from u_l(tau).
    // w_{lb,lf}(iom) = 1/beta \sum_{iOm} u^b_{lb}(iOm) * u^f_{lf}(iom+iOm)
    // which is equivalent to
    // w_{lb,lf}(iom) = Integrate[ u^b_{lb}(beta-t) * u^f_{lf}(t) e^{i*om*t}, {t,0,beta}]
    int Nl2 = lmax_bose*lmax;
#ifdef _MPI
    int ipr_proc = (Nl2 % mpi_size==0) ? Nl2/mpi_size : Nl2/mpi_size+1;
    int iistart = ipr_proc*my_rank;
    int iiend   = ipr_proc*(my_rank+1);
    if (iistart>Nl2) iistart=Nl2;
    if (iiend  >Nl2) iiend=Nl2;
    int _Nl2_ = max(Nl2, mpi_size*ipr_proc);
    if (BRISI) cout<<"istart="<<iistart<<" iend="<<iiend<<" pr_proc="<<ipr_proc<<" pr_proc*mpi_size="<<ipr_proc*mpi_size<<" Nl2="<<Nl2<<endl;
#else
    int iistart = 0;
    int iiend = Nl2;
    int ipr_proc = Nl2;
    int _Nl2_ = Nl2;
#endif    
    function2D<complex<double> > w_iom(_Nl2_, iom.size());
    w_iom = 0.0;
    for (int ii=iistart; ii<iiend; ii++){
      // ii = lb * lmax + lf
      int lb = ii / lmax;
      int lf = ii % lmax;
      spline1D<double> wt(Nt);
      for (int it=0; it<tau.size(); it++)
	wt[it] = fU_bose[lb][Nt-1-it]*fU[lf][it]; // wt[t] = u^b_{lb}(beta-t) * u^f_{lf}(t)
      
      int n = tau.size();
      double x1 = 0.5*(tau[1]+tau[0]);
      double df1 = (wt[1]-wt[0])/(tau[1]-tau[0]);
      double x2 = 0.5*(tau[2]+tau[1]);
      double df2 = (wt[2]-wt[1])/(tau[2]-tau[1]);
      double df0 = df1 + (df2-df1)*(0.0-x1)/(x2-x1);
      x1 = 0.5*(tau[n-1]+tau[n-2]);
      df1 = (wt[n-1]-wt[n-2])/(tau[n-1]-tau[n-2]);
      x2 = 0.5*(tau[n-2]+tau[n-3]);
      df2 = (wt[n-2]-wt[n-3])/(tau[n-2]-tau[n-3]);
      double dfn = df1 + (df2-df1)*(beta-x1)/(x2-x1);
      wt.splineIt(tau, df0, dfn);
      // w(iom) = Fourier[ wt[t] ]
      for (int in=0; in<iom.size(); in++){
	double om = (2*iom[in]+1.)*M_PI/beta;
	//w_om(lb,lf,in) = wt.Fourier(om, tau);
	w_iom(ii,in) = wt.Fourier(om, tau);
      }
    }
#ifdef _MPI
    if (mpi_size>1)
      MPI_Allgather(MPI_IN_PLACE, ipr_proc*iom.size()*2, MPI_DOUBLE, w_iom.MemPt(), ipr_proc*iom.size()*2, MPI_DOUBLE, MPI_COMM_WORLD);
#endif
    if (BRISI) cout<<"# Starting last part"<<endl;
    int Nl = (lmax*lmax*lmax_bose);
    int offset=2;
    spline1D<double> ftmp(ntail+offset);
    mesh1D wom(ntail+offset);
    int ifirst = iom.size()-ntail-offset;
    for (int in=ifirst; in<iom.size(); in++) wom[in-ifirst] = iom[in];
    wom.SetUp(0);

    deque<pair<int,int> > icase;
    for (int i=0; i<Nl; i++){
      for (int j=i; j<Nl; j++){
	vector<int> li = get_lll(i);
	vector<int> lj = get_lll(j);
	if ( (li[0] + li[1] + li[2] + lj[0] + lj[1] + lj[2])%2 ) continue; // this vanishes due to symmetry
	icase.push_back( make_pair(i,j) );
      }
    }

#ifdef _MPI
    int pr_proc = (icase.size() % mpi_size==0) ? icase.size()/mpi_size : icase.size()/mpi_size+1;
    int istart = pr_proc*my_rank;
    int iend   = pr_proc*(my_rank+1);
    if (istart>icase.size()) istart=icase.size();
    if (iend  >icase.size()) iend=icase.size();
    if (BRISI) cout<<"istart="<<istart<<" iend="<<iend<<" pr_proc="<<pr_proc<<" pr_proc*mpi_size="<<pr_proc*mpi_size<<" icase.size="<<icase.size()<<endl;
#else
    int istart = 0;
    int iend = icase.size();
    int pr_proc = icase.size();
#endif    
    function1D<double> Cl0(pr_proc);
    Cl0=0;
    // This loop is parallelized
    for (int ii=istart; ii<iend; ii++){
      int i = icase[ii].first;
      int j = icase[ii].second;
      vector<int> li = get_lll(i);
      vector<int> lj = get_lll(j);
      if ( (li[0] + li[1] + li[2] + lj[0] + lj[1] + lj[2])%2 ) continue; // this vanishes due to symmetry

      // lj[2],lj[1],lj[0]  :  bose,fermi,fermi
      //                       lb * lmax + lf
      //complex<double> * w_1 = w_iom[ lj[2]*lmax + li[1] ].MemPt();
      //complex<double> * w_2 = w_iom[ li[2]*lmax + lj[1] ].MemPt();
      //complex<double> * u_1 = u_om[ li[0] ].MemPt();
      //complex<double> * u_2 = u_om[ lj[0] ].MemPt();
      complex<double> * w_1 = w_iom[ lj[2]*lmax + li[0] ].MemPt();
      complex<double> * w_2 = w_iom[ li[2]*lmax + lj[0] ].MemPt();
      complex<double> * u_1 = u_om[ li[1] ].MemPt();
      complex<double> * u_2 = u_om[ lj[1] ].MemPt();

      
      double wsum = 0.0;
      for (int in=0; in<iom.size()-ntail; in++)
	wsum += (w_1[in] * conj(w_2[in]) * u_1[in] * conj(u_2[in])).real();
      // fill int the 1D-interpolation for the tail
      for (int in=ifirst; in<iom.size(); in++)
	ftmp[in-ifirst] = (w_1[in] * conj(w_2[in]) * u_1[in] * conj(u_2[in])).real();
      // intepolate the tail
      int n=ftmp.size();
      ftmp.splineIt(wom, (ftmp[1]-ftmp[0])/(wom[1]-wom[0]), (ftmp[n-1]-ftmp[n-2])/(wom[n-1]-wom[n-2]) );
      // Now anlytically evaluate the sum over all points, but loop only over a few points in the tail.
      // The sum of the spline in the interval [n1,n2] can be analytically evaluated:
      //  1/2*(f[n2]+f[n1])*(n2-n1)  - 1/24*(n2-n1)*((n2-n1)^2-1) * (f2[n2]+f2[n1])
      //  where f2[n1], f2[n2] is the second derivative at n1 and n2 points.
      //  Note that the first and the last point come in with the weight 1/2, hence we need to correct for that.
      wsum += 0.5*ftmp[iom.size()-ntail-ifirst] + 0.5*ftmp[iom.size()-1-ifirst];
      for (int in=iom.size()-ntail; in<iom.size()-1; in++){
	double dh = iom[in+1]-iom[in];
	wsum += 0.5*(ftmp[in-ifirst]+ftmp[in+1-ifirst])*dh - dh/24.0*(dh*dh - 1.0)*(ftmp.f2[in+1-ifirst]+ftmp.f2[in-ifirst]);
      }
      /* //The above few lines are equivalent to
	 int iom_last = iom[iom.size()-1]; tint ia=0;
	 for (int iw=iom.size()-ntail; iw<iom_last; iw++){
	 double ww = ftmp( wom.InterpLeft(iw, ia) );
	 wsum += ww;
	 }
      */
      /*
      Cl(i,j) = 2*(wsum)/beta;
      Cl(j,i) = Cl(i,j);
      */
      Cl0[ii-istart] = 2*(wsum)/beta;
      //cout<<setw(4)<<i<<" "<<setw(4)<<j<<"     "<<setw(3)<<li[2]<<" "<<setw(3)<<li[1]<<" "<<setw(3)<<li[0]<<"; "<<setw(3)<<lj[2]<<" "<<setw(3)<<lj[1]<<" "<<setw(3)<<lj[0]<<"    "<<setw(10)<<left<<Cl(i,j)<<right<<endl;
    }
#ifdef _MPI
    if (mpi_size>1){
      if (BRISI){
	cout<<"MPI_Gather"<<endl;
	cout.flush();
      }
      function1D<double> cCl0;
      cCl0.resize(pr_proc*mpi_size);
      MPI_Gather(Cl0.MemPt(), pr_proc, MPI_DOUBLE, cCl0.MemPt(), pr_proc, MPI_DOUBLE, Master, MPI_COMM_WORLD);
      if (my_rank==Master){
	Cl0.resize(icase.size());
	for (int i=0; i<icase.size(); i++) Cl0[i] = cCl0[i];
      }
    }
    if (my_rank == Master){
#endif
      Cl.resize(Nl,Nl);
      Cl=0.0;
      for (int ii=0; ii<icase.size(); ii++){
	int i = icase[ii].first;
	int j = icase[ii].second;
	Cl(i,j) = Cl0[ii];
	Cl(j,i) = Cl0[ii];
      }
#ifdef _MPI    
    }
#endif  
  }
};

/*
double fGs(double tau, double beta){
  // On real axis, this corresponds to
  // A(w) = 1/2*(delta(w-x0)+delta(w+x0))
  // with x0=1.
  return -exp(-beta/2.)*cosh(beta/2.-tau);
}

double chis(double tau, double beta){
// On real axis this corresponds to
// chi''(x) = pi*[delta(w-x0)-delta(w+x0)]
// with x0=1/beta 
  const double bx2 = 0.5;
  return cosh(bx2*(1-2.*tau/beta))/sinh(bx2);
}

class fGspl{
public:
  spline1D<double>& fg;
  mesh1D t;
  fGspl(spline1D<double>& fg_, mesh1D& t_) : fg(fg_), t(t_) {}
  double operator()(double x, double beta) const{
    return fg(t.Interp(x));
  }
};
  
int main(){
  double x0=0.005;
  double L=10.;
  int Nw = 500;
  //double beta=9.35887;
  double beta=100;
  int lmax = 40;
  int k_ovlp=5;
  int Ntau = static_cast<int>(5*beta+100);
  bool readGFile=false;//true;
  
  spline1D<double> ginp;
  mesh1D tc;
  if (readGFile){
    ifstream  inp("../Gtau.dat_4");
    deque<double> ct, gt;
    while (inp){
      double t, g;
      inp >> t >> g;
      //cout<<t<<" "<<g<<endl;
      if (!inp) break;
      ct.push_back(t);
      gt.push_back(g);
      inp.ignore(numeric_limits<streamsize>::max(), '\n');
    }

    beta = ct.back()+ct[0];
    cout<<"beta="<<beta<<endl;
    ginp.resize(gt.size());
    tc.resize(ct.size());
    {
      for (int i=0; i<ct.size(); i++){
	tc[i] = ct[i];
	ginp[i] = gt[i];
      }
      tc.SetUp(0);
      double df0 = (ginp[1]-ginp[0])/(tc[1]-tc[0]);
      int n = tc.size();
      double dfn = (ginp[n-1]-ginp[n-2])/(tc[n-1]-tc[n-2]);
      cout<<"df0="<<df0<<" dfn="<<dfn<<endl;
      ginp.splineIt(tc, df0, dfn);
    }
  }
  fGspl Fg(ginp, tc);

  clog<<"beta="<<beta<<" Ntau="<<Ntau<<endl;
  
  SVDFunc svdf;
  
  svdf.cmp("both", beta, lmax, Ntau, x0, L, Nw, clog);

  // Print SVD functions
  svdf.Print("uls.dat",beta);
  svdf.Print("ulb.dat",beta, "bose");
  
  // Now we check how well is the orthogonality obeyed after reorthogonalization
  clog<<"Final overlap is smaller than "<<svdf.CheckOverlap(beta)<<endl;
  
  ofstream fo("cc.dat");
  mesh1D& tau = svdf.tau;
  function1D<double> Gts(tau.size());
  for (int i=0; i<tau.size(); i++){
    if (readGFile) Gts[i] = Fg(tau[i],beta);
    else Gts[i] = fGs(tau[i], beta);
    fo<<tau[i]<<" "<<Gts[i]<<" "<<svdf[0][i]<<" "<<svdf[2][i]<<endl;
  }
  fo.close();

  function1D<double> gl(lmax);
  if (readGFile) svdf.ExpandAnalyticFunction(Fg, gl);
  else svdf.ExpandAnalyticFunction(fGs, gl);


  function1D<double> gl2(lmax);
  svdf.SetUpFastInterp();
  {
    k_ovlp=5;
    gl2.resize(lmax);
    gl2 = 0.0;
    int M = pow(2,k_ovlp)+1;    // Number of points in Romberg integration routine
    function2D<double> utu_T(M,lmax); // storage for Romberg
    double beta = tau[tau.size()-1];
    for (int i=0; i<tau.size()-1; i++){
      utu_T=0;
      double a = tau[i];  // integral only between t_i and t_{i+1} points
      double b = tau[i+1];
      tint ia=0;
      for (int j=0; j<M; j++){
	double t = a + (b-a)*j/(M-1.0);
	intpar p = tau.InterpLeft( t, ia);
	double cst;
	if (readGFile) cst = Fg(t,beta);
	else cst = fGs(t, beta);
	///// Slow equivalent ////
	//for (int l=0; l<svdf.lmax; l++) utu_T(j,l) += cst * svdf.fU[l](p);
	svdf.FastInterp(utu_T[j], p, cst );
      }
      function1D<double> utu(M);
      for (int l=0; l<lmax; l++){
	for (int j=0; j<M; j++) utu[j] = utu_T(j,l);
	gl2[l] += romberg(utu, b-a);
      }
    }
  }
  cout<<"Should be equal:"<<endl;
  for (int l=0; l<lmax; l++)
    cout<<l<<" "<<gl[l]<<" "<<gl2[l]<<" "<<gl[l]-gl2[l]<<endl;
  
  spline1D<double> gf;
  svdf.CreateSplineFromCoeff(gf, gl);
  int nmax=5000;
  function1D<complex<double> > giom;
  svdf.MatsubaraFrequency(giom, gf, nmax);
  
  ofstream fig("giom.dat");
  fig.precision(16);
  for (int in=0; in<nmax; in++){
    double iom = (2*in+1)*M_PI/beta;
    double gi_exact = -iom/(iom*iom+1.);
    fig<<iom<<" "<<giom[in]<<" "<<gi_exact<<endl;
  }

  ofstream fog("gl.dat");
  fog.precision(16);
  for (int l=0; l<svdf.lmax; l++)
    if (fabs(gl[l])>1e-10)
      fog<<l<<" "<<gl[l]<<endl;
  fog.close();

  ofstream foG("Gapp.dat");
  foG.precision(16);
  for (int i=0; i<tau.size(); i++){
    double dval=0;
    for (int l=0; l<svdf.lmax; l++) dval += svdf[l][i]*gl[l];
    foG<<tau[i]<<" "<<dval<<" "<<Gts[i]<<" "<<fabs(dval-Gts[i])<<endl;
  }
  foG.close();


  //int lc = svdf.l_critical(gl);
  //cout<<lc<<endl;

  function2D<double> Aw;
  ofstream fAw("Aw.out");
  svdf.GiveRealAxis(Aw, gl, 2.0);
  for (int i=0; i<svdf.om.size(); i++){
    fAw << svdf.om[i]<<" ";
    for (int l=0; l<Aw.size_N(); l++){
      fAw << Aw(l,i) <<" ";
    }
    fAw << endl;
  }


  svdf.ExpandAnalyticFunction(chis, gl, "bose");
  svdf.CreateSplineFromCoeff(gf, gl, "bose");
  ofstream foc("Chid.dat");
  foc.precision(16);
  for (int i=0; i<svdf.tau.size(); i++){
    double val = chis(svdf.tau[i],beta);
    foc<<svdf.tau[i]<<" "<<gf[i]<<" "<<val<<" "<<fabs(val-gf[i])<<endl;
  }
  foc.close();
  
  svdf.MatsubaraFrequency(giom, gf, nmax, "bose");
  ofstream sig("siom.dat");
  sig.precision(16);
  for (int in=0; in<nmax; in++){
    double iom = (2*in+0.0)*M_PI/beta;
    double x0 = 1./beta;
    double gi_exact = 2*x0/(x0*x0+iom*iom);
    sig<<iom<<" "<<giom[in]<<" "<<gi_exact<<" "<<giom[in]-gi_exact<<endl;
  }


  cout<<"bose coefficients are"<<gl<<endl;
  
  svdf.GiveRealAxis(Aw, gl, 10.0, "bose");
  ofstream sAw("sw.out");
  for (int i=0; i<svdf.om.size(); i++){
    sAw << svdf.om[i]<<" ";
    for (int l=0; l<Aw.size_N(); l++) sAw << Aw(l,i) <<" ";
    sAw << endl;
  }

  svd.Print("usb.dat",common::beta,"bose");
  for (int l=0; l<svdf.lmax_bose; l++) svdf.MatsubaraFrequency(bUw[l], svdf.fU_bose[l], N_W);
  
  return 0;
}
*/
