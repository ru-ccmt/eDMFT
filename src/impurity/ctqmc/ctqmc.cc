// @Copyright 2007 Kristjan Haule
// 
//  WHAT WAS DONE in March 2017:
//  Sampling in SVD basis rather than frequency
//  Nmax is now not needed anymore, as the arrays are enlarged dynamically.
//  Acceptance step was speed up substantially by taking into account proportionality of the state_evolutions.
//  Implemented HB2 for general case
//  allows now for different modes "[S|G][H|M]" with S=self-energy sampling, G=Green's function sampling, H=Hubbard I tail, M=high frequency moment tails
//  Implemented Exchange_Two_Intervals for general case, i.e., new exchange move that reduces autocorrelation time
//  Allows for some percentage local moves in full
//  Prints two particle vertex in terms of (l1,l2,l3)
//  The S-mode for general case was finally implemented within svd scheme.
//
// Current problems: - It seems Exchange move for non-segment case has some problems with detail-balance. It is switched-off by default. You should have PMove=1 for non-segment case.
/// WHAT SHOULD BE DONE
//
//  Implement generalized (off-diagonal) susceptibility
//  Add two kinks like Patrick
//  Frequency dependent interaction
//
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <set>
#include <vector>
#include <deque>
#include <list>
#include <map>
#include <algorithm>
#include <list>
#include <fstream>
#include "logging.h"
#include "smesh.h"
#include "sfunction.h"
#include "random.h"
#include "number.h"
#include "mpi.h"
#include "common.h"
#include "intervals.h"
#include "svdfunc.h"
#include "matrixm.h"
#include "state.h"
#include "inout.h"
#include "bcast.h"
#include "operators.h"
#include "local.h"
#include "stateim.h"
#include "hb1.h"
#include "bathp.h"
#include "segment.h"

using namespace std;

typedef list<pair<double,double> > seg_typ;
ostream* xout, *yout;

template <class Rand>
class CTQMC{
  Rand& rand;                     // random number generator
  ClusterData& cluster;           // information about the atomic states of the particular problem
  BathProb BathP;          // What is the probability to create a kink in any of the baths
  vector<nIntervals> intervals;   // Stores information to compute the bath determinant and green's function: times and types corresponding to each bath sorted in the proper time order
  NOperators Operators;           // Stores information to compute local trace: all times and types in the proper time order.
public:
  function1D<NState> npraStates;  // list of all atomic states
private:
  function2D<NState> state_evolution_left, state_evolution_right; // evolution of atomic states for fast update of atomic trace
  function1D<NState> state_evolution_copy; // for fast update, we sometimes need to store previous copy
  vector<function2D<double> > MD; // matrix of the inverse of delta which is needed to compute local Green's function
  vector<MatrixM> TMD;            // class which manipulates matrix of deltas and local green's function
public:  
  vector<function2D<spline1D<double>* > > Deltat; // Delta in imaginary time but in matrix form (needed for non-diagonal cases)
  mesh1D tau;                     // imaginary time mesh
private:  
  Number matrix_element;          // current matrix element due to local part of the system
  function1D<Number> Trace;       // current atomic trace (for computing average)
  vector<double> detD;            // bath determinants
  
  mesh1D iom;                     // small imaginary frequency mesh on which local Green's function is manipulated
  mesh1D biom;                    // small bosonic imaginary frequency mesh for susceptibilities
  SVDFunc& svd;                    // SVD functions for storing the Green's function
    
  function1D<double> histogram;   // histogram
public:
  int Naver;                      // Number of measurements
  vector<function2D<dcomplex> > Gaver;     // Current approximation for the local Green's function
  function2D<dcomplex> Faver;     // Two particle F function, which can be used to obtain Self-energy
  
  vector<function2D<double> > Gsvd;        // Current approximation for the local Green's function in SVD mode
  vector<function2D<double> > Fsvd;        // Two particle F function, which can be used to obtain Self-energy
  function2D<double> Fsvd1;
  vector<function2D<double> > Gsvd_s;
  vector<function2D<double> > Gsvd_e;
  vector<function3D<double> > Gsvd_es;
private:
  function1D<double> observables; // some zero time quantities needed for calculation of occupancies and susceptibilities
  function1D<double> aver_observables; // the same zero-time quantities but the average in simulation rather than the current state
  function2D<double> Prob;        // Probability for each atomic state in current configuration
public:
  function2D<double> AProb;       // Average Probability
private:
  function1D<double> kaver;       // average perturbation order for each bath - measures kinetic energy
  function2D<double> P_transition;// Transition probability from state to state
public:  
  function2D<double> AP_transition;// Avergae Transition probability from state to state
  function2D<dcomplex> asuscg;       // sum/average of generalized susceptibility
private:
  int nsusc;
  function2D<double> susc;           // current frequency dependent susceptibility
  function2D<double> aver_susc;      // sampled frequency dependent susceptibility
  function2D<dcomplex> dsusc;        // temporary variable for frequency dependent susceptibility
  function2D<dcomplex> dsusc1;       // temporary variable for frequency dependent susceptibility
  function2D<dcomplex> suscg;        // generalized susceptibility
  vector<function2D<double> > Gtau;  // Green's function in imaginary time (not so precise)
  int NGta;                          // number of measurements of Gtau
  function2D<dcomplex> dexp;         // differences e^{iom*tau_{l+1}}-e^{iom*tau_l} needed for frequency dependent response functions
  //function2D<dcomplex> bexp;       // temporary for e^{iom*tau}
  int successful;                    // number of successful moves
  function1D<int> successfullC, successfullM;
  
  function1D<double> dtau;            // temporary for computeaverage
  function2D<double> tMD;             // just temporary used in global flip
  //function2D<double> Cf;              // C(Nmax,Nmax), mD(Nmax,Nmax)  : temporary space for GlobalFlip
  nIntervals tinterval;               // just temporary used in global flip

  int nomv, nOm;
  vector<function2D<dcomplex> > Mv; // M(om1,om2) matrix in imaginary frequency used to compute two particle vertex
public:
  function5D<dcomplex> VertexH, VertexF;
  function5D<double> VH;
  function2D<double> tMsMe;
  int NGtv;                        // number of measurements of vertex
  int gsign;
  double asign, asign_fast;
  double ratio_minM, ratio_minD;
  vector<deque<double> > Variation;
  //function1D<double> Njc;
  //vector<function2D<double> > Njcs;
  vector<function1D<double> > Njcs1;
  deque<pair<int,int> > g_exchange;
  static const bool Vertex_subtract = true;
private:
  Timer t_trial1, t_trial2, t_trial3, t_accept;
  Timer t_nj, t_fgf, t_vh, t_f, t_g;
  Timer t_mv, t_add, t_rm;
  //Timer t_vh1, t_vh2, t_vh3, t_vh4, t_vh5;
public:
  CTQMC(Rand& rand_, ClusterData& cluster_, BathProb& BathP_, int Nmax, const mesh1D& iom_large, const function2D<dcomplex>& Deltaw, SVDFunc& svd_,
	const vector<pair<double,double> >& ah, int nom, int nomv, int nomb, int nOm, int Ntau, double minM, double minD, bool IamMaster);
  ~CTQMC();
  
  double sample(long long max_steps); // main function which does the Monte Carlo sampling
 
  void RecomputeCix(); // tried to find better atomic base for the next iteration
  void ComputeFinalAverage(function1D<double>& mom, double& nf, double& TrSigma_G, double& Epot, double& Ekin); // From atomic probabilities computes various zero time quantities
  
  bool SaveStatus(int my_rank);     // saves all kinks to disc for shorter warmup in the next step
  bool RetrieveStatus(int my_rank, long long istep); // retrieves all kinks from disc
  void Resize(int newNmax);
  
  const mesh1D& ioms() const { return iom;}
  const mesh1D& iomb() const { return biom;}
  function1D<double>& Histogram() {return histogram;}
  function2D<double>& AverageProbability() {return AProb;}
  function1D<double>& Observables() {return aver_observables;}
  function1D<double>& k_aver() {return kaver;}
  vector<function2D<dcomplex> >& G_aver() { return Gaver;}
  const vector<function2D<double> >& gtau() const { return Gtau;}
  function2D<double>& Susc() {return aver_susc;}
  void EvaluateVertex(long long i);
  void WriteVertex(const function2D<dcomplex>& Gd, bool print_vertex_xxx);
  void WriteVertexSVD(const function2D<dcomplex>& Gd, const function2D<double>& Gd_, int my_rank, int mpi_size, int Master, bool print_vertex_xxx);
  void ReleaseSomeTempMemory();
private:
  void Add_One_Kink(long long istep, int ifl);
  void Remove_One_Kink(long long istep, int ifl);
  void Move_A_Kink(long long istep, int ifl);
  void Add_Two_Kinks(long long istep, int ifl);
  void Remove_Two_Kinks(long long istep, int ifl);
  void ComputeAverage(long long istep);
  void CleanUpdate(long long istep);
  void StoreCurrentState(long long istep);
  void StoreCurrentStateFast(function1D<double>& aver, long long istep);
  void PrintInfo(long long i, const function1D<double>& aver);
  void GlobalFlip(long long istep);
  void GlobalFlipFull(long long istep);
  void Segment_Add_One_Kink(long long istep, int ifl);
  void Segment_Remove_One_Kink(long long istep, int ifl);
  void Segment_Move_A_Kink(long long istep, int ifl);
  void Segment_Exchange_Two_Intervals(long long istep);
  void Exchange_Two_Intervals(long long istep);
  void Exchange_Two_Intervals_old(long long istep);
  void write_times(ostream& out);
  void CheckTimes(int istep);

public:
  enum AlgorithmType_t {OLD, WERNER, NEW};
  AlgorithmType_t AlgorithmType;
};

template <class Rand>
CTQMC<Rand>::CTQMC(Rand& rand_, ClusterData& cluster_, BathProb& BathP_, int Nmax, const mesh1D& iom_large, const function2D<dcomplex>& Deltaw, SVDFunc& svd_,
		   const vector<pair<double,double> >& ah, int nom, int nomv_, int nomb, int nOm_, int Ntau, double minM, double minD, bool IamMaster) :
  rand(rand_), cluster(cluster_), BathP(BathP_), intervals(common::N_ifl), Operators(Nmax, cluster),
  npraStates(cluster.nsize), state_evolution_left(cluster.nsize,Nmax), state_evolution_right(cluster.nsize,Nmax),state_evolution_copy(Nmax),
  MD(common::N_ifl), TMD(common::N_ifl), Deltat(cluster.N_ifl), Trace(cluster.nsize), detD(common::N_ifl),
  svd(svd_), histogram(Nmax),
  Naver(0), observables(4), aver_observables(4),
  Prob(cluster.nsize,common::max_size), AProb(cluster.nsize,common::max_size),
  kaver(common::N_ifl), nsusc(2), NGta(0), 
  successful(0), dtau(Nmax+1),
  nomv(nomv_), nOm(nOm_),
  NGtv(0), gsign(1), asign(0), asign_fast(0), ratio_minM(minM), ratio_minD(minD), successfullC(common::N_ifl), successfullM(common::N_ifl)
{

  //cerr<<"Zatom="<<cluster.Zatom<<" log(Zatom)="<<cluster.Zatom.exp_dbl()<<endl;
  //(*xout)<<"log(Zatom)="<<cluster.Zatom.exp_dbl()<<endl;
  
  // all atomic states are generated
  for (int j=0; j<npraStates.size(); j++) npraStates[j].SetPraState(j,cluster);

  successfullC=0; successfullM=0;
  histogram=0;

  {
    vector<spline1D<double> > Delta(cluster.N_unique_fl);// Delta in imaginary time
    //vector<function1D<double> > Delta(cluster.N_unique_fl);// Delta in imaginary time
    DeltaFourier(Ntau, common::beta, tau, Delta, iom_large, Deltaw, ah);
  
    int ibath=0;
    int ntau=tau.size();
    for (int ifl=0; ifl<cluster.N_ifl; ifl++){
      Deltat[ifl].resize(cluster.ifl_dim[ifl],cluster.ifl_dim[ifl]);
      for (int i1=0; i1<cluster.ifl_dim[ifl]; i1++){
	for (int i2=0; i2<cluster.ifl_dim[ifl]; i2++){
	  Deltat[ifl][i1][i2] = new spline1D<double>(Ntau);
	  int ib = cluster.tfl_index[ifl][i1][i2];
	  int fl = cluster.bfl_index[ifl][ib];
	  if (i1==i2){
	    // This is new change! Checks for diagonal Deltas to be causal! changed 18.10.2008
	    // Delta(tau) should be concave. Hence, if it becomes very small at two points, it
	    // should remain zero in all points between the two points. This is very important in
	    // insulators, because Delta(tau) can overshoot to positive values multiple times
	    // and kinks can be trapped in the range between the two crossing points, where Delta(tau)
	    // is causal, but should be zero.
	    const int Nx=100; // how many points between each to taut points to add
	    int first=tau.size();
	    for (int it=0; it<tau.size()-1; it++){
	      for (int ix=0; ix<Nx; ix++){
		//double x = tau[it]+(tau[it+1]-tau[it])*ix/(Nx-1.0);
		if ( Delta[fl](intpar(it, ix/(Nx-1.0)))>-common::minDeltat){
		  first=it;
		  goto break_me_first;
		}
	      }
	    }
	    break_me_first:
	    int last=-1;
	    for (int it=tau.size()-1; it>0; it--){
	      for (int ix=0; ix<Nx; ix++){
		//double x = tau[it]+(tau[it-1]-tau[it])*ix/(Nx-1.0);
		if ( Delta[fl](intpar(it-1, 1.0-ix/(Nx-1.0)))>-common::minDeltat){
		  last=it;
		  goto break_me_last;
		}
	      }
	    }
	    break_me_last:
	    bool FoundPositive = (first<=last);
	    if (FoundPositive){
	      for (int it=first; it<=last; it++) Delta[fl][it] = -common::minDeltat; // changed 18.10.2008
	      // we should avoid ringing of the spine, hence make the drop gradual.
	      Delta[fl][first] = (Delta[fl][first-1]+Delta[fl][first])/3.;
	      if (first+1<last) Delta[fl][first+1] = (Delta[fl][first-1]+Delta[fl][first])/21.;
	      if (first+2<last) Delta[fl][first+2] = (Delta[fl][first-1]+Delta[fl][first])/150.;
	      Delta[fl][last] = (Delta[fl][last]+Delta[fl][last+1])/3.;
	      if (last-1>first) Delta[fl][last-1] = (Delta[fl][last]+Delta[fl][last+1])/21.;
	      if (last-2>first) Delta[fl][last-2] = (Delta[fl][last]+Delta[fl][last+1])/150.;
	    }
	    for (int it=0; it<tau.size(); it++){ // check ones more
	      if (Delta[fl][it]>-common::minDeltat){
		FoundPositive=true;
		Delta[fl][it] = -common::minDeltat;
	      }
	    }
	    //if (FoundPositive) cout<<"Found it "<<ifl<<" "<<FoundPositive<<endl;
	    // Here you should probably smoothen
	  }
	  for (int it=0; it<tau.size(); it++){
	    // Finally set the sline array.
	    double Dt = Delta[fl][it];
	    if (cluster.conjg[ifl][ib]) Dt = -Delta[fl][ntau-it-1];
	    Dt *= cluster.sign[ifl][ib];
	    (*Deltat[ifl][i1][i2])[it] = Dt;
	  }
	  int n = tau.size();
	  /*
	    double df0 = ((*Deltat[ifl][i1][i2])[1]-(*Deltat[ifl][i1][i2])[0])/(tau[1]-tau[0]);
	    double dfn = ((*Deltat[ifl][i1][i2])[n-1]-(*Deltat[ifl][i1][i2])[n-2])/(tau[n-1]-tau[n-2]);
	    (*Deltat[ifl][i1][i2]).splineIt(tau,df0,dfn);
	  */
	  spline1D<double>& Dc =  *Deltat[ifl](i1,i2);
	  double x1 = 0.5*(tau[1]+tau[0]);
	  double df1 = (Dc[1]-Dc[0])/(tau[1]-tau[0]); // derivative in the first midpoint
	  double x2 = 0.5*(tau[2]+tau[1]);
	  double df2 = (Dc[2]-Dc[1])/(tau[2]-tau[1]); // derivative in the second midpoint
	  double df0 = df1 + (df2-df1)*(0.0-x1)/(x2-x1); // extrapolated derivative at 0
	  x1 = 0.5*(tau[n-1]+tau[n-2]);
	  df1 = (Dc[n-1]-Dc[n-2])/(tau[n-1]-tau[n-2]); // derivative at the last mindpoint
	  x2 = 0.5*(tau[n-2]+tau[n-3]);
	  df2 = (Dc[n-2]-Dc[n-3])/(tau[n-2]-tau[n-3]); // derivative at the point before the last point
	  double dfn = df1 + (df2-df1)*(common::beta-x1)/(x2-x1);  // extrapolated derivative
	  Dc.splineIt(tau,df0,dfn);

	  ibath++;
	}
      }
    }
  }
  
  if (IamMaster & common::fastFilesystem){
    for (int ifl=0; ifl<cluster.N_ifl; ifl++){
      for (int i1=0; i1<cluster.ifl_dim[ifl]; i1++){
	for (int i2=0; i2<cluster.ifl_dim[ifl]; i2++){
	  ofstream out(NameOfFile_("Delta.tau",ifl,i1*2+i2).c_str());
	  const int Ncheck_points=10000;
	  for (int it=0; it<Ncheck_points; it++){
	    double t = it*common::beta/(Ncheck_points-1.);
	    double Dt = (*Deltat[ifl][i1][i2])(tau.Interp(t));
	    out<<t<<" "<<Dt<<endl;
	  }
	  ofstream out2(NameOfFile_("rDelta.tau",ifl,i1*2+i2).c_str());
	  for (int it=0; it<tau.size(); it++){
	    out2<<tau[it]<<" "<<(*Deltat[ifl][i1][i2])[it]<<endl;
	  }
	}
      }
    }
  }
  if (!common::Qsvd){
    // Matsubara Frequency meshes
    iom.resize(nom);
    for (int i=0; i<iom.size(); i++) iom[i] = iom_large[i];
    biom.resize(nomb);
    for (int i=0; i<nomb; i++) biom[i] = 2*i*M_PI/common::beta;
  }
  
  // propagators are set-up
  for (size_t ifl=0; ifl<intervals.size(); ifl++) intervals[ifl].SetUp(Nmax,iom,biom,common::beta,cluster.ifl_dim[ifl]);
  
  tinterval.SetUp(Nmax,iom,biom,common::beta,0);
  // matrix M is set up
  for (size_t ifl=0; ifl<MD.size(); ifl++) MD[ifl].resize(Nmax,Nmax);
  for (size_t ifl=0; ifl<TMD.size(); ifl++) TMD[ifl].SetUp(Nmax, tau, Deltat[ifl], cluster.tfl_index[ifl], iom, common::beta, cluster.v2fl_index[ifl], nomv, svd.lmax);

  // set-up local part of the system
  StartStateEvolution(state_evolution_left, state_evolution_right, npraStates, Operators, cluster, Trace, Prob);

  if (! common::Qsvd){
    // approximation for Glocal
    Gaver.resize(cluster.N_ifl);
    for (int ifl=0; ifl<cluster.N_ifl; ifl++){
      Gaver[ifl].resize(cluster.N_baths(ifl),iom.size());
      Gaver[ifl] = 0;
    }  
    if (common::QHB2){
      //Njc.resize(cluster.Nvfl);
      Faver.resize(cluster.Nvfl, iom.size());
      Faver=0.0;
    }
  }else{
    Gsvd.resize(cluster.N_ifl);
    for (int ifl=0; ifl<cluster.N_ifl; ifl++){
      Gsvd[ifl].resize(cluster.N_baths(ifl),svd.lmax);
      Gsvd[ifl] = 0;
    }
    if (common::QHB2){
      /*
      if (common::Segment){
	Njcs.resize(cluster.N_ifl);
	Fsvd.resize(cluster.Nvfl);
	for (int fl=0; fl<cluster.Nvfl; fl++){
	  Fsvd[fl].resize(cluster.Nvfl,svd.lmax);
	  Fsvd[fl] = 0.0;
	}
      }else{
      */
	Njcs1.resize(cluster.N_ifl);
	Fsvd1.resize(cluster.Nvfl,svd.lmax);
	Fsvd1=0;
	/*}*/
    }
  }
  
  for (size_t ifl=0; ifl<detD.size(); ifl++) detD[ifl] = 1;
  
  for (int j=0; j<cluster.nsize; j++){
    Prob[j].resize(cluster.msize(j+1));
    AProb[j].resize(cluster.msize(j+1));
  }
  AProb = 0.0;
  
  Gtau.resize(cluster.N_ifl);
  if (common::SampleGtau>0){
    for (size_t ifl=0; ifl<Gtau.size(); ifl++){
      Gtau[ifl].resize(cluster.N_baths(ifl),Ntau);
      Gtau[ifl]=0;
    }
  }

  kaver = 0;

  if (common::SampleSusc && !common::Qsvd){
    susc.resize(nsusc,nomb);
    aver_susc.resize(nsusc,nomb);
    dsusc.resize(nomb,nsusc);
    dsusc1.resize(nsusc*cluster.nsize, nomb);
    dexp.resize(Nmax+1,nom);
    aver_susc = 0;
    suscg.resize(cluster.DOsize,biom.size());
    asuscg.resize(cluster.DOsize,biom.size());
    asuscg = 0;
  }else{
    susc.resize(nsusc,1);
  }
  

  if (common::cmp_vertex){
    if (common::Qsvd){
      tMsMe.resize(svd.lmax, svd.lmax);
      VH.resize(cluster.Nvfl,cluster.Nvfl,svd.lmax_bose,svd.lmax,svd.lmax);
      VH=0.0;
      Gsvd_s.resize(cluster.Nvfl);
      Gsvd_e.resize(cluster.Nvfl);
      Gsvd_es.resize(cluster.Nvfl);
    }else{
      VertexH.resize(cluster.Nvfl,cluster.Nvfl,2*nOm-1,2*nomv,2*nomv);
      VertexF.resize(cluster.Nvfl,cluster.Nvfl,2*nOm-1,2*nomv,2*nomv);
      VertexH=dcomplex(0,0);
      VertexF=dcomplex(0,0);
      Mv.resize(cluster.Nvfl);
      for (unsigned ifl=0; ifl<Mv.size(); ifl++) Mv[ifl].resize(2*nomv,2*nomv);//(-nomv,nomv,-nomv,nomv);
    }
  }

  // Segment : 1 , Algorithm=NEW,    flips only global flip
  //           2 , Algorithm=NEW,    flips all equivalent baths from cix file
  //           3 , Algorithm=WERNER, flips only global flip
  //           4 , Algorithm=WERNER, flips all equivalent baths from cix file
  //           5 , Algorithm=OLD,    flips only global flip
  //           6 , Algorithm=OLD,    flips all equivalent baths from cix file
  if ((common::Segment+1)/2 == 1)
    AlgorithmType = NEW;
  else if ((common::Segment+1)/2 == 2)
    AlgorithmType = WERNER;
  else if ((common::Segment+1)/2 == 3)
    AlgorithmType = OLD;
  if (common::Segment){
    for (int iu=0; iu<cluster.gflip_ifl.size(); iu++){
      int ifa = cluster.gflip_ifl[iu].first;
      int ifb = cluster.gflip_ifl[iu].second;
      int dim=cluster.ifl_dim[ifa];
      for (int b1=0; b1<dim; b1++){
	if (common::Segment % 2 == 1){
	  int fla = cluster.fl_from_ifl[ifa][b1];
	  int flb = cluster.fl_from_ifl[ifb][b1];
	  g_exchange.push_back(make_pair(fla,flb));
	  g_exchange.push_back(make_pair(flb,fla));
	}else{
	  for (int b2=0; b2<dim; b2++){
	    if (cluster.bfl_index[ifa][b1+dim*b1]==cluster.bfl_index[ifb][b2+dim*b2]){
	      int fla = cluster.fl_from_ifl[ifa][b1];
	      int flb = cluster.fl_from_ifl[ifb][b2];
	      g_exchange.push_back(make_pair(fla,flb));
	      g_exchange.push_back(make_pair(flb,fla));
	    }
	  }
	}
      }
    }
  }else{
    for (int iu=0; iu<cluster.gflip_ifl.size(); iu++){
      int ifa = cluster.gflip_ifl[iu].first;
      int ifb = cluster.gflip_ifl[iu].second;
      int dim=cluster.ifl_dim[ifa];
      for (int b1=0; b1<dim; b1++){
	int fla = cluster.fl_from_ifl[ifa][b1];
	int flb = cluster.fl_from_ifl[ifb][b1];
	g_exchange.push_back(make_pair(fla,flb));
	g_exchange.push_back(make_pair(flb,fla));
      }
    }
    /*
    for (int i=0; i<g_exchange.size(); i++){
      cout<<"pair "<<g_exchange[i].first<<" "<<g_exchange[i].second<<endl;
    }
    */
  }
  
  if (common::SampleTransitionP){
    P_transition.resize(cluster.nsize,2*cluster.N_flavors);
    AP_transition.resize(cluster.nsize,2*cluster.N_flavors);
    AP_transition = 0.0;
  }

  //Cf.resize(Nmax,Nmax);//  temporary space for GlobalFlip
  tMD.resize(Nmax,Nmax);//  temporary space for GlobalFlip

}

template <class Rand>
CTQMC<Rand>::~CTQMC()
{
  for (size_t ifl=0; ifl<Deltat.size(); ifl++)
    for (int i1=0; i1<Deltat[ifl].size_N(); i1++)
      for (int i2=0; i2<Deltat[ifl].size_Nd(); i2++)
	delete Deltat[ifl][i1][i2];
}

template <class Rand>
void CTQMC<Rand>::ReleaseSomeTempMemory()
{// This routine releases some memory at the end, so that the high frequency
  // calculation would not run out of memory.
  state_evolution_left.ReleaseMemory();
  state_evolution_right.ReleaseMemory();
  state_evolution_copy.ReleaseMemory();
  for (size_t ifl=0; ifl<MD.size(); ifl++) MD[ifl].ReleaseMemory();

  if (common::Qsvd){
    if (common::QHB2){
      Fsvd1.ReleaseMemory();
      for (int ifl=0; ifl<cluster.N_ifl; ifl++) Njcs1[ifl].ReleaseMemory();
    }
  }
  
  if (common::cmp_vertex){
    if (common::Qsvd){
      tMsMe.ReleaseMemory();
      for (int ifl=0; ifl<cluster.Nvfl; ifl++){
	Gsvd_s[ifl].ReleaseMemory();
	Gsvd_e[ifl].ReleaseMemory();
	//Gsvd_es[ifl].~function3D();
      }
    }else{
      //VertexH.resize(cluster.Nvfl,cluster.Nvfl,2*nOm-1,2*nomv,2*nomv);
      //VertexF.resize(cluster.Nvfl,cluster.Nvfl,2*nOm-1,2*nomv,2*nomv);
      for (int ifl=0; ifl<cluster.Nvfl; ifl++) Mv[ifl].ReleaseMemory();
    }
  }
  tMD.ReleaseMemory();
  dexp.ReleaseMemory();
  dtau.ReleaseMemory();

  //Operators.ReleaseSomeTempMemory();
  //for (size_t ifl=0; ifl<intervals.size(); ifl++) intervals[ifl].ReleaseSomeTempMemory();
  for (size_t ifl=0; ifl<TMD.size(); ifl++) TMD[ifl].ReleaseSomeTempMemory();
}


template <class Rand>
void CTQMC<Rand>::Resize(int newNmax)
{
  /*
  if (common::my_rank==1 || true){
    cout<<"Before Resize left="<<endl;
    cout<<state_evolution_left.size_N()<<" "<<state_evolution_left.size_Nd()<<" "<<state_evolution_left.fullsize_N()<<" "<<state_evolution_left.fullsize_Nd()<<endl;
    for (int i=0; i<cluster.size(); i++){
      cout<<"i="<<i<<" size="<<state_evolution_left[i].size()<<" : ";
      for (int j=0; j<state_evolution_left[i].size(); j++){
	cout<<state_evolution_left(i,j).istate<<" ";
      }
      cout<<endl;
    }
    cout<<"right="<<endl;
    cout<<state_evolution_right.size_N()<<" "<<state_evolution_right.size_Nd()<<" "<<state_evolution_right.fullsize_N()<<" "<<state_evolution_right.fullsize_Nd()<<endl;
    for (int i=0; i<cluster.size(); i++){
      cout<<"i="<<i<<" size="<<state_evolution_right[i].size()<<" : ";
      for (int j=0; j<state_evolution_right[i].size(); j++){
	cout<<state_evolution_right(i,j).istate<<" ";
      }
      cout<<endl;
    }
  }
  */
  int Nmax = histogram.size();
  Operators.Resize(newNmax);
  for (size_t ifl=0; ifl<intervals.size(); ifl++) intervals[ifl].Resize(newNmax);
  tinterval.Resize(newNmax);

  for (size_t ifl=0; ifl<MD.size(); ifl++) MD[ifl].resize_preserve(newNmax,newNmax);
  for (size_t ifl=0; ifl<TMD.size(); ifl++) TMD[ifl].Resize(newNmax);
  
  state_evolution_left.resize_preserve(cluster.nsize,newNmax);
  state_evolution_right.resize_preserve(cluster.nsize,newNmax);
  state_evolution_copy.resize_preserve(newNmax);
  // simple acumulators
  histogram.resize_preserve(newNmax);
  histogram.resize_virtual(newNmax);
  for (int i=Nmax; i<newNmax; i++) histogram[i]=0;
  
  dtau.resize_preserve(newNmax+1);
  // temporary storage
  if (common::SampleSusc && !common::Qsvd) dexp.resize(newNmax+1,iom.size());
  tMD.resize(newNmax,newNmax);
  (*yout)<<"Nmax resized from "<<Nmax<<" to "<<newNmax<<endl;

  /*
  if (common::my_rank==1 || true){
    cout<<"After Resize left="<<endl;
    cout<<state_evolution_left.size_N()<<" "<<state_evolution_left.size_Nd()<<" "<<state_evolution_left.fullsize_N()<<" "<<state_evolution_left.fullsize_Nd()<<endl;
    for (int i=0; i<cluster.size(); i++){
      cout<<"i="<<i<<" size="<<state_evolution_left[i].size()<<" : ";
      for (int j=0; j<state_evolution_left[i].size(); j++){
	cout<<state_evolution_left(i,j).istate<<" ";
      }
      cout<<endl;
    }
    cout<<"right="<<endl;
    cout<<state_evolution_right.size_N()<<" "<<state_evolution_right.size_Nd()<<" "<<state_evolution_right.fullsize_N()<<" "<<state_evolution_right.fullsize_Nd()<<endl;
    for (int i=0; i<cluster.size(); i++){
      cout<<"i="<<i<<" size="<<state_evolution_right[i].size()<<" : ";
      for (int j=0; j<state_evolution_right[i].size(); j++){
	cout<<state_evolution_right(i,j).istate<<" ";
      }
      cout<<endl;
    }
  }
  */
}

template <class Rand>
void CTQMC<Rand>::Add_One_Kink(long long istep, int ifl)
{
  t_add.start();
  if (Operators.full()){
    int newNmax = static_cast<int>(Operators.max_size()*1.2+20);
    Resize(newNmax); // always resize for at least 20 percent + 20 more entries
    //return;
  }

  int dim = cluster.ifl_dim[ifl];
  int kn = intervals[ifl].size()/2;

  double t_start, t_end;
  int bfls, bfle;
  double P_to_add;
  
  if (1-drand48() > common::PlocalMoves){ // non-local moves
    bfls = static_cast<int>(dim*rand());
    bfle = static_cast<int>(dim*rand());    
    t_start = common::beta*rand();
    t_end = common::beta*rand();
    P_to_add = sqr(common::beta*dim)/sqr(kn+1.);
  }else{ // local moves
    int bfl  = static_cast<int>(dim*rand());  // type of creation/anhilation operators we attempt to add
    bfls = bfle = bfl;   // here we are always adding operators of the same type
    int knb = intervals[ifl].Nbtyp_s[bfl];    // how many operators of this type are there?
    int flse = cluster.fl_from_ifl[ifl][bfl]; // combined index for bath (ifl,bfl)
    if (knb==0){
      t_start = common::beta*rand();         // new random time in the existsing interval
      t_end = common::beta*rand();         // new random time in the existsing interval
      P_to_add = sqr(common::beta);
    }else{  
      double r1=rand();
      double r2=rand();
      if (r1==0 || r1==1 || r2==0 || r2==1) return;   //  We do now want two times to be exactly the same. The probability should vanish anyway.
      int type = static_cast<int>(rand()*2);          // either c or c^+
      int iop_1 = static_cast<int>(knb*rand());       // which operator from ordered list
      int i1 = iop_1;                                 // index in the multidimensional bath might be different in general
      if (dim>1) i1 = intervals[ifl].FindWhichOneIs(type, bfl, iop_1); // this is the index in multidimensional interval
      double t1 = intervals[ifl].time(type,i1);       // the beginning of the interval
      pair<int,int> i_iop = intervals[ifl].NextTimeIndex_(1-type, t1, bfl); // index of the next operator of the opposite type
      int i2 = i_iop.first;                           // index in the multidimensional bath index
      int iop_2 = i_iop.second;                       // index when counting only this bfl type
      double t2 = intervals[ifl].time(1-type,i2);     // the end of the interval
      double seg_length = t2-t1;                      // length of the segment
      if (seg_length<0) seg_length += common::beta;
      double t1_new = t1 + seg_length*r1;             // new random time in the existsing interval
      double t2_new = t1 + seg_length*r2;             // new random time in the existsing interval
      if (type==nIntervals::c){
	t_start = min(t1_new,t2_new);
	t_end = max(t1_new,t2_new);
      }else{
	t_end = min(t1_new,t2_new);
	t_start = max(t1_new,t2_new);
      }
      if (t_start>common::beta) t_start -= common::beta;
      if (t_end>common::beta) t_end -= common::beta;
      P_to_add = (sqr(seg_length)*knb)/(knb+1.)/2.;
    }
  }
  int fls = cluster.fl_from_ifl[ifl][bfls];
  int fle = cluster.fl_from_ifl[ifl][bfle];

  t_trial1.start();
  bool ssc1;

  ssc1 = Operators.Try_Add_Cd_C_(fls, t_start, fle, t_end, state_evolution_left, state_evolution_right);
  t_trial1.stop();
  t_add.stop();
  
  if (!ssc1 || t_start==t_end) return;
  
  t_trial2.start();
  t_add.start();
  
  pair<int,int> opt = Operators.Try_Add_Cd_C(fls, t_start, fle, t_end);

  int op_i = min(opt.first, opt.second), op_f = max(opt.first, opt.second);
  double ratioD = TMD[ifl].AddDetRatio(MD[ifl], kn, t_start, bfls, t_end, bfle, intervals[ifl]);
  t_trial2.stop();
  
  // First flip the coint, and than see if we accept the step.
  double Prand = 1-rand();    // 1-rand() because don't ever want to accept move with zero probability
  double P_part = P_to_add * fabs(ratioD); // all of acceptance prob. except trace ratio
  double P_r = Prand/P_part;

  t_trial3.start();
  double ms=1;
  bool Accept;
  if (common::LazyTrace){
    list< pair<int, double> > exp_sums;
    ComputeExpTrace(exp_sums, state_evolution_left, state_evolution_right, op_i, op_f, Operators, npraStates, cluster, 2);
    Accept = LazyComputeTrace(P_r, matrix_element, ms, exp_sums, state_evolution_left, state_evolution_right, op_i, op_f, Operators, npraStates, cluster, 2);
  }else{
    Number matrix_element_new = ComputeTrace(state_evolution_left, state_evolution_right, op_i, op_f, Operators, npraStates, cluster, 2);
    ms = fabs(divide(matrix_element_new, matrix_element));
    Accept = P_r<ms;
  }
  t_trial3.stop();
  
  if (ms>ratio_minM && Accept){
    t_accept.start();
    int is, ie;
    intervals[ifl].Find_is_ie(t_start, is, t_end, ie); // SIGN: HERE WE HAVE is,ie
    ratioD = TMD[ifl].AddUpdateMatrix(MD[ifl], kn, t_start, is, bfls, t_end, ie, bfle, intervals[ifl],istep);
    pair<int,int> ipl = intervals[ifl].InsertExponents(t_start, is, bfls, t_end, ie, bfle);
    if (!common::Qsvd) TMD[ifl].AddUpdate_Gf(intervals[ifl], Mv,istep);
    
    pair<int,int> opt = Operators.Add_Cd_C(fls, t_start, fle, t_end,   IntervalIndex(ifl, nIntervals::cd, ipl.first), IntervalIndex(ifl, nIntervals::c, ipl.second));
    int op_i_ = min(opt.first, opt.second), op_f_ = max(opt.first, opt.second); // Changed March 2013

    // And updates state evolution for fast computation of Trace
    Number new_matrix_element;
    new_matrix_element = UpdateStateEvolution(state_evolution_left, state_evolution_right, state_evolution_copy, op_i_, op_f_, Operators, npraStates, cluster, Trace, Prob, 2, istep); // Changed March 2013
    
    // Operators.sign() :  gives number of crossings + number of forward propagators (sign for the matrix element to bring it to original form)
    // sign_isie: gives permutation of columns and rows of the determinant to bring it to original form
    //
    int sign_isie = 1-2*((is+ie)%2);
    int from_orig = Operators.sign()*sign_isie;
    double ms1 = divide(new_matrix_element, matrix_element);
    double Ps = ratioD*ms1*from_orig;
    gsign *= sigP<0>(Ps);
    detD[ifl] *= ratioD;
    matrix_element = new_matrix_element;
    successful++;
    successfullC[ifl]++;
    ComputeAverage(istep);
    t_accept.stop();
  }
  t_add.stop();
}


template <class Rand>
void CTQMC<Rand>::Remove_One_Kink(long long istep, int ifl)
{
  int kn = intervals[ifl].size()/2;
  if (kn==0) return;
  int dim = cluster.ifl_dim[ifl];

  int ie, is, bfle, bfls;
  int iop_s, iop_e;
  double P_to_remove;
  
  if (1-drand48() > common::PlocalMoves){
    ie = static_cast<int>(kn*rand());
    is = static_cast<int>(kn*rand());
    bfle = intervals[ifl].btype_e(ie);
    bfls = intervals[ifl].btype_s(is);
    iop_s = intervals[ifl].FindSuccessive(0, is, bfls);
    iop_e = intervals[ifl].FindSuccessive(1, ie, bfle);
    P_to_remove = sqr(kn)/sqr(common::beta*dim);
  }else{
    int bfl = static_cast<int>(dim*rand());   // type of the multidimensional bath we want to remove
    bfle = bfls = bfl;
    int knb = intervals[ifl].Nbtyp_s[bfl];    // how many segments of this type are there?
    if (knb==0) return;
    int flse = cluster.fl_from_ifl[ifl][bfl];
    if (knb==1){
      P_to_remove = 1./sqr(common::beta);
      iop_s=0;
      is = iop_s;                            // index in the multidimensional bath might be different in general
      if (dim>1) is = intervals[ifl].FindWhichOneIs(nIntervals::cd, bfl, iop_s); // this is the index in multidimensional interval
      iop_e=0;
      ie = iop_e;
      if (dim>1) ie = intervals[ifl].FindWhichOneIs(nIntervals::c, bfl, iop_e); // this is the index in multidimensional interval
    }else{
      int type = static_cast<int>(rand()*2);          // either c or c^+
      int iop_1 = static_cast<int>(knb*rand());       // which operator from ordered list
      int i1 = iop_1;                                 // index in the multidimensional bath might be different in general
      if (dim>1) i1 = intervals[ifl].FindWhichOneIs(type, bfl, iop_1); // this is the index in multidimensional interval
      double t1 = intervals[ifl].time(type,i1);       // beginning of the interval
      pair<int,int> i_iop = intervals[ifl].NextTimeIndex_(1-type, t1, bfl); // index of the next operator of the opposite type
      int iop_2 = i_iop.second;                       // index when counting only this bfl type
      int i2 = i_iop.first;                           // index in the multidimensional bath index
      double t2 = intervals[ifl].time(1-type,i2);       // beginning of the interval
      double t1_previous = intervals[ifl].PreviousTime(1-type, i2, true);
      double t2_next = intervals[ifl].NextTime(type, i1, true);
      double seg_length = t2_next-t1_previous;         // length of the interval
      if (t2<t1) seg_length -= common::beta;           // jumped around twice when calculating t1_previous and t2_next
      if (type==nIntervals::cd){
	is=i1;  // if we started with cd, is=i1
	iop_s=iop_1;
	ie=i2;
	iop_e=iop_2;
      }else{
	ie=i1;
	iop_e=iop_1;
	is=i2;
	iop_s=iop_2;
      }
      P_to_remove = 2.*knb/(knb-1.)/sqr(seg_length);
    }
  }

  t_rm.start();
  t_trial1.start();

  int fle = cluster.fl_from_ifl[ifl][bfle];
  int fls = cluster.fl_from_ifl[ifl][bfls];
  int ipe, ips;
  bool ssc1;
  ssc1 = Operators.Try_Remove_C_Cd_(fle, iop_e, fls, iop_s, ipe, ips, state_evolution_left, state_evolution_right, istep);
  
  t_trial1.stop();
  t_rm.stop();
  
  if (!ssc1) return;

  t_rm.start();
  t_trial2.start();    
  int op_i = min(ipe, ips);
  int op_f = ipe>ips ? ipe-2 : ips-1;
  Operators.Try_Remove_C_Cd(ipe, ips);
  double ratioD = TMD[ifl].RemoveDetRatio(MD[ifl], kn, is, ie);
  t_trial2.stop();
  
  // First flip the coint, and than see if we accept the step.
  double Prand = 1-rand();  // 1-rand() because don't ever want to accept move with zero probability
  double P_part = fabs(ratioD)*P_to_remove;
  double P_r = Prand/P_part;

  t_trial3.start();
  double ms=1;
  bool Accept;

  if (common::LazyTrace){
    list< pair<int, double> > exp_sums;
    ComputeExpTrace(exp_sums, state_evolution_left, state_evolution_right, op_i, op_f, Operators, npraStates, cluster, -2);
    Accept = LazyComputeTrace(P_r, matrix_element, ms, exp_sums, state_evolution_left, state_evolution_right, op_i, op_f, Operators, npraStates, cluster, -2);
  }else{
    Number matrix_element_new = ComputeTrace(state_evolution_left, state_evolution_right, op_i, op_f, Operators, npraStates, cluster, -2);
    ms = fabs(divide(matrix_element_new,matrix_element));
    Accept = P_r<fabs(ms);
  }	
  t_trial3.stop();

  if (ms>ratio_minM && Accept){
    t_accept.start();
    Operators.Remove_C_Cd(ipe, ips);
    // And updates state evolution for fast computation of Trace
    Number new_matrix_element;
    new_matrix_element = UpdateStateEvolution(state_evolution_left, state_evolution_right, state_evolution_copy, op_i, op_f, Operators, npraStates, cluster, Trace, Prob, -2, istep);
    if (!common::Qsvd) TMD[ifl].RemoveUpdate_Gf(MD[ifl], kn, is, ie, intervals[ifl], Mv);
    int sign_isie = 1-2*((is+ie)%2);
    int from_orig = Operators.sign()*sign_isie;
    double ms1 = divide(new_matrix_element, matrix_element);
    double Ps = ratioD*ms1*from_orig;
    gsign *=sigP<0>(Ps);
    double ratioDn = TMD[ifl].RemoveUpdateMatrix(MD[ifl], kn, is, ie, intervals[ifl]);
    intervals[ifl].RemoveExponents(is,ie);
    detD[ifl] *= ratioDn;
    matrix_element = new_matrix_element;
    successful++;
    successfullC[ifl]++;

    ComputeAverage(istep);
    t_accept.stop();

  }
  t_rm.stop();
}

template <class Rand>
void CTQMC<Rand>::Move_A_Kink(long long istep, int ifl)
{
  int kn = intervals[ifl].size()/2;
  if (kn==0) return;
  t_mv.start();
  int dim = cluster.ifl_dim[ifl];
  int type = static_cast<int>(rand()*2);
  int to_move = static_cast<int>(kn*rand());
  int bfl = intervals[ifl].btype(type, to_move);
  int fl = cluster.fl_from_ifl[ifl][bfl];
  int opera = 2*fl+type;
  int to_move_o = intervals[ifl].FindSuccessive(type, to_move, bfl);

  double t_old = intervals[ifl].time(type,to_move);
  double t_prev = intervals[ifl].PreviousTime(type,to_move);
  double t_next = intervals[ifl].NextTime(type, to_move);

  double t_new;
  if (t_prev<t_old && t_next>t_old)
    t_new = t_prev + (t_next-t_prev)*rand();
  else if (t_prev>t_old)
    t_new = t_prev-common::beta + (t_next-t_prev+common::beta)*rand(); //t_prev = 0;//-= common::beta; SHOULD BE CORRECTED
  else if (t_next<t_old)
    t_new = t_prev + (t_next+common::beta-t_prev)*rand();//common::beta;//+= common::beta; SHOULD BE CORRECTED
  else if (fabs(t_next-t_prev)<1e-10){
    t_new = common::beta*rand();
  } else{ // have only one time
    cerr<<"DID NOT EXPECT TO HAPPEN"<<endl;
    (*yout)<<"t_old="<<t_old<<" t_next="<<t_next<<" t_prev="<<t_prev<<endl;
    t_new = t_old;
  }

  int iws=0;
  while (t_new<0) {t_new += common::beta; iws++;}
  while (t_new>common::beta) {t_new -= common::beta; iws++;}

  t_trial1.start();
  bool ssc1;
  ssc1 = Operators.Try_Move_(opera, t_old, to_move_o, t_new, state_evolution_left, state_evolution_right);
  t_trial1.stop();
  t_mv.stop();
  
  if (!ssc1) return;

  t_mv.start();
  t_trial2.start();    
  double ratioD;
  if (type==0) ratioD = TMD[ifl].Move_start_DetRatio(MD[ifl], kn, t_old, t_new, bfl, to_move, intervals[ifl]);
  else ratioD = TMD[ifl].Move_end_DetRatio(MD[ifl], kn, t_old, t_new, bfl, to_move, intervals[ifl]);
  int ip_old, ipo, ip_new;
  Operators.TryMove(opera, t_old, to_move_o, ip_old, ipo, t_new, ip_new);
  int op_i = min(ipo, ip_new), op_f = max(ipo, ip_new);
  t_trial2.stop();    

  double Prand = 1-rand();
  double P_part = fabs(ratioD);
  double P_r = Prand/P_part;

  t_trial3.start();
  double ms=1;
  bool Accept;
  if (common::LazyTrace){
    list< pair<int, double> > exp_sums;
    ComputeExpTrace(exp_sums, state_evolution_left, state_evolution_right, op_i, op_f, Operators, npraStates, cluster, 0);
    Accept = LazyComputeTrace(P_r, matrix_element, ms, exp_sums, state_evolution_left, state_evolution_right, op_i, op_f, Operators, npraStates, cluster, 0);
  }else{
    Number matrix_element_new = ComputeTrace(state_evolution_left, state_evolution_right, op_i, op_f, Operators, npraStates, cluster, 0);
    ms = fabs(divide(matrix_element_new, matrix_element));
    Accept = P_r<ms;
  }
  t_trial3.stop();
  
  if (ms>ratio_minM && Accept){
    t_accept.start();
    Operators.Move(ip_old, ip_new);

    Number new_matrix_element;
    new_matrix_element = UpdateStateEvolution(state_evolution_left, state_evolution_right, state_evolution_copy, op_i, op_f, Operators, npraStates, cluster, Trace, Prob, 0, istep);
    
    int i_new = intervals[ifl].FindIndex(t_new, t_old, type);
    if (type==0) TMD[ifl].Move_start_UpdateMatrix(MD[ifl], kn, t_old, t_new, bfl, to_move, i_new, intervals[ifl]);
    else TMD[ifl].Move_end_UpdateMatrix(MD[ifl], kn, t_old, t_new, bfl, to_move, i_new, intervals[ifl]);
    intervals[ifl].MoveExponents(type, t_old, to_move, t_new, i_new);
    if (!common::Qsvd) TMD[ifl].MoveUpdate_Gf(intervals[ifl], Mv, istep);
    ratioD *= 1-2*((i_new+to_move+iws)%2);
    int sign_isie = 1-2*((i_new+to_move+iws)%2);
    int from_orig = Operators.sign()*sign_isie;
    double ms1 = divide(new_matrix_element, matrix_element);
    double Ps = ratioD*ms1*from_orig;
    gsign *= sigP<0>(Ps);
    detD[ifl] *= ratioD;
    matrix_element = new_matrix_element;
    successful++;
    successfullM[ifl]++;
    ComputeAverage(istep);
    t_accept.stop();
  }
  t_mv.stop();
}


template <class Rand>
void CTQMC<Rand>::Exchange_Two_Intervals(long long istep)
{
  int ir = static_cast<int>(rand()*g_exchange.size());
  int fla = g_exchange[ir].first;
  int flb = g_exchange[ir].second;
  int ifa  = cluster.ifl_from_fl[fla];
  int bfla = cluster.bfl_from_fl[fla];
  int ifb  = cluster.ifl_from_fl[flb];
  int bflb = cluster.bfl_from_fl[flb];
  int kna = intervals[ifa].size()/2;
  int knb = intervals[ifb].size()/2;
  if (kna==0 || knb==0) return;
  int ka = intervals[ifa].Nbtyp_s[bfla];    // how many segments of first type are there?
  int kb = intervals[ifb].Nbtyp_s[bflb];    // how many segments of second type are there?
  if (ka==0 || kb==0) return;               // no segments, can't exchange...

  t_mv.start();
  
  int dim = cluster.ifl_dim[ifa];           // dimension of this bath
  int typea = static_cast<int>(rand()*2);   // either c or c+ for first segment
  int typeb = static_cast<int>(rand()*2);   // either c or c+ for second segment
  //int typeb = 1-typea;
  //int typeb=typea;
  int iop_a = static_cast<int>(ka*rand());       // which operator from ordered list
  int ia = iop_a;                                 // index in the multidimensional bath might be different in general
  if (dim>1) ia = intervals[ifa].FindWhichOneIs(typea, bfla, iop_a); // this is the index in multidimensional interval
  double t_a = intervals[ifa].time(typea,ia); // time in the first orbital

  pair<double,double> tv = intervals[ifa].Closest_Times(1-typea, bfla, t_a);
  pair<double,double> t_same;
  pair<int,int> ind_same;
  pair<double,double> t_oppo;
  pair<int,int> ind_oppo;
  intervals[ifb].Closest_Times(t_same, ind_same, 1-typeb, bflb, t_a);
  intervals[ifb].Closest_Times(t_oppo, ind_oppo, typeb, bflb, t_a);

  double t_b;
  int ib;
  int op_a = 2*fla+typea;
  int op_b = 2*flb+typeb;
  t_mv.stop();
  // Currently this move has problems with detail balance. I could not make it work for non-segment case.
  //#define _JUST_TRY
#ifdef _JUST_TRY
  vector<double> _times_(4);
  vector<int> _inds_(4);
  _times_[0] = t_same.first;
  _times_[1] = t_same.second;
  _times_[2] = t_oppo.first;
  _times_[3] = t_oppo.second;
  _inds_[0] = ind_same.first;
  _inds_[1] = ind_same.second;
  _inds_[2] = ind_oppo.first;
  _inds_[3] = ind_oppo.second;
  int im=0;
  for (int i=1; i<4; i++){
    if (fabs(_times_[i]-t_a) < fabs(_times_[im]-t_a)) im=i;
  }
  t_b = _times_[im];
  ib = _inds_[im];
  if (im<2){
    typeb = 1-typeb;
    op_b = 2*flb + typeb;
  }
#else
  if (typea==nIntervals::c){
    if (t_same.first<t_oppo.first){ //    [tv.1]-----------------[t]    [tv.2]----------------
      if (t_oppo.first>tv.first){   //    ----[ts.1]     [to.1]------------[ts.2]  [to.2]-----
	t_b = t_oppo.first;
	ib = ind_oppo.first;
      } else return;
    }else{
      if (t_oppo.second<tv.second){//  [tv.1]--------[t]            [tv.2]----------------
	t_b = t_oppo.second;       //[to.1]----[ts.1]       [to.2]------------[ts.2]
	ib = ind_oppo.second;
      } else return;
    }
  }else{
    if (t_same.second>t_oppo.second){ // ----[tv.1]     [t]-----------[tv.2]
      if (tv.second>t_oppo.second){   // [to.1] [ts.1]-----[to.2]  [ts.2]-----
	t_b = t_oppo.second;
	ib = ind_oppo.second;
      } else return;
    }else{
      if (tv.first<t_oppo.first){ //----[tv.1]        [t]---------------[tv.2]  
	t_b = t_oppo.first;	  //[ts.1]---[to.1]         [ts.2]-----[to.2]   -----
	ib = ind_oppo.first;
      }else return;
    }
  }
#endif
  
  int iws=0;
  if (t_b<0) {t_b += common::beta; iws=1;}
  if (t_b>=common::beta) {t_b -= common::beta; iws=1;}
  t_b = intervals[ifb].time(typeb,ib);
  if (t_a==t_b) return;

  t_mv.start();
  
  int ip_a, ip_b;
  if (t_a<t_b){
    ip_a = Operators.FindIndex(op_a, t_a, 0, istep);
    ip_b = Operators.FindIndex(op_b, t_b, ip_a, istep);
  }else{
    ip_b = Operators.FindIndex(op_b, t_b, 0, istep);
    ip_a = Operators.FindIndex(op_a, t_a, ip_b, istep);
  }

  Number matrix_element_new = ComputeTryalForExchange(ip_a, ip_b,  state_evolution_left,  state_evolution_right, Operators, npraStates, cluster, istep);
  double ms = fabs(divide(matrix_element_new, matrix_element));
  
  t_trial2.start();    
  double ratioDa, ratioDb;
  if (typea==nIntervals::c)
    ratioDa = TMD[ifa].Move_end_DetRatio(  MD[ifa], kna, t_a, t_b, bfla, ia, intervals[ifa]);
  else
    ratioDa = TMD[ifa].Move_start_DetRatio(MD[ifa], kna, t_a, t_b, bfla, ia, intervals[ifa]);
  if (typeb==nIntervals::c)
    ratioDb = TMD[ifb].Move_end_DetRatio(  MD[ifb], knb, t_b, t_a, bflb, ib, intervals[ifb]);
  else
    ratioDb = TMD[ifb].Move_start_DetRatio(MD[ifb], knb, t_b, t_a, bflb, ib, intervals[ifb]);
  
  t_trial2.stop();
  double P_to_accept = ms*fabs(ratioDa*ratioDb);
  double det_a = fabs(ratioDa);
  double det_b = fabs(ratioDb);
  bool Accept = 1-rand() < P_to_accept;
  
  if (Accept && ms>ratio_minM && det_a>ratio_minD && det_b>ratio_minD){
    t_accept.start();
    /*
    Number new_matrix_element = UpdateForExchange(ip_a, ip_b, state_evolution_left, state_evolution_right, Trace, Operators, npraStates, cluster, Prob, istep);
    if (fabs(divide(new_matrix_element,matrix_element_new)-1.)>1e-4){
      cerr<<"at istep "<<istep<<" ERROR Different matrix elements m_new="<<new_matrix_element.exp_dbl()<<" m_trial="<<matrix_element_new.exp_dbl()<<endl;
    }
    */
    Operators.ExchangeTwo(ip_a, ip_b);
    Number new_matrix_element = UpdateForExchange2(ip_a, ip_b, state_evolution_left, state_evolution_right, Trace, Operators, npraStates, cluster, Prob, istep);
    if (fabs(divide(new_matrix_element,matrix_element_new)-1.)>1e-4){
      cerr<<"at istep "<<istep<<" ERROR Different matrix elements m_new="<<new_matrix_element.exp_dbl()<<" m_trial="<<matrix_element_new.exp_dbl()<<endl;
    }
    
    int ia_new = intervals[ifa].FindIndex_nearby(t_b, typea, ia);
    int ib_new = intervals[ifb].FindIndex_nearby(t_a, typeb, ib);
    
    if (typea==nIntervals::c)
      TMD[ifa].Move_end_UpdateMatrix(  MD[ifa], kna, t_a, t_b, bfla, ia, ia_new, intervals[ifa]);
    else
      TMD[ifa].Move_start_UpdateMatrix(MD[ifa], kna, t_a, t_b, bfla, ia, ia_new, intervals[ifa]);
    if (typeb==nIntervals::c)
      TMD[ifb].Move_end_UpdateMatrix(  MD[ifb], knb, t_b, t_a, bflb, ib, ib_new, intervals[ifb]);
    else
      TMD[ifb].Move_start_UpdateMatrix(MD[ifb], knb, t_b, t_a, bflb, ib, ib_new, intervals[ifb]);
    IntervalsExchangeExponents(typea, intervals[ifa], ia, ia_new, typeb, intervals[ifb], ib, ib_new);
    if (!common::Qsvd) TMD[ifa].MoveUpdate_Gf(intervals[ifa], Mv, istep);
    if (!common::Qsvd) TMD[ifb].MoveUpdate_Gf(intervals[ifb], Mv, istep);
    int sign_a_isie = 1-2*((ia_new+ia+iws)%2);
    int sign_b_isie = 1-2*((ib_new+ib+iws)%2);
    ratioDa *= sign_a_isie;
    ratioDb *= sign_b_isie;
    int from_orig = -1*sign_a_isie*sign_b_isie;
    double ms1 = divide(new_matrix_element, matrix_element);
    double Ps = ratioDa*ratioDb*ms1*from_orig;
    gsign *= sigP<0>(Ps);
    detD[ifa] *= ratioDa;
    detD[ifb] *= ratioDb;
    matrix_element = new_matrix_element;
    successful++;
    successfullM[ifa]++;
    successfullM[ifb]++;
    ComputeAverage(istep);
    t_accept.stop();
  }
  t_mv.stop();
  return;
}


template <class Rand>
void CTQMC<Rand>::Segment_Add_One_Kink(long long istep, int ifl)
{
  if (Operators.full()){
    int newNmax = static_cast<int>(Operators.max_size()*1.2+20);
    //cout<<"Increasing Nmax to "<<newNmax<<" at step "<<istep<<endl;
    Resize(newNmax); // always resize for at least 20 percent + 20 more entries
    //return;
  }
  int dim = cluster.ifl_dim[ifl];           // dimensionality of this bath
  int kn = intervals[ifl].size()/2;         // total number of kinks/2 in the multidimensional bath
  int bfl  = static_cast<int>(dim*rand());  // type of creation/anhilation operators we attempt to add
  int knb = intervals[ifl].Nbtyp_s[bfl];    // how many operators of this type are there?
  int flse = cluster.fl_from_ifl[ifl][bfl]; // combined index for bath (ifl,bfl)
  double t_start, t_end;
  double P_to_add;

  t_trial1.start();
  if (AlgorithmType==OLD){
    t_start = common::beta*rand();
    t_end = common::beta*rand();
    if (!Segment_Try_Add_Cd_C_(bfl, t_start, bfl, t_end, intervals[ifl], istep)) return;
    P_to_add = sqr(common::beta/(knb+1.));
  }else if(AlgorithmType==WERNER) {
    double seg_length;
    double r = rand();
    if (r==0 || r==1) return; // We do now want two times to be exactly the same. The probability should vanish anyway.
    if (rand()<0.5){   // Adding a segment
      t_start = common::beta*rand();     // new creation operator
      double next_t_start, next_t_end, previous_t_start, previous_t_end;
      bool onseg = intervals[ifl].CheckOnSegmentHere(bfl,t_start,next_t_start,next_t_end,previous_t_start,previous_t_end,istep); // On segment, we can not add segment
      if (onseg || t_start==next_t_end) return; // On segment, we can not add segment
      if (t_start==previous_t_end) return; // The new and the old time are exactly the same
      seg_length = next_t_start-t_start;
      t_end = t_start + seg_length*r;
      if (t_end>common::beta) t_end -= common::beta;
    }else{             // Adding antisegment
      t_end = common::beta*rand();       // new anhilation operator
      double next_t_start, next_t_end, previous_t_start, previous_t_end;
      bool onseg = intervals[ifl].CheckOnSegmentHere(bfl,t_end,next_t_start,next_t_end,previous_t_start,previous_t_end,istep); // Can add only on segment
      if (!onseg) return; // Can add anti-segment only on segment
      seg_length = next_t_end-t_end;
      t_start = t_end + seg_length*r;
      if (t_start>common::beta) t_start -= common::beta;
    }
    P_to_add = (common::beta*seg_length)/(knb+1.);
  }else if(AlgorithmType==NEW){
    if (knb==0){
      t_start = common::beta*rand();         // new random time in the existsing interval
      t_end = common::beta*rand();         // new random time in the existsing interval
      P_to_add = sqr(common::beta);
    }else{
      double r1=rand();
      double r2=rand();
      if (r1==0 || r1==1 || r2==0 || r2==1) return;   //  We do now want two times to be exactly the same. The probability should vanish anyway.
      int type = static_cast<int>(rand()*2);          // either c or c^+
      int iop_1 = static_cast<int>(knb*rand());       // which operator from ordered list
      int i1 = iop_1;                                 // index in the multidimensional bath might be different in general
      if (dim>1) i1 = intervals[ifl].FindWhichOneIs(type, bfl, iop_1); // this is the index in multidimensional interval
      double t1 = intervals[ifl].time(type,i1);       // the beginning of the interval
      pair<int,int> i_iop = intervals[ifl].NextTimeIndex_(1-type, t1, bfl); // index of the next operator of the opposite type
      int i2 = i_iop.first;                           // index in the multidimensional bath index
      int iop_2 = i_iop.second;                       // index when counting only this bfl type
      double t2 = intervals[ifl].time(1-type,i2);     // the end of the interval
      double seg_length = t2-t1;                      // length of the segment
      if (seg_length<0) seg_length += common::beta;
      double t1_new = t1 + seg_length*r1;             // new random time in the existsing interval
      double t2_new = t1 + seg_length*r2;             // new random time in the existsing interval
      if (type==nIntervals::c){
	t_start = min(t1_new,t2_new);
	t_end = max(t1_new,t2_new);
      }else{
	t_end = min(t1_new,t2_new);
	t_start = max(t1_new,t2_new);
      }
      if (t_start>common::beta) t_start -= common::beta;
      if (t_end>common::beta) t_end -= common::beta;
      P_to_add = (sqr(seg_length)*knb)/(knb+1.)/2.;
    }
  }
  t_trial1.stop();
  
  if (t_start==t_end || Operators.FindExactTime(t_start)>=0) return; // There is already exactly the same time in the list

  t_add.start();
  
  t_trial2.start();
  pair<int,int> opt = Operators.Try_Add_Cd_C(flse, t_start, flse, t_end); 

  int op_i = min(opt.first, opt.second), op_f = max(opt.first, opt.second);
  double ratioD = TMD[ifl].AddDetRatio(MD[ifl], kn, t_start, bfl, t_end, bfl, intervals[ifl]);

  t_trial2.stop();
  t_trial3.start();
  double exp_sum = Segment_ComputeTryalExponentSum(state_evolution_left,op_i,Operators, npraStates, cluster, 2,istep);
  double ms = fabs(divide( Number(1.0, exp_sum), matrix_element));
  bool Accept = (1-rand()) < ms*fabs(ratioD)*P_to_add;
  //double det = fabs(detD[ifl]*ratioD);
  double det = fabs(ratioD);
  t_trial3.stop();

  if (Accept && ms>ratio_minM && det>ratio_minD){
    t_accept.start();
    int is, ie;
    intervals[ifl].Find_is_ie(t_start, is, t_end, ie); // SIGN: HERE WE HAVE is,ie
    ratioD = TMD[ifl].AddUpdateMatrix(MD[ifl], kn, t_start, is, bfl, t_end, ie, bfl, intervals[ifl],istep);
    pair<int,int> ipl = intervals[ifl].InsertExponents(t_start, is, bfl, t_end, ie, bfl);
    if (!common::Qsvd) TMD[ifl].AddUpdate_Gf(intervals[ifl], Mv,istep);
    pair<int,int> opt = Operators.Add_Cd_C(flse, t_start, flse, t_end,   IntervalIndex(ifl, nIntervals::cd, ipl.first), IntervalIndex(ifl, nIntervals::c, ipl.second));
    int op_i_ = min(opt.first, opt.second), op_f_ = max(opt.first, opt.second); // Changed March 2013
    // And updates state evolution for fast computation of Trace
    Number new_matrix_element = Segment_UpdateStateEvolution(state_evolution_left,op_i,Operators,npraStates, cluster, Trace, Prob, istep);
    // Operators.sign() :  gives number of crossings + number of forward propagators (sign for the matrix element to bring it to original form)
    // sign_isie: gives permutation of columns and rows of the determinant to bring it to original form
    //
    int sign_isie = 1-2*((is+ie)%2);
    int from_orig = Operators.sign()*sign_isie;
    double ms1 = divide(new_matrix_element, matrix_element);
    double Ps = ratioD*ms1*from_orig;
    gsign *= sigP<0>(Ps);
    detD[ifl] *= ratioD;
    matrix_element = new_matrix_element;
    successful++;
    successfullC[ifl]++;
    ComputeAverage(istep);
    t_accept.stop();
  }
  t_add.stop();
}
template <class Rand>
void CTQMC<Rand>::Segment_Remove_One_Kink(long long istep, int ifl)
{
  int kn = intervals[ifl].size()/2;
  if (kn==0) return;
  int dim = cluster.ifl_dim[ifl];           // dimension of this bath
  int bfl = static_cast<int>(dim*rand());   // type of the multidimensional bath we want to remove
  int knb = intervals[ifl].Nbtyp_s[bfl];    // how many segments of this type are there?
  if (knb==0) return;
  int flse = cluster.fl_from_ifl[ifl][bfl];
  int is, ie, iop_s, iop_e;
  double P_to_remove;

  t_trial1.start();
  if (AlgorithmType==OLD){
    P_to_remove = sqr(knb/common::beta);
    iop_s = static_cast<int>(knb*rand());  // removing this creation operator
    iop_e = static_cast<int>(knb*rand());  // removing this anhilation operator
    is = iop_s;                            // index in the multidimensional bath might be different in general
    if (dim>1) is = intervals[ifl].FindWhichOneIs(nIntervals::cd, bfl, iop_s); // this is the index in multidimensional interval
    ie = iop_e;
    if (dim>1) ie = intervals[ifl].FindWhichOneIs(nIntervals::c, bfl, iop_e); // this is the index in multidimensional interval
    if (!Segment_Try_Remove_C_Cd_(bfl, ie, bfl, is, intervals[ifl], istep)) return;
  }else if (AlgorithmType==WERNER){
    double t_start, t_end;
    double seg_length;
    if (rand()<0.5){   // Removing a segment
      iop_s = static_cast<int>(knb*rand());  // removing this creation operator
      is = iop_s;                            // index in the multidimensional bath might be different in general
      if (dim>1) is = intervals[ifl].FindWhichOneIs(nIntervals::cd, bfl, iop_s); // this is the index in multidimensional interval
      t_start = intervals[ifl].time_s(is);
      pair<int,int> i_iop = intervals[ifl].NextTimeIndex_(nIntervals::c, t_start, bfl); // next anhilation operator
      iop_e = i_iop.second;
      ie    = i_iop.first;
      double t_start_next = intervals[ifl].NextTime(nIntervals::cd, is, true);
      seg_length = t_start_next-t_start;
    }else{             // Removing an antisegment
      if (intervals[ifl].Nbtyp_s[bfl]<2) return; // Can not remove an anti-kink if there is only single kink
      iop_e = static_cast<int>(knb*rand());  // removing this destruction operator
      ie = iop_e;                            // index in the multidimensional bath might be different in general
      if (dim>1) ie = intervals[ifl].FindWhichOneIs(nIntervals::c, bfl, iop_e); // this is the index in multidimensional interval
      t_end = intervals[ifl].time_e(ie);
      pair<int,int> i_iop = intervals[ifl].NextTimeIndex_(nIntervals::cd, t_end, bfl); // next creation operator
      is    = i_iop.first;
      iop_s = i_iop.second;
      double t_end_next = intervals[ifl].NextTime(nIntervals::c, ie, true);
      seg_length = t_end_next-t_end;
    }
    P_to_remove = knb/(common::beta*seg_length);
  }else if (AlgorithmType==NEW){
    if (knb==1){
      P_to_remove = 1./sqr(common::beta);
      iop_s=0;
      is = iop_s;                            // index in the multidimensional bath might be different in general
      if (dim>1) is = intervals[ifl].FindWhichOneIs(nIntervals::cd, bfl, iop_s); // this is the index in multidimensional interval
      iop_e=0;
      ie = iop_e;
      if (dim>1) ie = intervals[ifl].FindWhichOneIs(nIntervals::c, bfl, iop_e); // this is the index in multidimensional interval
    }else{
      int type = static_cast<int>(rand()*2);          // either c or c^+
      int iop_1 = static_cast<int>(knb*rand());       // which operator from ordered list
      int i1 = iop_1;                                 // index in the multidimensional bath might be different in general
      if (dim>1) i1 = intervals[ifl].FindWhichOneIs(type, bfl, iop_1); // this is the index in multidimensional interval
      double t1 = intervals[ifl].time(type,i1);       // beginning of the interval
      pair<int,int> i_iop = intervals[ifl].NextTimeIndex_(1-type, t1, bfl); // index of the next operator of the opposite type
      int iop_2 = i_iop.second;                       // index when counting only this bfl type
      int i2 = i_iop.first;                           // index in the multidimensional bath index
      double t2 = intervals[ifl].time(1-type,i2);       // beginning of the interval
      double t1_previous = intervals[ifl].PreviousTime(1-type, i2, true);
      double t2_next = intervals[ifl].NextTime(type, i1, true);
      double seg_length = t2_next-t1_previous;         // length of the interval
      if (t2<t1) seg_length -= common::beta;           // jumped around twice when calculating t1_previous and t2_next
      if (type==nIntervals::cd){
	is=i1;  // if we started with cd, is=i1
	iop_s=iop_1;
	ie=i2;
	iop_e=iop_2;
      }else{
	ie=i1;
	iop_e=iop_1;
	is=i2;
	iop_s=iop_2;
      }
      P_to_remove = 2.*knb/(knb-1.)/sqr(seg_length);
    }
  }
  t_rm.start();
  t_trial1.stop();
  int ipe, ips;
  t_trial2.start();    
  ipe = Operators.FindSuccessive(2*flse+1, iop_e);
  ips = Operators.FindSuccessive(2*flse,   iop_s);
  if (ips>ipe) ips--;
  int op_i = min(ipe, ips);
  int op_f = ipe>ips ? ipe-2 : ips-1;
  Operators.Try_Remove_C_Cd(ipe, ips);
  double ratioD = TMD[ifl].RemoveDetRatio(MD[ifl], kn, is, ie);
  t_trial2.stop();

  t_trial3.start();
  double exp_sum = Segment_ComputeTryalExponentSum(state_evolution_left,op_i,Operators, npraStates, cluster, -2,istep);
  double ms = fabs(divide( Number(1.0, exp_sum), matrix_element));
  bool Accept = 1-rand() < ms*fabs(ratioD)*P_to_remove;
  //double det = fabs(detD[ifl]*ratioD);
  double det = fabs(ratioD);
  t_trial3.stop();

  if (Accept && ms>ratio_minM && det>ratio_minD){
    t_accept.start();
    Operators.Remove_C_Cd(ipe, ips);
    // And updates state evolution for fast computation of Trace
    Number new_matrix_element = Segment_UpdateStateEvolution(state_evolution_left, op_i, Operators,npraStates, cluster, Trace, Prob, istep);
    if (!common::Qsvd) TMD[ifl].RemoveUpdate_Gf(MD[ifl], kn, is, ie, intervals[ifl], Mv);
    int sign_isie = 1-2*((is+ie)%2);
    int from_orig = Operators.sign()*sign_isie;
    double ms1 = divide(new_matrix_element, matrix_element);
    double Ps = ratioD*ms1*from_orig;
    gsign *=sigP<0>(Ps);
    double ratioDn = TMD[ifl].RemoveUpdateMatrix(MD[ifl], kn, is, ie, intervals[ifl]);
    intervals[ifl].RemoveExponents(is,ie);
    detD[ifl] *= ratioDn;
    matrix_element = new_matrix_element;
    successful++;
    successfullC[ifl]++;
    ComputeAverage(istep);
    t_accept.stop();
  }
  t_rm.stop();
}

template <class Rand>
void CTQMC<Rand>::Segment_Move_A_Kink(long long istep, int ifl)
{
  int kn = intervals[ifl].size()/2;
  if (kn==0) return;
  t_trial1.start();
  int dim = cluster.ifl_dim[ifl];
  int type = static_cast<int>(rand()*2);
  int to_move = static_cast<int>(kn*rand());
  int bfl = intervals[ifl].btype(type, to_move);
  int fl = cluster.fl_from_ifl[ifl][bfl];
  int opera = 2*fl+type;
  int to_move_o = intervals[ifl].FindSuccessive(type, to_move, bfl);

  double t_old = intervals[ifl].time(type,to_move);
  pair<double,double> t_pn = intervals[ifl].Closest_Times(1-type, bfl, t_old);
  
  double t_new = t_pn.first + (t_pn.second-t_pn.first)*rand();
  if (t_new==t_pn.first || t_new==t_pn.second) return; // We just do not want to have two times exactly equal
  int iws=0;
  if (t_new<0) {t_new += common::beta; iws=1;}
  if (t_new>common::beta) {t_new -= common::beta; iws=1;}
  t_trial1.stop();

  t_mv.start();
  t_trial2.start();    
  double ratioD;
  if (type==0) ratioD = TMD[ifl].Move_start_DetRatio(MD[ifl], kn, t_old, t_new, bfl, to_move, intervals[ifl]);
  else ratioD = TMD[ifl].Move_end_DetRatio(MD[ifl], kn, t_old, t_new, bfl, to_move, intervals[ifl]);
  int ip_old, ipo, ip_new;
  Operators.TryMove(opera, t_old, to_move_o, ip_old, ipo, t_new, ip_new);
  int op_i = min(ipo, ip_new), op_f = max(ipo, ip_new);
  t_trial2.stop();    

  t_trial3.start();
  double P_part = fabs(ratioD);
  double exp_sum = Segment_ComputeTryalExponentSum(state_evolution_left,op_i,Operators, npraStates, cluster, 0,istep);
  double ms = fabs(divide( Number(1.0, exp_sum), matrix_element));
  bool Accept = 1-rand() < ms*P_part;
  //double det = fabs(detD[ifl]*ratioD);
  double det = fabs(ratioD);
  t_trial3.stop();
  
  if (Accept && ms>ratio_minM && det>ratio_minD){
    t_accept.start();
    Operators.Move(ip_old, ip_new);
    Number new_matrix_element = Segment_UpdateStateEvolution(state_evolution_left, op_i, Operators,npraStates, cluster, Trace, Prob, istep);
    int i_new = intervals[ifl].FindIndex(t_new, t_old, type);
    if (type==0) TMD[ifl].Move_start_UpdateMatrix(MD[ifl], kn, t_old, t_new, bfl, to_move, i_new, intervals[ifl]);
    else TMD[ifl].Move_end_UpdateMatrix(MD[ifl], kn, t_old, t_new, bfl, to_move, i_new, intervals[ifl]);
    intervals[ifl].MoveExponents(type, t_old, to_move, t_new, i_new);
    if (!common::Qsvd) TMD[ifl].MoveUpdate_Gf(intervals[ifl], Mv, istep);
    ratioD *= 1-2*((i_new+to_move+iws)%2);
    
    int sign_isie = 1-2*((i_new+to_move+iws)%2);
    int from_orig = Operators.sign()*sign_isie;
    double ms1 = divide(new_matrix_element, matrix_element);
    double Ps = ratioD*ms1*from_orig;
    gsign *= sigP<0>(Ps);
    
    detD[ifl] *= ratioD;
    matrix_element = new_matrix_element;
    successful++;
    successfullM[ifl]++;
    ComputeAverage(istep);
    t_accept.stop();
  }
  t_mv.stop();
}

template <class Rand>
void CTQMC<Rand>::Segment_Exchange_Two_Intervals(long long istep)
{
  int ir = static_cast<int>(rand()*g_exchange.size());
  int fla = g_exchange[ir].first;
  int flb = g_exchange[ir].second;
  int ifa  = cluster.ifl_from_fl[fla];
  int bfla = cluster.bfl_from_fl[fla];
  int ifb  = cluster.ifl_from_fl[flb];
  int bflb = cluster.bfl_from_fl[flb];
  int kna = intervals[ifa].size()/2;
  int knb = intervals[ifb].size()/2;
  if (kna==0 || knb==0) return;
  int ka = intervals[ifa].Nbtyp_s[bfla];    // how many segments of first type are there?
  int kb = intervals[ifb].Nbtyp_s[bflb];    // how many segments of second type are there?
  if (ka==0 || kb==0) return;               // no segments, can't exchange...
  int dim = cluster.ifl_dim[ifa];           // dimension of this bath
  int typea = static_cast<int>(rand()*2);   // either c or c+ for first segment
  int typeb = static_cast<int>(rand()*2);   // either c or c+ for second segment
  //int typeb = 1-typea;
  //int typeb=typea;
  int iop_a = static_cast<int>(ka*rand());       // which operator from ordered list
  int ia = iop_a;                                 // index in the multidimensional bath might be different in general
  if (dim>1) ia = intervals[ifa].FindWhichOneIs(typea, bfla, iop_a); // this is the index in multidimensional interval
  double t_a = intervals[ifa].time(typea,ia); // time in the first orbital

  pair<double,double> tv = intervals[ifa].Closest_Times(1-typea, bfla, t_a);
  pair<double,double> t_same;
  pair<int,int> ind_same;
  pair<double,double> t_oppo;
  pair<int,int> ind_oppo;
  intervals[ifb].Closest_Times(t_same, ind_same, 1-typeb, bflb, t_a);
  intervals[ifb].Closest_Times(t_oppo, ind_oppo, typeb, bflb, t_a);

  double t_b;
  int ib;
  int op_a = 2*fla+typea;
  int op_b = 2*flb+typeb;
  
  if (typea==nIntervals::c){
    if (t_same.first<t_oppo.first){ //    [tv.1]-----------------[t]    [tv.2]----------------
      if (t_oppo.first>tv.first){   //    ----[ts.1]     [to.1]------------[ts.2]  [to.2]-----
	t_b = t_oppo.first;
	ib = ind_oppo.first;
      } else return;
    }else{
      if (t_oppo.second<tv.second){//  [tv.1]--------[t]            [tv.2]----------------
	t_b = t_oppo.second;       //[to.1]----[ts.1]       [to.2]------------[ts.2]
	ib = ind_oppo.second;
      } else return;
    }
  }else{
    if (t_same.second>t_oppo.second){ // ----[tv.1]     [t]-----------[tv.2]
      if (tv.second>t_oppo.second){   // [to.1] [ts.1]-----[to.2]  [ts.2]-----
	t_b = t_oppo.second;
	ib = ind_oppo.second;
      } else return;
    }else{
      if (tv.first<t_oppo.first){ //----[tv.1]        [t]---------------[tv.2]  
	t_b = t_oppo.first;	  //[ts.1]---[to.1]         [ts.2]-----[to.2]   -----
	ib = ind_oppo.first;
      }else return;
    }
  }
  
  int iws=0;
  if (t_b<0) {t_b += common::beta; iws=1;}
  if (t_b>=common::beta) {t_b -= common::beta; iws=1;}
  t_b = intervals[ifb].time(typeb,ib); // BUG corrected 2019. Need very precise t_b, not approximate where beta is added and subtracted.
  if (t_a==t_b) return;

  t_mv.start();
  int ip_a, ip_b;
  if (t_a<t_b){
    ip_a = Operators.FindIndex(op_a, t_a, 0, istep);
    ip_b = Operators.FindIndex(op_b, t_b, ip_a, istep);
  }else{
    ip_b = Operators.FindIndex(op_b, t_b, 0, istep);
    ip_a = Operators.FindIndex(op_a, t_a, ip_b, istep);
  }

  double exp_sum = Segment_ComputeTryalExponentSumForExchange(ip_a,ip_b,Trace,state_evolution_left,Operators,npraStates,cluster,istep);
  double ms = fabs(divide( Number(1.0, exp_sum), matrix_element));
  
  t_trial2.start();    
  double ratioDa, ratioDb;
  if (typea==nIntervals::c)
    ratioDa = TMD[ifa].Move_end_DetRatio(  MD[ifa], kna, t_a, t_b, bfla, ia, intervals[ifa]);
  else
    ratioDa = TMD[ifa].Move_start_DetRatio(MD[ifa], kna, t_a, t_b, bfla, ia, intervals[ifa]);
  if (typeb==nIntervals::c)
    ratioDb = TMD[ifb].Move_end_DetRatio(  MD[ifb], knb, t_b, t_a, bflb, ib, intervals[ifb]);
  else
    ratioDb = TMD[ifb].Move_start_DetRatio(MD[ifb], knb, t_b, t_a, bflb, ib, intervals[ifb]);
  
  t_trial2.stop();
  double P_to_accept = ms*fabs(ratioDa*ratioDb);
  //double det_a = fabs(detD[ifa]*ratioDa);
  //double det_b = fabs(detD[ifb]*ratioDb);
  double det_a = fabs(ratioDa);
  double det_b = fabs(ratioDb);
  bool Accept = 1-rand() < P_to_accept;
  
  if (Accept && ms>ratio_minM && det_a>ratio_minD && det_b>ratio_minD){
    t_accept.start();
    Operators.ExchangeTwo(ip_a, ip_b);
    Number new_matrix_element = Segment_UpdateStateEvolution(state_evolution_left, min(ip_a,ip_b), Operators,npraStates, cluster, Trace, Prob, istep);
    if (fabs(fabs(divide(new_matrix_element,Number(1,exp_sum)))-1)>1e-5){
      cerr<<"at istep "<<istep<<" ERROR Different matrix elements m_new="<<new_matrix_element.exp_dbl()<<" m_old="<<exp_sum<<endl;
    }
    int ia_new=ia;
    int ib_new=ib;
    if (dim!=1 || iws!=0){
      ia_new = intervals[ifa].FindIndex(t_b, t_a, typea);
      ib_new = intervals[ifb].FindIndex(t_a, t_b, typeb);
    }
    if (typea==nIntervals::c)
      TMD[ifa].Move_end_UpdateMatrix(  MD[ifa], kna, t_a, t_b, bfla, ia, ia_new, intervals[ifa]);
    else
      TMD[ifa].Move_start_UpdateMatrix(MD[ifa], kna, t_a, t_b, bfla, ia, ia_new, intervals[ifa]);
    if (typeb==nIntervals::c)
      TMD[ifb].Move_end_UpdateMatrix(  MD[ifb], knb, t_b, t_a, bflb, ib, ib_new, intervals[ifb]);
    else
      TMD[ifb].Move_start_UpdateMatrix(MD[ifb], knb, t_b, t_a, bflb, ib, ib_new, intervals[ifb]);
    IntervalsExchangeExponents(typea, intervals[ifa], ia, ia_new, typeb, intervals[ifb], ib, ib_new);
    if (!common::Qsvd) TMD[ifa].MoveUpdate_Gf(intervals[ifa], Mv, istep);
    if (!common::Qsvd) TMD[ifb].MoveUpdate_Gf(intervals[ifb], Mv, istep);
    int sign_a_isie = 1-2*((ia_new+ia+iws)%2);
    int sign_b_isie = 1-2*((ib_new+ib+iws)%2);
    ratioDa *= sign_a_isie;
    ratioDb *= sign_b_isie;
    int from_orig = -1*sign_a_isie*sign_b_isie;
    double ms1 = divide(new_matrix_element, matrix_element);
    double Ps = ratioDa*ratioDb*ms1*from_orig;
    gsign *= sigP<0>(Ps);
    detD[ifa] *= ratioDa;
    detD[ifb] *= ratioDb;
    matrix_element = new_matrix_element;
    successful++;
    successfullM[ifa]++;
    successfullM[ifb]++;
    ComputeAverage(istep);
    t_accept.stop();
  }
  t_mv.stop();
  return;
}


template <class Rand>
void CTQMC<Rand>::Add_Two_Kinks(long long istep, int ifl)
{
  if (Operators.full()){
    int newNmax = static_cast<int>(Operators.max_size()*1.2+20);
    //cout<<"Increasing Nmax to "<<newNmax<<endl;
    Resize(newNmax); // always resize for at least 20 percent + 20 more entries
    //return;
  }
  int dim = cluster.ifl_dim[ifl];
  int kn = intervals[ifl].size()/2;
  int bfls1 = static_cast<int>(dim*rand());
  int bfle1 = static_cast<int>(dim*rand());
  int bfls2 = static_cast<int>(dim*rand());
  int bfle2 = static_cast<int>(dim*rand());
  int fls1 = cluster.fl_from_ifl[ifl][bfls1];
  int fle1 = cluster.fl_from_ifl[ifl][bfle1];
  int fls2 = cluster.fl_from_ifl[ifl][bfls2];
  int fle2 = cluster.fl_from_ifl[ifl][bfle2];

  double t_start1 = common::beta*rand();
  double t_end1 = common::beta*rand();
  double t_start2 = common::beta*rand();
  double t_end2 = common::beta*rand();
	  
  bool ssc1 = Operators.Try_Add_2_Cd_2_C_(fls1, t_start1, fle1, t_end1, fls2, t_start2, fle2, t_end2, state_evolution_left, state_evolution_right);

  if (t_start1==t_end1 || t_start2==t_end2 || t_start1==t_start2 || t_end1==t_end2) return;
  
  if (ssc1){
    double ratioD = TMD[ifl].AddDetRatio(MD[ifl], kn, bfls1, t_start1, bfle1, t_end1, bfls2, t_start2, bfle2, t_end2, intervals[ifl]);

    pair<int,int> opt = Operators.Try_Add_2_Cd_2_C(fls1, t_start1, fle1, t_end1, fls2, t_start2, fle2, t_end2);
    int op_i = min(opt.first, opt.second), op_f = max(opt.first, opt.second);
    Number matrix_element_new = ComputeTrace(state_evolution_left, state_evolution_right, op_i, op_f, Operators, npraStates, cluster, 4);
    double ms = fabs(divide(matrix_element_new, matrix_element));
    double P = sqr(common::beta*dim/(kn+1.))*sqr(common::beta*dim/(kn+2.))*(fabs(ratioD)*ms);
    if (1-rand()<P && ms>ratio_minM){
      int is1, ie1;
      intervals[ifl].Find_is_ie(t_start1, is1, t_end1, ie1);
      double ratioD1 = TMD[ifl].AddUpdateMatrix(MD[ifl], kn, t_start1, is1, bfls1, t_end1, ie1, bfle1, intervals[ifl],istep);
      pair<int,int> ipl1 = intervals[ifl].InsertExponents(t_start1, is1, bfls1, t_end1, ie1, bfle1);
      if (!common::Qsvd) TMD[ifl].AddUpdate_Gf(intervals[ifl], Mv,istep);
      
      pair<int,int> opt1 = Operators.Add_Cd_C(fls1, t_start1, fle1, t_end1, IntervalIndex(ifl, nIntervals::cd, ipl1.first), IntervalIndex(ifl, nIntervals::c, ipl1.second));
      
      int op_i1 = min(opt1.first, opt1.second);
      int op_f1 = max(opt1.first, opt1.second);
      Number new_matrix_element1 = UpdateStateEvolution(state_evolution_left, state_evolution_right, state_evolution_copy, op_i1, op_f1, Operators, npraStates, cluster, Trace, Prob, 2, istep);
      
      int is2, ie2;
      intervals[ifl].Find_is_ie(t_start2, is2, t_end2, ie2);
      double ratioD2 = TMD[ifl].AddUpdateMatrix(MD[ifl], kn+1, t_start2, is2, bfls2, t_end2, ie2, bfle2, intervals[ifl],istep);
      pair<int,int> ipl2 = intervals[ifl].InsertExponents(t_start2, is2, bfls2, t_end2, ie2, bfle2);
      if (!common::Qsvd) TMD[ifl].AddUpdate_Gf(intervals[ifl], Mv,istep);
      
      pair<int,int> opt2 = Operators.Add_Cd_C(fls2, t_start2, fle2, t_end2, IntervalIndex(ifl, nIntervals::cd, ipl2.first), IntervalIndex(ifl, nIntervals::c, ipl2.second));
	
      int op_i2 = min(opt2.first, opt2.second);
      int op_f2 = max(opt2.first, opt2.second);
      Number new_matrix_element2 = UpdateStateEvolution(state_evolution_left, state_evolution_right, state_evolution_copy, op_i2, op_f2, Operators, npraStates, cluster, Trace, Prob, 2, istep);
      
      double ms_new = divide(new_matrix_element2,matrix_element);
      if (fabs(fabs(ms)-fabs(ms_new))>1e-6) cerr<<"Matrix element not the same a2!"<<endl;
      double ratioD_new = ratioD1*ratioD2;
      if (fabs(fabs(ratioD_new)-fabs(ratioD))>1e-6) cerr<<"ratioD not the same!"<<endl;
	      
      detD[ifl] *= ratioD1*ratioD2;
      matrix_element = new_matrix_element2;
      successful++;
      successfullC[ifl]++;
	    
      ComputeAverage(istep);
    }
  }
}

template <class Rand>
void CTQMC<Rand>::Remove_Two_Kinks(long long istep, int ifl)
{
  int kn = intervals[ifl].size()/2;
  if (kn<=1) return;
  int dim = cluster.ifl_dim[ifl];
  int ie1 = static_cast<int>(kn*rand());
  int is1 = static_cast<int>(kn*rand());
  int bfle1 = intervals[ifl].btype_e(ie1);
  int bfls1 = intervals[ifl].btype_s(is1);
  int ie2, bfle2;
  do{
    ie2 = static_cast<int>(kn*rand());
    bfle2 = intervals[ifl].btype_e(ie2);
  } while(ie2==ie1);
  int is2, bfls2;
  do{
    is2 = static_cast<int>(kn*rand());
    bfls2 = intervals[ifl].btype_s(is2);
  } while (is2==is1);
  
  int fle1 = cluster.fl_from_ifl[ifl][bfle1];
  int fls1 = cluster.fl_from_ifl[ifl][bfls1];
  int fle2 = cluster.fl_from_ifl[ifl][bfle2];
  int fls2 = cluster.fl_from_ifl[ifl][bfls2];
	  
  int iop_s1 = intervals[ifl].FindSuccessive(0, is1, bfls1);
  int iop_e1 = intervals[ifl].FindSuccessive(1, ie1, bfle1);
  int iop_s2 = intervals[ifl].FindSuccessive(0, is2, bfls2);
  int iop_e2 = intervals[ifl].FindSuccessive(1, ie2, bfle2);
  
  int ipe1, ips1, ipe2, ips2;
  bool ssc1 = Operators.Try_Remove_2_C_2_Cd_(fle1, iop_e1, fls1, iop_s1, fle2, iop_e2, fls2, iop_s2, ipe1, ips1, ipe2, ips2, state_evolution_left, state_evolution_right);
    
  if (ssc1){
    double ratioD = TMD[ifl].RemoveDetRatio(MD[ifl], kn, is1, ie1, is2, ie2);
    int op_i = Operators.mini();
    int op_f = Operators.maxi()-2;
    Operators.Try_Remove_C_Cd(ipe1, ips1, ipe2, ips2);
    Number matrix_element_new = ComputeTrace(state_evolution_left, state_evolution_right, op_i, op_f, Operators, npraStates, cluster, -4);
    double ms = fabs(divide(matrix_element_new,matrix_element));
    double P = (fabs(ratioD)*ms)*sqr(kn/(common::beta*dim))*sqr((kn-1)/(common::beta*dim));
    //	    cout<<"r "<<setw(5)<<i<<" "<<setw(10)<<P<<" "<<setw(10)<<ms<<" "<<setw(10)<<ratioD<<endl;
    if (1-rand()<P && ms>ratio_minM){
      //      cout<<"r2 "<<setw(5)<<i<<" "<<setw(10)<<P<<" "<<setw(10)<<ms<<" "<<setw(10)<<ratioD<<" "<<kn<<endl;
      Operators.Remove_C_Cd(ipe1, ips1);
      int op_i = min(ipe1, ips1);
      int op_f = ipe1>ips1 ? ipe1-2 : ips1-1;
      Number new_matrix_element1 = UpdateStateEvolution(state_evolution_left, state_evolution_right, state_evolution_copy, op_i, op_f, Operators, npraStates, cluster, Trace, Prob, -2, istep);
      if (!common::Qsvd) TMD[ifl].RemoveUpdate_Gf(MD[ifl], kn, is1, ie1, intervals[ifl], Mv);
      double ratioD1 = TMD[ifl].RemoveUpdateMatrix(MD[ifl], kn, is1, ie1, intervals[ifl]);
      intervals[ifl].RemoveExponents(is1,ie1);

      Operators.Remove_C_Cd(ipe2, ips2);
      op_i = min(ipe2, ips2);
      op_f = ipe2>ips2 ? ipe2-2 : ips2-1;
      Number new_matrix_element2 = UpdateStateEvolution(state_evolution_left, state_evolution_right, state_evolution_copy, op_i, op_f, Operators, npraStates, cluster, Trace, Prob, -2, istep);
      if (is1<is2) is2--;
      if (ie1<ie2) ie2--;
      if (!common::Qsvd) TMD[ifl].RemoveUpdate_Gf(MD[ifl], kn-1, is2, ie2, intervals[ifl], Mv);
      double ratioD2 = TMD[ifl].RemoveUpdateMatrix(MD[ifl], kn-1, is2, ie2, intervals[ifl]);
      intervals[ifl].RemoveExponents(is2,ie2);
		
      double ms_new = divide(new_matrix_element2,matrix_element);
      if (fabs(fabs(ms)-fabs(ms_new))>1e-6) cerr<<"Matrix element not the same r2!"<<endl;
      double ratioD_new = ratioD1*ratioD2;
      if (fabs(fabs(ratioD_new)-fabs(ratioD))>1e-6) cerr<<"ratioD not the same!"<<endl;
		
      detD[ifl] *= ratioD1*ratioD2;
      matrix_element = new_matrix_element2;
      successful++;
      successfullC[ifl]++;
	    
      ComputeAverage(istep);
    }
  }
}

template <class Rand>
void CTQMC<Rand>::ComputeFinalAverage(function1D<double>& mom0, double& nf, double& TrSigma_G, double& Epot, double& Ekin)
{
  double Eimp_na = 0;
  (*xout)<<"mom="<<left;
  mom0.resize(cluster.N_unique_fl*cluster.HF_M.size());
  for (int op=0; op<cluster.HF_M.size(); op++){
    for (int fl=0; fl<cluster.N_unique_fl; fl++){
      double sum=0;
      for (int ist=0; ist<npraStates.size(); ist++)
	for (int m=0; m<cluster.msize(ist+1); m++)
	  sum += cluster.HF_M[op][fl][ist+1][m][m]*AProb[ist][m];
      int ind = fl+op*cluster.N_unique_fl;
      mom0[ind] = sum; // This is density matrix !
      (*xout)<<setw(8)<<mom0[ind]<<" ";
      // The first operator should be N_{bath}
      // We compute N_{bath}*E_{bath}
      if (op==0) Eimp_na += sum * cluster.epsk[fl];
    }
  }
  (*xout)<<endl;

  Epot=0;
  nf=0;
  for (int ist=0; ist<npraStates.size(); ist++){
    double dmuN = common::mu*cluster.Ns[ist+1]; // correction for N_{bath}*E_{bath}
    double Nm = cluster.Ns[ist+1];
    for (int m=0; m<cluster.msize(ist+1); m++){
      //Epot += (cluster.Ene[ist+1][m]+dmuN)*AProb[ist][m]; // Ene contains E[i]_{input}-Ns[i]*mu+U*Ns[i]*(Ns[i]-1)/2.
      Epot += (cluster.Enes[cluster.Eind(ist+1,m)]+dmuN)*AProb[ist][m]; // Ene contains E[i]_{input}-Ns[i]*mu+U*Ns[i]*(Ns[i]-1)/2.
      nf += Nm*AProb[ist][m];
    }
  }
  (*xout)<<left<<"Epot="<<setw(10)<<Epot<<" "; // This contains  E[i]_{input}+U*Ns[i]*(Ns[i]-1)/2.

  double kt=0;
  for (int ifl=0; ifl<kaver.size(); ifl++) kt += kaver[ifl];
  Ekin = -kt/common::beta;
  (*xout)<<left<<" Ekin="<<setw(10)<<Ekin<<" ";
  
  TrSigma_G = Epot - Eimp_na;  // Only TrSigmaG is free from single-particle (crystal-field splitting) terms, and contains only interaction terms, i.e., true potential energy
  (*xout)<<left<<"Tr(Sigma*G)="<<setw(10)<<TrSigma_G<<" ";
}


template <int boson_ferm>
void ComputeFrequencyPart(const vector<nIntervals>& intervals, const NOperators& Operators, int omsize, function2D<dcomplex>& dexp)
{ // This function computes the following factor:
  //            exp(iom*tau_{l+1})-exp(iom(tau_l)
  //            for frequencies iom sampled and all kink times tau_l.
  int N = Operators.size();
  dexp.resize(N+1, omsize);
  double eom_beta = e_om_beta<boson_ferm>();
  
  if (N==0){
    for (int im=0; im<omsize; im++)  dexp(0,im) = eom_beta-1;// e^{iom*beta}-e^{iom*0}
    return;
  }
  
  {// first kink
    const funProxy<dcomplex>& exp1 = find_Op_in_intervals<boson_ferm>(0, intervals, Operators);
    funProxy<dcomplex>& dexp_ = dexp[0];
    for (int im=0; im<omsize; im++)  dexp_[im] = exp1[im]-1; // exp(iom*0)==1
  }
  for (int ip=1; ip<N; ip++){// most of the kinks
    const funProxy<dcomplex>& exp0 = find_Op_in_intervals<boson_ferm>(ip-1, intervals, Operators);
    const funProxy<dcomplex>& exp1 = find_Op_in_intervals<boson_ferm>(ip,   intervals, Operators);
    funProxy<dcomplex>& dexp_ = dexp[ip];
    for (int im=0; im<omsize; im++)  dexp_[im] = exp1[im]-exp0[im];
  }
  {// after last kink
    const funProxy<dcomplex>& exp0 = find_Op_in_intervals<boson_ferm>(N-1, intervals, Operators);
    funProxy<dcomplex>& dexp_ = dexp[N];
    for (int im=0; im<omsize; im++)  dexp_[im] = eom_beta-exp0[im]; // exp(iom*beta)==+-1
  }
}


void ComputeFrequencyPart_Slow(const NOperators& Operators, const mesh1D& iom, function2D<dcomplex>& bexp, function2D<dcomplex>& dexp)
{ // This function computes the following factor:
  //            exp(iom*tau_{l+1})-exp(iom(tau_l)
  //            for frequencies iom sampled and all kink times tau_l.
  int N = Operators.size();
  int nom = iom.size();
  dexp.resize(N+1, iom.size());

  double eom_beta = cos(iom[0]*common::beta);
    
  if (N==0){
    dexp=0; // e^{iom*beta}-e^{iom*0} ==0
    return;
  }
  // This could be done at every insertion or removal of kink
  // But here it is done for all kinks every time
  for (int ip=0; ip<N; ip++){
    double tau = Operators.t(ip);
    for (int im=0; im<nom; im++){
      double phase = iom[im]*tau;
      bexp(ip,im).Set(cos(phase),sin(phase));
    }
  }

  const funProxy<dcomplex>& exp1 = bexp[0];
  funProxy<dcomplex>& dexp_ = dexp[0];
  for (int im=0; im<nom; im++)  dexp_[im] = exp1[im]-1; // exp(iom*0)==1
  
  for (int ip=1; ip<N; ip++){
    const funProxy<dcomplex>& exp0 = bexp[ip-1];
    const funProxy<dcomplex>& exp1 = bexp[ip];
    funProxy<dcomplex>& dexp_ = dexp[ip];
    for (int im=0; im<nom; im++)  dexp_[im] = exp1[im]-exp0[im];
  }

  {
    const funProxy<dcomplex>& exp0 = bexp[N-1];
    funProxy<dcomplex>& dexp_ = dexp[N];
    for (int im=0; im<nom; im++)  dexp_[im] = eom_beta-exp0[im]; // exp(iom*beta)==1
  }
}
template <class Rand>
void CTQMC<Rand>::ComputeAverage(long long istep)
{
  // components of observables:
  // 0 : nf
  // 1 : chi_{zz} = <Sz(t) Sz>
  // 2 : chi(density-density) = <1/2 n(t) 1/2 n>
  // 3 : <Sz>
  if (observables.size()!=4) cerr<<"Resize observables!"<<endl;
  observables=0;

  if (common::SampleSusc){
    susc=0;
    //    ComputeFrequencyPart_Slow(Operators, biom, bexp, dexp);
    ComputeFrequencyPart<0>(intervals, Operators, biom.size(), dexp);
  }

  if (common::SampleTransitionP) P_transition = 0;
  
  // time differences can be precomputed
  int N = Operators.size();
  dtau.resize(N+1);
  if (N>0){
    dtau[0] = Operators.t(0);
    for (int ip=1; ip<N; ip++) dtau[ip] = Operators.t(ip)-Operators.t(ip-1);
    dtau[N] = common::beta-Operators.t(N-1);
  } else dtau[0] = common::beta;


  deque<int> contrib;
  for (int ist=0; ist<npraStates.size(); ist++){
    // Checking which states do not contribute to trace
    // This was here before, but now I removed state_evolution!!
    //if (state_evolution_left[ist].size()<=Operators.size() || state_evolution_left[ist].last().istate!=npraStates[ist].istate) continue;
    if (npraStates[ist].empty()) continue;
    if (Trace[ist].mantisa==0) continue;
    
    double P_i = fabs(divide(Trace[ist],matrix_element)); // probability for this state
    if (fabs(P_i)<1e-10) continue; // Probability negligible
    contrib.push_back(ist);
    int istate = npraStates[ist].istate;    
    double Ns = cluster.Ns[istate];
    double Sz = cluster.Sz[istate];
    double Ns2 = Ns/2.;
    
    double dnf = Ns*dtau[0];
    double dmg = Sz*dtau[0];
 
    if (common::SampleSusc){
      for (int im=1; im<dexp.size_Nd(); im++){
	dsusc(im,0) = dexp[0][im] * Sz;
	dsusc(im,1) = dexp[0][im] * Ns2;
      }
    }

    double fact = P_i/common::beta;
    //int c_state=istate;
    for (int ip=0; ip<Operators.size(); ip++){
      int jstate = state_evolution_left(ist,ip).istate;
      double Ns = cluster.Ns[jstate];
      double Sz = cluster.Sz[jstate];
      double Ns2 = Ns/2.;
      
      dnf += Ns*dtau[ip+1];
      dmg += Sz*dtau[ip+1];

      //if (common::SampleTransitionP) P_transition(c_state-1, Operators.typ(ip) ) += fact;
      //c_state = jstate;
      
      if (common::SampleSusc){
	const funProxy<dcomplex>& dexp_ = dexp[ip+1];
	for (int im=1; im<dexp.size_Nd(); im++){
	  dsusc(im,0) += dexp_[im] * Sz;
	  dsusc(im,1) += dexp_[im] * Ns2;
	}
      }
    }
    if (common::SampleTransitionP){
      int c_state=istate;
      for (int ip=0; ip<Operators.size(); ip++){
	int jstate = state_evolution_left(ist,ip).istate;
	P_transition(c_state-1, Operators.typ(ip) ) += fact;
	c_state = jstate;
      }
    }
    
    observables[0] += dnf*fact;
    observables[1] += sqr(dmg)*fact;
    observables[2] += sqr(0.5*dnf)*fact;
    observables[3] += dmg*fact;
   
    if (common::SampleSusc){
      for (int im=1; im<susc.size_Nd(); im++){
	susc(0,im) += norm(dsusc(im,0))*fact/sqr(biom[im]);
	susc(1,im) += norm(dsusc(im,1))*fact/sqr(biom[im]);
      }
    }
  }
  susc(0,0) = observables[1];
  susc(1,0) = observables[2];
}


template <class Rand>
void CTQMC<Rand>::CleanUpdate(long long istep)
{
  static function2D<dcomplex> backup_Gf(4,iom.size());
  static function2D<double> backup_MD(MD[0].size_N(),MD[0].size_Nd());

  
  for (int ifl=0; ifl<common::N_ifl; ifl++){
    backup_MD.resize(MD[ifl].size_N(),MD[ifl].size_Nd());
    backup_MD = MD[ifl];
	
    TMD[ifl].CleanUpdateMatrix(MD[ifl], intervals[ifl].size()/2, intervals[ifl], tMD);

    for (int il=0; il<MD[ifl].size_N(); il++)
      for (int jl=0; jl<MD[ifl].size_Nd(); jl++)
	if (fabs(backup_MD[il][jl]-MD[ifl][il][jl])>1e-5) cerr<<"Difference in MD at "<<istep<<" "<<il<<" "<<jl<<" "<<backup_MD[il][jl]-MD[ifl][il][jl]<<endl;
	
    backup_Gf = TMD[ifl].gf();
	
    if (!common::Qsvd) TMD[ifl].CleanUpdateGf(MD[ifl], intervals[ifl].size()/2, intervals[ifl], Mv,istep);
	
    for (int il=0; il<backup_Gf.size_N(); il++)
      for (int jl=0; jl<backup_Gf.size_Nd(); jl++)
	if (norm(backup_Gf[il][jl]-TMD[ifl].gf()[il][jl])>1e-5) cerr<<"Difference in Gf at "<<istep<<" "<<il<<" "<<jl<<" "<<backup_Gf[il][jl]-TMD[ifl].gf()[il][jl]<<endl;
  }
}

template <class Rand>
void CTQMC<Rand>::StoreCurrentStateFast(function1D<double>& aver_observables, long long istep)
{

  int trusign = gsign;
  AProb.AddPart(Prob, trusign);
  aver_observables.AddPart(observables, trusign);
  if (common::SampleTransitionP) AP_transition.AddPart(P_transition, trusign);
  for (int ifl=0; ifl<kaver.size(); ifl++) kaver[ifl] += intervals[ifl].size()/2*trusign;
  if (common::SampleSusc) aver_susc.AddPart(susc, trusign);
  asign_fast += trusign;
  Naver++;  

  /*
  double dsum=0.0;
  for (int i=0; i<Prob.size_N(); i++){
    for (int m=0; m<Prob.size_Nd(); m++){
      dsum += Prob[i][m];
    }
  }
  if (fabs(dsum-1.0)>1e-3) {
    (*yout)<<"WARNING at istep="<<istep<<" sum(Probability)="<<dsum<<endl;
  }
  dsum=0.0;
  for (int i=0; i<Prob.size_N(); i++){
    for (int m=0; m<Prob.size_Nd(); m++){
      dsum += AProb[i][m];
    }
  }
  dsum *= 1./asign_fast;
  if (fabs(dsum-1.0)>1e-3) {
    (*yout)<<"WARNING2 at istep="<<istep<<" sum(AProbability)="<<dsum<<endl;
  }
  */
}
 
template <class Rand>
void CTQMC<Rand>::StoreCurrentState(long long istep)
{
  int trusign = gsign;
  static double  n_brisi[2]={0,0};
  static double nn=0;
  //AProb.AddPart(Prob, trusign);
  //aver_observables.AddPart(observables, trusign);
  //if (common::SampleTransitionP) AP_transition.AddPart(P_transition, trusign);
  //for (int ifl=0; ifl<kaver.size(); ifl++) kaver[ifl] += intervals[ifl].size()/2*trusign;
  
  asign += trusign;

  if (!common::Qsvd){
    
    for (size_t ifl=0; ifl<Gaver.size(); ifl++) Gaver[ifl].AddPart(TMD[ifl].gf(), trusign);
    
    if (common::QHB2){
      int N_ifl = cluster.N_ifl;
      for (int ip=0; ip<Operators.size(); ip++){// We go over all operators of all baths
	int op = Operators.typ(ip);
	if (op%2 == nIntervals::cd){
	  //int fls = op/2;
	  IntervalIndex p_ifl = Operators.iifl(ip);
	  int ifl = p_ifl.ifl; // type of bath
	  int iis = intervals[ifl].index_s_1[p_ifl.in]; // which contructor
	  if (fabs(intervals[ifl].time_s(iis)-Operators.t(ip))>1e-10) (*yout)<<"Times are different in StoreCurrentState t_interval="<<intervals[ifl].time_s(iis)<<" t_operators="<<Operators.t(ip)<<endl;
	  double Njcs1 = Get_N_at_Operator3(ip,Trace,matrix_element,state_evolution_left,state_evolution_right,Operators,npraStates,cluster, istep);
	  double nj_sign = Njcs1*trusign;
	  
	  //for (int ind=0; ind<cluster.N_baths(ifl); ind++){
	  //  dcomplex* __restrict__ _Ff = TMD[ifl].Ff[ind][iis].MemPt();
	  //  int fl_ai = cluster.v2fl_index[ifl][ind]; 
	  dcomplex* __restrict__ _Ff = TMD[ifl].Ff[iis].MemPt(); // needs to be generalized as above
	  int fl_ai = cluster.v2fl_index[ifl][0];                // needs to be generalized as above
	  for (int im=0; im<iom.size(); im++) Faver(fl_ai,im) += _Ff[im]*nj_sign;
	  //}
	}
      }
    }

    if (common::SampleSusc && cluster.DOsize>0){
      suscg=0.0;
      GeneralizedSusceptibility(suscg, Trace, matrix_element, biom, intervals, state_evolution_left, state_evolution_right, Operators, npraStates, cluster, istep);
      asuscg.AddPart(suscg,trusign);
    }
  }else{
    if (!common::QHB2 && !common::cmp_vertex){
      // We calculate only the Green's function
      for (size_t ifl=0; ifl<Gsvd.size(); ifl++) TMD[ifl].ComputeG_svd(Gsvd[ifl], trusign, svd, MD[ifl], intervals[ifl]);
      NGta += trusign;
    }else{
      int N_ifl = cluster.N_ifl;
      if (common::QHB2){// Need Nj before psi^+ operators
	t_nj.start();
	for (size_t ifl=0; ifl<N_ifl; ifl++){
	  // resizing the containers
	  int Nk = intervals[ifl].size()/2;// number of psi^\dagger kinks
	  Njcs1[ifl].resize(Nk);
	}
	for (int ip=0; ip<Operators.size(); ip++){// We go over all operators of all baths to calculate N at psi^dagger(tau)
	  int op = Operators.typ(ip);
	  if (op%2 == nIntervals::cd){
	    IntervalIndex p_ifl = Operators.iifl(ip);
	    int ifl = p_ifl.ifl; // type of bath
	    int iis = intervals[ifl].index_s_1[p_ifl.in]; // which contructor
	    if (fabs(intervals[ifl].time_s(iis)-Operators.t(ip))>1e-10) (*yout)<<"Times are different in StoreCurrentState t_interval="<<intervals[ifl].time_s(iis)<<" t_operators="<<Operators.t(ip)<<endl;
	    Njcs1[ifl][iis] = Get_N_at_Operator3(ip,Trace,matrix_element,state_evolution_left,state_evolution_right,Operators,npraStates,cluster, istep);
	  }
	}
	t_nj.stop();
      }
      
      int ifl_dim2_max = cluster.N_baths(0);
      for (int ifl=1; ifl<N_ifl; ifl++)
	if (ifl_dim2_max<cluster.N_baths(ifl))
	  ifl_dim2_max = cluster.N_baths(ifl);

      if (common::cmp_vertex){
	int lmax_lmax = svd.lmax*svd.lmax;
	for (size_t ifl=0; ifl<cluster.N_ifl; ifl++){
	  int Nk = intervals[ifl].size()/2;
	  t_fgf.start();
	  for (int ind=0; ind<cluster.v2fl_index[ifl].size(); ind++){
	    int ind2 = cluster.v2fl_index[ifl][ind];
	    Gsvd_s[ind2].resize(Nk,svd.lmax);
	    Gsvd_s[ind2] = 0.0;
	    Gsvd_e[ind2].resize(Nk,svd.lmax);
	    Gsvd_e[ind2] = 0.0;
	    Gsvd_es[ind2].resize(Nk,Nk,svd.lmax);
	    Gsvd_es[ind2] = 0.0;
	  }
	  
	  TMD[ifl].ComputeFGF_svd(Gsvd_s, Gsvd_e, Gsvd_es, Vertex_subtract, trusign, svd, MD[ifl], intervals[ifl]);  // the expensive part. Note that both Gsvd_s and Gsvd_e contain fermionic sign
	  t_fgf.stop();

	  t_g.start();
	  // The actual Green's function is just a straighforward sum over all starting times
	  for (int ind=0; ind<cluster.N_baths(ifl); ind++){
	    int fl_ai = cluster.v2fl_index[ifl][ind];
	    // *** slover equivalent ****
	    //for (int is=0; is<Nk; is++) Gsvd[ifl][ind][:] += Gsvd_s[fl_ai][is][:];
	    //double*       __restrict__ _Gsvd_   = Gsvd[ifl][ind].MemPt();
	    for (int is=0; is<Nk; is++){
	      //const double* __restrict__ _Gsvd_s_ = Gsvd_s[fl_ai][is].MemPt();
	      //for (int l=0; l<svd.lmax; l++) _Gsvd_[l] += _Gsvd_s_[l];
	      Gsvd[ifl][ind] += Gsvd_s[fl_ai][is];
	    }
	  }
	  t_g.stop();
	  t_f.start();
	  if (common::QHB2){
	    for (int ind=0; ind<cluster.N_baths(ifl); ind++){
	      int fl_ai = cluster.v2fl_index[ifl][ind];
	      if (Nk<10){
		for (int is=0; is<Nk; is++){
		  double njcs = Njcs1[ifl][is];
		  for (int l=0; l<svd.lmax; l++) Fsvd1(fl_ai,l) += Gsvd_s[fl_ai](is,l) * njcs;
		}
	      }else{
		Fsvd1[fl_ai].xgemv("T", Gsvd_s[fl_ai], Njcs1[ifl]);
	      }
	    }
	  }
	  t_f.stop();
	}
	
	t_vh.start();
	int lmax_x_lmax = svd.lmax*svd.lmax;
	function1D<double> G_e(svd.lmax), G_s(svd.lmax), bl(svd.lmax_bose);
	function2D<double> GG(svd.lmax,svd.lmax);
	for (int i0=0; i0<cluster.Nvfl; i0++){
	  int ifl0 = cluster.vfli_index[i0].ifl;
	  for (int i1=0; i1<cluster.Nvfl; i1++){
	    int ifl1 = cluster.vfli_index[i1].ifl;
	    for (int is=0; is<Gsvd_s[i0].size_N(); is++){
	      double t_s = intervals[ifl0].time_s(is);
	      for (int ie=0; ie<Gsvd_e[i1].size_N(); ie++){
		double t_e = intervals[ifl1].time_e(ie);
		double dt = t_s-t_e;
		if (dt<0) dt += common::beta;
		intpar ip = svd.tau.Interp(dt);
		//     _M0s_ = Gsvd_s[i0][is].MemPt(); // _M0s_[t_s,l] = \sum_{t_e} M[i0](t_e,t_s)*u[l](t_e-t_s) (-1)^{current_sign}
		//     _M1e_ = Gsvd_e[i1][ie].MemPt(); // _M1e_[t_e,l] = \sum_{t_s} M[i1](t_e,t_s)*u[l](t_e-t_s) (-1)^{current_sign}
		//     VH[i0,i1][lc,l1,l2] += _M0s_[t_s,l1] * _M1e_[t_e,l2] * uO[lc](t_s-t_e) (-1)^{current_sign}
		// hence
		// VH[i0,i1][lc,l1,l2] = \sum_{t_s1,t_e2} uO[lc](t_s1-t_e2) \sum_{t_e1,t_s2} M[i0](t_e1,t_s1)*u[l1](t_e1-t_s1) * M[i1](t_e2,t_s2)*u[l2](t_e2-t_s2)
		//                     = \sum_{t_s1,t_e1,t_s2,t_e2} M[i0](t_e1,t_s1) * M[i1](t_e2,t_s2) * u[l1](t_e1-t_s1) * u[l2](t_e2-t_s2) * uO[lc](t_s1-t_e2)
		//
		if (Vertex_subtract && i0==i1){
		  // This canceles the term when we intent to remove twice the same column and the same row. The Fock term should automatically remove this term,
		  // can be optimized by restrict pointers....
		  //t_vh1.start();
		  const double* __restrict__ _Gsvd_s_ = &Gsvd_s[i0](is,0);
		  const double* __restrict__ _Gsvd_e_ = &Gsvd_e[i0](ie,0);
		  const double* __restrict__ _Gsvd_es_ = &Gsvd_es[i0](ie,is,0);
		  double* __restrict__ _G_s_ = &G_s[0];
		  double* __restrict__ _G_e_ = &G_e[0];
		  for (int l=0; l<svd.lmax; l++){
		    _G_s_[l] = _Gsvd_s_[l] - _Gsvd_es_[l]; 
		    _G_e_[l] = _Gsvd_e_[l] - _Gsvd_es_[l];
		  }
		  //t_vh1.stop();
		  //t_vh2.start();
		  GG=0;
		  dger(GG.MemPt(), trusign, G_s, G_e);  // We need the sign again, because both Gsvd_s and Gsvd_e contain sign, hence the product is always +1.
		  //t_vh2.stop();
		  //t_vh3.start();
		  for (int lc=0; lc<svd.lmax_bose; lc++){
		    double blc = svd.fU_bose[lc](ip); 
		    //dger( &VH(i0,i1,lc,0,0),  blc, G_s,  G_e);
		    double* __restrict__ vh = &VH(i0,i1,lc,0,0);
		    const double* __restrict__ gg = GG.MemPt();
		    for (int ll=0; ll<lmax_x_lmax; ll++) vh[ll] += gg[ll]*blc;
		  }
		  //t_vh3.stop();
		}else{
		  //t_vh4.start();
		  GG=0;
		  dger(GG.MemPt(), trusign, Gsvd_s[i0][is], Gsvd_e[i1][ie]); // We need the sign again, because both Gsvd_s and Gsvd_e contain sign, hence the product is always +1.
		  //t_vh4.stop();
		  //t_vh5.start();
		  for (int lc=0; lc<svd.lmax_bose; lc++){
		    double blc = svd.fU_bose[lc](ip); 
		    //dger( &VH(i0,i1,lc,0,0),  blc, Gsvd_s[i0][is], Gsvd_e[i1][ie]);
		    double* __restrict__ vh = &VH(i0,i1,lc,0,0);
		    const double* __restrict__ gg = GG.MemPt();
		    for (int ll=0; ll<lmax_x_lmax; ll++) vh[ll] += gg[ll]*blc;
		  }
		  //t_vh5.stop();
		}
	      }
	    }
	  }
	}
	t_vh.stop();
      }else{ // only QHB2 but not cmp_vertex
	Gsvd_s.resize(ifl_dim2_max);
	for (size_t ifl=0; ifl<N_ifl; ifl++){
	  int Nk = intervals[ifl].size()/2;
	  for (int ind=0; ind<cluster.N_baths(ifl); ind++){
	    Gsvd_s[ind].resize(Nk,svd.lmax);
	    Gsvd_s[ind] = 0.0;
	  }

	  t_fgf.start();
	  TMD[ifl].ComputeFG_svd(Gsvd_s, trusign, svd, MD[ifl], intervals[ifl]); // the expensive part
	  t_fgf.stop();
	  
	  t_g.start();
	  // The actual Green's function is just a straighforward sum over all starting times
	  for (int ind=0; ind<cluster.N_baths(ifl); ind++)
	    for (int is=0; is<Nk; is++) Gsvd[ifl][ind] += Gsvd_s[ind][is];
	  t_g.stop();
	  
	  t_f.start();
	  for (int ind=0; ind<cluster.N_baths(ifl); ind++){
	    int fl_ai = cluster.v2fl_index[ifl][ind];
	    if (Nk<10){
	      for (int is=0; is<Nk; is++){
		double njcs = Njcs1[ifl][is];
		for (int l=0; l<svd.lmax; l++) Fsvd1(fl_ai,l) += Gsvd_s[ind](is,l) * njcs;
	      }
	    }else{
	      Fsvd1[fl_ai].xgemv("T", Gsvd_s[ind], Njcs1[ifl]);
	    }
	  }
	  t_f.stop();
	}
      }
    }
  }
}


template <class Rand>
void CTQMC<Rand>::PrintInfo(long long i, const function1D<double>& aver_observables)
{
  bool PRINT_DEBUG=false;
  (*yout)<<setw(9)<<i+1<<" ";
  double asignd = asign_fast/(Naver+0.0);
  double nf = aver_observables[0]/cluster.Nm/asign_fast;
  double mf = aver_observables[3]/cluster.Nm/asign_fast;

  yout->precision(4);
  (*yout)<<left<<setw(3)<<common::my_rank<<"  ";
  (*yout)<<left<<" nf="<<setw(7)<<nf<<"  ";
  (*yout)<<left<<" <Sz>="<<setw(8)<<mf<<"  ";
  (*yout)<<left<<" <s>="<<setw(8)<<asignd<<"  ";
  double chiS = aver_observables[1]/cluster.Nm/asign_fast;
  double chiC = aver_observables[2]/cluster.Nm/asign_fast;
  (*yout)<<left<<" chiS0="<<setw(12)<<(chiS)<<" ";
  //cout<<left<<" chiS="<<setw(12)<<(chiS-sqr(mf)*common::beta)<<"  ";
  (*yout)<<left<<" chiD="<<setw(12)<<(chiC-sqr(0.5*nf)*common::beta)<<"  ";
  
  {
    double Epot=0;
    for (int ist=0; ist<npraStates.size(); ist++){
      double dmuN = common::mu*cluster.Ns[ist+1];
      for (int m=0; m<cluster.msize(ist+1); m++)
	//Epot += (cluster.Ene[ist+1][m]+dmuN)*AProb[ist][m];
	Epot += (cluster.Enes[cluster.Eind(ist+1,m)]+dmuN)*AProb[ist][m];
    }
    Epot/= asign_fast;
    (*yout)<<left<<" Epot="<<setw(10)<<Epot<<" ";

    double kt=0;
    for (int ifl=0; ifl<kaver.size(); ifl++) kt += kaver[ifl];
    kt/= asign_fast;
    double Ekin = -kt/common::beta;
    (*yout)<<left<<" Ekin="<<setw(10)<<Ekin<<" ";
  }
  
  for (int op=0; op<cluster.HF_M.size(); op++){
    for (int fl=0; fl<cluster.N_unique_fl; fl++){
      
      double sum=0;
      for (int ist=0; ist<npraStates.size(); ist++)
	for (int m=0; m<cluster.msize(ist+1); m++)
	  sum += cluster.HF_M[op][fl][ist+1][m][m]*AProb[ist][m];
      sum/= asign_fast;

      (*yout)<<setw(8)<<sum<<" ";
      if (sum<0) PRINT_DEBUG=true;
    }
  }
  (*yout)<<right<<endl;

  if (PRINT_DEBUG){
    cout<<"asign="<<asign_fast<<endl;
    cout<<"Probabilities are exponent: ";
    for (int ist=0; ist<npraStates.size(); ist++)
      for (int m=0; m<cluster.msize(ist+1); m++)
	cout<<ist<<" "<<m<<" "<<AProb[ist][m]<<"  "<<Operators.exp_(0)[cluster.Eind(ist+1,m)]<<endl;
    
  }
  
  if (common::fastFilesystem && (i+1)%common::Naver==0){
    ofstream tout(NameOfFile_("Gaver",common::my_rank,i+1).c_str());
    tout.precision(16);
    for (int im=0; im<iom.size(); im++){
      tout<<setw(20)<<iom[im]<<" ";
      for (int ifl=0; ifl<common::N_ifl; ifl++)
	for (int b=0; b<Gaver[ifl].size_N(); b++)
	  tout<<setw(20)<<Gaver[ifl][b][im]/(asign)<<" ";
      tout<<endl;
    }
    ofstream pout(NameOfFile_("Probability",common::my_rank,i+1).c_str());
    pout.precision(16);
    for (int j=0; j<AProb.size_N(); j++)
      for (int k=0; k<AProb[j].size(); k++)
	pout<<setw(3)<<j+1<<" "<<setw(3)<<k<<" "<<setw(20)<<AProb[j][k]/(asign_fast)<<endl;

    if (common::SampleGtau){
      ofstream tout(NameOfFile_("Gatau",common::my_rank,i+1).c_str());
      tout.precision(16);
    
      for (int it=0; it<tau.size(); it++){
	tout<<tau[it]<<" ";
	for (int ifl=0; ifl<cluster.N_ifl; ifl++){
	  for (int ib=0; ib<cluster.N_baths(ifl); ib++){
	    double gf = Gtau[ifl](ib,it);
	    if (cluster.conjg[ifl][ib]) gf = -Gtau[ifl](ib,tau.size()-it-1);
	    gf *= cluster.sign[ifl][ib];
	    tout<<setw(20)<<gf<<" ";
	  }
	}
	tout<<endl;
      }
    }
  }
}

template<class container>
double variation(const container& gt, int dt=2)
{
  double gmax=0;
  //double gaver=0;
  int istart = max(static_cast<int>(0.05*gt.size()), dt);
  int iend = min(static_cast<int>(0.95*gt.size()), gt.size()-dt-1);
  //cout<<"istart="<<istart<<" iend="<<iend<<endl;
  if (iend>istart){
    for (int it=istart; it<iend; it++){
      double sm=0;
      for (int j=it-dt; j<it; j++) sm += gt[j];
      for (int j=it+1; j<it+dt+1; j++) sm += gt[j];    
      double ga = sm/(2*dt);
      if (fabs(ga)>1e-7 && fabs(gt[it])>1e-4) gmax = max(gmax, fabs(gt[it]/ga-1.));
      //if (fabs(ga)>1e-5) gaver += abs(gt[it]/ga-1.);
      //ga=(gt[it-1]+gt[it+1])/2.;
      //if (fabs(ga)>1e-5) gaver += gt[it]/ga;
    }
    //gaver *= 1./(iend-istart);
  }
  return gmax;
}


template <class Rand>
void CTQMC<Rand>::EvaluateVertex(long long i)
{
  int trusign = gsign;
  // We know that during sampling, the following was filled in
  //     Mv[(ifl,be,be)][w1,w2] = -1/beta \sum_{t_e,t_s} e^{i*w1*t_e-i*w2*t_s} M(t_e,t_s)
  // but only for positive w1. w2 was allowed to be positive and negative.
  // This part of the code fills in also negative w1 by using the relation
  //     Mv[w1,w2] = Mv[-w1,-w2]^*
  for (unsigned i0=0; i0<Mv.size(); i0++){
    for (int im1=0; im1<nomv; im1++){
      int nomv2 = 2*nomv-1;
      funProxy<dcomplex>& _Mvd = Mv[i0][im1];
      const funProxy<dcomplex>& _Mvs = Mv[i0][nomv2-im1];
      for (int im2=0; im2<2*nomv; im2++)
	_Mvd[im2] = _Mvs[nomv2-im2].conj();
    }
  }
      
  // Summary of Matsubara frequencies:
  //   we use the notation : nomv == N_w; nOm == N_W
  //   The fermionic frequencies are : om = (2*(iw-N_w)+1)*pi/beta   with iw=[0,....,2*N_w-1]  hence om=[ -(2*N_w-1)*pi/beta,....,(2*N_w-1)*pi/beta]
  //   The bosonic   frequencies are : Om = (2*(iW-N_W)+2)*pi/beta   with iW=[0,....,2*N_W-2]  hence Om=[ -2*(N_W-1)*pi/beta,....,2*(N_W-1)*pi/beta]
  //
  //   Example:
  //     iom-iOme can be obtained by im+dOm, where dOm is computed below:
  //     iom = (2*(im-N_w)+1)*pi/beta  and iOme = (2*(iOm-N_W)+2)*pi/beta
  //     so that iom-iOme = pi/beta*( 2*(im+ [N_W-iOm-1] -N_w)+1  )
  for (int i0=0; i0<VertexH.N0; i0++){
    for (int i1=0; i1<VertexH.N1; i1++){
      const function2D<dcomplex>& Mv0 = Mv[i0]; 
      const function2D<dcomplex>& Mv1 = Mv[i1]; 
      for (int iOm=0; iOm<2*nOm-1; iOm++){
	int dOm = nOm-1-iOm; // we get iw-iOm by im+dOm
	int sm1 = max(0, -dOm);
	int em1 = min(2*nomv, 2*nomv-dOm);
	for (int im1=sm1; im1<em1; im1++){
	  int sm2 = max(0, -dOm);
	  int em2 = min(2*nomv, 2*nomv-dOm);
	  dcomplex* VH = &VertexH(i0,i1,iOm,im1,0);
	  dcomplex* VF = &VertexF(i0,i1,iOm,im1,0);
	  dcomplex mv0 = Mv0(im1, im1+dOm);         // Mv[i0](iw1,iw1-iOm)
	  const funProxy<dcomplex>& _Mv0 = Mv0[im1];
	  for (int im2=sm2; im2<em2; im2++){
	    VH[im2] += mv0 * Mv1(im2+dOm, im2) * trusign;       // Mv[i0](iw1,iw1-iOm) * Mv[i1](iw2-iOm,iw2)
	    VF[im2] += _Mv0[im2]*Mv1(im2+dOm,im1+dOm) * trusign;// Mv[i0](iw1,iw2) * Mv[i1](iw2-iOm,iw1-iOm)
	  }
	}
      }
    }
  }
  NGtv += trusign;
}

template <class Rand>
void CTQMC<Rand>::write_times(ostream& out)
{
  out.precision(16);
  out<<"Operators contain"<<endl;
  for (int i=0; i<Operators.size(); i++){
    out<<i<<" "<<Operators.t(i)<<" "<<Operators.typ(i)<<endl;
  }
  out<<endl;
  for (int ifl=0; ifl<cluster.N_ifl; ifl++){
    out<<"interval "<<ifl<<" contains"<<endl;
    for (int m=0; m<intervals[ifl].size()/2; m++)
      out<<m<<" "<<intervals[ifl].time_s(m)<<" "<<intervals[ifl].time_e(m)<<endl;
    out<<endl;
  }
}
template <class Rand>
void CTQMC<Rand>::CheckTimes(int istep)
{
  cout.precision(16);
  for (int ifl=0; ifl<cluster.N_ifl; ifl++){
    int iis=0;
    int iie=0;
    for (int m=0; m<intervals[ifl].size()/2; m++){
      double ts = intervals[ifl].time_s(m);
      double te = intervals[ifl].time_e(m);
      while (iis<Operators.size() && Operators.t(iis)<ts) ++iis;
      if (Operators.t(iis)!=ts){
	cout<<"At "<<istep<<" have ts="<<ts<<" and can not find it in Operators ";
	if (iis>0) cout<<Operators.t(iis-1)<<",";
	cout<<Operators.t(iis)<<",";
	if (iis<Operators.size()-1) cout<<Operators.t(iis+1);
	cout<<endl;
      }
      while (iie<Operators.size() && Operators.t(iie)<te) ++iie;
      if (Operators.t(iie)!=te){
	cout<<"At "<<istep<<" have te="<<te<<" and can not find it in Operators ";
	if (iie>=2) cout<<Operators.t(iie-2)<<",";
	if (iie>=1) cout<<Operators.t(iie-1)<<",";
	cout<<Operators.t(iie)<<",";
	if (iie<Operators.size()-1) cout<<Operators.t(iie+1);
	cout<<endl;
	cout<<"The differences are"<<endl;
	if (iie>=2) cout<<Operators.t(iie-2)-te<<",";
	if (iie>=1) cout<<Operators.t(iie-1)-te<<",";
	cout<<Operators.t(iie)-te<<",";
	if (iie<Operators.size()-1) cout<<Operators.t(iie+1)-te;
	cout<<endl;
      }
    }
  }
}

template <class Rand>
double CTQMC<Rand>::sample(long long max_steps)
{
  Number Zloc = 0;
  for (int i=0; i<Trace.size(); i++) Zloc += Trace[i];
  
  matrix_element = Zloc;
  ComputeAverage(0);
  aver_observables=0;
  Timer t_measure;
  (*yout)<<"common::Segment="<<common::Segment<<endl;
  
  for (long long i=0; i<max_steps; i++){
    if (common::PChangeOrder>rand() || i<common::warmup/2){// change order, in warmup one should reach average order
      int ifl = BathP.ifl(rand());
      if (rand()>0.5){
	LOG(TRACE, i<<" Trying to add a kink: "<<ifl);
	if (common::Segment) Segment_Add_One_Kink(i, ifl);
	else Add_One_Kink(i, ifl);
      } else {
	LOG(TRACE, i<<" Trying to rem a kink: "<<ifl);
	if (common::Segment) Segment_Remove_One_Kink(i,ifl);
	else Remove_One_Kink(i, ifl);
      }
    }else{
      int ifl = BathP.ifl(rand());
      LOG(TRACE, i<<" Trying to move a kink: "<<ifl);
      if (common::Segment){
	if (common::PMove>1-rand())
	  Segment_Move_A_Kink(i, ifl);
	else
	  Segment_Exchange_Two_Intervals(i);
      } else{
	if (common::PMove>1-drand48())//1-rand()) // to keep the same status of RNG
	  Move_A_Kink(i, ifl);
	else
	  Exchange_Two_Intervals(i);
      }
    }
    
    if ((i+1)%common::CleanUpdate==0) CleanUpdate(i);
    if (common::GlobalFlip>0 && (i+1)%common::GlobalFlip==0) GlobalFlipFull(i);
    //if (common::GlobalFlip>0 && (i+1)%common::GlobalFlip==0) GlobalFlip(i);

    if (i>common::warmup){
      histogram[Operators.size()/2]++;
      if ((i+1)%common::tsampleFast==0){
	t_measure.start();
	StoreCurrentStateFast(aver_observables,i); // measuaring of probabiity and nf...
	t_measure.stop();
      }
      if ((i+1)%common::tsample==0){
	t_measure.start();
	StoreCurrentState(i); // most measuaring, G, susc, ....
	t_measure.stop();
      }
      if (common::SampleGtau>0 && (i+1)%common::SampleGtau==0){ // some tau sampling, but not necessary
	t_measure.start();
	// this is only for debuggig! It does not take into account properly the sign. It should depend on the sign, but we neglect that here.
	for (size_t ifl=0; ifl<TMD.size(); ifl++) TMD[ifl].ComputeGtau(MD[ifl], intervals[ifl], Gtau[ifl]);
	NGta++;
	t_measure.stop();
      }
      if (common::cmp_vertex && (!common::Qsvd) && (i+1)%common::SampleVertex==0){
	t_measure.start();
	EvaluateVertex(i);
	t_measure.stop();
      }
      if ((i+1)%common::Ncout==0 && i>2*common::warmup){
	PrintInfo(i, aver_observables);
#ifdef _TIME
	(*yout)<<"t_add="<<t_add.elapsed()<<" t_rm="<<t_rm.elapsed()<<" t_mv="<<t_mv.elapsed();
	(*yout)<<" t_accept="<<t_accept.elapsed()<<" t_measure="<<t_measure.elapsed();
	(*yout)<<" t_trial1="<<t_trial1.elapsed()<<" t_trial2="<<t_trial2.elapsed()<<" t_trial3="<<t_trial3.elapsed();
	//<<" t_evolve="<<NState::t_evolve.elapsed()<<" t_apply="<<NState::t_apply.elapsed()<<"; ";
	(*yout)<<" t_g="<<t_g.elapsed();
	if (common::QHB2) (*yout)<<" t_f="<<t_f.elapsed()<<" t_nj="<<t_nj.elapsed()<<" t_fgf="<<t_fgf.elapsed();
	if (common::SampleVertex>0) (*yout)<<" t_vh="<<t_vh.elapsed();
	//if (common::SampleVertex>0) (*yout)<<" t_vh1="<<t_vh1.elapsed()<<" t_vh2="<<t_vh2.elapsed()<<" t_vh3="<<t_vh3.elapsed()<<" t_vh4="<<t_vh4.elapsed()<<" t_vh5="<<t_vh5.elapsed()<<endl;
	(*yout)<<endl;
#endif  
      }
    }
  }

  if (common::Qsvd){
    for (size_t ifl=0; ifl<Gsvd.size(); ifl++) Gsvd[ifl] *= 1./asign;
    if (common::QHB2){
      /*
      if (common::Segment)
	for (size_t fl2=0; fl2<cluster.Nvfl; fl2++) Fsvd[fl2] *= 1./asign;
      else
      */
	for (size_t fl2=0; fl2<cluster.Nvfl; fl2++) Fsvd1[fl2] *= 1./asign;
    }
  }else{
    for (size_t ifl=0; ifl<Gaver.size(); ifl++) Gaver[ifl] *= 1./asign;
    if (common::QHB2)
      //for (int fl2=0; fl2<cluster.Nvfl; fl2++) Faver[fl2] *= 1./asign;
      Faver *= 1./asign;
    if (common::SampleSusc && cluster.DOsize>0) asuscg *= 1./asign;
  }
  
  kaver *= 1./asign_fast;
  aver_observables *= 1./(asign_fast*cluster.Nm);
  aver_observables[2] -= sqr(0.5*aver_observables[0])*common::beta;
  
  if (common::SampleSusc){
    aver_susc *= 1./(asign_fast*cluster.Nm);
    aver_susc[0][0] = aver_observables[1];
    aver_susc[1][0] = aver_observables[2];
  }
  if (common::cmp_vertex){
    if (common::Qsvd){
      VH *= 1./asign;
    }else{
      VertexH *= (1.0/NGtv);
      VertexF *= (1.0/NGtv);
    }
  }
  if (common::SampleGtau>0){
    double fct = (Gtau[0].size_Nd()/common::beta)/(NGta+0.0);
    for (size_t ifl=0; ifl<TMD.size(); ifl++) Gtau[ifl] *= fct;
  }


  //ofstream outv(NameOfFile_("Variation",common::my_rank,0).c_str());
  Variation.resize(cluster.N_ifl);
  if (common::SampleGtau>0){
    ostream &outv = cout;
    outv<<common::my_rank<<":variation=";
    for (int ifl=0; ifl<cluster.N_ifl; ifl++)
      for (int ib=0; ib<cluster.N_baths(ifl); ib++){
	double var = variation(Gtau[ifl][ib]);
	Variation[ifl].push_back( var );
	outv<<var<<" ";
      }
    outv<<endl;
  }else{
    for (int ifl=0; ifl<cluster.N_ifl; ifl++)
      for (int ib=0; ib<cluster.N_baths(ifl); ib++)
	Variation[ifl].push_back(0.0);
  }
  //outv.close();

  double succM=0, succAR;
  for (int ifl=0; ifl<common::N_ifl; ifl++){
    succM += successfullM[ifl];
    succAR += successfullC[ifl];
  }

  if (common::my_rank==0){
    (*xout)<<"Acceptance ratio: add/rm="<<succAR/(max_steps*common::PChangeOrder)<<"  mv="<<succM/(max_steps*(1.-common::PChangeOrder+1e-20))<<" number of accepted steps="<<successful<<" sign at finish "<<gsign<<endl;
    for (int ifl=0; ifl<common::N_ifl; ifl++){
      (*xout)<<"acceptance add/rm for orbital ["<<ifl<<"]="<<setw(7)<<std::left<<successfullC[ifl]/static_cast<double>(successful)<<" "<<std::right;
      (*xout)<<"acceptance  move  for orbital ["<<ifl<<"]="<<setw(7)<<std::left<<successfullM[ifl]/static_cast<double>(successful)<<" "<<std::right<<endl;
    }
#ifdef _TIME
    (*xout)<<"t_add="<<t_add.elapsed()<<" t_rm="<<t_rm.elapsed()<<" t_mv="<<t_mv.elapsed();
    (*xout)<<" t_accept="<<t_accept.elapsed()<<" t_measure="<<t_measure.elapsed();
    (*xout)<<" t_trial1="<<t_trial1.elapsed()<<" t_trial2="<<t_trial2.elapsed()<<" t_trial3="<<t_trial3.elapsed();
    (*xout)<<" t_g="<<t_g.elapsed();
    if (common::QHB2) (*xout)<<" t_f="<<t_f.elapsed()<<" t_nj="<<t_nj.elapsed()<<" t_fgf="<<t_fgf.elapsed();
    if (common::SampleVertex>0) (*xout)<<" t_vh="<<t_vh.elapsed();
    (*xout)<<endl;
#endif  
  }
  // HERE YOU SHOULD CHECK GTAU VARIATION
  return aver_observables[0];
}

template <class Rand>
bool CTQMC<Rand>::SaveStatus(int my_rank)
{
  ofstream out(NameOfFile("status",my_rank).c_str());
  out.precision(16);
  if (!out){cerr<<"Could not save status!"<<endl; return false;}
  out<<"# Status file for ctqmc. Used to save all times t_s, t_e."<<endl;
  out<<intervals.size()<<endl;
  for (size_t ifl=0; ifl<intervals.size(); ifl++){
    out<<intervals[ifl].size()/2<<endl;
    for (int it=0; it<intervals[ifl].size()/2; it++){
      out<<intervals[ifl].time_s(it)<<" "<<intervals[ifl].time_e(it)<<" "<<intervals[ifl].btype_s(it)<<" "<<intervals[ifl].btype_e(it)<<endl;
    }
  }
  return true;
}

int Nmax_from_StatusFiles(int my_rank)
{
  string str;
  ifstream inp(NameOfFile("status",my_rank).c_str());
  if (!inp) return false;

  int cNmax = 0;
  getline(inp,str); // comment
  int intervals_size;
  inp>>intervals_size;
  getline(inp,str);
  if (!inp) return false;
  for (size_t ifl=0; ifl<intervals_size; ifl++){
    int nt;
    inp>>nt;
    getline(inp,str);
    if (!inp) {clog<<"Only partly restored file status."<<my_rank<<endl; break;}
    if (nt<0){cerr<<"Wrong input times. Skipping...."<<endl; continue;}
    cNmax += nt;
    for (int it=0; it<nt; it++) getline(inp,str);// skipping all times
  }
  return cNmax*2+6;
}


template <class Rand>
bool CTQMC<Rand>::RetrieveStatus(int my_rank, long long istep)
{
  string str;
  ifstream inp(NameOfFile("status",my_rank).c_str());
  if (!inp) return false;
  getline(inp,str); // comment
  int tN_ifl;
  inp>>tN_ifl;
  getline(inp,str);
  if (!inp) return false;

  for (size_t ifl=0; ifl<intervals.size(); ifl++){
    int dim = cluster.ifl_dim[ifl];
    int nt;
    inp>>nt;
    getline(inp,str);
    if (!inp) {clog<<"Only partly restored file status."<<my_rank<<endl; break;}
    if (nt<0 || nt>=intervals[ifl].fullsize()/2){cerr<<"Too many input times. Skipping...."<<endl; continue;}
    for (int it=0; it<nt; it++){
      double t_start, t_end;
      inp>>t_start>>t_end;
      if (t_start<0 || t_start>common::beta || t_end<0 || t_end>common::beta) {cerr<<"The times from starting file are wrong. Skipping..."<<endl; continue;}
      int bfls, bfle;
      inp>>bfls>>bfle;
      if (bfls<0 || bfls>dim || bfle<0 || bfle>dim) {cerr<<"The baths from starting file are wrong. Skipping..."<<endl; continue;}
      getline(inp,str);
      
      int kn = intervals[ifl].size()/2;// should be same as it!
      int fls = cluster.fl_from_ifl[ifl][bfls];
      int fle = cluster.fl_from_ifl[ifl][bfle];
      
      int is, ie;      
      intervals[ifl].Find_is_ie(t_start, is, t_end, ie);// shoule be same as it
      double ratioD = TMD[ifl].AddUpdateMatrix(MD[ifl], kn, t_start, is, bfls, t_end, ie, bfle, intervals[ifl],-1);

      //if (isnan(ratioD)){
      //	cout<<"ratioD="<<ratioD<<" t_start"<<t_start<<" t_end="<<t_end<<" bls="<<bfls<<" bfle="<<bfle<<" is="<<is<<" ie="<<ie<<" kn="<<kn<<endl;
      //}
      
      pair<int,int> ipl = intervals[ifl].InsertExponents(t_start, is, bfls, t_end, ie, bfle);
      
      if (!common::Qsvd) TMD[ifl].AddUpdate_Gf(intervals[ifl], Mv,-1);
      
      pair<int,int> opt = Operators.Add_Cd_C(fls, t_start, fle, t_end, IntervalIndex(ifl, nIntervals::cd, ipl.first), IntervalIndex(ifl, nIntervals::c, ipl.second));
      
      int op_i = min(opt.first, opt.second);
      int op_f = max(opt.first, opt.second);

      Number new_matrix_element = UpdateStateEvolution(state_evolution_left, state_evolution_right, state_evolution_copy, op_i, op_f, Operators, npraStates, cluster, Trace, Prob, 2, istep, false);

      gsign *= sigP<0>(ratioD*Operators.sign());
      
      detD[ifl] *= ratioD;
      matrix_element = new_matrix_element;
    }
  }
  gsign *= sigP<0>(matrix_element.mantisa);
  
  (*yout)<<"Old status retrieved with sign="<<gsign<<" and determinants ";
  for (int ifl=0; ifl<common::N_ifl; ifl++){
    (*yout)<<ifl<<":"<<detD[ifl]<<" ";
    if (isnan(detD[ifl])) detD[ifl]=1.0;
    // It can happen that detD is nan, because we set up the matrix through a series of steps adding two kinks at a time. If we have a configuration, which can not be reached in
    // such steps (non-ergodic configuration for two kinks addition), than we will get ratio=0 at some step, and finally ratio=nan at the next step. Here we reset it, because
    // for MC evolution, it is irrelevant.
  }
  (*yout)<<endl;
  ComputeAverage(0);
  CleanUpdate(0);
  return true;
}

template <class Rand>
void CTQMC<Rand>::RecomputeCix()
{  cluster.RecomputeCix(AProb,asign_fast,common::treshold,npraStates);}

template <class Rand>
void CTQMC<Rand>::GlobalFlip(long long istep)
{
  for (int iu=0; iu<cluster.gflip_ifl.size(); iu++){
    int ifa = cluster.gflip_ifl[iu].first;
    int ifb = cluster.gflip_ifl[iu].second;
    int sizea = cluster.ifl_dim[ifa];
    int sizeb = cluster.ifl_dim[ifb];
    if (sizea!=sizeb) cerr<<"Sizes of symmetric baths are different! This is not supported!"<<endl;
    //    int ka = intervals[ifa].size()/2;
    //    int kb = intervals[ifb].size()/2;

    double ratioDb=1, ratioDa=1;
    if (!cluster.ifl_equivalent[ifa][ifb]){
      ratioDb = TMD[ifa].GlobalFlipDetRatio(intervals[ifa], MD[ifa], Deltat[ifb], /*Cf,*/ tMD);
      ratioDa = TMD[ifb].GlobalFlipDetRatio(intervals[ifb], MD[ifb], Deltat[ifa], /*Cf,*/ tMD);
    }
    Number matrix_element_new = ComputeGlobalTrace(Operators, npraStates, cluster, cluster.gflip_fl[iu]);
    
    double ms = divide(matrix_element_new,matrix_element);
    double P = fabs(ms*ratioDa*ratioDb);

    //cout<<"Global flip: ms="<<ms<<" ratioD="<<ratioDa*ratioDb<<endl;
    
    if (P>1-rand()){
      tinterval.copy_data(intervals[ifa]); intervals[ifa].copy_data(intervals[ifb]); intervals[ifb].copy_data(tinterval);//swap
      
      if (cluster.ifl_equivalent[ifa][ifb]){
	tMD = MD[ifa]; MD[ifa] = MD[ifb]; MD[ifb] = tMD;
	if (!common::Qsvd) swapGf(TMD[ifa], TMD[ifb], Mv);
      }else{
	TMD[ifa].CleanUpdateMatrix(MD[ifa], intervals[ifa].size()/2, intervals[ifa], tMD);
	TMD[ifb].CleanUpdateMatrix(MD[ifb], intervals[ifb].size()/2, intervals[ifb], tMD);
	if (!common::Qsvd) TMD[ifa].CleanUpdateGf(MD[ifa], intervals[ifa].size()/2, intervals[ifa], Mv,istep);
	if (!common::Qsvd) TMD[ifb].CleanUpdateGf(MD[ifb], intervals[ifb].size()/2, intervals[ifb], Mv,istep);
      }

      Operators.GlobalFlipChange(cluster.gflip_fl[iu], ifa, ifb);
      
      Number new_matrix_element = UpdateStateEvolution(state_evolution_left, state_evolution_right, state_evolution_copy, 0, Operators.size()-1, Operators, npraStates, cluster, Trace, Prob, 0, istep);
      if (fabs(fabs(divide(matrix_element_new,new_matrix_element))-1)>1e-6) {cerr<<"Matrix elements are not the same in global flip "<<matrix_element_new<<" "<<new_matrix_element<<endl;}
      
      swap(detD[ifa],detD[ifb]);
      detD[ifa] *= ratioDa;
      detD[ifb] *= ratioDb;
      
      matrix_element = new_matrix_element;
      successful++;
      ComputeAverage(0);
    }
  }
}

template <class Rand>
void CTQMC<Rand>::GlobalFlipFull(long long istep)
{
  if (cluster.gflip_ifl.size()==0) return; // Nothing to flip.

  //cout<<"Trying global flip"<<endl;
  vector<pair<double,double> > Dratios(cluster.gflip_ifl.size());
  double ratioD=1;
  for (int iu=0; iu<cluster.gflip_ifl.size(); iu++){
    int ifa = cluster.gflip_ifl[iu].first;
    int ifb = cluster.gflip_ifl[iu].second;
    int sizea = cluster.ifl_dim[ifa];
    int sizeb = cluster.ifl_dim[ifb];
    if (sizea!=sizeb) cerr<<"Sizes of symmetric baths are different! This is not supported!"<<endl;    
    double ratioDb=1, ratioDa=1;
    if (!cluster.ifl_equivalent[ifa][ifb]){
      ratioDb = TMD[ifa].GlobalFlipDetRatio(intervals[ifa], MD[ifa], Deltat[ifb], /*Cf,*/ tMD);
      ratioDa = TMD[ifb].GlobalFlipDetRatio(intervals[ifb], MD[ifb], Deltat[ifa], /*Cf,*/ tMD);
    }
    Dratios[iu].first = ratioDa;
    Dratios[iu].second = ratioDb;
    ratioD *= ratioDa*ratioDb;
  }

  //cout<<"RatioD="<<ratioD<<endl;
  //cout<<"size="<<cluster.gflip_fl.size_N()<<endl;
  
  Number matrix_element_new = ComputeGlobalTrace(Operators, npraStates, cluster, cluster.gflip_fl[cluster.gflip_fl.size_N()-1]);
  double ms = divide(matrix_element_new,matrix_element);
  double P = fabs(ms*ratioD);

  //cout<<"Global flip on all pairs of baths  ratioD="<<ratioD<<" ms="<<ms<<endl;//" check="<<divide(mm_brisi,matrix_element)<<endl;
  
  if (P>1-rand()){
    (*yout)<<"Global flip succeeded: P="<<P<<" ratioD="<<ratioD<<" ms="<<ms<<endl;
    for (int iu=0; iu<cluster.gflip_ifl.size(); iu++){
      int ifa = cluster.gflip_ifl[iu].first;
      int ifb = cluster.gflip_ifl[iu].second;
      tinterval.copy_data(intervals[ifa]); intervals[ifa].copy_data(intervals[ifb]); intervals[ifb].copy_data(tinterval);//swap
      if (cluster.ifl_equivalent[ifa][ifb]){
	tMD = MD[ifa]; MD[ifa] = MD[ifb]; MD[ifb] = tMD;
	if (!common::Qsvd) swapGf(TMD[ifa], TMD[ifb], Mv);
      }else{
	TMD[ifa].CleanUpdateMatrix(MD[ifa], intervals[ifa].size()/2, intervals[ifa], tMD);
	TMD[ifb].CleanUpdateMatrix(MD[ifb], intervals[ifb].size()/2, intervals[ifb], tMD);
	if (!common::Qsvd) TMD[ifa].CleanUpdateGf(MD[ifa], intervals[ifa].size()/2, intervals[ifa], Mv,istep);
	if (!common::Qsvd) TMD[ifb].CleanUpdateGf(MD[ifb], intervals[ifb].size()/2, intervals[ifb], Mv,istep);
      }
    }

    Operators.GlobalFlipChange(cluster.gflip_fl[cluster.gflip_fl.size_N()-1], cluster.gflip_ifl);
      
    Number new_matrix_element = UpdateStateEvolution(state_evolution_left, state_evolution_right, state_evolution_copy, 0, Operators.size()-1, Operators, npraStates, cluster, Trace, Prob, 0, istep);

    if (fabs(fabs(divide(matrix_element_new,new_matrix_element))-1)>1e-6) {cerr<<"Matrix elements are not the same in global flip "<<matrix_element_new<<" "<<new_matrix_element<<endl;}

    for (int iu=0; iu<cluster.gflip_ifl.size(); iu++){
      int ifa = cluster.gflip_ifl[iu].first;
      int ifb = cluster.gflip_ifl[iu].second;
      double ratioDa = Dratios[iu].first;
      double ratioDb = Dratios[iu].second;
      swap(detD[ifa],detD[ifb]);
      detD[ifa] *= ratioDa;
      detD[ifb] *= ratioDb;
    }
    double mms = divide(new_matrix_element,matrix_element);
    matrix_element = new_matrix_element;
    successful++;
    ComputeAverage(0);
  }
}

bool Qnonzero(const functionb<dcomplex>& Gf){
  // Check that the function is nonzero
  for (int im=0; im<min(10,Gf.size()); im++)
    if (abs(Gf[im])<1e-10) return false;
  double wsum=0;
  for (int im=0; im<min(10,Gf.size()); im++) wsum += abs(Gf[im]);
  if (wsum<1e-5) return false;
  return true;
}

struct CustomLess{
  function2D<int>* pF_i;
  int ist;
  CustomLess(function2D<int>* pF_i_, int ist_) : pF_i(pF_i_), ist(ist_){}
  bool operator()(int op1, int op2)
  {   
    return pF_i->operator()(op1,ist) < pF_i->operator()(op2,ist);
  }
};


template <class Rand>
void CTQMC<Rand>::WriteVertexSVD(const function2D<dcomplex>& Gd, const function2D<double>& Gd_, int my_rank, int mpi_size, int Master, bool print_vertex_xxx)
{
  //function2D<double> Gd_(cluster.N_unique_fl,svd.lmax);
  // Summary of Matsubara frequencies:
  //   we use the notation : nomv == N_w; nOm == N_W
  //   The fermionic frequencies are : om = (2*(iw-N_w)+1)*pi/beta   with iw=[0,....,2*N_w-1]  hence om=[ -(2*N_w-1)*pi/beta,....,(2*N_w-1)*pi/beta]
  //   The bosonic   frequencies are : Om = (2*(iW-N_W)+2)*pi/beta   with iW=[0,....,2*N_W-2]  hence Om=[ -2*(N_W-1)*pi/beta,....,2*(N_W-1)*pi/beta]
  //
  //   Example:
  //     iom-iOme can be obtained by im+dOm, where dOm is computed below:
  //     iom = (2*(im-N_w)+1)*pi/beta  and iOme = (2*(iOm-N_W)+2)*pi/beta
  //     so that iom-iOme = pi/beta*( 2*(im+ [N_W-iOm-1] -N_w)+1  )

  /*
  function2D<double> Cll;
  if (false){
    svd.ComputeClCoefficients(Cll, my_rank, mpi_size, Master);
  }
  */
  
  int Nl = svd.lmax_bose*svd.lmax*svd.lmax;

  ofstream outv("svd_vertex.dat"); outv.precision(12);
  outv<<"# SVD parameters:  lmax  lmax_bose  beta  Ntau  Nw  L  x0"<<endl;
  outv<<svd.lmax<<" "<<svd.lmax_bose<<" "<<common::beta<<" "<<tau.size()<<" "<< (svd.om.size()/2) <<" "<<svd.L<<" "<<svd.x0<<endl;
  outv<<" Nvfl  N_unique_fl  Nl "<<endl;
  outv<<cluster.Nvfl<<" "<<cluster.N_unique_fl<<" "<<Nl<<endl;

  /*
  outv<<"# fl  l    G_l :"<<endl;
  for (int i=0; i<cluster.N_unique_fl; i++){
    for (int l=0; l<svd.lmax; l++)
      outv<<setw(3)<<i<<" "<<setw(3)<<l<<"   "<<Gd_(i,l)<<endl;
  }
  */

  outv<<"# fl  l    G_l :"<<endl;
  for (int l=0; l<svd.lmax; l++){
    outv<<setw(3)<<l<<" ";
    for (int i=0; i<cluster.Nvfl; i++){
      int ifl = cluster.vfli_index[i].ifl;
      int bfl = cluster.vfli_index[i].bind;
      outv<<setw(20)<<Gsvd[ifl](bfl,l)<<" ";
    }
    outv<<endl;
  }
  
  outv<<"# i0  i1  l2  l1  l0   VH :"<<endl;
  for (int i0=0; i0<cluster.Nvfl; i0++){
    for (int i1=0; i1<cluster.Nvfl; i1++){
      for (int i=0; i<Nl; i++){
	vector<int> li = svd.get_lll(i);
	outv<<setw(3)<<i0<<" "<<setw(3)<<i1<<" "<<setw(3)<<li[2]<<" "<<setw(3)<<li[1]<<" "<<setw(3)<<li[0]<<"    "<<setw(15)<<VH(i0,i1,li[2],li[1],li[0])<<endl;
      }
    }
  }
  
  if (!print_vertex_xxx) return;
  return ;


    function2D<dcomplex> ul_om;
    function2D<dcomplex> ul_Om;
    int N_W = max(nOm,2*nomv);
    int N_w = nomv; 

    N_w = nomv+nOm-1;  // increase it, so that we do not need to cut-off frequencies in Fock term
    ul_om.resize(2*N_w,svd.lmax);
    ul_Om.resize(2*N_W-1,svd.lmax_bose);
    {
      function2D<dcomplex> fUw(svd.lmax,N_w), bUw(svd.lmax_bose,N_W);
      for (int l=0; l<svd.lmax; l++) svd.MatsubaraFrequency(fUw[l], svd.fU[l], N_w);
      for (int l=0; l<svd.lmax_bose; l++) svd.MatsubaraFrequency(bUw[l], svd.fU_bose[l], N_W, "bose");

      for (int im=0; im<2*N_w; im++)
	for (int l=0; l<svd.lmax; l++)
	  ul_om(im,l ) = (im >=N_w)  ? fUw(l,im-N_w) : fUw(l,N_w-im-1).conj();
      
      for (int iOm=0; iOm<2*N_W-1; iOm++)
	for (int l=0; l<svd.lmax_bose; l++)
	  ul_Om(iOm,l) = (iOm>=N_W-1) ? bUw(l,iOm-(N_W-1)) : bUw(l,N_W-1-iOm).conj();
    }
    
  
    function2D<dcomplex> BubbleH(2*nomv,2*nomv), BubbleF(2*nOm-1,2*nomv);
    int nom = iom.size();

  
    function3D<double> VFX(svd.lmax_bose, svd.lmax, svd.lmax);

    outv<<"# i0  i1  l2  l1  l0   VF: "<<endl;
      
    (*yout)<<"Printing vertex"<<endl;
    for (int i0=0; i0<cluster.Nvfl; i0++){
      TBath ind0 = cluster.vfli_index[i0];
      int ind0t = cluster.bfl_index[ind0.ifl][ind0.bind];      
      for (int i1=0; i1<cluster.Nvfl; i1++){
	TBath ind1 = cluster.vfli_index[i1];
	int ind1t = cluster.bfl_index[ind1.ifl][ind1.bind];
	const funProxy<dcomplex>& _Gd0 = Gd[ind0t];
	const funProxy<dcomplex>& _Gd1 = Gd[ind1t];

	// BubbleH
	for (int im1=0; im1<2*nomv; im1++){
	  dcomplex gd0 = (im1>=nomv) ? _Gd0[im1-nomv] : _Gd0[nomv-im1-1].conj();
	  for (int im2=0; im2<2*nomv; im2++){
	    dcomplex gd1 = (im2>=nomv) ? _Gd1[im2-nomv] : _Gd1[nomv-im2-1].conj();
	    BubbleH(im1,im2) = gd0*gd1; // BH[iw_1,iw_2] = G[i0][iw_1] * G[i1][iw_2]
	  }
	}
      
	// BubbleF
	for (int iOm=0; iOm<2*nOm-1; iOm++){
	  // iom-iOme can be obtained by im+dOm
	  // iom = (2*(im-N_w)+1)*pi/beta  and iOme = (2*(iOm-N_W)+2)*pi/beta
	  // so that iom-iOme = pi/beta*( 2*(im+ [N_W-iOm-1] -N_w)+1  )
	  int dOm = nOm-1-iOm; 
	  int sm1 = max(0, -dOm);
	  int em1 = min(2*nomv, 2*nomv-dOm);
	  for (int im1=sm1; im1<em1; im1++){
	    dcomplex gd0 = (im1>=nomv) ? _Gd0[im1-nomv] : _Gd0[nomv-im1-1].conj();	    
	    int im2 = im1+dOm; // this is iw-iOm
	    if (im2<0 || im2>=2*nomv) continue;
	    dcomplex gd1 = (im2>=nomv) ? _Gd1[im2-nomv] : _Gd1[nomv-im2-1].conj();
	    BubbleF(iOm,im1) = gd0*gd1; // BF[iOm,iw_1] = G[i0][iw_1] * G[i1][iw_2-iOm]
	  }
	}

	/*
	if (false){
	  for (int i=0; i<Nl; i++){
	    vector<int> li = svd.get_lll(i);
	    double _VF_ = 0.0;
	    for (int j=0; j<Nl; j++){
	      vector<int> lj = svd.get_lll(j);
	      if ( (li[0] + li[1] + li[2] + lj[0] + lj[1] + lj[2])%2 ) continue; // this vanishes due to symmetry
	      _VF_ += Cll(i,j) * VH(i0,i1,lj[2],lj[1],lj[0]);
	    }
	    VFX(li[2],li[1],li[0]) = _VF_;
	    outv<<setw(3)<<i0<<" "<<setw(3)<<i1<<" "<<setw(3)<<li[2]<<" "<<setw(3)<<li[1]<<" "<<setw(3)<<li[0]<<"    "<<setw(15)<<VFX(li[2],li[1],li[0])<<endl;
	  }
	}
	*/
	  
	function3D<dcomplex> Cl; 
	Cl.resize(2*N_w,2*N_w,svd.lmax_bose);
	Cl=0.0;
	function2D<dcomplex> ul_VH(2*N_w,svd.lmax);
	function2D<dcomplex> ul_VH_ul(2*N_w,2*N_w);
	function2D<dcomplex> _VH_(VH.N3, VH.N4);
	for (int lc=0; lc<svd.lmax_bose; lc++){
	  // ul_VH(im1,l2) = sum_{l1} ul_om(im1,l1)*VH(i0,i1,lc,l1,l2);
	  // ul_VH_ul(im1,im2) = sum_{l2} ul_VH(im1,l2)*ul_om(im2,l2)
	  // Cl(im1,im2,lc) = ul_VH_ul(im1,im2) 
	  for (int i3=0; i3<VH.N3; i3++)
	    for (int i4=0; i4<VH.N4; i4++) 
	      _VH_(i3,i4) = VH(i0,i1,lc,i3,i4);
	  ul_VH.Product(   "N", "N", ul_om, _VH_ , 1.0, 0.0);
	  ul_VH_ul.Product("N", "T", ul_VH, ul_om, 1.0, 0.0);
	  for (int i=0; i<2*N_w; i++)
	    for (int j=0; j<2*N_w; j++) 
	      Cl(i,j,lc) = ul_VH_ul(i,j);
	}

	function3D<dcomplex> Dl; 
	Dl.resize(2*N_w,2*N_w,svd.lmax_bose);
	Dl=0.0;
	for (int lc=0; lc<svd.lmax_bose; lc++){
	  // ul_VH(im1,l2) = sum_{l1} ul_om(im1,l1)*VH(i0,i1,lc,l1,l2);
	  // ul_VH_ul(im1,im2) = sum_{l2} ul_VH(im1,l2)*ul_om(im2,l2)
	  // Dl(im1,im2,lc) = ul_VH_ul(im1,im2) 
	  for (int i3=0; i3<VH.N3; i3++)
	    for (int i4=0; i4<VH.N4; i4++) 
	      _VH_(i3,i4) = VFX(lc,i3,i4);
	  ul_VH.Product(   "N", "N", ul_om, _VH_ , 1.0, 0.0);
	  ul_VH_ul.Product("N", "T", ul_VH, ul_om, 1.0, 0.0);
	  for (int i=0; i<2*N_w; i++)
	    for (int j=0; j<2*N_w; j++) 
	      Dl(i,j,lc) = ul_VH_ul(i,j);
	}
      

      
	for (int iOm=0; iOm<2*nOm-1; iOm++){
	  double rOmega = 2*(iOm-nOm+1)*M_PI/common::beta;       // actual center of mass frequency
	  int dOm = nOm-1-iOm;
	  for (int im1=0; im1<2*nomv; im1++){
	    stringstream strc;
	    strc<<"Cvertex."<<i0<<"."<<i1<<"."<<(iOm-nOm+1)<<"."<<(2*(im1-nomv)+1);
	    ofstream vout(strc.str().c_str());
	    for (int im2=0; im2<2*nomv; im2++){
	      dcomplex bH = (iOm == nOm-1) ? BubbleH(im1,im2) : 0.0; // we have BH = delta(iOm)*BH[iw_1,iw_2]
	      dcomplex bF = (im1 == im2) ? BubbleF(iOm,im1) : 0.0;   // we have BF = delta(iw_1-iw_2)*BF[iOm,iw_1]
	      double om2 = (2*(im2-nomv)+1)*M_PI/common::beta;       // actual frequency

	      dcomplex VH=0, VF=0;
	      // computing the frequency dependent vertex
	      int NOm_off = N_W-nOm;
	      int nw_off = N_w-nomv;
	      dcomplex* __restrict__ _Cl_    = &Cl(im1+nw_off,im2+nw_off,0);        // Cl(om1,om2,lc=0)
	      dcomplex* __restrict__ _Dl_    = &Dl(im1+nw_off,im2+nw_off,0);        // Dl(om1,om2,lc=0)
	      dcomplex* __restrict__ _ul_Om_ = &ul_Om(iOm+NOm_off,0);               // ul_Om(Om,lc=0)
	      for (int lc=0; lc<svd.lmax_bose; lc++) VH += _Cl_[lc] * _ul_Om_[lc];                                         // Cl(om1,om2) * ul_Om(Om,lc)
	      for (int lc=0; lc<svd.lmax_bose; lc++) VF += _Dl_[lc] * _ul_Om_[lc];                                         // Cl(om1,om2) * ul_Om(Om,lc)
	    
	      vout <<om2<<" "<<VH-bH<<" "<<bH<<" "<<VF-bF<<" "<<bF<<endl;
	    }
	  }
	}
      }
    }
}


template <class Rand>
void CTQMC<Rand>::WriteVertex(const function2D<dcomplex>& Gd, bool print_vertex_xxx)
{
  // Summary of Matsubara frequencies:
  //   we use the notation : nomv == N_w; nOm == N_W
  //   The fermionic frequencies are : om = (2*(iw-N_w)+1)*pi/beta   with iw=[0,....,2*N_w-1]  hence om=[ -(2*N_w-1)*pi/beta,....,(2*N_w-1)*pi/beta]
  //   The bosonic   frequencies are : Om = (2*(iW-N_W)+2)*pi/beta   with iW=[0,....,2*N_W-2]  hence Om=[ -2*(N_W-1)*pi/beta,....,2*(N_W-1)*pi/beta]
  //
  //   Example:
  //     iom-iOme can be obtained by im+dOm, where dOm is computed below:
  //     iom = (2*(im-N_w)+1)*pi/beta  and iOme = (2*(iOm-N_W)+2)*pi/beta
  //     so that iom-iOme = pi/beta*( 2*(im+ [N_W-iOm-1] -N_w)+1  )
  function2D<dcomplex> ul_om;
  function2D<dcomplex> ul_Om;
  int N_W = max(nOm,2*nomv);
  int N_w = nomv; 
  if (common::Qsvd){
    N_w = nomv+nOm-1;  // increase it, so that we do not need to cut-off frequencies in Fock term
    ul_om.resize(2*N_w,svd.lmax);
    ul_Om.resize(2*N_W-1,svd.lmax_bose);
    {
      function2D<dcomplex> fUw(svd.lmax,N_w), bUw(svd.lmax_bose,N_W);
      for (int l=0; l<svd.lmax; l++) svd.MatsubaraFrequency(fUw[l], svd.fU[l], N_w);
      for (int l=0; l<svd.lmax_bose; l++) svd.MatsubaraFrequency(bUw[l], svd.fU_bose[l], N_W, "bose");

      for (int im=0; im<2*N_w; im++)
	for (int l=0; l<svd.lmax; l++)
	  ul_om(im,l ) = (im >=N_w)  ? fUw(l,im-N_w) : fUw(l,N_w-im-1).conj();
      
      for (int iOm=0; iOm<2*N_W-1; iOm++)
	for (int l=0; l<svd.lmax_bose; l++)
	  ul_Om(iOm,l) = (iOm>=N_W-1) ? bUw(l,iOm-(N_W-1)) : bUw(l,N_W-1-iOm).conj();
    }
  }

  function2D<dcomplex> BubbleH(2*nomv,2*nomv), BubbleF(2*nOm-1,2*nomv);
  ofstream cvout("tvertex.dat");
  cvout<<"# beta, Nvfl, nomv, nOm nom"<<endl;

  int nom = iom.size();
  cvout<<common::beta<<" "<<cluster.Nvfl<<" "<<nomv<<" "<<nOm<<" "<<nom<<endl;
  for (int i0=0; i0<cluster.Nvfl; i0++){
    TBath ind0 = cluster.vfli_index[i0];
    int ind0t = cluster.bfl_index[ind0.ifl][ind0.bind];      
    cvout<<i0<<" "<<ind0t<<endl;
  }
  cvout<<"# b0 b1 Om om1"<<endl;
    
  (*yout)<<"Printing vertex"<<endl;
  for (int i0=0; i0<cluster.Nvfl; i0++){
    TBath ind0 = cluster.vfli_index[i0];
    int ind0t = cluster.bfl_index[ind0.ifl][ind0.bind];      
    for (int i1=0; i1<cluster.Nvfl; i1++){
      TBath ind1 = cluster.vfli_index[i1];
      int ind1t = cluster.bfl_index[ind1.ifl][ind1.bind];
      const funProxy<dcomplex>& _Gd0 = Gd[ind0t];
      const funProxy<dcomplex>& _Gd1 = Gd[ind1t];

      // BubbleH
      for (int im1=0; im1<2*nomv; im1++){
	dcomplex gd0 = (im1>=nomv) ? _Gd0[im1-nomv] : _Gd0[nomv-im1-1].conj();
	for (int im2=0; im2<2*nomv; im2++){
	  dcomplex gd1 = (im2>=nomv) ? _Gd1[im2-nomv] : _Gd1[nomv-im2-1].conj();
	  BubbleH(im1,im2) = gd0*gd1; // BH[iw_1,iw_2] = G[i0][iw_1] * G[i1][iw_2]
	}
      }
      
      // BubbleF
      for (int iOm=0; iOm<2*nOm-1; iOm++){
	// iom-iOme can be obtained by im+dOm
	// iom = (2*(im-N_w)+1)*pi/beta  and iOme = (2*(iOm-N_W)+2)*pi/beta
	// so that iom-iOme = pi/beta*( 2*(im+ [N_W-iOm-1] -N_w)+1  )
	int dOm = nOm-1-iOm; 
	int sm1 = max(0, -dOm);
	int em1 = min(2*nomv, 2*nomv-dOm);
	for (int im1=sm1; im1<em1; im1++){
	  dcomplex gd0 = (im1>=nomv) ? _Gd0[im1-nomv] : _Gd0[nomv-im1-1].conj();	    
	  int im2 = im1+dOm; // this is iw-iOm
	  if (im2<0 || im2>=2*nomv) continue;
	  dcomplex gd1 = (im2>=nomv) ? _Gd1[im2-nomv] : _Gd1[nomv-im2-1].conj();
	  BubbleF(iOm,im1) = gd0*gd1; // BF[iOm,iw_1] = G[i0][iw_1] * G[i1][iw_2-iOm]
	}
      }

      function3D<dcomplex> Cl; 
      if (common::Qsvd){
	Cl.resize(2*N_w,2*N_w,svd.lmax_bose);
	Cl=0.0;
	function2D<dcomplex> ul_VH(2*N_w,svd.lmax);
	function2D<dcomplex> ul_VH_ul(2*N_w,2*N_w);
	function2D<dcomplex> _VH_(VH.N3, VH.N4);
	for (int lc=0; lc<svd.lmax_bose; lc++){
	  // ul_VH(im1,l2) = sum_{l1} ul_om(im1,l1)*VH(i0,i1,lc,l1,l2);
	  // ul_VH_ul(im1,im2) = sum_{l2} ul_VH(im1,l2)*ul_om(im2,l2)
	  // Cl(im1,im2,lc) = ul_VH_ul(im1,im2) 
	  for (int i3=0; i3<VH.N3; i3++)
	    for (int i4=0; i4<VH.N4; i4++) 
	      _VH_(i3,i4) = VH(i0,i1,lc,i3,i4);
	  ul_VH.Product(   "N", "N", ul_om, _VH_ , 1.0, 0.0);
	  ul_VH_ul.Product("N", "T", ul_VH, ul_om, 1.0, 0.0);
	  for (int i=0; i<2*N_w; i++)
	    for (int j=0; j<2*N_w; j++) 
	      Cl(i,j,lc) = ul_VH_ul(i,j);
	}
      }
      for (int iOm=0; iOm<2*nOm-1; iOm++){
	double rOmega = 2*(iOm-nOm+1)*M_PI/common::beta;       // actual center of mass frequency
	int dOm = nOm-1-iOm;

	int sm1 = 0;
	int em1 = 2*nomv;
	if (!common::Qsvd){ // Here we do not have all frequencies available
	  sm1 = max(0, -dOm);
	  em1 = min(2*nomv, 2*nomv-dOm);
	}
	for (int im1=sm1; im1<em1; im1++){
	  int sm2 = 0;
	  int em2 = 2*nomv;
	  if (!common::Qsvd){ // Here we do not have all frequencies available
	    sm2 = max(0, -dOm);
	    em2 = min(2*nomv, 2*nomv-dOm);
	  }
	  
	  ofstream vout;
	  if (print_vertex_xxx){
	    stringstream strc;
	    strc<<"vertex."<<i0<<"."<<i1<<"."<<(iOm-nOm+1)<<"."<<(2*(im1-nomv)+1);
	    vout.open(strc.str().c_str());
	  }
	  
	  cvout<<setw(3)<<i0<<" "<<setw(3)<<i1<<" "<<setw(14)<<2*(iOm-nOm+1)*M_PI/common::beta<<" "<<setw(14)<<(2*(im1-nomv)+1)*M_PI/common::beta<<endl;
	    
	  for (int im2=sm2; im2<em2; im2++){
	    dcomplex bH = (iOm == nOm-1) ? BubbleH(im1,im2) : 0.0; // we have BH = delta(iOm)*BH[iw_1,iw_2]
	    dcomplex bF = (im1 == im2) ? BubbleF(iOm,im1) : 0.0;   // we have BF = delta(iw_1-iw_2)*BF[iOm,iw_1]
	    double om2 = (2*(im2-nomv)+1)*M_PI/common::beta;       // actual frequency

	    dcomplex VH=0, VF=0;
	    if (common::Qsvd){ // computing the frequency dependent vertex
	      int NOm_off = N_W-nOm;
	      int nw_off = N_w-nomv;
	      dcomplex* __restrict__ _Cl_    = &Cl(im1+nw_off,im2+nw_off,0);        // Cl(om1,om2,lc=0)
	      dcomplex* __restrict__ _ul_Om_ = &ul_Om(iOm+NOm_off,0);               // ul_Om(Om,lc=0)
	      for (int lc=0; lc<svd.lmax_bose; lc++) VH += _Cl_[lc] * _ul_Om_[lc];                                         // Cl(om1,om2) * ul_Om(Om,lc)
	      int NW_off = N_W-1;
	      _Cl_    = &Cl(im1+nw_off,im1+dOm+nw_off,0);                           // Cl(om1,om1-Om,lc=0)
	      _ul_Om_ = &ul_Om(im1-im2+NW_off,0);                                   //  ul_Om(om1-om2,lc=0)
	      for (int lc=0; lc<svd.lmax_bose; lc++) VF += _Cl_[lc] * _ul_Om_[lc];  // Cl(om1,om1-Om) * ul_Om(om1-om2,lc)
	    }else{
	      VH = VertexH(i0,i1,iOm,im1,im2);
	      VF = VertexF(i0,i1,iOm,im1,im2);
	    }
	    if (print_vertex_xxx) vout <<om2<<" "<<VH-bH<<" "<<bH<<" "<<VF-bF<<" "<<bF<<endl;
	    cvout<<setw(10)<<om2<<" "<<setw(12)<<VH-bH<<"   "<<setw(12)<<VF-bF<<endl;
	  }
	}
      }
      
    }
  }
}

void SaveProbability(function2D<double>& AProb, double asign_fast)
{
  ofstream out("Probability.dat"); out.precision(16);
  out<<"# asign="<<asign_fast<<endl;
  for (int j=0; j<AProb.size_N(); j++)
    for (int k=0; k<AProb[j].size(); k++)
      out<<setw(3)<<j+1<<" "<<setw(3)<<k<<" "<<setw(20)<<AProb[j][k]<<endl;
}

class fDspl{
public:
  const spline1D<double>& fg;
  const mesh1D t;
  fDspl(const spline1D<double>& fg_, const mesh1D& t_) : fg(fg_), t(t_) {}
  double operator()(double x, double beta) const{
    return fg(t.Interp(x));
  }
};
double ComputeEkin_svd(const function2D<double>& Gd, const vector<function2D<spline1D<double>* > >& Deltat, const mesh1D& tau, const ClusterData& cluster, const SVDFunc& svd)
{// We compute kinetic energy by Tr(Delta*G) in SVD basis, where it takes the form sum_ g_l * d_l
  function1D<double> conjl(svd.lmax);
  for (int l=0; l<svd.lmax; l++){
    bool even = (svd(l,0)*svd(l,-1)>0.0);
    conjl[l] = even ? -1 : 1; // (-1)^{l+1}
  }
  function1D<double> Dl(svd.lmax);
  double Ekin=0;
  for (int ifl=0; ifl<cluster.N_ifl; ifl++){
    for (int i1=0; i1<cluster.ifl_dim[ifl]; i1++){
      for (int i2=0; i2<cluster.ifl_dim[ifl]; i2++){
	fDspl fD(*Deltat[ifl][i1][i2], tau);
	svd.ExpandAnalyticFunction(fD, Dl);
	int ind = cluster.tfl_index[ifl](i1,i2);
	int i_unique = cluster.bfl_index[ifl][ind];
	const funProxy<double>& gl = Gd[i_unique];
	double sgn = cluster.sign[ifl][ind];
	double dsum=0;
	/*
	cout<<"ifl="<<ifl<<" i1="<<i1<<" i2="<<i2<<" sign="<<sgn<<endl;
	if (cluster.conjg[ifl][ind]){
	  for (int l=0; l<svd.lmax; l++){
	    cout<<"    l="<<l<<" gl="<<gl[l]*conjl[l]*sgn<<" Dl="<<Dl[l]<<"  conjl="<<conjl[l]<<endl;
	  }
	}else{
	  for (int l=0; l<svd.lmax; l++){
	    cout<<"    l="<<l<<" gl="<<gl[l]*sgn<<" Dl="<<Dl[l]<<"  conjl="<<conjl[l]<<endl;
	  }
	}
	*/
	if (cluster.conjg[ifl][ind]){
	  for (int l=0; l<svd.lmax; l++) dsum += gl[l]*Dl[l]*conjl[l]*sgn;
	}else{
	  for (int l=0; l<svd.lmax; l++) dsum += gl[l]*Dl[l]*sgn;
	}
	Ekin -= dsum;
      }
    }
  }
  return Ekin;
}


int main(int argc, char *argv[])
{
  setvbuf (stdout, NULL, _IONBF, BUFSIZ);
  int my_rank, mpi_size, Master;
  MPI_Init(argc, argv, my_rank, mpi_size, Master);

  if (mpi_size>1){ // we have parallel run
    string filename = NameOfFile("nohup_imp.out",my_rank); // each processor has his own log-file
#ifdef _DEBUG
    yout = new ofstream(filename.c_str());  // in debug mode all processor print out where they are
#else
    if (my_rank==Master){                   // in non-debug mode, we print only on master
      yout = new ofstream(filename.c_str());
    }else{
      yout = new ofstream("/dev/null");    // the rest of processors are redirected to /dev/null. Otherwise it slows down the code...
    }
#endif
  }else{
    yout = &cout;
  }
  (*yout) << "Hello World! I am " << my_rank << " of " << mpi_size << std::endl;
  
  if (my_rank==Master){
    xout = new ofstream("ctqmc.log");
  }else{
    xout = &cout;
  }
  
#ifdef _LOGGING
  ofstream otrace(NameOfFile("traces_log",my_rank).c_str());  
  slog[TRACE]=&otrace;
  ofstream ostats(NameOfFile("stats_log",my_rank).c_str());  
  slog[STATS] = &ostats;
#endif
  string mode   = "";          //default "GH" where the options are [G|S][H|M], i.e., first letter can be G or S, the second can be H or M. 
  // The first letter being G meand we are sampling the Green's function, while S stand for sampling the self-energy from the two particle vertex
  // The second letter beging H means the tail is computed by Hubbard I, while M computes the high frequency tails from evaluating the high frequency moments.
  string inputf = "PARAMS";      // filename of the parameter file
  string fDelta   = "Delta.inp"; // input bath function Delta(iom)
  string fcix     = "imp.cix";   // input file with atom eigenvalues and eigenvectors (exact diagonalization results)
  string fGf      = "Gf.out";    // output file for green's function G(iom)
  string fSig     = "Sig.out";   // output file for self-energy Sig(iom)
  double mu       = 0;           // chemical potential
  double beta     = 100;         // inverse temperature
  double U        = 0;           // Coulomb repulsion (should be used only in single site DMFT)
  long long M     = 1000000;     // number of all qmc steps
  int    Ntau     = 0;           // number of imaginary time points for Delta(tau). Will be set to 5*beta+100, if not specified in the PARAMS file.
  int    Nmax     = 50;          // maximum number of kinks
  int    nom      = 0;           // number of frequency points sample in qmc
  int    nomD     = 0;           // number of frequency points for self-energy, obtained by the Dyson equation.
                                 // The rest will be computed by two particle correlation function (the direct trick)
  int    nomb     = 100;         // number of bosonic frequencies for susceptibilities
  double PChangeOrder = 0.9;     // probability to change the order (as compared to move a kink)
  double PMove = 0.3;            // Probability to move a kink, versus permute times.
  double PlocalMoves = 0.0;      // When using non-segment scheme, how much local moves should be tryed. It is advisable to use mostly non-local moves.
  int tsampleFast= 10;          // how often to measure probabilities and nf
  int    tsample  = 200;         // how often to record measurements (each tsample steps)
  int    warmup   = 0;           // how many qmc steps should be ignored
  int CleanUpdate = 100000000;   // clean update is sometimes useful
  double minF     = 1e-10;       // after applying F or F^+ and getting all matrix elements smaller than minF, we can assume that Pauli principle does not allow such step
  double minM     = 1e-15;       // trace shuld always be larger than this minimal value when accepting the step
  double minD     = 0.0;         // determinant of hybrodization should always be larger than this number when accepting the step
  int    Ncout    = 500000;      // screen output for basic information every Ncout steps
  long long  Naver= 50000000000; // much more information about simulation every Naver steps
  double TwoKinks = 0.;          // probability to add two kinks
  int    GlobalFlip= -1;         // global flip is tried after GlobalFlip qmc steps
  double treshold = 1e-10;       // energy to throw away atomic states when creating new.cix for next iteration
  int  SampleGtau = -1;          // How often to sample for Gtau (if at all)
  int SampleVertex =-1;          // How often to sample two particle vertex
  int  ReduceComm = 0;           // only part of Gtau is transfered between nodes (if slow communication)
  int    Ncorrect = -1;          // baths with index higher than Ncorrect should not be added a high-frequency tail (if -1, all tails are added)
  int         aom = 1;           // number of frequency points to find Sigma(om_last_sampled)
  int         som = 3;           // number of frequency points to find Susc(om_last_sampled)
  int    PreciseP = 1;           // computes probabilities more precisely
  double   sderiv = 0.01;        // maximum mismatch when matching high and low frequency of imaginary part of Sigmas
  double minDeltat= 1e-7;        // Delta(tau) is sometimes not causal due to numerical error. In this case we set it to small value minDeltat.
  bool SampleSusc = false;       // If spin and charge dynamic susceptibility should be sampled during simulation
  int        nomv = 50;          // number of fermionic frequencies for computing vertex
  int         nOm = 1;           // number of bosonic frequencies for vertex calculation
  double maxNoise = 1e100;        // Maximum amount of noise allowed in Green's function
  vector<double> BathProbability(100); // If bath is choosen randomly or depending on number of successful kinks in the bath
  bool  LazyTrace = true;        // This will speed up the code for Full-U calculation. It computes trace only for a few superstate, and neglects some small contributions to trace. While this seems to work very well, it is not exact.
  //double MinTraceRatio=1e-6;     // For LazyTrace we need a cutoff for how precise we want to compute trace.
  int     AllRead = 0;           // Should all processes read cix file, or just Master?
  int fastFilesystem = 0;        // Should we write-out Delta.tau and Gaver...
  int iseed_start = 0;           // starting seed
  int Segment = -1;              // Assuming superstates/blocks are all one dimensional -- Segment picture
  int Nhigh   = 40;              // The last Nhigh number of Matsubara points will be used to determine the tails
  bool SampleTransitionP = false;// Should we sample transition probability?
  bool saveStatus=1;             // Should always save status
  ifstream inp(inputf.c_str());
  int svd_lmax = 0;              // for SVD functions, the number of SVD functions used to store G(tau) data.
  double svd_L = 0;              // for SVD functions, the real axis high energy cutoff (10 by default)
  double svd_x0 = 0;             // for SVD functions, the real axis low energy cutoff (0.005 by default)
  int svd_Nw = 0;                // for SVD functions, the real axis mesh size  (500 by default)
  int svd_Ntau = 0;              // for SVD functions, the number of points on imaginary time axis (by default 4*beta+50)
  int print_vertex_xxx = 0;       // Do not print vertex.a.b.Om.om1
  if (!inp.good()) {cerr<<"Missing input file!... exiting "<<endl; return 1;}
  for (int i=0; i<BathProbability.size(); i++) BathProbability[i]=1.0;
  
  while (inp){
    string str;
    inp>>str;
    if (str[0]=='#') inp.ignore(2000,'\n');
    if (str=="mode")        inp>>mode;
    if (str=="Delta")       inp>>fDelta;
    if (str=="cix")         inp>>fcix;
    if (str=="mu")          inp>>mu;
    if (str=="beta")        inp>>beta;
    if (str=="U")           inp>>U;
    if (str=="M"){
      double dM;
      inp>>dM;
      M = static_cast<long long>(dM);
    }
    if (str=="Ntau")        inp>>Ntau;
    if (str=="Nmax")        inp>>Nmax;
    if (str=="nom")         inp>>nom;
    if (str=="iseed")       inp>>iseed_start;
    if (str=="nomD")        inp>>nomD;
    if (str=="nomv")        inp>>nomv;
    if (str=="nomb")        inp>>nomb;
    if (str=="nOm")         inp>>nOm;
    if (str=="tsample")     inp>>tsample;
    if (str=="tsampleFast")inp>>tsampleFast;
    if (str=="warmup")      inp>>warmup;
    if (str=="CleanUpdate") inp>>CleanUpdate;
    if (str=="minF")        inp>>minF;
    if (str=="minM")        inp>>minM;
    if (str=="minD")        inp>>minD;
    if (str=="PChangeOrder")inp>>PChangeOrder;
    if (str=="PMove")       inp>>PMove;
    if (str=="Ncout")       inp>>Ncout;
    if (str=="Naver")       inp>>Naver;
    if (str=="Gf")          inp>>fGf;
    if (str=="Sig")         inp>>fSig;
    if (str=="TwoKinks")    inp>>TwoKinks;
    if (str=="GlobalFlip")  inp>>GlobalFlip;
    if (str=="treshold")    inp>>treshold;
    if (str=="SampleGtau")  inp>>SampleGtau;
    if (str=="ReduceComm")  inp>>ReduceComm;
    if (str=="Ncorrect")    inp>>Ncorrect;
    if (str=="aom")         inp>>aom;
    if (str=="som")         inp>>som;
    if (str=="PreciseP")    inp>>PreciseP;
    if (str=="sderiv")      inp>>sderiv;
    if (str=="minDeltat")   inp>>minDeltat;
    if (str=="SampleSusc")  {int iSampleSusc; inp>>iSampleSusc; SampleSusc=iSampleSusc;}
    if (str=="SampleVertex")inp>>SampleVertex;
    if (str=="maxNoise")    inp>>maxNoise;
    if (str=="LazyTrace")   inp>>LazyTrace;
    if (str=="Segment")     inp>>Segment;
    //if (str=="MinTraceRatio")inp>>MinTraceRatio;
    if (str=="Nhigh")       inp>>Nhigh;
    if (str=="SampleTransitionP")  inp>>SampleTransitionP;
    if (str=="svd_lmax")    inp>>svd_lmax;
    if (str=="svd_L")       inp>>svd_L;
    if (str=="svd_x0")      inp>>svd_x0;
    if (str=="svd_Nw")      inp>>svd_Nw;
    if (str=="svd_Ntau")    inp>>svd_Ntau;
    if (str=="saveStatus")  inp>>saveStatus;
    if (str=="PlocalMoves") inp>>PlocalMoves;
    if (str=="print_vertex_xxx") inp>>print_vertex_xxx;
    if (str=="BathProbability"){
      char chr[1001];
      inp.getline(chr, 1000, '\n');
      //cout<<"BathProbability="<<chr<<endl;
      stringstream str(chr);
      int ii=0;
      while(str){
	double x;
	str>>x;
	if (!str) break;
	BathProbability[ii++]=x;
	//cout<<x<<endl;
      }
    }
    if (str=="AllRead")     inp>>AllRead;
    if (str=="fastFilesystem") inp>>fastFilesystem;
  }
  //  common::SaveStatus=true;  !?? SHOULD SET TO TRUE ALWAYS
  common::SaveStatus=saveStatus;
  if (nom==0) nom = static_cast<int>(2.0*beta);
  if (nomv>nom) nomv = nom;
  if (nOm>nom) nOm = nom;
  if (SampleVertex<=0){nomv=0; nOm=0; }
  if (Ntau==0) Ntau = 5*beta+100;
  if (svd_lmax>0){    // reasonable choise : svd_lmax = 40
    if (svd_L==0)    svd_L    = 10.0;
    if (svd_x0==0)   svd_x0   = 0.005;
    if (svd_Nw==0)   svd_Nw   = 500;
    if (svd_Ntau==0) svd_Ntau = 5*beta+100;
  }
  (*yout)<<"Input parameters are:"<<endl;
  (*yout)<<"fDelta="<<fDelta<<endl;
  (*yout)<<"fcix="<<fcix<<endl;
  (*yout)<<"mu="<<mu<<endl;
  (*yout)<<"beta="<<beta<<endl;
  (*yout)<<"U="<<U<<endl;
  (*yout)<<"M="<<M<<endl;
  (*yout)<<"Ntau="<<Ntau<<endl;
  (*yout)<<"Nmax="<<Nmax<<endl;
  (*yout)<<"nom="<<nom<<endl;
  (*yout)<<"nomD="<<nomD<<endl;
  (*yout)<<"nomv="<<nomv<<endl;
  (*yout)<<"nomb="<<nomb<<endl;
  (*yout)<<"nOm="<<nOm<<endl;
  (*yout)<<"PChangeOrder="<<PChangeOrder<<endl;
  (*yout)<<"PMove="<<PMove<<endl;
  (*yout)<<"PlocalMoves="<<PlocalMoves<<endl;
  (*yout)<<"tsample="<<tsample<<endl;
  (*yout)<<"tsampleFast="<<tsampleFast<<endl;
  (*yout)<<"warmup="<<warmup<<endl;
  (*yout)<<"CleanUpdate="<<CleanUpdate<<endl;
  (*yout)<<"minF="<<minF<<endl;
  (*yout)<<"minM="<<minM<<endl;
  (*yout)<<"minD="<<minD<<endl;
  (*yout)<<"Ncout="<<Ncout<<endl;
  (*yout)<<"Naver="<<Naver<<endl;
  (*yout)<<"Gf="<<fGf<<endl;
  (*yout)<<"Sig="<<fSig<<endl;
  (*yout)<<"TwoKinks="<<TwoKinks<<endl;
  (*yout)<<"GlobalFlip="<<GlobalFlip<<endl;
  (*yout)<<"treshold="<<treshold<<endl;
  (*yout)<<"SampleGtau="<<SampleGtau<<endl;
  (*yout)<<"ReduceComm="<<ReduceComm<<endl;
  (*yout)<<"Ncorrect="<<Ncorrect<<endl;
  (*yout)<<"aom="<<aom<<endl;
  (*yout)<<"som="<<som<<endl;
  (*yout)<<"PreciseP="<<PreciseP<<endl;
  (*yout)<<"sderiv="<<sderiv<<endl;
  (*yout)<<"minDeltat="<<minDeltat<<endl;
  (*yout)<<"SampleSusc="<<SampleSusc<<endl;
  (*yout)<<"SampleVertex="<<SampleVertex<<endl;
  (*yout)<<"maxNoise="<<maxNoise<<endl;
  (*yout)<<"AllRead="<<AllRead<<endl;
  (*yout)<<"LazyTrace="<<LazyTrace<<endl;
  //(*yout)<<"MinTraceRatio"<<MinTraceRatio<<endl;
  (*yout)<<"fastFilesystem="<<fastFilesystem<<endl;
  (*yout)<<"SampleTransitionP="<<SampleTransitionP<<endl;
  (*yout)<<"SaveStatus="<<common::SaveStatus<<endl;
  ClusterData cluster;
  ifstream incix(fcix.c_str());
  
  if (AllRead){
    cluster.ParsInput(incix, *yout);
  }else{
    if (my_rank==Master) cluster.ParsInput(incix, *yout);
    cluster.BcastClusterData(my_rank, Master);
  }

  if (Segment<0){    // Segment was not yet initialized
    if (cluster.max_size==1) Segment=1; //Segment can be used
    else Segment=0;                     // otherwise it can not be used
  }
  if (cluster.max_size>1) Segment=0; // Segment can not be used when we have matrices.
  (*yout)<<"Segment="<<Segment<<endl;
  
  if (Segment==0){
    PMove = 1.0; // It seems ExchangeTwoIntervals does not work in general.
    (*yout)<<"PMove corrected to "<<PMove<<endl;
  }
  if (Segment==0 && mode[0]=='S' && svd_lmax==0){
    (*yout)<<"ERROR: In non-segment case, and with svd_lmax==0 the S mode is not implemented. Either choose svd_lmax>0 or change S mode to G mode."<<endl;
    return 0;
  }
  if (mode==""){ // by default, mode="" and we set it either to "SM" or "GM". It seems "SM" does not work for general (non-segment) case.
    if (cluster.max_size==1) mode = "SM";
    else mode = "GM";
  }
  
  // It just read cix file, so that it knows if the tensor matrix U is available.
  if (mode[0]=='S' && !common::QHB2){
    (*yout)<<"WARNING: You have mode S to sample the self-energy, but did not provide the Coulomb interaction tensor in the cix file. Will not sample Sigma!"<<endl;
    mode[0]='G';
  }
  if (mode[1]=='M' && !common::QHB2){
    (*yout)<<"WARNING: You have mode M to add high-frequency moments, but did not provide the Coulomb interaction tensor in the cix file. Will use HB1 instead!"<<endl;
    mode[1]='H';
  }
  if (cluster.QHB1 && mode[1]!='M') mode[1]='H';
  if (mode[0]=='G') common::QHB2=false;

  (*yout)<<"mode="<<mode<<endl;

  if (!common::QHB2) nomD = nom;
  
  BathProb BathP(BathProbability, cluster);
  cluster.EvaluateMa(beta,mu,U);


  if (SampleSusc && svd_lmax){
    (*yout)<<"Currently susceptibiliy can be measured only when sampling with matsubara frequency. We will turn off svd functions! "<<endl;
    svd_lmax=0;
  }
  if (svd_lmax>0){
    (*yout)<<"sampling with projection to SVD functions is turned on"<<endl;
    (*yout)<<"svd_lmax="<<svd_lmax<<endl;
    (*yout)<<"svd_L="<<svd_L<<endl;
    (*yout)<<"svd_x0="<<svd_x0<<endl;
    (*yout)<<"svd_Nw="<<svd_Nw<<endl;
    (*yout)<<"svd_Ntau="<<svd_Ntau<<endl;
  }
  
  common::SetParameters(my_rank,mpi_size,mu,U,beta,cluster.max_size,cluster.N_flavors,cluster.N_ifl,tsample,tsampleFast,warmup,
			CleanUpdate,minF,PChangeOrder,PMove,PlocalMoves,Ncout,Naver,TwoKinks,GlobalFlip,treshold,SampleGtau,PreciseP,
			minDeltat,SampleSusc, SampleVertex, maxNoise,LazyTrace,Segment,fastFilesystem,SampleTransitionP, (svd_lmax>0));


  mesh1D iom_large;
  function2D<dcomplex> Deltaw;
  vector<pair<double,double> > ah;
  
  if (AllRead){
    ifstream input(fDelta.c_str());
    int n, m;
    if (!CheckStream(input, n, m)) { cerr<<"Something wrong in reading input Delta="<<fDelta<<endl; return 1;}
    if (!ReadDelta(cluster.N_unique_fl, input, n, m, iom_large, Deltaw, ah, common::beta,Nhigh)){ cerr<<"Something wrong in reading input Delta="<<fDelta<<endl; return 1;}
  }else{
    if (my_rank==Master){
      ifstream input(fDelta.c_str());
      int n, m;
      if (!CheckStream(input, n, m)) { cerr<<"Something wrong in reading input Delta="<<fDelta<<endl; return 1;}
      if (!ReadDelta(cluster.N_unique_fl, input, n, m, iom_large, Deltaw, ah, common::beta,Nhigh)){ cerr<<"Something wrong in reading input Delta="<<fDelta<<endl; return 1;}
    }
    BCastDelta(my_rank, Master, iom_large, Deltaw, ah);
  }
  
  if (nom>iom_large.size()) nom = iom_large.size();
  
  //int iseed = time(0)+my_rank*10;
  int iseed = my_rank*10 + iseed_start;
  //int iseed=0;
  RanGSL rand(iseed);
  (*yout)<<"iseed="<<iseed<<" in rank "<<my_rank<<endl;

  SVDFunc svd;
  if (svd_lmax>0){
    if (common::cmp_vertex)
      svd.cmp("both", common::beta, svd_lmax, svd_Ntau, svd_x0, svd_L, svd_Nw, *yout);
    else
      svd.cmp("fermi", common::beta, svd_lmax, svd_Ntau, svd_x0, svd_L, svd_Nw, *yout);
    svd.SetUpFastInterp();

    
    bool Print_SVD_Functions = common::cmp_vertex && common::Qsvd;
    if (my_rank==Master && Print_SVD_Functions){
      // Print SVD functions
      svd.Print("ulf.dat","fermi");
      svd.Print("ulb.dat","bose");
      //if (common::cmp_vertex) svd.Print("ulb.dat","bose");
      // Now we check how well is the orthogonality obeyed after reorthogonalization
      clog<<"Final overlap is smaller than "<<svd.CheckOverlap()<<endl;
    }
  }

  if (common::QHB2){
    cluster.Compute_Nmatrices2(common::U, false);
  }
  
  int cNmax = Nmax_from_StatusFiles(my_rank);
  cNmax = max(Nmax,cNmax);
  
  if (my_rank==Master) (*xout)<<"log(Zatom)="<<cluster.Zatom.exp_dbl()<<endl;
  
  CTQMC<RanGSL> ctqmc(rand, cluster, BathP, cNmax, iom_large, Deltaw, svd, ah, nom, nomv, nomb, nOm, Ntau, minM, minD, my_rank==Master);

  // Getting previous status from status.xxx files
  ctqmc.RetrieveStatus(my_rank,-1);
  
  // Main part of the code : sampling
  double nf = ctqmc.sample(M);

  // Saving the current status to status.xxx files
  if (common::SaveStatus) ctqmc.SaveStatus(my_rank);

  // Green's function is compressed to a few components only
  function2D<dcomplex> Gd;                     // needed for !Qsvd
  function2D<double> Gd_;                      // needed for Qsvd
  function1D<int> Gd_deg(cluster.N_unique_fl); // needed by both
  function2D<dcomplex> Sd;                     // needed for !Qsvd
  function2D<double> Ft;                       // needed for Qsvd
  function2D<dcomplex> Fd;
  Gd_deg=0;
  if (!common::Qsvd){
    Gd.resize(cluster.N_unique_fl,ctqmc.ioms().size());    
    Gd=0;
    for (int ifl=0; ifl<cluster.N_ifl; ifl++){
      for (int ib=0; ib<cluster.N_baths(ifl); ib++){
	//cout<<"variation="<<ctqmc.Variation[ifl][ib]<<" and maxvaraition="<<common::maxNoise<<endl;
	if (ctqmc.Variation[ifl][ib]>common::maxNoise) continue;
	int i_unique = cluster.bfl_index[ifl][ib];
	for (int im=0; im<ctqmc.ioms().size(); im++){
	  dcomplex gf = ctqmc.G_aver()[ifl](ib,im);
	  if (cluster.conjg[ifl][ib]) gf = gf.conj();
	  gf *= cluster.sign[ifl][ib];
	  Gd(i_unique,im) += gf;
	}
	Gd_deg[i_unique]++;
      }
    }
    //Fd.resize(cluster.N_unique_fl, ctqmc.ioms().size());
    if (common::QHB2){
      //Multiply_F_with_U(Fd, ctqmc.Faver, U, cluster); // computes  Fd_{a,b} = <psi_a psi_i^+ psi_j^+ psi_k> U_{i,j,k,b}
      Fd.resize(cluster.N_unique_fl, ctqmc.ioms().size());
      Fd=0.0;
      for (int ifl=0; ifl<cluster.N_ifl; ifl++){
	for (int ib=0; ib<cluster.N_baths(ifl); ib++){
	  int fl = cluster.v2fl_index[ifl][ib];
	  int i_unique = cluster.bfl_index[ifl][ib];
	  Fd[i_unique] += ctqmc.Faver[fl];
	}
      }
    }
  }else{
    
    Gd_.resize(cluster.N_unique_fl,svd.lmax);    
    Gd_=0;
    for (int l=0; l<svd.lmax; l++){
      bool even = (svd(l,0)*svd(l,-1)>0.0);
      for (int ifl=0; ifl<cluster.N_ifl; ifl++){
	for (int ib=0; ib<cluster.N_baths(ifl); ib++){
	  int i_unique = cluster.bfl_index[ifl][ib];
	  double gl = ctqmc.Gsvd[ifl](ib,l);
	  if (cluster.conjg[ifl][ib] && even) gl *= -1.; // G(t)->-G(-t), hence even coefficients get (-1)^{l+1}.
	  gl *= cluster.sign[ifl][ib];
	  Gd_(i_unique,l) += gl;
	  Gd_deg[i_unique]++;
	}
      }
    }
    for (int fl=0; fl<cluster.N_unique_fl; fl++) Gd_deg[fl]=Gd_deg[fl]/svd.lmax;

    if (common::QHB2){
      /*
      if (common::Segment){
	Multiply_F_with_U(Ft, ctqmc.Fsvd, U, cluster);
      }else{
      */
	Ft.resize(cluster.N_unique_fl,svd.lmax);
	Ft=0.0;
	for (int ifl=0; ifl<cluster.N_ifl; ifl++){
	  for (int ib=0; ib<cluster.N_baths(ifl); ib++){
	    int fl = cluster.v2fl_index[ifl][ib];
	    int i_unique = cluster.bfl_index[ifl][ib];
	    //for (int l=0; l<svd.lmax; l++) Ft[i_unique][l] += ctqmc.Fsvd1[fl][l];
	    Ft[i_unique] += ctqmc.Fsvd1[fl];
	  }
	}
	/*}*/
    }
  }
  

  // G(tau) is also compressed to few components only
  // the number of time slices is reduces to first and last 3 in case ReduceComm is set to true
  // In this case, full Gtau will not be exchaned or printed, but only 6 times will be printed
  // This is done to prevent SIGFAULT in parallel execution due to low communication pipeline
  int ntau = ctqmc.gtau()[0].size_Nd();
  function2D<double> Gtau;
  function1D<int> tau_ind;
  if (!ReduceComm || ntau<6){
    Gtau.resize(cluster.N_unique_fl,ntau);
    tau_ind.resize(ntau);
    for (int it=0; it<tau_ind.size(); it++) tau_ind[it]=it;
  }  else {
    Gtau.resize(cluster.N_unique_fl,6);
    tau_ind.resize(6);
    for (int it=0; it<3; it++) tau_ind[it] = it;
    for (int it=0; it<3; it++) tau_ind[6-it-1] = ntau-it-1;
  }
  if (common::SampleGtau>0){
    //    double dt = common::beta/ntau;
    Gtau=0;
    for (int iti=0; iti<tau_ind.size(); iti++){
      int it = tau_ind[iti];
      for (int ifl=0; ifl<cluster.N_ifl; ifl++){
	for (int ib=0; ib<cluster.N_baths(ifl); ib++){
	  double gf = ctqmc.gtau()[ifl](ib,it);
	  if (cluster.conjg[ifl][ib]) gf = -ctqmc.gtau()[ifl](ib,ntau-it-1);
	  gf *= cluster.sign[ifl][ib];
	  Gtau(cluster.bfl_index[ifl][ib],iti) += gf;
	}
      }
      for (int fl=0; fl<cluster.N_unique_fl; fl++) Gtau(fl,iti) /= cluster.fl_deg[fl];
    }
  }

  if (!common::Qsvd)
    Reduce(my_rank, Master, mpi_size, ctqmc.Histogram(), Gd, Fd, ctqmc.AverageProbability(), ctqmc.asign, ctqmc.asign_fast, ctqmc.Observables(),
	   ctqmc.k_aver(), ctqmc.Susc(), ctqmc.asuscg, cluster.DOsize, Gtau, ctqmc.VertexH, ctqmc.VertexF, Gd_deg, ctqmc.AP_transition,
	   common::cmp_vertex, common::QHB2, common::SampleSusc, common::SampleTransitionP);
  else
    ReduceS(my_rank, Master, mpi_size, ctqmc.Histogram(), Gd_, Ft, ctqmc.AverageProbability(), ctqmc.asign, ctqmc.asign_fast, ctqmc.Observables(),
	    ctqmc.k_aver(), Gtau, ctqmc.VH, Gd_deg, ctqmc.AP_transition, ctqmc.Gsvd,
	    common::cmp_vertex, common::QHB2, common::SampleTransitionP, *yout);

  if (my_rank==Master){
    (*xout)<<"Number of successfully sampled Green's functions is ";
    for (int i=0; i<Gd_deg.size(); i++) (*xout)<<Gd_deg[i]<<" ";
    (*xout)<<endl;

    if (!common::Qsvd){
      // Renormalizing Green's function
      for (int fl=0; fl<cluster.N_unique_fl; fl++) Gd[fl] *= 1.0/Gd_deg[fl];
    
      if (common::QHB2){
	for (int fl=0; fl<cluster.N_unique_fl; fl++) Fd[fl] *= 1.0/Gd_deg[fl];
	Sd.resize(cluster.N_unique_fl,ctqmc.ioms().size());
      }
    }else{
      Gd.resize(cluster.N_unique_fl, iom_large.size());
      vector<spline1D<double> > gf(cluster.N_unique_fl);
      for (int fl=0; fl<cluster.N_unique_fl; fl++){
	Gd_[fl] *= 1.0/Gd_deg[fl];
	svd.CreateSplineFromCoeff(gf[fl], Gd_[fl]);
	svd.MatsubaraFrequency(Gd[fl], gf[fl], iom_large.size());
      }
      if (common::QHB2){
	Fd.resize(cluster.N_unique_fl, iom_large.size());
	for (int fl=0; fl<cluster.N_unique_fl; fl++){
	  Ft[fl] *= 1.0/Gd_deg[fl];
	  spline1D<double> fh;
	  svd.CreateSplineFromCoeff(fh, Ft[fl]);
	  svd.MatsubaraFrequency(Fd[fl], fh, iom_large.size());
	}
	Sd.resize(cluster.N_unique_fl, iom_large.size());
      }
      {
	ofstream out("Gt.dat");
	out.precision(13);
	//vector<spline1D<double> > gf(cluster.N_unique_fl);
	for (int i=0; i<svd.tau.size(); i++){
	  out<<svd.tau[i]<<" ";
	  for (int fl=0; fl<cluster.N_unique_fl; fl++) out<<gf[fl][i]<<" ";
	  out<<endl;
	}
	out.close();
      }
    }

    // Common part for Qsvd and !Qsvd
    if (common::QHB2){
      // Finally computing self-energy by Sg = Ff * Gf^{-1}
      function1D<dcomplex> Gf(cluster.N_unique_fl), Ff(cluster.N_unique_fl), Sg(cluster.N_unique_fl);
      for (int im=0; im<Sd.size_Nd(); im++){
	for (int fl=0; fl<cluster.N_unique_fl; fl++) {
	  Gf[fl] = Gd(fl, im);
	  Ff[fl] = Fd(fl, im);
	}
	F_times_Inverse_G(Sg, Ff, Gf, cluster); // It does the following: Sg = Ff * Gf^{-1}
	for (int fl=0; fl<cluster.N_unique_fl; fl++) Sd(fl,im) = Sg[fl];
      }
    }
    
    // some printing which might be turned off at some point
    if (common::Qsvd){
      {
	ofstream ou("Gcoeff.dat");
	ou.precision(16);
	for (int l=0; l<svd.lmax; l++){
	  ou<<setw(2)<<l<<" ";
	  for (int fl=0; fl<cluster.N_unique_fl; fl++) ou<<setw(20)<<Gd_(fl,l)<<" ";
	  ou<<endl;
	}
	ou.close();
      }
      {
	ofstream ou("Gw.dat");
	ou.precision(13);
	for (int im=0; im<iom_large.size(); im++){
	  //double iom = (2*im+1)*M_PI/common::beta;
	  ou<<iom_large[im]<<" ";
	  for (int fl=0; fl<cluster.N_unique_fl; fl++) ou<<Gd(fl,im)<<"  ";
	  ou<<endl;
	}
	ou.close();
      }
      if (common::QHB2){
	ofstream ou("Sw.dat");
	ou.precision(13);
	for (int im=0; im<Sd.size_Nd(); im++){
	  double iom = (2*im+1)*M_PI/common::beta;
	  ou<<iom<<" ";
	  for (int fl=0; fl<cluster.N_unique_fl; fl++) ou<<Sd(fl,im)<<"  ";
	  ou<<endl;
	}
	ou.close();
      }
      {
	vector<function2D<double> > Aw(cluster.N_unique_fl);
	for (int fl=0; fl<cluster.N_unique_fl; fl++) svd.GiveRealAxis(Aw[fl], Gd_[fl], 2.0);
	for (int fl=0; fl<cluster.N_unique_fl; fl++){
	  ofstream fAw(NameOfFile("Aw.out",fl).c_str());
	  for (int i=0; i<svd.om.size(); i++){
	    fAw << svd.om[i]<<" ";
	    int lmax=Aw[fl].size_N();
	    for (int l=lmax-1; l>=0; --l) fAw << Aw[fl](l,i) <<" ";
	    fAw << endl;
	  }
	}
      }
    }
    
    for (int i=0; i<cluster.nsize; i++)
      for (int j=0; j<cluster.msize(i+1); j++)
	ctqmc.AProb[i][j]/=ctqmc.asign_fast;
    SaveProbability(ctqmc.AProb,ctqmc.asign_fast);
    
    function1D<double> mom;
    double TrSigmaG=0, Epot=0, Ekin=0;
    ctqmc.ComputeFinalAverage(mom, nf, TrSigmaG, Epot, Ekin);
    if (common::Qsvd){
      double Ekin_svd = ComputeEkin_svd(Gd_, ctqmc.Deltat, ctqmc.tau, cluster, svd);
      (*xout)<<" Ekins="<<Ekin_svd<<" ";
    }
    (*xout)<<endl;
    double asign=ctqmc.asign_fast/(ctqmc.Naver);
    (*xout)<<"nf="<<ctqmc.Observables()[0]<<" chiS0="<<ctqmc.Observables()[1]<<" chiD="<<ctqmc.Observables()[2]<<" <<Sz>>="<<ctqmc.Observables()[3]<<" TrSigmaG="<<TrSigmaG<<" <s>="<<asign<<endl;
    function2D<dcomplex> Gf(cluster.N_unique_fl,iom_large.size());
    function2D<dcomplex> Sd0(cluster.N_unique_fl,nom);
    function2D<dcomplex> Gf_1(cluster.N_unique_fl,iom_large.size());
    function1D<bool> nonzero(cluster.N_ifl), nonzerou(cluster.N_unique_fl);
    nonzero = true; nonzerou = true;
    for (int ifl=0; ifl<cluster.N_ifl; ifl++) // NO DOT INVERSE IF GF=0
      nonzero[ifl] = Qnonzero(Gd[cluster.bfl_index[ifl][0]]);
    for (int fl=0; fl<cluster.N_unique_fl; fl++)
      nonzerou[fl] = Qnonzero(Gd[fl]);
    
    Sd0=0;
    // Computes the high frequency part of Green's function and self-energy 
    Gf=0;
    Gf_1=0;
    for (int im=0; im<nom; im++){
      dcomplex iomega(0,iom_large[im]);
      for (int fl=0; fl<cluster.N_unique_fl; fl++) Gf(fl,im) = Gd(fl, im);
      for (int ifl=0; ifl<cluster.N_ifl; ifl++){
	if (nonzero[ifl]) // NO DOT INVERSE IF GF=0
	  Inverse_Gf(cluster.ifl_dim[ifl], cluster.N_unique_fl, cluster.bfl_index[ifl], cluster.sign[ifl], cluster.conjg[ifl], im, Gf, Gf_1);
      }
      for (int fl=0; fl<cluster.N_unique_fl; fl++) Gf_1(fl,im) /= cluster.fl_deg[fl];
      for (int fl=0; fl<cluster.N_unique_fl; fl++){
	if (nonzerou[fl])
	  Sd0(fl,im) = (iomega+mu)*cluster.Id[fl]-cluster.epsk[fl]-Gf_1(fl,im)-Deltaw(fl,im);
      }
    }
    
    function2D<double> nt(2,cluster.N_unique_fl);
    nt=0;
    double dt = common::beta/ntau;
    if (common::SampleGtau>0){
      int nti = tau_ind.size();
      for (int fl=0; fl<cluster.N_unique_fl; fl++){
	if (nonzerou[fl]){
	  nt[0][fl] = expin(0.0,make_pair(0.5*dt,Gtau[fl][0]),make_pair(1.5*dt,Gtau[fl][1]));
	  nt[1][fl] = expin(common::beta,make_pair((ntau-0.5)*dt,Gtau[fl][nti-1]),make_pair((ntau-1.5)*dt,Gtau[fl][nti-2]));
	  (*xout)<<"nt[0]["<<fl<<"]="<<nt[0][fl]<<endl;
	  (*xout)<<"nt[1]["<<fl<<"]="<<nt[1][fl]<<endl;
	}
      }
    }

    function2D<dcomplex> Sigma(cluster.N_unique_fl,iom_large.size());
    Sigma=0;
    for (int fl=0; fl<cluster.N_unique_fl; fl++){
      for (int im=0; im<min(nomD,nom); im++) Sigma(fl,im) = Sd0(fl,im); // Sigma from Dyson
      for (int im=min(nomD,nom); im<nom; im++) Sigma(fl,im) = Sd(fl,im); // Sigma by S-trick
    }
    
    // This routine releases most of the memory, so that the high frequency HB2 calculation would not run out of memory.
    ctqmc.ReleaseSomeTempMemory();
  
    cluster.Read_for_HB1(incix, mu, U, *yout);
    
    switch(mode[1]){
    case 'M':
      cluster.HB3(beta, mu, U, iom_large, ctqmc.AverageProbability(), nom, aom, Deltaw, Sigma, nonzero, nonzerou, *yout, sderiv);
      break;
    case 'H':
      cluster.HB1(beta, mu, iom_large, ctqmc.AverageProbability(), nom, aom, Deltaw, Sigma, nonzero, nonzerou, *yout, sderiv);
      break;
    default:
      ComputeMoments(cluster, iom_large, nom, nf, mom, ctqmc.k_aver(), nt, Sigma, true, Ncorrect, aom, sderiv);
    }

    Gf=0;
    Gf_1=0;
    for (int im=0; im<iom_large.size(); im++){
      dcomplex iomega(0,iom_large[im]);
      for (int fl=0; fl<cluster.N_unique_fl; fl++){
	if (nonzerou[fl])
	  Gf_1(fl,im) = (iomega+mu)*cluster.Id[fl]-cluster.epsk[fl]-Sigma(fl,im)-Deltaw(fl,im);
      }
      for (int ifl=0; ifl<cluster.N_ifl; ifl++){
	if (nonzero[ifl])
	  Inverse_Gf(cluster.ifl_dim[ifl], cluster.N_unique_fl, cluster.bfl_index[ifl], cluster.sign[ifl], cluster.conjg[ifl], im, Gf_1, Gf);
      }
      for (int fl=0; fl<cluster.N_unique_fl; fl++) Gf(fl,im) /= cluster.fl_deg[fl];
    }
    
    ofstream gout(fGf.c_str()); gout.precision(16);
    gout<<"# nf="<<nf<<" mu="<<mu<<" T="<<1/beta<<" TrSigmaG="<<TrSigmaG<<" Epot="<<Epot<<" Ekin="<<Ekin;
    gout<<" mom=[";
    for (int i=0; i<mom.size(); i++) gout<<mom[i]<<",";
    gout<<"]"<<endl;
    for (int im=0; im<iom_large.size(); im++){
      gout<<iom_large[im]<<" ";
      for (int fl=0; fl<Gf.size_N(); fl++) gout<<Gf(fl,im)<<" ";
      gout<<endl;
    }
    ofstream sout(fSig.c_str()); sout.precision(16);
    sout<<"# nf="<<nf<<" mu="<<mu<<" T="<<1/beta<<" TrSigmaG="<<TrSigmaG<<" Epot="<<Epot<<" Ekin="<<Ekin;
    sout<<" mom=[";
    for (int i=0; i<mom.size(); i++) sout<<mom[i]<<",";
    sout<<"]"<<endl;
    for (int im=0; im<iom_large.size(); im++){
      sout<<iom_large[im]<<" ";
      for (int fl=0; fl<Sigma.size_N(); fl++) sout<<Sigma(fl,im)<<" ";
      sout<<endl;
    }
    
    double sum=0;
    for (int i=0; i<ctqmc.Histogram().size(); i++) sum += ctqmc.Histogram()[i];
    for (int i=0; i<ctqmc.Histogram().size(); i++) ctqmc.Histogram()[i]/=sum;
    ofstream hout("histogram.dat");
    hout<<"# nf="<<nf<<" TrSigmaG="<<TrSigmaG<<" mom=[";
    for (int i=0; i<mom.size(); i++) hout<<mom[i]<<",";
    hout<<"]"<<endl;
    for (int i=0; i<ctqmc.Histogram().size()/2; i++)
      hout<<setw(3)<<i<<"  "<<setw(10)<<ctqmc.Histogram()[i]<<endl;

    if (common::QHB2){
      ofstream sout((fSig+"D").c_str()); sout.precision(16);
      sout<<"# nf="<<nf<<" mu="<<mu<<" T="<<1/beta<<" TrSigmaG="<<TrSigmaG<<" mom=[";
      for (int i=0; i<mom.size(); i++) sout<<mom[i]<<",";
      sout<<"]"<<endl;
      for (int im=0; im<nom; im++){
	sout<<iom_large[im]<<" ";
	for (int fl=0; fl<Sd0.size_N(); fl++) sout<<Sd0(fl,im)<<" ";
	sout<<endl;
      }
      sout.close();
      ofstream fout((fSig+"B").c_str()); fout.precision(16);
      fout<<"# nf="<<nf<<" mu="<<mu<<" T="<<1/beta<<" TrSigmaG="<<TrSigmaG<<" mom=[";
      for (int i=0; i<mom.size(); i++) fout<<mom[i]<<",";
      fout<<"]"<<endl;
      for (int im=0; im<nom; im++){
	fout<<iom_large[im]<<" ";
	for (int fl=0; fl<Sd.size_N(); fl++){
	  fout<<Sd(fl,im)<<" ";
	}
	fout<<endl;
      }
      fout.close();
    }

    if (common::SampleGtau>0){
      ofstream out("Gtau.dat"); out.precision(16);
      double dt = common::beta/ntau;
      for (int iti=0; iti<Gtau.size_Nd(); iti++){
	int it = tau_ind[iti];
	out<<setw(20)<<dt*(it+0.5)<<" ";
	for (int fl=0; fl<cluster.N_unique_fl; fl++){
	  out<<setw(25)<<Gtau[fl][iti]<<" ";
	}
	out<<endl;
      }
    }

    if (common::SampleSusc){
      ofstream out("susc.dat"); out.precision(16);
      int nsusc = ctqmc.Susc().size_N();
      int nom_s = ctqmc.iomb().size();
      int nom_l = iom_large.size();
      out<<0.0<<" "<<ctqmc.Observables()[1]<<" "<<ctqmc.Observables()[2]<<endl;
      for (int im=0; im<nom_s; im++){
	out<<ctqmc.iomb()[im]<<" ";
	for (int l=0; l<nsusc; l++) out<<ctqmc.Susc()(l,im)<<" ";
	out<<endl;
      }
      function1D<double> alpha_susc(nsusc);
      alpha_susc=0;
      for (int im=nom_s-1; im>nom_s-1-som; im--){
	for (int l=0; l<nsusc; l++) alpha_susc[l] += ctqmc.Susc()(l,im);
      }
      double om2 = sqr(ctqmc.iomb()[nom_s-1-som/2]);
      alpha_susc *= om2/som;
      for (int im=nom_s; im<nom_l; im++){
	double omega = 2*im*M_PI/common::beta;
	out<<omega<<" ";
	for (int l=0; l<nsusc; l++) out<<alpha_susc[l]/sqr(omega)<<" ";
	out<<endl;
      }
      if (cluster.DOsize>0){
	ofstream out("suscg.dat"); out.precision(16);
	int nom_s = ctqmc.iomb().size();
	for (int im=0; im<nom_s; im++){
	  out<<ctqmc.iomb()[im]<<" ";
	  for (int l=0; l<cluster.DOsize; l++) out<<ctqmc.asuscg(l,im)<<" ";
	  out<<endl;
	}
      }
    }
    if (common::SampleTransitionP){
      ofstream out("TProbability.dat"); out.precision(16);
      out<<"# asign="<<ctqmc.asign_fast<<endl;
      for (int ist=1; ist<=ctqmc.AP_transition.size_N(); ist++){
	// For convenience we sort the operators such that final states are sorted in the output
	vector<int> index(2*common::N_flavors);
	for (int op=0; op<2*common::N_flavors; op++) index[op]=op;
	CustomLess customLess(&cluster.F_i,ist);
	std::sort(index.begin(), index.end(), customLess);
	
	for (int op=0; op<2*common::N_flavors; op++){
	  int op_sorted = index[op];
	  int inew = cluster.F_i(op_sorted,ist);
	  if (inew>0)
	    out<<setw(3)<<ist<<" "<<setw(3)<<inew<<" "<<setw(20)<<ctqmc.AP_transition[ist-1][op_sorted]/ctqmc.asign_fast<<endl;
	}
      }
    }
  }

  if (common::cmp_vertex && my_rank==Master){  
    if (common::Qsvd){
      ctqmc.WriteVertexSVD(Gd, Gd_, my_rank, mpi_size, Master, print_vertex_xxx);
    }
    ctqmc.WriteVertex(Gd, print_vertex_xxx);
  }
  
  if (my_rank==Master){
    ctqmc.RecomputeCix();
    delete xout;
  }
  
  MPI_finalize();
  return 0;
}
double common::beta;
double common::U;
double common::mu;
int common::max_size;
int common::N_flavors;
double common::minF;
int common::tsample;
int common::tsampleFast;
int common::warmup;
int common::CleanUpdate;
double common::PChangeOrder;
double common::PMove;
int common::Ncout;
long long common::Naver;
const int nIntervals::c;
const int nIntervals::cd;
int common::my_rank;
int common::mpi_size;
int common::N_ifl;
double common::TwoKinks;
double common::treshold;
int common::GlobalFlip;
int common::SampleGtau;
int common::PreciseP;
double common::minDeltat;
bool common::SampleSusc;
bool common::cmp_vertex;
int common::SampleVertex;
double common::maxNoise;
int common::fastFilesystem;
bool common::LazyTrace;
bool common::QHB2;
Timer NState::t_evolve;
Timer NState::t_apply;
double common::smallest_dt;
int common::Segment;
bool common::SampleTransitionP;
bool common::SaveStatus;
bool common::Qsvd;
double common::PlocalMoves;
