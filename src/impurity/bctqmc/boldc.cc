//#pragma GCC diagnostic ignored "-Wwrite-strings"
#include <sstream>
#include <vector>
#include <deque>
#include <map>
#include <limits>
#include <algorithm>
#include <ctime>
#include <cmath>   // H1: std::isfinite
#include <iostream>
#include <fstream>
#include <string>
#include <set>
#include <blitz/array.h>
#include <numeric>
#include <regex>
#include "basics.h"
#include "svd_func.h"
#define _PCG
//#define _GSL
#include "Random2.h"
#include "wmpi.h"
MPIcls *mpi;
#include "read_cix.h"
#define FMT_HEADER_ONLY
#include "fmt/core.h"
#include "fmt/format.h"
#include "fmt/ostream.h"
#include "readinput.h"
#include "json.hpp"
std::ostream* xout;
#include "skeletons.h"
//#include "whybrid.h"
#include "wsign.h"
#include "wfile.h"
//#include "generatediags.h"
//#include "unique_diags.h"
#include "dynamic_mesh.h"
using json = nlohmann::json;
#include "Gexp.h"
#include "wbroad.h"
#include "any_axis.h"
#include "real_axis.h"
#include "frequency_loops_real_axis.h"
const int max_order=30;
//#include "wevaluate.h"

const bool DEBUG=false;
#include "imag_axis.h"
#include "topologymove.h"

#define TIMING
Array<int,1> Skeletons::Freq(max_order);

static void MergeKnownKeys(json& pars, const json& patch)
{
  for (auto it = patch.begin(); it != patch.end(); ++it){
    if (pars.contains(it.key())){
      pars[it.key()] = it.value();
    }
  }
}

template<class Axis1, typename T>
void EvaluateSingleDiag(/*Array<T,1>& Phi,*/double& Phi, Array<int,2>& ppstates, Array<int,1>& ops, Axis1& ra, 
			const Array<int,2>& diag, const Array<int,1>& dirspin,
			const Array<double,1>& times, const Cix& cix)
{
  int Norder = diag.extent(0);
  int PPs = PPconfigs2(diag,dirspin,cix);
  Phi = ra.Topology_trial_move(PPs,0,diag,dirspin,times,/*ppstates,ops,*/cout);
  ra.Accept_topology_move(0,PPs,Norder);
}

template<class Axis1, bool ImagAxis, typename T, class Axis0, class RNG>
void CurrentMC6(Array<double,1>& Phi_tot, Array<double,2>& Sm, Array<T,2>& G00, Array<T,2>& F00, Array<T,2>& F01,
		const Axis0& ra0, const Cix& cix, const Array<double,1>& xn, Array<T,2>& Gm, Array<T,2>& gm, RNG& drand, json& p,
		const TinyVector<int,3>& recalc_weight_t, bool cmpPhysG, double Phi_lowest_order, MPIcls& mpi)
{
  double weight_power = 0.0;   // extra power with beta to tune weight
  double some_factor = 1.0;    // extra factor to tune jumping
  int warmup=p["warmup"], Ncout=p["Ncout"], tsample=p["tmeassure"];
  double beta=p["beta"];
  int xNt = p["Nt"], xNd = p["Nd"];
  double rescale=p["rescale"];
  int N_steps = p["Ms"];
  int minorder=p["N_min_order"], maxorder=p["Norder"];
  double min_weight = p["min_weight"];
  bool cOCA_needed = (ImagAxis && minorder==2 && ra0.OCA_cut_sum<10);
  /*
  bool split_OCA = false;
  double OCA_cut1 = 1.0, OCA_cut2 = 1.0;
  if (minorder==2){
    if (pars["OCA_cut1"].get<double>()<1.0 or pars["OCA_cut2"].get<double>()<1.0){
      split_OCA = true;
      OCA_cut1 = pars["OCA_cut1"];
      OCA_cut2 = pars["OCA_cut2"];
    }
  }
  */
  
  double prob_exchange = p["prob_exchange"];
  double prob_dirspin = p["prob_dirspin"];
  double pko = p["prob_keep_order"];
  TinyVector<double,3> prob_move(0.5*(1-pko),pko,0.5*(1-pko));// Probabilities for [remove, swap, add]
  if (mpi.rank==mpi.master)
    (*mpi.out)<<"prob_move="<<prob_move<<endl;

  //int Nbath=Nds.extent(0), Npseudo=Nds.extent(1);
  //int Npseudo_unique = deg_pp.extent(0);
  //int Nbath_unique =  deg_bath.extent(0);
  
  ChangeHybridization<ImagAxis> ch(beta,p["cutoffw"],p["distr_width"]);
  
  //ImagAxis1new ra(ra0,xn,p,Npseudo_unique,Nbath,ibth_ind,ipsd_ind,2,Fm,mpi); 
  Axis1 ra(ra0,xn,Gm,cix,p,cmpPhysG,mpi);
  // We will spline outside in 
  //ra.SplineG(Gm,gm); // relevant only for imaginary axis
  
  //// This is for using default C++ settings in printing.
  mpi.out->unsetf(std::ios::floatfield); 
  mpi.out->precision(6);

  // currently sqrbetalist is just unity. Might want to tune....
  Array<double,1> sqrtbetalist(maxorder-minorder+1);
  for (int i = 0; i < sqrtbetalist.size(); ++i) sqrtbetalist(i) = std::pow(beta, weight_power * i);
  // we precompute the factor for accepting step. pow is slow in C++.
  double JumpF;
  if (ImagAxis)
    JumpF= cix.Nbath*beta*beta*pow(beta,weight_power)/some_factor;
  else 
    JumpF = 4*cix.Nbath*(2*p["cutoffw"].get<double>())*pow(beta,weight_power)/some_factor;
  //cout<<"JumpF="<<JumpF<<endl;
  // We need factorial for acceptance step probability. We can evaluate it up to max_order using double numbers....
  Array<double,1> factorial(maxorder+1);
  factorial(0)=1;
  for (int i=1; i<factorial.size(); ++i) factorial(i) = factorial(i-1)*i;
  // this still needs tuning.... Currently just take constant.
  Array<double,1> desired_distribution(maxorder-minorder+1);
  desired_distribution = 1.0;
  desired_distribution *= 1.0/sum(desired_distribution);
  
  // reweighting is connected to desired distribution.
  Array<double,1> weight(maxorder-minorder+1);
  weight = 1.0/weight.size();
  //
  //  Prepares splines for several quantities needed on the real axis:
  //   femx (==f(-xn)), lfemx (==log(f(-xn))), fep  (==f(om) ), lfep  (==log(f(om)))
  //   lmfx (combination of f(-x) and log(f(-x)) )
  //   Gms ( G(x) and log(Im(G(x))) so that we can compute Im(G(x))/f(-x) at large negative x)
  ra.Create_Fermi_Splines(Gm);
  ra.InitializeStorageMC();//Npseudo_unique,Nbath_unique,2); // SHOULD BE CHANGED

  string stmfr = ImagAxis ? "times" : "freq";
  
  // starting diagram of the middle order
  int Norder = (minorder+maxorder)/2;
  //int Norder=maxorder;
  Array<int,2> cdiag(Norder,2);
  for (int i=0; i<cdiag.extent(0); ++i){
    cdiag(i,0)=2*i; // setup NCA type diagram, and later find skeleton
    cdiag(i,1)=2*i+1;
  }
  Skeletons skl(cdiag);
  Array<int,1> cdirspin(Norder);
  int PPs=0;
  Array<double,1> ctimes(ra.number_of2(Norder));;
  if (not ra.ReadStatusFile(cdiag,cdirspin,ctimes,/**/"status."s+to_string(mpi.rank),*mpi.out)){
    (*mpi.out)<<"Note: proper status file status."<<mpi.rank<<" not read. Starting from random configuration." << endl;
    ra.SetRandom_time_freq(ctimes,Norder,drand);
    //
    bool Reducible = skl.IsReducible();
    //(*mpi.out)<<cdiag<<" is Reducible="<<Reducible<<endl;
    if (Reducible){
      int sgn_dlt=1;
      FindNextSkeleton(sgn_dlt, cdiag, skl, drand);
    }
    //(*mpi.out)<<cdiag<<" should be skeleton"<<endl;
    PPs=0;
    while (PPs==0){
      for (int i=0; i < Norder; ++i) cdirspin(i) = static_cast<int>(drand()*cix.Nbath*2);
      PPs = PPconfigs2(cdiag,cdirspin,cix);
    }
    if (ImagAxis)
      (*mpi.out)<<"Starting with cdiag="<<cdiag<<" dirspin="<<cdirspin<<" "<<stmfr<<"="<<ctimes<<endl;
    else
      (*mpi.out)<<"Starting with cdiag="<<cdiag<<" dirspin="<<cdirspin<<" "<<stmfr<<"="<<ctimes<<endl;
  }else{
    Norder = cdiag.extent(0);
    PPs = PPconfigs2(cdiag,cdirspin,cix);
    (*mpi.out)<<"Read cdiag="<<cdiag<<" dirspin="<<cdirspin<<" "<<stmfr<<"="<<ctimes<<endl;
  }
  if (DEBUG) cout<<stmfr<<"_0="<<ctimes<<endl;
  (*mpi.out)<<"wNbin_lmax="<<ra.wNbin_lmax<<endl;
  
  T cPhi = ra.Topology_trial_move(PPs,0,cdiag,cdirspin,ctimes,*mpi.out);
  ra.Accept_topology_move(0,PPs,Norder);
  (*mpi.out)<<"Phi0="<<cPhi<<" diag="<<cdiag<<" dirspin="<<cdirspin<<" "<<stmfr<<"="<<ctimes<<endl;
  //(*mpi.out)<<"  Gm="<<ra.cGm(Range(0,PPs-1),Range(0,2*Norder-1))<<" Dl="<<ra.cDelta(Range(0,Norder-1))<<endl;
  
  // Stval are sample pseudo self-energies. On real axis, we have wNbin and on imag axis lmax
  // dG00 and dF00 are for physical greens function.
  Array<double,1> Phi_accumulate(maxorder-minorder+1); Phi_accumulate=0;  // This Phi is stored in sampled_data from all iterations
  Array<double,1> Phi_accumulatef(maxorder-minorder+1); Phi_accumulatef=0;// This Phi is stored only when cmpPhysG is on, and is used to get correct normalization for G00 and F00
  Array<double,1> N_measure(maxorder-minorder+1); N_measure=0; // How many times we meassured
  Array<T,3> Stval(maxorder-minorder+1,cix.Npseudo_unique,ra.wNbin_lmax); Stval=0;
  //Array<double,3> Strel(maxorder-minorder+1,cix.Npseudo_unique,ra.wNbin_lmax); Strel=0;
  Array<T,3> dG00(maxorder-minorder+1,cix.Nbath_unique,ra.wNbin_lmax); dG00=0;
  Array<T,3> dF00(maxorder-minorder+1,cix.Nbath_unique,ra.wNbin_lmax); dF00=0;
  Array<T,3> dF01(maxorder-minorder+1,cix.Nbath_unique,ra.wNbin_lmax); dF01=0;
  //// Added by coxed
  Array<T,2> Stval_cOCA;
  double Phi_cOCA=0, Phi_cOCAf=0;
  Array<T,2> G00_cOCA, F00_cOCA, F01_cOCA;
  if (cOCA_needed){
    Stval_cOCA.resize(cix.Npseudo_unique,ra.wNbin_lmax);
    Stval_cOCA=0;
    G00_cOCA.resize(cix.Nbath_unique,ra.wNbin_lmax); G00_cOCA=0;
    F00_cOCA.resize(cix.Nbath_unique,ra.wNbin_lmax); F00_cOCA=0;
    F01_cOCA.resize(cix.Nbath_unique,ra.wNbin_lmax); F01_cOCA=0;
    // REMOVE
    //ra._F00_OCA.resize(cix.Nbath_unique,ra.wNbin_lmax); ra._F00_OCA=0;
    //ra._F01_OCA.resize(cix.Nbath_unique,ra.wNbin_lmax); ra._F01_OCA=0;
    //ra._F02_OCA.resize(cix.Nbath_unique,ra.wNbin_lmax); ra._F02_OCA=0;
    //ra._F03_OCA.resize(cix.Nbath_unique,ra.wNbin_lmax); ra._F03_OCA=0;
  }
  //(*mpi.out)<<"wNbin_lmax="<<ra.wNbin_lmax<<endl;
  
  // "sampled_data" from previous run contains reasonable data from previous iteration, which is a good starting point
  // It is being renormalized, and reused the data contains : Phi_accumulate,N_measure,dG00,dF00
  if (ra.Read_sampled_data(Phi_accumulate,Phi_accumulatef,N_measure,weight,Stval,dG00,dF00,dF01,
			   Phi_cOCA,Phi_cOCAf,Stval_cOCA,G00_cOCA,F00_cOCA,F01_cOCA,
			   "sampled_data", mpi)){
    double measure_before = sum(N_measure);
    double measure_this_time = N_steps/tsample;
    double rescale_factor = rescale*measure_this_time/measure_before;
    (*mpi.out)<<"Note: 'sampled_data' was read and rescaled by "<<rescale_factor;
    (*mpi.out)<<" : (measure_this_time/measure_before="<<(measure_this_time/measure_before)<<" * rescale="<<rescale<<")"<<endl;
    N_measure         *= rescale_factor;
    (*mpi.out)<<"   weight="<<weight<<endl;
    (*mpi.out)<<"   N_measure"<<N_measure<<endl;
    Phi_accumulate    *= rescale_factor; // Previous contribution is reduced, so that new contribution is more important.
    Phi_accumulatef   *= rescale_factor;
    Stval             *= rescale_factor;
    dG00              *= rescale_factor;
    dF00              *= rescale_factor;
    dF01              *= rescale_factor;
    if (cOCA_needed){
      Phi_cOCA          *= rescale_factor;
      Phi_cOCAf         *= rescale_factor;
      Stval_cOCA        *= rescale_factor;
      G00_cOCA          *= rescale_factor;
      F00_cOCA          *= rescale_factor;
      F01_cOCA          *= rescale_factor;
    }
  }else{
    mpi.Print("Note: sampled_data file was not read, starting from scratch\n");
    N_measure=0;
    Phi_accumulate=0;
    Phi_accumulatef=0;
    Stval=0;
    dG00=0;
    dF00=0;
    dF01=0;
    weight = 1.0/weight.size();
  }
  //(*mpi.out)<<"wNbin_lmax="<<ra.wNbin_lmax<<endl;
  if (not cmpPhysG){ // We don't want to polute these quantities from many steps ago.
    Phi_accumulatef=0;
    dG00=0;
    dF00=0;
    dF01=0;
  }
  Array<double,1> Nall_t(maxorder-minorder+1); Nall_t=0;
  Array<double,1> Nacc_t(maxorder-minorder+1); Nacc_t=0;
  Array<double,1> Nall_d(maxorder-minorder+1); Nall_d=0;
  Array<double,1> Nacc_d(maxorder-minorder+1); Nacc_d=0;
  
  double Q = ra.recmp_Q(Gm,cix);
  mpi.Print("Q recomputed to Q="s+to_string(Q)+" using Npseudo_unique="+to_string(cix.Npseudo_unique));
#ifdef TIMING
  double t_topot, t_topo, t_mtime, t_acc, t_rej, t_measure;
  t_topot=t_topo=t_mtime=t_acc=t_rej=t_measure=0;
#endif

  int N_d=0, N_t=0;
  int topo_update = 1;
  Range rangetf(0,ra.number_of2(Norder)-1);
  (*mpi.out)<<"Starting with diag="<<cdiag<<" dirspin="<<cdirspin<<" time_freq="<<ctimes(rangetf)<<endl;

  int tPPs=PPs, ircw=0;
  Array<int,2> tdiag = cdiag.copy();
  Array<int,1> tdirspin = cdirspin.copy();
  Array<double,1> ttimes = ctimes.copy();
  //Array<int,1> tops = ops.copy();
  //Array<int,2> tpstates;
  //Array<T,1> tPhi(PPs);
  T tPhi=0;
  /////////////////////////////////////////////////////////
  // Start of the long loop over MC steps /////////////////
  /////////////////////////////////////////////////////////

  (*mpi.out)<<"minorder="<<minorder<<" maxorder="<<maxorder<<endl;
  
  for (int itt=0; itt<N_steps; ++itt){
    if (N_d>xNd){
      topo_update = 0;
      N_d=N_t=0;
    }
    if (N_t>xNt){
      topo_update = 1;
      N_t=N_d=0;
    }
    int new_different = 0; // can be 0,1
                           // meaning 0:same diagram&time, 1:different time or diagram 2:only dirspin changed
#ifdef TIMING    
    auto t0 = get_cpu_time();
#endif

    //cout<<"itt="<<itt<<" topo_update="<<topo_update<<" N_d="<<N_d<<" N_t="<<N_t<<endl;
    if (topo_update) {// change diagram, including diagram order
      N_d += 1;
      tPPs=PPs;
      
      if (maxorder>1 and drand()>=prob_dirspin){
	topo_update=1;
	// new_different here can be 0 = (topology unchanged) or 1 = (topology changed)
	new_different = Next_valid_config_bold4(itt,tdiag,tdirspin,ttimes,tPPs,ch,
						cix,prob_move,some_factor,prob_exchange,minorder,maxorder,drand);
	if (new_different and allclose(cdiag,tdiag) and allclose(ctimes,ttimes)) topo_update = 2; // 2 means only tdirspin!=cdirspin
      }else{
	topo_update=2;
	new_different=1;
	// changing only dirspin. But we don't set ttimes. Be careful at accepting step to not set times=ttimes
	tdiag=cdiag;
	tPPs=0;
	while (tPPs==0) {
	  // Propose a new random configuration
	  for (int k=0; k < cdiag.extent(0); ++k) tdirspin(k) = static_cast<int>(drand()*cix.Nbath*2);
	  tPPs = PPconfigs2(cdiag, tdirspin, cix);
	}
      }
      if (new_different){
	int nn = tdiag.extent(0);
	Nall_d(nn-minorder) += 1;
	//dirspin2ops(tops,
	//	    tdiag, tdirspin);
	//PPstates1(tpstates,tPPs,
	//	  tdiag,tops,Fi);
	//if (tpstates.extent(0)!=tPPs){
	//  cout<<"ERROR itt="<<itt<<" tpstates.size="<<tpstates.extent(0)<<" tPPs="<<tPPs<<endl;
	//}
#ifdef TIMING    
	auto t1 = get_cpu_time();
	t_topot += t1-t0;
#endif
	if (topo_update==1)
	  // Now we evaluate all necessary quantities for topology change of a diagram
	  tPhi = ra.Topology_trial_move(tPPs,itt,tdiag,tdirspin,ttimes,*mpi.out);
	else if (topo_update==2){
	  tPhi = ra.Dirspin_trial_move(tPPs,itt,cdiag,tdirspin,cdirspin,ctimes,*mpi.out);
	}else
	  cerr<<"Should not happen! topo_update="<<topo_update<<endl;
	//if (DEBUG) (*mpi.out)<<"itt="<<itt<<" trial: tPhi="<<tPhi<<" tGm="<<ra.tGm(Range(0,tPPs-1),Range(0,2*nn-1))<<" tDl="<<ra.tDelta(Range(0,nn-1))<<" tpstates="<<tpstates<<endl;
#ifdef TIMING	
	auto t2 = get_cpu_time();
	t_topo += t2-t0;
#endif	
      }
    } else {
      N_t += 1;
      int n = cdiag.extent(0);
      Nall_t(n-minorder) += 1;
      // changing one of the frequencies
      double xi2=drand();
      int it_loop = static_cast<int>(drand()*ra.number_of(n)); // this is the time index
      if (ra.random_tw(ttimes,
		       it_loop, drand(),n,itt)){ // instead of fabs(tw_new)<ra.cutoff()
	new_different = 1;
	if (DEBUG) cout<<"itt="<<itt<<" ntimes="<<ttimes<<" tchange_ind="<<it_loop<<" times="<<ctimes<<" PPs="<<PPs<<" diag="<<cdiag<<" dirspin="<<cdirspin<<endl;
	// Now update all Green's functions and Phi for the trial step
	tPhi = ra.FrequencyTime_trial_move(itt,it_loop,cdiag,cdirspin,ttimes,ctimes,*mpi.out);
	//if (DEBUG) cout << "itt="<<itt<<" trial: dPhi="<<tPhi<<" tGm="<<ra.tGm(Range(0,PPs-1),Range(0,2*n-1))<<" tDl="<<ra.tDelta(Range(0,n-1))<<" cGm="<<ra.cGm(Range(0,PPs-1),Range(0,2*n-1))<<endl;
#ifdef TIMING
	auto t4 = get_cpu_time();
	t_mtime += t4-t0;
#endif	
      }
    }
    if (new_different){
      double ratio = abs(tPhi)/abs(cPhi);
      int n = cdiag.extent(0);
      int nn=n;
      if (topo_update==1){
	nn = tdiag.extent(0);
	if (nn!=n){
	  double cc = std::pow( JumpF, nn-n ) * factorial(n)/factorial(nn);
	  if (not ImagAxis)
	    cc *= sqr(factorial(nn-1)/factorial(n-1));
	  ratio *= cc*weight(n-minorder)/weight(nn-minorder); // here is old/new. if one of the orders is sampled too much, we can tune this ratio to make sure it is sampled evenly.
	}
      }
      double xi1=drand();
      bool accept = xi1<ratio;
      //if (DEBUG) (*mpi.out)<< "itt="<<itt<<" accept="<<accept<<" ratio="<<ratio<<" xi1="<<xi1<<" tPhi="<<tPhi<<" cPhi="<<cPhi<<" topo_update="<<topo_update<<endl;
#ifdef TIMING      
      auto t0 = get_cpu_time();
#endif      
      if (accept){ // Now we are actually changing the configuration
	if (topo_update==1){
	  Nacc_d(nn-minorder) += 1;
	  //if (PPs!=tPPs) // here we need resizing of cPhi,sign,diag,dirspin,times,ppstates. Maybe we should keep some of fixed size.
	  //cPhi.resize(tPPs);
	  if (nn!=n){
	    cdiag.resize(nn,2);
	    cdirspin.resize(nn);
	    ctimes.resize(ra.number_of2(nn));
	    //ops.resize(2*nn);
	  }
	  //if (PPs!=tPPs or nn!=n){
	  //  ppstates.resize(tPPs,2*nn);
	  //}
	  cdiag = tdiag;
	  cdirspin = tdirspin;
	  //ppstates = tpstates;
	  //ops=tops;
	  PPs=tPPs;
	  n=nn;
	  ctimes = ttimes;
	  ra.Accept_topology_move(itt,PPs,n);
	}else if (topo_update==0){
	  Nacc_t(n-minorder) += 1;
	  ctimes = ttimes;
	  ra.Accept_frequencyTime_move(itt,PPs,n);
	}else if (topo_update==2){
	  Nacc_d(n-minorder) += 1;
	  ra.Accept_dirspin_move(itt,PPs,n,tdirspin,cdirspin/*,tpstates,ppstates*/);
	  //if (PPs!=tPPs){
	    //cPhi.resize(tPPs);
	    //ppstates.resize(tPPs,2*nn);
	  //}
	  cdirspin = tdirspin;
	  //ppstates = tpstates;
	  //ops=tops;
	  PPs=tPPs;
	  //topo_update=1;
	}
	//if (DEBUG){
	//  cout<<"itt="<<itt<<" accept: cPhi="<<cPhi<<" cdiag="<<cdiag<<" cdirspin="<<cdirspin<<" ppstates="<<ppstates;
	//  cout<<" ctimes="<<ctimes<<" cGm="<<ra.cGm(Range(0,PPs-1),Range(0,2*n-1))<<endl;
	//}
	//if (cPhi.size()!=tPhi.size() or tPhi.size()!=tPPs) cout<<"ERROR cPhi.size="<<cPhi.size()<<" tPhi.size="<<tPhi.size()<<" tPPs="<<tPPs<<endl;
	cPhi = tPhi;
#ifdef TIMING	
	auto t1 = get_cpu_time();
	t_acc += t1-t0;
#endif	
      }else{// backtrack
	if (topo_update==1){
	  if (nn!=n){
	    tdiag.resize(cdiag.shape());
	    tdirspin.resize(cdirspin.size());
	    ttimes.resize(ctimes.size());
	    //tops.resize(ops.size());
	  }
	  //if (tPPs!=PPs){
	  //tPhi.resize(PPs);
	  //}
	  //if (nn!=n or tPPs!=PPs){
	  //  tpstates.resize(ppstates.shape());
	  //}
	  tdiag = cdiag; // copies data only, but truncates when not enough reserved memory
	  tdirspin = cdirspin;
	  tPPs=PPs;
	  ttimes = ctimes;
	}else if (topo_update==0)
	  ttimes = ctimes;
	else if (topo_update==2){
	  //if (tPPs!=PPs){
	    //tPhi.resize(PPs);
	    //tpstates.resize(ppstates.shape());
	  //}
	  if (tdirspin.size()!=cdirspin.size()){
	    cerr<<"at itt="<<itt<<" it should not happen tdirspin.size="<<tdirspin.size();
	    cerr<<" cdirspin.size="<<cdirspin.size()<<" cdiag="<<cdiag<<" tdiag="<<tdiag<<endl;
	  }
	  tdirspin = cdirspin;
	  tPPs=PPs;
	  //topo_update=1;
	}
	//
	n=cdiag.extent(0);
	//if (DEBUG) cout<<"itt="<<itt<<" backtracking diag="<<cdiag<<" dirspin="<<cdirspin<<" cGm="<<ra.cGm(Range(0,PPs-1),Range(0,2*n-1))<<endl;
#ifdef TIMING	
	auto t1 = get_cpu_time();
	t_rej += t1-t0;
#endif	
      }
      //if (DEBUG) (*mpi.out)<<"itt="<<itt<<" cGm="<<ra.cGm(Range(0,PPs-1),Range(0,2*n-1))<<" PPs="<<PPs<<" ppstates="<<ppstates<<endl;
      //cout<<"itt="<<itt<<" cdiag="<<cdiag<<" times="<<ctimes<<endl;
    }else{
      if (DEBUG) cout<<"itt="<<itt<<" diagram failed to change"<<endl;
    }
    
    //// For debugging, we can recompute all quantities from scratch and check their accuracy
    //ra.CheckQuantities(cPhi,itt,cdiag,cdirspin,ctimes,cix,cout);
    ///////////////////////////// print ////////////////////////////////////
    if (itt>0 and itt%Ncout==0){
      Array<double,1> smn = N_measure.copy();
      smn/=(sum(smn)+1e-200);
      //cout<<"smn="<<smn<<endl;
      fmt::print(*mpi.out, "itt={} visits=", itt);
      for (int iord=0; iord<maxorder-minorder+1; ++iord)
	fmt::print(*mpi.out, "{}:{:.4f},", iord+minorder, smn(iord));
      Array<double,1> smv = Phi_accumulate.copy();
      //cout<<"smv="<<smv<<endl;
      smv /= (sum(smv)+1e-200);
      fmt::print(*mpi.out, " aPhi=");
      for (int iord=0; iord<maxorder-minorder+1; ++iord)
	fmt::print(*mpi.out, "{}:{:.4f},", iord+minorder, smv(iord));
      //fmt::print(*mpi.out, " pp_chng={:.4f}", ra.PP_prop_changed/ra.PP_N_changed);
      /*
      if (smn(0)<0.0001){
	int n = cdiag.extent(0);
	Range rangetf(0,ra.number_of2(n)-1);
	(*mpi.out)<<"   diag="<<cdiag<<" dirspin="<<cdirspin<<" time_freq="<<ctimes(rangetf)<<endl;
      }
      */
      //(*mpi.out)<<"cd="<<cdiag<<" td="<<tdiag;
      fmt::print(*mpi.out, "\n");
      // now for combined results on master
      mpi.Reduce_sum_array(smn,N_measure);
      mpi.Reduce_sum_array(smv,Phi_accumulate);
      if (mpi.rank==mpi.master){
	smn /= (sum(smn)+1e-200);
	smv /= (sum(smv)+1e-200);
	fmt::print(cout, "itt={} visits=", itt);
	for (int iord=0; iord<maxorder-minorder+1; ++iord)
	  fmt::print(cout, "{}:{:.4f},", iord+minorder, smn(iord));
	fmt::print(cout, " aPhi=");
	for (int iord=0; iord<maxorder-minorder+1; ++iord)
	  fmt::print(cout, "{}:{:.4f},", iord+minorder, smv(iord));
	fmt::print(cout, "\n");
      }
      //(*mpi.out)<<"wNbin_lmax="<<ra.wNbin_lmax<<endl;
    }
    ///////////////////////////// meassure ////////////////////////////////////
    if (itt>warmup && itt%tsample==0){
#ifdef TIMING	
      auto t0 = get_cpu_time();
#endif      
      if (DEBUG) (*mpi.out)<<"itt="<<itt<<" Meassuring c_Phi="<<cPhi<<endl;
      ra.Meassure1(Phi_accumulate,Phi_accumulatef,Stval,dG00,dF00,dF01,
		   Phi_cOCA,Phi_cOCAf,Stval_cOCA,G00_cOCA,F00_cOCA,F01_cOCA,
		   N_measure,itt,cPhi,cdiag,cdirspin,ctimes,PPs,weight,sqrtbetalist,drand,*mpi.out);
#ifdef TIMING	      
      auto t1 = get_cpu_time();
      t_measure += t1-t0;
#endif      
    }
    if (bl::any(recalc_weight_t == TinyVector<int,3>(itt,itt,itt))){//  itt==recalc_weight_t[0] or itt==recalc_weight_t[1] or itt==recalc_weight_t[2]){
      ircw += 1;
      //mpi.Reduce_sum_array(N_measure);
      Array<double,1> N_measure_combined = N_measure.copy();
      mpi.Reduce_sum_array(N_measure_combined);
      if (mpi.rank==mpi.master){
	double Nall_meassure = sum(N_measure_combined);
	// Note: sometimes if the punishment of high order is very bad it is possible that the Markov chain never visits some of the orders, and the reweighting will be screwed.
	// so we add a safety check here. if the order is never visited, we add a very small value.
	for (int iord=0; iord<maxorder-minorder+1; ++iord){
	  double msr = (itt+0.0)*mpi.size/tsample;
	  if (N_measure_combined(iord) < 1e-3*msr) N_measure_combined(iord)=1e-3*msr;
	}
	//(*mpi.out)<<"previous weight="<<weight<<endl;
	//(*mpi.out)<<"N_measure_combined="<<N_measure_combined<<" desired_distribution="<<desired_distribution<<endl;
	weight *= N_measure_combined/desired_distribution;
	weight *= 1./sum(weight);
	if (min(weight)<min_weight){
	  for (int iord=0; iord<weight.size(); ++iord)
	    if (weight(iord)<min_weight) weight(iord) = min_weight;
	  weight *= 1./sum(weight);
	}
	//(*mpi.out)<<"weight="<<weight<<endl;
      }
      // broadcast the weight to all ranks
      mpi.Bcast_simple_array(weight);
      if (mpi.rank==mpi.master){
	N_measure_combined /= mpi.size;
	(*mpi.out) << "itt="<<itt<<" new weight="<<weight<<" N_measure="<<N_measure_combined<<endl;
      }else{
	(*mpi.out) << "itt="<<itt<<" new weight="<<weight<<endl;
      }
      N_measure *= rescale;
      Phi_accumulate *= rescale;
      Stval *= rescale;
      if (cOCA_needed){
	Phi_cOCA *= rescale;
	Phi_cOCAf *= rescale;
	Stval_cOCA *= rescale;
      }
      if (cmpPhysG){
	Phi_accumulatef *= rescale;
	dG00 *= rescale;
	dF00 *= rescale;
	dF01 *= rescale;
	if (cOCA_needed){
	  G00_cOCA *= rescale;
	  F00_cOCA *= rescale;
	  F01_cOCA *= rescale;
	}
      }
    }
  }
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Here the long MC loop is finished.
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  {// We save last data on each process status.rank
    std::string status_file = "status."s + std::to_string(mpi.rank);
    ofstream out(status_file); out.precision(16);
    out << cdiag.extent(0)<<endl;
    for (int i=0; i<cdiag.extent(0); ++i) out<<cdiag(i,0)<<" "<<cdiag(i,1)<<"  ";
    out<<endl;
    for (int i=0; i<cdirspin.extent(0); ++i) out<<cdirspin(i)<<" ";
    out<<endl;
    for (int i=0; i<ra.number_of2(cdiag.extent(0)); ++i) out<<ctimes(i)<<" ";
    out<<endl;
  }
  //
  mpi.Reduce_sum_array(Phi_accumulate);
  mpi.Reduce_sum_array(N_measure);
  mpi.Reduce_sum_array(Stval);
  if (cOCA_needed){
    mpi.Reduce_sum_array(Stval_cOCA);
    mpi.Reduce_sum_number(Phi_cOCA);
    if (cmpPhysG) mpi.Reduce_sum_number(Phi_cOCAf);
  }
  if (cmpPhysG){
    mpi.Reduce_sum_array(Phi_accumulatef);
    mpi.Reduce_sum_array(dG00);
    mpi.Reduce_sum_array(dF00);
    mpi.Reduce_sum_array(dF01);
    if (cOCA_needed){
      mpi.Reduce_sum_array(G00_cOCA);
      mpi.Reduce_sum_array(F00_cOCA);
      mpi.Reduce_sum_array(F01_cOCA);
    }
  }
  if (mpi.rank==mpi.master){
    std::string sampled_file = "sampled_data";
    ofstream out2(sampled_file); out2.precision(16);
    for (int iorder=0; iorder<=maxorder-minorder; ++iorder)
      out2<<iorder<<" "<<Phi_accumulate(iorder)/mpi.size<<" "<<Phi_accumulatef(iorder)/mpi.size<<" "<<N_measure(iorder)/mpi.size<<" "<<weight(iorder)<<endl;
    out2<<ra.wNbin_lmax<<endl;
    for (int iorder=0; iorder<Stval.extent(0); ++iorder){
      for (int atomic_state=0; atomic_state<Stval.extent(1); ++atomic_state){
	out2<<iorder<<" "<<atomic_state<<" ";
	for (int i=0; i<Stval.extent(2); ++i) 
	  out2<< to_str(Stval(iorder,atomic_state,i)/(mpi.size+0.0)) <<" ";
	out2<<endl;
      }
    }
    if (cmpPhysG){
      out2<<1<<endl;
      for (int iorder=0; iorder<dG00.extent(0); ++iorder){
	for (int ibath=0; ibath<dG00.extent(1); ++ibath){
	  out2<<iorder<<" "<<ibath<<" ";
	  for (int l=0; l<dG00.extent(2); ++l){ 
	    out2<< to_str(dG00(iorder,ibath,l)/(mpi.size+0.0))  <<" "
		<< to_str(dF00(iorder,ibath,l)/(mpi.size+0.0))  <<" "
		<< to_str(dF01(iorder,ibath,l)/(mpi.size+0.0))  <<" ";
	  }
	  out2<<endl;
	}
      }
    }else{
      out2<<0<<endl;
    }
    if (cOCA_needed){
      out2<<2<<endl;
      out2<< Phi_cOCA/mpi.size<<" "<<Phi_cOCAf/mpi.size<<endl;
      for (int atomic_state=0; atomic_state<Stval_cOCA.extent(0); ++atomic_state){
	out2<<atomic_state<<" ";
	for (int i=0; i<Stval_cOCA.extent(1); ++i) 
	  out2<< to_str(Stval_cOCA(atomic_state,i)/(mpi.size+0.0)) <<" ";
	out2<<endl;
      }
      if (cmpPhysG){
	out2<<1<<endl;
	for (int ibath=0; ibath<G00_cOCA.extent(0); ++ibath){
	  out2<<ibath<<" ";
	  for (int l=0; l<G00_cOCA.extent(1); ++l){
	    out2<< to_str(G00_cOCA(ibath,l)/(mpi.size+0.0))  <<" "
		<< to_str(F00_cOCA(ibath,l)/(mpi.size+0.0))  <<" "
		<< to_str(F01_cOCA(ibath,l)/(mpi.size+0.0))  <<" ";
	  }
	  out2<<endl;
	}
      }else{
	out2<<0<<endl;
      }
    }else{
      out2<<0<<endl;
    }
  }

  if (mpi.rank==mpi.master){
    double factor = Phi_lowest_order*Q/Phi_accumulate(0);   // normalization from NCA or OCA
    double factorf = cmpPhysG ? Phi_lowest_order*Q/Phi_accumulatef(0) : 0.0;   // normalization from NCA
    if (cOCA_needed){
      factor = Phi_lowest_order*Q/Phi_cOCA;   // normalization from NCA or OCA
      factorf = cmpPhysG ? Phi_lowest_order*Q/Phi_cOCAf : 0.0;
    }
    // H1: a non-finite normalization factor means this iteration's entire output would be
    // garbage (NaN propagates into PPSigma/Sig files and every later iteration). Die loudly
    // NOW, MPI-safely: mpi.Abort() -> MPI_Abort(MPI_COMM_WORLD,1) kills ALL ranks (same
    // idiom as the Delta-read failure path in this file).
    if (!std::isfinite(factor) or !std::isfinite(factorf)){
      std::cerr << "FATAL: non-finite normalization factor -- aborting before writing output."
                << " factor=" << factor << " factorf=" << factorf
                << " Phi_lowest_order=" << Phi_lowest_order << " Q=" << Q
                << " Phi_accumulate(0)=" << Phi_accumulate(0)
                << " Phi_accumulatef(0)=" << Phi_accumulatef(0)
                << " Phi_cOCA=" << Phi_cOCA << " Phi_cOCAf=" << Phi_cOCAf
                << " cOCA_needed=" << cOCA_needed << " cmpPhysG=" << cmpPhysG << std::endl;
      (*mpi.out) << "FATAL: non-finite normalization factor -- aborting before writing output."
                 << " factor=" << factor << " factorf=" << factorf
                 << " Phi_lowest_order=" << Phi_lowest_order << " Q=" << Q
                 << " Phi_accumulate(0)=" << Phi_accumulate(0)
                 << " Phi_accumulatef(0)=" << Phi_accumulatef(0)
                 << " Phi_cOCA=" << Phi_cOCA << " Phi_cOCAf=" << Phi_cOCAf
                 << " cOCA_needed=" << cOCA_needed << " cmpPhysG=" << cmpPhysG << std::endl;
      (*mpi.out).flush();
      mpi.Abort();
      exit(1);
    }
    double discrete_factor = ra.discrete_factor();
    Phi_accumulate *= factor/Q;
    Stval *= factor*discrete_factor;
    if (cOCA_needed){ // added by codex
      Phi_cOCA *= factor/Q;
      Stval_cOCA *= factor*discrete_factor;
    }
    (*mpi.out)<<" Phi_lowest_order="<<Phi_lowest_order<<" factor="<<factor<<" factorf="<<factorf<<endl;
    if (cmpPhysG){
      Phi_accumulatef *= factorf/Q;
      dG00  *= factorf*discrete_factor/Q;
      dF00  *= factorf*discrete_factor/Q;
      dF01  *= factorf*discrete_factor/Q;
      if (cOCA_needed){
	G00_cOCA  *= factorf*discrete_factor/Q;
	F00_cOCA  *= factorf*discrete_factor/Q;
	F01_cOCA  *= factorf*discrete_factor/Q;
      }
    }
    (*mpi.out) << " Phi=" << Phi_accumulate << ", " << Phi_accumulatef;
    if (cOCA_needed) (*mpi.out) << ", Phi_cOCA=" << Phi_cOCA << ", " << Phi_cOCAf;
    (*mpi.out)<<endl;

    for (int pp=0; pp<cix.Npseudo_unique; ++pp)
      Stval(Range::all(),pp,Range::all()) *= 1.0/cix.deg_pp(pp);
    if (cOCA_needed){ 
      for (int pp=0; pp<cix.Npseudo_unique; ++pp)
	Stval_cOCA(pp,Range::all()) *= 1.0/cix.deg_pp(pp);
      int iorder2 = 2-minorder;
      ra.PrintStval_OCA(Stval, Stval_cOCA, 2-minorder);
      /*
      if (iorder2>=0 and iorder2<Stval.extent(0)){
	
	double ref_norm = sum(abs(Stval(iorder2,Range::all(),Range::all())));
	double diff_norm = sum(abs(Stval(iorder2,Range::all(),Range::all())-Stval_cOCA));
	(*mpi.out)<<"Matrix-MC order-2 validator: Phi_exact="<<Phi_accumulate(iorder2)
		  <<" Phi_mc_ref="<<Phi_cOCA
		  <<" rel_diff(Stval)="<<(diff_norm/(ref_norm+1e-200))<<endl;

	ra.PrintStval_OCA(Stval, Stval_cOCA, iorder2);
      }
      */
    }

    if (cmpPhysG){
      Array<bool,1> bath_SC(cix.Nbath_unique); bath_SC=false;
      for (int iorb=0; iorb<cix.Nbath; ++iorb) bath_SC(cix.bath_sigind(iorb)) = cix.SC_bath(iorb);
      for (int ibath=0; ibath<cix.Nbath_unique; ++ibath){
	dG00(Range::all(),ibath,Range::all()) *= 1.0/cix.deg_bth(ibath);
	double extra_2 = bath_SC(ibath) ? 2.0 : 1.0;
	dF00(Range::all(),ibath,Range::all()) *= extra_2/cix.deg_bth(ibath);
	dF01(Range::all(),ibath,Range::all()) *= extra_2/cix.deg_bth(ibath);
	if (cOCA_needed){
	  G00_cOCA(ibath,Range::all()) *= 1.0/cix.deg_bth(ibath);
	  F00_cOCA(ibath,Range::all()) *= extra_2/cix.deg_bth(ibath);
	  F01_cOCA(ibath,Range::all()) *= extra_2/cix.deg_bth(ibath);
	}
      }
      if (cOCA_needed) ra.PrintG00_OCA(G00_cOCA, F00_cOCA, F01_cOCA, dG00, dF00, dF01, 2-minorder);
    }

    double broad=p["broad"];
    if (broad>0)
      ra.Broaden_Quantities(Stval,dG00,dF00,dF01,cmpPhysG,/**/p,mpi);

    ra.InterpolateOnRightMesh(Sm,G00,F00,F01,/**/Stval,dG00,dF00,dF01,Stval_cOCA,G00_cOCA,F00_cOCA,F01_cOCA,cix,cmpPhysG,p,mpi);
  }
  double _Nacc_t=sum(Nacc_t), _Nall_t=sum(Nall_t), _Nacc_d=sum(Nacc_d), _Nall_d=sum(Nall_d);
  (*mpi.out) << "acceptance rates: time/frequency="<<(_Nacc_t/_Nall_t)<<" topology="<<(_Nacc_d/_Nall_d)<<" topology/time="<<(_Nall_d/_Nall_t)<<endl;
  {
    Array<double,1> tmp(Nall_t.size());
    tmp = Nall_t/_Nall_t;
    (*mpi.out) << " freq-trial  rate per order="<<tmp<<endl;
    tmp = Nacc_t/_Nacc_t;
    (*mpi.out)<<  " freq-accept rate per order="<<tmp<<endl;
    tmp = Nall_d/_Nall_d;
    (*mpi.out) << " topology-trial  rate per-order="<<tmp<<endl;
    tmp = Nacc_d/_Nacc_d;
    (*mpi.out) << " topology-accept rate per order="<< tmp<<endl;
  }
#ifdef TIMING	
  fmt::print(*mpi.out,"t_topo={:.4f}s (t_topo_search={:.4f},t_topo_eval={:.4f})\n", t_topo,t_topot,t_topo-t_topot);
  fmt::print(*mpi.out,"t_move={:.4f}s t_acc={:.4f}s t_rej={:.4f}s t_meas={:.4f}s\n", t_mtime, t_acc, t_rej, t_measure);
#endif
  ra.PrintTimings(*mpi.out);
  
  mpi.Bcast_simple_array(Phi_accumulate);
  mpi.Bcast_simple_array(Phi_accumulatef);
  Phi_tot.resize(Phi_accumulate.size());
  Phi_tot = Phi_accumulate;
  if (mpi.rank!=mpi.master)
    ra.Resize_Output(Sm,G00,F00,F01,cix,cmpPhysG);
  (*mpi.out)<<"Sm.size="<<Sm.size()<<" G00.size="<<G00.size()<<" F00.size="<<F00.size()<<" F01.size="<<F01.size()<<endl;
  mpi.Bcast_simple_array(Sm);
  if (cmpPhysG){
    mpi.Bcast_simple_array(G00);
    mpi.Bcast_simple_array(F00);
    mpi.Bcast_simple_array(F01);
  }
}

void RereadParameters(json& pars){
  // Reread PARAMS file in case something was updated
  ifstream inp("PARAMS"); // requires file "PARAMS"
  stringstream inp2;
  RemoveComments(inp2, inp); // we want to use comments in the input PARAMS
  if (!(inp2.str().empty())){
    json patch = json::parse(inp2);
    pars.merge_patch(patch);
  }
}


template<class Axis0, class Axis1, typename T, bool ImagAxis, class RNG>
void ComputeNew(json& pars, RNG& drand, Cix& cix, MPIcls& mpi){
  double beta = pars["beta"];
  int N_max_order=pars["Norder"], N_min_order=pars["N_min_order"];
  {
    int two_2_i = 1;
    for (int i=0; i<Skeletons::Freq.extent(0); ++i){
      Skeletons::Freq(i) = two_2_i;
      two_2_i *= 2;
    }
  }
  ///////////////////////////////////////////////////////////////

  Axis0 ra(beta,cix,pars,mpi);
  
  Array<double,1> om;
  Array<complex<double>,2> Deltaw;
  try {
    auto result = ReadDelta_mpi(pars["Delta"], cix.Nbath_unique_unique, mpi);
    om.resize(result.om.shape());
    om = result.om;
    Deltaw.resize(result.Deltaw.shape());
    Deltaw = result.Deltaw;
    ra.PrepareDelta(Deltaw,om,pars,mpi);
  } catch (const std::exception& e) {
    if (mpi.rank == mpi.master) {
      std::cerr << "Error while reading Delta file: " << e.what() << std::endl;
    }
    // In MPI code you usually want to abort all ranks explicitly:
    mpi.Abort();
    exit(0);
  }
  //auto [om,Deltaw] = ReadDelta_mpi(pars["Delta"],cix.Ebath.size(), mpi);
  
  string Sfilename = pars["PPSig"].get<string>()+to_string(0);
  auto [xn,Sm,Sinitialized] = ra.ReadPseudoSig00(Sfilename,mpi);
  if (!Sinitialized)
    ra.InitializeSig00(xn,Sm,om,cix,pars["fastFilesystem"],mpi);
  
  //cout<<"1) xn(0)="<<xn(0)<<endl;
  
  Array<T,2> Gm(cix.Npseudo_unique,xn.size()); Gm=0;// Pseudo spectral functions
  Array<T,2> gm(cix.Nblock,xn.size()); gm=0; // largest eigenvalues in each block, something it's not used at the moment, but could be used later...
  Array<T,2> G00, dG00; // local Green's function
  Array<T,2> F00, F01, dF00, dF01; // local two-particle function
  Array<complex<double>,2> Sigm(cix.Nbath_unique,om.size()); // physical self-energy
  
  Array<T,2> Sm_new(cix.Npseudo_unique,xn.size());
  Array<double,2> Smt(cix.Npseudo_unique,xn.size()); // The correction to pseudo self-energy
  Array<complex<double>,2> Gc;
  Array<double,1> density_matrix; Array<double,1> density_matrix0(cix.Nbath_unique_unique); double Ekin=0;
  Array<double,1> Prob(cix.Npseudo_unique); Prob=0;
  double nf=0; double mu=pars["mu"];
  
  ra.InitializeStorage(cix.Npseudo_unique,xn);
  
  int itt=0;
  double mix=pars["mix"];

  Array<double,1> Phi_n(max_order+1);
  
  //std::vector<std::vector<int>> blocks;   // blocks[b] = list of orbitals in block b
  //std::vector<bool> is_block_sc;          // is_sc_block[b] = true if block b has SC (Nambu) terms

  BathBlocks bb = cix.FindBlocks();
  (*mpi.out) << bb << endl;
  
  double t_nca=0, t_oca=0, t_mc=0;
  bool cmpPhysG=false;
  bool next_break= (1==pars["Nitt"].get<int>()); // if Nitt==1, we break after first iteration
  if (next_break)
    (*mpi.out)<<itt<<") This should be the last iteration"<<endl;
  double lambda0_fixed=0;
  int followPeak_fixed=pars["followPeak"].get<int>();
  if (pars["lambda0_fixed"].get<double>()!=0) lambda0_fixed = pars["lambda0_fixed"];
  int include_oca=0; bool cOCA_needed=false;
  Array<double,1> Phi_tot(N_max_order); Phi_tot=0;
  for (itt=0; itt<pars["Nitt"]; itt++){
    RereadParameters(pars);
    if (next_break) cmpPhysG=true;
    
    auto [nd,Q,lambda0,Fimp,Epot,TrSigmaG,followPeak,density_matrix0,Prob]=ra.DeterminePseudos(Gm,gm,Sm,xn,
											  cix,lambda0_fixed,followPeak_fixed,itt,pars,mpi);
    if (mpi.rank==mpi.master){
      int ist = next_break ? 0 : itt+1;
      ra.Print_Pseudo(pars["PPG"].get<string>()+to_string(ist), Gm, pars);
    }
    auto _tst_ = get_cpu_time();

    Phi_tot(0) = ra.NCA(Smt, G00, F00, F01,
			cmpPhysG, Q, xn, Gm, cix, mpi);

    auto _tst1_ = get_cpu_time();
    t_nca += _tst1_ - _tst_;
    //if (mpi.rank==mpi.master) ra.Print_Pseudo("Smt_NCA."+to_string(itt+1), Smt, pars);

    if (mpi.rank==mpi.master){
      ra.Print_Pseudo("aStval.1", Smt, pars);
      ra.Print_Physical("ag00.1", G00, F00, F01);
    }
    Sm_new.resize(cix.Npseudo_unique,xn.size()); Sm_new=0;
    ra.Update_Sm(Sm_new,Smt,xn,cix,mpi);

    if (N_max_order>N_min_order or (N_min_order>=2 and N_max_order>=2)) ra.SplineG(xn,Gm,gm); // relevant only for imaginary axis
    
    if (N_min_order>=2 and N_max_order>=2){
      std::array<double,2> OCA_cut({pars["OCA_cut1"].get<double>(),pars["OCA_cut2"].get<double>()});
      if ((OCA_cut[0]-1.0)<-1e-5 or (OCA_cut[1]-1)<-1e-5){
	// Here we compute G_importance, which will allow us to speed up analytic OCA calculation
	ra.G_importance = ra.ImportanceOfPseudo(xn,Gm,OCA_cut, false);
	ra.OCA_cut_sum = pars["OCA_cut_sum"].get<int>();
	if (mpi.rank==mpi.master){
	  (*mpi.out)<<"G_importance="<< ra.G_importance << endl;
	  (*mpi.out)<<"OCA_cut_sum="<< ra.OCA_cut_sum << endl;
	}
	cOCA_needed = ImagAxis && N_min_order==2 && (N_max_order>2 or pars["OCA_MC"].get<bool>()); // only imaginary-axis minorder==2 uses cOCA semantically for now
      }else{
	ra.OCA_cut_sum = 10;
	cOCA_needed = false;
      }
      auto _tst_ = get_cpu_time();
      Phi_tot(1) = ra.OCA(Smt,dG00,dF00,dF01,
			  cmpPhysG, Q, xn, Gm, cix, pars, mpi);
      auto _tst1_ = get_cpu_time();
      t_oca += _tst1_ - _tst_;
      if (mpi.rank==mpi.master){
	ra.Print_Pseudo("aStval.2", Smt, pars);
	if (cmpPhysG) ra.Print_Physical("ag00.2", dG00, dF00, dF01);
      }
      ra.Update_Sm(Sm_new,Smt,xn,cix,mpi);
      if (cmpPhysG){
	G00 += dG00;
	F00 += dF00;
	F01 += dF01;
      }
      include_oca=1;
    }
    
    (*mpi.out)<<setw(3)<<itt<<") Phi="<<Phi_tot(Range(0,include_oca))<<"  t(NCA)="<<t_nca<<" t(OCA)="<<t_oca<<endl;mpi.out->flush();
    
    if (N_max_order>N_min_order or (N_max_order==2 and cOCA_needed)){ // even when N_max_order==2 but OCA is not complete
      TinyVector<int,3> recalc_weight_t;
      if (pars["recalc_weight"].get<double>()>0){
	double recalc_weight = pars["recalc_weight"];
	double N_steps = pars["Ms"];
	for (int i=1; i<=recalc_weight_t.length(); ++i) recalc_weight_t[i-1] = max(500000*i, int(N_steps*recalc_weight*i));
	(*mpi.out)<<"recalc_weight="<<recalc_weight_t<<endl;
      }else
	recalc_weight_t=pars["Ms"].get<double>()*pars["multiples_last_step"].get<double>()+1;

      if (next_break){pars["Ms"] = pars["Ms"].get<double>()*pars["multiples_last_step"].get<double>();}

      Array<double,1> Phi_n;
      auto _tst_ = get_cpu_time();
      CurrentMC6<Axis1,ImagAxis>(Phi_n,Smt,dG00,dF00,dF01,
				 ra, cix, xn,Gm,gm,drand, pars, recalc_weight_t, cmpPhysG, Phi_tot(N_min_order-1), mpi);
      
      auto _tst1_ = get_cpu_time();
      t_mc += _tst1_ - _tst_;
      Phi_tot(Range(N_min_order-1,N_max_order-1)) = Phi_n;
      ra.Update_Sm(Sm_new,Smt,xn,cix,mpi);
      if (cmpPhysG){
	G00+=dG00;
	F00+=dF00;
	F01+=dF01;
      }
    }

    if (cmpPhysG) ra.PhysicalSelfEnergy(Sigm, Gc, density_matrix, Ekin, nf,
					G00, F00, F01, cix, bb, pars, itt+1, mpi, "Sdyson.");
    
    double diff = ra.Difference(Sm_new,Sm,Prob);   // sum(abs(Sm_new-Sm))/Sm_new.size();
    ra.StopAdaptingMesh(diff<pars["diff_adapt_mesh"].get<double>());
    if (diff<pars["diff_adapt_mesh"].get<double>()) lambda0_fixed=lambda0; // this fixes lambda0 as well
    if (diff<pars["diff_followPeak"].get<double>()) followPeak_fixed=followPeak;
    Sm = Sm*(1.0-mix) + mix*Sm_new;
    
    (*mpi.out)<<setw(3)<<itt<<") Phi=["<<Phi_tot<<"]"<<endl;
    fmt::print(*mpi.out, "{:3d}) Difference between steps: {:10.6f}\n", itt, diff);

    
    if (mpi.rank==mpi.master){
      int ist = next_break ? 0 : itt+1;
      ra.Print_Pseudo(pars["PPSig"].get<string>()+to_string(ist), Sm, pars);
      if (cmpPhysG){
	string ending = (not next_break) ? string(".")+to_string(ist) : "";
	ra.Print_Physical2(pars["fSig"].get<string>()+ending, Sigm, density_matrix0, nd, mu, Fimp, Epot, TrSigmaG, Ekin, nf, density_matrix, cix.deg_bath);
	ra.Print_Physical2(pars["fGf"].get<string>()+ending, Gc, density_matrix0, nd, mu, Fimp, Epot, TrSigmaG, Ekin, nf, density_matrix, cix.deg_bath);
      }
    }
    
    if (next_break) break;
    if (diff<pars["diff_stop"].get<double>() or itt+2==pars["Nitt"].get<int>()){
      (*mpi.out)<<itt<<") Next should be the last iteration"<<endl;
      next_break=true;
    }
    //return ;

    (*mpi.out).flush();
  }
}

//#include <iostream>
//#include <sstream>
//#include <regex>
//#include <cmath>
//#include <string>
//#include <nlohmann/json.hpp>

void ReportJsonParseError(const nlohmann::json::parse_error& e, const std::string& txt_params, MPIcls* mpi)
{
  // Print the raw error message
  std::cerr << "JSON parse error: " << e.what() << "\n";
  // Try to extract line/column from the error text
  // Typically looks like "... parse error at line 4, column 5 ..."
  std::string msg = e.what();
  std::regex re(R"(line\s+(\d+),\s+column\s+(\d+))");
  std::istringstream iss(txt_params);
  std::smatch match;
  if (std::regex_search(msg, match, re)){
    int line_num = std::stoi(match[1]);
    int col_num  = std::stoi(match[2]);
    std::cerr << "Error at line " << line_num << ", column " << col_num << ":\n";
    // Print exactly that line from the input
    std::string line;
    for (int i = 1; std::getline(iss, line); ++i){
      //cout<<"i="<<i<<" line="<<line<<endl;
      if (i>line_num-2) std::cerr << i<<": "<< line << "\n";
      // Optionally print an indicator pointing to the column
      if (i==line_num){
	if (col_num > 0 && col_num <= (int)line.size()){
	  // Spaces up to col_num - 1, then a caret
	  std::cerr << std::string(col_num+2+std::floor(std::log10(line_num)), ' ') << "^\n";
	}
	break;
      }
    }
  }
  delete mpi;
  exit(1);
}

int main(int argc, char *argv[])
{
  //ios::sync_with_stdio(false);
  //cin.tie(nullptr);
  
  setvbuf (stdout, NULL, _IONBF, BUFSIZ);
  
  mpi = new MPIcls(argc, argv);
  
  ifstream inp("PARAMS"); // requires file "PARAMS"
  stringstream inp2;
  RemoveComments(inp2, inp); // we want to use comments in the input PARAMS
  
  json pars, pars_r, pars_i;
  pars["U"]               = 0.0;  // Coulomb repulsion, required from input
  pars["mu"]              = 0.0;  // chemical potential, required from input
  pars["beta"]            = 0.0;  // inverse temperature, required from input
  pars["Norder"]          = 5;    // at which order do we stop with MC
  pars["N_min_order"]     = 1;    // at which order do we start with MC
  pars["Ms"]              = 1e8;
  pars["cix"]             = "bctqmc.cix";
  pars["Nitt"]            = 100;  // max number of self-consistent steps for bold calculation
  pars["recalc_weight"]   = 0;    // we recalculate weight for equal distribution per order after recalc_weight*Ms MC steps have passed.
  pars["Nt"]              = 150;  // number of MC steps changing time (instead of changing diagram)
  pars["Nd"]              = 1;    // number of MC steps changing diagram&topology&order instead of just time
  pars["prob_keep_order"] = 0.75; // when topological update, how often to keep order versus jump order
  pars["prob_dirspin"]    = 0.25; // when topological update, how often to only change dirspin, but keep diagram topology&times
  pars["prob_exchange"]   = 0.25; // when topological update but not dirspin change, and don't change order, how often to exchange two operators
  pars["tmeassure"]       = 99;   // how often to meassure (it used to be called tsample)
  pars["diff_stop"]       = 0.001;// difference between successive steps in pseudo self-energy that causes the algorithm to finish
  pars["diff_adapt_mesh"] = 0.01; // difference between successive steps in pseudo self-energy that stops adapting mesh
  pars["diff_followPeak"] = 0.1;  // difference when we freeze into which peak to follow to determine lambda0
  pars["followPeak"]= -1;         // if we want to fix in which sector to compute lambda0
  pars["Nbath"]           = 1;    // magnetic vs paramagnetic calculation
  pars["mix"]             = 0.3;  // mixing for pseudo-self energy
  pars["iseed"]           = 0;    // starting iseed for first processor. 0 is special.
  pars["multiples_last_step"]=3;  // if this is last step and computes observables, how much more sampling should we do?
  pars["rescale"]         = 0.25; // observables are rescaled once they are read from previous step.
  pars["Ncout"]    = 1000000;
  pars["warmup"]       = 1;    // minimal number of steps for warmup
  pars["cmpPhysG"] = true;
  pars["Delta"]    = "Delta.inp"; // input file for hybridization
  pars["fSig"]     = "Sig.out";   // output file for self-energy Sig(iom)
  pars["PPSig"]    = "PPSigma.";
  pars["fGf"]      = "Gf.out";    // output file for green's function G(iom)
  pars["PPG"]      = "PPG.";
  pars["fastFilesystem"]=false;
  pars["converge_OCA"]=false;
  pars["subExprElim"]=true;
  pars["min_weight"]=1e-2;         // when weight for any order becomes smaller than that, dont reduce it further
  pars["lambda0_fixed"] = 0.0;     // means lambda0 should be recomputed
  pars["peak_pos"]      = 0.0;     // shift pseudo-particle threshold by peak_pos*T when determining lambda0
  pars["cutoffw"]      = 30.;      // real_axis_only: cutoff frequency on the real axis, the largest energy scale. Should not be in pars_r, because needed above
  pars["distr_width"]  = 1.5;      // real_axis_only: when adding new frequency, what is the width of the distribution. Should not be in pars_r, because needed above
  pars["broad"]        = 2.0;      // real_axis_only: how much to broad results compared to estimated error. Should not be in pars_r, because needed above
  pars["OCA_cut1"] = 1.0;
  pars["OCA_cut2"] = 1.0;
  pars["OCA_cut_sum"] = 3;
  pars["OCA_MC"] = true;          // When analytic OCA is the last order, should we add MC corrections or not.
  pars["Sdyson_Bulla_blend"] = true;  // if true, Sig.out is Sdyson at low energy smoothly blended to Bulla at high energy.
  // added by codex
  //pars["matrix_mc_samples"] = 0;
  //pars["matrix_mc_burnin"] = 32;
  // real axis calculation
  pars_r["Wtypical"]     = 1.0;     // typical frequency scale
  pars_r["wNbin"]        = 400;     // # bins for observables on real axis
  pars_r["gauss_x_width"]= 0.5;     // step in changing frequency
  pars_r["wNg"]=100;
  pars_r["x0"]=0.05;
  // adjustable mesh on the real axis
  pars_r["AdaptMesh"] = true;    // should we dynamically adapt mesh on the real axis
  pars_r["AdaptPhysicalMesh"] = true; // should we dynamically adapt mesh for physical greens function on the real axis
  pars_r["Nw_1"]      = 150;  // number of frequency points in the linear mesh at high-frequencies [-infty,Nw_La]&[Nw_Lb,infty]
  pars_r["Np_1"]      = 150;  // number of frequency points in the linear mesh at high-frequencies [-infty,Nw_La]&[Nw_Lb,infty]
  pars_r["Nw_2"]      = 150;  // number of frequency points in the adjustable mesh at low frequency between [Nw_La,Nw_Lb]
  pars_r["Np_2"]      = 150;  // number of frequency points in the adjustable mesh at low frequency between [Nw_La,Nw_Lb]
  pars_r["Nw_max"]    = 1000; // maximum number of frequency points after checking precision of the integral
  pars_r["Np_max"]    = 1000; // maximum number of frequency points after checking precision of the integral
  pars_r["Nw_itt"]    = 20;   // how many recursive steps in adding frequency points at the midpoint before the integral reaches desired precision
  pars_r["Np_itt"]    = 20;   // how many recursive steps in adding frequency points at the midpoint before the integral reaches desired precision
  pars_r["Nw_La"]     = -1;   // start of adjustable mesh for PP
  pars_r["Np_La"]     = -2;   // start of adjustable mesh for physical G
  pars_r["Nw_Lb"]     =  1;   // end of adjustable mesh for PP
  pars_r["Np_Lb"]     =  2;   // end of adjustable mesh for physical G
  pars_r["Nw_small"]  = 1e-5; // how precise should the integral be to stop adding points in the mid-point for PP G
  pars_r["Np_small"]  = 1e-5; // how precise should the integral be to stop adding points in the mid-point for physical G
  pars_r["Nw_alpha"]  = 0.4;  // for adaptive mesg to transform spectra for mesh we exponentiate with alpha, i.e., pow(A(w)*dw,alpha)
  pars_r["brd_exclude"] = 1.0;// removing noise for frequency above that energy
  pars_r["brd_max_der"] = 10.0;// removing when derivative due to noise is bigger than that 
  pars_r["brd_beta"] = 2.0;   // artificial cutoff beta for fermi function which cuts-off high frequency noise
  pars_r["Sdyson_Bulla_L_low"] = 5.0;   // real-axis energy below which Sdyson is always used; above L-L_low Bulla is used.
  // imaginary axis calculation
  pars_i["diff_adapt_mesh"]=0; // so that we don't adapt mesh or fix lambda0
  pars_i["Ntau"]     = 512;
  pars_i["svd_x0"]   = 0.01; // low energy cutoff for Kernel in SVD
  pars_i["svd_lmax"] = 20;   // number of SVD functions to use in the expansion
  pars_i["svd_L"]    = 20.;  // high energy cutoff for Kernel in SVD 
  pars_i["svd_Nw"]   = 500;  // number of points in SVD
  pars_i["broad"]    = 0.0;  // no broadening on imag axis
  pars_i["S_sigma_pts"]=7;   // gaussian width for attenuation ratio
  pars_i["S_ratio"] = 1.5;   // how much larger can smaller order be before we cut it off. If |G[iorder]|>|G[iorder-1]|*S_ratio, we start to cut
  pars_i["S_kappa"] = 4.0;   // how aggresively we cut once the treshold for cutting is achieved
  pars_i["S_min_order"]=2;   // order at which we start to correct with attenuation
  pars_i["r_dyson_laplace"]=2;// how many times larger uniform mesh we built (2^r_dyson_laplace points) for solving Dyson-Laplace
  pars_i["volterra_scheduler_auto"] = true;          // if true, choose scheduler thresholds from mpi.size automatically; if false, use the manual values below verbatim
  pars_i["volterra_cmp2_small_cutoff"] = 4;          // block-size cutoff for local solver choice: N<=cutoff uses cmp2_small, larger blocks use cmp2_large
  pars_i["volterra_no_inside_parallel_below"] = 64;  // blocks smaller than this are never parallelized inside the block; they stay block-parallel over ib
  pars_i["volterra_prefer_inside_parallel_at_or_above"] = 94; // blocks at or above this size are first candidates for dedicated inside-block MPI teams
  pars_i["volterra_block_parallel_rank_fraction"] = 0.50;     // fraction of MPI ranks reserved for ordinary block-parallel work when large dedicated teams are created
  pars_i["volterra_min_team_size"] = 2;              // smallest MPI team allowed for inside-block parallelization of one large block
  pars_i["volterra_max_team_size"] = 32;             // absolute upper bound on MPI team size for one block (e.g. one node with 32 cores)
  pars_i["volterra_practical_team_size_cap"] = 8;    // practical cap used by scheduler; larger teams usually do not help because Volterra scaling saturates
  pars_i["volterra_max_dedicated_blocks"] = 5;       // maximum number of largest blocks that may receive dedicated inside-parallel teams
  pars_i["volterra_promote_medium_blocks"] = true;   // if true, medium-large blocks between the two size thresholds may also get dedicated teams if ranks remain
  pars_i["volterra_print_diagnostics"] = true;       // print scheduler decisions (chosen teams, block sizes, rank split) for tuning and debugging
  pars_i["volterra_lambda_rescue"] = true;            // if Q becomes too large, scan a nearby lambda0 in the culprit pseudo-particle block and recompute propagators
  pars_i["volterra_lambda_rescue_Q"] = 500.0;         // trigger rescue when raw sum_i G_i(beta) exceeds this threshold
  pars_i["volterra_lambda_rescue_target_Q"] = 50.0;   // accept the first nearby lambda shift with culprit block max endpoint below this value
  pars_i["volterra_lambda_rescue_negative_tol"] = 0.5; // also rescue if a diagonal endpoint weight is more negative than this
  pars_i["volterra_lambda_rescue_max_delta"] = 0.05;  // largest |delta lambda0| tested around the static lambda0
  pars_i["volterra_lambda_rescue_step"] = 0.001;      // lambda scan step used only when rescue is triggered
  pars_i["r_OCA_mesh"] = 1.0;                        // OCA uses equidistant mesh, and r_OCA_mesh makes equidistant mesh larger for Nt->Nt*r_OCA_mesh.
  pars_i["Sdyson_Bulla_alpha"] = 1.0;                // alpha<=1 moves the blend to lower frequency; reduce when MC noise is present.
  pars_i["Sdyson_Bulla_cutoff_fraction"] = 0.6;      // largest clean blend end is alpha*cutoff_fraction*omega_Nyquist.
  pars_i["Sdyson_Bulla_start_ratio"] = 0.6;          // blend starts at start_ratio*w_end.
  // 2) Read PARAMS file (with comments removed) and parse once
  if (inp2.str().empty()) {
    std::cerr << "PARAMS not present or empty\n";
    delete mpi;
    std::exit(1);
  }
  const std::string txt_params = inp2.str();
  json patch;
  try {
    patch = json::parse(txt_params);
  } catch (const nlohmann::json::parse_error& e) {
    ReportJsonParseError(e, txt_params, mpi);
  }
  // 3) Decide RealAxisQ using Delta from patch if present, else from defaults
  const std::string delta_file = patch.value("Delta", pars["Delta"].get<std::string>());
  const bool RealAxisQ = IsRealAxis(delta_file);
  if (RealAxisQ)
    pars.merge_patch(pars_r);
  else
    pars.merge_patch(pars_i);

  MergeKnownKeys(pars, patch);

  if (pars["fastFilesystem"] and mpi->rank!=mpi->master){
    delete mpi->out;
    mpi->out = new ofstream("debug."s+to_string(mpi->rank));
  }
  (*mpi->out)<<"I am rank="<<mpi->rank<<" of size="<<mpi->size<<endl;
  
  
  if (mpi->rank==mpi->master){ (*mpi->out) << pars.dump(4) << std::endl;
    //for (auto& element : pars.items()) std::cout << element.key() << " = " << element.value() << '\n';
  }
  
  
#ifdef _PCG
  RanPCG drand(0u);
  int base_seed = pars["iseed"].get<int>();
  //uint64_t real_iseed = drand.seed_repro(static_cast<uint64_t>(base_seed), static_cast<uint64_t>(mpi->rank));
  int real_iseed = base_seed + 1000*mpi->rank;
  drand.set_iseed_DistanceFromZero(real_iseed);
#endif
#ifdef _GSL
  int base_seed;
  if (pars["iseed"].get<int>()==0)
    base_seed = time(0);
  else
    base_seed = pars["iseed"].get<int>();
  int real_iseed = base_seed + mpi->rank*1000;
  RanGSL drand(real_iseed);
#endif  
  (*mpi->out)<<"iseed="<<real_iseed<<" in rank "<<mpi->rank<<endl;
  
  double beta   = pars["beta"];

  Cix cix;
  cix.Read_cix(pars["cix"],pars["U"],pars["mu"],beta);

  Diag2Matrix diag2matrix(pars["Norder"], cix);

  {
    double mu = pars["mu"];
    double U  = pars["U"];
    if (cix.Nsite!=1){ // In single-site DMFT we can change Slater F0 by simply adding U to double-occupied sites. We can't do that in cluster-DMFT.
      U=0;
    }
    for (int ib=0; ib<cix.Nblock; ++ib){
      Array<int,1>& indices = cix.Block[ib];
      for (int j=0; j<indices.size(); ++j){
	int N = cix.iN(indices(j)-1);
	double dE = U*N*(N-1)/2. - mu*N;
	cix.Ham[ib](j,j) += dE;
	//cout<<"Correcting Ham["<<ib<<"]("<<j<<","<<j<<")+="<<dE<<" H="<<cix.Ham[ib](j,j)<<endl;
      }
      //cout<<"ib="<<ib<<" N="<<N<<" dE="<<dE<<" H[i,i]=";
      //for (int i=0; i<indices.size(); ++i) cout<<cix.Ham[ib](i,i)<<", ";
      //cout<<endl;
    }
    Array<bool,1> bath_diag(cix.Nbath_unique_unique);
    for (int i=0; i<cix.Nbath; ++i){
      // Psi_a Delta_{ab} Psi_b have to correspond to the same orbital, and should not be
      // both creation or anhilation.
      bath_diag(abs(cix.bath_sigind0(i))) = not cix.bath_SC(i) and (cix.Fi_indx(2*i)/2 == cix.Fi_indx(2*i+1)/2);
    }
    //cout<<"bath_SC="<<bath_SC<<endl;
    for (int i=0; i<cix.Nbath_unique_unique; ++i)
      if (bath_diag(i)) cix.Ebath(i) -= mu;
    //cout<<"corrected cix.Ebath="<<cix.Ebath<<endl;
    if (mpi->rank==mpi->master) *mpi->out<<"F_atom="<<cix.F_atom<<" log(Zatom)="<<-cix.F_atom*beta<<endl;
  }
  {
    int two_2_i = 1;
    for (int i=0; i<Skeletons::Freq.extent(0); ++i){
      Skeletons::Freq(i) = two_2_i;
      two_2_i *= 2;
    }
  }

  {// Here we check if Sindex exists. If yes, we just read it. If no, we generate it and bcast to all nodes.
    bool sindex_exists = false;
    if (mpi->rank==mpi->master){
      ifstream fSindex_in("Sindex.dat");
      sindex_exists = fSindex_in.good();
      if (sindex_exists){
        *mpi->out<<"Reading frozen symmetry from Sindex.dat"<<endl;
      }else{
        *mpi->out<<"Sindex.dat not found. Determining symmetry and freezing it."<<endl;
        Determine_Symmetry(cix, diag2matrix);
        ofstream fSindex("Sindex.dat");
        cix.PrintSymmetry(fSindex,false);
        cix.PrintSymmetry(*mpi->out);
      }
    }
    mpi->bcast_one_boolean(sindex_exists);
    if (sindex_exists){
      ReadSindex("Sindex.dat", cix);
      if (mpi->rank==mpi->master) *mpi->out<<"Sindex.dat successfully read."<<endl;
    }else{
      Bcast_Sindex(cix, mpi);
    }
  }
  
  if (RealAxisQ)
    ComputeNew<RealAxis0,RealAxis1new,complex<double>,false>(pars,drand,cix,*mpi);
  else
    ComputeNew<ImagAxis0,ImagAxis1new,double,true>(pars,drand,cix,*mpi);
  
  delete mpi;
  return 0;
}
