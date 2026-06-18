#ifndef Real_axis
#define Real_axis
#include <vector>
#include <unordered_set>
#include "nca_oca_real_axis.h"
#include "real_axis_spectra.h"
#include "frequency_loops_real_axis.h"
#include "sub_expression_elimination.h"
#include "Eigen/Dense"
using Eigen::MatrixXcd;  // complex<double> matrix
#define _TIMEIT_

// conclusion 11 seems best, but 1 is probably good to keep sometimes.
// 11 -> sum_fp * product |G|
// 1  -> Im(G*Sigma)*(exp(-beta*x)+1)
// // Next we need to speed up measuring, which is very very slow. 10-times slower than move or 100-times slower than topology.

// 2 best for Stval, actually with 2 site is worse for Stval
// 11 best for Physical G, for 2 sites best
// 15 not bad at all

// two sites is best with 1 and 11. 2 is pretty bad for Stval
// 15 and 2 are bad for Sig.out
// 11 seems best for 2 sites both for Stval and G00

const int New_form_Phi=11;
// 1  -> Im(G*Sigma)*(exp(-beta*x)+1)
// 11 -> sum_fp * product |G|
//
// 2  -> widetilde{Sigma}*Phi^R
// 10 -> sum_fp*f(x)* product |G|
// 15 -> sum_fp * product |Gt|
//
// 1,2 -> 1,2
// 3,5 -> 10,15
// >=10   sum_fp
// ==10   sum_fp  contains f(x)
// >=15   Go_ and Go_t

const bool check_equivalence=false;

class RealAxis0 : public AnyAxis{
public:
  bool AdaptMesh;
  Array<double,2> Deltai; // Deltai(om) = -1/pi*Im(iDeltaw(om))
  Array<double,2> mDelta; // mDelta(om) = -1/pi*Im(iDeltaw(-om))
  Array<double,1> om;     // frequency for Delta
  Array<double,1> feo;    // fermi function on om mesh, which corresponds to Delta
  Array<double,1> xn_;
  Array<double,1> dh, logo;
  Array<double,1> Smi, Smr;
  Array<double,1> few;
  Array<double,2> PhiR;
  Array<complex<double>,2> iDeltaw;
  const Cix& cix;
  vector<int> G_importance; int OCA_cut_sum;
  //using AnyAxis::AnyAxis;  // <-- inherit constructors: All constructors of AnyAxis are also constructors of RealAxis0
  RealAxis0(double _beta_, const Cix& _cix_, json& pars, MPIcls& mpi) :
    AnyAxis(_beta_), cix(_cix_), AdaptMesh(pars["AdaptMesh"]),
    G_importance(cix.Npseudo_unique,0), OCA_cut_sum(10) {};
  int Coeff_size() const {
    return 0;
  }
  const Array<double,1>& freq() const override { return om; }
  void PrepareDelta(Array<complex<double>,2>& Deltaw, Array<double,1>& _om_, const json& pars, MPIcls& mpi)
  {
    if (cix.Nbath_unique_unique != Deltaw.extent(0))
      cout<<"ERROR Nbath_unique_unique="<<cix.Nbath_unique_unique<<" while Delta.extent(0)="<<Deltaw.extent(0)<<endl;
    om.reference(_om_);
    iDeltaw.reference(Deltaw);

    Array<int,1> sigind0(cix.Nbath_unique);
    for (int ib=0; ib<cix.Nbath; ++ib) sigind0(cix.bath_sigind(ib)) = cix.bath_sigind0(ib);
    
    //cout<<"sigind0="<<sigind0<<endl;
    Deltai.resize(cix.Nbath_unique,Deltaw.extent(1));
    mDelta.resize(cix.Nbath_unique,Deltaw.extent(1));
    for (int ib=0; ib<cix.Nbath_unique; ++ib){
      double ib_sg = (sigind0(ib) < 0 ? -1 : 1);
      int iib   = std::abs(sigind0(ib));
      int pos=om.size()-1;
      for (int i=0; i<Deltaw.extent(1); ++i){
	// some baths can have negative bath_sigind0 index, which is here eliminated, and we create two copies for positive and negative baths
	Deltai(ib,i) = -Deltaw(iib,i).imag()*ib_sg/M_PI;            // This is spectral function for Deltaw(om)
	intpar ip = InterpRight(-om(i),pos,om);
	mDelta(ib,i) = -linInterp(Deltaw,iib,ip).imag()*ib_sg/M_PI; // This is spectral function for Deltaw(-om)
      }
      //cout<<"Setting Delta["<<ib<<"] = Deltaw["<<iib<<"]*"<<ib_sg<<endl;
    }
    feo.resize(om.size());
    for (int i=0; i<om.size(); ++i) feo(i)=fermi(beta*om(i));

    if (true){
      ofstream out("debug_Delta.dat");
      for (int iw=0; iw<om.size(); ++iw){
	out<<om(iw)<<" ";
	for (int ib=0; ib<cix.Nbath_unique; ++ib) out<<Deltai(ib,iw)<<" "<<mDelta(ib,iw)<<" ";
	out<<endl;
      }
    }
  }
  auto ReadPseudoSig00(string Sfilename, MPIcls& mpi)
  { return _ReadPseudoSig00_<complex<double>>(Sfilename, cix.Npseudo_unique, mpi);}
  void InitializeSig00(Array<double,1>& xn, Array<complex<double>,2>& Sm,
		       const Array<double,1>& om, const Cix& cix, bool fastFilesystem, MPIcls& mpi)
  {
    Sm.resize(cix.Npseudo_unique,om.size()); 
    
    CreateSigma000(Sm,
		   om,Deltai,feo,cix,beta,fastFilesystem,mpi);
    xn.resize(om.size()); xn=om;
  }
  void InitializeStorage(int Npseudo_unique, const Array<double,1>& xn)
  {
    //// The current pseudo self-energy twidle
    //Smt.resize(Npseudo_unique,xn.size());    // I don't like Sm to be defined inside, while in imag axis it is outside.          
    //
    dh.resize(xn.size()); logo.resize(xn.size());
    Smi.resize(xn.size()); Smr.resize(xn.size());
  }

  auto DeterminePseudos(Array<complex<double>,2>& Gm, Array<complex<double>,2>& gm, Array<complex<double>,2>& Sm, Array<double,1>& xn, 
			const Cix& cix, double lambda0_fixed, int followPeak_fixed, int itt, json& pars, MPIcls& mpi)
  {
    if (lambda0_fixed!=0 or followPeak_fixed>0){
      if (lambda0_fixed!=0) (*mpi.out)<<"lambda0_fixed="<<lambda0_fixed<<" ";
      if (followPeak_fixed>0) (*mpi.out)<<"followPeak_fixed="<<followPeak_fixed<<" ";
      (*mpi.out)<<endl;
    }
    
    auto [lambda0,followPeak,_few_] = DeterminSpectralFunctions(Gm,Sm,xn,
								cix,1./beta,AdaptMesh,pars,mpi,pars["peak_pos"],lambda0_fixed,followPeak_fixed);
    //FindLargestEigenvalues(gm,
    //			   Gm,tau,cix,itt,mpi);
    double mu = pars["mu"];
    auto r = DeterminProbabilities(xn,Gm,cix,_few_,1./beta,mu,lambda0,followPeak,itt,mpi);
    
    xn_.reference(xn);
    dh.resize(xn.size());
    dh = GiveDh(xn);
    logo.resize(xn.size());
    logo = CalcLogoOnMesh(xn);
    
    few.resize(_few_.size());
    few=_few_;
    
    PhiR.resize(cix.Npseudo_unique,xn.size());
    PhiR = DeterminPhiR(Gm,xn,few,dh,logo,beta);

    if (followPeak>0){
      ofstream out("GroundState."+to_string(itt+1));
      int ib=followPeak-1;
      int N = cix.Block[ib].size(); // block size

      std::vector<int> iis;
      std::unordered_set<int> seen;
      for (int i=0; i<N; ++i){
	int _ii_ = cix.Sgm_symmetry[ib](i,i);
	int ii = abs(_ii_)-1;
	if (seen.insert(ii).second){ // true if inserted
	  iis.push_back(ii);
        }
      }
      for (int i=0; i<N; ++i){
	for (int j=0; j<N; ++j){
	  if (i==j) continue;
	  int _ii_ = cix.Sgm_symmetry[ib](i,j);
	  int ii = abs(_ii_)-1;
	  if (ii>=0){
	    if (seen.insert(ii).second){ // true if inserted
	      iis.push_back(ii);
	    }
	  }
	}
      }
      out<<"#  followPeak="<<followPeak<<" iis="<<iis<<endl;
      for (int ix=0; ix<xn.size(); ++ix){
	out<<xn(ix)<<" ";
	for (int i=0; i<iis.size(); ++i){
	  int ii = iis[i];
	  out << Gm(ii,ix).real()<<" "<<Gm(ii,ix).imag()<<" ";
	}
	out<<endl;
      }
    }
    
    return r;
  }
  vector<int> ImportanceOfPseudo(const Array<double,1>& xn, const Array<complex<double>,2>& Gm, std::array<double,2> cut={0.1,0.5}, bool debug=false)
  {
    Array<double,1> Gwidth(cix.Npseudo_unique);
    Array<double,1> dxh(xn.size()); dxh = GiveDh(xn);
    double x0=0.2;
    Array<double,1> chb(xn.size());
    for (int ix=0; ix<xn.size(); ++ix) chb(ix) = 1.0/cosh(0.5*beta*xn(ix));
    for (int ipp=0; ipp<cix.Npseudo_unique; ++ipp){
      double gi=0;
      for (int ix=0; ix<xn.size(); ++ix){
	double g = fabs(Gm(ipp,ix).imag())*dxh(ix);
	gi += g*chb(ix);
      }
      double sgm=1/gi;
      Gwidth(ipp) = sgm;
      //cout<<"width["<<ipp+1<<"]=" << sgm << endl;
    }
    std::vector<int> idx(cix.Npseudo_unique);
    std::iota(idx.begin(), idx.end(), 0); // 0..N-1
    std::sort(idx.begin(), idx.end(), [&](int a, int b){ return Gwidth(a) < Gwidth(b); });

    vector<int> _G_importance_(cix.Npseudo_unique,0);
    for (int k = 0; k < (int)idx.size(); ++k) {
      int ipp = idx[k];
      if (k>cut[0]*cix.Npseudo_unique) _G_importance_[ipp]=1;
      if (k>cut[1]*cix.Npseudo_unique) _G_importance_[ipp]=2;
    }
    if (debug){
      for (int k = 0; k < (int)idx.size(); ++k) {
	int ipp = idx[k];
	fmt::print(std::cout, "{:<3d} {:<4d} {:3d} {:10.3g}\n", k, ipp+1, _G_importance_[ipp], Gwidth(ipp));
      }
    }
    return _G_importance_;
  }  
  double NCA(Array<double,2>& Smt_new, Array<complex<double>,2>& G00, Array<complex<double>,2>& F00, Array<complex<double>,2>& F01,
	     bool cmpPhysG, double Q, const Array<double,1>& xn, const Array<complex<double>,2>& Gm, Cix& cix, MPIcls& mpi){
    double phi = R_NCA(Smt_new,G00,F00,F01,
		       cmpPhysG,Q,xn,Gm,few,PhiR,om,Deltai,mDelta,feo,cix,1./beta,mpi);
    return phi;
  }
  double OCA(Array<double,2>& Smt, Array<complex<double>,2>& G00, Array<complex<double>,2>& F00, Array<complex<double>,2>& F01,
	     bool cmpPhysG, double Q, const Array<double,1>& xn, const Array<complex<double>,2>& Gm, Cix& cix, json& pars, MPIcls& mpi){
    /*
    std::array<double,2> OCA_cut({pars["OCA_cut1"].get<double>(),pars["OCA_cut2"].get<double>()});
    if ((OCA_cut[0]-1.0)<-1e-5 or (OCA_cut[1]-1)<-1e-5){
      G_importance = ImportanceOfPseudo(xn, Gm, OCA_cut);
      OCA_cut_sum = pars["OCA_cut_sum"].get<int>();
    }else{
      OCA_cut_sum = 10;
    }
    */
    return R_OCA(Smt,G00,F00,F01,
		 cmpPhysG,Q,xn,Gm,cix,1./beta,few,PhiR,om,Deltai,mDelta,feo,pars,G_importance,OCA_cut_sum,mpi);
  }

  void SplineG(const Array<double,1>& , const Array<complex<double>,2>&,const Array<complex<double>,2>&) {};

  void StopAdaptingMesh(bool stop_adapting){if (stop_adapting) AdaptMesh=false;}
  void StartAdaptingMesh() {AdaptMesh=true;}
  void Update_Sm(Array<complex<double>,2>& Sm, Array<double,2>& Smt, const Array<double,1>& xn, const Cix& cix, MPIcls& mpi, bool check_causality=false)
  {
    // First we check for causality. We just force all diagonal components to be causal.
    if (check_causality){
      for (int ib=0; ib<cix.Nblock; ++ib){
	int N = cix.Block[ib].size();              // dimension of the block
	for (int i=0; i<N; ++i){
	  int ii = abs(cix.Sgm_symmetry[ib](i,i))-1;
	  for (int iw=0; iw<xn.size(); ++iw){
	    if (Smt(ii,iw)<0){
	      Smt(ii,iw)=1e-100;
	      (*mpi.out)<<" WARN Smt noncausal at ib="<<ib+1<<",i="<<i<<" om="<<xn(iw)<<endl;
	    }
	  }
	}
      }
    }
    Smi.resize(xn.size());
    Smr.resize(xn.size());
    for (int l=0; l<Smt.extent(0); l++){
      for (int i=0; i<Smt.extent(1); i++) Smi(i) = -Smt(l,i)*(1-few(i));
      Smr = KramarsKronig(Smi, xn, dh, logo);
      for (int i=0; i<xn.size(); i++) Sm(l,i) += complex<double>(Smr(i), Smt(l,i));
    }
  }
  double Difference(const Array<complex<double>,2>& Sm_new, const Array<complex<double>,2>& Sm, const Array<double,1>& Prob){
    //double diff = sum(abs(Sm_new-Sm)/(1.0+abs(Sm)))/Sm.size();
    double diff=0;
    double nrm = sum(abs(Prob));
    int nom = Sm_new.extent(1);
    for (int ii=0; ii<cix.Npseudo_unique; ++ii)
      diff += sum(abs(Sm_new(ii,Range::all())-Sm(ii,Range::all())))/nom*abs(Prob(ii))/nrm;
    return diff;
  }
  void PhysicalSelfEnergy(Array<complex<double>,2>& Sigm, Array<complex<double>,2>& Gc, Array<double,1>& density_matrix, double& Ekin, double& nf,
			  Array<complex<double>,2>& G00, const Array<complex<double>,2>& F00, const Array<complex<double>,2>& F01,
			  const Cix& cix, const BathBlocks& bb, const json& pars, int ist, MPIcls& mpi, string dysonFname)
  {
    bool debug = true;
    const double small=1e-10;
    Array<int,1> sigind0(cix.Nbath_unique), deg_sigind(cix.Nbath_unique_unique); sigind0=0; deg_sigind=0;
    for (int ib=0; ib<cix.Nbath; ++ib) sigind0(cix.bath_sigind(ib)) = cix.bath_sigind0(ib);
    for (int ib=0; ib<cix.Nbath_unique; ++ib) deg_sigind(abs(sigind0(ib)))+=1;
    //cout<<"sigind0="<<sigind0<<" deg_sigind="<<deg_sigind<<endl;
    // First average over components that are opposite in value, if any.
    Array<complex<double>,2> Gimp(cix.Nbath_unique_unique,om.extent(0)); Gimp=0;
    Array<complex<double>,2> Fimp1(cix.Nbath_unique_unique,om.extent(0)); Fimp1=0;
    Array<complex<double>,2> Fimp2(cix.Nbath_unique_unique,om.extent(0)); Fimp2=0;
    //cout<<"om.extent="<<om.extent(0)<<" Gimp.shape="<<Gimp.shape()<<" G00.shape="<<G00.shape()<<" Fimp1.shape="<<Fimp1.shape()<<" F00.shape="<<F00.shape()<<" F01.shape="<<F01.shape()<<endl;
    for (int ib=0; ib<cix.Nbath_unique; ++ib){
      double ib_sg = (sigind0(ib) < 0 ? -1 : 1);
      int iib   = std::abs(sigind0(ib));
      //cout << "setting Gim[" << iib << "]=G00[" << ib << "]*" << (ib_sg/deg_sigind(iib)) << endl;
      Gimp(iib,Range::all())  += G00(ib,Range::all())*(ib_sg/deg_sigind(iib));
      Fimp1(iib,Range::all()) += F00(ib,Range::all())*(ib_sg/deg_sigind(iib));
      Fimp2(iib,Range::all()) += F01(ib,Range::all())*(ib_sg/deg_sigind(iib));
    }
    if (false){
      ofstream out1("Gimp__.dat");
      ofstream out2("Fimp1_.dat");
      ofstream out3("Fimp2_.dat");
      for (int iw=0; iw<om.size(); ++iw){
	out1<<om(iw)<<" ";
	out2<<om(iw)<<" ";
	out3<<om(iw)<<" ";
	for (int ib=0; ib<cix.Nbath_unique_unique; ++ib){
	  out1 << Gimp(ib,iw).real()<<" "<<Gimp(ib,iw).imag()<<" ";
	  out2 << Fimp1(ib,iw).real()<<" "<<Fimp1(ib,iw).imag()<<" ";
	  out3 << Fimp2(ib,iw).real()<<" "<<Fimp2(ib,iw).imag()<<" ";
	}
	out1<<endl;
	out2<<endl;
	out3<<endl;
      }
    }
    Sigm.resize(cix.Nbath_unique_unique,om.extent(0));
    Gc.resize(cix.Nbath_unique_unique,om.extent(0));
    Array<complex<double>,2> SigD(cix.Nbath_unique_unique,om.extent(0));
    density_matrix.resize(cix.Nbath_unique_unique); density_matrix=0; // still need to implement
    Ekin=0; nf=0;// still need to implement
    
    ComputePhysicalSelfEnergyOnMesh(Sigm, SigD, Gc,
    				    om, Gimp, Fimp1, Fimp2, iDeltaw, bb, cix, *mpi.out, small, true);

    if (pars["AdaptPhysicalMesh"].get<bool>()){

      if (mpi.rank==mpi.master){
	/*
	int Nps=0;
	for (int iorb=0; iorb<cix.Nbath; ++iorb)
	  if (bb.orb_ind[iorb][1]==bb.orb_ind[iorb][2]) Nps += 1;
	int nps=0;
	for (int iorb=0; iorb<cix.Nbath; ++iorb){
	  if (bb.orb_ind[iorb][1]==bb.orb_ind[iorb][2]){
	    int iib = abs(cix.bath_sigind0(iorb));
	    App(nps,Range::all()) += bl::fabs(bl::imag(Gc(iib, bl::Range::all()))) / M_PI;
	    nps+=1;
	  }
	}
	*/
	int Nps = cix.Nbath_unique_unique;
	Array<double,2> App(Nps, om.extent(0)); App=0;
	for (int iorb=0; iorb<Nps; ++iorb)
	  App(iorb,Range::all()) = bl::fabs(bl::imag(Gc(iorb, bl::Range::all()))) / M_PI;
	
	Array<double,1> dh = GiveDh(om);
	Array<double,1> integral0(Nps);
	for (int i=0; i<Nps; ++i)
	  integral0(i) = sum(App(i,Range::all())*dh(Range::all()));
	
	double small_mesh=pars["Np_small"], alpha=pars["Nw_alpha"], La=pars["Np_La"], Lb=pars["Np_Lb"];
	int Nw_1=pars["Np_1"], Nw_2=pars["Np_2"], Nw_max=pars["Np_max"], Nw_itt=pars["Np_itt"];
	
	// mesh is constructed with the use of fw0 = \sum_p G_p(om)*dom/Integrate(G(om)*dom)
	// so that new points are such that each interval contributes the same weight
	bl::Array<double,1> fw2 = BuildFw2FromApp(App, dh, alpha, 1e-10);
	// Build an adaptive frequency mesh xn from a global weight fw2(ω) on an existing grid om.
	Array<double,1> om_new = BuildAdaptiveMeshFromFw2(om, fw2, App, La, Lb, Nw_1, Nw_2, small_mesh, Nw_max, Nw_itt, *mpi.out, debug);
	{
	  Array<complex<double>,2> Gimp_new(cix.Nbath_unique_unique,om_new.extent(0));
	  ResampleComplex2D(om, om_new, Gimp, Gimp_new);
	  Gimp.resize(cix.Nbath_unique_unique,om_new.extent(0));
	  Gimp = Gimp_new;
	}
	{
	  Array<complex<double>,2> Fimp1_new(cix.Nbath_unique_unique,om_new.extent(0));
	  ResampleComplex2D(om, om_new, Fimp1, Fimp1_new);
	  Fimp1.resize(cix.Nbath_unique_unique,om_new.extent(0));
	  Fimp1 = Fimp1_new;
	}
	{
	  Array<complex<double>,2> Fimp2_new(cix.Nbath_unique_unique,om_new.extent(0));
	  ResampleComplex2D(om, om_new, Fimp2, Fimp2_new);
	  Fimp2.resize(cix.Nbath_unique_unique,om_new.extent(0));
	  Fimp2 = Fimp2_new;
	}
	{
	  Array<complex<double>,2> iDeltaw_new(cix.Nbath_unique_unique,om_new.extent(0));
	  ResampleComplex2D(om, om_new, iDeltaw, iDeltaw_new);
	  iDeltaw.resize(cix.Nbath_unique_unique,om_new.extent(0));
	  iDeltaw = iDeltaw_new;
	}
	om.resize(om_new.extent(0));
	om = om_new;

	Sigm.resize(cix.Nbath_unique_unique,om_new.extent(0));
	Gc.resize(cix.Nbath_unique_unique,om_new.extent(0));
	SigD.resize(cix.Nbath_unique_unique,om_new.extent(0));
	ComputePhysicalSelfEnergyOnMesh(Sigm, SigD, Gc,
					om, Gimp, Fimp1, Fimp2, iDeltaw, bb, cix, *mpi.out, small, false);
    
	dh.resize(om.size());
	dh = GiveDh(om);
	Array<double,1> integral1(Nps); integral1=0;
	/*
	nps=0;
	for (int iorb=0; iorb<cix.Nbath; ++iorb){
	  if (bb.orb_ind[iorb][1]==bb.orb_ind[iorb][2]){
	    int iib = abs(cix.bath_sigind0(iorb));
	    integral1(nps) = sum( bl::fabs(bl::imag(Gc(iib, bl::Range::all()))) * dh ) / M_PI;
	    nps+=1;
	  }
	}
	*/
	for (int iorb=0; iorb<Nps; ++iorb)
	  integral1(iorb) =  sum( bl::fabs(bl::imag(Gc(iorb, bl::Range::all())))*dh ) / M_PI;
	
	for (int ib=0; ib<Nps; ++ib)
	  (*mpi.out)<<"# Final pintegral["<<ib<<"]="<< integral1(ib) <<" and original "<< integral0(ib) << " with diff="<< (integral1(ib)-integral0(ib)) << endl;

	Print_Physical1(dysonFname+"out", SigD);
      }

      mpi.Bcast_New_mesh(om, Sigm, Gc, iDeltaw);
      
      // Next we update Delta's that we use in MC, which have just been changed.
      Deltai.resize(cix.Nbath_unique,om.extent(0));
      mDelta.resize(cix.Nbath_unique,om.extent(0));
      for (int ib=0; ib<cix.Nbath_unique; ++ib){
	double ib_sg = (sigind0(ib) < 0 ? -1 : 1);
	int iib   = std::abs(sigind0(ib));
	int pos=om.size()-1;
	for (int i=0; i<om.extent(0); ++i){
	  // some baths can have negative bath_sigind0 index, which is here eliminated, and we create two copies for positive and negative baths
	  Deltai(ib,i) = -iDeltaw(iib,i).imag()*ib_sg/M_PI;            // This is spectral function for Deltaw(om)
	  intpar ip = InterpRight(-om(i),pos,om);
	  mDelta(ib,i) = -linInterp(iDeltaw,iib,ip).imag()*ib_sg/M_PI; // This is spectral function for Deltaw(-om)
	}
      }
      feo.resize(om.size());
      for (int i=0; i<om.size(); ++i) feo(i)=fermi(beta*om(i));
    }
	
    /*
    vector<bool> component_set(cix.Nbath_unique_unique,false);
    
    for (int ib=0; ib<(int)bb.blocks.size(); ++ib){
      int N = bb.blocks[ib].size();
      if (N==1){ // no matrix needed
	for (int iorb=0; iorb<cix.Nbath; ++iorb){
	  if (bb.orb_ind[iorb][0]==ib){
	    int iib0 = cix.bath_sigind0(iorb);
	    double ib_sg = (iib0 < 0 ? -1 : 1);
	    int iib   = std::abs(iib0);
	    if (not component_set[iib]){
	      double Eimp = cix.Ebath(iib);
	      Sigm(iib,Range::all())  = Fimp1(iib,Range::all())/Gimp(iib,Range::all());
	      SigD(iib,Range::all()) = -1.0/Gimp(iib,Range::all())*ib_sg-iDeltaw(iib,Range::all())*ib_sg-Eimp*ib_sg+om(Range::all());
	      Gc  (iib,Range::all()) =  1.0/(om(Range::all())-Eimp*ib_sg-Sigm(iib,Range::all())*ib_sg-iDeltaw(iib,Range::all())*ib_sg);
	      component_set[iib]=true;
	      //
	      for (int iw=0; iw<om.size(); ++iw){
		if (Sigm(iib,iw).imag()>0){ // causality violation
		  //std::cout << "WARNING Sg["<<iorb<<",om="<<om(iw)<<"] is not causal "<< Sigm(iib,iw).imag() << std::endl;
		  Sigm(iib,iw) = std::complex<double>(Sigm(iib,iw).real(),-small);
		}
	      }
	    }
	  }
	}
      }else{
	bool SC = bb.is_block_sc[ib];
	MatrixXcd Gm(N, N), Fm(N, N), Dm(N, N), om_Eimp(N, N); // Eigen matrices for linear algebra
	ofstream out1("Sm.dat");
	ofstream out2("Gmc.dat");
	ofstream out3("Dm.dat");
	ofstream out4("Gm.dat");
	ofstream out5("Fm.dat");
	ofstream out6("om_Eimp.dat");
	for (int iw=0; iw<om.size(); ++iw){
	  intpar ip;
	  if (SC) ip = Interp(-om(iw), om);
	  Gm.setZero();
	  Fm.setZero();
	  Dm.setZero();
	  om_Eimp.setZero();
	  // Diagonal: ω * I
	  for (int i=0; i<N; ++i) om_Eimp(i,i) = om(iw);
	  // Fill Gm, Fm, Dm, om_Eimp from bath data
	  for (int iorb=0; iorb<cix.Nbath; ++iorb){
	    if (bb.orb_ind[iorb][0]==ib){
	      int j=bb.orb_ind[iorb][1], k=bb.orb_ind[iorb][2];
	      int iib0 = cix.bath_sigind0(iorb);
	      double ib_sg = (iib0 < 0 ? -1 : 1);
	      int iib   = std::abs(iib0);
	      if (SC && j>=N/2 && k>=N/2){
		if (iib0>=0){
		  // SC/Nambu “conjugated” sector, namely, G_{22}(om+i*delta) = -G^H_{11}(-om+i*delta)
		  Gm(k,j) = -std::conj(linInterp(Gimp, iib, ip));
		  Dm(k,j) = -std::conj(linInterp(iDeltaw, iib, ip));
		  om_Eimp(k,j) += cix.Ebath(iib);
		}else{
		  // SC/Nambu “conjugated” sector, namely, G_{22}(om+i*delta) = -G^H_{11}(-om+i*delta)
		  Gm(k,j) = std::conj(linInterp(Gimp, iib, ip));
		  Dm(k,j) = std::conj(linInterp(iDeltaw, iib, ip));
		  om_Eimp(k,j) -= cix.Ebath(iib);
		}
	      }else{
		if (iib0>=0){
		  // “normal” state or “normal” sector
		  Gm(j,k) = Gimp(iib,iw);
		  Dm(j,k) = iDeltaw(iib,iw);
		  om_Eimp(j,k) -= cix.Ebath(iib);
		}else{
		  // “normal” state or “normal” sector
		  Gm(j,k) = -Gimp(iib,iw);
		  Dm(j,k) = -iDeltaw(iib,iw);
		  om_Eimp(j,k) += cix.Ebath(iib);
		}
	      }
	      if (SC){
		if (k>=N/2){
		  if (j>=N/2) // F_{22} has an extra minus sign, and takes down, which is stored in F2
		    Fm(k,j) = ib_sg*std::conj(linInterp(Fimp2, iib, ip));
		  else        // F_{12} has also an extra minus sign, and it destroys cooper pair, hence it is stored in F2
		    Fm(j,k) = -ib_sg*Fimp2(iib,iw);
		}else{
		  Fm(j,k) = ib_sg*Fimp1(iib,iw);
		}
	      }else{
		//Fm(j,k) = ib_sg*0.5*(Fimp1(iib,iw)+Fimp2(iib,iw));
		Fm(j,k) = ib_sg*Fimp1(iib,iw);
	      }
	    }
	  }
	  MatrixXcd Gminv = Gm.inverse();          // Gm^{-1}
	  MatrixXcd Sm    = Gminv * Fm;            // Sm = Gm^{-1} * Fm
	  MatrixXcd SmD   = om_Eimp - Gminv - Dm;  // Dyson: Σ = (ω - E_imp) - G^{-1} - Δ
	  if (SC){ // we need to symmetrize anomalous components
	    const int N2 = N/2;
	    MatrixXcd S12 = Sm.block(0, N2, N2, N2);
	    MatrixXcd S21 = Sm.block(N2, 0, N2, N2);
	    MatrixXcd A = 0.5*(S12 + S21);
	    Sm.block(0, N2, N2, N2) = A;
	    Sm.block(N2, 0, N2, N2) = A;
	  }
	  {
	    Eigen::MatrixXcd A = complex<double>(0.0, 0.5) * (Sm - Sm.adjoint());
	    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es(A);
	    if (es.info() != Eigen::Success){
	      (*mpi.out)<<"WARNING Eigenvalues of self-energy could not be computed!"<<endl;
	    }else{
	      bool causality_violated = false;
	      Eigen::VectorXd evals = es.eigenvalues();
	      for (int i = 0; i < evals.size(); ++i){
		if (evals(i)<0){
		  causality_violated = true;
		  (*mpi.out)<<"WARNING correcting causality S[om="<<om(iw)<<",ib="<<ib<<"] eval="<<evals(i)<<endl;
		}
		evals(i) = std::max(evals(i), small);
	      }
	      const Eigen::MatrixXcd& V = es.eigenvectors();
	      Eigen::MatrixXcd Apsd = V * evals.asDiagonal() * V.adjoint();
	      // Clean up Hermitian symmetry again
	      Apsd = 0.5 * (Apsd + Apsd.adjoint());
	      if (causality_violated){
		Sm = complex<double>(0.5,0) * (Sm + Sm.adjoint()) - complex<double>(0,1.)*Apsd;
	      }
	    }
	  }
	  MatrixXcd Ginv   = om_Eimp - Sm - Dm;    // G^{-1} = (ω - E_imp - Σ - Δ)
	  MatrixXcd Gmc  = Ginv.inverse();         // G = (ω - E_imp - Σ - Δ)^{-1}
	  // Storin Σ back into Blitz arrays
	  vector<bool> component_set(cix.Nbath_unique_unique,false);
	  for (int iorb=0; iorb<cix.Nbath; ++iorb){
	    if (bb.orb_ind[iorb][0]==ib){
	      int iib0 = cix.bath_sigind0(iorb);
	      double ib_sg = (iib0 < 0 ? -1 : 1);
	      int iib   = std::abs(iib0);
	      if (not component_set[iib]){
		int j=bb.orb_ind[iorb][1], k=bb.orb_ind[iorb][2];
		Sigm(iib,iw) = Sm(j,k)*ib_sg;  // just find the first occurence of this component.
		SigD(iib,iw) = SmD(j,k)*ib_sg;
		Gc(iib,iw) = Gmc(j,k)*ib_sg;
		component_set[iib] = true;
		// There might be a bit of a problem when G_{dn} is different than G_{up} in SC state, i.e.
		// both SC and magnetism are present. In this case we would need to write into iw that
		// correspond to -omega, which is not readily available.
		// We will take care of that specific case when needed.
	      }
	    }
	  }
	  out1<<om(iw)<<" ";
	  out2<<om(iw)<<" ";
	  out3<<om(iw)<<" ";
	  out4<<om(iw)<<" ";
	  out5<<om(iw)<<" ";
	  out6<<om(iw)<<" ";
	  for (int i=0; i<N; ++i)
	    for (int j=0; j<N; ++j){
	      out1 << Sm(i,j).real() << " " << Sm(i,j).imag()<< " ";
	      out2 << Gmc(i,j).real()<< " " << Gmc(i,j).imag()<< " ";
	      out3 << Dm(i,j).real() << " " << Dm(i,j).imag()<< " ";
	      out4 << Gm(i,j).real() << " " << Gm(i,j).imag()<< " ";
	      out5 << Fm(i,j).real() << " " << Fm(i,j).imag()<< " ";
	      out6 << om_Eimp(i,j).real()<< " " << om_Eimp(i,j).imag()<< " ";
	    }
	  out1<<endl;
	  out2<<endl;
	  out3<<endl;
	  out4<<endl;
	  out5<<endl;
	  out6<<endl;
	}
      }
    }
    */
  }
  Array<double,2> imagt(Array<complex<double>,2>& Gt){return imag(Gt);}
  //
  template<class T>
  void Print_Pseudo(const string& filename, const Array<T,2>& Smt_new, json& pars){
    ofstream dout(filename); dout.precision(12);
    dout<<"# beta="<<beta<<" U="<<pars["U"]<<" Norder="<<pars["Norder"]<<endl;
    for (int i=0; i<xn_.size(); ++i){
      dout<<xn_(i)<<" ";
      for (int ip=0; ip<Smt_new.extent(0); ++ip) Prints(dout,Smt_new(ip,i));
      dout<<endl;
    }
  }
  void Print_Physical(const string& filename, const Array<complex<double>,2>& G00, const Array<complex<double>,2>& F00, const Array<complex<double>,2>& F01)
  {
    ofstream sout(filename); sout.precision(12);
    for (int i=0; i<om.size(); ++i){
      sout<<om(i)<<" ";
      for (int ip=0; ip<G00.extent(0); ++ip) Prints(sout,G00(ip,i)); 
      for (int ip=0; ip<F00.extent(0); ++ip) Prints(sout,F00(ip,i));
      for (int ip=0; ip<F01.extent(0); ++ip) Prints(sout,F01(ip,i));
      sout<<endl;
    }
  }
};


void Give_GtimesFMatrix(Array<complex<double>,2>& GGi,
		  const Token& t0, const vector<mGexp>& Gmx, const vector<Array<double,1>>& Fi, vector<Array<complex<double>,2>>& H,
		  const vector<Array<int,2>>& alw, Array<complex<double>,2>& storage_GGi,
		  const Array<nFexp,1>& frms_, int twoN)
{
  if (t0.kind==Kind::R || t0.kind==Kind::C || t0.kind==Kind::I){
    const int l   = t0.idx;
    const int lp1 = (l+1) % twoN;                 // assumes twoN in scope
    const int N1  = alw[l].extent(0);
    const int N2  = alw[lp1].extent(0);
    // Bind GGi to caller-provided storage (no allocation)
    GGi.reference( storage_GGi(Range(0,N1-1), Range(0,N2-1)) );
    const auto& gc = Gmx[l].c;
    assert(gc.extent(0)==N1 && gc.extent(1)==N2);
    if (t0.kind==Kind::R) {
      // Copy ONLY the real part into complex buffer (imag = 0) — no temp 'real(A)'
      for (int i=0;i<N1;++i)
        for (int j=0;j<N2;++j)
          GGi(i,j) = complex<double>(gc(i,j).real(), 0.0);
    } else if (t0.kind==Kind::C) { // Kind::C
      // Copy complex matrix into pre-bound storage (no rebind to source)
      for (int i=0;i<N1;++i)
        for (int j=0;j<N2;++j)
          GGi(i,j) = gc(i,j);
    } else{ // if (t0.kind==Kind::I) {
      // Fill GGi in-place (no temp)
      Gmx[l].Set_Gt_ferm(GGi, frms_(l).x);
    }
    tiny::dot_diag(GGi, Fi[lp1]);
  }else{
    // Kind::H — just return a header referencing an existing matrix (zero-copy)
    GGi.reference(H[t0.idx]);    // shallow header copy; no data movement
  }
}

void Give_GMatrix(Array<complex<double>,2>& GGi,
		  const Token& t0, const vector<mGexp>& Gmx, vector<Array<complex<double>,2>>& H,
		  const vector<Array<int,2>>& alw, Array<complex<double>,2>& storage_GGi,
		  const Array<nFexp,1>& frms_, int twoN)
{
  if (t0.kind==Kind::R || t0.kind==Kind::C || t0.kind==Kind::I){
    const int l   = t0.idx;
    const int lp1 = (l+1) % twoN;                 // assumes twoN in scope
    const int N1  = alw[l].extent(0);
    const int N2  = alw[lp1].extent(0);
    // Bind GGi to caller-provided storage (no allocation)
    GGi.reference( storage_GGi(Range(0,N1-1), Range(0,N2-1)) );
    const auto& gc = Gmx[l].c;
    assert(gc.extent(0)==N1 && gc.extent(1)==N2);
    if (t0.kind==Kind::R) {
      // Copy ONLY the real part into complex buffer (imag = 0) — no temp 'real(A)'
      for (int i=0;i<N1;++i)
        for (int j=0;j<N2;++j)
          GGi(i,j) = complex<double>(gc(i,j).real(), 0.0);
    } else if (t0.kind==Kind::C) { // Kind::C
      // Copy complex matrix into pre-bound storage (no rebind to source)
      for (int i=0;i<N1;++i)
        for (int j=0;j<N2;++j)
          GGi(i,j) = gc(i,j);
    } else{ // if (t0.kind==Kind::I) {
      // Fill GGi in-place (no temp)
      Gmx[l].Set_Gt_ferm(GGi, frms_(l).x);
    }
  }else{
    // Kind::H — just return a header referencing an existing matrix (zero-copy)
    GGi.reference(H[t0.idx]);    // shallow header copy; no data movement
  }
}

void Give_GMatrix_only(Array<complex<double>,2>& GGi,
		       const Token& t0, const mGexp& Gt, vector<Array<complex<double>,2>>& H, Array<complex<double>,2>& storage_GGi, const Array<nFexp,1>& frms_, int twoN)
{
  if (t0.kind==Kind::R || t0.kind==Kind::C || t0.kind==Kind::I){
    const int N1  = Gt.c.extent(0);
    const int N2  = Gt.c.extent(1);
    // Bind GGi to caller-provided storage (no allocation)
    GGi.reference( storage_GGi(Range(0,N1-1), Range(0,N2-1)) );
    const auto& gc = Gt.c;
    if (t0.kind==Kind::R) {
      for (int i=0;i<N1;++i)
        for (int j=0;j<N2;++j)
          GGi(i,j) = complex<double>(gc(i,j).real(), 0.0);
    } else if (t0.kind==Kind::C) { // Kind::C
      // Copy complex matrix into pre-bound storage (no rebind to source)
      for (int i=0;i<N1;++i)
        for (int j=0;j<N2;++j)
          GGi(i,j) = gc(i,j);
    } else{ // if (t0.kind==Kind::I) {
      // Fill GGi in-place (no temp)
      const int l   = t0.idx;
      Gt.Set_Gt_ferm(GGi, frms_(l).x);
    }
  }else{
    // Kind::H — just return a header referencing an existing matrix (zero-copy)
    GGi.reference(H[t0.idx]);    // shallow header copy; no data movement
  }
}


class RealAxis1new{
public:
  const double& beta;
  const Array<double,1>& dlom;   // mesh for Delta
  const Array<double,1>& xn;     // mesh for PP
  const Array<double,2>& Deltai; // hybridization on the real axis
  const Array<double,1>& few;    // f(xn) from RealAxis0
  const Array<complex<double>,2>& Gm;// (Re(G_pp), widetilde{G}_pp)
  const Array<double,2>& PhiR;   // quantity to calculate Phi from self-energy, computed in RealAxis0
  const Cix& cix;
  const std::vector<int>& G_importance;
  const int& OCA_cut_sum;
  bool Offdiagonal;
  int N_min_order, N_max_order;//, Nsize, Nbath;
  Diag2Matrix tdiag2matrix, cdiag2matrix;
  int wNbin_lmax;                      // pars["wNbin"]
  double cutoffw, gauss_x_width;
  int current_iloop;
  Array<int,1> n2two;                  // 2^i
  // computed by Create_Fermi_Splines
  Spline1DLin<double> femx, lfemx;     // splines for f(-xn) and log(f(-xn))
  Spline1DLin<double> fep, lfep;       // splines for f(dlom), log(f(dlom))
  Array<nFexp,1> lmfx;                 // combination of f(-x) and log(f(-x))
  Array<Gexp,2> Gms;                  // stores G(x) and log(Im(G(x))), so that we can compute Im(G(x))/f(-x) at large negative x.
  // quantities to be updated for each diagram
  double sum_fp, sum_fp_t;             // combination of fermi factors for fake Phi
  int sgn_top, sgn_top_t;              // sign due to crossing or propagators
  Array<int,2> wtb, wtb_t;             // frequency to loop connection
  Array<double,1> sgn_tot, sgn_tot_t;      // sgn_tot(PPs) total sign of a diagram
  Array<double,1> pp_freq, pp_freq_t;  // pp_freq(2*Norder) contains all pp-frequencies
  Array<double,1> Dlty, Dlty_t;        // Dlty(Norder) and trial equivalent
  double cPhi_dd, tPhi_dd;             // product of all hybridization functions
  Array<complex<double>,1> cPhiv, tPhiv; // diagonal Phi for each block
  Array<intpar,1> pxs, pxs_t;          // pxs_(2*Norder)=Interp(pp_freq(2*Norder), xn), pp-frequency for all PP.
  Array<nFexp,1> fry_, fry_t;          // fry_(Norder) contains f(y) where y is frequency for hybridization function.
  Array<double,1> lfpt, lfpt_t;        // contains log(f(y1)*f(-y2)..*f(y_n)), the product needed for computing Sigma_j of pseudoparticle j.
  Array<complex<double>,2> Gmo_, Gmo_t;
  Array<complex<double>,1> mGppt, mGppt_t; // mGppt(PPs) is product of Gmt_
  //Array<complex<double>,2> Gmt_, Gmt_t;// Gmt_(PPs,2*Norder) contains (Re(G), widetilde{G}).
  Array<Gexp,2> Gmt_, Gmt_t;           // Gmt_(PPs,2*Norder) contains (Re(G), widetilde{G}).
  Array<double,2> phiR, phiR_t;        //phiR(tPPs,2*Norder)
  Array<int,1> xfsgn, xfsgn_t;         // contains signs for fermi functions in a diagram
  Array<nFexp,1> frms_, frms_t;        // frms_(2*Norder) contains f(pp_freq(2*Norder))
  //
  Array<complex<double>,2> Gmc_;       // Gmc_(PPs,2*Norder) contains (Re(G), Im(G)=-\widetilde{G}f(-x))
  Array<int,1> xftype_conj;            // which propagator should be conjugated
  Array<int,2> xftype_real;            // which propagator should be real
  Array<complex<double>,2> GG;
  bool topology_change;
  Array<bool,1> frequency_changed;
  bool cmpPhysG, subExpr, subExprElim;
  //Array<complex<double>,2> _F00_OCA, _F01_OCA, _F02_OCA, _F03_OCA;
  int sgns;
  vector<mGexp> Gmx;
  vector<Array<double,1>> Fi;
  vector<mGexp> Lft, Rgh;
  Array<double,2> dS;
  //
  Array<complex<double>,2> storage_GGi, storage_GGj, storage_GGl;
  Array<double,2> storage_GGk;
  mGexp Gtm;
  mGexp Gtp;
  Plan plan;
  vector<Plan> fplan;
  vector<Array<complex<double>,2>> H, gH, fH;
#ifdef _TIMEIT_         
  double t_lft, t_rgh, t_dphi, t_stval, t_interp, t_GF1, t_GF2, t_plan, t_pre_plan, t_post_plan, t_H, t_Hg, t_Hf, t_GGin, t_GGln, t_GGi, t_GGl;
#endif    
  RealAxis1new(const RealAxis0& ra0, const Array<double,1>& _xn_, const Array<complex<double>,2>& _Gm_,
	       const Cix& _cix_, json& pars, bool _cmpPhysG_,
	       MPIcls& mpi) :
    beta(ra0.beta), dlom(ra0.om), xn(_xn_), Gm(_Gm_), Deltai(ra0.Deltai),
    few(ra0.few), PhiR(ra0.PhiR), cix(_cix_), G_importance(ra0.G_importance), OCA_cut_sum(ra0.OCA_cut_sum),
    //ibth_ind(_ibth_ind_), ipsd_ind(_ipsd_ind_), Fm(_Fm_), max_pp_dim(_max_pp_dim_), 
    N_min_order(pars["N_min_order"]), N_max_order(pars["Norder"]), tdiag2matrix(N_max_order,_cix_), cdiag2matrix(N_max_order,_cix_), cmpPhysG(_cmpPhysG_)
    //Nsize(N_max_order+1), Nbath(_Nbath_)
  {
    wNbin_lmax = pars["wNbin"]; // how many slices of binned data to store in MC
    cutoffw = pars["cutoffw"];
    gauss_x_width = pars["gauss_x_width"]; // frequency is changed by this amount
    topology_change=true;
    frequency_changed.resize(2*N_max_order); frequency_changed=false;
    //xt.resize(wNbin_lmax); 
    //for (int i=0; i<wNbin_lmax; ++i) xt(i) = cutoffw*( -1.0 + (2*i+1.0)/wNbin_lmax );
    if (cmpPhysG){
      storage_GGi.resize(cix.Nmaxsize,cix.Nmaxsize);
      storage_GGj.resize(cix.Nmaxsize,cix.Nmaxsize);
      storage_GGl.resize(cix.Nmaxsize,cix.Nmaxsize);
      storage_GGk.resize(cix.Nmaxsize,cix.Nmaxsize);
      //Gtm.resize(cix.Nmaxsize,cix.Nmaxsize);
      //Gtp.resize(cix.Nmaxsize,cix.Nmaxsize);
    }
    Offdiagonal=false;
    for (int ib=0; ib<cix.Block.size(); ++ib)
      if (cix.Block[ib].size()>1){ Offdiagonal=true; break;}
    subExpr=false;
    subExprElim = pars["subExprElim"];
#ifdef _TIMEIT_    
    t_lft=t_rgh=t_dphi=t_stval=t_interp=t_GF1=t_GF2=0;
    t_plan=t_pre_plan=t_post_plan=0;
    t_H=t_Hg=t_Hf=t_GGin=t_GGln=t_GGi=t_GGl=0;
#endif    
  };
  int number_of(int Norder){return Norder+1;}
  int number_of2(int Norder){return Norder+1;}
  bool random_tw(Array<double,1>& freq, int iloop, double xi, int Norder, int itt)
  {
    // This function modifies the frequency array at index ind.
    // We called this Next_tarray in Python
    double x_new = freq(iloop) + (xi-0.5)*gauss_x_width; // change of w
    if (fabs(x_new)>cutoffw) return false;
    //if (iloop==0 and Norder==1 and fabs(x_new)>10.){
    //  cout<<itt<<" It happens with x="<<x_new<<endl;
    //}
    freq(iloop) = x_new;
    current_iloop = iloop;
    return true;
  }
  //double cutoff(){return cutoffw;} // frequency should be smaller than cutoffw
  double recmp_Q(const Array<complex<double>,2>& Gm, const Cix& cix)
  {
    Array<double,1> dh(xn.size());
    dh = GiveDh(xn);
    double Q=0;
    for (int ib=0; ib<cix.Nblock; ++ib){
      for (int i=0; i<cix.Block[ib].size(); ++i){
	int _ii_ = cix.Sgm_symmetry[ib](i,i);
	double dQ=0;
	if (_ii_!=0){
	  int ii = abs(_ii_)-1;
	  for (int i=0; i<Gm.extent(1); i++) dQ += Gm(ii,i).imag()*few(i)*dh(i);
	}
	Q += dQ/M_PI;
      }
    }
    return Q;
  }
  void Create_Fermi_Splines(const Array<complex<double>,2>& Gm)
  {
    // Here we fill the following arrays:
    //
    //   Spline1DLin<double> femx, lfemx;  // f(-xn), lof(f(-xn))
    //   Spline1DLin<double> fep, lfep;    // f(dlom), log(f(dlom))
    //   Array<nFexp,1> lmfx;              // combination of f(-x) and log(f(-x))
    //   Array<nGexp,2> Gms;               // stores G(x) and log(Im(G(x))), so that we can compute Im(G(x))/f(-x) at large negative x.
    //
    // We create special splines, which use combination of linear and cubic interpolation.
    // For fermi functions and their logs, the linear interpolation is far more accurate at large x, while
    // cubic interpolation is good only close to zero frequency. Hence we evaluate the function at the mid-point and
    // compare the quality of linear and cubic interpolation. We them remember the limits at which the cubic interpolation
    // becomes more precise.
    //
    //#define _CHECK_NANS
    femx.resize(xn.size()); lfemx.resize(xn.size());  // f(-xn), lof(f(-xn))
    for (int i=0; i<xn.size(); ++i){
      femx[i] = fermi(-xn(i)*beta);
      lfemx[i] = logfermi(-xn(i)*beta);
    }
    {
      Array<double,1> fem_mid(xn.size()-1), lfem_mid(xn.size()-1); // evaluation the same functions at the midpoint
      for (int i=0; i<xn.size()-1; ++i){
	double x = 0.5*(xn(i)+xn(i+1));
	fem_mid(i) = fermi(-x*beta);
	lfem_mid(i) = logfermi(-x*beta);
      }
      femx.splineIt(xn,fem_mid); // now using the mid-point knowledge, we can see at which point is linear interpolation becomes more accurate
      lfemx.splineIt(xn,lfem_mid,0,-1);
    }
    
    fep.resize(dlom.size()); lfep.resize(dlom.size());   // f(dlom), log(f(dlom))
    for (int i=0; i<dlom.size(); ++i){
      fep[i] = fermi(dlom(i)*beta);
      lfep[i] = logfermi(dlom(i)*beta);
    }
    {
      Array<double,1> fep_mid(dlom.size()-1), lfep_mid(dlom.size()-1); // here we also evaluate mid-points
      for (int i=0; i<dlom.size()-1; ++i){
	double x = 0.5*(dlom(i)+dlom(i+1));
	fep_mid(i) = fermi(x*beta);
	lfep_mid(i) = logfermi(x*beta);
      }
      fep.splineIt(dlom,fep_mid); // so that we can use cubic spline only in the regions where it is more accurate
      lfep.splineIt(dlom,lfep_mid,0,-1);
    }

    // This class uses logaritm of the fermi function lmfx -> f(-x) and log(f(-x))
    lmfx.resize(xn.size()); // Array<nFexp,1> lmfx(xn.size());
    for (int i=0; i<xn.size(); ++i) lmfx(i) = nFexp( femx[i], lfemx[i], -xn(i)*beta );

    int Npseudo_unique = Gm.extent(0);
    Array<double,2> common_a(Npseudo_unique,xn.size());
    for (int ix=0; ix<xn.size(); ++ix){
      // Computes average Im(G) in every block, and across equivalent blocks, so that
      // we take average over all blocks that share a common component
      for (int i=0; i<cix.unique_blocks.size(); ++i){
	double average_Gim=0;
	int how_many = 0;
	for (int j=0; j<cix.unique_blocks[i].size(); ++j){
	  int ib = cix.unique_blocks[i][j];
	  int Ni = cix.Block[ib].size();
	  for (int a=0; a<Ni; ++a)
	    if (cix.Sgm_symmetry[ib](a,a)!=0){
	      int ii = abs(cix.Sgm_symmetry[ib](a,a))-1;
	      average_Gim += fabs(Gm(ii,ix).imag());
	      how_many += 1;
	    }
	}
	average_Gim /= how_many;
	// pretending that |Gtc.imag|>=1e-16 because otherwise we get nan's when taking log
	double av = (average_Gim>1e-16) ? log(average_Gim) : log(1e-16);
	// now we can compute a of components in the block and across the blocks.
	// They must have the same exponent a if they share one common matrix element.
	// DEBUG
	//set<int> iis;
	for (int j=0; j<cix.unique_blocks[i].size(); ++j){
	  int ib = cix.unique_blocks[i][j];
	  int Ni = cix.Block[ib].size();
	  for (int a=0; a<Ni; ++a){
	    for (int b=0; b<Ni; ++b){
	      int ii = abs(cix.Sgm_symmetry[ib](a,b))-1;
	      if (ii>=0){
		common_a(ii,ix) = av;
		// DEBUG
		//iis.insert(ii);
	      }
	    }
	  }
	}
	//if (ix==xn.size()/2) cout<<"average_Gim["<<ix<<","<<cix.unique_blocks[i]<<"]="<<average_Gim<<" common_a["<<iis<<"]="<< av <<endl;
      }
    }
    // This class stores the Green's function values and the logarith of the imaginary part, so that we can devide the imaginary
    // part by f(-w), i.e., Im(G(x))/f(-x) can be computed when needed.
    Gms.resize(cix.Npseudo_unique,xn.size()); //Array<nGexp,2> Gms(Npseudo_unique,xn.size());
    for (int ii=0; ii<cix.Npseudo_unique; ++ii){
      for (int ix=0; ix<xn.size(); ++ix){
	Gms(ii,ix).set(Gm(ii,ix),common_a(ii,ix),lmfx(ix));
	//if (ix==xn.size()/2) cout<<"Gms(ii="<<ii<<",ix="<<ix<<")="<<Gms(ii,ix)<<" common_a="<<common_a(ii,ix)<<" Gm="<<Gm(ii,ix)<<" lmfx="<<lmfx(ix)<<endl;
      }
    }
  }
  void InitializeStorageMC(){//int Npseudo_unique, int Nbath_unique, int max_pp_dim=2){
    n2two.resize(2*N_max_order+1);
    for (int i=0; i<=2*N_max_order; ++i) n2two(i) = (1<<i);

    fry_.resize(N_max_order); fry_t.resize(N_max_order);
    Dlty.resize (N_max_order); Dlty=0; Dlty_t.resize (N_max_order); Dlty_t=0;
    
    lfpt.resize(2*N_max_order); lfpt_t.resize(2*N_max_order);
    pp_freq.resize(2*N_max_order); pp_freq_t.resize(2*N_max_order);
    pxs.resize(2*N_max_order); pxs_t.resize(2*N_max_order);
    frms_.resize(2*N_max_order); frms_t.resize(2*N_max_order);
  }
  void DeleteStorageMC(){
    n2two.resize(0);
    pp_freq.resize(0); pp_freq_t.resize(0);
    Dlty.resize (0); Dlty_t.resize (0);
    pxs.resize(0); pxs_t.resize(0);
    fry_.resize(0); fry_t.resize(0);
    lfpt.resize(0); lfpt_t.resize(0);
    //
    frms_.resize(0); frms_t.resize(0);
  }
  //void Zero_time_freq(Array<double,1>& freq1, int Norder) {freq1 = 0;}
  
  template<class RNG>
  void SetRandom_time_freq(Array<double,1>& freq, int Norder, RNG& drand)
  {
    for (int i=0; i<Norder+1; ++i) freq(i) = (0.5-drand())*gauss_x_width;
  }
  bool ReadStatusFile(Array<int,2>& cdiag, Array<int,1>& cdirspin, Array<double,1>& freq1,
		      const string& status_file, ostream& log)
  {
    Array<int,2> ndiag;
    Array<int,1> ndirspin;
    Array<double,1> nfreq;
    
    ifstream inp(status_file);
    int Norder=0;
    if (!(inp >> Norder)) return false;
    if (Norder>N_max_order or Norder<N_min_order){
      log<<"Read diagram is of too high/low order not("<<N_min_order<<"<="<<Norder<<"<="<<N_max_order<<") hence ignoring file "<<status_file<<endl;
      return false;
    }
    ndiag.resize(Norder,2);
    for (int i=0; i<Norder; ++i) inp >> ndiag(i,0) >> ndiag(i,1);
    if (!inp) return false;
    ndirspin.resize(Norder);
    for (int i=0; i<Norder; ++i) inp >> ndirspin(i);
    if (!inp) return false;
    nfreq.resize(number_of2(Norder));
    for (int i=0; i<nfreq.extent(0); ++i) inp >> nfreq(i);
    if (!inp) return false;
    for (int i=0; i<nfreq.extent(0); ++i){
      if (nfreq(i)<-cutoffw or nfreq(i)>cutoffw){
	log<<"Note: status file contains frequency "<<nfreq(i)<<", which is outside cutoff="<<cutoffw<<endl;
	return false;
      }
    }
    cdiag.resize(ndiag.shape()); cdiag=ndiag;
    cdirspin.resize(ndirspin.shape()); cdirspin=ndirspin;
    freq1.resize(nfreq.shape()); freq1=nfreq;
    log<<"Note: status file " << status_file << " succesfully read. Starting from configuration read." << endl;
    return true;
  }
  complex<double> Topology_trial_move(int& tPPs, int itt, const Array<int,2>& diag, const Array<int,1>& dirspin, const Array<double,1>& freq, ostream& cout)
  { // Input Arguments:
    //	   Array<int,2>    diag
    //	   Array<int,1>    dirspin
    //     Array<double,1> freq                 // current frequencies for all Norder+1 loops
    //  The quantities from this class that are used, but not modified:
    //
    //	   Cix                      cix
    //	   Array<double,1>          dlom        // frequency on which hybridizations delta is specified
    //	   Array<double,2>          Deltai      // imaginary part of hybridization Delta, as needed on real axis
    //     Array<double,1>          xn          // frequency on which PP are specified
    //     Array<complex<double>,2> Gm          // PP-green's function values
    //	   Spline1DLin<double> fep, lfep        // fermi function quantities on hybridization mesh
    // The quantities from this class that are modified:
    //     Array<int,2>             wtb_t       // The matrix of arguments, which uniquely defines frequency loops
    //     Array<int,1>             xfsgn_t     // contains signs for fermi functions in a given diagrams
    //     Array<nFexp,1>           fry_t       // contains fermi functions of hybridization frequencies f(y). Note always f(+y).
    //     Array<double,1>          Dlty_t(n)   // hybridizations in this diagram
    //     complex<double>          tPhi_dd     // Product of Deltas for this diagram
    //
    //     Array<double,1>          pp_freq_t   // contains pseudo-particle frequencies
    //     Array<intpar,1>          pxs_t       // interpolations for all pp-propagators values
    //     Array<double,1>          lfpt_t      // contains log(f(y1)*f(-y2)..*f(y_n)), the product needed for computing Sigma_j of pseudoparticle j.
    //     double                   sum_fp_t    // combination of fermi functions as needed for fake Phi
    //     Array<complex<double>,2> Gmt_t(tPPs,2*Norder) // all necessary pp propagators in the complex form (Re(G), widetilde{G})
    //     Array<complex<double>,1> tPhiv(iblck) // Resulting fake Phi functional being sampled, and only for each block
    //
    // Note: if we change dirspin only, we can update wtb to -1*wtb_old for those wtb where dire has been changed.
    // The following quantities do not change:
    //    xfsgn_t, fry_t[:], lfpt_t[:], sum_fp_t, sgn_top
    //
    // Quantities that change:
    //    Dlty_t for those iy on which dirspin has changed ( either wtb or spin).
    //    tPhi_dd has to be updated.
    //    Gmt_t is updated if pseudo is changed. Should have new and old pseudo
    //
    bool const DEBUG=false;
    int Norder = diag.extent(0);
    DetermineLoops(wtb_t,/**/diag, dirspin);
    //
    DetermineFermiSigns_fsgn_only(xfsgn_t,/**/wtb_t);
    tPhi_dd = 1;
    for (int iy=0; iy<Norder; ++iy){ // 0:Norder for all possible configurations
      int dire = dirspin(iy) % 2;    // direction of hybridization propagator
      int iorb = dirspin(iy) / 2;    // its spin 
      int s = cix.bath_sigind(iorb); // actual entry for hybridization
      double y = freq(1+iy);         // frequency on hybridization
      intpar jpy = Interp(y, dlom);  // 
      fry_t(iy).set(fep(jpy), lfep(jpy), beta*y ); // fermi function for f(y) stored as logarithm log(f(y))
      if (wtb_t(iy,iy+1)>0){
	Dlty_t(iy) = linInterp(Deltai, s, jpy); // actual interpolation and saving the interpolated hybridization
      }else{
	intpar jmpy = Interp(-y, dlom);
	Dlty_t(iy) = linInterp(Deltai, s, jmpy); // actual interpolation and saving the interpolated hybridization
      }
      tPhi_dd *= Dlty_t(iy);
    }
    //if (DEBUG) cout<<" wtb="<<wtb_t<<" xfsgn="<<to_binary(xfsgn_t)<<" Phi_dd="<<tPhi_dd<<" Dlty="<<Dlty_t(Range(0,Norder-1))<<" fry="<<fry_t(Range(0,Norder-1))<<endl;
    auto& tallowed = tdiag2matrix.Cmp_diagonal3(diag, dirspin, cix);      // map<int,Array<int,2>>
    for (int j=0; j<2*Norder; ++j){ // going over all pseudo-particles
      // the G[j] argument can be omega+y_0-y_2+y_3+.... Such combination would correspond to array wtb[Norder+j,:]=[1,1,0,-1,1,..]
      Range all_freq(0,Norder);
      double x = sum(freq(all_freq)*wtb_t(Norder+j,all_freq)); // This is the frequency argument for this PP, namely, sum(wtb[Norder+j,:]*freq[:])
      pp_freq_t(j) = x;                           // remember PP-frequency
      intpar px = Interp(x, xn);                  // interpolation of x in xn
      pxs_t(j) = px;
      //
      frms_t(j).set( femx(px), lfemx(px), -beta*x); // fermi function f(-x) with argument x of pseudo propagator j.
      //if (DEBUG) cout<<"  x["<<j<<"]="<<x<<" px="<<px<<" f(-x)="<<frms_t(j)<<endl;
    }
    /*
    for (int j=0; j<2*Norder; ++j){ // going over the ring
      // Note: lfpt(ipseudo) = log(  \product_{iy} f(\pm_{ipseudo,iy} y)  )
      double lfp=0; // This is the product of fermi functions needed for Self-energy Sigma_j of pseudoparticle j.
      for (int iy=0; iy<Norder; ++iy){
	if (xfsgn_t(j) & (1<<iy)){ // xfsgn is represented as binary array in an integer. 
	  lfp += fry_t(iy).plf(); // log(f(y))
	}else{
	  lfp += fry_t(iy).mlf(); // log(f(-y))
	}
      }
      lfpt_t(j)=lfp; // This is the product of fermi functions needed for Self-energy Sigma_j of pseudoparticle j.
    }
    */
    {// This is new implementation of lfpt, which is based on realization that lfpt(ip)=lfpt(0)-beta*om_0+beta*om_{ip}
      // this is because Phi cntribution is invariant with respect to ip, and exp(-beta*om_{ip})*fpt_{ip}=const.
      // Note: lfpt(ipseudo) = log(  \product_{iy} f(\pm_{ipseudo,iy} y)  )
      double lfp=0; // This is the product of fermi functions needed for Self-energy Sigma_j of pseudoparticle j.
      for (int iy=0; iy<Norder; ++iy){
	if (xfsgn_t(0) & (1<<iy)){ // xfsgn is represented as binary array in an integer. 
	  lfp += fry_t(iy).plf(); // log(f(y))
	}else{
	  lfp += fry_t(iy).mlf(); // log(f(-y))
	}
      }
      lfpt_t(0)=lfp; // This is the product of fermi functions needed for Self-energy Sigma_j of pseudoparticle j.
      for (int j=1; j<2*Norder; ++j) lfpt_t(j) = lfpt_t(0)+frms_t(0).x-frms_t(j).x;
      //if (DEBUG) cout<<"lfpt="<< lfpt_t <<endl;
    }
    tPPs=0;
    for (const auto& [ib, alw] : tallowed)  tPPs += alw.extent(0);
    Gmt_t.resize(tPPs,2*Norder); //Gmt_t = 0;
    if (New_form_Phi==2) phiR_t.resize(tPPs,2*Norder);
    if (New_form_Phi>=15) {mGppt_t.resize(tPPs);mGppt_t(Range(0,tPPs-1))=1.0;Gmo_t.resize(tPPs,2*Norder);}
    int idx=0;
    for (const auto& [ib, alw] : tallowed){ // over all blocks that contribute on diagonal
      for (int j=0; j<alw.extent(0); ++j){  // in principle more than one term can contribute per block
	for (int i=0; i<2*Norder; ++i){     // over the ring of PPs
	  //int _ab_ = alw(j,i);
	  //int _a_ = abs(_ab_);
	  //int ab = cix.inBlock(_a_);      // _ab_ is given in terms of atomic state index, but ab is index within the block. cix.inBlock translates this.
	  //int _ib_ = cix.iBlck(_a_)-1;
	  //int ii = cix.Sgm_symmetry[_ib_](ab,ab)-1; // Green's function index for propagator G(a0,b0)==G[|_ii0_|-1]
	  int _ii_ = alw(j,i);
	  int sg = isign(_ii_);
	  int ii = sg*_ii_ - 1;
	  //int ii = abs(_ii_)-1;
	  if (ii>=0){
	    //double sg = sign(_ab_);
	    //double sg = sign(_ii_);
	    //Gmt_t(idx,i) = linInterp(Gm, ii, pxs_t(i)) * sg;  // now actual interpolation of Gm and setting value of Gm*F. This is not just Gm, but multiplied by F.
	    Gmt_t(idx,i).linInterp(Gms, ii, pxs_t(i), sg);
	    //if (DEBUG) cout<<"  Gmt["<<ib<<","<<j<<","<<i<<"]="<<Gmt_t(idx,i)<<" Gms[i]="<<Gms(ii,pxs_t(i).i)<<" Gms[i+1]"<<Gms(ii,pxs_t(i).i+1)<<endl;
	    if (New_form_Phi==2) phiR_t(idx,i) = linInterp(PhiR, ii, pxs_t(i));
	    if (New_form_Phi>=15) Gmo_t (idx,i) = linInterp(Gm,   ii, pxs_t(i))*(sg+0.0);  // like for single orbital
	  }else{
	    Gmt_t(idx,i).set2zero();
	    if (New_form_Phi==2) phiR_t(idx,i) = 0;
	    if (New_form_Phi>=15) Gmo_t(idx,i) = 0;
	  }
	  if (New_form_Phi>=15){
	    mGppt_t(idx) *= Gmo_t(idx,i);
	    //cout << "mGppt["<<idx<<"]="<<mGppt_t(idx)<<" Gmo_t["<<idx<<","<<i<<"]="<<Gmo_t(idx,i)<<endl;
	  }
	  //if (DEBUG)  cout<<"   ib="<<ib<<" j="<<j<<" i="<<i<<" ii="<<ii<<" Gmt_t["<<idx<<","<<i<<"]="<<Gmt_t(idx,i)<<" pxs="<<pxs_t(i)<<" _ii_="<<_ii_<<endl;
	}
	idx += 1;
      }
    }
    //if (DEBUG) cout<<"Gmt="<<Gmt_t<<endl;
    if (New_form_Phi>=10){
      sum_fp_t=0;
      for (int iy=0; iy<Norder; ++iy){ // Now compute the product of fermi functions entering the physical green's function
	// This looks like a strange algorithm, but the idea is the following:
	// It is expensive to compute a lot of exponents, and the reults would be exp(x_1)+exp(x_2)+....., where x_i vary in size a lot.
	// For MC weight, we don't need a very precise value of Psi. It is sufficient to redefine
	// Psi where only the biggest contribution enters instead of the sum, i.e., exp(max(x_i)) as an estimation.
	// 
	double max_fp=-10*beta;
	for (int ip=0; ip<2*Norder; ++ip){
	  // sum_fp[ic] = \sum_{iy_ext}\sum_{ip}\product_{iy\ne iy_ext} f(\pm_{ip,iy} y)
	  // sum_fp[ic] = \sum_{ip}\sum_{iy_ext} 1/f(\pm_{ip,iy_ext} y_ext) \product_{iy} f(\pm_{ip,iy} y)
	  double lfp = lfpt_t(ip); // lfp(ipseudo) = log(  \product_{iy} f(\pm_{ipseudo,iy} y)  )
	  if (xfsgn_t(ip) & (1<<iy)){
	    lfp -= fry_t(iy).plf(); // now dividing : + log( 1/f(\pm_{ip,iy_ext} y_ext) )
	  }else{
	    lfp -= fry_t(iy).mlf(); // now dividing : + log( 1/f(\pm_{ip,iy_ext} y_ext) )
	  }
	  if (New_form_Phi==10) lfp += frms_t(ip).mlf(); // we are multiplying by f(+x_ip)
	  if (lfp>max_fp) max_fp=lfp; // we would need to take exponent of lfp, but we will just find the largest, and take exponent only once.
	  //double crnt= (xfsgn_t(ip) & (1<<iy)) ? fry_t(iy).plf(): fry_t(iy).mlf();
	  //cout<<".... for iy="<<iy<<" ip="<<ip<<" lfp="<<lfp<<" total="<<lfpt_t(ip)<<" current="<<crnt<<" max_fp="<<max_fp<<endl;
	}
	// instead of computing \sum_{ip} exp(lfp_{iy_ext,ip}) we just compute exp(max(lfp_{iy_ext,ip})), which is
	// a reasonable proxy, good enough for MC walk.
	sum_fp_t += exp(max_fp);
      }
      //if (DEBUG) cout<<"sum_fp="<<sum_fp_t<<endl;
    }
    //sgn_top_t = SignSort(diag);    // sign due to topological crossing in the diagram
    tPhiv.resize(tallowed.size());
    complex<double> Phi_pp=0;
    idx=0; int _idx_=0;
    for (const auto& [ib, alw] : tallowed){ // over all blocks that contribute on diagonal
      complex<double> phi_pp=0;
      for (int j=0; j<alw.extent(0); ++j){ // in principle more than one term can contribute per block
	Gexp pd; pd.Init2Unity();
	for (int ip=0; ip<2*Norder; ++ip){
	  pd *= Gmt_t(idx,ip);
	  //if (DEBUG) cout<<"pd["<<ib<<","<<j<<","<<ip<<"]="<<pd<<" Gmt="<<Gmt_t(idx,ip)<<endl;
	}
	//
	if (New_form_Phi==1){
	  for (int ip=0; ip<2*Norder; ++ip)
	    phi_pp -= pd.get0(lfpt_t(ip), frms_t(ip))/M_PI;
	}else if (New_form_Phi==2){
	  for (int ip=0; ip<2*Norder; ++ip){
	    Gexp Phis_over_Gm(pd);
	    Phis_over_Gm.divide_by_positive(Gmt_t(idx,ip)); // We should divide by diagonal Gms, but Gmt_t stores G*F, and F is either +-1, hence we
	    Phis_over_Gm.matmul_lferm(lfpt_t(ip));
	    phi_pp += Phis_over_Gm.getImt(frms_t(ip).lf)*phiR_t(idx,ip);
	  }
	}else if (New_form_Phi==10 or New_form_Phi==11){
	  phi_pp += abs(pd.c) * sum_fp_t * (2*Norder);
	}else if (New_form_Phi==15){
	  phi_pp += mGppt_t(idx) * sum_fp_t * (2.0*Norder);
	}
	//if (DEBUG) cout<<" phi_pp["<<ib<<","<<j<<"]="<<phi_pp<<" pd="<<pd<<" sum_fp="<<sum_fp_t<<endl;
	idx += 1;
      }
      phi_pp /= (2.*Norder);
      Phi_pp += phi_pp;
      tPhiv(_idx_) = phi_pp*tPhi_dd;
      _idx_ += 1;
    }
    complex<double> ret = Phi_pp*tPhi_dd;
    //if (DEBUG) cout<<"ret="<<ret<<" Phi="<<tPhiv<<" Phi_pp="<<Phi_pp<<" Phi_dd="<<tPhi_dd<<" sum_fp="<<sum_fp_t<<endl;
    /*
    if (false){
      //if (abs(ret)==0){
      cout<<"ERROR ret="<<ret<<" Phi_pp="<<Phi_pp<<" tPhi_dd="<<tPhi_dd<<" tallowed.size="<<tallowed.size()<<" tPPs="<<tPPs<<endl;
      idx=0; int _idx_=0;
      for (const auto& [ib, alw] : tallowed){ // over all blocks that contribute on diagonal
	complex<double> phi_pp=0;
	for (int j=0; j<alw.extent(0); ++j){ // in principle more than one term can contribute per block
	  Gexp pd; pd.Init2Unity();
	  for (int ip=0; ip<2*Norder; ++ip) pd *= Gmt_t(idx,ip);
	  cout<<"idx="<<idx<<" pd="<<pd<<" Gmt.shape="<<Gmt_t.shape()<<" Gmt_t="<<Gmt_t(idx,Range::all())<<endl;
	  //
	  if (New_form_Phi==1){
	    for (int ip=0; ip<2*Norder; ++ip)
	      phi_pp -= pd.get0(lfpt_t(ip), frms_t(ip))/M_PI;
	  }else if (New_form_Phi==2){
	    for (int ip=0; ip<2*Norder; ++ip){
	      Gexp Phis_over_Gm(pd);
	      Phis_over_Gm.divide_by_positive(Gmt_t(idx,ip)); // We should divide by diagonal Gms, but Gmt_t stores G*F, and F is either +-1, hence we
	      cout<<"Phis_over_Gm="<<Phis_over_Gm<<" lfpt_t["<<ip<<"]="<<lfpt_t(ip)<<endl;
	      Phis_over_Gm.matmul_lferm(lfpt_t(ip));
	      phi_pp += Phis_over_Gm.getImt(frms_t(ip).lf)*phiR_t(idx,ip);
	      cout<<"Phis_over_Gm.getImt="<<Phis_over_Gm.getImt(frms_t(ip).lf)<<" phiR="<<phiR_t(idx,ip)<<" phi_pp="<<phi_pp<<endl;
	    }
	  }else if (New_form_Phi==10 or New_form_Phi==11){
	    double pd_im = pd.p*exp(pd.a), pd_abs = abs(pd.c);
	    phi_pp += pd_abs * sum_fp_t * (2*Norder);
	    //phi_pp += pd_im * sum_fp_t * (2*Norder);
	    //for (int ip=0; ip<2*Norder; ++ip)
	    //  phi_pp -= pd_im*exp(frms_t(ip).x)*sum_fp_t/M_PI;
	  }else if (New_form_Phi==15){
	    phi_pp += mGppt_t(idx) * sum_fp_t * (2.0*Norder);
	    //cout<<" mGpp["<<idx<<"]="<< mGppt_t(idx) <<" sum_fp_t="<<sum_fp_t<<endl;
	  }
	  idx += 1;
	}
	//phi_pp /= (2.*Norder);
	//Phi_pp += phi_pp;
	//tPhiv(_idx_) = phi_pp*tPhi_dd;
	_idx_ += 1;
      }
    }
    */
    return ret;
  }
  complex<double> Dirspin_trial_move(int&tPPs, int ii, const Array<int,2>& diag, const Array<int,1>& dirspin,
				     const Array<int,1>& old_dirspin, const Array<double,1>& freq, ostream& xout)
  { //  The quantities from this class that are used, but not modified:
    //	   Array<int,2> diag
    //	   Array<int,1> dirspin, old_dirspin
    //     Array<double,1> freq
    //
    //     Cix          cix
    //	   Array<double,1> dlom                 // frequency on which hybridizations delta is specified
    //	   Array<double,2> Deltai               // imaginary part of hybridization Delta, as needed on real axis
    //     Array<double,1> xn                   // frequency on which PP are specified
    //     Array<complex<double>,2> Gm          // PP-green's function values
    //	   Spline1DLin<double> fep, lfep        // fermi function quantities on hybridization mesh
    // The quantities from this class that are modified:
    //     Array<int,2> wtb_t                   // The matrix of arguments, which uniquely defines frequency loops
    //     Array<double,1> Dlty_t(n)            // hybridization, changed only when dirspin for that propagator changed.
    //     tPhi_dd                              // contains the product of all hybridizations for particular diagram
    //
    //     Array<complex<double>,2> Gmt_t(tPPs,2*Norder) // all necessary pp propagators in the complex form (Re(G), widetilde{G})
    //
    // Note: if we change dirspin only, we can update wtb to -1*wtb_old for those wtb where dire has been changed.
    // The following quantities do not change:
    //    xfsgn_, fry_[:], lfpt_[:], sum_fp_, pp_freq_, pxs_, frms_
    //
    // Quantities that change:
    //    wtb_t
    //    Dlty_t for those iy on which dirspin has changed ( either wtb or spin).
    //    tPhi_dd has to be updated.
    //    tPhiv  has to be updated.
    //    Gmt_t is updated if pseudo is changed. Can see from old and new pseudo array
    //
    //const bool DEBUG=false;
    int Norder = diag.extent(0);
    //if (DEBUG) cout<<"Dirspin_trial_move for diag="<<diag<<" dirspin="<<dirspin<<endl;
    // Determins wtb(3*Norder,Nloops), which contains the signs for arguments of all propagators
    //DetermineLoops(wtb_t, diag, dirspin);
    wtb_t.resize(wtb.shape());
    Range range_loops(0,Norder);
    Range range_props(0,3*Norder-1);
    // first copy old wtb and than slightly update
    wtb_t(range_props,range_loops) = wtb(range_props,range_loops);
    for (int i=0; i<Norder; ++i){
      // is direction of hybridization changed?
      if ((dirspin(i)%2) != (old_dirspin(i)%2)) wtb_t(i,i+1) *= -1; // yes, hybridization direction changed.
    }
    //
    tPhi_dd = 1;
    for (int iy=0; iy<Norder; ++iy){ // 0:Norder for all possible configurations
      if (dirspin(iy)!=old_dirspin(iy)){
	int dire = dirspin(iy) % 2;    // direction of hybridization propagator
	int iorb = dirspin(iy) / 2;    // its spin
	int s = cix.bath_sigind(iorb); // actual entry for hybridization
	double y = freq(1+iy);         // frequency on hybridization
	if (wtb_t(iy,iy+1)>0){
	  intpar jpy = Interp(y, dlom);
	  Dlty_t(iy) = linInterp(Deltai, s, jpy); // actual interpolation and saving the interpolated hybridization
	}else{
	  intpar jmpy = Interp(-y, dlom);
	  Dlty_t(iy) = linInterp(Deltai, s, jmpy); // actual interpolation and saving the interpolated hybridization
	}
	tPhi_dd *= Dlty_t(iy);
      }else{
	tPhi_dd *= Dlty(iy);
      }
    }
    //idx_stored(ib,j) = idx 
    // Here is the problem why I did not implement to copy Gmt_ to Gmt_t. We would need to store index idx for previous
    // type diagram, because for the same ib and j Gmt_(idx,i) is not compatible with Gmt_t(idx,i) but with some Gmt_(idxp,i).
    // We would need to write Gmt_( idx_stored(ib,j), i) and we would need to store idx_stored to be compatible with stored Gmt_.
    // For now we don't store idx_stored, so we need to recompute every time.
    
    auto& tallowed = tdiag2matrix.Cmp_diagonal3(diag, dirspin, cix);      // map<int,Array<int,2>>
    tPPs=0;
    for (const auto& [ib, alw] : tallowed)  tPPs += alw.extent(0);
    Gmt_t.resize(tPPs,2*Norder); //Gmt_t = 0;
    if (New_form_Phi==2) phiR_t.resize(tPPs,2*Norder);
    if (New_form_Phi>=15) {mGppt_t.resize(tPPs);mGppt_t=1;Gmo_t.resize(tPPs,2*Norder);}
    int idx=0;
    for (const auto& [ib, alw] : tallowed){ // over all blocks that contribute on diagonal
      for (int j=0; j<alw.extent(0); ++j){  // in principle more than one term can contribute per block
	for (int i=0; i<2*Norder; ++i){ // going over all pseudo-particles
	  int _ii_ = alw(j,i);
	  int sg = isign(_ii_);
	  int ii = sg*_ii_ - 1;
	  //int ii = abs(_ii_)-1;
	  //int _ab_ = alw(i,j);
	  //int _a_ = abs(_ab_);
	  //int ab = cix.inBlock(_a_);      // _ab_ is given in terms of atomic state index, but ab is index within the block. cix.inBlock translates this.
	  //int _ib_ = cix.iBlck(_a_)-1;
	  //int ii = cix.Sgm_symmetry[_ib_](ab,ab)-1; // Green's function index for propagator G(a0,b0)==G[|_ii0_|-1]
	  if (ii>=0){
	    //double sg = sign(_ab_);
	    //double sg = sign(_ii_);
	    //Gmt_t(idx,i) = linInterp(Gm, ii, pxs(i)) * sg;  // now actual interpolation of Gm and seeting value of Gm*F. This is not just Gm, but multiplied by F.
	    Gmt_t(idx,i).linInterp(Gms, ii, pxs(i), sg);  // now actual interpolation of Gm and seeting value of Gm*F. This is not just Gm, but multiplied by F.
	    if (New_form_Phi==2) phiR_t(idx,i) = linInterp(PhiR, ii, pxs(i));
	    if (New_form_Phi>=15) Gmo_t (idx,i) = linInterp(Gm,   ii, pxs(i))*(sg+0.0);  // like for single orbital
	  }else{
	    Gmt_t(idx,i).set2zero();
	    if (New_form_Phi==2) phiR_t(idx,i) = 0;
	    if (New_form_Phi>=15) Gmo_t(idx,i) = 0;
	  }
	  if (New_form_Phi>=15) mGppt_t(idx) *= Gmo_t(idx,i);
	}
	idx += 1;
      }
    }
    tPhiv.resize(tallowed.size());
    complex<double> Phi_pp=0;
    idx=0; int _idx_=0;
    for (const auto& [ib, alw] : tallowed){ // over all blocks that contribute on diagonal
      complex<double> phi_pp=0;
      for (int j=0; j<alw.extent(0); ++j){ // in principle more than one term can contribute per block
	Gexp pd; pd.Init2Unity();
	for (int ip=0; ip<2*Norder; ++ip) pd *= Gmt_t(idx,ip);
	//
	if (New_form_Phi==1){
	  for (int ip=0; ip<2*Norder; ++ip)
	    phi_pp -= pd.get0(lfpt(ip), frms_(ip))/M_PI;
	}else if (New_form_Phi==2){
	  for (int ip=0; ip<2*Norder; ++ip){
	    Gexp Phis_over_Gm(pd);
	    Phis_over_Gm.divide_by_positive(Gmt_t(idx,ip));
	    Phis_over_Gm.matmul_lferm(lfpt(ip));
	    // divide by (-1)^p *G, where p is such that imaginary part of G is negative. That should work for diagonal components of G.
	    phi_pp += Phis_over_Gm.getImt(frms_(ip).lf)*phiR_t(idx,ip);
	  }
	}else if (New_form_Phi==10 or New_form_Phi==11){
	  double pd_im = pd.p*exp(pd.a), pd_abs = abs(pd.c);
	  phi_pp += pd_abs * sum_fp * (2*Norder);
	  //phi_pp += pd_im * sum_fp * (2*Norder);
	}else if (New_form_Phi==15){
	  phi_pp += mGppt_t(idx) * sum_fp * (2.0*Norder);
	}
	idx += 1;
      }
      phi_pp /= (2.*Norder);
      Phi_pp += phi_pp;
      tPhiv(_idx_) = phi_pp*tPhi_dd;
      _idx_ += 1;
    }
    return Phi_pp*tPhi_dd;
  }
  complex<double> FrequencyTime_trial_move(int itt, int iloop, const Array<int,2>& diag, const Array<int,1>& dirspin,
					   const Array<double,1>& tfreq, const Array<double,1>& freq, ostream& xout)
  { // Input Arguments:
    //	   Array<int,2>    diag
    //	   Array<int,1>    dirspin
    //     Array<double,1> tfreq, freq        // current and trial frequencies for all Norder+1 loops
    // Routine that computes the value of Phi (within certain approximation), which is used for MC sampling
    //
    // Updates the \widetilde{G} and Delta arrays when frequency changes.
    // Only modifies entries that are changed, to improve performance.
    // The quantities from this class that are used, but not modified:
    //  Cix                      cix
    //  map<int,Array<int,2>>    cdiag2matrix.diagonal   // which PP-indices give contributoon to Phi
    //  Array<int,2>             wtb                     // The matrix of arguments, which uniquely defines frequency loops for all propagators
    //  Array<complex<double>,2> Gm                      // all necessary pp propagators in the complex form (Re(G), widetilde{G})
    //  Array<double,2>          Deltai                  // imaginary part of hybridization Delta, as needed on real axis
    //  Array<double,1>          Dlty                    // hybridizations in this diagram
    //  Array<double,1>          pp_freq, xn             // contains pseudo-particle frequencies and pp-mesh
    //        
    // The quantities from this class that are modified:
    //
    //  Array<nFexp,1>           fry_t(iloop-1)          // contains fermi functions of hybridization frequencies f(y). Note always f(+y).
    //  Array<double,1>          Dlty_t(iloop-1)         // hybridizations in this diagram
    //  Array<complex<double>,2> Gmt_t  , but only if wtb(Norder+i,iloop)!=0
    //  Array<double,1>          lfpt_t                  // contains log(f(y1)*f(-y2)..*f(y_n)), the product needed for computing Sigma_j of pseudoparticle j.
    //  double                   sum_fp_t                // combination of fermi functions as needed for fake Phi
    //  Array<complex<double>,1> tPhiv(iblck)             // Resulting fake Phi functional being sampled, and only for each block
    //  Array<double,1>          pp_freq_t               // pseudo-particle frequencies 
    //  Array<intpar,1>          pxs_t                   // interpolation of pp propagators.
    //
    //int PPs = ppstates.extent(0);
    //tPhiv.resize(PPs);
    
    if (Gmt_t.shape()!=Gmt_.shape()){ Gmt_t.resize(Gmt_.shape()); } // We need Gmt_t to be of right size.
    if (New_form_Phi==2 and (phiR_t.shape()!=phiR.shape())){ phiR_t.resize(phiR.shape()); }
    if (New_form_Phi>=15 and (Gmo_t.shape()!=Gmo_.shape())){ Gmo_t.resize(Gmo_.shape()); }
    int Norder = diag.extent(0);
    current_iloop = iloop;
    double w_old = freq(iloop), w_new = tfreq(iloop);
    tPhi_dd = cPhi_dd;
    if (iloop>0){ // which hybridizations need to be updated. If iloop=0 it is omega loop,
      //                                                   that contains no hybridization
      int iy=iloop-1;               // which hybridization loop?
      int iorb = dirspin(iy) / 2;
      int s = cix.bath_sigind(iorb); // actual entry for hybridization
      intpar jpy = Interp(w_new, dlom);   // need intpar for f(w_new) or Delta(w_new)
      fry_t(iy).set(fep(jpy), lfep(jpy), beta*w_new);  // fermi function for f(y) stored as logarithm log(f(y))
      if (wtb(iy,iy+1)>0){
	Dlty_t(iy) = linInterp(Deltai, s, jpy); // actual interpolation and saving the interpolated hybridization
      }else{
	intpar jmpy = Interp(-w_new, dlom);
	Dlty_t(iy) = linInterp(Deltai, s, jmpy); // actual interpolation and saving the interpolated hybridization
      }
      tPhi_dd *= Dlty_t(iy)/Dlty(iy);
    }
    ////// computig fermi functions for pseudo frequencies f(pp_freq(i)) but here pp_freq(i) with updated frequencies
    //mGppt_t(Range(0,PPs-1)) = mGppt(Range(0,PPs-1));

    if (New_form_Phi>=15) { mGppt_t.resize(mGppt.shape());mGppt_t=mGppt;}
    //cout<<"previous mGppt="<<mGppt(Range(0,PPs-1))<<endl;
    for (int i=0; i<2*Norder; ++i){ // going over all pseudo-particles
      if (wtb(Norder+i,iloop)!=0){  // iloop is change, and this shows if PP(i) is affected. What is the frequency argument for pseudoparticle i?
	// the G[i] argument can be omega+y_0-y_2+y_3+.... Such combination would correspond to array wtb[Norder+i,:]=[1,1,0,-1,1,..]
	double x_new = pp_freq(i) + (w_new-w_old)*wtb(Norder+i,iloop);
	pp_freq_t(i) = x_new;
	intpar px = Interp(x_new, xn); // intpar for trial pseudo frequency
	pxs_t(i) = px;
	//
	frms_t(i).set( femx(px), lfemx(px), -beta*x_new); // fermi function f(-x) with argument x of pseudo propagator j.
	//
	int idx=0;
	for (const auto& [ib, alw] : cdiag2matrix.diagonal){ // over all blocks that contribute on diagonal
	  for (int j=0; j<alw.extent(0); ++j){  // in principle more than one term can contribute per block
	    int _ii_ = alw(j,i);
	    int sg = isign(_ii_);
	    int ii = sg*_ii_ - 1;
	    if (ii>=0){
	      Gmt_t(idx,i).linInterp(Gms, ii, px, sg);  // now actual interpolation of Gm and seeting value of Gm*F. This is not just Gm, but multiplied by F.
	      if (New_form_Phi==2) phiR_t(idx,i) = linInterp(PhiR, ii, px);
	      if (New_form_Phi>=15) Gmo_t (idx,i) = linInterp(Gm,  ii, px)*(sg+0.0);
	    }else{
	      Gmt_t(idx,i).set2zero();
	      if (New_form_Phi==2) phiR_t(idx,i) = 0;
	      if (New_form_Phi>=15) Gmo_t(idx,i) = 0;
	    }
	    if (New_form_Phi>=15){
	      mGppt_t(idx) *= Gmo_t(idx,i)/Gmo_(idx,i);
	      //cout<<" mGppt["<<idx<<"]="<<mGppt_t(idx)<<" Gmo_t["<<idx<<","<<i<<"]="<< Gmo_t(idx,i) << " Gmo="<<Gmo_(idx,i)<<endl;
	    }
	    idx += 1;
	  }
	}
      }
    }

    // This is the product of fermi functions needed for Self-energy Sigma_j of pseudoparticle j.
    lfpt_t(Range(0,2*Norder-1)) = lfpt(Range(0,2*Norder-1));
    // lfpt(ipseudo) = log(  \product_{iy} f(\pm_{ipseudo,iy} y)  )
    // here we changed fry_t(iy=iloop-1), hence we correct for f(\pm y_{iloop-1})/f(\pm y_{iloop-1})
    if (iloop>0){
      for (int j=0; j<2*Norder; ++j){
    	int iy=iloop-1;
    	if (xfsgn(j) & (1<<iy)){
    	  lfpt_t(j) += fry_t(iy).plf() - fry_(iy).plf(); //
    	}else{
    	  lfpt_t(j) += fry_t(iy).mlf() - fry_(iy).mlf(); //
    	}
      }
    }
    //cout<<"lfpt_t="<<lfpt_t(Range(0,2*Norder-1))<<" lfpt="<<lfpt(Range(0,2*Norder-1))<<endl;
    if (New_form_Phi>=10){
      sum_fp_t = 0;
      double lfp = 0;
      // sum_fp[ic] = \sum_{iy_ext}\sum_{ip}\product_{iy\ne iy_ext} f(\pm(ic,ip,iy) y)
      for (int iy=0; iy<Norder; ++iy){ // Now compute the product of fermi functions entering the physical green's function
	const nFexp& fry_iy = (iy==iloop-1) ? fry_t(iy) : fry_(iy);
	double max_fp=-10*beta;
	for (int ip=0; ip<2*Norder; ++ip){
	  if (xfsgn(ip) & (1<<iy)){
	    lfp = lfpt_t(ip) - fry_iy.plf() ;
	  }else{
	    lfp = lfpt_t(ip) - fry_iy.mlf() ;
	  }
	  Array<nFexp,1>& frms = (wtb(Norder+ip,iloop)!=0) ? frms_t : frms_;
	  if (New_form_Phi==10) lfp += frms(ip).mlf(); // f(+x_p)
	  if (lfp>max_fp) max_fp=lfp;
	}
	// instead of computing \sum_{ip} exp(lfp_{iy_ext,ip}) we just compute exp(max(lfp_{iy_ext,ip})), which is
	// a reasonable proxy, good enough for MC walk.
	sum_fp_t += exp(max_fp);
      }
    }
    //cout<<"sum_fp="<<sum_fp<<" mGppt_t="<<mGppt_t<<" tPhi_dd="<<tPhi_dd<<endl;

    tPhiv.resize(cdiag2matrix.diagonal.size());
    complex<double> Phi_pp=0;
    int idx=0; int _idx_=0;
    for (const auto& [ib, alw] : cdiag2matrix.diagonal){ // over all blocks that contribute on diagonal
      complex<double> phi_pp=0;
      for (int j=0; j<alw.extent(0); ++j){ // in principle more than one term can contribute per block
	Gexp pd; pd.Init2Unity();
	for (int ip=0; ip<2*Norder; ++ip) pd *= (wtb(Norder+ip,iloop)!=0) ? Gmt_t(idx,ip) : Gmt_(idx,ip);
	//
	//cout<<itt<<")  pd_ftm["<<_idx_<<"]="<< pd << endl;
	//
	if (New_form_Phi==1){
	  for (int ip=0; ip<2*Norder; ++ip){
	    Array<nFexp,1>& frms = (wtb(Norder+ip,iloop)!=0) ? frms_t : frms_;
	    phi_pp -= pd.get0(lfpt_t(ip), frms(ip))/M_PI;
	  }
	}else if (New_form_Phi==2){
	  for (int ip=0; ip<2*Norder; ++ip){
	    Gexp Phis_over_Gm(pd);
	    if (wtb(Norder+ip,iloop)!=0){
	      Phis_over_Gm.divide_by_positive( Gmt_t(idx,ip) );
	      Phis_over_Gm.matmul_lferm(lfpt_t(ip));
	      phi_pp += Phis_over_Gm.getImt(frms_t(ip).lf)*phiR_t(idx,ip);
	    }else{
	      Phis_over_Gm.divide_by_positive( Gmt_(idx,ip) );
	      Phis_over_Gm.matmul_lferm(lfpt_t(ip));
	      phi_pp += Phis_over_Gm.getImt(frms_(ip).lf)*phiR(idx,ip);
	    }
	  }
	}else if (New_form_Phi==10 or New_form_Phi==11){
	  double pd_im = pd.p*exp(pd.a), pd_abs = abs(pd.c);
	  phi_pp += pd_abs * sum_fp_t * (2*Norder);
	  //phi_pp += pd_im * sum_fp_t * (2*Norder);
	}else if (New_form_Phi==15){
	  phi_pp += mGppt_t(idx) * sum_fp_t * (2.0*Norder);
	  //cout<<" phi_pp["<<idx<<"]="<< phi_pp<< " mGppt="<<mGppt_t(idx)<<" sum_fp="<< sum_fp_t <<endl;
	}
	//cout<<itt<<")  phi_pp["<<_idx_<<"]="<< phi_pp << endl;
	idx+=1;
      }
      phi_pp /= (2.*Norder);
      Phi_pp += phi_pp;
      tPhiv(_idx_++) = phi_pp*tPhi_dd;
    }
    complex<double> ret = Phi_pp*tPhi_dd;
    return ret;
  }
  
  void Accept_topology_move(int ii, int PPs, int Norder)
  {
    wtb.resize(wtb_t.shape()); wtb = wtb_t;
    xfsgn.resize(xfsgn_t.size()); xfsgn = xfsgn_t;
    
    cPhi_dd = tPhi_dd;
    
    Range rngh(0,Norder-1);
    fry_(rngh) = fry_t(rngh);
    Dlty(rngh) = Dlty_t(rngh);
    
    Range rngp(0,2*Norder-1);
    lfpt(rngp) = lfpt_t(rngp);
    pp_freq(rngp) = pp_freq_t(rngp);
    frms_(rngp) = frms_t(rngp);
    pxs(rngp) = pxs_t(rngp);

    Gmt_.resize(Gmt_t.shape()); Gmt_ = Gmt_t;
    if (New_form_Phi==2) {phiR.resize(phiR_t.shape()); phiR = phiR_t;}
    else if (New_form_Phi>=10){ sum_fp = sum_fp_t;}
    if (New_form_Phi>=15){
      Gmo_.resize(Gmo_t.shape()); Gmo_(Range::all(),Range::all()) = Gmo_t(Range::all(),Range::all());
      mGppt.resize(mGppt_t.shape()); mGppt = mGppt_t;
    }
    cPhiv.resize(tPhiv.shape()); cPhiv = tPhiv;
    cdiag2matrix = tdiag2matrix;
    topology_change = true;
  }
  void Accept_dirspin_move(int ii, int PPs, int Norder, const Array<int,1>& dirspin, const Array<int,1>& old_dirspin)
  {
    if (wtb.shape()!=wtb_t.shape()){ cout<<"ERROR in Accept_dirspin_move wtb.shape=" << wtb.shape() << " while wbt_t.shape=" << wtb_t.shape() << endl; }
    wtb = wtb_t; // I think reshape is not needed, because wtp_t and wtb are of equal shape for Dirspin_trial_move
    
    cPhi_dd = tPhi_dd;
    
    for (int iy=0; iy<Norder; ++iy)
      if (dirspin(iy)!=old_dirspin(iy))
	Dlty(iy) = Dlty_t(iy);
    
    Gmt_.resize(Gmt_t.shape()); Gmt_ = Gmt_t;
    if (New_form_Phi==2) {phiR.resize(phiR_t.shape()); phiR = phiR_t;}
    if (New_form_Phi>=15){
      Gmo_.resize(Gmo_t.shape()); Gmo_ = Gmo_t;
      mGppt.resize(mGppt_t.shape()); mGppt = mGppt_t;
    }
    cPhiv.resize(tPhiv.shape()); cPhiv = tPhiv;
    cdiag2matrix = tdiag2matrix;
    topology_change = true;
  }
  void Accept_frequencyTime_move(int ii, int PPs, int Norder)
  {
    if (New_form_Phi>=10){ sum_fp = sum_fp_t;}
    if (current_iloop>0){
      fry_(current_iloop-1) = fry_t(current_iloop-1);
      Dlty(current_iloop-1)=Dlty_t(current_iloop-1);
      cPhi_dd = tPhi_dd;
      lfpt(Range(0,2*Norder-1)) = lfpt_t(Range(0,2*Norder-1));
    }
    //int PPs=0;
    //for (const auto& [ib, alw] : cdiag2matrix.diagonal) // over all blocks that contribute on diagonal
    //  PPs += alw.extent(1);  // in principle more than one term can contribute per block
    Range rngs(0,PPs-1);
    for (int i=0; i<2*Norder; ++i){
      if (wtb(Norder+i,current_iloop)!=0){
	pp_freq(i) = pp_freq_t(i);
	pxs(i) = pxs_t(i);
	frms_(i) = frms_t(i);
	Gmt_(rngs,i) = Gmt_t(rngs,i);
	if (New_form_Phi==2) phiR(rngs,i) = phiR_t(rngs,i);
	if (New_form_Phi>=15) Gmo_(rngs,i) = Gmo_t(rngs,i);
	frequency_changed(i)=true;
      }
    }
    if (New_form_Phi>=15){ mGppt(rngs) = mGppt_t(rngs); }
    cPhiv(Range::all()) = tPhiv(Range::all());
  }
  auto C(int k){ return Token{Kind::C, k}; };
  auto R(int k){ return Token{Kind::R, k}; };
  auto I(int k){ return Token{Kind::I, k}; };
  
  template<class RNG>
  void Meassure1(Array<double,1>& Phi_accumulate, Array<double,1>& Phi_accumulatef, Array<complex<double>,3>& Stval,
		 Array<complex<double>,3>& G00, Array<complex<double>,3>& F00, Array<complex<double>,3>& F01,
		 double& Phi_cOCA, double& Phi_cOCAf, Array<complex<double>,2>& Stval_cOCA,
		 Array<complex<double>,2>& G00_cOCA, Array<complex<double>,2>& F00_cOCA, Array<complex<double>,2>& F01_cOCA,
		 Array<double,1>& N_measure,
		 int itt, const complex<double>& cPhi, const Array<int,2>& diag,
		 const Array<int,1>& dirspin, const Array<double,1>& freq, int PPs, const Array<double,1>& weight,
		 const Array<double,1>& sqrtbetalist, RNG&, ostream& cout)
  {
    // The quantities from this class that are used, but not modified:
    /////// These should be updated at every move, and should be ready to be used here
    //
    // Array<int,2> wtb(Norder,Nloops)     -- frequency to loop connection
    // Array<double,1> sgn_tot(PPs)        -- total sign of a diagram
    // Array<double,1> pp_freq(2*Norder)   -- frequency on all pseudo propagators
    // double Phi_dd                           -- product of all hybridizations
    // Array<double,1> Dlty(Norder)        -- hybridizations for the current diagram
    // Array<intpar,1> pxs(2*Norder)       -- interpolation for pp_freq(j) for all PP propagators
    // Array<nFexp,1> fry_(Norder)         -- contains fermi functions of hybridization frequencies f(y). Note always f(+y).
    // Array<complex,2> Gmt_(PPs,2*Norder) -- all necessary PP propagators in the complex form (Re(G), \widetilde{G})
    // Array<double,1> lfpt(2*Norder)      -- contains log(f(y1)*f(-y2)..*f(y_n)), the product needed for computing Sigma_pp of pseudoparticle pp.
    // Array<nFexp,1> frms_(2*Norder)      -- contains fermi functions for pp-frequencies
    // Array<double,2> phiR(PPs,2*Norder)        -- the quantity for computing exact BK functional phi, not MC approximation
    // Array<int,1> n2two                  -- 2^n
    // double cutoffw
    // int wNbin
    //
    // The quantities computed here for measurement:
    //
    // Array<int,1> xfsgn(2*Norder)              -- contains signs for diagrams. This one should already be available, but we recompute it.
    // Array<int,1> xftype_conj(2*Norder)        -- contains conjugation information for physical green's function
    // Array<int,2> xftype_real(Norder,2*Norder) -- contains info if pp-propagator need to be real in green's function
    // Array<nGexp,2> uG_(PPs,2*Norder)          -- all necessary pp propagators (Re(G),e^log(-Im(G))) but in logarithmic form
    // Array<nGexp,1> uGpp(PPs)                  -- contains the product of all uG_ pp-propagators of this diagram in logarithmic representation
    // Array<complex,2> Gmc_(PPs,2*Norder)       -- all necessary pp propagators in the complex form (Re(G), -widetilde{G}f(-x))
    // Array<complex<double>,1> mGppc(PPs)       -- contains the product of pp-propagators Gmc_ in complex form (Re(G),Im{G})
    //
    const bool debug=false;
    int Norder = diag.extent(0);
    const int twoN = 2*Norder;
    int n_mino = Norder-N_min_order;
    N_measure(n_mino) += 1;
    string ccstr=" ";

    double phi_weight = weight(n_mino)/sqrtbetalist(n_mino);
    double W = abs(cPhi)/phi_weight;
    if (W==0){cout<<itt<<") ERROR cPhi="<<cPhi<<" W="<<W<<endl; exit(1);}
    
    if (topology_change){
      // computes which pseudo propagators are real, conjugated, and what are signs for fermi functions.
      DetermineFermiSigns(xfsgn,xftype_conj,xftype_real,wtb); // xfsgn was already known, but the rest are here computed for the first time.
      
      auto& allowem = cdiag2matrix.Cmp_matrix(diag, dirspin, cix);
      Lft.resize(2*Norder); // Array<Array<double,2>,2> Lft(:,:)(N1,N2)
      Rgh.resize(2*Norder); // Array<Array<double,2>,2> Rgh(:,:)(N1,N1)
      Gmx.resize(2*Norder);
      //for (int i=0; i<2*Norder; ++i) Gmx[i].resize(cix.Nmaxsize,cix.Nmaxsize);
      Fi.resize(2*Norder); // vector<Array<double,1>>  Fi[2*Norder](N1)
      //dS.resize(cix.Npseudo_unique,2*Norder);
      //
      sgns = SignSort(diag);    // sign due to topological crossing in the diagram
      if (cmpPhysG){
	if (not (subExprElim and Offdiagonal and Norder>2)){
	  subExpr=false;
	}else{
	  subExpr = true;
#ifdef _TIMEIT_     
	  auto t0 = get_cpu_time();
#endif	  
	  // products of (2*N-2) propagators like I(0)*R(1)*C(2)..., which means Im(G_0) * Re(G_1) * G_2
	  // we curve out two propagators, which are next to creation operator, which needs Njmm operator rather than F operator.
	  vector<vector<Token>> products(twoN*Norder, vector<Token>(twoN-2));// cproducts(twoN*Norder, vector<Token>(2));
	  // the two propagators adjacent to creation propagator iy fproduct[iy], which will need to be treated differently by inserting Njmm operator.
	  vector<vector<vector<Token>>> fproducts(Norder, vector<vector<Token>>(twoN, vector<Token>(2)));
	  for (int iy=0; iy<Norder; ++iy){ // over all hybridization lines
	    int dire = dirspin(iy)%2;      // direction of the hybridization
	    double y = (wtb(iy,iy+1)>0) ? freq(iy+1) : -freq(iy+1); // Delta(y) can enter with + or - into PP propagators.
	    int ibin = static_cast<int>((y/(2*cutoffw)+0.5)*wNbin_lmax); // which bin does this y fall into
	    if (ibin<0 or ibin>=wNbin_lmax) continue;
	    int ic = diag(iy,1-dire); // the end of hybridization line ie, we will use ie for vertex corrections. We could decide to use is, and equations would be slightly modified.
	    int ic_p_1 = (ic+1     )%twoN;
	    int ic_m_1 = (ic-1+twoN)%twoN;
	    for (int ip=0; ip<twoN; ++ip){ // this loop is over 2N different terms that we generate on the real axis. On Imaginary axis we would have only one term.
	      int X = xftype_real(iy,ip);  // bits of X tell us which Green's functions should be taken principle part only (real part) and which should be taken with imaginary part also
	      // First we will compute the common part, which is equal for G00 and F00. We need to exclude vertex ei and the two propagators comming toegether at this vertex ei.
	      // Below is plan to multiply (G[ic+1]*F[ic+2])* (G[ic+2]*F[ic+3])* ... (G[ic-2]*F[ic-1])
	      int l = ic_p_1;
	      for (int k=2; k<2*Norder; ++k){ // product[iy*twoN+ip][k-2] contains (G[l]*F[l+1]), which can be Im(), Re() or complex
		int l = (ic+k-1)%twoN;
		int lp1 = (l+1)%twoN;
		if (l==ip) products[iy*twoN+ip][k-2] = I(l); // // Im(G_l(x))*exp(-beta*x)
		else       products[iy*twoN+ip][k-2] = ( X & n2two(l) ) ? R(l) : C(l); // either Re(G_l(x)) or G_l(x)
	      }
	      {// Here we have only two propagators, ie-1 and ie, which are adjacent to the anhilation operator for hybridization iy
		l=ic_m_1;
		if (l==ip){
		  Token tok = I(l);
		  fproducts[iy][ip][0] = tok;
		}else{
		  Token tok = ( X & n2two(l) ) ? R(l) : C(l);
		  fproducts[iy][ip][0] = tok;//( X & n2two(l) ) ? R(l) : C(l);
		}
		l = ic;
		if (l==ip){
		  Token tok = I(l);
		  fproducts[iy][ip][1] = tok;
		}else{
		  Token tok = ( X & n2two(l) ) ? R(l) : C(l);
		  fproducts[iy][ip][1] = tok;//( X & n2two(l) ) ? R(l) : C(l);
		}
	      }
	    }
	  }
	  if (debug){
	    print_products(cout, products, to_string(itt)+") Original products");
	    for (int iy=0; iy<Norder; ++iy)
	      print_products(cout, fproducts[iy], to_string(itt)+") products iy="+to_string(iy));
	  }
	  int all_matrix_multiply=products.size()*(products[0].size()-1);
	  int second_multiply = Norder*twoN;
	  if (debug) cout<<"# matrix multiplications="<< all_matrix_multiply <<"  + "<< second_multiply << endl;
#ifdef _TIMEIT_     
	  auto t1 = get_cpu_time();
	  t_pre_plan += t1-t0;
#endif	  
	  CSEPlanner planner{products};
	  plan = planner.run();
	  
	  fplan.resize(Norder);
	  for (int iy=0; iy<Norder; ++iy){
	    CSEPlanner fplanner{fproducts[iy]};
	    fplan[iy] = fplanner.run();
	  }
#ifdef _TIMEIT_     
	  auto t2 = get_cpu_time();
	  t_plan += t2-t1;
#endif	  
	  if (debug){
	    cout << to_string(itt)+") Temporary sub-products found:\n";
	    int matrix_multiply=0;
	    for (auto& t : plan.temps){
	      cout << "H" << t.id << " = ";
	      for (size_t i=0;i<t.sig.size();++i){
		if (i) cout << " · ";
		cout << tok2str(t.sig[i]);
	      }
	      cout << "   used at: ";
	      for (auto& oc : t.occs) cout << "[P" << oc.prod << ":" << oc.start << "-" << oc.end << "] ";
	      cout << "\n";
	      matrix_multiply += t.sig.size()-1;
	    }
	    cout<<"matrix_multiplications for temps="<< matrix_multiply << endl;
	    // Print reduced products
	    cout << to_string(itt)<<") Reduced products:\n";
	    for (size_t p=0;p<plan.reduced.size();++p){
	      const auto& seq = plan.reduced[p];
	      cout << "P" << p << " = ";
	      for (size_t i=0;i<seq.size();++i){
		if (i) cout << " · ";
		cout << tok2str(seq[i]);
	      }
	      cout << "\n";
	      matrix_multiply += seq.size()-1;
	    }
	    int matrix_multiply2=0;
	    for (int iy=0; iy<Norder; ++iy){
	      cout << to_string(itt)+") sub-products iy="<<iy<<":\n";
	      for (auto& t : fplan[iy].temps){
		cout << "H" <<t.id << " = ";
		for (size_t i=0;i<t.sig.size();++i){
		  if (i) cout << " · ";
		  cout << tok2str(t.sig[i]);
		}
		cout << "   used at: ";
		for (auto& oc : t.occs) cout << "[P" << oc.prod << ":" << oc.start << "-" << oc.end << "] ";
		cout << "\n";
		matrix_multiply2 += t.sig.size()-1;
	      }
	      cout<<"matrix_multiplications for temps="<< matrix_multiply2 << endl;
	      // Print reduced products
	      cout << to_string(itt)<<") Second reduced products iy="<<iy<<":\n";
	      for (size_t p=0;p<fplan[iy].reduced.size();++p){
		const auto& seq = fplan[iy].reduced[p];
		cout << "P" << p << " = ";
		for (size_t i=0;i<seq.size();++i){
		  if (i) cout << " · ";
		  cout << tok2str(seq[i]);
		}
		cout << "\n";
		matrix_multiply2 += seq.size()-1;
	      }
	    }
	  }
	  //cout<<"total number of matrix multiplications ="<< matrix_multiply << " + "<< matrix_multiply2;
	  //cout<<" with max speedup="<< (second_multiply+all_matrix_multiply)/(matrix_multiply+matrix_multiply2+0.0) << endl;
#ifdef _TIMEIT_     
	  auto t3 = get_cpu_time();
	  t_post_plan += t3-t2;
#endif	  
	}
      }
    }
    // 1. Simplified code for not Offdiagonal
    // 2. Use products and fproducts for Offdiagonal, which should speed it up.
    //
    //
    Array<double,1> partial_phi(cdiag2matrix.allowed.size()); partial_phi=0;
    // Additional factors (appart from PP) for self-energy and Phi, i.e., sign and hybridization part, and MC-fake-weight and sign of the diagram.
    double extra_factors = (cPhi_dd/W)*sgns; // chi_dd==Hybridization part, W==Markov-chain part, sgns due to crossing propagators.
    int idx=0;
    double sum_phi=0;
    for (const auto& [ib, alw] : cdiag2matrix.allowed){
#ifdef _TIMEIT_     
      auto t0 = get_cpu_time();
#endif
      double dphi=0;
      if (! Offdiagonal){
#ifdef _TIMEIT_     
	auto t2 = get_cpu_time();
#endif	
	Gexp Phis_PP;
	Phis_PP.Init2Unity();
	for (int i=0; i<2*Norder; ++i){
	  int ip1 = (i+1)%(2*Norder);
	  Gmx[i].resize(1,1,true); // We could use Gmx[i].setActiveShape(N1,N2) and above use Gmx[i].resize(cix.Nmaxsize,cix.Nmaxsize); Will see if it is worth it.
	  Fi[i].resize(1);          // We store Fi for later.
	  int ib_a = cix.iBlck(abs(alw[i](0,1)))-1; // index of block, which should be =cix.iBlck(alw[ip1](j2,0))-1
	  int _a_ = alw[i](0,1);    // first index of the block G_{a,b}
	  double mm0 = sign(_a_);    // matrix element of op (is F or F^+) operator <alw[i](j1,1)|op|alw[i](j1,0)> is stored in the sign of op acting on index,
	  Fi[i](0)=mm0;             // Now remember matrix elements of F or F^+ for this operator at i. We choose this to be a diagonal matrix for convenience.
	  int ii = cix.Sgm_symmetry[ib_a](0,0)-1; // Green's function index for propagator G(a,b)==G[|_ii0_|-1]. This is due to our compact storage.
	  Gmx[i].linInterpolate(0,0, Gms, ii, pxs(i), 1); // Gmx[i](j1,j2) = ii_sign0*linInterpolate(Gms[ii](pxs(i))) : this is compatible with uG_(ipp,i)
	  Phis_PP *= Gmx[i];
	  Phis_PP *= mm0;
	}
	{
	  int ist=2*Norder-1;
	  dphi = -Phis_PP.get1(lfpt(ist), frms_(ist))/M_PI; // Gives exp(-beta*x) * exp(lfpt[ist]) * Im(_phi_pp_)
	}
#ifdef _TIMEIT_     
	auto t3 = get_cpu_time();
	t_dphi += t3-t2;
#endif
	for (int ist=0; ist<2*Norder; ++ist){
	  double x = pp_freq(ist);                 // what is frequency argument of this pseudo-particle
	  int ibin = static_cast<int>((x/(2.*cutoffw)+0.5)*wNbin_lmax); // Which bin in frequency?
	  int ip1 = (ist+1)%(2*Norder);
	  int ib_a = cix.iBlck(abs(alw[ist](0,1)))-1; // block index
	  Gexp Phis_over_Gm = Phis_PP;
	  Phis_over_Gm /= Gmx[ist]; // This is matrix equivalent of Phi/Gmx[i]
	  // Re(Phi_over_Gm_complex)*exp(lfpt(ist)), Im(Phi_over_Gm_complex)/f(-om)*exp(lfpt(ist))
	  complex<double> Phis_over_Gm_complex = Phis_over_Gm.getc1(lfpt(ist), frms_(ist).lf);
	  int ii = cix.Sgm_symmetry[ib_a](0,0)-1;
	  if (ibin>=0 and ibin<wNbin_lmax) Stval(n_mino,ii,ibin) += Phis_over_Gm_complex*extra_factors;
	}
#ifdef _TIMEIT_     
	auto t4 = get_cpu_time();//high_resolution_clock::now();
	t_stval += t4-t3;
#endif
      }else{
	int N0 = alw[0].extent(0);
	Lft[0].resize(N0,N0,true);
	for (int j=0; j<N0; ++j)  // Lft(:,0) is just set to F.
	  Lft[0].Set(j,j, sign(alw[0](j,1)));
	for (int i=0; i<2*Norder; ++i){
	  int ip1 = (i+1)%(2*Norder);
	  int N1 = alw[i].extent(0), N2 = alw[ip1].extent(0);
	  mGexp GxF(N1,N2);          // This will contain G_i @ F_{i+1}
	  Gmx[i].resize(N1,N2,true); // We could use Gmx[i].setActiveShape(N1,N2) and above use Gmx[i].resize(cix.Nmaxsize,cix.Nmaxsize); Will see if it is worth it.
	  Fi[i].resize(N1);          // We store Fi for later.
	  int ib_a = cix.iBlck(abs(alw[i](0,1)))-1; // index of block, which should be =cix.iBlck(alw[ip1](j2,0))-1
	  for (int j1=0; j1<N1; ++j1){ // over the matrix for this block
	    int _a_ = alw[i](j1,1);    // first index of the block G_{a,b}
	    int a = cix.inBlock(abs(_a_)); // alw[i] index (_a_) is given in terms of atomic state index, but a is index within the block. cix.inBlock translates this.
	    double mm0 = sign(_a_);    // matrix element of op (is F or F^+) operator <alw[i](j1,1)|op|alw[i](j1,0)> is stored in the sign of op acting on index,
	    // hence sign of alw[i](j1,1) contains op matrix element.
	    Fi[i](j1)=mm0;             // Now remember matrix elements of F or F^+ for this operator at i. We choose this to be a diagonal matrix for convenience.
	    // In our implementation G is full matrix and F is diagonal (with matrix elements either +1,-1 or 0). 
	    for (int j2=0; j2<N2; ++j2){// over the second index of the matrix for this block
	      int _b_ = alw[ip1](j2,0);// second index of the block G_{a,b}
	      int b = cix.inBlock(abs(_b_));// again needs to be translated to in-block index from global index.
	      //int ibb = cix.iBlck(_b_)-1;// should be equal to ib_a.
	      int _ii0_ = cix.Sgm_symmetry[ib_a](a,b); // Green's function index for propagator G(a,b)==G[|_ii0_|-1]. This is due to our compact storage.
	      if (_ii0_==0) continue;     // if index vanishes, G(a,b) exactly vanishes.
	      int ii_sign0 = sign(_ii0_); // Green's function value for off-diagonal propagators could be negative or opposite, but otherwise equal in magnitude. 
	      int ii = abs(_ii0_)-1;      // index of where G(a,b) is stored.
	      Gmx[i].linInterpolate(j1,j2, Gms, ii, pxs(i), ii_sign0); // Gmx[i](j1,j2) = ii_sign0*linInterpolate(Gms[ii](pxs(i))) : this is compatible with uG_(ipp,i)
	      if (i<2*Norder-1){                                       // this is equivalent to Gmx[i] @ diag(F[i+1])
		int mm1 = sign(alw[ip1](j2,1)); // matrix element: <alw[i+1](j2,1)|op|alw[i+1](j2,0)>. It should be equivalent to F[(i+1)%(2N)]
		GxF.Set(j1,j2, Gmx[i], mm1);    // G_i*F_{i+1}, i.e., GxF(j1,j2) = Gmx[i](j1,j2)*mm1
	      }
	    }
	  }
	  if (i<2*Norder-1){ // Lft[i] = F[0] @ Gmx[0]@F[1] ... Gmx[i-1]@F[i]
	    // Now we compute Lft[i+1] = Lft[i] @ G[i]@F[i+1]
	    Lft[i+1].matmul( Lft[i], GxF);// this is similar to uGpp(ipp)
	  }
	}
#ifdef _TIMEIT_     
	auto t1 = get_cpu_time();
	t_lft += t1-t0;
#endif
	int N1 = alw[2*Norder-1].extent(0);
	// Rgh[i] = F[2N-i] @ Gmx[2N-i] @ .... @ Fi[2N-1] @ Gmx[2N-1]
	Rgh[1].resize(N1,N0);
	Rgh[1].dot_diag(Fi[2*Norder-1], Gmx[2*Norder-1]); // Rgh(idx,1) = diag(Fi[2*n-1]) @ Gmx[2*n-1]
	for (int i=1; i<2*Norder-1; ++i){
	  int N1 = alw[2*Norder-1-i].extent(0), N2 = alw[2*Norder-i].extent(0);
	  mGexp FGx(N1,N2); // This will contain F_i @ G_i
	  FGx.dot_diag(Fi[2*Norder-1-i], Gmx[2*Norder-1-i]); // FGx = diag(Fi[2*n-1-i]) @ Gmx[2*n-1-i]
	  Rgh[i+1].matmul(FGx, Rgh[i]); // Rgh[i+1] = FGx @ Rgh[i];
	}
#ifdef _TIMEIT_     
	auto t2 = get_cpu_time();
	t_rgh += t2-t1;
#endif
	// Here we do the trace of multiplication from the left.
	// There are 2*N possible ways to compute that and they should be equal when Compatible_interp=true.
	// When Compatible_interp=false, they should be close, but not identical. I wonder if we should use
	// at least two ways to compute it from the left+right?
	dphi=0;
	if (not check_equivalence){
	  int ist=2*Norder-1;
	  Gexp _phi_pp_ = trace_dot_Gexp(Lft[ist], Gmx[ist]);// _phi_pp_=Tr(F[0]@Gmx[0]@F[1] ... Gmx[2N-2]@F[2N-1]@Gmx[2N-1])
	  dphi = -_phi_pp_.get1(lfpt(ist), frms_(ist))/M_PI; // Gives exp(-beta*x) * exp(lfpt[ist]) * Im(_phi_pp_)
	  // where exp(-beta*x)  is computed from frms_[ist]=f(-x)
	}
#ifdef _TIMEIT_     
	auto t3 = get_cpu_time();
	t_dphi += t3-t2;
#endif
	//cout<<itt<<") "<<idx<<": check=[";
	for (int ist=0; ist<2*Norder; ++ist){
	  double x = pp_freq(ist);                 // what is frequency argument of this pseudo-particle
	  int ibin = static_cast<int>((x/(2.*cutoffw)+0.5)*wNbin_lmax); // Which bin in frequency?
	  int ip1 = (ist+1)%(2*Norder);
	  int N1 = alw[ist].extent(0), N2 = alw[ip1].extent(0); // dimension of G(a,b)
	  int ib_a = cix.iBlck(abs(alw[ist](0,1)))-1; // block index
	  mGexp Phis_over_Gm(N2,N1); // This is matrix equivalent of Phi/Gmx[i]
	  if (ist<2*Norder-1)
	    // Phi/Gmx[i] = Rgh[2N-1-i]@Lft[i] = F[i+1] @ Gmx[i+1] @ .... @ Fi[2N-1] @ Gmx[2N-1] @ F[0] @ Gmx[0]@F[1] ... Gmx[i-1]@F[i]
	    Phis_over_Gm.matmul(Rgh[2*Norder-1-ist], Lft[ist] ); // Right*Left gives the self-energy contribution, i.e., roughly Phi/G
	  else{
	    // Phi/Gmx[2N-1] = Lft[2N-1] = F[0] @ Gmx[0]@F[1] ... Gmx[2N-2]@F[2N-1]
	    Phis_over_Gm = Lft[ist];
	  }
	  if (check_equivalence){
	    // It turns out that phi_pp_ should be independent of ist in exact numerics. Because of linear interpolation it might not be.... Check!!!
	    // I just proved that exp(lfpt(ist)+frms_(ist).x)=const, hence we don't need lfpt and frms_ at the same time.
	    // We could write lfpt(i) = lfpt(0)+frms_(0).x-frms_(i).x
	    Gexp _phi_pp_ = trace_dot_Gexp(Phis_over_Gm, Gmx[ist]);
	    // PP-contribution to Phi
	    dphi += -_phi_pp_.get1(lfpt(ist), frms_(ist))/(M_PI*2*Norder); // Gives exp(-beta*x) * exp(lfpt[ist]) * Im(_phi_pp_),
	    // where exp(-beta*x)  is computed from frms_[ist]=f(-x)
	  }
	  Array<complex<double>,2> Phis_over_Gm_complex(N2,N1); // extract Phis_over_Gm_complex = Phis_over_Gm*exp(lfpt[ist])
	  // But actually it gives widetilde{Sigma} in imaginary part. More specifically, it computes
	  // Re(Phis_over_Gm_complex) = Re(Phis_over_Gm)*exp(lfpt[ist])
	  // Im(Phis_over_Gm_complex) = -Im(Phis_over_Gm)*exp(lfpt[ist])/f(-x), where f(-x)~frms_[ist]
	  Phis_over_Gm.getc1(Phis_over_Gm_complex, lfpt(ist), frms_(ist).lf); // The imaginary part is this the tilde component
	  Phis_over_Gm_complex *= extra_factors; 
	  for (int j1=0; j1<N1; ++j1){ // Now we store it into Stval. Over the first index of the matrix
	    int _a_ = alw[ist](j1,1);  // component a of Sigma_{b,a}
	    int a = cix.inBlock(abs(_a_)); // alw[i] index is given in terms of global atomic index, but a is index within the block.
	    for (int j2=0; j2<N2; ++j2){// over the second index of the matrix Sigma_{b,a}
	      int _b_ = alw[ip1](j2,0); // second index b in Sigma_{b,a} or G_{a,b}
	      int b = cix.inBlock(abs(_b_)); // alw[i+1] is also given in terms of global atomic index, now translate to within the block
	      int _ii0_ = cix.Sgm_symmetry[ib_a](b,a); // Green's function index for propagator Sigma(b,a)==Sigma[|_ii0_|-1] ~ G(a,b)
	      if (_ii0_==0) continue;
	      double ii_sign0 = sign(_ii0_); // Green's function value for off-diagonal propagators could be negative/opposite, hence sign
	      int ii = abs(_ii0_)-1;         // actual index in the G/Sigma array
	      complex<double> dSt = ii_sign0*Phis_over_Gm_complex(j2,j1); // need to store opposite sign when the component is opposite
	      if (ibin>=0 and ibin<wNbin_lmax){
		Stval(n_mino,ii,ibin) += dSt; // ii is PP index, ibin is pp-frequency index, and n_mino is order-minimal order.
	      }
	    }
	  }
	}
#ifdef _TIMEIT_     
	auto t4 = get_cpu_time();
	t_stval += t4-t3;
#endif
      } // Finished Offdiagonal
      //cout<<"]"<<endl;
      dphi *= extra_factors;
      sum_phi += dphi;
      partial_phi(idx) = dphi;
      Phi_accumulate(n_mino)  += dphi;
      if (cmpPhysG) Phi_accumulatef(n_mino) += dphi;
      
      if (cmpPhysG){
	Array<complex<double>,2> GGin; // brisi: temporary to comparison
	Array<complex<double>,2> GGi, GGj, GGl; // views on storage_GG
	Array<complex<double>,2> A, B, C; // also views on other storage
	Array<double,2> GGk;
	Array<bool,1> needed(cix.Nmaxsize);
	if (subExpr){
#ifdef _TIMEIT_     
	  auto t0 = get_cpu_time();
#endif
	  // We divide backbone into two parts, 2N-2 part encoded in plan, and additional 2 propagators
	  // that depend on incoming hybridization fplan[iy]
	  // First we compute all temporary products needed below for 2N-2 part of the backbone.
	  // They are stored into plan.temps[:].sig[:]
	  H.resize(plan.temps.size());
	  for (const auto& t : plan.temps) {
	    if (debug){
	      cout <<itt<< ") ib="<<ib<<" H" << t.id << " = ";
	      for (size_t i=0;i<t.sig.size();++i){
		if (i) cout << " · ";
		cout << tok2str(t.sig[i]);
	      }
	      cout << endl;
	    }
	    const int L = static_cast<int>(t.sig.size());
	    if (L == 0) continue;
	    Give_GtimesFMatrix(A,  // A into storage_GGi or H[previous], it prepares t.sig[0] component from Gmx[l]*Fi
			       t.sig[0], Gmx, Fi, H, alw, storage_GGi, frms_, twoN); 
	    if (L == 1) { // We are done, no actual multiplication, just store into t.id
	      H[t.id].resize(A.extent(0), A.extent(1));
	      H[t.id] = A;
	      if (debug) cout<<"  H=A="<<A<<endl;
	      continue;
	    }
	    if (debug) cout<<"  "<<tok2str(t.sig[0])<<"="<<A<<endl;
	    //cout<<"  note that corresponding C["<<t.sig[0].idx<<"]="<<Gmx[t.sig[0].idx].c<<endl;
	    //  Here is the plan of how we are going to compute the chain product:
	    // il=1         il=2        il=3         il=4
	    // A -> GGi|H   A->C->GGl   A->C->GGi    A->C->GGl
	    // flip-false   flip=true   flip=false   flip=true
	    // B -> GGj|H   B -> GGj|H  B -> GGj|H   B -> GGj|H
	    // C -> GGl	    C -> GGi    C -> GGl     C -> GGi
	    // C = A @ B    C = A @ B   C = A @ B    C = A @ B
	    //               
	    bool flip = false;
	    for (int il = 1; il < L; ++il) { // Below we chain the remaining factors
	      Give_GtimesFMatrix(B,  // B saved into storage_GGj or take H[previous]
			   t.sig[il], Gmx, Fi, H, alw, storage_GGj, frms_, twoN);
	      const bool is_last = (il == L - 1);
	      const int m = A.extent(0);
	      const int n = B.extent(1);
	      if (is_last) {
		H[t.id].resize(m, n); // Write final result into H[t.id]
		C.reference(H[t.id]); // final result is written into C but actually into H[t.id]
	      } else { // if not last, write into the ping-pong scratch buffer
		Array<complex<double>,2>& target = flip ? storage_GGi : storage_GGl;
		C.reference(target(Range(0, m-1), Range(0, n-1)));
	      }
	      tiny::dot(C, A, B); // C = A @ B
	      if (debug){
		cout<<"  A="<<A<<endl;
		cout<<"  "<<tok2str(t.sig[il])<<"="<<B<<endl;
		cout<<"  C="<<C<<endl;
	      }
	      // Prepare for next step: A becomes C.
	      A.reference(C); // no real copying, just taking reference
	      flip = !flip;
	    }
	    if (debug) cout<<"  H["<<t.id<<"]="<< H[t.id] << endl;
	  }
#ifdef _TIMEIT_     
	  auto t1 = get_cpu_time();
	  t_H += t1-t0;
#endif
	}
	// We precomputed common products for all hybridizations. Now we cut one of the hybridizations iy
	for (int iy=0; iy<Norder; ++iy){ // over all hybridization lines
#ifdef _TIMEIT_     
	  auto t0 = get_cpu_time();
#endif
	  int dire = dirspin(iy)%2;
	  int iorb = dirspin(iy)/2;
	  int s = cix.bath_sigind(iorb); // actual entry for hybridization
	  //
	  double y = (wtb(iy,iy+1)>0) ? freq(iy+1) : -freq(iy+1);
	  int ibin = static_cast<int>((y/(2*cutoffw)+0.5)*wNbin_lmax); // which bin does this y fall into
	  if (ibin<0 or ibin>=wNbin_lmax) continue;
	  //
	  double Dda = cPhi_dd/Dlty(iy);      // product of all hybridizations divided by this hybridization
	  //
	  //                --->---
	  //              /	    \
	  //  pseud[ie]  |           | 
	  //      ic==ie |----<------| is
	  //             |           | pseud[is]
	  //  pseud[ie-1] \         /  
	  //                ---<---
	  //             F_R = psi^\dagger_s(tau_{ie}) n_{1-s}(\tau_{ie}^{-}) -> n_{1-s}(pseud[ie-1]) * <pseud[ie]|psi^\dagger_s|pseud[ie-1]>
	  //             F_L = psi_s(\tau_{is}) n_{1-s}(\tau_{is}^{+})        -> n_{1-s}(pseud[is])   * <pseud[is]|psi_s|pseu[is-1]>
	  int is = diag(iy,dire);   // the begining of hybridization line
	  int ie = diag(iy,1-dire); // the end of hybridization line
	  int ic = ie;              // we will use ie for vertex corrections. We could decide to use is, and equations would be slightly modified.
	  int ic_m_1 = (ic-1+twoN)%twoN;
	  int ic_p_1 = (ic+1     )%twoN;
	  int ibm = cix.iBlck(abs(alw[ie](0,0)))-1;// a bit complicated way to find which block corresponds to this ist
	  int ibp = cix.iBlck(abs(alw[ie](0,1)))-1;// a bit complicated way to find which block corresponds to this ist
	  Array<double,2> njmm;           // matrix when njm is not proportional to F
	  //int njmc = cix.Njm_c(ibm,iorb); // index to identify whether njm is proportional to F or not.
	  NjmResult nx = cix.Njm_ind(ibm,2*iorb);
	  double njmr = nx.njm_r;         // constant when njm is proportional to F
	  //cout<<"njmc="<<njmc<<" at iy="<<iy<<" ibm="<<ibm<<" iorb="<<iorb<<endl;
	  //if (njmc==1){                   // njm is proportional to F, hence simple
	  //  njmr = cix.Njm_r.at(cix.encode_block_bath(ibm,iorb));
	  //}else
	  if (nx.njm_c==2){             // njm is not proportional
	    // Need to compute Trace(Phi/(Gm[ie-1]*F*Gm[ie]) * Gtm * njmm * Gtp)
	    // where the sizes are:
	    //     Gtm(alw[ic_m_1].extent(0),cix.Block[ibm].size());
	    //     njmm(cix.Block[ibm].size(),cix.Block[ibp].size());
	    //     Gtp(cix.Block[ibp].size(),alw[ic_p_1].extent(0));
	    // The fundamental issue is that Fi zeros out many entries when computing PP self-energy or physical Green's function.
	    // Only those components that are in alw[] are allowed, and not the entire block.
	    // However, Njm matrix, which is used for calculation of F00, is nonzero in the entire block, hence we need
	    // to extent PP-Green's functions beyond what is normally needed for self-enerhy and physical G calculation.
	    // This interpolation is carried out below.
	    //
	    //njmm.reference(cix.Njm[ibm][iorb]);
	    // WARNING: This could be done withouth copying njmm. Of not transposed, just take reference. If transposed, use dot with transposed.
	    if (not nx.transposed){
	      njmm.resize(nx.njm->extent(0),nx.njm->extent(1));
	      njmm = *nx.njm;
	    }else{
	      njmm.resize(nx.njm->extent(1),nx.njm->extent(0));
	      for (int n1=0; n1<njmm.extent(0); ++n1)
		for (int n2=0; n2<njmm.extent(1); ++n2)
		  njmm(n1,n2) = (*nx.njm)(n2,n1);
	    }
 	    //cout<<" ibm="<<ibm<<" iorb="<<iorb<<" Njm="<<njmm<<"="<<cix.Njm[ibm][iorb]<<endl;
	    {
	      // Below we interpolate Gtm, which is G-PP[ie-1] before hybridizaton ends, and is needed for vertex F00.
	      int im=ic_m_1, ip=ic; // This propagator is between the vertices with indices im,ip. The index of propagator is im by our choice.
	      int N1 = alw[im].extent(0), N2 = alw[ip].extent(0); // sizes are of the im and ip vertex size. Note that the index arrays are alw[im](:,1) and alw[ip](:,0)
	      // which is obvious from the vertex perspective, first index is the second index of the previous vertex, and second index is the first index of the next vertex.
	      int Nm = cix.Block[ibm].size();
	      if (N2==Nm){// No need to create larger matrix and interpolate, alw index contains all components
		Gtm.reference(Gmx[im]); // note that reference only points Gtm.c&p to Gmx[imp].c&p but Gtm.storage_c&p remain independent.
	      }else { // Need to increase Gtm and add new interpolates
		Gtm.resize(N1,Nm,true); // Need to resize and set it to unity
		Gtm.a = Gmx[im].a; // setting exponent directly to Gmx[ie-1]. The entire block has a single exponent
		needed(Range(0,Nm-1))=true; // we will copy and remember which were copied and dont need to be interpolated.
		for (int j2=0; j2<N2; ++j2){ // Those values are available in Gmx[im]. Just copy them from c,p arrays.
		  int b = cix.inBlock(abs(alw[ip](j2,0))); // j2 is index in alw[ip], and alw[ip](j2,0) is ending index of propagator, but specified in global atomic list
		  // because alw[ip](j2,0) is given in global settings, cix.inBlock gives index in this block, which turns out to be block ibm.
		  Gtm.c(Range::all(),b) = Gmx[im].c(Range::all(),j2);
		  Gtm.p(Range::all(),b) = Gmx[im].p(Range::all(),j2);
		  needed(b)=false;
		}
		for (int b=0; b<Nm; ++b){ // now we go over all indices in the block ibm, and interpolate those that have not been interpolated yet.
		  if (needed(b)){  // those were not copied, hence need to interpolate
		    for (int j1=0; j1<N1; ++j1){
		      int a = cix.inBlock(abs(alw[im](j1,1)));
		      int _ii0_ = cix.Sgm_symmetry[ibm](a,b);
		      if (_ii0_!=0){
			int ii_sign0 = sign(_ii0_);
			int ii = abs(_ii0_)-1;
			Gtm.linInterpolate(j1,b, Gms, ii, pxs(im), ii_sign0); // Gtm(j1,b) = ii_sign0*linInterpolate(Gms[ii](pxs(im))) : this is compatible with uG_(ipp,i)
		      }
		    }
		  }
		}
	      }
	    }
	    {
	      // Below we interpolate Gtp=G-PP[ic], similarly to Gtm, it is G-PP just after hybridization ends.
	      int im = ic, ip=ic_p_1; // This propagator is between the vertices im,ip
	      int N1 = alw[im].extent(0), N2 = alw[ip].extent(0); // again the sizes of blocks, which are needed for pseudo self-energy. Note than Nmm needs the entire block, not
	      // just what is allowed by Fi and used in self-energy calculation.
	      int Np = cix.Block[ibp].size();
	      if (Np==N1){// No need to create larger matrix and interpolate, because alw contains the entire block
		Gtp.reference(Gmx[im]);
	      }else{ // Need to increase Gtp and add new interpolates
		Gtp.resize(Np,N2,true);
		Gtp.a = Gmx[im].a; // should have the same exponent as Gmx[ie]
		needed(Range(0,Np-1))=true; // those that can't be copied
		for (int j1=0; j1<N1; ++j1){
		  int a = cix.inBlock(abs(alw[im](j1,1)));
		  Gtp.c(a,Range::all()) = Gmx[im].c(j1,Range::all()); // copy what can be copied
		  Gtp.p(a,Range::all()) = Gmx[im].p(j1,Range::all()); 
		  needed(a)=false; // once copied, no need to interpolate
		}
		for (int a=0; a<Np; ++a){
		  if (needed(a)){ // those were not copied, hence need to interpolate
		    for (int j2=0; j2<N2; ++j2){
		      int b = cix.inBlock(abs(alw[ip](j2,0)));
		      int _ii0_ = cix.Sgm_symmetry[ibp](a,b);
		      if (_ii0_!=0){
			int ii_sign0 = sign(_ii0_);
			int ii = abs(_ii0_)-1;
			Gtp.linInterpolate(a,j2, Gms, ii, pxs(im), ii_sign0);
		      }
		    }
		  }
		}
	      }
	    }
	  }
#ifdef _TIMEIT_     
	  auto t1 = get_cpu_time();
	  t_interp += t1-t0;
#endif
	  // We now have Gtm and Gtp covering the entire block around vertex ie of hybridization iy.
	  // Next we compute all temporary products from common for Gmx[ie-1]*F*Gmx[ie], which are
	  // stored in gH, and also Gtm*Njm*Gtp products, which are stored in fH.
	  // These products are listed in fplan[iy].temps
	  /////////////////NEW//////////////////////
	  if (subExpr and fplan[iy].temps.size()>0){
#ifdef _TIMEIT_     
	    auto t0 = get_cpu_time();
#endif
	    gH.resize(fplan[iy].temps.size()); // This one is for G_physical but not for F_physical, unless Njm is proportional for F
	    for (const auto& t : fplan[iy].temps) { // t is a temporary product stored in gH[t.id]
	      if (debug){
		cout <<itt<< ") ib="<<ib<<" iy="<<iy<<" gH" << t.id << " = ";
		for (size_t i=0;i<t.sig.size();++i){
		  if (i) cout << " · ";
		  cout << tok2str(t.sig[i]);
		}
		cout << endl;
	      }
	      const int L = static_cast<int>(t.sig.size()); // t has L term-product.
	      if (L != 2){ cout<<"ERROR expect exactly two matrices to multiply in this temp. Have L="<<L<<endl;}
	      Give_GtimesFMatrix(A, // A in GGi_storage or gH[previous]: it computes Gmx[ie-1]*Fi[ie] or gH[previous]
				 t.sig[0], Gmx, Fi, gH, alw, storage_GGi, frms_, twoN); 
	      Give_GMatrix(B, // B stored in GGj_storage of gH[previous]: it computes Gmx[ie] but does not contain Fi[ie+1]
			   t.sig[1], Gmx, gH, alw, storage_GGj, frms_, twoN); // must WRITE INTO storage_GGj or H
	      const int m = A.extent(0);
	      const int n = B.extent(1);
	      gH[t.id].resize(m, n); // Write final result directly into H[t.id] (no extra copy)
	      tiny::dot(gH[t.id], A, B); // gH[t.id] = A @ B
	      if (debug){
		cout<<"  "<<tok2str(t.sig[0])<<"="<<A<<endl;
		cout<<"  "<<tok2str(t.sig[1])<<"="<<B<<endl;
		cout<<"  gH["<<t.id<<"]="<< gH[t.id] << endl;
	      }
	    }
#ifdef _TIMEIT_     
	    auto t1 = get_cpu_time();
	    t_Hg += t1-t0;
#endif
	    if (nx.njm_c==2){ // Need also fH[] = Gm*Njm*Gp
	      fH.resize(fplan[iy].temps.size());
	      for (const auto& t : fplan[iy].temps) { // the same products needed, but using the entire blocks around Njm, instead of F operators.
		if (debug){
		  cout <<itt<< ") ib="<<ib<<" fH" << t.id << " = ";
		  for (size_t i=0;i<t.sig.size();++i){
		    if (i) cout << " · ";
		    cout << tok2str(t.sig[i]);
		  }
		  cout << endl;
		}
		const int L = static_cast<int>(t.sig.size());
		if (L != 2){ cout<<"ERROR expect exactly two matrices to multiply in this temp. Have L="<<L<<endl;}
		Give_GMatrix_only(A, // A stored into storage_GGi: contains Gtm, or Re(Gtm) or Im(Gtm)*exp(-beta*x)
				  t.sig[0], Gtm, fH, storage_GGi, frms_, twoN); // must WRITE INTO storage_GGi or H
		Give_GMatrix_only(B, // B stored into storage_GGj: contains Gtp or Re(Gtp) or Im(Gtp)*exp(-beta*x)
				  t.sig[1], Gtp, fH, storage_GGj, frms_, twoN); // must WRITE INTO storage_GGj or H
		// WARNING: here we should use transposed instead of copy njmm
		const int m = A.extent(0);
		const int k = njmm.extent(1);
		const int n = B.extent(1);
		C.reference( storage_GGl(Range(0,m-1),Range(0,k-1)) ); // C->storage_GGl
		tiny::dot(C,A,njmm);       // C = A @ njmm = Gtm @ njmm
		fH[t.id].resize(m, n);     // Write final result directly into H[t.id] (no extra copy)
		tiny::dot(fH[t.id], C, B); // fH[] = C @ B = Gtm @ njmm @ Gtp
		if (debug){
		  cout<<"  "<<tok2str(t.sig[0])<<"="<<A<<endl;
		  cout<<"  njmm="<<njmm<<endl;
		  cout<<"  "<<tok2str(t.sig[1])<<"="<<B<<endl;
		  cout<<"  fH["<<t.id<<"]="<< fH[t.id] << endl;
		}
	      }
	    }
#ifdef _TIMEIT_     
	    auto t2 = get_cpu_time();
	    t_Hf += t2-t1;
#endif
	  }// Now we have fH and gH temporary products, below we used them for all 2N products
	  /////////////////NEW//////////////////////
	  complex<double> Ac=0, Af=0;
	  for (int ip=0; ip<2*Norder; ++ip){ // this loop is over 2N different terms that we generate on the real axis. On Imaginary axis we would have only one term.
#ifdef _TIMEIT_     
	    auto ta = get_cpu_time();
#endif
	    int X = xftype_real(iy,ip);      // bits of X tell us which Green's functions should be taken principle part only (real part) and which should be taken with imaginary part also
	    complex<double> ggi=0, ffi=0;    // contribution to G00 and F00
	    //complex<double> _ggi_=0, _ffi_=0;
	    // First we compute the common part, which is equal for G00 and F00. We need to exclude vertex ei and the two propagators comming toegether at this vertex ei.
	    // Below is multiplication of F[ic+1]* (G[ic+1]*F[ic+2])* (G[ic+2]*F[ic+3])* ... (G[ic-2]*F[ic-1])
	    int l = ic_p_1;// l=(ic+1)%twoN;
	    if (not Offdiagonal){ // Ising or single band model type
	      if (not (nx.njm_c==1 or nx.njm_c==0)) cout<<"ERROR for Diagonal system we should have njmc==1 or 0, but njmc="<<nx.njm_c<<endl;
	      ggi=1.0;
	      for (int l=0; l<twoN; ++l){// G[0]*F[1] * G[1]*F[2] * G[2]*F[3] ... G[2N-1]*F[0]
		int lp1 = (l+1)%twoN;
		if (l==ip) ggi *= Gmx[l].Set_Gt_ferm_1D(frms_(l).x);     // gives -Im(G_l(x))*exp(-beta*x)=p*exp(a-beta*x)=p*exp(a+mf.x)
		else       ggi *= ( X & n2two(l) ) ?  real(Gmx[l].c(0,0)) : Gmx[l].c(0,0); // we need to take principle (real) part of G.
		ggi *= Fi[lp1](0);
	      }
	      ffi = ggi*njmr;
	    }else{ // Offdiagonal with matrix form.
	      if (subExpr){ // use above temporary products to compute all products
#ifdef _TIMEIT_     
		auto t0 = get_cpu_time();
#endif
		Array<complex<double>,2> A, B, C;
		const auto& seq = plan.reduced[iy*twoN+ip]; // all products when using temporary products stored in H
		int L = static_cast<int>(seq.size());
		if (debug){
		  cout <<itt<< ") P[" << iy*twoN+ip << "] = ";
		  for (size_t i=0;i<seq.size();++i){
		    if (i) cout << " · ";
		    cout << tok2str(seq[i]);
		  }
		  cout<<endl;
		}
		if (L == 0) continue;
		Give_GtimesFMatrix(A, // A=seq[0] stored into storage_GGi, Gmx[l]*Fi or H[], depending on seq[0]
			     seq[0], Gmx, Fi, H, alw, storage_GGi, frms_, twoN); // must WRITE INTO storage_GGi or H
		if (L == 1) {
		  GGi.reference(A); // GGi is result, and can point to storage_GGi or H[], previously used by A.
		  if (debug) cout<<"  GGi=A="<<GGi<<endl;
		  continue;
		}
		if (debug) cout<<"  "<<tok2str(seq[0])<<"="<<A<<endl;
		//
		// il=1         il=2        il=3         il=4
		// A -> GGi|H   A->C->GGl   A->C->GGi    A->C->GGl
		// flip-false	  flip=true   flip=false   flip=true
		// B -> GGj|H	  B -> GGj|H  B -> GGj|H   B -> GGj|H
		// C -> GGl	  C -> GGi    C -> GGl     C -> GGi
		// C = A @ B    C = A @ B   C = A @ B    C = A @ B
		//               
		bool flip = false;
		// Computing the remaining factors in the product seq[1]...seq[L-1]
		for (int il = 1; il < L; ++il) {
		  Give_GtimesFMatrix(B, // B=seq[il] stored into storage_GGj, Gmx[l]*Fi or H[], depending on seq[il]
			       seq[il], Gmx, Fi, H, alw, storage_GGj, frms_, twoN); // must WRITE INTO storage_GGj or H
		  const bool is_last = (il == L - 1);
		  const int m = A.extent(0);
		  const int n = B.extent(1);
		  Array<complex<double>,2>& target = flip ? storage_GGi : storage_GGl;
		  C.reference(target(Range(0, m-1), Range(0, n-1))); // write either into storage_GGi or storage_GGl
		  tiny::dot(C, A, B); // C = A @ B
		  if (debug){
		    cout<<"  A="<<A<<endl;
		    cout<<"  "<<tok2str(seq[il])<<"="<<B<<endl;
		    cout<<"  C="<<C<<endl;
		  }
		  if (is_last){
		    GGi.reference(storage_GGi(Range(0,m-1),Range(0,n-1))); // last always goes into storage_GGi and is called GGi
		    if (not flip) GGi = C; // copying data when C points to storage_GGl, and not storage_GGi.
		  }else{
		    // 4) Prepare for next step: A becomes C.
		    A.reference(C);
		  }
		  flip = !flip;
		}
		if (debug) cout<<"  GGi["<<iy<<","<<ip<<"]="<<GGi<<endl;
		tiny::dot_diag(Fi[ic_p_1],GGi); // Now we add F[ic+1] at the very beginning.
		// temporary storing into GGin because GGi is likely rewriten with temporary use below.
		GGin.resize(GGi.extent(0), GGi.extent(1));
		GGin = GGi;
		if (debug) cout<<"  F*GGi=GGin["<<iy<<","<<ip<<"]="<<GGin<<endl;
		//
		// Next we compute Gmx[ie-1]*F*Gmx[ie] products using previously prepared gH temporary products.
		const auto& ser = fplan[iy].reduced[ip]; // simplified products using tempe gH.
		L = static_cast<int>(ser.size()); // how long is this product?
		if (!(L==1 or L==2) ){ cout<<"ERROR expected two matrices or single matrix. Have L="<<L<<endl;}
		if (debug){
		  cout <<itt<< ") Px[" << ip << "] = ";
		  for (size_t i=0;i<ser.size();++i){
		    if (i) cout << " · ";
		    cout << tok2str(ser[i]);
		  }
		  cout<<endl;
		}
#ifdef _TIMEIT_     
		auto t1 = get_cpu_time();
		t_GGin += t1-t0;
#endif
		//cout<<"survived L="<<L<<endl; cout.flush();
		if (L==1){ // this product is just equal to temporary product stored in gH[tok.idx]
		  Token tok = ser[0];
		  ggi = tiny::trace_dot(gH[tok.idx],GGin);            // ggi = Tr( GGl @ GGin )
		  if (debug) cout<<"using gH["<<tok.idx<<"]="<<gH[tok.idx]<<endl;
		  //// This part should be removed when finished, because its only used for comparison below.
		  //GGl1.resize(gH[tok.idx].extent(0),gH[tok.idx].extent(1));
		  //GGl1 = gH[tok.idx];
		  //// end of temporary code
		  if (nx.njm_c==1){
		    ffi = ggi*njmr;
		  }else if (nx.njm_c==2){
		    ffi = tiny::trace_dot(fH[tok.idx],GGin);            // ffi = Tr( GGl @ GGi )
		  }else{
		    ffi = 0;
		  }
		}else if (L==2){ // We need to compute the product Gmx[ie-1]@Fi@Gmx[ie] here.
		  Give_GtimesFMatrix(A, // A=ser[0] from Gmx[ser[0].idx]*Fi stored into storage_GGi or gH (can ruing GGi hence we needed GGin)
				     ser[0], Gmx, Fi, gH, alw, storage_GGi, frms_, twoN); // must WRITE INTO storage_GGi or H
		  Give_GMatrix_only(B, // B=ser[1] from Gmx[ser[1].idx] stored into storage_GGj
				    ser[1], Gmx[ser[1].idx], gH, storage_GGj, frms_, twoN); // must WRITE INTO storage_GGj or H
		  int m = A.extent(0);
		  int n = B.extent(1);
		  GGl.reference(storage_GGl(Range(0,m-1),Range(0,n-1))); // result stored in GGl -> storage_GGl
		  tiny::dot(GGl, A, B); // GGl = A @ B = Gmx[ie-1]*Fi[ie]*Gmx[ie]
		  ggi = tiny::trace_dot(GGl, GGin);         // ggi = Tr( GGl @ GGin )
		  if (debug){
		    cout<<"  "<<tok2str(ser[0])<<"="<<A<<endl;
		    cout<<"  "<<tok2str(ser[1])<<"="<<B<<endl;
		    cout<<"  Px["<<ip<<"]="<< GGl << endl;
		    cout.flush();
		  }
		  if (nx.njm_c==1){
		    ffi = ggi*njmr;
		  } else if (nx.njm_c==2){ // Njm is not proportional to Fi, hence compute Gtm@njm@Gtp
		    Give_GMatrix_only(A, // A stored into storage_GGl, contains Gtm or Re(Gtm) or Im(Gtm)*exp(-beta*x)
				      ser[0], Gtm, fH, storage_GGl, frms_, twoN); // must WRITE INTO storage_GGi or H
		    Give_GMatrix_only(B, // B stored into storage_GGj, contains Gtp or Re(Gtp) or Im(Gtp)*exp(-beta*x)
				      ser[1], Gtp, fH, storage_GGj, frms_, twoN); // must WRITE INTO storage_GGj or H
		    m = A.extent(0);
		    int k = njmm.extent(1);
		    n = B.extent(1);
		    C.reference( storage_GGi(Range(0,m-1),Range(0,k-1)) ); // Result writen into storage_GGi
		    tiny::dot(C,A,njmm); // C = A @ njmm = Gtm @ njm -> storage_GGi
		    if (debug){
		      cout<<"  "<<tok2str(ser[0])<<"="<<A<<endl;
		      cout<<"  njmm="<<njmm<<endl;
		      cout<<"  "<<tok2str(ser[0])<<"*njmm="<<C<<endl;
		      cout<<"  "<<tok2str(ser[1])<<"="<<B<<endl;
		    }
		    GGl.reference(storage_GGl(Range(0,m-1),Range(0,n-1))); // Now result goes in GGl -> storage_GGl
		    tiny::dot(GGl, C, B); // GGl = C @ B = Gtm @ njm @ Gtp
		    ffi = tiny::trace_dot(GGl, GGin);        // ffi = Tr( GGl @ GGin )
		    if (debug){
		      cout<<"  Px["<<ip<<"]="<< GGl << endl;
		      cout.flush();
		    }
		  } else{
		    ffi = 0;
		  }
		}else{
		  cout<<"ERROR unexpected L="<<L<<endl;
		}
#ifdef _TIMEIT_     
		auto t2 = get_cpu_time();
		t_GGln += t2-t1;
#endif
	      }else{
		//if (subExpr and debug){
		//cout<<itt<<") ib="<<ib<<" old calculation Norder="<<Norder<<endl;
		//}
#ifdef _TIMEIT_     
		auto t3 = get_cpu_time();
#endif
		
		int N0 = alw[l].extent(0);
		for (int k=2; k<2*Norder; ++k){
		  // Here we do GGi <- GGi * (G[ic+k-1]*F[ic+k]). at k==2 we start by GGi=identity
		  int l = (ic+k-1)%twoN;
		  int lp1 = (l+1)%twoN;
		  int N1 = alw[l].extent(0), N2 = alw[lp1].extent(0);
		  // We compute GGj <- GGi * { Im(Gmx[l].c)*exp(-beta*x) | Re(Gmx[l].c) | Gmx[l].c }
		  GGj.reference( storage_GGj(Range(0,N0-1),Range(0,N2-1)) );
		  if (l==ip){
		    // One Green's function has the form Im(G(x))*exp(-beta*x), which means opposite sign of the fermi function f(x)*\widetilde{G}(x)=Amt(x)
		    // where Im(G)= -f(-x)*\widetilde{G}(x)
		    if (k>2){
		      GGk.reference( storage_GGk(Range(0,N1-1),Range(0,N2-1)) );		
		      Gmx[l].Set_Gt_ferm(GGk, frms_(l).x);     // gives -Im(G_l(x))*exp(-beta*x)=p*exp(a-beta*x)=p*exp(a+mf.x)
		      tiny::dot(GGj, GGi, GGk);                // GGj = GGi @ GGk
		      if (subExpr and debug){
			cout<<"  A="<<GGi<<endl;
			cout<<"  I["<<l<<"]="<<GGk<<endl;
			cout<<"  C="<<GGj<<endl;
		      }
		    }else{
		      Gmx[l].Set_Gt_ferm(GGj, frms_(l).x); // GGj = GGk = -Im(G_l(x))*exp(-beta*x)
		      if (subExpr and debug) cout<<"  I["<<l<<"]="<<GGj<<endl;
		    }
		  }else{
		    if ( X & n2two(l) ){ // we need to take principle (real) part of G.
		      if (k>2){
			GGk.reference( storage_GGk(Range(0,N1-1),Range(0,N2-1)) );		
			GGk = real(Gmx[l].c);     // GGk = Re(G_l(x))
			tiny::dot(GGj, GGi, GGk); // GGj = GGi @ GGk
			if (subExpr and debug){
			  cout<<"  A="<<GGi<<endl;
			  cout<<"  R["<<l<<"]="<<GGk<<endl;
			  //cout<<"  note that corresponding C["<<l<<"]="<<Gmx[l].c<<endl;
			  cout<<"  C="<<GGj<<endl;
			}
		      } else{
			GGj = real(Gmx[l].c);// GGj = real(G_l(x))
			if (subExpr and debug){
			  cout<<"  R["<<l<<"]="<<GGj<<endl;
			}
		      }
		    }else{
		      if (k>2){
			tiny::dot(GGj, GGi, Gmx[l].c);// GGj = GGi @ G_l(x)
			if (subExpr and debug){
			  cout<<"  A="<<GGi<<endl;
			  cout<<"  C["<<l<<"]="<<Gmx[l].c<<endl;
			  cout<<"  C="<<GGj<<endl;
			}
		      }
		      else{
			GGj = Gmx[l].c;                   // GGj = G_l(x)
			if (subExpr and debug) cout<<"  C["<<l<<"]="<<GGj<<endl;
		      }
		    }
		  }
		  // GGi <- GGj
		  GGi.reference( storage_GGi(Range(0,N0-1),Range(0,N2-1)) );
		  GGi = GGj;
		  tiny::dot_diag(GGi,Fi[lp1]); // GGj *= Fi[l+1], or, GGj *= F[ic+k]. Note, the last is Fi[ic-1]
		  if (subExpr and debug) cout<<"GGi["<<iy<<","<<ip<<"]*F="<<GGi<<endl;
		}
		if (Norder>1){
		  //cout<<"GGi.shape="<<GGi.shape()<<" Fi["<<ic_p_1<<"].shape="<<Fi[ic_p_1].shape()<<" Fi="<<Fi[ic_p_1]<<endl;
		  tiny::dot_diag(Fi[ic_p_1],GGi); // Now we add F[ic+1] at the very beginning.
		}else{ // in this case GGi is empty
		  GGi.reference( storage_GGi(Range(0,N0-1),Range(0,N0-1)) ); // Final result is not in GGi
		  tiny::cdiag(GGi,Fi[ic_p_1]);
		  if (subExpr and debug) cout<<"F*GGi["<<iy<<","<<ip<<"]="<<GGi<<endl;
		}
#ifdef _TIMEIT_     
		auto t4 = get_cpu_time();
		t_GGi += t4-t3;
#endif
		// GGi now contains: F[ic+1]*G[ic+1]*F[ic+2]*G[ic+2]*F[ic+3] ... G[ic-2]*F[ic-1]
		//
		// Finally we need to add GGl==G[ic-1]*F[ic]*G[ic], which is different for G00 and F00.
		// We continue at ie-1 to get G[ie-1]*F[ie]*G[ie]. For F00 we need G[ie-1]*Njmm*G[ie] and Njmm has larger dimension than F[ie].
		GGl.reference( storage_GGl(Range(0,alw[ic_m_1].extent(0)-1),Range(0,alw[ic_p_1].extent(0)-1)) );
		l=ic_m_1;
		int N1 = alw[ic_m_1].extent(0), N2 = alw[ic].extent(0);
		GGj.reference( storage_GGj(Range(0,N1-1),Range(0,N2-1)) );
		if (l==ip){
		  Gmx[l].Set_Gt_ferm(GGj, frms_(l).x);      // GGj = -Im(G_l(x))*exp(-beta*x)
		}else{
		  if ( X & n2two(l) ) GGj = real(Gmx[l].c);// GGj = Re(G_l(x))
		  else GGj = Gmx[l].c;                     // GGj = G_l(x)
		}
		tiny::dot_diag(GGj,Fi[ic]);                // GGj = GGj @ Fi[ic]
		l = ic;
		N1 = alw[ic].extent(0); N2 = alw[ic_p_1].extent(0);
		if (l==ip){
		  GGk.reference( storage_GGk(Range(0,N1-1),Range(0,N2-1)) );		
		  Gmx[l].Set_Gt_ferm(GGk, frms_(l).x);     // GGk = -Im(G_l(x))*exp(-beta*x)
		  tiny::dot(GGl, GGj, GGk);                // GGl = GGj @ GGk
		}else{
		  if ( X & n2two(l) ){
		    GGk.reference( storage_GGk(Range(0,N1-1),Range(0,N2-1)) );		
		    GGk = real(Gmx[l].c);                  // GGk = Re(G_l(x))
		    tiny::dot(GGl, GGj, GGk);              // GGl = GGj @ GGk
		  }else tiny::dot(GGl, GGj, Gmx[l].c);     // GGl = GGj @ G_l(x)
		}
		ggi = tiny::trace_dot(GGl,GGi);            // ggi = Tr( GGl @ GGi )	      
		if (nx.njm_c==1){
		  ffi = ggi*njmr;
		}else if (nx.njm_c==2){
		  // Now we calculate GGl==Gtm*Njmm[ic]*Gtp, where Njmm is full matrix.
		  int l=ic_m_1; // back to l=(ic-1)%twoN
		  GGl.reference( storage_GGl(Range(0,Gtm.N1-1),Range(0,Gtm.N2-1)) );
		  if (l==ip){
		    Gtm.Set_Gt_ferm(GGl, frms_(l).x);      // GGl = -Im(Gtm(x))*exp(-beta*x)
		  }else{
		    if (X & n2two(l)) GGl = real(Gtm.c);   // GGl = Re(Gtm(x))
		    else GGl = Gtm.c;                      // GGl = Gtm(x)
		  }
		  GGj.reference( storage_GGj(Range(0,Gtm.N1-1),Range(0,njmm.extent(1)-1)) );
		  tiny::dot(GGj,GGl,njmm);                 // GGj = GGl @ njmm
		  GGl.reference( storage_GGl(Range(0,Gtm.N1-1),Range(0,Gtp.N2-1)) );
		  l = ic;
		  if (l==ip){
		    GGk.reference( storage_GGk(Range(0,Gtp.N1-1),Range(0,Gtp.N2-1)) );
		    Gtp.Set_Gt_ferm(GGk, frms_(l).x);      // GGk = -Im(Gtp(x))*exp(-beta*x)
		    tiny::dot(GGl, GGj, GGk);              // GGl = GGj @ GGk
		  }else{
		    if ( X & n2two(l) ){
		      GGk.reference( storage_GGk(Range(0,Gtp.N1-1),Range(0,Gtp.N2-1)) );		
		      GGk = real(Gtp.c);                   // GGk = Re(Gtp(x))
		      tiny::dot(GGl, GGj, GGk);            // GGl = GGj @ GGk
		    }else tiny::dot(GGl, GGj, Gtp.c);      // GGk = GGj @ Gtp(x)
		  }
		  ffi = tiny::trace_dot(GGl,GGi);          // ffi = Tr( GGl @ GGi )
		}
#ifdef _TIMEIT_     
		auto t5 = get_cpu_time();
		t_GGl += t5-t4;
#endif
	      }
	    }
	    // in this case the entire term must be conjugated. This is the rule for the real axis.
	    if ( xftype_conj(iy) & n2two(ip) ){
	      ggi = -conj(ggi);
	      ffi = -conj(ffi);
	    }
	    // Fermi functions except the one that we took functional derivative
	    double lfp = lfpt(ip);    // all fermi functions for ip self-energy
	    if (xfsgn(ip) & (1<<iy))
	      lfp -= fry_(iy).plf();  // we divide by one fermi function, which corresponds to the external frequency
	    else
	      lfp -= fry_(iy).mlf();  // we divide by one fermi function, which corresponds to the external frequency
	    double _fp_ = exp(lfp);   // we were dealing with log(f....)
	    Ac += ggi*_fp_;
	    Af += ffi*_fp_;
#ifdef _TIMEIT_     
	    auto tb = get_cpu_time();
	    t_GF1 += tb-ta;
#endif
	  }
	  bool destroys_Cooper_pair = cix.Fi_indx(2*iorb)%2==1 and cix.Fi_indx(2*iorb+1)%2==1;
	  bool normal_down = cix.SC_bath(iorb) and ((cix.Fi_indx(2*iorb)+cix.Fi_indx(2*iorb+1))%2==1 and cix.bath_spin(iorb)<0);
	  bool destroys_Cooper_pair_or_normal_down = destroys_Cooper_pair or normal_down;
	  
	  G00(n_mino,s,ibin)  += Ac*(Dda*sgns/M_PI)/W; // contribution to physical Greens function
	  complex<double> Fval = Af*(Dda*sgns/M_PI)/W; // contribution to physical Greens function
	  if (destroys_Cooper_pair_or_normal_down)
	    F01(n_mino,s,ibin) += Fval; // contribution to F vertex function
	  else
	    F00(n_mino,s,ibin) += Fval; // contribution to F vertex function

#ifdef _TIMEIT_     
	  auto t6 = get_cpu_time();
	  t_GF2 += t6-t1;
#endif
	}
      }
      idx += 1; // Finally we proceed to the next block
    } // end of block [ib,alw] loop
    frequency_changed=false; // now all frequencies are updated.
    topology_change=false;
  }
  void Resize_Output(Array<double,2>& Sm, Array<complex<double>,2>& G00, Array<complex<double>,2>& F00, Array<complex<double>,2>& F01,
		     const Cix& cix, bool cmpPhysG)
  {
    //int Npseudo_unique = deg_pp.extent(0);
    //int Nbath_unique =  deg_bath.extent(0);
    Sm.resize(cix.Npseudo_unique,xn.size());
    if (cmpPhysG){
      G00.resize(cix.Nbath_unique,dlom.size());
      F00.resize(cix.Nbath_unique,dlom.size());
      F01.resize(cix.Nbath_unique,dlom.size());
    }
  }
  int size(int Norder){ return Norder+1;}
  
  double discrete_factor(){return wNbin_lmax/(2*cutoffw);}

  void Broaden_Quantities(Array<complex<double>,3>& Stval, Array<complex<double>,3>& G00, Array<complex<double>,3>& F00, Array<complex<double>,3>& F01,
			  bool cmpPhysG, const json& pars, MPIcls& mpi)
  {
    Array<double,1> xt(wNbin_lmax); 
    for (int i=0; i<wNbin_lmax; ++i) xt(i) = cutoffw*( -1.0 + (2*i+1.0)/wNbin_lmax );
    //(*mpi.out)<<"wNbin="<<wNbin_lmax<<" xt.size="<<xt.size()<<" Stval.shape="<<Stval.shape()<<endl;
    double Wtypical = pars["Wtypical"];
    double broad = pars["broad"];
    bool fastFilesystem = pars["fastFilesystem"];
    int Norders = N_max_order-N_min_order;
    //(*mpi.out)<<"Wtypical="<<Wtypical<<" broad="<<broad <<" fastFilesystem="<<fastFilesystem<<" Norders="<<Norders<<endl;
    
    //(*mpi.out)<<"Before broaden sum(Stval)="<<sum(Stval)<<endl;
    if (fastFilesystem and mpi.rank==mpi.master){
      for (int iorder=0; iorder<Stval.extent(0); ++iorder){
	ofstream hout("Stval_."s+to_string(iorder+N_min_order)); hout.precision(12);
	ofstream rout("Strel_."s+to_string(iorder+N_min_order)); hout.precision(12);
	for (int i=0; i<xt.size(); ++i){
	  hout<<xt(i)<<" ";
	  rout<<xt(i)<<" ";
	  for (int ip=0; ip<Stval.extent(1); ++ip) hout<<Stval(iorder,ip,i).imag()<<" ";
	  for (int ip=0; ip<Stval.extent(1); ++ip) rout<<Stval(iorder,ip,i).real()<<" ";
	  hout<<endl;
	  rout<<endl;
	}
	if (cmpPhysG){
	  ofstream gout("g00_."s+to_string(iorder+N_min_order)); gout.precision(12);
	  for (int i=0; i<xt.size(); ++i){
	    gout<<xt(i)<<" ";
	    for (int ip=0; ip<G00.extent(1); ++ip) gout<<G00(iorder,ip,i).real()<<" "<<G00(iorder,ip,i).imag()<<" ";
	    for (int ip=0; ip<G00.extent(1); ++ip) gout<<F00(iorder,ip,i).real()<<" "<<F00(iorder,ip,i).imag()<<" ";
	    for (int ip=0; ip<G00.extent(1); ++ip) gout<<F01(iorder,ip,i).real()<<" "<<F01(iorder,ip,i).imag()<<" ";
	    gout<<endl;
	  }
	}
      }
    }
    double exclude = pars["brd_exclude"];
    double max_der = pars["brd_max_der"];
    double artificial_beta = pars["brd_beta"];
    Array<double,2> Sm(Stval.extent(1),Stval.extent(2));
    for (int iorder=0; iorder<Stval.extent(0); ++iorder){
      Sm(Range::all(),Range::all()) = imag(Stval(iorder,Range::all(),Range::all()));
      double md = FindMaxDeriv(xt, Sm, Wtypical);
      //if (md<0.2) md=0.2;
      int n = iorder+N_min_order;
      double ky = 0.015*md*broad*n*n/16.;  // might need some finetuning, but the width of gaussian will increase as sqrt(kwidth)*om
      double kwidth= ky/(1.+ky);  // it should never be >1 because integral does not converge 
      //double kwidth=0.005*FindMaxDeriv(xt, Sm, Wtypical)*broad; // WARNING: here I reduced from 0.01 to 0.005
      double _exclude_ = exclude;
      if (iorder<=1) _exclude_ = 1000; // the lowest two orders should not be cut off at all.
      Array<double,2> x_remove = FindRangeToRemove(xt, Sm, _exclude_, max_der);
      (*mpi.out)<<" broadening kwidth["<<N_min_order+iorder<<"]="<<kwidth<<" remove=";
      for (int j=0; j<Stval.extent(1); ++j){
	(*mpi.out)<<"(";
	if (x_remove(j,0)>xt(0)) (*mpi.out)<<x_remove(j,0);
	(*mpi.out)<<",";
	if (x_remove(j,1)<xt(xt.extent(0)-1)) (*mpi.out)<<x_remove(j,1);
	(*mpi.out)<<") ";
      }
      (*mpi.out)<<endl;
      //
      //Stval(iorder,Range::all(),Range::all()) = BroadenFunctions(xt,Stval(iorder,Range::all(),Range::all()),kwidth,*mpi.out);
      // Make a contiguous slab
      Array<complex<double>,2> slab(Stval.extent(1), Stval.extent(2));
      slab = Stval(iorder, Range::all(), Range::all());
      slab = BroadenFunctions(xt, slab, kwidth, *mpi.out);
      // Write it back to the view
      Stval(iorder, Range::all(), Range::all()) = slab;
      //
      //(*mpi.out)<<"During broaden1 iorder="<<iorder<<"broaden sum(Stval)="<<sum(Stval(iorder,Range::all(),Range::all()))<<" sum="<<sum(Stval)<<endl;
      //
      for (int j=0; j<Stval.extent(1); ++j){
	if (x_remove(j,1)<xt(xt.extent(0)-1)){
	  double x_start = x_remove(j,1)-1.0/artificial_beta;
	  for (int i=0; i<xt.extent(0); ++i)
	    Stval(iorder,j,i) *= fermi((xt(i)-x_start)*artificial_beta);
	}
	if (x_remove(j,0)>xt(0)){
	  double x_start = x_remove(j,0)+1.0/artificial_beta;
	  for (int i=0; i<xt.extent(0); ++i)
	    Stval(iorder,j,i) *= 1-fermi((xt(i)-x_start)*artificial_beta);
	}
      }
      /*
      for (int j=0; j<Stval.extent(1); ++j)
	for (int i=0; i<xt.extent(0); ++i)
	  if (std::isnan(Stval(iorder,j,i))){
	    (*mpi.out)<<"WARNING At iorder="<<iorder<<" ip="<<j<<" xt["<<i<<"]="<<xt(i)<<" Stval is nan. Correcting to 0."<<endl;
	    Stval(iorder,j,i)=0;
	  }
      (*mpi.out)<<"During broaden2 iorder="<<iorder<<"broaden sum(Stval)="<<sum(Stval(iorder,Range::all(),Range::all()))<<" sum="<<sum(Stval)<<endl;
      */
    }
    //(*mpi.out)<<"After broaden sum(Stval)="<<sum(Stval)<<endl;
    if (cmpPhysG){
      int Nbath_unique = G00.extent(1);
      Array<double,2> zg00(6*Nbath_unique,xt.size());
      for (int iorder=0; iorder<=Norders; ++iorder){
	zg00(Range(             0,  Nbath_unique-1),Range::all()) = real(G00(iorder,Range::all(),Range::all()));
	zg00(Range(  Nbath_unique,2*Nbath_unique-1),Range::all()) = imag(G00(iorder,Range::all(),Range::all()));
	zg00(Range(2*Nbath_unique,3*Nbath_unique-1),Range::all()) = real(F00(iorder,Range::all(),Range::all()));
	zg00(Range(3*Nbath_unique,4*Nbath_unique-1),Range::all()) = imag(F00(iorder,Range::all(),Range::all()));
	zg00(Range(4*Nbath_unique,5*Nbath_unique-1),Range::all()) = real(F01(iorder,Range::all(),Range::all()));
	zg00(Range(5*Nbath_unique,6*Nbath_unique-1),Range::all()) = imag(F01(iorder,Range::all(),Range::all()));
	//double kwidth=0.01*FindMaxDeriv(xt, zg00, Wtypical)*broad;
	double md = FindMaxDeriv(xt, zg00, Wtypical);
	int n = iorder+N_min_order;
	double ky = 0.015*md*broad*n*n/16.;
	double kwidth = ky/(1.+ky);
	double _exclude_ = exclude;
	if (iorder<=0) _exclude_ = 1000; // the lowest orders should not be cut off at all.
	Array<double,2> x_remove = FindRangeToRemove(xt, zg00, _exclude_, max_der);
	Array<double,2> y_remove(Nbath_unique,2);
	for (int ib=0; ib<Nbath_unique; ++ib){
	  y_remove(ib,0) = -1e100;
	  y_remove(ib,1) = 1e100;
	  for (int itype=0; itype<6; ++itype){// itype=Re(G00),Im(G00),Re(F00),Im(F00)
	    int row = ib + Nbath_unique*itype;
	    y_remove(ib,0) = std::max(y_remove(ib,0), x_remove(row,0));
	    y_remove(ib,1) = std::min(y_remove(ib,1), x_remove(row,1));
	  }
	  //y_remove(ib,0) = max(x_remove(Range(4*ib,4*ib+3),0));
	  //y_remove(ib,1) = min(x_remove(Range(4*ib,4*ib+3),1));
	}
	(*mpi.out)<<" broadening kwidth["<<N_min_order+iorder<<"]="<<kwidth<<" remove=";
	for (int j=0; j<G00.extent(1); ++j){
	  (*mpi.out)<<"(";
	  if (y_remove(j,0)>xt(0)) (*mpi.out)<<y_remove(j,0);
	  (*mpi.out)<<",";
	  if (y_remove(j,1)<xt(xt.extent(0)-1)) (*mpi.out)<<y_remove(j,1);
	  (*mpi.out)<<") ";
	}
	(*mpi.out)<<endl;
	//(*mpi.out)<<" broadening kwidth["<<iorder<<"]="<<kwidth<<endl;
	zg00 = BroadenFunctions(xt,zg00,kwidth,*mpi.out);
	//
	for (int j=0; j<Nbath_unique; ++j){
	  int r0 = j + 0*Nbath_unique; // Re(G00)[j]
	  int r1 = j + 1*Nbath_unique; // Im(G00)[j]
	  int r2 = j + 2*Nbath_unique; // Re(F00)[j]
	  int r3 = j + 3*Nbath_unique; // Im(F00)[j]
	  int r4 = j + 4*Nbath_unique; // Re(F01)[j]
	  int r5 = j + 5*Nbath_unique; // Im(F01)[j]
	  if (y_remove(j,1)<xt(xt.extent(0)-1)){
	    double x_start = y_remove(j,1)-(1.0+0.15*y_remove(j,1))/artificial_beta; // because of broadening we should start cutting off a bit earlier than at y_remove, and when y_remove is larger spectra is broader, hence should cutoff earlier
	    for (int i=0; i<xt.extent(0); ++i){
	      double w = fermi((xt(i) - x_start)*artificial_beta);
	      zg00(r0,i) *= w;
	      zg00(r1,i) *= w;
	      zg00(r2,i) *= w;
	      zg00(r3,i) *= w;
	      zg00(r4,i) *= w;
	      zg00(r5,i) *= w;
	    }
	  }
	  if (y_remove(j,0)>xt(0)){
	    double x_start = y_remove(j,0)+(1.0+0.15*fabs(y_remove(j,0)))/artificial_beta; // because of broadening we should start cutting off a bit earlier than at y_remove, and when y_remove is larger spectra is broader, hence should cutoff earlier
	    for (int i=0; i<xt.extent(0); ++i){
	      double w = 1.0 - fermi((xt(i) - x_start)*artificial_beta);
	      zg00(r0,i) *= w;
	      zg00(r1,i) *= w;
	      zg00(r2,i) *= w;
	      zg00(r3,i) *= w;
	      zg00(r4,i) *= w;
	      zg00(r5,i) *= w;
	    }
	  }
	}
	//
	complex<double> I(0.0,1.0);
	G00(iorder,Range::all(),Range::all()) = zg00(Range(0,               Nbath_unique-1),Range::all()) + zg00(Range(1*Nbath_unique,2*Nbath_unique-1),Range::all())*I;
	F00(iorder,Range::all(),Range::all()) = zg00(Range(2*Nbath_unique,3*Nbath_unique-1),Range::all()) + zg00(Range(3*Nbath_unique,4*Nbath_unique-1),Range::all())*I;
	F01(iorder,Range::all(),Range::all()) = zg00(Range(4*Nbath_unique,5*Nbath_unique-1),Range::all()) + zg00(Range(5*Nbath_unique,6*Nbath_unique-1),Range::all())*I;
      }
    }
  }

  void PrintStval_OCA(Array<complex<double>,3>& Stval, Array<complex<double>,2>& Stval_matrix_mc, int iorder2){};
  void PrintG00_OCA(Array<complex<double>,2>& dG00_OCA, Array<complex<double>,2>& dF00_OCA, Array<complex<double>,2>& dF01_OCA, Array<complex<double>,3>& dG00, Array<complex<double>,3>& dF00, Array<complex<double>,3>& dF01, int iorder2){};
  
  void InterpolateOnRightMesh(Array<double,2>& _Sm_, Array<complex<double>,2>& G00, Array<complex<double>,2>& F00, Array<complex<double>,2>& F01, 
			      const Array<complex<double>,3>& Stval, 
			      const Array<complex<double>,3>& dG00, const Array<complex<double>,3>& dF00, const Array<complex<double>,3>& dF01, 
			      const Cix& cix, bool cmpPhysG,
			      const json& pars, const MPIcls& mpi)
  {
    Array<double,1> xt(wNbin_lmax); 
    for (int i=0; i<wNbin_lmax; ++i) xt(i) = cutoffw*( -1.0 + (2*i+1.0)/wNbin_lmax );
    //int Npseudo = Stval.extent(1);
    //int Nbath = dG00.extent(1);
    //int Npseudo_unique = deg_pp.extent(0);
    //int Nbath_unique =  deg_bath.extent(0);
    int Norders = Stval.extent(0);
    //(*mpi.out)<<"Before InterpolateOnRightMesh sum(Stval)="<<sum(Stval)<<endl;
    Array<complex<double>,2> Sm(cix.Npseudo_unique,xn.size());
    //Sm.resize(cix.Npseudo_unique,xn.size());
    Sm=0;
    for (int ii=0; ii<cix.Npseudo_unique; ++ii){
      Array<complex<double>,1> Sml(xt.size()); Sml=0; // Sml is a sum of Stval over all orders being sampled.
      // WARNING: here we remove iorder=0 because it is added analytically!
      for (int iorder=1; iorder<Norders; ++iorder) Sml(Range::all()) += Stval(iorder,ii,Range::all());
      int posx=0;
      for (int i=0; i<xn.size(); ++i){
	if (fabs(xn(i))<=cutoffw){
	  intpar px  = InterpLeft(xn(i),posx,xt);
	  Sm(ii,i) += linInterp(Sml,px);
	}
      }
    }
    // We still do Kramers-Kronig, because real axis turns out to be less stable than imaginary axis.
    // Maybe should consider to use both the real and imaginary part in future.
    _Sm_.resize(Sm.shape());
    _Sm_ = imag(Sm);
    
    if (mpi.rank==mpi.master and pars["fastFilesystem"].get<bool>()){
      for (int iorder=0; iorder<Norders; ++iorder){
	ofstream hout("Stval."s+to_string(iorder+N_min_order));
	ofstream rout("Strel."s+to_string(iorder+N_min_order));
	for (int i=0; i<xt.size(); ++i){
	  hout<<xt(i)<<" ";
	  rout<<xt(i)<<" ";
	  for (int ii=0; ii<cix.Npseudo_unique; ++ii) hout<<Stval(iorder,ii,i).imag()<<" ";
	  for (int ii=0; ii<cix.Npseudo_unique; ++ii) rout<<Stval(iorder,ii,i).real()<<" ";
	  hout<<endl;
	  rout<<endl;
	}
      }
    }
    if (cmpPhysG){
      G00.resize(cix.Nbath_unique,dlom.size()); G00=0;
      F00.resize(cix.Nbath_unique,dlom.size()); F00=0;
      F01.resize(cix.Nbath_unique,dlom.size()); F01=0;
      for (int ib=0; ib<cix.Nbath_unique; ++ib){
	Array<complex<double>,1> G00l(xt.size()), F00l(xt.size()), F01l(xt.size()); G00l=0; F00l=0; F01l = 0; 
	// WARNING iorder==0 is computed analytically, hence removed here.
	for (int iorder=1; iorder<Norders; ++iorder){
	  G00l(Range::all()) += dG00(iorder,ib,Range::all());
	  F00l(Range::all()) += dF00(iorder,ib,Range::all());
	  F01l(Range::all()) += dF01(iorder,ib,Range::all());
	}
	int posx=0;
	for (int i=0; i<dlom.size(); ++i){
	  if (fabs(dlom(i))<=cutoffw){
	    intpar px  = InterpLeft(dlom(i),posx,xt);
	    G00(ib,i) += linInterp(G00l,px);
	    F00(ib,i) += linInterp(F00l,px);
	    F01(ib,i) += linInterp(F01l,px); 
	  }
	}
      }
      //for (int ib=0; ib<cix.Nbath_unique; ++ib){
      //G00(ib,Range::all()) *= 1./cix.deg_bth(ib);
      //F00(ib,Range::all()) *= 1./cix.deg_bth(ib);
      //}
      if (mpi.rank==mpi.master and pars["fastFilesystem"].get<bool>()){
	for (int iorder=0; iorder<Norders; ++iorder){
	  ofstream gout("g00."s+to_string(iorder+N_min_order));
	  for (int i=0; i<xt.size(); ++i){
	    gout<<xt(i)<<" ";
	    for (int ip=0; ip<cix.Nbath_unique; ++ip) gout<<dG00(iorder,ip,i).real()<<" "<<dG00(iorder,ip,i).imag()<<" ";
	    for (int ip=0; ip<cix.Nbath_unique; ++ip) gout<<dF00(iorder,ip,i).real()<<" "<<dF00(iorder,ip,i).imag()<<" ";
	    for (int ip=0; ip<cix.Nbath_unique; ++ip) gout<<dF01(iorder,ip,i).real()<<" "<<dF01(iorder,ip,i).imag()<<" ";
	    gout<<endl;
	  }
	}
	ofstream tout("g00.t");
	for (int i=0; i<dlom.size(); ++i){
	  tout<<dlom(i)<<" ";
	  for (int ip=0; ip<cix.Nbath_unique; ++ip) tout<<G00(ip,i).real()<<" "<<G00(ip,i).imag()<<" ";
	  for (int ip=0; ip<cix.Nbath_unique; ++ip) tout<<F00(ip,i).real()<<" "<<F00(ip,i).imag()<<" ";
	  for (int ip=0; ip<cix.Nbath_unique; ++ip) tout<<F01(ip,i).real()<<" "<<F01(ip,i).imag()<<" "; 
	  tout<<endl;
	}
      }
    }
  }
  void InterpolateOnRightMesh(Array<double,2>& _Sm_, Array<complex<double>,2>& G00, Array<complex<double>,2>& F00, Array<complex<double>,2>& F01,
			      Array<complex<double>,3>& Stval,
			      Array<complex<double>,3>& dG00, Array<complex<double>,3>& dF00, Array<complex<double>,3>& dF01,
			      const Array<complex<double>,2>& Stval_cOCA, const Array<complex<double>,2>& G00_cOCA, const Array<complex<double>,2>& F00_cOCA, const Array<complex<double>,2>& F01_cOCA,
			      const Cix& cix, bool cmpPhysG,
			      const json& pars, const MPIcls& mpi)
  {
    InterpolateOnRightMesh(_Sm_, G00, F00, F01, static_cast<const Array<complex<double>,3>&>(Stval),
			   static_cast<const Array<complex<double>,3>&>(dG00), static_cast<const Array<complex<double>,3>&>(dF00), static_cast<const Array<complex<double>,3>&>(dF01),
			   cix, cmpPhysG, pars, mpi);
  }
  bool Read_sampled_data(Array<double,1>& Phi_accumulate, Array<double,1>& Phi_accumulatef,
			 Array<double,1>& N_measure, Array<double,1>& weight, Array<complex<double>,3>& Stval,
			 Array<complex<double>,3>& dG00, Array<complex<double>,3>& dF00, Array<complex<double>,3>& dF01,
			 const string& fsampled_data, MPIcls& mpi)
  {
    int Norders = Stval.extent(0);
    int Npseudo = Stval.extent(1);
    int _wNbin_lmax = Stval.extent(2);
    int Nbaths = dG00.extent(1);
    int _wNw_lmax = dG00.extent(2);
    int _wNbin_lmax_=0;
    bool previous_cmpPhysG=false;
    bool read_succeeded=true;
    if (mpi.rank==mpi.master){
      std::ifstream in(fsampled_data);
      if (!in) {
        (*mpi.out) << "Could not open file " << fsampled_data << std::endl;
	read_succeeded=false;
      }
      if (read_succeeded){
	int iorder, atomic_state, ibath, l;
	double tmp;
	// Read Phi_accumulate, N_measure, weight
	int iord;
	for (iorder = 0; iorder < Norders; ++iorder) {
	  in >> iord >> Phi_accumulate(iorder) >> Phi_accumulatef(iorder) >> N_measure(iorder) >> weight(iorder);
	  if (iord!=iorder or (not in)){
	    read_succeeded=false;
	    break;
	  }
	}
	if (not in or in.eof()){
	  (*mpi.out)<<"Note1 "<<fsampled_data<<" was not properly read."<<endl;
	  read_succeeded=false;
	}
      }
      if (read_succeeded){
	// Read wNbin_lmax
	in >> _wNbin_lmax_;
	if ( _wNbin_lmax_ != _wNbin_lmax){
	  read_succeeded=false;
	  (*mpi.out)<<"Note2 "<<fsampled_data<<" was not properly read."<<endl;
	}
      }
      if (read_succeeded){
	// Read Stval: shape (norder, n_atomic_states, n_tau)
	for (int iorder = 0; iorder < Norders; ++iorder) {
	  for (int atomic_state = 0; atomic_state < Npseudo; ++atomic_state) {
	    int idx1, idx2;
	    in >> idx1 >> idx2;
	    if (idx1!=iorder or idx2!=atomic_state or (not in)){
	      read_succeeded=false;
	      (*mpi.out)<<"Note3 "<<fsampled_data<<" was not properly read."<<endl;
	      break;
	    }
	    double _re_, _im_;
	    for (int i = 0; i < _wNbin_lmax; ++i){
	      in >> _re_ >> _im_;
	      Stval(iorder, atomic_state, i) = complex<double>(_re_,_im_);
	    }
	    if (not in or in.eof()){
	      read_succeeded=false;
	      (*mpi.out)<<"Note4 "<<fsampled_data<<" was not properly read."<<endl;
	      break;
	    }
	  }
	  if (not read_succeeded) break;
	}
      }
      if (read_succeeded){
	int _previous_cmpPhysG_;
	in >> _previous_cmpPhysG_;
	previous_cmpPhysG = static_cast<bool>(_previous_cmpPhysG_);
	if (previous_cmpPhysG){
	  (*mpi.out)<<"previous_cmpPhysG exists.... reading it"<<endl;
	  // Read dG00 and dF00, dF01: both (norder, nbath, lmax)
	  for (int iorder = 0; iorder < Norders; ++iorder) {
	    for (int ibath = 0; ibath < Nbaths; ++ibath) {
	      int idx1, idx2;
	      in >> idx1 >> idx2;
	      if (idx1!=iorder or idx2!=ibath or (not in) or in.eof()){
		read_succeeded=false;
		(*mpi.out)<<"Note5 "<<fsampled_data<<" was not properly read."<<endl;
		break;
	      }
	      double _re_, _im_;
	      for (int l=0; l<_wNw_lmax; ++l){
		// dG00
		in >> _re_ >> _im_;
		dG00(iorder,ibath,l) = complex<double>(_re_,_im_);
		// dF00
		in >> _re_ >> _im_;
		dF00(iorder,ibath,l) = complex<double>(_re_,_im_);
		// dF01                                        // NEW
		in >> _re_ >> _im_;
		dF01(iorder,ibath,l) = complex<double>(_re_,_im_);
	      }
	      if (!in || in.eof()){
		read_succeeded=false;
		(*mpi.out)<<"Note6 "<<fsampled_data<<" was not properly read."<<endl;
		break;
	      }
	    }
	  }
	}else{
	  dG00=0;
	  dF00=0;
	  dF01=0;
	  Phi_accumulatef=0; // Phi_accumulatef should only be from calculations with cmpPhysG turned on.
	}
      }
      in.close();
    }
    mpi.bcast_two_boolean(read_succeeded,previous_cmpPhysG);
    (*mpi.out)<<"read_succeeded="<<read_succeeded<<" cmpPhysG="<<previous_cmpPhysG<<endl;
    if (read_succeeded){
      mpi.Bcast_simple_array(Phi_accumulate);
      mpi.Bcast_simple_array(Phi_accumulatef);
      mpi.Bcast_simple_array(N_measure);
      mpi.Bcast_simple_array(weight);
      mpi.Bcast_simple_array(Stval);
      if (previous_cmpPhysG){
	mpi.Bcast_simple_array(dG00);
	mpi.Bcast_simple_array(dF00);
	mpi.Bcast_simple_array(dF01);
      }else{
	dG00=0;
	dF00=0;
	dF01=0;
      }
    }else{
      Phi_accumulate=0;
      Phi_accumulatef=0;
      N_measure=0;
      Stval=0;
      dG00=0;
      dF00=0;
      dF01=0;
      weight=1.0/weight.size();
    }
    return read_succeeded;
  }
  bool Read_sampled_data(Array<double,1>& Phi_accumulate, Array<double,1>& Phi_accumulatef,
			 Array<double,1>& N_measure, Array<double,1>& weight, Array<complex<double>,3>& Stval,
			 Array<complex<double>,3>& dG00, Array<complex<double>,3>& dF00, Array<complex<double>,3>& dF01,
			 double& Phi_cOCA, double& Phi_cOCAf, Array<complex<double>,2>& Stval_cOCA, Array<complex<double>,2>& G00_cOCA, Array<complex<double>,2>& F00_cOCA, Array<complex<double>,2>& F01_cOCA,
			 const string& fsampled_data, MPIcls& mpi)
  {
    bool read_succeeded = Read_sampled_data(Phi_accumulate, Phi_accumulatef, N_measure, weight, Stval, dG00, dF00, dF01, fsampled_data, mpi);
    Phi_cOCA = 0;
    Phi_cOCAf = 0;
    if (Stval_cOCA.size()>0) Stval_cOCA = complex<double>(0,0);
    if (G00_cOCA.size()>0) G00_cOCA = complex<double>(0,0);
    if (F00_cOCA.size()>0) F00_cOCA = complex<double>(0,0);
    if (F01_cOCA.size()>0) F01_cOCA = complex<double>(0,0);
    return read_succeeded;
  }
  void CheckQuantities(const complex<double>& cPhi, int itt, const Array<int,2>& diag, const Array<int,1>& dirspin, const Array<double,1>& freq, const Cix& cix, ostream& cout)
  {
    int Norder = diag.extent(0);
    DetermineLoops(wtb_t,/**/diag, dirspin);

    if ( sum(abs(wtb-wtb_t)) > 0){
      cout<<"ERROR wtb="<<wtb<<" while wtb_t="<<wtb_t<<endl;
    }
    //
    DetermineFermiSigns_fsgn_only(xfsgn_t,/**/wtb_t);
    if ( sum(abs(xfsgn_t-xfsgn)) > 0){
      cout<<"ERROR xfsgn="<< xfsgn << " xfsgn_t="<<xfsgn_t<<endl;
    }
    
    tPhi_dd = 1;
    for (int iy=0; iy<Norder; ++iy){ // 0:Norder for all possible configurations
      int dire = dirspin(iy) % 2;    // direction of hybridization propagator
      int iorb = dirspin(iy) / 2;    // its spin 
      int s = cix.bath_sigind(iorb); // actual entry for hybridization
      double y = freq(1+iy);         // frequency on hybridization
      intpar jpy = Interp(y, dlom);  // 
      fry_t(iy).set(fep(jpy), lfep(jpy), beta*y ); // fermi function for f(y) stored as logarithm log(f(y))
      double t_Dlty;
      if (wtb_t(iy,iy+1)>0){
	t_Dlty = linInterp(Deltai, s, jpy); // actual interpolation and saving the interpolated hybridization
      }else{
	intpar jmpy = Interp(-y, dlom);
	t_Dlty = linInterp(Deltai, s, jmpy); // actual interpolation and saving the interpolated hybridization
      }
      if (abs(t_Dlty-Dlty(iy))>1e-8){
	cout<<"ERROR Dlty["<<iy<<"]="<<Dlty(iy)<<" while Dlty_t="<<t_Dlty<<endl;
      }
      tPhi_dd *= t_Dlty;
    }
    if (abs(cPhi_dd-tPhi_dd)>1e-8){
      cout<<"ERROR cPhi_dd="<<cPhi_dd<<" tPhi_dd="<<tPhi_dd<<endl;
    }
    
    auto& tallowed = tdiag2matrix.Cmp_diagonal3(diag, dirspin, cix);      // map<int,Array<int,2>>

    for (const auto& [ib, alw] : tallowed){
      auto& alwt = cdiag2matrix.diagonal[ib];
      if (sum(abs(alw-alwt))!=0){
	cout<<"ERROR:  alw["<<ib<<"]="<<alw<<" while alwt="<<alwt<<endl;
      }
    }
    
    for (int j=0; j<2*Norder; ++j){ // going over all pseudo-particles
      // the G[j] argument can be omega+y_0-y_2+y_3+.... Such combination would correspond to array wtb[Norder+j,:]=[1,1,0,-1,1,..]
      Range all_freq(0,Norder);
      double x = sum(freq(all_freq)*wtb(Norder+j,all_freq)); // This is the frequency argument for this PP, namely, sum(wtb[Norder+j,:]*freq[:])
      if (abs(pp_freq(j)-x)>1e-8){
	cout<<"ERROR pp_freq="<<pp_freq(j)<<" while x="<<x<<endl;
      }
      intpar px = Interp(x, xn);                  // interpolation of x in xn
      if (pxs(j).i != px.i or fabs(pxs(j).p-px.p)>1e-8) cout<<"ERROR pxs["<<j<<"]="<<pxs(j)<<" while px="<<px<<endl;

      nFexp xfrms( femx(px), lfemx(px), -beta*x); // fermi function f(-x) with argument x of pseudo propagator j.
      if (fabs(frms_(j).f-xfrms.f)>1e-8 or fabs(frms_(j).lf-xfrms.lf)>1e-8 or  fabs(frms_(j).x-xfrms.x)>1e-8 ){
	cout<<"ERROR frms_["<<j<<"]="<<frms_(j)<<" while xfrms="<<xfrms<<endl;
      }
    }
    
    
    for (int j=0; j<2*Norder; ++j){ // going over the ring
      // Note: lfpt(ipseudo) = log(  \product_{iy} f(\pm_{ipseudo,iy} y)  )
      double lfp=0; // This is the product of fermi functions needed for Self-energy Sigma_j of pseudoparticle j.
      for (int iy=0; iy<Norder; ++iy){
	if (xfsgn(j) & (1<<iy)){ // xfsgn is represented as binary array in an integer. 
	  lfp += fry_(iy).plf(); // log(f(y))
	}else{
	  lfp += fry_(iy).mlf(); // log(f(-y))
	}
      }
      //lfpt_t(j)=lfp; // This is the product of fermi functions needed for Self-energy Sigma_j of pseudoparticle j.
      if (abs(lfpt(j)-lfp)>1e-8){
	cout<<"ERROR lfpt["<<j<<"]="<<lfpt(j)<<" while lfp="<<lfp<<endl;
      }
      //fmt::print("fpt[{:d},{:d}]={:g},{:f}\n", ic,j, fpt(ic,j).f, fpt(ic,j).lf);
      //if (DEBUG) cout<<" pp_fr["<<j<<"]="<<pp_freq_t(j)<<" px="<<px<<" Gmt="<<Gmt_t(Range(0,PPs-1),j)<<" mGppt="<<mGppt_t(Range(0,PPs-1))<<endl;
    }

    //cout<<" pp_freq="<<pp_freq_t(Range(0,2*Norder-1))<<" pxs="<<pxs_t(Range(0,2*Norder-1))<<" lfpt="<<lfpt_t(Range(0,2*Norder-1))<<endl;
    
    //tPPs=0;
    //for (const auto& [ib, alw] : tallowed)  tPPs += alw.extent(1);
    //Gmt_t.resize(tPPs,2*Norder); Gmt_t = 0;
    int idx=0;
    for (const auto& [ib, alw] : tallowed){ // over all blocks that contribute on diagonal
      for (int j=0; j<alw.extent(0); ++j){  // in principle more than one term can contribute per block
	complex<double> _mGpp_ = 1.0;
	for (int i=0; i<2*Norder; ++i){     // over the ring of PPs
	  //int _ab_ = alw(j,i);
	  //int _a_ = abs(_ab_);
	  //int ab = cix.inBlock(_a_);      // _ab_ is given in terms of atomic state index, but ab is index within the block. cix.inBlock translates this.
	  //int _ib_ = cix.iBlck(_a_)-1;
	  //int ii = cix.Sgm_symmetry[_ib_](ab,ab)-1; // Green's function index for propagator G(a0,b0)==G[|_ii0_|-1]
	  int _ii_ = alw(j,i);
	  int ii = abs(_ii_)-1;
	  if (ii>=0){
	    //double sg = sign(_ab_);
	    double sg = sign(_ii_);
	    //complex<double> Gmt_t = linInterp(Gm, ii, pxs(i)) * sg;  // now actual interpolation of Gm and setting value of Gm*F. This is not just Gm, but multiplied by F.
	    //if (abs(Gmt_(idx,i)-Gmt_t)>1e-8) cout<<"ERROR Gmt["<<idx<<","<<i<<"]="<<Gmt_(idx,i)<<" while Gmt_t="<<Gmt_t<<endl;
	    Gexp Gmt_t; Gmt_t.linInterp(Gms, ii, pxs(i), sg);
	    if (abs(Gmt_(idx,i).c-Gmt_t.c)>1e-8 or fabs(Gmt_(idx,i).p-Gmt_t.p)>1e-8 or fabs(Gmt_(idx,i).a-Gmt_t.a)>1e-8 )
	      cout<<"ERROR Gmt["<<idx<<","<<i<<"]="<<Gmt_(idx,i)<<" while Gmt_t="<<Gmt_t<<endl;
	    if (New_form_Phi==2){
	      double phiR_t = linInterp(PhiR, ii, pxs(i));
	      if (fabs(phiR_t-phiR(idx,i))>1e-8){
		cout<<"ERROR phiR["<<idx<<","<<i<<"]="<< phiR(idx,i)<<" phiR_t="<<phiR_t<<endl;
	      }
	    }
	    if (New_form_Phi>=15){
	      complex<double> _Gmo_ = linInterp(Gm,ii, pxs(i))*sg;
	      if (abs(Gmo_(idx,i)-_Gmo_)>1e-6){
		cout<<itt<<") ERROR Gmo_["<<idx<<","<<i<<"]="<<Gmo_(idx,i)<<" _Gmo_="<<_Gmo_<<endl;
	      }
	    }
	  }
	  if (New_form_Phi>=15) _mGpp_ *= Gmo_(idx,i);
	}
	if (New_form_Phi>=15 and abs(mGppt(idx)-_mGpp_)>1e-6){
	  cout<<itt<<") ERROR mGppt["<<idx<<"]="<<mGppt(idx)<<" while _mGpp_="<<_mGpp_<<endl;
	}
	idx += 1;
      }
    }
    if (New_form_Phi>=10){
      sum_fp_t=0;
      for (int iy=0; iy<Norder; ++iy){ // Now compute the product of fermi functions entering the physical green's function
	double max_fp=-10*beta;
	for (int ip=0; ip<2*Norder; ++ip){
	  // sum_fp[ic] = \sum_{iy_ext}\sum_{ip}\product_{iy\ne iy_ext} f(\pm_{ip,iy} y)
	  // sum_fp[ic] = \sum_{ip}\sum_{iy_ext} 1/f(\pm_{ip,iy_ext} y_ext) \product_{iy} f(\pm_{ip,iy} y)
	  double lfp = lfpt(ip); // lfp(ipseudo) = log(  \product_{iy} f(\pm_{ipseudo,iy} y)  )
	  if (xfsgn(ip) & (1<<iy)){
	    lfp -= fry_(iy).plf(); // now dividing : + log( 1/f(\pm_{ip,iy_ext} y_ext) )
	  }else{
	    lfp -= fry_(iy).mlf(); // now dividing : + log( 1/f(\pm_{ip,iy_ext} y_ext) )
	  }
	  if (New_form_Phi==10) lfp += frms_(ip).mlf(); // we are multiplying by f(+x_ip)
	  if (lfp>max_fp) max_fp=lfp; // we would need to take exponent of lfp, but we will just find the largest, and take exponent only once.
	  //double crnt= (xfsgn_t(ip) & (1<<iy)) ? fry_t(iy).plf(): fry_t(iy).mlf();
	  //cout<<".... for iy="<<iy<<" ip="<<ip<<" lfp="<<lfp<<" total="<<lfpt_t(ip)<<" current="<<crnt<<" max_fp="<<max_fp<<endl;
	}
	// instead of computing \sum_{ip} exp(lfp_{iy_ext,ip}) we just compute exp(max(lfp_{iy_ext,ip})), which is
	// a reasonable proxy, good enough for MC walk.
	sum_fp_t += exp(max_fp);
      }
      if (fabs(sum_fp-sum_fp_t)>1e-8) cout<<"ERROR sum_fp="<<sum_fp<<" while sum_fp_t="<<sum_fp_t<<endl;
    }
    
    complex<double> Phi_pp=0;
    idx=0; int _idx_=0;
    for (const auto& [ib, alw] : tallowed){ // over all blocks that contribute on diagonal
      complex<double> phi_pp=0;
      for (int j=0; j<alw.extent(0); ++j){ // in principle more than one term can contribute per block
	Gexp pd; pd.Init2Unity();
	for (int ip=0; ip<2*Norder; ++ip) pd *= Gmt_(idx,ip);
	if (New_form_Phi==1){
	  for (int ip=0; ip<2*Norder; ++ip)
	    phi_pp -= pd.get0(lfpt(ip), frms_(ip))/M_PI;
	}else if (New_form_Phi==2){
	  for (int ip=0; ip<2*Norder; ++ip){
	    Gexp Phis_over_Gm(pd);
	    Phis_over_Gm.divide_by_positive(Gmt_(idx,ip));
	    Phis_over_Gm.matmul_lferm(lfpt(ip));
	    phi_pp += Phis_over_Gm.getImt(frms_(ip).lf)*phiR(idx,ip);
	  }
	}else if (New_form_Phi==10 or New_form_Phi==11){
	  double pd_im = pd.p*exp(pd.a), pd_abs = abs(pd.c);
	  phi_pp += pd_abs * sum_fp_t * (2*Norder);
	  //phi_pp += pd_im * sum_fp_t * (2*Norder);
	  //for (int ip=0; ip<2*Norder; ++ip)
	  //  phi_pp -= pd_im*exp(frms_t(ip).x)*sum_fp_t/M_PI;
	}else if (New_form_Phi==15){
	  phi_pp += mGppt(idx) * sum_fp * (2.0*Norder);
	  //cout<<" mGpp["<<idx<<"]="<< mGppt_t(idx) <<" sum_fp_t="<<sum_fp_t<<endl;
	}
	idx += 1;
      }
      phi_pp /= (2.*Norder);
      Phi_pp += phi_pp;
      complex<double> tPhiv = phi_pp*cPhi_dd;
      if (abs(cPhiv(_idx_) - tPhiv)>1e-8){
	cout<<itt<<") ERROR cPhiv["<<_idx_<<"]="<<cPhiv(_idx_)<<" tPhiv="<<tPhiv<<endl;
      }
      _idx_ += 1;
    }
    complex<double> ret = Phi_pp*cPhi_dd;
    if ( abs(ret)==0 or abs(1.0-ret/cPhi)>1e-6 ){
      cout<<"ERROR cPhi="<<cPhi<<" ret="<<ret<<endl;
      cout<<"Phi_pp="<<Phi_pp<<" cPhi_dd="<<cPhi_dd<<" tPhi_dd="<<tPhi_dd<<endl;
    }
  }
  void PrintTimings(ostream& log)
  {
#ifdef _TIMEIT_
    //fsec t_lft, t_rgh, t_dphi, t_stval;
    fmt::print(log,"[t_lft={:.4f}, t_rgh={:.4f}, t_dphi={:.4f}, t_stval={:.4f}] [t_plan={:.4f}, t_pre_plan={:.4f}, t_post_plan={:.4f}] [t_H={:.4f}, t_interp={:.4f}, t_GF2={:.4f}]\n",
	       t_lft, t_rgh , t_dphi , t_stval , t_plan , t_pre_plan , t_post_plan , t_H , t_interp , t_GF2 );

    fmt::print(log," t_GF2=[t_GF1={:.4f} + t_Hf={:.4f} + t_Hg={:.4f}] t_GF1=[t_GGin={:.4f} + t_GGln={:.4f} + t_GGi={:.4f} + t_GGl={:.4f}]\n",
	       t_GF1 ,t_Hf ,t_Hg , t_GGin ,t_GGln ,t_GGi ,t_GGl );
    /////// [t_H,t_interp,t_GF2]
    /////// t_GF2=t_GF1+t_Hf+t_Hg]
    /////// t_GF1=t_GGin+t_GGln+t_GGi+t_GGl
#endif    
  }
};
#endif // Real_axis
