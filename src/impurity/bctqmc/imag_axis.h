#ifndef Imag_axis
#define Imag_axis

#include "nca_oca_imag_axis.h"
#include "imag_axis_spectra.h"
#include "topologymove.h"
#include "tca_paths.h"
#include "aaa.h"
#include "dlr.h"
#include "imag_axis_dlr.h"
#include "oca_dlr.h"
#include "tca_dlr.h"

#include <chrono>
#include <memory>

#define _TIMEIT_

double Product(const Array<double,1>& A)
{
  double res=1.0;
  for (int i=0; i<A.size(); ++i) res *= A(i);
  return res;
}

class ImagAxis0 : public AnyAxis{
public:
  Array<double,1> tau, tau_dh;
  Array<double,2> Delta;
  //vector<Spline1D<double>> Deltas;
  Spline1DMany Deltas;
  Spline1DMany Gms, gms;  
  Array<double,1> iom;
  Array<int,1> iom_n;       // fermionic Matsubara indices corresponding to iom
  Array<double,1> w_iom;    // weights which reproduce the full Matsubara sum from the compressed mesh
  Array<complex<double>,2> iDeltaw;
  SVDFunc svdf;
  const Cix& cix;
  vector<int> G_importance;  int OCA_cut_sum;
  // Compact-diagram state.  These adapters are intentionally optional: the
  // legacy spline/uniform-mesh OCA remains the production result until the
  // comparisons below are accepted for all target problems.
  bool aaa_adapter_enabled=false, dlr_adapter_enabled=false;

  // AAA objects have distinct roles. aaa_result is the barycentric fit used
  // to discover poles; aaa_hybridization is the real-pole partial-fraction
  // expansion actually used by the diagram code.
  edmft::aaa::Result aaa_result;
  edmft::aaa::PoleExpansion aaa_hybridization;
  // Optional lower-rank hybridization used only by approximate TCA.  OCA and
  // the exact TCA audit always use aaa_hybridization above.
  edmft::aaa::PoleExpansion tca_active_hybridization;

  // The DLR basis depends only on beta, spectral width, and tolerance.  The
  // workspace combines it with the current AAA hybridization and current
  // pseudo propagators.  Plans contain only topology/sparse-index data and
  // may therefore be cached across self-consistency iterations.
  std::unique_ptr<edmft::dlr::FermionicBasis> dlr_basis;
  std::unique_ptr<edmft::imag_dlr::Workspace> diagram_dlr_workspace;
  std::unique_ptr<edmft::imag_dlr::Workspace> tca_active_dlr_workspace;
  std::unique_ptr<edmft::oca_dlr::P3Plan> oca_dlr_p3_plan;
  std::unique_ptr<edmft::oca_dlr::SelfEnergyPlan> oca_dlr_self_plan;
  std::unique_ptr<edmft::oca_dlr::PhysicalPlan> oca_dlr_physical_plan;
  // True when the cached plans contain only this MPI rank's reduced bath-pair
  // shard.  The configured strings are retained so rereading PARAMS with a
  // different distribution policy invalidates the old topology immediately.
  bool oca_dlr_plans_distributed=false;
  int oca_dlr_topology_shard_count=1;
  std::string oca_dlr_cached_topology_setting;
  double oca_dlr_cached_replication_budget=-1.0;
  // Legacy integration needs the complete (and potentially very large)
  // enumeration.  Retain it only for fixed-order runs.  Genuine higher-order
  // Monte Carlo needs the memory more than it needs to avoid re-enumeration.
  std::unique_ptr<OCAEnumeration> oca_legacy_enumeration;
  Eigen::MatrixXd dlr_G_coefficients, dlr_g_coefficients;
  // The active pseudo mask is chosen once from the OCA-converged propagators.
  // The frozen correction is stored at DLR nodes, so its memory cost is only
  // O(r*Npseudo), independent of the dense tau mesh.
  std::vector<unsigned char> tca_active_pseudo_mask;
  Eigen::MatrixXd tca_frozen_sigma_tail;
  std::array<Eigen::MatrixXd,2> tca_frozen_sigma_tail_by_topology;
  std::array<double,2> tca_frozen_phi_tail_by_topology{{0.0,0.0}};
  bool tca_frozen_tail_ready=false;
  double tca_cached_keep_weight=1.0;
  double aaa_relative_rms_error=0.0, aaa_max_error=0.0;
  double dlr_G_relative_error=0.0, dlr_g_relative_error=0.0;
  //using AnyAxis::AnyAxis; // <-- inherit constructors: All constructors of AnyAxis are also constructors of ImagAxis0
  //
  ImagAxis0(double _beta_, const Cix& _cix_, json& pars, MPIcls& mpi) :
    AnyAxis(_beta_), cix(_cix_), G_importance(cix.Npseudo_unique,0), OCA_cut_sum(10)
  {
    //cout<<"In ImagAxis0 constructor beta="<<beta<<endl; cout.flush();
    pars["svd_x0"] = min(1./beta, pars["svd_x0"].get<double>());
    //cout<<"svd_x0="<< pars["svd_x0"].get<double>()<<endl; cout.flush();
    svdf.cmp("fermi", beta, pars["svd_lmax"], pars["Ntau"], pars["svd_x0"], pars["svd_L"], pars["svd_Nw"], *mpi.out, mpi.rank==mpi.master);
    // Print SVD functions
    //cout<<"cmp done"<<endl; cout.flush();
    if (mpi.rank==mpi.master and pars["fastFilesystem"]) svdf.Print("uls.dat");
    svdf.SetUpFastInterp();
    // Now we check how well is the orthogonality obeyed after reorthogonalization
    if (mpi.rank==mpi.master) (*mpi.out)<<"Final overlap is smaller than "<<svdf.CheckOverlap()<<endl;

    // OCA_solver is the user-facing production switch.  Selecting "dlr"
    // implicitly enables both compact representations, so a PARAMS file needs
    // only one setting.  The two *_adapter flags remain useful for building
    // and diagnosing the representations while legacy OCA is still selected.
    const std::string oca_solver=pars["OCA_solver"].get<std::string>();
    if (oca_solver!="legacy" && oca_solver!="dlr")
      throw std::invalid_argument(
        "OCA_solver must be either \"legacy\" or \"dlr\", received \""+
        oca_solver+"\"");
    const std::string tca_solver=pars["TCA_solver"].get<std::string>();
    if (tca_solver!="dlr" && tca_solver!="mc")
      throw std::invalid_argument(
        "TCA_solver must be either \"dlr\" or \"mc\", received \""+
        tca_solver+"\"");
    // Analytic TCA uses exactly the same AAA and DLR representations as OCA,
    // even when a diagnostic run deliberately keeps legacy OCA selected.
    const bool analytic_tca_requested=tca_solver=="dlr" &&
      pars["N_min_order"].get<int>()==3 && pars["Norder"].get<int>()==3;
    const bool production_dlr=oca_solver=="dlr" || analytic_tca_requested;
    aaa_adapter_enabled = pars["aaa_adapter"].get<bool>() || production_dlr;
    dlr_adapter_enabled = pars["dlr_adapter"].get<bool>() || production_dlr;
    if (dlr_adapter_enabled){
      // A nonpositive dlr_omega_max means "follow the solver cutoff".  The
      // dimensionless difficulty printed below is Lambda=beta*omega_max;
      // panel_order samples the temporary construction grid, while rank is
      // the much smaller number of selected frequencies/time nodes.
      double omega_max=pars["dlr_omega_max"].get<double>();
      if (omega_max<=0.0) omega_max=pars["cutoffw"].get<double>();
      edmft::dlr::Options opts;
      opts.beta=beta;
      opts.omega_max=omega_max;
      opts.eps=pars["dlr_tolerance"].get<double>();
      opts.panel_order=pars["dlr_panel_order"].get<int>();
      dlr_basis=std::make_unique<edmft::dlr::FermionicBasis>(opts);
      diagram_dlr_workspace=std::make_unique<edmft::imag_dlr::Workspace>();
      if (mpi.rank==mpi.master)
        (*mpi.out) << "DLR adapter enabled: beta=" << beta
                   << " omega_max=" << omega_max
                   << " Lambda=" << dlr_basis->lambda()
                   << " tolerance=" << dlr_basis->tolerance()
                   << " rank=" << dlr_basis->rank() << std::endl;
    }
  };
  //void debug1(MPIcls& mpi){};
  int Coeff_size() const {return svdf.lmax;}
  const Array<double,1>& freq() const override { return iom; } // AnyAxis has now access to freq.
  void BuildSparseMatsubaraWeights(MPIcls& mpi)
  {
    // The imaginary-axis mesh iom contains true fermionic Matsubara frequencies
    // iom = (2*n+1)*pi/beta, but the tail can be compressed: only some integer
    // Matsubara indices n are kept explicitly.
    //
    // For quantities that decay sufficiently fast (at least as 1/iom^2) we
    // want to replace the exact full Matsubara sum
    //
    //   sum_{n on full mesh} f_splined(n)
    //
    // by a weighted sum on the compressed mesh
    //
    //   sum_{k on kept mesh} w_iom(k) * f(k)
    //
    // with weights chosen so that the two are *exactly equivalent* to the
    // cubic spline generated by Spline1D::splineIt.
    //
    // Because spline construction, spline evaluation, and summation are all
    // linear in the sampled values f(k), the whole "spline then sum" map is a
    // linear functional.  The weights are therefore obtained by feeding the
    // spline routine basis vectors e^(k)=(0,...,1,...,0), evaluating the
    // reconstructed basis function on all integer Matsubara points, and
    // summing the result.  This produces the k-th weight w_iom(k).
    const int N = iom.extent(0);
    iom_n.resize(N);
    w_iom.resize(N);
    w_iom = 0.0;

    if (N == 0) return;

    Array<double,1> nmesh_d(N);
    for (int i=0; i<N; ++i){
      double ni = 0.5*(iom(i)*beta/M_PI - 1.0);
      int n_round = static_cast<int>(std::round(ni));
      if (std::fabs(ni-n_round)>1e-6 && mpi.rank==mpi.master){
        (*mpi.out) << "WARN: iom(" << i << ")=" << iom(i)
                   << " does not correspond to an integer fermionic Matsubara index."
                   << " inferred n=" << ni << " rounded to " << n_round << std::endl;
      }
      iom_n(i) = n_round;
      nmesh_d(i) = n_round;
    }

    const int n0 = iom_n(0);
    const int n_last = iom_n(N-1);

    // Very common fast path: if the stored Matsubara mesh already contains
    // every consecutive fermionic Matsubara index from n0 to n_last, then the
    // "compressed" mesh is in fact the full mesh and no reconstruction is
    // needed.  In that case the exact full Matsubara sum is simply the plain
    // sum over the kept points, hence all weights are unity.
    bool mesh_is_dense = (N == n_last - n0 + 1);
    for (int i=1; mesh_is_dense && i<N; ++i)
      mesh_is_dense = (iom_n(i) == iom_n(i-1) + 1);
    if (mesh_is_dense){
      w_iom = 1.0;
      if (mpi.rank==mpi.master)
        (*mpi.out) << "Sparse Matsubara weights: mesh is dense, setting w_iom=1" << std::endl;
      return;
    }

    if (N == 1){
      w_iom(0) = 1.0;
      return;
    }

    if (N == 2){
      // Spline1D::splineIt uses boundary extrapolation that needs at least
      // three knots.  For the degenerate two-point case we fall back to the
      // unique linear interpolation and sum it exactly on the missing integer
      // Matsubara points.
      const double denom = nmesh_d(1) - nmesh_d(0);
      for (int n=n0; n<=n_last; ++n){
        double t = (n - nmesh_d(0))/denom;
        w_iom(0) += 1.0 - t;
        w_iom(1) += t;
      }
      return;
    }

    // The scalar version of this routine constructed one Spline1D object for
    // each basis vector e^(k).  Here we do the exact same mathematics with
    // Spline1DMany: every column k of basis_splines is one basis vector
    // e^(k), and splineIt solves all right-hand-sides at once.
    Spline1DMany basis_splines(N, N);
    for (int k=0; k<N; ++k){
      for (int i=0; i<N; ++i)
        basis_splines(i,k) = (i==k ? 1.0 : 0.0);
    }

    // This single call computes the cubic spline for all basis vectors at the
    // same time.  Column k is therefore the spline response to e^(k), i.e. the
    // k-th column of the linear map M that sends sampled values to values on
    // the full Matsubara mesh.
    basis_splines.splineIt(nmesh_d);

    for (int n=n0; n<=n_last; ++n){
      // Evaluate the whole basis at this integer Matsubara index.  Summing
      // over n is equivalent to summing over the first index of the matrix M.
      // The final result is w_iom(k) = sum_n M(n,k).
      intpar p = Interp(static_cast<double>(n), nmesh_d);
      for (int k=0; k<N; ++k)
        w_iom(k) += basis_splines(p,k);
    }

    if (mpi.rank==mpi.master){
      (*mpi.out) << "Sparse Matsubara indices iom_n=" << iom_n << std::endl;
      (*mpi.out) << "Sparse Matsubara weights w_iom=" << w_iom << std::endl;
    }
  }
  void PrepareDelta(Array<complex<double>,2>& Deltaw, Array<double,1>& om, const json& pars, MPIcls& mpi){
    if (cix.Nbath_unique_unique != Deltaw.extent(0))
      cout<<"ERROR Nbath_unique_unique="<<cix.Nbath_unique_unique<<" while Delta.extent(0)="<<Deltaw.extent(0)<<endl;
    InterpolateIfNeeded(om, Deltaw, beta, *mpi.out); // if Matsubara mesh is incompatible with current temperature, we interpolate.
    iom.reference(om);
    BuildSparseMatsubaraWeights(mpi);
    iDeltaw.reference(Deltaw);
    if (aaa_adapter_enabled) PrepareAAAHybridization(Deltaw,om,pars,mpi);
    //
    tau.resize(pars["Ntau"]);
    GiveDoubleExpMesh(tau, tau_dh, beta, pars["Ntau"]);
    //MakeEquidistantMesh(tau,pars["Ntau"],0,beta);
    // Notice that now we build larger array of size cix.Nbath_unique not cix.Nbath_unique_unique.
    Delta.resize(cix.Nbath_unique,tau.size()); 
    {
      Array<double,2> Deltat(cix.Nbath_unique_unique,tau.size());
      for (int ib=0; ib<cix.Nbath_unique_unique; ++ib){
	Array<double,1> delta(Deltat(ib,Range::all())); // taking a slice of Deltat
	FourierOfDelta(om,Deltaw(ib,Range::all()), tau, delta);
      }
      Array<int,1> sigind0(cix.Nbath_unique);
      for (int ib=0; ib<cix.Nbath; ++ib) sigind0(cix.bath_sigind(ib)) = cix.bath_sigind0(ib);
      for (int ib=0; ib<cix.Nbath_unique; ++ib){
	// Here we create larger array of size cix.Nbath_unique.
	// bath_sigind0 might have index with negative numbers, which means that we need to take the opposite sig.
	// bath_sigind is just regular
	double ib_sg = (sigind0(ib) < 0 ? -1 : 1);
	int iib   = std::abs(sigind0(ib));
	Delta(ib,Range::all()) = ib_sg*Deltat(iib,Range::all());
	//cout<<"Taking Delta["<<ib<<"] = Deltat["<<iib<<"]*"<<ib_sg<<endl; 
      }
      if (pars["fastFilesystem"] and mpi.rank==mpi.master){
	ofstream out("Deltat.inp");
	for (int it=0; it<tau.extent(0); ++it){
	  out<<tau(it)<<" ";
	  for (int ib=0; ib<cix.Nbath_unique; ++ib) out<<Delta(ib,it)<<" ";
	  out<<endl;
	}
      }
    }
    /*
    Deltas.resize(cix.Nbath_unique);
    for (int ib=0; ib<cix.Nbath_unique; ++ib){
      Deltas[ib].resize(tau.size());
      for (int it=0; it<tau.size(); ++it) Deltas[ib][it] = Delta(ib,it);
      Deltas[ib].splineIt(tau);
    }
    */
    Deltas.resize(TinyVector<int,2>(tau.size(),cix.Nbath_unique));
    for (int ib=0; ib<cix.Nbath_unique; ++ib)
      for (int it=0; it<tau.size(); ++it)
	Deltas(it,ib) = Delta(ib,it);
    Deltas.splineIt(tau);
    /*
    ofstream dout2("SD_debug1.dat"); dout2.precision(14);
    for (int it=0; it<tau.size(); ++it){
      dout2 << tau(it) << " ";
      for (int b=0; b<cix.Nbath_unique_unique; ++b) dout2 << Deltas(it,b) << " ";
      dout2 << endl;
    }
    */
  }
  // The CIX file stores independent hybridization channels, not necessarily
  // a dense orbital matrix. Put those channels on the diagonal of an
  // auxiliary matrix so matrix AAA chooses one common support/pole set
  // without changing or guessing the physical CIX mapping.
  void PrepareAAAHybridization(const Array<complex<double>,2>& Deltaw,
                               const Array<double,1>& om,
                               const json& pars, MPIcls& mpi)
  {
    using edmft::aaa::Complex;
    using edmft::aaa::MatrixXcd;
    // The solver stores independent scalar CIX channels.  Embedding them as
    // diagonal matrices makes matrix AAA choose one common pole set, which is
    // essential for batching channels later.  No off-diagonal physics is
    // invented by this auxiliary matrix representation.
    std::vector<Complex> z(om.extent(0));
    std::vector<MatrixXcd> values;
    values.reserve(om.extent(0));
    for (int n=0; n<om.extent(0); ++n){
      z[n]=Complex(0.0,om(n));
      MatrixXcd fn=MatrixXcd::Zero(Deltaw.extent(0),Deltaw.extent(0));
      for (int b=0; b<Deltaw.extent(0); ++b) fn(b,b)=Deltaw(b,n);
      values.push_back(std::move(fn));
    }

    edmft::aaa::Options opts;
    opts.tolerance=pars["aaa_tolerance"].get<double>();
    opts.max_support_pairs=pars["aaa_max_support_pairs"].get<int>();
    opts.relative_tolerance=true;
    const std::string requested_mode=pars["aaa_fit_mode"].get<std::string>();
    if (requested_mode!="auto" && requested_mode!="joint" &&
        requested_mode!="balanced" && requested_mode!="channelwise")
      throw std::invalid_argument(
        "aaa_fit_mode must be \"auto\", \"joint\", \"balanced\", or "
        "\"channelwise\"");

    // The expensive SVD/eigensystem/L-BFGS work is deterministic and has no
    // useful MPI decomposition.  Performing it independently on every rank
    // merely multiplies startup time and may introduce last-bit differences
    // between the pole tables.  Rank zero constructs one expansion below;
    // the compact poles and diagonal residues are broadcast afterwards.
    int selected_kind=0; // 1=joint, 2=scalar-pole union, 3=balanced joint
    int support_count=0, raw_pole_count=0, refinement_iterations=0;
    double barycentric_max_error=0.0;
    std::string root_error;

    // A channel whose norm is less than this fraction of the largest channel
    // is treated as numerically zero.  This is deliberately an amplitude
    // threshold (its square is used below for accumulated norms): a small but
    // physical anomalous/Nambu bath is still normalized and protected, while
    // roundoff in an exactly zero channel is not amplified into a fit target.
    constexpr double channel_scale_floor_ratio=1.e-12;

    const auto expansion_error=[&](const edmft::aaa::PoleExpansion& expansion) {
      double error2=0.0,scale2=0.0,max_error=0.0;
      for (int n=0; n<static_cast<int>(z.size()); ++n) {
        const MatrixXcd d=expansion.evaluate(z[n])-values[n];
        error2+=d.squaredNorm();
        scale2+=values[n].squaredNorm();
        max_error=std::max(max_error,d.norm());
      }
      return std::make_pair(
        std::sqrt(error2/std::max(scale2,std::numeric_limits<double>::min())),
        max_error);
    };

    // A global Frobenius error favors the largest normal-state bath elements.
    // Measure every independent CIX channel in its own norm as well, so a
    // smaller anomalous hybridization can veto a poor joint approximation.
    // A negative result marks a channel that is zero at the scale above and
    // is intentionally excluded from the worst-channel comparison.
    const auto channel_relative_errors=[&](
        const edmft::aaa::PoleExpansion& expansion) {
      const int channel_count=Deltaw.extent(0);
      std::vector<double> error2(channel_count,0.0);
      std::vector<double> scale2(channel_count,0.0);
      for (int n=0; n<static_cast<int>(z.size()); ++n) {
        const MatrixXcd approximation=expansion.evaluate(z[n]);
        for (int b=0; b<channel_count; ++b) {
          const Complex d=approximation(b,b)-values[n](b,b);
          error2[b]+=std::norm(d);
          scale2[b]+=std::norm(values[n](b,b));
        }
      }
      const double largest_scale2=*std::max_element(
        scale2.begin(),scale2.end());
      const double scale2_floor=largest_scale2*
        channel_scale_floor_ratio*channel_scale_floor_ratio;
      std::vector<double> relative(channel_count,-1.0);
      for (int b=0; b<channel_count; ++b)
        if (scale2[b]>std::max(scale2_floor,
                              std::numeric_limits<double>::min()))
          relative[b]=std::sqrt(error2[b]/scale2[b]);
      return relative;
    };

    const auto worst_channel_error=[](const std::vector<double>& errors) {
      double worst=0.0;
      bool found=false;
      for (double error:errors)
        if (error>=0.0) {
          worst=std::max(worst,error);
          found=true;
        }
      return found ? worst : std::numeric_limits<double>::infinity();
    };

    if (mpi.rank==mpi.master) {
      try {
        const double pole_imag_tolerance=
          pars["aaa_pole_imag_tolerance"].get<double>();
        const bool refine=pars["aaa_refine_poles"].get<bool>();
        double bound=pars["aaa_omega_bound"].get<double>();
        if (bound<=0.0) bound=pars["cutoffw"].get<double>();

        const auto finish_expansion=[&](
            const std::vector<MatrixXcd>& fit_values,
            const Eigen::VectorXd& initial_poles,
            int& iterations) {
          if (initial_poles.size()==0) return edmft::aaa::PoleExpansion{};
          if (!refine)
            return edmft::aaa::fit_residues(z,fit_values,initial_poles);
          const double local_bound=std::max(
            bound,1.05*initial_poles.cwiseAbs().maxCoeff());
          auto report=edmft::aaa::refine_poles(
            z,fit_values,initial_poles,-local_bound,local_bound,
            pars["aaa_refine_max_iterations"].get<int>(),
            pars["aaa_refine_gradient_tolerance"].get<double>());
          iterations+=report.iterations;
          return std::move(report.expansion);
        };

        edmft::aaa::PoleExpansion joint_expansion;
        double joint_relative_error=std::numeric_limits<double>::infinity();
        if (requested_mode!="balanced") {
          aaa_result=edmft::aaa::fit(z,values,opts);
          support_count=static_cast<int>(aaa_result.support_points.size());
          raw_pole_count=aaa_result.poles.size();
          barycentric_max_error=aaa_result.max_error_history.empty()
            ? 0.0 : aaa_result.max_error_history.back();
          Eigen::VectorXd poles=aaa_result.real_poles(pole_imag_tolerance);
          int joint_iterations=0;
          joint_expansion=finish_expansion(values,poles,joint_iterations);
          refinement_iterations=joint_iterations;
          if (joint_expansion.poles.size()>0)
            joint_relative_error=expansion_error(joint_expansion).first;
        }

        // Balanced matrix AAA keeps a single denominator and hence a compact,
        // well-conditioned real-pole set.  Channel normalization is used only
        // while discovering/refining those poles.  The final residues are
        // always fitted to the original, unscaled hybridization.
        const bool try_balanced=requested_mode=="balanced" ||
          (requested_mode=="auto" &&
           joint_relative_error>pars["aaa_channelwise_trigger"].get<double>());
        edmft::aaa::PoleExpansion balanced_expansion;
        int balanced_support_count=0,balanced_raw_poles=0;
        int balanced_iterations=0;
        double balanced_barycentric_error=0.0;
        if (try_balanced) {
          std::vector<double> channel_scale(Deltaw.extent(0),0.0);
          for (const auto& value:values)
            for (int b=0; b<Deltaw.extent(0); ++b)
              channel_scale[b]+=std::norm(value(b,b));
          for (double& scale:channel_scale)
            scale=std::sqrt(scale/std::max<size_t>(values.size(),1));
          const double largest_channel_scale=*std::max_element(
            channel_scale.begin(),channel_scale.end());
          const double channel_scale_floor=std::max(
            channel_scale_floor_ratio*largest_channel_scale,
            std::numeric_limits<double>::min());
          std::vector<MatrixXcd> balanced_values;
          balanced_values.reserve(values.size());
          for (const auto& value:values) {
            MatrixXcd balanced=MatrixXcd::Zero(
              Deltaw.extent(0),Deltaw.extent(0));
            for (int b=0; b<Deltaw.extent(0); ++b)
              balanced(b,b)=value(b,b)/
                std::max(channel_scale[b],channel_scale_floor);
            balanced_values.push_back(std::move(balanced));
          }
          auto balanced_result=edmft::aaa::fit(z,balanced_values,opts);
          balanced_support_count=balanced_result.support_points.size();
          balanced_raw_poles=balanced_result.poles.size();
          balanced_barycentric_error=
            balanced_result.max_error_history.empty()
            ? 0.0 : balanced_result.max_error_history.back();
          const Eigen::VectorXd poles=
            balanced_result.real_poles(pole_imag_tolerance);
          const auto scaled_expansion=finish_expansion(
            balanced_values,poles,balanced_iterations);
          if (scaled_expansion.poles.size()>0)
            balanced_expansion=edmft::aaa::fit_residues(
              z,values,scaled_expansion.poles);
        }

        // The scalar-union construction is intentionally available only when
        // requested explicitly.  A large union can interpolate Matsubara data
        // extremely well while being ill conditioned between sample points;
        // it is useful for diagnosis but unsafe as an automatic fallback.
        const bool try_channelwise=requested_mode=="channelwise";
        edmft::aaa::PoleExpansion channel_expansion;
        int channel_support_count=0,channel_raw_poles=0,channel_iterations=0;
        double channel_barycentric_error=0.0;
        if (try_channelwise) {
          std::vector<double> channel_poles;
          // Fit each physical scalar channel in its natural one-dimensional
          // space.  This prevents a large normal component from controlling
          // support selection for a much smaller anomalous component.
          for (int b=0; b<Deltaw.extent(0); ++b) {
            std::vector<MatrixXcd> scalar_values;
            scalar_values.reserve(values.size());
            for (const auto& value:values) {
              MatrixXcd scalar(1,1);
              scalar(0,0)=value(b,b);
              scalar_values.push_back(std::move(scalar));
            }
            auto scalar_result=edmft::aaa::fit(z,scalar_values,opts);
            channel_support_count+=scalar_result.support_points.size();
            channel_raw_poles+=scalar_result.poles.size();
            if (!scalar_result.max_error_history.empty())
              channel_barycentric_error=std::max(
                channel_barycentric_error,
                scalar_result.max_error_history.back());
            Eigen::VectorXd poles=scalar_result.real_poles(pole_imag_tolerance);
            int scalar_iterations=0;
            const auto scalar_expansion=
              finish_expansion(scalar_values,poles,scalar_iterations);
            channel_iterations+=scalar_iterations;
            for (int l=0; l<scalar_expansion.poles.size(); ++l)
              channel_poles.push_back(scalar_expansion.poles(l));
          }

          // Equivalent channels should produce identical poles analytically,
          // but their independent finite-precision fits can differ slightly.
          // Merge such copies before the final simultaneous residue fit.
          std::sort(channel_poles.begin(),channel_poles.end());
          std::vector<double> merged;
          std::vector<int> multiplicity;
          const double merge_tolerance=
            pars["aaa_pole_merge_tolerance"].get<double>();
          for (double pole:channel_poles) {
            if (merged.empty() ||
                std::abs(pole-merged.back())>
                merge_tolerance*std::max({1.0,std::abs(pole),
                                           std::abs(merged.back())})) {
              merged.push_back(pole);
              multiplicity.push_back(1);
            } else {
              const int count=++multiplicity.back();
              merged.back()+= (pole-merged.back())/count;
            }
          }
          Eigen::VectorXd union_poles(merged.size());
          for (int l=0; l<union_poles.size(); ++l) union_poles(l)=merged[l];
          if (union_poles.size()>0)
            channel_expansion=edmft::aaa::fit_residues(z,values,union_poles);
        }

        const double balanced_relative_error=balanced_expansion.poles.size()>0
          ? expansion_error(balanced_expansion).first
          : std::numeric_limits<double>::infinity();
        const double channel_relative_error=channel_expansion.poles.size()>0
          ? expansion_error(channel_expansion).first
          : std::numeric_limits<double>::infinity();
        const std::vector<double> joint_channel_errors=
          joint_expansion.poles.size()>0
          ? channel_relative_errors(joint_expansion) : std::vector<double>{};
        const std::vector<double> balanced_channel_errors=
          balanced_expansion.poles.size()>0
          ? channel_relative_errors(balanced_expansion)
          : std::vector<double>{};
        const double joint_worst_channel=
          joint_channel_errors.empty()
          ? std::numeric_limits<double>::infinity()
          : worst_channel_error(joint_channel_errors);
        const double balanced_worst_channel=
          balanced_channel_errors.empty()
          ? std::numeric_limits<double>::infinity()
          : worst_channel_error(balanced_channel_errors);

        // Preserve the existing global-error decision, but add a minimax
        // safeguard for small physical channels.  Requiring a 5% reduction
        // avoids changing pole sets because of last-bit differences between
        // nearly equivalent fits.  In particular, auto may retain the compact
        // joint fit only when it is also competitive channel by channel.
        const bool balanced_improves_global=
          balanced_relative_error<joint_relative_error;
        const bool balanced_improves_worst_channel=
          balanced_worst_channel<0.95*joint_worst_channel;
        const bool balanced_auto_preferred=requested_mode=="auto" &&
          (balanced_improves_global || balanced_improves_worst_channel);
        const bool channelwise_safe=try_channelwise &&
          channel_relative_error<joint_relative_error &&
          channel_expansion.poles.size()<=
            std::max(2*joint_expansion.poles.size(),
                     joint_expansion.poles.size()+8);
        if (channelwise_safe) {
          aaa_hybridization=std::move(channel_expansion);
          selected_kind=2;
          support_count=channel_support_count;
          raw_pole_count=channel_raw_poles;
          refinement_iterations=channel_iterations;
          barycentric_max_error=channel_barycentric_error;
        } else if (try_balanced &&
                   (requested_mode=="balanced" ||
                    balanced_auto_preferred)) {
          aaa_hybridization=std::move(balanced_expansion);
          selected_kind=3;
          support_count=balanced_support_count;
          raw_pole_count=balanced_raw_poles;
          refinement_iterations=balanced_iterations;
          barycentric_max_error=balanced_barycentric_error;
        } else {
          aaa_hybridization=std::move(joint_expansion);
          selected_kind=1;
          if (try_channelwise)
            (*mpi.out) << "WARN: rejected ill-conditioned channelwise AAA "
                       << "candidate (joint_poles="
                       << aaa_hybridization.poles.size()
                       << ", channel_union_poles="
                       << channel_expansion.poles.size()
                       << ", joint_relative_rms=" << joint_relative_error
                       << ", channel_relative_rms=" << channel_relative_error
                       << ")" << std::endl;
        }
        if (requested_mode=="auto" && try_balanced) {
          (*mpi.out) << "AAA auto per-channel relative RMS:";
          const int diagnostic_channels=static_cast<int>(std::max(
            joint_channel_errors.size(),balanced_channel_errors.size()));
          for (int b=0; b<diagnostic_channels; ++b) {
            const double joint_error=
              b<static_cast<int>(joint_channel_errors.size())
              ? joint_channel_errors[b]
              : std::numeric_limits<double>::infinity();
            const double balanced_error=
              b<static_cast<int>(balanced_channel_errors.size())
              ? balanced_channel_errors[b]
              : std::numeric_limits<double>::infinity();
            (*mpi.out) << " ch" << b << "=(joint="
                       << joint_error << ",balanced="
                       << balanced_error << ")";
          }
          (*mpi.out) << " worst=(joint=" << joint_worst_channel
                     << ",balanced=" << balanced_worst_channel << ")"
                     << " selected="
                     << (selected_kind==3 ? "balanced" : "joint")
                     << " reason="
                     << (balanced_improves_global ? "global-error" :
                         balanced_improves_worst_channel
                           ? "worst-channel-safeguard" :
                             "joint-remains-competitive")
                     << std::endl;
        }
        if (aaa_hybridization.poles.size()>0) {
          const auto errors=expansion_error(aaa_hybridization);
          aaa_relative_rms_error=errors.first;
          aaa_max_error=errors.second;
        }
      } catch (const std::exception& e) {
        root_error=e.what();
        aaa_hybridization=edmft::aaa::PoleExpansion{};
      }
    }

    int header[7]={root_error.empty() ? 1 : -1,
                   static_cast<int>(aaa_hybridization.poles.size()),
                   Deltaw.extent(0),selected_kind,support_count,
                   raw_pole_count,refinement_iterations};
    double diagnostics[3]={barycentric_max_error,
                           aaa_relative_rms_error,aaa_max_error};
#ifdef _MPI
    MPI_Bcast(header,7,MPI_INT,mpi.master,MPI_COMM_WORLD);
    MPI_Bcast(diagnostics,3,MPI_DOUBLE,mpi.master,MPI_COMM_WORLD);
#endif
    if (header[0]<0) {
      if (mpi.rank==mpi.master)
        (*mpi.out) << "AAA construction failed: " << root_error << std::endl;
      throw std::runtime_error("AAA construction failed on the master MPI rank");
    }
    selected_kind=header[3]; support_count=header[4];
    raw_pole_count=header[5]; refinement_iterations=header[6];
    barycentric_max_error=diagnostics[0];
    aaa_relative_rms_error=diagnostics[1]; aaa_max_error=diagnostics[2];

    // Broadcast only the physically used diagonal residue entries.  The
    // auxiliary dense matrices are reconstructed locally and contain exact
    // zero off-diagonal elements.
    const int pole_count=header[1],channel_count=header[2];
    std::vector<double> pole_buffer(pole_count);
    std::vector<Complex> residue_buffer(
      static_cast<size_t>(pole_count)*channel_count);
    if (mpi.rank==mpi.master) {
      for (int l=0; l<pole_count; ++l) {
        pole_buffer[l]=aaa_hybridization.poles(l);
        for (int b=0; b<channel_count; ++b)
          residue_buffer[static_cast<size_t>(l)*channel_count+b]=
            aaa_hybridization.residues[l](b,b);
      }
    }
#ifdef _MPI
    MPI_Bcast(pole_buffer.data(),pole_count,MPI_DOUBLE,mpi.master,MPI_COMM_WORLD);
    MPI_Bcast(residue_buffer.data(),pole_count*channel_count,
              MPI_C_DOUBLE_COMPLEX,mpi.master,MPI_COMM_WORLD);
#endif
    if (mpi.rank!=mpi.master) {
      aaa_hybridization.poles.resize(pole_count);
      aaa_hybridization.residues.assign(
        pole_count,MatrixXcd::Zero(channel_count,channel_count));
      for (int l=0; l<pole_count; ++l) {
        aaa_hybridization.poles(l)=pole_buffer[l];
        for (int b=0; b<channel_count; ++b)
          aaa_hybridization.residues[l](b,b)=
            residue_buffer[static_cast<size_t>(l)*channel_count+b];
      }
    }

    if (pole_count==0) {
      if (mpi.rank==mpi.master)
        (*mpi.out) << "WARN: AAA adapter found no usable real poles; "
                   << "the existing dense hybridization remains in use."
                   << std::endl;
      return;
    }

    // Build the optional TCA-only reduced pole table after the accurate AAA
    // expansion has been agreed by every MPI rank.  Greedy deletion/refitting
    // runs only on rank zero; broadcasting the result avoids both duplicated
    // startup work and rank-dependent last-bit choices.  A target of zero, or
    // a target no smaller than the accurate table, is an exact copy.
    const int requested_tca_poles=pars["tca_active_aaa_poles"].get<int>();
    const double requested_keep_weight=
      pars["tca_pseudo_keep_weight"].get<double>();
    if (requested_tca_poles<0)
      throw std::invalid_argument("tca_active_aaa_poles must be nonnegative");
    if (!(requested_keep_weight>0.0 && requested_keep_weight<=1.0))
      throw std::invalid_argument(
        "tca_pseudo_keep_weight must lie in (0,1]");
    const int effective_tca_poles=requested_tca_poles==0
      ? pole_count : std::min(requested_tca_poles,pole_count);
    double tca_compression_diagnostics[2]={aaa_relative_rms_error,0.0};
    int tca_compression_ok=1;
    if (mpi.rank==mpi.master) {
      try {
        if (effective_tca_poles>0 && effective_tca_poles<pole_count) {
          auto compressed=edmft::aaa::compress_poles_greedy(
            z,values,aaa_hybridization,effective_tca_poles);
          tca_active_hybridization=std::move(compressed.expansion);
          tca_compression_diagnostics[0]=compressed.relative_rms_error;
          tca_compression_diagnostics[1]=
            compressed.worst_channel_relative_rms;
        } else {
          tca_active_hybridization=aaa_hybridization;
          const auto channel_errors=channel_relative_errors(
            tca_active_hybridization);
          tca_compression_diagnostics[1]=
            worst_channel_error(channel_errors);
        }
      } catch (const std::exception& e) {
        tca_compression_ok=0;
        (*mpi.out) << "TCA active AAA compression failed: " << e.what()
                   << std::endl;
      }
    }
#ifdef _MPI
    MPI_Bcast(&tca_compression_ok,1,MPI_INT,mpi.master,MPI_COMM_WORLD);
    MPI_Bcast(tca_compression_diagnostics,2,MPI_DOUBLE,mpi.master,
              MPI_COMM_WORLD);
#endif
    if (!tca_compression_ok)
      throw std::runtime_error(
        "TCA active AAA compression failed on the master MPI rank");
    const int tca_pole_count=effective_tca_poles;
    std::vector<double> tca_pole_buffer(tca_pole_count);
    std::vector<Complex> tca_residue_buffer(
      static_cast<size_t>(tca_pole_count)*channel_count);
    if (mpi.rank==mpi.master)
      for (int l=0;l<tca_pole_count;++l) {
        tca_pole_buffer[l]=tca_active_hybridization.poles(l);
        for (int b=0;b<channel_count;++b)
          tca_residue_buffer[static_cast<size_t>(l)*channel_count+b]=
            tca_active_hybridization.residues[l](b,b);
      }
#ifdef _MPI
    MPI_Bcast(tca_pole_buffer.data(),tca_pole_count,MPI_DOUBLE,mpi.master,
              MPI_COMM_WORLD);
    MPI_Bcast(tca_residue_buffer.data(),tca_pole_count*channel_count,
              MPI_C_DOUBLE_COMPLEX,mpi.master,MPI_COMM_WORLD);
#endif
    if (mpi.rank!=mpi.master) {
      tca_active_hybridization.poles.resize(tca_pole_count);
      tca_active_hybridization.residues.assign(
        tca_pole_count,MatrixXcd::Zero(channel_count,channel_count));
      for (int l=0;l<tca_pole_count;++l) {
        tca_active_hybridization.poles(l)=tca_pole_buffer[l];
        for (int b=0;b<channel_count;++b)
          tca_active_hybridization.residues[l](b,b)=
            tca_residue_buffer[static_cast<size_t>(l)*channel_count+b];
      }
    }

    if (mpi.rank==mpi.master)
      (*mpi.out) << "AAA hybridization adapter: channels=" << channel_count
                 << " mode=" << (selected_kind==2 ? "channelwise" :
                                  selected_kind==3 ? "balanced" : "joint")
                 << " support_points=" << support_count
                 << " barycentric_max_error=" << barycentric_max_error
                 << " raw_poles=" << raw_pole_count
                 << " real_poles=" << pole_count
                 << " refinement_iterations=" << refinement_iterations
                 << " pole_fit_relative_rms=" << aaa_relative_rms_error
                 << " pole_fit_max_error=" << aaa_max_error << std::endl;
    if (diagram_dlr_workspace && aaa_hybridization.poles.size()>0)
      diagram_dlr_workspace->prepare_hybridization(*dlr_basis,aaa_hybridization,
                                                   cix.Nbath_unique_unique);
    const bool active_workspace_requested=
      tca_pole_count<pole_count || requested_keep_weight<1.0;
    if (active_workspace_requested) {
      tca_active_dlr_workspace=
        std::make_unique<edmft::imag_dlr::Workspace>();
      tca_active_dlr_workspace->prepare_hybridization(
        *dlr_basis,tca_active_hybridization,cix.Nbath_unique_unique);
      if (mpi.rank==mpi.master)
        (*mpi.out) << "TCA active workspace: poles=" << tca_pole_count
                   << "/" << pole_count
                   << " pseudo_keep_weight=" << requested_keep_weight
                   << " pole_fit_relative_rms="
                   << tca_compression_diagnostics[0]
                   << " worst_channel_relative_rms="
                   << tca_compression_diagnostics[1] << std::endl;
    }
  }

  Eigen::MatrixXd PseudoToDLRCoefficients(const Array<double,1>& dense_tau,
                                           const Array<double,2>& values) const
  {
    if (!dlr_basis) throw std::logic_error("DLR adapter is not enabled");
    if (dense_tau.extent(0)!=values.extent(1))
      throw std::invalid_argument("PseudoToDLRCoefficients: tau/value size mismatch");
    // Input pseudo functions remain on the solver's excellent dense spline
    // mesh.  This least-squares projection is the boundary between that
    // representation and compact diagram evaluation; Dyson is untouched.
    Eigen::VectorXd etau(dense_tau.extent(0));
    Eigen::MatrixXd evals(dense_tau.extent(0),values.extent(0));
    for (int it=0; it<dense_tau.extent(0); ++it){
      etau(it)=dense_tau(it);
      for (int p=0; p<values.extent(0); ++p) evals(it,p)=values(p,it);
    }
    return dlr_basis->fit_values(etau,evals);
  }

  void DLRToPseudo(const Eigen::MatrixXd& coefficients,
                   const Array<double,1>& dense_tau,
                   Array<double,2>& values) const
  {
    if (!dlr_basis) throw std::logic_error("DLR adapter is not enabled");
    if (coefficients.rows()!=dlr_basis->rank())
      throw std::invalid_argument("DLRToPseudo: coefficient count differs from DLR rank");
    values.resize(coefficients.cols(),dense_tau.extent(0));
    for (int it=0; it<dense_tau.extent(0); ++it){
      const auto row=dlr_basis->evaluate(dense_tau(it),coefficients);
      for (int p=0; p<coefficients.cols(); ++p) values(p,it)=row(p);
    }
  }

  void UpdateDLRPseudos(const Array<double,1>& dense_tau,
                        const Array<double,2>& Gm,
                        const Array<double,2>& gm,
                        MPIcls& mpi)
  {
    if (!dlr_adapter_enabled) return;
    dlr_G_coefficients=PseudoToDLRCoefficients(dense_tau,Gm);
    dlr_g_coefficients=PseudoToDLRCoefficients(dense_tau,gm);
    if (diagram_dlr_workspace && diagram_dlr_workspace->ready())
      diagram_dlr_workspace->update_pseudos(dlr_G_coefficients);
    // Round-trip both G and the shifted/scaled g used elsewhere in the solver
    // even though OCA currently consumes G only.  Reporting this error at
    // every update separates DLR representation error from diagram algebra
    // or legacy quadrature error.
    Array<double,2> Gcheck,gcheck;
    DLRToPseudo(dlr_G_coefficients,dense_tau,Gcheck);
    DLRToPseudo(dlr_g_coefficients,dense_tau,gcheck);
    double Gerr2=0.0,Gscale2=0.0,gerr2=0.0,gscale2=0.0;
    for (int p=0; p<Gm.extent(0); ++p)
      for (int it=0; it<Gm.extent(1); ++it){
        Gerr2+=std::pow(Gcheck(p,it)-Gm(p,it),2);
        Gscale2+=std::pow(Gm(p,it),2);
      }
    for (int p=0; p<gm.extent(0); ++p)
      for (int it=0; it<gm.extent(1); ++it){
        gerr2+=std::pow(gcheck(p,it)-gm(p,it),2);
        gscale2+=std::pow(gm(p,it),2);
      }
    dlr_G_relative_error=std::sqrt(Gerr2/std::max(Gscale2,std::numeric_limits<double>::min()));
    dlr_g_relative_error=std::sqrt(gerr2/std::max(gscale2,std::numeric_limits<double>::min()));
    if (mpi.rank==mpi.master)
      (*mpi.out) << "DLR pseudo adapter: rank=" << dlr_basis->rank()
                 << " G_relative_rms=" << dlr_G_relative_error
                 << " g_relative_rms=" << dlr_g_relative_error << std::endl;
  }

  bool HasAAAHybridization() const { return aaa_hybridization.poles.size()>0; }
  bool HasDLRPseudos() const { return dlr_G_coefficients.size()>0; }
  const edmft::aaa::PoleExpansion& AAAHybridization() const { return aaa_hybridization; }
  const edmft::dlr::FermionicBasis& DLRBasis() const {
    if (!dlr_basis) throw std::logic_error("DLR adapter is not enabled");
    return *dlr_basis;
  }

  auto ReadPseudoSig00(string Sfilename, MPIcls& mpi){
    bool debug = true;
    auto retval = _ReadPseudoSig00_<double>(Sfilename, cix.Npseudo_unique, mpi);
    if (retval.Sinitialized){
      double beta_S00 = retval.xn(retval.xn.size()-1);
      double beta = tau(tau.size()-1);
      double stretch = 1.0;
      if (fabs(beta_S00-beta)>1e-5){
	(*mpi.out) <<"WARN: The read "<<Sfilename<<" has different temperature."<<endl;
	stretch = beta_S00/beta;
      }
      double diff = 0;
      if (retval.xn.size()!=tau.size()) diff = 1.0;
      else diff = sum(abs(retval.xn-tau)); // Array<double,1> xn;
      if (diff > 1e-5){
	(*mpi.out) <<"WARN: The read "<<Sfilename<<" has different mesh than Deltat.dat. We need to interpolate."<<endl;
	Spline1DMany fSig(retval.xn.extent(0),retval.Sm.extent(0));
	for (int k=0; k<retval.Sm.extent(0); ++k)
	  for (int i=0; i<retval.xn.extent(0); ++i) fSig(i,k) = retval.Sm(k,i);
	fSig.splineIt(retval.xn);
	Array<double,2> Sm_new(retval.Sm.extent(0),tau.extent(0));
	int pos=0;
	for (int i=0; i<tau.extent(0); ++i){
	  intpar ips = InterpLeft(tau(i)*stretch, pos,retval.xn);
	  for (int k=0; k<retval.Sm.extent(0); ++k)
	    Sm_new(k,i) = fSig(ips,k);
	}
	retval.xn.resize(tau.size());
	retval.xn = tau;
	retval.Sm.resize(Sm_new.shape());
	retval.Sm = Sm_new;
	if (debug){
	  ofstream out("PPSigma.0.interpolated"); out.precision(12);
	  for (int i=0; i<tau.extent(0); ++i){
	    out << tau(i) << " ";
	    for (int k=0; k<retval.Sm.extent(0); ++k)
	      out << retval.Sm(k,i) << " ";
	    out << endl;
	  }
	}
      }
    }
    return retval;
  }
  void InitializeSig00(Array<double,1>& _tau_, Array<double,2>& Sm,
		       const Array<double,1>& om, const Cix& cix, bool fastFilesystem, MPIcls& mpi)
  {
    CreateSigma000(Sm,
		   tau,Deltas,cix,beta,fastFilesystem,mpi);
    // Note that tau was created in PrepareDelta and is non-uniform doubleExp mesh.
    _tau_.resize(tau.size()); _tau_ = tau;
  }
  
  void InitializeStorage(int Npseudo_unique, const Array<double,1>& xn){};
  auto DeterminePseudos(Array<double,2>& Gm, Array<double,2>& gm,
			const Array<double,2>& Sm, const Array<double,1>& _tau_, const Cix& cix, double lambda0_fixed, int followPeak_fixed,
			int itt, const json& pars, MPIcls& mpi)
  {
    if (lambda0_fixed!=0 or followPeak_fixed>0){
      if (lambda0_fixed!=0) (*mpi.out)<<"lambda0_fixed="<<lambda0_fixed<<" ";
      if (followPeak_fixed>0) (*mpi.out)<<"followPeak_fixed="<<followPeak_fixed<<" ";
      (*mpi.out)<<endl;
    }
    
    auto [lambda0,followPeak] = DeterminePseudoPropagators(Gm,/**/tau,tau_dh,Sm,cix,beta,lambda0_fixed,followPeak_fixed,itt,mpi,pars,pars["peak_pos"]);
    
    FindLargestEigenvalues(gm,
			   Gm,tau,cix,itt,mpi);
    double mu = pars["mu"];
    auto r = DeterminProbabilitiesI(tau,Gm,cix,beta,mu,lambda0,followPeak,itt,mpi);
    return r;
  }
  
  vector<int> ImportanceOfPseudo(const Array<double,1>& pseudo_tau,
                                 const Array<double,2>& pseudo_G,
                                 std::array<double,2> cut={0.1,0.5},
                                 bool debug=false)
  {
    if (pseudo_G.extent(0)!=cix.Npseudo_unique ||
        pseudo_G.extent(1)!=pseudo_tau.extent(0))
      throw std::invalid_argument(
        "ImportanceOfPseudo: pseudo Green function/tau dimensions mismatch");
  
    Array<double,1> Gwidth(cix.Npseudo_unique);
    const int middle=pseudo_tau.extent(0)/2;
    // Use the original dense pseudo mesh directly.  The mesh is symmetric and
    // its middle point represents beta/2, so constructing Gms splines merely
    // to rank |G_p(beta/2)| was unnecessary.  This keeps DLR-only OCA entirely
    // independent of the legacy pseudo spline representation.
    for (int ipp=0; ipp<cix.Npseudo_unique; ++ipp)
      Gwidth(ipp)=fabs(pseudo_G(ipp,middle));
  
    std::vector<int> idx(cix.Npseudo_unique);
    std::iota(idx.begin(), idx.end(), 0); // 0..N-1
    std::sort(idx.begin(), idx.end(), [&](int a, int b){ return Gwidth(a) > Gwidth(b); });

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
  double NCA(Array<double,2>& Sm_new, Array<double,2>& G00t, Array<double,2>& F00t, Array<double,2>& F01t,
	     double cmpPhysG, double Q, const Array<double,1>& tau, const Array<double,2>& Gm, Cix& cix, MPIcls& mpi){
    double phi = I_NCA(Sm_new,G00t,F00t,F01t,
		       cmpPhysG,Q,beta,tau,tau_dh,Gm,Delta,cix,mpi);
    return phi;
  }
  double OCA(Array<double,2>& dSm, Array<double,2>& G00t, Array<double,2>& F00t, Array<double,2>& F01t,
	     bool cmpPhysG, double Q, const Array<double,1>& tau, const Array<double,2>& Gm, Cix& cix, json& pars, MPIcls& mpi){
    /* Algorithm visual summary:

       Current dense pseudo G
             │
             ├── fit pseudo G directly to DLR coefficients
             ├── compute G_importance directly on the dense mesh
             │
             ▼
       EnumerateOCAPaths (unless a fixed-order cache is reusable)
             │
             ├── pass 1: count unique CIX paths
             ├── pass 2: accumulate matrix-element weights
             └── one shared, mesh-independent enumeration
                         │
               ┌─────────┴─────────┐
               ▼                   ▼
       legacy requested       DLR requested
       I_OCA(enumeration)     build SelfEnergyPlan + PhysicalPlan
       uniform mesh +                    │
       spline integration                ▼
                               one pole-major loop:
                               shared pair convolutions
                                  ├── four cuts → Sigma_OCA
                                  ├── two cuts  → G,F00,F01
                                  └── Sigma*G identity → Phi_OCA

                               MPI execution:
                               poles → independent teams
                               pair columns → team Allgatherv
                               cut/physical blocks → team ranks
                               compact DLR nodes → world Allreduce

       Diagnostic mode follows both branches; OCA_solver determines which
       result is copied back to the self-consistency loop.
    */
    const std::string oca_solver=pars["OCA_solver"].get<std::string>();
    if (oca_solver!="legacy" && oca_solver!="dlr")
      throw std::invalid_argument(
        "OCA_solver must be either \"legacy\" or \"dlr\", received \""+
        oca_solver+"\"");
    const bool use_dlr=oca_solver=="dlr";

    // Diagnostic mode sends the shared enumeration through both numerical
    // branches and compares them on precisely the same external tau mesh.
    // With legacy
    // selected, diagnostics never overwrite the legacy output.  With DLR
    // selected, the compact arrays become the returned production output;
    // enabling diagnostics merely adds the legacy calculation and comparison.
    const bool dlr_diagnostic=pars["oca_dlr_diagnostic"].get<bool>();
    const bool evaluate_dlr=use_dlr || dlr_diagnostic;
    // In production DLR mode a diagnostic request means "compute both and
    // compare, but return DLR".  Without diagnostics the legacy integrator is
    // not called at all; only the shared mesh-independent enumerator runs.
    const bool run_legacy=!use_dlr || dlr_diagnostic;
    const bool dlr_cut_diagnostic=dlr_diagnostic && pars["oca_dlr_cut_diagnostic"].get<bool>();

    // pars["Norder"] is N_max_order in the driver.  Equal minimum and maximum
    // orders describe a fixed-order calculation.  In particular, order two
    // may still run a small complementary-OCA Monte Carlo, but retaining OCA
    // topology there is the intended speed/memory compromise.  When a range
    // of orders is sampled, all persistent OCA storage is released so the
    // substantially larger higher-order MC has the memory available.
    const bool retain_oca_topology=
      pars["N_min_order"].get<int>()==pars["Norder"].get<int>();

    const bool topology_mpi_enabled=
      pars["oca_dlr_mpi"].get<bool>() && mpi.size>1;
    const std::string requested_topology_setting=
      pars["oca_dlr_mpi_topology"].get<std::string>();
    // A serial evaluator cannot consume rank-local topology: forcing the
    // replicated path when oca_dlr_mpi is disabled keeps that switch a true
    // numerical-communication switch rather than silently dropping diagrams.
    const std::string topology_setting=topology_mpi_enabled
      ? requested_topology_setting : "replicated";
    const double replication_budget=
      pars["oca_dlr_max_replicated_mb"].get<double>();
    if (requested_topology_setting!="auto" &&
        requested_topology_setting!="replicated" &&
        requested_topology_setting!="distributed")
      throw std::invalid_argument(
        "oca_dlr_mpi_topology must be \"auto\", \"replicated\", or "
        "\"distributed\"");
    if (!(replication_budget>0.0))
      throw std::invalid_argument("oca_dlr_max_replicated_mb must be positive");
    if (oca_dlr_cached_topology_setting!=topology_setting ||
        oca_dlr_cached_replication_budget!=replication_budget) {
      oca_dlr_p3_plan.reset();
      oca_dlr_self_plan.reset();
      oca_dlr_physical_plan.reset();
      oca_dlr_topology_shard_count=1;
      oca_dlr_cached_topology_setting=topology_setting;
      oca_dlr_cached_replication_budget=replication_budget;
    }

    // Drop objects that the selected numerical path cannot use.  This also
    // handles changing OCA_solver while parameters are reread between
    // iterations.
    if (!run_legacy) oca_legacy_enumeration.reset();
    if (!evaluate_dlr) {
      oca_dlr_p3_plan.reset();
      oca_dlr_self_plan.reset();
      oca_dlr_physical_plan.reset();
    }
    if (!dlr_cut_diagnostic) oca_dlr_p3_plan.reset();
    if (!retain_oca_topology) {
      oca_legacy_enumeration.reset();
      oca_dlr_p3_plan.reset();
      oca_dlr_self_plan.reset();
      oca_dlr_physical_plan.reset();
    }

    // Invalid objects are destroyed before enumeration starts.  For a large
    // cluster this matters: otherwise an obsolete plan/map and its
    // replacement can coexist at the peak and nearly double memory use.
    const auto self_plan_matches=[&]() {
      return oca_dlr_self_plan &&
        oca_dlr_self_plan->matches(
          cix.Npseudo_unique,cix.Nbath_unique_unique,
          G_importance,OCA_cut_sum);
    };
    const auto physical_plan_matches=[&]() {
      return oca_dlr_physical_plan &&
        oca_dlr_physical_plan->matches(
          cix.Npseudo_unique,cix.Nbath_unique_unique,
          G_importance,OCA_cut_sum);
    };
    const auto p3_plan_matches=[&]() {
      return oca_dlr_p3_plan &&
        oca_dlr_p3_plan->matches(
          cix.Npseudo_unique,cix.Nbath_unique_unique,
          G_importance,OCA_cut_sum);
    };
    if (evaluate_dlr && !self_plan_matches()) oca_dlr_self_plan.reset();
    if (evaluate_dlr && !physical_plan_matches()) oca_dlr_physical_plan.reset();
    if (dlr_cut_diagnostic && !p3_plan_matches()) oca_dlr_p3_plan.reset();
    if (run_legacy && oca_legacy_enumeration &&
        !oca_legacy_enumeration->matches(
          cmpPhysG,cix,G_importance,OCA_cut_sum))
      oca_legacy_enumeration.reset();

    const bool self_ready=evaluate_dlr && self_plan_matches();
    const bool physical_ready=!cmpPhysG ||
      (physical_plan_matches() && self_ready &&
       edmft::oca_dlr::CombinedEvaluator::compatible(
         *oca_dlr_self_plan,*oca_dlr_physical_plan));
    const bool p3_ready=!dlr_cut_diagnostic || p3_plan_matches();
    const bool need_dlr_paths=evaluate_dlr &&
      (!self_ready || !physical_ready || !p3_ready);
    const bool legacy_ready=run_legacy && oca_legacy_enumeration;
    const bool need_enumeration=need_dlr_paths ||
      (run_legacy && !legacy_ready);
    if (!need_enumeration && mpi.rank==mpi.master)
      (*mpi.out) << "OCA topology cache hit: retained="
                 << (run_legacy ? "legacy-enumeration" : "dlr-plans")
                 << std::endl;
    if (dlr_diagnostic && mpi.rank==mpi.master && diagram_dlr_workspace &&
        diagram_dlr_workspace->ready()) {
      const auto& nodes=diagram_dlr_workspace->basis().time_nodes();
      const auto& delta_plus=diagram_dlr_workspace->delta_plus();
      const auto& delta_minus=diagram_dlr_workspace->delta_minus();
      for (int s=0; s<cix.Nbath_unique_unique; ++s) {
        double error2=0.0,scale2=0.0,max_error=0.0;
        for (int k=0; k<nodes.size(); ++k) {
          int pos=0;
          const intpar ip=InterpLeft(nodes(k),pos,tau);
          pos=0;
          const intpar im=InterpLeft(beta-nodes(k),pos,tau);
          const double reference_plus=Deltas(ip,s);
          const double reference_minus=Deltas(im,s);
          // FourierOfDelta stores the opposite sign from the standard
          // fermionic Lehmann kernel used by AAA/DLR.  OCA contains two
          // hybridization lines, so this convention cancels in the diagram.
          const double dp=delta_plus(k,s)+reference_plus;
          const double dm=delta_minus(k,s)+reference_minus;
          error2+=dp*dp+dm*dm;
          scale2+=reference_plus*reference_plus+reference_minus*reference_minus;
          max_error=std::max({max_error,std::abs(dp),std::abs(dm)});
        }
        (*mpi.out) << "AAA hybridization time diagnostic: channel=" << s
                   << " relative_rms="
                   << std::sqrt(error2/std::max(scale2,std::numeric_limits<double>::min()))
                   << " max_abs=" << max_error << std::endl;
      }
    }
    // Enumeration is independent of the numerical integration method.  A
    // fixed-order legacy calculation retains the full tables.  DLR retains
    // only its compiled plans: the raw hash maps and path records are merely
    // construction material and are never cached for DLR-only production.
    std::unique_ptr<OCAEnumeration> temporary_enumeration;
    OCAEnumeration* legacy_enumeration=nullptr;
    std::vector<OCAPathRecord> paths;
    if (need_enumeration) {
      auto fresh=std::make_unique<OCAEnumeration>();
      const bool use_disk_cache=need_dlr_paths && !run_legacy &&
        retain_oca_topology && pars["oca_dlr_disk_cache"].get<bool>();
      const std::string cache_base=
        pars["oca_dlr_cache_file"].get<std::string>();
      // Self-energy-only and measurement iterations contain different path
      // weights/pair requirements; keep both files instead of making them
      // overwrite and invalidate one another on alternating iterations.
      const std::string cache_file=cache_base.empty() ? std::string{} :
        cache_base+(cmpPhysG ? ".physical" : ".self");
      std::uint64_t cix_hash=0;
      // Avoid a metadata/read storm on a shared filesystem: only the master
      // hashes the CIX contents, then distributes the cache signature.
      if (use_disk_cache && mpi.rank==mpi.master)
        cix_hash=OCAFileHash(pars["cix"].get<std::string>());
#ifdef _MPI
      if (use_disk_cache)
        MPI_Bcast(&cix_hash,1,MPI_UINT64_T,mpi.master,MPI_COMM_WORLD);
#endif
      bool loaded_cache=false;
      if (use_disk_cache)
        loaded_cache=LoadOCAPathCache(
          cache_file,cix_hash,cmpPhysG,cix,G_importance,OCA_cut_sum,
          topology_setting,replication_budget,mpi,*fresh);
      if (!loaded_cache) {
        *fresh=EnumerateOCAPaths(
          cmpPhysG,cix,G_importance,OCA_cut_sum,run_legacy,need_dlr_paths,
          topology_setting,replication_budget,mpi);
        if (use_disk_cache) WriteOCAPathCache(cache_file,cix_hash,*fresh,mpi);
      }
      oca_dlr_plans_distributed=fresh->paths_distributed;
      oca_dlr_topology_shard_count=fresh->topology_shard_count;
      if (need_dlr_paths) paths=std::move(fresh->paths);
      if (run_legacy && retain_oca_topology) {
        oca_legacy_enumeration=std::move(fresh);
        legacy_enumeration=oca_legacy_enumeration.get();
      } else if (run_legacy) {
        temporary_enumeration=std::move(fresh);
        legacy_enumeration=temporary_enumeration.get();
      }
      // In DLR-only mode `fresh` dies here after its compact path records have
      // been moved out, releasing the much larger legacy maps before plan
      // construction begins.
    } else if (run_legacy) {
      legacy_enumeration=oca_legacy_enumeration.get();
    }

    std::array<Array<double,2>,4> cut_references;
    double legacy_phi=0.0;
    if (run_legacy) {
      if (!legacy_enumeration)
        throw std::logic_error("OCA legacy evaluation has no enumeration");
      legacy_phi=I_OCA(*legacy_enumeration,dSm,G00t,F00t,F01t,
                       cmpPhysG,Q,beta,tau,Gms,Deltas,cix,pars,mpi,
                       dlr_cut_diagnostic ? &cut_references : nullptr);
    } else {
      // DLR will overwrite the requested outputs below.  Zero-initialize all
      // arrays now so an iteration without physical measurements also has a
      // well-defined state, without entering any legacy spline code.
      dSm.resize(cix.Npseudo_unique,tau.extent(0));
      dSm=0.0;
      G00t.resize(cix.Nbath_unique,tau.extent(0));
      F00t.resize(cix.Nbath_unique,tau.extent(0));
      F01t.resize(cix.Nbath_unique,tau.extent(0));
      G00t=0.0;
      F00t=0.0;
      F01t=0.0;
      if (mpi.rank==mpi.master)
        (*mpi.out) << "OCA DLR: topology ready; "
                   << "legacy integrator not called" << std::endl;
    }

    // A non-retained legacy enumeration is no longer needed after I_OCA.
    // Release it before allocating DLR pair catalogs and contraction matrices.
    temporary_enumeration.reset();
    double selected_phi=legacy_phi;
    if (evaluate_dlr){
      if (!diagram_dlr_workspace || !diagram_dlr_workspace->ready() ||
          !diagram_dlr_workspace->has_pseudos()) {
        throw std::logic_error(
          "OCA DLR evaluation requested without prepared AAA hybridization "
          "and DLR pseudo-particle data");
      } else {
        // Legacy self-energy cuts divide by the degeneracy of their output
        // pseudo state.  Physical cuts instead receive Q/bath normalization
        // below, matching the two separate legacy code paths.
        std::vector<double> degeneracy(cix.Npseudo_unique);
        for (int p=0; p<cix.Npseudo_unique; ++p) degeneracy[p]=cix.deg_pp(p);
        if (dlr_cut_diagnostic){
          if (!oca_dlr_p3_plan ||
              !oca_dlr_p3_plan->matches(cix.Npseudo_unique,cix.Nbath_unique_unique,
                                        G_importance,OCA_cut_sum)) {
            oca_dlr_p3_plan=std::make_unique<edmft::oca_dlr::P3Plan>();
            oca_dlr_p3_plan->build(paths,cix.Npseudo_unique,cix.Nbath_unique_unique,
                                   G_importance,OCA_cut_sum);
          }
          const Eigen::MatrixXd compact=oca_dlr_p3_plan->evaluate(
            *diagram_dlr_workspace,degeneracy,tau);
          double error2=0.0,scale2=0.0,max_error=0.0;
          for (int p=0; p<cix.Npseudo_unique; ++p)
            for (int it=0; it<tau.extent(0); ++it) {
              const double d=compact(p,it)-cut_references[3](p,it);
              error2+=d*d;
              scale2+=cut_references[3](p,it)*cut_references[3](p,it);
              max_error=std::max(max_error,std::abs(d));
            }
          if (mpi.rank==mpi.master)
            (*mpi.out) << "OCA DLR p3 diagnostic: paths=" << oca_dlr_p3_plan->nonzero_path_count()
                       << " inner_pairs=" << oca_dlr_p3_plan->inner_pair_count()
                       << " outer_pairs=" << oca_dlr_p3_plan->outer_pair_count()
                       << " relative_rms="
                       << std::sqrt(error2/std::max(scale2,std::numeric_limits<double>::min()))
                       << " max_abs=" << max_error << std::endl;
        }

        // Building a plan performs only hash deduplication and sparse-matrix
        // assembly.  When physical quantities are requested, build both plans
        // before numerical evaluation so one pole-major loop can feed every
        // self-energy and physical cut from the same pair convolutions.
        const auto plan_start=std::chrono::steady_clock::now();
        if (!oca_dlr_self_plan ||
            !oca_dlr_self_plan->matches(cix.Npseudo_unique,cix.Nbath_unique_unique,
                                        G_importance,OCA_cut_sum)) {
          oca_dlr_self_plan=std::make_unique<edmft::oca_dlr::SelfEnergyPlan>();
          oca_dlr_self_plan->build(paths,cix.Npseudo_unique,cix.Nbath_unique_unique,
                                   G_importance,OCA_cut_sum);
        }
        if (cmpPhysG) {
          if (!oca_dlr_physical_plan ||
              !oca_dlr_physical_plan->matches(
                cix.Npseudo_unique,cix.Nbath_unique_unique,
                G_importance,OCA_cut_sum)) {
            oca_dlr_physical_plan=
              std::make_unique<edmft::oca_dlr::PhysicalPlan>();
            oca_dlr_physical_plan->build(
              paths,cix.Npseudo_unique,cix.Nbath_unique_unique,
              G_importance,OCA_cut_sum);
          }

          // cmpPhysG can be false on inexpensive intermediate iterations and
          // true only on measurement iterations.  A cached self plan built in
          // the former case does not contain physical-only Bulla pairs.  On
          // the first measurement iteration rebuild it from the current path
          // list so both sparse plans use exactly the same pair numbering.
          if (!edmft::oca_dlr::CombinedEvaluator::compatible(
                *oca_dlr_self_plan,*oca_dlr_physical_plan)) {
            oca_dlr_self_plan->build(
              paths,cix.Npseudo_unique,cix.Nbath_unique_unique,
              G_importance,OCA_cut_sum);
          }
        }
        const auto plan_end=std::chrono::steady_clock::now();
        // Every required plan now owns its compact topology.  Path records are
        // construction-only data; release their capacity before the numerical
        // evaluator allocates its potentially large intermediate matrices.
        std::vector<OCAPathRecord>().swap(paths);
        std::array<Eigen::MatrixXd,4> compact_cuts;
        Eigen::MatrixXd compact_all;
        edmft::oca_dlr::PhysicalResult compact_physical;
        const auto eval_start=std::chrono::steady_clock::now();
        // Both measurement and inexpensive self-energy-only iterations use
        // the same pole-range/additive-accumulator path.  Passing a null
        // physical plan simply omits G/F00/F01 storage and contractions.  A
        // MPI execution below assigns pole ranges to teams, subdivides a pole
        // when a team contains several ranks, and reduces compact nodes before
        // final interpolation.  The serial and parallel paths therefore share
        // one numerical implementation.
        const edmft::oca_dlr::PhysicalPlan* physical_for_iteration=
          cmpPhysG ? oca_dlr_physical_plan.get() : nullptr;
        edmft::oca_dlr::CombinedResult combined;
        edmft::oca_dlr::ParallelLayout parallel_layout;
#ifdef _MPI
        const bool parallel_dlr=
          pars["oca_dlr_mpi"].get<bool>() && mpi.size>1;
        if (parallel_dlr) {
          if (oca_dlr_plans_distributed)
            combined=edmft::oca_dlr::CombinedEvaluator::evaluate_sharded_mpi(
              *oca_dlr_self_plan,physical_for_iteration,
              *diagram_dlr_workspace,degeneracy,tau,MPI_COMM_WORLD,
              oca_dlr_topology_shard_count,
              pars["oca_dlr_mpi_teams"].get<int>(),
              &parallel_layout,
              dlr_cut_diagnostic ? &compact_cuts : nullptr,
              pars["oca_dlr_whole_poles"].get<bool>(),
              pars["oca_dlr_profile"].get<bool>(),mpi.out);
          else
            combined=edmft::oca_dlr::CombinedEvaluator::evaluate_mpi(
              *oca_dlr_self_plan,physical_for_iteration,
              *diagram_dlr_workspace,degeneracy,tau,
              MPI_COMM_WORLD,pars["oca_dlr_mpi_teams"].get<int>(),
              &parallel_layout,
              dlr_cut_diagnostic ? &compact_cuts : nullptr,
              pars["oca_dlr_whole_poles"].get<bool>(),
              pars["oca_dlr_profile"].get<bool>(),mpi.out);
        } else
#endif
        {
          combined=edmft::oca_dlr::CombinedEvaluator::evaluate(
            *oca_dlr_self_plan,physical_for_iteration,
            *diagram_dlr_workspace,degeneracy,tau,
            dlr_cut_diagnostic ? &compact_cuts : nullptr);
          parallel_layout.world_size=mpi.size;
          parallel_layout.owned_poles=diagram_dlr_workspace->pole_count();
        }
        compact_all=std::move(combined.sigma);
        if (cmpPhysG)
          compact_physical=std::move(combined.physical);
        const auto eval_end=std::chrono::steady_clock::now();

        // Phi is recovered from its skeleton homogeneity identity.  Every OCA
        // Phi diagram has four pseudo lines, so cutting all of them gives
        //
        //   Phi_OCA = (1/4Q) sum_p deg_p int_0^beta
        //                        Sigma_p(tau) G_p(beta-tau) d tau.
        //
        // compact_all already contains the sum of all four cuts and includes
        // the 1/deg_p convention of the solver self-energy.  Multiplying the
        // degeneracy back here reproduces the legacy Phi normalization.  This
        // final one-dimensional integral uses the accurate existing dense
        // pseudo mesh; the costly nested OCA integration remains eliminated.
        double compact_phi=0.0;
        if (tau_dh.extent(0)!=tau.extent(0))
          throw std::logic_error("OCA DLR Phi quadrature mesh mismatch");
        for (int p=0; p<cix.Npseudo_unique; ++p)
          for (int it=0; it<tau.extent(0); ++it)
            compact_phi += 0.25*degeneracy[p]*tau_dh(it)*
              compact_all(p,it)*Gm(p,tau.extent(0)-1-it);
        compact_phi/=Q;
        if (use_dlr) selected_phi=compact_phi;

        if (dlr_cut_diagnostic && mpi.rank==mpi.master){
          for (int cut=0; cut<4; ++cut){
            double cut_error2=0.0,cut_scale2=0.0,cut_max=0.0;
            for (int p=0; p<cix.Npseudo_unique; ++p)
              for (int it=0; it<tau.extent(0); ++it){
                const double d=compact_cuts[cut](p,it)-cut_references[cut](p,it);
                cut_error2+=d*d;
                cut_scale2+=cut_references[cut](p,it)*cut_references[cut](p,it);
                cut_max=std::max(cut_max,std::abs(d));
              }
            (*mpi.out) << "OCA DLR cut p" << cut << ": relative_rms="
                       << std::sqrt(cut_error2/std::max(cut_scale2,std::numeric_limits<double>::min()))
                       << " max_abs=" << cut_max << std::endl;
          }
        }
        unsigned long long plan_counts[3]={
          static_cast<unsigned long long>(oca_dlr_self_plan->nonzero_path_count()),
          static_cast<unsigned long long>(oca_dlr_self_plan->forward_pair_count()),
          static_cast<unsigned long long>(oca_dlr_self_plan->backward_pair_count())};
        double plan_seconds=
          std::chrono::duration<double>(plan_end-plan_start).count();
        double evaluate_seconds=
          std::chrono::duration<double>(eval_end-eval_start).count();
#ifdef _MPI
        if (mpi.size>1) {
          double local_seconds[2]={plan_seconds,evaluate_seconds};
          double maximum_seconds[2]={};
          MPI_Allreduce(local_seconds,maximum_seconds,2,MPI_DOUBLE,MPI_MAX,
                        MPI_COMM_WORLD);
          plan_seconds=maximum_seconds[0];
          evaluate_seconds=maximum_seconds[1];
          if (oca_dlr_plans_distributed) {
            // A hybrid shard is replicated among its pole workers.  Count each
            // compiled plan only on the cyclic-group leader (world ranks
            // 0..topology_shard_count-1), otherwise diagnostics would multiply
            // global paths and pair catalogs by the workers per shard.
            if (mpi.rank>=oca_dlr_topology_shard_count)
              std::fill(std::begin(plan_counts),std::end(plan_counts),0ULL);
            unsigned long long global_counts[3]={};
            MPI_Allreduce(plan_counts,global_counts,3,MPI_UNSIGNED_LONG_LONG,
                          MPI_SUM,MPI_COMM_WORLD);
            std::copy(std::begin(global_counts),std::end(global_counts),
                      std::begin(plan_counts));
          }
        }
#endif
        if (mpi.rank==mpi.master) {
          (*mpi.out) << "OCA DLR self-energy: selected="
                     << (use_dlr ? "dlr" : "legacy")
                     << " evaluator="
                     << (cmpPhysG ? "combined-pole-loop" : "self-only")
                     << " topology="
                     << (oca_dlr_plans_distributed
                           ? "hybrid-distributed" : "replicated")
                     << " paths=" << plan_counts[0]
                     << " pair_instances=[" << plan_counts[1]
                     << "," << plan_counts[2] << "]"
                     << " plan_seconds=" << plan_seconds
                     << " evaluate_seconds=" << evaluate_seconds
                     << " topology_shards="
                     << parallel_layout.topology_shard_count
                     << " workers_per_shard_range=["
                     << parallel_layout.min_topology_group_size << ","
                     << parallel_layout.max_topology_group_size << "]"
                     << " pole_teams_per_shard=" << parallel_layout.team_count
                     << " pole_team_size_range=[" << parallel_layout.min_team_size
                     << "," << parallel_layout.max_team_size << "]"
                     << " Phi=" << compact_phi << std::endl;
        }
        if (dlr_diagnostic) {
          double error2=0.0,scale2=0.0,max_error=0.0;
          for (int p=0; p<cix.Npseudo_unique; ++p)
            for (int it=0; it<tau.extent(0); ++it) {
              const double d=compact_all(p,it)-dSm(p,it);
              error2+=d*d;
              scale2+=dSm(p,it)*dSm(p,it);
              max_error=std::max(max_error,std::abs(d));
            }
          if (mpi.rank==mpi.master)
            (*mpi.out) << "OCA DLR all-cuts diagnostic: relative_rms="
                       << std::sqrt(error2/std::max(
                            scale2,std::numeric_limits<double>::min()))
                       << " max_abs=" << max_error
                       << " legacy_Phi=" << legacy_phi
                       << " dlr_Phi=" << compact_phi << std::endl;
        }

        if (use_dlr) {
          // Convert Eigen's [pseudo,time] result to the Blitz array consumed
          // by Update_Sm.  This is the actual production switch: from this
          // point onward the self-consistency loop sees the compact result.
          dSm.resize(cix.Npseudo_unique,tau.extent(0));
          for (int p=0; p<cix.Npseudo_unique; ++p)
            for (int it=0; it<tau.extent(0); ++it)
              dSm(p,it)=compact_all(p,it);
        }

        // The physical Green function and Bulla numerators are only evaluated
        // on iterations on which the legacy solver computed their reference
        // values.  The path list then also contains records whose ordinary
        // OCA matrix element is zero but whose N_jm (Bulla) matrix element is
        // nonzero; PhysicalPlan must retain those records for F00/F01.
        if (cmpPhysG) {
          // CombinedEvaluator produced these values during the same pole loop
          // as compact_all; no tilted propagators or pair convolutions are
          // recomputed here.
          const edmft::oca_dlr::PhysicalResult& physical=compact_physical;

          // Reproduce exactly the normalization and channel reconstruction at
          // the end of legacy I_OCA.  The compact plan intentionally returns
          // the unnormalized independent channels so it remains independent
          // of CIX storage conventions.
          Array<double,1> bath_degeneracy(cix.Nbath_unique_unique);
          bath_degeneracy=0.0;
          for (int iorb=0; iorb<cix.Nbath; ++iorb)
            bath_degeneracy(std::abs(cix.bath_sigind0(iorb)))+=1.0;
          Array<int,1> sigind(cix.Nbath_unique);
          Array<bool,1> superconducting(cix.Nbath_unique);
          superconducting=false;
          for (int iorb=0; iorb<cix.Nbath; ++iorb) {
            const int ib=cix.bath_sigind(iorb);
            sigind(ib)=cix.bath_sigind0(iorb);
            superconducting(ib)=cix.SC_bath(iorb);
          }

          std::array<double,3> observable_error2{0.0,0.0,0.0};
          std::array<double,3> observable_scale2{0.0,0.0,0.0};
          std::array<double,3> observable_max{0.0,0.0,0.0};
          for (int ib=0; ib<cix.Nbath_unique; ++ib) {
            const int s=std::abs(sigind(ib));
            const double bath_sign=sigind(ib)<0 ? -1.0 : 1.0;
            const double normalizer=1.0/(bath_degeneracy(s)*Q);
            const double bulla_factor=superconducting(ib) ? 2.0 : 1.0;
            for (int it=0; it<tau.extent(0); ++it) {
              const std::array<double,3> compact_value{
                physical.g(s,it)*normalizer*bath_sign,
                physical.f00(s,it)*normalizer*bulla_factor*bath_sign,
                physical.f01(s,it)*normalizer*bulla_factor*bath_sign};
              // Read all references before production DLR mode overwrites the
              // Blitz arrays.  In enumeration-only mode diagnostics are off,
              // so the zero placeholders returned by I_OCA are never treated
              // as references.
              const std::array<double,3> reference_value{
                G00t(ib,it),F00t(ib,it),F01t(ib,it)};
              if (dlr_diagnostic) {
                for (int observable=0; observable<3; ++observable) {
                  const double d=compact_value[observable]-reference_value[observable];
                  observable_error2[observable]+=d*d;
                  observable_scale2[observable]+=
                    reference_value[observable]*reference_value[observable];
                  observable_max[observable]=std::max(
                    observable_max[observable],std::abs(d));
                }
              }
              if (use_dlr) {
                G00t(ib,it)=compact_value[0];
                F00t(ib,it)=compact_value[1];
                F01t(ib,it)=compact_value[2];
              }
            }
          }
          if (mpi.rank==mpi.master) {
            constexpr const char* names[3]={"G","F00","F01"};
            (*mpi.out) << "OCA DLR physical: selected="
                       << (use_dlr ? "dlr" : "legacy")
                       << " records="
                       << oca_dlr_physical_plan->nonzero_record_count()
                       << " forward_pairs="
                       << oca_dlr_physical_plan->forward_pair_count()
                       << " backward_pairs="
                       << oca_dlr_physical_plan->backward_pair_count()
                       << " evaluation=shared-pole-loop"
                       << std::endl;
            if (dlr_diagnostic)
              for (int observable=0; observable<3; ++observable)
                (*mpi.out) << "OCA DLR physical " << names[observable]
                           << ": relative_rms="
                           << std::sqrt(observable_error2[observable]/
                                std::max(observable_scale2[observable],
                                         std::numeric_limits<double>::min()))
                           << " max_abs=" << observable_max[observable]
                           << std::endl;
          }
        }
      }
    }
    if (!retain_oca_topology) {
      // Higher-order MC follows this call.  Do not let OCA topology compete
      // with its diagram and measurement storage.  Fixed-order calculations
      // retain these plans (or the legacy enumeration) for the next iteration.
      oca_dlr_p3_plan.reset();
      oca_dlr_self_plan.reset();
      oca_dlr_physical_plan.reset();
      oca_legacy_enumeration.reset();
    }
    return selected_phi;
  }

  // Complete analytic third order (TCA): both Phi topologies, all six
  // pseudo-particle cuts, and—on measurement iterations—all three physical
  // hybridization cuts.  No scalar diagram catalog is materialized; every
  // fixed bath/orientation sector is contracted in atomic block matrices.
  double TCA(Array<double,2>& dSm,
             Array<double,2>& G00t,
             Array<double,2>& F00t,
             Array<double,2>& F01t,
             bool cmpPhysG,
             double Q,
             const Array<double,1>& output_tau,
             const Array<double,2>& Gm,
             Cix& cix,
             json& pars,
             MPIcls& mpi)
  {
    if (pars["TCA_solver"].get<std::string>()!="dlr")
      throw std::logic_error("analytic TCA called while TCA_solver is not dlr");
    if (!diagram_dlr_workspace || !diagram_dlr_workspace->ready() ||
        !diagram_dlr_workspace->has_pseudos())
      throw std::logic_error(
        "analytic TCA requires prepared AAA hybridization and DLR pseudos");
    if (!(Q>0.0)) throw std::invalid_argument("TCA requires positive Q");

    std::vector<double> degeneracy(cix.Npseudo_unique);
    for (int p=0;p<cix.Npseudo_unique;++p) degeneracy[p]=cix.deg_pp(p);
    const double progress_seconds=
      pars["tca_progress_seconds"].get<double>();
    if (progress_seconds<0.0)
      throw std::invalid_argument("tca_progress_seconds must be nonnegative");
    const auto begin=std::chrono::steady_clock::now();
    edmft::tca_dlr::SelfEnergyNodeResult compact;
    const int root_sum_min_baths=
      pars["tca_root_sum_min_baths"].get<int>();
    if (root_sum_min_baths<1)
      throw std::invalid_argument("tca_root_sum_min_baths must be positive");
    const bool root_sum=pars["tca_root_sum"].get<bool>() &&
                        cix.Nbath>=root_sum_min_baths;
    const bool store_self_diagnostics=
      pars["tca_store_self_diagnostics"].get<bool>();
    const bool topology_b_self_c3=pars["tca_self_c3"].get<bool>();
    const std::string tca_mpi_scheduler=
      pars["tca_dlr_mpi_scheduler"].get<std::string>();
    if (tca_mpi_scheduler!="dynamic" && tca_mpi_scheduler!="static")
      throw std::invalid_argument(
        "tca_dlr_mpi_scheduler must be 'dynamic' or 'static'");
    const bool dynamic_tca_mpi=tca_mpi_scheduler=="dynamic";
    const int tca_mpi_micro_batch=
      pars["tca_dlr_mpi_micro_batch"].get<int>();
    if (tca_mpi_micro_batch<0)
      throw std::invalid_argument(
        "tca_dlr_mpi_micro_batch must be nonnegative");
    if (mpi.rank==mpi.master && pars["tca_root_sum"].get<bool>() &&
        !root_sum && mpi.out)
      (*mpi.out) << "TCA DLR bath summation: using fixed-bath plans because "
                 << "Nbath=" << cix.Nbath << " < tca_root_sum_min_baths="
                 << root_sum_min_baths << std::endl;
    // Keep the exact evaluator in one callable so the same code is used for
    // the accurate baseline/audit and for the cheaper active workspace.  The
    // physical Green/Bulla functions are requested only from the accurate
    // workspace; they are evaluated once, on the final measurement step.
    const auto evaluate_workspace=[&](
        const edmft::imag_dlr::Workspace& workspace,
        bool evaluate_physical) {
      edmft::tca_dlr::SelfEnergyNodeResult result;
      if (root_sum) {
#ifdef _MPI
        if (pars["tca_dlr_mpi"].get<bool>() && mpi.size>1)
          result=edmft::tca_dlr::FullSummedRootSelfEnergyEvaluator::evaluate_mpi(
            cix,workspace,degeneracy,MPI_COMM_WORLD,
            store_self_diagnostics,
            mpi.out,progress_seconds,topology_b_self_c3,
            dynamic_tca_mpi,tca_mpi_micro_batch);
        else
#endif
          result=edmft::tca_dlr::FullSummedRootSelfEnergyEvaluator::evaluate_serial(
            cix,workspace,degeneracy,store_self_diagnostics,
            mpi.out,progress_seconds,topology_b_self_c3);

        // The physical-line grouping is a separate contraction: the cut line
        // becomes the measured external channel rather than an unexpanded
        // Delta.  The caller decides whether the final measurement uses the
        // full workspace (an exact audit) or the same pole-capped workspace
        // used during self-consistency.
        if (evaluate_physical) {
          edmft::tca_dlr::SelfEnergyNodeResult physical;
#ifdef _MPI
          if (pars["tca_dlr_mpi"].get<bool>() && mpi.size>1)
            physical=edmft::tca_dlr::FullSummedPhysicalEvaluator::evaluate_mpi(
              cix,workspace,MPI_COMM_WORLD,
              mpi.out,progress_seconds,dynamic_tca_mpi,
              tca_mpi_micro_batch);
          else
#endif
            physical=edmft::tca_dlr::FullSummedPhysicalEvaluator::evaluate_serial(
              cix,workspace,mpi.out,progress_seconds);
          result.g=std::move(physical.g);
          result.f00=std::move(physical.f00);
          result.f01=std::move(physical.f01);
          result.physical_enabled=true;
        }
      } else {
        if (topology_b_self_c3 && mpi.rank==mpi.master && mpi.out)
          (*mpi.out) << "TCA DLR topology-B self C3 candidate is available "
                     << "only with tca_root_sum=true at this stage; using the "
                     << "six-cut fixed-bath reference" << std::endl;
#ifdef _MPI
        if (pars["tca_dlr_mpi"].get<bool>() && mpi.size>1)
          result=edmft::tca_dlr::FullSelfEnergyEvaluator::evaluate_mpi(
            cix,workspace,degeneracy,MPI_COMM_WORLD,
            store_self_diagnostics,evaluate_physical,
            mpi.out,progress_seconds);
        else
#endif
          result=edmft::tca_dlr::FullSelfEnergyEvaluator::evaluate_serial(
            cix,workspace,degeneracy,
            store_self_diagnostics,evaluate_physical,
            mpi.out,progress_seconds);
      }
      return result;
    };

    const double keep_weight=pars["tca_pseudo_keep_weight"].get<double>();
    if (!(keep_weight>0.0 && keep_weight<=1.0))
      throw std::invalid_argument(
        "tca_pseudo_keep_weight must lie in (0,1]");
    const bool approximate_workspace_available=
      tca_active_dlr_workspace && tca_active_dlr_workspace->ready();
    // Topology/cut diagnostics are intended as exact algebra tests.  Running
    // them through a frozen correction would require storing twelve extra
    // tail matrices and obscures their purpose, so diagnostic runs stay exact.
    const bool use_active=approximate_workspace_available &&
                          !store_self_diagnostics;
    edmft::tca_dlr::SelfEnergyNodeResult active_compact;
    bool have_active_compact=false;
    bool using_frozen_tail=false;
    bool exact_audit=false;
    bool had_frozen_tail=false;

    if (use_active) {
      // Select the active pseudo set once from the OCA-converged propagators.
      // deg_p is included because it is restored in the Phi trace.  Keeping
      // the mask fixed across TCA iterations is essential: otherwise changes
      // caused by a moving cutoff would contaminate the self-consistency map.
      if (tca_active_pseudo_mask.empty() ||
          keep_weight!=tca_cached_keep_weight) {
        std::vector<std::pair<double,int>> importance;
        importance.reserve(cix.Npseudo_unique);
        double total=0.0;
        for (int p=0;p<cix.Npseudo_unique;++p) {
          double weight=0.0;
          for (int it=0;it<output_tau.extent(0);++it)
            weight+=tau_dh(it)*std::abs(Gm(p,it));
          weight*=degeneracy[p];
          total+=weight;
          importance.emplace_back(weight,p);
        }
        std::sort(importance.begin(),importance.end(),
          [](const auto& left,const auto& right) {
            return left.first>right.first;
          });
        tca_active_pseudo_mask.assign(cix.Npseudo_unique,0);
        double retained=0.0;
        int active_count=0;
        for (const auto& [weight,pseudo]:importance) {
          if (active_count>0 && retained>=keep_weight*total) break;
          tca_active_pseudo_mask[pseudo]=1;
          retained+=weight;
          ++active_count;
        }
        tca_cached_keep_weight=keep_weight;
        tca_frozen_tail_ready=false;
        int active_blocks=0;
        for (int block=0;block<cix.Nblock;++block) {
          bool block_active=false;
          const int n=cix.Block[block].size();
          for (int i=0;i<n && !block_active;++i)
            for (int j=0;j<n && !block_active;++j) {
              const int label=cix.Sgm_symmetry[block](i,j);
              block_active=label!=0 &&
                tca_active_pseudo_mask[std::abs(label)-1];
            }
          active_blocks+=block_active ? 1 : 0;
        }
        if (mpi.rank==mpi.master && mpi.out)
          (*mpi.out) << "TCA active pseudo space: active=" << active_count
                     << "/" << cix.Npseudo_unique
                     << " active_blocks=" << active_blocks
                     << "/" << cix.Nblock
                     << " requested_weight=" << keep_weight
                     << " retained_weight="
                     << (total>0.0 ? retained/total : 1.0) << std::endl;
      }

      Eigen::MatrixXd active_coefficients=dlr_G_coefficients;
      for (int p=0;p<cix.Npseudo_unique;++p)
        if (!tca_active_pseudo_mask[p])
          active_coefficients.col(p).setZero();
      tca_active_dlr_workspace->update_pseudos(active_coefficients);
      if (mpi.rank==mpi.master && mpi.out)
        (*mpi.out) << "TCA active self evaluation: poles="
                   << tca_active_dlr_workspace->pole_count() << "/"
                   << diagram_dlr_workspace->pole_count() << std::endl;
      const bool freeze_tail=pars["tca_freeze_weak_tail"].get<bool>();
      const bool final_exact_audit=
        pars["tca_final_exact_audit"].get<bool>();
      exact_audit=(cmpPhysG && final_exact_audit) ||
                  (freeze_tail && !tca_frozen_tail_ready);
      // When the final exact audit is disabled, evaluate the physical
      // contractions in the active TCA workspace as well.  OCA is unaffected:
      // it always keeps the full hybridization pole expansion.
      active_compact=evaluate_workspace(
        *tca_active_dlr_workspace,cmpPhysG && !exact_audit);
      have_active_compact=true;
      // The cut propagator is absent from the internal chain.  Explicitly
      // zeroing its inactive output column completes the restriction that all
      // six pseudo-particle lines of the parent Phi diagram are active.
      for (int p=0;p<cix.Npseudo_unique;++p)
        if (!tca_active_pseudo_mask[p]) {
          active_compact.sigma.col(p).setZero();
          for (auto& value:active_compact.sigma_by_topology)
            value.col(p).setZero();
        }

      had_frozen_tail=tca_frozen_tail_ready;
      if (exact_audit) {
        if (mpi.rank==mpi.master && mpi.out)
          (*mpi.out) << "TCA exact "
                     << (tca_frozen_tail_ready ? "final audit" :
                                                   "tail calibration")
                     << ": evaluating full self-energy"
                     << (cmpPhysG ? " and physical functions" : "")
                     << std::endl;
        compact=evaluate_workspace(*diagram_dlr_workspace,cmpPhysG);
        const Eigen::MatrixXd new_tail=compact.sigma-active_compact.sigma;
        if (tca_frozen_tail_ready && mpi.rank==mpi.master && mpi.out) {
          const double full_scale=std::max(
            compact.sigma.norm(),std::numeric_limits<double>::min());
          const double tail_scale=std::max(
            tca_frozen_sigma_tail.norm(),
            std::numeric_limits<double>::min());
          const double relative_full=
            (new_tail-tca_frozen_sigma_tail).norm()/full_scale;
          const double relative_tail=
            (new_tail-tca_frozen_sigma_tail).norm()/tail_scale;
          (*mpi.out) << "TCA frozen-tail audit: change/full="
                     << relative_full << " change/old_tail="
                     << relative_tail << " tolerance="
                     << pars["tca_approx_audit_tolerance"].get<double>()
                     << (relative_full>
                           pars["tca_approx_audit_tolerance"].get<double>()
                           ? " WARNING" : " accepted") << std::endl;
        }
        tca_frozen_sigma_tail=new_tail;
        for (int topology=0;topology<2;++topology)
          tca_frozen_sigma_tail_by_topology[topology]=
            compact.sigma_by_topology[topology]-
            active_compact.sigma_by_topology[topology];
      } else if (freeze_tail) {
        compact=active_compact;
        compact.sigma+=tca_frozen_sigma_tail;
        for (int topology=0;topology<2;++topology)
          compact.sigma_by_topology[topology]+=
            tca_frozen_sigma_tail_by_topology[topology];
        using_frozen_tail=true;
      } else {
        compact=active_compact;
      }
    } else {
      compact=evaluate_workspace(*diagram_dlr_workspace,cmpPhysG);
    }
    const auto end=std::chrono::steady_clock::now();

    const auto& basis=diagram_dlr_workspace->basis();
    const auto to_dense=[&](const Eigen::MatrixXd& nodes) {
      const Eigen::MatrixXd coefficients=basis.values_to_coefficients(nodes);
      Eigen::MatrixXd dense(nodes.cols(),output_tau.extent(0));
      for (int it=0;it<output_tau.extent(0);++it)
        dense.col(it)=basis.evaluate(output_tau(it),coefficients).transpose();
      return dense;
    };
    const Eigen::MatrixXd sigma=to_dense(compact.sigma);
    dSm.resize(cix.Npseudo_unique,output_tau.extent(0));
    for (int p=0;p<cix.Npseudo_unique;++p)
      for (int it=0;it<output_tau.extent(0);++it)
        dSm(p,it)=sigma(p,it);

    // A third-order Phi ring contains six pseudo propagators.  Summing its
    // six functional cuts and closing once with G therefore counts Phi six
    // times.  The dense symmetric solver mesh is retained only for this final
    // inexpensive one-dimensional normalization integral.
    if (tau_dh.extent(0)!=output_tau.extent(0) ||
        Gm.extent(1)!=output_tau.extent(0))
      throw std::logic_error("TCA Phi quadrature mesh mismatch");
    double phi=0.0;
    std::array<double,2> phi_by_topology{{0.0,0.0}};
    std::array<Eigen::MatrixXd,2> sigma_by_topology{{
      to_dense(compact.sigma_by_topology[0]),
      to_dense(compact.sigma_by_topology[1])}};
    std::array<std::array<Eigen::MatrixXd,6>,2> sigma_cut_by_topology;
    if (store_self_diagnostics)
      for (int topology=0;topology<2;++topology)
        for (int cut=0;cut<6;++cut)
          sigma_cut_by_topology[topology][cut]=
            to_dense(compact.cuts_by_topology[topology][cut]);
    for (int p=0;p<cix.Npseudo_unique;++p)
      for (int it=0;it<output_tau.extent(0);++it)
        for (int topology=0;topology<2;++topology)
          phi_by_topology[topology]+=(1.0/6.0)*degeneracy[p]*tau_dh(it)*
            sigma_by_topology[topology](p,it)*
            Gm(p,output_tau.extent(0)-1-it);
    phi=phi_by_topology[0]+phi_by_topology[1];
    phi/=Q;
    for (double& value:phi_by_topology) value/=Q;

    // Sigma and Phi are related by closing a cut with the *current* G.  If a
    // weak Sigma tail is frozen, recomputing that closure would nevertheless
    // make its Phi contribution drift with G.  Freeze the reported Phi tail
    // explicitly as well.  This is a controlled numerical approximation, not
    // a claim that the frozen finite-step map remains exactly Phi-derivable;
    // the final full audit restores the exact diagrammatic result.
    if (have_active_compact) {
      std::array<double,2> active_phi{{0.0,0.0}};
      for (int topology=0;topology<2;++topology) {
        const Eigen::MatrixXd active_sigma=
          to_dense(active_compact.sigma_by_topology[topology]);
        for (int p=0;p<cix.Npseudo_unique;++p)
          for (int it=0;it<output_tau.extent(0);++it)
            active_phi[topology]+=(1.0/6.0)*degeneracy[p]*tau_dh(it)*
              active_sigma(p,it)*
              Gm(p,output_tau.extent(0)-1-it)/Q;
      }
      if (exact_audit && pars["tca_freeze_weak_tail"].get<bool>()) {
        std::array<double,2> new_phi_tail{{
          phi_by_topology[0]-active_phi[0],
          phi_by_topology[1]-active_phi[1]}};
        if (had_frozen_tail && mpi.rank==mpi.master && mpi.out) {
          const double old_tail_norm=std::hypot(
            tca_frozen_phi_tail_by_topology[0],
            tca_frozen_phi_tail_by_topology[1]);
          const double change=std::hypot(
            new_phi_tail[0]-tca_frozen_phi_tail_by_topology[0],
            new_phi_tail[1]-tca_frozen_phi_tail_by_topology[1]);
          (*mpi.out) << "TCA frozen-Phi audit: change/old_tail="
                     << change/std::max(
                          old_tail_norm,std::numeric_limits<double>::min())
                     << std::endl;
        }
        tca_frozen_phi_tail_by_topology=new_phi_tail;
        tca_frozen_tail_ready=true;
      } else if (using_frozen_tail) {
        for (int topology=0;topology<2;++topology)
          phi_by_topology[topology]=active_phi[topology]+
            tca_frozen_phi_tail_by_topology[topology];
        phi=phi_by_topology[0]+phi_by_topology[1];
      }
    }

    // The diagnostic files deliberately use the same dense tau mesh and
    // pseudo ordering as aStval.3.  They can therefore be compared directly
    // with a saved six-cut or reduced run without reconstructing DLR data.
    std::array<std::array<double,6>,2> phi_by_topology_cut{};
    if (store_self_diagnostics) {
      for (int p=0;p<cix.Npseudo_unique;++p)
        for (int it=0;it<output_tau.extent(0);++it)
          for (int topology=0;topology<2;++topology)
            for (int cut=0;cut<6;++cut)
              phi_by_topology_cut[topology][cut]+=
                (1.0/6.0)*degeneracy[p]*tau_dh(it)*
                sigma_cut_by_topology[topology][cut](p,it)*
                Gm(p,output_tau.extent(0)-1-it)/Q;

      if (mpi.rank==mpi.master) {
        const auto save_dense=[&](const std::string& filename,
                                  const Eigen::MatrixXd& dense) {
          Array<double,2> values(
            cix.Npseudo_unique,output_tau.extent(0));
          for (int p=0;p<cix.Npseudo_unique;++p)
            for (int it=0;it<output_tau.extent(0);++it)
              values(p,it)=dense(p,it);
          Print_Pseudo(filename,values,pars);
        };
        constexpr const char* topology_name[2]={"A","B"};
        for (int topology=0;topology<2;++topology) {
          save_dense("aStval.3.topology"s+topology_name[topology],
                     sigma_by_topology[topology]);
          Eigen::MatrixXd reconstructed=Eigen::MatrixXd::Zero(
            compact.sigma.rows(),compact.sigma.cols());
          for (int cut=0;cut<6;++cut) {
            reconstructed+=compact.cuts_by_topology[topology][cut];
            save_dense("aStval.3.topology"s+topology_name[topology]+
                         ".cut"+std::to_string(cut),
                       sigma_cut_by_topology[topology][cut]);
          }
          const double scale=std::max(
            1.e-300,compact.sigma_by_topology[topology].norm());
          (*mpi.out) << "TCA DLR self diagnostic topology="
                     << topology_name[topology]
                     << " cut_sum_relative="
                     << (reconstructed-
                         compact.sigma_by_topology[topology]).norm()/scale
                     << " Phi_by_cut=[";
          for (int cut=0;cut<6;++cut)
            (*mpi.out) << (cut ? "," : "")
                       << phi_by_topology_cut[topology][cut];
          (*mpi.out) << "] Phi_cut_sum="
                     << std::accumulate(
                          phi_by_topology_cut[topology].begin(),
                          phi_by_topology_cut[topology].end(),0.0)
                     << std::endl;
        }
      }
    }

    G00t.resize(cix.Nbath_unique,output_tau.extent(0)); G00t=0.0;
    F00t.resize(cix.Nbath_unique,output_tau.extent(0)); F00t=0.0;
    F01t.resize(cix.Nbath_unique,output_tau.extent(0)); F01t=0.0;
    if (cmpPhysG) {
      const Eigen::MatrixXd g=to_dense(compact.g);
      const Eigen::MatrixXd f00=to_dense(compact.f00);
      const Eigen::MatrixXd f01=to_dense(compact.f01);
      Array<double,1> bath_degeneracy(cix.Nbath_unique_unique);
      bath_degeneracy=0.0;
      for (int orbital=0;orbital<cix.Nbath;++orbital)
        bath_degeneracy(std::abs(cix.bath_sigind0(orbital)))+=1.0;
      Array<int,1> sigind(cix.Nbath_unique);
      Array<bool,1> superconducting(cix.Nbath_unique);
      superconducting=false;
      for (int orbital=0;orbital<cix.Nbath;++orbital) {
        const int ib=cix.bath_sigind(orbital);
        sigind(ib)=cix.bath_sigind0(orbital);
        superconducting(ib)=cix.SC_bath(orbital);
      }
      for (int ib=0;ib<cix.Nbath_unique;++ib) {
        const int s=std::abs(sigind(ib));
        const double sign=sigind(ib)<0 ? -1.0 : 1.0;
        const double normalizer=sign/(bath_degeneracy(s)*Q);
        const double bulla=superconducting(ib) ? 2.0 : 1.0;
        for (int it=0;it<output_tau.extent(0);++it) {
          G00t(ib,it)=normalizer*g(s,it);
          F00t(ib,it)=normalizer*bulla*f00(s,it);
          F01t(ib,it)=normalizer*bulla*f01(s,it);
        }
      }
    }

    double seconds=std::chrono::duration<double>(end-begin).count();
#ifdef _MPI
    if (mpi.size>1) {
      double maximum=0.0;
      MPI_Allreduce(&seconds,&maximum,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
      seconds=maximum;
    }
#endif
    if (mpi.rank==mpi.master)
      (*mpi.out) << "TCA DLR: " << (root_sum ? "self_groups=" : "plans=")
                 << compact.candidate_plans
                 << " nonempty_topologies=["
                 << compact.nonempty_plans_by_topology[0] << ","
                 << compact.nonempty_plans_by_topology[1] << "]"
                 << " physical=" << (cmpPhysG ? "yes" : "no")
                 << " mpi=" << (pars["tca_dlr_mpi"].get<bool>() ? "on" : "off")
                 << " self_mode="
                 << (!use_active ? "exact" :
                     exact_audit ? "exact-audit" :
                     using_frozen_tail ? "active+frozen-tail" :
                                         "active-only")
                 << " active_pseudos="
                 << (use_active
                       ? std::count(tca_active_pseudo_mask.begin(),
                                    tca_active_pseudo_mask.end(),
                                    static_cast<unsigned char>(1))
                       : cix.Npseudo_unique)
                 << "/" << cix.Npseudo_unique
                 << " poles="
                 << (use_active ? tca_active_dlr_workspace->pole_count() :
                                  diagram_dlr_workspace->pole_count())
                 << "/" << diagram_dlr_workspace->pole_count()
                 << " evaluate_seconds=" << seconds
                 << " Phi_by_topology=[" << phi_by_topology[0] << ","
                 << phi_by_topology[1] << "]"
                 << " Phi=" << phi << std::endl;
    return phi;
  }

  void SplineG(const Array<double,1>& pptau, Array<double,2>& Gm, Array<double,2>& gm)
  {
    Gms.resize(TinyVector<int,2>(pptau.extent(0),Gm.extent(0)));
    for (int l=0; l<Gm.extent(0); ++l)
      for (int it=0; it<pptau.extent(0); ++it) Gms(it,l) = Gm(l,it);
    Gms.splineIt(pptau);
    gms.resize(TinyVector<int,2>(pptau.extent(0),gm.extent(0)));
    for (int l=0; l<gm.extent(0); ++l)
      for (int it=0; it<pptau.extent(0); ++it) gms(it,l) = gm(l,it);
    gms.splineIt(pptau);
    /*
    ofstream dout1("SG_debug1.dat");
    for (int it=0; it<pptau.extent(0); ++it){
      dout1 << pptau(it) << " ";
      for (int a=0; a<cix.Npseudo_unique; ++a) dout1 << Gms(it,a) << " ";
      dout1 << endl;
    }
    cout<<"Gms splined!"<<endl;
    */
  }
  //  void ImportanceOfPseudo(vector<int>& G_importance,
  //			  const Array<double,1>& tau, const Array<double,2>& Gm, std::array<double,2> cut={0.1,0.5}){};
  
  void StopAdaptingMesh(bool ){};
  void StartAdaptingMesh() {};
  void Update_Sm(Array<double,2>& Sm, const Array<double,2>& Smt, const Array<double,1>& tau, const Cix& cix, MPIcls& mpi, bool check_causality=false) {
    Sm += Smt;
  }
  double Difference(const Array<double,2>& Sm_new, const Array<double,2>& Sm, const Array<double,1>& Prob){
    //return sum(abs(Sm_new-Sm))/Sm_new.size();
    double diff=0;
    double nrm = sum(abs(Prob));
    int ntau = Sm_new.extent(1);
    for (int ii=0; ii<cix.Npseudo_unique; ++ii)
      diff += sum(abs(Sm_new(ii,Range::all())-Sm(ii,Range::all())))/ntau*abs(Prob(ii))/nrm;
    return diff;
  }
  void PhysicalSelfEnergy(Array<complex<double>,2>& Sigm, Array<complex<double>,2>& Gc, Array<double,1>& density_matrix, double& Ekin, double& nf,
			  const Array<double,2>& G00, const Array<double,2>& F00, const Array<double,2>& F01,
			  const Cix& cix, const BathBlocks& bb, const json& pars, int ist, MPIcls& mpi, string dysonFname)
  {
    bool debug = false;
    Array<int,1> sigind0(cix.Nbath_unique), deg_sigind(cix.Nbath_unique_unique); sigind0=0; deg_sigind=0;
    for (int ib=0; ib<cix.Nbath; ++ib) sigind0(cix.bath_sigind(ib)) = cix.bath_sigind0(ib);
    for (int ib=0; ib<cix.Nbath_unique; ++ib) deg_sigind(abs(sigind0(ib)))+=1;
    // First we combine components that are opposite in value, if any.
    Array<double,2> Gimp(cix.Nbath_unique_unique,G00.extent(1)); Gimp=0;
    Array<double,2> Fimp1(cix.Nbath_unique_unique,F00.extent(1)); Fimp1=0;
    Array<double,2> Fimp2(cix.Nbath_unique_unique,F01.extent(1)); Fimp2=0;
    for (int ib=0; ib<cix.Nbath_unique; ++ib){
      double ib_sg = (sigind0(ib) < 0 ? -1 : 1);
      int iib   = std::abs(sigind0(ib));
      Gimp(iib,Range::all())  += G00(ib,Range::all())*(ib_sg/deg_sigind(iib));
      Fimp1(iib,Range::all()) += F00(ib,Range::all())*(ib_sg/deg_sigind(iib));
      Fimp2(iib,Range::all()) += F01(ib,Range::all())*(ib_sg/deg_sigind(iib));
    }
    Gc.resize(cix.Nbath_unique_unique,iom.extent(0));
    Array<complex<double>,2> Fc1(cix.Nbath_unique_unique,iom.extent(0)), Fc2(cix.Nbath_unique_unique,iom.extent(0));
    {
      /*
      Spline1D<double> gf(tau.extent(0)), ff1(tau.extent(0)), ff2(tau.extent(0));
      for (int ib=0; ib<cix.Nbath_unique_unique; ++ib){
	for (int i=0; i<tau.size(); ++i){
	  gf[i]=-Gimp(ib,i);  // We choose the positive sign for G(tau)
	  ff1[i]=-Fimp1(ib,i);
	  ff2[i]=-Fimp2(ib,i);
	}
	gf.splineIt(tau);
	ff1.splineIt(tau);
	ff2.splineIt(tau);
	for (int im=0; im<iom.extent(0); ++im){
	  Gc(ib,im) = gf.Fourier(iom(im), tau);
	  Fc1(ib,im) = ff1.Fourier(iom(im), tau);
	  Fc2(ib,im) = ff2.Fourier(iom(im), tau);
	}
      }
      */
      Spline1DMany gf(tau.extent(0),cix.Nbath_unique_unique*3);
      for (int i=0; i<tau.size(); ++i){
	for (int ib=0; ib<cix.Nbath_unique_unique; ++ib){
	  gf(i,3*ib)   = -Gimp(ib,i);  // We choose the positive sign for G(tau)
	  gf(i,3*ib+1) = -Fimp1(ib,i); // and F1(tau) and F2(tau)
	  gf(i,3*ib+2) = -Fimp2(ib,i);
	}
      }
      gf.splineIt(tau);
      Array<complex<double>,1> fiom(cix.Nbath_unique_unique*3);
      for (int im=0; im<iom.extent(0); ++im){
	gf.Fourier(fiom, iom(im));
	for (int ib=0; ib<cix.Nbath_unique_unique; ++ib){
	  Gc(ib,im)  = fiom(3*ib);
	  Fc1(ib,im) = fiom(3*ib+1);
	  Fc2(ib,im) = fiom(3*ib+2);
	}
      }
    }
    Sigm.resize(cix.Nbath_unique_unique,iom.extent(0)); Sigm=0;
    Array<complex<double>,2> SigB(cix.Nbath_unique_unique,iom.extent(0)); SigB=0;
    Array<complex<double>,2> SigD(cix.Nbath_unique_unique,iom.extent(0)); SigD=0;
    Array<complex<double>,2> Gcn(cix.Nbath_unique_unique,iom.extent(0)); Gcn=0;
    density_matrix.resize(cix.Nbath_unique_unique); density_matrix=0;
    Ekin = 0; nf=0;

    const bool blend_sdyson_bulla = pars["Sdyson_Bulla_blend"].get<bool>();
    const double blend_alpha = pars["Sdyson_Bulla_alpha"].get<double>();
    const double blend_cutoff_fraction = pars["Sdyson_Bulla_cutoff_fraction"].get<double>();
    const double blend_start_ratio = pars["Sdyson_Bulla_start_ratio"].get<double>();
    const double r_OCA_mesh = pars["r_OCA_mesh"].get<double>();
    const double Ntau = pars["Ntau"].get<int>();
    const double omega_N = M_PI*Ntau*r_OCA_mesh/beta;
    const double w_end = blend_alpha*blend_cutoff_fraction*omega_N;
    const double w_start = blend_start_ratio*w_end;
    auto smoothstep = [](double x){
      x = std::max(0.0, std::min(1.0, x));
      return x*x*x*(10.0 - 15.0*x + 6.0*x*x);
    };
    auto bulla_weight = [&](double w){
      if (!blend_sdyson_bulla) return 1.0;
      if (blend_alpha <= 0.0) return 0.0;
      if (w_end <= w_start) return 1.0;
      return smoothstep((std::abs(w)-w_start)/(w_end-w_start));
    };
    
    vector<bool> component_set(cix.Nbath_unique_unique,false);
    complex<double> j(0.,1.);
    Array<double,1> _Ekin(cix.Nbath_unique_unique), _nf(cix.Nbath_unique_unique);
    _Ekin=0; _nf=0;
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
	      int ilast = iom.extent(0)-1;
	      double E_inf_scalar;
	      complex<double> om_Eimp, tail_sum=0;
	      double _Ekin_ = 0, _nf_ = 0;
	      //SigB(iib,Range::all())  = 0.5*(Fc1(iib,Range::all())+Fc2(iib,Range::all()))/Gc(iib,Range::all());
	      SigB(iib,Range::all())  = Fc1(iib,Range::all())/Gc(iib,Range::all());
	      E_inf_scalar = Eimp*ib_sg + std::real(SigB(iib,ilast))*ib_sg;
	      SigD(iib,Range::all()) = -1.0/Gc(iib,Range::all())*ib_sg-iDeltaw(iib,Range::all())*ib_sg-Eimp*ib_sg+iom(Range::all())*j;
	      for (int iw=0; iw<iom.size(); ++iw){
		const double W = bulla_weight(iom(iw));
		Sigm(iib,iw) = (1.0-W)*SigD(iib,iw) + W*SigB(iib,iw);
	      }
	      Gcn (iib,Range::all()) =  1.0/(iom(Range::all())*j-Eimp*ib_sg-Sigm(iib,Range::all())*ib_sg-iDeltaw(iib,Range::all())*ib_sg);
	      for (int iw=0; iw<iom.size(); ++iw){
		om_Eimp = iom(iw)*j - E_inf_scalar;
		// w_iom(iw) replaces the naive sparse-mesh sum by the exact sum over
		// all integer Matsubara points, while remaining equivalent to the
		// spline reconstruction defined above.
		tail_sum += w_iom(iw) * (Gcn(iib,iw) - 1.0/om_Eimp);
		_Ekin_ += w_iom(iw) * std::real(Gcn(iib,iw)*iDeltaw(iib,iw)*ib_sg);
		_nf_ += w_iom(iw) * std::real(Gcn(iib,iw) - 1.0/om_Eimp);
	      }
	      _Ekin(iib) = _Ekin_*2.0/beta;
	      _nf(iib) = _nf_*2.0/beta + fermi(beta*E_inf_scalar);
	      density_matrix(iib) = std::real(tail_sum)*2.0/beta + fermi(beta*E_inf_scalar);
	      Ekin += _Ekin(iib);
	      nf += _nf(iib);
	      component_set[iib]=true;
	    }else{
	      Ekin += _Ekin(iib);
	      nf += _nf(iib);
	    }
	  }
	}
      }else{
	bool SC = bb.is_block_sc[ib];
	MatrixXcd Gm(N, N), Fm(N, N), Dm(N, N), om_Eimp(N, N), _density_matrix_(N,N); // Eigen matrices for linear algebra
	_density_matrix_.setZero();
	double _Ekin_=0, _nf_=0;
	
	ofstream out1, out2, out3, out4, out5, out6, out7;
	if (mpi.rank==mpi.master){
	  out1.open("Sm."+to_string(ib)+".dat");
	  out2.open("Gmc."+to_string(ib)+".dat");
	  out3.open("Dm."+to_string(ib)+".dat");
	  out4.open("Gm."+to_string(ib)+".dat");
	  out5.open("Fm."+to_string(ib)+".dat");
	  out6.open("Sms."+to_string(ib)+".dat");
	  out7.open("Smd."+to_string(ib)+".dat");
	}

	Array<double,2> E_inf(N,N); E_inf=0;
	{// next few lines determines Sigma_\infty and high frequency limit of G(iom) as 1/(iom-E_{infty})
	  int ilast = iom.extent(0)-1;
	  MatrixXcd Gm_inf(N, N), Fm_inf(N, N);
	  for (int iorb=0; iorb<cix.Nbath; ++iorb){
	    if (bb.orb_ind[iorb][0]==ib){
	      int j=bb.orb_ind[iorb][1], k=bb.orb_ind[iorb][2];
	      int iib0 = cix.bath_sigind0(iorb);
	      double ib_sg = (iib0 < 0 ? -1 : 1);
	      int iib   = std::abs(iib0);
	      if (SC && j>=N/2 && k>=N/2){ // This is hole-hole (or 22) part of Nambu matrix
		// SC/Nambu “conjugated” sector, namely, G_{22}(i*om) = -G^T_{11}(-i*om)
		Gm_inf(k,j) = -ib_sg*std::conj(Gc(iib,ilast));
	      }else{
		// “normal” state or “normal” sector
		Gm_inf(j,k) = ib_sg*Gc(iib,ilast);
	      }
	      if (SC){
		if (k>=N/2){
		  if (j>=N/2) // G_22 has an extra minus sign
		    Fm_inf(k,j) = ib_sg*std::conj(Fc2(iib,ilast));
		  else        // G_{12} has also an extra minus sign
		    Fm_inf(j,k) = -ib_sg*Fc2(iib,ilast);
		}else{
		  Fm_inf(j,k) = ib_sg*Fc1(iib,ilast);
		}
	      }else{
		Fm_inf(j,k) = ib_sg*Fc1(iib,ilast); // for normal state only Fc1 is used
	      }
	    }
	  }
	  MatrixXcd Gminv = Gm_inf.inverse();     // Gm^{-1}
	  MatrixXcd Sm_inf= Gminv * Fm_inf;       // Sm = Gm^{-1} * Fm
	  if (SC){ // we know that anomalous self-energy should vanish in infinity.
	    const int N2 = N/2;
	    Sm_inf.block(0,  N2, N2, N2).setZero();
	    Sm_inf.block(N2, 0,  N2, N2).setZero();
	  }
	  for (int iorb=0; iorb<cix.Nbath; ++iorb){
	    if (bb.orb_ind[iorb][0]==ib){
	      int j=bb.orb_ind[iorb][1], k=bb.orb_ind[iorb][2];
	      int iib0 = cix.bath_sigind0(iorb);
	      double ib_sg = (iib0 < 0 ? -1 : 1);
	      int iib   = std::abs(iib0);
	      if (SC && j>=N/2 && k>=N/2){ // This is hole-hole (or 22) part of Nambu matrix
		// SC/Nambu “conjugated” sector, namely, G_{22}(i*om) = -G^T_{11}(-i*om)
		E_inf(k,j) += -ib_sg*cix.Ebath(iib);
	      }else{
		// “normal” state or “normal” sector
		E_inf(j,k) += ib_sg*cix.Ebath(iib);
	      }
	      E_inf(j,k) += std::real(Sm_inf(j,k));
	    }
	  }
	}
	
	
	for (int iw=0; iw<iom.size(); ++iw){
	  Gm.setZero();
	  Fm.setZero();
	  Dm.setZero();
	  om_Eimp.setZero();
	  // Diagonal: ω * I
	  for (int i=0; i<N; ++i) om_Eimp(i,i) = iom(iw)*j;
	  // Fill Gm, Fm, Dm, om_Eimp from bath data
	  for (int iorb=0; iorb<cix.Nbath; ++iorb){
	    if (bb.orb_ind[iorb][0]==ib){
	      int j=bb.orb_ind[iorb][1], k=bb.orb_ind[iorb][2];
	      int iib0 = cix.bath_sigind0(iorb);
	      double ib_sg = (iib0 < 0 ? -1 : 1);
	      int iib   = std::abs(iib0);
	      if (SC && j>=N/2 && k>=N/2){ // This is hole-hole (or 22) part of Nambu matrix
		if (iib0>=0){
		  // SC/Nambu “conjugated” sector, namely, G_{22}(i*om) = -G^T_{11}(-i*om)
		  Gm(k,j) = -std::conj(Gc(iib,iw));
		  Dm(k,j) = -std::conj(iDeltaw(iib,iw));
		  om_Eimp(k,j) += cix.Ebath(iib);
		}else{
		  // SC/Nambu “conjugated” sector, namely, G_{22}(i*om) = -G^T_{11}(-i*om)
		  // Here bath has -sign as compared to the stored bath index, hence one more minus sign. For example, d-wave SC for y-direction
		  Gm(k,j) = std::conj(Gc(iib,iw));
		  Dm(k,j) = std::conj(iDeltaw(iib,iw));
		  om_Eimp(k,j) -= cix.Ebath(iib);
		}
	      }else{
		if (iib0>=0){
		  // “normal” state or “normal” sector
		  Gm(j,k) = Gc(iib,iw);
		  Dm(j,k) = iDeltaw(iib,iw);
		  om_Eimp(j,k) -= cix.Ebath(iib);
		}else{
		  // “normal” state or “normal” sector
		  Gm(j,k) = -Gc(iib,iw);
		  Dm(j,k) = -iDeltaw(iib,iw);
		  om_Eimp(j,k) += cix.Ebath(iib);
		}
	      }
	      if (SC){
		if (k>=N/2){
		  if (j>=N/2) // G_22 has an extra minus sign
		    Fm(k,j) = ib_sg*std::conj(Fc2(iib,iw));
		  else        // G_{12} has also an extra minus sign
		    Fm(j,k) = -ib_sg*Fc2(iib,iw);
		}else{
		  Fm(j,k) = ib_sg*Fc1(iib,iw);
		}
	      }else{
		Fm(j,k) = ib_sg*Fc1(iib,iw); // for normal state only Fc1 is used
	      }
	    }
	  }
	  MatrixXcd Gminv = Gm.inverse();         // Gm^{-1}
	  MatrixXcd SmB   = Gminv * Fm;           // Bulla/EOM estimator: SmB = Gm^{-1} * Fm
	  MatrixXcd SmD   = om_Eimp - Gminv - Dm; // Dyson: Σ = (ω - E_imp) - G^{-1} - Δ
	  
	  if (mpi.rank==mpi.master){
	    out1<<iom(iw)<<" ";
	    out3<<iom(iw)<<" ";
	    out4<<iom(iw)<<" ";
	    out5<<iom(iw)<<" ";
	    out7<<iom(iw)<<" ";
	    for (int i=0; i<N; ++i)
	      for (int j=0; j<N; ++j){
		out1 << SmB(i,j).real() << " " << SmB(i,j).imag()<< " ";
		out3 << Dm(i,j).real() << " " << Dm(i,j).imag()<< " ";
		out4 << Gm(i,j).real() << " " << Gm(i,j).imag()<< " ";
		out5 << Fm(i,j).real() << " " << Fm(i,j).imag()<< " ";
		out7 << SmD(i,j).real() << " " << SmD(i,j).imag()<< " ";
	      }
	    out1<<endl;
	    out3<<endl;
	    out4<<endl;
	    out5<<endl;
	    out7<<endl;
	  }
	  
	  if (SC){ // we need to symmetrize anomalous components
	    const int N2 = N/2;
	    MatrixXcd S12 = SmB.block(0, N2, N2, N2);
	    MatrixXcd S21 = SmB.block(N2, 0, N2, N2);
	    MatrixXcd A = 0.5*(S12 + S21);
	    SmB.block(0, N2, N2, N2) = A;
	    SmB.block(N2, 0, N2, N2) = A;
	  }
	  const double W = bulla_weight(iom(iw));
	  MatrixXcd Sm = (1.0-W)*SmD + W*SmB;     // NEW after realizing Sdyson is needed at low w: Combined self-energy used by DMFT: Sdyson at low w, Bulla at high w.
	  MatrixXcd Ginv  = om_Eimp - Sm - Dm;    // G^{-1} = (ω - E_imp - Σ - Δ)
	  MatrixXcd Gmc  = Ginv.inverse();        // G = (ω - E_imp - Σ - Δ)^{-1}
	  
	  for (int a=0; a<N; ++a){
	    // density_matrix(a,b) = \sum_{iom} G_{a,b}(iom) - delta(a,b)/(iom-Eimp)
	    for (int b=0; b<N; ++b){
	      _density_matrix_(a,b) += w_iom(iw) * Gmc(a,b);
	      _Ekin_ += w_iom(iw) * std::real(Gmc(a,b)*Dm(b,a));
	    }
	    _density_matrix_(a,a) -= w_iom(iw) * 1.0/(iom(iw)*j-E_inf(a,a));
	    _nf_ += w_iom(iw) * std::real( Gmc(a,a) - 1.0/(iom(iw)*j-E_inf(a,a)) );
	  }
	  
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
		SigB(iib,iw) = SmB(j,k)*ib_sg;
		SigD(iib,iw) = SmD(j,k)*ib_sg;
		Gcn(iib,iw) = Gmc(j,k)*ib_sg;
		//cout<<"Sigm["<<iib<<","<<iw<<"]= Sm["<<j<<","<<k<<"]*"<<ib_sg<<endl;
		component_set[iib] = true;
		// There might be a bit of a problem when G_{dn} is different than G_{up} in SC state, i.e.
		// both SC and magnetism are present. In this case we would need to write into iw that
		// correspond to -omega, which is not readily available.
		// We will take care of that specific case when needed.
	      }
	    }
	  }
	  
	  if (mpi.rank==mpi.master){
	    out2<<iom(iw)<<" ";
	    out6<<iom(iw)<<" ";
	    for (int i=0; i<N; ++i)
	      for (int j=0; j<N; ++j){
		out2 << Gmc(i,j).real()<< " " << Gmc(i,j).imag()<< " ";
		out6 << Sm(i,j).real() << " " << Sm(i,j).imag()<< " ";
	      }
	    out2<<endl;
	    out6<<endl;	    
	  }
	}
	
	_density_matrix_ *= 2.0/beta; // 2 because sum from -infty->infty
	_Ekin_ *= 2.0/beta;
	_nf_ *= 2.0/beta;
	Array<int,1> deg_dm(cix.Nbath_unique_unique); deg_dm=0;
	for (int iorb=0; iorb<cix.Nbath; ++iorb){
	  if (bb.orb_ind[iorb][0]==ib){
	    int j=bb.orb_ind[iorb][1], k=bb.orb_ind[iorb][2];
	    int iib0 = cix.bath_sigind0(iorb);
	    double ib_sg = (iib0 < 0 ? -1 : 1);
	    int iib   = std::abs(iib0);
	    // density_matrix(iib) is a reduced quantity: several physical
	    // orbitals can map to the same reduced bath index iib.  To keep the
	    // same convention as density_matrix0 from DeterminePseudos, we
	    // average over all such equivalent representatives.
	    density_matrix(iib) += std::real(_density_matrix_(j,k))*ib_sg;
	    deg_dm(iib) += 1;
	    if (k==j){
	      // The Fermi-function add-back is different for density_matrix and
	      // nf on purpose:
	      //
	      // 1) density_matrix(iib) is an averaged reduced-component object,
	      //    hence we add one representative Fermi term only *after*
	      //    averaging over equivalent orbitals.
	      //
	      // 2) nf is the total physical occupancy, hence every physical
	      //    orbital contributes and we must sum the Fermi term over all
	      //    diagonal orbital representatives.
	      //_nf_ += fermi(beta*cix.Ebath(iib));
	      _nf_ += fermi(beta*E_inf(k,k));
	      density_matrix(iib) += fermi(beta*E_inf(k,k));
	    }
	  }
	}
	for (int iib=0; iib<cix.Nbath_unique_unique; ++iib){
	  if (deg_dm(iib)>0){
	    density_matrix(iib) /= deg_dm(iib);
	    //density_matrix(iib) += fermi(beta*cix.Ebath(iib));
	  }
	}
	Ekin += _Ekin_;
	nf += _nf_;
      }
    }

    if (mpi.rank==mpi.master){
      Print_Physical1("Gimpw.out", Gc);
    }
    Gc = Gcn; // Be careful. Gc becomes 1/(om-Eimp-Sigma-Delta)
    if (mpi.rank==mpi.master){
      Print_Physical1("Sbulla.out", SigB);
      Print_Physical1(dysonFname+"out", SigD);
    }
    /*
    if (mpi.rank==mpi.master){
      ofstream out4("Sgm.debug"); out4.precision(12);
      for (int iw=0; iw<iom.extent(0); ++iw){
	out4<<iom(iw)<<" ";
	for (int ip=0; ip<Sigm.extent(0); ++ip) out4 << Sigm(ip,iw).real()<<" "<< Sigm(ip,iw).imag()<<" "; 
	out4 << endl;
      }
      ofstream out5("Ggm.debug"); out5.precision(12);
      for (int iw=0; iw<iom.extent(0); ++iw){
	out5<<iom(iw)<<" ";
	for (int ip=0; ip<Sigm.extent(0); ++ip) out5 << Gc(ip,iw).real()<<" "<< Gc(ip,iw).imag()<<" "; 
	out5 << endl;
      }
    }
    */
  }
  Array<double,2>& imagt(Array<double,2>& Gt){return Gt;}
  //
  void Print_Pseudo(const string& filename, const Array<double,2>& Smt_new, json& pars){
    ofstream dout(filename); dout.precision(12);
    dout<<"# beta="<<beta<<" U="<<pars["U"]<<" Norder="<<pars["Norder"]<<endl;
    for (int i=0; i<tau.size(); ++i){
      dout<<tau(i)<<" ";
      for (int ip=0; ip<Smt_new.extent(0); ++ip) Prints(dout,Smt_new(ip,i));
      dout<<endl;
    }
  }
  void Print_Physical(const string& filename, const Array<double,2>& G00, const Array<double,2>& F00, const Array<double,2>& F01)
  {
    ofstream sout(filename); sout.precision(12);
    for (int i=0; i<tau.size(); ++i){
      sout<<tau(i)<<" ";
      for (int ip=0; ip<G00.extent(0); ++ip) Prints(sout,G00(ip,i)); 
      for (int ip=0; ip<F00.extent(0); ++ip) Prints(sout,F00(ip,i));
      for (int ip=0; ip<F01.extent(0); ++ip) Prints(sout,F01(ip,i));
      sout<<endl;
    }
  }
};


class ImagAxis1new{
public:
  const bool Phi_diagonal=true;
  const double& beta;
  const Array<double,1>& dltau;
  const Array<double,1>& pptau;
  //const vector<Spline1D<double>>& Deltas;
  const Spline1DMany& Deltas;
  const SVDFunc& svdf;
  const Cix& cix;
  const std::vector<int>& G_importance;
  const int& OCA_cut_sum;
  //int matrix_mc_samples, matrix_mc_burnin; // added by codex
  int N_min_order, N_max_order;
  Diag2Matrix tdiag2matrix, cdiag2matrix;
  const Spline1DMany& Gms, gms;
  int wNbin_lmax;
  Array<double,2> cGm;   // current Gm(max_pp_dim,2*N_max_order) 
  Array<double,2> tGm;   //   trial Gm(max_pp_dim,2*N_max_order)
  vector<unordered_map<int,double>> tGmi; // similar to tGm, but in case we compute Phi as diagonal values.
  vector<unordered_map<int,double>> cGmi; // similar to tGm, but in case we compute Phi as diagonal values.
  Array<double,1> cDelta;// current Delta(N_max_order)
  Array<double,1> tDelta;//   trial Delta(N_max_order)
  Array<double,1> one_over_Phi;
  Array<double,1> cPhi, tPhi;
  Array<int,1> cops, tops;
  Array<int,2> ppstates, tpstates;
  //
  Array<Array<double,2>,2> Lft, Rgh;
  // only one of those two is needed.....
  Array<Array<double,2>,2> Gmx;
  Array<Array<int,2>,2> remember_ii;
  vector<Array<double,2>> Gx;
  //
  vector<Array<double,1>> Fi;
  bool topology_change;
  Array<bool,1> times_changed, times_changed0;
  TinyVector<int,2> _tind_;
  Array<double,2> dS, _dS_;
  Array<intpar,1> icps; // used for interpolate t_{i+1}-t_i.
  Array<double,2> utli, utl;
  int sgns;
  bool F_use_ie, cmpPhysG;
  Array<Array<double,2>,2> Gtm;
  Array<Array<double,2>,2> Gtp;
  Array<Array<int,2>,2> importance_m;
  Array<Array<int,2>,2> importance_p;
#ifdef _TIMEIT_  
  double t_intrp, t_allowm, t_Gx, t_Lft, t_Rgh, t_stval, t_g00, t_f00, t_dS;
  double t_sdot, t_gintrp, t_gcpy, t_fdot;
  double t_cOCA1, t_cOCA2;
#endif
  //
  //Array<double,2> _F00_OCA, _F01_OCA, _F02_OCA, _F03_OCA;
  //
  //
  ImagAxis1new(const ImagAxis0& ra0, const Array<double,1>& _pptau_, const Array<double,2>& _Gm_,
	       const Cix& _cix_, json& pars, bool _cmpPhysG_, MPIcls& mpi) :
    beta(ra0.beta),dltau(ra0.tau), pptau(_pptau_), Deltas(ra0.Deltas),svdf(ra0.svdf), cix(_cix_), Gms(ra0.Gms), gms(ra0.gms),
    G_importance(ra0.G_importance), OCA_cut_sum(ra0.OCA_cut_sum),
    //matrix_mc_samples(pars["matrix_mc_samples"]), matrix_mc_burnin(pars["matrix_mc_burnin"]),
    N_min_order(pars["N_min_order"]), N_max_order(pars["Norder"]), tdiag2matrix(N_max_order,_cix_), cdiag2matrix(N_max_order,_cix_), cmpPhysG(_cmpPhysG_)
  {
    wNbin_lmax = svdf.lmax; // how many components of physical quantities we need to store
    topology_change=true;
    times_changed.resize(2*N_max_order);
    times_changed=false;
    times_changed0.resize(2*N_max_order);
    times_changed0=false;
#ifdef _TIMEIT_
    t_intrp=t_allowm=t_Gx=t_Lft=t_Rgh=t_stval=0;
    t_g00=t_f00=t_dS=t_sdot=t_gintrp=t_gcpy=t_fdot=0;
    t_cOCA1=0, t_cOCA2=0;
#endif    
  };
  void PrintTimings(ostream& log)
  {
#ifdef _TIMEIT_    
    fmt::print(log,"t_intrp={:.4f}s t_allowm={:.4f}s t_Gx={:.4f}s t_Lft={:.4f}s t_Rgh={:.4f}s t_stval={:.4f}s t_f00={:.4f}s t_g00={:.4f}s t_dS={:.4f}\n",
	       t_intrp , t_allowm , t_Gx , t_Lft , t_Rgh , t_stval , t_f00 , t_g00 , t_dS );
    fmt::print(log,"t_fdot={:.4f} t_sdot={:.4f} t_gcpy={:.4f} t_gintrp={:.4f} t_cOCA={:.4f}+{:.4f}\n", t_fdot , t_sdot , t_gcpy , t_gintrp, t_cOCA1, t_cOCA2 );
#endif    
  }
  int number_of(int Norder){return 2*Norder;}
  int number_of2(int Norder){return 2*Norder+1;}
  bool random_tw(Array<double,1>& times, int ind, double xi, int Norder, int itt)
  {
    // This function modifies the time array at index ind.
    // We called this Next_tarray in Python
    int size_2n = times.extent(0) - 1;  // Equivalent to len(times) - 1 in Python
    if (ind > 0){
      times(ind) = times(ind-1)*(1-xi) + xi*times(ind+1);
    } else if (ind == 0) {
      double t0 = xi * times(1) + (1-xi)*(times(size_2n-1)-beta);
      for (int i = 1; i<size_2n; ++i) times(i) -= t0;
    }
    return true;
  }
  /*
  void SplineG(Array<double,2>& Gm, Array<double,2>& gm)
  {
    Gms.resize(Gm.extent(0));
    for (int l=0; l<Gm.extent(0); ++l){
      //Gms[l].resize(pptau.size());
      //for (int it=0; it<pptau.size(); ++it) Gms[l][it] = Gm(l,it);
      Array<double,1> row;
      row.reference(Gm(l, blitz::Range::all()));
      Gms[l].Init(row);
      Gms[l].splineIt(pptau);
    }
    gms.resize(gm.extent(0));
    for (int l=0; l<gm.extent(0); ++l){
      Array<double,1> row;
      row.reference(gm(l, blitz::Range::all()));
      gms[l].Init(row);
      gms[l].splineIt(pptau);
    }
  } 
  */ 
  double recmp_Q(const Array<double,2>& Gm, const Cix& cix)
  {
    int Nt = pptau.size();
    double Q=0;
    for (int ib=0; ib<cix.Nblock; ++ib){
      for (int i=0; i<cix.Block[ib].size(); ++i){
	int _ii_ = cix.Sgm_symmetry[ib](i,i);
	if (_ii_!=0) Q += Gm(abs(_ii_)-1,Nt-1);
      }
    }
    return Q;
  }
  void Create_Fermi_Splines(const Array<double,2>& Gm)
  {}
  void InitializeStorageMC(){//int Npseudo_unique, int Nbath_unique, int max_pp_dim=2){
    cDelta.resize(N_max_order); cDelta=0;
    tDelta.resize(N_max_order); tDelta=0;
  }
  void DeleteStorageMC(){
    cDelta.resize(0);
    tDelta.resize(0);
  }
  //void Zero_time_freq(Array<double,1>& times1, int Norder)
  //{
  //  times1 = 0;
  //  times1(2*Norder) = beta;
  //}
  template<class RNG>
  void SetRandom_time_freq(Array<double,1>& times1, int Norder, RNG& drand)
  {
    times1(0)=0; times1(2*Norder)=beta;
    for (int it=1; it<2*Norder; ++it) times1(it) = drand()*beta;
    sort(times1.data(),times1.data()+2*Norder);
  }
  bool ReadStatusFile(Array<int,2>& cdiag, Array<int,1>& cdirspin, Array<double,1>& ctimes,
		      const string& status_file, ostream& log)
  {
    Array<int,2> ndiag;
    Array<int,1> ndirspin;
    Array<double,1> ntimes;
    
    ifstream inp(status_file);
    int Norder=0;
    inp>>Norder;
    if (not inp or inp.eof()) return false;
    if (Norder>N_max_order or Norder<N_min_order){
      log<<"Read diagram is of too high/low order not("<<N_min_order<<"<="<<Norder<<"<="<<N_max_order<<") hence ignoring file "<<status_file<<endl;
      return false;
    }
    ndiag.resize(Norder,2);
    for (int i=0; i<Norder; ++i) inp >> ndiag(i,0) >> ndiag(i,1);
    if (not inp) return false;
    ndirspin.resize(Norder);
    for (int i=0; i<Norder; ++i) inp >> ndirspin(i);
    if (not inp) return false;
    ntimes.resize(number_of2(Norder));
    for (int i=0; i<number_of2(Norder); ++i) inp >> ntimes(i);
    if (not inp) return false;
    for (int i=0; i<2*Norder; ++i){
      if (ntimes(i)<0 or ntimes(i)>beta){
	log<<"Note: status file contains frequency/time outside cutoff="<<beta<<", hence aborting reading" << endl;
	return false;
      }
    }
    if (ntimes(2*Norder)-ntimes(0)!=beta){
      log<<"Note correcting: time(-1)-beta="<<ntimes(2*Norder)-beta<<"!= time(0)="<<ntimes(0)<<endl;
      ntimes(2*Norder)=ntimes(0)+beta;
    }
    cdiag.resize(ndiag.shape()); cdiag=ndiag;
    cdirspin.resize(ndirspin.shape()); cdirspin=ndirspin;
    ctimes.resize(ntimes.shape()); ctimes=ntimes;
    log<<"Note: status file " << status_file << " succesfully read. Starting from configuration read." << endl;
    return true;
  }
  double Topology_trial_move(int& tPPs, int ii, const Array<int,2>& tdiag, const Array<int,1>& tdirspin, const Array<double,1>& ttimes, ostream& xout)
  {
    int n = tdiag.extent(0);
    // Compute Delta for all n propagators
    for (int i=0; i<n; ++i) {
      int dire = tdirspin(i) % 2;
      int iorb = tdirspin(i) / 2;
      int s = cix.bath_sigind(iorb);  // actual entry for hybridization
      double del_t = ((ttimes(tdiag(i,1))-ttimes(tdiag(i,0)))*(1-2*dire) + beta*dire);
      intpar ip = Interp(del_t, dltau);
      tDelta(i) = Deltas[s](ip);
    }
    double Phi_dd = Product(tDelta(Range(0,n-1)));
    
    Array<intpar,1> ips(2*n);
    for (int i=0; i<2*n; ++i){
      int ip1 = (i+1)%(2*n);
      double dt = ttimes(i+1)-ttimes(i);
      ips(i) = Interp(dt,pptau);
    }

    if (not Phi_diagonal){
      dirspin2ops(tops,/**/tdiag, tdirspin);
      PPstates1(tpstates,/**/tPPs,tdiag,tops,cix);
      
      tGm.resize(tpstates.extent(0),2*n);
      // Here we compute approximation for Phi as the product of largest eigenvalue in each block, which is  a proxy for matrix product.
      // It is often too large, because creation/anhilation operator kills off many terms. But it is extremely fast.
      tPhi.resize(tPPs);
      double Phi_pp=0;
      for (int ipp=0; ipp<tpstates.extent(0); ++ipp){
	int ib = tpstates(ipp,2*n-1);
	double cphi=1.0;
	for (int i=0; i<2*n; ++i){ 
	  int ii = tpstates(ipp,i);
	  double Gx = gms[ii](ips(i));
	  tGm(ipp,i) = Gx;
	  cphi *= Gx;
	}
	tPhi(ipp) = cphi*Phi_dd;
	Phi_pp += cphi;
      }
      //if (debug) cout<<"Final Phi_pp="<<Phi_pp<<endl;
      return Phi_pp*Phi_dd;
    }else{
      // Here we compute approximation for Phi by taking only diagonal PP-propagators into account.
      // This is often too small compared to the actual matrix product, but it is fairly close, and of course many times faster.
      //tGmi.resize(2*n);
      //auto& tallowed = tdiag2matrix.Cmp_diagonal(tdiag, tdirspin, cix);      // map<int,Array<int,2>>
      auto& tallowed = tdiag2matrix.Cmp_diagonal3(tdiag, tdirspin, cix);      // map<int,Array<int,2>>

      tPPs=0;
      //for (const auto& [ib, alw] : tallowed)  tPPs += alw.extent(1);
      for (const auto& [ib, alw] : tallowed)  tPPs += alw.extent(0);

      tGm.resize(tPPs,2*n); tGm=0;
      // Compute PP propagators for diagonal case
      int idx=0;
      for (const auto& [ib, alw] : tallowed){ // over all blocks that contribute on diagonal
	//cout<<"   alw["<<ib<<"]="<<alw<<endl;
	//for (int j=0; j<alw.extent(1); ++j){ // in principle more than one term can contribute per block
	for (int j=0; j<alw.extent(0); ++j){ // in principle more than one term can contribute per block
	  for (int i=0; i<2*n; ++i){             // over the ring of PPs
	    //int _ab_ = alw(i,j);
	    //int _a_ = abs(_ab_);
	    //int ab = cix.inBlock(_a_);        // _ab_ is given in terms of atomic state index, but ab is index within the block. cix.inBlock translates this.
	    //int _ib_ = cix.iBlck(_a_)-1;
	    //int ii = cix.Sgm_symmetry[_ib_](ab,ab)-1; // Green's function index for propagator G(a0,b0)==G[|_ii0_|-1]
	    int _ii_ = alw(j,i);
	    int sg = isign(_ii_);
	    int ii = sg*_ii_ - 1;
	    if (ii>=0){
	      //int sg = sign(_ab_);
	      tGm(idx,i) = Gms[ii](ips(i)) * sg; // this is actually G*F
	    }
	  }
	  idx += 1;
	}
      }
      //int sgns = SignSort(tdiag);    // sign due to topological crossing in the diagram
      tPhi.resize(tallowed.size());
      double Phi_pp=0;
      idx=0; int _idx_=0;
      for (const auto& [ib, alw] : tallowed){ // over all blocks that contribute on diagonal
	double phi=0;
	//for (int j=0; j<alw.extent(1); ++j) // in principle more than one term can contribute per block
	for (int j=0; j<alw.extent(0); ++j) // in principle more than one term can contribute per block
	  phi += Product(tGm(idx++,Range(0,2*n-1)));  // over the ring of PPs
	tPhi(_idx_++) = phi*Phi_dd;
	Phi_pp += phi;
      }
      //if (debug) cout<<"Final Phi_pp="<<Phi_pp<<endl;
      return Phi_pp*Phi_dd;
    }
  }
  
  double Dirspin_trial_move(int tPPs, int ii, const Array<int,2>& tdiag, const Array<int,1>& tdirspin, const Array<int,1>& old_dirspin, const Array<double,1>& ttimes, ostream& cout)
  { return Topology_trial_move(tPPs, ii, tdiag, tdirspin, ttimes, cout);}
  
  
  double FrequencyTime_trial_move(int itt, int tind, const Array<int,2>& diag, const Array<int,1>& dirspin, const Array<double,1>& ttimes,
				  const Array<double,1>& ctimes, ostream& cout)
  { // Routine that computes the value of Phi (within certain approximation), which is used for MC sampling
    //
    // Updates the \widetilde{G} and Delta arrays when time points change.
    // Only modifies entries that are changed, to improve performance.
    // The quantities from this class that are used, but not modified:
    //        vector<Spline1D<double>> Gms,
    //        vector<Spline1D<double>> Deltas
    //        Array<double,2> cGm
    //        Array<double,1> cDelta
    // The quantities from this class that are modified:
    //        Array<double,2> tGm
    //        Array<double,1> tDelta
    //
    int n = diag.extent(0);  // Number of hybridizations
    // First copy previous configuration
    tDelta(Range(0,n-1))=cDelta(Range(0,n-1));

    // Find the hybridization index corresponding to `tind`
    int hybind = -1;
    for (int i = 0; i<n; ++i)
      if (diag(i,0)==tind or diag(i,1)==tind){
	hybind=i;
	break;
      }
    if (hybind != -1) { // Ensure valid index
      int dir = dirspin(hybind) % 2;
      int iorb = dirspin(hybind) / 2;
      int s = cix.bath_sigind(iorb);
      double del_t = ((ttimes(diag(hybind,1))-ttimes(diag(hybind,0)))*(1-2*dir) + beta*dir);
      intpar ip = Interp(del_t, dltau);
      tDelta(hybind) = Deltas[s](ip);
    }else{
      cout<<"ERROR hybind=-1 for tind="<<tind<<" and 2*n="<<2*n<<" diag="<<diag<<endl;
    }
    double Phi_dd = Product(tDelta(Range(0,n-1)));      
    
    tGm.resize(cGm.shape());
    tGm=cGm;
    
    // Update the modified Gm values for two PP propagators,  the one before and the one after.
    TinyVector<intpar,2> ip;
    //TinyVector<int,2> _tind_;
    ip[0] = Interp(ttimes(tind+1)-ttimes(tind), pptau);
    _tind_[0] = tind;
    int tti = (tind!=0) ? tind : 2*n;
    ip[1] = Interp(ttimes(tti)-ttimes(tti-1), pptau);
    _tind_[1] = tti-1;
    
    if (not Phi_diagonal){
      // Here we compute approximation for Phi as the product of largest eigenvalue in each block, which is  a proxy for matrix product.
      // It is often too large, because creation/anhilation operator kills off many terms. But it is extremely fast.
      int PPs = ppstates.extent(0);
      tPhi.resize(PPs);
      for (int ipp=0; ipp<PPs; ++ipp){
	for (int ici=0; ici<2; ++ici){
	  int ii = ppstates(ipp,_tind_[ici]);
	  tGm(ipp,_tind_[ici]) = gms[ii](ip[ici]);
	}
      }
      double Phi_pp=0;
      for (int ipp=0; ipp<PPs; ++ipp){
	double cphi = Product(tGm(ipp,Range(0,2*n-1)));  // over the ring of PPs
	tPhi(ipp) = cphi*Phi_dd;
	Phi_pp += cphi;
      }
      //cout<<"Final Phi_pp="<<Phi_pp<<endl;
      return Phi_pp*Phi_dd;
    }else{
      auto& allowed = cdiag2matrix.diagonal;
      int idx=0;
      for (const auto& [ib, alw] : allowed){ // over all blocks that contribute on diagonal
	//for (int j=0; j<alw.extent(1); ++j){ // in principle more than one term can contribute per block
	for (int j=0; j<alw.extent(0); ++j){ // in principle more than one term can contribute per block
	  for (int ici=0; ici<2; ++ici){
	    int i = _tind_[ici];
	    //int _ab_ = alw(i,j);
	    //int _a_ = abs(_ab_);
	    //int ab = cix.inBlock(_a_);        // _ab_ is given in terms of atomic state index, but ab is index within the block. cix.inBlock translates this.
	    //int _ib_ = cix.iBlck(_a_)-1;
	    //int ii = cix.Sgm_symmetry[_ib_](ab,ab)-1; // Green's function index for propagator G(a0,b0)==G[|_ii0_|-1]
	    int _ii_ = alw(j,i);
	    int sg = isign(_ii_);
	    int ii = sg*_ii_ - 1;
	    if (ii>=0){
	      //int sg = sign(_ab_);
	      tGm(idx,i) = Gms[ii](ip[ici]) * sg; // this is actually G*F
	    }else{
	      tGm(idx,i) = 0;
	    }
	  }
	  idx += 1;
	}
      }
      tPhi.resize(allowed.size());
      double Phi_pp=0;
      idx=0; int _idx_=0;
      for (const auto& [ib, alw] : allowed){ // over all blocks that contribute on diagonal
	double phi=0;
	//for (int j=0; j<alw.extent(1); ++j) // in principle more than one term can contribute per block
	for (int j=0; j<alw.extent(0); ++j) // in principle more than one term can contribute per block
	  phi += Product(tGm(idx++,Range(0,2*n-1)));  // over the ring of PPs
	tPhi(_idx_++) = phi*Phi_dd;
	Phi_pp += phi;
      }
      //cout<<"Final Phi_pp="<<Phi_pp<<endl;
      return Phi_pp*Phi_dd;
    }
  }
  
  void Accept_topology_move(int ii, int PPs, int Norder)
  {
    if (not Phi_diagonal){
      ppstates.resize(tpstates.shape());  // Resizes ppstates to match tpstates
      ppstates = tpstates;                // Copies the data
      cops.resize(tops.shape());
      cops = tops;
    }else{
      cdiag2matrix = tdiag2matrix;
    }
    cPhi.resize(tPhi.shape());
    cPhi = tPhi;
    cGm.resize(tGm.shape());
    cGm = tGm;
    //cDelta(Range(0,Norder-1)) = tDelta(Range(0,Norder-1));
    cDelta.resize(tDelta.shape());
    cDelta = tDelta;
    topology_change = true;
  }
  void Accept_frequencyTime_move(int ii, int PPs, int Norder)
  {
    cGm.resize(tGm.shape());
    cGm = tGm;
    //cGm(Range(0,PPs-1),Range(0,2*Norder-1)) = tGm(Range(0,PPs-1),Range(0,2*Norder-1));
    cDelta.resize(tDelta.shape());
    cDelta = tDelta;
    //cDelta(Range(0,Norder-1)) = tDelta(Range(0,Norder-1));
    cPhi.resize(tPhi.shape());
    cPhi = tPhi;
    times_changed(_tind_[0])=true;
    times_changed(_tind_[1])=true;
    times_changed0(_tind_[0])=true;
  }
  void Accept_dirspin_move(int ii, int PPs, int Norder, const Array<int,1>& dirspin, const Array<int,1>& old_dirspin)
  {
    Accept_topology_move(ii, PPs, Norder);
  }
  void Get_utl_Interpolations(Array<double,2>& utli, Array<double,2>& utl,
			      const SVDFunc& svdf, const Array<int,2>& diag, const Array<int,1>& dirspin, const Array<double,1>& times, double beta)
  {
    int n = diag.extent(0);
    if (topology_change) {
      utli.resize(2*n,svdf.lmax);
      utli = 0.0;
    }
    // Compute utli values
    for (int ist = 0; ist < 2*n; ++ist){
      if (topology_change or times_changed(ist)){
	double dtau = beta - (times(ist+1)-times(ist));
	for (int l=0; l<svdf.lmax; ++l)
	  utli(ist,l) = svdf.fU[l](Interp(dtau, svdf.tau));
      }
    }
    // Compute utl values
    if (topology_change) {
      utl.resize(n,svdf.lmax);
      utl = 0.0;
    }
    for (int ihyb = 0; ihyb < n; ++ihyb) {
      int i1 = diag(ihyb,1), i0 = diag(ihyb,0);
      if (topology_change or times_changed0(i1) or times_changed0(i0)){
	int dire = dirspin(ihyb) % 2;
	int spin = dirspin(ihyb) / 2;
	double dtau = (times(i1)-times(i0))*(1-2*dire) + dire*beta;
	double d_tau = beta - dtau;
	for (int l=0; l<svdf.lmax; ++l)
	  utl(ihyb, l) = svdf.fU[l](Interp(d_tau, svdf.tau));
      }
    }
  }

  template<class RNG>
  void Meassure1(Array<double,1>& Phi_accumulate, Array<double,1>& Phi_accumulatef, Array<double,3>& Stval,
		 Array<double,3>& G00, Array<double,3>& F00,  Array<double,3>& F01,
		 double& Phi_cOCA, double& Phi_cOCAf, Array<double,2>& Stval_cOCA,  // added by codex
		 Array<double,2>& G00_cOCA, Array<double,2>& F00_cOCA, Array<double,2>& F01_cOCA,
		 Array<double,1>& N_measure,
		 int itt, double sumPhi, const Array<int,2>& diag, const Array<int,1>& dirspin,
		 const Array<double,1>& times, int PPs, const Array<double,1>& weight,
		 const Array<double,1>& sqrtbetalist, RNG& drand, ostream& cout)
  {
    const bool debug=false;
    int n = diag.extent(0);
    int minorder = N_min_order;
    int n_mino = n-minorder;
    N_measure(n_mino) += 1;
    if (DEBUG) cout<<"   physical_visits["<<n_mino<<"]="<<N_measure(n_mino)<<endl;
    string ccstr=" ";
    bool cOCA_needed = (n==2 and OCA_cut_sum<10);
	
    double Phi_dd = Product(cDelta(Range(0,n-1)));
#ifdef _TIMEIT_     
    auto t0 = get_cpu_time();
#endif
    //Array<double,2> utli, utl;
    Get_utl_Interpolations(utli, utl, // here we are interpolating (beta-(t_{i+1}-t_i)), and t'-t for hybridization.
			   svdf, diag, dirspin, times, beta);
    // here we interpolate t_{i+1}-t_i. This is related to beta-(t_{i+1}-t_i) but not directly related...
    if (topology_change) icps.resize(2*n);  //Array<intpar,1> icps(2*n);
    for (int i=0; i<2*n; ++i){
      if (topology_change or times_changed(i)){
	icps(i) = Interp(times(i+1)-times(i),pptau);
      }
    }
    if (DEBUG){
      cout<<"   utl="<<utl<<endl;
      cout<<"   utli="<<utli<<endl;
      flush(cout);
    }
    double phi_weight = weight(n_mino)/sqrtbetalist(n_mino);
    double W = abs(sumPhi) / phi_weight;

#ifdef _TIMEIT_ 
    auto t1 = get_cpu_time();
    t_intrp += t1-t0;
#endif
    //map<int,vector<Array<int,2>>>* allowem;
    if (topology_change){
      auto& allowem = cdiag2matrix.Cmp_matrix(diag, dirspin, cix);
      Lft.resize(allowem.size(),2*n); // Array<Array<double,2>,2> Lft(:,:)(N1,N2)
      Rgh.resize(allowem.size(),2*n); // Array<Array<double,2>,2> Rgh(:,:)(N1,N1)
      Gmx.resize(allowem.size(),2*n); // Array<Array<double,2>,2>  Gx(:,2*Norder)(N1,N2)
      if (cOCA_needed) remember_ii.resize(allowem.size(),2*n);
      int idx=0;
      for (const auto& [ib, alw] : allowem){
	int N0 = alw[0].extent(0);
	//Lft(idx,0).resize(N0,N0);
	for (int i=0; i<2*n; ++i){
	  int ip1 = (i+1)%(2*n);
	  int N1 = alw[i].extent(0), N2 = alw[ip1].extent(0);
	  Gmx(idx,i).resize(N1,N2);
	  Lft(idx,i).resize(N0,N1);
	  // remember_ii stores the Green-function index corresponding to Gmx(idx,i)(j1,j2)
	  if (cOCA_needed) remember_ii(idx,i).resize(N1,N2);
	}
	for (int i=1; i<2*n; ++i) Rgh(idx,i).resize(alw[2*n-i].extent(0),N0);
	idx += 1;
      }
      //
      Fi.resize(2*n); // vector<Array<double,1>>  Fi[2*Norder](N1)
      dS.resize(cix.Npseudo_unique,2*n);
      //
      sgns = SignSort(diag);    // sign due to topological crossing in the diagram
      //
      F_use_ie=true; // F_use_ie is true when vertex 0 does not have incoming hybridization. Needed for computing F.
      for (int ihyb = 0; ihyb<n; ++ihyb) {
	int dire = dirspin(ihyb) % 2;
	int ie = diag(ihyb,1-dire);
	if (ie==0) F_use_ie=false;
      }
      // WAS BUG Oct 16 2025
      Gtm.resize(allowem.size(),n);
      Gtp.resize(allowem.size(),n);
      if (cOCA_needed){
	_dS_.resize(cix.Npseudo_unique,4);
	importance_m.resize(allowem.size(),2);
	importance_p.resize(allowem.size(),2);
      }
    }
#ifdef _TIMEIT_
    auto t2 = get_cpu_time();
    t_allowm += t2-t1;
#endif    
    dS = 0.0;
    if (cOCA_needed) _dS_=0;
    Array<double,2> tmp, tmp1, tmp_storage(cix.Nmaxsize,cix.Nmaxsize), tmp1_storage(cix.Nmaxsize,cix.Nmaxsize);
    Array<double,2> Phis_over_Gm;
    Array<double,1> Phi_F(n); Phi_F=0;
    Array<double,1> Phi_F_OCA(n); Phi_F_OCA=0;
    Array<double,2> FGx, GxF;
    double tot_Phi=0, tot_Phi_OCA=0;
    int idx=0;
    for (const auto& [ib, alw] : cdiag2matrix.allowed){
      double phiL=0, phiR=0, phi_diag=0;
      // We would like to know how is cPhi(ii) related to phiL and phiR through cdiag2matrix
      //// Here we want to update only Gx[i] and Fi[i] that were changed recently.
      //// The information of which time needs to be updated is contained in times_changed(0...2*n-1)
      //// If there was topological change, than we recompute eveything (topology_change).
      //// If only few times are changed, times_changed(i) are true
#ifdef _TIMEIT_
      auto t3 = get_cpu_time();
#endif	
      for (int i=0; i<2*n; ++i){
	int N1 = alw[i].extent(0);
	Fi[i].resize(N1);
	for (int j1=0; j1<N1; ++j1) Fi[i](j1)=sign(alw[i](j1,1)); // matrix element of operator op(0) <alw[0](j0,1)|op(0)|alw[0](j0,0)> is stored in the sign of op(0) acting on index, hence sign of alw[0](j0,1) contains op(0) matrix element.
	if (topology_change or times_changed(i)){
	  int ip1 = (i+1)%(2*n);
	  intpar ip = icps(i);
	  int N2 = alw[ip1].extent(0);
	  if (topology_change){
	    //Gmx(idx,i).resize(N1,N2);
	    Gmx(idx,i)=0;
	  }
	  int ib_a = cix.iBlck(abs(alw[i](0,1)))-1;
	  for (int j1=0; j1<N1; ++j1){
	    int _a_ = alw[i](j1,1);
	    int a = cix.inBlock(abs(_a_)); // alw[0] index is given in terms of atomic state index, but a0 is index within the block. cix.inBlock translates this.
	    //double mm0 = sign(_a_);         // matrix element of operator op(0) <alw[0](j0,1)|op(0)|alw[0](j0,0)> is stored in the sign of op(0) acting on index, hence sign of alw[0](j0,1) contains op(0) matrix element.
	    //Fi[i](j1)=mm0;
	    for (int j2=0; j2<N2; ++j2){
	      int _b_ = alw[ip1](j2,0);
	      int b = cix.inBlock(abs(_b_));
	      int _ii0_ = cix.Sgm_symmetry[ib_a](a,b); // Green's function index for propagator G(a0,b0)==G[|_ii0_|-1]
	      if (cOCA_needed) remember_ii(idx,i)(j1,j2) = _ii0_;
	      if (_ii0_==0) continue;
	      int ii_sign0 = sign(_ii0_);                 // Green's function value for off-diagonal propagators could be negative/opposite, hence sign
	      int ii = abs(_ii0_)-1;
	      Gmx(idx,i)(j1,j2) = ii_sign0*Gms[ii](ip);
	      //G_importance[ii]; 
	    }
	  }
	}
      }
#ifdef _TIMEIT_
      auto t4 = get_cpu_time();	  
      t_Gx += t4-t3;
#endif
      int N0 = alw[0].extent(0);
      bool left_i_changed = topology_change;
      if (left_i_changed){
	//Lft(idx,0).resize(N0,N0);
	Lft(idx,0)=0;
	for (int j=0; j<N0; ++j)  // Lft(:,0) is just set to F
	  Lft(idx,0)(j,j) = sign(alw[0](j,1));
      }
      
#ifdef _TIMEIT_
      auto t5 = get_cpu_time();	  
#endif	
      for (int i=0; i<2*n; ++i){
	left_i_changed = left_i_changed or times_changed(i);
	if (i<2*n-1){
	  if (left_i_changed){
	    int ip1 = (i+1)%(2*n);
	    int N1 = alw[i].extent(0), N2 = alw[ip1].extent(0);
	    GxF.resize(N1,N2);
	    tiny::dot_diag(GxF, Gmx(idx,i), Fi[ip1]);  // GxF = Gx[i] @ diag(Fi[ip1])
	    tiny::dot(Lft(idx,i+1), Lft(idx,i), GxF);   //Lft(idx,i+1) = matmul(Lft(idx,i),GxF);??????
	  }
	}else{
	  phiL = tiny::trace_dot(Lft(idx,i),Gmx(idx,i));
	}	
      }
#ifdef _TIMEIT_
      auto t6 = get_cpu_time();	
      t_Lft += t6-t5;
#endif	
      bool right_i_changed = topology_change or times_changed(2*n-1);
      if (right_i_changed){
	int N1 = alw[2*n-1].extent(0);
	//Rgh(idx,1).resize(N1,N0);
	tiny::dot_diag(Rgh(idx,1), Fi[2*n-1], Gmx(idx,2*n-1)); // Rgh(idx,1) = diag(Fi[2*n-1]) @ Gx[2*n-1]
      }
      for (int i=1; i<2*n; ++i){
	right_i_changed = right_i_changed or times_changed(2*n-1-i);
	if (right_i_changed or i==2*n-1){
	  int N1 = alw[2*n-1-i].extent(0), N2 = alw[2*n-i].extent(0);
	  FGx.resize(N1,N2);
	  tiny::dot_diag(FGx, Fi[2*n-1-i], Gmx(idx,2*n-1-i)); // FGx = diag(Fi[2*n-1-i]) @ Gx[2*n-1-i]
	  if (i<2*n-1){
	    //Rgh(idx,i+1).resize(N1,N0);
	    tiny::dot(Rgh(idx,i+1), FGx, Rgh(idx,i)); // Rgh(idx,i+1) = matmul(FGx,Rgh(idx,i));
	  }else{
	    phiR = tiny::trace_dot(FGx,Rgh(idx,i));
	  }
	}
      }
#ifdef _TIMEIT_
      auto t7 = get_cpu_time();
      t_Rgh += t7-t6;
#endif      
      if (debug){
	if (not is_close(phiL,phiR)) { cout << "ERROR phiL=" << phiL << " while phiR=" << phiR << endl;}
	cout<<ib<<"   phiL="<<phiL<<" phiR="<<phiR<<" phi_diag="<<phi_diag<<endl;
	if (phi_diag!=0 and phi_diag<1e-8) ccstr="*";
	else ccstr=" ";
      }
      double Phi_this_block = 0.5*(phiL+phiR)*(Phi_dd/W)*sgns;
      tot_Phi += Phi_this_block;
      double Phi_this_block_OCA=0;
      if (cOCA_needed and Phi_this_block!=0){
#ifdef _TIMEIT_ 	
	auto t1 = get_cpu_time();
#endif
	const TinyVector<int,4> Nm(alw[0].extent(0), alw[1].extent(0), alw[2].extent(0), alw[3].extent(0));
	if (true){
	  for (int j0=0; j0<Nm[0]; ++j0){
	    for (int j1=0; j1<Nm[1]; ++j1){
	      int _ii0_ = remember_ii(idx,0)(j0,j1);
	      if (_ii0_==0) continue;
	      int ii0 = abs(_ii0_)-1;
	      int imp_ii0 = G_importance[ii0];
	      double g0 = Gmx(idx,0)(j0,j1);
	      double _f0g0_ = Fi[0](j0)*g0 * (Phi_dd/W)*sgns;
	      double sgn_ii0_g0 = sign(_ii0_)/g0;
	      for (int j2=0; j2<Nm[2]; ++j2){
		int _ii1_ = remember_ii(idx,1)(j1,j2);
		if (_ii1_==0) continue;
		int ii1 = abs(_ii1_)-1;
		int imp_ii1 = G_importance[ii1];
		if (imp_ii0 + imp_ii1 > OCA_cut_sum) continue;
		double g1 = Gmx(idx,1)(j1,j2);
		double f1g1 = Fi[1](j1)*g1;
		double _f0g0_f1g1_ = _f0g0_*f1g1;
		double sgn_ii1_g1 = sign(_ii1_)/g1;
		for (int j3=0; j3<Nm[3]; ++j3){
		  int _ii2_ = remember_ii(idx,2)(j2,j3);
		  if (_ii2_==0) continue;
		  int ii2 = abs(_ii2_)-1;
		  int imp_ii2 = G_importance[ii2];
		  int _ii3_ = remember_ii(idx,3)(j3,j0);
		  if (_ii3_==0) continue;
		  int ii3 = abs(_ii3_)-1;
		  int importance = imp_ii0 + imp_ii1 + imp_ii2 + G_importance[ii3];
		  if (importance <= OCA_cut_sum ){
		    double g2 = Gmx(idx,2)(j2,j3);
		    double f2g2 = Fi[2](j2)*g2;
		    double g3 = Gmx(idx,3)(j3,j0);
		    double f3g3 = Fi[3](j3)*g3;
		    //double phiLR = Fi[0](j0)*Gmx(idx,0)(j0,j1)*Fi[1](j1)*Gmx(idx,1)(j1,j2)*Fi[2](j2)*Gmx(idx,2)(j2,j3)*Fi[3](j3)*Gmx(idx,3)(j3,j0);
		    //double Phi_this_block = phiLR*(Phi_dd/W)*sgns;
		    double Phi_this_block = _f0g0_f1g1_*f2g2*f3g3;
		    Phi_this_block_OCA += Phi_this_block;
		    _dS_(ii0,0) += sgn_ii0_g0*Phi_this_block;
		    _dS_(ii1,1) += sgn_ii1_g1*Phi_this_block;
		    _dS_(ii2,2) += sign(_ii2_)/g2*Phi_this_block;
		    _dS_(ii3,3) += sign(_ii3_)/g3*Phi_this_block;
		  }
		}
	      }
	    }
	  }
	  //Phi_cOCA += Phi_this_block_OCA;
	  tot_Phi_OCA += Phi_this_block_OCA;
	}
#ifdef _TIMEIT_  
	auto t2 = get_cpu_time();
	t_cOCA1 += t2-t1;
#endif
      }
      // End inserted by codex
      // computing self-energy and writing it into compact array.
      for (int ist=0; ist<2*n; ++ist){
	int ip1 = (ist+1)%(2*n);
	int N1 = alw[ist].extent(0), N2 = alw[ip1].extent(0);
	int ib_a = cix.iBlck(abs(alw[ist](0,1)))-1;
	if (N1==1 and N2==1){
	  if (Gmx(idx,ist)(0,0)!=0){
	    double Phi_over_Gm = Phi_this_block/Gmx(idx,ist)(0,0);
	    int a = cix.inBlock(abs(alw[ist](0,1))); // alw[0] index is given in terms of atomic state index, but a0 is index within the block. cix.inBlock translates this.
	    int b = cix.inBlock(abs(alw[ip1](0,0)));
	    int _ii0_ = cix.Sgm_symmetry[ib_a](b,a); // Green's function index for propagator G(a0,b0)==G[|_ii0_|-1]
	    if (_ii0_!=0){
	      int ii_sign0 = sign(_ii0_);            // Green's function value for off-diagonal propagators could be negative/opposite, hence sign
	      int ii = abs(_ii0_)-1;                 // correction to which PP self-energy
	      if (debug) cout<<"   "<<ccstr<<"P(0,0)=S("<<b<<","<<a<<")=S("<<_ii0_<<")="<<Phi_over_Gm<<endl;
	      dS(ii,ist) += ii_sign0*Phi_over_Gm;    // ii is PP index, ist is the time-difference index.
	      //Stval(n_mino,ii,Range::all()) += utli(ist,Range::all())*ii_sign0*Phi_over_Gm;
	    }
	  }
	}else{
	  if (ist<2*n-1){
	    Phis_over_Gm.resize(Rgh(idx,2*n-1-ist).extent(0), Lft(idx,ist).extent(1));
	    tiny::dot(Phis_over_Gm, Rgh(idx,2*n-1-ist), Lft(idx,ist) ); // Right*Left gives the self-energy contribution, i.e., roughly Phi/G
	  }else{
	    Phis_over_Gm.resize(Lft(idx,ist).shape());
	    Phis_over_Gm = Lft(idx,ist);
	  }
	  Phis_over_Gm *= (Phi_dd/W)*sgns; // also need sign and hybridization part, and MC-fake-weight.
	  if (debug){
	    cout<<" P.shape="<<Phis_over_Gm.shape()<<" Rgh("<<idx<<","<<2*n-1-ist<<")="<<Rgh(idx,2*n-1-ist).shape()<<" Lft("<<idx<<","<<ist<<")="<<Lft(idx,ist).shape();
	    cout<<" N1="<<N1<<" N2="<<N2<<" ib="<<ib_a+1<<endl;
	  }
	  for (int j1=0; j1<N1; ++j1){
	    int _a_ = alw[ist](j1,1);
	    int a = cix.inBlock(abs(_a_)); // alw[0] index is given in terms of atomic state index, but a0 is index within the block. cix.inBlock translates this.
	    for (int j2=0; j2<N2; ++j2){
	      int _b_ = alw[ip1](j2,0);
	      int b = cix.inBlock(abs(_b_));
	      int _ii0_ = cix.Sgm_symmetry[ib_a](b,a); // Green's function index for propagator G(a0,b0)==G[|_ii0_|-1]
	      if (_ii0_==0) continue;
	      int ii_sign0 = sign(_ii0_);              // Green's function value for off-diagonal propagators could be negative/opposite, hence sign
	      int ii = abs(_ii0_)-1;                   // correction to which PP self-energy
	      if (debug) cout<<"   "<<ccstr<<"P("<<j2<<","<<j1<<")=S("<<_b_<<","<<_a_<<")=S("<<_ii0_<<")="<<Phis_over_Gm(j2,j1)<<endl;
	      dS(ii,ist) += ii_sign0*Phis_over_Gm(j2,j1); // ii is PP index, ist is the time-difference index.
	      //Stval(n_mino,ii,Range::all()) += utli(ist,Range::all())*ii_sign0*Phis_over_Gm(j2,j1);
	    }
	  }
	}
      }
#ifdef _TIMEIT_
      auto t8 = get_cpu_time();      
      t_stval += t8-t7;
#endif
      if (cmpPhysG){
	// Below we compute F vertex function for the self-energy calculation Sigma(w)*G(w) = F(w)
	int N_2=2*n;
	Array<bool,1> needed(cix.Nmaxsize); 
	if (topology_change){
	  for (int ihyb = 0; ihyb<n; ++ihyb){
	    int dire = dirspin(ihyb) % 2;
	    int iorb = dirspin(ihyb) / 2;
	    int ic, ibm, ibp, ibc;
	    if (F_use_ie){ // replacing operator ie with Njm, because Left&Right allow to do for every hybridization
	      int ie = diag(ihyb,1-dire);
	      ic = ie;
	      ibm = cix.iBlck(abs(alw[ie](0,0)))-1;// a bit complicated way to find which block corresponds to this ist
	      ibp = cix.iBlck(abs(alw[ie](0,1)))-1;// a bit complicated way to find which block corresponds to this ist
	      ibc = ibm;
	    }else{ // replacing operator is with Njm, because Left&Right allow to do for every hybridization
	      int is = diag(ihyb,dire);
	      ic = is;
	      ibm = cix.iBlck(abs(alw[is](0,0)))-1;// a bit complicated way to find which block corresponds to this ist
	      ibp = cix.iBlck(abs(alw[is](0,1)))-1;// a bit complicated way to find which block corresponds to this ist
	      ibc = ibp;
	    }
	    int ic_m_1 = (ic-1+N_2)%N_2;
	    int ic_p_1 = (ic+1)%N_2;
	    //int njmc = cix.Njm_c(ibc,iorb); // index to identify whether njm is proportional to F or not.
	    NjmResult nx = F_use_ie ? cix.Njm_ind(ibm,2*iorb) : cix.Njm_ind(ibm,2*iorb+1);
	    if (nx.njm_c==2){                   // njm is proportional to F, hence simple
	      Gtm(idx,ihyb).resize(alw[ic_m_1].extent(0),cix.Block[ibm].size());
	      Gtm(idx,ihyb)=0;
	      Gtp(idx,ihyb).resize(cix.Block[ibp].size(),alw[ic_p_1].extent(0));
	      Gtp(idx,ihyb)=0;
	      if (cOCA_needed){
		importance_m(idx,ihyb).resize(alw[ic_m_1].extent(0),cix.Block[ibm].size());
		importance_p(idx,ihyb).resize(cix.Block[ibp].size(),alw[ic_p_1].extent(0));
		importance_m(idx,ihyb)=0;
		importance_p(idx,ihyb)=0;
	      }
	    }
	  }
	  //cout<<"Gtp and Gtm resized"<<endl;
	}
	// Need to compute Trace(Phis_over_Gm * Gtm * njmm * Gtp)
	//Phis_over_Gm.resize(2*n-1, Rgh(idx,2*n-1-ic).extent(0), Lft(idx,ic-1).extent(1) );
	//tmp.resize(alw[ic_m_1].extent(0),cix.Block[ibp].size());
	//tmp1.resize(alw[ic_m_1].extent(0),alw[(ic+1)%(2*n)].extent(0));
	//tmp = Gtm * njmm;
	//tmp1= tmp * Gtp;
	
	for (int ihyb = 0; ihyb<n; ++ihyb) {
	  int dire = dirspin(ihyb) % 2;
	  int iorb = dirspin(ihyb) / 2;
	  int s = cix.bath_sigind(iorb);
	  // if dt>0 : G00(-t) has a minus sign while Delta(t) is positive
	  // if dt<0 : Delta(t) has minus sign, while G00(-t) is positive.
	  // We did not add the minus sign in Delta, hence we need to account for both minus signs here,
	  // which leads to overal (-1)
	  //
	  //                --->---
	  //              /         \
	  //  pseud[ic]  | ibp       | 
	  //    ie or is |-----------| is or is
	  //             | ibm       | pseud[is]
	  //  pseud[ic-1] \         /  
	  //                ---<---
	  //   Notice that for ic==ie, we have Njmm==Njm[ibm][iorb], while for
	  //                   ic==is, we have Njmm==Njm[ibp][iobr]^T.
	  //   This is because psi^+_a should be replaced by -psi^+_i psi^+_j psi_k * U_{i j k a}   while
	  //                   psi_a   should be replaced by -psi^+_k psi_j psi_i   * U_{i j k a}
	  //             which requires the matrix Njm to be transposed (in general hermitian conjugated)
	  int is = diag(ihyb,dire);
	  int ie = diag(ihyb,1-dire);
	  int ic, ibc, ibm, ibp;
	  if (F_use_ie){ // replacing operator ie with Njm, because Left&Right allow to do for every hybridization
	    ic = ie;
	    ibm = cix.iBlck(abs(alw[ie](0,0)))-1;// a bit complicated way to find which block corresponds to this ist
	    ibp = cix.iBlck(abs(alw[ie](0,1)))-1;// a bit complicated way to find which block corresponds to this ist
	    ibc = ibm;
	  }else{ // replacing operator is with Njm, because Left&Right allow to do for every hybridization
	    ic = is;
	    ibm = cix.iBlck(abs(alw[is](0,0)))-1;// a bit complicated way to find which block corresponds to this ist
	    ibp = cix.iBlck(abs(alw[is](0,1)))-1;// a bit complicated way to find which block corresponds to this ist
	    ibc = ibp;
	  }
	  //int njmc = cix.Njm_c(ibc,iorb); // index to identify whether njm is proportional to F or not.
	  NjmResult nx = F_use_ie ? cix.Njm_ind(ibm,2*iorb) : cix.Njm_ind(ibm,2*iorb+1);
	  if (nx.njm_c==1){                   // njm is proportional to F, hence simple
	    //double njmr = cix.Njm_r.at(cix.encode_block_bath(ibc,iorb));
	    //Phi_F(ihyb) += Phi_this_block*njmr;
	    Phi_F(ihyb) += Phi_this_block*nx.njm_r;
	  }else if (nx.njm_c==2){              // njm is not proportional
	    // Need to compute Trace(Phis_over_Gm * Gtm * njmm * Gtp)
	    // where Phis_over_Gm = Rgh(idx,2*n-1-ic) * Lft(idx,ic-1)
	    // and Gtm == Gm(ic-1) and Gp == Gm(ic) and ic is either is or ie
#ifdef _TIMEIT_
	    auto t0 = get_cpu_time();      
#endif	    
	    Array<double,2> njmm;          // matrix when njm is not proportional to F
	    njmm.reference(*nx.njm);
	    bool trans_njmm = nx.transposed;
	    //cout << "Setting trans_njmm="<<trans_njmm<<" for ib=" << ib << endl;
	    int ic_m_1 = ic>0 ? ic-1 : 2*n-1;
	    if (ic==2*n-1){
	      Phis_over_Gm.resize(Lft(idx,ic_m_1).shape());
	      Phis_over_Gm = Lft(idx,ic_m_1);
	      //cout<<" ic==2n-1 and after dot Phis_over_Gm.shape="<<Phis_over_Gm.shape()<<" ibm="<<ibm+1<<" ibp="<<ibp+1<<" ibc="<<ibc+1;
	      //cout<<" alw["<<ic_m_1<<"]="<<alw[ic_m_1]<<" alw["<<(ic+1)%(2*n)<<"]="<<alw[(ic+1)%(2*n)]<<endl;
	    }else if (ic==0)
	      // we don't have a good formula
	      cout<<"ERROR F_use_ie="<<F_use_ie<<" but ic==0 at ihyb="<<ihyb<<" with diag="<<diag<<" and dirspin="<<dirspin<<endl;
	    else{
	      //cout<<"dot at ic="<<ic<<" Rgh["<<idx<<","<<2*n-1-ic<<"].shape="<<Rgh(idx,2*n-1-ic).shape()<<" Lft["<<idx<<","<<ic-1<<"].shape="<<Lft(idx,ic-1).shape()<<" ibm="<<ibm+1<<" ibp="<<ibp+1<<" ibc="<<ibc+1<<endl;
	      Phis_over_Gm.resize(Rgh(idx,2*n-1-ic).extent(0),Lft(idx,ic-1).extent(1));
	      tiny::dot(Phis_over_Gm, Rgh(idx,2*n-1-ic), Lft(idx,ic-1));
	      //cout<<" after dot Phis_over_Gm.shape="<<Phis_over_Gm.shape()<<endl;
	    }
	    //cout<<"Phis_over_Gm.shape="<<Phis_over_Gm.shape()<<endl;
	    if (topology_change or times_changed(ic_m_1)){
#ifdef _TIMEIT_
	      auto t1 = get_cpu_time();
	      t_fdot += t1-t0;
#endif	    
	      int im=ic_m_1, ip=ic;
	      int N1 = alw[im].extent(0), N2 = alw[ip].extent(0);
	      int Nm = cix.Block[ibm].size();
	      needed(Range(0,Nm-1))=true;
	      for (int j2=0; j2<N2; ++j2){
		int b = cix.inBlock(abs(alw[ip](j2,0)));
		Gtm(idx,ihyb)(Range::all(),b) = Gmx(idx,im)(Range::all(),j2);
		needed(b)=false;
		//
		if (cOCA_needed){
		  for (int j1=0; j1<N1; ++j1){
		    int _ii_ = remember_ii(idx,im)(j1,j2);
		    if (_ii_!=0) importance_m(idx,ihyb)(j1,b)  = G_importance[abs(_ii_)-1];
		  }
		}
	      }
#ifdef _TIMEIT_
	      auto t2 = get_cpu_time();
	      t_gcpy += t2-t1;
#endif
	      for (int j1=0; j1<N1; ++j1){
		int a = cix.inBlock(abs(alw[im](j1,1)));
		for (int b=0; b<Nm; ++b){
		  if (needed(b)){
		    int _ii0_ = cix.Sgm_symmetry[ibm](a,b);
		    if (_ii0_!=0){
		      int ii_sign0 = sign(_ii0_);
		      int ii = abs(_ii0_)-1;
		      Gtm(idx,ihyb)(j1,b) = ii_sign0*Gms[ii](icps(im));
		      if (cOCA_needed) importance_m(idx,ihyb)(j1,b) = G_importance[ii];
		    }
		  }
		}
	      }
#ifdef _TIMEIT_
	      auto t3 = get_cpu_time();
	      t_gintrp += t3-t2;
#endif
	    }
	    if (topology_change or times_changed(ic)){
#ifdef _TIMEIT_
	      auto t3 = get_cpu_time();
#endif	      
	      int im = ic, ip=(ic+1)%(2*n);
	      int N1 = alw[im].extent(0), N2 = alw[ip].extent(0);
	      int Np = cix.Block[ibp].size();
	      //Gtp.resize(Np,N2); Gtp=0;
	      //needed.resize(Np);
	      needed(Range(0,Np-1))=true;
	      for (int j1=0; j1<N1; ++j1){
		int a = cix.inBlock(abs(alw[im](j1,1)));
		Gtp(idx,ihyb)(a,Range::all()) = Gmx(idx,im)(j1,Range::all()); // HERE IS BUG
		needed(a)=false;
		if (cOCA_needed){
		  for (int j2=0; j2<N2; ++j2){
		    int _ii_ = remember_ii(idx,im)(j1,j2);
		    if (_ii_!=0)  importance_p(idx,ihyb)(a,j2) = G_importance[abs(_ii_)-1];
		  }
		}
	      }
#ifdef _TIMEIT_
	      auto t4 = get_cpu_time();
	      t_gcpy += t4-t3;
#endif	    
	      for (int a=0; a<Np; ++a){
		if (needed(a)){
		  for (int j2=0; j2<N2; ++j2){
		    int b = cix.inBlock(abs(alw[ip](j2,0)));
		    int _ii0_ = cix.Sgm_symmetry[ibp](a,b);
		    if (_ii0_!=0){
		      int ii_sign0 = sign(_ii0_);
		      int ii = abs(_ii0_)-1;
		      Gtp(idx,ihyb)(a,j2) = ii_sign0*Gms[ii](icps(im));
		      if (cOCA_needed) importance_p(idx,ihyb)(a,j2) = G_importance[ii];
		    }
		  }
		}
	      }
#ifdef _TIMEIT_
	      auto t5 = get_cpu_time();
	      t_gintrp += t5-t4;
#endif
	    }
	    //cout<<"Gtm&Gtp set, shape(Gtm)="<<Gtm(idx,ihyb).shape() << " shape(Gtp)="<<Gtp(idx,ihyb).shape() << " shape(njmm)="<< njmm.shape() << " trans_njmm="<< trans_njmm << endl;
#ifdef _TIMEIT_
	    auto t6 = get_cpu_time();
#endif
	    //cout<<" before dot(Gtm,njmm) shape(Gtm)="<< Gtm(idx,ihyb).shape() <<" shape(njmm)="<< njmm.shape() <<" trans_njmm="<< trans_njmm <<endl;
	    tmp.reference(tmp_storage);
	    tiny::dot(tmp, Gtm(idx,ihyb), njmm, false, trans_njmm);
	    //cout<<"dot(Gtm,njmm) done shape(tmp)="<< tmp.shape()<<" shape(Gtp)="<< Gtp(idx,ihyb).shape() << endl;
	    tmp1.reference(tmp1_storage);
	    tiny::dot(tmp1, tmp, Gtp(idx,ihyb));
	    //cout<<"dot(Gtm*njm,Gtp) done shape(tmp1)="<< tmp1.shape()<<" shape(Phis_over_Gm)="<< Phis_over_Gm.shape() << endl;
	    Phi_F(ihyb) += tiny::trace_dot(Phis_over_Gm, tmp1) * (Phi_dd/W)*sgns;
	    //cout<<"Phi_F extracted="<< Phi_F(ihyb) << endl;
#ifdef _TIMEIT_
	    auto t7 = get_cpu_time();
	    t_sdot += t7-t6;
#endif	    
	  }
	  // Now handle the case of approximate analytic OCA
	  if (cOCA_needed and Phi_this_block!=0){
#ifdef _TIMEIT_
	    auto t1 = get_cpu_time();
#endif	    
	    if (nx.njm_c==1){                   // njm is proportional to F, hence simple
	      Phi_F_OCA(ihyb) += Phi_this_block_OCA*nx.njm_r;
	    }else if (nx.njm_c==2){              // njm is not proportional
	      const TinyVector<int,4> Nm(alw[0].extent(0), alw[1].extent(0), alw[2].extent(0), alw[3].extent(0));
	      if (true){
		// things that are already defined above:
		// is, ie, F_use_ie, ic, ibc, ibm, ibp, nx
		Array<double,2> njmm;          // matrix when njm is not proportional to F
		njmm.reference(*nx.njm);
		bool trans_njmm = nx.transposed;
		int ic_m_1 = ic>0 ? ic-1 : 2*n-1;
		int ic_p_1 = (ic+1)%(2*n);
		int N1 = alw[ic_m_1].extent(0), N2 = alw[ic].extent(0), N3 = alw[(ic+1)%(2*n)].extent(0);
		int Nm = cix.Block[ibm].size();
		int Np = cix.Block[ibp].size();
		int ic_0 = (ic+2)%(2*n);
		int N0 = alw[ic_0].extent(0);
		Array<double,3> tmp1(N1,Np,OCA_cut_sum+1), tmp2(N1,N3,OCA_cut_sum+1), tmp3(N3,N1,OCA_cut_sum+1);
		Array<double,1> tmp0(OCA_cut_sum+1), rtmp(N0);
		for (int i=0; i<N1; ++i){
		  for (int j=0; j<Np; ++j){
		    tmp0=0;
		    if (trans_njmm){
		      for (int k=0; k<Nm; ++k){
			int imp_ii = importance_m(idx,ihyb)(i,k);
			tmp0(imp_ii) += Gtm(idx,ihyb)(i,k)*njmm(j,k);
		      }
		    }else{
		      for (int k=0; k<Nm; ++k){
			int imp_ii = importance_m(idx,ihyb)(i,k);
			tmp0(imp_ii) += Gtm(idx,ihyb)(i,k)*njmm(k,j);
		      }
		    }
		    tmp1(i,j,Range::all()) = tmp0(Range::all());
		  }
		}
		tmp2=0;
		for (int i=0; i<N1; ++i){
		  for (int j=0; j<N3; ++j){
		    tmp0=0;
		    for (int k=0; k<Np; ++k){
		      int imp_ii1 = importance_p(idx,ihyb)(k,j);
		      for (int imp_ii0=0; imp_ii0<=OCA_cut_sum-imp_ii1; ++imp_ii0)
			tmp0(imp_ii1+imp_ii0) += tmp1(i,k,imp_ii0)*Gtp(idx,ihyb)(k,j);
		    }
		    tmp2(i,j,Range::all()) = tmp0(Range::all());
		  }
		}
		tmp3 = 0;
		for (int j3=0; j3<N3; ++j3){
		  for (int k=0; k<N0; ++k)
		    rtmp(k) = Gmx(idx,ic_p_1)(j3,k) * Fi[ic_0](k);
		  for (int j1=0; j1<N1; ++j1){
		    tmp0=0;
		    for (int k=0; k<N0; ++k){
		      int _ii3_ = remember_ii(idx,ic_p_1)(j3,k);
		      int _ii0_ = remember_ii(idx,ic_0)(k,j1);
		      if (_ii3_==0 or _ii0_==0) continue;
		      int ii3 = abs(_ii3_)-1;
		      int ii0 = abs(_ii0_)-1;
		      int imp_ii = G_importance[ii3]+G_importance[ii0];
		      if (imp_ii>OCA_cut_sum) continue;
		      tmp0(imp_ii) += rtmp(k) * Gmx(idx,ic_0)(k,j1);
		    }
		    tmp3(j3,j1,Range::all()) += Fi[ic_p_1](j3)*tmp0(Range::all())*Fi[ic_m_1](j1);
		  }
		}
		double Phi_F_oca=0;
		for (int i=0; i<N1; ++i){
		  for (int j=0; j<N3; ++j){
		    for (int imp_ii2=0; imp_ii2<=OCA_cut_sum; ++imp_ii2)
		      for (int imp_ii3=0; imp_ii3<=OCA_cut_sum; ++imp_ii3)
			if (imp_ii2+imp_ii3<=OCA_cut_sum)
			  Phi_F_oca += tmp2(i,j,imp_ii2)*tmp3(j,i,imp_ii3);
		  }
		}
		Phi_F_OCA(ihyb) += Phi_F_oca*(Phi_dd/W)*sgns;
	      }
	    }
#ifdef _TIMEIT_
	    auto t2 = get_cpu_time();
	    t_cOCA2 += t2-t1;
#endif	    
	  } // cOCA_needed
	} // over ihyb
      } // end cmpPhysG
#ifdef _TIMEIT_
      auto t9 = get_cpu_time();
      t_f00 += t9-t8;
#endif      
      idx+=1;
    }

#ifdef _TIMEIT_
    auto t10 = get_cpu_time();
#endif
    for (int ii=0; ii<cix.Npseudo_unique; ++ii){
      for (int ist=0; ist<2*n; ++ist){
	if (dS(ii,ist)!=0 and dS(ii,ist)<1e100 and dS(ii,ist)>-1e100) // it seems sometimes nan can be added.
	  Stval(n_mino,ii,Range::all()) += utli(ist,Range::all())*dS(ii,ist);
	if (cOCA_needed and (_dS_(ii,ist)!=0 and _dS_(ii,ist)<1e100 and _dS_(ii,ist)>-1e100))
	  Stval_cOCA(ii,Range::all()) += utli(ist,Range::all())*_dS_(ii,ist);
      }
    }

#ifdef _TIMEIT_
    auto t11 = get_cpu_time();
    t_dS += t11-t10;
#endif    
    Phi_accumulate(n_mino) += tot_Phi;
    if (cOCA_needed){
      Phi_cOCA += tot_Phi_OCA;
      if (cmpPhysG) Phi_cOCAf += tot_Phi_OCA;
    }
    if (cmpPhysG){
      Phi_accumulatef(n_mino) += tot_Phi;
      for (int ihyb = 0; ihyb<n; ++ihyb){ // Loop over hybridization lines for physical G & F
	int dire = dirspin(ihyb) % 2;
	int iorb = dirspin(ihyb) / 2;
	int s = cix.bath_sigind(iorb);
	G00(n_mino,s,Range::all()) += utl(ihyb,Range::all()) * tot_Phi     / cDelta(ihyb);
	// I am not entirely sure this statement is correct. Needs extensive testing to understand its correctness!  
	bool destroys_Cooper_pair = cix.Fi_indx(2*iorb)%2==1 and cix.Fi_indx(2*iorb+1)%2==1;
	bool normal_down = cix.SC_bath(iorb) and ((cix.Fi_indx(2*iorb)+cix.Fi_indx(2*iorb+1))%2==1 and cix.bath_spin(iorb)<0);
	bool destroys_Cooper_pair_or_normal_down = destroys_Cooper_pair or normal_down;
	if (cix.SC_bath(iorb)){
	  // When F_use_ie=true we substitude the end of bath (i_e or psi_{i_e}^\dagger or end of bath in SC state) for [H,\psi_{i_e}^\dagger]
	  // and when F_use_ie=false we substitude the begining of bath (i_s or psi_{i_s} or beginning of bath operator in SC state) for [H,\psi_{i_s}].
	  // The first gives G*\Sigma and the second gives \Sigma*G.
	  // In normal state G and \Sigma are symmetric therefore it does not matter in which F we put. We just use F00 only.
	  // In SC state we split them:
	  //          if F_use_ie=true and we deal with F11 or F21 term in Nambu, it means we are doing equation of motion for creation operator.
	  //          if F_use_ie=true and we deal with F21 or F22 (destroys_Cooper_pair or normal_down), it means we are doing EOM for anhilation operator.
	  //          if F_use_ie=false and we deal with F11 or F21 term in Nambu, it means we are doing EOM on anhilation for F11 and creation for F21
	  //          if F_use_ie=false and we deal with F21 or F22, it means we are doing EOM for anhilation for F21 and creation for F22
	  if (F_use_ie){
	    if (destroys_Cooper_pair_or_normal_down)
	      F01(n_mino,s,Range::all()) += utl(ihyb,Range::all()) * Phi_F(ihyb) / cDelta(ihyb);
	    else
	      F00(n_mino,s,Range::all()) += utl(ihyb,Range::all()) * Phi_F(ihyb) / cDelta(ihyb);
	  }else{
	    if (destroys_Cooper_pair_or_normal_down)
	      F00(n_mino,s,Range::all()) += utl(ihyb,Range::all()) * Phi_F(ihyb) / cDelta(ihyb);
	    else
	      F01(n_mino,s,Range::all()) += utl(ihyb,Range::all()) * Phi_F(ihyb) / cDelta(ihyb);
	  }
	}else{ // If superconductivity is absent, all goes into F00.
	  F00(n_mino,s,Range::all()) += utl(ihyb,Range::all()) * Phi_F(ihyb) / cDelta(ihyb);
	}
	if (cOCA_needed){
	  G00_cOCA(s,Range::all()) += utl(ihyb,Range::all()) * tot_Phi_OCA/ cDelta(ihyb);
	  if (cix.SC_bath(iorb)){
	    if (F_use_ie){
	      if (destroys_Cooper_pair_or_normal_down){
		F01_cOCA(s,Range::all()) += utl(ihyb,Range::all()) * Phi_F_OCA(ihyb) / cDelta(ihyb);
		//_F01_cOCA(s,Range::all()) += utl(ihyb,Range::all()) * Phi_F_OCA(ihyb) / cDelta(ihyb);
	      }else{
		F00_cOCA(s,Range::all()) += utl(ihyb,Range::all()) * Phi_F_OCA(ihyb) / cDelta(ihyb);
		//_F00_cOCA(s,Range::all()) += utl(ihyb,Range::all()) * Phi_F_OCA(ihyb) / cDelta(ihyb);
	      }
	    }else{
	      if (destroys_Cooper_pair_or_normal_down){
		F00_cOCA(s,Range::all()) += utl(ihyb,Range::all()) * Phi_F_OCA(ihyb) / cDelta(ihyb);
		//_F02_OCA(s,Range::all()) += utl(ihyb,Range::all()) * Phi_F_OCA(ihyb) / cDelta(ihyb);
	      }else{
		F01_cOCA(s,Range::all()) += utl(ihyb,Range::all()) * Phi_F_OCA(ihyb) / cDelta(ihyb);
		//_F03_OCA(s,Range::all()) += utl(ihyb,Range::all()) * Phi_F_OCA(ihyb) / cDelta(ihyb);
	      }
	    }
	  }else{
	    F00_cOCA(s,Range::all()) += utl(ihyb,Range::all()) * Phi_F_OCA(ihyb) / cDelta(ihyb);
	  }
	}
      }
    }
#ifdef _TIMEIT_
    auto t12 = get_cpu_time();
    t_g00 += t12-t11;
#endif    
    if (DEBUG){
      //cout<<"   ppstates="<<ppstates<<" dirspin="<<dirspin<<" cPhi="<<cPhi<<endl;
      cout<<"   Vph["<<n_mino<<"]="<<Phi_accumulate(n_mino)<<endl;
      cout<<"   cSm["<<n_mino<<"]="<<Stval(n_mino,Range::all(),Range::all())<<endl;
      cout<<"   c_G["<<n_mino<<"]="<<G00(n_mino,Range::all(),Range::all())<<endl;
      flush(cout);
    }
    topology_change=false;
    times_changed=false;
    times_changed0=false;

    //if (cOCA_needed){
    //  if ( sum(abs(Phi_F(Range::all())-Phi_F_OCA(Range::all()))) > 1e-6 )
    //	for (int ihyb=0; ihyb<2; ++ihyb)
    //	  cout << " "<<setw(4)<<itt<<" Phi_F-Phi_F_OCA="<<abs(Phi_F(ihyb)-Phi_F_OCA(ihyb))<<" , ["<<Phi_F(ihyb)<<" "<<Phi_F_OCA(ihyb)<<"]"<<endl;
    //}
    //if (cOCA_needed){
    //  if (abs(Phi_accumulatef(n_mino)-Phi_cOCA)/abs(Phi_accumulatef(n_mino))>1e-4)
    //	cout << " "<<setw(4)<<itt<<" phi_exact(2)="<<Phi_accumulatef(n_mino)<<" phi_mc(2)="<<Phi_cOCA<<endl;
    //}
  }
  void Resize_Output(Array<double,2>& Sm, Array<double,2>& G00, Array<double,2>& F00, Array<double,2>& F01,
		     const Cix& cix, bool cmpPhysG)
  {
    Sm.resize(cix.Npseudo_unique,pptau.extent(0));
    G00.resize(cix.Nbath_unique,dltau.extent(0));
    F00.resize(cix.Nbath_unique,dltau.extent(0));
    F01.resize(cix.Nbath_unique,dltau.extent(0));
  }
  int size(int Norder){ return 2*Norder+1;}

  double discrete_factor(){return 1.0;}
  
  void Broaden_Quantities(Array<double,3>& Stval, Array<double,3>& G00, Array<double,3>& F00, Array<double,3>& F01,
			  bool cmpPhysG, const json& pars, MPIcls& mpi){};

  void PrintStval_OCA(const Array<double,3>& Stval, const Array<double,2>& Stval_cOCA, int iorder2)
  {
    ofstream outSo("Stval_OCA.orig");
    ofstream outSm("Stval_OCA.mc");
    
    Array<double,2> dSmo, dSmm;
    Interpolate_SVD_OnNewMesh(dSmo,svdf,pptau,Stval(iorder2,Range::all(),Range::all()),beta);
    Interpolate_SVD_OnNewMesh(dSmm,svdf,pptau,Stval_cOCA,beta);
    
    for (int it=0; it<pptau.size(); it++){
      outSo << pptau(it) << " ";
      for (int ip=0; ip<cix.Npseudo_unique; ++ip) outSo << dSmo(ip,it)<<" ";
      outSo << endl;
      outSm << pptau(it) << " ";
      for (int ip=0; ip<cix.Npseudo_unique; ++ip) outSm << dSmm(ip,it)<<" ";
      outSm << endl;
    }
  }
  void PrintG00_OCA(Array<double,2>& dG00_cOCA, Array<double,2>& dF00_cOCA, Array<double,2>& dF01_cOCA, Array<double,3>& dG00, Array<double,3>& dF00, Array<double,3>& dF01, int iorder2)
  {
    ofstream gout("dG00_OCA.orig");
    ofstream gOut("dG00_OCA.mc");

    int Nt = dltau.extent(0);
    int Nbth = cix.Nbath_unique;

    Array<double,2> G00(Nbth,Nt); G00=0;
    Array<double,2> F00(Nbth,Nt); F00=0;
    Array<double,2> F01(Nbth,Nt); F01=0;
    // precompute svd functions on G00 mesh
    Array<double,2> u_t(Nt,svdf.lmax);
    for (int it=0; it<Nt; it++){
      intpar p = Interp( dltau(it),  svdf.tau);
      for (int l=0; l<svdf.lmax; l++) u_t(it,l) = svdf.fU[l](p);
    }

    Array<double,2> dG00t(Nbth,Nt), dF00t(Nbth,Nt), dF01t(Nbth,Nt), dG00o(Nbth,Nt), dF00o(Nbth,Nt), dF01o(Nbth,Nt);
    for (int ik=0; ik<Nbth; ++ik){
      for (int it=0; it<Nt; ++it){
	double dg=0, df00=0, df01=0, dgO=0, df00O=0, df01O=0;
	for (int l=0; l<svdf.lmax; ++l){
	  dg += dG00(iorder2,ik,l)*u_t(it,l);
	  df00 += dF00(iorder2,ik,l)*u_t(it,l);
	  df01 += dF01(iorder2,ik,l)*u_t(it,l);
	  dgO += dG00_cOCA(ik,l)*u_t(it,l);
	  df00O += dF00_cOCA(ik,l)*u_t(it,l);
	  df01O += dF01_cOCA(ik,l)*u_t(it,l);
	}
	dG00t(ik,it) = dg;
	dF00t(ik,it) = df00;
	dF01t(ik,it) = df01;
	dG00o(ik,it) = dgO;
	dF00o(ik,it) = df00O;
	dF01o(ik,it) = df01O;
      }
    }
    for (int it=0; it<Nt; it++){
      gout<<dltau(it)<<" ";
      for (int ik=0; ik<Nbth; ++ik) gout<< dG00t(ik,it)<<" ";
      for (int ik=0; ik<Nbth; ++ik) gout<< dF00t(ik,it)<<" ";
      for (int ik=0; ik<Nbth; ++ik) gout<< dF01t(ik,it)<<" "; 
      gout<<endl;
      
      gOut<<dltau(it)<<" ";
      for (int ik=0; ik<Nbth; ++ik) gOut<< dG00o(ik,it)<<" ";
      for (int ik=0; ik<Nbth; ++ik) gOut<< dF00o(ik,it)<<" ";
      for (int ik=0; ik<Nbth; ++ik) gOut<< dF01o(ik,it)<<" "; 
      gOut<<endl;
    }
    // REMOVE
    /*
    ofstream f00_out("dF00_OCA.mc");
    ofstream f01_out("dF01_OCA.mc");
    ofstream f02_out("dF02_OCA.mc");
    ofstream f03_out("dF03_OCA.mc");
    Array<double,2> dF02t(Nbth,Nt), dF03t(Nbth,Nt);
    for (int ik=0; ik<Nbth; ++ik){
      for (int it=0; it<Nt; ++it){
	double df00=0, df01=0, df02=0, df03=0;
	for (int l=0; l<svdf.lmax; ++l){
	  df00 += _F00_OCA(ik,l)*u_t(it,l);
	  df01 += _F01_OCA(ik,l)*u_t(it,l);
	  df02 += _F02_OCA(ik,l)*u_t(it,l);
	  df03 += _F03_OCA(ik,l)*u_t(it,l);
	}
	dF00t(ik,it) = df00;
	dF01t(ik,it) = df01;
	dF02t(ik,it) = df02;
	dF03t(ik,it) = df03;
      }
    }
    for (int it=0; it<Nt; it++){
      f00_out<<dltau(it)<<" ";
      f01_out<<dltau(it)<<" ";
      f02_out<<dltau(it)<<" ";
      f03_out<<dltau(it)<<" ";
      for (int ik=0; ik<Nbth; ++ik) f00_out << dF00t(ik,it)<<" ";
      for (int ik=0; ik<Nbth; ++ik) f01_out << dF01t(ik,it)<<" ";
      for (int ik=0; ik<Nbth; ++ik) f02_out << dF02t(ik,it)<<" "; 
      for (int ik=0; ik<Nbth; ++ik) f03_out << dF03t(ik,it)<<" ";
      f00_out <<endl;
      f01_out <<endl;
      f02_out <<endl;
      f03_out <<endl;
    }
    */
  }
  
  void InterpolateOnRightMesh(Array<double,2>& Sm, Array<double,2>& G00, Array<double,2>& F00, Array<double,2>& F01,
			      Array<double,3>& Stval, Array<double,3>& dG00l, Array<double,3>& dF00l, Array<double,3>& dF01l,
			      const Array<double,2>& Stval_cOCA, const Array<double,2>& G00_cOCA, const Array<double,2>& F00_cOCA, const Array<double,2>& F01_cOCA,
			      const Cix& cix, bool cmpPhysG, const json& pars, const MPIcls& mpi)
  {
    int Norders = Stval.extent(0);
    bool cOCA_needed = (OCA_cut_sum<10 and N_min_order==2);
    if (cOCA_needed){
      int iorder=0;
      // the lowest order is OCA, but there is missing part from OCA neglected in analytic calculation. We add it here.
      Stval(iorder,Range::all(),Range::all()) -= Stval_cOCA(Range::all(),Range::all());
      if (cmpPhysG){
	dG00l(iorder,Range::all(),Range::all()) -= G00_cOCA(Range::all(),Range::all());
	dF00l(iorder,Range::all(),Range::all()) -= F00_cOCA(Range::all(),Range::all());
	dF01l(iorder,Range::all(),Range::all()) -= F01_cOCA(Range::all(),Range::all());
      }
    }
    Sm.resize(cix.Npseudo_unique,pptau.extent(0)); Sm=0;
    Array<double,2> dSm;
    for (int iorder=0; iorder<Norders; ++iorder){
      Interpolate_SVD_OnNewMesh(dSm,svdf,pptau,Stval(iorder,Range::all(),Range::all()),beta);
      if (iorder>0 or (iorder==0 and cOCA_needed)) Sm += dSm; // notice that the lowest order is computed analytically, in order to normalize
      if (mpi.rank==mpi.master and pars["fastFilesystem"].get<bool>()){
	ofstream bout("Stval."s+to_string(iorder+N_min_order));
	for (int it=0; it<pptau.size(); it++){
	  bout<<pptau(it)<<" ";
	  for (int ip=0; ip<cix.Npseudo_unique; ++ip) bout<<dSm(ip,it)<<" ";
	  bout<<endl;
	}
      }
    }
    if (cmpPhysG){
      int Nt = dltau.extent(0);
      int Nbth = cix.Nbath_unique;
      G00.resize(Nbth,Nt); G00=0;
      F00.resize(Nbth,Nt); F00=0;
      F01.resize(Nbth,Nt); F01=0;
      Array<double,3> mg00_raw(Norders,Nbth,Nt), mf00_raw(Norders,Nbth,Nt), mf01_raw(Norders,Nbth,Nt);          // added by codex
      Array<double,3> mg00_smooth(Norders,Nbth,Nt), mf00_smooth(Norders,Nbth,Nt), mf01_smooth(Norders,Nbth,Nt); // added by codex
      mg00_raw = 0; mf00_raw = 0; mf01_raw = 0;
      mg00_smooth = 0; mf00_smooth = 0; mf01_smooth = 0;
      // precompute svd functions on G00 mesh
      Array<double,2> u_t(Nt,svdf.lmax);
      for (int it=0; it<Nt; it++){
	intpar p = Interp( dltau(it),  svdf.tau);
	for (int l=0; l<svdf.lmax; l++) u_t(it,l) = svdf.fU[l](p);
      }

      SmoothCorrectOrders smco(Nt,pars["S_sigma_pts"].get<double>());
      
      for (int ik=0; ik<Nbth; ++ik){
	Array<double,2> dG00t(Norders,Nt), dF00t(Norders,Nt), dF01t(Norders,Nt);	
	for (int iorder=0; iorder<Norders; ++iorder){
	  for (int it=0; it<Nt; ++it){
	    double dg=0, df00=0, df01=0;
	    for (int l=0; l<svdf.lmax; ++l){
	      dg += dG00l(iorder,ik,l)*u_t(it,l);
	      df00 += dF00l(iorder,ik,l)*u_t(it,l);
	      df01 += dF01l(iorder,ik,l)*u_t(it,l);
	    }
	    dG00t(iorder,it) = dg;
	    dF00t(iorder,it) = df00;
	    dF01t(iorder,it) = df01;
	  }
	}
	for (int iorder=0; iorder<Norders; ++iorder){
	  mg00_raw(iorder,ik,Range::all()) = dG00t(iorder,Range::all());
	  mf00_raw(iorder,ik,Range::all()) = dF00t(iorder,Range::all());
	  mf01_raw(iorder,ik,Range::all()) = dF01t(iorder,Range::all());
	}
	if (mpi.rank==mpi.master and pars["fastFilesystem"].get<bool>()){
	  ofstream gout("g00_."s+to_string(ik));
	  for (int it=0; it<Nt; it++){
	    gout<<dltau(it)<<" ";
	    for (int iorder=0; iorder<Norders; ++iorder) gout<< dG00t(iorder,it)<<" ";
	    for (int iorder=0; iorder<Norders; ++iorder) gout<< dF00t(iorder,it)<<" ";
	    for (int iorder=0; iorder<Norders; ++iorder) gout<< dF01t(iorder,it)<<" "; 
	    gout<<endl;
	  }
	}
	//
	for (int iorder=0; iorder<Norders; ++iorder){
	  if (iorder+1>=pars["S_min_order"].get<int>()){
	    Array<double,1> a(Nt);
	    a = smco.cmp(dG00t(iorder,Range::all()), dG00t(iorder-1,Range::all()), pars["S_ratio"].get<double>(), pars["S_kappa"].get<double>(), true, true);
	    double attenuation_ratio = sum(a)/a.extent(0);
	    if (attenuation_ratio<0.999) (*mpi.out)<<"<attenuation_ratio>[k="<<ik<<",order="<<iorder+1<<"]=" << attenuation_ratio << endl;
	    dG00t(iorder,Range::all()) *= a(Range::all());
	    dF00t(iorder,Range::all()) *= a(Range::all());
	    dF01t(iorder,Range::all()) *= a(Range::all()); 
	  }
	  //
	  if (iorder>0 or (iorder==0 and cOCA_needed)){  // normally we skip order==0 because it was determined analytically, however, when OCA is computed approximately, we do add the difference.
	    G00(ik,Range::all()) += dG00t(iorder,Range::all());
	    F00(ik,Range::all()) += dF00t(iorder,Range::all());
	    F01(ik,Range::all()) += dF01t(iorder,Range::all());
	  }
	}
	for (int iorder=0; iorder<Norders; ++iorder){
	  mg00_smooth(iorder,ik,Range::all()) = dG00t(iorder,Range::all());
	  mf00_smooth(iorder,ik,Range::all()) = dF00t(iorder,Range::all());
	  mf01_smooth(iorder,ik,Range::all()) = dF01t(iorder,Range::all());
	}
	if (mpi.rank==mpi.master and pars["fastFilesystem"].get<bool>()){
	  ofstream gout("g00."s+to_string(ik));
	  for (int it=0; it<Nt; it++){
	    gout<<dltau(it)<<" ";
	    for (int iorder=0; iorder<Norders; ++iorder) gout<< dG00t(iorder,it)<<" ";
	    for (int iorder=0; iorder<Norders; ++iorder) gout<< dF00t(iorder,it)<<" ";
	    for (int iorder=0; iorder<Norders; ++iorder) gout<< dF01t(iorder,it)<<" ";
	    gout<<endl;
	  }
	}
      }
      if (mpi.rank==mpi.master and pars["fastFilesystem"].get<bool>()){
	for (int iorder=0; iorder<Norders; ++iorder){
	  ofstream gout("mg00_."s+to_string(iorder+N_min_order));
	  for (int it=0; it<Nt; ++it){
	    gout<<dltau(it)<<" ";
	    for (int ik=0; ik<Nbth; ++ik) gout<< mg00_raw(iorder,ik,it)<<" ";
	    for (int ik=0; ik<Nbth; ++ik) gout<< mf00_raw(iorder,ik,it)<<" ";
	    for (int ik=0; ik<Nbth; ++ik) gout<< mf01_raw(iorder,ik,it)<<" ";
	    gout<<endl;
	  }
	}
      }
      if (mpi.rank==mpi.master and pars["fastFilesystem"].get<bool>()){
	for (int iorder=0; iorder<Norders; ++iorder){
	  ofstream gout("mg00."s+to_string(iorder+N_min_order));
	  for (int it=0; it<Nt; ++it){
	    gout<<dltau(it)<<" ";
	    for (int ik=0; ik<Nbth; ++ik) gout<< mg00_smooth(iorder,ik,it)<<" ";
	    for (int ik=0; ik<Nbth; ++ik) gout<< mf00_smooth(iorder,ik,it)<<" ";
	    for (int ik=0; ik<Nbth; ++ik) gout<< mf01_smooth(iorder,ik,it)<<" ";
	    gout<<endl;
	  }
	}
      }
    }
  }
  
  bool Read_sampled_data(Array<double,1>& Phi_accumulate, Array<double,1>& Phi_accumulatef,
			 Array<double,1>& N_measure, Array<double,1>& weight, Array<double,3>& Stval,
			 Array<double,3>& dG00, Array<double,3>& dF00, Array<double,3>& dF01,
			 double& Phi_cOCA, double& Phi_cOCAf, Array<double,2>& Stval_cOCA, Array<double,2>& G00_cOCA, Array<double,2>& F00_cOCA, Array<double,2>& F01_cOCA,
			 const string& sampled_file, MPIcls& mpi)
  {
    int Norders = Stval.extent(0);
    int Npseudo = Stval.extent(1);
    int _wNbin_lmax = Stval.extent(2);
    int Nbaths = dG00.extent(1);
    int _wNw_lmax = dG00.extent(2);
    bool previous_cmpPhysG=false, previous_cOCA=false, previous_cOCAf_exists=false;
    bool read_succeeded=true;
    if (mpi.rank==mpi.master){
      ifstream inp2(sampled_file);
      Phi_cOCA = 0;
      Phi_cOCAf = 0;
      if (Stval_cOCA.size()>0) Stval_cOCA = 0;
      if (G00_cOCA.size()>0) G00_cOCA = 0;
      if (F00_cOCA.size()>0) F00_cOCA = 0;
      if (F01_cOCA.size()>0) F01_cOCA = 0;

      int _iorder_;
      for (int iorder=0; iorder<Norders; ++iorder){
	inp2>>_iorder_>>Phi_accumulate(iorder)>>Phi_accumulatef(iorder)>>N_measure(iorder)>>weight(iorder);
	if (_iorder_!=iorder or (not inp2)){
	  (*mpi.out)<<"Note1 "<<sampled_file<<" was not properly read."<<endl;
	  read_succeeded=false;
	  break;
	}
      }
      if (not inp2 or inp2.eof()){
	if (read_succeeded) (*mpi.out)<<"Note1 "<<sampled_file<<" was not properly read."<<endl;
	read_succeeded=false;
      }
      if (read_succeeded){
	int lmax=0;
	inp2>>lmax;
	if (lmax!=_wNbin_lmax){
	  (*mpi.out)<<"Note2 "<<sampled_file<<" was not properly read."<<endl;
	  read_succeeded=false;
	}
      }
      if (read_succeeded){
	int _iorder2_, _atomic_state_;
	for (int iorder=0; iorder<Norders; ++iorder){
	  for (int atomic_state=0; atomic_state<Npseudo; ++atomic_state){
	    inp2>>_iorder2_>>_atomic_state_;
	    for (int l=0; l<_wNbin_lmax; ++l) inp2 >> Stval(iorder,atomic_state,l);
	    if (_iorder2_!=iorder or _atomic_state_!=atomic_state){
	      read_succeeded=false;
	      (*mpi.out)<<"Note3 "<<sampled_file<<" was not properly read."<<endl;
	      break;
	    }
	    if (not inp2 or inp2.eof()){
	      read_succeeded=false;
	      (*mpi.out)<<"Note4 "<<sampled_file<<" was not properly read."<<endl;
	      break;
	    }
	  }
	  if (not read_succeeded) break;
	}
      }
      if (read_succeeded){
	int _previous_cmpPhysG_=0;
	inp2>>_previous_cmpPhysG_;
	previous_cmpPhysG = static_cast<bool>(_previous_cmpPhysG_);
	if (previous_cmpPhysG){
	  int _iorder2_, _ibath_;
	  for (int iorder=0; iorder<Norders; ++iorder){
	    for (int ibath=0; ibath<Nbaths; ++ibath){
	      inp2 >> _iorder2_ >> _ibath_;
	      for (int l=0; l<_wNw_lmax; ++l)
		inp2 >> dG00(iorder,ibath,l) >> dF00(iorder,ibath,l) >> dF01(iorder,ibath,l);
	      if (_iorder2_!=iorder or _ibath_!=ibath){
		read_succeeded=false;
		(*mpi.out)<<"Note5 "<<sampled_file<<" was not properly read."<<endl;
		break;
	      }
	      if (not inp2 or inp2.eof()){
		read_succeeded=false;
		(*mpi.out)<<"Note6 "<<sampled_file<<" was not properly read."<<endl;
		break;
	      }
	    }
	  }
	}else{
	  dG00=0;
	  dF00=0;
	  dF01=0;
	  Phi_accumulatef=0;
	}
      }
      if (read_succeeded){
	int coca_marker=0;
	if (inp2 >> coca_marker){
	  if (coca_marker==0){
	    previous_cOCA = false;
	  }else if (coca_marker!=2){
	    read_succeeded=false;
	    (*mpi.out)<<"Note7 "<<sampled_file<<" was not properly read."<<endl;
	  }else{
	    previous_cOCA = true;
	    inp2 >> Phi_cOCA >> Phi_cOCAf;
	    if (not inp2){
	      read_succeeded=false;
	      (*mpi.out)<<"Note8 "<<sampled_file<<" was not properly read."<<endl;
	    }
	  }
	}else{
	  inp2.clear();
	}
      }
      if (read_succeeded and previous_cOCA){
	int Npseudo_cOCA = (Stval_cOCA.size()>0) ? Stval_cOCA.extent(0) : Npseudo;
	int lmax_cOCA = (Stval_cOCA.size()>0) ? Stval_cOCA.extent(1) : _wNbin_lmax;
	int _atomic_state_;
	double tmp_val;
	for (int atomic_state=0; atomic_state<Npseudo_cOCA; ++atomic_state){
	  inp2 >> _atomic_state_;
	  for (int l=0; l<lmax_cOCA; ++l){
	    inp2 >> tmp_val;
	    if (Stval_cOCA.size()>0) Stval_cOCA(atomic_state,l) = tmp_val;
	  }
	  if (_atomic_state_!=atomic_state or (not inp2)){
	    read_succeeded=false;
	    (*mpi.out)<<"Note9 "<<sampled_file<<" was not properly read."<<endl;
	    break;
	  }
	}
      }
      if (read_succeeded and previous_cOCA){
	int _previous_cOCAf_exists_=0;
	inp2 >> _previous_cOCAf_exists_;
	previous_cOCAf_exists = static_cast<bool>(_previous_cOCAf_exists_);
	if (previous_cOCAf_exists){
	  int Nbaths_cOCA = (G00_cOCA.size()>0) ? G00_cOCA.extent(0) : Nbaths;
	  int lmax_cOCA = (G00_cOCA.size()>0) ? G00_cOCA.extent(1) : _wNw_lmax;
	  int _ibath_;
	  double gtmp, f0tmp, f1tmp;
	  for (int ibath=0; ibath<Nbaths_cOCA; ++ibath){
	    inp2 >> _ibath_;
	    for (int l=0; l<lmax_cOCA; ++l){
	      inp2 >> gtmp >> f0tmp >> f1tmp;
	      if (G00_cOCA.size()>0) G00_cOCA(ibath,l) = gtmp;
	      if (F00_cOCA.size()>0) F00_cOCA(ibath,l) = f0tmp;
	      if (F01_cOCA.size()>0) F01_cOCA(ibath,l) = f1tmp;
	    }
	    if (_ibath_!=ibath or (not inp2)){
	      read_succeeded=false;
	      (*mpi.out)<<"Note10 "<<sampled_file<<" was not properly read."<<endl;
	      break;
	    }
	  }
	}
      }
    }

    mpi.bcast_two_boolean(read_succeeded,previous_cmpPhysG);
    mpi.bcast_two_boolean(previous_cOCA,previous_cOCAf_exists);
    (*mpi.out)<<"read_succeeded="<<read_succeeded<<" cmpPhysG="<<previous_cmpPhysG<<" cOCA="<<previous_cOCA<<" cOCAf="<<previous_cOCAf_exists<<endl;
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
      if (previous_cOCA){
	Array<double,1> Phi_cOCA_arr(2);
	if (mpi.rank==mpi.master){
	  Phi_cOCA_arr(0) = Phi_cOCA;
	  Phi_cOCA_arr(1) = Phi_cOCAf;
	}
	mpi.Bcast_simple_array(Phi_cOCA_arr);
	Phi_cOCA = Phi_cOCA_arr(0);
	Phi_cOCAf = Phi_cOCA_arr(1);
	if (Stval_cOCA.size()>0) mpi.Bcast_simple_array(Stval_cOCA);
	if (previous_cOCAf_exists){
	  if (G00_cOCA.size()>0) mpi.Bcast_simple_array(G00_cOCA);
	  if (F00_cOCA.size()>0) mpi.Bcast_simple_array(F00_cOCA);
	  if (F01_cOCA.size()>0) mpi.Bcast_simple_array(F01_cOCA);
	}else{
	  if (G00_cOCA.size()>0) G00_cOCA = 0;
	  if (F00_cOCA.size()>0) F00_cOCA = 0;
	  if (F01_cOCA.size()>0) F01_cOCA = 0;
	}
      }else{
	Phi_cOCA = 0;
	Phi_cOCAf = 0;
	if (Stval_cOCA.size()>0) Stval_cOCA = 0;
	if (G00_cOCA.size()>0) G00_cOCA = 0;
	if (F00_cOCA.size()>0) F00_cOCA = 0;
	if (F01_cOCA.size()>0) F01_cOCA = 0;
      }
    }else{
      Phi_accumulate=0;
      Phi_accumulatef=0;
      N_measure=0;
      Stval=0;
      dG00=0;
      dF00=0;
      dF01=0;
      Phi_cOCA = 0;
      Phi_cOCAf = 0;
      if (Stval_cOCA.size()>0) Stval_cOCA = 0;
      if (G00_cOCA.size()>0) G00_cOCA = 0;
      if (F00_cOCA.size()>0) F00_cOCA = 0;
      if (F01_cOCA.size()>0) F01_cOCA = 0;
      weight=1.0/weight.size();
    }
    return read_succeeded;
  }

  //double distribution_x_jump(){ return 1.0;}//1./beta;}
  void CheckQuantities(double cPhi, int ii, const Array<int,2>& diag, const Array<int,1>& dirspin,
		       const Array<double,1>& times, const Cix& cix, ostream& cout)
  {
  }

};
#endif
