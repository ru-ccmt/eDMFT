// Arbitrary-pole fourth-order reference evaluator.
//
// No atomic-path enumeration is used. Rotationally identical rooted cuts
// can share one representative, with exact orbit multiplicity. A switch
// retains the unreduced reference. All arrows and all cuts remain included.
// The numerical block operations are the SAME tested operations used by TCA.
// In particular, actual atomic blocks are retained until the final Sindex
// projection; equality of symmetry labels never identifies atomic paths.
#ifndef EDMFT_FCA_DLR_H
#define EDMFT_FCA_DLR_H
#include "fca_topology.h"
#include "tca_dlr.h"
#include <functional>
#include <climits>
#include <memory>
#include <limits>
#include <list>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <type_traits>
#include <set>

namespace edmft::fca_dlr {
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Engine=tca_dlr::SummedRootSelfEnergyEvaluator;
using Physical=tca_dlr::SummedPhysicalEvaluator;
using Function=Engine::SparseFunction;
using Workspace=imag_dlr::Workspace;

struct Options {
  bool reduce_cuts=true;
  bool physical_reuse=true;
  // Keep intermediate physical arcs in DLR coefficients. Constant sparse
  // F operators and bath sums act on coefficients directly; only the final
  // convolution returns node values for the G/F boundary contraction.
  // False is the previous evaluator. The basis also gates tight tolerances.
  bool physical_coeff_chain=true;
  // Look up the deepest reusable arc prefix BEFORE doing any convolution.
  // In particular, a cached bath-summed prefix can survive even when the
  // much larger branched intermediate before it did not fit the cache.
  bool physical_cache_resume=true;
  bool physical_cache_retain=true;
  // E/F sharing is independent of the existing q=4 A--D family. It uses
  // the existing physical_shared_bytes allowance, never a bath-cubed table.
  bool physical_ef_shared=true;
  bool physical_ef_presum=true; // combine closed E inputs and share E/F boundary weights
  bool physical_ef_boundary=false; // optional: small families may not amortize weight preparation
  bool self_reuse=true;
  // Retain ordinary self-energy chains in DLR coefficients between vertices.
  // The root and final boundary still use time-node values. A basis-level
  // stability gate falls back to the previous node path when needed.
  bool self_coeff_chain=true;
  bool self_adaptive_anchor=true;
  bool reuse_prepared=true;
  bool root_plans=true;
  bool indexed_states=true;
  bool convolution_scratch=true;
  bool self_pole_sum=false;
  bool self_direction_sum=true; // sum BOTH arrows at closure; no symmetry assumption
  int self_pole_sum_min_poles=3;
  bool shared_branch_transforms=true;
  // Experimental extension to shared families and physical arcs. Small
  // cluster A/B tests did not establish a speedup, so leave production's
  // previous paths intact. The older ordinary-self optimization above
  // remains enabled independently.
  bool extended_branch_transforms=false;
  // Extra retained dressed-propagator storage per MPI rank. Zero keeps the
  // previous task-local preparation; this is NOT a total solver-memory cap.
  std::size_t prepared_cache_bytes=128ULL*1024*1024;
  // Experimental family evaluation. False is the previous production path.
  bool shared_prefix=false;
  std::size_t shared_prefix_bytes=256ULL*1024*1024;
  bool physical_shared=false; // opt-in cross-diagram physical arc sharing
  std::size_t physical_shared_bytes=128ULL*1024*1024;
  bool shared_boundary=true;
  // Sum compatible, multiplicity-weighted B/C right arcs BEFORE their
  // common final convolution. Only the topology-total API may use this:
  // individual cut diagnostics must retain the separate member evaluation.
  bool physical_final_sum=true;
  std::size_t boundary_cache_bytes=64ULL*1024*1024;
  std::size_t self_prefix_bytes=64ULL*1024*1024;
  std::size_t physical_cache_bytes=64ULL*1024*1024;
  bool profile=false;
  // Stream per-task cost and gross contribution norms to rank zero while
  // the calculation is running. This is diagnostic only: no diagrams are cut.
  bool live_profile=false;
  // Opt-in bath-resolved screening pilot for ordinary FCA pseudo-self-energy
  // chains. It measures the streamed bath before that bath is summed into the
  // common suffix; no contribution is removed from the calculation.
  bool bath_sector_profile=false;
  // Diagnostic only: inspect one shared-self task without changing its result.
  // The selected worker writes singular spectra at vertices 2, 3, and 4.
  int tt_probe_task=-1;
  int tt_probe_bath=-1; // -1 selects the first nonzero streamed bath
  int tt_probe_pole=-1; // -1 selects the first nonzero streamed pole
};
struct PhysicalProfile {
  double convolution_functions=0,convolution_batches=0,cache_hits=0;
  double shared_branch_inputs=0;
  double resumed_vertices=0,cache_evictions=0;
  double cross_diagram_hits=0;
  double final_convolutions_saved=0,boundary_hits=0,boundary_builds=0;
  std::size_t cross_diagram_peak_bytes=0;
  std::size_t cache_peak_bytes=0;
};
struct SelfProfile {
  struct BathSector {
    double visits=0,seconds=0,prescore=0,intermediate_norm=0;
  };
  double convolution_functions=0,convolution_batches=0,root_functions=0;
  double prepared_unique=0,prepared_reused=0,prefix_fallbacks=0;
  int anchor=-1;
  std::size_t prefix_bytes=0;
  double root_seconds=0,vertex_seconds=0,convolution_seconds=0;
  double pack_seconds=0,transform_seconds=0,product_seconds=0,unpack_seconds=0;
  double pole_groups=0,shared_branch_inputs=0,shared_transform_seconds=0;
  std::vector<BathSector> bath_sectors;
  // Flat [bath][source block][target block] scores, allocated only for the
  // opt-in pilot and only when its dense table fits the fixed safety cap.
  std::vector<double> bath_block_prefix,bath_block_intermediate;
};

struct PreparedCacheProfile {
  std::uint64_t hits=0,misses=0,evictions=0,bypasses=0;
  std::size_t reserved_bytes=0,peak_bytes=0,entries=0;
};

struct SharedPrefixProfile {
  double prefix_convolutions=0,tail_convolutions=0,root_functions=0;
  double shared_branch_inputs=0;
  std::size_t peak_bytes=0;
  int fallbacks=0;
  std::vector<SelfProfile::BathSector> bath_sectors;
  std::vector<double> bath_block_prefix,bath_block_intermediate;
};

// The SC cluster uses 12 baths and 30 atomic blocks: 10800 entries per
// family, < 0.2 MiB for both scores. Large f-shell spaces skip this optional
// block breakdown rather than multiplying that storage by hundreds of ranks.
inline constexpr std::size_t bath_block_score_limit=100000;
inline std::size_t bath_block_entries(int baths,int blocks) {
  if (baths<=0 || blocks<=0) return 0;
  const long double count=static_cast<long double>(baths)*blocks*blocks;
  return count<=bath_block_score_limit ? static_cast<std::size_t>(count) : 0;
}
struct PhysicalSharedProfile {
  double convolution_functions=0,cache_hits=0,left_reuses=0,cache_bypasses=0,cache_evictions=0;
  double shared_branch_inputs=0,branch_transform_hits=0;
  std::size_t transform_peak_bytes=0;
  double boundary_builds=0,boundary_hits=0,boundary_fallbacks=0;
  double boundary_prepare_seconds=0,boundary_contract_seconds=0,boundary_reference_seconds=0;
  double final_sum_groups=0,final_convolutions_saved=0;
  std::size_t boundary_peak_bytes=0;
  std::size_t peak_bytes=0;
  int fallbacks=0;
};

inline fca::Cut physical_family_cut(const fca::CutOrbit& orbit) {
  std::array<int,4> arrows{};
  arrows[orbit.representative/2]=orbit.representative%2;
  auto cut=fca::physical_cut(orbit.topology,orbit.representative/2,arrows);
  std::sort(cut.split.begin(),cut.split.end(),[](const auto& a,const auto& b){ return a.a<b.a; });
  return cut;
}
inline bool physical_family_member(const fca::CutOrbit& orbit) {
  const auto cut=physical_family_cut(orbit);
  if (cut.q!=4) return false;
  for (int j=0;j<3;++j) if (cut.split[j].a!=j+1 || cut.split[j].b<5) return false;
  return true;
}

inline std::pair<int,unsigned> physical_ef_group(const fca::CutOrbit& orbit) {
  const auto cut=physical_family_cut(orbit);
  unsigned mask=0;
  for (int j=0;j<3;++j) if (cut.split[j].a<cut.q && cut.split[j].b>cut.q) mask|=1U<<j;
  if (orbit.topology<4 || (cut.q!=3 && cut.q!=5)) return {-1,0};
  return {cut.q,mask};
}

// Canonical labels mean FIRST, SECOND, THIRD opening, not the original
// hybridization numbering. This permutation is essential when sharing cuts
// from different topologies. Directions and poles passed to the family
// evaluator use these canonical labels; fallback maps them back explicitly.
inline fca::Cut prefix_family_cut(const fca::CutOrbit& orbit) {
  auto cut=fca::self_cut(orbit.topology,orbit.representative,{{0,0,0,0}});
  std::sort(cut.split.begin(),cut.split.end(),[](const auto& a,const auto& b){ return a.a<b.a; });
  return cut;
}
inline bool prefix_family_member(const fca::CutOrbit& orbit) {
  const auto cut=prefix_family_cut(orbit);
  if (cut.q!=4) return false;
  for (int j=0;j<3;++j) if (cut.split[j].a!=j+1 || cut.split[j].b<5) return false;
  return true;
}

// Direct-address storage for the at-most-two VARIABLE bath labels of one
// chain stage. The streamed label is fixed and is checked, not indexed.
// Matrices exist only for active entries; the address table is at most Nb^2.
// Each stage has one common open-label mask. A mask mismatch is an error,
// never an invitation to allocate an Nb^3 table or merge different baths.
class IndexedStates {
  using State=std::array<int,3>;
  int nb,excluded,mask=-1,fixed=-1;
  std::vector<int> slot;
  std::vector<std::pair<State,Function>> active;
public:
  IndexedStates(int baths,int streamed=-1):nb(baths),excluded(streamed) {}
  auto begin() { return active.begin(); }
  auto end() { return active.end(); }
  auto begin() const { return active.begin(); }
  auto end() const { return active.end(); }
  bool empty() const { return active.empty(); }
  std::size_t size() const { return active.size(); }
  std::size_t count(const State& key) const {
    for (const auto& entry:active) if (entry.first==key) return 1;
    return 0;
  }
  const Function& at(const State& key) const {
    // Only used for the one-entry closed physical boundary, not in hot sums.
    for (const auto& entry:active) if (entry.first==key) return entry.second;
    throw std::out_of_range("FCA missing closed state");
  }
  std::size_t index_bytes() const { return slot.capacity()*sizeof(int)+
    active.capacity()*sizeof(std::pair<State,Function>); }
  void clear() { active.clear(); mask=-1; slot.clear(); }
  Function& operator[](const State& key) {
    int current=0,count=0;
    std::size_t index=0,capacity=1;
    for (int j=0;j<3;++j) if (j!=excluded && key[j]>=0) {
      if (key[j]>=nb) throw std::logic_error("FCA indexed bath out of range");
      current|=1<<j; ++count;
      index=index*nb+key[j]; capacity*=nb;
    }
    if (count>2) throw std::logic_error("FCA indexed states would be cubic");
    if (mask<0) {
      mask=current; fixed=excluded>=0 ? key[excluded] : -1;
      slot.assign(capacity,-1); active.reserve(capacity);
    }
    if (mask!=current || (excluded>=0 && key[excluded]!=fixed))
      throw std::logic_error("FCA indexed state layout mismatch");
    int& position=slot[index];
    if (position<0) { position=active.size(); active.emplace_back(key,Function{}); }
    return active[position].second;
  }
  void emplace(const State& key,Function value) { (*this)[key]=std::move(value); }
};

// A diagnostic Schmidt spectrum of the open-bath histories against EVERYTHING
// carried by their matrix-valued pseudo-propagator functions.  This does not
// truncate or otherwise change the FCA evaluator.  Form the small row Gram
// matrix instead of materializing a potentially enormous dense history-by-
// (block,time,matrix-entry) tensor.  In the SC cluster a fixed streamed bath
// leaves at most Nbath^2=144 rows, so this is bounded even when the functions
// themselves contain many atomic blocks.
//
// The metric is physical L2 imaginary time, NOT Euclidean DLR coefficients:
// the dyadic quadrature resolves the sharp endpoint kernels, and the cached
// DLR node solve maps the node-valued shared-family functions into that
// quadrature.  Every missing sparse block is treated as an exact zero.
template<class StateMap>
void write_tt_probe_spectrum(std::ostream& out,const StateMap& states,
                             const dlr::FermionicBasis& basis,
                             int vertex,int arrow,int pole,int bath,
                             double stage_seconds)
{
  const auto begin=std::chrono::steady_clock::now();
  std::vector<const Function*> rows;
  std::set<std::pair<int,int>> blocks;
  for (const auto& [state,function]:states) if (!function.empty()) {
    rows.push_back(&function);
    for (const auto& [key,values]:function) if (!values.empty()) blocks.insert(key);
  }
  const int n=static_cast<int>(rows.size()),r=basis.rank();
  if (n==0) {
    out<<vertex<<','<<arrow<<','<<pole<<','<<bath
       <<",0,0,0,0,0,0,0,"<<stage_seconds<<",0\n";
    out.flush();
    return;
  }

  // Reuse the same endpoint-resolving quadrature construction as the DLR
  // basis.  sqrt(beta) converts its dimensionless weights to d(tau).
  const auto fine=dlr::detail::make_fine_grid(basis.lambda(),24);
  MatrixXd evaluation(fine.time_relative.size(),r);
  for (int q=0;q<evaluation.rows();++q)
    for (int l=0;l<r;++l)
      evaluation(q,l)=std::sqrt(basis.beta())*fine.sqrt_time_weights(q)*
        dlr::kernel_relative(fine.time_relative(q),basis.dimensionless_frequencies()(l));
  // self_prefix_family's states are values at the representative DLR nodes.
  const MatrixXd identity=MatrixXd::Identity(r,r);
  const MatrixXd inverse_nodes=basis.values_to_coefficients(identity);
  const MatrixXd physical=evaluation*inverse_nodes;
  const MatrixXd metric=physical.transpose()*physical;
  Eigen::SelfAdjointEigenSolver<MatrixXd> metric_eig(metric);
  if (metric_eig.info()!=Eigen::Success) throw std::runtime_error("FCA TT probe time metric eigensolve failed");
  MatrixXd metric_sqrt=metric_eig.eigenvectors().transpose();
  for (int l=0;l<r;++l)
    metric_sqrt.row(l)*=std::sqrt(std::max(0.0,metric_eig.eigenvalues()(l)));

  MatrixXd gram=MatrixXd::Zero(n,n);
  // Tile matrix elements so neither a dense atomic-space tensor nor a
  // second full copy of all sparse functions is ever allocated.
  constexpr int tile_size=8;
  for (const auto& key:blocks) {
    int nr=-1,nc=-1;
    for (const Function* row:rows) {
      const auto it=row->find(key);
      if (it==row->end()) continue;
      if (static_cast<int>(it->second.size())!=r)
        throw std::logic_error("FCA TT probe function has wrong DLR rank");
      nr=it->second[0].rows(); nc=it->second[0].cols(); break;
    }
    if (nr<0 || nc<0) continue;
    for (int first=0;first<nr*nc;first+=tile_size) {
      const int tile=std::min(tile_size,nr*nc-first);
      MatrixXd transformed=MatrixXd::Zero(n,r*tile);
      for (int u=0;u<n;++u) {
        const auto it=rows[u]->find(key);
        if (it==rows[u]->end()) continue;
        const auto& values=it->second;
        MatrixXd flat(r,tile);
        for (int t=0;t<r;++t) {
          if (values[t].rows()!=nr || values[t].cols()!=nc)
            throw std::logic_error("FCA TT probe sparse block shape mismatch");
          for (int e=0;e<tile;++e)
            flat(t,e)=values[t]((first+e)/nc,(first+e)%nc);
        }
        const MatrixXd weighted=metric_sqrt*flat;
        for (int t=0;t<r;++t)
          for (int e=0;e<tile;++e) transformed(u,t*tile+e)=weighted(t,e);
      }
      gram.noalias()+=transformed*transformed.transpose();
    }
  }
  Eigen::SelfAdjointEigenSolver<MatrixXd> eig(gram);
  if (eig.info()!=Eigen::Success) throw std::runtime_error("FCA TT probe Gram eigensolve failed");
  // The row-normalized spectrum reveals whether a small raw rank comes from
  // genuine shared histories or merely from a few much larger amplitudes.
  MatrixXd unit=gram;
  for (int u=0;u<n;++u) for (int v=0;v<n;++v) {
    const double denom=std::sqrt(std::max(0.0,gram(u,u))*std::max(0.0,gram(v,v)));
    unit(u,v)=denom>0.0 ? gram(u,v)/denom : 0.0;
  }
  Eigen::SelfAdjointEigenSolver<MatrixXd> unit_eig(unit);
  if (unit_eig.info()!=Eigen::Success) throw std::runtime_error("FCA TT probe unit-Gram eigensolve failed");
  const double total=std::max(0.0,eig.eigenvalues().sum());
  const double unit_total=std::max(0.0,unit_eig.eigenvalues().sum());
  double remainder=total,unit_remainder=unit_total;
  const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count();
  for (int j=0;j<n;++j) {
    const double value=std::max(0.0,eig.eigenvalues()(n-1-j));
    const double unit_value=std::max(0.0,unit_eig.eigenvalues()(n-1-j));
    remainder=std::max(0.0,remainder-value);
    unit_remainder=std::max(0.0,unit_remainder-unit_value);
    out<<vertex<<','<<arrow<<','<<pole<<','<<bath<<','<<n<<','<<blocks.size()<<','<<(j+1)<<','
       <<std::sqrt(value)<<','<<(total>0.0 ? std::sqrt(remainder/total) : 0.0)<<','
       <<std::sqrt(unit_value)<<','<<(unit_total>0.0 ? std::sqrt(unit_remainder/unit_total) : 0.0)<<','
       <<stage_seconds<<','<<seconds<<'\n';
  }
  out.flush(); // useful even if a long SC run is stopped after the probe
}

// Lifetime: one E/F family and ONE pair of crossing-bath labels. The key is
// the complete numeric/operation history of a prefix, not a topology hash.
// Values keep full atomic blocks. Canonical opening order fixes the meaning
// of each open-bath State slot across members. Matrix storage is immutable.
template<class StateMap> struct PhysicalArcCache {
  using Key=std::vector<std::uint64_t>;
  std::map<Key,std::shared_ptr<const StateMap>> entries;
  std::size_t bytes=0;
};

// Ordinary per-member cache, optionally kept alive while a family streams
// crossing-bath tuples. The family partitions its combined cache allowance
// between ordinary and cross-diagram storage, never Nmember full allowances.
template<class StateMap> struct PhysicalOrdinaryCache {
  using Key=std::array<int,4>;
  using Ptr=std::shared_ptr<const StateMap>;
  std::array<std::map<Key,Ptr>,2> entries;
  std::array<std::map<Key,std::pair<std::size_t,std::uint64_t>>,2> metadata;
  std::size_t bytes=0;
  std::uint64_t clock=0;
};

struct SectionTimer {
  double* sum;
  std::chrono::steady_clock::time_point start;
  explicit SectionTimer(double* output):sum(output) {
    if (sum) start=std::chrono::steady_clock::now();
  }
  ~SectionTimer() { if (sum) *sum+=std::chrono::duration<double>(
    std::chrono::steady_clock::now()-start).count(); }
};

struct Decorations {
  std::array<VectorXd,8> vertex,edge;
  double inverse_denominator=1.0;
};

// These scalar identities are independent of diagram order. The denominator
// ALWAYS uses a nonnegative effective frequency, so |K(0)| >= 1/2. Negative
// physical times use the solver's beta-t convention, not an extra guessed
// antiperiodic sign. The one -R per hybridization is applied at bath closure.
inline Decorations decorate(const Workspace& w,const fca::Cut& cut,
                            const std::array<int,3>& poles,bool physical)
{
  Decorations d;
  for (auto& v:d.vertex) v=VectorXd::Ones(w.rank());
  for (auto& v:d.edge) v=VectorXd::Ones(w.rank());
  for (int j=0;j<3;++j) {
    const auto& l=cut.split[j];
    const int p=poles[j],s=l.d==0 ? 1 : -1;
    if (p<0 || p>=w.pole_count()) throw std::out_of_range("FCA pole index");
    const double effective=s*w.poles()(p);
    const auto kernel=[&](int sign)->VectorXd {
      return sign>0 ? VectorXd(w.kernel_plus().col(p))
                    : VectorXd(w.kernel_minus().col(p));
    };
    const auto zero=[&](int sign) {
      return tca_dlr::KernelZero(w,p,sign);
    };
    if (!physical && effective<=0.0) {
      d.vertex[l.a].array()*=kernel(-s).array();
      d.vertex[l.b].array()*=kernel(s).array();
      d.inverse_denominator/=zero(-s);
    } else {
      // For physical cuts the complementary arc has duration beta-(tb-ta).
      // It replaces negative-frequency endpoint factors by positive-frequency
      // edge factors, allowing two independent convolution chains.
      const bool complement=physical && effective<0.0;
      const int frequency_sign=complement ? -s : s;
      int count=0;
      for (int e=0;e<(physical ? 8 : 7);++e) {
        const bool inside=l.a<=e && e<l.b;
        if (inside==complement) continue;
        d.edge[e].array()*=kernel(frequency_sign).array();
        ++count;
      }
      d.inverse_denominator/=std::pow(zero(frequency_sign),count-1);
    }
  }
  return d;
}

struct Result {
  MatrixXd sigma,g,f00,f01;
  std::array<MatrixXd,6> sigma_topology;
  std::array<std::array<MatrixXd,8>,6> sigma_cut;
  std::array<Physical::Result,6> physical_topology;
};

class Evaluator {
  const Cix& cix;
  const Workspace& w;
  Engine::PseudoBlocks pseudos;
  std::vector<Engine::OperatorTranspose> operators,bulla;
  std::vector<Engine::StoredRootPair> roots;
  std::vector<double> degeneracy;
  VectorXd ones;
  std::vector<VectorXd> root_weights;
  struct RootConnection { int root,left,right; Engine::BlockKey output; };
  // Connections only: no dense superoperator, no scalar pseudo-path catalog.
  // Populated lazily for encountered block pairs and destroyed with evaluator.
  mutable std::map<Engine::BlockKey,std::vector<RootConnection>> root_plan;
  mutable std::size_t root_plan_bytes=0;
  static constexpr std::size_t root_plan_limit=16ULL*1024*1024;

  // Geometry-only inverse routing for Tr[C L A R]. For a right block R(c,d),
  // C must have source d and A target c. The remaining block is L(C.target,
  // A.source). These references point into this evaluator's immutable
  // operators; no sampled values survive a bath tuple or Dyson iteration.
  struct BoundaryConnection {
    int bath,creation,annihilation,bulla;
    Engine::BlockKey left;
  };
  mutable std::map<Engine::BlockKey,std::vector<BoundaryConnection>> boundary_plan;
  mutable std::size_t boundary_plan_bytes=0;
  static constexpr std::size_t boundary_plan_limit=16ULL*1024*1024;

  struct BoundaryChannel {
    int channel=0;
    std::vector<MatrixXd> g,f00,f01;
  };
  using BoundaryBlock=std::vector<BoundaryChannel>;
  using BoundaryCache=std::map<Engine::BlockKey,BoundaryBlock>;

  const std::vector<BoundaryConnection>& boundary_connections(
      const Engine::BlockKey& key,std::vector<BoundaryConnection>& temporary) const {
    const auto found=boundary_plan.find(key);
    if (found!=boundary_plan.end()) return found->second;
    for (int bath=0;bath<cix.Nbath;++bath) {
      const auto& creation=operators[2*bath];
      const auto& annihilation=operators[2*bath+1];
      if (key.second<0 || key.second>=static_cast<int>(creation.by_source.size()) ||
          key.first<0 || key.first>=static_cast<int>(annihilation.by_target.size())) continue;
      const int ci=creation.by_source[key.second];
      if (ci<0) continue;
      const auto& ce=creation.entries[ci];
      const int bi=bulla[bath].by_source[ce.source];
      if (bi>=0 && bulla[bath].entries[bi].target!=ce.target)
        throw std::logic_error("FCA boundary Bulla connector changed target block");
      for (int ai:annihilation.by_target[key.first])
        temporary.push_back({bath,ci,ai,bi,{ce.target,annihilation.entries[ai].source}});
    }
    const std::size_t bytes=sizeof(key)+sizeof(temporary)+128+
      temporary.capacity()*sizeof(BoundaryConnection);
    if (bytes<=boundary_plan_limit-boundary_plan_bytes) {
      boundary_plan_bytes+=bytes;
      return boundary_plan.emplace(key,std::move(temporary)).first->second;
    }
    return temporary;
  }

  // Store (C L A)^T, so joining with R is a contiguous Frobenius product.
  // The explicit exact nonzero lists work for BOTH ordinary signed maps and
  // general sparse Bulla connectors. Do not use dense GEMM for Bulla here:
  // sparse scatter is O(nnz(C)*nnz(A)), not a new cubic matrix multiplication.
  static void add_boundary_weight(MatrixXd& destination,const Engine::OperatorEntry& creation,
      const MatrixXd& left,const Engine::OperatorEntry& annihilation,double sign) {
    for (const auto& c:creation.nonzero) for (const auto& a:annihilation.nonzero)
      destination(a.column,c.row)+=sign*c.value*a.value*left(c.column,a.row);
  }

  bool prepare_boundary_block(const Function& left,const Engine::BlockKey& key,
      int rows,int cols,std::size_t available,BoundaryBlock& result,std::size_t& bytes,
      bool fixed_right=false) const {
    std::vector<BoundaryConnection> temporary;
    // Mirrored q=5 E/F families share R instead of L. In that case prepare
    // (A R C)^T for contraction with L. Do NOT transpose the physical output,
    // change the channel sign or swap the Bulla operator: only reassociate
    // the same trace C L A R. The connector's `left` field denotes the fixed
    // arc's block in either mode.
    if (fixed_right) for (int bath=0;bath<cix.Nbath;++bath) {
      const auto& creation=operators[2*bath];
      const auto& annihilation=operators[2*bath+1];
      if (key.first<0 || key.first>=int(creation.by_target.size()) ||
          key.second<0 || key.second>=int(annihilation.by_source.size())) continue;
      const int ai=annihilation.by_source[key.second];
      if (ai<0) continue;
      for (int ci:creation.by_target[key.first]) {
        const auto& ce=creation.entries[ci];
        const int bi=bulla[bath].by_source[ce.source];
        if (bi>=0 && bulla[bath].entries[bi].target!=ce.target)
          throw std::logic_error("FCA mirrored Bulla connector changed target block");
        temporary.push_back({bath,ci,ai,bi,{annihilation.entries[ai].target,ce.source}});
      }
    }
    const auto& plan=fixed_right ? temporary : boundary_connections(key,temporary);
    // Determine required channels before allocating matrix data. Different
    // external baths with the same physical output channel are summed inside
    // these weights, including their individual channel signs below.
    std::map<int,int> masks;
    for (const auto& link:plan) if (left.count(link.left)) {
      const int channel=std::abs(cix.bath_sigind0(link.bath));
      masks[channel]|=1;
      if (link.bulla>=0) masks[channel]|=Physical::IsF01Channel(cix,link.bath) ? 4 : 2;
    }
    long double estimate=sizeof(BoundaryBlock)+sizeof(key)+128;
    for (const auto& [channel,mask]:masks) {
      const int count=1+bool(mask&2)+bool(mask&4);
      estimate+=sizeof(BoundaryChannel)+64+
        count*w.rank()*(sizeof(MatrixXd)+64+sizeof(double)*static_cast<long double>(rows)*cols);
    }
    if (estimate>available) return false;
    bytes=static_cast<std::size_t>(estimate);
    result.reserve(masks.size());
    for (const auto& [channel,mask]:masks) {
      result.push_back({channel,{},{},{}});
      auto& entry=result.back();
      entry.g.assign(w.rank(),MatrixXd::Zero(rows,cols));
      if (mask&2) entry.f00.assign(w.rank(),MatrixXd::Zero(rows,cols));
      if (mask&4) entry.f01.assign(w.rank(),MatrixXd::Zero(rows,cols));
    }
    for (const auto& link:plan) {
      const auto found=left.find(link.left);
      if (found==left.end()) continue;
      const int signed_channel=cix.bath_sigind0(link.bath),channel=std::abs(signed_channel);
      auto& target=*std::find_if(result.begin(),result.end(),[&](const auto& e){return e.channel==channel;});
      const auto& ce=operators[2*link.bath].entries[link.creation];
      const auto& ae=operators[2*link.bath+1].entries[link.annihilation];
      if (!fixed_right && (ae.matrix.cols()!=rows || ce.matrix.rows()!=cols ||
          found->second[0].rows()!=ce.matrix.cols() || found->second[0].cols()!=ae.matrix.rows())
        ) throw std::logic_error("FCA prepared boundary dimension mismatch");
      if (fixed_right && (ce.matrix.cols()!=rows || ae.matrix.rows()!=cols ||
          found->second[0].rows()!=ae.matrix.cols() || found->second[0].cols()!=ce.matrix.rows()))
        throw std::logic_error("FCA mirrored boundary dimension mismatch");
      const double sign=tca_dlr::SignedIndexSign(signed_channel);
      const auto scatter=[&](MatrixXd& out,const Engine::OperatorEntry& c,const MatrixXd& fixed) {
        if (!fixed_right) add_boundary_weight(out,c,fixed,ae,sign);
        else for (const auto& cv:c.nonzero) for (const auto& av:ae.nonzero)
          out(cv.column,av.row)+=sign*cv.value*av.value*fixed(av.column,cv.row);
      };
      for (int node=0;node<w.rank();++node) {
        scatter(target.g[node],ce,found->second[node]);
        if (link.bulla>=0) {
          auto& estimator=Physical::IsF01Channel(cix,link.bath) ? target.f01 : target.f00;
          scatter(estimator[node],bulla[link.bath].entries[link.bulla],found->second[node]);
        }
      }
    }
    return true;
  }

  // Cache only blocks actually requested by a right arc. This avoids paying
  // for all geometrically possible closures when most right blocks vanish.
  // Admission completes BEFORE adding any output: on a budget failure the
  // caller can safely evaluate the whole member through the old path once.
  bool accumulate_shared_boundary(const Function& left,const Function& right,
      double factor,Physical::Result& output,BoundaryCache& cache,std::size_t& used,
      std::size_t limit,PhysicalSharedProfile* profile,bool fixed_right=false) const {
    if (!limit) return false;
    {
      SectionTimer timer(profile ? &profile->boundary_prepare_seconds : nullptr);
      for (const auto& [key,values]:right) {
        if (cache.count(key)) { if (profile) ++profile->boundary_hits; continue; }
        BoundaryBlock block; std::size_t bytes=0;
        if (!prepare_boundary_block(left,key,values[0].rows(),values[0].cols(),limit-used,block,bytes,fixed_right))
          return false;
        cache.emplace(key,std::move(block)); used+=bytes;
        if (profile) {
          ++profile->boundary_builds;
          profile->boundary_peak_bytes=std::max(profile->boundary_peak_bytes,used);
        }
      }
    }
    SectionTimer timer(profile ? &profile->boundary_contract_seconds : nullptr);
    for (const auto& [key,values]:right) for (const auto& entry:cache.at(key))
      for (int node=0;node<w.rank();++node) {
        output.g(node,entry.channel)+=factor*(entry.g[node].array()*values[node].array()).sum();
        if (!entry.f00.empty())
          output.f00(node,entry.channel)+=factor*(entry.f00[node].array()*values[node].array()).sum();
        if (!entry.f01.empty())
          output.f01(node,entry.channel)+=factor*(entry.f01[node].array()*values[node].array()).sum();
      }
    return true;
  }

  // One evaluator belongs to ONE immutable pseudo-propagator workspace and
  // lives through self and physical evaluation of that iteration only.
  // Never move this cache to a static/global object: Dyson changes G even
  // when all pole/kernel keys remain identical.
  struct PreparedCacheEntry {
    std::uint64_t hash;
    VectorXd scalar;
    std::shared_ptr<Engine::PreparedPseudoCoefficients> value;
    std::size_t reserved;
  };
  mutable std::list<PreparedCacheEntry> prepared_cache; // MRU first
  mutable PreparedCacheProfile prepared_cache_stats;

  // Reserve for FULL lazy block expansion at insertion time, not just the
  // compact unique-Sindex matrices currently allocated. Otherwise a small
  // apparent cache could silently grow far beyond its cap during convolution.
  // Include both coefficient/Hilbert arrays, sampled key, matrix/container
  // storage, and conservative allocator overhead. This bounds accounted
  // retained storage, not process RSS or uncached task-local temporaries.
  std::size_t prepared_entry_reservation() const {
    long double bytes=sizeof(PreparedCacheEntry)+sizeof(Engine::PreparedPseudoCoefficients)+512;
    bytes+=8.L*w.rank()+16.L*w.rank()*w.pseudo_count();
    bytes+=cix.Nblock*(2*sizeof(tca_dlr::MatrixFunction)+2+256);
    for (int b=0;b<cix.Nblock;++b) {
      const long double n=cix.Block[b].size();
      bytes+=2.L*w.rank()*(8.L*n*n+sizeof(MatrixXd)+64);
    }
    if (bytes>=std::numeric_limits<std::size_t>::max())
      return std::numeric_limits<std::size_t>::max();
    return static_cast<std::size_t>(bytes);
  }

  std::shared_ptr<Engine::PreparedPseudoCoefficients> prepared_edge(
      const VectorXd& scalar,std::size_t limit) const {
    // Evict only objects not held by an active chain. Removing a pinned
    // object from accounting would disguise still-live cache memory.
    const auto evict_one=[&]() {
      auto it=prepared_cache.end();
      while (it!=prepared_cache.begin()) {
        --it;
        if (it->value.use_count()!=1) continue;
        prepared_cache_stats.reserved_bytes-=it->reserved;
        prepared_cache.erase(it);
        ++prepared_cache_stats.evictions;
        prepared_cache_stats.entries=prepared_cache.size();
        return true;
      }
      return false;
    };
    while (prepared_cache_stats.reserved_bytes>limit && evict_one()) {}
    const auto build=[&]() {
      return std::make_shared<Engine::PreparedPseudoCoefficients>(
        Engine::PreparePseudoCoefficients(w,scalar));
    };
    if (!limit) { ++prepared_cache_stats.bypasses; return build(); }

    // The key is the EXACT sampled decoration, not topology, bath, or a
    // rounded frequency sum. Hash collisions are checked byte-for-byte;
    // kernels multiplied in a different floating-point order are not merged
    // unless their samples really agree. No new approximation is introduced.
    std::uint64_t hash=14695981039346656037ULL;
    const auto* raw=reinterpret_cast<const unsigned char*>(scalar.data());
    const std::size_t length=sizeof(double)*scalar.size();
    for (std::size_t i=0;i<length;++i) { hash^=raw[i]; hash*=1099511628211ULL; }
    for (auto it=prepared_cache.begin();it!=prepared_cache.end();++it)
      if (it->hash==hash && it->scalar.size()==scalar.size() &&
          std::memcmp(it->scalar.data(),scalar.data(),length)==0) {
        auto value=it->value;
        prepared_cache.splice(prepared_cache.begin(),prepared_cache,it);
        ++prepared_cache_stats.hits;
        return value;
      }
    ++prepared_cache_stats.misses;
    const std::size_t reservation=prepared_entry_reservation();
    if (reservation>limit) { ++prepared_cache_stats.bypasses; return build(); }
    while (prepared_cache_stats.reserved_bytes>limit-reservation)
      if (!evict_one()) { ++prepared_cache_stats.bypasses; return build(); }
    auto value=build();
    prepared_cache.push_front({hash,scalar,value,reservation});
    prepared_cache_stats.reserved_bytes+=reservation;
    prepared_cache_stats.peak_bytes=std::max(prepared_cache_stats.peak_bytes,
                                          prepared_cache_stats.reserved_bytes);
    prepared_cache_stats.entries=prepared_cache.size();
    return value;
  }

  double residue(int bath,int pole) const {
    const int channel=cix.bath_sigind0(bath);
    return -tca_dlr::SignedIndexSign(channel)*w.residues()(std::abs(channel),pole);
  }
  using State=std::array<int,3>;
  using States=std::map<State,Function>;

  static std::size_t function_bytes(const Function& function) {
    std::size_t bytes=sizeof(Function)+128;
    for (const auto& [key,matrices]:function) {
      bytes+=sizeof(key)+sizeof(matrices)+128;
      for (const auto& matrix:matrices)
        bytes+=sizeof(MatrixXd)+sizeof(double)*matrix.size()+64;
    }
    return bytes;
  }

  // Estimate real compatible operator branches, rather than nominal Nbath.
  // A transform of a large input is not a saving when most F_b blocks are
  // absent or when they select much smaller output blocks. The direction
  // here is the actual right-multiplying operator, including Opposite at
  // opening vertices. This test never merges physical blocks by Sindex.
  std::pair<double,double> branch_entries(const Function& input,int direction,int pole) const {
    double incoming=0,outgoing=0;
    for (const auto& [key,values]:input) {
      incoming+=values[0].size();
      for (int bath=0;bath<cix.Nbath;++bath) {
        if (residue(bath,pole)==0.) continue;
        const auto& op=operators[2*bath+direction];
        const int index=op.by_source[key.second];
        if (index>=0) outgoing+=double(values[0].rows())*op.entries[index].matrix.cols();
      }
    }
    return {incoming,outgoing};
  }

  template<class StateMap>
  static std::size_t storage_bytes(const StateMap& states) {
    std::size_t bytes=sizeof(StateMap)+128;
    if constexpr (std::is_same_v<StateMap,IndexedStates>) bytes+=states.index_bytes();
    for (const auto& [state,function]:states) {
      bytes+=sizeof(State)+sizeof(Function)+128;
      for (const auto& [key,matrices]:function) {
        bytes+=sizeof(key)+sizeof(matrices)+128;
        for (const auto& matrix:matrices)
          bytes+=sizeof(MatrixXd)+sizeof(double)*matrix.size()+64;
      }
    }
    return bytes;
  }

  // Batch all blocks of ONE common incoming function in a multiple-RHS
  // transform. A separate tiny solve for each atomic block is expensive in
  // the normal state. The caller's branch-memory guard also bounds these
  // three rank-by-entry buffers. No transformed prefix survives a task.
  std::pair<Function,Function> prepare_left(const Function& input,
                                          const VectorXd& scalar) const {
    Eigen::Index entries=0;
    for (const auto& [key,values]:input) entries+=values[0].size();
    MatrixXd flat(w.rank(),entries);
    Eigen::Index offset=0;
    for (const auto& [key,values]:input) {
      const int nr=values[0].rows(),nc=values[0].cols();
      for (int i=0;i<nr;++i) for (int j=0;j<nc;++j)
        for (int l=0;l<w.rank();++l)
          flat(l,offset+i*nc+j)=scalar(l)*values[l](i,j);
      offset+=nr*nc;
    }
    const MatrixXd cf=w.basis().values_to_coefficients(flat);
    const MatrixXd hf=w.basis().convolution_hilbert(cf);
    std::pair<Function,Function> result;
    offset=0;
    for (const auto& [key,values]:input) {
      const int nr=values[0].rows(),nc=values[0].cols();
      auto& first=result.first[key]; auto& second=result.second[key];
      first.resize(w.rank()); second.resize(w.rank());
      for (int l=0;l<w.rank();++l) {
        first[l].resize(nr,nc); second[l].resize(nr,nc);
        for (int i=0;i<nr;++i) for (int j=0;j<nc;++j) {
          first[l](i,j)=cf(l,offset+i*nc+j);
          second[l](i,j)=hf(l,offset+i*nc+j);
        }
      }
      offset+=nr*nc;
    }
    return result;
  }

  // Apply the same small DLR-time linear map to every atomic matrix entry.
  // This is a multiple-RHS GEMM, not one transform per block or bath state.
  // In coefficient-chain mode the map B^{-1} diag(vertex) B dresses a vertex
  // while staying in coefficients; B itself is used only at the root, whose
  // hybridization is sampled at the selected time nodes.
  Function transform_time(const Function& input,const MatrixXd& map) const {
    Eigen::Index entries=0;
    for (const auto& [key,values]:input) entries+=values[0].size();
    MatrixXd flat(w.rank(),entries);
    Eigen::Index offset=0;
    for (const auto& [key,values]:input) {
      const int nr=values[0].rows(),nc=values[0].cols();
      for (int i=0;i<nr;++i) for (int j=0;j<nc;++j)
        for (int l=0;l<w.rank();++l)
          flat(l,offset+i*nc+j)=values[l](i,j);
      offset+=nr*nc;
    }
    const MatrixXd transformed=map*flat;
    Function output;
    offset=0;
    for (const auto& [key,values]:input) {
      const int nr=values[0].rows(),nc=values[0].cols();
      auto& result=output[key]; result.resize(w.rank());
      for (int l=0;l<w.rank();++l) {
        result[l].resize(nr,nc);
        for (int i=0;i<nr;++i) for (int j=0;j<nc;++j)
          result[l](i,j)=transformed(l,offset+i*nc+j);
      }
      offset+=nr*nc;
    }
    return output;
  }

  // A prepared edge also owns lazy block caches. Sharing the object, rather
  // than copying its matrices, shares those caches as well. Exact equality of
  // the sampled decoration is sufficient; no approximate frequency grouping
  // is permitted. The bounded evaluator cache extends reuse across tasks,
  // cuts and self/physical phases, but cannot outlive a PP iteration.
  using Prepared=std::array<std::shared_ptr<Engine::PreparedPseudoCoefficients>,8>;
  bool coefficient_chain(const Options& options) const {
    return options.physical_reuse && options.physical_coeff_chain &&
      w.basis().coefficient_convolution_available();
  }
  Function physical_seed(const VectorXd& scalar,const Options& options,bool coefficients) const {
    if (!coefficients) return Engine::PseudoValues(pseudos,scalar);
    // Fit only the unique Sindex columns, then expand the signed atomic
    // blocks. Never perform a separate interpolation solve for each copy
    // of the same pseudo-propagator. This shares the dressed-edge cache.
    const auto prepared=prepared_edge(scalar,options.reuse_prepared ? options.prepared_cache_bytes : 0);
    // The undressed workspace can contain coefficients from an overdetermined
    // fit, whereas the old LEFT arc starts from its reconstructed node data.
    // Re-solve those unique columns once at the arc entrance so both routes
    // use the same interpolation convention even in an ill-conditioned basis.
    // Dressed inputs already went through this solve in prepared_edge().
    MatrixXd canonical;
    if (prepared->use_workspace_coefficients)
      canonical=w.basis().values_to_coefficients(w.pseudo_values());
    const auto& unique=prepared->use_workspace_coefficients ? canonical : prepared->dressed;
    Function result;
    for (int block=0;block<cix.Nblock;++block) {
      auto values=Engine::PseudoCoefficientBlock(cix,block,unique);
      bool nonzero=false;
      for (const auto& m:values) nonzero=nonzero || m.squaredNorm()!=0.;
      if (nonzero) result.emplace(Engine::BlockKey{block,block},std::move(values));
    }
    return result;
  }
  Prepared prepare_edges(const Decorations& d,int end,int skip,const Options& options,
                         SelfProfile* profile=nullptr) const {
    Prepared result{};
    for (int e=1;e<end;++e) if (e!=skip) {
      if (options.reuse_prepared) for (int previous=1;previous<e;++previous)
        if (result[previous] && (d.edge[e].array()==d.edge[previous].array()).all()) {
          result[e]=result[previous]; break;
        }
      if (result[e]) { if (profile) ++profile->prepared_reused; }
      else {
        result[e]=prepared_edge(d.edge[e],options.reuse_prepared ? options.prepared_cache_bytes : 0);
        if (profile) ++profile->prepared_unique;
      }
    }
    return result;
  }

  // One root contraction is shared by all split-bath assignments reaching
  // this stage. Both root arrows are present in the full-orientation domain.
  Function root_contract(const Function& h) const {
    return Engine::ContractRoot(h,roots,w,ones);
  }

  Function planned_root(const Function& h) const {
    Function result;
    for (const auto& [key,values]:h) {
      auto found=root_plan.find(key);
      std::vector<RootConnection> temporary;
      if (found==root_plan.end()) {
        for (int r=0;r<static_cast<int>(roots.size());++r) {
          const auto& root=roots[r];
          if (key.first<0 || key.first>=static_cast<int>(root.zero.by_target.size()) ||
              key.second<0 || key.second>=static_cast<int>(root.endpoint.by_source.size())) continue;
          const int right=root.endpoint.by_source[key.second];
          if (right<0) continue;
          for (int left:root.zero.by_target[key.first])
            temporary.push_back({r,left,right,{root.zero.entries[left].source,root.endpoint.entries[right].target}});
        }
        const std::size_t bytes=128+sizeof(key)+sizeof(temporary)+temporary.capacity()*sizeof(RootConnection);
        if (bytes<=root_plan_limit-root_plan_bytes) {
          root_plan_bytes+=bytes;
          found=root_plan.emplace(key,std::move(temporary)).first;
        }
      }
      const auto& connections=found==root_plan.end() ? temporary : found->second;
      for (const auto& link:connections) {
        const auto& left=roots[link.root].zero.entries[link.left];
        const auto& right=roots[link.root].endpoint.entries[link.right];
        auto [out,inserted]=result.try_emplace(link.output);
        if (inserted) {
          out->second.resize(values.size());
          for (auto& node:out->second) node=MatrixXd::Zero(left.matrix.rows(),right.matrix.cols());
        }
        for (std::size_t node=0;node<values.size();++node)
          Engine::AccumulateBothOperators(out->second[node],left,values[node],right,root_weights[link.root](node));
      }
    }
    return result;
  }

public:
  PreparedCacheProfile prepared_cache_profile() const { return prepared_cache_stats; }
  std::vector<fca::CutOrbit> select_prefix_family(const std::vector<fca::CutOrbit>& orbits,
                                                const Options& options) const {
    // Geometry-only admission estimate, followed by the actual live-byte
    // guard in the evaluator. A first-opening line closing at 5 retains a
    // quadratic continuation accumulator; closures at 6/7 retain only a
    // linear/scalar one. Excluding those two expensive consumers keeps the
    // SC family useful without committing several quadratic tensors.
    long double base=512;
    for (int b=0;b<cix.Nblock;++b) {
      const long double n=cix.Block[b].size();
      base+=256+w.rank()*(8*n*n+sizeof(MatrixXd)+64);
    }
    const long double nb=cix.Nbath;
    std::vector<fca::CutOrbit> family;
    if (base*(nb*nb+3*nb+4)>options.shared_prefix_bytes) return family;
    const bool include_early=base*(4*nb*nb+2*nb+4)<=options.shared_prefix_bytes;
    for (const auto& orbit:orbits)
      if (prefix_family_member(orbit) && (include_early || prefix_family_cut(orbit).split[0].b>=6))
        family.push_back(orbit);
    if (family.size()<2) family.clear();
    return family;
  }
  Evaluator(const Cix& c,const Workspace& workspace):cix(c),w(workspace),
      pseudos(Engine::PreparePseudoBlocks(c,w)),ones(VectorXd::Ones(w.rank()))
  {
    if (!w.ready() || !w.has_pseudos())
      throw std::invalid_argument("FCA requires hybridization and pseudo DLR data");
    degeneracy.resize(cix.Npseudo_unique);
    for (int p=0;p<cix.Npseudo_unique;++p) degeneracy[p]=cix.deg_pp(p);
    for (int op=0;op<2*cix.Nbath;++op)
      operators.push_back(Engine::MakeOperatorTranspose(cix,op));
    for (int bath=0;bath<cix.Nbath;++bath) {
      bulla.push_back(Physical::MakeBullaTranspose(cix,bath));
      const int channel=cix.bath_sigind0(bath);
      for (int direction=0;direction<2;++direction) {
        const int op=2*bath+direction;
        roots.push_back({operators[Engine::Opposite(op)],operators[op],
          std::abs(channel),direction,
          static_cast<double>(tca_dlr::SignedIndexSign(channel))});
      }
    }
    for (const auto& root:roots) {
      const auto& delta=root.direction==0 ? w.delta_plus() : w.delta_minus();
      root_weights.push_back((-root.bath_sign*delta.col(root.channel)).eval());
    }
  }

  // pole[] refers to the three NON-root lines in the rotated cut. The root
  // bath/direction is summed locally; its hybridization comes from precisely
  // this workspace (never the more accurate OCA workspace).
  MatrixXd self_nodes_reference(int topology,int cut_edge,int directions,
                      const std::array<int,3>& pole,SelfProfile* profile=nullptr) const
  {
    auto cut=fca::self_cut(topology,cut_edge,{{0,0,0,0}});
    for (int j=0;j<3;++j) cut.split[j].d=(directions>>j)&1;
    const Decorations d=decorate(w,cut,pole,false);
    std::array<Engine::PreparedPseudoCoefficients,7> prepared;
    for (int e=1;e<7;++e)
      prepared[e]=Engine::PreparePseudoCoefficients(w,d.edge[e]);
    if (profile) profile->prepared_unique+=6;
    MatrixXd output=MatrixXd::Zero(w.rank(),w.pseudo_count());

    // Stream the bath of the latest-closing split line. This bounds live
    // bath-state storage by Nbath^2 rather than Nbath^3, while coefficient
    // preparation stays shared across all anchor baths. Other labels are
    // summed immediately at their second endpoint. In topology E, a label
    // closes BEFORE the third one opens, so no cubic intermediate is formed.
    int anchor=0;
    for (int j=1;j<3;++j)
      if (cut.split[j].b>cut.split[anchor].b) anchor=j;
    if (profile) profile->anchor=anchor;
    for (int anchor_bath=0;anchor_bath<cix.Nbath;++anchor_bath) {
      States states;
      states[{{-1,-1,-1}}]=Engine::PseudoValues(pseudos,d.edge[0]);
      for (int vertex=1;vertex<8 && !states.empty();++vertex) {
        States next;
        if (vertex==cut.q) {
          for (auto& [state,h]:states) {
            if (profile) ++profile->root_functions;
            auto value=root_contract(h);
            if (!value.empty()) next.emplace(state,std::move(value));
          }
        } else {
          int j=0;
          while (j<3 && cut.split[j].a!=vertex && cut.split[j].b!=vertex) ++j;
          if (j==3) throw std::logic_error("FCA unmatched vertex");
          const auto& line=cut.split[j];
          const bool opening=line.a==vertex;
          for (auto& [state,h]:states) {
            const int first=opening ? (j==anchor ? anchor_bath : 0) : state[j];
            const int end=opening ? (j==anchor ? anchor_bath+1 : cix.Nbath) : first+1;
            for (int bath=first;bath<end;++bath) {
              if (bath<0) throw std::logic_error("FCA bath closed before opening");
              const double weight=opening ? 1.0 : residue(bath,pole[j]);
              if (weight==0.0) continue;
              const int op=2*bath+line.d;
              auto value=Engine::MultiplyRight(h,
                operators[opening ? Engine::Opposite(op) : op],
                (weight*d.vertex[vertex]).eval());
              if (value.empty()) continue;
              State target=state;
              target[j]=opening ? bath : -1;
              Engine::AddFunction(next[target],std::move(value));
            }
          }
        }
        states=std::move(next);
        if (vertex<7 && !states.empty()) {
          std::vector<State> keys;
          std::vector<Function> batch;
          for (auto& [key,h]:states) { keys.push_back(key); batch.push_back(std::move(h)); }
          if (profile) { profile->convolution_functions+=batch.size(); ++profile->convolution_batches; }
          Engine::ConvolveRightPseudoPreparedBatch(batch,cix,w,prepared[vertex]);
          states.clear();
          for (std::size_t k=0;k<keys.size();++k)
            if (!batch[k].empty()) states.emplace(keys[k],std::move(batch[k]));
        }
      }
      for (const auto& [state,h]:states) {
        if (state!=State{{-1,-1,-1}}) throw std::logic_error("FCA unclosed bath");
        Engine::AccumulateCut(h,cix,
          fca::symmetry[topology]*fca::crossing_sign[topology]*d.inverse_denominator,
          degeneracy,output);
      }
    }
    return output;
  }

  // Share everything outside the streamed line's two endpoints. The prefix
  // retains OTHER open bath labels. At closure we add only states with equal
  // remaining labels, then perform the closing-edge convolution and suffix
  // once. This is distributivity of linear operators, not diagram pruning.
  MatrixXd self_nodes(int topology,int cut_edge,int directions,
                      const std::array<int,3>& pole,const Options& options=Options{},
                      SelfProfile* profile=nullptr) const
  {
    if (options.indexed_states)
      return self_nodes_impl<IndexedStates>(topology,cut_edge,directions,pole,options,profile);
    return self_nodes_impl<States>(topology,cut_edge,directions,pole,options,profile);
  }

  // A task holds two poles fixed and streams the third pole TOGETHER with
  // its bath label. Outside that line's endpoints both labels are absent.
  // This retains the existing at-most-Nbath^2 bath storage; no pole tensor
  // is allocated. Choose the line using the same live-label work model as
  // bath streaming, now charging its interval once per pole.
  int pole_anchor(int topology,int cut_edge,int repeats=0) const {
    const auto cut=fca::self_cut(topology,cut_edge,{{0,0,0,0}});
    int best=0; double best_cost=std::numeric_limits<double>::infinity();
    for (int candidate=0;candidate<3;++candidate) {
      double cost=0;
      for (int v=1;v<8;++v) {
        int before=0,after=0;
        for (int j=0;j<3;++j) if (j!=candidate) {
          before+=cut.split[j].a<v && v<=cut.split[j].b;
          after+=cut.split[j].a<=v && v<cut.split[j].b;
        }
        const double nb=std::pow(double(cix.Nbath),before);
        const double na=std::pow(double(cix.Nbath),after);
        const bool inside=cut.split[candidate].a<=v && v<=cut.split[candidate].b;
        const double repeat=inside ? double(cix.Nbath)*(repeats ? repeats : w.pole_count()) : 1.;
        cost+=repeat*(v==cut.q ? 2*cix.Nbath*nb : std::max(nb,na));
        if (v<7) cost+=na*(inside && v<cut.split[candidate].b ? repeat : 1.);
      }
      if (cost<best_cost) { best=candidate; best_cost=cost; }
    }
    return best;
  }

  // fixed_directions contains just the TWO other arrows, in increasing
  // split-line order. Both arrows of anchor are explicitly evaluated, with
  // their own operators, endpoint factors and normalization, then added at
  // closure. This is not identifying rotated diagrams or assuming equality.
  MatrixXd self_nodes_direction_sum(int topology,int cut_edge,int fixed_directions,
      const std::array<int,3>& input_poles,bool sum_poles,const Options& options,
      SelfProfile* profile=nullptr) const {
    const int baseline_anchor=bath_anchor(topology,cut_edge,options);
    const int anchor=sum_poles ? pole_anchor(topology,cut_edge) :
      (baseline_anchor>=0 ? baseline_anchor : pole_anchor(topology,cut_edge,2));
    int directions=0; auto poles=input_poles;
    for (int j=0,k=0;j<3;++j) if (j!=anchor) {
      directions|=((fixed_directions>>k)&1)<<j;
      if (sum_poles) poles[j]=input_poles[k];
      ++k;
    }
    if (sum_poles) poles[anchor]=0;
    if (profile && sum_poles) ++profile->pole_groups;
    // Never change a previously unstreamed bath plan merely to group arrows.
    // That can increase convolution/batch counts in SC blocks. Both arrows
    // still belong to this task, but retain their previous evaluator there.
    if (!sum_poles && baseline_anchor<0) {
      MatrixXd result=self_nodes(topology,cut_edge,directions,poles,options,profile);
      result+=self_nodes(topology,cut_edge,directions^(1<<anchor),poles,options,profile);
      return result;
    }
    if (options.indexed_states)
      return self_nodes_impl<IndexedStates>(topology,cut_edge,directions,poles,options,profile,
                                            sum_poles ? anchor : -1,anchor);
    return self_nodes_impl<States>(topology,cut_edge,directions,poles,options,profile,
                                   sum_poles ? anchor : -1,anchor);
  }

  MatrixXd self_nodes_pole_sum(int topology,int cut_edge,int directions,
                      const std::array<int,2>& fixed,const Options& options=Options{},
                      SelfProfile* profile=nullptr) const {
    const int anchor=pole_anchor(topology,cut_edge);
    std::array<int,3> poles{};
    for (int j=0,k=0;j<3;++j) if (j!=anchor) poles[j]=fixed[k++];
    if (!options.self_reuse) {
      MatrixXd result=MatrixXd::Zero(w.rank(),w.pseudo_count());
      for (int l=0;l<w.pole_count();++l) {
        poles[anchor]=l;
        result+=self_nodes(topology,cut_edge,directions,poles,options,profile);
      }
      return result;
    }
    if (profile) ++profile->pole_groups;
    if (options.indexed_states)
      return self_nodes_impl<IndexedStates>(topology,cut_edge,directions,poles,options,profile,anchor);
    return self_nodes_impl<States>(topology,cut_edge,directions,poles,options,profile,anchor);
  }

  // Shared convolution prefix for the q=4, three-consecutive-openings
  // family. One task fixes the last two poles and sums the FIRST opening's
  // pole and bath. At most two bath labels vary in a live prefix: there is
  // no Nbath^3 stored tensor and no stored pole-triple tensor.
  //
  // All members have the SAME decorated edges/vertices through vertex 4.
  // Their inverse kernel normalizations and diagram signs are NOT shared.
  // Each continuation merges the streamed bath/pole AT ITS OWN CLOSURE,
  // before that closing-edge convolution. Remaining tails are done once,
  // preserving the earlier early-sum optimization rather than undoing it.
  std::vector<MatrixXd> self_prefix_family(const std::vector<fca::CutOrbit>& family,
                      int directions,const std::array<int,2>& fixed,
                      const Options& options,SharedPrefixProfile* profile=nullptr,
                      bool sum_directions=false,std::ostream* tt_probe=nullptr) const {
    if (family.empty()) return {};
    for (const auto& member:family)
      if (!prefix_family_member(member)) throw std::invalid_argument("incompatible FCA prefix family");
    const auto zeros=[&]() {
      return std::vector<MatrixXd>(family.size(),MatrixXd::Zero(w.rank(),w.pseudo_count()));
    };
    const auto fallback=[&]() {
      if (profile) ++profile->fallbacks;
      if (tt_probe) {
        (*tt_probe)<<"# shared-prefix memory/reuse fallback; no checkpoint spectrum for this task\n";
        tt_probe->flush();
      }
      auto output=zeros();
      for (std::size_t k=0;k<family.size();++k) {
        const auto original=fca::self_cut(family[k].topology,family[k].representative,{{0,0,0,0}});
        const auto canonical=prefix_family_cut(family[k]);
        int arrows=0; std::array<int,3> permutation{};
        for (int j=0;j<3;++j) for (int old=0;old<3;++old)
          if (canonical.split[j].original==original.split[old].original) {
            permutation[j]=old; arrows|=((directions>>j)&1)<<old;
          }
        // This task sums canonical pole zero, which need not be the old
        // adaptive pole anchor. Use the existing grouped implementation with
        // that anchor to give exactly this slice while preserving early sums.
        std::array<int,3> p{};
        p[permutation[1]]=fixed[0]; p[permutation[2]]=fixed[1];
        if (options.self_reuse) {
          if (options.indexed_states)
            output[k]=self_nodes_impl<IndexedStates>(family[k].topology,family[k].representative,
                          arrows,p,options,nullptr,permutation[0],sum_directions ? permutation[0] : -1);
          else output[k]=self_nodes_impl<States>(family[k].topology,family[k].representative,
                          arrows,p,options,nullptr,permutation[0],sum_directions ? permutation[0] : -1);
          continue;
        }
        for (int l=0;l<w.pole_count();++l) {
          p[permutation[0]]=l;
          for (int a=0;a<(sum_directions ? 2 : 1);++a)
            output[k]+=self_nodes(family[k].topology,family[k].representative,
                                  arrows^(a<<permutation[0]),p,options);
        }
      }
      return output;
    };
    if (!options.shared_prefix_bytes || !options.self_reuse) return fallback();
    struct BudgetExceeded {};
    // Put allocations inside the try scope: a budget fallback must RELEASE
    // the prefix/branch accumulators before starting the old evaluator.
    try {
      if (profile && options.bath_sector_profile)
        profile->bath_sectors.resize(cix.Nbath);
      if (profile && options.bath_sector_profile && bath_block_entries(cix.Nbath,cix.Nblock)) {
        profile->bath_block_prefix.assign(bath_block_entries(cix.Nbath,cix.Nblock),0.);
        profile->bath_block_intermediate.assign(bath_block_entries(cix.Nbath,cix.Nblock),0.);
      }
      std::vector<fca::Cut> cuts;
      for (const auto& member:family) {
        auto cut=prefix_family_cut(member);
        for (int j=0;j<3;++j) cut.split[j].d=(directions>>j)&1;
        cuts.push_back(cut);
      }
      std::vector<States> merged(family.size());
      std::vector<Decorations> decorations(family.size());
      std::vector<Prepared> prepared(family.size());
      Engine::ConvolutionScratch scratch;
      States pending_hilbert;
      bool pending_coefficients=false;
      const auto guard=[&](const States* prefix,const States* current=nullptr) {
        long double bytes=0;
        if (pending_coefficients) bytes+=storage_bytes(pending_hilbert);
        for (const auto& m:merged) bytes+=storage_bytes(m);
        if (prefix) bytes+=storage_bytes(*prefix);
        if (current && current!=prefix) bytes+=storage_bytes(*current);
        const auto accounted=static_cast<std::size_t>(std::min<long double>(
          bytes,std::numeric_limits<std::size_t>::max()));
        if (profile) profile->peak_bytes=std::max(profile->peak_bytes,accounted);
        if (bytes>options.shared_prefix_bytes) throw BudgetExceeded{};
      };
      const auto convolve=[&](States& states,int edge,const Prepared& edges,bool prefix) {
        if (edge>=7 || states.empty()) return;
        std::vector<State> keys; std::vector<Function> batch;
        std::vector<Function> hilbert;
        for (auto& [key,h]:states) {
          keys.push_back(key); batch.push_back(std::move(h));
          if (pending_coefficients) hilbert.push_back(std::move(pending_hilbert.at(key)));
        }
        if (profile) (prefix ? profile->prefix_convolutions : profile->tail_convolutions)+=batch.size();
        Engine::ConvolveRightPseudoPreparedBatch(batch,cix,w,*edges[edge],false,
                                                 options.convolution_scratch ? &scratch : nullptr,
                                                 pending_coefficients ? &hilbert : nullptr);
        pending_hilbert.clear(); pending_coefficients=false;
        states.clear();
        for (std::size_t i=0;i<keys.size();++i) if (!batch[i].empty())
          states.emplace(keys[i],std::move(batch[i]));
      };
      const auto advance=[&](const States& states,int vertex,const fca::Cut& cut,
                            const Decorations& d,const std::array<int,3>& pole,int bath0) {
        States next;
        pending_hilbert.clear(); pending_coefficients=false;
        if (vertex==4) {
          for (const auto& [state,h]:states) {
            if (profile) ++profile->root_functions;
            auto value=options.root_plans ? planned_root(h) : root_contract(h);
            if (!value.empty()) next.emplace(state,std::move(value));
          }
          return next;
        }
        int j=0;
        while (j<3 && cut.split[j].a!=vertex && cut.split[j].b!=vertex) ++j;
        if (j==3) throw std::logic_error("unmatched FCA prefix vertex");
        const bool opening=cut.split[j].a==vertex;
        // Only unbound openings branch here. The streamed first bath is
        // fixed; closures must first sum their NODE values, so neither is
        // converted by this path. Dress BEFORE transforming, then commute
        // the time-independent F through BOTH coefficient transforms.
        bool share=false;
        if (options.shared_branch_transforms && options.extended_branch_transforms &&
            opening && j!=0 && cix.Nbath>1) {
          double input=0,output=0;
          for (const auto& [state,h]:states) {
            const auto count=branch_entries(h,Engine::Opposite(cut.split[j].d),pole[j]);
            input+=count.first; output+=count.second;
          }
          const double temporary=5.*storage_bytes(states)+32.*w.rank()*output;
          share=output>1.5*input && temporary<=options.self_prefix_bytes;
        }
        pending_coefficients=share;
        for (const auto& [state,h]:states) {
          std::pair<Function,Function> transformed;
          if (share) {
            transformed=prepare_left(h,d.vertex[vertex]);
            if (profile) ++profile->shared_branch_inputs;
          }
          const int first=opening ? (j==0 ? bath0 : 0) : state[j];
          const int end=opening ? (j==0 ? bath0+1 : cix.Nbath) : first+1;
          for (int bath=first;bath<end;++bath) {
            if (bath<0) throw std::logic_error("FCA prefix closes an unbound bath");
            const double weight=residue(bath,pole[j]);
            if (weight==0.) continue;
            const int op=2*bath+cut.split[j].d;
            const auto& matrix=operators[opening ? Engine::Opposite(op) : op];
            auto value=share ? Engine::MultiplyRight(transformed.first,matrix,ones) :
              Engine::MultiplyRight(h,matrix,((opening ? 1. : weight)*d.vertex[vertex]).eval());
            if (value.empty()) continue;
            State target=state; target[j]=opening ? bath : -1;
            Engine::AddFunction(next[target],std::move(value));
            if (share) Engine::AddFunction(pending_hilbert[target],
              Engine::MultiplyRight(transformed.second,matrix,ones));
          }
        }
        return next;
      };
      std::array<int,3> pole{{0,fixed[0],fixed[1]}};
      int selected_probe_pole=options.tt_probe_pole;
      int selected_probe_bath=options.tt_probe_bath;
      int recorded_probe_slices=0;
      for (int a=0;a<(sum_directions ? 2 : 1);++a) {
      for (auto& cut:cuts) cut.split[0].d=(directions&1)^a;
      for (int l=0;l<w.pole_count();++l) {
        pole[0]=l;
        for (std::size_t k=0;k<family.size();++k) {
          decorations[k]=decorate(w,cuts[k],pole,false);
          prepared[k]=prepare_edges(decorations[k],7,-1,options);
          // An explicit numerical guard protects this fast path if the
          // decoration convention ever changes. Factors outside the shared
          // prefix (including inverse_denominator) intentionally differ.
          for (int v=0;v<=4;++v)
            if (!(decorations[k].edge[v].array()==decorations[0].edge[v].array()).all() ||
                !(decorations[k].vertex[v].array()==decorations[0].vertex[v].array()).all())
              throw std::logic_error("FCA shared prefix decorations differ");
        }
        for (int bath=0;bath<cix.Nbath;++bath) {
          if (residue(bath,l)==0.) continue;
          // Auto-select once and retain the same slice for both arrows.
          if (tt_probe && selected_probe_pole<0) selected_probe_pole=l;
          if (tt_probe && l==selected_probe_pole && selected_probe_bath<0)
            selected_probe_bath=bath;
          const bool record_probe=tt_probe && l==selected_probe_pole && bath==selected_probe_bath;
          if (record_probe) ++recorded_probe_slices;
          double prescore=0,intermediate_squared=0;
          std::chrono::steady_clock::time_point sector_start;
          States prefix;
          prefix[{{-1,-1,-1}}]=Engine::PseudoValues(pseudos,decorations[0].edge[0]);
          if (profile && options.bath_sector_profile) {
            const auto& op=operators[Engine::Opposite(2*bath+cuts[0].split[0].d)];
            for (const auto& [state,h]:prefix) for (const auto& [key,values]:h) {
              if (key.second<0 || key.second>=int(op.by_source.size())) continue;
              const int index=op.by_source[key.second];
              if (index<0) continue;
              double prefix_squared=0;
              for (const auto& value:values) prefix_squared+=value.squaredNorm();
              const double score=std::sqrt(prefix_squared)*op.entries[index].matrix.norm();
              prescore+=score;
              if (!profile->bath_block_prefix.empty() && key.first>=0 &&
                  key.first<cix.Nblock && key.second<cix.Nblock)
                profile->bath_block_prefix[(std::size_t(bath)*cix.Nblock+key.first)*cix.Nblock+
                  key.second]+=std::abs(residue(bath,l))*score;
            }
            prescore*=std::abs(residue(bath,l));
            sector_start=std::chrono::steady_clock::now();
          }
          guard(&prefix);
          auto previous_vertex=std::chrono::steady_clock::now();
          int last_vertex=0;
          for (int v=1;v<=4 && !prefix.empty();++v) {
            prefix=advance(prefix,v,cuts[0],decorations[0],pole,bath);
            guard(&prefix); convolve(prefix,v,prepared[0],true); guard(&prefix);
            last_vertex=v;
            const auto finished_vertex=std::chrono::steady_clock::now();
            if (record_probe && v>=2)
              write_tt_probe_spectrum(*tt_probe,prefix,w.basis(),v,a,l,bath,
                std::chrono::duration<double>(finished_vertex-previous_vertex).count());
            // Exclude the diagnostic SVD from the next vertex's timing.
            previous_vertex=std::chrono::steady_clock::now();
          }
          if (record_probe && prefix.empty() && last_vertex<4) {
            (*tt_probe)<<"# selected prefix became empty before vertex 4; later checkpoint absent\n";
            tt_probe->flush();
          }
          for (std::size_t k=0;k<family.size();++k) {
            States tail;
            const int closing=cuts[k].split[0].b;
            for (int v=5;v<=closing;++v) {
              tail=advance(v==5 ? prefix : tail,v,cuts[k],decorations[k],pole,bath);
              guard(&prefix,&tail);
              if (v<closing) convolve(tail,v,prepared[k],false);
              guard(&prefix,&tail);
              if (tail.empty()) break;
            }
            // Normalization belongs to THIS diagram/pole triple and enters
            // BEFORE the pole sum, just as in the independent grouped path.
            for (auto& [state,h]:tail) {
              for (auto& [key,values]:h) for (auto& matrix:values)
                matrix*=decorations[k].inverse_denominator;
              if (profile && options.bath_sector_profile)
                for (const auto& [key,values]:h) {
                  double block_squared=0;
                  for (const auto& value:values) block_squared+=value.squaredNorm();
                  intermediate_squared+=block_squared;
                  if (!profile->bath_block_intermediate.empty() && key.first>=0 &&
                      key.first<cix.Nblock && key.second>=0 && key.second<cix.Nblock)
                    profile->bath_block_intermediate[
                      (std::size_t(bath)*cix.Nblock+key.first)*cix.Nblock+key.second]+=
                        std::sqrt(block_squared);
                }
              Engine::AddFunction(merged[k][state],std::move(h));
            }
            tail.clear(); guard(&prefix);
          }
          if (profile && options.bath_sector_profile) {
            auto& sector=profile->bath_sectors[bath];
            ++sector.visits;
            sector.seconds+=std::chrono::duration<double>(
              std::chrono::steady_clock::now()-sector_start).count();
            sector.prescore+=prescore;
            sector.intermediate_norm+=std::sqrt(intermediate_squared);
          }
        }
      }
      } // both anchor arrows have closed before the shared suffix
      if (tt_probe && recorded_probe_slices==0) {
        (*tt_probe)<<"# selected bath/pole has zero residue in this task\n";
        tt_probe->flush();
      }
      auto output=zeros();
      for (std::size_t k=0;k<family.size();++k) {
        auto states=std::move(merged[k]);
        const int closing=cuts[k].split[0].b;
        convolve(states,closing,prepared[k],false); guard(&states);
        // Pole zero has closed: later decorations no longer depend on it.
        for (int v=closing+1;v<8 && !states.empty();++v) {
          states=advance(states,v,cuts[k],decorations[k],pole,-1);
          guard(&states); convolve(states,v,prepared[k],false); guard(&states);
        }
        for (const auto& [state,h]:states) {
          if (state!=State{{-1,-1,-1}}) throw std::logic_error("unclosed FCA shared-prefix bath");
          const int t=family[k].topology;
          Engine::AccumulateCut(h,cix,fca::symmetry[t]*fca::crossing_sign[t],degeneracy,output[k]);
        }
      }
      return output;
    } catch (const BudgetExceeded&) {
      // The abandoned partial family must not be counted as completed work.
      if (profile) {
        profile->bath_sectors.clear();
        profile->bath_block_prefix.clear();
        profile->bath_block_intermediate.clear();
      }
      return fallback();
    }
  }
  // One common geometry-only planner for old and arrow-grouped evaluation.
  // Keeping this choice identical prevents an algebraic saving in the tail
  // from being outweighed by a more expensive streamed middle.
  int bath_anchor(int topology,int cut_edge,const Options& options) const {
    const auto cut=fca::self_cut(topology,cut_edge,{{0,0,0,0}});
    int anchor=0;
    for (int j=1;j<3;++j) if (cut.split[j].b>cut.split[anchor].b) anchor=j;
    if (options.self_adaptive_anchor) {
      // A geometry-only work estimate: count bath states processed by each
      // event/convolution, including repeated streamed intervals. Root
      // contractions include the sum over 2*Nbath operators. This is a
      // heuristic, not a timing prediction: block sparsity varies by bath.
      // Reject an unstreamed chain if it would retain THREE open bath labels.
      // All accepted plans keep state-count scaling at most Nbath^2.
      double best=std::numeric_limits<double>::infinity();
      double best_peak=std::numeric_limits<double>::infinity();
      double best_prefix=std::numeric_limits<double>::infinity();
      for (int candidate=-1;candidate<3;++candidate) {
        double cost=0,peak=0; bool valid=true;
        for (int v=1;v<8;++v) {
          int before=0,after=0;
          for (int j=0;j<3;++j) if (j!=candidate) {
            before+=cut.split[j].a<v && v<=cut.split[j].b;
            after+=cut.split[j].a<=v && v<cut.split[j].b;
          }
          if (before>2 || after>2) { valid=false; break; }
          const double nb=std::pow(static_cast<double>(cix.Nbath),before);
          const double na=std::pow(static_cast<double>(cix.Nbath),after);
          const bool inside=candidate>=0 && cut.split[candidate].a<=v && v<=cut.split[candidate].b;
          const double repeat=inside ? cix.Nbath : 1;
          cost+=repeat*(v==cut.q ? 2*cix.Nbath*nb : std::max(nb,na));
          // The closing convolution is AFTER the streamed-bath sum.
          if (v<7) cost+=na*(inside && v<cut.split[candidate].b ? cix.Nbath : 1);
          peak=std::max(peak,std::max(nb,na));
        }
        // Streaming different labels often has identical arithmetic count.
        // Break that tie by live-state count, then retained-prefix count;
        // otherwise choosing the first label can retain an Nbath^2 prefix
        // needlessly, especially expensive in the superconducting blocks.
        int prefix_open=0;
        if (candidate>=0) for (int j=0;j<3;++j)
          prefix_open+=cut.split[j].a<cut.split[candidate].a &&
                       cut.split[candidate].a<cut.split[j].b;
        const double prefix=candidate<0 ? 0 : std::pow(static_cast<double>(cix.Nbath),prefix_open);
        if (valid && (cost<best || (cost==best &&
            (peak<best_peak || (peak==best_peak && prefix<best_prefix))))) {
          anchor=candidate; best=cost; best_peak=peak; best_prefix=prefix;
        }
      }
    }
    return anchor;
  }

  template<class StateMap>
  MatrixXd self_nodes_impl(int topology,int cut_edge,int directions,
                      const std::array<int,3>& initial_pole,const Options& options,
                      SelfProfile* profile,int summed_pole=-1,int summed_direction=-1) const {
    auto pole=initial_pole;
    const auto independent=[&]() {
      MatrixXd result=MatrixXd::Zero(w.rank(),w.pseudo_count());
      auto p=initial_pole;
      for (int a=0;a<(summed_direction>=0 ? 2 : 1);++a)
        for (int l=0;l<(summed_pole>=0 ? w.pole_count() : 1);++l) {
          if (summed_pole>=0) p[summed_pole]=l;
          const int arrows=directions^(summed_direction>=0 ? a<<summed_direction : 0);
          result+=self_nodes(topology,cut_edge,arrows,p,options,profile);
        }
      return result;
    };
    if (!options.self_reuse && (summed_pole>=0 || summed_direction>=0)) return independent();
    if (!options.self_reuse) return self_nodes_reference(topology,cut_edge,directions,pole,profile);
    auto cut=fca::self_cut(topology,cut_edge,{{0,0,0,0}});
    for (int j=0;j<3;++j) cut.split[j].d=(directions>>j)&1;
    Decorations d=decorate(w,cut,pole,false);
    // A coefficient chain is meaningful only when the small DLR solve used
    // for coefficient output passed the basis's stability gate. All retained
    // states have the same representation after a given backbone edge.
    const bool chain=options.self_coeff_chain &&
      w.basis().coefficient_convolution_available();
    std::array<MatrixXd,8> vertex_coefficient_map;
    std::array<bool,8> vertex_has_factor{};
    const auto refresh_vertex_maps=[&]() {
      if (!chain) return;
      const auto& B=w.basis().coefficients_to_values_matrix();
      for (int v=2;v<7;++v) if (v!=cut.q) {
        vertex_has_factor[v]=!(d.vertex[v].array()==1.).all();
        if (vertex_has_factor[v])
          vertex_coefficient_map[v]=w.basis().values_to_coefficients(
            (d.vertex[v].asDiagonal()*B).eval());
      }
    };
    refresh_vertex_maps();
    int anchor=bath_anchor(topology,cut_edge,options);
    if (summed_pole>=0) anchor=summed_pole;
    if (summed_direction>=0) anchor=summed_direction;
    if (profile) profile->anchor=anchor;
    auto prepared=prepare_edges(d,7,-1,options,profile);
    const auto make_states=[&]() {
      if constexpr (std::is_same_v<StateMap,IndexedStates>) return StateMap(cix.Nbath,anchor);
      else return StateMap{};
    };
    Engine::ConvolutionScratch scratch;
    scratch.dlr_buffers.profile=profile!=nullptr;
    if (profile && options.bath_sector_profile)
      profile->bath_sectors.resize(cix.Nbath);
    if (profile && options.bath_sector_profile && bath_block_entries(cix.Nbath,cix.Nblock)) {
      profile->bath_block_prefix.assign(bath_block_entries(cix.Nbath,cix.Nblock),0.);
      profile->bath_block_intermediate.assign(bath_block_entries(cix.Nbath,cix.Nblock),0.);
    }
    std::vector<State> reusable_keys;
    std::vector<Function> reusable_batch;
    auto pending_hilbert=make_states();
    bool pending_coefficients=false;
    auto anchor_coefficients=make_states(),anchor_hilbert=make_states();
    bool anchor_prepared=false;
    const auto share_worthwhile=[&](const StateMap& states,int j) {
      if (!options.shared_branch_transforms || cix.Nbath<=1) return false;
      double input_entries=0,output_entries=0,output_bytes=0;
      for (const auto& [state,h]:states) for (const auto& [key,values]:h) {
        input_entries+=values[0].size();
        for (int bath=0;bath<cix.Nbath;++bath) {
          if (residue(bath,pole[j])==0.) continue;
          const auto& op=operators[Engine::Opposite(2*bath+cut.split[j].d)];
          const int index=op.by_source[key.second];
          if (index<0) continue;
          const double entries=double(values[0].rows())*op.entries[index].matrix.cols();
          output_entries+=entries;
          output_bytes+=w.rank()*(8*entries+sizeof(MatrixXd)+64)+256;
        }
      }
      // Count ACTUAL compatible branches, not nominal Nbath. Include target
      // block sizes and the retained input/transform buffers in the cap.
      return output_entries>1.5*input_entries &&
        5.*storage_bytes(states)+2.*output_bytes+
          storage_bytes(anchor_coefficients)+storage_bytes(anchor_hilbert)<=options.self_prefix_bytes;
    };
    const auto convolve=[&](StateMap& states,int edge) {
      if (edge>=7 || states.empty()) return;
      SectionTimer timer(profile ? &profile->convolution_seconds : nullptr);
      std::vector<State> local_keys;
      std::vector<Function> local_batch;
      auto& keys=options.convolution_scratch ? reusable_keys : local_keys;
      auto& batch=options.convolution_scratch ? reusable_batch : local_batch;
      std::vector<Function> left_hilbert;
      keys.clear(); batch.clear();
      for (auto& [key,h]:states) {
        keys.push_back(key); batch.push_back(std::move(h));
        if (pending_coefficients) left_hilbert.push_back(std::move(pending_hilbert[key]));
      }
      if (profile) { profile->convolution_functions+=batch.size(); ++profile->convolution_batches; }
      Engine::ConvolveRightPseudoPreparedBatch(batch,cix,w,*prepared[edge],false,
                                               options.convolution_scratch ? &scratch : nullptr,
                                               pending_coefficients ? &left_hilbert : nullptr,
                                               chain && edge>1 && edge!=cut.q,
                                               chain && edge<6);
      pending_hilbert.clear(); pending_coefficients=false;
      states.clear();
      for (std::size_t k=0;k<keys.size();++k)
        if (!batch[k].empty()) states.emplace(keys[k],std::move(batch[k]));
    };
    // Read-only input lets every streamed bath use the same prefix without
    // deep-copying its block matrices. Only newly produced functions move.
    const auto advance=[&](const StateMap& states,int vertex,int anchor_bath) {
      SectionTimer timer(profile ? (vertex==cut.q ? &profile->root_seconds : &profile->vertex_seconds) : nullptr);
      auto next=make_states();
      pending_hilbert.clear(); pending_coefficients=false;
      if (vertex==cut.q) {
        for (const auto& [state,h]:states) {
          if (profile) ++profile->root_functions;
          // The root bath is a pointwise function of time. Convert its input
          // back to representative time values before the root contraction.
          Function node_input;
          const Function* input=&h;
          if (chain && vertex>1) {
            node_input=transform_time(h,w.basis().coefficients_to_values_matrix());
            input=&node_input;
          }
          auto value=options.root_plans ? planned_root(*input) : root_contract(*input);
          if (!value.empty()) next.emplace(state,std::move(value));
        }
      } else {
        int j=0;
        while (j<3 && cut.split[j].a!=vertex && cut.split[j].b!=vertex) ++j;
        if (j==3) throw std::logic_error("FCA unmatched vertex");
        const auto& line=cut.split[j];
        const bool opening=line.a==vertex;
        // Branches differ by time-INDEPENDENT sparse F matrices. Transform
        // the common scalar-dressed input once, then apply F to coefficients
        // and H*coefficients. The next convolution accepts both directly.
        // Closing vertices deliberately keep the old node-level bath sum.
        // Cap the extra coefficient/Hilbert storage conservatively; large
        // atomic spaces fall back to the original batched node transforms.
        const bool share=!chain && opening && vertex<7 &&
          (j==anchor ? anchor_prepared : share_worthwhile(states,j));
        pending_coefficients=share;
        for (const auto& [state,h]:states) {
          // B^{-1} diag(d_v) B performs precisely the node dressing followed
          // by interpolation that the old path would do at the next
          // convolution. The bath operator is time independent and therefore
          // commutes with this transform; do it before bath branching.
          Function dressed;
          const Function* coefficient_input=&h;
          if (chain && vertex>1 && vertex<7 && vertex_has_factor[vertex]) {
            dressed=transform_time(h,vertex_coefficient_map[vertex]);
            coefficient_input=&dressed;
          }
          Function coefficients,hilbert;
          const Function* cf=nullptr; const Function* hf=nullptr;
          if (share && j==anchor) {
            cf=&anchor_coefficients.at(state); hf=&anchor_hilbert.at(state);
          } else if (share) {
            const auto begin=std::chrono::steady_clock::now();
            auto transformed=prepare_left(h,d.vertex[vertex]);
            coefficients=std::move(transformed.first); hilbert=std::move(transformed.second);
            cf=&coefficients; hf=&hilbert;
            if (profile) {
              ++profile->shared_branch_inputs;
              profile->shared_transform_seconds+=std::chrono::duration<double>(
                std::chrono::steady_clock::now()-begin).count();
            }
          }
          const int first=opening ? (j==anchor ? anchor_bath : 0) : state[j];
          const int end=opening ? (j==anchor ? anchor_bath+1 : cix.Nbath) : first+1;
          for (int bath=first;bath<end;++bath) {
            if (bath<0) throw std::logic_error("FCA bath closed before opening");
            // Zero residues are safe to skip already at opening. No weight
            // is applied twice: the actual residue still enters at closure.
            if (residue(bath,pole[j])==0.0) continue;
            const double weight=opening ? 1.0 : residue(bath,pole[j]);
            const int op=2*bath+line.d;
            const auto& op_matrix=operators[opening ? Engine::Opposite(op) : op];
            auto value=share ? Engine::MultiplyRight(*cf,op_matrix,ones) :
              (chain && vertex>1 && vertex<7) ?
                Engine::MultiplyRight(*coefficient_input,op_matrix,(weight*ones).eval()) :
                Engine::MultiplyRight(h,op_matrix,(weight*d.vertex[vertex]).eval());
            if (value.empty()) continue;
            State target=state; target[j]=opening ? bath : -1;
            Engine::AddFunction(next[target],std::move(value));
            if (share) Engine::AddFunction(pending_hilbert[target],
              Engine::MultiplyRight(*hf,op_matrix,ones));
          }
        }
      }
      return next;
    };
    auto states=make_states();
    states[{{-1,-1,-1}}]=Engine::PseudoValues(pseudos,d.edge[0]);
    if (anchor<0) {
      // Some chord geometries never have three simultaneously open labels;
      // they can sum every bath locally without streaming an anchor at all.
      for (int v=1;v<8 && !states.empty();++v) {
        states=advance(states,v,-1); convolve(states,v);
      }
    } else {
      const int opening=cut.split[anchor].a,closing=cut.split[anchor].b;
      for (int v=1;v<opening && !states.empty();++v) {
        states=advance(states,v,-1); convolve(states,v);
      }
      const std::size_t bytes=storage_bytes(states);
      if (profile) profile->prefix_bytes=bytes;
      if (bytes>options.self_prefix_bytes) {
        // This limit covers the EXTRA retained prefix, not total evaluator
        // memory. Large atomic spaces fall back to the original streaming
        // evaluator rather than keeping a large prefix alongside its work.
        states.clear();
        prepared={};
        if (profile) ++profile->prefix_fallbacks;
        if (summed_pole<0 && summed_direction<0)
          return self_nodes_reference(topology,cut_edge,directions,pole,profile);
        // Revert to independent triples with their ORIGINAL adaptive bath
        // plan if the grouped prefix is too large. Jumping straight to the
        // Cartesian reference here would discard previously validated reuse
        // and can be very expensive for SC/f-shell blocks. self_nodes has no
        // grouped-pole argument, so this fallback cannot recurse into itself.
        return independent();
      }
      auto merged=make_states();
      const int pole_count=summed_pole>=0 ? w.pole_count() : 1;
      for (int a=0;a<(summed_direction>=0 ? 2 : 1);++a) {
      if (summed_direction>=0) cut.split[anchor].d=((directions>>anchor)&1)^a;
      for (int l=0;l<pole_count;++l) {
        if (summed_pole>=0) pole[summed_pole]=l;
        if (summed_pole>=0 || summed_direction>=0) {
          d=decorate(w,cut,pole,false);
          refresh_vertex_maps();
          prepared=prepare_edges(d,7,-1,options,profile);
        }
        // The streamed bath is also a branch: keep its transformed prefix
        // once per pole rather than solving it anew in every bath iteration.
        // Endpoint dressing changes with this pole, so do NOT cache across l.
        anchor_prepared=!chain && share_worthwhile(states,anchor);
        if (anchor_prepared) {
          const auto begin=std::chrono::steady_clock::now();
          for (const auto& [state,h]:states) {
            auto transformed=prepare_left(h,d.vertex[opening]);
            anchor_coefficients.emplace(state,std::move(transformed.first));
            anchor_hilbert.emplace(state,std::move(transformed.second));
            if (profile) ++profile->shared_branch_inputs;
          }
          if (profile) {
            const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count();
            profile->shared_transform_seconds+=seconds;
            profile->vertex_seconds+=seconds;
          }
        }
        for (int bath=0;bath<cix.Nbath && !states.empty();++bath) {
          if (residue(bath,pole[anchor])==0.0) continue;
          // Cheap screening proxy: residue times the norms of geometrically
          // compatible prefix blocks and the opening operator. This is NOT
          // the final diagram contribution and is not an error bound.
          double prescore=0;
          std::chrono::steady_clock::time_point sector_start;
          if (profile && options.bath_sector_profile) {
            sector_start=std::chrono::steady_clock::now();
            const auto& op=operators[Engine::Opposite(2*bath+cut.split[anchor].d)];
            for (const auto& [state,h]:states) for (const auto& [key,values]:h) {
              if (key.second<0 || key.second>=int(op.by_source.size())) continue;
              const int index=op.by_source[key.second];
              if (index<0) continue;
              double prefix_squared=0;
              for (const auto& value:values) prefix_squared+=value.squaredNorm();
              const double score=std::sqrt(prefix_squared)*op.entries[index].matrix.norm();
              prescore+=score;
              if (!profile->bath_block_prefix.empty() && key.first>=0 &&
                  key.first<cix.Nblock && key.second<cix.Nblock)
                profile->bath_block_prefix[(std::size_t(bath)*cix.Nblock+key.first)*cix.Nblock+
                  key.second]+=std::abs(residue(bath,pole[anchor]))*score;
            }
            prescore*=std::abs(residue(bath,pole[anchor]));
          }
          auto middle=make_states();
          for (int v=opening;v<=closing;++v) {
            middle=advance(v==opening ? states : middle,v,bath);
            if (v<closing) convolve(middle,v);
            if (middle.empty()) break;
          }
          // The pole normalization used to be applied after the whole chain.
          // It MUST enter before summing different poles at this endpoint.
          double intermediate_squared=0;
          for (auto& [state,h]:middle) {
            if (summed_pole>=0 || summed_direction>=0) for (auto& [key,values]:h)
              for (auto& value:values) value*=d.inverse_denominator;
            if (profile && options.bath_sector_profile)
              for (const auto& [key,values]:h) {
                double block_squared=0;
                for (const auto& value:values) block_squared+=value.squaredNorm();
                intermediate_squared+=block_squared;
                if (!profile->bath_block_intermediate.empty() && key.first>=0 &&
                    key.first<cix.Nblock && key.second>=0 && key.second<cix.Nblock)
                  profile->bath_block_intermediate[
                    (std::size_t(bath)*cix.Nblock+key.first)*cix.Nblock+key.second]+=
                      std::sqrt(block_squared);
              }
            Engine::AddFunction(merged[state],std::move(h));
          }
          if (profile && options.bath_sector_profile) {
            auto& sector=profile->bath_sectors[bath];
            ++sector.visits;
            sector.seconds+=std::chrono::duration<double>(
              std::chrono::steady_clock::now()-sector_start).count();
            sector.prescore+=prescore;
            sector.intermediate_norm+=std::sqrt(intermediate_squared);
          }
        }
        anchor_coefficients.clear(); anchor_hilbert.clear(); anchor_prepared=false;
      }
      } // normalize each arrow BEFORE merging; the suffix has no anchor label
      states=std::move(merged); // releases prefix before the suffix
      convolve(states,closing);
      for (int v=closing+1;v<8 && !states.empty();++v) {
        states=advance(states,v,-1); convolve(states,v);
      }
    }
    MatrixXd output=MatrixXd::Zero(w.rank(),w.pseudo_count());
    for (const auto& [state,h]:states) {
      if (state!=State{{-1,-1,-1}}) throw std::logic_error("FCA unclosed bath");
      Engine::AccumulateCut(h,cix,fca::symmetry[topology]*fca::crossing_sign[topology]*
                            (summed_pole>=0 || summed_direction>=0 ? 1. : d.inverse_denominator),
                            degeneracy,output);
    }
    if (profile) {
      profile->pack_seconds+=scratch.dlr_buffers.pack_seconds;
      profile->transform_seconds+=scratch.dlr_buffers.transform_seconds;
      profile->product_seconds+=scratch.dlr_buffers.product_seconds;
      profile->unpack_seconds+=scratch.dlr_buffers.unpack_seconds;
    }
    return output;
  }

  struct PhysicalArcPair {
    std::shared_ptr<const Function> left,right;
    std::shared_ptr<Engine::PreparedPseudoCoefficients> final_prepared;
    VectorXd final_edge;
    double factor=0;
  };

  // E/F family: retain just the two-stage prefixes for one crossing-bath
  // pair, evaluate all matching members, then release those matrices. This
  // loop order makes reuse immediate (rather than hoping the next MPI task
  // will use a prefix left over from an entirely different bath tuple).
  // Return topology totals, INCLUDING each representative's orbit count.
  std::vector<Physical::Result> physical_ef_family(
      const std::vector<fca::CutOrbit>& family,int directions,
      const std::array<int,3>& pole,const Options& options,
      PhysicalProfile* profile=nullptr) const {
    if (options.indexed_states)
      return physical_ef_family_impl<IndexedStates>(family,directions,pole,options,profile);
    return physical_ef_family_impl<States>(family,directions,pole,options,profile);
  }
  template<class StateMap>
  std::vector<Physical::Result> physical_ef_family_impl(
      const std::vector<fca::CutOrbit>& family,int directions,
      const std::array<int,3>& pole,const Options& options,PhysicalProfile* profile) const {
    const bool cf_chain=coefficient_chain(options);
    const auto zero=[&]() { return MatrixXd::Zero(w.rank(),w.channel_count()).eval(); };
    std::vector<Physical::Result> output(6);
    for (auto& value:output) value={zero(),zero(),zero()};
    if (family.empty()) return output;
    const auto geometry=physical_ef_group(family.front());
    if (geometry.first<0) throw std::invalid_argument("invalid E/F family");
    std::vector<int> arrows;
    std::vector<std::array<int,3>> mapped;
    for (const auto& orbit:family) {
      if (physical_ef_group(orbit)!=geometry) throw std::invalid_argument("incompatible E/F family");
      const auto canonical=physical_family_cut(orbit);
      std::array<int,4> bits{}; bits[orbit.representative/2]=orbit.representative%2;
      const auto original=fca::physical_cut(orbit.topology,orbit.representative/2,bits);
      std::array<int,3> poles{};
      for (int j=0;j<3;++j) {
        bits[canonical.split[j].original]=((directions>>j)&1)^canonical.split[j].d;
        for (int a=0;a<3;++a) if (original.split[a].original==canonical.split[j].original) poles[a]=pole[j];
      }
      int packed=0; for (int l=0;l<4;++l) packed|=bits[l]<<l;
      arrows.push_back(packed); mapped.push_back(poles);
    }
    const auto accumulate=[&](std::size_t k,const Physical::Result& value) {
      const auto& orbit=family[k]; const double count=orbit.members.size();
      auto& target=output[orbit.topology];
      target.g+=count*value.g; target.f00+=count*value.f00; target.f01+=count*value.f01;
    };
    if (!options.physical_ef_shared || !options.physical_reuse || !options.physical_shared_bytes) {
      for (std::size_t k=0;k<family.size();++k)
        accumulate(k,physical_nodes(family[k].topology,family[k].representative/2,
                                   arrows[k],mapped[k],options,profile));
      return output;
    }
    std::vector<int> cross;
    for (int j=0;j<3;++j) if (geometry.second&(1U<<j)) cross.push_back(j);
    if (cross.size()!=2) throw std::logic_error("E/F family requires two crossing baths");
    auto streamed_options=options;
    // Preserve within-member reuse across tuples too. Otherwise cross-
    // diagram sharing could replace one repeated computation by another:
    // e.g. a long-arc prefix depending on only ONE crossing bath would be
    // recomputed for every value of the other bath. These separate caches
    // together share a bounded allowance. Reserve half the existing shared
    // budget for these member caches: a three-member SC family needs room for
    // both a bath-independent branched seed and a later bath-dependent prefix
    // in each member. A naive 64/3 MiB split cannot hold both and regresses.
    // Total retained storage still cannot exceed physical_cache_bytes plus
    // physical_shared_bytes. Explicitly disabling the ordinary cache is kept.
    const std::size_t member_reserve=options.physical_cache_bytes ? options.physical_shared_bytes/2 : 0;
    streamed_options.physical_cache_bytes=options.physical_cache_bytes/family.size()+member_reserve/family.size();
    streamed_options.physical_shared_bytes=options.physical_shared_bytes-member_reserve;
    std::vector<PhysicalOrdinaryCache<StateMap>> member_cache(family.size());
    // The short arc must be identical BEFORE sharing its external-bath
    // functional. Check its complete operation/decorations, not Sindex or
    // topology names. On future geometries fall back to ordinary joins.
    // The mirrored q=5 group offers boundary sharing only, not a same-
    // topology final sum. With that option off use the exact previous path,
    // avoiding even the capture/signature overhead of an inactive feature.
    bool presum=options.physical_ef_presum && (geometry.first==3 ||
      (options.physical_ef_boundary && options.shared_boundary && options.boundary_cache_bytes));
    std::vector<std::uint64_t> first_signature;
    for (std::size_t k=0;k<family.size() && presum;++k) {
      auto cut=physical_family_cut(family[k]);
      for (int j=0;j<3;++j) cut.split[j].d=(directions>>j)&1;
      const auto d=decorate(w,cut,pole,true);
      std::vector<std::uint64_t> signature;
      const int begin=geometry.first==3 ? 0 : 5,end=geometry.first==3 ? 3 : 8;
      for (int v=begin;v<end;++v) {
        for (int n=0;n<w.rank();++n) {
          std::uint64_t bits; double x=d.edge[v](n); std::memcpy(&bits,&x,sizeof(bits));
          signature.push_back(bits);
        }
        if (v>begin) for (int j=0;j<3;++j) {
          const auto& line=cut.split[j];
          if (v==line.a || v==line.b) {
            signature.push_back(j); signature.push_back(line.d); signature.push_back(v==line.a);
          }
        }
      }
      if (k==0) first_signature=signature;
      else presum=signature==first_signature;
    }
    // Keep the first crossing label encountered by the LONG arc outside.
    // For q=3 that label is the second canonical crossing line; in the
    // mirrored q=5 family it is the first. Ignoring this reversal causes
    // the large single-label prefixes to cycle through all baths and get
    // evicted before reuse, despite cross-diagram hits within each pair.
    const int slow=geometry.first==5 ? cross[0] : cross[1];
    const int fast=geometry.first==5 ? cross[1] : cross[0];
    // Reserve, rather than add, storage for the closed E accumulator. This
    // is one bath-pair's function, never a table over crossing bath pairs.
    const std::size_t final_budget=presum && geometry.first==3 ? streamed_options.physical_shared_bytes/4 : 0;
    streamed_options.physical_shared_bytes-=final_budget;
    const auto function_bytes=[](const Function& h) {
      std::size_t bytes=sizeof(Function)+128;
      for (const auto& [key,values]:h) {
        bytes+=sizeof(key)+sizeof(values)+128;
        for (const auto& m:values) bytes+=sizeof(MatrixXd)+64+sizeof(double)*m.size();
      }
      return bytes;
    };
    for (int b1=0;b1<cix.Nbath;++b1) for (int b0=0;b0<cix.Nbath;++b0) {
      State bath{{-1,-1,-1}}; bath[fast]=b0; bath[slow]=b1;
      if (residue(b0,pole[fast])==0 || residue(b1,pole[slow])==0) continue;
      PhysicalArcCache<StateMap> cache;
      BoundaryCache boundary; std::size_t boundary_bytes=0;
      PhysicalSharedProfile boundary_profile;
      PhysicalArcPair pending; Function combined; int combined_members=0;
      const auto finish=[&](const PhysicalArcPair& pair,const Function& right,int topology,double factor) {
        const bool mirrored=geometry.first==5;
        if (options.shared_boundary && options.physical_ef_boundary && accumulate_shared_boundary(
            mirrored ? right : *pair.left,mirrored ? *pair.left : right,
            factor,output[topology],boundary,boundary_bytes,options.boundary_cache_bytes,
            &boundary_profile,mirrored)) return;
        for (int external=0;external<cix.Nbath;++external) {
          const int channel=cix.bath_sigind0(external);
          Physical::AccumulateBoundaryEstimators(*pair.left,right,
            operators[2*external],operators[2*external+1],bulla[external],
            factor*tca_dlr::SignedIndexSign(channel),std::abs(channel),
            Physical::IsF01Channel(cix,external),output[topology]);
        }
      };
      const auto flush=[&]() {
        if (combined.empty()) return;
        std::vector<Function> batch; batch.push_back(std::move(combined));
        Engine::ConvolutionScratch scratch;
        Engine::ConvolveRightPseudoPreparedBatch(batch,cix,w,*pending.final_prepared,true,
                                                 options.convolution_scratch ? &scratch : nullptr,
                                                 nullptr,cf_chain,false);
        if (profile) {
          ++profile->convolution_functions; ++profile->convolution_batches;
          profile->final_convolutions_saved+=combined_members-1;
        }
        finish(pending,batch.front(),4,1.);
        combined.clear(); combined_members=0; pending={};
      };
      for (std::size_t k=0;k<family.size();++k) {
        const int t=family[k].topology;
        const bool deferred=presum && geometry.first==3 && t==4;
        PhysicalArcPair pair;
        const auto value=physical_nodes_impl<StateMap>(t,family[k].representative/2,
          arrows[k],mapped[k],streamed_options,profile,&cache,&bath,&member_cache[k],
          presum ? &pair : nullptr,deferred);
        if (!presum) { accumulate(k,value); continue; }
        if (!pair.left || !pair.right) continue;
        const double factor=pair.factor*family[k].members.size();
        if (!deferred) { finish(pair,*pair.right,t,factor); continue; }
        // E members have the same short arc and final propagator. Include
        // EACH member's denominator, crossing residues and orbit count in
        // X before evaluating (sum X)*G. F stays separate: topology-resolved
        // diagnostics must not be mixed even if its last edge also matches.
        if (pending.left && !(pending.final_edge.array()==pair.final_edge.array()).all()) flush();
        // Worst-case disjoint block sets bound the merged result and the
        // temporary weighted input. A tiny allowance reverts to separate
        // final convolutions without dropping a contribution or boundary.
        if (function_bytes(*pair.left)+2*function_bytes(combined)+
            3*function_bytes(*pair.right)>final_budget) {
          flush();
          std::vector<Function> batch{*pair.right};
          Engine::ConvolutionScratch scratch;
          Engine::ConvolveRightPseudoPreparedBatch(batch,cix,w,*pair.final_prepared,true,
                                                   options.convolution_scratch ? &scratch : nullptr,
                                                   nullptr,cf_chain,false);
          if (profile) { ++profile->convolution_functions; ++profile->convolution_batches; }
          finish(pair,batch.front(),t,factor);
          continue;
        }
        if (!pending.left) pending=pair;
        Function weighted=*pair.right;
        for (auto& [key,values]:weighted) for (auto& m:values) m*=factor;
        Engine::AddFunction(combined,std::move(weighted)); ++combined_members;
      }
      flush();
      if (profile) {
        profile->boundary_hits+=boundary_profile.boundary_hits;
        profile->boundary_builds+=boundary_profile.boundary_builds;
      }
      if (profile) {
        std::size_t total=0; for (const auto& c:member_cache) total+=c.bytes;
        profile->cache_peak_bytes=std::max(profile->cache_peak_bytes,total);
      }
    }
    return output;
  }

  // Physical q=4 family: all three internal lines cross between the arcs.
  // Canonical labels follow the openings 1,2,3. For a fixed pole triple and
  // canonical directions, the entire LEFT chain is identical across these
  // cuts. Stream bath triples, compute that chain once, and immediately join
  // it with each distinct right chain. Never retain a bath-cubed tensor.
  //
  // Short prefixes of BOTH arcs are cached across family members, but only
  // before all three bath labels have appeared. Cache IDs compare complete
  // operator-label histories and sampled decorations, not topology names or
  // Sindex labels. Full atomic block connectivity remains in every Function.
  std::vector<Physical::Result> physical_family(
      const std::vector<fca::CutOrbit>& family,int directions,
      const std::array<int,3>& pole,const Options& options,
      PhysicalSharedProfile* profile=nullptr,bool topology_totals=false) const {
    const bool cf_chain=coefficient_chain(options);
    const auto zero=[&]() { return MatrixXd::Zero(w.rank(),w.channel_count()).eval(); };
    std::vector<Physical::Result> output;
    for (const auto& orbit:family) {
      if (!physical_family_member(orbit)) throw std::invalid_argument("invalid FCA physical family");
    }
    // These are two deliberately distinct result contracts. The default
    // returns one unweighted result per member, as before. Topology totals
    // return six entries A..F with orbit multiplicities ALREADY included.
    for (std::size_t k=0;k<(topology_totals ? 6 : family.size());++k)
      output.push_back({zero(),zero(),zero()});
    if (family.empty()) return output;
    if (!options.physical_reuse || !options.physical_shared_bytes) {
      if (profile) ++profile->fallbacks;
      for (std::size_t k=0;k<family.size();++k) {
        const auto& orbit=family[k];
        const int external=orbit.representative/2;
        std::array<int,4> arrows{}; arrows[external]=orbit.representative%2;
        const auto original=fca::physical_cut(orbit.topology,external,arrows);
        const auto canonical=physical_family_cut(orbit);
        std::array<int,3> mapped{};
        for (int j=0;j<3;++j) {
          // Rotation can reverse a line. Convert canonical directions BACK
          // to the original four-line convention expected by physical_nodes.
          arrows[canonical.split[j].original]=((directions>>j)&1)^canonical.split[j].d;
          for (int a=0;a<3;++a) if (original.split[a].original==canonical.split[j].original)
            mapped[a]=pole[j];
        }
        int bits=0; for (int l=0;l<4;++l) bits|=arrows[l]<<l;
        const auto value=physical_nodes(orbit.topology,external,bits,mapped,options);
        auto& target=output[topology_totals ? orbit.topology : k];
        const double count=topology_totals ? orbit.members.size() : 1.;
        target.g+=count*value.g; target.f00+=count*value.f00; target.f01+=count*value.f01;
      }
      return output;
    }
    std::vector<fca::Cut> cuts;
    std::vector<Decorations> decor;
    std::vector<Prepared> prepared;
    for (const auto& orbit:family) {
      auto cut=physical_family_cut(orbit);
      for (int j=0;j<3;++j) cut.split[j].d=(directions>>j)&1;
      cuts.push_back(cut); decor.push_back(decorate(w,cut,pole,true));
      prepared.push_back(prepare_edges(decor.back(),8,4,options));
      for (int e=0;e<4;++e)
        if (!(decor.back().edge[e].array()==decor.front().edge[e].array()).all())
          throw std::logic_error("FCA physical shared left decorations differ");
    }
    // Group only within a topology, and compare the ACTUAL final-edge
    // decoration. In the q=4 family edge 7 is outside every internal line:
    // positive effective poles do not dress it; negative ones dress it via
    // the complementary interval. Thus it is independent of closing order.
    // Keep the comparison anyway, so a future decoration change falls back
    // to singleton groups instead of silently violating linearity.
    std::vector<std::vector<std::size_t>> final_groups;
    const auto member_scale=[&](std::size_t k) {
      return decor[k].inverse_denominator*(topology_totals ? family[k].members.size() : 1.);
    };
    for (std::size_t k=0;k<family.size();++k) {
      bool joined=false;
      if (topology_totals && options.physical_final_sum)
        for (auto& group:final_groups) {
          const auto lead=group.front();
          const double scale=member_scale(lead);
          if (family[k].topology==family[lead].topology && scale!=0. &&
              std::isfinite(scale) && std::isfinite(member_scale(k)/scale) &&
              (decor[k].edge[7].array()==decor[lead].edge[7].array()).all()) {
            group.push_back(k); joined=true; break;
          }
        }
      if (!joined) final_groups.push_back({k});
    }
    // Plan IDs are local to this pole/direction task. A prefix includes its
    // initial propagator and every following operator/convolution. Reflection
    // is only done at the final right edge, which is deliberately NOT cached.
    struct PrefixPlan { int arc,depth; std::array<int,3> order; std::size_t member; };
    std::vector<PrefixPlan> plans;
    using ArcOrder=std::array<std::array<int,3>,2>;
    std::vector<ArcOrder> order(family.size()),ids(family.size());
    for (std::size_t k=0;k<family.size();++k) for (int arc=0;arc<2;++arc) {
      const int begin=arc*4;
      for (int j=0;j<3;++j)
        order[k][arc][(arc ? cuts[k].split[j].b : cuts[k].split[j].a)-begin-1]=j;
      for (int depth=1;depth<=2;++depth) {
        int id=-1;
        for (std::size_t p=0;p<plans.size() && id<0;++p) {
          const auto& plan=plans[p];
          if (plan.arc!=arc || plan.depth!=depth) continue;
          bool equal=true;
          for (int v=0;v<depth;++v) equal&=plan.order[v]==order[k][arc][v];
          for (int e=begin;e<=begin+depth && equal;++e)
            equal=(decor[k].edge[e].array()==decor[plan.member].edge[e].array()).all();
          if (equal) id=p;
        }
        if (id<0) { id=plans.size(); plans.push_back({arc,depth,order[k][arc],k}); }
        ids[k][arc][depth-1]=id;
      }
    }
    using Ptr=std::shared_ptr<const Function>;
    using Key=std::array<int,4>; // prefix plan, then only encountered bath labels
    struct CacheEntry { Ptr value; std::size_t bytes; std::list<Key>::iterator recency; };
    std::list<Key> recency; // current bath prefixes stay hot; evict old tuples
    std::map<Key,CacheEntry> cache;
    std::size_t retained=0;
    const auto bytes=[](const Function& value) {
      std::size_t size=sizeof(Function)+sizeof(Key)+192;
      for (const auto& [key,matrices]:value) {
        size+=sizeof(key)+sizeof(matrices)+128;
        for (const auto& m:matrices) size+=sizeof(MatrixXd)+sizeof(double)*m.size()+64;
      }
      return size;
    };
    std::vector<std::array<Ptr,2>> initial(family.size());
    for (std::size_t k=0;k<family.size();++k) for (int arc=0;arc<2;++arc)
      initial[k][arc]=std::make_shared<Function>(physical_seed(decor[k].edge[4*arc],options,cf_chain));
    Engine::ConvolutionScratch scratch;
    // Node-prefix reuse and transform reuse are distinct: siblings have
    // different F_b and therefore different node outputs, but the SAME
    // incoming transform. Keep immutable owners with pointer keys so an
    // evicted/reallocated node prefix can never produce a false cache hit.
    // This extra cache is bounded by physical_cache_bytes (normally 64 MiB),
    // independently of the physical_shared_bytes node-prefix cache. Charge
    // input owners as well as coefficients/Hilbert data, conservatively even
    // when the node cache also owns the input. No entry survives this task.
    using LeftPair=std::pair<Function,Function>;
    struct LeftEntry {
      Ptr input;
      std::shared_ptr<const LeftPair> value;
      std::size_t bytes;
      std::list<const Function*>::iterator recency;
    };
    std::map<const Function*,LeftEntry> left_cache;
    std::list<const Function*> left_recency;
    std::size_t left_bytes=0;
    const auto transformed_left=[&](const Ptr& input,int direction,int pole_index)
        ->std::shared_ptr<const LeftPair> {
      if (cf_chain || !options.shared_branch_transforms || !options.extended_branch_transforms || input->empty()) return {};
      const auto hit=left_cache.find(input.get());
      if (hit!=left_cache.end()) {
        left_recency.splice(left_recency.begin(),left_recency,hit->second.recency);
        if (profile) ++profile->branch_transform_hits;
        return hit->second.value;
      }
      const auto count=branch_entries(*input,direction,pole_index);
      const auto size=3*function_bytes(*input)+sizeof(LeftEntry)+128;
      // Five input-size units cover retained input, two transforms and
      // packing work approximately; this is not a process RSS limit.
      if (count.second<=1.5*count.first ||
          5.*function_bytes(*input)+512>options.physical_cache_bytes) return {};
      while (size>options.physical_cache_bytes-left_bytes && !left_recency.empty()) {
        const auto old=left_cache.find(left_recency.back());
        left_bytes-=old->second.bytes;
        left_cache.erase(old); left_recency.pop_back();
      }
      auto result=std::make_shared<const LeftPair>(prepare_left(*input,ones));
      left_recency.push_front(input.get());
      left_cache.emplace(input.get(),LeftEntry{input,result,size,left_recency.begin()});
      left_bytes+=size;
      if (profile) {
        ++profile->shared_branch_inputs;
        profile->transform_peak_bytes=std::max(profile->transform_peak_bytes,left_bytes);
      }
      return result;
    };
    const auto chain=[&](std::size_t k,int arc,const std::array<int,3>& bath,
                         bool before_final=false) {
      Ptr current=initial[k][arc];
      Key key{{-1,-1,-1,-1}};
      for (int depth=1;depth<=3 && !current->empty();++depth) {
        const int j=order[k][arc][depth-1],v=4*arc+depth;
        key[j+1]=bath[j];
        if (depth<3) {
          key[0]=ids[k][arc][depth-1];
          const auto hit=cache.find(key);
          if (hit!=cache.end()) {
            current=hit->second.value;
            recency.splice(recency.begin(),recency,hit->second.recency);
            if (profile) ++profile->cache_hits;
            continue;
          }
        }
        const int op=2*bath[j]+cuts[k].split[j].d;
        const auto& matrix=operators[arc ? op : Engine::Opposite(op)];
        // b2 is the inner bath loop. Only its siblings are visited together;
        // retaining transforms for the more distant b0/b1 siblings thrashes
        // a bounded cache and can cost more than the solves it eliminates.
        // Node-prefix caching already handles those outer-label outputs.
        // The fused tail and all siblings have the SAME representation:
        // coefficients with chaining, nodes with the fallback. Its final
        // convolution below explicitly consumes that representation.
        // The older branch-transform mode must still stop at node values.
        const bool stop=before_final && arc==1 && depth==3;
        const auto left=j==2 && !stop ? transformed_left(current,arc ? cuts[k].split[j].d :
                                        Engine::Opposite(cuts[k].split[j].d),pole[j]) : nullptr;
        auto next=Engine::MultiplyRight(left ? left->first : *current,matrix,ones);
        if (stop) return std::make_shared<const Function>(std::move(next));
        if (!next.empty()) {
          std::vector<Function> batch; batch.push_back(std::move(next));
          std::vector<Function> hilbert;
          if (left) hilbert.push_back(Engine::MultiplyRight(left->second,matrix,ones));
          Engine::ConvolveRightPseudoPreparedBatch(batch,cix,w,*prepared[k][v],arc==1 && depth==3,
                                                    options.convolution_scratch ? &scratch : nullptr,
                                                    left ? &hilbert : nullptr,
                                                    cf_chain,cf_chain && depth<3);
          next=std::move(batch[0]);
          if (profile) ++profile->convolution_functions;
        }
        current=std::make_shared<Function>(std::move(next));
        if (depth<3) {
          const auto size=bytes(*current);
          // This cap bounds RETAINED cache entries, not live chain/scratch
          // storage or RSS. LRU eviction prevents early bath tuples from
          // monopolizing the cache. An evicted object still referenced by a
          // live chain stays alive through its shared_ptr, outside cache
          // accounting. Oversized entries are recomputed without caching.
          if (size<=options.physical_shared_bytes) {
            while (size>options.physical_shared_bytes-retained && !recency.empty()) {
              const auto old=cache.find(recency.back());
              retained-=old->second.bytes;
              cache.erase(old); recency.pop_back();
              if (profile) ++profile->cache_evictions;
            }
            recency.push_front(key);
            cache.emplace(key,CacheEntry{current,size,recency.begin()}); retained+=size;
            if (profile) profile->peak_bytes=std::max(profile->peak_bytes,retained);
          } else if (profile) ++profile->cache_bypasses;
        }
      }
      return current;
    };
    for (int b0=0;b0<cix.Nbath;++b0) for (int b1=0;b1<cix.Nbath;++b1)
      for (int b2=0;b2<cix.Nbath;++b2) {
        const std::array<int,3> bath{{b0,b1,b2}};
        const double weight=residue(b0,pole[0])*residue(b1,pole[1])*residue(b2,pole[2]);
        if (weight==0.) continue;
        const auto left=chain(0,0,bath);
        if (left->empty()) continue;
        // All family members use this exact left arc. Retain its boundary
        // functional only while streaming their right arcs for THIS tuple;
        // never cache a bath-cubed collection or retain six complete rights.
        BoundaryCache boundaries;
        std::size_t boundary_bytes=0;
        for (const auto& group:final_groups) {
          const auto k=group.front();
          Ptr right;
          if (group.size()==1) right=chain(k,1,bath);
          else {
            // Linearity is used AFTER each member's last sparse operator,
            // never across unclosed bath labels. Only this fixed bath tuple
            // is live. Normalization and orbit weights enter the sum; common
            // residues, crossing sign and topology symmetry stay outside.
            // Scaling relative to the leader avoids multiplying every node
            // by the potentially large inverse kernel denominator.
            Function combined;
            int nonempty=0;
            for (auto member:group) {
              const auto input=chain(member,1,bath,true);
              if (input->empty()) continue;
              ++nonempty;
              const double ratio=member_scale(member)/member_scale(k);
              for (const auto& [key,values]:*input) {
                auto& target=combined[key];
                if (target.empty()) {
                  target=values;
                  if (ratio!=1.) for (auto& matrix:target) matrix*=ratio;
                } else for (std::size_t l=0;l<values.size();++l)
                  target[l].noalias()+=ratio*values[l];
              }
            }
            if (profile) ++profile->final_sum_groups;
            if (!combined.empty()) {
              std::vector<Function> batch; batch.push_back(std::move(combined));
              Engine::ConvolveRightPseudoPreparedBatch(batch,cix,w,*prepared[k][7],true,
                  options.convolution_scratch ? &scratch : nullptr,nullptr,cf_chain,false);
              combined=std::move(batch[0]);
              if (profile) {
                ++profile->convolution_functions;
                profile->final_convolutions_saved+=std::max(0,nonempty-1);
              }
            }
            right=std::make_shared<const Function>(std::move(combined));
          }
          if (right->empty()) continue;
          if (profile && k>0) ++profile->left_reuses;
          const int t=family[k].topology;
          auto& target=output[topology_totals ? t : k];
          const double factor=weight*fca::symmetry[t]*fca::crossing_sign[t]*member_scale(k);
          if (options.shared_boundary && accumulate_shared_boundary(*left,*right,factor,
                  target,boundaries,boundary_bytes,options.boundary_cache_bytes,profile)) continue;
          if (options.shared_boundary && profile) ++profile->boundary_fallbacks;
          SectionTimer boundary_timer(profile ? &profile->boundary_reference_seconds : nullptr);
          // Preserve all external bath channels and the creation boundary at
          // zero, including anomalous-channel signs and both Bulla estimators.
          for (int external=0;external<cix.Nbath;++external) {
            const int channel=cix.bath_sigind0(external);
            Physical::AccumulateBoundaryEstimators(*left,*right,operators[2*external],
              operators[2*external+1],bulla[external],factor*tca_dlr::SignedIndexSign(channel),
              std::abs(channel),Physical::IsF01Channel(cix,external),target);
          }
        }
      }
    return output;
  }

  // Production consumes topology totals, not individual oriented cuts.
  // Keep this named entry point explicit so callers cannot accidentally
  // multiply orbit weights twice. physical_family() remains the per-cut
  // reference/diagnostic API and does NOT fuse its member outputs.
  std::vector<Physical::Result> physical_family_topologies(
      const std::vector<fca::CutOrbit>& family,int directions,
      const std::array<int,3>& pole,const Options& options,
      PhysicalSharedProfile* profile=nullptr) const {
    return physical_family(family,directions,pole,options,profile,true);
  }

  // Unmodified Cartesian reference, kept for exact A/B comparisons.
  Physical::Result physical_nodes_reference(int topology,int cut_line,int directions,
                                  const std::array<int,3>& pole,
                                  PhysicalProfile* profile=nullptr) const
  {
    std::array<int,4> arrows{};
    for (int l=0;l<4;++l) arrows[l]=(directions>>l)&1;
    const auto cut=fca::physical_cut(topology,cut_line,arrows);
    const Decorations d=decorate(w,cut,pole,true);
    std::array<Engine::PreparedPseudoCoefficients,8> prepared;
    for (int e=1;e<8;++e) if (e!=cut.q)
      prepared[e]=Engine::PreparePseudoCoefficients(w,d.edge[e]);
    Physical::Result output{MatrixXd::Zero(w.rank(),w.channel_count()),
      MatrixXd::Zero(w.rank(),w.channel_count()),MatrixXd::Zero(w.rank(),w.channel_count())};
    // Reference physical path: bounded-memory bath triples, with all pole
    // dressing/solves outside the bath loop. The external bath is NOT an
    // additional expensive chain: all external G/F connectors share L,R.
    for (int b0=0;b0<cix.Nbath;++b0)
      for (int b1=0;b1<cix.Nbath;++b1)
        for (int b2=0;b2<cix.Nbath;++b2) {
          const std::array<int,3> bath{{b0,b1,b2}};
          double factor=fca::symmetry[topology]*fca::crossing_sign[topology]*d.inverse_denominator;
          std::array<int,8> ops{};
          for (int j=0;j<3;++j) {
            factor*=residue(bath[j],pole[j]);
            const int op=2*bath[j]+cut.split[j].d;
            ops[cut.split[j].a]=Engine::Opposite(op);
            ops[cut.split[j].b]=op;
          }
          if (factor==0.0) continue;
          const auto chain=[&](int begin,int end) {
            Function h=Engine::PseudoValues(pseudos,d.edge[begin]);
            for (int v=begin+1;v<end && !h.empty();++v) {
              h=Engine::MultiplyRight(h,operators[ops[v]],ones);
              if (!h.empty()) {
                if (profile) { ++profile->convolution_functions; ++profile->convolution_batches; }
                h=Engine::ConvolveRightPseudoPrepared(h,cix,w,prepared[v]);
              }
            }
            return h;
          };
          auto first=chain(0,cut.q);
          if (first.empty()) continue;
          auto second=chain(cut.q,8);
          if (second.empty()) continue;
          Physical::Reflect(second,w.basis());
          for (int external=0;external<cix.Nbath;++external) {
            const int channel=cix.bath_sigind0(external);
            Physical::AccumulateBoundaryEstimators(first,second,
              operators[2*external],operators[2*external+1],bulla[external],
              factor*tca_dlr::SignedIndexSign(channel),std::abs(channel),
              Physical::IsF01Channel(cix,external),output);
          }
        }
    return output;
  }

  // Build the two physical arcs while remembering only bath labels that
  // are still needed. A line wholly inside an arc is summed at its closing
  // vertex, BEFORE subsequent convolutions. A line crossing between arcs
  // remains fixed until the two arcs are joined: summing those labels
  // independently would lose their correlation and change the diagram.
  Physical::Result physical_nodes(int topology,int cut_line,int directions,
                                  const std::array<int,3>& pole,
                                  const Options& options=Options{},
                                  PhysicalProfile* profile=nullptr) const
  {
    if (options.indexed_states)
      return physical_nodes_impl<IndexedStates>(topology,cut_line,directions,pole,options,profile);
    return physical_nodes_impl<States>(topology,cut_line,directions,pole,options,profile);
  }
  template<class StateMap>
  Physical::Result physical_nodes_impl(int topology,int cut_line,int directions,
                                  const std::array<int,3>& pole,
                                  const Options& options,PhysicalProfile* profile,
                                  PhysicalArcCache<StateMap>* shared=nullptr,
                                  const State* selected_bath=nullptr,
                                  PhysicalOrdinaryCache<StateMap>* retained=nullptr,
                                  PhysicalArcPair* captured=nullptr,bool before_final=false) const
  {
    const bool cf_chain=coefficient_chain(options);
    if (captured && !selected_bath) throw std::logic_error("FCA arc capture needs one bath tuple");
    if (!options.physical_reuse)
      return physical_nodes_reference(topology,cut_line,directions,pole,profile);
    std::array<int,4> arrows{};
    for (int l=0;l<4;++l) arrows[l]=(directions>>l)&1;
    auto cut=fca::physical_cut(topology,cut_line,arrows);
    auto active_pole=pole;
    if (shared) {
      // Family directions and poles were mapped back into the ordinary API
      // by its caller. Restore canonical opening order for shared StateMap
      // slots. All other ordinary calls retain their original sum order.
      const auto original=cut;
      std::sort(cut.split.begin(),cut.split.end(),[](const auto& a,const auto& b){return a.a<b.a;});
      for (int j=0;j<3;++j) for (int k=0;k<3;++k)
        if (cut.split[j].original==original.split[k].original) active_pole[j]=pole[k];
    }
    const Decorations d=decorate(w,cut,active_pole,true);
    const auto prepared=prepare_edges(d,8,cut.q,options);
    if (before_final && (cut.q!=3 || topology!=4 || !captured))
      throw std::logic_error("FCA deferred final convolution requires q=3 E capture");
    std::array<int,8> event;
    event.fill(-1);
    std::array<bool,3> crossing{};
    std::vector<int> cross_lines;
    for (int j=0;j<3;++j) {
      const auto& line=cut.split[j];
      event[line.a]=event[line.b]=j;
      crossing[j]=line.a<cut.q && line.b>cut.q;
      if (crossing[j]) cross_lines.push_back(j);
    }

    const auto make_states=[&]() {
      if constexpr (std::is_same_v<StateMap,IndexedStates>) return StateMap(cix.Nbath);
      else return StateMap{};
    };
    Engine::ConvolutionScratch scratch;
    using StatePtr=std::shared_ptr<const StateMap>;
    std::vector<State> reusable_keys;
    std::vector<Function> reusable_batch;
    using CacheKey=std::array<int,4>; // last vertex, then encountered cross baths
    PhysicalOrdinaryCache<StateMap> task_cache;
    auto& storage=retained ? *retained : task_cache;
    auto& cache=storage.entries;
    auto& cache_meta=storage.metadata;
    auto& cache_bytes=storage.bytes;
    auto& cache_clock=storage.clock;
    // Dependencies are topological: internal baths are carried in StateMap
    // and summed on closure; only crossing baths label the outer cache key.
    std::array<unsigned,8> dependency{};
    for (int arc=0;arc<2;++arc) {
      unsigned mask=0;
      for (int v=(arc==0 ? 0 : cut.q)+1;v<(arc==0 ? cut.q : 8);++v) {
        if (event[v]<0) throw std::logic_error("FCA physical unmatched vertex");
        if (crossing[event[v]]) mask|=1U<<event[v];
        dependency[v]=mask;
      }
    }
    unsigned all_crossing=0;
    for (int j:cross_lines) all_crossing|=1U<<j;
    const auto prefix_key=[&](int vertex,const State& bath) {
      CacheKey key{{vertex,-1,-1,-1}};
      for (int j=0;j<3;++j) if (dependency[vertex]&(1U<<j)) key[j+1]=bath[j];
      return key;
    };
    const auto cached_bytes=[](const StateMap& states) {
      // Conservative accounting of cached matrix data and container nodes.
      // Live convolution buffers are separately bounded by the TCA batcher.
      std::size_t bytes=sizeof(StateMap)+128;
      if constexpr (std::is_same_v<StateMap,IndexedStates>) bytes+=states.index_bytes();
      for (const auto& [state,function]:states) {
        bytes+=sizeof(State)+sizeof(Function)+128;
        for (const auto& [key,matrices]:function) {
          bytes+=sizeof(key)+sizeof(matrices)+128;
          for (const auto& matrix:matrices)
            bytes+=sizeof(MatrixXd)+sizeof(double)*matrix.size()+64;
        }
      }
      return bytes;
    };
    std::array<StatePtr,2> initial;
    for (int arc=0;arc<2;++arc) {
      auto value=std::make_shared<StateMap>(make_states());
      const int begin=arc==0 ? 0 : cut.q,end=arc==0 ? cut.q : 8;
      // A one-edge arc has no convolution and must remain at nodes.
      (*value)[{{-1,-1,-1}}]=physical_seed(d.edge[begin],options,cf_chain && end>begin+1);
      initial[arc]=value;
    }

    const auto chain=[&](int arc,const State& bath)->StatePtr {
      const int begin=arc==0 ? 0 : cut.q,end=arc==0 ? cut.q : 8;
      StatePtr current=initial[arc];
      int first_vertex=begin+1;
      typename PhysicalArcCache<StateMap>::Key shared_key;
      const int shared_vertex=begin+2;
      if (shared && shared_vertex<end) {
        // Compare the first TWO stages, including the internal bath branch
        // and all positive/negative-pole decorations. A future endpoint is
        // deliberately absent: only operations already executed determine
        // this value. Residues of crossing lines and the global denominator
        // stay outside the cache and are applied independently per member.
        const auto append_edge=[&](int v) {
          for (int i=0;i<d.edge[v].size();++i) {
            std::uint64_t bits; const double value=d.edge[v](i);
            std::memcpy(&bits,&value,sizeof(bits)); shared_key.push_back(bits);
          }
        };
        append_edge(begin);
        for (int v=begin+1;v<=shared_vertex;++v) {
          const int j=event[v]; const auto& line=cut.split[j];
          const bool opening=v==line.a;
          shared_key.push_back(j); // canonical State slot, not original line
          shared_key.push_back(opening); shared_key.push_back(line.d);
          shared_key.push_back(crossing[j]);
          shared_key.push_back(crossing[j] ? bath[j] : active_pole[j]);
          shared_key.push_back(arc==1 && v==end-1); // beta-tau output matters
          shared_key.push_back(cf_chain && v<end-1); // coefficients vs final node values
          append_edge(v);
        }
        const auto hit=shared->entries.find(shared_key);
        if (hit!=shared->entries.end()) {
          current=hit->second; first_vertex=shared_vertex+1;
          if (profile) ++profile->cross_diagram_hits;
        }
      }
      if (options.physical_cache_resume) {
        for (int v=end-1;v>begin;--v) {
          if (v<first_vertex || dependency[v]==all_crossing) continue;
          const auto key=prefix_key(v,bath);
          const auto found=cache[arc].find(key);
          if (found==cache[arc].end()) continue;
          current=found->second;
          first_vertex=v+1;
          cache_meta[arc].at(key).second=++cache_clock;
          if (profile) { ++profile->cache_hits; profile->resumed_vertices+=v-begin; }
          break;
        }
      }
      for (int vertex=first_vertex;vertex<end && !current->empty();++vertex) {
        const int j=event[vertex];
        if (j<0) throw std::logic_error("FCA physical unmatched vertex");
        const auto& line=cut.split[j];
        const bool opening=vertex==line.a;
        const auto key=prefix_key(vertex,bath);
        // Once all crossing baths have appeared, this prefix cannot recur
        // for another tuple. Never cache it: this excludes N_bath^3 storage.
        const bool reusable=dependency[vertex]!=all_crossing;
        if (reusable) {
          const auto found=cache[arc].find(key);
          if (found!=cache[arc].end()) {
            current=found->second;
            cache_meta[arc].at(key).second=++cache_clock;
            if (profile) ++profile->cache_hits;
            continue;
          }
        }
        auto next=std::make_shared<StateMap>(make_states());
        auto next_hilbert=make_states();
        bool share=false;
        // Internal lines branch over baths within this arc. Apply the same
        // prepare-before-branch identity here (including E/F), but leave
        // closing sums at nodes. Cross-arc baths are fixed by the outer tuple
        // in this ordinary evaluator and do not branch in this loop.
        if (!cf_chain && options.shared_branch_transforms && options.extended_branch_transforms &&
            opening && !crossing[j] && cix.Nbath>1) {
          double input=0,output=0;
          for (const auto& [state,h]:*current) {
            const auto count=branch_entries(h,Engine::Opposite(line.d),active_pole[j]);
            input+=count.first; output+=count.second;
          }
          share=output>1.5*input &&
            5.*cached_bytes(*current)+32.*w.rank()*output<=options.physical_cache_bytes;
        }
        for (const auto& [state,h]:*current) {
          std::pair<Function,Function> transformed;
          if (share) {
            transformed=prepare_left(h,ones);
            if (profile) ++profile->shared_branch_inputs;
          }
          const int first=crossing[j] ? bath[j] : (opening ? 0 : state[j]);
          const int last=crossing[j] || !opening ? first+1 : cix.Nbath;
          for (int b=first;b<last;++b) {
            if (b<0) throw std::logic_error("FCA physical bath closed before opening");
            // Cross-line residues are applied once, outside both chains.
            // Internal residues enter at closure, allowing an immediate sum.
            if (!crossing[j] && residue(b,active_pole[j])==0.0) continue;
            const double weight=!crossing[j] && !opening ? residue(b,active_pole[j]) : 1.0;
            const int op=2*b+line.d;
            const auto& matrix=operators[opening ? Engine::Opposite(op) : op];
            auto value=Engine::MultiplyRight(share ? transformed.first : h,matrix,
                                             (weight*ones).eval());
            if (value.empty()) continue;
            State target=state;
            if (!crossing[j]) target[j]=opening ? b : -1;
            Engine::AddFunction((*next)[target],std::move(value));
            if (share) Engine::AddFunction(next_hilbert[target],
              Engine::MultiplyRight(transformed.second,matrix,ones));
          }
        }
        // This state has closed its internal bath but has NOT yet passed
        // through edge 7. Never insert it into a cache of convolved prefixes.
        // Only q=3 E uses this path; its shared two-stage prefix ends earlier.
        if (before_final && arc==1 && vertex==7) return next;
        // Batch all remaining internal bath states through the same DLR
        // transform. The last convolution of the second arc returns values
        // at beta-tau directly, avoiding a separate reflection solve.
        std::vector<State> local_keys;
        std::vector<Function> local_batch;
        auto& keys=options.convolution_scratch ? reusable_keys : local_keys;
        auto& batch=options.convolution_scratch ? reusable_batch : local_batch;
        std::vector<Function> hilbert;
        keys.clear(); batch.clear();
        for (auto& [state,h]:*next) {
          keys.push_back(state); batch.push_back(std::move(h));
          if (share) hilbert.push_back(std::move(next_hilbert[state]));
        }
        if (!batch.empty()) {
          if (profile) { profile->convolution_functions+=batch.size(); ++profile->convolution_batches; }
          Engine::ConvolveRightPseudoPreparedBatch(batch,cix,w,*prepared[vertex],
                                                   arc==1 && vertex==end-1,
                                                   options.convolution_scratch ? &scratch : nullptr,
                                                   share ? &hilbert : nullptr,
                                                   cf_chain,cf_chain && vertex<end-1);
        }
        next->clear();
        for (std::size_t k=0;k<keys.size();++k)
          if (!batch[k].empty()) next->emplace(keys[k],std::move(batch[k]));
        current=next;
        if (shared && vertex==shared_vertex && !shared_key.empty()) {
          const std::size_t bytes=cached_bytes(*next)+sizeof(shared_key)+
            shared_key.capacity()*sizeof(std::uint64_t)+128;
          if (bytes<=options.physical_shared_bytes-shared->bytes &&
              shared->entries.emplace(shared_key,current).second) {
            shared->bytes+=bytes;
            if (profile) profile->cross_diagram_peak_bytes=
              std::max(profile->cross_diagram_peak_bytes,shared->bytes);
          }
        }
        if (reusable && options.physical_cache_bytes) {
          const std::size_t bytes=cached_bytes(*next);
          const auto erase_prefix=[&](int a,const CacheKey& victim) {
            cache_bytes-=cache_meta[a].at(victim).first;
            cache_meta[a].erase(victim); cache[a].erase(victim);
            if (profile) ++profile->cache_evictions;
          };
          // A later prefix with exactly the same crossing-bath dependencies
          // supersedes earlier ones. Often it is also dramatically smaller:
          // an internal closing vertex has already summed its bath branches.
          // Never discard a fallback unless the replacement itself fits.
          const bool retain=options.physical_cache_resume && options.physical_cache_retain;
          if (retain && bytes<=options.physical_cache_bytes) {
            std::vector<CacheKey> superseded;
            std::size_t reclaimable=0;
            const bool terminal=vertex+1==end || dependency[vertex+1]!=dependency[vertex];
            for (const auto& [old,entry]:cache[arc]) {
              bool same=old[0]<vertex && dependency[old[0]]==dependency[vertex];
              for (int l=1;l<4;++l) same=same && old[l]==key[l];
              if (same) {
                superseded.push_back(old);
                // Nonzero-dependency entries are counted by the eviction
                // scan below when terminal is true; do not count twice.
                if (!terminal || dependency[old[0]]==0) reclaimable+=cache_meta[arc].at(old).first;
              }
            }
            if (terminal) for (int a=0;a<2;++a) for (const auto& [old,meta]:cache_meta[a])
              if (dependency[old[0]]!=0) reclaimable+=meta.first;
            const bool can_admit=bytes<=options.physical_cache_bytes-cache_bytes+reclaimable;
            if (can_admit) for (const auto& old:superseded) erase_prefix(arc,old);
            // Only a terminal prefix of this dependency interval may evict
            // other tuples. Bulky temporary branches may use spare capacity
            // but must not push out useful bath-summed prefixes. Protect the
            // bath-independent seed, reused by every tuple, from this churn.
            while (can_admit && terminal && bytes>options.physical_cache_bytes-cache_bytes) {
              int victim_arc=-1; CacheKey victim{};
              std::uint64_t oldest=std::numeric_limits<std::uint64_t>::max();
              for (int a=0;a<2;++a) for (const auto& [old,meta]:cache_meta[a])
                if (dependency[old[0]]!=0 && meta.second<oldest) {
                  oldest=meta.second; victim_arc=a; victim=old;
                }
              if (victim_arc<0) break;
              erase_prefix(victim_arc,victim);
            }
          }
          if (bytes<=options.physical_cache_bytes-cache_bytes) {
            cache[arc].emplace(key,current);
            cache_meta[arc].emplace(key,std::make_pair(bytes,++cache_clock));
            cache_bytes+=bytes;
            if (profile) profile->cache_peak_bytes=std::max(profile->cache_peak_bytes,cache_bytes);
          }
        }
      }
      // Handle a length-one reflected arc, which has no convolution.
      if (arc==1 && end==begin+1) {
        auto reflected=std::make_shared<StateMap>(*current);
        for (auto& [state,h]:*reflected) Physical::Reflect(h,w.basis());
        current=reflected;
      }
      return current;
    };

    Physical::Result output{MatrixXd::Zero(w.rank(),w.channel_count()),
      MatrixXd::Zero(w.rank(),w.channel_count()),MatrixXd::Zero(w.rank(),w.channel_count())};
    std::uint64_t tuples=1;
    if (!selected_bath) for (std::size_t j=0;j<cross_lines.size();++j) tuples*=cix.Nbath;
    for (std::uint64_t tuple=0;tuple<tuples;++tuple) {
      auto index=tuple;
      State bath{{-1,-1,-1}};
      double factor=fca::symmetry[topology]*fca::crossing_sign[topology]*d.inverse_denominator;
      for (int j:cross_lines) {
        bath[j]=selected_bath ? (*selected_bath)[j] : index%cix.Nbath; index/=cix.Nbath;
        factor*=residue(bath[j],active_pole[j]);
      }
      if (factor==0.0) continue;
      const auto first=chain(0,bath);
      if (first->empty()) continue;
      const auto second=chain(1,bath);
      if (second->empty()) continue;
      const State closed{{-1,-1,-1}};
      if (first->size()!=1 || second->size()!=1 || !first->count(closed) || !second->count(closed))
        throw std::logic_error("FCA physical unclosed internal bath");
      if (captured) {
        // Aliasing pointers keep the immutable arc alive without copying its
        // matrices, including when the same arc is owned by the prefix cache.
        captured->left=std::shared_ptr<const Function>(first,&first->at(closed));
        captured->right=std::shared_ptr<const Function>(second,&second->at(closed));
        captured->factor=factor; captured->final_edge=d.edge[7];
        captured->final_prepared=prepared[7];
        continue;
      }
      for (int external=0;external<cix.Nbath;++external) {
        const int channel=cix.bath_sigind0(external);
        Physical::AccumulateBoundaryEstimators(first->at(closed),second->at(closed),
          operators[2*external],operators[2*external+1],bulla[external],
          factor*tca_dlr::SignedIndexSign(channel),std::abs(channel),
          Physical::IsF01Channel(cix,external),output);
      }
    }
    return output;
  }
};

// One completed task's diagnostic payload. The norms are sums of Frobenius
// norms BEFORE cancellation between tasks, not estimates of the error from
// dropping a sector. Keeping this fixed-size also avoids a second MPI message.
struct LiveTaskSample {
  std::uint64_t id=0;
  int group=-1,cut=-1,directions=-1;
  int pole[3]={-1,-1,-1};
  double worker_seconds=0,sigma_norm=0,g_norm=0,f00_norm=0,f01_norm=0;
  double convolved_functions=0,root_functions=0;
};
static_assert(std::is_trivially_copyable_v<LiveTaskSample>);

// A task is one (topology, cut, arrow sector, pole triple), with all its
// bath/backbone paths summed locally. Completion already requires a message
// to the dispatcher; attach a tiny timing/norm payload to that same message.
// Rank zero periodically flushes a task CSV and a compact log summary, so
// interrupted long runs still leave useful profiling data. Shared-family
// time is kept separate rather than spuriously assigned to one topology.
template<class Work>
void schedule(std::uint64_t tasks,int rank,int size,bool mpi_enabled,
              std::ostream* log,double interval,const char* phase,Work work,
              bool reduced=false,bool live=false,std::uint64_t evaluation=0)
{
  const auto begin=std::chrono::steady_clock::now();
  auto last=begin;
  constexpr const char* names[]={"A","B","C","D","E","F","shared-self","shared-physical","shared-EF"};
  struct Summary {
    std::uint64_t count=0;
    double seconds=0,sigma=0,g=0,f00=0,f01=0;
  };
  std::array<Summary,9> summary{};
  std::ofstream live_file;
  const std::string live_path=std::string("fca_live_")+phase+
    (evaluation ? "."+std::to_string(evaluation) : std::string{})+".csv";
  if (rank==0 && live) {
    live_file.open(live_path,std::ios::out|std::ios::trunc);
    if (live_file) {
      live_file<<"completed_elapsed_s,task,group,cut,directions,pole0,pole1,pole2,worker_s,"
                 "sigma_gross_norm,g_gross_norm,f00_gross_norm,f01_gross_norm,"
                 "convolved_functions,root_functions\n";
      live_file<<std::setprecision(17);
    } else if (log) {
      (*log)<<"FCA DLR "<<phase<<" live diagnostic unavailable: cannot open "<<live_path<<std::endl;
    }
  }
  const auto record=[&](const LiveTaskSample& sample) {
    if (rank!=0 || !live || sample.group<0 || sample.group>=9) return;
    auto& s=summary[sample.group];
    ++s.count; s.seconds+=sample.worker_seconds; s.sigma+=sample.sigma_norm;
    s.g+=sample.g_norm; s.f00+=sample.f00_norm; s.f01+=sample.f01_norm;
    if (live_file) {
      const double elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count();
      live_file<<elapsed<<','<<sample.id<<','<<names[sample.group]<<','<<sample.cut<<','
        <<sample.directions<<','<<sample.pole[0]<<','<<sample.pole[1]<<','<<sample.pole[2]<<','
        <<sample.worker_seconds<<','<<sample.sigma_norm<<','<<sample.g_norm<<','
        <<sample.f00_norm<<','<<sample.f01_norm<<','<<sample.convolved_functions<<','
        <<sample.root_functions<<'\n';
    }
  };
  // Standalone algebra tests also use this scheduler with void callbacks.
  const auto invoke=[&](std::uint64_t task) -> LiveTaskSample {
    if constexpr (std::is_void_v<std::invoke_result_t<Work,std::uint64_t>>) {
      work(task);
      return {};
    } else return work(task);
  };
  const auto report=[&](std::uint64_t done) {
    const auto now=std::chrono::steady_clock::now();
    if (rank!=0 || (done!=tasks && (interval<=0.0 ||
        std::chrono::duration<double>(now-last).count()<interval))) return;
    const double elapsed=std::chrono::duration<double>(now-begin).count();
    if (log && interval>0.0)
      (*log)<<"FCA DLR "<<phase<<" progress: completed="<<done<<"/"<<tasks
            <<" elapsed_seconds="<<elapsed<<" remaining_seconds="
            <<(done ? elapsed*(tasks-done)/done : 0.0)<<std::endl;
    if (live && log) {
      (*log)<<"FCA DLR "<<phase<<" live: file="<<live_path<<" groups=[";
      bool first=true;
      for (int g=0;g<9;++g) if (summary[g].count) {
        if (!first) (*log)<<",";
        first=false;
        const auto& s=summary[g];
        (*log)<<names[g]<<":tasks="<<s.count<<":worker_s="<<s.seconds;
        if (phase[0]=='s') (*log)<<":gross_sigma_norm="<<s.sigma;
        else (*log)<<":gross_g_f00_f01_norm="<<s.g<<"/"<<s.f00<<"/"<<s.f01;
      }
      (*log)<<"]"<<std::endl;
    }
    if (live_file) live_file.flush();
    last=now;
  };
  if (rank==0 && log)
    (*log)<<"FCA DLR "<<phase<<" start: tasks="<<tasks<<" ranks="<<size
          <<" scheduler="<<(mpi_enabled && size>1 ? "dynamic-dispatch" : "serial")
          <<" cuts="<<(reduced ? "rotation-reduced" : "unreduced")
          <<" orientations=all live_profile="<<(live ? "yes" : "no")<<std::endl;
#ifdef _MPI
  if (mpi_enabled && size>1) {
    constexpr int assignment=27810,completion=27811;
    if (rank==0) {
      std::uint64_t next=0,done=0;
      int active=0;
      for (int worker=1;worker<size;++worker) {
        const std::uint64_t task=next<tasks ? next++ : tasks;
        MPI_Send(&task,1,MPI_UINT64_T,worker,assignment,MPI_COMM_WORLD);
        if (task<tasks) ++active;
      }
      while (active) {
        LiveTaskSample finished;
        MPI_Status status;
        MPI_Recv(&finished,sizeof(finished),MPI_BYTE,MPI_ANY_SOURCE,completion,MPI_COMM_WORLD,&status);
        record(finished);
        ++done;
        const std::uint64_t task=next<tasks ? next++ : tasks;
        MPI_Send(&task,1,MPI_UINT64_T,status.MPI_SOURCE,assignment,MPI_COMM_WORLD);
        if (task==tasks) --active;
        report(done);
      }
      if (done!=tasks) throw std::logic_error("FCA scheduler lost tasks");
    } else {
      for (;;) {
        std::uint64_t task;
        MPI_Recv(&task,1,MPI_UINT64_T,0,assignment,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        if (task==tasks) break;
        const auto task_begin=std::chrono::steady_clock::now();
        LiveTaskSample finished=invoke(task);
        finished.id=task;
        finished.worker_seconds=std::chrono::duration<double>(
          std::chrono::steady_clock::now()-task_begin).count();
        MPI_Send(&finished,sizeof(finished),MPI_BYTE,0,completion,MPI_COMM_WORLD);
      }
    }
    return;
  }
#endif
  if (rank==0) for (std::uint64_t task=0;task<tasks;++task) {
    const auto task_begin=std::chrono::steady_clock::now();
    LiveTaskSample finished=invoke(task);
    finished.id=task;
    finished.worker_seconds=std::chrono::duration<double>(
      std::chrono::steady_clock::now()-task_begin).count();
    record(finished);
    report(task+1);
  }
}

inline Result evaluate(const Cix& cix,const Workspace& w,bool physical,
                       bool diagnostics,int rank,int size,bool mpi_enabled,
                       std::ostream* log,double interval,
                       const std::function<void(const Result&)>& self_complete={},
                       const Options& options=Options{})
{
  if (!(interval>=0.0)) throw std::invalid_argument("invalid FCA progress interval");
  const std::uint64_t p=w.pole_count();
  if (p==0 || p>10000) throw std::invalid_argument("invalid FCA pole count");
  const std::uint64_t triples=p*p*p;
  if (options.self_pole_sum_min_poles<1)
    throw std::invalid_argument("FCA pole-sum minimum must be positive");
  // Two-pole SC samples can lose more from changing the bath-streaming
  // anchor than they gain from pole reuse. Keep the previous schedule in
  // that regime by default; a minimum of one explicitly forces grouping.
  const bool sum_poles=options.self_pole_sum &&
    p>=static_cast<std::uint64_t>(options.self_pole_sum_min_poles);
  // Preserve separate live files for each FCA self-consistency iteration.
  // Every MPI rank enters evaluate the same number of times, but only rank
  // zero opens the file. A new solver process starts again at iteration 0.
  static std::uint64_t live_evaluation_counter=0;
  const std::uint64_t live_evaluation=(options.live_profile || options.bath_sector_profile) ?
    live_evaluation_counter++ : 0;
  static std::uint64_t tt_probe_evaluation_counter=0;
  const std::uint64_t tt_probe_evaluation=tt_probe_evaluation_counter++;
  Result r;
  const auto sigma_zero=[&] { return MatrixXd::Zero(w.rank(),w.pseudo_count()).eval(); };
  const auto physical_zero=[&] { return MatrixXd::Zero(w.rank(),w.channel_count()).eval(); };
  r.sigma=sigma_zero(); r.g=physical_zero(); r.f00=physical_zero(); r.f01=physical_zero();
  for (int t=0;t<6;++t) {
    r.sigma_topology[t]=sigma_zero();
    if (diagnostics) for (auto& cut:r.sigma_cut[t]) cut=sigma_zero();
    r.physical_topology[t]={physical_zero(),physical_zero(),physical_zero()};
  }
  Evaluator evaluator(cix,w);
  const auto self_orbits=fca::cut_orbits(false,options.reduce_cuts);
  const auto physical_orbits=fca::cut_orbits(true,options.reduce_cuts);
  std::vector<fca::CutOrbit> shared_family,ordinary_orbits;
  if (options.shared_prefix && options.self_reuse && options.reduce_cuts && options.shared_prefix_bytes)
    shared_family=evaluator.select_prefix_family(self_orbits,options);
  if (shared_family.size()<2) shared_family.clear();
  for (const auto& orbit:self_orbits)
    if (std::none_of(shared_family.begin(),shared_family.end(),[&](const auto& member) {
          return member.topology==orbit.topology && member.representative==orbit.representative;
        })) ordinary_orbits.push_back(orbit);
  MatrixXd family_stats=MatrixXd::Zero(1,7);
  double family_peak=0;
  // Sum of worker seconds, task counts, and convolution counts by topology.
  // These are not wall times; they identify expensive diagrams separately
  // from load balancing. Cache peak is the largest per-task retained cache.
  MatrixXd profile=MatrixXd::Zero(6,27);
  // One small reduction after the self phase; no per-bath MPI messages and
  // no path-sized storage. Rows are (topology, streamed bath).
  MatrixXd bath_sector_stats=MatrixXd::Zero(7*cix.Nbath,4);
  const std::size_t block_entries=options.bath_sector_profile ?
    bath_block_entries(cix.Nbath,cix.Nblock) : 0;
  MatrixXd bath_block_stats=MatrixXd::Zero(7*block_entries,2);
  const auto add_block_scores=[&](int group,const auto& counters) {
    if (!block_entries || counters.bath_block_prefix.size()!=block_entries ||
        counters.bath_block_intermediate.size()!=block_entries) return;
    for (std::size_t i=0;i<block_entries;++i) {
      bath_block_stats(group*block_entries+i,0)+=counters.bath_block_prefix[i];
      bath_block_stats(group*block_entries+i,1)+=counters.bath_block_intermediate[i];
    }
  };
  double cache_peak=0;
  const auto poles=[&](std::uint64_t id) {
    return std::array<int,3>{{static_cast<int>(id%p),static_cast<int>((id/p)%p),
                              static_cast<int>((id/(p*p))%p)}};
  };
  if (rank==0 && log)
    (*log)<<"FCA evaluator: self_cut_representatives="<<self_orbits.size()
          <<"/48 physical_oriented_representatives="<<physical_orbits.size()
          <<"/48 physical_reuse="<<(options.physical_reuse ? "yes" : "no")
          <<" physical_coeff_chain="<<(options.physical_reuse && options.physical_coeff_chain &&
              w.basis().coefficient_convolution_available() ? "yes" : "no")
          <<" physical_coeff_chain_requested="<<(options.physical_coeff_chain ? "yes" : "no")
          <<" physical_cache_resume="<<(options.physical_cache_resume ? "yes" : "no")
          <<" physical_cache_retain="<<(options.physical_cache_resume && options.physical_cache_retain ? "yes" : "no")
          <<" self_reuse="<<(options.self_reuse ? "yes" : "no")
          <<" self_coeff_chain="<<(options.self_coeff_chain &&
              w.basis().coefficient_convolution_available() ? "yes" : "no")
          <<" self_coeff_chain_requested="<<(options.self_coeff_chain ? "yes" : "no")
          <<" adaptive_anchor="<<(options.self_adaptive_anchor ? "yes" : "no")
          <<" reuse_prepared="<<(options.reuse_prepared ? "yes" : "no")
          <<" root_plans="<<(options.root_plans ? "yes" : "no")
          <<" indexed_states="<<(options.indexed_states ? "yes" : "no")
          <<" convolution_scratch="<<(options.convolution_scratch ? "yes" : "no")
          <<" self_pole_sum="<<(sum_poles ? "yes" : "no")
          <<" self_direction_sum="<<(options.self_direction_sum ? "yes" : "no")
          <<" pole_sum_min_poles="<<options.self_pole_sum_min_poles
          <<" shared_branch_transforms="<<(options.shared_branch_transforms ? "yes" : "no")
          <<" extended_branch_transforms="<<(options.shared_branch_transforms && options.extended_branch_transforms ? "yes" : "no")
          <<" prepared_cache_mb="<<options.prepared_cache_bytes/(1024.*1024.)
          <<" shared_prefix_members="<<shared_family.size()
          <<" shared_prefix_mb="<<options.shared_prefix_bytes/(1024.*1024.)
          <<" physical_cache_mb="<<options.physical_cache_bytes/(1024.*1024.)<<std::endl;
  if (rank==0 && log && options.physical_coeff_chain && options.physical_reuse &&
      !w.basis().coefficient_convolution_available())
    (*log)<<"FCA coefficient-chain fallback: using node-based physical arcs; dlr_tolerance="
          <<w.basis().tolerance()<<" (requires tolerance>=1e-7 and a finite, accurate B*M=T solve)"<<std::endl;
  const std::uint64_t self_pole_tasks=sum_poles ? p*p : triples;
  const int direction_tasks=options.self_direction_sum ? 4 : 8;
  const std::uint64_t ordinary_tasks=ordinary_orbits.size()*direction_tasks*self_pole_tasks;
  const std::uint64_t family_tasks=shared_family.empty() ? 0 : direction_tasks*p*p;
  if (options.tt_probe_task>=0 &&
      static_cast<std::uint64_t>(options.tt_probe_task)>=family_tasks)
    throw std::invalid_argument("fca_tt_probe_task requires an enabled shared-self task");
  if (options.tt_probe_bath < -1 || options.tt_probe_pole < -1 ||
      options.tt_probe_bath>=cix.Nbath || options.tt_probe_pole>=static_cast<int>(p))
    throw std::invalid_argument("FCA TT probe bath or pole out of range");
  if (rank==0 && log && options.tt_probe_task>=0)
    (*log)<<"FCA TT probe enabled: shared_task="<<options.tt_probe_task
          <<" streamed_bath="<<options.tt_probe_bath
          <<" streamed_pole="<<options.tt_probe_pole
          <<" output=fca_tt_probe.self.iter"<<tt_probe_evaluation
          <<".task"<<options.tt_probe_task<<".rank<worker>.csv"
          <<" (diagnostic only; exact evaluation continues)"<<std::endl;
  schedule(ordinary_tasks+family_tasks,rank,size,mpi_enabled,log,interval,"self",[&](std::uint64_t id) {
    const auto start=std::chrono::steady_clock::now();
    LiveTaskSample sample;
    // Dispatch the larger shared-family tasks first. The remaining ordinary
    // tasks then fill workers as families finish, rather than leaving a
    // small, expensive family-only tail on a many-rank machine. When the
    // feature is disabled family_tasks is zero and the old order is exact.
    if (id<family_tasks) {
      const std::array<int,2> fixed{{static_cast<int>(id%p),static_cast<int>((id/p)%p)}};
      const int directions=(id/(p*p))*(options.self_direction_sum ? 2 : 1);
      sample.group=6; sample.directions=directions;
      sample.pole[0]=fixed[0]; sample.pole[1]=fixed[1];
      SharedPrefixProfile counters;
      std::ofstream tt_file;
      std::ostream* tt_output=nullptr;
      if (options.tt_probe_task>=0 && id==static_cast<std::uint64_t>(options.tt_probe_task)) {
        const std::string name="fca_tt_probe.self.iter"+std::to_string(tt_probe_evaluation)+
          ".task"+std::to_string(id)+".rank"+std::to_string(rank)+".csv";
        tt_file.open(name,std::ios::out|std::ios::trunc);
        if (!tt_file) throw std::runtime_error("cannot open FCA TT probe file "+name);
        tt_file<<std::setprecision(17)
          <<"vertex,arrow_branch,streamed_pole,streamed_bath,active_rows,block_keys,"
             "retained_rank,singular_value,relative_tail,unit_row_singular_value,"
             "unit_row_relative_tail,vertex_seconds,spectrum_seconds\n";
        tt_output=&tt_file;
      }
      const auto values=evaluator.self_prefix_family(shared_family,directions,fixed,options,&counters,
                                                     options.self_direction_sum,tt_output);
      if (options.bath_sector_profile)
        for (int bath=0;bath<int(counters.bath_sectors.size());++bath) {
          const auto& sector=counters.bath_sectors[bath];
          const int row=6*cix.Nbath+bath;
          bath_sector_stats(row,0)+=sector.visits;
          bath_sector_stats(row,1)+=sector.seconds;
          bath_sector_stats(row,2)+=sector.prescore;
          bath_sector_stats(row,3)+=sector.intermediate_norm;
        }
      if (options.bath_sector_profile) add_block_scores(6,counters);
      for (std::size_t k=0;k<shared_family.size();++k) {
        const auto& orbit=shared_family[k];
        r.sigma_topology[orbit.topology]+=double(orbit.members.size())*values[k];
        if (diagnostics) for (int cut:orbit.members) r.sigma_cut[orbit.topology][cut]+=values[k];
        if (options.live_profile) sample.sigma_norm+=double(orbit.members.size())*values[k].norm();
      }
      family_stats(0,0)+=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
      family_stats(0,1)+=1; family_stats(0,2)+=counters.fallbacks;
      family_stats(0,3)+=counters.prefix_convolutions;
      family_stats(0,4)+=counters.tail_convolutions;
      family_stats(0,5)+=counters.root_functions;
      family_stats(0,6)+=counters.shared_branch_inputs;
      family_peak=std::max(family_peak,double(counters.peak_bytes));
      sample.convolved_functions=counters.prefix_convolutions+counters.tail_convolutions;
      sample.root_functions=counters.root_functions;
      return sample;
    }
    id-=family_tasks;
    const auto pole=poles(id);
    id/=self_pole_tasks;
    const int directions=id%direction_tasks; id/=direction_tasks;
    const auto& orbit=ordinary_orbits[id];
    const int topology=orbit.topology;
    sample.group=topology; sample.cut=orbit.representative; sample.directions=directions;
    for (int j=0;j<3;++j) sample.pole[j]=pole[j];
    if (sum_poles) sample.pole[2]=-1; // one pole is internally summed, not fixed by this task
    SelfProfile counters;
    const auto value=options.self_direction_sum ?
      evaluator.self_nodes_direction_sum(topology,orbit.representative,directions,pole,sum_poles,options,
                                         (options.profile || options.live_profile || options.bath_sector_profile) ? &counters : nullptr) : sum_poles ?
      evaluator.self_nodes_pole_sum(topology,orbit.representative,directions,{{pole[0],pole[1]}},options,
                                   (options.profile || options.live_profile || options.bath_sector_profile) ? &counters : nullptr) :
      evaluator.self_nodes(topology,orbit.representative,directions,pole,options,
                           (options.profile || options.live_profile || options.bath_sector_profile) ? &counters : nullptr);
    if (options.bath_sector_profile)
      for (int bath=0;bath<int(counters.bath_sectors.size());++bath) {
        const auto& sector=counters.bath_sectors[bath];
        const int row=topology*cix.Nbath+bath;
        bath_sector_stats(row,0)+=sector.visits;
        bath_sector_stats(row,1)+=sector.seconds;
        bath_sector_stats(row,2)+=sector.prescore;
        bath_sector_stats(row,3)+=sector.intermediate_norm;
      }
    if (options.bath_sector_profile) add_block_scores(topology,counters);
    r.sigma_topology[topology]+=static_cast<double>(orbit.members.size())*value;
    if (options.live_profile) sample.sigma_norm=double(orbit.members.size())*value.norm();
    // Preserve the existing per-cut output contract. Equality holds after
    // the complete internal-arrow/pole sum, not necessarily for one task.
    if (diagnostics) for (int cut:orbit.members) r.sigma_cut[topology][cut]+=value;
    if (options.profile) {
      profile(topology,0)+=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
      profile(topology,1)+=1;
      profile(topology,7)+=counters.convolution_functions;
      profile(topology,8)+=counters.convolution_batches;
      profile(topology,9)+=counters.root_functions;
      profile(topology,10)+=counters.prepared_unique;
      profile(topology,11)+=counters.prepared_reused;
      profile(topology,12)+=counters.prefix_fallbacks;
      profile(topology,13)+=counters.anchor<0;
      profile(topology,14)+=counters.root_seconds;
      profile(topology,15)+=counters.vertex_seconds;
      profile(topology,16)+=counters.convolution_seconds;
      profile(topology,17)+=counters.pack_seconds;
      profile(topology,18)+=counters.transform_seconds;
      profile(topology,19)+=counters.product_seconds;
      profile(topology,20)+=counters.unpack_seconds;
      profile(topology,21)+=counters.pole_groups;
      profile(topology,22)+=counters.shared_branch_inputs;
      profile(topology,23)+=counters.shared_transform_seconds;
    }
    sample.convolved_functions=counters.convolution_functions;
    sample.root_functions=counters.root_functions;
    return sample;
  },options.reduce_cuts,options.live_profile,live_evaluation);
  const auto reduce=[&](MatrixXd& x) {
#ifdef _MPI
    if (size>1) for (Eigen::Index i=0;i<x.size();) {
      const int count=static_cast<int>(std::min<Eigen::Index>(x.size()-i,INT_MAX));
      MPI_Allreduce(MPI_IN_PLACE,x.data()+i,count,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      i+=count;
    }
#endif
  };
  // Cumulative counters over workers, not just the dispatcher. Memory is
  // the largest per-rank reservation, including future lazy block expansion.
  const auto report_prepared_cache=[&](const char* phase) {
    const auto s=evaluator.prepared_cache_profile();
    MatrixXd counts(1,4); counts<<double(s.hits),double(s.misses),double(s.evictions),double(s.bypasses);
    reduce(counts);
    double peak=double(s.peak_bytes);
#ifdef _MPI
    if (size>1) MPI_Allreduce(MPI_IN_PLACE,&peak,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
#endif
    if (rank==0 && log)
      (*log)<<"FCA prepared cache: phase="<<phase<<" counters=cumulative-all-ranks"
        <<" hits="<<counts(0,0)<<" misses="<<counts(0,1)
        <<" evictions="<<counts(0,2)<<" bypasses="<<counts(0,3)
        <<" max_reserved_mb="<<peak/(1024.*1024.)<<std::endl;
  };
  for (int t=0;t<6;++t) {
    reduce(r.sigma_topology[t]); r.sigma+=r.sigma_topology[t];
    if (diagnostics) for (auto& cut:r.sigma_cut[t]) reduce(cut);
  }
  if (options.bath_sector_profile) {
    reduce(bath_sector_stats);
    if (block_entries) reduce(bath_block_stats);
    if (rank==0) {
      const std::string path="fca_bath_self."+std::to_string(live_evaluation)+".csv";
      std::ofstream out(path,std::ios::out|std::ios::trunc);
      if (out) {
        out<<"topology,bath,visits,worker_s,prefix_score_gross,intermediate_norm_gross\n";
        out<<std::setprecision(17);
        for (int t=0;t<7;++t) for (int bath=0;bath<cix.Nbath;++bath) {
          const int row=t*cix.Nbath+bath;
          if (t<6) out<<fca::names[t]; else out<<"shared-self";
          out<<','<<bath;
          for (int j=0;j<4;++j) out<<','<<bath_sector_stats(row,j);
          out<<'\n';
        }
      }
      if (log) (*log)<<"FCA bath self diagnostic: file="<<path
        <<" scope=streamed-self-chains score=proxy,not-error-bound"
        <<" shared-family=separate-row reference-fallbacks=not-covered"
        <<" status="<<(out ? "written" : "write-failed")<<std::endl;
      if (block_entries) {
        const std::string block_path="fca_bath_backbone."+
          std::to_string(live_evaluation)+".csv";
        std::ofstream block_out(block_path,std::ios::out|std::ios::trunc);
        if (block_out) {
          block_out<<"topology,bath,block_first,block_second,opening_proxy_gross,closing_norm_gross\n";
          block_out<<std::setprecision(17);
          for (int t=0;t<7;++t) for (int bath=0;bath<cix.Nbath;++bath)
            for (int first=0;first<cix.Nblock;++first)
              for (int second=0;second<cix.Nblock;++second) {
                const std::size_t row=t*block_entries+
                  (std::size_t(bath)*cix.Nblock+first)*cix.Nblock+second;
                if (bath_block_stats(row,0)==0. && bath_block_stats(row,1)==0.) continue;
                if (t<6) block_out<<fca::names[t]; else block_out<<"shared-self";
                block_out<<','<<bath<<','<<first<<','<<second<<','
                  <<bath_block_stats(row,0)<<','<<bath_block_stats(row,1)<<'\n';
              }
        }
        if (log) (*log)<<"FCA bath/block diagnostic: file="<<block_path
          <<" status="<<(block_out ? "written" : "write-failed")
          <<" block_key=intermediate-matrix-pair,not-complete-backbone"<<std::endl;
      } else if (log) (*log)<<"FCA bath/block diagnostic skipped: Nbath*Nblock^2 exceeds "
        <<bath_block_score_limit<<" entries per family"<<std::endl;
    }
  }
  // At this boundary every rank has the COMPLETE self-energy, not just its
  // local tasks. The production adapter can now transform and save it before
  // any expensive physical contraction starts. Invoke on every rank so the
  // adapter can also populate its usual distributed return arrays. Only the
  // master writes files. Exceptions must be handled by the MPI-aware caller.
  if (self_complete) self_complete(r);
  report_prepared_cache("self");
  if (!shared_family.empty()) {
    reduce(family_stats);
#ifdef _MPI
    if (size>1) MPI_Allreduce(MPI_IN_PLACE,&family_peak,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
#endif
    if (rank==0 && log)
      (*log)<<"FCA shared prefix: members="<<shared_family.size()
        <<" worker_seconds="<<family_stats(0,0)<<" tasks="<<family_stats(0,1)
        <<" budget_fallbacks="<<family_stats(0,2)
        <<" prefix_convolved_functions="<<family_stats(0,3)
        <<" tail_convolved_functions="<<family_stats(0,4)
        <<" root_functions="<<family_stats(0,5)
        <<" shared_branch_inputs="<<family_stats(0,6)
        <<" max_accounted_mb="<<family_peak/(1024.*1024.)
        <<" (separate from ordinary per-topology profiles; peak includes rejected families)"<<std::endl;
  }
  // Print self profiling immediately: a user stopping after the saved Sigma
  // must not lose these measurements while waiting for physical diagrams.
  if (options.profile) {
    MatrixXd self_profile=profile; reduce(self_profile);
    if (rank==0 && log) for (int t=0;t<6;++t)
      (*log)<<"FCA self profile: topology="<<fca::names[t]
        <<" worker_seconds="<<self_profile(t,0)<<" tasks="<<self_profile(t,1)
        <<" convolved_functions="<<self_profile(t,7)<<" convolution_batches="<<self_profile(t,8)
        <<" root_functions="<<self_profile(t,9)<<" prepared_unique="<<self_profile(t,10)
        <<" prepared_reused="<<self_profile(t,11)<<" prefix_fallbacks="<<self_profile(t,12)
        <<" unstreamed_tasks="<<self_profile(t,13)
        <<" pole_groups="<<self_profile(t,21)<<" shared_branch_inputs="<<self_profile(t,22)
        <<" shared_transform_seconds="<<self_profile(t,23)
        <<" seconds_root_vertex_convolution=["<<self_profile(t,14)<<","<<self_profile(t,15)<<","<<self_profile(t,16)<<"]"
        <<" convolution_seconds_pack_transform_products_unpack=["<<self_profile(t,17)<<","<<self_profile(t,18)
        <<","<<self_profile(t,19)<<","<<self_profile(t,20)<<"]"<<std::endl;
  }
  std::vector<fca::CutOrbit> physical_family,ordinary_physical;
  for (const auto& orbit:physical_orbits) {
    if (options.physical_shared && options.physical_reuse && options.reduce_cuts &&
        options.physical_shared_bytes && physical_family_member(orbit)) physical_family.push_back(orbit);
    else ordinary_physical.push_back(orbit);
  }
  if (physical_family.size()<2) { physical_family.clear(); ordinary_physical=physical_orbits; }
  std::vector<std::vector<fca::CutOrbit>> ef_families;
  if (options.physical_ef_shared && options.physical_reuse && options.reduce_cuts && options.physical_shared_bytes) {
    std::map<std::pair<int,unsigned>,std::vector<fca::CutOrbit>> candidates;
    std::vector<fca::CutOrbit> remaining;
    for (const auto& orbit:ordinary_physical) {
      const auto group=physical_ef_group(orbit);
      if (group.first<0) remaining.push_back(orbit);
      else candidates[group].push_back(orbit);
    }
    for (const auto& [key,family]:candidates) {
      if (family.size()>1) ef_families.push_back(family);
      else remaining.push_back(family.front());
    }
    ordinary_physical=std::move(remaining);
  }
  const std::uint64_t ef_tasks=ef_families.size()*8*triples;
  MatrixXd ef_stats=MatrixXd::Zero(1,8);
  double ef_peak=0;
  const std::uint64_t physical_family_tasks=physical_family.empty() ? 0 : 8*triples;
  MatrixXd physical_family_stats=MatrixXd::Zero(1,18);
  double physical_family_peak=0;
  double physical_transform_peak=0;
  double physical_boundary_peak=0;
  if (physical && rank==0 && log)
    (*log)<<"FCA physical sharing: members="<<physical_family.size()
          <<" cache_mb="<<options.physical_shared_bytes/(1024.*1024.)
          <<" shared_boundary="<<(options.shared_boundary ? "yes" : "no")
          <<" final_sum="<<(options.physical_final_sum ? "yes" : "no")
          <<" boundary_cache_mb="<<options.boundary_cache_bytes/(1024.*1024.)
          <<" family_tasks="<<physical_family_tasks<<std::endl;
  if (physical && rank==0 && log)
    (*log)<<"FCA E/F physical sharing: families="<<ef_families.size()<<" tasks="<<ef_tasks
          <<" presum="<<(options.physical_ef_presum ? "yes" : "no")
          <<" shared_boundary="<<(options.physical_ef_boundary && options.shared_boundary ? "yes" : "no")
          <<" cache_mb="<<options.physical_shared_bytes/(1024.*1024.)
          <<" lifetime=one-crossing-bath-pair"<<std::endl;
  if (physical)
    schedule(ordinary_physical.size()*8*triples+physical_family_tasks+ef_tasks,rank,size,mpi_enabled,log,interval,"physical",[&](std::uint64_t id) {
      const auto start=std::chrono::steady_clock::now();
      LiveTaskSample sample;
      // Larger families first; ordinary work fills workers as they finish.
      if (id<physical_family_tasks) {
        sample.group=7; sample.directions=id/triples;
        const auto pole=poles(id);
        for (int j=0;j<3;++j) sample.pole[j]=pole[j];
        PhysicalSharedProfile counters;
        const auto values=evaluator.physical_family_topologies(physical_family,id/triples,poles(id),options,&counters);
        for (int t=0;t<6;++t) {
          auto& target=r.physical_topology[t];
          target.g+=values[t].g; target.f00+=values[t].f00; target.f01+=values[t].f01;
          if (options.live_profile) {
            sample.g_norm+=values[t].g.norm();
            sample.f00_norm+=values[t].f00.norm();
            sample.f01_norm+=values[t].f01.norm();
          }
        }
        physical_family_stats(0,0)+=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        physical_family_stats(0,1)+=1; physical_family_stats(0,2)+=counters.fallbacks;
        physical_family_stats(0,3)+=counters.convolution_functions;
        physical_family_stats(0,4)+=counters.cache_hits;
        physical_family_stats(0,5)+=counters.left_reuses;
        physical_family_stats(0,6)+=counters.cache_bypasses;
        physical_family_stats(0,7)+=counters.cache_evictions;
        physical_family_stats(0,8)+=counters.shared_branch_inputs;
        physical_family_stats(0,9)+=counters.branch_transform_hits;
        physical_family_stats(0,10)+=counters.boundary_builds;
        physical_family_stats(0,11)+=counters.boundary_hits;
        physical_family_stats(0,12)+=counters.boundary_fallbacks;
        physical_family_stats(0,13)+=counters.boundary_prepare_seconds;
        physical_family_stats(0,14)+=counters.boundary_contract_seconds;
        physical_family_stats(0,15)+=counters.boundary_reference_seconds;
        physical_family_stats(0,16)+=counters.final_sum_groups;
        physical_family_stats(0,17)+=counters.final_convolutions_saved;
        physical_boundary_peak=std::max(physical_boundary_peak,double(counters.boundary_peak_bytes));
        physical_transform_peak=std::max(physical_transform_peak,double(counters.transform_peak_bytes));
        physical_family_peak=std::max(physical_family_peak,double(counters.peak_bytes));
        sample.convolved_functions=counters.convolution_functions;
        return sample;
      }
      id-=physical_family_tasks;
      if (id<ef_tasks) {
        sample.group=8; sample.cut=id/(8*triples); sample.directions=(id/triples)%8;
        const auto pole=poles(id);
        for (int j=0;j<3;++j) sample.pole[j]=pole[j];
        PhysicalProfile counters;
        const auto values=evaluator.physical_ef_family(ef_families[id/(8*triples)],
          (id/triples)%8,poles(id),options,&counters);
        for (int t=0;t<6;++t) {
          r.physical_topology[t].g+=values[t].g;
          r.physical_topology[t].f00+=values[t].f00;
          r.physical_topology[t].f01+=values[t].f01;
          if (options.live_profile) {
            sample.g_norm+=values[t].g.norm();
            sample.f00_norm+=values[t].f00.norm();
            sample.f01_norm+=values[t].f01.norm();
          }
        }
        ef_stats(0,0)+=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        ef_stats(0,1)+=1; ef_stats(0,2)+=counters.convolution_functions;
        ef_stats(0,3)+=counters.convolution_batches; ef_stats(0,4)+=counters.cross_diagram_hits;
        ef_stats(0,5)+=counters.final_convolutions_saved;
        ef_stats(0,6)+=counters.boundary_builds; ef_stats(0,7)+=counters.boundary_hits;
        ef_peak=std::max(ef_peak,double(counters.cross_diagram_peak_bytes));
        sample.convolved_functions=counters.convolution_functions;
        return sample;
      }
      id-=ef_tasks;
      const auto pole=poles(id);
      id/=triples;
      const int internal_directions=id%8; id/=8;
      const auto& orbit=ordinary_physical[id];
      const int topology=orbit.topology,cut=orbit.representative/2;
      int directions=(orbit.representative%2)<<cut,bit=0;
      for (int l=0;l<4;++l) if (l!=cut) directions|=((internal_directions>>bit++)&1)<<l;
      sample.group=topology; sample.cut=cut; sample.directions=directions;
      for (int j=0;j<3;++j) sample.pole[j]=pole[j];
      PhysicalProfile counters;
      const auto value=evaluator.physical_nodes(topology,cut,directions,pole,options,
                                                (options.profile || options.live_profile) ? &counters : nullptr);
      const double multiplicity=orbit.members.size();
      r.physical_topology[topology].g+=multiplicity*value.g;
      r.physical_topology[topology].f00+=multiplicity*value.f00;
      r.physical_topology[topology].f01+=multiplicity*value.f01;
      if (options.live_profile) {
        sample.g_norm=multiplicity*value.g.norm();
        sample.f00_norm=multiplicity*value.f00.norm();
        sample.f01_norm=multiplicity*value.f01.norm();
      }
      if (options.profile) {
        profile(topology,2)+=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        profile(topology,3)+=1;
        profile(topology,4)+=counters.convolution_functions;
        profile(topology,5)+=counters.convolution_batches;
        profile(topology,6)+=counters.cache_hits;
        profile(topology,24)+=counters.shared_branch_inputs;
        profile(topology,25)+=counters.resumed_vertices;
        profile(topology,26)+=counters.cache_evictions;
        cache_peak=std::max(cache_peak,static_cast<double>(counters.cache_peak_bytes));
      }
      sample.convolved_functions=counters.convolution_functions;
      return sample;
    },options.reduce_cuts,options.live_profile,live_evaluation);
  if (physical && !ef_families.empty()) {
    reduce(ef_stats);
#ifdef _MPI
    if (size>1) MPI_Allreduce(MPI_IN_PLACE,&ef_peak,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
#endif
    if (rank==0 && log)
      (*log)<<"FCA E/F shared physical: worker_seconds="<<ef_stats(0,0)
        <<" tasks="<<ef_stats(0,1)<<" convolved_functions="<<ef_stats(0,2)
        <<" convolution_batches="<<ef_stats(0,3)<<" cross_diagram_hits="<<ef_stats(0,4)
        <<" final_convolutions_saved="<<ef_stats(0,5)
        <<" boundary_builds="<<ef_stats(0,6)<<" boundary_hits="<<ef_stats(0,7)
        <<" max_cache_mb="<<ef_peak/(1024.*1024.)<<std::endl;
  }
  if (physical && !physical_family.empty()) {
    reduce(physical_family_stats);
#ifdef _MPI
    if (size>1) MPI_Allreduce(MPI_IN_PLACE,&physical_family_peak,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    if (size>1) MPI_Allreduce(MPI_IN_PLACE,&physical_transform_peak,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    if (size>1) MPI_Allreduce(MPI_IN_PLACE,&physical_boundary_peak,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
#endif
    if (rank==0 && log)
      (*log)<<"FCA shared physical: members="<<physical_family.size()
            <<" worker_seconds="<<physical_family_stats(0,0)<<" tasks="<<physical_family_stats(0,1)
            <<" fallbacks="<<physical_family_stats(0,2)<<" convolved_functions="<<physical_family_stats(0,3)
            <<" prefix_hits="<<physical_family_stats(0,4)<<" complete_left_reuses="<<physical_family_stats(0,5)
            <<" cache_bypasses="<<physical_family_stats(0,6)<<" max_cache_mb="<<physical_family_peak/(1024.*1024.)
            <<" cache_evictions="<<physical_family_stats(0,7)
            <<" shared_branch_inputs="<<physical_family_stats(0,8)
            <<" branch_transform_hits="<<physical_family_stats(0,9)
            <<" transform_cache_peak_mb="<<physical_transform_peak/(1024.*1024.)
            <<" boundary_builds="<<physical_family_stats(0,10)
            <<" boundary_hits="<<physical_family_stats(0,11)
            <<" boundary_fallbacks="<<physical_family_stats(0,12)
            <<" boundary_seconds_prepare_contract_reference=["<<physical_family_stats(0,13)
            <<","<<physical_family_stats(0,14)<<","<<physical_family_stats(0,15)<<"]"
            <<" boundary_cache_peak_mb="<<physical_boundary_peak/(1024.*1024.)
            <<" final_sum_groups="<<physical_family_stats(0,16)
            <<" final_convolutions_saved="<<physical_family_stats(0,17)
            <<" (separate from ordinary topology profiles; cache cap excludes live chains/scratch)"<<std::endl;
  }
  for (int t=0;t<6;++t) {
    if (physical) {
      auto& v=r.physical_topology[t];
      reduce(v.g); reduce(v.f00); reduce(v.f01);
      r.g+=v.g; r.f00+=v.f00; r.f01+=v.f01;
    }
  }
  if (options.profile) {
    reduce(profile);
#ifdef _MPI
    if (size>1) MPI_Allreduce(MPI_IN_PLACE,&cache_peak,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
#endif
    if (rank==0 && log) {
      for (int t=0;t<6;++t)
        (*log)<<"FCA profile: topology="<<fca::names[t]
          <<" self_worker_seconds="<<profile(t,0)<<" self_tasks="<<profile(t,1)
          <<" physical_worker_seconds="<<profile(t,2)<<" physical_tasks="<<profile(t,3)
          <<" physical_convolved_functions="<<profile(t,4)
          <<" physical_convolution_batches="<<profile(t,5)
          <<" physical_prefix_hits="<<profile(t,6)
          <<" physical_shared_branch_inputs="<<profile(t,24)
          <<" physical_resumed_vertices="<<profile(t,25)
          <<" physical_cache_evictions="<<profile(t,26)<<std::endl;
      (*log)<<"FCA profile: physical_prefix_cache_peak_mb="<<cache_peak/(1024.*1024.)<<std::endl;
    }
  }
  if (physical) report_prepared_cache("self+physical");
  return r;
}
} // namespace edmft::fca_dlr
#endif
