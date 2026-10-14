// Discrete third-order (TCA) skeleton-diagram enumeration.
//
// This file deliberately contains no time mesh, DLR basis, AAA fit, or
// numerical integration.  It performs only the CIX/operator walk and returns
// the six pseudo propagators, three hybridization labels, and local matrix
// element of every cyclically unique third-order Phi path.  Keeping this
// stage independent will let the future TCA integrator reorganize the same
// path stream for all six pseudo-self-energy cuts.
#ifndef EDMFT_TCA_PATHS_H
#define EDMFT_TCA_PATHS_H

#include "tca_topology.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdlib>
#include <cstdint>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#ifdef _MPI
#include <mpi.h>
#endif

// A key is exactly the ordered pseudo-propagator ring requested for TCA:
// (p0,p1,p2,p3,p4,p5).  Topologies are kept in separate hash tables, so a
// topology tag is intentionally absent from this key.
struct TCAPseudoKey {
  // This is the full 12-component unsymmetrized identity
  // (p0a,p0b,...,p5a,p5b), packed as six directed pairs.  Sindex labels must
  // never be used here: distinct matrix elements may share one propagated
  // pseudo function but are still distinct diagrams before cyclic reduction.
  std::array<std::uint64_t,6> pseudo{};

  bool operator==(const TCAPseudoKey& other) const noexcept
  {
    return pseudo==other.pseudo;
  }
  bool operator<(const TCAPseudoKey& other) const noexcept
  {
    return pseudo<other.pseudo;
  }
};

struct TCAPseudoKeyHash {
  std::size_t operator()(const TCAPseudoKey& key) const noexcept
  {
    // Hash-combine all bits of the six packed (a,b) pairs.
    std::uint64_t hash=1469598103934665603ULL;
    for (const std::uint64_t value:key.pseudo) {
      hash^=value;
      hash*=1099511628211ULL;
      hash^=(value>>32);
      hash*=1099511628211ULL;
    }
    return static_cast<std::size_t>(hash);
  }
};

// One retained cyclic equivalence class.  dirspin uses the existing solver
// convention: 2*orbital+direction.  It is stored for later hybridization
// evaluation, but it is not part of diagram identity: adjacent pseudo states
// determine a unique operator, hence a unique consistently rotated dirspin.
struct TCAPathRecord {
  std::array<int,6> pseudo{};
  std::array<int,3> dirspin{};
  double matrix_element=0.0;
  std::uint32_t cyclic_copies=1;
  std::uint8_t topology=0;
};

struct TCAEnumeration {
  // Paths remain distributed over MPI ranks.  Numerical contraction should
  // also remain distributed; gathering a third-order catalog on every rank
  // would defeat the memory scaling required for the large f-shell case.
  std::vector<TCAPathRecord> paths;
  unsigned long long global_unique_paths=0;
  unsigned long long global_cyclic_copies=0;
  unsigned long long global_pauli_rejections=0;
  unsigned long long global_candidate_paths=0;
  std::array<unsigned long long,2> global_unique_paths_by_topology{{0,0}};
  std::array<unsigned long long,2> global_cyclic_copies_by_topology{{0,0}};
};

namespace tca_detail {

using Endpoints=edmft::tca::Endpoints;
inline constexpr auto topology_endpoints=edmft::tca::topology_endpoints;

inline constexpr std::array<int,6> all_rotations{{0,1,2,3,4,5}};
inline constexpr std::array<int,2> half_turn_rotations{{0,3}};

template<typename T>
inline std::array<T,6> RotatePseudos(const std::array<T,6>& pseudo,
                                    int rotation)
{
  std::array<T,6> result{};
  for (int vertex=0;vertex<6;++vertex)
    result[vertex]=pseudo[(vertex-rotation+6)%6];
  return result;
}

struct RotatedLines {
  Endpoints endpoints{};
  std::array<int,3> dirspin{};
};

inline RotatedLines RotateLines(int topology,
                                const std::array<int,3>& dirspin,
                                int rotation)
{
  if (topology<0 || topology>=2)
    throw std::out_of_range("TCA topology index must be zero or one");

  struct LabeledLine { int first=0,second=0,dirspin=0; };
  std::array<LabeledLine,3> lines{};
  for (int line=0;line<3;++line) {
    int first=(topology_endpoints[topology][line][0]+rotation)%6;
    int second=(topology_endpoints[topology][line][1]+rotation)%6;
    int ds=dirspin[line];
    if (first>second) {
      std::swap(first,second);
      // dirspin denotes the operator at the second endpoint.  Swapping the
      // endpoints therefore exchanges creation and annihilation.
      ds^=1;
    }
    lines[line]={first,second,ds};
  }
  RotatedLines result;
  std::array<bool,3> used{{false,false,false}};
  bool same_topology=true;
  for (int target=0;target<3;++target) {
    int source=-1;
    for (int line=0;line<3;++line)
      if (!used[line] &&
          std::array<int,2>{{lines[line].first,lines[line].second}}==
            topology_endpoints[topology][target]) {
        source=line;
        break;
      }
    if (source<0) {
      same_topology=false;
      break;
    }
    used[source]=true;
    result.endpoints[target]=topology_endpoints[topology][target];
    result.dirspin[target]=lines[source].dirspin;
  }

  if (!same_topology) {
    // RotationPreservesTopology needs a deterministic nonmatching result for
    // rotations outside the topology's automorphism group.
    std::sort(lines.begin(),lines.end(),[](const auto& a,const auto& b) {
      return std::tie(a.first,a.second)<std::tie(b.first,b.second);
    });
    for (int line=0;line<3;++line) {
      result.endpoints[line]={lines[line].first,lines[line].second};
      result.dirspin[line]=lines[line].dirspin;
    }
  }
  return result;
}

inline bool RotationPreservesTopology(int topology,int rotation)
{
  const std::array<int,3> dummy_dirspin{{0,0,0}};
  return RotateLines(topology,dummy_dirspin,rotation).endpoints==
    topology_endpoints[topology];
}

struct CanonicalPath {
  TCAPseudoKey key;
  std::array<int,3> dirspin{};
  int rotation=0;
};

inline CanonicalPath Canonicalize(int topology,
                                  const std::array<std::uint64_t,6>& pseudo,
                                  const std::array<int,3>& dirspin)
{
  CanonicalPath best{{pseudo},dirspin,0};
  const auto consider=[&](int rotation) {
    if (!RotationPreservesTopology(topology,rotation))
      throw std::logic_error("invalid rotation in TCA topology orbit");
    const auto rotated_pseudo=RotatePseudos(pseudo,rotation);
    const auto rotated_lines=RotateLines(topology,dirspin,rotation);
    if (rotated_pseudo<best.key.pseudo) {
      best.key.pseudo=rotated_pseudo;
      best.dirspin=rotated_lines.dirspin;
      best.rotation=rotation;
    } else if (rotated_pseudo==best.key.pseudo &&
               rotated_lines.dirspin!=best.dirspin) {
      // Equal adjacent pseudo transitions must imply equal operators.  A
      // violation means either a CIX ambiguity or an error in rotation
      // bookkeeping; silently identifying the paths would be unsafe.
      throw std::logic_error(
        "cyclically identical TCA pseudo ring has inconsistent dirspin");
    }
  };

  if (topology==0) {
    consider(3);
  } else {
    for (int rotation=1;rotation<6;++rotation) consider(rotation);
  }
  return best;
}

// Outer MPI scheduling must keep rotated dirspin copies on one rank, because
// their pseudo keys can become identical.  This key is used only to assign
// work; diagram uniqueness itself is still decided solely by TCAPseudoKey.
struct DirspinOrbitKey {
  int topology=0;
  std::array<int,3> dirspin{};
  bool operator<(const DirspinOrbitKey& other) const noexcept
  {
    return std::tie(topology,dirspin)<
      std::tie(other.topology,other.dirspin);
  }
};

inline DirspinOrbitKey CanonicalDirspinOrbit(
    int topology,const std::array<int,3>& dirspin)
{
  DirspinOrbitKey best{topology,dirspin};
  if (topology==0) {
    best.dirspin=std::min(best.dirspin,RotateLines(0,dirspin,3).dirspin);
  } else {
    for (int rotation=1;rotation<6;++rotation)
      best.dirspin=std::min(
        best.dirspin,RotateLines(1,dirspin,rotation).dirspin);
  }
  return best;
}

inline int IntegerSign(int value)
{
  // Zero is the ordinary first bath channel, not a vanishing sign.  Negative
  // values encode symmetry-related channels with a relative minus.
  return value>=0 ? +1 : -1;
}

inline std::uint64_t RawPseudoIdentity(int global_a,int global_b)
{
  return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(global_a))
          <<32) |
         static_cast<std::uint32_t>(global_b);
}

} // namespace tca_detail

// Enumerate the two third-order skeleton topologies.  This is intentionally
// an independent entry point: merely compiling it has no runtime cost, and it
// will be connected to the future TCA evaluator only after its path counts and
// matrix elements have been checked against the archived MC references.
inline TCAEnumeration EnumerateTCAPaths(const Cix& cix,
                                        const std::vector<int>& G_importance,
                                        int TCA_cut_sum,
                                        MPIcls& mpi)
{
  using namespace tca_detail;
  using Clock=std::chrono::steady_clock;

  if (static_cast<int>(G_importance.size())<cix.Npseudo_unique)
    throw std::invalid_argument(
      "EnumerateTCAPaths requires one importance value per pseudo propagator");

  Diag2Matrix diag2matrix(6,cix);
  blitz::Array<int,2> diag(3,2);
  blitz::Array<int,1> dirspin_array(3);

  struct Task {
    int topology=0;
    std::array<int,3> dirspin{};
    long double cost=0.0L;
    std::size_t orbit=0;
  };
  struct Orbit {
    long double cost=0.0L;
    std::vector<std::size_t> tasks;
  };

  const auto set_topology=[&](int topology) {
    for (int line=0;line<3;++line) {
      diag(line,0)=topology_endpoints[topology][line][0];
      diag(line,1)=topology_endpoints[topology][line][1];
    }
  };
  const auto set_dirspin=[&](const std::array<int,3>& ds) {
    for (int line=0;line<3;++line) dirspin_array(line)=ds[line];
  };

  const auto scheduling_start=Clock::now();
  std::vector<Task> tasks;
  std::vector<Orbit> orbits;
  std::map<DirspinOrbitKey,std::size_t> orbit_lookup;
  unsigned long long pauli_rejections=0;

  // The requested external order is (i,k,j), not (i,j,k).  Each variable is
  // a full dirspin label in [0,2*Nbath), including its direction bit.
  for (int topology=0;topology<2;++topology) {
    set_topology(topology);
    for (int i=0;i<2*cix.Nbath;++i)
      for (int j=0;j<2*cix.Nbath;++j)
        for (int k=0;k<2*cix.Nbath;++k) {
          const std::array<int,3> ds{{i,k,j}};
          set_dirspin(ds);
          if (PPconfigs2(diag,dirspin_array,cix)==0) {
            ++pauli_rejections;
            continue;
          }

          auto& allowed=diag2matrix.Cmp_matrix(diag,dirspin_array,cix);
          long double cost=0.0L;
          for (const auto& [ib,alw]:allowed) {
            long double block_cost=1.0L;
            for (int vertex=0;vertex<6;++vertex)
              block_cost*=alw[vertex].extent(0);
            cost+=block_cost;
          }
          if (cost==0.0L) continue;

          const auto orbit_key=CanonicalDirspinOrbit(topology,ds);
          const auto [it,inserted]=
            orbit_lookup.try_emplace(orbit_key,orbits.size());
          if (inserted) orbits.push_back({});
          const std::size_t task_index=tasks.size();
          tasks.push_back({topology,ds,cost,it->second});
          orbits[it->second].cost+=cost;
          orbits[it->second].tasks.push_back(task_index);
        }
  }

  // Greedy whole-orbit scheduling mirrors OCA's bath-pair scheduler, but the
  // indivisible unit is now a cyclic dirspin orbit.  Equivalent rotations can
  // therefore never reach different rank-local uniqueness tables.
  std::vector<int> task_owner(tasks.size(),0);
  std::vector<long double> rank_load(std::max(mpi.size,1),0.0L);
  std::vector<std::size_t> sorted_orbits(orbits.size());
  std::iota(sorted_orbits.begin(),sorted_orbits.end(),0);
  std::stable_sort(sorted_orbits.begin(),sorted_orbits.end(),
                   [&](std::size_t a,std::size_t b) {
                     return orbits[a].cost>orbits[b].cost;
                   });
  for (const std::size_t orbit_index:sorted_orbits) {
    const int owner=static_cast<int>(std::min_element(
      rank_load.begin(),rank_load.end())-rank_load.begin());
    for (const std::size_t task_index:orbits[orbit_index].tasks)
      task_owner[task_index]=owner;
    rank_load[owner]+=orbits[orbit_index].cost;
  }
  const auto scheduling_end=Clock::now();
  const long double estimated_total=
    std::accumulate(rank_load.begin(),rank_load.end(),0.0L);
  const long double estimated_largest_rank=
    rank_load.empty() ? 0.0L :
    *std::max_element(rank_load.begin(),rank_load.end());

  if (mpi.rank==mpi.master) {
    const auto [minimum,maximum]=
      std::minmax_element(rank_load.begin(),rank_load.end());
    (*mpi.out) << "TCA enumeration scheduler: tasks=" << tasks.size()
               << " cyclic_dirspin_orbits=" << orbits.size()
               << " ranks=" << mpi.size
               << " estimated_load_range=[" << static_cast<double>(*minimum)
               << "," << static_cast<double>(*maximum) << "]"
               << " max_over_mean="
               << static_cast<double>(estimated_total>0 ?
                    *maximum/(estimated_total/std::max(mpi.size,1)) : 0)
               << " scheduling_seconds="
               << std::chrono::duration<double>(
                    scheduling_end-scheduling_start).count()
               << std::endl;
  }

  // An explicit catalog stores both a 48-byte raw key (in the uniqueness
  // hash) and a symmetry-reduced evaluation record for every retained path.
  // It is therefore unsafe for cluster TCA: the SC four-site test expands to
  // billions of scalar paths even though its atomic blocks are small.  When
  // importance cutting is inactive, reject such jobs before allocating the
  // hash tables.  Large TCA must instead contract block matrices first and
  // apply DLR operations to the resulting matrix-valued functions.
  const int maximum_importance=G_importance.empty() ? 0 :
    *std::max_element(G_importance.begin(),G_importance.end());
  const bool importance_cut_inactive=
    TCA_cut_sum>=6*maximum_importance;
  constexpr long double safe_paths_per_rank=2000000.0L;
  if (importance_cut_inactive &&
      estimated_largest_rank>safe_paths_per_rank) {
    throw std::length_error(
      "EnumerateTCAPaths: explicit raw-path materialization is unsafe "
      "(estimated largest-rank Cartesian work="+
      std::to_string(static_cast<unsigned long long>(
        estimated_largest_rank))+ "). Use the block-matrix TCA contraction; "
      "Sindex labels may reduce propagated functions but must not define "
      "raw cyclic diagram identity.");
  }

  const auto enumeration_start=Clock::now();
  TCAEnumeration result;
  std::array<std::unordered_map<TCAPseudoKey,std::size_t,TCAPseudoKeyHash>,2>
    seen;
  unsigned long long local_candidates=0;
  unsigned long long local_copies=0;
  std::array<unsigned long long,2> local_unique_by_topology{{0,0}};
  std::array<unsigned long long,2> local_copies_by_topology{{0,0}};

  for (std::size_t task_index=0;task_index<tasks.size();++task_index) {
    if (task_owner[task_index]!=mpi.rank) continue;
    const Task& task=tasks[task_index];
    set_topology(task.topology);
    set_dirspin(task.dirspin);
    auto& allowed=diag2matrix.Cmp_matrix(diag,dirspin_array,cix);

    const int bath_sign=
      IntegerSign(cix.bath_sigind0(task.dirspin[0]/2))*
      IntegerSign(cix.bath_sigind0(task.dirspin[1]/2))*
      IntegerSign(cix.bath_sigind0(task.dirspin[2]/2));

    for (const auto& [ib,alw]:allowed) {
      bool empty=false;
      for (int vertex=0;vertex<6;++vertex)
        empty=empty || alw[vertex].extent(0)==0;
      if (empty) continue;

      // Block immediately after each operator.  ib is both the block before
      // op(0) and after op(5), closing the pseudo-particle trace.
      std::array<int,6> ibs{};
      ibs[0]=cix.Fbc(diag2matrix.op(0),ib);
      for (int vertex=1;vertex<5;++vertex)
        ibs[vertex]=cix.Fbc(diag2matrix.op(vertex),ibs[vertex-1]);
      ibs[5]=ib;

      for (int j0=0;j0<alw[0].extent(0);++j0) {
        const int ga0=std::abs(alw[0](j0,1));
        const int gb5=std::abs(alw[0](j0,0));
        for (int j1=0;j1<alw[1].extent(0);++j1) {
          const int gb0=std::abs(alw[1](j1,0));
          // A zero Sindex entry means that this propagator matrix element is
          // forbidden.  This is only early viability pruning; the nonzero
          // symmetry label is deliberately not used for diagram identity.
          if (cix.Sgm_symmetry[ibs[0]](
                cix.inBlock(ga0),cix.inBlock(gb0))==0) continue;
          const int ga1=std::abs(alw[1](j1,1));
          for (int j2=0;j2<alw[2].extent(0);++j2) {
            const int gb1=std::abs(alw[2](j2,0));
            if (cix.Sgm_symmetry[ibs[1]](
                  cix.inBlock(ga1),cix.inBlock(gb1))==0) continue;
            const int ga2=std::abs(alw[2](j2,1));
            for (int j3=0;j3<alw[3].extent(0);++j3) {
              const int gb2=std::abs(alw[3](j3,0));
              if (cix.Sgm_symmetry[ibs[2]](
                    cix.inBlock(ga2),cix.inBlock(gb2))==0) continue;
              const int ga3=std::abs(alw[3](j3,1));
              for (int j4=0;j4<alw[4].extent(0);++j4) {
                const int gb3=std::abs(alw[4](j4,0));
                if (cix.Sgm_symmetry[ibs[3]](
                      cix.inBlock(ga3),cix.inBlock(gb3))==0) continue;
                const int ga4=std::abs(alw[4](j4,1));
                for (int j5=0;j5<alw[5].extent(0);++j5) {
                  const int gb4=std::abs(alw[5](j5,0));
                  if (cix.Sgm_symmetry[ibs[4]](
                        cix.inBlock(ga4),cix.inBlock(gb4))==0) continue;
                  const int ga5=std::abs(alw[5](j5,1));
                  if (cix.Sgm_symmetry[ibs[5]](
                        cix.inBlock(ga5),cix.inBlock(gb5))==0) continue;

                  ++local_candidates;
                  const std::array<std::uint64_t,6> raw_pseudo{{
                    RawPseudoIdentity(ga0,gb0),RawPseudoIdentity(ga1,gb1),
                    RawPseudoIdentity(ga2,gb2),RawPseudoIdentity(ga3,gb3),
                    RawPseudoIdentity(ga4,gb4),RawPseudoIdentity(ga5,gb5)}};

                  // Resolve cyclic identity with all twelve unsymmetrized
                  // state indices.  Only a genuinely new raw diagram may be
                  // reduced to the six propagated Sindex labels below.
                  const CanonicalPath canonical=
                    Canonicalize(task.topology,raw_pseudo,task.dirspin);
                  auto& topology_seen=seen[task.topology];
                  const auto where=topology_seen.find(canonical.key);
                  if (where!=topology_seen.end()) {
                    TCAPathRecord& previous=result.paths[where->second];
                    if (previous.dirspin!=canonical.dirspin)
                      throw std::logic_error(
                        "identical TCA pseudo key produced different dirspin");
                    ++previous.cyclic_copies;
                    ++local_copies;
                    ++local_copies_by_topology[task.topology];
                    continue;
                  }

                  const int q0=cix.Sgm_symmetry[ibs[0]](
                    cix.inBlock(ga0),cix.inBlock(gb0));
                  const int q1=cix.Sgm_symmetry[ibs[1]](
                    cix.inBlock(ga1),cix.inBlock(gb1));
                  const int q2=cix.Sgm_symmetry[ibs[2]](
                    cix.inBlock(ga2),cix.inBlock(gb2));
                  const int q3=cix.Sgm_symmetry[ibs[3]](
                    cix.inBlock(ga3),cix.inBlock(gb3));
                  const int q4=cix.Sgm_symmetry[ibs[4]](
                    cix.inBlock(ga4),cix.inBlock(gb4));
                  const int q5=cix.Sgm_symmetry[ibs[5]](
                    cix.inBlock(ga5),cix.inBlock(gb5));
                  const std::array<int,6> sym_pseudo{{
                    std::abs(q0)-1,std::abs(q1)-1,std::abs(q2)-1,
                    std::abs(q3)-1,std::abs(q4)-1,std::abs(q5)-1}};

                  int importance=0;
                  for (const int p:sym_pseudo) importance+=G_importance[p];
                  if (importance>TCA_cut_sum) continue;

                  topology_seen.emplace(canonical.key,result.paths.size());

                  const int matrix_sign=
                    IntegerSign(alw[0](j0,1))*
                    IntegerSign(alw[1](j1,1))*
                    IntegerSign(alw[2](j2,1))*
                    IntegerSign(alw[3](j3,1))*
                    IntegerSign(alw[4](j4,1))*
                    IntegerSign(alw[5](j5,1));
                  const int propagator_sign=
                    IntegerSign(q0)*IntegerSign(q1)*IntegerSign(q2)*
                    IntegerSign(q3)*IntegerSign(q4)*IntegerSign(q5);

                  TCAPathRecord record;
                  record.pseudo=RotatePseudos(sym_pseudo,canonical.rotation);
                  record.dirspin=canonical.dirspin;
                  record.topology=static_cast<std::uint8_t>(task.topology);
                  record.matrix_element=
                    static_cast<double>(bath_sign*matrix_sign*propagator_sign);
                  result.paths.push_back(record);
                  ++local_unique_by_topology[task.topology];
                }
              }
            }
          }
        }
      }
    }
  }

  unsigned long long local_unique=result.paths.size();
  unsigned long long global_values[4]={
    local_unique,local_copies,pauli_rejections,local_candidates};
#ifdef _MPI
  if (mpi.size>1) {
    unsigned long long reduced[4]{};
    MPI_Allreduce(global_values,reduced,4,MPI_UNSIGNED_LONG_LONG,MPI_SUM,
                  MPI_COMM_WORLD);
    std::copy(std::begin(reduced),std::end(reduced),global_values);
    // Every rank performed the same inexpensive task-discovery pass, so its
    // Pauli rejection count must not be summed over ranks.
    global_values[2]/=mpi.size;
  }
#endif
  result.global_unique_paths=global_values[0];
  result.global_cyclic_copies=global_values[1];
  result.global_pauli_rejections=global_values[2];
  result.global_candidate_paths=global_values[3];

  unsigned long long topology_values[4]{
    local_unique_by_topology[0],local_unique_by_topology[1],
    local_copies_by_topology[0],local_copies_by_topology[1]};
#ifdef _MPI
  if (mpi.size>1) {
    unsigned long long reduced[4]{};
    MPI_Allreduce(topology_values,reduced,4,MPI_UNSIGNED_LONG_LONG,MPI_SUM,
                  MPI_COMM_WORLD);
    std::copy(std::begin(reduced),std::end(reduced),topology_values);
  }
#endif
  result.global_unique_paths_by_topology=
    {{topology_values[0],topology_values[1]}};
  result.global_cyclic_copies_by_topology=
    {{topology_values[2],topology_values[3]}};

  const double enumeration_seconds=
    std::chrono::duration<double>(Clock::now()-enumeration_start).count();
  double maximum_seconds=enumeration_seconds;
#ifdef _MPI
  if (mpi.size>1)
    MPI_Allreduce(&enumeration_seconds,&maximum_seconds,1,MPI_DOUBLE,MPI_MAX,
                  MPI_COMM_WORLD);
#endif
  if (mpi.rank==mpi.master)
    (*mpi.out) << "TCA enumeration MPI: global_unique_paths="
               << result.global_unique_paths
               << " cyclic_copies_removed=" << result.global_cyclic_copies
               << " Pauli_rejected_dirspin="
               << result.global_pauli_rejections
               << " candidate_paths=" << result.global_candidate_paths
               << " unique_by_topology=["
               << result.global_unique_paths_by_topology[0] << ","
               << result.global_unique_paths_by_topology[1] << "]"
               << " copies_by_topology=["
               << result.global_cyclic_copies_by_topology[0] << ","
               << result.global_cyclic_copies_by_topology[1] << "]"
               << " enumeration_max_seconds=" << maximum_seconds
               << std::endl;

  return result;
}

#endif // EDMFT_TCA_PATHS_H
