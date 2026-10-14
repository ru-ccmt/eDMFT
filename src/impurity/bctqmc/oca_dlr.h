// Compact DLR/AAA evaluation plans for imaginary-axis OCA diagrams.
#ifndef EDMFT_OCA_DLR_H
#define EDMFT_OCA_DLR_H

#include "imag_axis_dlr.h"
#include "oca_dlr_types.h"

#include <blitz/array.h>
#include "Eigen/Sparse"

#include <array>
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

#ifdef _MPI
#include <mpi.h>
#endif

namespace edmft::oca_dlr {

using Eigen::MatrixXd;
using Eigen::SparseMatrix;
using Eigen::Triplet;
using Eigen::VectorXd;

inline std::uint64_t pair_key(int first, int second)
{
  // Preserve ordering while packing two signed 32-bit indices into a hash
  // key.  Pair (a,b) and pair (b,a) are generally different convolutions.
  return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(first))<<32) |
         static_cast<std::uint32_t>(second);
}

// First compact prototype: evaluate only the p3 self-energy cut.  It is kept
// as a small readable reference implementation and as a diagnostic against
// the general all-four-cut SelfEnergyPlan below.
class P3Plan {
public:
  void build(const std::vector<OCAPathRecord>& paths,
             int pseudo_count,
             int channel_count,
             const std::vector<int>& importance,
             int cut_sum)
  {
    pseudo_count_=pseudo_count;
    channel_count_=channel_count;
    importance_=importance;
    cut_sum_=cut_sum;
    inner_pairs_.clear(); outer_pairs_.clear();
    std::unordered_map<std::uint64_t,int> inner_ids,outer_ids;

    // For the p3 cut the remaining pseudo ring separates naturally into an
    // inner pair (p0,p1) and an outer pair (p2,p3).  Deduplicating these pairs
    // is the main cluster-DMFT optimization: thousands of atomic paths reuse
    // the same numerical convolution.
    for (const auto& path:paths) {
      if (path.weights[0]==0.0) continue;
      const auto ikey=pair_key(path.pseudo[0],path.pseudo[1]);
      const auto okey=pair_key(path.pseudo[2],path.pseudo[3]);
      if (!inner_ids.contains(ikey)) {
        const int id=inner_pairs_.size();
        inner_ids.emplace(ikey,id);
        inner_pairs_.push_back({path.pseudo[0],path.pseudo[1]});
      }
      if (!outer_ids.contains(okey)) {
        const int id=outer_pairs_.size();
        outer_ids.emplace(okey,id);
        outer_pairs_.push_back({path.pseudo[2],path.pseudo[3]});
      }
    }

    // C_(s0,s1)[outer,inner] performs the complete local-state and operator
    // matrix-element sum after the numerical functions have been computed.
    std::vector<std::vector<Triplet<double>>> triplets(channel_count*channel_count);
    nonzero_path_count_=0;
    for (const auto& path:paths) {
      if (path.weights[0]==0.0) continue;
      if (path.s0<0 || path.s0>=channel_count || path.s1<0 || path.s1>=channel_count)
        throw std::out_of_range("OCA bath channel in DLR plan");
      const int col=inner_ids.at(pair_key(path.pseudo[0],path.pseudo[1]));
      const int row=outer_ids.at(pair_key(path.pseudo[2],path.pseudo[3]));
      // The explicit minus is the crossing sign used by the reference OCA.
      triplets[path.s0*channel_count+path.s1].emplace_back(row,col,-path.weights[0]);
      ++nonzero_path_count_;
    }

    contractions_.clear();
    contractions_.resize(channel_count*channel_count);
    for (int i=0; i<static_cast<int>(contractions_.size()); ++i) {
      contractions_[i].resize(outer_pairs_.size(),inner_pairs_.size());
      contractions_[i].setFromTriplets(triplets[i].begin(),triplets[i].end());
      contractions_[i].makeCompressed();
    }

    // Batched convolution shares its left operand.  Grouping ordered pairs by
    // the second pseudo index gives precisely that common operand with the
    // convolution convention used by FermionicBasis.
    inner_by_second_.assign(pseudo_count,{});
    for (int id=0; id<static_cast<int>(inner_pairs_.size()); ++id)
      inner_by_second_.at(inner_pairs_[id][1]).push_back(id);
    outer_by_first_.assign(pseudo_count,{});
    for (int id=0; id<static_cast<int>(outer_pairs_.size()); ++id)
      outer_by_first_.at(outer_pairs_[id][0]).push_back(id);

    all_pseudos_.resize(pseudo_count);
    for (int p=0; p<pseudo_count; ++p) all_pseudos_[p]=p;
  }

  bool matches(int pseudo_count,
               int channel_count,
               const std::vector<int>& importance,
               int cut_sum) const
  {
    return pseudo_count_==pseudo_count && channel_count_==channel_count &&
           importance_==importance && cut_sum_==cut_sum && !inner_pairs_.empty();
  }

  int inner_pair_count() const { return inner_pairs_.size(); }
  int outer_pair_count() const { return outer_pairs_.size(); }
  std::size_t nonzero_path_count() const { return nonzero_path_count_; }

  MatrixXd evaluate(const imag_dlr::Workspace& workspace,
                    const std::vector<double>& pseudo_degeneracy,
                    const blitz::Array<double,1>& output_tau) const
  {
    if (!workspace.has_pseudos()) throw std::logic_error("OCA DLR evaluation requires pseudo-particle coefficients");
    if (workspace.pseudo_count()!=pseudo_count_ || workspace.channel_count()!=channel_count_)
      throw std::invalid_argument("OCA DLR plan/workspace dimension mismatch");
    if (static_cast<int>(pseudo_degeneracy.size())!=pseudo_count_)
      throw std::invalid_argument("OCA DLR pseudo degeneracy size mismatch");

    const auto& basis=workspace.basis();
    const int r=workspace.rank();
    MatrixXd sigma_nodes=MatrixXd::Zero(r,pseudo_count_);

    // Sum the expanded hybridization pole-by-pole.  At no point do we sample
    // or integrate on the legacy uniform OCA mesh.
    for (int l=0; l<workspace.pole_count(); ++l) {
      const bool nonpositive=workspace.poles()(l)<=0.0;
      const int kernel_sign=nonpositive ? -1 : +1;
      const MatrixXd tilted=workspace.tilted_pseudo_coefficients(all_pseudos_,l,kernel_sign);

      MatrixXd b_values(r,inner_pairs_.size());
      for (int second=0; second<pseudo_count_; ++second) {
        const auto& ids=inner_by_second_[second];
        if (ids.empty()) continue;
        MatrixXd rhs(r,ids.size());
        VectorXd lhs;
        if (nonpositive) {
          lhs=workspace.pseudo_coefficients().col(second);
          for (int j=0; j<static_cast<int>(ids.size()); ++j)
            rhs.col(j)=tilted.col(inner_pairs_[ids[j]][0]);
        } else {
          lhs=tilted.col(second);
          for (int j=0; j<static_cast<int>(ids.size()); ++j)
            rhs.col(j)=workspace.pseudo_coefficients().col(inner_pairs_[ids[j]][0]);
        }
        const MatrixXd values=basis.time_ordered_convolution_values(lhs,rhs);
        for (int j=0; j<static_cast<int>(ids.size()); ++j) b_values.col(ids[j])=values.col(j);
      }

      MatrixXd h_transposed=MatrixXd::Zero(outer_pairs_.size(),r);
      const MatrixXd b_transposed=b_values.transpose();
      for (int s0=0; s0<channel_count_; ++s0) {
        const double residue=workspace.residues()(s0,l);
        if (residue==0.0) continue;
        for (int s1=0; s1<channel_count_; ++s1) {
          const auto& contraction=contractions_[s0*channel_count_+s1];
          if (contraction.nonZeros()==0) continue;
          MatrixXd contribution=contraction*b_transposed;
          contribution.array().rowwise() *= workspace.delta_plus().col(s1).transpose().array();
          h_transposed.noalias() += residue*contribution;
        }
      }
      const MatrixXd h_values=h_transposed.transpose();
      const MatrixXd h_coefficients=basis.values_to_coefficients(h_values);

      VectorXd external(r);
      if (nonpositive)
        external=workspace.kernel_plus().col(l)/workspace.kernel_zero_minus()(l);
      else
        external.setConstant(1.0/workspace.kernel_zero_plus()(l));
      for (int outer=0; outer<pseudo_count_; ++outer) {
        const auto& ids=outer_by_first_[outer];
        if (ids.empty()) continue;
        MatrixXd rhs(r,ids.size());
        for (int j=0; j<static_cast<int>(ids.size()); ++j) rhs.col(j)=h_coefficients.col(ids[j]);
        const VectorXd lhs=nonpositive
          ? VectorXd(workspace.pseudo_coefficients().col(outer))
          : VectorXd(tilted.col(outer));
        const MatrixXd values=basis.time_ordered_convolution_values(lhs,rhs);
        for (int j=0; j<static_cast<int>(ids.size()); ++j) {
          const int output=outer_pairs_[ids[j]][1];
          sigma_nodes.col(output).array() +=
            external.array()*values.col(j).array()/pseudo_degeneracy[output];
        }
      }
    }

    const MatrixXd sigma_coefficients=basis.values_to_coefficients(sigma_nodes);
    MatrixXd result(pseudo_count_,output_tau.extent(0));
    for (int it=0; it<output_tau.extent(0); ++it)
      result.col(it)=basis.evaluate(output_tau(it),sigma_coefficients).transpose();
    return result;
  }

private:
  int pseudo_count_=0,channel_count_=0,cut_sum_=0;
  std::size_t nonzero_path_count_=0;
  std::vector<int> importance_,all_pseudos_;
  std::vector<std::array<int,2>> inner_pairs_,outer_pairs_;
  std::vector<std::vector<int>> inner_by_second_,outer_by_first_;
  std::vector<SparseMatrix<double>> contractions_;
};

// General compact OCA pseudo self-energy.  Cutting each of the four pseudo
// lines of Phi gives Sigma_p0 ... Sigma_p3.  The topology is fixed, while all
// model dependence is isolated in sparse contraction matrices built from the
// legacy enumerator's OCAPathRecords.
class SelfEnergyPlan {
  struct PairCatalog {
    // Unique ordered pseudo pairs G_first(x)G_second(t-x).
    std::vector<std::array<int,2>> pairs;
    // Pair ids grouped by their second member for batched convolutions.
    std::vector<std::vector<int>> by_second;
  };
  struct CutPlan {
    int sigma=1;                 // expanded line sees sigma*omega_l
    int local_orientation=1;     // unexpanded Delta(t) (+1) or Delta(beta-t) (-1)
    // (w,q): w is convolved with the already contracted inner object; q is
    // the external/output pseudo self-energy index.
    std::vector<std::array<int,2>> outer_pairs;
    std::vector<std::vector<int>> outer_by_first;
    std::vector<SparseMatrix<double>> contractions;
  };

public:
  void build(const std::vector<OCAPathRecord>& paths,
             int pseudo_count,
             int channel_count,
             const std::vector<int>& importance,
             int cut_sum)
  {
    pseudo_count_=pseudo_count;
    channel_count_=channel_count;
    importance_=importance;
    cut_sum_=cut_sum;
    nonzero_path_count_=0;
    all_pseudos_.resize(pseudo_count);
    for (int p=0; p<pseudo_count; ++p) all_pseudos_[p]=p;

    // One table row per cut.  Columns are
    //   q,u,v,w, expanded-bath selector, local-bath selector,
    //   effective pole sign sigma, local-Delta orientation.
    // q is the cut/output line; (u,v) is the cached inner convolution; w is
    // the remaining propagator in the outer convolution.  Bath selectors 0
    // and 1 mean path.s0 and path.s1.  Written out, the four cuts are:
    //
    //   cut p0: inner(p1,p2), outer(p3,p0), expand s1, local s0(t)
    //   cut p1: inner(p2,p3), outer(p0,p1), expand s0, local s1(beta-t)
    //   cut p2: inner(p3,p0), outer(p1,p2), expand s1, local s0(beta-t)
    //   cut p3: inner(p0,p1), outer(p2,p3), expand s0, local s1(t).
    //
    // Centralizing this table is safer than maintaining four nearly equal
    // formulas with subtly permuted indices.
    constexpr int spec[4][8]={
      {0,1,2,3,1,0,-1,+1},
      {1,2,3,0,0,1,-1,-1},
      {2,3,0,1,1,0,+1,-1},
      {3,0,1,2,0,1,+1,+1}
    };

    // Only two numerical pair families are necessary.  They correspond to
    // effective +omega and -omega flow and are shared between opposite cuts.
    //
    // Include physical-only Bulla paths in these catalogs even though their
    // ordinary OCA weight is zero.  This makes the pair numbering identical
    // to PhysicalPlan whenever physical observables are requested.  A single
    // pole-major evaluator can then compute each pair value once and feed both
    // plans without a permutation or a second convolution pass.
    std::array<std::unordered_map<std::uint64_t,int>,2> pair_ids;
    for (auto& catalog:pair_catalogs_) { catalog.pairs.clear(); catalog.by_second.clear(); }
    const auto insert_pair=[&](int catalog_index,int first,int second) {
      const auto key=pair_key(first,second);
      auto [it,inserted]=pair_ids[catalog_index].emplace(
        key,static_cast<int>(pair_catalogs_[catalog_index].pairs.size()));
      if (inserted) pair_catalogs_[catalog_index].pairs.push_back({first,second});
      return it->second;
    };
    for (const auto& path:paths) {
      const bool need_s1=path.weights[0]!=0.0 || path.weights[1]!=0.0 ||
                         path.weights[3]!=0.0;
      const bool need_s0=path.weights[0]!=0.0 || path.weights[2]!=0.0 ||
                         path.weights[4]!=0.0;
      // Keep precisely the same insertion order as PhysicalPlan::build.
      // Equal insertion order means equal pair ids, not merely equal sets.
      if (need_s1) {
        insert_pair(1,path.pseudo[2],path.pseudo[3]);
        insert_pair(0,path.pseudo[0],path.pseudo[1]);
      }
      if (need_s0) {
        insert_pair(0,path.pseudo[3],path.pseudo[0]);
        insert_pair(1,path.pseudo[1],path.pseudo[2]);
      }
      if (path.weights[0]!=0.0) ++nonzero_path_count_;
    }
    for (int ci=0; ci<2; ++ci) {
      auto& catalog=pair_catalogs_[ci];
      // Group equal second pseudo indices into contiguous pair-id ranges.
      // Today this improves locality in the gather/scatter loops.  More
      // importantly, an intra-pole MPI team can later assign whole ranges to
      // ranks and assemble the resulting pair matrix with one Allgatherv.
      std::sort(catalog.pairs.begin(),catalog.pairs.end(),
        [](const auto& a,const auto& b) {
          return a[1]<b[1] || (a[1]==b[1] && a[0]<b[0]);
        });
      pair_ids[ci].clear();
      pair_ids[ci].reserve(catalog.pairs.size());
      for (int id=0; id<static_cast<int>(catalog.pairs.size()); ++id)
        pair_ids[ci].emplace(
          pair_key(catalog.pairs[id][0],catalog.pairs[id][1]),id);
      catalog.by_second.assign(pseudo_count,{});
      for (int id=0; id<static_cast<int>(catalog.pairs.size()); ++id)
        catalog.by_second.at(catalog.pairs[id][1]).push_back(id);
    }

    for (int cut=0; cut<4; ++cut) {
      auto& cp=cuts_[cut];
      cp.sigma=spec[cut][6];
      cp.local_orientation=spec[cut][7];
      cp.outer_pairs.clear();
      std::unordered_map<std::uint64_t,int> outer_ids;
      for (const auto& path:paths) {
        if (path.weights[0]==0.0) continue;
        const int w=path.pseudo[spec[cut][3]],q=path.pseudo[spec[cut][0]];
        const auto key=pair_key(w,q);
        if (!outer_ids.contains(key)) {
          const int id=cp.outer_pairs.size();
          outer_ids.emplace(key,id);
          cp.outer_pairs.push_back({w,q});
        }
      }
      // Outer self-energy work is naturally batched by its first propagator
      // w.  Contiguous w ranges are independent work blocks for a future
      // (pole,cut,w-block) decomposition and retain the existing batched
      // final convolution.
      std::sort(cp.outer_pairs.begin(),cp.outer_pairs.end(),
        [](const auto& a,const auto& b) {
          return a[0]<b[0] || (a[0]==b[0] && a[1]<b[1]);
        });
      outer_ids.clear();
      outer_ids.reserve(cp.outer_pairs.size());
      for (int id=0; id<static_cast<int>(cp.outer_pairs.size()); ++id)
        outer_ids.emplace(
          pair_key(cp.outer_pairs[id][0],cp.outer_pairs[id][1]),id);
      // Each sparse matrix maps cached inner-pair values to the cut-specific
      // outer-pair rows for one (expanded bath, local bath) combination.
      std::vector<std::vector<Triplet<double>>> triplets(channel_count*channel_count);
      const int catalog_index=cp.sigma>0 ? 0 : 1;
      for (const auto& path:paths) {
        if (path.weights[0]==0.0) continue;
        const int expanded=spec[cut][4]==0 ? path.s0 : path.s1;
        const int local=spec[cut][5]==0 ? path.s0 : path.s1;
        const int u=path.pseudo[spec[cut][1]],v=path.pseudo[spec[cut][2]];
        const int w=path.pseudo[spec[cut][3]],q=path.pseudo[spec[cut][0]];
        const int col=pair_ids[catalog_index].at(pair_key(u,v));
        const int row=outer_ids.at(pair_key(w,q));
        triplets[expanded*channel_count+local].emplace_back(row,col,-path.weights[0]);
      }
      cp.contractions.clear();
      cp.contractions.resize(channel_count*channel_count);
      for (int i=0; i<static_cast<int>(cp.contractions.size()); ++i) {
        cp.contractions[i].resize(cp.outer_pairs.size(),pair_catalogs_[catalog_index].pairs.size());
        cp.contractions[i].setFromTriplets(triplets[i].begin(),triplets[i].end());
        cp.contractions[i].makeCompressed();
      }
      cp.outer_by_first.assign(pseudo_count,{});
      for (int id=0; id<static_cast<int>(cp.outer_pairs.size()); ++id)
        cp.outer_by_first.at(cp.outer_pairs[id][0]).push_back(id);
    }
  }

  bool matches(int pseudo_count,
               int channel_count,
               const std::vector<int>& importance,
               int cut_sum) const
  {
    return pseudo_count_==pseudo_count && channel_count_==channel_count &&
           importance_==importance && cut_sum_==cut_sum &&
           !pair_catalogs_[0].pairs.empty();
  }

  std::size_t nonzero_path_count() const { return nonzero_path_count_; }
  int forward_pair_count() const { return pair_catalogs_[0].pairs.size(); }
  int backward_pair_count() const { return pair_catalogs_[1].pairs.size(); }

  MatrixXd evaluate(const imag_dlr::Workspace& workspace,
                    const std::vector<double>& pseudo_degeneracy,
                    const blitz::Array<double,1>& output_tau,
                    std::array<MatrixXd,4>* cut_results=nullptr) const
  {
    if (!workspace.has_pseudos()) throw std::logic_error("OCA DLR evaluation requires pseudo-particle coefficients");
    if (workspace.pseudo_count()!=pseudo_count_ || workspace.channel_count()!=channel_count_)
      throw std::invalid_argument("OCA DLR plan/workspace dimension mismatch");
    const auto& basis=workspace.basis();
    const int r=workspace.rank();
    std::array<MatrixXd,4> cut_nodes;
    for (auto& nodes:cut_nodes) nodes=MatrixXd::Zero(r,pseudo_count_);

    // For a fixed AAA pole compute each unique inner pair once, then reuse it
    // in all four cuts and all bath-channel contractions.
    for (int l=0; l<workspace.pole_count(); ++l) {
      // Both orientations use the same bounded kernel dressing for a fixed
      // original pole: K(+omega_l) for omega_l>0 and K(-omega_l) otherwise.
      const int tilt_sign=workspace.poles()(l)>0.0 ? +1 : -1;
      const MatrixXd tilted=workspace.tilted_pseudo_coefficients(all_pseudos_,l,tilt_sign);
      std::array<MatrixXd,2> pair_values;
      for (int ci=0; ci<2; ++ci) {
        const int sigma=ci==0 ? +1 : -1;
        const bool nonpositive=sigma*workspace.poles()(l)<=0.0;
        const auto& catalog=pair_catalogs_[ci];
        pair_values[ci].resize(r,catalog.pairs.size());
        for (int second=0; second<pseudo_count_; ++second) {
          const auto& ids=catalog.by_second[second];
          if (ids.empty()) continue;
          MatrixXd rhs(r,ids.size());
          VectorXd lhs;
          if (nonpositive) {
            lhs=workspace.pseudo_coefficients().col(second);
            for (int j=0; j<static_cast<int>(ids.size()); ++j)
              rhs.col(j)=tilted.col(catalog.pairs[ids[j]][0]);
          } else {
            lhs=tilted.col(second);
            for (int j=0; j<static_cast<int>(ids.size()); ++j)
              rhs.col(j)=workspace.pseudo_coefficients().col(catalog.pairs[ids[j]][0]);
          }
          const MatrixXd values=basis.time_ordered_convolution_values(lhs,rhs);
          for (int j=0; j<static_cast<int>(ids.size()); ++j)
            pair_values[ci].col(ids[j])=values.col(j);
        }
      }

      for (int cut=0; cut<4; ++cut) {
        const auto& cp=cuts_[cut];
        const int ci=cp.sigma>0 ? 0 : 1;
        const bool nonpositive=cp.sigma*workspace.poles()(l)<=0.0;
        const MatrixXd& local_delta=cp.local_orientation>0
          ? workspace.delta_plus() : workspace.delta_minus();
        // h(row,t_k) is the sum over every atomic path, inner pair, and local
        // hybridization channel that feeds a particular outer pair.  Sparse
        // matrix times dense matrix performs that entire internal-index sum.
        MatrixXd h_transposed=MatrixXd::Zero(cp.outer_pairs.size(),r);
        const MatrixXd b_transposed=pair_values[ci].transpose();
        for (int expanded=0; expanded<channel_count_; ++expanded) {
          const double residue=workspace.residues()(expanded,l);
          if (residue==0.0) continue;
          for (int local=0; local<channel_count_; ++local) {
            const auto& contraction=cp.contractions[expanded*channel_count_+local];
            if (contraction.nonZeros()==0) continue;
            MatrixXd contribution=contraction*b_transposed;
            contribution.array().rowwise() *= local_delta.col(local).transpose().array();
            h_transposed.noalias() += residue*contribution;
          }
        }
        const MatrixXd h_values=h_transposed.transpose();
        const MatrixXd h_coefficients=basis.values_to_coefficients(h_values);

        // Restore the kernel factor depending on the external time after the
        // stable pole factorization.  The two branches below are algebraic
        // rearrangements chosen so only bounded K(t,+|w|) values appear.
        VectorXd external(r);
        if (nonpositive) {
          if (cp.sigma>0)
            external=workspace.kernel_plus().col(l)/workspace.kernel_zero_minus()(l);
          else
            external=workspace.kernel_minus().col(l)/workspace.kernel_zero_plus()(l);
        } else {
          const double k0=cp.sigma>0
            ? workspace.kernel_zero_plus()(l) : workspace.kernel_zero_minus()(l);
          external.setConstant(1.0/k0);
        }
        // Final time-ordered convolution with G_w.  Group outer pairs by w so
        // all q outputs sharing that propagator form one batched call.
        for (int outer=0; outer<pseudo_count_; ++outer) {
          const auto& ids=cp.outer_by_first[outer];
          if (ids.empty()) continue;
          MatrixXd rhs(r,ids.size());
          for (int j=0; j<static_cast<int>(ids.size()); ++j)
            rhs.col(j)=h_coefficients.col(ids[j]);
          const VectorXd lhs=nonpositive
            ? VectorXd(workspace.pseudo_coefficients().col(outer))
            : VectorXd(tilted.col(outer));
          const MatrixXd values=basis.time_ordered_convolution_values(lhs,rhs);
          for (int j=0; j<static_cast<int>(ids.size()); ++j) {
            const int output=cp.outer_pairs[ids[j]][1];
            cut_nodes[cut].col(output).array() +=
              external.array()*values.col(j).array()/pseudo_degeneracy[output];
          }
        }
      }
    }

    MatrixXd sigma_nodes=MatrixXd::Zero(r,pseudo_count_);
    for (const auto& nodes:cut_nodes) sigma_nodes+=nodes;
    const MatrixXd sigma_coefficients=basis.values_to_coefficients(sigma_nodes);
    MatrixXd result(pseudo_count_,output_tau.extent(0));
    for (int it=0; it<output_tau.extent(0); ++it)
      result.col(it)=basis.evaluate(output_tau(it),sigma_coefficients).transpose();
    if (cut_results){
      for (int cut=0; cut<4; ++cut){
        const MatrixXd coefficients=basis.values_to_coefficients(cut_nodes[cut]);
        (*cut_results)[cut].resize(pseudo_count_,output_tau.extent(0));
        for (int it=0; it<output_tau.extent(0); ++it)
          (*cut_results)[cut].col(it)=basis.evaluate(output_tau(it),coefficients).transpose();
      }
    }
    return result;
  }

private:
  friend class CombinedEvaluator;
  int pseudo_count_=0,channel_count_=0,cut_sum_=0;
  std::size_t nonzero_path_count_=0;
  std::vector<int> importance_,all_pseudos_;
  std::array<PairCatalog,2> pair_catalogs_;
  std::array<CutPlan,4> cuts_;
};

// Results of cutting either hybridization line in the OCA Phi diagram.
//
// The matrices contain the independent bath channels only.  They have not
// yet been divided by Q or by the bath degeneracy, and superconducting bath
// copies have not yet been reconstructed.  Keeping those CIX-specific steps
// outside this class makes the contraction below depend only on diagram
// topology and therefore reusable by higher-order implementations.
struct PhysicalResult {
  MatrixXd g;
  MatrixXd f00;
  MatrixXd f01;
};

// Compact evaluator for the physical Green function and the two Bulla
// numerators generated by the OCA Phi diagram.
//
// Cutting s1 leaves the s0 hybridization in the diagram.  With T denoting the
// time of the cut line, x=t1 and y=t3-T, its scalar integral is
//
//   int_0^T dx G2(x) G3(T-x)
//   int_0^(beta-T) dy G0(y) G1(beta-T-y) Delta_s0(beta-T-y+x).
//
// Cutting s0 similarly gives
//
//   int_0^T dx G3(x) G0(T-x)
//   int_0^(beta-T) dy G1(y) G2(beta-T-y) Delta_s1(T+y-x).
//
// After expanding the remaining Delta in real AAA poles, each double integral
// factorizes into a product of two one-dimensional DLR convolutions.  Thus no
// numerical integration over T, x, or y remains.  The many atomic matrix
// elements are summed by sparse matrices after equal pseudo-particle pairs
// have been combined.
class PhysicalPlan {
  struct PairCatalog {
    // An ordered pair (a,b) represents G_a(x) G_b(t-x).  Ordering matters
    // because the AAA kernel may dress either x or t-x.
    std::vector<std::array<int,2>> pairs;

    // Pair ids grouped by b.  A whole group shares the same left operand in
    // the batched convolution, reducing many small Eigen operations to one.
    std::vector<std::vector<int>> by_second;
  };

  struct BranchPlan {
    // Catalog holding the pair evaluated at T and the pair evaluated at
    // beta-T.  Catalog 0 is the forward family; catalog 1 is the backward
    // family, with the same convention as SelfEnergyPlan.
    int left_catalog=0;
    int right_catalog=0;

    // One sparse contraction for every (expanded bath, cut/output bath)
    // combination and for each observable.  Rows index the left pair and
    // columns index the right pair.
    std::vector<SparseMatrix<double>> g;
    std::vector<SparseMatrix<double>> f00;
    std::vector<SparseMatrix<double>> f01;
  };

public:
  void build(const std::vector<OCAPathRecord>& paths,
             int pseudo_count,
             int channel_count,
             const std::vector<int>& importance,
             int cut_sum)
  {
    pseudo_count_=pseudo_count;
    channel_count_=channel_count;
    importance_=importance;
    cut_sum_=cut_sum;
    nonzero_record_count_=0;
    all_pseudos_.resize(pseudo_count_);
    for (int p=0; p<pseudo_count_; ++p) all_pseudos_[p]=p;

    for (auto& catalog:pair_catalogs_) {
      catalog.pairs.clear();
      catalog.by_second.clear();
    }
    std::array<std::unordered_map<std::uint64_t,int>,2> pair_ids;

    // Insert a pair only once.  The two physical branches use exactly the
    // same forward/backward pair families that occur in the four self-energy
    // cuts:
    //
    //   cut s1: (p2,p3)_backward at T, (p0,p1)_forward at beta-T
    //   cut s0: (p3,p0)_forward  at T, (p1,p2)_backward at beta-T.
    const auto insert_pair=[&](int catalog_index,int first,int second) {
      const auto key=pair_key(first,second);
      auto [it,inserted]=pair_ids[catalog_index].emplace(
        key,static_cast<int>(pair_catalogs_[catalog_index].pairs.size()));
      if (inserted) pair_catalogs_[catalog_index].pairs.push_back({first,second});
      return it->second;
    };

    for (const auto& path:paths) {
      const bool need_s1=path.weights[0]!=0.0 || path.weights[1]!=0.0 ||
                         path.weights[3]!=0.0;
      const bool need_s0=path.weights[0]!=0.0 || path.weights[2]!=0.0 ||
                         path.weights[4]!=0.0;
      if (need_s1) {
        insert_pair(1,path.pseudo[2],path.pseudo[3]);
        insert_pair(0,path.pseudo[0],path.pseudo[1]);
      }
      if (need_s0) {
        insert_pair(0,path.pseudo[3],path.pseudo[0]);
        insert_pair(1,path.pseudo[1],path.pseudo[2]);
      }
      if (need_s1 || need_s0) ++nonzero_record_count_;
    }

    for (int ci=0; ci<2; ++ci) {
      auto& catalog=pair_catalogs_[ci];
      // Use the same canonical order as SelfEnergyPlan: all pairs sharing a
      // second pseudo index form one contiguous range.  Identical ordering is
      // required by the shared pole cache and is also the unit of future
      // distributed pair construction.
      std::sort(catalog.pairs.begin(),catalog.pairs.end(),
        [](const auto& a,const auto& b) {
          return a[1]<b[1] || (a[1]==b[1] && a[0]<b[0]);
        });
      pair_ids[ci].clear();
      pair_ids[ci].reserve(catalog.pairs.size());
      for (int id=0; id<static_cast<int>(catalog.pairs.size()); ++id)
        pair_ids[ci].emplace(
          pair_key(catalog.pairs[id][0],catalog.pairs[id][1]),id);
      catalog.by_second.assign(pseudo_count_,{});
      for (int id=0; id<static_cast<int>(catalog.pairs.size()); ++id)
        catalog.by_second.at(catalog.pairs[id][1]).push_back(id);
    }

    // Branch zero means that s1 was cut; branch one means that s0 was cut.
    // This naming avoids confusing the branch number with the output bath.
    branches_[0].left_catalog=1;
    branches_[0].right_catalog=0;
    branches_[1].left_catalog=0;
    branches_[1].right_catalog=1;

    for (int branch=0; branch<2; ++branch) {
      auto& bp=branches_[branch];
      const int matrix_count=channel_count_*channel_count_;
      std::vector<std::vector<Triplet<double>>> tg(matrix_count),
                                                   tf00(matrix_count),
                                                   tf01(matrix_count);
      for (const auto& path:paths) {
        const bool cut_s1=branch==0;
        const int expanded=cut_s1 ? path.s0 : path.s1;
        const int output=cut_s1 ? path.s1 : path.s0;
        if (expanded<0 || expanded>=channel_count_ ||
            output<0 || output>=channel_count_)
          throw std::out_of_range("OCA physical bath channel in DLR plan");

        const double wg=path.weights[0];
        const double wf00=cut_s1 ? path.weights[1] : path.weights[2];
        const double wf01=cut_s1 ? path.weights[3] : path.weights[4];
        if (wg==0.0 && wf00==0.0 && wf01==0.0) continue;

        const int lf=cut_s1 ? path.pseudo[2] : path.pseudo[3];
        const int ls=cut_s1 ? path.pseudo[3] : path.pseudo[0];
        const int rf=cut_s1 ? path.pseudo[0] : path.pseudo[1];
        const int rs=cut_s1 ? path.pseudo[1] : path.pseudo[2];
        const int row=pair_ids[bp.left_catalog].at(pair_key(lf,ls));
        const int col=pair_ids[bp.right_catalog].at(pair_key(rf,rs));
        const int matrix_id=expanded*channel_count_+output;

        // The minus sign is the crossing sign used in the legacy evaluator.
        // Different weights correspond to replacing the cut F operator by
        // the ordinary operator, Nambu-diagonal Bulla operator, or anomalous
        // Bulla operator, respectively.
        if (wg!=0.0) tg[matrix_id].emplace_back(row,col,-wg);
        if (wf00!=0.0) tf00[matrix_id].emplace_back(row,col,-wf00);
        if (wf01!=0.0) tf01[matrix_id].emplace_back(row,col,-wf01);
      }

      const int rows=pair_catalogs_[bp.left_catalog].pairs.size();
      const int cols=pair_catalogs_[bp.right_catalog].pairs.size();
      const auto finish=[&](std::vector<SparseMatrix<double>>& matrices,
                            const std::vector<std::vector<Triplet<double>>>& trips) {
        matrices.resize(trips.size());
        for (int i=0; i<static_cast<int>(trips.size()); ++i) {
          matrices[i].resize(rows,cols);
          // setFromTriplets deliberately sums duplicate atomic paths that
          // reduce to the same two pair ids and bath channels.
          matrices[i].setFromTriplets(trips[i].begin(),trips[i].end());
          matrices[i].makeCompressed();
        }
      };
      finish(bp.g,tg);
      finish(bp.f00,tf00);
      finish(bp.f01,tf01);
    }
  }

  bool matches(int pseudo_count,
               int channel_count,
               const std::vector<int>& importance,
               int cut_sum) const
  {
    return pseudo_count_==pseudo_count && channel_count_==channel_count &&
           importance_==importance && cut_sum_==cut_sum &&
           !pair_catalogs_[0].pairs.empty() && !pair_catalogs_[1].pairs.empty();
  }

  std::size_t nonzero_record_count() const { return nonzero_record_count_; }
  int forward_pair_count() const { return pair_catalogs_[0].pairs.size(); }
  int backward_pair_count() const { return pair_catalogs_[1].pairs.size(); }

  PhysicalResult evaluate(const imag_dlr::Workspace& workspace,
                          const blitz::Array<double,1>& output_tau) const
  {
    if (!workspace.has_pseudos())
      throw std::logic_error("OCA DLR physical evaluation requires pseudo-particle coefficients");
    if (workspace.pseudo_count()!=pseudo_count_ ||
        workspace.channel_count()!=channel_count_)
      throw std::invalid_argument("OCA DLR physical plan/workspace dimension mismatch");

    const auto& basis=workspace.basis();
    const int r=workspace.rank();
    MatrixXd g_nodes=MatrixXd::Zero(r,channel_count_);
    MatrixXd f00_nodes=MatrixXd::Zero(r,channel_count_);
    MatrixXd f01_nodes=MatrixXd::Zero(r,channel_count_);

    // Convert values at the DLR nodes to values at beta-tau at those same
    // nodes.  Physical diagrams naturally contain both intervals T and
    // beta-T; doing this on only r representative points avoids returning to
    // the original 512-point integration grid.
    const auto reflected_values=[&](const MatrixXd& node_values) {
      const MatrixXd coefficients=basis.values_to_coefficients(node_values);
      MatrixXd reflected(r,node_values.cols());
      for (int k=0; k<r; ++k)
        reflected.row(k)=basis.evaluate(
          basis.beta()-basis.time_nodes()(k),coefficients);
      return reflected;
    };

    // Evaluate all ordered pairs in one catalog with the AAA kernel attached
    // either to the first factor G_a(x), or to the second factor G_b(t-x).
    const auto dressed_pair_values=[&](const PairCatalog& catalog,
                                       const MatrixXd& tilted,
                                       bool dress_first) {
      MatrixXd values(r,catalog.pairs.size());
      for (int second=0; second<pseudo_count_; ++second) {
        const auto& ids=catalog.by_second[second];
        if (ids.empty()) continue;
        MatrixXd rhs(r,ids.size());
        VectorXd lhs;
        if (dress_first) {
          // Integral G_first(x) K(x,w) G_second(t-x): x is the
          // integration variable (the right operand of our convolution).
          lhs=workspace.pseudo_coefficients().col(second);
          for (int j=0; j<static_cast<int>(ids.size()); ++j)
            rhs.col(j)=tilted.col(catalog.pairs[ids[j]][0]);
        } else {
          // Integral G_first(x) G_second(t-x) K(t-x,w): the common
          // second pseudo-particle is the dressed left operand.
          lhs=tilted.col(second);
          for (int j=0; j<static_cast<int>(ids.size()); ++j)
            rhs.col(j)=workspace.pseudo_coefficients().col(catalog.pairs[ids[j]][0]);
        }
        const MatrixXd batch=basis.time_ordered_convolution_values(lhs,rhs);
        for (int j=0; j<static_cast<int>(ids.size()); ++j)
          values.col(ids[j])=batch.col(j);
      }
      return values;
    };

    // Apply all atomic matrix-element sums for one physical branch.  C is a
    // sparse map from right-pair ids to left-pair ids.  For every DLR time
    // node, sum_left L(left,T) [C R(beta-T)](left).
    const auto contract=[&](const std::vector<SparseMatrix<double>>& matrices,
                            const MatrixXd& left,
                            const MatrixXd& right_reflected,
                            int pole,
                            MatrixXd& output) {
      const MatrixXd right_transposed=right_reflected.transpose();
      const MatrixXd left_transposed=left.transpose();
      for (int expanded=0; expanded<channel_count_; ++expanded) {
        // FourierOfDelta in the legacy imaginary-axis solver stores
        // Delta_code(tau) = -sum_l R_l K(tau,omega_l), whereas Workspace uses
        // the standard Lehmann kernel sum_l R_l K.  A self-energy diagram has
        // two Delta lines and the two convention signs cancel.  A physical
        // cut leaves exactly one Delta line, hence the explicit minus here.
        const double prefactor=-workspace.residues()(expanded,pole)/
          (workspace.poles()(pole)>=0.0
             ? workspace.kernel_zero_plus()(pole)
             : workspace.kernel_zero_minus()(pole));
        if (prefactor==0.0) continue;
        for (int out=0; out<channel_count_; ++out) {
          const auto& c=matrices[expanded*channel_count_+out];
          if (c.nonZeros()==0) continue;
          const MatrixXd mapped=c*right_transposed;
          output.col(out).array() += prefactor*
            (left_transposed.array()*mapped.array()).colwise().sum().transpose();
        }
      }
    };

    for (int l=0; l<workspace.pole_count(); ++l) {
      // The stable factorization always uses a positive-frequency kernel:
      // K(+w) for w>=0 and K(-w) for w<0.  Its value at zero is O(1), so the
      // division in contract() never introduces exp(beta*|w|) overflow.
      const bool positive=workspace.poles()(l)>=0.0;
      const int tilt_sign=positive ? +1 : -1;
      const MatrixXd tilted=workspace.tilted_pseudo_coefficients(
        all_pseudos_,l,tilt_sign);

      std::array<MatrixXd,2> first,second,first_reflected,second_reflected;
      for (int ci=0; ci<2; ++ci) {
        first[ci]=dressed_pair_values(pair_catalogs_[ci],tilted,true);
        second[ci]=dressed_pair_values(pair_catalogs_[ci],tilted,false);
        first_reflected[ci]=reflected_values(first[ci]);
        second_reflected[ci]=reflected_values(second[ci]);
      }

      // Exact stable kernel identities used here:
      //
      // cut s1, w>=0: K(B-y+x,w)=K(x,w)K(B-y,w)/K(0,w)
      // cut s1, w< 0: K(B-y+x,w)=K(T-x,-w)K(y,-w)/K(0,-w)
      // cut s0, w>=0: K(T+y-x,w)=K(T-x,w)K(y,w)/K(0,w)
      // cut s0, w< 0: K(T+y-x,w)=K(x,-w)K(B-y,-w)/K(0,-w)
      //
      // 'first' means the kernel dresses x; 'second' means it dresses
      // t-x.  Reflected matrices are evaluated at B=beta-T.
      const MatrixXd& s1_left=positive ? first[1] : second[1];
      const MatrixXd& s1_right=positive ? second_reflected[0]
                                        : first_reflected[0];
      const MatrixXd& s0_left=positive ? second[0] : first[0];
      const MatrixXd& s0_right=positive ? first_reflected[1]
                                        : second_reflected[1];

      contract(branches_[0].g,s1_left,s1_right,l,g_nodes);
      contract(branches_[0].f00,s1_left,s1_right,l,f00_nodes);
      contract(branches_[0].f01,s1_left,s1_right,l,f01_nodes);
      contract(branches_[1].g,s0_left,s0_right,l,g_nodes);
      contract(branches_[1].f00,s0_left,s0_right,l,f00_nodes);
      contract(branches_[1].f01,s0_left,s0_right,l,f01_nodes);
    }

    const auto to_output_mesh=[&](const MatrixXd& nodes) {
      const MatrixXd coefficients=basis.values_to_coefficients(nodes);
      MatrixXd result(channel_count_,output_tau.extent(0));
      for (int it=0; it<output_tau.extent(0); ++it)
        result.col(it)=basis.evaluate(output_tau(it),coefficients).transpose();
      return result;
    };
    return {to_output_mesh(g_nodes),to_output_mesh(f00_nodes),
            to_output_mesh(f01_nodes)};
  }

private:
  friend class CombinedEvaluator;
  int pseudo_count_=0,channel_count_=0,cut_sum_=0;
  std::size_t nonzero_record_count_=0;
  std::vector<int> importance_,all_pseudos_;
  std::array<PairCatalog,2> pair_catalogs_;
  std::array<BranchPlan,2> branches_;
};

// Result of the production OCA evaluator when physical cuts are requested.
// Keeping both objects together makes it explicit that they were accumulated
// during the same AAA-pole pass and therefore share all expensive pair data.
struct CombinedResult {
  MatrixXd sigma;
  PhysicalResult physical;
};

// Description of the communicator layout selected by evaluate_mpi().  It is
// returned for diagnostics only; numerical ownership is internal to the
// evaluator and every world rank receives the same final result.
struct ParallelLayout {
  int world_size=1;
  int topology_shard_count=1;
  int topology_shard_id=0;
  int topology_group_rank=0;
  int topology_group_size=1;
  int min_topology_group_size=1;
  int max_topology_group_size=1;
  int team_count=1;
  int team_id=0;
  int team_rank=0;
  int team_size=1;
  int min_team_size=1;
  int max_team_size=1;
  int owned_poles=0;
};

// Additive result on the representative DLR time nodes.
//
// This is deliberately the boundary between expensive diagram work and the
// inexpensive final interpolation to the solver's dense tau mesh.  Every AAA
// pole contributes linearly to these matrices, so independent pole ranges can
// be accumulated on different MPI teams and combined by addition.  A future
// MPI implementation therefore needs to reduce only O(r*(Np+3*Nb)) numbers,
// rather than dense pair tables or dense-tau functions.
//
// The four individual self-energy cuts are diagnostic data.  Production runs
// normally keep only their sum in sigma_nodes, which avoids four additional
// r-by-Np matrices in every rank.  Likewise physical matrices are absent on
// inexpensive iterations that do not request G/F00/F01.
struct NodeAccumulator {
  MatrixXd sigma_nodes;
  MatrixXd g_nodes;
  MatrixXd f00_nodes;
  MatrixXd f01_nodes;
  std::array<MatrixXd,4> cut_nodes;
  bool physical_enabled=false;
  bool cuts_enabled=false;

  NodeAccumulator()=default;

  NodeAccumulator(int rank,
                  int pseudo_count,
                  int channel_count,
                  bool enable_physical,
                  bool keep_individual_cuts)
    : sigma_nodes(MatrixXd::Zero(rank,pseudo_count)),
      physical_enabled(enable_physical),
      cuts_enabled(keep_individual_cuts)
  {
    if (enable_physical) {
      g_nodes=MatrixXd::Zero(rank,channel_count);
      f00_nodes=MatrixXd::Zero(rank,channel_count);
      f01_nodes=MatrixXd::Zero(rank,channel_count);
    }
    if (keep_individual_cuts)
      for (auto& nodes:cut_nodes)
        nodes=MatrixXd::Zero(rank,pseudo_count);
  }

  // Combine independently evaluated work ranges.  This is also the exact
  // algebra that a later MPI_Allreduce performs element by element.
  void add(const NodeAccumulator& other)
  {
    const auto same_shape=[](const MatrixXd& a,const MatrixXd& b) {
      return a.rows()==b.rows() && a.cols()==b.cols();
    };
    if (physical_enabled!=other.physical_enabled ||
        cuts_enabled!=other.cuts_enabled ||
        !same_shape(sigma_nodes,other.sigma_nodes))
      throw std::invalid_argument("incompatible OCA DLR node accumulators");
    sigma_nodes+=other.sigma_nodes;
    if (physical_enabled) {
      if (!same_shape(g_nodes,other.g_nodes) ||
          !same_shape(f00_nodes,other.f00_nodes) ||
          !same_shape(f01_nodes,other.f01_nodes))
        throw std::invalid_argument(
          "incompatible OCA DLR physical node accumulators");
      g_nodes+=other.g_nodes;
      f00_nodes+=other.f00_nodes;
      f01_nodes+=other.f01_nodes;
    }
    if (cuts_enabled)
      for (int cut=0; cut<4; ++cut) {
        if (!same_shape(cut_nodes[cut],other.cut_nodes[cut]))
          throw std::invalid_argument(
            "incompatible OCA DLR cut node accumulators");
        cut_nodes[cut]+=other.cut_nodes[cut];
      }
  }
};

// Joint pseudo-self-energy and physical-cut evaluator.
//
// SelfEnergyPlan::evaluate and PhysicalPlan::evaluate remain useful as small,
// independent reference paths and for focused tests.  Production iterations
// that request G/F00/F01 use this class instead.  The important reordering is
//
//   old: all poles -> self energy, then all poles -> physical observables
//   new: pole l -> shared pair values -> every self and physical cut.
//
// For a fixed pole only two pair families are mathematically needed.  If the
// pole is positive they are {catalog 0, kernel on second factor} and
// {catalog 1, kernel on first factor}; a negative pole reverses those roles.
// Those are exactly the two families already needed by the four pseudo cuts.
// Physical cuts use the same values, one directly at T and one reflected to
// beta-T.  Thus this loop replaces six pair-family evaluations by two.
class CombinedEvaluator {
  struct TeamContext {
    int rank=0;
    int size=1;
#ifdef _MPI
    MPI_Comm communicator=MPI_COMM_NULL;
#endif
  };

  struct PhysicalTask {
    int branch=0;
    int observable=0; // 0=G, 1=F00, 2=F01
    int matrix_id=0;  // expanded*Nb+output
    double cost=0.0;
  };

public:
  // The two plans are built independently so they can still be used alone.
  // Joint evaluation requires identical pair ids, not merely identical pair
  // sets, because every sparse row/column index refers to those ids.
  static bool compatible(const SelfEnergyPlan& self,
                         const PhysicalPlan& physical)
  {
    if (self.pseudo_count_!=physical.pseudo_count_ ||
        self.channel_count_!=physical.channel_count_)
      return false;
    for (int ci=0; ci<2; ++ci)
      if (self.pair_catalogs_[ci].pairs!=physical.pair_catalogs_[ci].pairs)
        return false;
    return true;
  }

  // Allocate only the node matrices needed by this iteration.  Keeping this
  // operation explicit makes ownership clear when each MPI team later owns a
  // private partial accumulator.
  static NodeAccumulator make_accumulator(
      const SelfEnergyPlan& self,
      const imag_dlr::Workspace& workspace,
      bool enable_physical,
      bool keep_individual_cuts=false)
  {
    if (!workspace.has_pseudos())
      throw std::logic_error(
        "combined OCA DLR evaluation requires pseudo-particle coefficients");
    if (workspace.pseudo_count()!=self.pseudo_count_ ||
        workspace.channel_count()!=self.channel_count_)
      throw std::invalid_argument(
        "combined OCA plan/workspace dimension mismatch");
    return NodeAccumulator(workspace.rank(),self.pseudo_count_,
                           self.channel_count_,enable_physical,
                           keep_individual_cuts);
  }

  // Add the half-open AAA-pole range [first_pole,last_pole) to accumulator.
  // No final interpolation occurs here.  Calls on disjoint pole ranges are
  // independent and their NodeAccumulators may be added in any order.
  //
  // The next parallel layer can subdivide one pole further.  Pair catalogs
  // are already stored in contiguous equal-second ranges, outer self-energy
  // pairs in contiguous equal-w ranges, and reflected pair rows are local.
  // Thus suitable intra-team jobs are (pole,catalog,second),
  // (pole,cut,w-range), and (pole,physical-branch,left-row-range).
  static void accumulate_poles(
      const SelfEnergyPlan& self,
      const PhysicalPlan* physical,
      const imag_dlr::Workspace& workspace,
      const std::vector<double>& pseudo_degeneracy,
      int first_pole,
      int last_pole,
      NodeAccumulator& accumulator)
  {
    if (first_pole<0 || last_pole<first_pole ||
        last_pole>workspace.pole_count())
      throw std::out_of_range("invalid OCA DLR AAA-pole range");
    validate(self,physical,workspace,pseudo_degeneracy,accumulator);
    TeamContext serial;
    accumulate_sequence(self,physical,workspace,pseudo_degeneracy,
                        first_pole,last_pole,1,serial,accumulator);
  }

  // Perform the one values-to-coefficients solve and dense-mesh evaluation
  // after all local ranges (and, later, the global MPI reduction) are done.
  static CombinedResult finalize(
      const SelfEnergyPlan& self,
      const imag_dlr::Workspace& workspace,
      const NodeAccumulator& accumulator,
      const blitz::Array<double,1>& output_tau,
      std::array<MatrixXd,4>* cut_results=nullptr)
  {
    if (accumulator.sigma_nodes.rows()!=workspace.rank() ||
        accumulator.sigma_nodes.cols()!=self.pseudo_count_)
      throw std::invalid_argument("OCA DLR accumulator dimension mismatch");
    if (cut_results && !accumulator.cuts_enabled)
      throw std::invalid_argument(
        "individual OCA cuts requested from an accumulator that did not keep them");

    const auto& basis=workspace.basis();
    const MatrixXd sigma_coefficients=
      basis.values_to_coefficients(accumulator.sigma_nodes);

    CombinedResult result;
    result.sigma.resize(self.pseudo_count_,output_tau.extent(0));
    for (int it=0; it<output_tau.extent(0); ++it)
      result.sigma.col(it)=
        basis.evaluate(output_tau(it),sigma_coefficients).transpose();

    if (cut_results) {
      for (int cut=0; cut<4; ++cut) {
        const MatrixXd coefficients=
          basis.values_to_coefficients(accumulator.cut_nodes[cut]);
        (*cut_results)[cut].resize(self.pseudo_count_,output_tau.extent(0));
        for (int it=0; it<output_tau.extent(0); ++it)
          (*cut_results)[cut].col(it)=
            basis.evaluate(output_tau(it),coefficients).transpose();
      }
    }

    if (accumulator.physical_enabled) {
      const auto to_output_mesh=[&](const MatrixXd& nodes) {
        const MatrixXd coefficients=basis.values_to_coefficients(nodes);
        MatrixXd values(self.channel_count_,output_tau.extent(0));
        for (int it=0; it<output_tau.extent(0); ++it)
          values.col(it)=
            basis.evaluate(output_tau(it),coefficients).transpose();
        return values;
      };
      result.physical={to_output_mesh(accumulator.g_nodes),
                       to_output_mesh(accumulator.f00_nodes),
                       to_output_mesh(accumulator.f01_nodes)};
    }
    return result;
  }

  // Serial convenience wrapper.  Its implementation intentionally uses the
  // same range API and additive accumulator that distributed execution will
  // use; there is no separate numerical algorithm to keep synchronized.
  static CombinedResult evaluate(
      const SelfEnergyPlan& self,
      const PhysicalPlan* physical,
      const imag_dlr::Workspace& workspace,
      const std::vector<double>& pseudo_degeneracy,
      const blitz::Array<double,1>& output_tau,
      std::array<MatrixXd,4>* cut_results=nullptr)
  {
    if (physical && !compatible(self,*physical))
      throw std::invalid_argument(
        "combined OCA evaluator requires identical self/physical pair catalogs");
    NodeAccumulator accumulator=make_accumulator(
      self,workspace,physical!=nullptr,cut_results!=nullptr);
    accumulate_poles(self,physical,workspace,pseudo_degeneracy,
                     0,workspace.pole_count(),accumulator);
    return finalize(self,workspace,accumulator,output_tau,cut_results);
  }

  static CombinedResult evaluate(
      const SelfEnergyPlan& self,
      const PhysicalPlan& physical,
      const imag_dlr::Workspace& workspace,
      const std::vector<double>& pseudo_degeneracy,
      const blitz::Array<double,1>& output_tau,
      std::array<MatrixXd,4>* cut_results=nullptr)
  {
    return evaluate(self,&physical,workspace,pseudo_degeneracy,output_tau,
                    cut_results);
  }

#ifdef _MPI
  // Two-dimensional memory-distributed entry point.  There are
  // `topology_shard_count` disjoint sparse plans globally.  Ranks with equal
  // world_rank%topology_shard_count hold the same plan and form its worker
  // group.  Different pole teams inside that group own disjoint AAA poles;
  // when workers outnumber poles, multiple ranks cooperate on a pole using the
  // existing pair-column and contraction-block decomposition.  Thus topology
  // memory remains bounded while ranks beyond the memory-required shard count
  // accelerate numerical integration instead of creating duplicate catalogs.
  static CombinedResult evaluate_sharded_mpi(
      const SelfEnergyPlan& local_self,
      const PhysicalPlan* local_physical,
      const imag_dlr::Workspace& workspace,
      const std::vector<double>& pseudo_degeneracy,
      const blitz::Array<double,1>& output_tau,
      MPI_Comm world,
      int topology_shard_count,
      int requested_pole_team_count=0,
      ParallelLayout* layout=nullptr,
      std::array<MatrixXd,4>* cut_results=nullptr)
  {
    int world_rank=0,world_size=1;
    MPI_Comm_rank(world,&world_rank);
    MPI_Comm_size(world,&world_size);
    if (local_physical && !compatible(local_self,*local_physical))
      throw std::invalid_argument(
        "sharded combined OCA evaluator requires locally compatible plans");
    const int shard_count=std::max(1,std::min(topology_shard_count,world_size));
    const int shard_id=world_rank%shard_count;
    MPI_Comm shard_comm=MPI_COMM_NULL;
    MPI_Comm_split(world,shard_id,world_rank,&shard_comm);
    int shard_rank=0,shard_size=1;
    MPI_Comm_rank(shard_comm,&shard_rank);
    MPI_Comm_size(shard_comm,&shard_size);

    const int pole_count=workspace.pole_count();
    const int automatic_pole_teams=std::max(1,std::min(pole_count,shard_size));
    const int pole_team_count=requested_pole_team_count<=0
      ? automatic_pole_teams
      : std::max(1,std::min({requested_pole_team_count,pole_count,shard_size}));
    const int pole_team_id=shard_rank%pole_team_count;
    MPI_Comm pole_team_comm=MPI_COMM_NULL;
    MPI_Comm_split(shard_comm,pole_team_id,shard_rank,&pole_team_comm);
    TeamContext pole_team;
    pole_team.communicator=pole_team_comm;
    MPI_Comm_rank(pole_team_comm,&pole_team.rank);
    MPI_Comm_size(pole_team_comm,&pole_team.size);

    if (layout) {
      *layout=ParallelLayout{};
      layout->world_size=world_size;
      layout->topology_shard_count=shard_count;
      layout->topology_shard_id=shard_id;
      layout->topology_group_rank=shard_rank;
      layout->topology_group_size=shard_size;
      layout->min_topology_group_size=world_size/shard_count;
      layout->max_topology_group_size=(world_size+shard_count-1)/shard_count;
      layout->team_count=pole_team_count;
      layout->team_id=pole_team_id;
      layout->team_rank=pole_team.rank;
      layout->team_size=pole_team.size;
      layout->min_team_size=shard_size/pole_team_count;
      layout->max_team_size=(shard_size+pole_team_count-1)/pole_team_count;
      layout->owned_poles=pole_team_id<pole_count
        ? 1+(pole_count-1-pole_team_id)/pole_team_count : 0;
    }
    NodeAccumulator accumulator=make_accumulator(
      local_self,workspace,local_physical!=nullptr,cut_results!=nullptr);
    try {
      validate(local_self,local_physical,workspace,pseudo_degeneracy,accumulator);
      accumulate_sequence(local_self,local_physical,workspace,
                          pseudo_degeneracy,pole_team_id,pole_count,
                          pole_team_count,pole_team,accumulator);
    } catch (...) {
      MPI_Comm_free(&pole_team_comm);
      MPI_Comm_free(&shard_comm);
      throw;
    }
    MPI_Comm_free(&pole_team_comm);
    MPI_Comm_free(&shard_comm);
    allreduce_accumulator(accumulator,world);
    return finalize(local_self,workspace,accumulator,output_tau,cut_results);
  }

  // Distributed production entry point.
  //
  // World ranks are split cyclically into `team_count` communicators.  Team
  // t owns poles t, t+team_count, ... .  When there are no more ranks than
  // poles every team has one rank and this is pure pole parallelism.  With
  // hundreds of ranks, multiple ranks cooperate on one pole by constructing
  // disjoint pair columns and independent self/physical contraction blocks.
  // All ranks finally sum only the compact node accumulator and therefore
  // leave this routine with identical dense-tau results.
  static CombinedResult evaluate_mpi(
      const SelfEnergyPlan& self,
      const PhysicalPlan* physical,
      const imag_dlr::Workspace& workspace,
      const std::vector<double>& pseudo_degeneracy,
      const blitz::Array<double,1>& output_tau,
      MPI_Comm world,
      int requested_team_count=0,
      ParallelLayout* layout=nullptr,
      std::array<MatrixXd,4>* cut_results=nullptr)
  {
    int world_rank=0,world_size=1;
    MPI_Comm_rank(world,&world_rank);
    MPI_Comm_size(world,&world_size);
    if (world_size==1) {
      if (layout) {
        *layout=ParallelLayout{};
        layout->world_size=1;
        layout->topology_group_size=1;
        layout->owned_poles=workspace.pole_count();
      }
      return evaluate(self,physical,workspace,pseudo_degeneracy,output_tau,
                      cut_results);
    }

    const int pole_count=workspace.pole_count();
    const int automatic_teams=std::max(1,std::min(pole_count,world_size));
    const int team_count=requested_team_count<=0
      ? automatic_teams
      : std::max(1,std::min({requested_team_count,pole_count,world_size}));
    const int team_id=world_rank%team_count;
    MPI_Comm team_comm=MPI_COMM_NULL;
    MPI_Comm_split(world,team_id,world_rank,&team_comm);

    TeamContext team;
    team.communicator=team_comm;
    MPI_Comm_rank(team_comm,&team.rank);
    MPI_Comm_size(team_comm,&team.size);

    if (layout) {
      layout->world_size=world_size;
      layout->topology_shard_count=1;
      layout->topology_shard_id=0;
      layout->topology_group_rank=world_rank;
      layout->topology_group_size=world_size;
      layout->min_topology_group_size=world_size;
      layout->max_topology_group_size=world_size;
      layout->team_count=team_count;
      layout->team_id=team_id;
      layout->team_rank=team.rank;
      layout->team_size=team.size;
      layout->min_team_size=world_size/team_count;
      layout->max_team_size=(world_size+team_count-1)/team_count;
      layout->owned_poles=team_id<pole_count
        ? 1+(pole_count-1-team_id)/team_count : 0;
    }

    NodeAccumulator accumulator=make_accumulator(
      self,workspace,physical!=nullptr,cut_results!=nullptr);
    try {
      validate(self,physical,workspace,pseudo_degeneracy,accumulator);
      accumulate_sequence(self,physical,workspace,pseudo_degeneracy,
                          team_id,pole_count,team_count,team,accumulator);
    } catch (...) {
      MPI_Comm_free(&team_comm);
      throw;
    }
    MPI_Comm_free(&team_comm);

    // No pair matrix crosses team boundaries.  The only world-wide
    // communication is this small sum over final representative-node data.
    allreduce_accumulator(accumulator,world);
    return finalize(self,workspace,accumulator,output_tau,cut_results);
  }
#endif

private:
  static void validate(
      const SelfEnergyPlan& self,
      const PhysicalPlan* physical,
      const imag_dlr::Workspace& workspace,
      const std::vector<double>& pseudo_degeneracy,
      const NodeAccumulator& accumulator)
  {
    if (!workspace.has_pseudos())
      throw std::logic_error(
        "combined OCA DLR evaluation requires pseudo-particle coefficients");
    if (workspace.pseudo_count()!=self.pseudo_count_ ||
        workspace.channel_count()!=self.channel_count_)
      throw std::invalid_argument(
        "combined OCA plan/workspace dimension mismatch");
    if (static_cast<int>(pseudo_degeneracy.size())!=self.pseudo_count_)
      throw std::invalid_argument("combined OCA pseudo degeneracy size mismatch");
    if ((physical!=nullptr)!=accumulator.physical_enabled)
      throw std::invalid_argument(
        "OCA DLR accumulator physical mode does not match evaluation");
    if (physical && !compatible(self,*physical))
      throw std::invalid_argument(
        "combined OCA evaluator requires identical self/physical pair catalogs");
    if (accumulator.sigma_nodes.rows()!=workspace.rank() ||
        accumulator.sigma_nodes.cols()!=self.pseudo_count_)
      throw std::invalid_argument("OCA DLR accumulator dimension mismatch");
    if (accumulator.physical_enabled &&
        (accumulator.g_nodes.rows()!=workspace.rank() ||
         accumulator.g_nodes.cols()!=self.channel_count_ ||
         accumulator.f00_nodes.rows()!=workspace.rank() ||
         accumulator.f00_nodes.cols()!=self.channel_count_ ||
         accumulator.f01_nodes.rows()!=workspace.rank() ||
         accumulator.f01_nodes.cols()!=self.channel_count_))
      throw std::invalid_argument(
        "OCA DLR physical accumulator dimension mismatch");
    if (accumulator.cuts_enabled)
      for (const auto& nodes:accumulator.cut_nodes)
        if (nodes.rows()!=workspace.rank() ||
            nodes.cols()!=self.pseudo_count_)
          throw std::invalid_argument("OCA DLR cut accumulator dimension mismatch");
  }

  // Deterministic longest-job-first assignment.  Every rank constructs this
  // tiny schedule independently, avoiding another broadcast.  Costs need not
  // be exact timings; pair counts and sparse nonzeros are good stable proxies.
  static std::vector<int> greedy_owners(const std::vector<double>& costs,
                                        int worker_count)
  {
    std::vector<int> order(costs.size()),owners(costs.size(),-1);
    std::iota(order.begin(),order.end(),0);
    std::stable_sort(order.begin(),order.end(),[&](int a,int b) {
      return costs[a]>costs[b];
    });
    std::vector<double> load(worker_count,0.0);
    for (int task:order) {
      if (!(costs[task]>0.0)) continue;
      const int worker=static_cast<int>(
        std::min_element(load.begin(),load.end())-load.begin());
      owners[task]=worker;
      load[worker]+=costs[task];
    }
    return owners;
  }

  static MatrixXd build_pair_transposed(
      const SelfEnergyPlan& self,
      const imag_dlr::Workspace& workspace,
      const MatrixXd& tilted,
      int pole,
      int catalog_index,
      const TeamContext& team)
  {
    const auto& basis=workspace.basis();
    const auto& catalog=self.pair_catalogs_[catalog_index];
    const int r=workspace.rank();
    const int pair_count=static_cast<int>(catalog.pairs.size());
    const int first=pair_count*team.rank/team.size;
    const int last=pair_count*(team.rank+1)/team.size;
    MatrixXd local_values(r,last-first);
    const int sigma=catalog_index==0 ? +1 : -1;
    const bool nonpositive=sigma*workspace.poles()(pole)<=0.0;

    for (int second=0; second<self.pseudo_count_; ++second) {
      const auto& ids=catalog.by_second[second];
      const auto begin=std::lower_bound(ids.begin(),ids.end(),first);
      const auto end=std::lower_bound(ids.begin(),ids.end(),last);
      const int count=static_cast<int>(end-begin);
      if (count==0) continue;
      MatrixXd rhs(r,count);
      VectorXd lhs;
      if (nonpositive) {
        lhs=workspace.pseudo_coefficients().col(second);
        for (int j=0; j<count; ++j)
          rhs.col(j)=tilted.col(catalog.pairs[*(begin+j)][0]);
      } else {
        lhs=tilted.col(second);
        for (int j=0; j<count; ++j)
          rhs.col(j)=workspace.pseudo_coefficients().col(
            catalog.pairs[*(begin+j)][0]);
      }
      const MatrixXd values=
        basis.time_ordered_convolution_values(lhs,rhs);
      for (int j=0; j<count; ++j)
        local_values.col(*(begin+j)-first)=values.col(j);
    }

    if (team.size==1) return local_values.transpose();
#ifdef _MPI
    // pair_values is [time,pair], so a contiguous interval of pair columns is
    // contiguous in Eigen's column-major storage and can be gathered without
    // packing.  The transpose happens only after the collective.
    MatrixXd pair_values(r,pair_count);
    std::vector<int> counts(team.size),displacements(team.size);
    for (int worker=0; worker<team.size; ++worker) {
      const int b=pair_count*worker/team.size;
      const int e=pair_count*(worker+1)/team.size;
      counts[worker]=(e-b)*r;
      displacements[worker]=b*r;
    }
    MPI_Allgatherv(local_values.data(),static_cast<int>(local_values.size()),
                   MPI_DOUBLE,pair_values.data(),counts.data(),
                   displacements.data(),MPI_DOUBLE,team.communicator);
    return pair_values.transpose();
#else
    throw std::logic_error("multi-rank OCA requested without MPI support");
#endif
  }

  static MatrixXd reflect_pair_transposed(
      const dlr::FermionicBasis& basis,
      const MatrixXd& pair_transposed,
      const TeamContext& team)
  {
    if (team.size==1)
      return basis.reflect_transposed_values(pair_transposed);
#ifdef _MPI
    const int pair_count=pair_transposed.rows();
    const int r=pair_transposed.cols();
    const int first=pair_count*team.rank/team.size;
    const int last=pair_count*(team.rank+1)/team.size;
    const MatrixXd local_reflected=
      basis.reflect_transposed_values(
        MatrixXd(pair_transposed.middleRows(first,last-first)));
    // Gather the transpose for the same reason as pair construction above:
    // each owned pair interval is then a contiguous set of columns.
    const MatrixXd local_values=local_reflected.transpose();
    MatrixXd reflected_values(r,pair_count);
    std::vector<int> counts(team.size),displacements(team.size);
    for (int worker=0; worker<team.size; ++worker) {
      const int b=pair_count*worker/team.size;
      const int e=pair_count*(worker+1)/team.size;
      counts[worker]=(e-b)*r;
      displacements[worker]=b*r;
    }
    MPI_Allgatherv(local_values.data(),static_cast<int>(local_values.size()),
                   MPI_DOUBLE,reflected_values.data(),counts.data(),
                   displacements.data(),MPI_DOUBLE,team.communicator);
    return reflected_values.transpose();
#else
    throw std::logic_error("multi-rank OCA requested without MPI support");
#endif
  }

  static void accumulate_sequence(
      const SelfEnergyPlan& self,
      const PhysicalPlan* physical,
      const imag_dlr::Workspace& workspace,
      const std::vector<double>& pseudo_degeneracy,
      int first_pole,
      int last_pole,
      int pole_stride,
      const TeamContext& team,
      NodeAccumulator& accumulator)
  {
    const auto& basis=workspace.basis();
    const int r=workspace.rank();

    // Intra-pole self-energy tasks are whole equal-w outer blocks.  A block
    // contains every q paired with one outer propagator, preserving the
    // batched final convolution while making blocks independent across ranks.
    std::vector<double> self_costs(4*self.pseudo_count_,0.0);
    if (team.size>1)
      for (int cut=0; cut<4; ++cut)
        for (int outer=0; outer<self.pseudo_count_; ++outer)
          self_costs[cut*self.pseudo_count_+outer]=
            self.cuts_[cut].outer_by_first[outer].size();
    const std::vector<int> self_owner=team.size>1
      ? greedy_owners(self_costs,team.size) : std::vector<int>{};

    // Physical contractions are independent sparse matrix products.  Their
    // nonzero counts give a much better balance than assigning only the two
    // diagram branches, particularly for large cluster CIX files.
    std::vector<PhysicalTask> physical_tasks;
    std::vector<int> physical_owner;
    if (physical) {
      for (int branch=0; branch<2; ++branch)
        for (int observable=0; observable<3; ++observable) {
          const auto& matrices=observable==0 ? physical->branches_[branch].g
            : observable==1 ? physical->branches_[branch].f00
                            : physical->branches_[branch].f01;
          for (int matrix_id=0;
               matrix_id<static_cast<int>(matrices.size()); ++matrix_id)
            if (matrices[matrix_id].nonZeros()!=0)
              physical_tasks.push_back(
                {branch,observable,matrix_id,
                 static_cast<double>(matrices[matrix_id].nonZeros())});
        }
      std::vector<double> costs;
      costs.reserve(physical_tasks.size());
      for (const auto& task:physical_tasks) costs.push_back(task.cost);
      physical_owner=greedy_owners(costs,team.size);
    }

    for (int l=first_pole; l<last_pole; l+=pole_stride) {
      // Stable kernel dressing is formed once on each cooperating rank.  It
      // is only r-by-Npseudo; distributing it would cost more than recomputing
      // it and would add a collective before both pair catalogs.
      const int tilt_sign=workspace.poles()(l)>0.0 ? +1 : -1;
      const MatrixXd tilted=workspace.tilted_pseudo_coefficients(
        self.all_pseudos_,l,tilt_sign);

      // Each rank computes a contiguous slice of pair columns.  Team-local
      // Allgatherv assembles the two complete [pair,time] tables needed by
      // sparse contractions; no pair data is communicated between teams.
      std::array<MatrixXd,2> pair_transposed;
      for (int ci=0; ci<2; ++ci)
        pair_transposed[ci]=build_pair_transposed(
          self,workspace,tilted,l,ci,team);

      for (int cut=0; cut<4; ++cut) {
        const auto& cp=self.cuts_[cut];
        const int ci=cp.sigma>0 ? 0 : 1;
        const bool nonpositive=cp.sigma*workspace.poles()(l)<=0.0;
        const MatrixXd& local_delta=cp.local_orientation>0
          ? workspace.delta_plus() : workspace.delta_minus();
        VectorXd external(r);
        if (nonpositive) {
          if (cp.sigma>0)
            external=workspace.kernel_plus().col(l)/
                     workspace.kernel_zero_minus()(l);
          else
            external=workspace.kernel_minus().col(l)/
                     workspace.kernel_zero_plus()(l);
        } else {
          const double k0=cp.sigma>0
            ? workspace.kernel_zero_plus()(l)
            : workspace.kernel_zero_minus()(l);
          external.setConstant(1.0/k0);
        }

        if (team.size==1) {
          // Preserve the large batched serial path when pole parallelism is
          // sufficient.  This is faster than many tiny outer-block solves.
          MatrixXd h_transposed=MatrixXd::Zero(cp.outer_pairs.size(),r);
          for (int expanded=0; expanded<self.channel_count_; ++expanded) {
            const double residue=workspace.residues()(expanded,l);
            if (residue==0.0) continue;
            for (int local=0; local<self.channel_count_; ++local) {
              const auto& contraction=
                cp.contractions[expanded*self.channel_count_+local];
              if (contraction.nonZeros()==0) continue;
              MatrixXd contribution=contraction*pair_transposed[ci];
              contribution.array().rowwise() *=
                local_delta.col(local).transpose().array();
              h_transposed.noalias() += residue*contribution;
            }
          }
          const MatrixXd h_values=h_transposed.transpose();
          const MatrixXd h_coefficients=
            basis.values_to_coefficients(h_values);
          for (int outer=0; outer<self.pseudo_count_; ++outer) {
            const auto& ids=cp.outer_by_first[outer];
            if (ids.empty()) continue;
            MatrixXd rhs(r,ids.size());
            for (int j=0; j<static_cast<int>(ids.size()); ++j)
              rhs.col(j)=h_coefficients.col(ids[j]);
            const VectorXd lhs=nonpositive
              ? VectorXd(workspace.pseudo_coefficients().col(outer))
              : VectorXd(tilted.col(outer));
            const MatrixXd values=
              basis.time_ordered_convolution_values(lhs,rhs);
            add_self_values(cp,cut,outer,ids.front(),values,external,
                            pseudo_degeneracy,accumulator);
          }
        } else {
          // Multi-rank teams own disjoint equal-w row blocks.  Because outer
          // pairs were canonically sorted, each block is a contiguous sparse
          // row interval and can be multiplied without copying a global h.
          for (int outer=0; outer<self.pseudo_count_; ++outer) {
            const auto& ids=cp.outer_by_first[outer];
            if (ids.empty() ||
                self_owner[cut*self.pseudo_count_+outer]!=team.rank)
              continue;
            const int row_first=ids.front();
            const int row_count=static_cast<int>(ids.size());
            MatrixXd h_transposed=MatrixXd::Zero(row_count,r);
            for (int expanded=0; expanded<self.channel_count_; ++expanded) {
              const double residue=workspace.residues()(expanded,l);
              if (residue==0.0) continue;
              for (int local=0; local<self.channel_count_; ++local) {
                const auto& contraction=
                  cp.contractions[expanded*self.channel_count_+local];
                if (contraction.nonZeros()==0) continue;
                MatrixXd contribution=
                  contraction.middleRows(row_first,row_count)*
                  pair_transposed[ci];
                contribution.array().rowwise() *=
                  local_delta.col(local).transpose().array();
                h_transposed.noalias() += residue*contribution;
              }
            }
            const MatrixXd h_values=h_transposed.transpose();
            const MatrixXd h_coefficients=
              basis.values_to_coefficients(h_values);
            const VectorXd lhs=nonpositive
              ? VectorXd(workspace.pseudo_coefficients().col(outer))
              : VectorXd(tilted.col(outer));
            const MatrixXd values=
              basis.time_ordered_convolution_values(lhs,h_coefficients);
            add_self_values(cp,cut,outer,row_first,values,external,
                            pseudo_degeneracy,accumulator);
          }
        }
      }

      if (physical) {
        // Reflection itself is split by pair rows, then team-gathered.  This
        // avoids repeating an r-by-r multiplication for every pair on every
        // cooperating rank.
        std::array<MatrixXd,2> reflected;
        reflected[0]=reflect_pair_transposed(
          basis,pair_transposed[0],team);
        reflected[1]=reflect_pair_transposed(
          basis,pair_transposed[1],team);

        for (int task_id=0;
             task_id<static_cast<int>(physical_tasks.size()); ++task_id) {
          if (physical_owner[task_id]!=team.rank) continue;
          const auto& task=physical_tasks[task_id];
          const auto& branch=physical->branches_[task.branch];
          const auto& matrices=task.observable==0 ? branch.g
            : task.observable==1 ? branch.f00 : branch.f01;
          const auto& contraction=matrices[task.matrix_id];
          const int expanded=task.matrix_id/self.channel_count_;
          const int output=task.matrix_id%self.channel_count_;
          const double prefactor=-workspace.residues()(expanded,l)/
            (workspace.poles()(l)>=0.0
               ? workspace.kernel_zero_plus()(l)
               : workspace.kernel_zero_minus()(l));
          if (prefactor==0.0) continue;
          const MatrixXd& left=task.branch==0
            ? pair_transposed[1] : pair_transposed[0];
          const MatrixXd& right=task.branch==0
            ? reflected[0] : reflected[1];
          const MatrixXd mapped=contraction*right;
          MatrixXd& target=task.observable==0 ? accumulator.g_nodes
            : task.observable==1 ? accumulator.f00_nodes
                                 : accumulator.f01_nodes;
          target.col(output).array() += prefactor*
            (left.array()*mapped.array()).colwise().sum().transpose();
        }
      }
    }
  }

  static void add_self_values(
      const SelfEnergyPlan::CutPlan& cut_plan,
      int cut,
      int outer,
      int row_first,
      const MatrixXd& values,
      const VectorXd& external,
      const std::vector<double>& pseudo_degeneracy,
      NodeAccumulator& accumulator)
  {
    for (int j=0; j<values.cols(); ++j) {
      const int row=row_first+j;
      if (cut_plan.outer_pairs[row][0]!=outer)
        throw std::logic_error("noncontiguous OCA outer-pair work block");
      const int output=cut_plan.outer_pairs[row][1];
      const VectorXd increment=(external.array()*values.col(j).array()/
        pseudo_degeneracy[output]).matrix();
      accumulator.sigma_nodes.col(output)+=increment;
      if (accumulator.cuts_enabled)
        accumulator.cut_nodes[cut].col(output)+=increment;
    }
  }

#ifdef _MPI
  static void allreduce_accumulator(NodeAccumulator& accumulator,
                                    MPI_Comm communicator)
  {
    const auto sum=[](MatrixXd& values,MPI_Comm comm) {
      if (values.size()==0) return;
      if (values.size()>std::numeric_limits<int>::max())
        throw std::overflow_error("OCA DLR MPI reduction exceeds int count");
      MPI_Allreduce(MPI_IN_PLACE,values.data(),
                    static_cast<int>(values.size()),MPI_DOUBLE,MPI_SUM,comm);
    };
    sum(accumulator.sigma_nodes,communicator);
    if (accumulator.physical_enabled) {
      sum(accumulator.g_nodes,communicator);
      sum(accumulator.f00_nodes,communicator);
      sum(accumulator.f01_nodes,communicator);
    }
    if (accumulator.cuts_enabled)
      for (auto& nodes:accumulator.cut_nodes) sum(nodes,communicator);
  }
#endif
};

} // namespace edmft::oca_dlr

#endif
