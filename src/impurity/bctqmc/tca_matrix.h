// Block-matrix contraction prototype for third-order (TCA) atomic traces.
//
// The explicit TCA enumerator visits six independent atomic matrix indices
// j0,...,j5.  For a fixed topology and fixed directed bath triple, its atomic
// part has the form
//
//   sum_{j0...j5} F0(j0) G0(j0,j1) ... F5(j5) G5(j5,j0).
//
// This is exactly
//
//   Tr[(diag(F0) G0) ... (diag(F5) G5)].
//
// Performing the contraction as six dense block products avoids constructing
// the enormous scalar path catalog.  This header intentionally stops at that
// algebraic kernel: hybridization functions, DLR convolutions, topology signs,
// root/orientation multiplicities, and final symmetry factors remain outside.
// Keeping those factors outside makes it possible to validate the atomic
// contraction before committing to a convention for the two TCA topologies.
#ifndef EDMFT_TCA_MATRIX_H
#define EDMFT_TCA_MATRIX_H

#include "Eigen/Dense"
#include "read_cix.h"
#include "tca_topology.h"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace edmft::tca_matrix {

using PseudoValues=std::array<Eigen::VectorXd,6>;

inline blitz::Array<int,2> MakeTopology(int topology)
{
  return edmft::tca::MakeTopology(topology);
}

// One pseudo-propagator segment G_v.  A matrix element stores the signed,
// one-based Sindex label used throughout bctqmc:
//
//   0       forbidden matrix element,
//   +k/-k   +/- pseudo_values[v](k-1).
//
// We preserve the sign because off-diagonal pseudo Green's functions can be
// symmetry-related with a minus sign.
struct PropagatorMap {
  Eigen::MatrixXi signed_label;
};

// One block contributing to the atomic trace.  vertex_sign[v](j) is the
// matrix element of the operator at vertex v.  CIX encodes that matrix
// element in the sign of allowed[v](j,1), exactly as in the existing
// Topology_trial_move matrix implementation in imag_axis.h.
struct BlockPlan {
  int initial_block=-1;
  // Atomic block in which the pseudo propagator following each vertex lives.
  // The incoming block of vertex v is propagated_block[(v+5)%6].
  std::array<int,6> propagated_block{{-1,-1,-1,-1,-1,-1}};
  std::array<Eigen::VectorXd,6> vertex_sign;
  // Local atomic-state indices on the two sides of every operator vertex.
  // They are retained because a physical Bulla cut replaces the ordinary
  // one-body operator by a generally dense N_jm matrix.  In the ordinary
  // trace only matching entries are used and vertex_sign is sufficient.
  std::array<Eigen::VectorXi,6> incoming_state;
  std::array<Eigen::VectorXi,6> outgoing_state;
  std::array<PropagatorMap,6> propagator;
};

class AtomicTracePlan {
public:
  AtomicTracePlan()=default;

  // Compile all Pauli-allowed atomic blocks for one topology and one directed
  // bath triple.  Diag2Matrix is passed in so callers scanning many triples
  // can reuse its internal storage rather than reallocating the object.
  static AtomicTracePlan Compile(const Cix& cix,
                                 Diag2Matrix& diag2matrix,
                                 int topology,
                                 const std::array<int,3>& dirspin)
  {
    AtomicTracePlan plan;
    plan.topology_=topology;
    plan.dirspin_=dirspin;

    blitz::Array<int,2> diag=MakeTopology(topology);
    blitz::Array<int,1> ds(3);
    for (int line=0;line<3;++line) {
      if (dirspin[line]<0 || dirspin[line]>=2*cix.Nbath)
        throw std::out_of_range("TCA dirspin is outside [0,2*Nbath)");
      ds(line)=dirspin[line];
    }

    const auto& allowed=diag2matrix.Cmp_matrix(diag,ds,cix);
    plan.blocks_.reserve(allowed.size());

    for (const auto& [initial_block,alw]:allowed) {
      if (alw.size()!=6)
        throw std::logic_error("third-order Cmp_matrix did not return six vertices");

      bool empty=false;
      for (int vertex=0;vertex<6;++vertex)
        empty=empty || alw[vertex].extent(0)==0;
      if (empty) continue;

      BlockPlan block;
      block.initial_block=initial_block;

      for (int vertex=0;vertex<6;++vertex) {
        const int next=(vertex+1)%6;
        const int rows=alw[vertex].extent(0);
        const int cols=alw[next].extent(0);

        block.vertex_sign[vertex].resize(rows);
        block.incoming_state[vertex].resize(rows);
        block.outgoing_state[vertex].resize(rows);
        for (int row=0;row<rows;++row) {
          block.vertex_sign[vertex](row)=
            alw[vertex](row,1)>=0 ? 1.0 : -1.0;
          block.incoming_state[vertex](row)=
            cix.inBlock(std::abs(alw[vertex](row,0)));
          block.outgoing_state[vertex](row)=
            cix.inBlock(std::abs(alw[vertex](row,1)));
        }

        auto& labels=block.propagator[vertex].signed_label;
        labels=Eigen::MatrixXi::Zero(rows,cols);

        // Both states belong to the block immediately after operator vertex.
        // iBlck is one-based, whereas Sgm_symmetry is indexed from zero.
        const int state_after=std::abs(alw[vertex](0,1));
        const int propagated_block=cix.iBlck(state_after)-1;
        if (propagated_block<0 || propagated_block>=cix.Nblock)
          throw std::logic_error("invalid atomic block returned by Cmp_matrix");
        block.propagated_block[vertex]=propagated_block;

        for (int row=0;row<rows;++row) {
          const int a=cix.inBlock(std::abs(alw[vertex](row,1)));
          for (int col=0;col<cols;++col) {
            const int b=cix.inBlock(std::abs(alw[next](col,0)));
            const int label=cix.Sgm_symmetry[propagated_block](a,b);
            labels(row,col)=label;
            plan.maximum_pseudo_label_=
              std::max(plan.maximum_pseudo_label_,std::abs(label));
          }
        }
      }

      // The six rectangular factors must form a closed ring.  Checking this
      // here turns an obscure Eigen assertion into a useful CIX diagnostic.
      for (int vertex=0;vertex<6;++vertex) {
        const int next=(vertex+1)%6;
        if (block.propagator[vertex].signed_label.cols()!=
            block.propagator[next].signed_label.rows())
          throw std::logic_error("TCA atomic block does not form a closed matrix ring");
      }

      plan.blocks_.push_back(std::move(block));
    }
    return plan;
  }

  bool empty() const noexcept { return blocks_.empty(); }
  int topology() const noexcept { return topology_; }
  const std::array<int,3>& dirspin() const noexcept { return dirspin_; }
  int maximum_pseudo_label() const noexcept { return maximum_pseudo_label_; }
  const std::vector<BlockPlan>& blocks() const noexcept { return blocks_; }

  // Number of scalar tuples that the reference implementation would visit
  // before zero-Sindex pruning.  Saturation avoids unsigned overflow in a
  // deliberately diagnostic quantity.
  std::uint64_t scalar_cartesian_terms() const noexcept
  {
    std::uint64_t total=0;
    for (const BlockPlan& block:blocks_) {
      std::uint64_t product=1;
      for (int vertex=0;vertex<6;++vertex) {
        const auto size=static_cast<std::uint64_t>(
          block.vertex_sign[vertex].size());
        if (size!=0 && product>std::numeric_limits<std::uint64_t>::max()/size) {
          product=std::numeric_limits<std::uint64_t>::max();
          break;
        }
        product*=size;
      }
      if (total>std::numeric_limits<std::uint64_t>::max()-product)
        return std::numeric_limits<std::uint64_t>::max();
      total+=product;
    }
    return total;
  }

  // Production candidate: evaluate the same scalar sum by six dense matrix
  // multiplications per contributing atomic block.  No scalar path list is
  // allocated.  Vertex signs are folded into the rows of each G matrix.
  double EvaluateMatrix(const PseudoValues& pseudo_values) const
  {
    ValidatePseudoValues(pseudo_values);
    double result=0.0;
    for (const BlockPlan& block:blocks_) {
      Eigen::MatrixXd product=Factor(block,0,pseudo_values[0]);
      for (int vertex=1;vertex<6;++vertex)
        product=(product*Factor(block,vertex,pseudo_values[vertex])).eval();
      result+=product.trace();
    }
    return result;
  }

  // Validation only: direct six-index summation of the identical expression.
  // The hard term limit prevents accidental use on the multi-billion-path SC
  // case that motivated this prototype.
  double EvaluateScalarReference(
    const PseudoValues& pseudo_values,
    std::uint64_t maximum_terms=1000000ULL) const
  {
    ValidatePseudoValues(pseudo_values);
    const std::uint64_t terms=scalar_cartesian_terms();
    if (terms>maximum_terms)
      throw std::length_error(
        "TCA scalar reference requires "+std::to_string(terms)+
        " Cartesian terms; increase the diagnostic limit explicitly");

    double result=0.0;
    for (const BlockPlan& block:blocks_) {
      std::array<int,6> index{{0,0,0,0,0,0}};
      std::function<void(int,double)> visit;
      visit=[&](int vertex,double weight) {
        if (vertex==6) {
          const int label=block.propagator[5].signed_label(index[5],index[0]);
          if (label!=0)
            result+=weight*SignedPseudoValue(pseudo_values[5],label);
          return;
        }

        const int count=static_cast<int>(block.vertex_sign[vertex].size());
        for (int j=0;j<count;++j) {
          index[vertex]=j;
          double next_weight=weight*block.vertex_sign[vertex](j);
          if (vertex>0) {
            const int label=block.propagator[vertex-1].signed_label(
              index[vertex-1],index[vertex]);
            if (label==0) continue;
            next_weight*=SignedPseudoValue(pseudo_values[vertex-1],label);
          }
          visit(vertex+1,next_weight);
        }
      };
      visit(0,1.0);
    }
    return result;
  }

private:
  static double SignedPseudoValue(const Eigen::VectorXd& values,int label)
  {
    const double value=values(std::abs(label)-1);
    return label>0 ? value : -value;
  }

  static Eigen::MatrixXd Factor(const BlockPlan& block,
                                int vertex,
                                const Eigen::VectorXd& values)
  {
    const Eigen::MatrixXi& labels=
      block.propagator[vertex].signed_label;
    Eigen::MatrixXd factor=Eigen::MatrixXd::Zero(labels.rows(),labels.cols());
    for (int row=0;row<labels.rows();++row) {
      const double vertex_value=block.vertex_sign[vertex](row);
      for (int col=0;col<labels.cols();++col) {
        const int label=labels(row,col);
        if (label!=0)
          factor(row,col)=vertex_value*SignedPseudoValue(values,label);
      }
    }
    return factor;
  }

  void ValidatePseudoValues(const PseudoValues& pseudo_values) const
  {
    for (int vertex=0;vertex<6;++vertex) {
      if (pseudo_values[vertex].size()<maximum_pseudo_label_)
        throw std::invalid_argument(
          "TCA pseudo-value vector is shorter than its largest Sindex label");
    }
  }

  int topology_=-1;
  std::array<int,3> dirspin_{{-1,-1,-1}};
  int maximum_pseudo_label_=0;
  std::vector<BlockPlan> blocks_;
};

} // namespace edmft::tca_matrix

#endif // EDMFT_TCA_MATRIX_H
