// Reusable DLR/AAA decomposition primitives for third-order diagrams.
//
// This header contains both the small, testable decoration tables of
// Eqs. (28)-(31) of SciPost Phys. 19, 121 (2025) and their production matrix-
// valued convolution engine.  It evaluates all six pseudo self-energy cuts
// and all three physical G/F cuts without materializing scalar atomic paths.
#ifndef EDMFT_TCA_DLR_H
#define EDMFT_TCA_DLR_H

#include "imag_axis_dlr.h"
#include "tca_matrix.h"
#include "tca_topology.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <map>
#include <numeric>
#include <ostream>
#include <set>
#include <stdexcept>
#include <vector>

#ifdef _MPI
#include <mpi.h>
#endif

namespace edmft::tca_dlr {

using Eigen::MatrixXd;
using Eigen::VectorXd;
using MatrixFunction=std::vector<MatrixXd>;

enum class PlacementKind {
  endpoint,
  segment
};

// A kernel sign is relative to the real AAA pole omega_l:
//   +1 -> K(t,+omega_l),  -1 -> K(t,-omega_l),  0 -> no factor.
//
// vertex_kernel[v] decorates the matrix/operator multiplication at vertex v.
// edge_kernel[e] decorates pseudo propagator G_e between vertices e and e+1.
// Multiple separated lines may decorate the same vertex or edge; execution
// then multiplies all their node-value factors before a DLR convolution.
struct LineDecoration {
  int topology=-1;
  int line=-1;
  int direction=0;
  int first_vertex=-1;
  int second_vertex=-1;
  PlacementKind kind=PlacementKind::segment;
  std::array<int,6> vertex_kernel{{0,0,0,0,0,0}};
  std::array<int,6> edge_kernel{{0,0,0,0,0,0}};
  int denominator_kernel=0;
  int denominator_power=0;

  int edge_count() const noexcept { return second_vertex-first_vertex; }
};

// Construct the stable decomposition for one of the two non-root lines.
// direction=0 means Delta(t_b-t_a); direction=1 means
// Delta(beta-t_b+t_a), matching dirspin=2*orbital+direction.
//
// With the solver's signed DLR kernel, every physical Delta_code line has the
// same single overall minus from Delta_code=-sum_l R_l K_l.  That minus is
// returned separately by code_pole_prefactor() and is NOT encoded as a
// direction-dependent sign here.
inline LineDecoration MakeRootedLineDecoration(int first_vertex,
                                               int second_vertex,
                                               int direction,
                                               double pole)
{
  if (direction!=0 && direction!=1)
    throw std::out_of_range("TCA hybridization direction must be zero or one");

  LineDecoration result;
  result.direction=direction;
  result.first_vertex=first_vertex;
  result.second_vertex=second_vertex;
  const int a=result.first_vertex;
  const int b=result.second_vertex;
  if (!(0<=a && a<b && b<6))
    throw std::logic_error("canonical TCA line does not satisfy 0 <= a < b < 6");

  // Forward/nonpositive and backward/positive poles use the two-endpoint
  // identity.  This avoids constructing exponentially small segment factors
  // followed by division by an equally small K(0).
  const bool endpoint=
    (direction==0 && pole<=0.0) || (direction==1 && pole>0.0);
  if (endpoint) {
    result.kind=PlacementKind::endpoint;
    if (direction==0) {
      // K(t_b-t_a,w)=K(t_b,w) K(t_a,-w)/K(0,-w).
      result.vertex_kernel[a]=-1;
      result.vertex_kernel[b]=+1;
      result.denominator_kernel=-1;
    } else {
      // K(beta-t_b+t_a,w)=K(t_b,-w) K(t_a,w)/K(0,w).
      result.vertex_kernel[a]=+1;
      result.vertex_kernel[b]=-1;
      result.denominator_kernel=+1;
    }
    result.denominator_power=1;
    return result;
  }

  // The complementary stable branch uses every backbone edge between the
  // endpoints.  If q=b-a edges are multiplied, q-1 joins are removed by the
  // denominator K(0)^(q-1), exactly as in Eq. (28).  A backward line uses
  // the reflected frequency -omega.
  result.kind=PlacementKind::segment;
  const int effective_frequency=direction==0 ? +1 : -1;
  for (int edge=a;edge<b;++edge)
    result.edge_kernel[edge]=effective_frequency;
  result.denominator_kernel=effective_frequency;
  result.denominator_power=(b-a)-1;
  return result;
}

inline LineDecoration MakeLineDecoration(int topology,
                                         int line,
                                         int direction,
                                         double pole)
{
  if (topology<0 || topology>=2)
    throw std::out_of_range("TCA topology index must be zero or one");
  if (line<=edmft::tca::root_line || line>=3)
    throw std::out_of_range("only non-root TCA lines may be separated");
  auto result=MakeRootedLineDecoration(
    edmft::tca::topology_endpoints[topology][line][0],
    edmft::tca::topology_endpoints[topology][line][1],direction,pole);
  result.topology=topology;
  result.line=line;
  return result;
}

// Stable scalar multiplying the decorated backbone for a single physical
// hybridization pole.  kernel_zero must be K(0,+/-omega_l), selected according
// to decoration.denominator_kernel.  The minus is the bctqmc Fourier
// convention and is independent of propagation direction.
inline double CodePolePrefactor(double residue,
                                double kernel_zero,
                                const LineDecoration& decoration)
{
  if (decoration.denominator_kernel==0)
    throw std::logic_error("TCA line decoration has no denominator kernel");
  return -residue/std::pow(kernel_zero,decoration.denominator_power);
}

// The four backbone types generated by two separated hybridization lines.
// Pole values, rather than only their signs, are accepted to make the zero
// convention unambiguous: zero belongs to the nonpositive branch.
struct PolePairDecoration {
  int topology=-1;
  std::array<int,3> direction{{0,0,0}};
  std::array<LineDecoration,2> split;
};

inline PolePairDecoration MakePolePairDecoration(
    int topology,int direction1,int direction2,double pole1,double pole2)
{
  PolePairDecoration result;
  result.topology=topology;
  result.direction={{0,direction1,direction2}};
  result.split={{
    MakeLineDecoration(topology,1,direction1,pole1),
    MakeLineDecoration(topology,2,direction2,pole2)
  }};
  return result;
}

// One self-energy cut written as an open, time-ordered backbone.  Local
// vertex zero is the vertex immediately after the removed pseudo edge;
// local vertex five is immediately before it.  original_vertex and
// original_edge map every local object back to AtomicTracePlan storage.
struct RootedCut {
  struct Line {
    int original_line=-1;
    int first_vertex=-1;
    int second_vertex=-1;
    int direction=0;
  };

  int topology=-1;
  int cut_edge=-1;
  std::array<int,6> original_vertex{{0,1,2,3,4,5}};
  std::array<int,6> original_edge{{0,1,2,3,4,5}};
  Line root;
  std::array<Line,2> split;
};

// Rotate an arbitrary pseudo cut to the common open-backbone convention.
// If an endpoint order reverses across the new time origin, its direction bit
// is flipped so the represented Delta is unchanged.
inline RootedCut MakeRootedCut(int topology,
                               const std::array<int,3>& direction,
                               int cut_edge)
{
  if (topology<0 || topology>=2)
    throw std::out_of_range("TCA topology index must be zero or one");
  if (cut_edge<0 || cut_edge>=6)
    throw std::out_of_range("TCA pseudo cut must be in [0,6)");
  for (const int d:direction)
    if (d!=0 && d!=1)
      throw std::out_of_range("TCA hybridization direction must be zero or one");

  RootedCut result;
  result.topology=topology;
  result.cut_edge=cut_edge;
  const int start=(cut_edge+1)%6;
  for (int local=0;local<6;++local) {
    result.original_vertex[local]=(start+local)%6;
    result.original_edge[local]=(start+local)%6;
  }

  std::array<RootedCut::Line,3> lines;
  int root_count=0,split_count=0;
  for (int line=0;line<3;++line) {
    int a=(edmft::tca::topology_endpoints[topology][line][0]-start+6)%6;
    int b=(edmft::tca::topology_endpoints[topology][line][1]-start+6)%6;
    int d=direction[line];
    if (a>b) {
      std::swap(a,b);
      d^=1;
    }
    const RootedCut::Line rotated{line,a,b,d};
    lines[line]=rotated;
    if (a==0) {
      result.root=rotated;
      ++root_count;
    } else {
      if (split_count>=2)
        throw std::logic_error("TCA rooted cut has more than two split lines");
      result.split[split_count++]=rotated;
    }
  }
  if (root_count!=1 || split_count!=2)
    throw std::logic_error(
      "TCA rooted cut must have one root and two separated lines");
  return result;
}

inline int SignedIndexSign(int value)
{
  return value>=0 ? +1 : -1;
}

inline double KernelZero(const imag_dlr::Workspace& workspace,
                         int pole,
                         int frequency_sign)
{
  return frequency_sign>=0
    ? workspace.kernel_zero_plus()(pole)
    : workspace.kernel_zero_minus()(pole);
}

// First executable TCA self-energy kernel.  It evaluates one fixed atomic
// block plan and any one of its six pseudo-propagator cuts entirely on the DLR
// grid.  It is intentionally serial and fine-grained; later plan compilation
// will batch equal matrix shapes and distribute (bath triple,pole pair,cut)
// work without changing this algebra.
class RootedSelfEnergyEvaluator {
public:
  using PreparedPlan=std::vector<std::array<MatrixFunction,6>>;

  static MatrixXd evaluate_cut_nodes(
      const tca_matrix::AtomicTracePlan& plan,
      const Cix& cix,
      const imag_dlr::Workspace& workspace,
      const std::vector<double>& pseudo_degeneracy,
      int cut_edge,
      const PreparedPlan* prepared_plan=nullptr)
  {
    Validate(plan,cix,workspace,pseudo_degeneracy);
    PreparedPlan local_prepared;
    if (!prepared_plan) {
      local_prepared=PreparePropagators(plan,workspace);
      prepared_plan=&local_prepared;
    }
    std::array<int,3> direction{};
    std::array<int,3> orbital{};
    std::array<int,3> channel{};
    double bath_sign=1.0;
    for (int line=0;line<3;++line) {
      direction[line]=plan.dirspin()[line]&1;
      orbital[line]=plan.dirspin()[line]/2;
      const int signed_channel=cix.bath_sigind0(orbital[line]);
      channel[line]=std::abs(signed_channel);
      bath_sign*=SignedIndexSign(signed_channel);
    }
    const RootedCut cut=MakeRootedCut(plan.topology(),direction,cut_edge);
    const auto& basis=workspace.basis();
    const int rank=workspace.rank();
    MatrixXd output=MatrixXd::Zero(rank,workspace.pseudo_count());

    // The two AAA sums are independent.  This is also the natural MPI work
    // unit once correctness is established.
    for (int l0=0;l0<workspace.pole_count();++l0)
      for (int l1=0;l1<workspace.pole_count();++l1) {
        const std::array<int,2> pole{{l0,l1}};
        std::array<LineDecoration,2> decoration;
        double scalar_prefactor=
          bath_sign*edmft::tca::topology_crossing_sign[plan.topology()]*
          edmft::tca::topology_symmetry_factor[plan.topology()];

        std::array<VectorXd,6> vertex_factor;
        std::array<VectorXd,5> edge_factor;
        for (auto& factor:vertex_factor)
          factor=VectorXd::Ones(rank);
        for (auto& factor:edge_factor)
          factor=VectorXd::Ones(rank);

        for (int split=0;split<2;++split) {
          const auto& line=cut.split[split];
          const int original=line.original_line;
          const double omega=workspace.poles()(pole[split]);
          decoration[split]=MakeRootedLineDecoration(
            line.first_vertex,line.second_vertex,line.direction,omega);
          decoration[split].topology=plan.topology();
          decoration[split].line=original;
          const double residue=workspace.residues()(
            channel[original],pole[split]);
          scalar_prefactor*=CodePolePrefactor(
            residue,KernelZero(workspace,pole[split],
                               decoration[split].denominator_kernel),
            decoration[split]);

          for (int vertex=0;vertex<6;++vertex) {
            const int sign=decoration[split].vertex_kernel[vertex];
            if (sign==0) continue;
            vertex_factor[vertex].array()*=
              (sign>0 ? workspace.kernel_plus().col(pole[split])
                      : workspace.kernel_minus().col(pole[split])).array();
          }
          for (int edge=0;edge<5;++edge) {
            const int sign=decoration[split].edge_kernel[edge];
            if (sign==0) continue;
            edge_factor[edge].array()*=
              (sign>0 ? workspace.kernel_plus().col(pole[split])
                      : workspace.kernel_minus().col(pole[split])).array();
          }
        }

        // The unseparated line is evaluated as a multiplication at its second
        // endpoint.  delta_plus/minus are standard kernel sums, hence the
        // explicit minus converting them to Delta_code.
        const int root_channel=channel[cut.root.original_line];
        const VectorXd root_delta=cut.root.direction==0
          ? VectorXd(workspace.delta_plus().col(root_channel))
          : VectorXd(workspace.delta_minus().col(root_channel));
        vertex_factor[cut.root.second_vertex].array()*=-root_delta.array();

        for (int block_index=0;
             block_index<static_cast<int>(plan.blocks().size());
             ++block_index) {
          const auto& block=plan.blocks()[block_index];
          std::array<MatrixFunction,5> edge_coefficients;
          // Edge zero starts the ordered product in value form.  Converting
          // it to coefficients only to immediately convert it back at the
          // interpolation nodes is an exact algebraic no-op and was a
          // substantial cost for large atomic blocks.
          MatrixFunction values=
            (*prepared_plan)[block_index][cut.original_edge[0]];
          for (int node=0;node<rank;++node)
            values[node]*=edge_factor[0](node);
          for (int local_edge=1;local_edge<5;++local_edge) {
            const int original_edge=cut.original_edge[local_edge];
            MatrixFunction edge_values=
              (*prepared_plan)[block_index][original_edge];
            for (int node=0;node<rank;++node)
              edge_values[node]*=edge_factor[local_edge](node);
            edge_coefficients[local_edge]=
              basis.matrix_values_to_coefficients(edge_values);
          }

          // Start with D_0 G_0 D_1.  Subsequent steps convolve with the next
          // interval propagator on the right and then apply D_v at the new
          // absolute-time vertex.  This preserves the atomic trace order.
          ApplyLeftVertex(values,
            block.vertex_sign[cut.original_vertex[0]],
            vertex_factor[0]);
          ApplyRightVertex(values,
            block.vertex_sign[cut.original_vertex[1]],vertex_factor[1]);
          MatrixFunction coefficients=
            basis.matrix_values_to_coefficients(values);

          for (int local_edge=1;local_edge<5;++local_edge) {
            values=basis.time_ordered_convolution_values(
              coefficients,edge_coefficients[local_edge]);
            ApplyRightVertex(values,
              block.vertex_sign[cut.original_vertex[local_edge+1]],
              vertex_factor[local_edge+1]);
            if (local_edge<4)
              coefficients=basis.matrix_values_to_coefficients(values);
          }

          AccumulateCut(values,
            block.propagator[cut.original_edge[5]].signed_label,
            scalar_prefactor,pseudo_degeneracy,output);
        }
      }
    return output;
  }

  static std::array<MatrixXd,6> evaluate_all_cut_nodes(
      const tca_matrix::AtomicTracePlan& plan,
      const Cix& cix,
      const imag_dlr::Workspace& workspace,
      const std::vector<double>& pseudo_degeneracy,
      std::array<double,6>* cut_seconds=nullptr,
      double* prepare_seconds=nullptr)
  {
    std::array<MatrixXd,6> cuts;
    const auto prepare_begin=std::chrono::steady_clock::now();
    const PreparedPlan prepared=PreparePropagators(plan,workspace);
    if (prepare_seconds)
      *prepare_seconds=std::chrono::duration<double>(
        std::chrono::steady_clock::now()-prepare_begin).count();
    if (cut_seconds) cut_seconds->fill(0.0);
    for (int cut=0;cut<6;++cut) {
      const auto cut_begin=std::chrono::steady_clock::now();
      cuts[cut]=evaluate_cut_nodes(
        plan,cix,workspace,pseudo_degeneracy,cut,&prepared);
      if (cut_seconds)
        (*cut_seconds)[cut]=std::chrono::duration<double>(
          std::chrono::steady_clock::now()-cut_begin).count();
    }
    return cuts;
  }

private:
  static PreparedPlan PreparePropagators(
      const tca_matrix::AtomicTracePlan& plan,
      const imag_dlr::Workspace& workspace)
  {
    PreparedPlan prepared(plan.blocks().size());
    for (int block=0;block<static_cast<int>(plan.blocks().size());++block)
      for (int edge=0;edge<6;++edge)
        prepared[block][edge]=PropagatorValues(
          plan.blocks()[block].propagator[edge].signed_label,workspace);
    return prepared;
  }

  static void Validate(const tca_matrix::AtomicTracePlan& plan,
                       const Cix& cix,
                       const imag_dlr::Workspace& workspace,
                       const std::vector<double>& pseudo_degeneracy)
  {
    if (!workspace.has_pseudos())
      throw std::logic_error("TCA DLR evaluation requires pseudo coefficients");
    if (workspace.pseudo_count()!=cix.Npseudo_unique ||
        workspace.channel_count()!=cix.Nbath_unique_unique)
      throw std::invalid_argument("TCA plan/workspace dimension mismatch");
    if (static_cast<int>(pseudo_degeneracy.size())!=cix.Npseudo_unique)
      throw std::invalid_argument("TCA pseudo degeneracy size mismatch");
    for (const double degeneracy:pseudo_degeneracy)
      if (!(degeneracy>0.0))
        throw std::invalid_argument("TCA pseudo degeneracy must be positive");
    if (plan.empty())
      throw std::invalid_argument("cannot evaluate an empty TCA atomic plan");
  }

  static MatrixFunction PropagatorValues(
      const Eigen::MatrixXi& signed_label,
      const imag_dlr::Workspace& workspace)
  {
    MatrixFunction values(
      workspace.rank(),MatrixXd::Zero(signed_label.rows(),signed_label.cols()));
    for (int i=0;i<signed_label.rows();++i)
      for (int j=0;j<signed_label.cols();++j) {
        const int label=signed_label(i,j);
        if (label==0) continue;
        const int pseudo=std::abs(label)-1;
        for (int node=0;node<workspace.rank();++node)
          values[node](i,j)=SignedIndexSign(label)*
            workspace.pseudo_values()(node,pseudo);
      }
    return values;
  }

  static void ApplyLeftVertex(MatrixFunction& values,
                              const VectorXd& atomic_sign,
                              const VectorXd& scalar)
  {
    for (int node=0;node<static_cast<int>(values.size());++node) {
      if (values[node].rows()!=atomic_sign.size())
        throw std::logic_error("TCA left vertex/matrix dimension mismatch");
      values[node].array().colwise()*=atomic_sign.array();
      values[node]*=scalar(node);
    }
  }

  static void ApplyRightVertex(MatrixFunction& values,
                               const VectorXd& atomic_sign,
                               const VectorXd& scalar)
  {
    for (int node=0;node<static_cast<int>(values.size());++node) {
      if (values[node].cols()!=atomic_sign.size())
        throw std::logic_error("TCA right vertex/matrix dimension mismatch");
      values[node].array().rowwise()*=atomic_sign.transpose().array();
      values[node]*=scalar(node);
    }
  }

  static void AccumulateCut(const MatrixFunction& values,
                            const Eigen::MatrixXi& cut_label,
                            double prefactor,
                            const std::vector<double>& pseudo_degeneracy,
                            MatrixXd& output)
  {
    if (values.empty() || cut_label.rows()!=values[0].cols() ||
        cut_label.cols()!=values[0].rows())
      throw std::logic_error("TCA cut closure has incompatible dimensions");
    for (int cut_row=0;cut_row<cut_label.rows();++cut_row)
      for (int cut_col=0;cut_col<cut_label.cols();++cut_col) {
        const int label=cut_label(cut_row,cut_col);
        if (label==0) continue;
        const int pseudo=std::abs(label)-1;
        const double weight=prefactor*SignedIndexSign(label)/
          pseudo_degeneracy[pseudo];
        for (int node=0;node<static_cast<int>(values.size());++node)
          output(node,pseudo)+=weight*values[node](cut_col,cut_row);
      }
  }
};

// Self-energy kernel with the bath carried by the unexpanded (root) line
// summed *inside* the atomic matrix contraction.
//
// The fixed-plan evaluator above applies the operator at local time zero
// before doing any convolution.  That makes the whole backbone depend on the
// root bath and forces an outer Nbath loop.  Algebraically that operator is a
// constant boundary matrix and can be deferred.  We therefore build
//
//   G0 O1^T * G1 ... * G_{b-1}
//
// with both root endpoints open.  At vertex b we form the correlated sum
//
//   sum_a O0(a)^T [backbone] Ob(a)^T [-Delta_a(t_b)] ,
//
// and only then continue along the common remainder of the backbone.  This
// is Step 2 of Sec. 4.3 of the SciPost algorithm.  The number of expensive
// DLR backbones is consequently proportional to Nbath^2 rather than
// Nbath^3.
//
// Atomic symmetry blocks are kept explicitly.  A SparseFunction block (i,j)
// maps atomic block j on its right boundary to atomic block i on its left
// boundary.  This is essential for the 256-state cluster example: its dense
// Hilbert space is large, while its largest symmetry block has dimension 22.
class SummedRootSelfEnergyEvaluator {
public:
  using PseudoBlocks=std::vector<MatrixFunction>;
  // A pseudo propagator multiplied by a pole kernel is first known at the
  // representative imaginary-time nodes.  Convolutions need its DLR
  // coefficients.  In the original bath-major loop that conversion was
  // repeated for every (bath0,bath1) pair even though the multiplier depends
  // only on the cut geometry and the two AAA poles.  This small holder lets a
  // pole-major caller prepare the conversion once and reuse it for every bath
  // pair in the current pole tile.  The undressed propagator aliases the
  // coefficients already owned by Workspace and therefore allocates nothing.
  struct PreparedPseudoCoefficients {
    bool use_workspace_coefficients=true;
    MatrixXd dressed;
    MatrixXd hilbert;
    // Atomic-block expansions are independent of the bath operators.  They
    // are filled lazily and then reused by every bath pair for this pole
    // pair.  The cache lifetime is one pole pair, so memory stays bounded
    // even for large f-shell calculations.
    mutable std::vector<MatrixFunction> coefficient_blocks;
    mutable std::vector<MatrixFunction> hilbert_blocks;
    mutable std::vector<unsigned char> block_ready;
    mutable std::vector<unsigned char> block_nonzero;
  };

  struct BathPairNodeResult {
    MatrixXd nodes;
    std::uint64_t nonempty_pairs=0;
  };
  struct SparseMatrixEntry {
    int row=0;
    int column=0;
    double value=0.0;
  };
  struct OperatorEntry {
    int source=-1;
    int target=-1;
    MatrixXd matrix; // O^T: source rows, target columns
    // Ordinary creation/annihilation operators are signed partial
    // permutations: every source state has at most one target state.  Store
    // that exact map so hot contractions can gather rows/columns in O(n^2)
    // instead of multiplying by a mostly-zero dense matrix in O(n^3).
    // Bulla connectors may be genuinely dense and leave this flag false.
    bool signed_partial_permutation=false;
    std::vector<int> target_column;
    std::vector<double> row_sign;
    // Exact nonzero entries are also retained for terminal trace
    // contractions.  In the superconducting cluster, even Bulla connectors
    // classified as general matrices contain only about 18% nonzeros.  A
    // sparse bilinear trace can therefore avoid constructing a dense matrix
    // product whose off-trace elements would immediately be discarded.
    std::vector<SparseMatrixEntry> nonzero;
  };
  struct OperatorTranspose {
    std::vector<OperatorEntry> entries;
    std::vector<int> by_source;
    std::vector<std::vector<int>> by_target;
  };
  struct StoredRootPair {
    OperatorTranspose zero;
    OperatorTranspose endpoint;
    int channel=0;
    int direction=0;
    double bath_sign=1.0;
  };

  // orbital[root.original_line] is ignored and may be -1.  The other two
  // entries specify the explicitly enumerated split-line baths.  direction
  // uses the original, unrotated topology convention (line zero fixed to
  // direction zero by the caller).
  static MatrixXd evaluate_cut_nodes(
      const Cix& cix,
      const imag_dlr::Workspace& workspace,
      const std::vector<double>& pseudo_degeneracy,
      int topology,
      int cut_edge,
      const std::array<int,3>& orbital,
      const std::array<int,3>& direction,
      bool sum_free_root_direction=false)
  {
    Validate(cix,workspace,pseudo_degeneracy,topology,cut_edge,
             orbital,direction);
    const RootedCut cut=MakeRootedCut(topology,direction,cut_edge);
    const int root_line=cut.root.original_line;
    const int root_vertex=cut.root.second_vertex;
    const int rank=workspace.rank();
    MatrixXd output=MatrixXd::Zero(rank,workspace.pseudo_count());
    const PseudoBlocks pseudo_blocks=PreparePseudoBlocks(cix,workspace);

    // Operators on split-line vertices do not depend on the AAA poles.
    // Build their transposes once.  O^T is the matrix appearing in the
    // left-to-right atomic product used by AtomicTracePlan.
    std::array<OperatorTranspose,6> vertex_operator;
    std::array<bool,6> has_operator{{false,false,false,false,false,false}};
    for (const auto& line:cut.split) {
      const int bath=orbital[line.original_line];
      const int ds=2*bath+line.direction;
      vertex_operator[line.first_vertex]=
        MakeOperatorTranspose(cix,Opposite(ds));
      vertex_operator[line.second_vertex]=
        MakeOperatorTranspose(cix,ds);
      has_operator[line.first_vertex]=true;
      has_operator[line.second_vertex]=true;
    }

    // Root operator pairs and their channel metadata are likewise static.
    //
    // The canonical Phi convention fixes the direction of original line
    // zero.  When line one or two becomes the unexpanded line after a pseudo
    // cut, however, both of its original directions must be included.  Those
    // two directions change only this local root contraction: neither split
    // line, its AAA decoration, nor any convolution along the backbone is
    // affected.  Production evaluation therefore puts both directions in
    // this list and evaluates the expensive pole-pair backbone only once.
    // The optional flag preserves the old one-direction behavior for focused
    // tests and other fine-grained callers.
    std::vector<StoredRootPair> roots;
    const int root_direction_count=
      sum_free_root_direction && root_line!=edmft::tca::root_line ? 2 : 1;
    roots.reserve(root_direction_count*cix.Nbath);
    for (int root_choice=0;root_choice<root_direction_count;++root_choice) {
      int local_root_direction=cut.root.direction;
      if (root_direction_count==2) {
        auto root_variant_direction=direction;
        root_variant_direction[root_line]=root_choice;
        const RootedCut root_variant=MakeRootedCut(
          topology,root_variant_direction,cut_edge);
        if (root_variant.root.original_line!=root_line ||
            root_variant.root.first_vertex!=cut.root.first_vertex ||
            root_variant.root.second_vertex!=cut.root.second_vertex)
          throw std::logic_error(
            "TCA root direction unexpectedly changed cut geometry");
        local_root_direction=root_variant.root.direction;
      }
      for (int bath=0;bath<cix.Nbath;++bath) {
        const int ds=2*bath+local_root_direction;
        const int signed_channel=cix.bath_sigind0(bath);
        roots.push_back({MakeOperatorTranspose(cix,Opposite(ds)),
                         MakeOperatorTranspose(cix,ds),
                         std::abs(signed_channel),
                         local_root_direction,
                         static_cast<double>(SignedIndexSign(signed_channel))});
      }
    }

    // The signs attached to the two explicitly enumerated bath entries do
    // not depend on pole number.  The root sign is included term-by-term in
    // ContractRoot because its bath is summed there.
    double fixed_bath_sign=1.0;
    for (const auto& line:cut.split)
      fixed_bath_sign*=SignedIndexSign(
        cix.bath_sigind0(orbital[line.original_line]));

    for (int l0=0;l0<workspace.pole_count();++l0)
      for (int l1=0;l1<workspace.pole_count();++l1) {
        const std::array<int,2> pole{{l0,l1}};
        std::array<VectorXd,6> vertex_factor;
        std::array<VectorXd,5> edge_factor;
        for (auto& x:vertex_factor) x=VectorXd::Ones(rank);
        for (auto& x:edge_factor) x=VectorXd::Ones(rank);
        double prefactor=
          fixed_bath_sign*edmft::tca::topology_crossing_sign[topology]*
          edmft::tca::topology_symmetry_factor[topology];

        for (int split=0;split<2;++split) {
          const auto& line=cut.split[split];
          const double omega=workspace.poles()(pole[split]);
          LineDecoration decoration=MakeRootedLineDecoration(
            line.first_vertex,line.second_vertex,line.direction,omega);
          const int signed_channel=
            cix.bath_sigind0(orbital[line.original_line]);
          prefactor*=CodePolePrefactor(
            workspace.residues()(std::abs(signed_channel),pole[split]),
            KernelZero(workspace,pole[split],
                       decoration.denominator_kernel),decoration);
          for (int vertex=0;vertex<6;++vertex) {
            const int sign=decoration.vertex_kernel[vertex];
            if (sign!=0)
              vertex_factor[vertex].array()*=
                (sign>0 ? workspace.kernel_plus().col(pole[split])
                        : workspace.kernel_minus().col(pole[split])).array();
          }
          for (int edge=0;edge<5;++edge) {
            const int sign=decoration.edge_kernel[edge];
            if (sign!=0)
              edge_factor[edge].array()*=
                (sign>0 ? workspace.kernel_plus().col(pole[split])
                        : workspace.kernel_minus().col(pole[split])).array();
          }
        }

        // Leave O0^T unapplied.  The diagonal G0 blocks include every
        // possible block reached by a root operator, and the later root
        // contraction selects only compatible operator paths.
        SparseFunction chain=PseudoValues(pseudo_blocks,edge_factor[0]);
        for (int vertex=1;vertex<6 && !chain.empty();++vertex) {
          if (vertex==root_vertex) {
            chain=ContractRoot(chain,roots,workspace,
                               vertex_factor[vertex]);
          } else {
            if (!has_operator[vertex])
              throw std::logic_error(
                "non-root TCA vertex has no split-line operator");
            chain=MultiplyRight(chain,vertex_operator[vertex],
                                vertex_factor[vertex]);
          }
          if (vertex<5 && !chain.empty())
            chain=ConvolveRightPseudo(chain,cix,workspace,
                                      edge_factor[vertex]);
        }
        AccumulateCut(chain,cix,prefactor,pseudo_degeneracy,output);
      }
    return output;
  }

  // Pole-major counterpart of evaluate_cut_nodes().  bath_pairs contains the
  // two baths in the same order as cut.split[].  All pairs share the rooted
  // geometry, pole kernels, root contraction, and (most importantly) the DLR
  // values-to-coefficients conversions.  Only the sparse impurity operators,
  // AAA residues, and signed bath channels vary between pairs.
  //
  // Keeping this as a separate entry point is intentional: the single-pair
  // routine above remains a simple numerical reference while this routine is
  // validated in production-sized normal and superconducting block spaces.
  static BathPairNodeResult evaluate_cut_bath_pairs_nodes(
      const Cix& cix,
      const imag_dlr::Workspace& workspace,
      const std::vector<double>& pseudo_degeneracy,
      int topology,
      int cut_edge,
      const std::array<int,3>& direction,
      const std::vector<std::array<int,2>>& bath_pairs,
      bool sum_free_root_direction=false,
      int first_pole_pair=0,
      int last_pole_pair=-1,
      const PseudoBlocks* shared_pseudo_blocks=nullptr,
      const std::vector<OperatorTranspose>* shared_operator_cache=nullptr)
  {
    BathPairNodeResult result{
      MatrixXd::Zero(workspace.rank(),workspace.pseudo_count()),0};
    if (bath_pairs.empty()) return result;

    // Validate with the first pair, then validate every remaining bath index
    // while constructing its compact execution plan below.
    const RootedCut cut=MakeRootedCut(topology,direction,cut_edge);
    const int root_line=cut.root.original_line;
    std::array<int,3> first_orbital{{-1,-1,-1}};
    first_orbital[cut.split[0].original_line]=bath_pairs[0][0];
    first_orbital[cut.split[1].original_line]=bath_pairs[0][1];
    Validate(cix,workspace,pseudo_degeneracy,topology,cut_edge,
             first_orbital,direction);

    const int root_vertex=cut.root.second_vertex;
    const int rank=workspace.rank();
    const int pole_count=workspace.pole_count();
    const int total_pole_pairs=pole_count*pole_count;
    if (last_pole_pair<0) last_pole_pair=total_pole_pairs;
    if (first_pole_pair<0 || first_pole_pair>last_pole_pair ||
        last_pole_pair>total_pole_pairs)
      throw std::out_of_range("summed-root TCA pole-pair range is invalid");
    PseudoBlocks owned_pseudo_blocks;
    if (shared_pseudo_blocks==nullptr)
      owned_pseudo_blocks=PreparePseudoBlocks(cix,workspace);
    const PseudoBlocks& pseudo_blocks=shared_pseudo_blocks!=nullptr
      ? *shared_pseudo_blocks : owned_pseudo_blocks;

    // Each physical one-body operator is built only once for this geometry.
    // A bath-pair plan below stores integer operator indices rather than
    // copying dense matrices four times per pair.
    std::vector<OperatorTranspose> owned_operator_cache;
    if (shared_operator_cache==nullptr) {
      owned_operator_cache.reserve(2*cix.Nbath);
      for (int op=0;op<2*cix.Nbath;++op)
        owned_operator_cache.push_back(MakeOperatorTranspose(cix,op));
    }
    const std::vector<OperatorTranspose>& operator_cache=
      shared_operator_cache!=nullptr
        ? *shared_operator_cache : owned_operator_cache;
    if (static_cast<int>(operator_cache.size())!=2*cix.Nbath)
      throw std::invalid_argument("TCA shared operator cache size mismatch");

    // Root operators and Delta metadata are independent of the two explicit
    // baths and are shared by the complete pole-major batch.
    std::vector<StoredRootPair> roots;
    const int root_direction_count=
      sum_free_root_direction && root_line!=edmft::tca::root_line ? 2 : 1;
    roots.reserve(root_direction_count*cix.Nbath);
    for (int root_choice=0;root_choice<root_direction_count;++root_choice) {
      int local_root_direction=cut.root.direction;
      if (root_direction_count==2) {
        auto root_variant_direction=direction;
        root_variant_direction[root_line]=root_choice;
        const RootedCut root_variant=MakeRootedCut(
          topology,root_variant_direction,cut_edge);
        if (root_variant.root.original_line!=root_line ||
            root_variant.root.first_vertex!=cut.root.first_vertex ||
            root_variant.root.second_vertex!=cut.root.second_vertex)
          throw std::logic_error(
            "TCA root direction unexpectedly changed cut geometry");
        local_root_direction=root_variant.root.direction;
      }
      for (int bath=0;bath<cix.Nbath;++bath) {
        const int ds=2*bath+local_root_direction;
        const int signed_channel=cix.bath_sigind0(bath);
        roots.push_back({operator_cache[Opposite(ds)],operator_cache[ds],
                         std::abs(signed_channel),local_root_direction,
                         static_cast<double>(
                           SignedIndexSign(signed_channel))});
      }
    }

    struct BathPlan {
      std::array<int,6> vertex_operator{{-1,-1,-1,-1,-1,-1}};
      std::array<int,2> channel{{0,0}};
      double fixed_prefactor=1.0;
      bool nonempty=false;
    };
    std::vector<BathPlan> plans(bath_pairs.size());
    for (std::size_t pair=0;pair<bath_pairs.size();++pair) {
      BathPlan& plan=plans[pair];
      plan.fixed_prefactor=
        edmft::tca::topology_crossing_sign[topology]*
        edmft::tca::topology_symmetry_factor[topology];
      for (int split=0;split<2;++split) {
        const int bath=bath_pairs[pair][split];
        if (bath<0 || bath>=cix.Nbath)
          throw std::out_of_range("summed-root TCA split bath is invalid");
        const auto& line=cut.split[split];
        const int ds=2*bath+line.direction;
        plan.vertex_operator[line.first_vertex]=Opposite(ds);
        plan.vertex_operator[line.second_vertex]=ds;
        const int signed_channel=cix.bath_sigind0(bath);
        plan.channel[split]=std::abs(signed_channel);
        plan.fixed_prefactor*=SignedIndexSign(signed_channel);
      }
    }

    // The production MPI path passes the complete Cartesian product of the
    // two split-line bath labels.  In that case there is no reason to carry a
    // bath label through the whole six-vertex backbone.  A label is born at
    // the first endpoint of its hybridization line and is dead immediately
    // after the second endpoint.  Since all later operations are linear, we
    // may sum a dead label before performing the next DLR convolution.
    //
    // This is the matrix-valued analogue of eliminating an internal index at
    // the earliest possible point in a tensor network.  It is especially
    // useful for disjoint lines (only Nbath live states); nested and crossing
    // lines need Nbath^2 states only between their overlapping endpoints.
    bool complete_cartesian=
      bath_pairs.size()==static_cast<std::size_t>(cix.Nbath*cix.Nbath);
    if (complete_cartesian)
      for (int pair=0;pair<cix.Nbath*cix.Nbath;++pair)
        if (bath_pairs[pair][0]!=pair/cix.Nbath ||
            bath_pairs[pair][1]!=pair%cix.Nbath) {
          complete_cartesian=false;
          break;
        }

    // Connectivity is independent of pole number and of all scalar kernel
    // factors.  Determine the old `nonempty_pairs` diagnostic once with block
    // keys alone.  This exactly mirrors MultiplyRight/ContractRoot without
    // allocating r matrix values for every bath pair.
    // The dynamic caller records connectivity only for the first pole tile of
    // a structure.  Later tiles must not repeat this pole-independent walk;
    // with fine SC tiles that needless work would otherwise be multiplied by
    // O(100) scheduler tasks per structure.
    if (complete_cartesian && first_pole_pair==0) {
      using KeySet=std::set<BlockKey>;
      for (const BathPlan& plan:plans) {
        KeySet keys;
        for (int block=0;block<cix.Nblock;++block)
          keys.insert({block,block});
        for (int vertex=1;vertex<6 && !keys.empty();++vertex) {
          KeySet next;
          if (vertex==root_vertex) {
            for (const auto& root:roots)
              for (const auto& key:keys) {
                if (key.first<0 || key.first>=
                      static_cast<int>(root.zero.by_target.size())) continue;
                if (key.second<0 || key.second>=
                      static_cast<int>(root.endpoint.by_source.size()))
                  continue;
                const int right_index=root.endpoint.by_source[key.second];
                if (right_index<0) continue;
                const auto& right=root.endpoint.entries[right_index];
                for (const int left_index:root.zero.by_target[key.first]) {
                  const auto& left=root.zero.entries[left_index];
                  next.insert({left.source,right.target});
                }
              }
          } else {
            const int op=plan.vertex_operator[vertex];
            if (op<0)
              throw std::logic_error(
                "non-root TCA vertex has no split-line operator");
            const auto& transpose=operator_cache[op];
            for (const auto& key:keys) {
              if (key.second<0 || key.second>=
                    static_cast<int>(transpose.by_source.size())) continue;
              const int index=transpose.by_source[key.second];
              if (index<0) continue;
              next.insert({key.first,transpose.entries[index].target});
            }
          }
          keys.swap(next);
        }
        if (!keys.empty()) ++result.nonempty_pairs;
      }
    }

    for (int pole_pair=first_pole_pair;
         pole_pair<last_pole_pair;++pole_pair) {
        const int l0=pole_pair/pole_count;
        const int l1=pole_pair%pole_count;
        const std::array<int,2> pole{{l0,l1}};
        std::array<VectorXd,6> vertex_factor;
        std::array<VectorXd,5> edge_factor;
        for (auto& x:vertex_factor) x=VectorXd::Ones(rank);
        for (auto& x:edge_factor) x=VectorXd::Ones(rank);
        std::array<double,2> pole_prefactor{{1.0,1.0}};

        for (int split=0;split<2;++split) {
          const auto& line=cut.split[split];
          const double omega=workspace.poles()(pole[split]);
          const LineDecoration decoration=MakeRootedLineDecoration(
            line.first_vertex,line.second_vertex,line.direction,omega);
          pole_prefactor[split]=CodePolePrefactor(
            1.0,KernelZero(workspace,pole[split],
                           decoration.denominator_kernel),decoration);
          for (int vertex=0;vertex<6;++vertex) {
            const int sign=decoration.vertex_kernel[vertex];
            if (sign!=0)
              vertex_factor[vertex].array()*=
                (sign>0 ? workspace.kernel_plus().col(pole[split])
                        : workspace.kernel_minus().col(pole[split])).array();
          }
          for (int edge=0;edge<5;++edge) {
            const int sign=decoration.edge_kernel[edge];
            if (sign!=0)
              edge_factor[edge].array()*=
                (sign>0 ? workspace.kernel_plus().col(pole[split])
                        : workspace.kernel_minus().col(pole[split])).array();
          }
        }

        // Only edges one through four enter a convolution.  Several often
        // carry exactly the same multiplier (especially the all-ones base
        // case), so deduplicate them before the relatively expensive DLR
        // solve.  The cache lives for one pole pair only, bounding memory by
        // a handful of r x Npseudo matrices even for the f-shell problem.
        std::vector<VectorXd> unique_scalars;
        std::vector<PreparedPseudoCoefficients> prepared;
        std::array<int,5> prepared_index{{-1,-1,-1,-1,-1}};
        for (int edge=1;edge<5;++edge) {
          int found=-1;
          for (int i=0;i<static_cast<int>(unique_scalars.size());++i)
            if ((unique_scalars[i].array()==edge_factor[edge].array()).all()) {
              found=i;
              break;
            }
          if (found<0) {
            found=static_cast<int>(unique_scalars.size());
            unique_scalars.push_back(edge_factor[edge]);
            prepared.push_back(
              PreparePseudoCoefficients(workspace,edge_factor[edge]));
          }
          prepared_index[edge]=found;
        }

        if (complete_cartesian) {
          // State index zero means neither split bath is live.  A live bath b
          // is encoded as b+1, so the complete state space is only
          // (Nbath+1)^2.  Most entries stay empty; at any vertex the number of
          // populated entries is 1, Nbath, or Nbath^2 according to the line
          // overlap pattern.
          const int state_extent=cix.Nbath+1;
          const int state_count=state_extent*state_extent;
          const auto encode_state=[&](int bath0,int bath1) {
            return (bath0+1)*state_extent+(bath1+1);
          };
          std::array<int,6> event_split{{-1,-1,-1,-1,-1,-1}};
          std::array<bool,6> event_open{{false,false,false,false,false,false}};
          for (int split=0;split<2;++split) {
            event_split[cut.split[split].first_vertex]=split;
            event_open[cut.split[split].first_vertex]=true;
            event_split[cut.split[split].second_vertex]=split;
          }

          std::array<std::vector<double>,2> line_weight;
          for (int split=0;split<2;++split) {
            line_weight[split].resize(cix.Nbath);
            for (int bath=0;bath<cix.Nbath;++bath) {
              const int signed_channel=cix.bath_sigind0(bath);
              line_weight[split][bath]=
                SignedIndexSign(signed_channel)*pole_prefactor[split]*
                workspace.residues()(std::abs(signed_channel),pole[split]);
            }
          }

          std::vector<SparseFunction> states(state_count);
          states[encode_state(-1,-1)]=
            PseudoValues(pseudo_blocks,edge_factor[0]);
          for (int vertex=1;vertex<6;++vertex) {
            std::vector<SparseFunction> next(state_count);
            for (int state=0;state<state_count;++state) {
              if (states[state].empty()) continue;
              std::array<int,2> live{{
                state/state_extent-1,state%state_extent-1}};
              if (vertex==root_vertex) {
                AddFunction(next[state],ContractRoot(
                  states[state],roots,workspace,vertex_factor[vertex]));
                continue;
              }

              const int split=event_split[vertex];
              if (split<0)
                throw std::logic_error(
                  "non-root TCA vertex has no split-line event");
              const auto& line=cut.split[split];
              if (event_open[vertex]) {
                if (live[split]>=0)
                  throw std::logic_error("TCA split bath opened twice");
                for (int bath=0;bath<cix.Nbath;++bath) {
                  // A zero residue makes the complete branch zero.  Dropping
                  // it at the opening vertex avoids every intervening matrix
                  // multiplication and convolution.
                  if (line_weight[split][bath]==0.0) continue;
                  auto target_live=live;
                  target_live[split]=bath;
                  const int ds=2*bath+line.direction;
                  AddFunction(next[encode_state(target_live[0],
                                                target_live[1])],
                    MultiplyRight(states[state],operator_cache[Opposite(ds)],
                                  vertex_factor[vertex]));
                }
              } else {
                const int bath=live[split];
                if (bath<0)
                  throw std::logic_error(
                    "TCA split bath closed before it was opened");
                auto target_live=live;
                target_live[split]=-1;
                const int ds=2*bath+line.direction;
                AddFunction(next[encode_state(target_live[0],
                                              target_live[1])],
                  MultiplyRight(states[state],operator_cache[ds],
                                vertex_factor[vertex],
                                line_weight[split][bath]));
              }
            }
            if (vertex<5)
              ConvolveRightPseudoPreparedBatch(
                next,cix,workspace,prepared[prepared_index[vertex]]);
            states.swap(next);
          }
          const double common_prefactor=
            edmft::tca::topology_crossing_sign[topology]*
            edmft::tca::topology_symmetry_factor[topology];
          AccumulateCut(states[encode_state(-1,-1)],cix,
                        common_prefactor,pseudo_degeneracy,result.nodes);
          continue;
        }

        // Reference/fallback path for an arbitrary subset or ordering of bath
        // pairs.  The static scheduler still uses this path; production's
        // dynamic pole-tile scheduler supplies the Cartesian domain above.
        for (std::size_t pair=0;pair<plans.size();++pair) {
          BathPlan& plan=plans[pair];
          double prefactor=plan.fixed_prefactor;
          for (int split=0;split<2;++split)
            prefactor*=pole_prefactor[split]*workspace.residues()(
              plan.channel[split],pole[split]);
          if (prefactor==0.0) continue;

          SparseFunction chain=PseudoValues(pseudo_blocks,edge_factor[0]);
          for (int vertex=1;vertex<6 && !chain.empty();++vertex) {
            if (vertex==root_vertex) {
              chain=ContractRoot(chain,roots,workspace,
                                 vertex_factor[vertex]);
            } else {
              const int op=plan.vertex_operator[vertex];
              if (op<0)
                throw std::logic_error(
                  "non-root TCA vertex has no split-line operator");
              chain=MultiplyRight(chain,operator_cache[op],
                                  vertex_factor[vertex]);
            }
            if (vertex<5 && !chain.empty())
              chain=ConvolveRightPseudoPrepared(
                chain,cix,workspace,prepared[prepared_index[vertex]]);
          }
          if (!chain.empty()) plan.nonempty=true;
          AccumulateCut(chain,cix,prefactor,pseudo_degeneracy,result.nodes);
        }
      }

    if (!complete_cartesian)
      for (const auto& plan:plans)
        if (plan.nonempty) ++result.nonempty_pairs;
    return result;
  }

public:
  using BlockKey=std::pair<int,int>;
  using SparseFunction=std::map<BlockKey,MatrixFunction>;

  static int Opposite(int op) noexcept
  {
    return 1-(op&1)+2*(op>>1);
  }

  static OperatorTranspose MakeOperatorTranspose(const Cix& cix,int op)
  {
    OperatorTranspose result;
    result.by_source.assign(cix.Nblock,-1);
    result.by_target.resize(cix.Nblock);
    for (int source=0;source<cix.Nblock;++source) {
      const int target=cix.Fbc(op,source);
      if (target<0) continue;
      MatrixXd matrix=MatrixXd::Zero(cix.Block[source].size(),
                                     cix.Block[target].size());
      std::vector<int> target_column(cix.Block[source].size(),-1);
      std::vector<double> row_sign(cix.Block[source].size(),0.0);
      for (int i=0;i<cix.Block[source].size();++i) {
        const int state=cix.Block[source](i);
        const int mapped=cix.Fi(op,state);
        if (mapped==0) continue;
        const int j=cix.inBlock(std::abs(mapped));
        const double sign=SignedIndexSign(mapped);
        matrix(i,j)=sign;
        target_column[i]=j;
        row_sign[i]=sign;
      }
      const int index=static_cast<int>(result.entries.size());
      OperatorEntry entry;
      entry.source=source;
      entry.target=target;
      entry.matrix=std::move(matrix);
      entry.signed_partial_permutation=true;
      entry.target_column=std::move(target_column);
      entry.row_sign=std::move(row_sign);
      for (int i=0;i<entry.matrix.rows();++i)
        for (int j=0;j<entry.matrix.cols();++j)
          if (entry.matrix(i,j)!=0.0)
            entry.nonzero.push_back({i,j,entry.matrix(i,j)});
      result.entries.push_back(std::move(entry));
      result.by_source[source]=index;
      result.by_target[target].push_back(index);
    }
    return result;
  }

  static PseudoBlocks PreparePseudoBlocks(
      const Cix& cix,const imag_dlr::Workspace& workspace)
  {
    PseudoBlocks result(cix.Nblock);
    for (int block=0;block<cix.Nblock;++block) {
      const auto& labels=cix.Sgm_symmetry[block];
      const int n=cix.Block[block].size();
      result[block].assign(workspace.rank(),MatrixXd::Zero(n,n));
      for (int i=0;i<n;++i)
        for (int j=0;j<n;++j) {
          const int label=labels(i,j);
          if (label==0) continue;
          const int pseudo=std::abs(label)-1;
          const double sign=SignedIndexSign(label);
          for (int node=0;node<workspace.rank();++node)
            result[block][node](i,j)=
              sign*workspace.pseudo_values()(node,pseudo);
        }
    }
    return result;
  }

  static SparseFunction PseudoValues(
      const PseudoBlocks& pseudo_blocks,const VectorXd& scalar)
  {
    SparseFunction result;
    for (int block=0;block<static_cast<int>(pseudo_blocks.size());++block) {
      MatrixFunction values=pseudo_blocks[block];
      bool nonzero=false;
      for (int node=0;node<static_cast<int>(values.size());++node)
        if (scalar(node)!=0.0) {
          values[node]*=scalar(node);
          nonzero=nonzero || values[node].squaredNorm()!=0.0;
        } else {
          values[node].setZero();
        }
      // In an active-space calculation an entire atomic symmetry block can
      // contain only inactive pseudo labels.  Do not carry its exact zero
      // matrices through five operator/convolution stages.
      if (nonzero)
        result.emplace(BlockKey{block,block},std::move(values));
    }
    return result;
  }

  // Scatter coefficients stored once per unique pseudo-particle label into
  // one atomic symmetry block.  Sgm_symmetry contains signed, one-based
  // labels; its sign is part of the matrix element.  Keeping the DLR solve in
  // the unique-label representation is substantially cheaper than solving
  // independently for every dense entry of every block (especially for the
  // superconducting cluster and large f-shell calculations).
  static MatrixFunction PseudoCoefficientBlock(
      const Cix& cix,int block,const MatrixXd& unique_coefficients)
  {
    const auto& labels=cix.Sgm_symmetry[block];
    const int n=cix.Block[block].size();
    MatrixFunction result(
      unique_coefficients.rows(),MatrixXd::Zero(n,n));
    for (int i=0;i<n;++i)
      for (int j=0;j<n;++j) {
        const int label=labels(i,j);
        if (label==0) continue;
        const int pseudo=std::abs(label)-1;
        const double sign=SignedIndexSign(label);
        for (int coefficient=0;
             coefficient<unique_coefficients.rows();++coefficient)
          result[coefficient](i,j)=
            sign*unique_coefficients(coefficient,pseudo);
      }
    return result;
  }

  static void Add(SparseFunction& destination,const BlockKey& key,
                  MatrixFunction values)
  {
    auto [it,inserted]=destination.try_emplace(key,std::move(values));
    if (inserted) return;
    if (it->second.size()!=values.size())
      throw std::logic_error("TCA sparse-function rank mismatch");
    for (int node=0;node<static_cast<int>(values.size());++node)
      it->second[node]+=values[node];
  }

  // Merge a complete sparse block function.  Taking the source by value lets
  // callers pass a temporary and moves newly encountered blocks without an
  // additional copy.  This is the reduction primitive used when a split-line
  // bath reaches its second endpoint and can be summed immediately.
  static void AddFunction(SparseFunction& destination,SparseFunction source)
  {
    for (auto& [key,values]:source)
      Add(destination,key,std::move(values));
  }

  static MatrixXd MultiplyRightOperator(
      const MatrixXd& values,const OperatorEntry& op)
  {
    if (!op.signed_partial_permutation)
      return (values*op.matrix).eval();
    if (values.cols()!=static_cast<int>(op.target_column.size()))
      throw std::logic_error("TCA right signed-map dimension mismatch");
    MatrixXd result=MatrixXd::Zero(values.rows(),op.matrix.cols());
    for (int source=0;source<static_cast<int>(op.target_column.size());
         ++source) {
      const int target=op.target_column[source];
      if (target<0) continue;
      result.col(target).noalias()+=op.row_sign[source]*values.col(source);
    }
    return result;
  }

  static MatrixXd MultiplyLeftOperator(
      const OperatorEntry& op,const MatrixXd& values)
  {
    if (!op.signed_partial_permutation)
      return (op.matrix*values).eval();
    if (values.rows()!=op.matrix.cols())
      throw std::logic_error("TCA left signed-map dimension mismatch");
    MatrixXd result=MatrixXd::Zero(op.matrix.rows(),values.cols());
    for (int source=0;source<static_cast<int>(op.target_column.size());
         ++source) {
      const int target=op.target_column[source];
      if (target<0) continue;
      result.row(source)=op.row_sign[source]*values.row(target);
    }
    return result;
  }

  static MatrixXd MultiplyBothOperators(
      const OperatorEntry& left,const MatrixXd& values,
      const OperatorEntry& right)
  {
    if (!left.signed_partial_permutation ||
        !right.signed_partial_permutation)
      return (left.matrix*values*right.matrix).eval();
    if (values.rows()!=left.matrix.cols() ||
        values.cols()!=right.matrix.rows())
      throw std::logic_error("TCA two-sided signed-map dimension mismatch");
    MatrixXd result=MatrixXd::Zero(left.matrix.rows(),right.matrix.cols());
    // Directly gather both indices.  Besides replacing two dense GEMMs, this
    // avoids allocating the intermediate values*right matrix in the root
    // contraction, which is repeated for every root bath and DLR node.
    for (int left_source=0;
         left_source<static_cast<int>(left.target_column.size());
         ++left_source) {
      const int left_target=left.target_column[left_source];
      if (left_target<0) continue;
      for (int right_source=0;
           right_source<static_cast<int>(right.target_column.size());
           ++right_source) {
        const int right_target=right.target_column[right_source];
        if (right_target<0) continue;
        result(left_source,right_target)+=
          left.row_sign[left_source]*right.row_sign[right_source]*
          values(left_target,right_source);
      }
    }
    return result;
  }

  // Add scalar*left*values*right directly to an existing destination block.
  //
  // ContractRoot used to form this product in a fresh MatrixXd for every
  // root bath and DLR node, store those matrices in a temporary
  // MatrixFunction, and then traverse the complete function once more to add
  // it to the sparse result.  Ordinary impurity operators are signed partial
  // permutations, so every output element can instead be scattered directly
  // into its final accumulator.  This removes both the per-contribution
  // allocation/zeroing and the subsequent full-matrix addition.
  static void AccumulateBothOperators(
      MatrixXd& destination,const OperatorEntry& left,
      const MatrixXd& values,const OperatorEntry& right,double scalar)
  {
    if (values.rows()!=left.matrix.cols() ||
        values.cols()!=right.matrix.rows() ||
        destination.rows()!=left.matrix.rows() ||
        destination.cols()!=right.matrix.cols())
      throw std::logic_error(
        "TCA fused two-sided accumulation dimension mismatch");

    // Keep a general-matrix fallback so this primitive remains correct if it
    // is reused with Bulla or other dense connectors.  Every root operator
    // currently takes the signed-map branch below.
    if (!left.signed_partial_permutation ||
        !right.signed_partial_permutation) {
      destination.noalias()+=scalar*left.matrix*values*right.matrix;
      return;
    }

    for (int left_source=0;
         left_source<static_cast<int>(left.target_column.size());
         ++left_source) {
      const int left_target=left.target_column[left_source];
      if (left_target<0) continue;
      const double left_scale=scalar*left.row_sign[left_source];
      for (int right_source=0;
           right_source<static_cast<int>(right.target_column.size());
           ++right_source) {
        const int right_target=right.target_column[right_source];
        if (right_target<0) continue;
        destination(left_source,right_target)+=
          left_scale*right.row_sign[right_source]*
          values(left_target,right_source);
      }
    }
  }

  static double TraceLeftOperator(
      const OperatorEntry& op,const MatrixXd& values)
  {
    if (!op.signed_partial_permutation)
      return (op.matrix.array()*values.transpose().array()).sum();
    if (values.rows()!=op.matrix.cols() || values.cols()!=op.matrix.rows())
      throw std::logic_error("TCA signed-map trace dimension mismatch");
    double result=0.0;
    for (int source=0;source<static_cast<int>(op.target_column.size());
         ++source) {
      const int target=op.target_column[source];
      if (target>=0) result+=op.row_sign[source]*values(target,source);
    }
    return result;
  }

  // Evaluate Tr[op * left * right] without materializing left*right.  The
  // matrix element selected by op(i,j) is the dot product of row j of left
  // and column i of right.  This costs O(nnz(op)*N), instead of O(N^3), and
  // applies both to ordinary signed-map operators and sparse Bulla matrices.
  static double TraceLeftProduct(
      const OperatorEntry& op,const MatrixXd& left,const MatrixXd& right)
  {
    if (left.rows()!=op.matrix.cols() ||
        right.cols()!=op.matrix.rows() || left.cols()!=right.rows())
      throw std::logic_error("TCA sparse product-trace dimension mismatch");
    double result=0.0;
    for (const auto& element:op.nonzero)
      result+=element.value*
        left.row(element.column).dot(right.col(element.row));
    return result;
  }

  static SparseFunction MultiplyRight(
      const SparseFunction& input,const OperatorTranspose& op,
      const VectorXd& scalar,double prefactor=1.0)
  {
    SparseFunction result;
    for (const auto& [key,values]:input) {
      const int middle=key.second;
      if (middle<0 || middle>=static_cast<int>(op.by_source.size())) continue;
      const int index=op.by_source[middle];
      if (index<0) continue;
      const auto& entry=op.entries[index];
      MatrixFunction product(values.size());
      for (int node=0;node<static_cast<int>(values.size());++node)
        product[node]=MultiplyRightOperator(values[node],entry)*
          (prefactor*scalar(node));
      Add(result,{key.first,entry.target},std::move(product));
    }
    return result;
  }

  static SparseFunction ContractRoot(
      const SparseFunction& input,const std::vector<StoredRootPair>& roots,
      const imag_dlr::Workspace& workspace,const VectorXd& vertex_scalar)
  {
    SparseFunction result;
    for (const auto& root:roots) {
      const auto& delta_matrix=root.direction==0
        ? workspace.delta_plus() : workspace.delta_minus();
      if (root.channel>=delta_matrix.cols())
        throw std::out_of_range("TCA root channel is outside Delta");
      // This scalar is common to every compatible atomic-block path for the
      // current root bath.  Form it once rather than copying Delta and
      // multiplying the two node factors inside every block contribution.
      const VectorXd root_scalar=
        (-root.bath_sign*delta_matrix.col(root.channel).array()*
         vertex_scalar.array()).matrix();
      for (const auto& [key,values]:input) {
        if (key.first<0 || key.first>=
              static_cast<int>(root.zero.by_target.size())) continue;
        if (key.second<0 || key.second>=
              static_cast<int>(root.endpoint.by_source.size())) continue;
        const int right_index=root.endpoint.by_source[key.second];
        if (right_index<0) continue;
        const auto& right=root.endpoint.entries[right_index];
        for (const int left_index:root.zero.by_target[key.first]) {
          const auto& left=root.zero.entries[left_index];
          const BlockKey output_key{left.source,right.target};
          auto [output,inserted]=result.try_emplace(output_key);
          if (inserted) {
            output->second.resize(values.size());
            for (auto& node_values:output->second)
              node_values=MatrixXd::Zero(
                left.matrix.rows(),right.matrix.cols());
          } else {
            if (output->second.size()!=values.size())
              throw std::logic_error("TCA root contraction rank mismatch");
            // All nodes belonging to one sparse block are constructed with
            // identical dimensions; checking the first one is sufficient.
            if (!output->second.empty() &&
                (output->second.front().rows()!=left.matrix.rows() ||
                 output->second.front().cols()!=right.matrix.cols()))
              throw std::logic_error(
                "TCA root contraction block dimension mismatch");
          }
          for (int node=0;node<static_cast<int>(values.size());++node) {
            AccumulateBothOperators(
              output->second[node],left,values[node],right,
              root_scalar(node));
          }
        }
      }
    }
    return result;
  }

  static SparseFunction ConvolveRightPseudo(
      const SparseFunction& input,const Cix& cix,
      const imag_dlr::Workspace& workspace,const VectorXd& scalar)
  {
    if (scalar.size()!=workspace.rank())
      throw std::invalid_argument(
        "TCA pseudo convolution scalar has wrong DLR rank");

    // Dressing by scalar is diagonal in the representative-time values.
    // Perform the values->coefficients solve once for all unique pseudo
    // functions, rather than once for every matrix entry in every block.
    // If there is no dressing, the coefficients already held by Workspace
    // can be used directly without either a solve or a copy.
    const PreparedPseudoCoefficients prepared=
      PreparePseudoCoefficients(workspace,scalar);
    return ConvolveRightPseudoPrepared(input,cix,workspace,prepared);
  }

  static PreparedPseudoCoefficients PreparePseudoCoefficients(
      const imag_dlr::Workspace& workspace,const VectorXd& scalar)
  {
    if (scalar.size()!=workspace.rank())
      throw std::invalid_argument(
        "TCA pseudo convolution scalar has wrong DLR rank");
    PreparedPseudoCoefficients result;
    if ((scalar.array()==1.0).all()) {
      result.hilbert=workspace.pseudo_hilbert_coefficients();
      return result;
    }
    result.use_workspace_coefficients=false;
    MatrixXd dressed_values=workspace.pseudo_values();
    dressed_values.array().colwise()*=scalar.array();
    result.dressed=
      workspace.basis().values_to_coefficients(dressed_values);
    result.hilbert=workspace.basis().convolution_hilbert(result.dressed);
    return result;
  }

  static SparseFunction ConvolveRightPseudoPrepared(
      const SparseFunction& input,const Cix& cix,
      const imag_dlr::Workspace& workspace,
      const PreparedPseudoCoefficients& prepared,
      bool reflect_output=false)
  {
    const MatrixXd& unique_coefficients=
      prepared.use_workspace_coefficients
        ? workspace.pseudo_coefficients() : prepared.dressed;

    SparseFunction result;
    if (prepared.block_ready.empty()) {
      prepared.coefficient_blocks.resize(cix.Nblock);
      prepared.hilbert_blocks.resize(cix.Nblock);
      prepared.block_ready.assign(cix.Nblock,0);
      prepared.block_nonzero.assign(cix.Nblock,0);
    } else if (static_cast<int>(prepared.block_ready.size())!=cix.Nblock) {
      throw std::logic_error("TCA prepared pseudo block-cache mismatch");
    }
    std::vector<BlockKey> keys;
    std::vector<const MatrixFunction*> value_batch;
    std::vector<const MatrixFunction*> coefficient_batch;
    std::vector<const MatrixFunction*> hilbert_batch;
    keys.reserve(input.size());
    value_batch.reserve(input.size());
    coefficient_batch.reserve(input.size());
    hilbert_batch.reserve(input.size());
    for (const auto& [key,values]:input) {
      const int block=key.second;
      if (!prepared.block_ready[block]) {
        prepared.coefficient_blocks[block]=
          PseudoCoefficientBlock(cix,block,unique_coefficients);
        prepared.hilbert_blocks[block]=
          PseudoCoefficientBlock(cix,block,prepared.hilbert);
        for (const auto& matrix:prepared.coefficient_blocks[block])
          if (matrix.squaredNorm()!=0.0) {
            prepared.block_nonzero[block]=1;
            break;
          }
        prepared.block_ready[block]=1;
      }
      if (!prepared.block_nonzero[block]) continue;
      keys.push_back(key);
      value_batch.push_back(&values);
      coefficient_batch.push_back(&prepared.coefficient_blocks[block]);
      hilbert_batch.push_back(&prepared.hilbert_blocks[block]);
    }
    auto convolved=
      workspace.basis().time_ordered_convolution_values_right_prepared_batch(
        value_batch,coefficient_batch,hilbert_batch,reflect_output);
    for (int item=0;item<static_cast<int>(keys.size());++item)
      Add(result,keys[item],std::move(convolved[item]));
    return result;
  }

  // Convolve every populated member of a bath-state array with the same
  // dressed pseudo propagator.  The scalar matrix entries of many symmetry
  // blocks and many live bath states are collected into moderately sized
  // multiple-RHS batches.  Thus each chunk pays for one DLR interpolation
  // solve and two basis GEMMs, while the cap keeps memory bounded for large
  // f-shell Hilbert spaces.
  static void ConvolveRightPseudoPreparedBatch(
      std::vector<SparseFunction>& functions,const Cix& cix,
      const imag_dlr::Workspace& workspace,
      const PreparedPseudoCoefficients& prepared,
      bool reflect_output=false)
  {
    const MatrixXd& unique_coefficients=
      prepared.use_workspace_coefficients
        ? workspace.pseudo_coefficients() : prepared.dressed;
    if (prepared.block_ready.empty()) {
      prepared.coefficient_blocks.resize(cix.Nblock);
      prepared.hilbert_blocks.resize(cix.Nblock);
      prepared.block_ready.assign(cix.Nblock,0);
      prepared.block_nonzero.assign(cix.Nblock,0);
    } else if (static_cast<int>(prepared.block_ready.size())!=cix.Nblock) {
      throw std::logic_error("TCA prepared pseudo block-cache mismatch");
    }

    struct Item {
      int state=0;
      BlockKey key;
      const MatrixFunction* values=nullptr;
      const MatrixFunction* coefficients=nullptr;
      const MatrixFunction* hilbert=nullptr;
      std::size_t flat_size=0;
    };
    std::vector<Item> items;
    for (int state=0;state<static_cast<int>(functions.size());++state)
      for (const auto& [key,values]:functions[state]) {
        const int block=key.second;
        if (!prepared.block_ready[block]) {
          prepared.coefficient_blocks[block]=
            PseudoCoefficientBlock(cix,block,unique_coefficients);
          prepared.hilbert_blocks[block]=
            PseudoCoefficientBlock(cix,block,prepared.hilbert);
          for (const auto& matrix:prepared.coefficient_blocks[block])
            if (matrix.squaredNorm()!=0.0) {
              prepared.block_nonzero[block]=1;
              break;
            }
          prepared.block_ready[block]=1;
        }
        if (!prepared.block_nonzero[block]) continue;
        const auto& right=prepared.coefficient_blocks[block];
        const std::size_t input_size=
          static_cast<std::size_t>(values[0].rows())*values[0].cols();
        const std::size_t output_size=
          static_cast<std::size_t>(values[0].rows())*right[0].cols();
        items.push_back({state,key,&values,&right,
                         &prepared.hilbert_blocks[block],
                         std::max(input_size,output_size)});
      }
    if (items.empty()) return;

    std::vector<SparseFunction> output(functions.size());
    // Six or seven rank-by-M temporaries coexist inside the DLR kernel.  A
    // cap of 32768 scalar matrix entries is about 50 MiB at rank 32 and also
    // leaves enough columns for efficient level-3 BLAS on small blocks.
    constexpr std::size_t max_flat_entries=32768;
    for (std::size_t first=0;first<items.size();) {
      std::size_t last=first,flat_entries=0;
      while (last<items.size() &&
             (last==first ||
              flat_entries+items[last].flat_size<=max_flat_entries)) {
        flat_entries+=items[last].flat_size;
        ++last;
      }
      std::vector<const MatrixFunction*> value_batch;
      std::vector<const MatrixFunction*> coefficient_batch;
      std::vector<const MatrixFunction*> hilbert_batch;
      value_batch.reserve(last-first);
      coefficient_batch.reserve(last-first);
      hilbert_batch.reserve(last-first);
      for (std::size_t item=first;item<last;++item) {
        value_batch.push_back(items[item].values);
        coefficient_batch.push_back(items[item].coefficients);
        hilbert_batch.push_back(items[item].hilbert);
      }
      auto convolved=
        workspace.basis().time_ordered_convolution_values_right_prepared_batch(
          value_batch,coefficient_batch,hilbert_batch,reflect_output);
      for (std::size_t item=first;item<last;++item)
        Add(output[items[item].state],items[item].key,
            std::move(convolved[item-first]));
      first=last;
    }
    functions.swap(output);
  }

  static void AccumulateCut(
      const SparseFunction& values,const Cix& cix,double prefactor,
      const std::vector<double>& pseudo_degeneracy,MatrixXd& output)
  {
    for (const auto& [key,function]:values) {
      if (key.first!=key.second) continue;
      const int block=key.first;
      const auto& labels=cix.Sgm_symmetry[block];
      const int n=cix.Block[block].size();
      for (int incoming=0;incoming<n;++incoming)
        for (int outgoing=0;outgoing<n;++outgoing) {
          const int label=labels(outgoing,incoming);
          if (label==0) continue;
          const int pseudo=std::abs(label)-1;
          const double weight=prefactor*SignedIndexSign(label)/
                              pseudo_degeneracy[pseudo];
          for (int node=0;node<static_cast<int>(function.size());++node)
            output(node,pseudo)+=
              weight*function[node](incoming,outgoing);
        }
    }
  }

  static void Validate(
      const Cix& cix,const imag_dlr::Workspace& workspace,
      const std::vector<double>& pseudo_degeneracy,int topology,int cut_edge,
      const std::array<int,3>& orbital,
      const std::array<int,3>& direction)
  {
    if (!workspace.has_pseudos())
      throw std::logic_error("summed-root TCA requires pseudo coefficients");
    if (workspace.pseudo_count()!=cix.Npseudo_unique ||
        workspace.channel_count()!=cix.Nbath_unique_unique)
      throw std::invalid_argument("summed-root TCA workspace mismatch");
    if (static_cast<int>(pseudo_degeneracy.size())!=cix.Npseudo_unique)
      throw std::invalid_argument("summed-root TCA degeneracy mismatch");
    if (topology<0 || topology>=2 || cut_edge<0 || cut_edge>=6)
      throw std::out_of_range("summed-root TCA topology/cut is invalid");
    const RootedCut cut=MakeRootedCut(topology,direction,cut_edge);
    for (const auto& line:cut.split)
      if (orbital[line.original_line]<0 ||
          orbital[line.original_line]>=cix.Nbath)
        throw std::out_of_range("summed-root TCA split bath is invalid");
  }
};

// Physical counterpart of SummedRootSelfEnergyEvaluator.  The cut
// hybridization is an external measurement line, so it contributes no Delta.
// Its creation and annihilation operators are nevertheless a correlated pair.
// We leave both boundary operators open while evaluating the T and beta-T
// backbones, then contract every external bath locally and accumulate into
// its physical channel.  The same contraction is repeated with N_jm at the
// creation vertex for the Bulla estimators.
class SummedPhysicalEvaluator {
public:
  using Engine=SummedRootSelfEnergyEvaluator;
  using SparseFunction=Engine::SparseFunction;
  using OperatorTranspose=Engine::OperatorTranspose;
  using PseudoBlocks=Engine::PseudoBlocks;

  struct Result {
    MatrixXd g;
    MatrixXd f00;
    MatrixXd f01;
  };

  struct BathPairResult {
    Result values;
    std::uint64_t nonempty_pairs=0;
  };

  // orbital[cut_line] is ignored.  The other entries are the two explicitly
  // enumerated hybridization lines.
  static Result evaluate_line_nodes(
      const Cix& cix,const imag_dlr::Workspace& workspace,
      int topology,int cut_line,const std::array<int,3>& orbital,
      const std::array<int,3>& direction,int selected_external_bath=-1)
  {
    if (!workspace.has_pseudos())
      throw std::logic_error("summed TCA physical evaluation needs pseudos");
    const PhysicalCut cut=MakePhysicalCut(topology,direction,cut_line);
    const int rank=workspace.rank();
    Result result{MatrixXd::Zero(rank,workspace.channel_count()),
                  MatrixXd::Zero(rank,workspace.channel_count()),
                  MatrixXd::Zero(rank,workspace.channel_count())};
    const PseudoBlocks pseudo_blocks=Engine::PreparePseudoBlocks(cix,workspace);

    std::array<OperatorTranspose,6> vertex_operator;
    std::array<bool,6> has_operator{{false,false,false,false,false,false}};
    for (const auto& line:cut.split) {
      const int bath=orbital[line.original_line];
      if (bath<0 || bath>=cix.Nbath)
        throw std::out_of_range("summed TCA physical split bath is invalid");
      const int ds=2*bath+line.direction;
      vertex_operator[line.first_vertex]=
        Engine::MakeOperatorTranspose(cix,Engine::Opposite(ds));
      vertex_operator[line.second_vertex]=
        Engine::MakeOperatorTranspose(cix,ds);
      has_operator[line.first_vertex]=true;
      has_operator[line.second_vertex]=true;
    }

    struct ExternalPair {
      OperatorTranspose creation;
      OperatorTranspose annihilation;
      OperatorTranspose bulla;
      int channel=0;
      double bath_sign=1.0;
      bool anomalous=false;
    };
    std::vector<ExternalPair> external;
    external.reserve(cix.Nbath);
    if (selected_external_bath<-1 || selected_external_bath>=cix.Nbath)
      throw std::out_of_range("selected TCA external bath is invalid");
    const int first_external=
      selected_external_bath<0 ? 0 : selected_external_bath;
    const int last_external=
      selected_external_bath<0 ? cix.Nbath : selected_external_bath+1;
    for (int bath=first_external;bath<last_external;++bath) {
      const int signed_channel=cix.bath_sigind0(bath);
      external.push_back({
        Engine::MakeOperatorTranspose(cix,2*bath),
        Engine::MakeOperatorTranspose(cix,2*bath+1),
        MakeBullaTranspose(cix,bath),
        std::abs(signed_channel),
        static_cast<double>(SignedIndexSign(signed_channel)),
        IsF01Channel(cix,bath)});
    }

    double fixed_bath_sign=1.0;
    for (const auto& line:cut.split)
      fixed_bath_sign*=SignedIndexSign(
        cix.bath_sigind0(orbital[line.original_line]));

    for (int l0=0;l0<workspace.pole_count();++l0)
      for (int l1=0;l1<workspace.pole_count();++l1) {
        const std::array<int,2> pole{{l0,l1}};
        std::array<VectorXd,6> edge_factor;
        for (auto& x:edge_factor) x=VectorXd::Ones(rank);
        double prefactor=
          fixed_bath_sign*edmft::tca::topology_crossing_sign[topology]*
          edmft::tca::topology_symmetry_factor[topology];

        for (int split=0;split<2;++split) {
          const auto& line=cut.split[split];
          const double omega=workspace.poles()(pole[split]);
          const int signed_channel=
            cix.bath_sigind0(orbital[line.original_line]);
          const bool direct=
            (line.direction==0 && omega>=0.0) ||
            (line.direction==1 && omega<0.0);
          const int frequency_sign=direct
            ? (line.direction==0 ? +1 : -1)
            : (line.direction==0 ? -1 : +1);
          int edge_count=0;
          for (int edge=0;edge<6;++edge) {
            const bool selected=direct
              ? (line.first_vertex<=edge && edge<line.second_vertex)
              : (edge<line.first_vertex || edge>=line.second_vertex);
            if (!selected) continue;
            ++edge_count;
            edge_factor[edge].array()*=
              (frequency_sign>0
                 ? workspace.kernel_plus().col(pole[split])
                 : workspace.kernel_minus().col(pole[split])).array();
          }
          prefactor*=-workspace.residues()(
            std::abs(signed_channel),pole[split])/
            std::pow(KernelZero(workspace,pole[split],frequency_sign),
                     edge_count-1);
        }

        // Backbone from the creation boundary at local zero to the open
        // annihilation boundary q.  Neither external operator is applied.
        SparseFunction first=Engine::PseudoValues(
          pseudo_blocks,edge_factor[0]);
        for (int vertex=1;vertex<cut.external_vertex && !first.empty();
             ++vertex) {
          if (!has_operator[vertex])
            throw std::logic_error("missing TCA physical split operator");
          first=Engine::MultiplyRight(
            first,vertex_operator[vertex],VectorXd::Ones(rank));
          first=Engine::ConvolveRightPseudo(
            first,cix,workspace,edge_factor[vertex]);
        }

        // Independent beta-T backbone.  It starts immediately after the
        // annihilation boundary and includes the closing edge five.
        SparseFunction second=Engine::PseudoValues(
          pseudo_blocks,edge_factor[cut.external_vertex]);
        for (int vertex=cut.external_vertex+1;
             vertex<6 && !second.empty();++vertex) {
          if (!has_operator[vertex])
            throw std::logic_error("missing TCA physical split operator");
          second=Engine::MultiplyRight(
            second,vertex_operator[vertex],VectorXd::Ones(rank));
          second=Engine::ConvolveRightPseudo(
            second,cix,workspace,edge_factor[vertex]);
        }
        Reflect(second,workspace.basis());

        for (const auto& pair:external) {
          // Both estimators differ only at the creation boundary.  Closing
          // the ordinary and Bulla chains separately would form
          //
          //   L G_left R G_right
          //
          // twice and then perform another dense product merely to take its
          // trace.  Associate the product instead as
          //
          //   X = G_left R G_right,
          //   G = Tr[L_ordinary X],  F = Tr[L_Bulla X].
          //
          // Rather than forming X, the final implementation evaluates each
          // trace directly from the nonzero entries of its left operator.
          // This is exactly the delayed creation-operator contraction
          // described for physical diagrams, but it needs neither a dense X
          // nor a matrix-valued accumulator over pole or bath pairs.
          AccumulateBoundaryEstimators(
            first,second,pair.creation,pair.annihilation,pair.bulla,
            prefactor*pair.bath_sign,pair.channel,pair.anomalous,result);
        }
      }
    return result;
  }

  // Pole-major physical evaluator.  As for the self-energy counterpart, the
  // two entries of every bath pair follow cut.split[] order.  Kernel-dressed
  // pseudo coefficients are prepared once per pole pair and shared across
  // all locally owned bath pairs.  External G/F closure still sums its bath
  // inside the atomic contraction and is therefore also reused.
  static BathPairResult evaluate_line_bath_pairs_nodes(
      const Cix& cix,const imag_dlr::Workspace& workspace,
      int topology,int cut_line,const std::array<int,3>& direction,
      const std::vector<std::array<int,2>>& bath_pairs,
      int selected_external_bath=-1,
      int first_pole_pair=0,
      int last_pole_pair=-1)
  {
    BathPairResult batch{
      {MatrixXd::Zero(workspace.rank(),workspace.channel_count()),
       MatrixXd::Zero(workspace.rank(),workspace.channel_count()),
       MatrixXd::Zero(workspace.rank(),workspace.channel_count())},0};
    if (bath_pairs.empty()) return batch;
    if (!workspace.has_pseudos())
      throw std::logic_error("summed TCA physical evaluation needs pseudos");
    const PhysicalCut cut=MakePhysicalCut(topology,direction,cut_line);
    const int rank=workspace.rank();
    const int pole_count=workspace.pole_count();
    const int total_pole_pairs=pole_count*pole_count;
    if (last_pole_pair<0) last_pole_pair=total_pole_pairs;
    if (first_pole_pair<0 || first_pole_pair>last_pole_pair ||
        last_pole_pair>total_pole_pairs)
      throw std::out_of_range("summed physical TCA pole-pair range is invalid");
    const PseudoBlocks pseudo_blocks=Engine::PreparePseudoBlocks(cix,workspace);
    const VectorXd ones=VectorXd::Ones(rank);

    std::vector<OperatorTranspose> operator_cache;
    operator_cache.reserve(2*cix.Nbath);
    for (int op=0;op<2*cix.Nbath;++op)
      operator_cache.push_back(Engine::MakeOperatorTranspose(cix,op));

    struct BathPlan {
      std::array<int,6> vertex_operator{{-1,-1,-1,-1,-1,-1}};
      std::array<int,2> channel{{0,0}};
      double fixed_prefactor=1.0;
      bool nonempty=false;
    };
    std::vector<BathPlan> plans(bath_pairs.size());
    for (std::size_t pair=0;pair<bath_pairs.size();++pair) {
      BathPlan& plan=plans[pair];
      plan.fixed_prefactor=
        edmft::tca::topology_crossing_sign[topology]*
        edmft::tca::topology_symmetry_factor[topology];
      for (int split=0;split<2;++split) {
        const int bath=bath_pairs[pair][split];
        if (bath<0 || bath>=cix.Nbath)
          throw std::out_of_range("summed TCA physical split bath is invalid");
        const auto& line=cut.split[split];
        const int ds=2*bath+line.direction;
        plan.vertex_operator[line.first_vertex]=Engine::Opposite(ds);
        plan.vertex_operator[line.second_vertex]=ds;
        const int signed_channel=cix.bath_sigind0(bath);
        plan.channel[split]=std::abs(signed_channel);
        plan.fixed_prefactor*=SignedIndexSign(signed_channel);
      }
    }

    struct ExternalPair {
      OperatorTranspose creation;
      OperatorTranspose annihilation;
      OperatorTranspose bulla;
      int channel=0;
      double bath_sign=1.0;
      bool anomalous=false;
    };
    std::vector<ExternalPair> external;
    external.reserve(cix.Nbath);
    if (selected_external_bath<-1 || selected_external_bath>=cix.Nbath)
      throw std::out_of_range("selected TCA external bath is invalid");
    const int first_external=
      selected_external_bath<0 ? 0 : selected_external_bath;
    const int last_external=
      selected_external_bath<0 ? cix.Nbath : selected_external_bath+1;
    for (int bath=first_external;bath<last_external;++bath) {
      const int signed_channel=cix.bath_sigind0(bath);
      external.push_back({operator_cache[2*bath],operator_cache[2*bath+1],
                          MakeBullaTranspose(cix,bath),
                          std::abs(signed_channel),
                          static_cast<double>(SignedIndexSign(signed_channel)),
                          IsF01Channel(cix,bath)});
    }

    for (int pole_pair=first_pole_pair;
         pole_pair<last_pole_pair;++pole_pair) {
        const int l0=pole_pair/pole_count;
        const int l1=pole_pair%pole_count;
        const std::array<int,2> pole{{l0,l1}};
        std::array<VectorXd,6> edge_factor;
        for (auto& x:edge_factor) x=VectorXd::Ones(rank);
        std::array<double,2> pole_prefactor{{1.0,1.0}};

        for (int split=0;split<2;++split) {
          const auto& line=cut.split[split];
          const double omega=workspace.poles()(pole[split]);
          const bool direct=
            (line.direction==0 && omega>=0.0) ||
            (line.direction==1 && omega<0.0);
          const int frequency_sign=direct
            ? (line.direction==0 ? +1 : -1)
            : (line.direction==0 ? -1 : +1);
          int edge_count=0;
          for (int edge=0;edge<6;++edge) {
            const bool selected=direct
              ? (line.first_vertex<=edge && edge<line.second_vertex)
              : (edge<line.first_vertex || edge>=line.second_vertex);
            if (!selected) continue;
            ++edge_count;
            edge_factor[edge].array()*=
              (frequency_sign>0
                 ? workspace.kernel_plus().col(pole[split])
                 : workspace.kernel_minus().col(pole[split])).array();
          }
          pole_prefactor[split]=
            -1.0/std::pow(KernelZero(workspace,pole[split],frequency_sign),
                          edge_count-1);
        }

        // The edge containing the external time is the boundary between the
        // two independent backbones and is never convolved.  Prepare and
        // deduplicate only the remaining convolution edges.
        std::vector<VectorXd> unique_scalars;
        std::vector<Engine::PreparedPseudoCoefficients> prepared;
        std::array<int,6> prepared_index{{-1,-1,-1,-1,-1,-1}};
        for (int edge=1;edge<6;++edge) {
          const bool convolved=edge<cut.external_vertex ||
                               edge>cut.external_vertex;
          if (!convolved) continue;
          int found=-1;
          for (int i=0;i<static_cast<int>(unique_scalars.size());++i)
            if ((unique_scalars[i].array()==edge_factor[edge].array()).all()) {
              found=i;
              break;
            }
          if (found<0) {
            found=static_cast<int>(unique_scalars.size());
            unique_scalars.push_back(edge_factor[edge]);
            prepared.push_back(
              Engine::PreparePseudoCoefficients(workspace,edge_factor[edge]));
          }
          prepared_index[edge]=found;
        }

        for (std::size_t pair_index=0;pair_index<plans.size();++pair_index) {
          BathPlan& plan=plans[pair_index];
          double prefactor=plan.fixed_prefactor;
          for (int split=0;split<2;++split)
            prefactor*=pole_prefactor[split]*workspace.residues()(
              plan.channel[split],pole[split]);
          if (prefactor==0.0) continue;

          SparseFunction first=Engine::PseudoValues(
            pseudo_blocks,edge_factor[0]);
          for (int vertex=1;vertex<cut.external_vertex && !first.empty();
               ++vertex) {
            const int op=plan.vertex_operator[vertex];
            if (op<0)
              throw std::logic_error("missing TCA physical split operator");
            first=Engine::MultiplyRight(first,operator_cache[op],ones);
            first=Engine::ConvolveRightPseudoPrepared(
              first,cix,workspace,prepared[prepared_index[vertex]]);
          }

          SparseFunction second=Engine::PseudoValues(
            pseudo_blocks,edge_factor[cut.external_vertex]);
          bool second_reflected=false;
          for (int vertex=cut.external_vertex+1;
               vertex<6 && !second.empty();++vertex) {
            const int op=plan.vertex_operator[vertex];
            if (op<0)
              throw std::logic_error("missing TCA physical split operator");
            second=Engine::MultiplyRight(second,operator_cache[op],ones);
            second=Engine::ConvolveRightPseudoPrepared(
              second,cix,workspace,prepared[prepared_index[vertex]],
              vertex==5);
            if (vertex==5) second_reflected=true;
          }
          // The last convolution (vertex five) directly returns values at
          // beta-tau, so no separate DLR reflection pass is needed here.
          if (!second.empty() && !second_reflected)
            Reflect(second,workspace.basis());
          if (!first.empty() && !second.empty()) plan.nonempty=true;

          for (const auto& pair:external)
            AccumulateBoundaryEstimators(
              first,second,pair.creation,pair.annihilation,pair.bulla,
              prefactor*pair.bath_sign,pair.channel,pair.anomalous,
              batch.values);
        }
      }

    for (const auto& plan:plans)
      if (plan.nonempty) ++batch.nonempty_pairs;
    return batch;
  }

private:
  struct PhysicalCut {
    std::array<RootedCut::Line,2> split;
    int external_vertex=0;
  };

  static PhysicalCut MakePhysicalCut(
      int topology,const std::array<int,3>& direction,int cut_line)
  {
    if (topology<0 || topology>=2 || cut_line<0 || cut_line>=3)
      throw std::out_of_range("invalid summed TCA physical cut");
    PhysicalCut result;
    const int a=edmft::tca::topology_endpoints[topology][cut_line][0];
    const int b=edmft::tca::topology_endpoints[topology][cut_line][1];
    const int start=direction[cut_line]==0 ? b : a;
    const int finish=direction[cut_line]==0 ? a : b;
    result.external_vertex=(finish-start+6)%6;
    int split=0;
    for (int line=0;line<3;++line) {
      if (line==cut_line) continue;
      int x=(edmft::tca::topology_endpoints[topology][line][0]-start+6)%6;
      int y=(edmft::tca::topology_endpoints[topology][line][1]-start+6)%6;
      int d=direction[line];
      if (x>y) { std::swap(x,y); d^=1; }
      result.split[split++]={line,x,y,d};
    }
    return result;
  }

  static bool IsF01Channel(const Cix& cix,int orbital)
  {
    const bool destroys_pair=
      cix.Fi_indx(2*orbital)%2==1 && cix.Fi_indx(2*orbital+1)%2==1;
    const bool normal_down=cix.SC_bath(orbital) &&
      ((cix.Fi_indx(2*orbital)+cix.Fi_indx(2*orbital+1))%2==1) &&
      cix.bath_spin(orbital)<0;
    return destroys_pair || normal_down;
  }

  static OperatorTranspose MakeBullaTranspose(const Cix& cix,int orbital)
  {
    OperatorTranspose result;
    result.by_source.assign(cix.Nblock,-1);
    result.by_target.resize(cix.Nblock);
    const OperatorTranspose ordinary=
      Engine::MakeOperatorTranspose(cix,2*orbital);
    for (int source=0;source<cix.Nblock;++source) {
      const int ordinary_index=ordinary.by_source[source];
      if (ordinary_index<0) continue;
      const auto& base=ordinary.entries[ordinary_index];
      const NjmResult nx=cix.Njm_ind(source,2*orbital);
      if (nx.njm_c==0) continue;
      MatrixXd matrix=MatrixXd::Zero(base.matrix.rows(),base.matrix.cols());
      if (nx.njm_c==1) {
        matrix=nx.njm_r*base.matrix;
      } else if (nx.njm_c==2 && nx.njm!=nullptr) {
        const int njm_rows=nx.transposed
          ? nx.njm->extent(1) : nx.njm->extent(0);
        const int njm_cols=nx.transposed
          ? nx.njm->extent(0) : nx.njm->extent(1);
        if (njm_rows!=matrix.rows() || njm_cols!=matrix.cols())
          throw std::logic_error(
            "summed TCA N_jm block dimensions differ from bath operator");
        for (int i=0;i<matrix.rows();++i)
          for (int j=0;j<matrix.cols();++j)
            // The fixed-plan product basis contains one index only for a
            // nonzero transition of the ordinary one-body operator.  Its
            // dense Bulla connector may mix those transition indices, but
            // it cannot introduce a source row or target column absent from
            // that support.  Projecting here makes the full block matrix
            // exactly equivalent to BullaConnector(BlockPlan,...).
            if (base.matrix.row(i).cwiseAbs().maxCoeff()!=0.0 &&
                base.matrix.col(j).cwiseAbs().maxCoeff()!=0.0)
              matrix(i,j)=nx.transposed
                ? (*nx.njm)(j,i) : (*nx.njm)(i,j);
      } else {
        throw std::logic_error("unsupported summed TCA N_jm representation");
      }
      const int index=static_cast<int>(result.entries.size());
      Engine::OperatorEntry entry;
      entry.source=source;
      entry.target=base.target;
      entry.matrix=std::move(matrix);
      // N_jm may mix transition indices, so it deliberately retains the
      // dense fallback rather than claiming partial-permutation structure.
      entry.signed_partial_permutation=false;
      for (int i=0;i<entry.matrix.rows();++i)
        for (int j=0;j<entry.matrix.cols();++j)
          if (entry.matrix(i,j)!=0.0)
            entry.nonzero.push_back({i,j,entry.matrix(i,j)});
      result.entries.push_back(std::move(entry));
      result.by_source[source]=index;
      result.by_target[base.target].push_back(index);
    }
    return result;
  }

  // Close the two physical backbones for one external bath and accumulate
  // both the ordinary Green's function and its Bulla estimator.  Sparse block
  // compatibility is resolved once, outside the DLR-node loop.  An ordinary
  // creation entry and its Bulla replacement have identical source/target
  // blocks by construction; the Bulla entry may simply be absent when N_jm
  // is zero for that source block.
  static void AccumulateBoundaryEstimators(
      const SparseFunction& first,const SparseFunction& second,
      const OperatorTranspose& creation,
      const OperatorTranspose& annihilation,
      const OperatorTranspose& bulla,
      double prefactor,int channel,bool anomalous,Result& output)
  {
    struct Closure {
      const Engine::OperatorEntry* creation=nullptr;
      const Engine::OperatorEntry* bulla=nullptr;
      const MatrixFunction* second=nullptr;
    };

    for (const auto& [key,values]:first) {
      if (key.first<0 || key.first>=
            static_cast<int>(creation.by_target.size())) continue;
      if (key.second<0 || key.second>=
            static_cast<int>(annihilation.by_source.size())) continue;
      const int right_index=annihilation.by_source[key.second];
      if (right_index<0) continue;
      const auto& right=annihilation.entries[right_index];

      std::vector<Closure> closures;
      closures.reserve(creation.by_target[key.first].size());
      for (const int left_index:creation.by_target[key.first]) {
        const auto& left=creation.entries[left_index];
        const auto reflected=second.find({right.target,left.source});
        if (reflected==second.end()) continue;

        const Engine::OperatorEntry* bulla_entry=nullptr;
        if (left.source>=0 &&
            left.source<static_cast<int>(bulla.by_source.size())) {
          const int bulla_index=bulla.by_source[left.source];
          if (bulla_index>=0) {
            const auto& replacement=bulla.entries[bulla_index];
            // MakeBullaTranspose retains the ordinary connector's target.
            // Keep the check here because using a mismatched block would
            // otherwise surface only as a cryptic Eigen dimension failure.
            if (replacement.target!=key.first)
              throw std::logic_error(
                "TCA Bulla connector changed creation target block");
            bulla_entry=&replacement;
          }
        }
        closures.push_back(
          Closure{&left,bulla_entry,&reflected->second});
      }
      if (closures.empty()) continue;

      // G_left*R is independent of the creation source block and is shared
      // by every compatible closure.  The remaining product is sampled only
      // at entries selected by the sparse creation/Bulla connector.
      for (int node=0;node<static_cast<int>(values.size());++node) {
        const MatrixXd left_annihilation=
          Engine::MultiplyRightOperator(values[node],right);
        for (const Closure& closure:closures) {
          const MatrixXd& reflected=(*closure.second)[node];
          output.g(node,channel)+=prefactor*
            Engine::TraceLeftProduct(
              *closure.creation,left_annihilation,reflected);
          if (closure.bulla) {
            MatrixXd& estimator=anomalous ? output.f01 : output.f00;
            estimator(node,channel)+=prefactor*
              Engine::TraceLeftProduct(
                *closure.bulla,left_annihilation,reflected);
          }
        }
      }
    }
  }

  static void Reflect(SparseFunction& function,
                      const dlr::FermionicBasis& basis)
  {
    for (auto& [key,values]:function) {
      if (values.empty()) continue;
      const int rows=values[0].rows(),cols=values[0].cols();
      MatrixXd flat(basis.rank(),rows*cols);
      for (int node=0;node<basis.rank();++node)
        for (int i=0;i<rows;++i)
          for (int j=0;j<cols;++j)
            flat(node,i*cols+j)=values[node](i,j);
      const MatrixXd reflected_flat=basis.reflect_values(flat);
      MatrixFunction reflected(
        basis.rank(),MatrixXd::Zero(rows,cols));
      for (int node=0;node<basis.rank();++node)
        for (int i=0;i<rows;++i)
          for (int j=0;j<cols;++j)
            reflected[node](i,j)=reflected_flat(node,i*cols+j);
      values=std::move(reflected);
    }
  }
};

struct SelfEnergyNodeResult {
  MatrixXd sigma;
  std::array<MatrixXd,2> sigma_by_topology;
  MatrixXd g;
  MatrixXd f00;
  MatrixXd f01;
  std::array<MatrixXd,6> cuts;
  // Optional diagnostic decomposition.  cuts[] is retained for existing
  // callers and contains the sum of both topologies.  This second array keeps
  // every topology/cut contribution separate, which is essential when
  // checking a proposed automorphism of an open pseudo-particle diagram.
  std::array<std::array<MatrixXd,6>,2> cuts_by_topology;
  std::array<std::uint64_t,2> nonempty_plans_by_topology{{0,0}};
  std::uint64_t candidate_plans=0;
  bool physical_enabled=false;
};

// Physical third-order diagrams are obtained by cutting one of the three
// hybridization lines of Phi.  The cut leaves all six atomic operators and
// all six pseudo propagators, but removes the selected Delta.  Its two
// endpoints become the physical creation vertex at time zero and annihilation
// vertex at time T.  The remaining four integration vertices separate into
// two independent ordered backbones, one of total length T and one of length
// beta-T, as in Appendix B of SciPost Phys. 19, 121 (2025).
class RootedPhysicalEvaluator {
public:
  struct Result {
    MatrixXd g;
    MatrixXd f00;
    MatrixXd f01;
  };

  static Result evaluate_all_lines_nodes(
      const tca_matrix::AtomicTracePlan& plan,
      const Cix& cix,
      const imag_dlr::Workspace& workspace,
      std::array<double,3>* line_seconds=nullptr,
      double* prepare_seconds=nullptr,
      int selected_cut_line=-1)
  {
    if (!workspace.has_pseudos())
      throw std::logic_error("TCA physical evaluation requires pseudos");
    const int rank=workspace.rank();
    Result result{MatrixXd::Zero(rank,workspace.channel_count()),
                  MatrixXd::Zero(rank,workspace.channel_count()),
                  MatrixXd::Zero(rank,workspace.channel_count())};

    std::array<int,3> direction{},orbital{},channel{};
    double bath_sign=1.0;
    for (int line=0;line<3;++line) {
      direction[line]=plan.dirspin()[line]&1;
      orbital[line]=plan.dirspin()[line]/2;
      const int signed_channel=cix.bath_sigind0(orbital[line]);
      channel[line]=std::abs(signed_channel);
      bath_sign*=SignedIndexSign(signed_channel);
    }
    const auto prepare_begin=std::chrono::steady_clock::now();
    std::vector<std::array<MatrixFunction,6>> prepared(plan.blocks().size());
    for (int block=0;block<static_cast<int>(plan.blocks().size());++block)
      for (int edge=0;edge<6;++edge)
        prepared[block][edge]=PropagatorValues(
          plan.blocks()[block].propagator[edge].signed_label,workspace);
    if (prepare_seconds)
      *prepare_seconds=std::chrono::duration<double>(
        std::chrono::steady_clock::now()-prepare_begin).count();
    if (line_seconds) line_seconds->fill(0.0);

    if (selected_cut_line<-1 || selected_cut_line>=3)
      throw std::out_of_range("selected TCA physical cut line is invalid");
    const int first_line=selected_cut_line<0 ? 0 : selected_cut_line;
    const int last_line=selected_cut_line<0 ? 3 : selected_cut_line+1;
    for (int cut_line=first_line;cut_line<last_line;++cut_line) {
      const auto line_begin=std::chrono::steady_clock::now();
      const auto rooted=MakePhysicalCut(plan.topology(),direction,cut_line);
      const int output_channel=channel[cut_line];
      const bool anomalous_output=IsF01Channel(cix,orbital[cut_line]);

      for (int l0=0;l0<workspace.pole_count();++l0)
        for (int l1=0;l1<workspace.pole_count();++l1) {
          const std::array<int,2> pole{{l0,l1}};
          std::array<VectorXd,6> edge_factor;
          for (auto& x:edge_factor) x=VectorXd::Ones(rank);
          double prefactor=
            bath_sign*edmft::tca::topology_crossing_sign[plan.topology()]*
            edmft::tca::topology_symmetry_factor[plan.topology()];

          for (int split=0;split<2;++split) {
            const auto& line=rooted.split[split];
            const double omega=workspace.poles()(pole[split]);
            const double residue=workspace.residues()(
              channel[line.original_line],pole[split]);
            // On a circular physical diagram every separated line has two
            // possible paths between its endpoints.  Choose the path whose
            // effective frequency has K(0)=O(1): this is precisely the four
            // stable cases in Appendix B, but written without topology-
            // specific pictures.  Direct d=0 / positive-pole and direct d=1
            // / negative-pole lines use edges [a,b); the other two cases use
            // the complementary edges [0,a) U [b,6).  Products on the two
            // sides of the external vertex automatically become the T and
            // beta-T backbones, with no dangerous K(T)/K(0) ratio.
            const bool direct=
              (line.direction==0 && omega>=0.0) ||
              (line.direction==1 && omega<0.0);
            const int frequency_sign=direct
              ? (line.direction==0 ? +1 : -1)
              : (line.direction==0 ? -1 : +1);
            int edge_count=0;
            for (int edge=0;edge<6;++edge) {
              const bool selected=direct
                ? (line.first_vertex<=edge && edge<line.second_vertex)
                : (edge<line.first_vertex || edge>=line.second_vertex);
              if (!selected) continue;
              ++edge_count;
              edge_factor[edge].array()*=
                (frequency_sign>0
                   ? workspace.kernel_plus().col(pole[split])
                   : workspace.kernel_minus().col(pole[split])).array();
            }
            if (edge_count<=0)
              throw std::logic_error("empty stable TCA physical line path");
            prefactor*=-residue/std::pow(
              KernelZero(workspace,pole[split],frequency_sign),edge_count-1);
          }

          for (int block_index=0;
               block_index<static_cast<int>(plan.blocks().size());
               ++block_index) {
            const auto& block=plan.blocks()[block_index];
            const MatrixFunction ordinary=EvaluateBlock(
              block,rooted,workspace,edge_factor,
              prepared[block_index],nullptr);
            AccumulateTrace(ordinary,prefactor,output_channel,result.g);

            const int original_vertex=rooted.original_vertex[0];
            const int incoming_block=
              block.propagated_block[(original_vertex+5)%6];
            const NjmResult nx=
              cix.Njm_ind(incoming_block,2*orbital[cut_line]);
            if (nx.njm_c==1) {
              // Here N_jm is exactly nx.njm_r times the ordinary connector.
              // Every later operation is linear, so the complete Bulla trace
              // is the same scalar times the already computed ordinary trace.
              AccumulateTrace(ordinary,prefactor*nx.njm_r,output_channel,
                              anomalous_output ? result.f01 : result.f00);
            } else if (nx.njm_c==2) {
              const MatrixXd bulla=BullaConnector(
                block,original_vertex,cix,orbital[cut_line]);
              if (bulla.size()==0 || bulla.cwiseAbs().maxCoeff()==0.0)
                continue;
              const MatrixFunction replaced=EvaluateBlock(
                block,rooted,workspace,edge_factor,
                prepared[block_index],&bulla);
              AccumulateTrace(replaced,prefactor,output_channel,
                              anomalous_output ? result.f01 : result.f00);
            }
          }
        }
      if (line_seconds)
        (*line_seconds)[cut_line]=std::chrono::duration<double>(
          std::chrono::steady_clock::now()-line_begin).count();
    }
    return result;
  }

private:
  struct PhysicalCut {
    std::array<int,6> original_vertex{{0,1,2,3,4,5}};
    std::array<int,6> original_edge{{0,1,2,3,4,5}};
    std::array<RootedCut::Line,2> split;
    int external_vertex=0;
  };

  static PhysicalCut MakePhysicalCut(
      int topology,const std::array<int,3>& direction,int cut_line)
  {
    PhysicalCut result;
    const int a=edmft::tca::topology_endpoints[topology][cut_line][0];
    const int b=edmft::tca::topology_endpoints[topology][cut_line][1];
    // op=2*orbital is the creation-side operator in the CIX convention.
    // For direction zero it is at b; reversing the line puts it at a.
    const int start=direction[cut_line]==0 ? b : a;
    const int finish=direction[cut_line]==0 ? a : b;
    result.external_vertex=(finish-start+6)%6;
    if (result.external_vertex<=0 || result.external_vertex>=6)
      throw std::logic_error("invalid TCA physical external separation");
    for (int local=0;local<6;++local) {
      result.original_vertex[local]=(start+local)%6;
      result.original_edge[local]=(start+local)%6;
    }
    int split=0;
    for (int line=0;line<3;++line) {
      if (line==cut_line) continue;
      int x=(edmft::tca::topology_endpoints[topology][line][0]-start+6)%6;
      int y=(edmft::tca::topology_endpoints[topology][line][1]-start+6)%6;
      int d=direction[line];
      if (x>y) { std::swap(x,y); d^=1; }
      result.split[split++]={line,x,y,d};
    }
    return result;
  }

  static bool IsF01Channel(const Cix& cix,int orbital)
  {
    const bool destroys_pair=
      cix.Fi_indx(2*orbital)%2==1 && cix.Fi_indx(2*orbital+1)%2==1;
    const bool normal_down=cix.SC_bath(orbital) &&
      ((cix.Fi_indx(2*orbital)+cix.Fi_indx(2*orbital+1))%2==1) &&
      cix.bath_spin(orbital)<0;
    return destroys_pair || normal_down;
  }

  static MatrixXd OrdinaryConnector(
      const tca_matrix::BlockPlan& block,int original_vertex)
  {
    return block.vertex_sign[original_vertex].asDiagonal();
  }

  static MatrixXd BullaConnector(const tca_matrix::BlockPlan& block,
                                 int original_vertex,
                                 const Cix& cix,
                                 int orbital)
  {
    const int incoming_block=
      block.propagated_block[(original_vertex+5)%6];
    const NjmResult nx=cix.Njm_ind(incoming_block,2*orbital);
    const int n=block.vertex_sign[original_vertex].size();
    MatrixXd connector=MatrixXd::Zero(n,n);
    if (nx.njm_c==0) return connector;
    if (nx.njm_c==1)
      return nx.njm_r*OrdinaryConnector(block,original_vertex);
    if (nx.njm_c!=2 || nx.njm==nullptr)
      throw std::logic_error("unsupported TCA N_jm representation");
    for (int incoming=0;incoming<n;++incoming)
      for (int outgoing=0;outgoing<n;++outgoing) {
        const int i=block.incoming_state[original_vertex](incoming);
        const int j=block.outgoing_state[original_vertex](outgoing);
        connector(incoming,outgoing)=nx.transposed
          ? (*nx.njm)(j,i) : (*nx.njm)(i,j);
      }
    return connector;
  }

  static MatrixFunction PropagatorValues(
      const Eigen::MatrixXi& labels,const imag_dlr::Workspace& workspace)
  {
    MatrixFunction values(
      workspace.rank(),MatrixXd::Zero(labels.rows(),labels.cols()));
    for (int i=0;i<labels.rows();++i)
      for (int j=0;j<labels.cols();++j) {
        const int label=labels(i,j);
        if (label==0) continue;
        for (int node=0;node<workspace.rank();++node)
          values[node](i,j)=SignedIndexSign(label)*
            workspace.pseudo_values()(node,std::abs(label)-1);
      }
    return values;
  }

  static void MultiplyLeft(MatrixFunction& values,const MatrixXd& matrix)
  {
    for (int node=0;node<static_cast<int>(values.size());++node)
      values[node]=(matrix*values[node]).eval();
  }

  static void MultiplyRight(MatrixFunction& values,const MatrixXd& matrix)
  {
    for (int node=0;node<static_cast<int>(values.size());++node)
      values[node]=(values[node]*matrix).eval();
  }

  static MatrixFunction ReflectValues(const MatrixFunction& values,
                                      const dlr::FermionicBasis& basis)
  {
    MatrixFunction reflected(
      basis.rank(),MatrixXd::Zero(values[0].rows(),values[0].cols()));
    const MatrixXd& reflection=basis.reflection_values_matrix();
    for (int i=0;i<basis.rank();++i)
      for (int j=0;j<basis.rank();++j)
        reflected[i].noalias()+=reflection(i,j)*values[j];
    return reflected;
  }

  // Return Tr[A(T) B(beta-T)] as a one-by-one matrix function.  The optional
  // replacement is the N_jm connector at local vertex zero.
  static MatrixFunction EvaluateBlock(
      const tca_matrix::BlockPlan& block,
      const PhysicalCut& cut,
      const imag_dlr::Workspace& workspace,
      const std::array<VectorXd,6>& edge_factor,
      const std::array<MatrixFunction,6>& prepared,
      const MatrixXd* replacement)
  {
    const auto& basis=workspace.basis();
    const int q=cut.external_vertex;
    std::array<MatrixXd,6> connector;
    std::array<MatrixFunction,6> edge_coefficients;
    for (int local=0;local<6;++local) {
      connector[local]=OrdinaryConnector(block,cut.original_vertex[local]);
      if (local==0 || local==q) continue;
      MatrixFunction edge=prepared[cut.original_edge[local]];
      for (int node=0;node<workspace.rank();++node)
        edge[node]*=edge_factor[local](node);
      edge_coefficients[local]=basis.matrix_values_to_coefficients(edge);
    }
    if (replacement) connector[0]=*replacement;

    // The first edge on each side is already represented at the DLR time
    // nodes.  Avoid an identity values->coefficients->values round trip for
    // local edges zero and q.
    MatrixFunction first=prepared[cut.original_edge[0]];
    for (int node=0;node<workspace.rank();++node)
      first[node]*=edge_factor[0](node);
    MultiplyLeft(first,connector[0]);
    MultiplyRight(first,connector[1]);
    for (int edge=1;edge<q;++edge) {
      first=basis.time_ordered_convolution_values(
        basis.matrix_values_to_coefficients(first),edge_coefficients[edge]);
      MultiplyRight(first,connector[edge+1]);
    }

    MatrixFunction second=prepared[cut.original_edge[q]];
    for (int node=0;node<workspace.rank();++node)
      second[node]*=edge_factor[q](node);
    for (int edge=q+1;edge<6;++edge) {
      MultiplyRight(second,connector[edge]);
      second=basis.time_ordered_convolution_values(
        basis.matrix_values_to_coefficients(second),edge_coefficients[edge]);
    }
    const MatrixFunction reflected=ReflectValues(second,basis);
    MatrixFunction trace(workspace.rank(),MatrixXd::Zero(1,1));
    for (int node=0;node<workspace.rank();++node)
      trace[node](0,0)=(first[node]*reflected[node]).trace();
    return trace;
  }

  static void AccumulateTrace(const MatrixFunction& trace,double prefactor,
                              int output_channel,MatrixXd& output)
  {
    for (int node=0;node<static_cast<int>(trace.size());++node)
      output(node,output_channel)+=prefactor*trace[node](0,0);
  }
};

// Correctness-first complete TCA self-energy driver.  It fixes the direction
// of canonical root (0,3), sums every bath carried by that root, and includes
// both directions of each remaining line.  Fixing the root removes the
// twofold arrow reversal.  The residual threefold automorphism of topology B
// is independent of that arrow choice and is applied by
// topology_symmetry_factor.
//
// The implementation intentionally keeps no catalog of AtomicTracePlan
// objects.  A plan is compiled, evaluated, and released immediately, which is
// important for the large SC/f-shell cases.  MPI distributes precisely these
// independent fixed-(topology,bath triple,orientation) work units and reduces
// only the compact DLR-node outputs.
struct TCAPlanProfile {
  double total_seconds=0.0;
  double compile_seconds=0.0;
  double self_prepare_seconds=0.0;
  std::array<double,6> self_cut_seconds{{0,0,0,0,0,0}};
  double physical_prepare_seconds=0.0;
  std::array<double,3> physical_line_seconds{{0,0,0}};
};

// A plain candidate%nranks assignment couples MPI ownership to the four
// direction pairs in the innermost loop.  In particular, with 12 ranks each
// rank receives only one direction pair.  SplitMix64 breaks that arithmetic
// correlation while remaining deterministic and requiring no communication.
inline std::uint64_t MixedTCACandidate(std::uint64_t value) noexcept
{
  value+=0x9e3779b97f4a7c15ULL;
  value=(value^(value>>30))*0xbf58476d1ce4e5b9ULL;
  value=(value^(value>>27))*0x94d049bb133111ebULL;
  return value^(value>>31);
}

inline int TCACandidateOwner(std::uint64_t candidate,int mpi_size) noexcept
{
  return static_cast<int>(
    MixedTCACandidate(candidate)%static_cast<std::uint64_t>(mpi_size));
}

inline std::uint64_t TCAAssignedCandidateCount(
    std::uint64_t total,int mpi_rank,int mpi_size) noexcept
{
  std::uint64_t count=0;
  for (std::uint64_t candidate=0;candidate<total;++candidate)
    if (TCACandidateOwner(candidate,mpi_size)==mpi_rank) ++count;
  return count;
}

// Choose a small pole-pair batch while retaining the complete bath sum in
// every task.  This ordering is important: pseudo-propagator preparation is
// now performed once per pole pair and reused for all bath pairs, instead of
// being repeated by every bath tile.  The automatic choice deliberately
// overdecomposes the phase into at least about 32 tasks per rank.  It also
// caps the work represented by one tile in units of
//
//                 pole pairs * Nbath^2.
//
// That second bound matters on a laptop: 34 pole pairs are a modest tile for
// a few baths but an opaque, several-minute task for the 12-bath SC cluster.
// The chosen cap gives 18 pole pairs for Nbath=8, 8 for Nbath=12, and 5 for
// Nbath=14.  At hundreds of ranks the task-count bound naturally drives the
// value to one.  A positive requested value remains an exact user override;
// zero selects this automatic policy.
inline int TCAChooseMicroBatchSize(std::uint64_t total_groups,
                                   int pole_pairs_per_structure,
                                   int mpi_size,int requested,int bath_count)
{
  if (pole_pairs_per_structure<1 || mpi_size<1 || requested<0 ||
      bath_count<1)
    throw std::invalid_argument("invalid TCA MPI micro-batch parameters");
  if (requested>0)
    return std::min(requested,pole_pairs_per_structure);
  constexpr std::uint64_t target_tasks_per_rank=32;
  const std::uint64_t target_tasks=
    target_tasks_per_rank*static_cast<std::uint64_t>(mpi_size);
  const std::uint64_t task_limited=std::max<std::uint64_t>(
    1,(total_groups+target_tasks-1)/target_tasks);
  constexpr std::uint64_t target_pole_bath_units_per_task=1152;
  const std::uint64_t bath_pairs=
    static_cast<std::uint64_t>(bath_count)*bath_count;
  const std::uint64_t bath_limited=std::max<std::uint64_t>(
    1,target_pole_bath_units_per_task/bath_pairs);
  const std::uint64_t automatic=std::min(task_limited,bath_limited);
  return static_cast<int>(std::min<std::uint64_t>(
    automatic,static_cast<std::uint64_t>(pole_pairs_per_structure)));
}

#ifdef _MPI
inline void TCAReportDynamicBalance(std::ostream* progress,const char* phase,
                                    std::uint64_t local_groups,
                                    double local_compute_seconds,
                                    double phase_seconds,int mpi_rank,
                                    int mpi_size,MPI_Comm communicator)
{
  std::uint64_t group_min=0,group_max=0,group_sum=0;
  double compute_min=0.0,compute_max=0.0,compute_sum=0.0;
  MPI_Reduce(&local_groups,&group_min,1,MPI_UINT64_T,MPI_MIN,0,communicator);
  MPI_Reduce(&local_groups,&group_max,1,MPI_UINT64_T,MPI_MAX,0,communicator);
  MPI_Reduce(&local_groups,&group_sum,1,MPI_UINT64_T,MPI_SUM,0,communicator);
  MPI_Reduce(&local_compute_seconds,&compute_min,1,MPI_DOUBLE,MPI_MIN,0,
             communicator);
  MPI_Reduce(&local_compute_seconds,&compute_max,1,MPI_DOUBLE,MPI_MAX,0,
             communicator);
  MPI_Reduce(&local_compute_seconds,&compute_sum,1,MPI_DOUBLE,MPI_SUM,0,
             communicator);
  if (mpi_rank==0 && progress) {
    const double average=compute_sum/std::max(1,mpi_size);
    const double efficiency=compute_sum/
      std::max(1.e-12,phase_seconds*static_cast<double>(mpi_size));
    (*progress) << "TCA DLR dynamic " << phase
                << " balance: work_units_per_rank=[" << group_min << ','
                << static_cast<double>(group_sum)/std::max(1,mpi_size) << ','
                << group_max << "] compute_seconds_per_rank=["
                << compute_min << ',' << average << ',' << compute_max
                << "] phase_seconds=" << phase_seconds
                << " parallel_efficiency=" << efficiency << std::endl;
  }
}

struct TCAHybridScheduleResult {
  std::uint64_t local_groups=0;
  double local_compute_seconds=0.0;
  double phase_seconds=0.0;
};

// Portable hybrid scheduler.  Three quarters of the deterministically mixed
// micro-batches are first distributed exactly evenly over every rank.  This
// keeps rank zero productive and needs no MPI progress while the expensive
// contractions run.  When rank zero finishes its bulk share it becomes a
// dispatcher for the reserved tail; workers request one tail micro-batch at a
// time, so early bulk finishers absorb work from slow ranks.  Only the tail
// sacrifices rank zero, limiting that cost to 1/(4*nranks), while avoiding
// passive-target RMA progress assumptions that are not portable across MPI
// implementations.
template <class Work>
TCAHybridScheduleResult RunTCAHybridSchedule(
    std::uint64_t total_tasks,std::uint64_t total_groups,
    int mpi_rank,int mpi_size,MPI_Comm communicator,
    std::ostream* progress,double progress_seconds,const char* phase,
    Work&& work)
{
  if (mpi_size<2 || total_tasks==0)
    throw std::invalid_argument("hybrid TCA scheduler needs MPI work");
  std::vector<std::uint64_t> order(total_tasks);
  std::iota(order.begin(),order.end(),0);
  std::sort(order.begin(),order.end(),[](std::uint64_t a,std::uint64_t b) {
    const auto ha=MixedTCACandidate(a);
    const auto hb=MixedTCACandidate(b);
    return ha<hb || (ha==hb && a<b);
  });
  // Round down to a complete rank-wide round so every rank receives exactly
  // the same number of bulk tasks.  At least one quarter remains available to
  // correct differences in actual matrix-contraction cost.
  std::uint64_t bulk_tasks=(3*total_tasks)/4;
  bulk_tasks-=bulk_tasks%static_cast<std::uint64_t>(mpi_size);
  if (bulk_tasks==0) bulk_tasks=std::min<std::uint64_t>(
    total_tasks,static_cast<std::uint64_t>(mpi_size));
  const auto begin=std::chrono::steady_clock::now();
  auto last_report=begin;
  TCAHybridScheduleResult summary;
  const auto run_one=[&](std::uint64_t task) {
    const auto task_begin=std::chrono::steady_clock::now();
    const std::uint64_t groups=work(task);
    summary.local_compute_seconds+=std::chrono::duration<double>(
      std::chrono::steady_clock::now()-task_begin).count();
    summary.local_groups+=groups;
    return groups;
  };

  // The dispatcher cannot report global progress until the dynamic tail
  // starts, because every rank works independently during this bulk stage.
  // Still emit a rank-zero heartbeat: all ranks receive the same number of
  // bulk tasks (apart from the defensive total_tasks<mpi_size case), so rank
  // zero's fraction is a useful indication of progress without adding any
  // synchronization to the expensive contraction loop.  The ETA is labelled
  // as bulk-only because the dynamically balanced tail has not started yet.
  const std::uint64_t local_bulk_tasks=
    bulk_tasks>static_cast<std::uint64_t>(mpi_rank) ?
      1+(bulk_tasks-1-static_cast<std::uint64_t>(mpi_rank))/
          static_cast<std::uint64_t>(mpi_size) : 0;
  std::uint64_t local_bulk_completed=0;
  for (std::uint64_t position=mpi_rank;position<bulk_tasks;
       position+=mpi_size) {
    run_one(order[position]);
    ++local_bulk_completed;
    if (mpi_rank==0 && progress && progress_seconds>0.0) {
      const auto now=std::chrono::steady_clock::now();
      const double since=
        std::chrono::duration<double>(now-last_report).count();
      if (since>=progress_seconds ||
          local_bulk_completed==local_bulk_tasks) {
        const double elapsed=
          std::chrono::duration<double>(now-begin).count();
        const double rate=local_bulk_completed/std::max(elapsed,1.e-12);
        (*progress) << "TCA DLR hybrid " << phase
                    << " bulk progress: rank0_bulk_completed="
                    << local_bulk_completed << "/" << local_bulk_tasks
                    << " ("
                    << 100.0*local_bulk_completed/
                         std::max<std::uint64_t>(1,local_bulk_tasks)
                    << "%) elapsed_seconds=" << elapsed
                    << " estimated_bulk_remaining_seconds="
                    << (local_bulk_tasks-local_bulk_completed)/rate
                    << std::endl;
        last_report=now;
      }
    }
  }

  constexpr int request_tag=27001;
  constexpr int assignment_tag=27002;
  const std::uint64_t stop_task=total_tasks;
  if (mpi_rank==0) {
    std::uint64_t next=bulk_tasks;
    std::uint64_t global_completed=summary.local_groups;
    int workers_remaining=mpi_size-1;
    while (workers_remaining>0) {
      std::uint64_t newly_completed=0;
      MPI_Status status;
      MPI_Recv(&newly_completed,1,MPI_UINT64_T,MPI_ANY_SOURCE,request_tag,
               communicator,&status);
      global_completed+=newly_completed;
      const std::uint64_t assignment=
        next<total_tasks ? order[next++] : stop_task;
      MPI_Send(&assignment,1,MPI_UINT64_T,status.MPI_SOURCE,assignment_tag,
               communicator);
      if (assignment==stop_task) --workers_remaining;

      if (progress && progress_seconds>0.0) {
        const auto now=std::chrono::steady_clock::now();
        const double since=
          std::chrono::duration<double>(now-last_report).count();
        if (since>=progress_seconds || global_completed==total_groups) {
          const double elapsed=
            std::chrono::duration<double>(now-begin).count();
          const double rate=global_completed/std::max(elapsed,1.e-12);
          (*progress) << "TCA DLR hybrid " << phase
                      << " progress: global_completed=" << global_completed
                      << "/" << total_groups << " ("
                      << 100.0*global_completed/
                           std::max<std::uint64_t>(1,total_groups)
                      << "%) rank0_completed=" << summary.local_groups
                      << " elapsed_seconds=" << elapsed
                      << " estimated_remaining_seconds="
                      << (total_groups-global_completed)/rate << std::endl;
          last_report=now;
        }
      }
    }
    if (global_completed!=total_groups)
      throw std::logic_error("hybrid TCA scheduler lost work units");
  } else {
    std::uint64_t newly_completed=summary.local_groups;
    for (;;) {
      MPI_Send(&newly_completed,1,MPI_UINT64_T,0,request_tag,communicator);
      std::uint64_t assignment=stop_task;
      MPI_Recv(&assignment,1,MPI_UINT64_T,0,assignment_tag,communicator,
               MPI_STATUS_IGNORE);
      if (assignment==stop_task) break;
      newly_completed=run_one(assignment);
    }
  }
  summary.phase_seconds=std::chrono::duration<double>(
    std::chrono::steady_clock::now()-begin).count();
  TCAReportDynamicBalance(progress,phase,summary.local_groups,
                          summary.local_compute_seconds,
                          summary.phase_seconds,mpi_rank,mpi_size,
                          communicator);
  return summary;
}
#endif

// Production driver for the root-summed pseudo-particle self-energy.  A work
// unit is (topology,cut,split-line directions,two split baths); the third bath
// and, when it is not canonical line zero, both directions of that root are
// contracted inside SummedRootSelfEnergyEvaluator.  There are 32*Nbath^2
// such units instead of 8*Nbath^3 fixed plans.  Units remain independent.
// MPI production runs overdecompose each structural unit into pole-pair
// tiles.  Every tile retains the complete bath-pair sum, so pole-dependent
// DLR preparation is never replicated merely to create more MPI work.  A
// balanced bulk plus dynamically dispatched tail distributes those tiles;
// the deterministic mixed-static scheduler remains as a reference.
class FullSummedRootSelfEnergyEvaluator {
public:
  static SelfEnergyNodeResult evaluate_serial(
      const Cix& cix,
      const imag_dlr::Workspace& workspace,
      const std::vector<double>& pseudo_degeneracy,
      bool keep_individual_cuts=false,
      std::ostream* progress=nullptr,
      double progress_seconds=0.0,
      bool topology_b_c3=false)
  {
    return Evaluate(cix,workspace,pseudo_degeneracy,0,1,
                    keep_individual_cuts,progress,progress_seconds,
                    topology_b_c3,false,0
#ifdef _MPI
                    ,MPI_COMM_NULL
#endif
                    );
  }

#ifdef _MPI
  static SelfEnergyNodeResult evaluate_mpi(
      const Cix& cix,
      const imag_dlr::Workspace& workspace,
      const std::vector<double>& pseudo_degeneracy,
      MPI_Comm communicator,
      bool keep_individual_cuts=false,
      std::ostream* progress=nullptr,
      double progress_seconds=0.0,
      bool topology_b_c3=false,
      bool dynamic_scheduler=true,
      int requested_micro_batch=0)
  {
    int rank=0,size=1;
    MPI_Comm_rank(communicator,&rank);
    MPI_Comm_size(communicator,&size);
    if (size==1)
      return evaluate_serial(cix,workspace,pseudo_degeneracy,
                             keep_individual_cuts,progress,progress_seconds,
                             topology_b_c3);
    return Evaluate(cix,workspace,pseudo_degeneracy,rank,size,
                    keep_individual_cuts,progress,progress_seconds,
                    topology_b_c3,dynamic_scheduler,
                    requested_micro_batch,
                    communicator);
  }
#endif

private:
  static SelfEnergyNodeResult MakeResult(
      const imag_dlr::Workspace& workspace,bool keep_individual_cuts)
  {
    SelfEnergyNodeResult result;
    result.sigma=MatrixXd::Zero(workspace.rank(),workspace.pseudo_count());
    for (auto& x:result.sigma_by_topology)
      x=MatrixXd::Zero(workspace.rank(),workspace.pseudo_count());
    if (keep_individual_cuts) {
      for (auto& x:result.cuts)
        x=MatrixXd::Zero(workspace.rank(),workspace.pseudo_count());
      for (auto& topology:result.cuts_by_topology)
        for (auto& x:topology)
          x=MatrixXd::Zero(workspace.rank(),workspace.pseudo_count());
    }
    return result;
  }

  static SelfEnergyNodeResult Evaluate(
      const Cix& cix,
      const imag_dlr::Workspace& workspace,
      const std::vector<double>& pseudo_degeneracy,
      int mpi_rank,int mpi_size,
      bool keep_individual_cuts,
      std::ostream* progress,double progress_seconds,
      bool topology_b_c3,bool dynamic_scheduler,
      int requested_micro_batch
#ifdef _MPI
      ,MPI_Comm communicator
#endif
      )
  {
    SelfEnergyNodeResult result=MakeResult(workspace,keep_individual_cuts);
    // Reference mode retains all six pseudo cuts in the fixed-arrow d0=0
    // convention.  The optional topology-B C3 candidate instead restores the
    // complete eight-direction domain for representative cuts zero and one.
    // That domain is closed under the two-vertex C3 rotation.  Multiplication
    // by three restores both three-cut orbits, while division by two removes
    // the global arrow-reversal duplication.  Thus the candidate multiplier
    // is 3/2, not the incorrect factor three previously applied inside d0=0.
    const std::uint64_t total_candidates=
      (topology_b_c3 ? 24ULL : 32ULL)*
      static_cast<std::uint64_t>(cix.Nbath)*
      static_cast<std::uint64_t>(cix.Nbath);
#ifdef _MPI
    if (mpi_size>1 && dynamic_scheduler)
      return EvaluateDynamic(cix,workspace,pseudo_degeneracy,mpi_rank,
                             mpi_size,keep_individual_cuts,progress,
                             progress_seconds,topology_b_c3,
                             requested_micro_batch,communicator);
#else
    (void)dynamic_scheduler;
    (void)requested_micro_batch;
#endif
    const std::uint64_t assigned=TCAAssignedCandidateCount(
      total_candidates,mpi_rank,mpi_size);
    // These expansions are identical for every topology, cut, direction and
    // pole tile in this self-energy evaluation.  Construct them once per
    // self-consistency iteration instead of once per scheduler work unit.
    const auto shared_pseudo_blocks=
      SummedRootSelfEnergyEvaluator::PreparePseudoBlocks(cix,workspace);
    std::vector<SummedRootSelfEnergyEvaluator::OperatorTranspose>
      shared_operator_cache;
    shared_operator_cache.reserve(2*cix.Nbath);
    for (int op=0;op<2*cix.Nbath;++op)
      shared_operator_cache.push_back(
        SummedRootSelfEnergyEvaluator::MakeOperatorTranspose(cix,op));
    const auto begin=std::chrono::steady_clock::now();
    auto last_report=begin;
    std::uint64_t candidate=0,completed=0;

    if (mpi_rank==0 && progress && progress_seconds>0.0)
      (*progress) << "TCA DLR root-summed self start: total_groups="
                  << total_candidates << " ranks=" << mpi_size
                  << " assigned_to_rank0="
                  << TCAAssignedCandidateCount(total_candidates,0,mpi_size)
                  << " scaling=" << (topology_b_c3 ? 24 : 32)
                  << "*Nbath^2 topology_B="
                  << (topology_b_c3 ?
                        "C3-all-directions" : "all-six-cuts")
                  << " evaluation=pole-major scheduler=mixed-static"
                  << std::endl;

    for (int topology=0;topology<2;++topology)
      for (int cut=0;
           cut<(topology_b_c3 && topology==1 ? 2 : 6);++cut)
        for (int d0=0;
             d0<(topology_b_c3 && topology==1 ? 2 : 1);++d0)
          for (int d1=0;d1<2;++d1)
            for (int d2=0;d2<2;++d2) {
            const std::array<int,3> direction{{d0,d1,d2}};
            const RootedCut rooted=MakeRootedCut(topology,direction,cut);
            const int root_line=rooted.root.original_line;
            // Direction zero of a noncanonical root is the representative
            // work unit.  evaluate_cut_nodes() sums original root directions
            // zero and one locally, so retaining the direction-one unit here
            // would evaluate the identical pole-pair backbone twice.
            if ((root_line==0 && d0!=0) ||
                (root_line==1 && d1!=0) ||
                (root_line==2 && d2!=0))
              continue;
            std::array<int,2> split_line{{-1,-1}};
            int ns=0;
            for (int line=0;line<3;++line)
              if (line!=root_line) split_line[ns++]=line;
            std::vector<std::array<int,2>> owned_bath_pairs;
            owned_bath_pairs.reserve(
              (cix.Nbath*cix.Nbath+mpi_size-1)/mpi_size);
            for (int bath0=0;bath0<cix.Nbath;++bath0)
              for (int bath1=0;bath1<cix.Nbath;
                   ++bath1,++candidate) {
                ++result.candidate_plans;
                if (TCACandidateOwner(candidate,mpi_size)!=mpi_rank) continue;
                owned_bath_pairs.push_back({{bath0,bath1}});
              }
            if (owned_bath_pairs.empty()) continue;

            const auto unit_begin=std::chrono::steady_clock::now();
            const auto batch=
              SummedRootSelfEnergyEvaluator::evaluate_cut_bath_pairs_nodes(
                cix,workspace,pseudo_degeneracy,topology,cut,direction,
                owned_bath_pairs,true,0,-1,&shared_pseudo_blocks,
                &shared_operator_cache);
            const double orbit_multiplicity=
              topology_b_c3 && topology==1 ? 1.5 : 1.0;
            result.sigma+=orbit_multiplicity*batch.nodes;
            result.sigma_by_topology[topology]+=
              orbit_multiplicity*batch.nodes;
            if (keep_individual_cuts) {
              result.cuts[cut]+=orbit_multiplicity*batch.nodes;
              result.cuts_by_topology[topology][cut]+=
                orbit_multiplicity*batch.nodes;
            }
            result.nonempty_plans_by_topology[topology]+=
              batch.nonempty_pairs;
            completed+=owned_bath_pairs.size();

            if (mpi_rank==0 && progress && progress_seconds>0.0) {
              const auto now=std::chrono::steady_clock::now();
              const double since=
                std::chrono::duration<double>(now-last_report).count();
              const std::uint64_t rank0_total=TCAAssignedCandidateCount(
                total_candidates,0,mpi_size);
              if (since>=progress_seconds || completed==rank0_total) {
                const double elapsed=
                  std::chrono::duration<double>(now-begin).count();
                const double rate=completed/std::max(elapsed,1.e-12);
                (*progress)
                  << "TCA DLR pole-major self progress: rank0_completed="
                  << completed << "/" << rank0_total << " ("
                  << 100.0*completed/std::max<std::uint64_t>(1,rank0_total)
                  << "%) elapsed_seconds=" << elapsed
                  << " estimated_remaining_seconds="
                  << (rank0_total-completed)/rate
                  << " last_tile=[topology=" << topology
                  << ",cut=" << cut << ",root_line=" << root_line
                  << ",split_lines=" << split_line[0] << ','
                  << split_line[1] << ",bath_pairs="
                  << owned_bath_pairs.size() << ",directions=" << d1 << ','
                  << d2 << "] last_seconds="
                  << std::chrono::duration<double>(now-unit_begin).count()
                  << std::endl;
                last_report=now;
              }
            }
          }

#ifdef _MPI
    if (mpi_size>1) {
      Allreduce(result.sigma,communicator);
      for (auto& x:result.sigma_by_topology) Allreduce(x,communicator);
      if (keep_individual_cuts) {
        for (auto& x:result.cuts) Allreduce(x,communicator);
        for (auto& topology:result.cuts_by_topology)
          for (auto& x:topology) Allreduce(x,communicator);
      }
      std::uint64_t local[3]={result.nonempty_plans_by_topology[0],
                              result.nonempty_plans_by_topology[1],completed};
      std::uint64_t global[3]={};
      MPI_Allreduce(local,global,3,MPI_UINT64_T,MPI_SUM,communicator);
      result.nonempty_plans_by_topology={{global[0],global[1]}};
      result.candidate_plans=total_candidates;
    }
#endif
    (void)assigned;
    return result;
  }

#ifdef _MPI
  struct DynamicStructure {
    int topology=0;
    int cut=0;
    int root_line=0;
    std::array<int,2> split_line{{-1,-1}};
    std::array<int,3> direction{{0,0,0}};
  };

  static SelfEnergyNodeResult EvaluateDynamic(
      const Cix& cix,
      const imag_dlr::Workspace& workspace,
      const std::vector<double>& pseudo_degeneracy,
      int mpi_rank,int mpi_size,
      bool keep_individual_cuts,
      std::ostream* progress,double progress_seconds,
      bool topology_b_c3,int requested_micro_batch,
      MPI_Comm communicator)
  {
    SelfEnergyNodeResult result=MakeResult(workspace,keep_individual_cuts);
    std::vector<DynamicStructure> structures;
    structures.reserve(topology_b_c3 ? 24 : 32);
    for (int topology=0;topology<2;++topology)
      for (int cut=0;
           cut<(topology_b_c3 && topology==1 ? 2 : 6);++cut)
        for (int d0=0;
             d0<(topology_b_c3 && topology==1 ? 2 : 1);++d0)
          for (int d1=0;d1<2;++d1)
            for (int d2=0;d2<2;++d2) {
              const std::array<int,3> direction{{d0,d1,d2}};
              const RootedCut rooted=MakeRootedCut(topology,direction,cut);
              const int root_line=rooted.root.original_line;
              if ((root_line==0 && d0!=0) ||
                  (root_line==1 && d1!=0) ||
                  (root_line==2 && d2!=0))
                continue;
              DynamicStructure item;
              item.topology=topology;
              item.cut=cut;
              item.root_line=root_line;
              item.direction=direction;
              int ns=0;
              for (int line=0;line<3;++line)
                if (line!=root_line) item.split_line[ns++]=line;
              structures.push_back(item);
            }

    const int pairs_per_structure=cix.Nbath*cix.Nbath;
    const std::uint64_t total_groups=
      static_cast<std::uint64_t>(structures.size())*pairs_per_structure;
    const std::uint64_t expected_groups=
      (topology_b_c3 ? 24ULL : 32ULL)*
      static_cast<std::uint64_t>(cix.Nbath)*cix.Nbath;
    if (total_groups!=expected_groups)
      throw std::logic_error("TCA dynamic self structure count mismatch");
    const int pole_count=workspace.pole_count();
    const int pole_pairs_per_structure=pole_count*pole_count;
    const std::uint64_t total_pole_bath_units=
      total_groups*static_cast<std::uint64_t>(pole_pairs_per_structure);
    // A task owns a tile of pole pairs and all bath pairs of one structural
    // diagram.  This order is essential: the kernel-dressed pseudo
    // coefficients are prepared once per pole pair and then reused by the
    // complete bath sum.  Bath-pair micro-batching repeated that preparation
    // once per tile (16 times for the 12-bath SC self-energy case).
    const int micro_batch=TCAChooseMicroBatchSize(
      static_cast<std::uint64_t>(structures.size())*
        pole_pairs_per_structure,
      pole_pairs_per_structure,mpi_size,requested_micro_batch,cix.Nbath);
    const std::uint64_t chunks_per_structure=
      (pole_pairs_per_structure+micro_batch-1)/micro_batch;
    const std::uint64_t total_tasks=
      static_cast<std::uint64_t>(structures.size())*chunks_per_structure;
    std::vector<std::array<int,2>> bath_pairs;
    bath_pairs.reserve(pairs_per_structure);
    for (int pair=0;pair<pairs_per_structure;++pair)
      bath_pairs.push_back({{pair/cix.Nbath,pair%cix.Nbath}});

    // Scheduler tiles differ only in pole range and cut geometry.  Pseudo
    // block scattering and the 2*Nbath one-body operator transposes are
    // immutable throughout the complete self-energy call, so share them
    // across every dynamically assigned tile on this MPI rank.
    const auto shared_pseudo_blocks=
      SummedRootSelfEnergyEvaluator::PreparePseudoBlocks(cix,workspace);
    std::vector<SummedRootSelfEnergyEvaluator::OperatorTranspose>
      shared_operator_cache;
    shared_operator_cache.reserve(2*cix.Nbath);
    for (int op=0;op<2*cix.Nbath;++op)
      shared_operator_cache.push_back(
        SummedRootSelfEnergyEvaluator::MakeOperatorTranspose(cix,op));

    if (mpi_rank==0 && progress && progress_seconds>0.0)
      (*progress) << "TCA DLR root-summed self start: total_groups="
                  << total_groups << " pole_pair_tiles=" << total_tasks
                  << " pole_pairs_per_tile=" << micro_batch
                  << " pole_bath_units=" << total_pole_bath_units
                  << " ranks=" << mpi_size
                  << " scaling=" << (topology_b_c3 ? 24 : 32)
                  << "*Nbath^2 topology_B="
                  << (topology_b_c3 ?
                        "C3-all-directions" : "all-six-cuts")
                  << " evaluation=pole-major scheduler=hybrid-dynamic"
                  << std::endl;

    const auto evaluate_task=[&](std::uint64_t task) {
      const std::uint64_t structure_index=task/chunks_per_structure;
      const std::uint64_t chunk=task%chunks_per_structure;
      const auto& item=structures[structure_index];
      const int first=static_cast<int>(chunk)*micro_batch;
      const int last=std::min(first+micro_batch,pole_pairs_per_structure);

      const auto batch=
        SummedRootSelfEnergyEvaluator::evaluate_cut_bath_pairs_nodes(
          cix,workspace,pseudo_degeneracy,item.topology,item.cut,
          item.direction,bath_pairs,true,first,last,&shared_pseudo_blocks,
          &shared_operator_cache);
      const double orbit_multiplicity=
        topology_b_c3 && item.topology==1 ? 1.5 : 1.0;
      result.sigma+=orbit_multiplicity*batch.nodes;
      result.sigma_by_topology[item.topology]+=
        orbit_multiplicity*batch.nodes;
      if (keep_individual_cuts) {
        result.cuts[item.cut]+=orbit_multiplicity*batch.nodes;
        result.cuts_by_topology[item.topology][item.cut]+=
          orbit_multiplicity*batch.nodes;
      }
      // Connectivity is independent of the pole value.  Count every bath
      // plan only in the first pole tile so the diagnostic remains a count
      // of structural bath pairs rather than pole-resolved work units.
      if (first==0)
        result.nonempty_plans_by_topology[item.topology]+=
          batch.nonempty_pairs;
      return static_cast<std::uint64_t>(last-first)*bath_pairs.size();
    };
    RunTCAHybridSchedule(total_tasks,total_pole_bath_units,mpi_rank,mpi_size,
                         communicator,progress,progress_seconds,"self",
                         evaluate_task);

    Allreduce(result.sigma,communicator);
    for (auto& x:result.sigma_by_topology) Allreduce(x,communicator);
    if (keep_individual_cuts) {
      for (auto& x:result.cuts) Allreduce(x,communicator);
      for (auto& topology:result.cuts_by_topology)
        for (auto& x:topology) Allreduce(x,communicator);
    }
    std::uint64_t local[2]={result.nonempty_plans_by_topology[0],
                            result.nonempty_plans_by_topology[1]};
    std::uint64_t global[2]={};
    MPI_Allreduce(local,global,2,MPI_UINT64_T,MPI_SUM,communicator);
    result.nonempty_plans_by_topology={{global[0],global[1]}};
    result.candidate_plans=total_groups;
    return result;
  }
#endif

#ifdef _MPI
  static void Allreduce(MatrixXd& matrix,MPI_Comm communicator)
  {
    MPI_Allreduce(MPI_IN_PLACE,matrix.data(),static_cast<int>(matrix.size()),
                  MPI_DOUBLE,MPI_SUM,communicator);
  }
#endif
};

// MPI driver for the bath-summed physical G/F cuts.  Topology A retains its
// three cut lines (12*Nbath^2 work groups), while the validated C3 reduction
// keeps one of topology B's three equivalent lines (4*Nbath^2), for a total
// of 16*Nbath^2 groups.  Dynamic MPI uses the same complete-bath pole-pair
// tiling and balanced-bulk plus dispatched-tail schedule as the self-energy
// driver.
class FullSummedPhysicalEvaluator {
public:
  static SelfEnergyNodeResult evaluate_serial(
      const Cix& cix,const imag_dlr::Workspace& workspace,
      std::ostream* progress=nullptr,double progress_seconds=0.0)
  {
    return Evaluate(cix,workspace,0,1,progress,progress_seconds,false,0
#ifdef _MPI
                    ,MPI_COMM_NULL
#endif
                    );
  }

#ifdef _MPI
  static SelfEnergyNodeResult evaluate_mpi(
      const Cix& cix,const imag_dlr::Workspace& workspace,
      MPI_Comm communicator,std::ostream* progress=nullptr,
      double progress_seconds=0.0,bool dynamic_scheduler=true,
      int requested_micro_batch=0)
  {
    int rank=0,size=1;
    MPI_Comm_rank(communicator,&rank);
    MPI_Comm_size(communicator,&size);
    if (size==1)
      return evaluate_serial(cix,workspace,progress,progress_seconds);
    return Evaluate(cix,workspace,rank,size,progress,progress_seconds,
                    dynamic_scheduler,requested_micro_batch,
                    communicator);
  }
#endif

private:
  static SelfEnergyNodeResult Evaluate(
      const Cix& cix,const imag_dlr::Workspace& workspace,
      int mpi_rank,int mpi_size,std::ostream* progress,double progress_seconds,
      bool dynamic_scheduler,int requested_micro_batch
#ifdef _MPI
      ,MPI_Comm communicator
#endif
      )
  {
    SelfEnergyNodeResult result;
    result.physical_enabled=true;
    result.g=MatrixXd::Zero(workspace.rank(),workspace.channel_count());
    result.f00=MatrixXd::Zero(workspace.rank(),workspace.channel_count());
    result.f01=MatrixXd::Zero(workspace.rank(),workspace.channel_count());
    // The three physical cuts of topology B are one C3 orbit.  Evaluate
    // line zero only and multiply by three, canceling the diagram's 1/3
    // symmetry factor.  Topology A retains all three inequivalent lines.
    const std::uint64_t total_candidates=
      16ULL*static_cast<std::uint64_t>(cix.Nbath)*
      static_cast<std::uint64_t>(cix.Nbath);
#ifdef _MPI
    if (mpi_size>1 && dynamic_scheduler)
      return EvaluateDynamic(cix,workspace,mpi_rank,mpi_size,progress,
                             progress_seconds,requested_micro_batch,
                             communicator);
#else
    (void)dynamic_scheduler;
    (void)requested_micro_batch;
#endif
    const std::uint64_t rank0_total=
      TCAAssignedCandidateCount(total_candidates,0,mpi_size);
    const auto begin=std::chrono::steady_clock::now();
    auto last_report=begin;
    std::uint64_t candidate=0,completed=0;
    if (mpi_rank==0 && progress && progress_seconds>0.0)
      (*progress) << "TCA DLR bath-summed physical start: total_groups="
                  << total_candidates << " ranks=" << mpi_size
                  << " assigned_to_rank0=" << rank0_total
                  << " scaling=16*Nbath^2 topology_B=C3-reduced"
                  << " evaluation=pole-major scheduler=mixed-static"
                  << std::endl;

    for (int topology=0;topology<2;++topology)
      for (int cut_line=0;cut_line<(topology==1 ? 1 : 3);++cut_line)
        for (int d1=0;d1<2;++d1)
          for (int d2=0;d2<2;++d2) {
            std::array<int,2> split_line{{-1,-1}};
            int ns=0;
            for (int line=0;line<3;++line)
              if (line!=cut_line) split_line[ns++]=line;
            std::vector<std::array<int,2>> owned_bath_pairs;
            owned_bath_pairs.reserve(
              (cix.Nbath*cix.Nbath+mpi_size-1)/mpi_size);
            for (int bath0=0;bath0<cix.Nbath;++bath0)
              for (int bath1=0;bath1<cix.Nbath;
                   ++bath1,++candidate) {
                ++result.candidate_plans;
                if (TCACandidateOwner(candidate,mpi_size)!=mpi_rank) continue;
                owned_bath_pairs.push_back({{bath0,bath1}});
              }
            if (owned_bath_pairs.empty()) continue;

            const auto unit_begin=std::chrono::steady_clock::now();
            const auto batch=
              SummedPhysicalEvaluator::evaluate_line_bath_pairs_nodes(
                cix,workspace,topology,cut_line,
                std::array<int,3>{{0,d1,d2}},owned_bath_pairs);
            const double orbit_multiplicity=topology==1 ? 3.0 : 1.0;
            result.g+=orbit_multiplicity*batch.values.g;
            result.f00+=orbit_multiplicity*batch.values.f00;
            result.f01+=orbit_multiplicity*batch.values.f01;
            result.nonempty_plans_by_topology[topology]+=
              batch.nonempty_pairs;
            completed+=owned_bath_pairs.size();

            if (mpi_rank==0 && progress && progress_seconds>0.0) {
              const auto now=std::chrono::steady_clock::now();
              const double since=
                std::chrono::duration<double>(now-last_report).count();
              if (since>=progress_seconds || completed==rank0_total) {
                const double elapsed=
                  std::chrono::duration<double>(now-begin).count();
                const double rate=completed/std::max(elapsed,1.e-12);
                (*progress)
                  << "TCA DLR pole-major physical progress: "
                  << "rank0_completed=" << completed << "/" << rank0_total
                  << " (" << 100.0*completed/
                       std::max<std::uint64_t>(1,rank0_total)
                  << "%) elapsed_seconds=" << elapsed
                  << " estimated_remaining_seconds="
                  << (rank0_total-completed)/rate
                  << " last_tile=[topology=" << topology
                  << ",cut_line=" << cut_line << ",split_lines="
                  << split_line[0] << ',' << split_line[1] << ",bath_pairs="
                  << owned_bath_pairs.size() << ",directions=" << d1 << ','
                  << d2 << "] last_seconds="
                  << std::chrono::duration<double>(now-unit_begin).count()
                  << std::endl;
                last_report=now;
              }
            }
          }

#ifdef _MPI
    if (mpi_size>1) {
      Allreduce(result.g,communicator);
      Allreduce(result.f00,communicator);
      Allreduce(result.f01,communicator);
      std::uint64_t local[2]={result.nonempty_plans_by_topology[0],
                              result.nonempty_plans_by_topology[1]};
      std::uint64_t global[2]={};
      MPI_Allreduce(local,global,2,MPI_UINT64_T,MPI_SUM,communicator);
      result.nonempty_plans_by_topology={{global[0],global[1]}};
      result.candidate_plans=total_candidates;
    }
#endif
    return result;
  }

#ifdef _MPI
  struct DynamicStructure {
    int topology=0;
    int cut_line=0;
    std::array<int,2> split_line{{-1,-1}};
    std::array<int,3> direction{{0,0,0}};
  };

  static SelfEnergyNodeResult EvaluateDynamic(
      const Cix& cix,const imag_dlr::Workspace& workspace,
      int mpi_rank,int mpi_size,std::ostream* progress,
      double progress_seconds,int requested_micro_batch,
      MPI_Comm communicator)
  {
    SelfEnergyNodeResult result;
    result.physical_enabled=true;
    result.g=MatrixXd::Zero(workspace.rank(),workspace.channel_count());
    result.f00=MatrixXd::Zero(workspace.rank(),workspace.channel_count());
    result.f01=MatrixXd::Zero(workspace.rank(),workspace.channel_count());
    std::vector<DynamicStructure> structures;
    structures.reserve(16);
    for (int topology=0;topology<2;++topology)
      for (int cut_line=0;cut_line<(topology==1 ? 1 : 3);++cut_line)
        for (int d1=0;d1<2;++d1)
          for (int d2=0;d2<2;++d2) {
            DynamicStructure item;
            item.topology=topology;
            item.cut_line=cut_line;
            item.direction={{0,d1,d2}};
            int ns=0;
            for (int line=0;line<3;++line)
              if (line!=cut_line) item.split_line[ns++]=line;
            structures.push_back(item);
          }

    const int pairs_per_structure=cix.Nbath*cix.Nbath;
    const std::uint64_t total_groups=
      static_cast<std::uint64_t>(structures.size())*pairs_per_structure;
    const std::uint64_t expected_groups=
      16ULL*static_cast<std::uint64_t>(cix.Nbath)*cix.Nbath;
    if (total_groups!=expected_groups)
      throw std::logic_error("TCA dynamic physical structure count mismatch");
    const int pole_count=workspace.pole_count();
    const int pole_pairs_per_structure=pole_count*pole_count;
    const std::uint64_t total_pole_bath_units=
      total_groups*static_cast<std::uint64_t>(pole_pairs_per_structure);
    const int micro_batch=TCAChooseMicroBatchSize(
      static_cast<std::uint64_t>(structures.size())*
        pole_pairs_per_structure,
      pole_pairs_per_structure,mpi_size,requested_micro_batch,cix.Nbath);
    const std::uint64_t chunks_per_structure=
      (pole_pairs_per_structure+micro_batch-1)/micro_batch;
    const std::uint64_t total_tasks=
      static_cast<std::uint64_t>(structures.size())*chunks_per_structure;
    std::vector<std::array<int,2>> bath_pairs;
    bath_pairs.reserve(pairs_per_structure);
    for (int pair=0;pair<pairs_per_structure;++pair)
      bath_pairs.push_back({{pair/cix.Nbath,pair%cix.Nbath}});

    if (mpi_rank==0 && progress && progress_seconds>0.0)
      (*progress) << "TCA DLR bath-summed physical start: total_groups="
                  << total_groups << " pole_pair_tiles=" << total_tasks
                  << " pole_pairs_per_tile=" << micro_batch
                  << " pole_bath_units=" << total_pole_bath_units
                  << " ranks=" << mpi_size
                  << " scaling=16*Nbath^2 topology_B=C3-reduced"
                  << " evaluation=pole-major scheduler=hybrid-dynamic"
                  << std::endl;

    const auto evaluate_task=[&](std::uint64_t task) {
      const std::uint64_t structure_index=task/chunks_per_structure;
      const std::uint64_t chunk=task%chunks_per_structure;
      const auto& item=structures[structure_index];
      const int first=static_cast<int>(chunk)*micro_batch;
      const int last=std::min(first+micro_batch,pole_pairs_per_structure);

      const auto batch=
        SummedPhysicalEvaluator::evaluate_line_bath_pairs_nodes(
          cix,workspace,item.topology,item.cut_line,item.direction,bath_pairs,
          -1,first,last);
      const double orbit_multiplicity=item.topology==1 ? 3.0 : 1.0;
      result.g+=orbit_multiplicity*batch.values.g;
      result.f00+=orbit_multiplicity*batch.values.f00;
      result.f01+=orbit_multiplicity*batch.values.f01;
      if (first==0)
        result.nonempty_plans_by_topology[item.topology]+=
          batch.nonempty_pairs;
      return static_cast<std::uint64_t>(last-first)*bath_pairs.size();
    };
    RunTCAHybridSchedule(total_tasks,total_pole_bath_units,
                         mpi_rank,mpi_size,
                         communicator,progress,progress_seconds,"physical",
                         evaluate_task);

    Allreduce(result.g,communicator);
    Allreduce(result.f00,communicator);
    Allreduce(result.f01,communicator);
    std::uint64_t local[2]={result.nonempty_plans_by_topology[0],
                            result.nonempty_plans_by_topology[1]};
    std::uint64_t global[2]={};
    MPI_Allreduce(local,global,2,MPI_UINT64_T,MPI_SUM,communicator);
    result.nonempty_plans_by_topology={{global[0],global[1]}};
    result.candidate_plans=total_groups;
    return result;
  }
#endif

#ifdef _MPI
  static void Allreduce(MatrixXd& matrix,MPI_Comm communicator)
  {
    MPI_Allreduce(MPI_IN_PLACE,matrix.data(),static_cast<int>(matrix.size()),
                  MPI_DOUBLE,MPI_SUM,communicator);
  }
#endif
};

class FullSelfEnergyEvaluator {
public:
  static SelfEnergyNodeResult evaluate_serial(
      const Cix& cix,
      const imag_dlr::Workspace& workspace,
      const std::vector<double>& pseudo_degeneracy,
      bool keep_individual_cuts=false,
      bool compute_physical=false,
      std::ostream* progress=nullptr,
      double progress_seconds=0.0,
      bool compute_self=true)
  {
    Diag2Matrix diag2matrix(3,cix);
    SelfEnergyNodeResult result=MakeResult(
      workspace,keep_individual_cuts,compute_physical);
    const std::uint64_t total_candidates=
      8ULL*static_cast<std::uint64_t>(cix.Nbath)*
      static_cast<std::uint64_t>(cix.Nbath)*
      static_cast<std::uint64_t>(cix.Nbath);
    const auto progress_begin=std::chrono::steady_clock::now();
    auto last_progress=progress_begin;
    const auto report_progress=[&](std::uint64_t completed,
                                   int topology,
                                   int root_orbital,
                                   int first_orbital,
                                   int second_orbital,
                                   int first_direction,
                                   int second_direction,
                                   const TCAPlanProfile& profile) {
      if (!progress || !(progress_seconds>0.0) || completed==0) return;
      const auto now=std::chrono::steady_clock::now();
      const double since_last=
        std::chrono::duration<double>(now-last_progress).count();
      if (since_last<progress_seconds && completed<total_candidates) return;
      const double elapsed=
        std::chrono::duration<double>(now-progress_begin).count();
      const double rate=completed/std::max(elapsed,1.e-12);
      const double remaining=(total_candidates-completed)/rate;
      (*progress) << "TCA DLR progress: completed=" << completed << "/"
                  << total_candidates << " ("
                  << 100.0*completed/std::max<std::uint64_t>(1,total_candidates)
                  << "%) elapsed_seconds=" << elapsed
                  << " rate_plans_per_second=" << rate
                  << " estimated_remaining_seconds=" << remaining
                  << " last_plan=[topology=" << topology
                  << ",baths=" << root_orbital << "," << first_orbital
                  << "," << second_orbital << ",directions="
                  << first_direction << "," << second_direction << "]"
                  << " last_seconds=[total=" << profile.total_seconds
                  << ",compile=" << profile.compile_seconds
                  << ",self_prepare=" << profile.self_prepare_seconds
                  << ",self_cuts=";
      for (int cut=0;cut<6;++cut)
        (*progress) << (cut ? "," : "") << profile.self_cut_seconds[cut];
      (*progress) << ",physical_prepare="
                  << profile.physical_prepare_seconds
                  << ",physical_lines=";
      for (int line=0;line<3;++line)
        (*progress) << (line ? "," : "")
                    << profile.physical_line_seconds[line];
      (*progress) << "]"
                  << std::endl;
      last_progress=now;
    };
    if (progress && progress_seconds>0.0)
      (*progress) << "TCA DLR start: total_candidates=" << total_candidates
                  << " ranks=1 assigned_to_rank0=" << total_candidates
                  << " physical=" << (compute_physical ? "yes" : "no")
                  << " progress_seconds=" << progress_seconds << std::endl;
    std::uint64_t candidate=0;
    for (int topology=0;topology<2;++topology)
      for (int root_orbital=0;root_orbital<cix.Nbath;++root_orbital)
        for (int first_orbital=0;first_orbital<cix.Nbath;++first_orbital)
          for (int second_orbital=0;second_orbital<cix.Nbath;
               ++second_orbital)
            for (int first_direction=0;first_direction<2;
                 ++first_direction)
              for (int second_direction=0;second_direction<2;
                   ++second_direction) {
                ++result.candidate_plans;
                ++candidate;
                const std::array<int,3> dirspin{{
                  2*root_orbital,
                  2*first_orbital+first_direction,
                  2*second_orbital+second_direction
                }};
                TCAPlanProfile profile;
                const auto plan_begin=std::chrono::steady_clock::now();
                const auto compile_begin=plan_begin;
                const auto plan=tca_matrix::AtomicTracePlan::Compile(
                  cix,diag2matrix,topology,dirspin);
                profile.compile_seconds=std::chrono::duration<double>(
                  std::chrono::steady_clock::now()-compile_begin).count();
                if (!plan.empty()) {
                  ++result.nonempty_plans_by_topology[topology];
                  if (compute_self) {
                    const auto cut_nodes=
                      RootedSelfEnergyEvaluator::evaluate_all_cut_nodes(
                        plan,cix,workspace,pseudo_degeneracy,
                        &profile.self_cut_seconds,
                        &profile.self_prepare_seconds);
                    for (int cut=0;cut<6;++cut) {
                      result.sigma+=cut_nodes[cut];
                      result.sigma_by_topology[topology]+=cut_nodes[cut];
                      if (keep_individual_cuts) {
                        result.cuts[cut]+=cut_nodes[cut];
                        result.cuts_by_topology[topology][cut]+=
                          cut_nodes[cut];
                      }
                    }
                  }
                  if (compute_physical) {
                    const auto physical=
                      RootedPhysicalEvaluator::evaluate_all_lines_nodes(
                        plan,cix,workspace,&profile.physical_line_seconds,
                        &profile.physical_prepare_seconds);
                    result.g+=physical.g;
                    result.f00+=physical.f00;
                    result.f01+=physical.f01;
                  }
                }
                profile.total_seconds=std::chrono::duration<double>(
                  std::chrono::steady_clock::now()-plan_begin).count();
                report_progress(candidate,topology,root_orbital,
                                first_orbital,second_orbital,
                                first_direction,second_direction,profile);
              }
    return result;
  }

#ifdef _MPI
  // Plans are the coarse MPI work unit.  For f shells and clusters there are
  // 8*Nbath^3 independent candidates, so this remains well populated on a few
  // hundred ranks without communicating atomic matrices or pole-pair data.
  // Only the compact r-by-(Np+3Nb) node arrays are summed at the end.
  static SelfEnergyNodeResult evaluate_mpi(
      const Cix& cix,
      const imag_dlr::Workspace& workspace,
      const std::vector<double>& pseudo_degeneracy,
      MPI_Comm communicator,
      bool keep_individual_cuts=false,
      bool compute_physical=false,
      std::ostream* progress=nullptr,
      double progress_seconds=0.0,
      bool compute_self=true)
  {
    int mpi_rank=0,mpi_size=1;
    MPI_Comm_rank(communicator,&mpi_rank);
    MPI_Comm_size(communicator,&mpi_size);
    if (mpi_size==1)
      return evaluate_serial(cix,workspace,pseudo_degeneracy,
                             keep_individual_cuts,compute_physical,
                             progress,progress_seconds,compute_self);
    Diag2Matrix diag2matrix(3,cix);
    SelfEnergyNodeResult result=MakeResult(
      workspace,keep_individual_cuts,compute_physical);
    const std::uint64_t total_candidates=
      8ULL*static_cast<std::uint64_t>(cix.Nbath)*
      static_cast<std::uint64_t>(cix.Nbath)*
      static_cast<std::uint64_t>(cix.Nbath);
    const std::uint64_t assigned_to_rank0=
      TCAAssignedCandidateCount(total_candidates,0,mpi_size);
    const auto progress_begin=std::chrono::steady_clock::now();
    auto last_progress=progress_begin;
    std::uint64_t owned_completed=0;
    const auto report_progress=[&](int topology,
                                   int root_orbital,
                                   int first_orbital,
                                   int second_orbital,
                                   int first_direction,
                                   int second_direction,
                                   const TCAPlanProfile& profile) {
      if (mpi_rank!=0 || !progress || !(progress_seconds>0.0) ||
          owned_completed==0) return;
      const auto now=std::chrono::steady_clock::now();
      const double since_last=
        std::chrono::duration<double>(now-last_progress).count();
      if (since_last<progress_seconds && owned_completed<assigned_to_rank0)
        return;
      const double elapsed=
        std::chrono::duration<double>(now-progress_begin).count();
      const double rate=owned_completed/std::max(elapsed,1.e-12);
      const double remaining=(assigned_to_rank0-owned_completed)/rate;
      (*progress) << "TCA DLR progress: rank0_completed="
                  << owned_completed << "/" << assigned_to_rank0 << " ("
                  << 100.0*owned_completed/
                       std::max<std::uint64_t>(1,assigned_to_rank0)
                  << "%) elapsed_seconds=" << elapsed
                  << " rate_plans_per_second=" << rate
                  << " estimated_remaining_seconds=" << remaining
                  << " rank0_nonempty=["
                  << result.nonempty_plans_by_topology[0] << ","
                  << result.nonempty_plans_by_topology[1] << "]"
                  << " last_plan=[topology=" << topology
                  << ",baths=" << root_orbital << "," << first_orbital
                  << "," << second_orbital << ",directions="
                  << first_direction << "," << second_direction << "]"
                  << " last_seconds=[total=" << profile.total_seconds
                  << ",compile=" << profile.compile_seconds
                  << ",self_prepare=" << profile.self_prepare_seconds
                  << ",self_cuts=";
      for (int cut=0;cut<6;++cut)
        (*progress) << (cut ? "," : "") << profile.self_cut_seconds[cut];
      (*progress) << ",physical_prepare="
                  << profile.physical_prepare_seconds
                  << ",physical_lines=";
      for (int line=0;line<3;++line)
        (*progress) << (line ? "," : "")
                    << profile.physical_line_seconds[line];
      (*progress) << "]"
                  << std::endl;
      last_progress=now;
    };
    if (mpi_rank==0 && progress && progress_seconds>0.0)
      (*progress) << "TCA DLR start: total_candidates=" << total_candidates
                  << " ranks=" << mpi_size
                  << " assigned_to_rank0=" << assigned_to_rank0
                  << " physical=" << (compute_physical ? "yes" : "no")
                  << " progress_seconds=" << progress_seconds
                  << " scheduler=mixed-static"
                  << " eta_basis=rank0_mixed_share" << std::endl;
    std::uint64_t candidate=0;
    for (int topology=0;topology<2;++topology)
      for (int root_orbital=0;root_orbital<cix.Nbath;++root_orbital)
        for (int first_orbital=0;first_orbital<cix.Nbath;++first_orbital)
          for (int second_orbital=0;second_orbital<cix.Nbath;
               ++second_orbital)
            for (int first_direction=0;first_direction<2;++first_direction)
              for (int second_direction=0;second_direction<2;
                   ++second_direction,++candidate) {
                ++result.candidate_plans;
                if (TCACandidateOwner(candidate,mpi_size)!=mpi_rank) continue;
                const std::array<int,3> dirspin{{
                  2*root_orbital,
                  2*first_orbital+first_direction,
                  2*second_orbital+second_direction
                }};
                TCAPlanProfile profile;
                const auto plan_begin=std::chrono::steady_clock::now();
                const auto compile_begin=plan_begin;
                const auto plan=tca_matrix::AtomicTracePlan::Compile(
                  cix,diag2matrix,topology,dirspin);
                profile.compile_seconds=std::chrono::duration<double>(
                  std::chrono::steady_clock::now()-compile_begin).count();
                if (!plan.empty()) {
                  ++result.nonempty_plans_by_topology[topology];
                  if (compute_self) {
                    const auto cut_nodes=
                      RootedSelfEnergyEvaluator::evaluate_all_cut_nodes(
                        plan,cix,workspace,pseudo_degeneracy,
                        mpi_rank==0 ? &profile.self_cut_seconds : nullptr,
                        mpi_rank==0 ? &profile.self_prepare_seconds : nullptr);
                    for (int cut=0;cut<6;++cut) {
                      result.sigma+=cut_nodes[cut];
                      result.sigma_by_topology[topology]+=cut_nodes[cut];
                      if (keep_individual_cuts) {
                        result.cuts[cut]+=cut_nodes[cut];
                        result.cuts_by_topology[topology][cut]+=
                          cut_nodes[cut];
                      }
                    }
                  }
                  if (compute_physical) {
                    const auto physical=
                      RootedPhysicalEvaluator::evaluate_all_lines_nodes(
                        plan,cix,workspace,
                        mpi_rank==0 ? &profile.physical_line_seconds : nullptr,
                        mpi_rank==0 ? &profile.physical_prepare_seconds : nullptr);
                    result.g+=physical.g;
                    result.f00+=physical.f00;
                    result.f01+=physical.f01;
                  }
                }
                profile.total_seconds=std::chrono::duration<double>(
                  std::chrono::steady_clock::now()-plan_begin).count();
                ++owned_completed;
                report_progress(topology,root_orbital,first_orbital,
                                second_orbital,first_direction,
                                second_direction,profile);
              }
    Allreduce(result.sigma,communicator);
    for (auto& topology:result.sigma_by_topology)
      Allreduce(topology,communicator);
    if (compute_physical) {
      Allreduce(result.g,communicator);
      Allreduce(result.f00,communicator);
      Allreduce(result.f01,communicator);
    }
    if (keep_individual_cuts) {
      for (auto& cut:result.cuts) Allreduce(cut,communicator);
      for (auto& topology:result.cuts_by_topology)
        for (auto& cut:topology) Allreduce(cut,communicator);
    }
    std::uint64_t counts[3]={result.candidate_plans,
                             result.nonempty_plans_by_topology[0],
                             result.nonempty_plans_by_topology[1]};
    std::uint64_t global[3]={};
    MPI_Allreduce(counts,global,3,MPI_UINT64_T,MPI_SUM,communicator);
    // candidate_plans was counted on every rank only for deterministic work
    // numbering; report the actual global candidate count once.
    result.candidate_plans=global[0]/mpi_size;
    result.nonempty_plans_by_topology={{global[1],global[2]}};
    return result;
  }
#endif

private:
  static SelfEnergyNodeResult MakeResult(
      const imag_dlr::Workspace& workspace,
      bool keep_individual_cuts,
      bool compute_physical)
  {
    if (!workspace.has_pseudos())
      throw std::logic_error("full TCA evaluation requires pseudo coefficients");
    SelfEnergyNodeResult result;
    result.sigma=MatrixXd::Zero(workspace.rank(),workspace.pseudo_count());
    for (auto& topology:result.sigma_by_topology)
      topology=MatrixXd::Zero(workspace.rank(),workspace.pseudo_count());
    result.physical_enabled=compute_physical;
    if (compute_physical) {
      result.g=MatrixXd::Zero(workspace.rank(),workspace.channel_count());
      result.f00=MatrixXd::Zero(workspace.rank(),workspace.channel_count());
      result.f01=MatrixXd::Zero(workspace.rank(),workspace.channel_count());
    }
    if (keep_individual_cuts) {
      for (auto& cut:result.cuts)
        cut=MatrixXd::Zero(workspace.rank(),workspace.pseudo_count());
      for (auto& topology:result.cuts_by_topology)
        for (auto& cut:topology)
          cut=MatrixXd::Zero(workspace.rank(),workspace.pseudo_count());
    }
    return result;
  }

#ifdef _MPI
  static void Allreduce(MatrixXd& matrix,MPI_Comm communicator)
  {
    MPI_Allreduce(MPI_IN_PLACE,matrix.data(),static_cast<int>(matrix.size()),
                  MPI_DOUBLE,MPI_SUM,communicator);
  }
#endif
};

} // namespace edmft::tca_dlr

#endif
