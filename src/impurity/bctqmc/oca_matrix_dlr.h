#ifndef EDMFT_OCA_MATRIX_DLR_H
#define EDMFT_OCA_MATRIX_DLR_H

#include "tca_dlr.h"
#include "oca_dlr.h"
#include <chrono>

namespace edmft::oca_matrix_dlr {

using Engine = tca_dlr::SummedRootSelfEnergyEvaluator;
using Physical = tca_dlr::SummedPhysicalEvaluator;
using Eigen::MatrixXd;
using Eigen::VectorXd;

// Enumeration-free OCA, in the very same product basis used by TCA/FCA.
// The Phi representative is (0,2),(1,3), both directed forwards. Cutting
// each of its four pseudo edges and rotating the open chain to zero gives
// exactly the four (root direction, split direction) combinations. There is
// NO extra orientation or symmetry factor. For a physical cut, put creation
// at zero and annihilation at two: the two cuts give split directions 0,1.
// We retain all cuts explicitly here to make comparison unambiguous.
//
// Root line: (0,2), unexpanded, bath summed by ContractRoot at vertex two.
// Split line: (1,3), expanded in AAA poles, bath held fixed in each task.
// Self chain: G0 O1 * G1, then O0 [chain] O2 Delta, then * G2 O3.
// Physical arcs: (G0 O1 * G1)(T), (G2 O3 * G3)(beta-T), closed by O0,O2.
// '*' means a time-ordered, MATRIX-valued DLR convolution. It sums internal
// atomic indices without ever materializing a four-propagator path catalog.
class Evaluator {
  const Cix& cix;
  const imag_dlr::Workspace& w;
  Engine::PseudoBlocks pseudo;
  std::vector<Engine::OperatorTranspose> operators, bulla;
  std::array<std::vector<Engine::StoredRootPair>,2> roots;
  std::vector<double> degeneracy;
  VectorXd one;

  // Unlike the fixed-path TCA prototype, OCA's dense Njm enumeration can
  // use every source/target element in the permitted atomic block. Do not
  // project dense Bulla entries onto the one-body operator's support.
  Engine::OperatorTranspose full_bulla(int bath) const {
    auto result=Physical::MakeBullaTranspose(cix,bath);
    for (auto& entry:result.entries) {
      const auto nx=cix.Njm_ind(entry.source,2*bath);
      if (nx.njm_c!=2) continue;
      for (int i=0;i<entry.matrix.rows();++i)
        for (int j=0;j<entry.matrix.cols();++j)
          entry.matrix(i,j)=nx.transposed ? (*nx.njm)(j,i) : (*nx.njm)(i,j);
      entry.nonzero.clear();
      for (int i=0;i<entry.matrix.rows();++i)
        for (int j=0;j<entry.matrix.cols();++j)
          if (entry.matrix(i,j)!=0.0)
            entry.nonzero.push_back({i,j,entry.matrix(i,j)});
    }
    return result;
  }

public:
  Evaluator(const Cix& input,const imag_dlr::Workspace& workspace)
    : cix(input),w(workspace),pseudo(Engine::PreparePseudoBlocks(cix,w)),
      degeneracy(cix.Npseudo_unique),one(VectorXd::Ones(w.rank())) {
    if (!w.has_pseudos() || w.pseudo_count()!=cix.Npseudo_unique ||
        w.channel_count()!=cix.Nbath_unique_unique)
      throw std::invalid_argument("OCA matrix workspace/CIX mismatch");
    for (int p=0;p<cix.Npseudo_unique;++p) degeneracy[p]=cix.deg_pp(p);
    for (int op=0;op<2*cix.Nbath;++op)
      operators.push_back(Engine::MakeOperatorTranspose(cix,op));
    for (int bath=0;bath<cix.Nbath;++bath) {
      bulla.push_back(full_bulla(bath));
      const int s=cix.bath_sigind0(bath);
      for (int d=0;d<2;++d)
        roots[d].push_back({operators[2*bath+1-d],operators[2*bath+d],
                           std::abs(s),d,s<0 ? -1.0 : 1.0});
    }
  }

  // One self cut / pole / split-bath task. Both signs of omega use the
  // established stable TCA endpoint/segment identity; no growing exponent
  // is introduced. Signed Sindex equivalences are used by the shared engine.
  MatrixXd self_nodes(int root_direction,int split_direction,int pole,int bath) const {
    const auto decoration=tca_dlr::MakeRootedLineDecoration(
      1,3,split_direction,w.poles()(pole));
    std::array<VectorXd,4> vertex;
    std::array<VectorXd,3> edge;
    for (auto& x:vertex) x=one;
    for (auto& x:edge) x=one;
    for (int i=0;i<4;++i) {
      const int s=decoration.vertex_kernel[i];
      if (s) vertex[i]=s>0 ? w.kernel_plus().col(pole) : w.kernel_minus().col(pole);
    }
    for (int i=0;i<3;++i) {
      const int s=decoration.edge_kernel[i];
      if (s) edge[i]=s>0 ? w.kernel_plus().col(pole) : w.kernel_minus().col(pole);
    }
    const int channel=cix.bath_sigind0(bath);
    // One crossing gives -1. Each Delta contributes the code's -residue
    // convention; ContractRoot supplies the unsplit Delta's minus itself.
    const double prefactor=-(channel<0 ? -1.0 : 1.0)*tca_dlr::CodePolePrefactor(
      w.residues()(std::abs(channel),pole),
      tca_dlr::KernelZero(w,pole,decoration.denominator_kernel),decoration);
    MatrixXd result=MatrixXd::Zero(w.rank(),w.pseudo_count());
    if (prefactor==0.0) return result;
    auto chain=Engine::PseudoValues(pseudo,edge[0]);
    chain=Engine::MultiplyRight(chain,operators[2*bath+1-split_direction],vertex[1]);
    chain=Engine::ConvolveRightPseudo(chain,cix,w,edge[1]);
    chain=Engine::ContractRoot(chain,roots[root_direction],w,vertex[2]);
    chain=Engine::ConvolveRightPseudo(chain,cix,w,edge[2]);
    chain=Engine::MultiplyRight(chain,operators[2*bath+split_direction],vertex[3]);
    Engine::AccumulateCut(chain,cix,prefactor,degeneracy,result);
    return result;
  }

  Physical::Result physical_nodes(int direction,int pole,int bath) const {
    Physical::Result result{MatrixXd::Zero(w.rank(),w.channel_count()),
      MatrixXd::Zero(w.rank(),w.channel_count()),MatrixXd::Zero(w.rank(),w.channel_count())};
    // On the closed ring the negative-frequency branch is placed on the
    // complementary arc. Both resulting kernels decay (no unstable division
    // by the exponentially small K(0) of a negative effective frequency).
    const bool direct=(direction==0 && w.poles()(pole)>=0.0) ||
                      (direction==1 && w.poles()(pole)<0.0);
    const int sign=direct ? (direction==0 ? 1 : -1) : (direction==0 ? -1 : 1);
    std::array<VectorXd,4> edge;
    for (int i=0;i<4;++i) {
      const bool inside=i>=1 && i<3;
      edge[i]=(inside==direct)
        ? VectorXd(sign>0 ? w.kernel_plus().col(pole) : w.kernel_minus().col(pole)) : one;
    }
    const int s=cix.bath_sigind0(bath);
    const double prefactor=(s<0 ? -1.0 : 1.0)*w.residues()(std::abs(s),pole)/
                          tca_dlr::KernelZero(w,pole,sign);
    if (prefactor==0.0) return result;
    auto first=Engine::PseudoValues(pseudo,edge[0]);
    first=Engine::MultiplyRight(first,operators[2*bath+1-direction],one);
    first=Engine::ConvolveRightPseudo(first,cix,w,edge[1]);
    auto second=Engine::PseudoValues(pseudo,edge[2]);
    second=Engine::MultiplyRight(second,operators[2*bath+direction],one);
    second=Engine::ConvolveRightPseudo(second,cix,w,edge[3]);
    Physical::Reflect(second,w.basis());
    for (int external=0;external<cix.Nbath;++external) {
      const int channel=cix.bath_sigind0(external);
      Physical::AccumulateBoundaryEstimators(first,second,
        operators[2*external],operators[2*external+1],bulla[external],
        prefactor*(channel<0 ? -1.0 : 1.0),std::abs(channel),
        Physical::IsF01Channel(cix,external),result);
    }
    return result;
  }
};

// Only compact node outputs are reduced. There is no path-table broadcast,
// pair Allgather, or collective inside a convolution. With MPI disabled all
// ranks evaluate the full answer independently, matching the existing API.
inline oca_dlr::CombinedResult evaluate(const Cix& cix,
    const imag_dlr::Workspace& w,const blitz::Array<double,1>& tau,
    bool physical,bool distributed,int rank,int ranks,std::ostream* log) {
  using Clock=std::chrono::steady_clock;
  const auto start=Clock::now();
  Evaluator evaluator(cix,w);
  MatrixXd sigma=MatrixXd::Zero(w.rank(),w.pseudo_count());
  Physical::Result phys{MatrixXd::Zero(w.rank(),w.channel_count()),
    MatrixXd::Zero(w.rank(),w.channel_count()),MatrixXd::Zero(w.rank(),w.channel_count())};
  const int workers=distributed ? ranks : 1, worker=distributed ? rank : 0;
  const int percut=w.pole_count()*cix.Nbath;
  const int total=(physical ? 6 : 4)*percut;
  int completed=0;
  double self_seconds=0.0,physical_seconds=0.0;
  auto last=start;
  if (rank==0 && log) *log<<"OCA matrix start: tasks="<<total<<" ranks="<<workers
    <<" poles="<<w.pole_count()<<" DLR_rank="<<w.rank()<<" physical="<<physical
    <<" enumeration=none"<<std::endl;
  for (int task=worker;task<total;task+=workers) {
    const int cut=task/percut,pole=(task%percut)/cix.Nbath,bath=task%cix.Nbath;
    const auto task_start=Clock::now();
    if (cut<4) sigma+=evaluator.self_nodes(cut/2,cut%2,pole,bath);
    else {
      const auto x=evaluator.physical_nodes(cut-4,pole,bath);
      phys.g+=x.g; phys.f00+=x.f00; phys.f01+=x.f01;
    }
    ++completed;
    const auto now=Clock::now();
    (cut<4 ? self_seconds : physical_seconds)+=
      std::chrono::duration<double>(now-task_start).count();
    if (rank==0 && log && std::chrono::duration<double>(now-last).count()>=30) {
      *log<<"OCA matrix progress: rank0_tasks="<<completed<<"/"
        <<(total+workers-1)/workers<<" elapsed_seconds="
        <<std::chrono::duration<double>(now-start).count()<<std::endl;
      last=now;
    }
  }
  const auto reduce_start=Clock::now();
  double work_min[2]={self_seconds,physical_seconds};
  double work_max[2]={self_seconds,physical_seconds};
#ifdef _MPI
  if (distributed && ranks>1) {
    for (auto* matrix:{&sigma,&phys.g,&phys.f00,&phys.f01}) {
      // Chunking avoids the signed-int MPI count limit for large pseudo sets.
      for (Eigen::Index offset=0;offset<matrix->size();) {
        const int n=static_cast<int>(std::min<Eigen::Index>(matrix->size()-offset,100000000));
        MPI_Allreduce(MPI_IN_PLACE,matrix->data()+offset,n,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        offset+=n;
      }
    }
    MPI_Allreduce(MPI_IN_PLACE,work_min,2,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,work_max,2,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  }
#endif
  const double reduction_seconds=std::chrono::duration<double>(Clock::now()-reduce_start).count();
  const auto dense=[&](const MatrixXd& nodes) {
    const MatrixXd coefficients=w.basis().values_to_coefficients(nodes);
    MatrixXd output(nodes.cols(),tau.extent(0));
    for (int it=0;it<tau.extent(0);++it)
      output.col(it)=w.basis().evaluate(tau(it),coefficients).transpose();
    return output;
  };
  oca_dlr::CombinedResult result;
  result.sigma=dense(sigma);
  if (physical) result.physical={dense(phys.g),dense(phys.f00),dense(phys.f01)};
  if (rank==0 && log) *log<<"OCA matrix done: evaluate_seconds="
    <<std::chrono::duration<double>(Clock::now()-start).count()
    <<" self_work_seconds_range=["<<work_min[0]<<","<<work_max[0]<<"]"
    <<" physical_work_seconds_range=["<<work_min[1]<<","<<work_max[1]<<"]"
    <<" rank0_reduce_wait_seconds="<<reduction_seconds<<std::endl;
  return result;
}
} // namespace edmft::oca_matrix_dlr
#endif
