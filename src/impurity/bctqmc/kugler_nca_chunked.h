#ifndef EDMFT_KUGLER_NCA_CHUNKED_H
#define EDMFT_KUGLER_NCA_CHUNKED_H
#include "tca_dlr.h"

namespace edmft::kugler_nca {
// Each time point is evaluated wholly on one rank, in the same atomic/bath
// summation order as the dense implementation. Reduce only the six small
// final correlators; never replicate full Nt x N x N propagator containers.
inline kugler::Correlators evaluate(const Cix& cix,const Array<double,2>& pp,
                                    MPIcls& mpi,int chunk_size=8) {
  if(chunk_size<=0)throw std::invalid_argument("Kugler NCA chunk size must be positive");
  using Physical=tca_dlr::SummedPhysicalEvaluator;using Engine=Physical::Engine;
  const int nt=pp.extent(1);kugler::Correlators dense;dense.ensure(cix.Nbath_unique_unique,nt);
  std::vector<Engine::OperatorTranspose> creation,annihilation;
  creation.reserve(cix.Nbath);annihilation.reserve(cix.Nbath);
  for(int b=0;b<cix.Nbath;++b){creation.push_back(Engine::MakeOperatorTranspose(cix,2*b));annihilation.push_back(Engine::MakeOperatorTranspose(cix,2*b+1));}
  int worker=0,workers=1;
#ifdef _MPI
  worker=mpi.rank;workers=mpi.size;
#endif
  for(int begin=worker*chunk_size;begin<nt;begin+=workers*chunk_size){
    const int ntime=std::min(chunk_size,nt-begin);Engine::SparseFunction first,second;
    for(int b=0;b<cix.Nblock;++b){const int n=cix.Block[b].size();auto& a=first[{b,b}];auto& r=second[{b,b}];a.assign(ntime,Eigen::MatrixXd::Zero(n,n));r=a;
      for(int i=0;i<n;++i)for(int j=0;j<n;++j){const int p=cix.Sgm_symmetry[b](i,j);if(!p)continue;const double sign=p<0?-1.:1.;
        for(int t=0;t<ntime;++t){a[t](i,j)=sign*pp(std::abs(p)-1,begin+t);r[t](i,j)=sign*pp(std::abs(p)-1,nt-1-begin-t);}
      }
    }
    kugler::Correlators local;local.ensure(ntime,cix.Nbath_unique_unique);
    for(int b=0;b<cix.Nbath;++b)Physical::AccumulateKugler(first,second,creation[b],annihilation[b],cix.bath_sigind0(b)<0?-1.:1.,std::abs(cix.bath_sigind0(b)),Physical::IsF01Channel(cix,b),local);
    for(int a=0;a<6;++a)dense.value[a].middleCols(begin,ntime)=local.value[a].transpose();
  }
#ifdef _MPI
  if(workers>1)for(auto& values:dense.value){
    for(Eigen::Index offset=0;offset<values.size();){const int n=static_cast<int>(std::min<Eigen::Index>(values.size()-offset,100000000));MPI_Allreduce(MPI_IN_PLACE,values.data()+offset,n,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);offset+=n;}
  }
#endif
  return dense;
}
}
#endif
