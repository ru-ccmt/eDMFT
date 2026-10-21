#ifndef EDMFT_OCA_DLR_BATCHED_H
#define EDMFT_OCA_DLR_BATCHED_H

#include "oca_path_spool.h"
#include "oca_batch_cache.h"
#include <cmath>
#include <memory>

namespace edmft::oca_dlr {

// Do not spool the four unused legacy inner_index fields from OCACompactPath.
// On the supported 64-bit platforms this saves 16 bytes per partial record.
struct OCAPhysicalDiskRecord {
  OCAKey key{};
  std::array<double,5> weights{};
};

using PhysicalBatchCache=OCABatchCache<OCAPhysicalDiskRecord>;

template<class T,int N>
inline void OCAHashArray(const blitz::Array<T,N>& a,uint64_t& hash) {
  for (int d=0;d<N;++d) {
    const int extent=a.extent(d);
    hash=PhysicalBatchCache::hash_bytes(&extent,sizeof(extent),hash);
  }
  for (auto it=a.begin();it!=a.end();++it) {
    const T value=*it;
    hash=PhysicalBatchCache::hash_bytes(&value,sizeof(value),hash);
  }
}

// Hash the actual loaded symmetry/operator/Bulla data, not just a filename or
// number of propagators. In particular Sindex.dat can change signs/mappings
// without changing Npseudo_unique. Numerical propagators and hybridizations
// do not enter this signature and may change freely between DMFT steps.
inline PhysicalBatchCache::Signature OCAPhysicalCacheSignature(
    const Cix& cix,const std::vector<int>& importance,int cut_sum,
    const std::string& cix_file) {
  const auto file_hash=OCAFileHash(cix_file);
  if (!file_hash) throw std::runtime_error("Cannot hash CIX for persistent OCA cache");
  uint64_t symmetry=PhysicalBatchCache::offset_basis;
  for (const auto& a:cix.Sgm_symmetry) OCAHashArray(a,symmetry);
  uint64_t operators=PhysicalBatchCache::offset_basis;
  OCAHashArray(cix.Fi,operators); OCAHashArray(cix.Fpsi,operators);
  OCAHashArray(cix.Fbc,operators); OCAHashArray(cix.Fbpsi,operators);
  OCAHashArray(cix.Fi_indx,operators);
  // inBlock is derived from Block below. Its unused zero-state slot is not
  // initialized by the CIX reader, so hashing its raw allocation would make
  // otherwise identical solver invocations spuriously miss the cache.
  OCAHashArray(cix.bath_sigind0,operators); OCAHashArray(cix.bath_sigind,operators);
  OCAHashArray(cix.bath_conj,operators); OCAHashArray(cix.bath_SC,operators);
  OCAHashArray(cix.SC_bath,operators);
  for (const auto& a:cix.Block) OCAHashArray(a,operators);
  OCAHashArray(cix.Njm_c,operators);
  // unordered_map traversal order is not a model property: hash sorted keys.
  const std::map<int,double> ratios(cix.Njm_r.begin(),cix.Njm_r.end());
  for (const auto& [key,value]:ratios) {
    operators=PhysicalBatchCache::hash_bytes(&key,sizeof(key),operators);
    operators=PhysicalBatchCache::hash_bytes(&value,sizeof(value),operators);
  }
  for (const auto& row:cix.Njm) for (const auto& a:row) OCAHashArray(a,operators);
  return {file_hash,symmetry,operators,OCAImportanceHash(importance),
          static_cast<uint64_t>(cix.Npseudo_unique),
          static_cast<uint64_t>(cix.Nbath_unique_unique),
          static_cast<uint64_t>(cut_sum),1}; // discrete convention revision
}

// Communication-light physical OCA driver. The ordinary enumerator assigns
// complete reduced bath pairs to ranks. Its bounded hook spills partial sums;
// the local spool combines complete keys, then feeds one small plan at a time
// through all poles. Only compact final node arrays cross MPI boundaries.
//
// The persistent cache contains merged pair streams, not partial spool files
// or full Eigen plans. It is read in bounded chunks and assigned afresh to the
// current MPI size. The old .self and full-plan .physical caches are unchanged.
inline CombinedResult EvaluateOCAPhysicalBatched(
    Cix& cix,const std::vector<int>& importance,int cut_sum,
    const imag_dlr::Workspace& workspace,const std::vector<double>& degeneracy,
    const blitz::Array<double,1>& tau,size_t max_paths,
    const std::string& scratch_directory,bool distributed,
    const SelfMemoryOptions& memory,MPIcls& mpi,
    const std::string& persistent_directory={},const std::string& cix_file={})
{
  using Clock=std::chrono::steady_clock;
  const auto start=Clock::now();
  auto last_report=start;
  const bool persistent=!persistent_directory.empty();
  const auto pair_of=[](const OCAKey& key) {
    const auto indices=unpack6(key);
    return (uint64_t(indices[0])<<32)|uint32_t(indices[1]);
  };
  // Only small startup/end metadata messages are added. Neither batch writes
  // nor cached batch reads perform MPI communication.
  const auto broadcast_text=[&](std::string& text) {
#ifdef _MPI
    int count=static_cast<int>(text.size());
    MPI_Bcast(&count,1,MPI_INT,mpi.master,MPI_COMM_WORLD);
    if (count<0 || count>16*1024*1024) throw std::runtime_error("Invalid OCA cache metadata size");
    text.resize(count);
    if (count) MPI_Bcast(text.data(),count,MPI_CHAR,mpi.master,MPI_COMM_WORLD);
#endif
  };
  PhysicalBatchCache::Signature signature{};
  if (persistent && mpi.rank==mpi.master)
    signature=OCAPhysicalCacheSignature(cix,importance,cut_sum,cix_file);
#ifdef _MPI
  if (persistent) MPI_Bcast(signature.data(),signature.size(),MPI_UINT64_T,mpi.master,MPI_COMM_WORLD);
#endif
  PhysicalBatchCache cache(persistent_directory,signature);
  std::string manifest,cache_reason;
  if (persistent) {
    if (mpi.rank==mpi.master) manifest=cache.load(cache_reason);
    broadcast_text(manifest);
  }
  const bool cache_hit=!manifest.empty();
  if (cache_hit) {
    cache.accept_manifest(manifest);
    if (mpi.rank==mpi.master)
      (*mpi.out)<<"OCA bounded cache hit: directory="<<persistent_directory
        <<" pairs="<<cache.entries().size()<<" unique_paths="<<cache.total_records()
        <<" batch_paths="<<max_paths<<" ranks="<<mpi.size<<std::endl;
  } else if (persistent) {
    std::string build;
    if (mpi.rank==mpi.master) {
      (*mpi.out)<<"OCA bounded cache miss: directory="<<persistent_directory
        <<" reason="<<cache_reason<<std::endl;
      build=cache.begin();
    }
    broadcast_text(build);
    if (mpi.rank!=mpi.master) cache.join(build);
  }
  using Spool=OCAPathSpool<OCAPhysicalDiskRecord,OCAKeyHash>;
  std::unique_ptr<Spool> spool;
  if (!cache_hit) spool=std::make_unique<Spool>(max_paths,scratch_directory,pair_of);
  // Every rank logs its private scratch path for recovery after SIGKILL. A
  // normal return or exception removes it automatically; a machine reboot
  // cannot execute destructors, so the printed directory can then be removed.
  (*mpi.out)<<"OCA bounded physical start: rank="<<mpi.rank
    <<" batch_paths="<<max_paths<<" poles="<<workspace.pole_count()
    <<" scratch="<<(spool ? spool->directory().string() : "not-needed")
    <<" mpi="<<(distributed ? "final-node-reduction" : "replicated")
    <<" cache="<<(cache_hit ? "persistent-hit" : persistent ? "persistent-build" : "temporary-local")<<std::endl;
  if (mpi.rank!=mpi.master && spool)
    std::cerr<<"OCA bounded scratch: rank="<<mpi.rank
             <<" directory="<<spool->directory().string()<<std::endl;
  if (!cache_hit) {
  OCABoundedEnumeration options;
  options.max_paths=max_paths;
  options.distributed=distributed;
  options.consume=[&](const std::vector<OCACompactPath>& records) {
    std::vector<OCAPhysicalDiskRecord> disk_records;
    disk_records.reserve(records.size());
    for (const auto& record:records)
      disk_records.push_back({record.key,record.weights});
    spool->append(disk_records);
    const auto now=Clock::now();
    if (std::chrono::duration<double>(now-last_report).count()>=60.0) {
      (*mpi.out)<<"OCA bounded enumeration progress: rank="<<mpi.rank
        <<" partial_records="<<spool->partial_records()
        <<" written_GiB="<<double(spool->bytes_written())/1073741824.0
        <<" elapsed_seconds="<<std::chrono::duration<double>(now-start).count()
        <<std::endl;
      last_report=now;
    }
  };
  EnumerateOCAPaths(true,cix,importance,cut_sum,false,true,
                    "distributed",1024.0,mpi,&options);
  const auto enumeration_end=Clock::now();
  (*mpi.out)<<"OCA bounded enumeration done: rank="<<mpi.rank
    <<" partial_records="<<spool->partial_records()
    <<" seconds="<<std::chrono::duration<double>(enumeration_end-start).count()
    <<std::endl;
  }
  NodeAccumulator accumulator(workspace.rank(),cix.Npseudo_unique,
                              cix.Nbath_unique_unique,true,false);
  uint64_t batches=0,paths=0;
  double plan_seconds=0.0,evaluate_seconds=0.0;
  const auto evaluate_batch=[&](std::vector<OCAPhysicalDiskRecord>& compact) {
    if (compact.empty()) return;
    const auto pair=pair_of(compact.front().key);
    for (const auto& path:compact) {
      const auto st=unpack6(path.key);
      if (pair_of(path.key)!=pair || st[0]>=cix.Nbath_unique_unique || st[1]>=cix.Nbath_unique_unique)
        throw std::runtime_error("Invalid OCA cached bath pair");
      for (int j=2;j<6;++j) if (st[j]>=cix.Npseudo_unique || st[j]<0)
        throw std::runtime_error("Invalid OCA cached pseudo index");
      for (auto w:path.weights) if (!std::isfinite(w))
        throw std::runtime_error("Nonfinite OCA cached diagram weight");
    }
    if (persistent && !cache_hit && (distributed || mpi.rank==mpi.master))
      cache.append(pair,compact);
    std::vector<OCAPathRecord> records;
    records.reserve(compact.size());
    for (const auto& path:compact) {
      const auto st=unpack6(path.key);
      OCAPathRecord record;
      record.pseudo={st[5],st[4],st[3],st[2]};
      record.s0=st[1]; record.s1=st[0]; record.weights=path.weights;
      records.push_back(record);
    }
    std::vector<OCAPhysicalDiskRecord>().swap(compact);
    paths+=records.size();
    const auto plan_start=Clock::now();
    SelfEnergyPlan self;
    PhysicalPlan physical;
    self.build(records,cix.Npseudo_unique,cix.Nbath_unique_unique,importance,cut_sum);
    physical.build(records,cix.Npseudo_unique,cix.Nbath_unique_unique,importance,cut_sum);
    if (!CombinedEvaluator::compatible(self,physical))
      throw std::logic_error("OCA batch self/physical pair catalogs disagree");
    std::vector<OCAPathRecord>().swap(records);
    const auto evaluation_start=Clock::now();
    plan_seconds+=std::chrono::duration<double>(evaluation_start-plan_start).count();
    CombinedEvaluator::accumulate_poles(self,&physical,workspace,degeneracy,
      0,workspace.pole_count(),accumulator,memory);
    const auto now=Clock::now();
    evaluate_seconds+=std::chrono::duration<double>(now-evaluation_start).count();
    ++batches;
    if (std::chrono::duration<double>(now-last_report).count()>=60.0) {
      (*mpi.out)<<"OCA bounded evaluation progress: rank="<<mpi.rank
        <<" batches="<<batches<<" unique_paths="<<paths
        <<" elapsed_seconds="<<std::chrono::duration<double>(now-start).count()
        <<std::endl;
      last_report=now;
    }
    // Both plans die here, before another batch's records are read.
  };
  if (cache_hit) {
    const auto owners=cache.owners(distributed ? mpi.size : 1);
    for (size_t i=0;i<cache.entries().size();++i)
      if (!distributed || owners[i]==mpi.rank) cache.read_pair(i,max_paths,evaluate_batch);
  } else spool->consume(evaluate_batch);
  if (persistent && !cache_hit) {
    if (distributed || mpi.rank==mpi.master) cache.finish_writer(distributed ? mpi.rank : 0);
#ifdef _MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    if (mpi.rank==mpi.master) {
      cache.publish(distributed ? mpi.size : 1);
      (*mpi.out)<<"OCA bounded cache published: directory="<<persistent_directory
        <<" pairs="<<cache.entries().size()<<" unique_paths="<<cache.total_records()
        <<" data_GiB="<<double(cache.total_records())*sizeof(OCAPhysicalDiskRecord)/1073741824.0
        <<std::endl;
    }
    // Publication errors go through the caller's MPI_Abort path, not a local
    // early return which would strand other ranks in collectives.
#ifdef _MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
  }
  const auto local_end=Clock::now();
  (*mpi.out)<<"OCA bounded local done: rank="<<mpi.rank
    <<" batches="<<batches<<" unique_paths="<<paths
    <<" peak_reduce_keys="<<(spool ? spool->peak_keys() : 0)
    <<" repartitions="<<(spool ? spool->partitions() : 0)
    <<" written_GiB="<<double(spool ? spool->bytes_written() : 0)/1073741824.0
    <<" plan_seconds="<<plan_seconds<<" evaluate_seconds="<<evaluate_seconds
    <<" local_seconds="<<std::chrono::duration<double>(local_end-start).count()
    <<std::endl;
#ifdef _MPI
  if (distributed && mpi.size>1)
    CombinedEvaluator::allreduce_accumulator(accumulator,MPI_COMM_WORLD);
#endif
  const auto reduced=Clock::now();
  // finalize needs only dimensions from the plan; an empty plan avoids keeping
  // the last batch alive just to interpolate the summed node arrays once.
  SelfEnergyPlan shape;
  shape.build({},cix.Npseudo_unique,cix.Nbath_unique_unique,importance,cut_sum);
  auto result=CombinedEvaluator::finalize(shape,workspace,accumulator,tau);
  if (mpi.rank==mpi.master)
    (*mpi.out)<<"OCA bounded physical done: reduction_and_wait_seconds="
      <<std::chrono::duration<double>(reduced-local_end).count()
      <<" total_seconds="<<std::chrono::duration<double>(Clock::now()-start).count()
      <<std::endl;
  return result;
}
} // namespace edmft::oca_dlr
#endif
