#ifndef EDMFT_OCA_SUBBLOCKS_H
#define EDMFT_OCA_SUBBLOCKS_H
#include "tca_dlr.h"
#include <array>
#include <deque>
#include <memory>
#include <numeric>

namespace edmft::oca_subblocks {
using Engine=tca_dlr::SummedRootSelfEnergyEvaluator;
using MatrixFunction=tca_dlr::MatrixFunction;
using Eigen::MatrixXd;
using Eigen::VectorXd;

// Exact right convolution. Pieces are rebuilt for the current pseudo workspace:
// no Sindex labels are removed and no small-value cutoff is used. Rectangular
// intermediates retain every row and every operator-induced inter-piece link.
class Convolution {
  const Cix& cix;
  const imag_dlr::Workspace& w;
  std::vector<std::vector<std::vector<int>>> parts;
  struct Piece { MatrixFunction g,hg; bool nonzero=false; };
  struct Slot {
    VectorXd scalar;
    Engine::PreparedPseudoCoefficients prepared;
    std::vector<std::vector<Piece>> blocks;
    std::vector<unsigned char> ready;
    unsigned long age=0;
  };
  mutable std::array<std::unique_ptr<Slot>,2> slots;
  mutable unsigned long clock=0;
  mutable dlr::FermionicBasis::ConvolutionScratch scratch;
  const std::size_t flat_cap;
public:
  struct Stats {
    std::uint64_t calls=0,preparations=0,tiles=0,packed_entries=0,dense_entries=0;
    std::size_t peak_batch_entries=0;
  };
  mutable Stats stats;
  explicit Convolution(const Cix& input,const imag_dlr::Workspace& workspace,
                       bool profile=false,std::size_t cap=16384,bool contiguous=true)
    :cix(input),w(workspace),parts(input.Nblock),flat_cap(std::max<std::size_t>(1,cap)) {
    scratch.profile=profile;scratch.contiguous_products=contiguous;
    std::vector<unsigned char> active(cix.Npseudo_unique,0);
    for (int p=0;p<cix.Npseudo_unique;++p)
      active[p]=(w.pseudo_values().col(p).array()!=0.0).any() ||
                (w.pseudo_coefficients().col(p).array()!=0.0).any() ||
                (w.pseudo_hilbert_coefficients().col(p).array()!=0.0).any();
    for (int b=0;b<cix.Nblock;++b) {
      const int n=cix.Block[b].size();
      std::vector<int> parent(n);std::iota(parent.begin(),parent.end(),0);
      const auto root=[&](int i){while(parent[i]!=i){parent[i]=parent[parent[i]];i=parent[i];}return i;};
      const auto& labels=cix.Sgm_symmetry[b];
      for(int i=0;i<n;++i)for(int j=0;j<i;++j){
        const int a=labels(i,j),d=labels(j,i);
        if((a && active[std::abs(a)-1]) || (d && active[std::abs(d)-1])){
          const int x=root(i),y=root(j);parent[std::max(x,y)]=std::min(x,y);
        }
      }
      std::vector<int> where(n,-1);
      for(int i=0;i<n;++i){const int r=root(i);if(where[r]<0){where[r]=parts[b].size();parts[b].emplace_back();}parts[b][where[r]].push_back(i);}
    }
  }
  const auto& partition() const {return parts;}
  void report(std::ostream& out) const {
    out<<"OCA convolution profile: calls="<<stats.calls<<" preparations="<<stats.preparations
       <<" contiguous="<<scratch.contiguous_products<<" tiles="<<stats.tiles<<" packed_entries="<<stats.packed_entries
       <<" dense_entries="<<stats.dense_entries<<" peak_batch_entries="<<stats.peak_batch_entries
       <<" pack_seconds="<<scratch.pack_seconds<<" transform_seconds="<<scratch.transform_seconds
       <<" product_seconds="<<scratch.product_seconds<<" unpack_seconds="<<scratch.unpack_seconds<<std::endl;
  }
private:
  Slot& prepare(const VectorXd& scalar) const {
    for(auto& slot:slots)if(slot && slot->scalar.size()==scalar.size() &&
      (slot->scalar.array()==scalar.array()).all()){slot->age=++clock;return *slot;}
    int k=0;if(slots[0] && (!slots[1] || slots[1]->age<slots[0]->age))k=1;
    slots[k]=std::make_unique<Slot>();auto& slot=*slots[k];slot.scalar=scalar;slot.age=++clock;
    slot.prepared=Engine::PreparePseudoCoefficients(w,scalar);
    slot.blocks.resize(cix.Nblock);slot.ready.assign(cix.Nblock,0);++stats.preparations;return slot;
  }
  void expand(Slot& slot,int block) const {
    if(slot.ready[block])return;
    const auto& unique=slot.prepared.use_workspace_coefficients ? w.pseudo_coefficients() : slot.prepared.dressed;
    auto& pieces=slot.blocks[block];pieces.resize(parts[block].size());
    for(std::size_t k=0;k<pieces.size();++k){
      const auto& ids=parts[block][k];auto& piece=pieces[k];const int n=ids.size();
      piece.g.assign(w.rank(),MatrixXd::Zero(n,n));piece.hg=piece.g;
      for(int i=0;i<n;++i)for(int j=0;j<n;++j){const int label=cix.Sgm_symmetry[block](ids[i],ids[j]);if(!label)continue;
        const int p=std::abs(label)-1;const double sign=label<0?-1.:1.;
        for(int l=0;l<w.rank();++l){
          piece.g[l](i,j)=sign*unique(l,p);piece.hg[l](i,j)=sign*slot.prepared.hilbert(l,p);
          piece.nonzero=piece.nonzero || piece.g[l](i,j)!=0.0 || piece.hg[l](i,j)!=0.0;
        }
      }
    }
    slot.ready[block]=1;
  }
public:
  Engine::SparseFunction operator()(const Engine::SparseFunction& input,
       const VectorXd& scalar,bool reflect_output=false) const {
    if(scalar.size()!=w.rank())throw std::invalid_argument("OCA subblock scalar rank mismatch");
    ++stats.calls;if(input.empty())return {};
    auto& slot=prepare(scalar);Engine::SparseFunction result;
    struct Job {Engine::BlockKey key;std::vector<int> rows;const std::vector<int>* columns;MatrixFunction values;const Piece* piece;int row_count;};
    std::deque<Job> jobs;std::size_t flat=0;
    auto flush=[&](){
      if(jobs.empty())return;
      std::vector<const MatrixFunction*> f,g,hg;f.reserve(jobs.size());g.reserve(jobs.size());hg.reserve(jobs.size());
      for(const auto& job:jobs){f.push_back(&job.values);g.push_back(&job.piece->g);hg.push_back(&job.piece->hg);}
      auto output=w.basis().time_ordered_convolution_values_right_prepared_batch(f,g,hg,reflect_output,&scratch);
      for(std::size_t k=0;k<jobs.size();++k){const auto& job=jobs[k];auto& dst=result[job.key];
        if(dst.empty())dst.assign(w.rank(),MatrixXd::Zero(job.row_count,cix.Block[job.key.second].size()));
        for(int l=0;l<w.rank();++l)for(std::size_t i=0;i<job.rows.size();++i)for(std::size_t j=0;j<job.columns->size();++j)
          dst[l](job.rows[i],(*job.columns)[j])=output[k][l](i,j);
      }
      stats.peak_batch_entries=std::max(stats.peak_batch_entries,flat);jobs.clear();flat=0;
    };
    for(const auto& [key,values]:input){
      if(static_cast<int>(values.size())!=w.rank())throw std::invalid_argument("OCA subblock input rank mismatch");
      expand(slot,key.second);const int nr=values.empty()?0:values.front().rows(),nc=cix.Block[key.second].size();
      for(const auto& matrix:values)if(matrix.rows()!=nr || matrix.cols()!=nc)throw std::invalid_argument("OCA subblock input dimensions mismatch");
      stats.dense_entries+=static_cast<std::uint64_t>(nr)*nc;
      for(std::size_t k=0;k<parts[key.second].size();++k){
        const auto& piece=slot.blocks[key.second][k];if(!piece.nonzero)continue;const auto& columns=parts[key.second][k];
        std::vector<int> rows;
        for(int i=0;i<nr;++i){bool nonzero=false;for(int l=0;l<w.rank() && !nonzero;++l)for(int j:columns)if(values[l](i,j)!=0.0){nonzero=true;break;}if(nonzero)rows.push_back(i);}
        const std::size_t row_tile=std::max<std::size_t>(1,flat_cap/columns.size());
        for(std::size_t first=0;first<rows.size();first+=row_tile){const auto last=std::min(rows.size(),first+row_tile);const std::size_t size=(last-first)*columns.size();
          if(flat && flat+size>flat_cap)flush();
          Job job;job.key=key;job.row_count=nr;job.rows.assign(rows.begin()+first,rows.begin()+last);job.columns=&columns;job.piece=&piece;
          job.values.resize(w.rank());for(int l=0;l<w.rank();++l){job.values[l].resize(job.rows.size(),columns.size());for(std::size_t i=0;i<job.rows.size();++i)for(std::size_t j=0;j<columns.size();++j)job.values[l](i,j)=values[l](job.rows[i],columns[j]);}
          jobs.push_back(std::move(job));flat+=size;++stats.tiles;stats.packed_entries+=size;
          if(flat>=flat_cap)flush();
        }
      }
    }
    flush();return result;
  }
};
}
#endif
