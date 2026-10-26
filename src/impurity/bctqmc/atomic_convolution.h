#ifndef EDMFT_ATOMIC_CONVOLUTION_H
#define EDMFT_ATOMIC_CONVOLUTION_H
#include "dlr.h"
#include <blitz/array.h>
#include <deque>
#include <map>
#include <numeric>

namespace edmft::atomic_convolution {
using Eigen::MatrixXd;
using MatrixFunction=std::vector<MatrixXd>;
using Key=std::pair<int,int>;
using Function=std::map<Key,MatrixFunction>;
struct Piece {std::vector<int> indices;MatrixFunction coefficients,hilbert;bool nonzero=false;};
struct Right {int dimension=0;std::vector<Piece> pieces;};
inline Right prepare(const blitz::Array<int,2>& labels,const MatrixXd& coefficients,
                     const MatrixXd& hilbert,const std::vector<unsigned char>& active) {
  Right result;const int n=labels.extent(0),rank=coefficients.rows();result.dimension=n;
  if(labels.extent(1)!=n || hilbert.rows()!=rank || hilbert.cols()!=coefficients.cols())throw std::invalid_argument("atomic convolution right dimensions differ");
  std::vector<int> parent(n);std::iota(parent.begin(),parent.end(),0);
  const auto root=[&](int i){while(parent[i]!=i){parent[i]=parent[parent[i]];i=parent[i];}return i;};
  for(int i=0;i<n;++i)for(int j=0;j<i;++j){const int a=labels(i,j),b=labels(j,i);
    if((a && active.at(std::abs(a)-1)) || (b && active.at(std::abs(b)-1))){const int x=root(i),y=root(j);parent[std::max(x,y)]=std::min(x,y);}}
  std::vector<int> location(n,-1);
  for(int i=0;i<n;++i){int r=root(i);if(location[r]<0){location[r]=result.pieces.size();result.pieces.emplace_back();}result.pieces[location[r]].indices.push_back(i);}
  for(auto& piece:result.pieces){int d=piece.indices.size();piece.coefficients.assign(rank,MatrixXd::Zero(d,d));piece.hilbert=piece.coefficients;
    for(int i=0;i<d;++i)for(int j=0;j<d;++j){int label=labels(piece.indices[i],piece.indices[j]);if(!label)continue;int p=std::abs(label)-1;double sign=label<0?-1.:1.;
      for(int l=0;l<rank;++l){double a=sign*coefficients(l,p),b=sign*hilbert(l,p);piece.coefficients[l](i,j)=a;piece.hilbert[l](i,j)=b;piece.nonzero=piece.nonzero || a!=0. || b!=0.;}}
  }
  return result;
}
struct Input {int state;Key key;const MatrixFunction* values;const MatrixFunction* left_hilbert;const Right* right;};
// Preserve complete rectangular boundary maps. Only the right propagator is
// block diagonal; its components do not identify or truncate operator paths.
inline std::vector<Function> evaluate(const dlr::FermionicBasis& basis,
    const std::vector<Input>& input,std::size_t states,bool reflect,
    dlr::FermionicBasis::ConvolutionScratch* scratch=nullptr,
    bool input_coefficients=false,bool output_coefficients=false,std::size_t cap=16384) {
  const int rank=basis.rank();cap=std::max<std::size_t>(1,cap);std::vector<Function> result(states);
  struct Job {int state;Key key;int rows,columns;std::vector<int> row;const Piece* piece;MatrixFunction values,hilbert;bool prepared_left;};
  std::deque<Job> jobs;std::size_t flat=0;bool prepared_left=false;
  const auto flush=[&](){if(jobs.empty())return;
    std::vector<const MatrixFunction*> f,g,hg,hl;
    for(const auto& job:jobs){f.push_back(&job.values);g.push_back(&job.piece->coefficients);hg.push_back(&job.piece->hilbert);if(prepared_left)hl.push_back(&job.hilbert);}
    auto output=basis.time_ordered_convolution_values_right_prepared_batch(f,g,hg,reflect,scratch,prepared_left?&hl:nullptr,input_coefficients,output_coefficients);
    for(std::size_t k=0;k<jobs.size();++k){const auto& job=jobs[k];auto& dst=result.at(job.state)[job.key];if(dst.empty())dst.assign(rank,MatrixXd::Zero(job.rows,job.columns));
      for(int l=0;l<rank;++l)for(std::size_t i=0;i<job.row.size();++i)for(std::size_t j=0;j<job.piece->indices.size();++j)dst[l](job.row[i],job.piece->indices[j])=output[k][l](i,j);
    }jobs.clear();flat=0;
  };
  for(const auto& item:input){
    if(!item.values || !item.right || static_cast<int>(item.values->size())!=rank)throw std::invalid_argument("atomic convolution input rank differs");
    const auto& values=*item.values;const int nr=values[0].rows(),nc=item.right->dimension;bool has_left=item.left_hilbert!=nullptr;
    if(has_left && static_cast<int>(item.left_hilbert->size())!=rank)throw std::invalid_argument("atomic convolution left Hilbert rank differs");
    for(int l=0;l<rank;++l)if(values[l].rows()!=nr || values[l].cols()!=nc || (has_left && ((*item.left_hilbert)[l].rows()!=nr || (*item.left_hilbert)[l].cols()!=nc)))throw std::invalid_argument("atomic convolution matrix dimensions differ");
    if(!jobs.empty() && has_left!=prepared_left)flush();prepared_left=has_left;
    for(const auto& piece:item.right->pieces){if(!piece.nonzero)continue;std::vector<int> rows;
      for(int i=0;i<nr;++i){bool nonzero=false;for(int l=0;l<rank && !nonzero;++l)for(int j:piece.indices)if(values[l](i,j)!=0. || (has_left && (*item.left_hilbert)[l](i,j)!=0.)){nonzero=true;break;}if(nonzero)rows.push_back(i);}
      const std::size_t row_tile=std::max<std::size_t>(1,cap/piece.indices.size());
      for(std::size_t first=0;first<rows.size();first+=row_tile){std::size_t last=std::min(rows.size(),first+row_tile),size=(last-first)*piece.indices.size();if(flat && flat+size>cap)flush();
        Job job;job.state=item.state;job.key=item.key;job.rows=nr;job.columns=nc;job.row.assign(rows.begin()+first,rows.begin()+last);job.piece=&piece;job.prepared_left=has_left;
        job.values.resize(rank);if(has_left)job.hilbert.resize(rank);
        for(int l=0;l<rank;++l){job.values[l].resize(job.row.size(),piece.indices.size());if(has_left)job.hilbert[l].resize(job.row.size(),piece.indices.size());
          for(std::size_t i=0;i<job.row.size();++i)for(std::size_t j=0;j<piece.indices.size();++j){job.values[l](i,j)=values[l](job.row[i],piece.indices[j]);if(has_left)job.hilbert[l](i,j)=(*item.left_hilbert)[l](job.row[i],piece.indices[j]);}}
        jobs.push_back(std::move(job));flat+=size;if(flat>=cap)flush();
      }
    }
  }
  flush();return result;
}
}
#endif
