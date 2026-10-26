#ifndef EDMFT_OCA_MEASUREMENTS_H
#define EDMFT_OCA_MEASUREMENTS_H
#include "tca_dlr.h"
#include <chrono>
#include <map>

namespace edmft::oca_measurements {

// Measurement-only fusion of the established boundary traces. Operator
// equality is exact; no small matrix elements or transitions are removed.
template<class CixType,class PhysicalType>
class Boundary {
  using Engine=typename PhysicalType::Engine;
  using Op=typename Engine::OperatorTranspose;
  using Entry=typename Engine::OperatorEntry;
  using Key=typename Engine::BlockKey;
  using Function=typename Engine::SparseFunction;
  using Result=typename PhysicalType::Result;
  using Matrix=Eigen::MatrixXd;
  using Clock=std::chrono::steady_clock;
  static constexpr bool has_kugler=requires(const CixType& c,Result& r) {
    c.kugler_estimators; r.kugler;
  };
  enum : unsigned { Green=1,Bulla=2,KLeft=4,KRight=8,KBoth=16 };
  struct Observation {
    const Op *creation,*bulla,*qcreation;
    int channel,id;
    bool anomalous;
    double sign;
  };
  struct Term {
    const Entry* left;
    int channel,observation;
    bool anomalous;
    double sign;
    unsigned mask;
  };
  struct Plan {
    const Entry *right=nullptr,*qright=nullptr;
    std::array<std::map<Key,std::vector<Term>>,2> closure;
  };
  struct Group {
    int physical;
    const Op *annihilation,*qannihilation;
    std::vector<Observation> observations;
    std::map<Key,Plan> plans;
  };
  const CixType& cix;
  bool kugler=false,dense,profile;
  std::vector<Group> groups;
  struct Stats {
    std::uint64_t calls=0,plans=0,ordinary_products=0,composite_products=0,
      dense_composite_products=0,dense_closures=0,sparse_closures=0,
      fused_bulla_kugler=0;
    double prepare_seconds=0,product_seconds=0,closure_seconds=0;
  } stats;

  static bool same_entry(const Entry& a,const Entry& b) {
    return a.source==b.source && a.target==b.target &&
      a.matrix.rows()==b.matrix.rows() && a.matrix.cols()==b.matrix.cols() &&
      (a.matrix.array()==b.matrix.array()).all();
  }
  static bool same_operator(const Op& a,const Op& b) {
    if(a.entries.size()!=b.entries.size() || a.by_source!=b.by_source ||
       a.by_target!=b.by_target)return false;
    for(std::size_t i=0;i<a.entries.size();++i)
      if(!same_entry(a.entries[i],b.entries[i]))return false;
    return true;
  }
  static const Entry* source(const Op* op,int block) {
    if(!op || block<0 || block>=static_cast<int>(op->by_source.size()))return nullptr;
    int index=op->by_source[block];return index<0?nullptr:&op->entries[index];
  }
  Plan make_plan(const Key& key,const Group& group) {
    Plan p;p.right=source(group.annihilation,key.second);
    if(!p.right)return p;
    p.qright=source(group.qannihilation,key.second);
    auto add=[&](int right_kind,const Entry& left,const Observation& o,unsigned mask) {
      const Entry* right=right_kind?p.qright:p.right;if(!right)return;
      Key second{right->target,left.source};auto& terms=p.closure[right_kind][second];
      // Reuse a scalar trace only for the SAME external bath observation.
      // Repeated equivalent bath representatives still contribute separately.
      if(mask==KLeft)for(auto& t:terms)
        if(t.observation==o.id && (t.mask&Bulla) && same_entry(*t.left,left)) {
          t.mask|=KLeft;++stats.fused_bulla_kugler;return;
        }
      terms.push_back({&left,o.channel,o.id,o.anomalous,o.sign,mask});
    };
    for(const auto& o:group.observations) {
      if(key.first<0 || key.first>=static_cast<int>(o.creation->by_target.size()))continue;
      for(int ci:o.creation->by_target[key.first]) {
        const auto& left=o.creation->entries[ci];add(0,left,o,Green);
        if(const Entry* b=source(o.bulla,left.source)) {
          if(b->target!=key.first)throw std::logic_error("OCA Bulla target differs from creation target");
          add(0,*b,o,Bulla);
        }
        if(kugler)add(1,left,o,KRight);
      }
      if(kugler && o.qcreation)
        for(int qi:o.qcreation->by_target.at(key.first)) {
          const auto& left=o.qcreation->entries[qi];
          add(0,left,o,KLeft);add(1,left,o,KBoth);
        }
    }
    return p;
  }
  void right_product(Matrix& out,const Matrix& values,const Entry& op,bool composite) {
    if(values.cols()!=op.matrix.rows())throw std::logic_error("OCA measurement product dimension mismatch");
    out.resize(values.rows(),op.matrix.cols());
    if(op.signed_partial_permutation) {
      out.setZero();
      for(int i=0;i<static_cast<int>(op.target_column.size());++i) {
        int j=op.target_column[i];if(j>=0)out.col(j).noalias()+=op.row_sign[i]*values.col(i);
      }
    } else if(dense && composite && op.matrix.size()>=256 &&
              op.nonzero.size()*5>=static_cast<std::size_t>(op.matrix.size())) {
      out.noalias()=values*op.matrix;++stats.dense_composite_products;
    } else {
      out.setZero();
      for(const auto& e:op.nonzero)out.col(e.column).noalias()+=e.value*values.col(e.row);
    }
  }
  void accumulate_term(const Term& term,double value,int node,double prefactor,Result& out) {
    const double v=prefactor*term.sign*value;
    if(term.mask&Green)out.g(node,term.channel)+=v;
    if(term.mask&Bulla)(term.anomalous?out.f01:out.f00)(node,term.channel)+=v;
    if constexpr(has_kugler) {
      if(term.mask&KLeft)out.kugler.value[int(term.anomalous)](node,term.channel)+=v;
      if(term.mask&KRight)out.kugler.value[2+int(term.anomalous)](node,term.channel)+=v;
      if(term.mask&KBoth)out.kugler.value[4+int(term.anomalous)](node,term.channel)+=v;
    }
  }
public:
  Boundary(const CixType& input,const std::vector<Op>& operators,
           const std::vector<Op>& bulla,bool use_dense=true,bool profiling=false)
    :cix(input),dense(use_dense),profile(profiling) {
    auto begin=Clock::now();
    if constexpr(has_kugler)kugler=cix.kugler_estimators;
    for(int b=0;b<cix.Nbath;++b) {
      const Op *qcreation=nullptr,*qannihilation=nullptr;
      if constexpr(has_kugler)if(kugler) {
        qcreation=&PhysicalType::InteractionOperator(operators[2*b]);
        qannihilation=&PhysicalType::InteractionOperator(operators[2*b+1]);
      }
      int physical=cix.Fi_indx(2*b+1),index=-1;
      for(int g=0;g<static_cast<int>(groups.size());++g)
        if(groups[g].physical==physical && same_operator(*groups[g].annihilation,operators[2*b+1]) &&
           (!kugler || same_operator(*groups[g].qannihilation,*qannihilation))) {index=g;break;}
      if(index<0) {index=groups.size();groups.push_back({physical,&operators[2*b+1],qannihilation,{},{}});}
      int s=cix.bath_sigind0(b);
      groups[index].observations.push_back({&operators[2*b],&bulla[b],qcreation,
        std::abs(s),b,PhysicalType::IsF01Channel(cix,b),s<0?-1.:1.});
    }
    stats.prepare_seconds=std::chrono::duration<double>(Clock::now()-begin).count();
  }
  void accumulate(const Function& first,const Function& second,double prefactor,Result& out) {
    ++stats.calls;if(first.empty() || second.empty())return;
    const int rank=first.begin()->second.size();
    if constexpr(has_kugler)if(kugler)out.kugler.ensure(rank,cix.Nbath_unique_unique);
    struct Active {const std::vector<Term>* terms;const tca_dlr::MatrixFunction* second;bool gemm;};
    Matrix product,closed;
    for(const auto& [key,values]:first)for(auto& group:groups) {
      auto found=group.plans.find(key);
      if(found==group.plans.end()){found=group.plans.emplace(key,make_plan(key,group)).first;++stats.plans;}
      const auto& plan=found->second;
      for(int kind=0;kind<(kugler?2:1);++kind) {
        const Entry* right=kind?plan.qright:plan.right;if(!right)continue;
        std::vector<Active> active;
        for(const auto& [second_key,terms]:plan.closure[kind]) {
          auto match=second.find(second_key);if(match==second.end() || match->second.empty())continue;
          std::size_t nnz=0;for(const auto& t:terms)nnz+=t.left->nonzero.size();
          std::size_t area=static_cast<std::size_t>(values[0].rows())*match->second[0].cols();
          bool gemm=dense && area>=256 && nnz*4>=area && right->matrix.cols()>=8;
          active.push_back({&terms,&match->second,gemm});
        }
        if(active.empty())continue;
        for(int node=0;node<rank;++node) {
          auto start=profile?Clock::now():Clock::time_point{};
          right_product(product,values[node],*right,kind!=0);
          if(kind)++stats.composite_products;else ++stats.ordinary_products;
          if(profile)stats.product_seconds+=std::chrono::duration<double>(Clock::now()-start).count();
          start=profile?Clock::now():Clock::time_point{};
          for(const auto& closure:active) {
            const auto& reflected=(*closure.second)[node];
            if(closure.gemm) {
              closed.noalias()=product*reflected;++stats.dense_closures;
              for(const auto& term:*closure.terms)
                accumulate_term(term,Engine::TraceLeftOperator(*term.left,closed),node,prefactor,out);
            } else {
              ++stats.sparse_closures;
              for(const auto& term:*closure.terms)
                accumulate_term(term,Engine::TraceLeftProduct(*term.left,product,reflected),node,prefactor,out);
            }
          }
          if(profile)stats.closure_seconds+=std::chrono::duration<double>(Clock::now()-start).count();
        }
      }
    }
  }
  void report(std::ostream& out) const {
    out<<"OCA measurement profile: calls="<<stats.calls<<" operator_groups="<<groups.size()
       <<" plans="<<stats.plans<<" ordinary_products="<<stats.ordinary_products
       <<" composite_products="<<stats.composite_products<<" dense_composite_products="<<stats.dense_composite_products
       <<" dense_closures="<<stats.dense_closures<<" sparse_closures="<<stats.sparse_closures
       <<" fused_bulla_kugler="<<stats.fused_bulla_kugler<<" prepare_seconds="<<stats.prepare_seconds
       <<" product_seconds="<<stats.product_seconds<<" closure_seconds="<<stats.closure_seconds<<std::endl;
  }
};
}
#endif
