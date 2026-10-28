#ifndef EDMFT_SPATIAL_SYMMETRY_H
#define EDMFT_SPATIAL_SYMMETRY_H
#include "time_reversal.h"
#include <functional>
namespace edmft::spatial_symmetry {
// Only signed, ordered orbital involutions are proposed. No point-group,
// spin, momentum, or Nambu convention is inferred from orbital numbering.
template<class CixType,class Workspace>
oca_symmetry::Verification discover(const CixType& c,const Workspace& w,int& attempts) {
  oca_symmetry::Verification best;best.reason="no verified ordered spatial involution";
  if(c.Npsi<2 || c.Npsi>63){best.reason="unsupported orbital count";return best;}
  using Signature=std::vector<std::array<int,4>>;
  std::vector<Signature> sig(c.Npsi);
  for(int b=0;b<c.Nbath;++b){int a=c.Fi_indx(2*b),d=c.Fi_indx(2*b+1);
    if(a<0 || d<0 || a/2>=c.Npsi || d/2>=c.Npsi){best.reason="invalid bath operator";return best;}
    if(a/2==d/2 && (a%2)!=(d%2))sig[a/2].push_back({c.bath_sigind0(b),int(c.bath_conj(b)),int(c.SC_bath(b)),a%2});
    for(int l=0;l<w.pole_count();++l)if(!std::isfinite(w.poles()(l)) || !std::isfinite(w.residues()(std::abs(c.bath_sigind0(b)),l))){best.reason="nonfinite bath";return best;}
  }
  std::map<Signature,std::vector<int>> grouped;
  for(int a=0;a<c.Npsi;++a){if(sig[a].empty()){best.reason="missing diagonal bath";return best;}std::sort(sig[a].begin(),sig[a].end());grouped[sig[a]].push_back(a);}
  std::vector<std::vector<int>> families;for(const auto& [s,ids]:grouped)families.push_back(ids);
  using Engine=tca_dlr::SummedRootSelfEnergyEvaluator;using Physical=tca_dlr::SummedPhysicalEvaluator;
  std::vector<typename Engine::OperatorTranspose> ops,bulla;
  for(int op=0;op<2*c.Nbath;++op)ops.push_back(Engine::MakeOperatorTranspose(c,op));
  for(int b=0;b<c.Nbath;++b){auto op=Physical::MakeBullaTranspose(c,b);
    // Verify the full dense measurement endpoint, including entries beyond
    // the sparse one-body support used by some higher-order backends.
    for(auto& e:op.entries){auto q=c.Njm_ind(e.source,2*b);if(q.njm_c==2)
      for(int i=0;i<e.matrix.rows();++i)for(int j=0;j<e.matrix.cols();++j)e.matrix(i,j)=q.transposed?(*q.njm)(j,i):(*q.njm)(i,j);}
    bulla.push_back(std::move(op));
  }
  bool anomalous=false;for(int b=0;b<c.Nbath;++b)anomalous=anomalous || c.bath_SC(b);
  if(!anomalous){
    // Keep the cheap, established normal-bath candidate. Exhaustive searches
    // over large f-shell product spaces would cost more than their saving.
    ++attempts;auto v=oca_symmetry::discover(c);
    oca_symmetry::verify_endpoints<Physical>(c,v,ops,bulla);return v;
  }
  int best_moved=0;std::vector<int> p(c.Npsi,-1);constexpr int budget=4096;
  std::function<void(int)> visit_family;
  visit_family=[&](int f){
    if(attempts>=budget)return;
    if(f==int(families.size())){
      bool moved=false;for(int a=0;a<c.Npsi;++a)moved=moved || p[a]!=a;if(!moved)return;
      ++attempts;auto v=oca_symmetry::candidate(c,p,true);if(!v.active)return;
      oca_symmetry::verify_endpoints<Physical>(c,v,ops,bulla);if(!v.active)return;
      int moved_baths=0;for(int b=0;b<c.Nbath;++b)moved_baths+=v.bath[b]!=b;
      if(moved_baths>best_moved){best_moved=moved_baths;best=std::move(v);}return;
    }
    const auto& ids=families[f];
    std::function<void()> pair;
    pair=[&](){if(attempts>=budget)return;int a=-1;for(int i:ids)if(p[i]<0){a=i;break;}
      if(a<0){visit_family(f+1);return;}
      p[a]=a;pair();p[a]=-1;
      for(int b:ids)if(b>a && p[b]<0){p[a]=b;p[b]=a;pair();p[a]=p[b]=-1;}
    };pair();
  };visit_family(0);
  if(best.active)best.reason="verified ordered signed spatial involution";
  else if(attempts>=budget)best.reason="spatial search budget exhausted";
  return best;
}

// Retain the time-reversal reconstruction while making its weighted row
// selector invariant under the spatial map. Conflicting signed components
// retain every row with weight one, rather than counting a symmetry twice.
template<class CixType>
time_reversal::Rows compatible_rows(const CixType& c,const oca_symmetry::Verification& tr,const oca_symmetry::Verification& spatial) {
  std::vector<int> sign(c.Natst+1),weight(c.Natst+1);
  for(int seed=1;seed<=c.Natst;++seed)if(!sign[seed]){
    std::vector<int> members{seed};sign[seed]=1;bool conflict=false;
    for(std::size_t k=0;k<members.size();++k){int i=members[k];
      for(auto [j,s]:{std::pair{tr.state[i],-sign[i]},std::pair{spatial.state[i],sign[i]}}){
        if(sign[j])conflict=conflict || sign[j]!=s;else{sign[j]=s;members.push_back(j);}}
    }
    for(int i:members)weight[i]=conflict?1:1+sign[i];
  }
  time_reversal::Rows rows;rows.index.resize(c.Nblock);rows.weight.resize(c.Nblock);
  for(int b=0;b<c.Nblock;++b)for(int i=0;i<c.Block[b].size();++i){int s=c.Block[b](i);++rows.full;
    if(weight[s]){rows.index[b].push_back(i);rows.weight[b].push_back(weight[s]);++rows.retained;}}
  return rows;
}

template<class CixType,class Workspace>
void configure(const CixType& c,Workspace& w,bool use_tr,bool enabled,std::ostream* log,const char* stage) {
  auto tr=time_reversal::configure(c,w,use_tr,log,stage);
  int attempts=0;oca_symmetry::Verification v;
  if(enabled)v=discover(c,w,attempts);
  if(v.active && tr.map.active)w.set_self_rows(compatible_rows(c,tr.map,v));
  w.set_spatial_symmetry(v);
  if(log){int reps=0;for(int b=0;b<c.Nbath;++b)reps+=!v.active || b<=v.bath[b];
    *log<<"Spatial symmetry: stage="<<stage<<" active="<<v.active<<" reason="<<v.reason<<" candidates="<<attempts
        <<" bath_representatives="<<reps<<"/"<<c.Nbath<<" self_rows="<<w.self_rows().retained<<"/"<<w.self_rows().full;
    if(v.active){*log<<" orbitals=";for(int a:v.orbital)*log<<a<<',';*log<<" phases=";for(int a:v.phase)*log<<a<<',';}
    *log<<std::endl;
  }
}
}
#endif
