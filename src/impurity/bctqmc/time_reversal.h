#ifndef EDMFT_TIME_REVERSAL_H
#define EDMFT_TIME_REVERSAL_H
#include "oca_symmetry.h"
#include "symmetry_rows.h"
#include "tca_dlr.h"
#include <functional>
#include <tuple>
namespace edmft::time_reversal {
struct Verification {
  oca_symmetry::Verification map;
  Rows rows;
  int candidates=0;
  double bath_defect=0,bath_bound=0;
};
template<class Physical,class CixType,class Op>
bool interaction_covariance(const CixType& c,oca_symmetry::Verification& v,
    const Op& a,const Op& b,int phase) {
  if constexpr(requires{Physical::InteractionOperator(a);})
    return oca_symmetry::covariant(c,v,Physical::InteractionOperator(a),
      Physical::InteractionOperator(b),phase,true,&v.endpoint_defect,&v.endpoint_roundoff_bound);
  return true;
}
// In this real product basis a Kramers antiunitary is represented by a signed
// orbital permutation followed by conjugation. The permutation is proposed
// from identical diagonal bath records; its meaning is NEVER inferred from
// orbital numbers, spin labels, j, cluster momenta, or the valence.
template<class CixType,class Workspace>
Verification candidate(const CixType& c,const Workspace& w,const std::vector<int>& permutation) {
  Verification result;auto& v=result.map;v.orbital=permutation;
  const auto reject=[&](const char* why){v.active=false;v.reason=why;return result;};
  int vacuum=0,maxN=0;
  for(int i=1;i<=c.Natst;++i){int n=c.iN(i-1);if(n<0||n>c.Npsi)return reject("non-product valence");if(n==0){if(vacuum)return reject("multiple vacua");vacuum=i;}maxN=std::max(maxN,n);}
  if(!vacuum)return reject("vacuum-connected product basis required");
  std::vector<std::uint64_t> occupation(c.Natst+1);
  v.state.assign(c.Natst+1,0);v.state_phase.assign(c.Natst+1,0);v.state[vacuum]=vacuum;v.state_phase[vacuum]=1;
  for(int n=1;n<=maxN;++n)for(int i=1;i<=c.Natst;++i)if(c.iN(i-1)==n){
    bool found=false;
    for(int a=0;a<c.Npsi;++a){int parent=c.Fpsi(2*a+1,i);if(!parent)continue;int p=std::abs(parent);
      if(p>c.Natst||c.iN(p-1)!=n-1||!v.state[p])return reject("non-product annihilation map");
      auto bit=std::uint64_t{1}<<a;if(occupation[p]&bit)return reject("non-fermionic occupation map");
      int original=c.Fpsi(2*a,p),mapped=c.Fpsi(2*permutation[a],v.state[p]);
      if(std::abs(original)!=i||!mapped||std::abs(mapped)>c.Natst)return reject("orbital map leaves atomic basis");
      int phase=v.state_phase[p]*(mapped<0?-1:1)*(original<0?-1:1);auto mask=occupation[p]|bit;
      if(found&&(v.state[i]!=std::abs(mapped)||v.state_phase[i]!=phase||occupation[i]!=mask))return reject("inconsistent fermionic paths");
      found=true;v.state[i]=std::abs(mapped);v.state_phase[i]=phase;occupation[i]=mask;
    }
    if(!found)return reject("disconnected product state");
  }
  for(int i=1;i<=c.Natst;++i){if(v.state[v.state[i]]!=i)return reject("state map not involutive");
    for(int a=0;a<c.Npsi;++a)for(int k=0;k<2;++k){int x=c.Fpsi(2*a+k,i),y=c.Fpsi(2*permutation[a]+k,v.state[i]);
      if(bool(x)!=bool(y))return reject("operator support differs");
      if(x&&(v.state[std::abs(x)]!=std::abs(y)||(x<0?-1:1)*v.state_phase[std::abs(x)]!=v.state_phase[i]*(y<0?-1:1)))return reject("operator signs differ");}}
  v.block.resize(c.Nblock);std::vector<int> seen(c.Nblock);
  for(int b=0;b<c.Nblock;++b){if(!c.Block[b].size())return reject("empty block");int target=c.iBlck(v.state[c.Block[b](0)])-1;
    if(target<0||target>=c.Nblock||c.Block[target].size()!=c.Block[b].size()||seen[target]++)return reject("block permutation differs");v.block[b]=target;
    for(int i=0;i<c.Block[b].size();++i)if(c.iBlck(v.state[c.Block[b](i)])-1!=target)return reject("permutation splits a block");}
  oca_symmetry::BinaryConstraints equations;
  for(int a=0;a<c.Npsi;++a)if(!equations.add((std::uint64_t{1}<<a)^(std::uint64_t{1}<<permutation[a]),true))return reject("Kramers phase constraint failed");
  const auto match=[&](auto x,auto y,int i,int j){if(x==0&&y==0)return true;int base=v.state_phase[i]*v.state_phase[j];
    if(x==base*y)return equations.add(occupation[i]^occupation[j],false);
    if(x==-base*y)return equations.add(occupation[i]^occupation[j],true);return false;};
  for(int b=0;b<c.Nblock;++b)for(int i=0;i<c.Block[b].size();++i)for(int j=0;j<c.Block[b].size();++j){
    int si=c.Block[b](i),sj=c.Block[b](j),tb=v.block[b],ti=c.inBlock(v.state[si]),tj=c.inBlock(v.state[sj]);
    double h=c.Ham[b](i,j),th=c.Ham[tb](ti,tj),u=c.SidecarHamiltonian(b)(i,j),tu=c.SidecarHamiltonian(tb)(ti,tj);
    if(!std::isfinite(h)||!std::isfinite(th)||!std::isfinite(u)||!std::isfinite(tu)||!match(h,th,si,sj)||!match(u,tu,si,sj))return reject("Hamiltonian covariance failed");
    if(!match(c.Sgm_symmetry[b](i,j),c.Sgm_symmetry[tb](ti,tj),si,sj))return reject("signed Sindex covariance failed");
    if(c.PseudoElementAllowed(b,i,j)!=c.PseudoElementAllowed(tb,ti,tj))return reject("sector projection differs");}
  // Compare the actual rational bath, not merely a few time samples. A
  // reversed anomalous endpoint order changes d and therefore E -> -E.
  // Normal residues remain unrestricted at opposite energies.
  using Key=std::tuple<int,int,double>;
  struct Coefficient{double value=0,scale=0;};std::map<Key,Coefficient> bath;
  for(int b=0;b<c.Nbath;++b)for(int d=0;d<2;++d)for(int l=0;l<w.pole_count();++l){int channel=c.bath_sigind0(b);
    double x=-(channel<0?-1.:1.)*w.residues()(std::abs(channel),l);
    if(!std::isfinite(x)||!std::isfinite(w.poles()(l)))return reject("nonfinite bath");
    auto& out=bath[{c.Fi_indx(2*b+1-d),c.Fi_indx(2*b+d),d ? -w.poles()(l):w.poles()(l)}];out.value+=x;out.scale+=std::abs(x);}
  for(const auto& [key,x]:bath){auto [a,b,e]=key;auto it=bath.find({2*permutation[a/2]+a%2,2*permutation[b/2]+b%2,e});
    if(it==bath.end()){if(x.value!=0.)return reject("missing bath pole partner");continue;}
    const auto& y=it->second;double bound=8*std::numeric_limits<double>::epsilon()*std::max(x.scale,y.scale);
    result.bath_bound=std::max(result.bath_bound,bound);
    bool same=std::abs(x.value-y.value)<=bound,opposite=std::abs(x.value+y.value)<=bound;
    result.bath_defect=std::max(result.bath_defect,std::min(std::abs(x.value-y.value),std::abs(x.value+y.value)));
    if(!same&&!opposite)return reject("rational bath covariance failed");
    if(same!=opposite&&!equations.add((std::uint64_t{1}<<(a/2))^(std::uint64_t{1}<<(b/2)),opposite))return reject("bath phase constraints failed");}
  auto solution=equations.solve();v.phase.resize(c.Npsi);
  for(int a=0;a<c.Npsi;++a)v.phase[a]=(solution&(std::uint64_t{1}<<a))?-1:1;
  for(int i=1;i<=c.Natst;++i)if(std::popcount(occupation[i]&solution)%2)v.state_phase[i]*=-1;
  for(int i=1;i<=c.Natst;++i)if(v.state_phase[i]*v.state_phase[v.state[i]]!=((c.iN(i-1)&1)?-1:1))return reject("many-body Kramers square failed");
  // Verify the actual operator matrices, including interaction endpoints.
  using Engine=tca_dlr::SummedRootSelfEnergyEvaluator;using Physical=tca_dlr::SummedPhysicalEvaluator;
  std::vector<Engine::OperatorTranspose> operators;std::map<int,int> lookup;
  for(int op=0;op<2*c.Nbath;++op){lookup[c.Fi_indx(op)]=op;operators.push_back(Engine::MakeOperatorTranspose(c,op));}
  for(int op=0;op<2*c.Nbath;++op){int raw=c.Fi_indx(op),target=2*permutation[raw/2]+raw%2;
    if(!lookup.count(target)||!oca_symmetry::covariant(c,v,operators[op],operators[lookup[target]],v.phase[raw/2]))return reject("operator matrix covariance failed");
    if(!interaction_covariance<Physical>(c,v,operators[op],operators[lookup[target]],v.phase[raw/2]))return reject("interaction endpoint covariance failed");}

  result.rows.index.resize(c.Nblock);result.rows.weight.resize(c.Nblock);
  for(int b=0;b<c.Nblock;++b)for(int i=0;i<c.Block[b].size();++i){int state=c.Block[b](i);++result.rows.full;
    if(state<=v.state[state]){result.rows.index[b].push_back(i);result.rows.weight[b].push_back(state==v.state[state]?1:2);++result.rows.retained;}}
  v.active=result.rows.active();v.reason=v.active?"verified Kramers antiunitary":"no row saving";return result;
}
template<class CixType,class Workspace>
Verification discover(const CixType& c,const Workspace& w) {
  Verification result;result.map.reason="unsupported real product basis";
  if(c.Npsi<2||c.Npsi>62||(c.Npsi&1))return result;
  // This first implementation accelerates anomalous baths. The existing
  // normal-bath task reduction remains the better normal-state schedule.
  bool sc=false;for(int b=0;b<c.Nbath;++b)sc=sc||c.bath_SC(b);
  if(!sc){result.map.reason="normal bath uses existing schedule";return result;}
  using Signature=std::vector<std::array<int,3>>;std::vector<Signature> sig(c.Npsi);
  for(int b=0;b<c.Nbath;++b){int a=c.Fi_indx(2*b),d=c.Fi_indx(2*b+1);
    if(a<0||d<0||a/2>=c.Npsi||d/2>=c.Npsi){result.map.reason="invalid bath operator";return result;}
    if(!c.bath_SC(b)&&!(a&1)&&(d&1)&&a/2==d/2)sig[a/2].push_back({c.bath_sigind0(b),int(c.bath_conj(b)),int(c.SC_bath(b))});}
  std::map<Signature,std::vector<int>> families;
  for(int a=0;a<c.Npsi;++a){if(sig[a].empty()){result.map.reason="missing normal diagonal bath";return result;}std::sort(sig[a].begin(),sig[a].end());families[sig[a]].push_back(a);}
  std::vector<std::vector<int>> family;for(const auto& [signature,ids]:families){if(ids.size()%2){result.map.reason="no Kramers pairing of diagonal baths";return result;}family.push_back(ids);}
  std::vector<int> permutation(c.Npsi,-1);int trials=0;bool stop=false;
  std::function<void(std::size_t)> next_family;
  std::function<void(std::size_t,const std::vector<int>&)> pair;
  pair=[&](std::size_t k,const std::vector<int>& remaining){if(stop)return;if(remaining.empty()){next_family(k+1);return;}
    int a=remaining.front();for(std::size_t j=1;j<remaining.size()&&!stop;++j){int b=remaining[j];permutation[a]=b;permutation[b]=a;
      std::vector<int> rest;for(std::size_t i=1;i<remaining.size();++i)if(i!=j)rest.push_back(remaining[i]);pair(k,rest);}};
  next_family=[&](std::size_t k){if(stop)return;if(k==family.size()){if(trials>=4096){stop=true;return;}++trials;
      result=candidate(c,w,permutation);if(result.map.active)stop=true;return;}pair(k,family[k]);};
  next_family(0);result.candidates=trials;if(!result.map.active&&trials>=4096)result.map.reason="candidate search limit; using full calculation";return result;
}
template<class CixType,class Workspace>
Verification configure(const CixType& c,Workspace& w,bool enabled,std::ostream* log,const char* stage) {
  Verification result;result.map.reason="disabled";w.set_self_rows({});
  if(enabled)result=discover(c,w);
  if(result.map.active)w.set_self_rows(result.rows);
  if(log){*log<<"Time reversal self rows: stage="<<stage<<" active="<<result.map.active<<" reason="<<result.map.reason
    <<" candidates="<<result.candidates<<" rows="<<result.rows.retained<<"/"<<result.rows.full
    <<" bath_defect="<<result.bath_defect<<" bath_roundoff_bound="<<result.bath_bound;
    if(result.map.active){*log<<" orbitals=[";for(int a:result.map.orbital)*log<<a<<",";*log<<"] phases=[";for(int a:result.map.phase)*log<<a<<",";*log<<"]";}*log<<std::endl;}
  return result;
}
}
#endif
