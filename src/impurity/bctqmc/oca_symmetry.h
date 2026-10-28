#ifndef EDMFT_OCA_SYMMETRY_H
#define EDMFT_OCA_SYMMETRY_H
#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <map>
#include <limits>
#include <numeric>
#include <string>
#include <vector>

namespace edmft::oca_symmetry {

// A verified signed permutation, not an assumption about spin or j. Labels
// must be invariant themselves: no compact output is silently relabelled.
struct Verification {
  bool active=false;
  std::string reason="disabled";
  std::vector<int> orbital,phase,state,state_phase,block,bath;
  double endpoint_defect=0.,endpoint_roundoff_bound=0.;
};

class BinaryConstraints {
  std::array<std::uint64_t,64> row{};
  std::array<bool,64> value{};
public:
  bool add(std::uint64_t bits,bool rhs) {
    for(int p=63;p>=0;--p)if(bits&(std::uint64_t{1}<<p)) {
      if(row[p]){bits^=row[p];rhs=rhs!=value[p];}
      else{row[p]=bits;value[p]=rhs;return true;}
    }
    return !rhs;
  }
  std::uint64_t solve() const {
    std::uint64_t result=0;
    for(int p=0;p<64;++p)if(row[p] &&
      (value[p] != bool(std::popcount(row[p]&result)%2)))
      result|=std::uint64_t{1}<<p;
    return result;
  }
};

template<class CixType>
Verification candidate(const CixType& c,const std::vector<int>& permutation,bool unitary=false) {
  Verification v;v.orbital=permutation;
  const auto reject=[&](const char* why){v.reason=why;return v;};
  if(c.Npsi<1 || c.Npsi>63 || int(permutation.size())!=c.Npsi)return reject("unsupported orbital count");
  for(int a=0;a<c.Npsi;++a)if(permutation[a]<0 || permutation[a]>=c.Npsi || permutation[permutation[a]]!=a)
    return reject("orbital map is not involutive");
  int vacuum=0,maxN=0;
  for(int i=1;i<=c.Natst;++i){
    const int n=c.iN(i-1);if(n<0 || n>c.Npsi)return reject("non-product valence");
    if(n==0){if(vacuum)return reject("multiple vacuum states");vacuum=i;}maxN=std::max(maxN,n);
  }
  if(!vacuum)return reject("vacuum-connected basis required");
  std::vector<std::uint64_t> occupation(c.Natst+1);
  v.state.assign(c.Natst+1,0);v.state_phase.assign(c.Natst+1,0);
  v.state[vacuum]=vacuum;v.state_phase[vacuum]=1;
  for(int n=1;n<=maxN;++n)for(int i=1;i<=c.Natst;++i)if(c.iN(i-1)==n){
    bool found=false;
    for(int a=0;a<c.Npsi;++a){
      const int parent=c.Fpsi(2*a+1,i);if(!parent)continue;
      const int p=std::abs(parent);
      if(p>c.Natst || c.iN(p-1)!=n-1 || !v.state[p])return reject("non-product annihilation map");
      const auto bit=std::uint64_t{1}<<a;
      if(occupation[p]&bit)return reject("non-fermionic occupation map");
      const int original=c.Fpsi(2*a,p),mapped=c.Fpsi(2*v.orbital[a],v.state[p]);
      if(std::abs(original)!=i || !mapped || std::abs(mapped)>c.Natst)return reject("permutation leaves atomic basis");
      const int sign=v.state_phase[p]*(mapped<0?-1:1)*(original<0?-1:1);
      const auto mask=occupation[p]|bit;
      if(found && (v.state[i]!=std::abs(mapped) || v.state_phase[i]!=sign || occupation[i]!=mask))return reject("inconsistent product-state paths");
      found=true;v.state[i]=std::abs(mapped);v.state_phase[i]=sign;occupation[i]=mask;
    }
    if(!found)return reject("disconnected product state");
  }
  for(int i=1;i<=c.Natst;++i){
    if(v.state[v.state[i]]!=i || c.iN(v.state[i]-1)!=c.iN(i-1))return reject("state permutation is not involutive");
    for(int a=0;a<c.Npsi;++a)for(int kind=0;kind<2;++kind){
      const int x=c.Fpsi(2*a+kind,i),y=c.Fpsi(2*v.orbital[a]+kind,v.state[i]);
      if(bool(x)!=bool(y))return reject("operator support differs");
      if(x && (v.state[std::abs(x)]!=std::abs(y) ||
        (x<0?-1:1)*v.state_phase[std::abs(x)]!=v.state_phase[i]*(y<0?-1:1)))return reject("operator signs differ");
    }
  }
  // An entire container must map to one container. This guarantees the same
  // block routing for all intermediate products and composite endpoints.
  v.block.resize(c.Nblock);std::vector<int> block_seen(c.Nblock,0);
  for(int b=0;b<c.Nblock;++b){
    if(!c.Block[b].size())return reject("empty atomic block");
    int target=c.iBlck(v.state[c.Block[b](0)])-1;
    if(target<0 || target>=c.Nblock || c.Block[target].size()!=c.Block[b].size() || block_seen[target]++)return reject("block permutation differs");
    v.block[b]=target;
    for(int i=0;i<c.Block[b].size();++i)if(c.iBlck(v.state[c.Block[b](i)])-1!=target)return reject("permutation splits an atomic block");
  }
  BinaryConstraints constraints;
  if(unitary)for(int a=0;a<c.Npsi;++a)
    if(!constraints.add((std::uint64_t{1}<<a)^(std::uint64_t{1}<<v.orbital[a]),false))
      return reject("unitary square constraint failed");
  const auto match=[&](auto x,auto y,int i,int j){
    if(x==0 && y==0)return true;
    const int base=v.state_phase[i]*v.state_phase[j];
    if(x==base*y)return constraints.add(occupation[i]^occupation[j],false);
    if(x==-base*y)return constraints.add(occupation[i]^occupation[j],true);
    return false;
  };
  for(int b=0;b<c.Nblock;++b){
    const int target=v.block[b];
    for(int a=0;a<c.Block[b].size();++a)for(int d=0;d<c.Block[b].size();++d){
      int i=c.Block[b](a),j=c.Block[b](d),ti=c.inBlock(v.state[i]),tj=c.inBlock(v.state[j]);
      double h=c.Ham[b](a,d),th=c.Ham[target](ti,tj);
      double u=c.SidecarHamiltonian(b)(a,d),tu=c.SidecarHamiltonian(target)(ti,tj);
      if(!std::isfinite(h)||!std::isfinite(th)||!std::isfinite(u)||!std::isfinite(tu)||
         !match(h,th,i,j)||!match(u,tu,i,j))return reject("Hamiltonian covariance failed");
      if(!match(c.Sgm_symmetry[b](a,d),c.Sgm_symmetry[target](ti,tj),i,j))return reject("signed Sindex covariance failed");
      if(c.PseudoElementAllowed(b,a,d)!=c.PseudoElementAllowed(target,ti,tj))return reject("sector projection differs");
    }
  }
  // A multiset preserves duplicate bath observations and signed channels.
  using Key=std::array<int,6>;
  std::map<Key,std::array<int,2>> baths;
  for(int b=0;b<c.Nbath;++b){int s=c.bath_sigind0(b);++baths[{c.Fi_indx(2*b),c.Fi_indx(2*b+1),std::abs(s),int(c.bath_conj(b)),int(c.SC_bath(b)),int(c.bath_SC(b))}][s<0];}
  for(const auto& [key,count]:baths){
    auto mapped=key;mapped[0]=2*v.orbital[key[0]/2]+key[0]%2;mapped[1]=2*v.orbital[key[1]/2]+key[1]%2;
    auto it=baths.find(mapped);if(it==baths.end())return reject("bath permutation differs");
    const bool same=count==it->second,opposite=count[0]==it->second[1]&&count[1]==it->second[0];
    if(!same&&!opposite)return reject("bath multiplicity differs");
    if(same!=opposite && !constraints.add((std::uint64_t{1}<<(key[0]/2))^(std::uint64_t{1}<<(key[1]/2)),opposite))return reject("bath phase constraints failed");
  }
  const auto solution=constraints.solve();v.phase.resize(c.Npsi);
  for(int a=0;a<c.Npsi;++a)v.phase[a]=(solution&(std::uint64_t{1}<<a))?-1:1;
  for(int i=1;i<=c.Natst;++i)if(std::popcount(occupation[i]&solution)%2)v.state_phase[i]*=-1;
  v.active=true;v.reason="verified";return v;
}

template<class CixType>
Verification discover(const CixType& c) {
  Verification v;
  const auto reject=[&](const char* why){v.reason=why;return v;};
  // The present automatic candidate finder covers real, normal product
  // bases with unambiguous diagonal bath partners. Other bases use the full
  // schedule. In particular, anomalous conjugation is never guessed here.
  if(c.Npsi<1 || c.Npsi>63)return reject("unsupported orbital count");
  using Diagonal=std::vector<std::array<int,3>>;
  std::vector<Diagonal> diagonal(c.Npsi);
  for(int b=0;b<c.Nbath;++b){
    const int a=c.Fi_indx(2*b),d=c.Fi_indx(2*b+1);
    if(c.bath_SC(b) || (a&1) || !(d&1))return reject("normal bath required");
    if(a/2>=c.Npsi || d/2>=c.Npsi)return reject("invalid bath operator");
    if(a/2==d/2)diagonal[a/2].push_back({c.bath_sigind0(b),int(c.bath_conj(b)),int(c.SC_bath(b))});
  }
  std::map<Diagonal,std::vector<int>> families;
  for(int a=0;a<c.Npsi;++a){
    if(diagonal[a].empty())return reject("missing diagonal bath");
    std::sort(diagonal[a].begin(),diagonal[a].end());families[diagonal[a]].push_back(a);
  }
  v.orbital.resize(c.Npsi);std::iota(v.orbital.begin(),v.orbital.end(),0);
  bool moved=false;
  for(const auto& [signature,ids]:families){
    if(ids.size()>2)return reject("ambiguous diagonal bath partners");
    if(ids.size()==2){v.orbital[ids[0]]=ids[1];v.orbital[ids[1]]=ids[0];moved=true;}
  }
  if(!moved)return reject("no diagonal bath partners");
  return candidate(c,v.orbital);
}


template<class CixType,class Op>
bool covariant(const CixType& c,const Verification& v,const Op& a,const Op& b,int phase,
 bool composite=false,double* defect=nullptr,double* bound=nullptr) {
  // H, labels, ordinary operators and bath covariance are exact above.
  // Composite operators are computed by floating-point sums/subtractions.
  // Permit only their arithmetic roundoff, never a user physics cutoff.
  double scale=1.;
  if(composite){for(const auto& x:a.entries)if(x.matrix.size())scale=std::max(scale,x.matrix.cwiseAbs().maxCoeff());for(const auto& x:b.entries)if(x.matrix.size())scale=std::max(scale,x.matrix.cwiseAbs().maxCoeff());}
  const double tolerance=composite?8.*std::numeric_limits<double>::epsilon()*scale:0.;
  if(bound)*bound=std::max(*bound,tolerance);
  for(int source=0;source<c.Nblock;++source){
    int ai=a.by_source[source],bi=b.by_source[v.block[source]];
    if(ai<0 || bi<0){
      if(ai>=0 && (!a.entries[ai].matrix.allFinite() || (a.entries[ai].matrix.array().abs()>tolerance).any()))return false;
      if(bi>=0 && (!b.entries[bi].matrix.allFinite() || (b.entries[bi].matrix.array().abs()>tolerance).any()))return false;
      continue;
    }
    const auto& x=a.entries[ai];const auto& y=b.entries[bi];
    if(v.block[x.target]!=y.target || x.matrix.rows()!=y.matrix.rows() || x.matrix.cols()!=y.matrix.cols())return false;
    for(int i=0;i<x.matrix.rows();++i)for(int j=0;j<x.matrix.cols();++j){
      int si=c.Block[source](i),sj=c.Block[x.target](j);
      double p=x.matrix(i,j),q=y.matrix(c.inBlock(v.state[si]),c.inBlock(v.state[sj]));
      const double error=std::abs(p-phase*v.state_phase[si]*v.state_phase[sj]*q);
      if(defect)*defect=std::max(*defect,error);
      if(!std::isfinite(p)||!std::isfinite(q)||error>tolerance)return false;
    }
  }
  return true;
}

// Explicit endpoint checks also protect against future custom composite
// operators that may not be derived from the local Hamiltonian.
template<class Physical,class CixType,class Op>
void verify_endpoints(const CixType& c,Verification& v,
 const std::vector<Op>& operators,const std::vector<Op>& bulla) {
  if(!v.active)return;
  using Key=std::array<int,7>;std::map<Key,std::vector<int>> lookup;
  for(int b=0;b<c.Nbath;++b){int s=c.bath_sigind0(b);lookup[{c.Fi_indx(2*b),c.Fi_indx(2*b+1),std::abs(s),s<0,int(c.bath_conj(b)),int(c.SC_bath(b)),int(c.bath_SC(b))}].push_back(b);}
  std::map<Key,std::size_t> used;v.bath.assign(c.Nbath,-1);
  for(int b=0;b<c.Nbath;++b){
    int a=c.Fi_indx(2*b)/2,d=c.Fi_indx(2*b+1)/2,s=c.bath_sigind0(b);
    Key key{2*v.orbital[a]+c.Fi_indx(2*b)%2,2*v.orbital[d]+c.Fi_indx(2*b+1)%2,std::abs(s),(s<0)!=(v.phase[a]*v.phase[d]<0),int(c.bath_conj(b)),int(c.SC_bath(b)),int(c.bath_SC(b))};
    auto found=lookup.find(key);auto offset=used[key]++;
    if(found==lookup.end() || offset>=found->second.size()){v.active=false;v.reason="bath endpoint pairing failed";return;}
    int target=found->second[offset];v.bath[b]=target;
    if(!covariant(c,v,operators[2*b],operators[2*target],v.phase[a]) ||
       !covariant(c,v,operators[2*b+1],operators[2*target+1],v.phase[d]) ||
       !covariant(c,v,bulla[b],bulla[target],v.phase[a],true,&v.endpoint_defect,&v.endpoint_roundoff_bound)) {
      v.active=false;v.reason="ordinary/Bulla covariance failed";return;
    }
    if constexpr(requires{c.kugler_estimators;Physical::InteractionOperator(operators[0]);}) {
      if(c.kugler_estimators &&
        (!covariant(c,v,Physical::InteractionOperator(operators[2*b]),Physical::InteractionOperator(operators[2*target]),v.phase[a],true,&v.endpoint_defect,&v.endpoint_roundoff_bound) ||
         !covariant(c,v,Physical::InteractionOperator(operators[2*b+1]),Physical::InteractionOperator(operators[2*target+1]),v.phase[d],true,&v.endpoint_defect,&v.endpoint_roundoff_bound))) {
        v.active=false;v.reason="Kugler covariance failed";return;
      }
    }
  }
  for(int b=0;b<c.Nbath;++b)if(v.bath[v.bath[b]]!=b){v.active=false;v.reason="bath map is not involutive";return;}
}
}
#endif
