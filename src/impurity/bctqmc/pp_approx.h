#ifndef BCTQMC_PP_APPROX_H
#define BCTQMC_PP_APPROX_H

// A reference partition for valence-selective pseudo-particle projection.
// Labels are derived from operator action, never from Jz or other CIX labels.
#include <cstdint>
#include <numeric>
#include <set>
#include "json.hpp"

namespace pp_approx {
inline std::vector<int> selection(const nlohmann::json& value)
{
  if (!value.is_array())
    throw std::invalid_argument("pp_approx_N must be an array of nonnegative integer valences");
  std::set<int> result;
  for (const auto& n : value) {
    if (!n.is_number_integer() || n.get<long long>() < 0 ||
        n.get<long long>() > std::numeric_limits<int>::max())
      throw std::invalid_argument("pp_approx_N must contain nonnegative integer valences");
    result.insert(n.get<int>());
  }
  return {result.begin(), result.end()};
}

inline std::string selector_text(const std::vector<int>& values)
{
  return nlohmann::json(values).dump();
}

struct Partition {
  std::vector<int> parent;
  explicit Partition(int n) : parent(n + 1) { std::iota(parent.begin(), parent.end(), 0); }
  int root(int i) { return parent[i] == i ? i : parent[i] = root(parent[i]); }
  bool join(int i, int j) {
    i = root(i); j = root(j);
    if (i == j) return false;
    parent[std::max(i,j)] = std::min(i,j);
    return true;
  }
};

// Correct the creation-operator interaction endpoint by [delta H,c^dagger].
// Njm is stored with source rows and target columns, as in Read_cix.
inline void correct_endpoints(Cix& cix)
{
  for (int ib=0; ib<cix.Nblock; ++ib) {
    for (int a=0; a<cix.Npsi; ++a) {
      const int target = cix.Fbpsi(2*a,ib);
      if (target < 0) continue;
      if (!cix.pp_unprojected_ham[ib].size() &&
          !cix.pp_unprojected_ham[target].size()) continue;
      const int ns=cix.Block[ib].size(), nt=cix.Block[target].size();
      Array<double,2> matrix(ns,nt); matrix=0.0;
      if (cix.Njm_c(ib,a)==2) matrix=cix.Njm[ib][a];
      else if (cix.Njm_c(ib,a)==1) {
        const double r=cix.Njm_r.at(cix.encode_block_bath(ib,a));
        for (int s=0; s<ns; ++s) {
          const int dest=cix.Fpsi(2*a,cix.Block[ib](s));
          if (dest) matrix(s,cix.inBlock(std::abs(dest)))=r*sign(dest);
        }
      }
      for (int s=0; s<ns; ++s) {
        const int dest=cix.Fpsi(2*a,cix.Block[ib](s));
        if (dest && cix.pp_unprojected_ham[target].size()) {
          const int t=cix.inBlock(std::abs(dest));
          for (int k=0; k<nt; ++k)
            matrix(s,k)+=sign(dest)*(cix.Ham[target](t,k)-cix.pp_unprojected_ham[target](t,k));
        }
        if (cix.pp_unprojected_ham[ib].size()) {
          for (int l=0; l<ns; ++l) {
            const int to=cix.Fpsi(2*a,cix.Block[ib](l));
            if (to) matrix(s,cix.inBlock(std::abs(to)))-=
              sign(to)*(cix.Ham[ib](s,l)-cix.pp_unprojected_ham[ib](s,l));
          }
        }
      }
      cix.Njm_c(ib,a)=2;
      cix.Njm[ib][a].resize(ns,nt);
      cix.Njm[ib][a]=matrix;
    }
  }
}

inline void configure(Cix& cix, const nlohmann::json& selector, double beta,
                      std::ostream& out)
{
  cix.pp_approx_N=selection(selector);
  if (cix.pp_approx_N.empty()) return;
  constexpr double tolerance=1e-9; // CIX Hamiltonians/levels are decimal text in eV.
  for (int n : cix.pp_approx_N) {
    bool present=false;
    for (int s=0; s<cix.Natst; ++s) present=present || cix.iN(s)==n;
    if (!present) throw std::invalid_argument("pp_approx_N selects an absent valence: "+std::to_string(n));
  }

  // Reconstruct off-diagonal h_ab c^dagger_a c_b from existing CIX bath
  // operators and impurity levels. The reference retains the full interaction.
  std::map<std::pair<int,int>,double> levels;
  for (int b=0; b<cix.Nbath; ++b) {
    if (cix.bath_SC(b))
      throw std::invalid_argument("pp_approx_N diagonal reference currently requires normal baths; an SC reference policy is not implemented");
    const int left=cix.Fi_indx(2*b), right=cix.Fi_indx(2*b+1);
    if (left/2==right/2) continue;
    const bool creation_first=(left%2)==0;
    const auto key=creation_first ? std::make_pair(left/2,right/2) : std::make_pair(right/2,left/2);
    const double value=(creation_first ? 1.0 : -1.0)*
      (cix.bath_sigind0(b)<0 ? -1.0 : 1.0)*cix.Ebath(std::abs(cix.bath_sigind0(b)));
    const auto [it,inserted]=levels.emplace(key,value);
    if (!inserted && std::abs(it->second-value)>tolerance)
      throw std::invalid_argument("conflicting CIX impurity levels for one operator pair");
  }
  for (const auto& [key,value] : levels) {
    auto reverse=levels.find({key.second,key.first});
    const double other=reverse==levels.end() ? 0.0 : reverse->second;
    if (!std::isfinite(value) || std::abs(value-other)>tolerance)
      throw std::invalid_argument("pp_approx_N requires a real symmetric one-body Hamiltonian");
  }
  std::vector<Array<double,2>> reference(cix.Nblock);
  for (int ib=0; ib<cix.Nblock; ++ib) reference[ib].reference(cix.Ham[ib].copy());
  for (const auto& [ab,value] : levels) {
    for (int s=1; s<=cix.Natst; ++s) {
      const int middle=cix.Fpsi(2*ab.second+1,s);
      if (!middle) continue;
      const int dest=cix.Fpsi(2*ab.first,std::abs(middle));
      if (!dest) continue;
      const int ib=cix.iBlck(s)-1;
      if (cix.iBlck(std::abs(dest))!=ib+1)
        throw std::invalid_argument("one-body operator leaves its CIX container block");
      reference[ib](cix.inBlock(std::abs(dest)),cix.inBlock(s))-=value*sign(middle)*sign(dest);
    }
  }

  Partition partition(cix.Natst);
  for (int ib=0; ib<cix.Nblock; ++ib) {
    for (int i=0; i<cix.Block[ib].size(); ++i) {
      for (int j=0; j<i; ++j) {
        // Remove roundoff left by subtracting separately rounded one-body terms.
        if (std::abs(reference[ib](i,j))<=tolerance && std::abs(reference[ib](j,i))<=tolerance)
          reference[ib](i,j)=reference[ib](j,i)=0.0;
        else partition.join(cix.Block[ib](i),cix.Block[ib](j));
      }
    }
  }
  // Diagonal hybridization transfers a coherence through the SAME physical
  // operator at both ends. Close this relation to a fixed point in all sectors.
  bool changed=true;
  int passes=0;
  while (changed) {
    changed=false; ++passes;
    for (int op=0; op<2*cix.Npsi; ++op) {
      std::vector<int> image(cix.Natst+1,0);
      for (int s=1; s<=cix.Natst; ++s) {
        const int dest=std::abs(cix.Fpsi(op,s));
        if (!dest) continue;
        const int root=partition.root(s);
        if (!image[root]) image[root]=dest;
        else {
          if (cix.iBlck(image[root])!=cix.iBlck(dest))
            throw std::invalid_argument("reference operator reaches multiple CIX containers");
          changed=partition.join(image[root],dest) || changed;
        }
      }
    }
  }
  cix.pp_reference_group.resize(cix.Natst+1);
  cix.pp_selected.resize(cix.Natst+1,false);
  cix.pp_unprojected_ham.resize(cix.Nblock);
  for (int s=1; s<=cix.Natst; ++s) {
    cix.pp_reference_group[s]=partition.root(s);
    cix.pp_selected[s]=std::binary_search(cix.pp_approx_N.begin(),cix.pp_approx_N.end(),cix.iN(s-1));
  }
  std::uint64_t hash=UINT64_C(14695981039346656037);
  auto feed=[&](int v) {
    for (int byte=0; byte<4; ++byte) {
      hash^=(static_cast<std::uint32_t>(v)>>(8*byte))&255;
      hash*=UINT64_C(1099511628211);
    }
  };
  feed(1); // Reference algorithm version, part of frozen-Sindex compatibility.
  for (int n : cix.pp_approx_N) feed(n);
  feed(-1); feed(cix.Nblock); feed(cix.Npsi); feed(cix.Nbath);
  for (int op=0; op<2*cix.Nbath; ++op) feed(cix.Fi_indx(op));
  for (int s=1; s<=cix.Natst; ++s) {
    feed(cix.iBlck(s)); feed(cix.inBlock(s)); feed(cix.iN(s-1));
    feed(cix.pp_reference_group[s]);
    for (int op=0; op<2*cix.Npsi; ++op) feed(cix.Fpsi(op,s));
  }
  std::ostringstream signature;
  signature<<std::hex<<hash;
  cix.pp_approx_signature=signature.str();
  out<<"Pseudo-particle approximation: pp_approx_N="<<selector_text(cix.pp_approx_N)
     <<" reference=diagonal_normal tolerance="<<tolerance<<" closure_passes="<<passes<<'\n';
  for (int ib=0; ib<cix.Nblock; ++ib) {
    bool selected=false;
    for (int i=0; i<cix.Block[ib].size(); ++i) selected=selected || cix.pp_selected[cix.Block[ib](i)];
    if (!selected) continue;
    cix.pp_unprojected_ham[ib].reference(cix.Ham[ib].copy());
    std::map<int,int> sizes;
    for (int i=0; i<cix.Block[ib].size(); ++i) ++sizes[cix.pp_reference_group[cix.Block[ib](i)]];
    int largest=0; long long allowed=0;
    for (const auto& [group,n] : sizes) { largest=std::max(largest,n); allowed+=static_cast<long long>(n)*n; }
    out<<"  CIX block "<<ib+1<<": dimension="<<cix.Block[ib].size()
       <<" reference_groups="<<sizes.size()<<" largest="<<largest<<" allowed_elements="<<allowed<<'\n';
    for (int i=0; i<cix.Block[ib].size(); ++i)
      for (int j=0; j<cix.Block[ib].size(); ++j)
        if (cix.pp_selected[cix.Block[ib](i)] || cix.pp_selected[cix.Block[ib](j)])
          cix.Ham[ib](i,j)=reference[ib](i,j);
  }
  correct_endpoints(cix);
  // Refresh the atomic energy origin and free energy for the projected H.
  std::vector<double> spectrum;
  for (int ib=0; ib<cix.Nblock; ++ib) {
    Array<double,2> h(cix.Ham[ib].copy());
    Array<double,1> e(h.extent(0)), work(3*h.extent(0));
    if (xsyev(h.extent(0),h.data(),h.extent(1),e.data(),work.data(),work.size(),"N"))
      throw std::runtime_error("pp_approx_N Hamiltonian diagonalization failed");
    for (int i=0; i<e.size(); ++i) spectrum.push_back(e(i));
  }
  cix.Emin=*std::min_element(spectrum.begin(),spectrum.end());
  cix.Emax=*std::max_element(spectrum.begin(),spectrum.end());
  double partition_sum=0.0;
  for (double e : spectrum) partition_sum+=std::exp(-beta*(e-cix.Emin));
  cix.F_atom=cix.Emin-std::log(partition_sum)/beta;
}

inline void check_restart_policy(const Cix& cix, const std::string& sindex,
                                  const std::string& restart)
{
  std::ifstream input(sindex);
  if (!input) {
    if (!cix.pp_approx_N.empty() && std::ifstream(restart).good())
      throw std::runtime_error("pp_approx_N requires a fresh restart when Sindex.dat is missing; archive "+restart+" first");
    return;
  }
  std::string line, policy="[]", signature;
  while (std::getline(input,line)) {
    if (line.rfind("# pp_approx_N ",0)==0) policy=line.substr(14);
    if (line.rfind("# pp_approx_signature ",0)==0) signature=line.substr(22);
  }
  if (policy!=selector_text(cix.pp_approx_N) ||
      (!cix.pp_approx_N.empty() && signature!=cix.pp_approx_signature))
    throw std::runtime_error("Sindex.dat belongs to a different pp_approx_N/reference partition; archive Sindex.dat and "+restart+" before changing the approximation");
}
} // namespace pp_approx
#endif
